import argparse
import hashlib
import fcntl
import json
import math
import os
import sys
import tempfile
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path

from .core import (CURVES, PROFILES, MAX_SKIP_DRAWS, MAX_KEYS_PER_SEED, attempt_at,
                   candidate, parse_time, stream_at, timestamp_count, total_candidates, validate_profile)
from .rngs import RNGS
from .targets import targets, valid_fingerprint

SCHEMA = 2


def atomic_json(path: Path, data: dict) -> None:
    descriptor, name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "w") as f:
            json.dump(data, f, indent=2, sort_keys=True)
            f.write("\n")
            f.flush()
            os.fsync(f.fileno())
        os.replace(temporary, path)
        directory = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        temporary.unlink(missing_ok=True)


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Offline public-key candidates for explicit timestamp/RNG/derivation hypotheses. Never outputs private keys.")
    sub = p.add_subparsers(dest="command", required=True)
    for name in ("estimate", "generate"):
        s = sub.add_parser(name)
        s.add_argument("--start", default="2000-01-01T00:00:00Z", help="inclusive integer epoch units or timezone-explicit ISO8601")
        s.add_argument("--end", default=datetime.now(timezone.utc).replace(microsecond=0).isoformat(), help="exclusive endpoint; default current UTC second")
        s.add_argument("--unit", choices=("seconds", "milliseconds"), default="seconds")
        s.add_argument("--step", type=int, default=1, help="positive number of selected units")
        s.add_argument("--rng", dest="rngs", choices=RNGS, nargs="+", default=["lrand48"])
        s.add_argument("--profile", dest="profiles", choices=PROFILES, nargs="+", default=["concat-be-reject"])
        s.add_argument("--curve", dest="curves", choices=CURVES, nargs="+", default=["p256"])
        s.add_argument("--match", action="append", help="local fingerprints, PEM public keys/certificates, or DER key/certificate; repeat to combine files")
        s.add_argument("--skip-draws", type=int, nargs="+", default=[0], help=f"initial RNG word counts to try, each 0..{MAX_SKIP_DRAWS}")
        s.add_argument("--keys-per-seed", type=int, default=1, help=f"successive candidate blocks per seed/profile/curve/skip, 1..{MAX_KEYS_PER_SEED}")
        if name == "generate":
            s.add_argument("--max-candidates", type=int, default=10000, help="attempt cap per invocation (not number of matches)")
            s.add_argument("--output", help="JSONL destination; default stdout; existing files are refused")
            s.add_argument("--checkpoint", help="atomic checkpoint path; requires --output")
            s.add_argument("--checkpoint-every", type=int, default=100, help="attempts per durable checkpoint")
            s.add_argument("--progress-every", type=float, default=5, help="seconds between JSON progress reports on stderr; 0 disables")
    r = sub.add_parser("resume")
    r.add_argument("checkpoint")
    r.add_argument("--max-candidates", type=int, default=10000)
    r.add_argument("--checkpoint-every", type=int, default=100)
    r.add_argument("--progress-every", type=float, default=5, help="seconds between JSON progress reports on stderr; 0 disables")
    return p


def make_config(args: argparse.Namespace) -> dict:
    config = {"start": parse_time(args.start, args.unit), "end": parse_time(args.end, args.unit),
              "unit": args.unit, "step": args.step, "rngs": list(dict.fromkeys(args.rngs)),
              "profiles": list(dict.fromkeys(args.profiles)), "curves": list(dict.fromkeys(args.curves)),
              "skip_draws": list(dict.fromkeys(args.skip_draws)), "keys_per_seed": args.keys_per_seed}
    if config["start"] < 0 or config["end"] <= config["start"] or config["step"] <= 0:
        raise ValueError("Require 0 <= start < end and step > 0")
    for rng in config["rngs"]:
        for profile in config["profiles"]:
            validate_profile(rng, profile)
    validate_consumption(config)
    selected = targets(args.match)
    config["match_fingerprints"] = sorted(selected) if selected is not None else None
    return config


def validate_consumption(config: dict) -> None:
    skips = config.get("skip_draws", [0])
    keys = config.get("keys_per_seed", 1)
    if (not isinstance(skips, list) or not skips or
            not all(type(n) is int and 0 <= n <= MAX_SKIP_DRAWS for n in skips)):
        raise ValueError(f"Skipped draws must be between 0 and {MAX_SKIP_DRAWS}")
    if len(set(skips)) != len(skips):
        raise ValueError("Duplicate skipped-draw counts")
    if type(keys) is not int or not 1 <= keys <= MAX_KEYS_PER_SEED:
        raise ValueError(f"Keys per seed must be between 1 and {MAX_KEYS_PER_SEED}")


def validate_checkpoint(data: dict) -> None:
    if not isinstance(data, dict) or type(data.get("schema_version")) is not int or data["schema_version"] not in (1, SCHEMA):
        raise ValueError("Unsupported checkpoint schema")
    required = {"config", "next_index", "emitted", "output", "output_offset", "output_sha256"}
    if not required <= data.keys() or not isinstance(data["config"], dict):
        raise ValueError("Incomplete checkpoint")
    c = data["config"]
    if data["schema_version"] == SCHEMA and not {"skip_draws", "keys_per_seed"} <= c.keys():
        raise ValueError("Missing checkpoint consumption settings")
    if data["schema_version"] == 1 and (c.get("skip_draws", [0]) != [0] or c.get("keys_per_seed", 1) != 1):
        raise ValueError("Schema 1 cannot contain extended consumption settings")
    validate_consumption(c)
    if c.get("unit") not in ("seconds", "milliseconds") or not all(type(c.get(k)) is int for k in ("start", "end", "step")):
        raise ValueError("Invalid checkpoint range")
    if c["start"] < 0 or c["end"] <= c["start"] or c["step"] <= 0:
        raise ValueError("Invalid checkpoint range")
    for key, allowed in (("rngs", RNGS), ("profiles", PROFILES), ("curves", CURVES)):
        choices = c.get(key)
        if not isinstance(choices, list) or not choices or not all(isinstance(x, str) and x in allowed for x in choices):
            raise ValueError("Unknown checkpoint algorithm")
        if len(set(choices)) != len(choices):
            raise ValueError("Duplicate checkpoint algorithm")
    for rng in c["rngs"]:
        for profile in c["profiles"]:
            validate_profile(rng, profile)
    if "match_fingerprints" not in c or (c["match_fingerprints"] is not None and
            (not isinstance(c["match_fingerprints"], list) or not all(valid_fingerprint(x) for x in c["match_fingerprints"]))):
        raise ValueError("Invalid checkpoint fingerprints")
    if type(data["next_index"]) is not int or not 0 <= data["next_index"] <= total_candidates(c):
        raise ValueError("Invalid checkpoint index")
    if type(data["emitted"]) is not int or not 0 <= data["emitted"] <= data["next_index"]:
        raise ValueError("Invalid checkpoint emitted count")
    if type(data["output_offset"]) is not int or data["output_offset"] < 0:
        raise ValueError("Invalid checkpoint output offset")
    if not valid_fingerprint(data["output_sha256"]):
        raise ValueError("Invalid checkpoint output hash")
    if not isinstance(data["output"], str) or "\0" in data["output"] or not Path(data["output"]).is_absolute():
        raise ValueError("Checkpoint output must be an absolute path")


@contextmanager
def checkpoint_lock(path: Path):
    # Keep this inode: unlinking a lock file allows waiters to lock different inodes.
    lock_path = path.with_name(path.name + ".lock")
    with os.fdopen(os.open(lock_path, os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW, 0o600), "r+b") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ValueError("Checkpoint is in use by another process") from exc
        yield lock_path


def run(args: argparse.Namespace) -> dict:
    checkpoint = getattr(args, "checkpoint", None)
    if checkpoint:
        cp = Path(checkpoint).resolve()
        with checkpoint_lock(cp) as lock_path:
            if args.command == "generate" and args.output and Path(args.output).resolve() in (cp, lock_path):
                raise ValueError("Output, checkpoint and checkpoint lock must differ")
            return run_locked(args)
    return run_locked(args)


def run_locked(args: argparse.Namespace) -> dict:
    if args.command != "estimate" and (args.max_candidates < 1 or args.checkpoint_every < 1):
        raise ValueError("Candidate and checkpoint limits must be positive")
    if args.command != "estimate" and (not math.isfinite(args.progress_every) or args.progress_every < 0):
        raise ValueError("Progress interval must be finite and nonnegative")
    if args.command == "resume":
        cp = Path(args.checkpoint).resolve()
        data = json.loads(cp.read_text())
        validate_checkpoint(data)
        config, index = dict(data["config"]), data["next_index"]
        config.setdefault("skip_draws", [0])
        config.setdefault("keys_per_seed", 1)
        output = Path(data["output"])
        if output.samefile(cp) or output.samefile(cp.with_name(cp.name + ".lock")):
            raise ValueError("Output, checkpoint and checkpoint lock must differ")
        # Acquire before verification/truncation so two resumptions cannot race.
        stream = output.open("r+b")
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
            if json.loads(cp.read_text()) != data:
                raise ValueError("Checkpoint changed before output lock; retry resume")
            remaining, digest = data["output_offset"], hashlib.sha256()
            while remaining:
                block = stream.read(min(1024 * 1024, remaining))
                if not block:
                    raise ValueError("Output is shorter than checkpoint")
                digest.update(block)
                remaining -= len(block)
            if digest.hexdigest() != data["output_sha256"]:
                raise ValueError("Checkpoint output prefix hash mismatch")
            stream.truncate(data["output_offset"])
            stream.seek(data["output_offset"])
        except BaseException:
            stream.close()
            raise
        emitted = data["emitted"]
    else:
        config = make_config(args)
        total = total_candidates(config)
        if args.command == "estimate":
            return {"config": config, "attempts": total, "timestamps": timestamp_count(config),
                    "note": "Attempts include rejected candidates; no RNG or EC generation performed. No universal timestamp key mapping."}
        index, emitted = 0, 0
        cp = Path(args.checkpoint).resolve() if args.checkpoint else None
        output = Path(args.output).resolve() if args.output else None
        if cp and not output:
            raise ValueError("--checkpoint requires --output")
        if cp and (cp == output or cp.exists()):
            raise ValueError("Checkpoint must be new and different from output")
        digest = hashlib.sha256()
        stream = os.fdopen(os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600), "wb") if output else sys.stdout.buffer
        if output:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
    total = total_candidates(config)
    selected = set(config["match_fingerprints"]) if config["match_fingerprints"] is not None else None
    stop = min(total, index + args.max_candidates)
    start_index = index
    rejected = nonmatching = 0
    started = last_report = time.monotonic()

    def statistics(now: float) -> dict:
        elapsed = max(0.0, now - started)
        attempted = index - start_index
        rate = attempted / elapsed if elapsed else 0.0
        return {"attempted_this_run": attempted, "next_index": index, "total_attempts": total,
                "emitted_total": emitted, "rejected_this_run": rejected, "nonmatching_this_run": nonmatching,
                "elapsed_seconds": elapsed, "attempts_per_second": rate, "complete": index == total}

    def save() -> None:
        stream.flush()
        if cp:
            os.fsync(stream.fileno())
            atomic_json(cp, {"schema_version": SCHEMA, "config": config, "next_index": index,
                             "emitted": emitted, "output": str(output), "output_offset": stream.tell(),
                             "output_sha256": digest.hexdigest()})

    try:
        previous_group, rng_stream = None, None
        while index < stop:
            timestamp, rng, profile, curve, skip, key_index = attempt_at(config, index)
            group = (timestamp, rng, profile, curve, skip)
            if group != previous_group:
                rng_stream = stream_at(rng, timestamp, profile, skip, key_index)
                previous_group = group
            record = candidate(timestamp, config["unit"], rng, profile, curve, skip, key_index, rng_stream)
            if record is None:
                rejected += 1
            elif selected is not None and record["spki_sha256"] not in selected:
                nonmatching += 1
            else:
                encoded = (json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n").encode()
                stream.write(encoded)
                digest.update(encoded)
                emitted += 1
            index += 1
            if (index - start_index) % args.checkpoint_every == 0:
                save()
            if args.progress_every:
                now = time.monotonic()
                if now - last_report >= args.progress_every:
                    progress = statistics(now)
                    rate = progress["attempts_per_second"]
                    progress.update(event="progress", remaining_this_run=stop - index,
                                    eta_seconds_this_run=(stop - index) / rate if rate else None)
                    print(json.dumps(progress, sort_keys=True), file=sys.stderr, flush=True)
                    last_report = now
        save()
    finally:
        if stream is not sys.stdout.buffer:
            stream.close()
    return statistics(time.monotonic())


def main() -> None:
    args = parser().parse_args()
    try:
        result = run(args)
    except (ValueError, OSError, KeyError, TypeError, OverflowError) as exc:
        print(f"timestampseeds: {exc}", file=sys.stderr)
        raise SystemExit(2)
    print(json.dumps(result, sort_keys=True), file=sys.stdout if args.command == "estimate" else sys.stderr)


if __name__ == "__main__":
    main()
