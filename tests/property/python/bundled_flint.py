"""Locate or build the bundled FLINT in a temporary cache, without installing it."""
import fcntl
import hashlib
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile


def bundled_flint_build():
    override=os.environ.get('SAGEMATH_AUDIT_FLINT_BUILD')
    if override:return Path(override)
    source=Path(__file__).resolve().parents[3]/'reference/flint'
    digest=hashlib.sha256(platform.machine().encode())
    for path in sorted(source.rglob('*')):
        if path.is_file() and path.suffix in ('.c','.h','.in','.ac'):
            digest.update(str(path.relative_to(source)).encode());digest.update(path.read_bytes())
    root=Path(tempfile.gettempdir())/('sage-bundled-flint-'+digest.hexdigest()[:16])
    root.mkdir(exist_ok=True)
    build=root/'src'
    with (root/'build.lock').open('w')as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if (root/'complete').exists():return build
        if not build.exists():shutil.copytree(source,build)
        prefix='/opt/homebrew' if Path('/opt/homebrew/include/gmp.h').exists() else '/usr/local'
        commands=[['./bootstrap.sh'],['./configure','--disable-static','--disable-assembly']]
        if Path(prefix+'/include/gmp.h').exists():commands[1]+=['--with-gmp='+prefix,'--with-mpfr='+prefix]
        commands.append(['make','-j4'])
        with (root/'build.log').open('w')as log:
            for command in commands:
                result=subprocess.run(command,cwd=build,stdout=log,stderr=subprocess.STDOUT)
                if result.returncode:
                    log.flush()
                    raise RuntimeError('Bundled FLINT build failed; see '+str(root/'build.log'))
        (root/'complete').touch()
    return build
