"""CPython 3.12/Unicode 15 oracle matching bundled Sage's Python dependency.

Set PYTHON312 to the executable if it is not discoverable as python3.12.
The text conversion stays in the original runtime; Sage 10.3 uses older Unicode data.
"""
import atexit
import json
import os
import shutil
import subprocess
import sys
from functools import lru_cache

_process = None

@lru_cache(maxsize=None)
def text_float(codes, as_bytes=False):
    global _process
    if _process is None:
        executable = os.environ.get('PYTHON312') or shutil.which('python3.12')
        if not executable:
            raise RuntimeError('CPython 3.12 is required for the bundled float oracle; set PYTHON312')
        _process = subprocess.Popen([executable, __file__, '--serve'], stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE, text=True)
        atexit.register(_process.terminate)
    _process.stdin.write(json.dumps([codes, as_bytes])+'\n')
    _process.stdin.flush()
    line = _process.stdout.readline()
    if not line:
        raise RuntimeError('CPython float oracle exited without a result')
    return json.loads(line)

if __name__ == '__main__' and '--serve' in sys.argv:
    import struct
    import unicodedata
    if sys.version_info[:2] != (3, 12) or unicodedata.unidata_version != '15.0.0':
        raise RuntimeError('float oracle requires CPython 3.12 with Unicode 15.0.0')
    for line in sys.stdin:
        codes, as_bytes = json.loads(line)
        value = bytes(codes) if as_bytes else ''.join(chr(c) for c in codes)
        try:
            result = ['ok', struct.pack('>d', float(value)).hex()]
        except Exception as error:
            result = ['error', type(error).__name__, str(error)]
        print(json.dumps(result), flush=True)
