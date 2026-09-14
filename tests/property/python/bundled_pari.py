"""Build the original PARI library in a temporary cache without installing it."""
import fcntl
import hashlib
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tempfile


def bundled_pari_build():
    override=os.environ.get('SAGEMATH_AUDIT_PARI_BUILD')
    if override:return Path(override)
    source=Path(__file__).resolve().parents[3]/'reference/pari'
    digest=hashlib.sha256(platform.machine().encode())
    for path in sorted(source.rglob('*')):
        if path.is_file() and path.suffix in ('.c','.h','.in','.sh'):
            digest.update(str(path.relative_to(source)).encode());digest.update(path.read_bytes())
    root=Path(tempfile.gettempdir())/('sage-bundled-pari-'+digest.hexdigest()[:16])
    root.mkdir(exist_ok=True)
    build=root/'src'
    with (root/'build.lock').open('w')as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if (root/'complete').exists():return build
        if not build.exists():shutil.copytree(source,build)
        command=['./Configure','--prefix='+str(root/'unused-install')]
        prefix='/opt/homebrew' if Path('/opt/homebrew/include/gmp.h').exists() else '/usr/local'
        if Path(prefix+'/include/gmp.h').exists():command+=['--with-gmp='+prefix]
        with (root/'build.log').open('w')as log:
            env=os.environ.copy()
            if Path('/opt/homebrew/opt/bison/bin/bison').exists():env['PATH']='/opt/homebrew/opt/bison/bin:'+env['PATH']
            subprocess.run(command,cwd=build,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
            subprocess.run(['make','-j4','all'],cwd=build,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        (root/'complete').touch()
    return build
