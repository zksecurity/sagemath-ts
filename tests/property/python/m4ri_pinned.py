"""Build the checksum-pinned original M4RI for isolated native boundary tests."""
from pathlib import Path
import fcntl,hashlib,os,subprocess,tarfile,tempfile

def executable():
    root=Path(__file__).resolve().parents[3]
    metadata=root/'reference/sage/build/pkgs/m4ri'
    version=(metadata/'package-version.txt').read_text().strip()
    fields=dict(line.split('=',1) for line in (metadata/'checksums.ini').read_text().splitlines() if '=' in line)
    checksum=fields['sha256'];url=fields['upstream_url'].replace('VERSION',version)
    cache=Path(tempfile.gettempdir())/'sagemath-ts-native'/('m4ri-'+version+'-'+checksum[:12]);cache.mkdir(parents=True,exist_ok=True)
    program=Path(__file__).with_name('m4ri_strassen_native.c')
    target=cache/('oracle-'+hashlib.sha256(program.read_bytes()).hexdigest()[:16])
    with (cache/'build.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        if target.exists():return target
        archive=cache/('m4ri-'+version+'.tar.gz')
        if not archive.exists():
            subprocess.run(['curl','--fail','--location','--silent','--show-error',url,'--output',str(archive)],check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if hashlib.sha256(archive.read_bytes()).hexdigest()!=checksum:raise RuntimeError('M4RI archive checksum mismatch')
        source=cache/('m4ri-'+version)
        if not source.exists():
            with tarfile.open(archive) as tar:
                for member in tar.getmembers():
                    path=Path(member.name)
                    if path.is_absolute() or '..' in path.parts:raise RuntimeError('invalid M4RI archive path')
                tar.extractall(cache)
        def run(command):
            completed=subprocess.run(command,cwd=source,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            if completed.returncode:raise RuntimeError('M4RI oracle build failed:\n'+completed.stdout)
        library=source/'.libs/libm4ri.a'
        if not library.exists():
            run(['./configure','--disable-shared','--enable-static','--disable-png'])
            run(['make','-j4'])
        run(['cc','-O2','-I'+str(source),str(program),str(library),'-lm','-o',str(target)])
    return target

def strassen(m,k,n,seed,cutoff,square):
    run=subprocess.run([str(executable()),*map(str,[m,k,n,seed,cutoff,square])],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if run.returncode:
        if run.returncode==-6:raise RuntimeError('Aborted')
        if run.returncode==-11:raise RuntimeError('Segmentation fault')
        raise RuntimeError('M4RI oracle failed: '+run.stderr)
    return [str(m),str(n),[str(int(line,16)) for line in run.stdout.splitlines()]]

def density(m,n,pattern,seed,resolution):
    run=subprocess.run([str(executable()),'density',*map(str,[m,n,pattern,seed,resolution])],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if run.returncode:
        if run.returncode==-6:raise RuntimeError('Aborted')
        if run.returncode==-11:raise RuntimeError('Segmentation fault')
        raise RuntimeError('M4RI oracle failed: '+run.stderr)
    return float.fromhex(run.stdout.strip())

def echelon(m,n,seed,algorithm,full,k):
    run=subprocess.run([str(executable()),'echelon',*map(str,[m,n,seed,algorithm,full,k])],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if run.returncode:
        if run.returncode==-6:raise RuntimeError('Aborted')
        if run.returncode==-11:raise RuntimeError('Segmentation fault')
        raise RuntimeError('M4RI oracle failed: '+run.stderr)
    lines=run.stdout.splitlines()
    return [int(lines[0]),[str(int(x,16)) for x in lines[1:]]]

def ple(m,n,seed,kind,k):
    run=subprocess.run([str(executable()),'ple',*map(str,[m,n,seed,kind,k])],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if run.returncode:
        if run.returncode==-11:raise RuntimeError('Segmentation fault')
        if run.returncode==-6:raise RuntimeError('Aborted')
        raise RuntimeError('M4RI PLE oracle failed: '+run.stderr)
    lines=run.stdout.splitlines()
    return [int(lines[0]),[int(x) for x in lines[1].split()],[int(x) for x in lines[2].split()],[str(int(x,16)) for x in lines[3:]]]

def trsm(m,n,seed,upper,cutoff):
    run=subprocess.run([str(executable()),'trsm',*map(str,[m,n,seed,upper,cutoff])],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if run.returncode:raise RuntimeError('M4RI TRSM oracle failed: '+run.stderr)
    return [str(int(x,16)) for x in run.stdout.splitlines()]

def empty_swap(m,n,axis,i,j):
    run=subprocess.run([str(executable()),'swap',*map(str,[m,n,1,axis,i,j])],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if run.returncode:raise RuntimeError('M4RI empty swap oracle failed: '+run.stderr)

def echelon_pattern(m,n,seed,pattern,algorithm,full,k):
    run=subprocess.run([str(executable()),'echelon_pattern',*map(str,[m,n,seed,pattern,algorithm,full,k])],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if run.returncode:raise RuntimeError('M4RI patterned echelon oracle failed: '+run.stderr)
    lines=run.stdout.splitlines();return [int(lines[0]),[str(int(x,16)) for x in lines[1:]]]

def factorization(m,n,seed,kind,k):
    run=subprocess.run([str(executable()),'factorization',*map(str,[m,n,seed,kind,k])],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=20)
    if run.returncode:
        if run.returncode==-6:raise RuntimeError('Aborted')
        if run.returncode==-11:raise RuntimeError('Segmentation fault')
        if run.returncode==-8:raise RuntimeError('Floating point exception')
        raise RuntimeError('M4RI factorization oracle failed: '+run.stderr)
    lines=run.stdout.splitlines()
    return [int(lines[0]),[int(x) for x in lines[1].split()],[int(x) for x in lines[2].split()],[str(int(x,16)) for x in lines[3:]]]

def inverse(n,seed,pattern,k):
    run=subprocess.run([str(executable()),'inverse',*map(str,[n,seed,pattern,k])],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if run.returncode:raise RuntimeError('M4RI inverse oracle failed: '+run.stderr)
    return [str(int(x,16)) for x in run.stdout.splitlines()]

def solve(m,n,bc,seed,pattern,cutoff,check):
    run=subprocess.run([str(executable()),'solve',*map(str,[m,n,bc,seed,pattern,cutoff,check])],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if run.returncode:
        if run.returncode==-11:raise RuntimeError('Segmentation fault')
        if run.returncode==-6:raise RuntimeError('Aborted')
        raise RuntimeError('M4RI solve oracle failed: '+run.stderr)
    lines=run.stdout.splitlines()
    return [int(lines[0]),[str(int(x,16)) for x in lines[1:m+1]],[str(int(x,16)) for x in lines[m+1:]]]

def kernel(m,n,seed,pattern,cutoff):
    run=subprocess.run([str(executable()),'kernel',*map(str,[m,n,seed,pattern,cutoff])],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if run.returncode:
        if run.returncode==-11:raise RuntimeError('Segmentation fault')
        if run.returncode==-6:raise RuntimeError('Aborted')
        raise RuntimeError('M4RI kernel oracle failed: '+run.stderr)
    lines=run.stdout.splitlines();width=int(lines[0])
    return [[str(int(x,16)) for x in lines[1:m+1]],None if width<0 else [width,[str(int(x,16)) for x in lines[m+1:]]]]

def basic(m,n,right,seed,op):
    run=subprocess.run([str(executable()),'basic',*map(str,[m,n,right,seed,op])],text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    if run.returncode:
        if run.returncode==-11:raise RuntimeError('Segmentation fault')
        if run.returncode==-6:raise RuntimeError('Aborted')
        raise RuntimeError('M4RI basic oracle failed: '+run.stderr)
    lines=run.stdout.splitlines();shape=[int(x) for x in lines[0].split()]
    return [*shape,[str(int(x,16)) for x in lines[1:]]]
