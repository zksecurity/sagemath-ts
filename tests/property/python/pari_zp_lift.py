"""Direct bundled native PARI scalar field comparisons."""
import atexit,json,hashlib,subprocess,tempfile
from pathlib import Path
from bundled_pari import bundled_pari_build
_process=None
class RangeError(Exception):pass
class PariError(Exception):pass

def call(values):
    global _process
    if _process is None:
        build=bundled_pari_build();obj=next(p.parent for p in build.glob('O*/pari.cfg'))
        source=Path(__file__).resolve().parents[1]/'native/pari_zp_lift.c'
        body=(Path(__file__).resolve().parents[3]/'reference/pari/src/basemath/FlxqE.c').read_text()
        parts=[]
        for name in ['RgX_circular_shallow','ZpXQ_frob_cyc','ZpXQ_frob','ZpXQ_sqrtnorm_pcyc']:
            start=body.rfind('static ',0,body.index('\n'+name+'('))
            brace=body.index('{',start);depth=1;end=brace+1
            while depth:
                depth+=(body[end]=='{')-(body[end]=='}');end+=1
            parts.append(body[start:end])
        helpers='\n'.join(parts)
        root_body=(Path(__file__).resolve().parents[3]/'reference/pari/src/basemath/Zp.c').read_text()
        root_parts=[]
        for name in ['mul2n','sqr2n','Fp_pow2n','Zp_sqrtnlift']:
            pos=root_body.index('\n'+name+'(')
            start=root_body.rfind('\nGEN' if name=='Zp_sqrtnlift' else 'static GEN',0,pos)
            brace=root_body.index('{',pos);depth=1;end=brace+1
            while depth:
                depth+=(root_body[end]=='{')-(root_body[end]=='}');end+=1
            root_parts.append(root_body[start:end])
        roots='\n'.join(root_parts).replace('Zp_sqrtnlift(', 'audit_Zp_sqrtnlift(',1)
        key=hashlib.sha256(source.read_bytes()+helpers.encode()+roots.encode()+str(build).encode()).hexdigest()[:16]
        folder=Path(tempfile.gettempdir())/('sage-bundled-pari-zp-lift-'+key);folder.mkdir(exist_ok=True)
        executable=folder/'oracle'
        if not executable.exists():
            (folder/'pari_padic_frobenius.h').write_text(helpers)
            (folder/'pari_guarded_root.h').write_text(roots)
            library=next(p for p in obj.glob('libpari*')if p.suffix in ('.dylib','.so'))
            compiled=subprocess.run(['cc','-O2','-I'+str(folder),'-I'+str(obj),'-I'+str(build/'src/headers'),str(source),str(library),
                            '-Wl,-rpath,'+str(obj),'-o',str(executable)],capture_output=True,text=True)
            if compiled.returncode:raise RuntimeError('PARI lifting oracle compilation failed: '+compiled.stderr[-4000:])
        _process=subprocess.Popen([str(executable)],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1)
        atexit.register(_process.terminate)
    compact=lambda v:'['+','.join(compact(c)for c in v)+']' if isinstance(v,list) else str(v)
    _process.stdin.write(' '.join(compact(v)for v in values)+'\n');_process.stdin.flush()
    line=_process.stdout.readline().rstrip('\n')
    if line.startswith('GUARD '):raise RangeError(line[6:])
    if line.startswith('ERROR '):raise PariError(json.loads(line[6:]))
    if not line.startswith('OK '):raise RuntimeError('bundled PARI scalar field oracle failed: '+line)
    return line[3:]

def native_pari_zp_lift(op,f,a,t,Q,p,e):
    return call([op,f,a,t,Q,p,e])
