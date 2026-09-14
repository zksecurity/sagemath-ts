"""Bundled NTL CRT, determinant bound and arbitrary-modulus determinant bodies."""
import atexit
import hashlib
import json
from pathlib import Path
import select
import subprocess
import tempfile
from sage.env import SAGE_LOCAL
_process=None
class Error(Exception):pass
class RangeError(Exception):pass

def _build():
 root=Path(__file__).resolve().parents[3]
 source=(root/'tests/property/native/ntl_integer_crt.cpp').read_text()
 paths=[root/'reference/ntl/src'/n for n in ['ZZ.cpp','mat_ZZ.cpp','mat_ZZ_p.cpp','lzz_p.cpp','ZZ_p.cpp','lip.cpp']]
 zz,mat,big,ctx,bigctx,lip=[p.read_text() for p in paths]
 cut=lambda s,start,end:s[s.index(start):s.index(end,s.index(start))]
 body=cut(zz,'long CRT(ZZ& gg, ZZ& a, long G, long p)','void sub(ZZ& x, long a, const ZZ& b)')
 body+=cut(mat,'static\nlong DetBound(const mat_ZZ& a)','void determinant(ZZ& rres')
 body+=cut(mat,'long CRT(mat_ZZ& gg, ZZ& a, const mat_zz_p& G)','void mul(mat_ZZ& X, const mat_ZZ& A, const ZZ& b_in)')
 body+=cut(big,'void determinant(ZZ_p& d, const mat_ZZ_p& M_in)','long IsIdent(const mat_ZZ_p& A, long n)')
 body+=cut(ctx,'zz_pInfoT::zz_pInfoT(long NewP, long maxroot)','zz_pInfoT::zz_pInfoT(INIT_FFT_TYPE')
 body+=cut(ctx,'void zz_p::init(long p, long maxroot)','void zz_p::FFTInit')
 body+=cut(ctx,'zz_pContext::zz_pContext(long p, long maxroot)','zz_pContext::zz_pContext(INIT_FFT_TYPE')
 body+=cut(bigctx,'ZZ_pInfoT::ZZ_pInfoT(const ZZ& NewP)','// we use a lazy strategy')
 source=source.replace('/*__NATIVE_BODIES__*/',body)
 local=Path(SAGE_LOCAL)
 headers=['ZZ.h','mat_ZZ.h','mat_ZZ_p.h','mat_lzz_p.h','ZZ_p.h','lzz_p.h','matrix.h','vec_ZZVec.h']
 for name in headers:
  bundled=(root/'reference/ntl/include/NTL'/name).read_text();installed=(local/'include/NTL'/name).read_text()
  # Full interfaces are part of the cache key; context layouts must agree exactly.
  if name in ['ZZ_p.h','lzz_p.h']:
   start='class '+('ZZ_pInfoT {' if name=='ZZ_p.h' else 'zz_pInfoT {')
   if cut(bundled,start,'extern ')!=cut(installed,start,'extern '):raise RuntimeError('native CRT modulus layout differs from bundled source')
  paths.append(root/'reference/ntl/include/NTL'/name)
 compiler=subprocess.check_output(['c++','--version']);flags=['-std=c++17','-O2']
 version=b''.join((local/'include/NTL'/n).read_bytes() for n in ['version.h','config.h','mach_desc.h',*headers])
 key=hashlib.sha256(source.encode()+b''.join(p.read_bytes() for p in paths)+compiler+version+str(local).encode()+repr(flags).encode()).hexdigest()[:16]
 folder=Path(tempfile.gettempdir())/('sage-bundled-ntl-integer-crt-'+key);folder.mkdir(exist_ok=True);executable=folder/'oracle'
 if not executable.exists():
  path=folder/'oracle.cpp';path.write_text(source)
  subprocess.run(['c++',*flags,'-I'+str(local/'include'),str(path),'-L'+str(local/'lib'),'-lntl','-lgmp','-Wl,-rpath,'+str(local/'lib'),'-o',str(executable)],check=True,capture_output=True)
 return executable

def ntl_integer_crt(*args):
 global _process
 if _process is None:
  _process=subprocess.Popen([str(_build())],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True,bufsize=1);atexit.register(_process.terminate)
 compact=lambda x:'['+' '.join(map(str,x))+']'if isinstance(x,list)else str(x)
 _process.stdin.write(' '.join(map(compact,args))+'\n');_process.stdin.flush()
 if not select.select([_process.stdout],[],[],30)[0]:
  _process.terminate();_process.wait();_process=None;raise RuntimeError('bundled NTL integer CRT oracle timed out; no comparison result')
 line=_process.stdout.readline().rstrip('\n')
 if line.startswith('ERROR '):raise Error(json.loads(line[6:]))
 if line.startswith('RANGE '):raise RangeError(json.loads(line[6:]))
 if not line:
  _process.wait();_process=None;raise RuntimeError('bundled NTL integer CRT oracle failed without a result')
 return json.dumps(json.loads(line),separators=(',',':'))
