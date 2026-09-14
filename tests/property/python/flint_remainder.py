import ctypes,hashlib,tempfile,subprocess,json,sys
from pathlib import Path
from bundled_flint import bundled_flint_build
_native=None

def native_rem(op,a,b,da,db):
 global _native
 if _native is None:
  build=bundled_flint_build();source=Path(__file__).resolve().parents[1] / 'native' / 'flint_remainder.c';key=hashlib.sha256(source.read_bytes()+str(build).encode()).hexdigest()[:16];folder=Path(tempfile.gettempdir())/('sage-bundled-remainder-'+key);folder.mkdir(exist_ok=True);lib=folder/'oracle.so'
  if not lib.exists():subprocess.run(['cc','-shared','-fPIC','-O2','-I'+str(build/'src'),str(source),'-L'+str(build),'-lflint','-Wl,-rpath,'+str(build),'-o',str(lib)],check=True,capture_output=True)
  _native=ctypes.CDLL(str(lib));_native.audit_remainder.argtypes=[ctypes.c_int,*[ctypes.c_char_p]*4];_native.audit_remainder.restype=ctypes.c_void_p;_native.audit_remainder_free.argtypes=[ctypes.c_void_p]
 encode=lambda a:(str(len(a))+'  '+' '.join(map(str,a))).encode()
 p=_native.audit_remainder(int(op),encode(a),encode(b),str(da).encode(),str(db).encode())
 try:q,r,d=ctypes.string_at(p).decode().splitlines()
 finally:_native.audit_remainder_free(p)
 values=lambda s:s.split()[1:]
 return json.dumps([values(q),values(r),d]if op<2 else[values(r),d],separators=(',',':'))
