#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/vec_ZZVec.h>
#define PAR_THRESH (40000.0)
#include <NTL/FFT.h>
#include <sys/wait.h>
#include <unistd.h>
#include <cerrno>
#include <cstdint>
#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/Lazy.h>
#include <NTL/fileio.h>
#include <NTL/SmartPtr.h>
#include <NTL/BasicThreadPool.h>
#include <sstream>
#include <iomanip>
#include <memory>
#include <vector>
#include <cstring>
#if defined(NTL_RANDOMSTREAM_AES256CTR) || defined(NTL_HAVE_AVX2) || defined(NTL_HAVE_SSSE3) || NTL_BITS_PER_LONG != 64 || NTL_BITS_PER_INT32 != 32
#error Unsupported native byte stream profile
#endif
static_assert(NTL_SP_NBITS==60 && NTL_ZZ_NBITS==64 && NTL_BITS_PER_NONCE==64, "unsupported prime word profile");
static_assert(NTL_FDOUBLE_PRECISION==double(1ULL<<52), "unsupported probability-bound profile");
namespace NTL {
NTL_IMPORT_FROM_STD
/*__NATIVE_RANDOM_BODY__*/
static RandomStream *active_stream;
RandomStream& GetCurrentRandomStream() { return *active_stream; }
static inline RandomStream& LocalGetCurrentRandomStream() { return *active_stream; }
// The explicit active stream is initialized before seeding; preserve assignment identity.
void SetSeed(const RandomStream& s) { *active_stream=s; }
/*__NATIVE_SEED_BODY__*/
/*__NATIVE_SAMPLING_BODY__*/
/*__NATIVE_PRIME_SEQUENCE__*/
/*__NATIVE_PRIME_BODY__*/
}

namespace NTL {
FFTTablesType FFTTables;
/*__NATIVE_FFT_PRIMES__*/
}

namespace NTL {
NTL_IMPORT_FROM_STD
/*__NATIVE_RECONSTRUCTION_BODY__*/
/*__NATIVE_CERTIFIED_GAUSS__*/
}
using namespace NTL;
static void error(const char*x){throw std::runtime_error(x);}
struct AdapterGuard:std::runtime_error{using std::runtime_error::runtime_error;};
static std::string hexbytes(const std::vector<unsigned char>&b){static const char*h="0123456789abcdef";std::string s;for(auto x:b){s+=h[x>>4];s+=h[x&15];}return s;}
static void scalar(std::ostream&out,long x){out<<'"'<<x<<'"';}
static void vector(std::ostream&out,const vec_long&v){out<<'[';for(long i=0;i<v.length();i++){if(i)out<<',';scalar(out,v[i]);}out<<']';}
static void reciprocal(std::ostream&out,double x){uint64_t bits;memcpy(&bits,&x,8);out<<'"'<<std::hex<<std::setfill('0')<<std::setw(16)<<bits<<std::dec<<'"';}
static void info(std::ostream&out,const FFTPrimeInfo&v){out<<'[';scalar(out,v.q);out<<',';reciprocal(out,v.qrecip);out<<',';vector(out,v.RootTable[0]);out<<',';vector(out,v.RootTable[1]);out<<',';vector(out,v.TwoInvTable);out<<']';}
static void emit(std::ostream&out,const ZZ& x){out<<'"'<<x<<'"';}
static void emit(std::ostream&out,const mat_ZZ& M){out<<'[';for(long i=0;i<M.NumRows();i++){if(i)out<<',';out<<'[';for(long j=0;j<M.NumCols();j++){if(j)out<<',';emit(out,M[i][j]);}out<<']';}out<<']';}
static void make(mat_ZZ&M,const vec_ZZ&v,long n,long m,bool ragged){
 if(ragged){Vec<vec_ZZ> rows;rows.SetLength(2);rows[0].SetLength(1);rows[1].SetLength(2);MakeMatrix(M,rows);return;}
 M.SetDims(n,m);if(v.length()!=n*m)throw std::runtime_error("invalid flat matrix length");
 for(long i=0,t=0;i<n;i++)for(long j=0;j<m;j++)M[i][j]=v[t++];
}
static void one(const std::string&line){try{
 std::istringstream in(line);vec_ZZ key,av,prev,cmd;vec_long dims;in>>key>>dims>>av>>prev>>cmd;
 if(!in||key.length()!=32||dims.length()!=6||cmd.length()%2)throw std::runtime_error("invalid NTL integer reconstruction input");
 unsigned char bytes[32];for(long i=0;i<32;i++)bytes[i]=conv<unsigned long>(key[i]);RandomStream stream(bytes);active_stream=&stream;
 mat_ZZ A,previous;make(A,av,dims[0],dims[1],dims[4]);make(previous,prev,dims[2],dims[3],dims[5]);
 std::ostringstream out;out<<'[';
 for(long i=0;i<cmd.length();i+=2){if(i)out<<',';long op=conv<long>(cmd[i]);const ZZ&a=cmd[i+1];std::ostringstream value;const char*kind=nullptr;std::string message;
 const zz_pInfoT*oldWord=zz_pInfo;const ZZ_pInfoT*oldBig=ZZ_pInfo;
 try{
  if(op==0||op==1){ZZ d;if(op==0)determinant(d,A);else determinant(d,A,conv<long>(a));emit(value,d);}
  else if(op==2){ZZ d;mat_ZZ X=previous;inv(d,X,A,conv<long>(a));value<<'[';emit(value,d);value<<',';emit(value,X);value<<']';}
  else if(op==3){mat_ZZ X=previous;inv(X,A);emit(value,X);}
  else if(op==4){UseFFTPrime(conv<long>(a));value<<"null";}
  else if(op==5){long q,w;NextFFTPrime(q,w,conv<long>(a));value<<'['<<'"'<<q<<"\",\""<<w<<"\"]";}
  else if(op==6){long n=conv<long>(a);std::vector<unsigned char>v(std::max(0L,n));stream.get(v.data(),n);value<<std::quoted(hexbytes(v));}
  else if(op==7){stream.set_nonce(conv<unsigned long>(a));value<<"null";}
  else if(op==8){zz_p::FFTInit(conv<long>(a));value<<"null";}
  else if(op==9){ZZ_p::init(a);value<<"null";}
  else if(op==10){ZZ P;GenPrime(P,conv<long>(a),90);emit(value,P);}
  else if(op==11){ZZ d;mat_ZZ R;gauss(d,R,A);value<<'[';emit(value,d);value<<',';emit(value,R);value<<']';}
  else if(op==12){scalar(value,GenPrime_long(conv<long>(a)));}
  else throw std::runtime_error("unknown integer reconstruction operation");
 }catch(const AdapterGuard&e){kind="RangeError";message=e.what();}catch(const std::exception&e){kind="Error";message=e.what();}
 if((op<=3||op==11)&&(zz_pInfo!=oldWord||ZZ_pInfo!=oldBig))throw std::runtime_error("native integer reconstruction leaked a modular context");
 out<<'[';if(kind)out<<std::quoted(kind)<<','<<std::quoted(message)<<",null";else out<<"null,null,"<<value.str();
 out<<",[";for(long j=0;j<FFTTables.length();j++){if(j)out<<',';info(out,*FFTTables[j]);}out<<']';
 RandomStream copy(stream);std::vector<unsigned char>tail(64);copy.get(tail.data(),64);out<<','<<std::quoted(hexbytes(tail))<<']';
 }
 out<<']';std::cout<<out.str()<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}
int main(){ErrorMsgCallback=error;SetNumThreads(1);std::string line;while(std::getline(std::cin,line)){
 // Each trace starts with the original global cache and static search state.
 // Keep the parent pristine; all operations in a trace share one child/cache.
 pid_t child=fork();if(child<0){std::cout<<"ERROR \"native integer reconstruction fork failed\""<<std::endl;continue;}
 if(child==0){one(line);_exit(0);}int status=0;pid_t waited;do{waited=waitpid(child,&status,0);}while(waited<0&&errno==EINTR);
 if(waited<0||!WIFEXITED(status)||WEXITSTATUS(status)!=0){std::cout<<"ERROR \"native integer reconstruction child failed\""<<std::endl;}
}}
