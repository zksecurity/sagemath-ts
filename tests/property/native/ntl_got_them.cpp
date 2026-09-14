#include <NTL/FFT.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/BasicThreadPool.h>
#include <sys/wait.h>
#include <unistd.h>
#include <cerrno>
#include <cstdint>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstring>
#include "ZZXFactoring.cpp"
#if defined(NTL_THREADS) || defined(NTL_RANDOM_AES256CTR) || defined(NTL_HAVE_AVX2) || defined(NTL_HAVE_SSSE3)
#error Unsupported native factor state profile
#endif
static_assert(NTL_BITS_PER_LONG==64 && NTL_SP_NBITS==60 && NTL_ZZ_NBITS==64 && NTL_BITS_PER_INT32==32 && NTL_BITS_PER_NONCE==64,"unsupported native word profile");
static_assert(NTL_FDOUBLE_PRECISION==double(1ULL<<52),"unsupported probability profile");
namespace NTL { void NextFFTPrime(long&, long&, long); }
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
static void emit(std::ostream&out,const ZZX&f){out<<'[';for(long i=0;i<=deg(f);i++){if(i)out<<',';emit(out,f[i]);}out<<']';}
static void emit(std::ostream&out,const vec_ZZX&f){out<<'[';for(long i=0;i<f.length();i++){if(i)out<<',';emit(out,f[i]);}out<<']';}
static std::string hextext(const std::string&s){return hexbytes(std::vector<unsigned char>(s.begin(),s.end()));}
static void one(const std::string&line){try{
 std::istringstream in(line);vec_ZZ key,flat,cmd;vec_long dims;vec_ZZX input,previous;ZZX f;ZZ bound,p;
 in>>key>>dims>>flat>>input>>f>>previous>>bound>>p>>cmd;
 if(!in||key.length()!=32||dims.length()!=2||cmd.length()%2)throw std::runtime_error("invalid GotThem input");
 unsigned char bytes[32];for(long i=0;i<32;i++)bytes[i]=conv<unsigned long>(key[i]);RandomStream initial(bytes);SetSeed(initial);RandomStream&stream=GetCurrentRandomStream();
 mat_ZZ B;B.SetDims(dims[0],dims[1]);if(flat.length()!=dims[0]*dims[1])throw std::runtime_error("invalid GotThem matrix length");for(long i=0,t=0;i<dims[0];i++)for(long j=0;j<dims[1];j++)B[i][j]=flat[t++];
 ZZ_p::init(p);vec_ZZ_pX W;W.SetLength(input.length());for(long i=0;i<input.length();i++)conv(W[i],input[i]);
 std::ostringstream out;out<<'[';
 for(long i=0;i<cmd.length();i+=2){if(i)out<<',';long op=conv<long>(cmd[i]);const ZZ&a=cmd[i+1];std::ostringstream value,diagnostics;const char*kind=nullptr;std::string message;
 const zz_pInfoT*oldWord=zz_pInfo;const ZZ_pInfoT*oldBig=ZZ_pInfo;auto* oldDiagnostics=std::cerr.rdbuf(diagnostics.rdbuf());
 try{
  if(op==0){if(bound<-power2_ZZ(63)||bound>=power2_ZZ(63))throw AdapterGuard("GotThem bound must be a signed native integer");vec_ZZX factors=previous;long status=GotThem(factors,B,W,f,conv<long>(bound),0);value<<'[';scalar(value,status);value<<',';emit(value,factors);value<<']';}
  else if(op==1){UseFFTPrime(conv<long>(a));value<<"null";}
  else if(op==2){stream.set_nonce(conv<unsigned long>(a));value<<"null";}
  else if(op==3){long n=conv<long>(a);std::vector<unsigned char>v(std::max(0L,n));stream.get(v.data(),n);value<<std::quoted(hexbytes(v));}
  else if(op==4){zz_p::FFTInit(conv<long>(a));value<<"null";}
  else throw std::runtime_error("unknown GotThem operation");
 }catch(const AdapterGuard&e){kind="RangeError";message=e.what();}catch(const std::exception&e){kind="Error";message=e.what();}
 std::cerr.rdbuf(oldDiagnostics);
 if(op==0&&(zz_pInfo!=oldWord||ZZ_pInfo!=oldBig))throw std::runtime_error("native GotThem leaked a modular context");
 out<<'[';if(kind)out<<std::quoted(kind)<<','<<std::quoted(message)<<",null";else out<<"null,null,"<<value.str();
 out<<",[";for(long j=0;j<FFTTables.length();j++){if(j)out<<',';info(out,*FFTTables[j]);}out<<']';
 RandomStream copy(stream);std::vector<unsigned char>tail(64);copy.get(tail.data(),64);out<<','<<std::quoted(hexbytes(tail))<<','<<std::quoted(hextext(diagnostics.str()))<<']';
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
