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
using namespace NTL;
static void error(const char*x){throw std::runtime_error(x);}
struct AdapterGuard:std::runtime_error{using std::runtime_error::runtime_error;};
static std::string hexbytes(const std::vector<unsigned char>&b){static const char*h="0123456789abcdef";std::string s;for(auto x:b){s+=h[x>>4];s+=h[x&15];}return s;}
static void scalar(std::ostream&out,long x){out<<'"'<<x<<'"';}
static void vector(std::ostream&out,const vec_long&v){out<<'[';for(long i=0;i<v.length();i++){if(i)out<<',';scalar(out,v[i]);}out<<']';}
static void reciprocal(std::ostream&out,double x){uint64_t bits;memcpy(&bits,&x,8);out<<'"'<<std::hex<<std::setfill('0')<<std::setw(16)<<bits<<std::dec<<'"';}
static void info(std::ostream&out,const FFTPrimeInfo&v){out<<'[';scalar(out,v.q);out<<',';reciprocal(out,v.qrecip);out<<',';vector(out,v.RootTable[0]);out<<',';vector(out,v.RootTable[1]);out<<',';vector(out,v.TwoInvTable);out<<']';}
static void one(const std::string&line){try{
 std::istringstream in(line);vec_ZZ key,cmd;in>>key>>cmd;if(!in||key.length()!=32||cmd.length()%4)throw std::runtime_error("invalid NTL FFT-prime input");
 unsigned char bytes[32];for(long i=0;i<32;i++)bytes[i]=conv<unsigned long>(key[i]);RandomStream stream(bytes);active_stream=&stream;
 std::ostringstream out;out<<'[';
 for(long i=0;i<cmd.length();i+=4){if(i)out<<',';long op=conv<long>(cmd[i]);const ZZ&a=cmd[i+1];const ZZ&b=cmd[i+2];std::ostringstream value;const char*kind=nullptr;std::string message;
 try{
  if(op==0||op==11){if(a<-(ZZ(1)<<63)||a>=(ZZ(1)<<63))throw AdapterGuard("IsFFTPrime: input must fit a signed 64-bit integer");long w=op==11?0:conv<long>(b);long yes=IsFFTPrime(conv<long>(a),w);value<<'['<<yes<<',';scalar(value,w);value<<']';}
  else if(op==1){if(a<=-(ZZ(1)<<63)||a>=(ZZ(1)<<63))throw AdapterGuard("CalcMaxRoot: p-1 must fit a signed 64-bit integer");if(a==1)throw std::runtime_error("CalcMaxRoot: p=1 has no terminating native result");value<<CalcMaxRoot(conv<long>(a));}
  else if(op==2){long q,w;NextFFTPrime(q,w,conv<long>(a));value<<'[';scalar(value,q);value<<',';scalar(value,w);value<<']';}
  else if(op==3){UseFFTPrime(conv<long>(a));value<<"null";}
  else if(op==4||op==5||op==6){long index=conv<long>(a);if(index<0||index>=FFTTables.length())throw AdapterGuard("FFTPrimeContext: prime index is not initialized");if(op==4)scalar(value,GetFFTPrime(index));else if(op==5)reciprocal(value,GetFFTPrimeRecip(index));else info(value,*FFTTables[index]);}
  else if(op==7){if(a<=1||a>=ZZ(1)<<60)throw AdapterGuard("InitFFTPrimeInfo: q must satisfy 1 < q < 2^60");if(b<0||b>=a)throw AdapterGuard("InitFFTPrimeInfo: root must satisfy 0 <= w < q");FFTPrimeInfo I;InitFFTPrimeInfo(I,conv<long>(a),conv<long>(b),-1);info(value,I);}
  else if(op==8)value<<FFTTables.length();
  else if(op==9){long n=conv<long>(a);std::vector<unsigned char>v(std::max(0L,n));stream.get(v.data(),n);value<<std::quoted(hexbytes(v));}
  else if(op==10){stream.set_nonce(conv<unsigned long>(a));value<<"null";}
  else throw std::runtime_error("unknown FFT-prime operation");
 }catch(const AdapterGuard&e){kind="RangeError";message=e.what();}catch(const std::exception&e){kind="Error";message=e.what();}
 out<<'[';if(kind)out<<std::quoted(kind)<<','<<std::quoted(message)<<",null";else out<<"null,null,"<<value.str();
 RandomStream copy(stream);std::vector<unsigned char>tail(64);copy.get(tail.data(),64);out<<','<<FFTTables.length()<<','<<std::quoted(hexbytes(tail))<<']';
 }
 out<<']';std::cout<<out.str()<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}
int main(){ErrorMsgCallback=error;SetNumThreads(1);std::string line;while(std::getline(std::cin,line)){
 // Each trace starts with the original global cache and static search state.
 // Keep the parent pristine; all operations in a trace share one child/cache.
 pid_t child=fork();if(child<0){std::cout<<"ERROR \"native FFT-prime fork failed\""<<std::endl;continue;}
 if(child==0){one(line);_exit(0);}int status=0;pid_t waited;do{waited=waitpid(child,&status,0);}while(waited<0&&errno==EINTR);
 if(waited<0||!WIFEXITED(status)||WEXITSTATUS(status)!=0){std::cout<<"ERROR \"native FFT-prime child failed\""<<std::endl;}
}}
