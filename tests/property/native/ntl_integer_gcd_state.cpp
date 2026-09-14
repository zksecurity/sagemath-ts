#include <NTL/FFT.h>
#include <NTL/FFT_impl.h>
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
#include <NTL/ZZXFactoring.h>
#include <NTL/lzz_pXFactoring.h>
#if !defined(__aarch64__)
#error This fixture profile requires AArch64 fused CRT arithmetic
#endif
#include <NTL/ZZX.h>
#if defined(NTL_THREADS) || defined(NTL_RANDOM_AES256CTR) || defined(NTL_HAVE_AVX2) || defined(NTL_HAVE_SSSE3)
#error Unsupported native factor state profile
#endif
static_assert(NTL_BITS_PER_LONG==64 && NTL_SP_NBITS==60 && NTL_ZZ_NBITS==64 && NTL_BITS_PER_INT32==32 && NTL_BITS_PER_NONCE==64,"unsupported native word profile");
static_assert(NTL_FDOUBLE_PRECISION==double(1ULL<<52),"unsupported probability profile");
namespace NTL { void NextFFTPrime(long&, long&, long); }
using namespace NTL;
using std::cerr;
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
static void emit(std::ostream&out,const zz_pX&f){out<<'[';for(long i=0;i<=deg(f);i++){if(i)out<<',';scalar(out,rep(f[i]));}out<<']';}
static void emit(std::ostream&out,const vec_pair_ZZX_long& a){out<<'[';for(long i=0;i<a.length();i++){if(i)out<<',';out<<'[';emit(out,a[i].a);out<<','<<a[i].b<<']';}out<<']';}
static void one(const std::string&line){try{
 std::istringstream in(line);vec_ZZ key,coeffs,cmd;vec_long lengths;long p,maxroot;in>>key>>p>>maxroot>>lengths>>coeffs>>cmd;
 if(!in||key.length()!=32||cmd.length()%2)throw std::runtime_error("invalid integer GCD state input");
 unsigned char bytes[32];for(long i=0;i<32;i++)bytes[i]=conv<unsigned long>(key[i]);SetSeed(RandomStream(bytes));RandomStream&stream=GetCurrentRandomStream();
 // p=0 selects no initial word context, for cold integer GCD calls.
 if(p!=0)zz_p::init(p,maxroot);
 std::vector<ZZX> raw;long offset=0;
 for(long i=0;i<lengths.length();i++){ZZX r;for(long j=0;j<lengths[i];j++)SetCoeff(r,j,coeffs[offset++]);raw.push_back(r);}
 if(offset!=coeffs.length())throw std::runtime_error("invalid polynomial lengths");
 std::ostringstream out;out<<'[';
 for(long i=0;i<cmd.length();i+=2){if(i)out<<',';long op=conv<long>(cmd[i]),index=conv<long>(cmd[i+1]);std::ostringstream value;const char*kind=nullptr;std::string message;
 try{
  const ZZX&r=raw.at(index);
  if(op==0){ZZX result;GCD(result,raw.at(conv<long>(coeff(r,0))),raw.at(conv<long>(coeff(r,1))));emit(value,result);}
  else if(op==1){vec_pair_ZZX_long result;SquareFreeDecomp(result,raw.at(conv<long>(coeff(r,0))));emit(value,result);}
  else if(op==2){ZZX result;mul(result,raw.at(conv<long>(coeff(r,0))),raw.at(conv<long>(coeff(r,1))));emit(value,result);}
  else throw std::runtime_error("unknown integer GCD state operation");
 }catch(const AdapterGuard&e){kind="RangeError";message=e.what();}catch(const std::exception&e){kind="Error";message=e.what();}
 out<<'[';if(kind)out<<std::quoted(kind)<<','<<std::quoted(message)<<",null";else out<<"null,null,"<<value.str();
 out<<",[";for(long j=0;j<FFTTables.length();j++){if(j)out<<',';info(out,*FFTTables[j]);}out<<']';
 RandomStream copy(stream);std::vector<unsigned char>tail(64);copy.get(tail.data(),64);out<<','<<std::quoted(hexbytes(tail))<<']';
 }
 out<<']';std::cout<<out.str()<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}
int main(){ErrorMsgCallback=error;SetNumThreads(1);std::string line;while(std::getline(std::cin,line)){
 // Each trace starts with the original global cache and static search state.
 // Keep the parent pristine; all operations in a trace share one child/cache.
 pid_t child=fork();if(child<0){std::cout<<"ERROR \"native integer GCD state fork failed\""<<std::endl;continue;}
 if(child==0){one(line);_exit(0);}int status=0;pid_t waited;do{waited=waitpid(child,&status,0);}while(waited<0&&errno==EINTR);
 if(waited<0||!WIFEXITED(status)||WEXITSTATUS(status)!=0){std::cout<<"ERROR \"native integer GCD state child failed\""<<std::endl;}
}}
