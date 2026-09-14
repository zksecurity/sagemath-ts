#include <NTL/FFT.h>
#include <NTL/FFT_impl.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_lZZ_p.h>
#include <NTL/BasicThreadPool.h>
#include <sys/wait.h>
#include <unistd.h>
#include <cerrno>
#include <cstdint>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <NTL/ZZ_pX.h>
#if !defined(__aarch64__)
#error This fixture profile requires AArch64 fused CRT arithmetic
#endif
#include <NTL/ZZX.h>
#if defined(NTL_THREADS) || defined(NTL_RANDOM_AES256CTR) || defined(NTL_HAVE_AVX2) || defined(NTL_HAVE_SSSE3)
#error Unsupported native factor state profile
#endif
static_assert(NTL_BITS_PER_LONG==64 && NTL_SP_NBITS==60 && NTL_ZZ_NBITS==64 && NTL_BITS_PER_INT32==32 && NTL_BITS_PER_NONCE==64,"unsupported native word profile");
static_assert(NTL_FDOUBLE_PRECISION==double(1ULL<<52),"unsupported probability profile");
namespace NTL { void FromModularRep(ZZ_p&,vec_long&,const ZZ_pFFTInfoT*,ZZ_pTmpSpaceT*); }
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
static void emit(std::ostream&out,const ZZ_pX&f){out<<'[';for(long i=0;i<=deg(f);i++){if(i)out<<',';emit(out,rep(f[i]));}out<<']';}
static void one(const std::string&line){try{
 std::istringstream in(line);vec_ZZ key,coeffs,cmd;vec_long lengths;ZZ p;in>>key>>p>>lengths>>coeffs>>cmd;
 if(!in||key.length()!=32||cmd.length()%2)throw std::runtime_error("invalid big quotient rebuild input");
 unsigned char bytes[32];for(long i=0;i<32;i++)bytes[i]=conv<unsigned long>(key[i]);SetSeed(RandomStream(bytes));RandomStream&stream=GetCurrentRandomStream();
 ZZ_p::init(p);std::vector<ZZ_pX> polynomials;std::vector<ZZX> raw;long offset=0;
 for(long i=0;i<lengths.length();i++){ZZ_pX f;ZZX r;for(long j=0;j<lengths[i];j++){SetCoeff(r,j,coeffs[offset]);SetCoeff(f,j,conv<ZZ_p>(coeffs[offset++]));}polynomials.push_back(f);raw.push_back(r);}
 if(offset!=coeffs.length())throw std::runtime_error("invalid polynomial lengths");
 ZZ_pXModulus F;ZZ_pXMultiplier B;std::ostringstream out;out<<'[';
 for(long i=0;i<cmd.length();i+=2){if(i)out<<',';long op=conv<long>(cmd[i]),index=conv<long>(cmd[i+1]);std::ostringstream value;const char*kind=nullptr;std::string message;
 try{
  const ZZ_pX&a=polynomials.at(index);ZZ_pX result;
  if(op==0){build(F,a);value<<"null";}
  else if(op==1){rem(result,a,F);emit(value,result);}
  else if(op==2){MulMod(result,a,a,F);emit(value,result);}
  else if(op==3){SqrMod(result,a,F);emit(value,result);}
  else if(op==4){build(B,a,F);value<<"null";}
  else if(op==5){MulMod(result,a,B,F);emit(value,result);}
  else if(op==6){value<<'[';emit(value,B.b);value<<',';scalar(value,B.UseFFT);value<<']';}
  else if(op==7){
   const ZZX&r=raw.at(index);long rows=conv<long>(coeff(r,0)),cols=conv<long>(coeff(r,1));
   const auto*fft=ZZ_p::GetFFTInfo();auto*tmp=ZZ_p::GetTmpSpace();
   if(rows!=fft->NumPrimes)throw AdapterGuard("FromModularRep: incorrect number of prime rows");
   if(cols<0||lengths[index]!=2+rows*cols)throw std::runtime_error("invalid CRT matrix encoding");
   value<<'[';
   for(long j=0;j<cols;j++){
    vec_long residues;residues.SetLength(rows);
    for(long i=0;i<rows;i++)residues[i]=rem(coeff(r,2+i*cols+j),GetFFTPrime(i));
    ZZ_p x;FromModularRep(x,residues,fft,tmp);
    if(j)value<<',';emit(value,rep(x));
   }value<<']';
  }
  else if(op==8){F=ZZ_pXModulus();value<<"null";}
  else if(op==9){emit(value,F.val());}
  else if(op==10){
   const ZZX&r=raw.at(index);ZZ_pPush restore(coeff(r,0));ZZ_pX g;
   for(long j=1;j<lengths[index];j++)SetCoeff(g,j-1,conv<ZZ_p>(coeff(r,j)));
   ZZ_pXModulus other(g);value<<"null";
  }
  else if(op==11){B=ZZ_pXMultiplier(a,F);value<<"null";}
  else if(op==12){F=ZZ_pXModulus(a);value<<"null";}
  else if(op==13){
   const ZZX&r=raw.at(index);ZZ_pPush restore(coeff(r,0));long length=conv<long>(coeff(r,1));ZZ_pX g,b;
   for(long j=0;j<length;j++)SetCoeff(g,j,conv<ZZ_p>(coeff(r,2+j)));
   for(long j=2+length;j<lengths[index];j++)SetCoeff(b,j-2-length,conv<ZZ_p>(coeff(r,j)));
   ZZ_pXModulus other(g);build(B,b,other);value<<"null";
  }
  else if(op==15){
   const ZZX&r=raw.at(index);long length=conv<long>(coeff(r,0));ZZ_pX a,b;
   for(long j=0;j<length;j++)SetCoeff(a,j,conv<ZZ_p>(coeff(r,1+j)));
   for(long j=1+length;j<lengths[index];j++)SetCoeff(b,j-1-length,conv<ZZ_p>(coeff(r,j)));
   MulMod(result,a,b,F);emit(value,result);
  }
  else if(op==16){ZZ_pX b=a;MulMod(result,a,b,F);emit(value,result);}
  else if(op==17){
   const ZZX&r=raw.at(index);long rows=conv<long>(coeff(r,0));const auto*fft=ZZ_p::GetFFTInfo();auto*tmp=ZZ_p::GetTmpSpace();
   if(rows!=fft->NumPrimes)throw AdapterGuard("FromModularRep: incorrect number of prime rows");
   long cols=conv<long>(coeff(r,1));
   for(long i=0;i<rows;i++)if(coeff(r,1+i)!=cols)throw AdapterGuard("FromModularRep: inconsistent coefficient counts");
   value<<'[';
   for(long j=0;j<cols;j++){vec_long residues;residues.SetLength(rows);for(long i=0;i<rows;i++)residues[i]=rem(coeff(r,1+rows+i*cols+j),GetFFTPrime(i));ZZ_p x;FromModularRep(x,residues,fft,tmp);if(j)value<<',';emit(value,rep(x));}value<<']';
  }
  else throw std::runtime_error("unknown quotient rebuild operation");
 }catch(const AdapterGuard&e){kind="RangeError";message=e.what();}catch(const std::exception&e){kind="Error";message=e.what();}
 out<<'[';if(kind)out<<std::quoted(kind)<<','<<std::quoted(message)<<",null";else out<<"null,null,"<<value.str();
 out<<',';scalar(out,F.n);out<<',';emit(out,F.f);out<<',';scalar(out,F.UseFFT);
 out<<",[";for(long j=0;j<FFTTables.length();j++){if(j)out<<',';info(out,*FFTTables[j]);}out<<']';
 RandomStream copy(stream);std::vector<unsigned char>tail(64);copy.get(tail.data(),64);out<<','<<std::quoted(hexbytes(tail))<<']';
 }
 out<<']';std::cout<<out.str()<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}
int main(){ErrorMsgCallback=error;SetNumThreads(1);std::string line;while(std::getline(std::cin,line)){
 // Each trace starts with the original global cache and static search state.
 // Keep the parent pristine; all operations in a trace share one child/cache.
 pid_t child=fork();if(child<0){std::cout<<"ERROR \"native big quotient rebuild fork failed\""<<std::endl;continue;}
 if(child==0){one(line);_exit(0);}int status=0;pid_t waited;do{waited=waitpid(child,&status,0);}while(waited<0&&errno==EINTR);
 if(waited<0||!WIFEXITED(status)||WEXITSTATUS(status)!=0){std::cout<<"ERROR \"native big quotient rebuild child failed\""<<std::endl;}
}}
