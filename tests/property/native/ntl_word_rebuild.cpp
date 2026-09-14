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
#include "lzz_pX.cpp"
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
static void one(const std::string&line){try{
 std::istringstream in(line);vec_ZZ key,coeffs,cmd;vec_long lengths;long p,maxroot;in>>key>>p>>maxroot>>lengths>>coeffs>>cmd;
 if(!in||key.length()!=32||cmd.length()%2)throw std::runtime_error("invalid word quotient rebuild input");
 unsigned char bytes[32];for(long i=0;i<32;i++)bytes[i]=conv<unsigned long>(key[i]);SetSeed(RandomStream(bytes));RandomStream&stream=GetCurrentRandomStream();
 zz_p::init(p,maxroot);std::vector<zz_pX> polynomials;std::vector<ZZX> raw;long offset=0;
 for(long i=0;i<lengths.length();i++){zz_pX f;ZZX r;for(long j=0;j<lengths[i];j++){SetCoeff(r,j,coeffs[offset]);SetCoeff(f,j,conv<zz_p>(coeffs[offset++]));}polynomials.push_back(f);raw.push_back(r);}
 if(offset!=coeffs.length())throw std::runtime_error("invalid polynomial lengths");
 zz_pXModulus F;zz_pXMultiplier B;std::ostringstream out;out<<'[';
 for(long i=0;i<cmd.length();i+=2){if(i)out<<',';long op=conv<long>(cmd[i]),index=conv<long>(cmd[i+1]);std::ostringstream value;const char*kind=nullptr;std::string message;
 try{
  const zz_pX&a=polynomials.at(index);zz_pX result;
  if(op==0){build(F,a);value<<"null";}
  else if(op==1){rem(result,a,F);emit(value,result);}
  else if(op==2){MulMod(result,a,a,F);emit(value,result);}
  else if(op==3){SqrMod(result,a,F);emit(value,result);}
  else if(op==4){build(B,a,F);value<<"null";}
  else if(op==5){MulMod(result,a,B,F);emit(value,result);}
  else if(op==6){value<<'[';emit(value,B.b);value<<',';scalar(value,B.UseFFT);value<<']';}
  else if(op==7||op==8){
   const ZZX&r=raw.at(index);long prime=conv<long>(coeff(r,0)),k=conv<long>(coeff(r,1)),yn=conv<long>(coeff(r,2)),xn=conv<long>(coeff(r,3));
   UseFFTPrime(prime);const auto&info=*FFTTables[prime];
   if(k<0||k>CalcMaxRoot(info.q))throw AdapterGuard("FFT transform: exponent exceeds the prime root table");
   long N=1L<<k;for(long length:{yn,op==7?xn:yn})if(length<1||length>N||FFTRoundUp(length,k)!=length)throw AdapterGuard("FFT transform: lengths must be admissible");
   if(lengths[index]-4<(op==7?xn:yn))throw AdapterGuard("FFT transform: input is shorter than its declared length");
   std::vector<long> input(N),output(N);for(long j=0;j<(op==7?xn:yn);j++)input[j]=rem(coeff(r,j+4),info.q);
   if(op==7)FFTFwd_trunc(output.data(),input.data(),k,info,yn,xn);else FFTRev1_trunc(output.data(),input.data(),k,info,yn);
   value<<'[';for(long j=0;j<yn;j++){if(j)value<<',';scalar(value,output[j]);}value<<']';
  }
  else if(op==9){
   const ZZX&r=raw.at(index);ZZ xn=coeff(r,0),k=coeff(r,1),limit=power2_ZZ(63);
   if(k<0||k>62)throw AdapterGuard("FFTRoundUp: exponent must be between 0 and 62");
   if(xn< -limit||xn>=limit)throw AdapterGuard("FFTRoundUp: length must be a signed native integer");
   if(xn>=limit-15)throw AdapterGuard("FFTRoundUp: length rounding overflows");
   scalar(value,FFTRoundUp(conv<long>(xn),conv<long>(k)));
  }
  else if(op==10){F=zz_pXModulus();value<<"null";}
  else if(op==11){
   // PlainUpdateMap writes x[n-deg(b)] first. Negative indices are undefined
   // in the bundled unchecked-vector build; degree exactly n remains defined.
   if(!B.UseFFT && deg(B.b)>F.n && a.rep.length()<=F.n)
    throw AdapterGuard("UpdateMap: multiplier degree exceeds native buffer bounds");
   vec_zz_p x;UpdateMap(x,a.rep,B,F);
   value<<'[';for(long j=0;j<x.length();j++){if(j)value<<',';scalar(value,rep(x[j]));}value<<']';
  }
  else if(op==13||op==14){
   const ZZX&r=raw.at(index);long prime=conv<long>(coeff(r,0)),k=conv<long>(coeff(r,1));
   UseFFTPrime(prime);const auto&info=*FFTTables[prime];
   if(k<0||k>CalcMaxRoot(info.q))throw AdapterGuard("FFT transform: exponent exceeds the prime root table");
   long N=1L<<k;
   if(lengths[index]-4<N)throw AdapterGuard("FFT transform: input is shorter than its declared length");
   std::vector<long> input(N),output(N);for(long j=0;j<N;j++)input[j]=rem(coeff(r,j+4),info.q);
   if(op==13)FFTFwd_trans(output.data(),input.data(),k,info);else FFTRev1_trans(output.data(),input.data(),k,info);
   value<<'[';for(long j=0;j<N;j++){if(j)value<<',';scalar(value,output[j]);}value<<']';
  }
  else if(op==12){
   const ZZX&r=raw.at(index);long rows=conv<long>(coeff(r,0)),cols=conv<long>(coeff(r,1));
   if(rows!=zz_pInfo->NumPrimes)throw AdapterGuard("FromModularRep: incorrect number of prime rows");
   if(cols<0||lengths[index]!=2+rows*cols)throw std::runtime_error("invalid CRT matrix encoding");
   fftRep R;R.SetSize(NextPowerOfTwo(std::max(1L,cols)));R.len=1L<<R.k;
   for(long i=0;i<rows;i++)for(long j=0;j<cols;j++)R.tbl[i][j]=rem(coeff(r,2+i*cols+j),GetFFTPrime(i));
   vec_zz_p out;out.SetLength(cols);FromModularRep(out.elts(),R,0,cols,zz_pInfo);
   value<<'[';for(long j=0;j<cols;j++){if(j)value<<',';scalar(value,rep(out[j]));}value<<']';
  }
  else if(op==15){
   const ZZX&r=raw.at(index);long p=conv<long>(coeff(r,0)),maxroot=conv<long>(coeff(r,1)),length=conv<long>(coeff(r,2));
   zz_pPush restore(p,maxroot);zz_pX g,b;
   for(long j=0;j<length;j++)SetCoeff(g,j,conv<zz_p>(coeff(r,3+j)));
   for(long j=3+length;j<lengths[index];j++)SetCoeff(b,j-3-length,conv<zz_p>(coeff(r,j)));
   zz_pXModulus other(g);build(B,b,other);value<<"null";
  }
  else if(op>=16 && op<=25){
   const ZZX&r=raw.at(index);long length=conv<long>(coeff(r,0));zz_pX aa,bb;
   for(long j=0;j<length;j++)SetCoeff(aa,j,conv<zz_p>(coeff(r,j+1)));
   for(long j=1+length;j<lengths[index];j++)SetCoeff(bb,j-1-length,conv<zz_p>(coeff(r,j)));
   if(op==16){zz_pX d,s,t;XGCD(d,s,t,aa,bb);value<<'[';emit(value,d);value<<',';emit(value,s);value<<',';emit(value,t);value<<']';}
   else if(op==17){mul(result,aa,bb);emit(value,result);}
   else if(op==18){sqr(result,aa);emit(value,result);}
   else if(op==19){mul(result,aa,aa);emit(value,result);}
   else if(op==20){InvMod(result,aa,bb);emit(value,result);}
   else if(op==21){long status=InvModStatus(result,aa,bb);value<<'['<<status<<',';emit(value,result);value<<']';}
   else if(op==22){PowerXMod(result,conv<ZZ>(coeff(r,1)),F);emit(value,result);}
   else if(op==23){PowerXPlusAMod(result,conv<zz_p>(coeff(r,1)),conv<ZZ>(coeff(r,2)),F);emit(value,result);}
   else if(op==24){PowerMod(result,aa,conv<ZZ>(coeff(r,1+length)),F);emit(value,result);}
   else {GCD(result,aa,bb);emit(value,result);}
  }
  else throw std::runtime_error("unknown quotient rebuild operation");
 }catch(const AdapterGuard&e){kind="RangeError";message=e.what();}catch(const std::exception&e){kind="Error";message=e.what();}
 out<<'[';if(kind)out<<std::quoted(kind)<<','<<std::quoted(message)<<",null";else out<<"null,null,"<<value.str();
 out<<',';scalar(out,F.n);out<<',';emit(out,F.f);
 out<<",[";for(long j=0;j<FFTTables.length();j++){if(j)out<<',';info(out,*FFTTables[j]);}out<<']';
 RandomStream copy(stream);std::vector<unsigned char>tail(64);copy.get(tail.data(),64);out<<','<<std::quoted(hexbytes(tail))<<']';
 }
 out<<']';std::cout<<out.str()<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}
int main(){ErrorMsgCallback=error;SetNumThreads(1);std::string line;while(std::getline(std::cin,line)){
 // Each trace starts with the original global cache and static search state.
 // Keep the parent pristine; all operations in a trace share one child/cache.
 pid_t child=fork();if(child<0){std::cout<<"ERROR \"native word quotient rebuild fork failed\""<<std::endl;continue;}
 if(child==0){one(line);_exit(0);}int status=0;pid_t waited;do{waited=waitpid(child,&status,0);}while(waited<0&&errno==EINTR);
 if(waited<0||!WIFEXITED(status)||WEXITSTATUS(status)!=0){std::cout<<"ERROR \"native word quotient rebuild child failed\""<<std::endl;}
}}
