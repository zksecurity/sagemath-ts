#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/SmartPtr.h>
#include <NTL/lzz_pX.h>
#include <NTL/BasicThreadPool.h>
#include <sstream>
#include <iomanip>
#include <memory>
#include <vector>
#include <cstring>
#if defined(NTL_RANDOMSTREAM_AES256CTR) || defined(NTL_HAVE_AVX2) || defined(NTL_HAVE_SSSE3) || NTL_BITS_PER_LONG != 64 || NTL_BITS_PER_INT32 != 32
#error Unsupported native byte stream profile
#endif
#include "__BUNDLED_FACTOR_SOURCE__"
#include "__BUNDLED_INTEGER_FACTOR_SOURCE__"
namespace NTL {
/*__NATIVE_CONTEXT_BODY__*/
/*__NATIVE_PRIME_SEQUENCE__*/
/*__NATIVE_RANDOM_BODY__*/
static RandomStream *active_stream;
RandomStream& GetCurrentRandomStream() { return *active_stream; }
static inline RandomStream& LocalGetCurrentRandomStream() { return *active_stream; }
/*__NATIVE_SAMPLING_BODY__*/
}
using namespace NTL;
static void error(const char*x){throw std::runtime_error(x);}
static std::string hex(const std::vector<unsigned char>&b){static const char*h="0123456789abcdef";std::string s;for(auto x:b){s+=h[x>>4];s+=h[x&15];}return s;}
static void poly(std::ostream&out,const zz_pX&x){out<<'[';for(long i=0;i<x.rep.length();i++){if(i)out<<',';out<<'"'<<rep(x.rep[i])<<'"';}out<<']';}
static void polys(std::ostream&out,const vec_zz_pX&v){out<<'[';for(long i=0;i<v.length();i++){if(i)out<<',';poly(out,v[i]);}out<<']';}
static zz_pX aspoly(const vec_ZZ&v){zz_pX f;for(long i=0;i<v.length();i++)SetCoeff(f,i,conv<zz_p>(v[i]));return f;}

static void init_info(LocalInfoT& I,long capacity){I.n=-1;I.NumPrimes=0;I.NumFactors=0;clear(I.PossibleDegrees);I.p.SetLength(capacity);for(long i=0;i<capacity;i++)I.p[i]=-17;I.p.SetLength(0);zz_pContext().restore();}
static void snapshot(std::ostream& out,const LocalInfoT&I){out<<'['<<I.n<<','<<I.NumPrimes<<','<<I.NumFactors<<",\""<<I.PossibleDegrees<<"\",[";for(long i=0;i<I.p.length();i++){if(i)out<<',';out<<'"'<<I.p[i]<<'"';}out<<"],[";for(long i=0;i<I.pattern.length();i++){if(i)out<<',';out<<'[';for(long j=0;j<I.pattern[i].length();j++){if(j)out<<',';out<<I.pattern[i][j];}out<<']';}out<<"],";if(zz_pInfo)out<<"[\""<<zz_p::modulus()<<"\","<<zz_pInfo->PrimeCnt<<','<<zz_pInfo->MaxRoot<<']';else out<<"null";out<<']';}
int main(){ErrorMsgCallback=error;SetNumThreads(1);unsigned char empty_key[32]={0};RandomStream fallback(empty_key);active_stream=&fallback;for(long i=0;i<4;i++)UseFFTPrime(i);std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);vec_ZZ commands,key;long capacity;in>>commands>>key>>capacity;if(!in||key.length()<32||capacity<0||capacity>10000)throw std::runtime_error("invalid NTL small-prime input");
 std::vector<unsigned char>kb(key.length());for(long i=0;i<key.length();i++)kb[i]=conv<unsigned long>(key[i]);RandomStream stream(kb.data());active_stream=&stream;
 auto I=std::make_unique<LocalInfoT>();init_info(*I,capacity);long pos=0;bool first=true;std::ostringstream out;out<<'[';
 auto scalar=[&]()->ZZ {if(pos>=commands.length())throw std::runtime_error("truncated small-prime command");return commands[pos++];};
 auto word=[&]()->long{return conv<long>(scalar());};
 auto vector=[&]()->vec_ZZ {long n=word();if(n<0||n>commands.length()-pos)throw std::runtime_error("invalid small-prime vector");vec_ZZ v;v.SetLength(n);for(long i=0;i<n;i++)v[i]=scalar();return v;};
 while(pos<commands.length()){
  long op=word();if(!first)out<<',';first=false;out<<'[';std::ostringstream value;
  try{
   if(op==0||op==9||op==10||op==11){long initial=7,maximum=50;if(op==0){initial=word();maximum=word();}else if(op==10)initial=word();else if(op==11)maximum=word();vec_ZZ v=vector();if(initial>=1&&initial<=10000&&maximum>=initial&&maximum<=10000&&initial>capacity)throw std::runtime_error("small-prime oracle requires initialized capacity");ZZX f;f.rep=v;f.normalize();ZZXFac_InitNumPrimes=initial;ZZXFac_MaxNumPrimes=maximum;std::unique_ptr<vec_zz_pX> W(SmallPrimeFactorization(*I,f,0));if(W)polys(value,*W);else value<<"null";}
   else if(op==1){I->s.reset(word());value<<"null";}
   else if(op==2){I->PossibleDegrees=scalar();value<<"null";}
   else if(op==3){vec_ZZ v=vector();vec_long w;w.SetLength(v.length());for(long i=0;i<v.length();i++)w[i]=conv<long>(v[i]);I->p=w;value<<"null";}
   else if(op==4){long n=word();vec_vec_long patterns;patterns.SetLength(n);for(long i=0;i<n;i++){vec_ZZ v=vector();patterns[i].SetLength(v.length());for(long j=0;j<v.length();j++)patterns[i][j]=conv<long>(v[j]);}I->pattern=patterns;value<<"null";}
   else if(op==5){I=std::make_unique<LocalInfoT>();init_info(*I,capacity);value<<"null";}
   else if(op==6){long n=word();value<<'[';for(long i=0;i<n;i++){if(i)value<<',';value<<'"'<<I->s.next()<<'"';}value<<']';}
   else if(op==7){long p=word(),maxroot=word();zz_p::init(p,maxroot);value<<"null";}
   else if(op==8){I->n=word();I->NumPrimes=word();I->NumFactors=word();value<<"null";}
   else throw std::runtime_error("unknown small-prime operation");
   out<<"null,null,"<<value.str();
  }catch(const std::exception&e){out<<"\"Error\","<<std::quoted(e.what())<<",null";}
  out<<',';snapshot(out,*I);RandomStream copy(stream);std::vector<unsigned char>tail(64);copy.get(tail.data(),64);out<<','<<std::quoted(hex(tail))<<']';
 }
 out<<']';std::cout<<out.str()<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
