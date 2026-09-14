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
#include <memory>
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

static void poly(std::ostream&out,const zz_pX&x){out<<'[';for(long i=0;i<x.rep.length();i++){if(i)out<<',';out<<'"'<<rep(x.rep[i])<<'"';}out<<']';}
static void polys(std::ostream&out,const vec_zz_pX&v){out<<'[';for(long i=0;i<v.length();i++){if(i)out<<',';poly(out,v[i]);}out<<']';}
static zz_pX aspoly(const vec_ZZ&v){zz_pX f;for(long i=0;i<v.length();i++)SetCoeff(f,i,conv<zz_p>(v[i]));return f;}

static void zpoly(std::ostream&out,const ZZX&f){out<<'[';for(long i=0;i<f.rep.length();i++){if(i)out<<',';out<<'"'<<f.rep[i]<<'"';}out<<']';}
static void zpolys(std::ostream&out,const vec_ZZX&v){out<<'[';for(long i=0;i<v.length();i++){if(i)out<<',';zpoly(out,v[i]);}out<<']';}
static void bigpolys(std::ostream&out,const vec_ZZ_pX&v){out<<'[';for(long i=0;i<v.length();i++){if(i)out<<',';out<<'[';for(long j=0;j<v[i].rep.length();j++){if(j)out<<',';out<<'"'<<rep(v[i].rep[j])<<'"';}out<<']';}out<<']';}
static void init_info(LocalInfoT& I,long capacity){I.n=-1;I.NumPrimes=0;I.NumFactors=0;clear(I.PossibleDegrees);I.p.SetLength(capacity);for(long i=0;i<capacity;i++)I.p[i]=-17;I.p.SetLength(0);zz_pContext().restore();}
static void snapshot(std::ostream& out,const LocalInfoT&I){out<<'['<<I.n<<','<<I.NumPrimes<<','<<I.NumFactors<<",\""<<I.PossibleDegrees<<"\",[";for(long i=0;i<I.p.length();i++){if(i)out<<',';out<<'"'<<I.p[i]<<'"';}out<<"],[";for(long i=0;i<I.pattern.length();i++){if(i)out<<',';out<<'[';for(long j=0;j<I.pattern[i].length();j++){if(j)out<<',';out<<I.pattern[i][j];}out<<']';}out<<"],";if(zz_pInfo)out<<"[\""<<zz_p::modulus()<<"\","<<zz_pInfo->PrimeCnt<<','<<zz_pInfo->MaxRoot<<']';else out<<"null";out<<']';}
static void one(const std::string& line){try{
 std::istringstream in(line);vec_ZZ commands,key;long capacity,warm;in>>commands>>key>>capacity>>warm;if(!in||key.length()<32||capacity<0||capacity>10000)throw std::runtime_error("invalid NTL small-prime input");
 std::vector<unsigned char>kb(key.length());for(long i=0;i<key.length();i++)kb[i]=conv<unsigned long>(key[i]);RandomStream initial(kb.data());SetSeed(initial);if(warm){unsigned char zero[32]={};SetSeed(RandomStream(zero));for(long i=0;i<warm;i++)UseFFTPrime(i);SetSeed(initial);}RandomStream& stream=GetCurrentRandomStream();
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
   else if(op==12||op==13||op==14||op==15){
     long vh=1,maxp=50;if(op==13){vh=word();maxp=word();}else if(op==14)vh=word();else if(op==15)maxp=word();
     ZZ modulus=scalar();long k=word();ZZ_p::init(modulus);
     long nw=word();vec_ZZ_pX W;W.SetLength(nw);for(long i=0;i<nw;i++){vec_ZZ v=vector();ZZX t;t.rep=v;t.normalize();conv(W[i],t);}
     long nf=word();vec_ZZX factors;factors.SetLength(nf);for(long i=0;i<nf;i++){factors[i].rep=vector();factors[i].normalize();}
     ZZX f;f.rep=vector();f.normalize();ZZXFac_van_Hoeij=vh;ZZXFac_MaxNumPrimes=maxp;
     bool changed=I->NumFactors<factors.length();vec_ZZ pdeg;UpdateLocalInfo(*I,pdeg,W,factors,f,k,0);
     if(changed){value<<'[';for(long i=0;i<pdeg.length();i++){if(i)value<<',';value<<'"'<<pdeg[i]<<'"';}value<<']';}else value<<"null";
   }
   else if(op==16||op==17||op==18){
     ZZ modulus=scalar();long k=word(),bound=word(),vh=word(),maxp=word(),pruning=word();ZZ_p::init(modulus);
     long nf=word();vec_ZZX factors;factors.SetLength(nf);for(long i=0;i<nf;i++){factors[i].rep=vector();factors[i].normalize();}
     ZZX f;f.rep=vector();f.normalize();long nw=word();vec_ZZX w;vec_ZZ_pX W;w.SetLength(nw);W.SetLength(nw);for(long i=0;i<nw;i++){w[i].rep=vector();w[i].normalize();conv(W[i],w[i]);}
     ZZXFac_van_Hoeij=vh;ZZXFac_MaxNumPrimes=maxp;ZZXFac_MaxPrune=pruning;
     // Explicit guards for undefined native indexing; these are labeled in shared cases.
     if(op==16&&k==0&&nw==0)throw std::runtime_error("CardinalitySearch: zero cardinality with no factors has no defined native result");
     if(op==17&&k>=2&&NumBits(k)<=30&&nw==0)throw std::runtime_error("CardinalitySearch1: empty factor vector has no defined native result");
     if(op==16)CardinalitySearch(factors,f,W,*I,k,bound,0);
     else if(op==17)CardinalitySearch1(factors,f,W,*I,k,bound,0);
     else{FindTrueFactors(factors,f,w,modulus,*I,0,bound);zpolys(value,factors);}
     if(op!=18){value<<'[';zpolys(value,factors);value<<',';zpoly(value,f);value<<',';bigpolys(value,W);value<<']';}
   }
   else if(op==20){UseFFTPrime(word());value<<"null";}
   else if(op==21){stream.set_nonce(conv<unsigned long>(scalar()));value<<"null";}
   else if(op==22){long n=word();std::vector<unsigned char> bytes(n);stream.get(bytes.data(),n);value<<std::quoted(hexbytes(bytes));}
   else if(op==23){ZZ p=scalar();vec_ZZ indices=vector();long nw=word();ZZ_p::init(p);vec_long selected;selected.SetLength(indices.length());for(long j=0;j<indices.length();j++)selected[j]=conv<long>(indices[j]);vec_ZZ_pX W;W.SetLength(nw);for(long j=0;j<nw;j++){ZZX t;t.rep=vector();t.normalize();conv(W[j],t);}ZZ_pX result;InvMul(result,W,selected);ZZX z;conv(z,result);zpoly(value,z);}
   else throw std::runtime_error("unknown small-prime operation");
   out<<"null,null,"<<value.str();
  }catch(const std::exception&e){out<<"\"Error\","<<std::quoted(e.what())<<",null";}
  out<<',';snapshot(out,*I);RandomStream copy(stream);std::vector<unsigned char>tail(64);copy.get(tail.data(),64);out<<','<<std::quoted(hexbytes(tail))<<",[";for(long j=0;j<FFTTables.length();j++){if(j)out<<',';info(out,*FFTTables[j]);}out<<"]]";
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
