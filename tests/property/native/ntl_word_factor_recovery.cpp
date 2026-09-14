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
namespace NTL {
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
int main(){ErrorMsgCallback=error;SetNumThreads(1);unsigned char empty_key[32]={0};RandomStream fallback(empty_key);std::string line;while(std::getline(std::cin,line)){try{
 active_stream=&fallback;
 std::istringstream in(line);long op,p,d;vec_ZZ ff,bb,roots,packed,key;in>>op>>p>>d>>ff>>bb>>roots>>packed>>key;if(!in||key.length()<32)throw std::runtime_error("invalid NTL factor recovery input");
 zz_p::init(p);zz_pX f=aspoly(ff),b=aspoly(bb);vec_zz_p r;r.SetLength(roots.length());for(long i=0;i<roots.length();i++)conv(r[i],roots[i]);
 vec_pair_zz_pX_long groups;for(long i=0;i<packed.length();){if(i+2>packed.length())throw std::runtime_error("invalid packed factor groups");long n=conv<long>(packed[i++]),degree=conv<long>(packed[i++]);if(n<0||n>packed.length()-i)throw std::runtime_error("invalid packed factor group length");zz_pX g;for(long j=0;j<n;j++)SetCoeff(g,j,conv<zz_p>(packed[i++]));append(groups,cons(g,degree));}
 std::vector<unsigned char>kb(key.length());for(long i=0;i<key.length();i++)kb[i]=conv<unsigned long>(key[i]);RandomStream stream(kb.data());active_stream=&stream;
 std::ostringstream out,value;out<<'[';
 try{
 if(op==0){vec_zz_p found;FindRoots(found,f);value<<'[';for(long i=0;i<found.length();i++){if(i)value<<',';value<<'"'<<rep(found[i])<<'"';}value<<']';}
 else if(op==1){zz_p root;FindRoot(root,f);value<<'"'<<rep(root)<<'"';}
 else if(op==8){zz_pX q,remainder;DivRem(q,remainder,f,b);value<<'[';poly(value,q);value<<',';poly(value,remainder);value<<']';}
 else{vec_zz_pX factors;
 if(op==2)RootEDF(factors,f,0);
 else if(op==3)FindFactors(factors,f,b,r);
 else if(op==4){if(d==0){zz_pXModulus F;build(F,f);throw std::runtime_error("EDFSplit: degree must be nonzero");}EDFSplit(factors,f,b,d);}
 else if(op==5){if(d==0&&IsOne(LeadCoeff(f)))throw std::runtime_error("EDF: degree must be nonzero");EDF(factors,f,b,d,0);}
 else if(op==6){long zero=-1;for(long i=0;i<groups.length();i++)if(groups[i].b==0){zero=i;break;}if(zero>=0){vec_pair_zz_pX_long prefix;prefix.SetLength(zero);for(long i=0;i<zero;i++)prefix[i]=groups[i];SFCanZass2(factors,prefix,b,0);throw std::runtime_error("SFCanZass2: factor degree must be nonzero");}SFCanZass2(factors,groups,b,0);}
 else if(op==7)SFCanZass(factors,f,0);
 else throw std::runtime_error("unknown factor recovery operation");polys(value,factors);}
 out<<"null,null,"<<value.str();
 }catch(const std::exception&e){out<<"\"Error\","<<std::quoted(e.what())<<",null";}
 std::vector<unsigned char>tail(64);stream.get(tail.data(),64);out<<','<<std::quoted(hex(tail))<<']';std::cout<<out.str()<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
