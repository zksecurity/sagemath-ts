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
namespace NTL {
void DoMinPolyMod(zz_pX&, const zz_pX&, const zz_pXModulus&, long, const vec_zz_p&);
/*__NATIVE_RANDOM_BODY__*/
static RandomStream *active_stream;
static inline RandomStream& LocalGetCurrentRandomStream() { return *active_stream; }
/*__NATIVE_SAMPLING_BODY__*/
}
using namespace NTL;
static void error(const char*x){throw std::runtime_error(x);}
static std::string hex(const std::vector<unsigned char>&b){static const char*h="0123456789abcdef";std::string s;for(auto x:b){s+=h[x>>4];s+=h[x&15];}return s;}
static void poly(std::ostream&out,const zz_pX&x){out<<'[';for(long i=0;i<x.rep.length();i++){if(i)out<<',';out<<'"'<<rep(x.rep[i])<<'"';}out<<']';}
static zz_pX aspoly(const vec_ZZ&v){zz_pX f;for(long i=0;i<v.length();i++)SetCoeff(f,i,conv<zz_p>(v[i]));return f;}
int main(){ErrorMsgCallback=error;SetNumThreads(1);unsigned char empty_key[32]={0};RandomStream fallback(empty_key);std::string line;while(std::getline(std::cin,line)){try{
 active_stream=&fallback;
 std::istringstream in(line);long op,p,m,ready;vec_ZZ f,g,r,key;in>>op>>p>>m>>ready>>f>>g>>r>>key;if(!in||key.length()<32)throw std::runtime_error("invalid NTL element minimum polynomial input");
 zz_p::init(p);zz_pXModulus F;if(ready)build(F,aspoly(f));zz_pX G=aspoly(g),h;vec_zz_p R;R.SetLength(r.length());for(long i=0;i<r.length();i++)conv(R[i],r[i]);
 std::vector<unsigned char>bytes(key.length());for(long i=0;i<key.length();i++)bytes[i]=conv<unsigned long>(key[i]);RandomStream stream(bytes.data());active_stream=&stream;
 std::ostringstream out;out<<'[';
 try{
 if(op==0)DoMinPolyMod(h,G,F,m,R);else if(op==1)IrredPolyMod(h,G,F,m);else if(op==2)ProbMinPolyMod(h,G,F,m);else if(op==3)MinPolyMod(h,G,F,m);else if(op==4)IrredPolyMod(h,G,F);else if(op==5)ProbMinPolyMod(h,G,F);else if(op==6)MinPolyMod(h,G,F);else throw std::runtime_error("unknown element minimum polynomial operation");
 out<<"null,null,";poly(out,h);
 }catch(const std::exception&e){out<<"\"Error\","<<std::quoted(e.what())<<",null";}
 std::vector<unsigned char>tail(64);stream.get(tail.data(),64);out<<','<<std::quoted(hex(tail))<<']';std::cout<<out.str()<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
