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
/*__NATIVE_CONTEXT_BODY__*/
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
 active_stream=&fallback;zz_pInfo=nullptr;
 std::istringstream in(line);long op,p,maxroot,d;vec_ZZ ff,gg,values,key;in>>op>>p>>maxroot>>d>>ff>>gg>>values>>key;if(!in||key.length()<32)throw std::runtime_error("invalid NTL word context input");
 zz_pInfoT native_context(p,maxroot);zz_pInfo=&native_context;zz_pX f=aspoly(ff),g=aspoly(gg);vec_zz_p v;v.SetLength(values.length());for(long i=0;i<values.length();i++)conv(v[i],values[i]);
 std::vector<unsigned char>kb(key.length());for(long i=0;i<key.length();i++)kb[i]=conv<unsigned long>(key[i]);RandomStream stream(kb.data());active_stream=&stream;
 std::ostringstream out,value;out<<"[["<<zz_pInfo->NumPrimes<<','<<zz_pInfo->PrimeCnt<<','<<zz_pInfo->MaxRoot<<','<<NTL_zz_pX_MOD_CROSSOVER<<','<<NTL_zz_pX_MUL_CROSSOVER<<','<<NTL_zz_pX_HalfGCD_CROSSOVER<<','<<NTL_zz_pX_GCD_CROSSOVER<<','<<NTL_zz_pX_BERMASS_CROSSOVER<<"],";
 try{
 if(op==0){zz_pXModulus F(f);value<<'['<<F.n<<','<<F.UseFFT<<']';}
 else if(op==1){zz_pX x;GCD(x,f,g);poly(value,x);}
 else if(op==2){zz_pX x;MinPolySeq(x,v,d);poly(value,x);}
 else if(op==3){zz_pX x;BuildFromRoots(x,v);poly(value,x);}
 else if(op==4){vec_zz_pX x;SFCanZass(x,f,0);polys(value,x);}
 else if(op==5){vec_zz_pX x;if(d==0&&IsOne(LeadCoeff(f)))throw std::runtime_error("EDF: degree must be nonzero");EDF(x,f,g,d,0);polys(value,x);}
 else if(op==6){zz_pX x;zz_pXModulus F(f);MinPolyMod(x,g,F,d);poly(value,x);}
 else if(op==7){vec_zz_pX x;FindFactors(x,f,g,v);polys(value,x);}
 else if(op==8||op==9){vec_pair_zz_pX_long groups;zz_pX h;if(op==8)SFCanZass1(groups,h,f,0);else{h=g;NewDDF(groups,f,h,0);}value<<"[[";for(long i=0;i<groups.length();i++){if(i)value<<',';value<<'[';poly(value,groups[i].a);value<<','<<groups[i].b<<']';}value<<"],";poly(value,h);value<<']';}
 else throw std::runtime_error("unknown NTL word context operation");
 out<<"null,null,"<<value.str();
 }catch(const std::exception&e){out<<"\"Error\","<<std::quoted(e.what())<<",null";}
 std::vector<unsigned char>tail(64);stream.get(tail.data(),64);out<<','<<std::quoted(hex(tail))<<']';std::cout<<out.str()<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
