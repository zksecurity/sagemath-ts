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
/*__NATIVE_RANDOM_BODY__*/
static RandomStream *active_stream;
RandomStream& GetCurrentRandomStream() { return *active_stream; }
static inline RandomStream& LocalGetCurrentRandomStream() { return *active_stream; }
/*__NATIVE_SAMPLING_BODY__*/
}
using namespace NTL;
static void error(const char*x){throw std::runtime_error(x);}
struct WordBoundError:std::runtime_error {using std::runtime_error::runtime_error;};
static long word(const ZZ&x){if(x<-(ZZ(1)<<63)||x>=(ZZ(1)<<63))throw WordBoundError("RandomBnd: word bound must fit a signed 64-bit integer");return conv<long>(x);}
static std::string hex(const std::vector<unsigned char>&b){static const char*h="0123456789abcdef";std::string s;for(auto x:b){s+=h[x>>4];s+=h[x&15];}return s;}
static void poly(std::ostream&out,const zz_pX&x){out<<'[';for(long i=0;i<x.rep.length();i++){if(i)out<<',';out<<'"'<<rep(x.rep[i])<<'"';}out<<']';}
static void state(std::ostream&out,const std::shared_ptr<RandomBndGenerator>(&g)[2]){out<<'[';for(long i=0;i<2;i++){if(i)out<<',';out<<"[\""<<g[i]->p<<"\",";if(g[i]->p)out<<g[i]->nb<<",\""<<g[i]->mask<<'"';else out<<"null,null";out<<']';}out<<']';}
int main(){ErrorMsgCallback=error;SetNumThreads(1);unsigned char empty_key[32]={0};RandomStream fallback(empty_key);std::string line;while(std::getline(std::cin,line)){try{
 active_stream=&fallback;
 std::istringstream in(line);long op,n;ZZ p;vec_ZZ data,key;in>>op>>p>>n>>data>>key;if(!in||key.length()<32)throw std::runtime_error("invalid NTL root prerequisite input");
 if(op<=1)zz_p::init(conv<long>(p));
 std::vector<unsigned char>kb(key.length());for(long i=0;i<key.length();i++)kb[i]=conv<unsigned long>(key[i]);RandomStream stream(kb.data()),secondary(kb.data());secondary.set_nonce(1);active_stream=&stream;
 std::ostringstream out,value;out<<'[';
 try{
 if(op==0){vec_zz_p roots;roots.SetLength(data.length());for(long i=0;i<data.length();i++)conv(roots[i],data[i]);zz_pX f;BuildFromRoots(f,roots);poly(value,f);}
 else if(op==1){zz_pX f;random(f,n);poly(value,f);}
 else if(op==2){long bound=word(p);std::vector<long>v(std::max(0L,n));VectorRandomBnd(n,v.data(),bound);value<<'[';for(long i=0;i<n;i++){if(i)value<<',';value<<'"'<<v[i]<<'"';}value<<']';}
 else if(op==3){std::shared_ptr<RandomBndGenerator>g[2];for(long j=0;j<2;j++){active_stream=(n==2&&j==1)?&secondary:&stream;g[j]=n?std::make_shared<RandomBndGenerator>(word(p)):std::make_shared<RandomBndGenerator>();}if(data.length()%3)throw std::runtime_error("invalid cached sampler commands");value<<'[';
 for(long i=0;i<data.length();i+=3){long action=conv<long>(data[i]),slot=conv<long>(data[i+1]);const ZZ&arg=data[i+2];if(i)value<<',';std::ostringstream result;value<<'[';try{
 if(action==0){active_stream=(n==2&&slot==1)?&secondary:&stream;g[slot]->build(word(arg));result<<"null";}
 else if(action==1){long count=conv<long>(arg);result<<'[';for(long j=0;j<count;j++){if(!g[slot]->p)throw std::runtime_error("RandomBndGenerator::next: uninitialized generator");if(j)result<<',';result<<'"'<<g[slot]->next()<<'"';}result<<']';}
 else if(action==2){g[slot]=std::make_shared<RandomBndGenerator>(*g[conv<long>(arg)]);result<<"null";}
 else if(action==3){g[slot]=g[conv<long>(arg)];result<<"null";}
 else if(action==4){g[slot]=std::make_shared<RandomBndGenerator>();result<<"null";}
 else if(action==5){long count=conv<long>(arg);std::vector<unsigned char>b(std::max(0L,count));stream.get(b.data(),count);result<<std::quoted(hex(b));}
 else if(action==6){stream.set_nonce(conv<unsigned long>(arg));result<<"null";}
 else if(action==7){*g[slot]=*g[conv<long>(arg)];result<<"null";}
 else throw std::runtime_error("unknown cached sampler command");
 value<<"null,null,"<<result.str();
 }catch(const WordBoundError&e){value<<"\"RangeError\","<<std::quoted(e.what())<<",null";}
 catch(const std::exception&e){value<<"\"Error\","<<std::quoted(e.what())<<",null";}
 value<<',';state(value,g);value<<']';}value<<']';}
 else throw std::runtime_error("unknown root prerequisite operation");
 out<<"null,null,"<<value.str();
 }catch(const WordBoundError&e){out<<"\"RangeError\","<<std::quoted(e.what())<<",null";}
 catch(const std::exception&e){out<<"\"Error\","<<std::quoted(e.what())<<",null";}
 std::vector<unsigned char>tail(64);stream.get(tail.data(),64);out<<',';if(op==3&&n==2){std::vector<unsigned char>other(64);secondary.get(other.data(),64);out<<'['<<std::quoted(hex(tail))<<','<<std::quoted(hex(other))<<']';}else out<<std::quoted(hex(tail));out<<']';std::cout<<out.str()<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
