#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/SmartPtr.h>
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
}
namespace reference_random {
using namespace NTL;
/*__NATIVE_RANDOM_CLASS__*/
}
using namespace NTL;
static void error(const char*x){throw std::runtime_error(x);}
static std::vector<unsigned char> bytes(const vec_ZZ& v){std::vector<unsigned char>b(v.length());for(long i=0;i<v.length();i++)b[i]=conv<unsigned long>(v[i]);return b;}
static std::string hex(const std::vector<unsigned char>&b){static const char*h="0123456789abcdef";std::string s;for(auto x:b){s+=h[x>>4];s+=h[x&15];}return s;}
static void words(std::ostream&out,const std::vector<_ntl_uint32>&v){out<<'[';for(long i=0;i<long(v.size());i++){if(i)out<<',';out<<'"'<<v[i]<<'"';}out<<']';}
int main(){ErrorMsgCallback=error;std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);long op,n;vec_ZZ av,bv,cmd;in>>op>>n>>av>>bv>>cmd;if(!in)throw std::runtime_error("invalid NTL random stream input");
 auto a=bytes(av),b=bytes(bv);std::ostringstream out;
 if(op<=2){std::vector<unsigned char> r(op==2?std::max(0L,n):32);if(op==0)sha256(a.data(),a.size(),r.data(),n);else if(op==1)hmac_sha256(b.data(),b.size(),a.data(),a.size(),r.data(),n);else DeriveKey(r.data(),n,a.data(),a.size());if(op<2)r.resize(std::max(0L,std::min(32L,n)));out<<std::quoted(hex(r));}
 else if(op==7){std::vector<unsigned char>r(32);long length=conv<long>(cmd[0]);if((static_cast<unsigned long>(length)&0xffffffffUL)>a.size())throw std::runtime_error("invalid virtual HMAC data length");hmac_sha256(b.data(),b.size(),a.data(),length,r.data(),n);r.resize(std::max(0L,std::min(32L,n)));out<<std::quoted(hex(r));}
 else if(op==3){
 if(a.size()<32)throw std::runtime_error("RandomStream: key must contain at least 32 bytes");
 std::shared_ptr<reference_random::RandomStream> streams[3];for(auto&s:streams)s=std::make_shared<reference_random::RandomStream>(a.data());out<<'[';
 if(cmd.length()%3)throw std::runtime_error("invalid stream commands");
 for(long i=0;i<cmd.length();i+=3){long action=conv<long>(cmd[i]),slot=conv<long>(cmd[i+1]);if(i)out<<',';try{
 if(action==0){long count=conv<long>(cmd[i+2]);std::vector<unsigned char>r(std::max(0L,count));streams[slot]->get(r.data(),count);out<<"[null,"<<std::quoted(hex(r))<<']';}
 else{if(action==1)streams[slot]->set_nonce(conv<unsigned long>(cmd[i+2]));else if(action==2)streams[slot]=std::make_shared<reference_random::RandomStream>(*streams[conv<long>(cmd[i+2])]);else if(action==3)streams[slot]=streams[conv<long>(cmd[i+2])];else if(action==4)*streams[slot]=*streams[conv<long>(cmd[i+2])];else throw std::runtime_error("unknown stream command");out<<"[null,null]";}
 }catch(const std::exception&e){out<<'['<<std::quoted(e.what())<<",null]";}}
 out<<']';}
 else if(op==4||op==5){if(av.length()<16)throw std::runtime_error(op==4?"salsa20_core: state must contain at least 16 words":"salsa20_apply: state must contain at least 16 words");std::vector<_ntl_uint32>state(av.length());for(long i=0;i<av.length();i++)state[i]=conv<unsigned long>(av[i]);if(op==4){salsa20_core(state.data());words(out,state);}else{std::vector<_ntl_uint32>block(16);salsa20_apply(state.data(),block.data());out<<'[';words(out,state);out<<',';words(out,block);out<<']';}}
 else if(op==6){if(a.size()<32)throw std::runtime_error("RandomStream: key must contain at least 32 bytes");std::vector<_ntl_uint32>state(16);salsa20_init(state.data(),a.data());words(out,state);}
 else throw std::runtime_error("unknown NTL random stream operation");
 std::cout<<out.str()<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
