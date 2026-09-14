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
static RandomStream *active_stream;
static inline RandomStream& LocalGetCurrentRandomStream() { return *active_stream; }
/*__NATIVE_SAMPLING_BODY__*/
}
using namespace NTL;
static void error(const char*x){throw std::runtime_error(x);}
struct WordBoundError:std::runtime_error {using std::runtime_error::runtime_error;};
static std::string hex(const std::vector<unsigned char>&b){static const char*h="0123456789abcdef";std::string s;for(auto x:b){s+=h[x>>4];s+=h[x&15];}return s;}
int main(){ErrorMsgCallback=error;std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);vec_ZZ key,cmd;in>>key>>cmd;if(!in||key.length()<32||cmd.length()%3)throw std::runtime_error("invalid NTL integer sampling input");
 std::vector<unsigned char>a(key.length());for(long i=0;i<key.length();i++)a[i]=conv<unsigned long>(key[i]);
 std::shared_ptr<reference_random::RandomStream>s[3];for(auto&x:s)x=std::make_shared<reference_random::RandomStream>(a.data());std::ostringstream out;out<<'[';
 for(long i=0;i<cmd.length();i+=3){long op=conv<long>(cmd[i]),slot=conv<long>(cmd[i+1]);const ZZ&arg=cmd[i+2];if(i)out<<',';std::ostringstream value;try{
 reference_random::active_stream=s[slot].get();
 if(op==0)value<<'"'<<reference_random::RandomWord()<<'"';
 else if(op==1){long n=conv<long>(arg);std::vector<unsigned long>v(std::max(0L,n));reference_random::VectorRandomWord(n,v.data());value<<'[';for(long j=0;j<n;j++){if(j)value<<',';value<<'"'<<v[j]<<'"';}value<<']';}
 else if(op==2)value<<'"'<<reference_random::RandomBits_long(conv<long>(arg))<<'"';
 else if(op==3)value<<'"'<<reference_random::RandomBits_ulong(conv<long>(arg))<<'"';
 else if(op==4)value<<'"'<<reference_random::RandomLen_long(conv<long>(arg))<<'"';
 else if(op==5||op==6||op==8||op==9){ZZ x;if(op==5||op==8)reference_random::RandomBits(x,conv<long>(arg));else reference_random::RandomLen(x,conv<long>(arg));value<<'"'<<x<<'"';}
 else if(op==7){ZZ x;reference_random::RandomBnd(x,arg);value<<'"'<<x<<'"';}
 else if(op==10){if(arg<-(ZZ(1)<<63)||arg>=(ZZ(1)<<63))throw WordBoundError("RandomBnd: word bound must fit a signed 64-bit integer");value<<'"'<<reference_random::RandomBnd(conv<long>(arg))<<'"';}
 else if(op==11){long n=conv<long>(arg);std::vector<unsigned char>b(std::max(0L,n));s[slot]->get(b.data(),n);value<<std::quoted(hex(b));}
 else if(op==12){s[slot]->set_nonce(conv<unsigned long>(arg));value<<"null";}
 else if(op==13){s[slot]=std::make_shared<reference_random::RandomStream>(*s[conv<long>(arg)]);value<<"null";}
 else if(op==14){*s[slot]=*s[conv<long>(arg)];value<<"null";}
 else if(op==15){s[slot]=s[conv<long>(arg)];value<<"null";}
 else throw std::runtime_error("unknown integer sampling operation");
 out<<"[null,null,"<<value.str()<<']';
 }catch(const WordBoundError&e){out<<"[\"RangeError\","<<std::quoted(e.what())<<",null]";}
 catch(const std::exception&e){out<<"[\"Error\","<<std::quoted(e.what())<<",null]";}}
 for(auto&x:s){std::vector<unsigned char>b(64);x->get(b.data(),64);if(cmd.length()||&x!=&s[0])out<<',';out<<"[null,null,"<<std::quoted(hex(b))<<']';}
 out<<']';std::cout<<out.str()<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
