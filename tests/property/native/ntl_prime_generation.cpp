#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/Lazy.h>
#include <NTL/fileio.h>
#include <NTL/SmartPtr.h>
#include <NTL/BasicThreadPool.h>
#include <sstream>
#include <iomanip>
#include <memory>
#include <vector>
#include <cstring>
#if defined(NTL_RANDOMSTREAM_AES256CTR) || defined(NTL_HAVE_AVX2) || defined(NTL_HAVE_SSSE3) || NTL_BITS_PER_LONG != 64 || NTL_BITS_PER_INT32 != 32
#error Unsupported native byte stream profile
#endif
static_assert(NTL_SP_NBITS==60 && NTL_ZZ_NBITS==64 && NTL_BITS_PER_NONCE==64, "unsupported prime word profile");
static_assert(NTL_FDOUBLE_PRECISION==double(1ULL<<52), "unsupported probability-bound profile");
namespace NTL {
NTL_IMPORT_FROM_STD
/*__NATIVE_RANDOM_BODY__*/
static RandomStream *active_stream;
RandomStream& GetCurrentRandomStream() { return *active_stream; }
static inline RandomStream& LocalGetCurrentRandomStream() { return *active_stream; }
// The explicit active stream is initialized before seeding; preserve assignment identity.
void SetSeed(const RandomStream& s) { *active_stream=s; }
/*__NATIVE_SEED_BODY__*/
/*__NATIVE_SAMPLING_BODY__*/
/*__NATIVE_PRIME_SEQUENCE__*/
/*__NATIVE_PRIME_BODY__*/
}
using namespace NTL;
static void error(const char*x){throw std::runtime_error(x);}
struct WordInputError:std::runtime_error {using std::runtime_error::runtime_error;};
static std::string hex(const std::vector<unsigned char>&b){static const char*h="0123456789abcdef";std::string s;for(auto x:b){s+=h[x>>4];s+=h[x&15];}return s;}
int main(){ErrorMsgCallback=error;SetNumThreads(1);std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);vec_ZZ key,cmd;in>>key>>cmd;if(!in||key.length()<32||cmd.length()%5)throw std::runtime_error("invalid NTL prime-generation input");
 std::vector<unsigned char>bytes(key.length());for(long i=0;i<key.length();i++)bytes[i]=conv<unsigned long>(key[i]);
 std::shared_ptr<RandomStream>s[3];for(auto&x:s)x=std::make_shared<RandomStream>(bytes.data());std::ostringstream out;out<<'[';
 for(long i=0;i<cmd.length();i+=5){long op=conv<long>(cmd[i]),slot=conv<long>(cmd[i+1]);const ZZ&a=cmd[i+2];long b=conv<long>(cmd[i+3]),c=conv<long>(cmd[i+4]);if(i)out<<',';std::ostringstream value;try{
 active_stream=s[slot].get();
 if(op==0){if(c){if(a<-(ZZ(1)<<63)||a>=(ZZ(1)<<63))throw WordInputError("ProbPrime: word input must fit a signed 64-bit integer");value<<ProbPrime(conv<long>(a),b);}else value<<ProbPrime(a,b);}
 else if(op==1)value<<ProbPrime(a);
 else if(op==2)value<<'"'<<ComputePrimeBound(conv<long>(a))<<'"';
 else if(op==3)value<<ErrBoundTest(conv<long>(a),b,c);
 else if(op>=4&&op<=8){long l=conv<long>(a);ZZ n;if(op==4){if(c)RandomPrime(n,l);else RandomPrime(n,l,b);}else if(op==5){OldRandomPrime(n,l,c?10:b);}else if(op==6){n=c?RandomPrime_long(l):RandomPrime_long(l,b);}else if(op==7){if(c)GenPrime(n,l);else GenPrime(n,l,b);}else{n=c?GenPrime_long(l):GenPrime_long(l,b);}value<<'"'<<n<<'"';}
 else if(op==9){long n=conv<long>(a);std::vector<unsigned char>v(std::max(0L,n));s[slot]->get(v.data(),n);value<<std::quoted(hex(v));}
 else if(op==10){s[slot]->set_nonce(conv<unsigned long>(a));value<<"null";}
 else if(op==11){s[slot]=std::make_shared<RandomStream>(*s[conv<long>(a)]);value<<"null";}
 else if(op==12){*s[slot]=*s[conv<long>(a)];value<<"null";}
 else if(op==13){s[slot]=s[conv<long>(a)];value<<"null";}
 else throw std::runtime_error("unknown prime-generation operation");
 out<<"[null,null,"<<value.str()<<']';
 }catch(const WordInputError&e){out<<"[\"RangeError\","<<std::quoted(e.what())<<",null]";}
 catch(const std::exception&e){out<<"[\"Error\","<<std::quoted(e.what())<<",null]";}}
 for(auto&x:s){std::vector<unsigned char>v(64);x->get(v.data(),64);if(cmd.length()||&x!=&s[0])out<<',';out<<"[null,null,"<<std::quoted(hex(v))<<']';}
 out<<']';std::cout<<out.str()<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
