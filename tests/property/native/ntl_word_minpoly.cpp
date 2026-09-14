#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <NTL/BasicThreadPool.h>
#include <sstream>
#include <iomanip>
using namespace NTL;
static void emit(const zz_pX& x){std::cout<<'[';for(long i=0;i<x.rep.length();i++){if(i)std::cout<<',';std::cout<<'"'<<rep(x.rep[i])<<'"';}std::cout<<']';}
static void error(const char*x){throw std::runtime_error(x);}
int main(){ErrorMsgCallback=error;SetNumThreads(1);std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);long p,m;vec_ZZ a;in>>p>>m>>a;if(!in)throw std::runtime_error("invalid NTL sequence input");
 zz_p::init(p);vec_zz_p seq;seq.SetLength(a.length());for(long i=0;i<a.length();i++)conv(seq[i],a[i]);zz_pX h;MinPolySeq(h,seq,m);emit(h);std::cout<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
