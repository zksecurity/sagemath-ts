#include <NTL/ZZX.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/BasicThreadPool.h>
#include <sstream>
#include <iomanip>
using namespace NTL;
namespace NTL { void SFCanZass1(vec_pair_zz_pX_long&, zz_pX&, const zz_pX&, long); }
static void emit(long x){std::cout<<'"'<<x<<'"';}
static void emit(zz_p x){emit(rep(x));}
static void emit(const zz_pX& x){std::cout<<'[';for(long i=0;i<x.rep.length();i++){if(i)std::cout<<',';emit(x.rep[i]);}std::cout<<']';}
static void emit(const vec_pair_zz_pX_long& x){std::cout<<'[';for(long i=0;i<x.length();i++){if(i)std::cout<<',';std::cout<<'[';emit(x[i].a);std::cout<<',';emit(x[i].b);std::cout<<']';}std::cout<<']';}
static void error(const char*x){throw std::runtime_error(x);}
int main(){ErrorMsgCallback=error;SetNumThreads(1);std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);long op,p;ZZX a,b;in>>op>>p>>a>>b;if(!in)throw std::runtime_error("invalid NTL distinct-degree input");
 zz_p::init(p);zz_pX f=conv<zz_pX>(a),h=conv<zz_pX>(b);vec_pair_zz_pX_long u;
 if(op==0){NewDDF(u,f,h);emit(u);}
 else if(op==1){SFCanZass1(u,h,f,0);std::cout<<'[';emit(u);std::cout<<',';emit(h);std::cout<<']';}
 else if(op==2){zz_pXModulus F(f);PowerXMod(h,p,F);NewDDF(u,f,h);emit(u);}
 else throw std::runtime_error("unknown NTL distinct-degree operation");
 std::cout<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
