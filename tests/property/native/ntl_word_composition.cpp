#include <NTL/ZZX.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/BasicThreadPool.h>
#include <sstream>
#include <iomanip>
using namespace NTL;
static void emit(long x){std::cout<<'"'<<x<<'"';}
static void emit(zz_p x){emit(rep(x));}
static void emit(const zz_pX& x);
template<class T>static void emit(const Vec<T>& x){std::cout<<'[';for(long i=0;i<x.length();i++){if(i)std::cout<<',';emit(x[i]);}std::cout<<']';}
static void emit(const zz_pX& x){emit(x.rep);}
static void emit(const mat_zz_p& x){std::cout<<'[';for(long i=0;i<x.NumRows();i++){if(i)std::cout<<',';emit(x[i]);}std::cout<<']';}
static void emit(const zz_pXNewArgument& x){std::cout<<'[';emit(x.mat);std::cout<<',';emit(x.poly);std::cout<<']';}
static void error(const char*x){throw std::runtime_error(x);}
int main(){ErrorMsgCallback=error;SetNumThreads(1);std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);long op,p;vec_long param;vec_ZZ packed;in>>op>>p>>param>>packed;if(!in)throw std::runtime_error("invalid NTL composition input");
 zz_p::init(p);vec_zz_pX w;w.SetLength(conv<long>(packed[0]));for(long i=0,k=1;i<w.length();i++){long n=conv<long>(packed[k++]);w[i].SetLength(n);for(long j=0;j<n;j++)conv(w[i][j],packed[k++]);w[i].normalize();}
 if(op==8){zz_pX x;GCD(x,w[2],w[3]);emit(x);std::cout<<std::endl;continue;}
 zz_pXModulus F;if(param[2])build(F,w[0]);zz_pXNewArgument H;zz_pX x,y,z;
 if(op==0){CompMod(x,w[2],w[1],F);emit(x);}
 else if(op==1){if(param[3])build(H,w[1],F,param[0]);CompMod(x,w[2],H,F);emit(x);}
 else if(op==2){build(H,w[1],F,param[0]);emit(H);}
 else if(op==3){if(param[3])build(H,w[1],F,param[0]);zz_pXModulus G(w[5]);reduce(H,G);CompMod(x,w[2],H,G);std::cout<<'[';emit(H);std::cout<<',';emit(x);std::cout<<']';}
 else if(op==4){Comp2Mod(x,y,w[2],w[3],w[1],F);std::cout<<'[';emit(x);std::cout<<',';emit(y);std::cout<<']';}
 else if(op==5){Comp3Mod(x,y,z,w[2],w[3],w[4],w[1],F);std::cout<<'[';emit(x);std::cout<<',';emit(y);std::cout<<',';emit(z);std::cout<<']';}
 else if(op==6){TraceMap(x,w[2],param[1],F,w[1]);emit(x);}
 else if(op==7){PowerCompose(x,w[1],param[1],F);emit(x);}
 else throw std::runtime_error("unknown NTL composition operation");
 std::cout<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
