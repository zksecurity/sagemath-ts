#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <NTL/BasicThreadPool.h>
#include <sstream>
#include <iomanip>
using namespace NTL;
static void emit(long x){std::cout<<'"'<<x<<'"';}
static void emit(zz_p x){emit(rep(x));}
template<class T>static void emit(const Vec<T>& x){std::cout<<'[';for(long i=0;i<x.length();i++){if(i)std::cout<<',';emit(x[i]);}std::cout<<']';}
static void emit(const zz_pX& x){emit(x.rep);}
static void emit(const zz_pXMultiplier& B){std::cout<<'[';emit(B.b);std::cout<<',';emit(B.UseFFT);std::cout<<',';emit(B.val());std::cout<<']';}
static void error(const char*x){throw std::runtime_error(x);}
int main(){ErrorMsgCallback=error;SetNumThreads(1);std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);long op,p;vec_long param;vec_ZZ packed;in>>op>>p>>param>>packed;if(!in)throw std::runtime_error("invalid NTL projection input");
 zz_p::init(p);Vec<vec_ZZ> raw;vec_zz_pX w;long count=conv<long>(packed[0]);raw.SetLength(count);w.SetLength(count);
 for(long i=0,k=1;i<count;i++){long n=conv<long>(packed[k++]);raw[i].SetLength(n);w[i].SetLength(n);for(long j=0;j<n;j++){raw[i][j]=packed[k++];conv(w[i][j],raw[i][j]);}w[i].normalize();}
 zz_pXModulus F;if(param[2])build(F,w[0]);zz_pXMultiplier B;zz_pXNewArgument H;zz_pX x;vec_zz_p a,y;conv(a,raw[2]);
 if(op<=2 || op==5){if(param[3])build(B,w[1],F);}
 if(op==0)emit(B);
 else if(op==1){MulMod(x,w[2],B,F);emit(x);}
 else if(op==2){UpdateMap(y,a,B,F);emit(y);}
 else if(op==3){ProjectPowers(y,a,param[1],w[3],F);emit(y);}
 else if(op==4){if(param[3])build(H,w[3],F,param[0]);
 // Adapter guard for the independently recorded native zero-row matrix crash.
 if(param[1]==0 && H.mat.NumRows()>0 && a.length()<=F.n)throw std::runtime_error("ProjectPowers: prepared argument requires a positive count");
 ProjectPowers(y,a,param[1],H,F);emit(y);}
 else if(op==5){std::string err;try{build(B,w[4],F);}catch(const std::exception&e){err=e.what();}
 std::ostringstream saved;std::streambuf* old=std::cout.rdbuf(saved.rdbuf());emit(B);std::cout.rdbuf(old);
 std::cout<<'[';if(err.empty())std::cout<<"null";else std::cout<<std::quoted(err);std::cout<<','<<saved.str();
 for(long which=0;which<2;which++){std::cout<<',';try{if(which==0){MulMod(x,w[2],B,F);std::cout<<"[null,";emit(x);std::cout<<']';}else{UpdateMap(y,a,B,F);std::cout<<"[null,";emit(y);std::cout<<']';}}catch(const std::exception&e){std::cout<<'['<<std::quoted(e.what())<<",null]";}}
 std::cout<<']';}
 else throw std::runtime_error("unknown NTL projection operation");
 std::cout<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
