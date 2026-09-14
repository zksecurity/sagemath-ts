#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <sstream>
#include <iomanip>
using namespace NTL;
static void emit(long x){std::cout<<'"'<<x<<'"';}
static void emit(zz_p x){emit(rep(x));}
template<class T>static void emit(const Vec<T>& x){std::cout<<'[';for(long i=0;i<x.length();i++){if(i)std::cout<<',';emit(x[i]);}std::cout<<']';}
static void emit(const zz_pX& x){emit(x.rep);}
static void error(const char*x){throw std::runtime_error(x);}
int main(){ErrorMsgCallback=error;std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);long op,p,initialized;ZZ exponent,scalar;vec_ZZ a,b,f;in>>op>>p>>initialized>>exponent>>scalar>>a>>b>>f;if(!in)throw std::runtime_error("invalid NTL word quotient input");
 zz_p::init(p);ZZX za,zb,zf;za.rep=a;za.normalize();zb.rep=b;zb.normalize();zf.rep=f;zf.normalize();zz_pX A,B,modulus,x;conv(A,za);conv(B,zb);conv(modulus,zf);
 if(op==4)MulByXMod(x,A,modulus);
 else if(op==5)InvMod(x,A,modulus);
 else if(op==6){long status=InvModStatus(x,A,modulus);std::cout<<'[';emit(status);std::cout<<',';emit(x);std::cout<<']'<<std::endl;continue;}
 else {
  zz_pXModulus F;if(initialized)build(F,modulus);
  if(op==0){std::cout<<'[';emit(F.n);std::cout<<',';emit(F.f);std::cout<<']'<<std::endl;continue;}
  else if(op==11){std::cout<<'[';emit(zz_pInfo->PrimeCnt);std::cout<<',';emit(NTL_zz_pX_MOD_CROSSOVER);std::cout<<',';emit(F.UseFFT);std::cout<<']'<<std::endl;continue;}
  else if(op==10){zz_pX before;rem(before,A,F);build(F,B);rem(x,A,F);std::cout<<'[';emit(before);std::cout<<',';emit(x);std::cout<<',';emit(F.f);std::cout<<',';emit(F.n);std::cout<<']'<<std::endl;continue;}
  else if(op==1)rem(x,A,F);
  else if(op==2)MulMod(x,A,B,F);
  else if(op==3)SqrMod(x,A,F);
  else if(op==7)PowerXMod(x,exponent,F);
  else if(op==8){zz_p c;conv(c,scalar);PowerXPlusAMod(x,c,exponent,F);}
  else if(op==9)PowerMod(x,A,exponent,F);
  else throw std::runtime_error("unknown NTL word quotient operation");
 }
 emit(x);std::cout<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
