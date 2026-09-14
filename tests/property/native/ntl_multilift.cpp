// Selected original ZZXFactoring.cpp bodies are inserted by the Python builder.
#include <NTL/ZZXFactoring.h>
#include <NTL/lzz_pXFactoring.h>
#include <iostream>
#include <sstream>
#include <iomanip>
using namespace NTL;
using std::cerr;
namespace audit_ntl {
// BUNDLED_NTL_LIFT_ROUTINES
}
static void emit(const ZZ& x) { std::cout << '"' << x << '"'; }
static void emit(const ZZX& x);
template<class T> static void emit(const Vec<T>& x) {
 std::cout<<'[';for(long i=0;i<x.length();i++){if(i)std::cout<<',';emit(x[i]);}std::cout<<']';
}
static void emit(const ZZX& x) {emit(x.rep);}
static void native_error(const char* s){throw std::runtime_error(s);}
int main(){
 ErrorMsgCallback=native_error;
 std::string line;
 while(std::getline(std::cin,line)){
  try{
   std::istringstream in(line);long op,p,e;vec_ZZ F,A,B;in>>op>>p>>e>>F>>A>>B;
   if(!in)throw std::runtime_error("invalid NTL lift input");
   ZZX f;f.rep=F;f.normalize();
   if(op==0||op==1){ZZX a,b,c;a.rep=A;a.normalize();b.rep=B;b.normalize();if(op==0)mul(c,a,b);else sqr(c,a);emit(c);}
   else if(op==2){
    zz_p::init(p);vec_zz_pX fac;fac.SetLength(conv<long>(A[0]));
    for(long i=0,k=1;i<fac.length();i++){long n=conv<long>(A[k++]);for(long j=0;j<n;j++)SetCoeff(fac[i],j,conv<zz_p>(A[k++]));}
    vec_ZZX out;audit_ntl::MultiLift(out,fac,f,e,0);emit(out);
   }else throw std::runtime_error("unknown NTL lift operation");
   std::cout<<std::endl;
  }catch(const std::exception& e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}
 }
}
