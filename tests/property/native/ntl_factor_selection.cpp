#include <NTL/ZZXFactoring.h>
#include <NTL/lzz_pXFactoring.h>
#include <sstream>
#include <iomanip>
using namespace NTL;
using std::cerr;
namespace audit_ntl {
using NTL::mul;
// BUNDLED_NTL_SELECTION_ROUTINES
}
static void emit(const ZZ& x){std::cout<<'"'<<x<<'"';}
static void emit(long x){std::cout<<'"'<<x<<'"';}
static void emit(const ZZ_p& x){emit(rep(x));}
static void emit(const ZZ_pX& x);
template<class T>static void emit(const Vec<T>& x){std::cout<<'[';for(long i=0;i<x.length();i++){if(i)std::cout<<',';emit(x[i]);}std::cout<<']';}
static void emit(const ZZ_pX& x){emit(x.rep);}
static void emit(const ZZX& x){emit(x.rep);}
static void error(const char*x){throw std::runtime_error(x);}
int main(){ErrorMsgCallback=error;std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);long op;ZZ p;vec_ZZ A,B,C,packed;in>>op>>p>>A>>B>>C>>packed;
 if(!in)throw std::runtime_error("invalid NTL selection input");
 ZZ_p::init(p);vec_ZZ_pX W;W.SetLength(conv<long>(packed[0]));for(long i=0,k=1;i<W.length();i++){long n=conv<long>(packed[k++]);W[i].SetLength(n);for(long j=0;j<n;j++)conv(W[i][j],packed[k++]);W[i].normalize();}
 vec_long a,b;conv(a,A);conv(b,B);
 if(op==0){ZZX f;f.rep=A;f.normalize();audit_ntl::inplace_rev(f);emit(f);}
 else if(op==1){zz_p::init(conv<long>(p));vec_pair_zz_pX_long fac;fac.SetLength(W.length());for(long i=0;i<W.length();i++){ZZX t;conv(t,W[i]);conv(fac[i].a,t);fac[i].b=b[i];}audit_ntl::RecordPattern(a,fac);emit(a);}
 else if(op==2)emit(audit_ntl::NumFactors(a));
 else if(op==3){ZZ pd;audit_ntl::CalcPossibleDegrees(pd,a);emit(pd);}
 else if(op==4){vec_ZZ S;audit_ntl::CalcPossibleDegrees(S,W,b[0]);emit(S);}
 else if(op==5){ZZ_p lc;conv(lc,B[1]);vec_ZZ_p prod;conv(prod,C);long len=b[2],ok=audit_ntl::ConstTermTest(W,a,B[0],lc,prod,len);std::cout<<'[';emit(ok);std::cout<<',';emit(prod);std::cout<<',';emit(len);std::cout<<']';}
 else if(op==6){ZZ_pX g;ZZX t;t.rep=A;t.normalize();conv(g,t);ZZX f;audit_ntl::BalCopy(f,g);emit(f);}
 else if(op==7){ZZ_pX g;audit_ntl::mul(g,W);emit(g);}
 else if(op==8){ZZ_pX g;audit_ntl::mul(g,W,a);emit(g);}
 else if(op==9){ZZ_pX g;audit_ntl::InvMul(g,W,a);emit(g);}
 else if(op==10){audit_ntl::RemoveFactors(W,a);emit(W);}
 else if(op==11){vec_long x;audit_ntl::unpack(x,A[0],b[0]);emit(x);}
 else if(op==12){audit_ntl::SubPattern(a,b);emit(a);}
 else throw std::runtime_error("unknown NTL selection operation");
 std::cout<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
