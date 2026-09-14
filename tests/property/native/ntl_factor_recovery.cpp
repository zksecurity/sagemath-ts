#include "LLL.cpp"
#include <NTL/ZZXFactoring.h>
#include <NTL/lzz_pXFactoring.h>
#include <sstream>
#include <iomanip>
using namespace NTL;
using std::cerr;
namespace audit_ntl {
#define MultiLift audit_MultiLift
// BUNDLED_NTL_RECOVERY_ROUTINES
#undef MultiLift
}
static void emit(const ZZ& x){std::cout<<'"'<<x<<'"';}
static void emit(long x){std::cout<<'"'<<x<<'"';}
static void emit(const ZZX& x);
template<class T>static void emit(const Vec<T>& x){std::cout<<'[';for(long i=0;i<x.length();i++){if(i)std::cout<<',';emit(x[i]);}std::cout<<']';}
static void emit(const ZZX& x){emit(x.rep);}
static void emit(const mat_ZZ& x){std::cout<<'[';for(long i=0;i<x.NumRows();i++){if(i)std::cout<<',';emit(x[i]);}std::cout<<']';}
static void error(const char*x){throw std::runtime_error(x);}
int main(){ErrorMsgCallback=error;std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);long op;vec_ZZ v,F,D,flat,W;in>>op>>v>>F>>D>>flat>>W;
 if(!in)throw std::runtime_error("invalid NTL recovery input");
 ZZX f;f.rep=F;f.normalize();
 if(op==0)emit(audit_ntl::PolyEval(f,v[0]));
 else if(op==1)emit(audit_ntl::RootBound(f));
 else if(op==2){long m=conv<long>(v[0]),n=conv<long>(v[1]),C=conv<long>(v[2]),r=conv<long>(v[3]),d=conv<long>(v[4]);mat_ZZ M,B;M.SetDims(m,n);for(long i=0,k=0;i<m;i++)for(long j=0;j<n;j++)M[i][j]=flat[k++];audit_ntl::CutAway(B,D,M,C,r,d);emit(B);}
 else if(op==3){ZZ P=v[0];long e=conv<long>(v[1]),p=conv<long>(v[2]),bound=conv<long>(v[3]),doubling=conv<long>(v[4]);vec_ZZX w;w.SetLength(conv<long>(W[0]));for(long i=0,k=1;i<w.length();i++){long n=conv<long>(W[k++]);w[i].SetLength(n);for(long j=0;j<n;j++)w[i][j]=W[k++];w[i].normalize();}audit_ntl::AdditionalLifting(P,e,w,p,bound,f,doubling,0);std::cout<<'[';emit(P);std::cout<<',';emit(e);std::cout<<',';emit(w);std::cout<<']';}
 else throw std::runtime_error("unknown NTL recovery operation");
 std::cout<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
