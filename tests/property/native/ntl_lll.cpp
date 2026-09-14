#include "LLL.cpp"
#include <sstream>
#include <iomanip>
using namespace NTL;
static void emit(const ZZ& x){std::cout<<'"'<<x<<'"';}
static void emit(long x){std::cout<<'"'<<x<<'"';}
template<class T>static void emit(const Vec<T>& x){std::cout<<'[';for(long i=0;i<x.length();i++){if(i)std::cout<<',';emit(x[i]);}std::cout<<']';}
static void emit(const mat_ZZ& x){std::cout<<'[';for(long i=0;i<x.NumRows();i++){if(i)std::cout<<',';emit(x[i]);}std::cout<<']';}
static void error(const char* x){throw std::runtime_error(x);}
int main(){ErrorMsgCallback=error;std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);long op,m,n,a,b,transform;vec_ZZ flat,y,initial;in>>op>>m>>n>>a>>b>>transform>>flat>>y>>initial;
 mat_ZZ M,U;M.SetDims(m,n);for(long i=0,k=0;i<m;i++)for(long j=0;j<n;j++)M[i][j]=flat[k++];
 if(op==4){ZZ d,s,t;XGCD(d,s,t,y[0],y[1]);std::cout<<'[';emit(s);std::cout<<',';emit(t);std::cout<<',';emit(d);std::cout<<']';}
 else if(op==3){vec_ZZ x=initial;long ok=LatticeSolve(x,M,y,a);std::cout<<'[';emit(ok);std::cout<<',';emit(x);std::cout<<']';}
 else{
  ZZ det;vec_ZZ D;long rank;
  if(op==0)rank=transform?LLL(det,M,U,a,b,0):LLL(det,M,a,b,0);
  else if(op==1)rank=transform?LLL_plus(D,M,U,a,b,0):LLL_plus(D,M,a,b,0);
  else rank=transform?image(det,M,U,0):image(det,M,0);
  std::cout<<'[';emit(rank);std::cout<<',';if(op==1)emit(D);else emit(det);std::cout<<',';emit(M);std::cout<<',';if(transform)emit(U);else std::cout<<"null";std::cout<<']';
 }
 std::cout<<std::endl;
 }catch(const std::exception& e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
