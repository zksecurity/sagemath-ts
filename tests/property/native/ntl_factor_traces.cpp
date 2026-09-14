// Selected original ZZXFactoring.cpp bodies are inserted by the Python builder.
#include <NTL/ZZXFactoring.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_vec_ZZ.h>
#include <iostream>
#include <sstream>
#include <iomanip>
using namespace NTL;
using std::cerr;
namespace audit_ntl {
// BUNDLED_NTL_TRACE_ROUTINES
}
static void emit(const ZZ& x) { std::cout << '"' << x << '"'; }
static void emit(long x) { std::cout << '"' << x << '"'; }
template<class T> static void emit(const Vec<T>& x) {
  std::cout << '[';
  for(long i=0;i<x.length();i++){ if(i)std::cout<<',';emit(x[i]); }
  std::cout<<']';
}
static void emit(const mat_ZZ& x) {
  std::cout << '[';
  for(long i=0;i<x.NumRows();i++){ if(i)std::cout<<',';emit(x[i]); }
  std::cout<<']';
}
static mat_ZZ unpack(const vec_ZZ& v,long columns) {
  mat_ZZ m;
  long rows=conv<long>(v[0]); m.SetDims(rows,columns);
  for(long i=0,k=1;i<rows;i++)for(long j=0;j<columns;j++)m[i][j]=v[k++];
  return m;
}
static void native_error(const char* message) { throw std::runtime_error(message); }
int main(){
  ErrorMsgCallback = native_error;
  std::string line;
  while(std::getline(std::cin,line)){
    try{
      std::istringstream in(line);
      long op;vec_ZZ v,F,Tr,C,pb,A,B;
      in>>op>>v>>F>>Tr>>C>>pb>>A>>B;
      if(!in)throw std::runtime_error("invalid native trace input");
      long p=conv<long>(v[0]),d=conv<long>(v[1]),d1=conv<long>(v[2]),n=conv<long>(v[3]),bit_delta=conv<long>(v[5]),delta=conv<long>(v[9]);
      const ZZ& root=v[4];const ZZ& lc=v[6];const ZZ& P=v[7];const ZZ& pd=v[8];
      ZZX f;f.rep=F;f.normalize();
      if(op==0){audit_ntl::ComputeTrace(Tr,f,d,P);emit(Tr);}
      else if(op==1){audit_ntl::ChopTraces(C,Tr,d,pb,pd,P,lc);emit(C);}
      else if(op==2){mat_ZZ m=unpack(A,d);audit_ntl::DenseChopTraces(C,Tr,d,d1,root,pd,P,lc,m);emit(C);}
      else if(op==3){vec_long b;conv(b,C);audit_ntl::Compute_pb(b,pb,p,d,root,n);std::cout<<'[';emit(b);std::cout<<',';emit(pb);std::cout<<']';}
      else if(op==4){ZZ z=pd;audit_ntl::Compute_pdelta(delta,z,p,bit_delta);std::cout<<'[';emit(delta);std::cout<<',';emit(z);std::cout<<']';}
      else if(op==5){mat_ZZ chop=unpack(A,d),base=unpack(B,n),M;vec_vec_ZZ cv;cv.SetLength(n);for(long i=0;i<n;i++)cv[i]=chop[i];long scale;audit_ntl::BuildReductionMatrix(M,scale,n,d,pd,cv,base,0);std::cout<<'[';emit(M);std::cout<<',';emit(scale);std::cout<<']';}
      else if(op==6){long b;ZZ z;audit_ntl::Compute_pb_eff(b,z,p,d,root,n,bit_delta);std::cout<<'[';emit(b);std::cout<<',';emit(z);std::cout<<']';}
      else if(op==7)emit(audit_ntl::d1_val(bit_delta,n,d1));
      else throw std::runtime_error("unknown native trace operation");
      std::cout<<std::endl;
    }catch(const std::exception& e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}
  }
}
