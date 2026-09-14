#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_lzz_p.h>
#include <NTL/vec_ZZVec.h>
#include <NTL/BasicThreadPool.h>
#include <sstream>
#include <iomanip>
#define PAR_THRESH (40000.0)
namespace NTL {
NTL_IMPORT_FROM_STD
/*__NATIVE_BODIES__*/
}
using namespace NTL;
static_assert(NTL_BITS_PER_LONG==64 && NTL_ZZ_NBITS==64 && NTL_SP_NBITS==60,"unsupported integer CRT profile");
static void emit(std::ostream&out,const ZZ& x){out<<'"'<<x<<'"';}
static void emit(std::ostream&out,const ZZ_p& x){emit(out,rep(x));}
template<class T>static void emit(std::ostream&out,const Vec<T>&x){out<<'[';for(long i=0;i<x.length();i++){if(i)out<<',';emit(out,x[i]);}out<<']';}
template<class T>static void emit(std::ostream&out,const Mat<T>&x){out<<'[';for(long i=0;i<x.NumRows();i++){if(i)out<<',';emit(out,x[i]);}out<<']';}
static void error(const char*x){throw std::runtime_error(x);}
struct RangeError:std::runtime_error{using std::runtime_error::runtime_error;};
template<class T>static void make(Mat<T>& M,const vec_ZZ& v,long n,long m,bool ragged){
 if(ragged){Vec<Vec<T>> rows;rows.SetLength(2);rows[0].SetLength(1);rows[1].SetLength(2);MakeMatrix(M,rows);return;}
 M.SetDims(n,m);if(v.length()!=n*m)throw std::runtime_error("invalid flat matrix length");
 for(long i=0,k=0;i<n;i++)for(long j=0;j<m;j++)conv(M[i][j],v[k++]);
}
int main(){ErrorMsgCallback=error;SetNumThreads(1);std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);long op;ZZ a,p;vec_long dims;vec_ZZ g,G;in>>op>>a>>p>>dims>>g>>G;
 if(!in||dims.length()!=6)throw std::runtime_error("invalid NTL integer CRT input");
 std::ostringstream out;
 if(op<=2){out<<'[';for(long i=0;i<g.length();i++){if(i)out<<',';try{
  if(op==0){out<<"[null,null,"<<CRTInRange(g[i],a)<<']';continue;}
  ZZ gg=g[i],aa=a,GG=G[i];
  if(op==2&&(p<-(ZZ(1)<<63)||p>=(ZZ(1)<<63)||GG<-(ZZ(1)<<63)||GG>=(ZZ(1)<<63)))throw RangeError("CRT: word inputs must fit signed 64-bit integers");
  if(a<=0||p<=1||GG<0||GG>=p)throw RangeError("CRT: require a > 0, p > 1 and 0 <= G < p");
  long changed=op==1?CRT(gg,aa,GG,p):CRT(gg,aa,conv<long>(GG),conv<long>(p));
  out<<"[null,null,["<<changed<<',';emit(out,gg);out<<',';emit(out,aa);out<<"]]";
 }catch(const RangeError&e){out<<"[\"RangeError\","<<std::quoted(e.what())<<",null]";}
 catch(const std::exception&e){out<<"[\"Error\","<<std::quoted(e.what())<<",null]";}}
 out<<']';
 }else if(op==3){
  if(p<-(ZZ(1)<<63)||p>=(ZZ(1)<<63))throw std::runtime_error(p<=1?"zz_pContext: p must be > 1":"zz_pContext: modulus too big");
  zz_p::init(conv<long>(p));mat_ZZ M;mat_zz_p R;make(M,g,dims[0],dims[1],dims[4]);make(R,G,dims[2],dims[3],dims[5]);
  if(M.NumRows()!=R.NumRows()||M.NumCols()!=R.NumCols())throw std::runtime_error("CRT: dimension mismatch");
  if(a<=0)throw RangeError("CRT: require a > 0");
  ZZ aa=a;long changed=CRT(M,aa,R);out<<'['<<changed<<',';emit(out,M);out<<',';emit(out,aa);out<<']';
 }else if(op==4){mat_ZZ M;make(M,g,dims[0],dims[1],dims[4]);out<<DetBound(M);}
 else if(op==5){ZZ_p::init(p);mat_ZZ_p M;make(M,g,dims[0],dims[1],dims[4]);ZZ_p d;determinant(d,M);emit(out,d);}
 else throw std::runtime_error("unknown integer CRT operation");
 std::cout<<out.str()<<std::endl;
 }catch(const RangeError&e){std::cout<<"RANGE "<<std::quoted(e.what())<<std::endl;}
 catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
