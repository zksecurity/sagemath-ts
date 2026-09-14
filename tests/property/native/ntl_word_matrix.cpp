#include <NTL/mat_lzz_p.h>
#include <NTL/BasicThreadPool.h>
#include <sstream>
#include <iomanip>
#include <cstdint>
using namespace NTL;
static void emit(zz_p x){std::cout<<'"'<<rep(x)<<'"';}
template<class T>static void emit(const Vec<T>& x){std::cout<<'[';for(long i=0;i<x.length();i++){if(i)std::cout<<',';emit(x[i]);}std::cout<<']';}
static void emit(const mat_zz_p& x){std::cout<<'[';for(long i=0;i<x.NumRows();i++){if(i)std::cout<<',';emit(x[i]);}std::cout<<']';}
static void error(const char*x){throw std::runtime_error(x);}
int main(){ErrorMsgCallback=error;SetNumThreads(1);std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);long op,p;unsigned long seed;vec_long dims;vec_ZZ a,b;in>>op>>p>>dims>>seed>>a>>b;if(!in)throw std::runtime_error("invalid NTL word matrix input");
 zz_p::init(p);long m=dims[0],k=dims[1],k2=dims[2],n=dims[3];mat_zz_p A,B,C;A.SetDims(m,k);B.SetDims(k2,n);
 uint64_t state=seed;auto next=[&](){state=state*UINT64_C(6364136223846793005)+UINT64_C(1442695040888963407);long value=((state^(state>>32))>>1)%uint64_t(p);return (state&1)?-value:value;};
 for(long i=0,t=0;i<m;i++)for(long j=0;j<k;j++){if(op)conv(A[i][j],next());else conv(A[i][j],a[t++]);}
 for(long i=0,t=0;i<k2;i++)for(long j=0;j<n;j++){if(op)conv(B[i][j],next());else conv(B[i][j],b[t++]);}
 mul(C,A,B);
 if(op){long width=(NumBits(p-1)+3)/4;std::cout<<"[\""<<C.NumRows()<<"\",\""<<C.NumCols()<<"\",\"";for(long i=0;i<C.NumRows();i++)for(long j=0;j<C.NumCols();j++)std::cout<<std::hex<<std::setfill('0')<<std::setw(width)<<rep(C[i][j]);std::cout<<std::dec<<"\"]";}
 else emit(C);
 std::cout<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
