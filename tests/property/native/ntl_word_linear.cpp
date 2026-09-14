#include <NTL/mat_lzz_p.h>
#include <NTL/BasicThreadPool.h>
#include <sstream>
#include <iomanip>
namespace NTL {
NTL_IMPORT_FROM_STD
/*__NATIVE_CONTEXT_BODY__*/
}
using namespace NTL;
static_assert(NTL_BITS_PER_LONG==64 && NTL_SP_NBITS==60,"unsupported word profile");
static void emit(std::ostream&out,zz_p x){out<<'"'<<rep(x)<<'"';}
template<class T>static void emit(std::ostream&out,const Vec<T>&x){out<<'[';for(long i=0;i<x.length();i++){if(i)out<<',';emit(out,x[i]);}out<<']';}
static void emit(std::ostream&out,const mat_zz_p&x){out<<'[';for(long i=0;i<x.NumRows();i++){if(i)out<<',';emit(out,x[i]);}out<<']';}
static void error(const char*x){throw std::runtime_error(x);}
int main(){ErrorMsgCallback=error;SetNumThreads(1);std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);ZZ p0;vec_long dims,cmd;vec_ZZ av,bv,prev,xprev;in>>p0>>dims>>av>>bv>>prev>>xprev>>cmd;
 if(!in||dims.length()!=4||cmd.length()%4)throw std::runtime_error("invalid NTL word linear-algebra input");
 // Labeled adapter guards only for values unrepresentable by native long.
 if(p0<-(ZZ(1)<<63)||p0>=(ZZ(1)<<63))throw std::runtime_error(p0<=1?"zz_pContext: p must be > 1":"zz_pContext: modulus too big");
 zz_p::init(conv<long>(p0));
 mat_zz_p A,previous;A.SetDims(dims[0],dims[1]);previous.SetDims(dims[2],dims[3]);
 if(av.length()!=dims[0]*dims[1]||prev.length()!=dims[2]*dims[3])throw std::runtime_error("invalid flat matrix length");
 for(long i=0,t=0;i<A.NumRows();i++)for(long j=0;j<A.NumCols();j++)conv(A[i][j],av[t++]);
 for(long i=0,t=0;i<previous.NumRows();i++)for(long j=0;j<previous.NumCols();j++)conv(previous[i][j],prev[t++]);
 vec_zz_p b,previous_x;conv(b,bv);conv(previous_x,xprev);std::ostringstream out;out<<'[';
 for(long i=0;i<cmd.length();i+=4){if(i)out<<',';long op=cmd[i],w=cmd[i+1];bool relax=cmd[i+2],left=cmd[i+3];std::ostringstream value;
 try{mat_zz_p X=previous;vec_zz_p x=previous_x;zz_p d;
 if(op==0){X=A;long rank=gauss(X,w);value<<'['<<rank<<',';emit(value,X);value<<']';}
 else if(op==1){X=A;long rank=gauss(X);value<<'['<<rank<<',';emit(value,X);value<<']';}
 else if(op==2){image(X,A);emit(value,X);}
 else if(op==3){kernel(X,A);emit(value,X);}
 else if(op==4){relaxed_inv(d,X,A,relax);value<<'[';emit(value,d);value<<',';emit(value,X);value<<']';}
 else if(op==5){X=inv(A);emit(value,X);}
 else if(op==6){relaxed_determinant(d,A,relax);emit(value,d);}
 else if(op==7){d=determinant(A);emit(value,d);}
 else if(op==8||op==9){if(op==8){if(left)relaxed_solve(d,x,A,b,relax);else relaxed_solve(d,A,x,b,relax);}else{if(left)solve(d,x,A,b);else solve(d,A,x,b);}value<<'[';emit(value,d);value<<',';emit(value,x);value<<']';}
 else if(op==10){relaxed_inv(d,X,A);value<<'[';emit(value,d);value<<',';emit(value,X);value<<']';}
 else if(op==11){d=relaxed_determinant(A);emit(value,d);}
 else if(op==12){relaxed_solve(d,A,x,b);value<<'[';emit(value,d);value<<',';emit(value,x);value<<']';}
 else if(op==13){Vec<vec_zz_p> rows;rows.SetLength(2);rows[0].SetLength(1);rows[1].SetLength(2);MakeMatrix(X,rows);emit(value,X);}
 else if(op==14){X.SetDims(0,-1);emit(value,X);}
 else throw std::runtime_error("unknown word linear-algebra operation");
 out<<"[null,null,"<<value.str()<<']';
 }catch(const std::exception&e){out<<"[\"Error\","<<std::quoted(e.what())<<",null]";}}
 if(cmd.length())out<<',';out<<"[null,null,";emit(out,A);out<<"]]";std::cout<<out.str()<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
