#include <NTL/ZZ.h>
#include <NTL/vec_ZZ.h>
#include <NTL/Lazy.h>
#include <NTL/SmartPtr.h>
#include <sstream>
#include <iomanip>
#include <memory>
namespace NTL {
/*__NATIVE_PRIME_SEQUENCE__*/
}
using namespace NTL;
static void error(const char*x){throw std::runtime_error(x);}
int main(){ErrorMsgCallback=error;std::string line;while(std::getline(std::cin,line)){try{
 std::istringstream in(line);vec_ZZ commands;in>>commands;if(!in||commands.length()%3)throw std::runtime_error("invalid NTL prime sequence input");
 std::unique_ptr<PrimeSeq> seq[2];for(auto& s:seq)s=std::make_unique<PrimeSeq>();std::cout<<'[';
 for(long i=0;i<commands.length();i+=3){long op=conv<long>(commands[i]),slot=conv<long>(commands[i+1]),arg=conv<long>(commands[i+2]);if(i)std::cout<<',';
 if(op==0){std::cout<<'[';for(long j=0;j<arg;j++){if(j)std::cout<<',';std::cout<<'"'<<seq[slot]->next()<<'"';}std::cout<<']';}
 else if(op==1){seq[slot]->reset(arg);std::cout<<"null";}
 else if(op==2){seq[slot]=std::make_unique<PrimeSeq>();std::cout<<"null";}
 else throw std::runtime_error("unknown NTL prime sequence operation");
 }std::cout<<']'<<std::endl;
 }catch(const std::exception&e){std::cout<<"ERROR "<<std::quoted(e.what())<<std::endl;}}}
