#include <NTL/ZZ_p.h>
#include <vector>

using namespace std;

void reduce( vector< vector<NTL::ZZ_p>* > M) {
    long m = M.size() - 2;
    long n = M[0]->size()-1;
    const NTL::ZZ_p zero((long)0);
    vector<NTL::ZZ_p> lambda(m, zero);

}



int main()
{
   NTL::ZZ p = (NTL::ZZ) 7;
   NTL::ZZ_p::init(p);

   NTL::ZZ_p x((long) 6);
   NTL::ZZ_p y((long) 1);

   cout << x+y << "\n";


}
