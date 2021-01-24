#include "MultiDimArray.cpp"
#include "NTL/ZZ_p.h"
#include <iostream>

using namespace std;

/******************************************************
*
*   TEST #2 (From Multidimensional paper, Arce et al.)
*
*******************************************************/

int main() {

    NTL::ZZ p(2);
    NTL::ZZ_p::init(p);
    typedef NTL::ZZ_p F;
    const int m = 2;
    blitz::Array<F,m> A(6,7);

    A = (F)0, (F)0, (F)1, (F)1, (F)0, (F)1, (F)0,
        (F)1, (F)0, (F)0, (F)0, (F)1, (F)1, (F)0,
        (F)0, (F)0, (F)0, (F)1, (F)1, (F)0, (F)1,
        (F)1, (F)1, (F)0, (F)1, (F)0, (F)0, (F)0,
        (F)0, (F)1, (F)0, (F)0, (F)0, (F)1, (F)1,
        (F)1, (F)0, (F)1, (F)0, (F)0, (F)0, (F)1;

    MultiDimArray<F,m> array(A);
    array.RST();
    cout << "Complexity: " << array.complexity() << endl;

    return 0;
}