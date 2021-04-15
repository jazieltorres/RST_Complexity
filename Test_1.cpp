//#include "MultiDimArray.cpp"
#include "NTL/ZZ_p.h"
#include "MultiDimArray_GF2.cpp"
#include "NTL/GF2.h"
#include <iostream>

using namespace std;

/******************************************************
*
*       TEST #1
*
*******************************************************/

int main() {
    vector<int> seq = {1,0,0,0,1,1,1,1,0,1,0,0,0,0,1,1,0,1,1,1,1,0,1,0,0,0,1,0,1,0,0};
    const unsigned int dim = 1;
    cout << "Length: " << seq.size() << endl;
    typedef NTL::GF2 F;
    F zero(0);
    F one(1);
    blitz::TinyVector<int,1> v(seq.size());
    MultiDimArray_GF2<dim> A(v);
    for (int i=0; i<seq.size(); i++) {
        blitz::TinyVector<int,1> position(i);
        if (seq[i] == 1)
            A.set_at(position, one);
        else
            A.set_at(position, zero);
    }

    A.RST_simple();

    cout << endl << endl;
    cout << A.complexity() << endl;

    return 0;
}