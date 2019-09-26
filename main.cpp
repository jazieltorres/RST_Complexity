#include <iostream>
#include "MultiDimArray.cpp"

using namespace blitz;



int main() {
    NTL::ZZ p = (NTL::ZZ) 11;
    NTL::ZZ_p::init(p);

    typedef NTL::ZZ_p F;
    const long m = 2;

    blitz::Array<F,m> A(2,2);
    A = (F)3, (F)10,
        (F)1, (F)8;
    MultiDimArray<F,m> array(A);
    array.RST();

    cout << "Compiles!" << endl;
    return 0;
}


// g++ -std=c++11 -pthread -march=native main.cpp  -lntl -lblitz -lgmp -lm
