#include "MultiDimArray.cpp"
#include "NTL/ZZ_p.h"

using namespace blitz;

int main() {
// TEST #1
    NTL::ZZ p = (NTL::ZZ) 11;
    NTL::ZZ_p::init(p);

    typedef NTL::ZZ_p F;
    const long m = 2;

    blitz::Array<F,m> A(2,2);
    A = (F)3, (F)10,
        (F)1, (F)8;
    MultiDimArray<F,m> array(A);
    array.RST();



// TEST #2
//    NTL::ZZ p = (NTL::ZZ) 2;
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const long m = 2;
//    blitz::Array<F,m> A(6,7);
//    A = (F)0, (F)0, (F)1, (F)1, (F)0, (F)1, (F)0,
//        (F)1, (F)0, (F)0, (F)0, (F)1, (F)1, (F)0,
//        (F)0, (F)0, (F)0, (F)1, (F)1, (F)0, (F)1,
//        (F)1, (F)1, (F)0, (F)1, (F)0, (F)0, (F)0,
//        (F)0, (F)1, (F)0, (F)0, (F)0, (F)1, (F)1,
//        (F)1, (F)0, (F)1, (F)0, (F)0, (F)0, (F)1;
//    MultiDimArray<F,m> array(A);
//    array.RST();


// TEST #3
//    NTL::ZZ p = (NTL::ZZ) 11;
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const long m = 2;
//    blitz::Array<F,m> A(5,5);
//    A = (F)4, (F)5, (F)7, (F)0, (F)8,
//        (F)8, (F)9, (F)0, (F)4, (F)1,
//        (F)0, (F)1, (F)3, (F)7, (F)0,
//        (F)8, (F)9, (F)0, (F)0, (F)0,
//        (F)10, (F)0, (F)0, (F)0, (F)0;
//    MultiDimArray<F,m> array(A);
//    array.RST();

    return 0;
}


// g++ -g -std=c++11 -pthread -march=native main.cpp  -lntl -lblitz -lgmp -lm
