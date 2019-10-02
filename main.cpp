#include "MultiDimArray.cpp"
#include "NTL/ZZ_p.h"

using namespace blitz;

int main() {
//    NTL::ZZ p = (NTL::ZZ) 11;
//    NTL::ZZ_p::init(p);
//
//    typedef NTL::ZZ_p F;
//    const long m = 2;
//
//    blitz::Array<F,m> A(2,2);
//    A = (F)3, (F)10,
//        (F)1, (F)8;

    NTL::ZZ p = (NTL::ZZ) 2;
    NTL::ZZ_p::init(p);
    typedef NTL::ZZ_p F;
    const long m = 2;
//    blitz::TinyVector<unsigned long, 2> v1, v2;
//    v1 = 0,2;
//    v2 = 1,1;
//    MExponent<2> e1(v1), e2(v2);
//    cout << e1 << " < " << e2 << " : " << e1.grlex_less(e2) << endl;
//    cout <<"SUM: " << sum(v1) << endl;
    blitz::Array<F,m> A(6,7);
    A = (F)0, (F)0, (F)1, (F)1, (F)0, (F)1, (F)0,
        (F)1, (F)0, (F)0, (F)0, (F)1, (F)1, (F)0,
        (F)0, (F)0, (F)0, (F)1, (F)1, (F)0, (F)1,
        (F)1, (F)1, (F)0, (F)1, (F)0, (F)0, (F)0,
        (F)0, (F)1, (F)0, (F)0, (F)0, (F)1, (F)1,
        (F)1, (F)0, (F)1, (F)0, (F)0, (F)0, (F)1;


    MultiDimArray<F,m> array(A);
    array.RST();

    return 0;
}


// g++ -g -std=c++11 -pthread -march=native main.cpp  -lntl -lblitz -lgmp -lm
