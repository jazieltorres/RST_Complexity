#include "MultiDimArray.cpp"
#include "NTL/ZZ_p.h"

using namespace blitz;


//bool LegendreSeq(unsigned long& i, unsigned long& p) {
//    if ((i % p) == 0)
//        return 0;
//    else
//        return ((unsigned long)(pow(i, (p-1)/2)+1) % p)/2;
//}

struct LegendreSeq {
    long mod;
    explicit LegendreSeq(long p) {
        mod = p;
    }
    NTL::ZZ_p operator () (const long& i) const {
        bool result;
        if ((i % mod) == 0)
            result = 0;
        else
            result = ((long)(pow(i, (mod-1)/2)+1) % mod)/2;
        return (NTL::ZZ_p) result;
    }
};

struct CostasSeq {
    long mod;
    long root;
    explicit CostasSeq(long m, long r){
        mod = m; root = r;
    }
    long operator () (const long& i) const {
        return (long)pow(root, i) % mod;
    }
};

template <typename Func>
void printSeq(Func func1){
    for (long i = 0; i < 6; i++) {
        cout << func1(i) << " ";
    }
    cout << endl;
}






int main() {
// TEST #1
//    NTL::ZZ p = (NTL::ZZ) 11;
//    NTL::ZZ_p::init(p);
//
//    typedef NTL::ZZ_p F;
//    const long m = 2;
//
//    blitz::Array<F,m> A(2,2);
//    A = (F)3, (F)10,
//        (F)1, (F)8;
//    MultiDimArray<F,m> array(A);
//    array.RST();



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

    NTL::ZZ p(2);
    NTL::ZZ_p::init(p);
    typedef NTL::ZZ_p F;
    const unsigned long dim = 2;
    MultiDimArray<F,dim> Array(CostasSeq(7,3), LegendreSeq(7), 6, 7);
    Array.RST();

    return 0;
}


// g++ -g -std=c++11 -pthread -march=native main.cpp  -lntl -lblitz -lgmp -lm
