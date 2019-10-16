#include "MultiDimArray.cpp"
//#include "Test.cpp"
#include "NTL/ZZ_p.h"
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>

using namespace blitz;
using namespace std;


bool isQuadratic(const long& n, const long& p){
    for (long i=0; i<p; i++){
        if ((i*i %p) == n) return true;
    }
    return false;
}

struct LegendreSeq {
    long mod;
    explicit LegendreSeq(long p) {
        mod = p;
    }
    NTL::ZZ_p operator () (const long& i) const {
        bool result;
        if ((i % mod) == 0)
            result = false;
        else
            result = isQuadratic(i,mod);
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

bool isPrime(long& n){
    if (n<3) return false;
    long bound = floor(sqrt(n));
    for (long i=2; i<=bound; i++){
        if (n%i == 0) return false;
    }
    return true;
}

bool isRoot(long& a, long& p) {
    long num = a*a;
    long ord = 1;
    while (num % p != a) {
        num = num*a % p;
        ord++;
    }
    if (ord == p-1) return true;
    return false;
}

long func(const long& i){
    return 0;
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
//    NTL::ZZ p(2);
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



// COMPOSITION METHOD TEST
//    long p_legendre = 17;
//    long p_costas = 11;
//    long root = 7;


//    long p_legendre = 29;
//    long p_costas = 3;
//    long root = 2;


    NTL::ZZ p(2);
    NTL::ZZ_p::init(p);
    typedef NTL::ZZ_p F;
    const unsigned long dim = 2;

// TO HARD CODE SEQUENCES
    long p_legendre = 5;
    long p_costas = 11;
    long root = 2;
    MultiDimArray<F,dim> Array(CostasSeq(p_costas,root), LegendreSeq(p_legendre), p_costas-1, p_legendre);
    Array.RST();
// ^ TO HARD CODE SEQUENCES ^

//    vector<long> primes;
//    for (long i=23; i<32; i++){
//        if (isPrime(i)) primes.push_back(i);
//    }
//
//    double testSatisfied = 0;
//    long numTest = 0;
//    srand(time(NULL));
//    ofstream outfile ("failCases.txt", ios::app) ;
//    for(long a=0; a<primes.size(); a++) {
//        long p_costas = primes[a];
//        vector<long> roots;
//        for(long c=2; c<p_costas; c++){
//            if(isRoot(c,p_costas)) roots.push_back(c);
//        }
//        for (long b=0; b<primes.size(); b++) {
//            long p_legendre = primes[b];
//
//// RANDOM PRIMES
////        long p_legendre = rand() % 3;
////        while (!isPrime(p_legendre)) {p_legendre++;}
////        cout << "p legendre: " << p_legendre << endl;
////        long p_costas = rand() % 30;
////        while (!isPrime(p_costas)) { p_costas++; }
////        cout << "p costas: " << p_costas << endl;
////        long root = rand() % p_costas;
////        while (!isRoot(root, p_costas)) { root = rand() % p_costas; }
////        cout << "root: " << root << endl;
//// ^ RANDOM PRIMES ^
//            for(auto root:roots) {
//                numTest++;
//                cout << "Test " << numTest << endl;
//                MultiDimArray<F, dim> Array(CostasSeq(p_costas, root), LegendreSeq(p_legendre), p_costas - 1,
//                                            p_legendre);
//                Array.RST();

                unsigned long d = Array.getDeltaSize();
                long n1 = p_costas - 1;
                long criteria = p_legendre % 8;
                bool satisfied = true;

                if (criteria == 1 && ((p_legendre - 1) / 2 * n1 != d)) {
                    cout << "p legendre: " << p_legendre << endl;
                    cout << "p costas: " << p_costas << endl;
                    cout << "root: " << root << endl;
                    cout << "Failed case 1\t" << (p_legendre - 1) / 2 * n1 << " : " << d << endl << endl;
//                    outfile << p_legendre << "\t" << p_costas << "\t" << root << '\t' << (p_legendre - 1) / 2 * n1 - d << "\n";
                    satisfied = false;
                }
                if (criteria == 3 && (n1 * (p_legendre - 1) + 1 != d)) {
                    cout << "p legendre: " << p_legendre << endl;
                    cout << "p costas: " << p_costas << endl;
                    cout << "root: " << root << endl;
                    cout << "Failed case 2\t" << n1 * (p_legendre - 1) + 1 << " : " << d << endl << endl;
//                    outfile << p_legendre << "\t" << p_costas << "\t" << root << '\t' << (n1 * (p_legendre - 1) + 1 - d) << "\n";
                    satisfied = false;
                }
                if (criteria == 5 && (n1 * (p_legendre - 1) != d)) {
                    cout << "p legendre: " << p_legendre << endl;
                    cout << "p costas: " << p_costas << endl;
                    cout << "root: " << root << endl;
                    cout << "Failed case 3\t" << n1 * (p_legendre - 1) << " : " << d << endl << endl;
//                    outfile << p_legendre << "\t" << p_costas << "\t" << root << '\t' << n1 * (p_legendre - 1) - d << "\n";
                    satisfied = false;
                }
                if (criteria == 7 && ((p_legendre - 1) / 2 * n1 + 1 != d)) {
                    cout << "p legendre: " << p_legendre << endl;
                    cout << "p costas: " << p_costas << endl;
                    cout << "root: " << root << endl;
                    cout << "Failed case 4\t" << (p_legendre - 1) / 2 * n1 + 1 << " : " << d << endl << endl;
//                    outfile << p_legendre << "\t" << p_costas << "\t" << root << '\t' << (p_legendre - 1) / 2 * n1 + 1 - d << "\n";
                    satisfied = false;
                }

//                if (satisfied) {
//                    testSatisfied = testSatisfied + 1;
//                }
            }
//        }
//    }

//    outfile.close();
//    cout << "\nProportion of successful tests: " << testSatisfied/numTest << endl;



    return 0;
}


// g++ -g -std=c++11 -pthread -march=native main.cpp  -lntl -lblitz -lgmp -lm
