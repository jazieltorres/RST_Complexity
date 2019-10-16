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

long powerMod(const long& x, const long& n, const long& mod){
    if (n == 0) return 1;
    long result = x;
    for (long i=1; i<n; i++){
        result = (result * x) % mod;
    }
    return result;
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
        return powerMod(root, i, mod);
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

long noShiftFunc(const long& i){
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



// HARD-CODE LEGENDRE AND COSTAS SEQUENCES
//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const unsigned long dim = 2;
//
//    long p_legendre = 5;
//    long p_costas = 19;
//    long root = 13;
//    cout << "Legendre:" << endl;
//    for (long i=0; i<p_legendre; i++) cout << LegendreSeq(p_legendre)(i) << " "; cout << endl;
//
//    cout << "Costas:" << endl;
//    for (long i=0; i<p_costas-1; i++) cout << CostasSeq(p_costas,root)(i) << " "; cout << "\n" << endl;
//    MultiDimArray<F,dim> Array(CostasSeq(p_costas,root), LegendreSeq(p_legendre), p_costas-1, p_legendre);
//    Array.RST();
//
//    unsigned long d = Array.getDeltaSize();
//    long n1 = p_costas - 1;
//    long criteria = p_legendre % 8;
//
//    if (criteria == 1 && ((p_legendre - 1) / 2 * n1 != d)) {
//        cout << "Failed case 1\t" << (p_legendre - 1) / 2 * n1 << " : " << d << endl;
//        cout << "p legendre: " << p_legendre << endl;
//        cout << "p costas: " << p_costas << endl;
//        cout << "root: " << root << endl;
//    }
//    if (criteria == 3 && (n1 * (p_legendre - 1) + 1 != d)) {
//        cout << "Failed case 2\t" << n1 * (p_legendre - 1) + 1 << " : " << d << endl;
//        cout << "p legendre: " << p_legendre << endl;
//        cout << "p costas: " << p_costas << endl;
//        cout << "root: " << root << endl;
//    }
//    if (criteria == 5 && (n1 * (p_legendre - 1) != d)) {
//        cout << "Failed case 3\t" << n1 * (p_legendre - 1) << " : " << d << endl;
//        cout << "p legendre: " << p_legendre << endl;
//        cout << "p costas: " << p_costas << endl;
//        cout << "root: " << root << endl;
//    }
//    if (criteria == 7 && ((p_legendre - 1) / 2 * n1 + 1 != d)) {
//        cout << "Failed case 4\t" << (p_legendre - 1) / 2 * n1 + 1 << " : " << d << endl;
//        cout << "p legendre: " << p_legendre << endl;
//        cout << "p costas: " << p_costas << endl;
//        cout << "root: " << root << endl;
//    }
// ^ HARD-CODE SEQUENCES ^



// LEGENDRE AND COSTAS SEQUENCES WITH RANDOM PRIMES
//    long numTest = 5;
//    long p_legendreMAX = 10;
//    long p_costasMAX = 10;
//
//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const unsigned long dim = 2;
//
//    double testSatisfied = 0;
//    srand(time(NULL));
//    for (long i=1; i<=numTest; i++) {
//        cout << "Test " << i << endl;
//
//        long p_legendre = rand() % p_legendreMAX;
//        while (!isPrime(p_legendre)) {p_legendre++;}
////        cout << "p legendre: " << p_legendre << endl;
//        long p_costas = rand() % p_costasMAX;
//        while (!isPrime(p_costas)) { p_costas++; }
////        cout << "p costas: " << p_costas << endl;
//        long root = rand() % p_costas;
//        while (!isRoot(root, p_costas)) { root = rand() % p_costas; }
////        cout << "root: " << root << endl;
//        MultiDimArray<F, dim> Array(CostasSeq(p_costas, root), LegendreSeq(p_legendre),
//                p_costas - 1, p_legendre);
//        Array.RST();
//
//        unsigned long d = Array.getDeltaSize();
//        long n1 = p_costas - 1;
//        long criteria = p_legendre % 8;
//        bool satisfied = true;
//
//        if (criteria == 1 && ((p_legendre - 1) / 2 * n1 != d)) {
//            cout << "p legendre: " << p_legendre << endl;
//            cout << "p costas: " << p_costas << endl;
//            cout << "root: " << root << endl;
//            cout << "Failed case 1\t" << (p_legendre - 1) / 2 * n1 << " : " << d << endl << endl;
//            satisfied = false;
//        }
//        if (criteria == 3 && (n1 * (p_legendre - 1) + 1 != d)) {
//            cout << "p legendre: " << p_legendre << endl;
//            cout << "p costas: " << p_costas << endl;
//            cout << "root: " << root << endl;
//            cout << "Failed case 2\t" << n1 * (p_legendre - 1) + 1 << " : " << d << endl << endl;
//            satisfied = false;
//        }
//        if (criteria == 5 && (n1 * (p_legendre - 1) != d)) {
//            cout << "p legendre: " << p_legendre << endl;
//            cout << "p costas: " << p_costas << endl;
//            cout << "root: " << root << endl;
//            cout << "Failed case 3\t" << n1 * (p_legendre - 1) << " : " << d << endl << endl;
//            satisfied = false;
//        }
//        if (criteria == 7 && ((p_legendre - 1) / 2 * n1 + 1 != d)) {
//            cout << "p legendre: " << p_legendre << endl;
//            cout << "p costas: " << p_costas << endl;
//            cout << "root: " << root << endl;
//            cout << "Failed case 4\t" << (p_legendre - 1) / 2 * n1 + 1 << " : " << d << endl << endl;
//            satisfied = false;
//        }
//
//        if (satisfied) {
//            testSatisfied = testSatisfied + 1;
//        }
//    } // end for
//    cout << "\nProportion of successful tests: " << testSatisfied/numTest << endl;
// ^ SEQUENCES WITH RANDOM PRIMES ^




// LEGENDRE AND COSTAS SEQUENCES WITH PRIMES LESS THAN primesUpTo
//    long primesUpTo = 20;
//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const unsigned long dim = 2;
//
//    vector<long> primes;
//    for (long i=3; i<=primesUpTo; i++){
//        if (isPrime(i)) primes.push_back(i);
//    }
//
//    double testSatisfied = 0;
//    long numTest = 0;
//    ofstream outfile ("failCases.txt", ios::app) ;
//    for(long p_costas : primes) {
//        vector<long> roots;
//        for(long c=2; c<p_costas; c++){
//            if(isRoot(c,p_costas)) roots.push_back(c);
//        }
//        for (long p_legendre : primes) {
//            for(auto root:roots) {
//                numTest++;
//                cout << "Test " << numTest << endl;
//                MultiDimArray<F, dim> Array(CostasSeq(p_costas, root), LegendreSeq(p_legendre),
//                        p_costas - 1, p_legendre);
//                Array.RST();
//
//                unsigned long d = Array.getDeltaSize();
//                long n1 = p_costas - 1;
//                long criteria = p_legendre % 8;
//                bool satisfied = true;
//
//                if (criteria == 1 && ((p_legendre - 1) / 2 * n1 != d)) {
//                    cout << "p legendre: " << p_legendre << endl;
//                    cout << "p costas: " << p_costas << endl;
//                    cout << "root: " << root << endl;
//                    cout << "Failed case 1\t" << (p_legendre - 1) / 2 * n1 << " : " << d << endl << endl;
//                    outfile << p_legendre << "\t" << p_costas << "\t" << root << '\t'
//                        << (p_legendre - 1) / 2 * n1 - d << "\n";
//                    satisfied = false;
//                }
//                if (criteria == 3 && (n1 * (p_legendre - 1) + 1 != d)) {
//                    cout << "p legendre: " << p_legendre << endl;
//                    cout << "p costas: " << p_costas << endl;
//                    cout << "root: " << root << endl;
//                    cout << "Failed case 2\t" << n1 * (p_legendre - 1) + 1 << " : " << d << endl << endl;
//                    outfile << p_legendre << "\t" << p_costas << "\t" << root << '\t'
//                        << (n1 * (p_legendre - 1) + 1 - d) << "\n";
//                    satisfied = false;
//                }
//                if (criteria == 5 && (n1 * (p_legendre - 1) != d)) {
//                    cout << "p legendre: " << p_legendre << endl;
//                    cout << "p costas: " << p_costas << endl;
//                    cout << "root: " << root << endl;
//                    cout << "Failed case 3\t" << n1 * (p_legendre - 1) << " : " << d << endl << endl;
//                    outfile << p_legendre << "\t" << p_costas << "\t" << root << '\t'
//                        << n1 * (p_legendre - 1) - d << "\n";
//                    satisfied = false;
//                }
//                if (criteria == 7 && ((p_legendre - 1) / 2 * n1 + 1 != d)) {
//                    cout << "p legendre: " << p_legendre << endl;
//                    cout << "p costas: " << p_costas << endl;
//                    cout << "root: " << root << endl;
//                    cout << "Failed case 4\t" << (p_legendre - 1) / 2 * n1 + 1 << " : " << d << endl << endl;
//                    outfile << p_legendre << "\t" << p_costas << "\t" << root << '\t'
//                        << (p_legendre - 1) / 2 * n1 + 1 - d << "\n";
//                    satisfied = false;
//                }
//
//                if (satisfied) {
//                    testSatisfied = testSatisfied + 1;
//                }
//            }
//        }
//    }
//
//    outfile.close();
//    cout << "\nProportion of successful tests: " << testSatisfied/numTest << endl;
// ^ SEQUENCES WITH PRIMES LESS THAN ^



    return 0;
}


// g++ -g -std=c++11 -pthread -march=native main.cpp  -lntl -lblitz -lgmp -lm
