#include "MultiDimArray.cpp"
//#include "Test.cpp"
#include "NTL/ZZ_p.h"
#include <cmath>
#include <chrono>
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

struct LogWelch {
    long mod;
    long root;
    vector<long> log;
    explicit LogWelch(long m, long r){
        mod = m;
        root = r;
        log.resize(m);
        for (long i=0; i<m-1; i++){
            log[powerMod(r,i,m)-1] = i;
        }
    }
    long operator () (const long& i) const {
            return log[i];
    }
};

long LogRoot(long i) {
    return i;
};

template <typename Func>
void printSeq(Func func1, long period){
    for (long i = 0; i < period; i++) {
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

long LegendreComplexity(long& p) {
    long criteria = p % 8;
    if(criteria == 1) return (p-1)/2;
    if(criteria == 3) return p;
    if(criteria == 5) return p-1;
    if(criteria == 7) return (p+1)/2;
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
//    long p_legendre =   7;
//    long p_costas =     7;
//    long root =         2;
//    cout << "Legendre: ";
//    for (long i=0; i<p_legendre; i++) cout << LegendreSeq(p_legendre)(i) << " "; cout << endl;
//
////    cout << "Log-Welch:" << endl;
////    for (long i=0; i<p_costas-1; i++) cout << LogWelch(p_costas,root)(i) << " "; cout << "\n" << endl;
//    MultiDimArray<F,dim> A(LogRoot, LegendreSeq(p_legendre), p_costas-1, p_legendre);
//    A.RST();
//
//    unsigned long d = A.getDeltaSize();
//    long n1 = A.period[0];
//    long criteria = p_legendre % 8;
//    long expected = 0;
//    bool satisfied = true;
//
//    if (criteria == 1 && ((p_legendre - 1) / 2 * n1 != d)) {
//        satisfied = false;
//        expected = (p_legendre - 1) / 2 * n1;
//    }
//    if (criteria == 3 && (n1 * (p_legendre - 1) + 1 != d)) {
//        satisfied = false;
//        expected = n1 * (p_legendre - 1) + 1;
//    }
//    if (criteria == 5 && (n1 * (p_legendre - 1) != d)) {
//        satisfied = false;
//        expected = n1 * (p_legendre - 1);
//    }
//    if (criteria == 7 && ((p_legendre - 1) / 2 * n1 + 1 != d)) {
//        satisfied = false;
//        expected = (p_legendre - 1) / 2 * n1 + 1;
//    }
//
//    if (!satisfied) {
//        cout << "Failed case " << criteria/2 + 1 << "\t" << expected << " : " << d << endl;
//        cout << "p legendre: " << p_legendre << endl;
//        cout << "p costas: " << p_costas << endl;
//        cout << "root: " << root << endl;
//
////        ofstream outfile ("failCases.txt", ios::app) ;
////        outfile << p_legendre << '\t' << p_costas << '\t' << root << '\t' << expected - d << "\t\t"
////            << "Case " << criteria/2 + 1 << '\t' << expected << " : " << d << '\n';
////        outfile.close();
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
    long primesFrom = 3;
    long primesUpTo = 37;
    NTL::ZZ p(2);
    NTL::ZZ_p::init(p);
    typedef NTL::ZZ_p F;
    const unsigned long dim = 2;

    vector<long> primes;
    for (long i=primesFrom; i<=primesUpTo; i++){
        if (isPrime(i)) primes.push_back(i);
    }

    double testSatisfied = 0;
    long numTest = 0;
    for(long p_costas : primes) {
        long root;
        for(long c=2; c<p_costas; c++){
            if(isRoot(c,p_costas)) root = c;
        }
        for (long p_legendre : primes) {
            numTest++;
            cout << "Test " << numTest << " of " << primes.size() * primes.size() << endl;
            MultiDimArray<F, dim> A(LogRoot, LegendreSeq(p_legendre),
                    p_costas - 1, p_legendre);
            A.RST();

            unsigned long d = A.getDeltaSize();
            long n1 = A.period[0];
            long criteria = p_legendre % 8;
            long expected = 0;
            bool satisfied = true;

            if (criteria == 1 && ((p_legendre - 1) / 2 * n1 != d)) {
                satisfied = false;
                expected = (p_legendre - 1) / 2 * n1;
            }
            if (criteria == 3 && (n1 * (p_legendre - 1) + 1 != d)) {
                satisfied = false;
                expected = n1 * (p_legendre - 1) + 1;
            }
            if (criteria == 5 && (n1 * (p_legendre - 1) != d)) {
                satisfied = false;
                expected = n1 * (p_legendre - 1);
            }
            if (criteria == 7 && ((p_legendre - 1) / 2 * n1 + 1 != d)) {
                satisfied = false;
                expected = (p_legendre - 1) / 2 * n1 + 1;
            }

            if (!satisfied) {
                cout << "Failed case " << criteria/2 + 1 << "\t" << expected << " : " << d << endl;
                cout << "p legendre: " << p_legendre << endl;
                cout << "p costas: " << p_costas << endl;
                cout << "root: " << root << endl;

                ofstream outfile ("failCases.txt", ios::app) ;
                outfile << p_legendre << '\t' << p_costas << '\t' << root << '\t' << expected - d << "\t\t"
                    << "Case " << criteria/2 + 1 << '\t' << expected << " : " << d << '\n';
                outfile.close();
            }
            else testSatisfied += 1;
        }
    }
    cout << "\nProportion of successful tests: " << testSatisfied/numTest << endl;
// ^ SEQUENCES WITH PRIMES LESS THAN ^



// COMPLEXITY OF LEGENDRE
//    long primesUpTo = 1000;
//
//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//
//    vector<long> primes;
//    for (long i=770; i<=primesUpTo; i++) {
//        if(isPrime(i)) primes.push_back(i);
//    }
//    double testSatisfied = 0;
//    int numTest = 0;
//    auto start = chrono::high_resolution_clock::now();
//    for(long p_legendre : primes) {
//        numTest++;
//        cout << "Test " << numTest << " of " << primes.size() << endl;
//        MultiDimArray<F,2> A(noShiftFunc, LegendreSeq(p_legendre), 1, p_legendre);
//        A.RST();
//        long d = A.getDeltaSize();
//        if(d != LegendreComplexity(p_legendre)) {
//            cout << "Failed with " << p_legendre << endl;
//        }
//        else testSatisfied = testSatisfied+1;
//    }
//    auto stop = chrono::high_resolution_clock::now();
//    auto duration = chrono::duration_cast<chrono::seconds>(stop - start);
//    cout << "\nProportion of successful tests: " << testSatisfied/numTest << endl;
//    cout << duration.count() << " ms" << endl;
// ^ COMPLEXITY OF LEGENDRE ^


// TEST FOR EACH ROOT
//    long primesFrom = 3;
//    long primesUpTo = 20;
//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const unsigned long dim = 2;
//
//    vector<long> primes;
//    for (long i=primesFrom; i<=primesUpTo; i++){
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
//
//        for (long p_legendre : primes) {
//            numTest++;
//            cout << "Test " << numTest << " of " << primes.size() * primes.size() << endl;
//
//            vector<bool> check_root;
//            for (long root : roots) {
//                MultiDimArray<F, dim> A(CostasSeq(p_costas, root), LegendreSeq(p_legendre),
//                                            p_costas - 1, p_legendre);
//                A.RST();
//
//                unsigned long d = A.getDeltaSize();
//                long n1 = p_costas - 1;
//                long criteria = p_legendre % 8;
//                bool satisfied = true;
//
//                if (criteria == 1 && ((p_legendre - 1) / 2 * n1 != d)) {
//                    satisfied = false;
//                }
//                if (criteria == 3 && (n1 * (p_legendre - 1) + 1 != d)) {
//                    satisfied = false;
//                }
//                if (criteria == 5 && (n1 * (p_legendre - 1) != d)) {
//                    satisfied = false;
//                }
//                if (criteria == 7 && ((p_legendre - 1) / 2 * n1 + 1 != d)) {
//                    satisfied = false;
//                }
//                check_root.push_back(satisfied);
//            }
//
//            bool check = true;
//            for (long i=1; i<check_root.size(); i++){
//                if (check_root[0] != check_root[i]) check = false;
//            }
//            if(check) testSatisfied = testSatisfied + 1;
//            else {
//                outfile << "p_legendre: " << p_legendre << endl;
//                outfile << "p_costas: " << p_costas << endl << endl;
//            }
//        }
//    }
//
//    outfile.close();
//    cout << "\nProportion of successful tests: " << testSatisfied/numTest << endl;
// ^ TEST FOR EACH ROOT ^




// PERMUTATIONS OF 1- -N1-2 FOR SHIFT
//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const unsigned long dim = 2;
//
//    long p_legendre =   5;
//    long p_costas =     7;
//    vector<long> shift(p_costas-1);
//    for (long i=0; i<p_costas; i++){
//        shift[i] = i;
//    }
//
//    long total=1, test=0;
//    for (long i=2; i<p_costas; i++){
//        total = total*i;
//    }
//    double testSatisfied = 0;
//
//    do {
//        test++;
//        cout << "Test " << test << " of " << total << endl;
//        MultiDimArray<F, dim> Array(shift, LegendreSeq(p_legendre), p_costas - 1, p_legendre);
//        Array.RST();
//
//        unsigned long d = Array.getDeltaSize();
//        long n1 = p_costas - 1;
//        long criteria = p_legendre % 8;
//        long expected = 0;
//        bool satisfied = true;
//
//        if (criteria == 1 && ((p_legendre - 1) / 2 * n1 != d)) {
//            satisfied = false;
//            expected = (p_legendre - 1) / 2 * n1;
//        }
//        if (criteria == 3 && (n1 * (p_legendre - 1) + 1 != d)) {
//            satisfied = false;
//            expected = n1 * (p_legendre - 1) + 1;
//        }
//        if (criteria == 5 && (n1 * (p_legendre - 1) != d)) {
//            satisfied = false;
//            expected = n1 * (p_legendre - 1);
//        }
//        if (criteria == 7 && ((p_legendre - 1) / 2 * n1 + 1 != d)) {
//            satisfied = false;
//            expected = (p_legendre - 1) / 2 * n1 + 1;
//        }
//
//        if(satisfied) testSatisfied += 1;
////        if (!satisfied) {
////            ofstream outfile("failCases.txt", ios::app);
////            outfile << expected - d << "\t\t" << "Case " << criteria/2 + 1 << '\t' << expected << " : " << d << '\n';
////            for (auto s : shift) outfile << s << " ";
////            outfile << endl;
////            outfile.close();
////        }
////        else {
////            testSatisfied += 1;
////        }
//    } while ( next_permutation(shift.begin(), shift.end()) );
//
//    cout << "Proportion " << testSatisfied/total;
// ^ PERMUTATIONS ^



    return 0;

}


// g++ -g -std=c++11 -pthread -march=native main.cpp  -lntl -lblitz -lgmp -lm
