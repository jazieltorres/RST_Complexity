#include "MultiDimArray.cpp"
//#include "Test.cpp"
#include "NTL/ZZ_p.h"
#include <cmath>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>


using namespace blitz;
using namespace std;

unsigned long readDim() {
    unsigned long dim;
    cin >> dim;
    return dim;
}

bool compare(vector<long>& v1, vector<long>& v2){
    return lexicographical_compare(v1.begin(), v1.end(), v2.begin(), v2.end());
}


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

struct M_Seq {
    vector<long> v;
    M_Seq() {
        v = {0,0,1,0,1,1,1};
    }
    NTL::ZZ_p operator () (const long& i) const {
        return (NTL::ZZ_p) v[i];
    }
};

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

// This is the polynomial (degree 3) that generates the PolynomialSeq
long f(const long& x, const long& mod, long& a, long& b) {
    return( a * powerMod(x, 3, mod) +
            b * powerMod(x, 2, mod) +
            x )
            % mod;
}

struct PolynomialSeq {
    long m;
    long alpha;
    long A;
    long B;
    explicit PolynomialSeq(long mod, long root, long a, long b) {
        m = mod; alpha = root; A = a; B = b;
    }
    long operator () (const long& i) {
        return f(powerMod(alpha, i, m), m, A, B);
    }
};


struct ExpQuadratic {
    long mod;
    long root;
    long a;
    explicit ExpQuadratic(long p, long r, long A) {
        mod = p; root = r; a = A;
    }
    long operator () (const long& i) {
        return (a * powerMod(root*root, i, mod) + powerMod(root, i, mod)) % mod;
    }
};

struct LogQuadratic {
    long mod;
    long root;
    vector<long> logTable;
    explicit LogQuadratic(long p, long r) {
        mod = p; root = r;
        logTable.resize(p);
        logTable[0] = -1;
        for (long i=0; i<mod-1; i++) {
            logTable[powerMod(root, i, mod)] = i;
        }
    }
    long operator () (const long& i) {
        long result = (4*powerMod(root*root, i, mod) + powerMod(root, i, mod)) % mod;
        return logTable[result];
    }
};


struct ExponentialSeq {
    long mod;
    long root;
    explicit ExponentialSeq(long m, long r){
        mod = m; root = r;
    }
    long operator () (const long& i) const {
        return powerMod(root, i, mod);
    }
};

struct OneSequence {
    long n;
    explicit OneSequence(long length) {
        n = length;
    }
    NTL::ZZ_p operator () (const long& i) const {
        if(i==0)    return (NTL::ZZ_p) 1;
        else        return (NTL::ZZ_p) 0;
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

bool isHere(vector<long> differences, long diff){
    for (auto d : differences){
        if (d==diff) return true;
    }
    return false;
}


bool isCostas(vector<long>& shift) {
    for (long i=1; i<shift.size()-1; i++){
        vector<long> differences;
        long first_index(0);
        while (first_index + i < shift.size()) {
            long diff = shift[first_index + i] - shift[first_index];
            if (isHere(differences, diff)) return false;
            else {
                differences.push_back(diff);
            }
            first_index++;
        }
    }
    return true;
}





int main() {
// TEST #1
    NTL::ZZ p(11);
    NTL::ZZ_p::init(p);

    typedef NTL::ZZ_p F;
    const long m = 2;

    blitz::Array<F,m> A(2,2);
//    A = (F)3, (F)10,
//        (F)1, (F)8;
//    MultiDimArray<F,m> array(A);
//    array.RST();

    return 0;



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



// HARD-CODE SEQUENCES
//    cout << "Started" << endl;
//
//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const unsigned long dim = 2;
//
//    long p_legendre =    11;
//    long p_shift =       11;
//    long root =          5;
//
////    For the vector of shifts
////    vector<long> shift({0, 1, 3, 1, 0, 3});
////    long p_shift = shift.size();
//
//    PolynomialSeq Poly(p_shift, root, 2, 3);
//
//    for (long i=0; i<p_shift-1; i++)
//        cout << Poly(i) << " ";
//    cout << endl;
////    return 0;
//
////    cout << "Legendre: ";
////    for (long i=0; i<p_legendre; i++) cout << LegendreSeq(p_legendre)(i) << " "; cout << endl << endl;
////    cout << "Log-Welch:" << endl;
////    for (long i=0; i<p_shift-1; i++) cout << LogWelch(p_shift,root)(i) << " "; cout << "\n" << endl;
////    cout << "Quadratic:" << endl;
////    for (long i=0; i<p_shift-1; i++) cout << ExpQuadratic(p_shift, root)(i) << " "; cout << endl;
////    cout << "Polynomial Sequence: ";
////    for (long i=0; i<p_shift-1; i++) cout << PolynomialSeq(p_shift, root)(i) << " "; cout << endl << endl;
//
//
//    MultiDimArray<F,dim> A(Poly, LegendreSeq(p_legendre), p_shift, p_legendre);
//
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
//        cout << "p costas: " << p_shift << endl;
//        cout << "root: " << root << endl;

//        ofstream outfile ("failCases.txt", ios::app) ;
//        outfile << p_legendre << '\t' << p_shift << '\t' << root << '\t' << expected - d << "\t\t"
//            << "Case " << criteria/2 + 1 << '\t' << expected << " : " << d << '\n';
//        vector<bool> check(shift.size(), false);
//        for(long i=0; i<shift.size(); i++){
//            if(!check[i]){
//                outfile << "(" << i;
//                check[i] = true;
//                long next = shift[i];
//                while (next != i){
//                    outfile << ", " << next;
//                    check[next] = true;
//                    next = shift[next];
//                }
//                outfile << ")";
//            }
//        }
//        outfile << endl;
//        outfile.close();
//    }
// ^ HARD-CODE SEQUENCES ^



// LEGENDRE AND COSTAS SEQUENCES WITH RANDOM PRIMES
//    long numTest = 5;
//    long p_legendreMAX = 10;
//    long p_shiftMAX = 10;
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
//        long p_shift = rand() % p_shiftMAX;
//        while (!isPrime(p_shift)) { p_shift++; }
////        cout << "p costas: " << p_shift << endl;
//        long root = rand() % p_shift;
//        while (!isRoot(root, p_shift)) { root = rand() % p_shift; }
////        cout << "root: " << root << endl;
//        MultiDimArray<F, dim> Array(CostasSeq(p_shift, root), LegendreSeq(p_legendre),
//                p_shift - 1, p_legendre);
//        Array.RST();
//
//        unsigned long d = Array.getDeltaSize();
//        long n1 = p_shift - 1;
//        long criteria = p_legendre % 8;
//        bool satisfied = true;
//
//        if (criteria == 1 && ((p_legendre - 1) / 2 * n1 != d)) {
//            cout << "p legendre: " << p_legendre << endl;
//            cout << "p costas: " << p_shift << endl;
//            cout << "root: " << root << endl;
//            cout << "Failed case 1\t" << (p_legendre - 1) / 2 * n1 << " : " << d << endl << endl;
//            satisfied = false;
//        }
//        if (criteria == 3 && (n1 * (p_legendre - 1) + 1 != d)) {
//            cout << "p legendre: " << p_legendre << endl;
//            cout << "p costas: " << p_shift << endl;
//            cout << "root: " << root << endl;
//            cout << "Failed case 2\t" << n1 * (p_legendre - 1) + 1 << " : " << d << endl << endl;
//            satisfied = false;
//        }
//        if (criteria == 5 && (n1 * (p_legendre - 1) != d)) {
//            cout << "p legendre: " << p_legendre << endl;
//            cout << "p costas: " << p_shift << endl;
//            cout << "root: " << root << endl;
//            cout << "Failed case 3\t" << n1 * (p_legendre - 1) << " : " << d << endl << endl;
//            satisfied = false;
//        }
//        if (criteria == 7 && ((p_legendre - 1) / 2 * n1 + 1 != d)) {
//            cout << "p legendre: " << p_legendre << endl;
//            cout << "p costas: " << p_shift << endl;
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
//    long primesFrom = 5;
//    long primesUpTo = 37;
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
//    for(long p_shift : primes) {
//        long root;
//        for(long c=2; c<p_shift; c++){
//            if(isRoot(c,p_shift)) root = c;
//        }
//        for (long p_legendre : primes) {
//            numTest++;
//            cout << "Test " << numTest << " of " << primes.size() * primes.size() << endl;
//            MultiDimArray<F, dim> A(ExpQuadratic(p_shift, root), LegendreSeq(p_legendre),
//                    p_shift - 1, p_legendre);
//            A.RST();
//
//            unsigned long d = A.getDeltaSize();
//            long n1 = A.period[0];
//            long criteria = p_legendre % 8;
//            long expected = 0;
//            bool satisfied = true;
//
//            if (criteria == 1 && ((p_legendre - 1) / 2 * n1 != d)) {
//                satisfied = false;
//                expected = (p_legendre - 1) / 2 * n1;
//            }
//            if (criteria == 3 && (n1 * (p_legendre - 1) + 1 != d)) {
//                satisfied = false;
//                expected = n1 * (p_legendre - 1) + 1;
//            }
//            if (criteria == 5 && (n1 * (p_legendre - 1) != d)) {
//                satisfied = false;
//                expected = n1 * (p_legendre - 1);
//            }
//            if (criteria == 7 && ((p_legendre - 1) / 2 * n1 + 1 != d)) {
//                satisfied = false;
//                expected = (p_legendre - 1) / 2 * n1 + 1;
//            }
//
//            if (!satisfied) {
//                cout << "Failed case " << criteria/2 + 1 << "\t" << expected << " : " << d << endl;
//                cout << "p legendre: " << p_legendre << endl;
//                cout << "p costas: " << p_shift << endl;
//                cout << "root: " << root << endl;
//
//                ofstream outfile ("failCases.txt", ios::app) ;
//                outfile << p_legendre << '\t' << p_shift << '\t' << root << '\t' << expected - d << "\t\t"
//                    << "Case " << criteria/2 + 1 << '\t' << expected << " : " << d << '\n';
//                outfile.close();
//            }
//            else testSatisfied += 1;
//        }
//    }
//    cout << "\nProportion of successful tests: " << testSatisfied/numTest << endl;
// ^ SEQUENCES WITH PRIMES LESS THAN ^



// COMPLEXITY OF LEGENDRE
//    long primesUpTo = 80;
//
//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//
//    vector<long> primes;
//    for (long i=3; i<=primesUpTo; i++) {
//        if(isPrime(i)) primes.push_back(i);
//    }
//    double testSatisfied = 0;
//    int numTest = 0;
//    auto start = chrono::high_resolution_clock::now();
//    for(long p_legendre : primes) {
//        numTest++;
////        cout << "Test " << numTest << " of " << primes.size() << endl;
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
//    for(long p_shift : primes) {
//        vector<long> roots;
//        for(long c=2; c<p_shift; c++){
//            if(isRoot(c,p_shift)) roots.push_back(c);
//        }
//
//        for (long p_legendre : primes) {
//            numTest++;
//            cout << "Test " << numTest << " of " << primes.size() * primes.size() << endl;
//
//            vector<bool> check_root;
//            for (long root : roots) {
//                MultiDimArray<F, dim> A(ExponentialSeq(p_shift, root), LegendreSeq(p_legendre),
//                                            p_shift - 1, p_legendre);
//                A.RST();
//
//                unsigned long d = A.getDeltaSize();
//                long n1 = p_shift - 1;
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
//                outfile << "p_shift: " << p_shift << endl << endl;
//            }
//        }
//    }
//
//    outfile.close();
//    cout << "\nProportion of successful tests: " << testSatisfied/numTest << endl;
// ^ TEST FOR EACH ROOT ^




// PERMUTATIONS OF NUMBERS FROM 0 TO N1-2 FOR SHIFT
//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const unsigned long dim = 2;
//
//    long p_legendre =   5;
//    long shift_len =    7;
//    vector<long> shift(shift_len);
//    for (long i=0; i<shift_len; i++){
//        shift[i] = i;
//    }
//
//    long total=1, numTest=0;
//    for (long i=2; i<shift_len; i++){
//        total = total*i;
//    }
//    double testSatisfied = 0;
//
//    vector< vector<long> > permutations;
//
//    do {
//        numTest++;
////        cout << "Test " << numTest << " of " << total << endl;
//        MultiDimArray<F, dim> Array(shift, LegendreSeq(p_legendre), shift_len, p_legendre);
//        Array.RST();
//
//        unsigned long d = Array.getDeltaSize();
//        long n1 = shift_len;
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
//        else {
//            permutations.push_back(shift);
////            cout << d << "\t" << expected << endl;
//        }
////        if (!satisfied) {
////            ofstream outfile("failCases.txt", ios::app);
////            outfile << expected - d << "\t\t" << "Case " << criteria/2 + 1 << '\t' << expected << " : " << d << "\t\t";
////            vector<bool> check(shift.size(), false);
////            for(long i=0; i<shift.size(); i++){
////                if(!check[i]){
////                    outfile << "(" << i;
////                    check[i] = true;
////                    long next = shift[i];
////                    while (next != i){
////                        outfile << ", " << next;
////                        check[next] = true;
////                        next = shift[next];
////                    }
////                    outfile << ")";
////                }
////            }
////            outfile << "\t\t";
////            for (auto i : shift) outfile << i << " ";
////
////            outfile << "\t\t";
////            for (long i=0; i<shift.size()-1; i++){
////                outfile << shift[i+1]-shift[i] << " ";
////            }
////            outfile << endl;
////            outfile.close();
////            permutations.push_back(shift);
////        }
////        else {
////            testSatisfied += 1;
////        }
//    } while ( next_permutation(shift.begin(), shift.end()) && numTest<20);
//
//    for(auto&& v : permutations){
//        while (v[0] != 0){
//            v.push_back(v[0]);
//            v.erase(v.begin());
//        }
//    }
//    sort(permutations.begin(), permutations.end(), compare);
//
//    for (auto v : permutations) {
//        for (auto s:v){
//            cout << s << ", ";
//        }
//        cout << endl;
//    }
//
//    cout << "Proportion " << testSatisfied << "/" << numTest << " = " << testSatisfied/numTest;
// ^ PERMUTATIONS ^



// EXPONENTIAL QUADRATIC (VALUE)
//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const unsigned long dim = 2;
//
//    vector<long> primes;
//    for (long i=11; i<=11; i++) {
//        if(isPrime(i)) primes.push_back(i);
//    }
//
//    for(long p_legendre : primes) {
//        cout << " - - - - - - - - - - - - - - - - - - - - -" << endl;
//        cout << "PRIME: " << p_legendre << endl << endl;
//        long p_shift = p_legendre;
//        long root(2);
//        while (!isRoot(root, p_shift)) root++;
//
//        cout << "Legendre: ";
//        for (long i = 0; i < p_legendre; i++)
//            cout << LegendreSeq(p_legendre)(i) << " ";
//        cout << endl << endl;
//        cout << "Shift Sequence:" << endl;
//        ExpQuadratic ShiftSeq(p_shift, root);
//        for (long i = 0; i < p_shift - 1; i++)
//            cout << ShiftSeq(i) << " ";
//        cout << endl << endl;
//
//        MultiDimArray<F, dim> A(ShiftSeq, LegendreSeq(p_legendre), p_shift - 1, p_legendre);
//        A.RST();
//
//        unsigned long d = A.getDeltaSize();
//        cout << "Linear complexity\t" << d << endl;
//        double size = p_legendre * (p_shift - 1);
//        cout << "Normalized complexity\t" << d / size << endl;
//
//        long n1 = A.period[0];
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
//        if (!satisfied) {
//            cout << "Failed case " << criteria / 2 + 1 << "\t" << expected << " : " << d << endl;
//            cout << "p legendre: " << p_legendre << endl;
//            cout << "p costas: " << p_shift << endl;
//            cout << "root: " << root << endl;
//        }
//    }
// ^^ EXPONENTIAL QUADRATIC (VALUE) ^^


// EXPONENTIAL QUADRATIC (LOG)
//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const unsigned long dim = 2;
//
//    vector<long> primes;
//    for (long i=7; i<=19; i++) {
//        if(isPrime(i)) primes.push_back(i);
//    }
//
//    for(long p_legendre : primes) {
//        cout << " - - - - - - - - - - - - - - - - - - - - -" << endl;
//        cout << "PRIME: " << p_legendre << endl << endl;
//
//        long p_shift = p_legendre;
//        long root(2);
//        while (!isRoot(root, p_shift)) root++;
//
//        cout << "Legendre: ";
//        for (long i = 0; i < p_legendre; i++)
//            cout << LegendreSeq(p_legendre)(i) << " ";
//        cout << endl << endl;
//        cout << "Shift Sequence (Root " << root << ")" << endl;
//        LogQuadratic ShiftSeq(p_shift, root);
//        for (long i = 0; i < p_shift - 1; i++)
//            cout << ShiftSeq(i) << " ";
//        cout << endl << endl;
//
//        MultiDimArray<F, dim> A(ShiftSeq, LegendreSeq(p_legendre), p_shift - 1, p_legendre);
//        A.RST();
//
//        unsigned long d = A.getDeltaSize();
//        cout << "Linear complexity\t" << d << endl;
//        double size = p_legendre * (p_shift - 1);
//        cout << "Normalized complexity\t" << d / size << endl;
//
//        long n1 = A.period[0];
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
//        if (!satisfied) {
//            cout << "Failed case " << criteria / 2 + 1 << "\t" << expected << " : " << d << endl;
//            cout << "p legendre: " << p_legendre << endl;
//            cout << "p costas: " << p_shift << endl;
//            cout << "root: " << root << endl;
//        }
//    }





// ALL COMBINATION OF COEFFICIENTS IN POLYNOMIAL SHIFT SEQUENCE
//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const unsigned long dim = 2;
//
//    long p_legendre =    23;
//    long p_shift =       23;
//    long root =          5;
//
//    long counter = 1;
//    for (long a = 4; a < p_legendre; a++) {
//        cout << "Test " << counter << endl;
//        counter++;
//        for (long b = 1; b < p_legendre; b++) {
//            MultiDimArray<F, dim> A(PolynomialSeq(p_shift, root, a, b), LegendreSeq(p_legendre), p_shift - 1, p_legendre);
//
//            A.RST();
//
//            unsigned long d = A.getDeltaSize();
////            cout << d << endl;
//
//
//            long n1 = A.period[0];
//            long criteria = p_legendre % 8;
//            long expected = 0;
//            bool satisfied = true;
//
//            if (criteria == 1 && ((p_legendre - 1) / 2 * n1 != d)) {
//                satisfied = false;
//                expected = (p_legendre - 1) / 2 * n1;
//            }
//            if (criteria == 3 && (n1 * (p_legendre - 1) + 1 != d)) {
//                satisfied = false;
//                expected = n1 * (p_legendre - 1) + 1;
//            }
//            if (criteria == 5 && (n1 * (p_legendre - 1) != d)) {
//                satisfied = false;
//                expected = n1 * (p_legendre - 1);
//            }
//            if (criteria == 7 && ((p_legendre - 1) / 2 * n1 + 1 != d)) {
//                satisfied = false;
//                expected = (p_legendre - 1) / 2 * n1 + 1;
//            }
//
//            if (!satisfied) {
//                cout << a << " " << b << endl;
////                for (int i=0; i<p_shift-1; i++)
////                    cout << PolynomialSeq(p_shift, root, a, b)(i) << " ";
////                cout << endl;
//            }
//        }
//    }




//GENERALIZED LEGENDRE FORM PAPER MULTI-DIM ARRAYS FOR WATERMARKING (SECTION IIB)
// Files are generated in the Sage program Generalized_Legendre

//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const unsigned long dim = 2;
//
//
//    string line;
//    vector<long> v(dim);
//    long size = 1;
//    for (int i=0; i<dim; i++) {
//        cin >> v[i];
//        size = size * v[i];
//    }
//
////    cout << "Size " << size << endl;
//
////    Removing character in buffer after cin
//    getline(cin, line);
//
//    int ctr = 0;
//
//    while(ctr < 30 && getline(cin, line)) {
//        ctr++;
//        // Skipping following two lines
//        getline(cin, line); getline(cin, line);
//
//        MultiDimArray<F, dim> A(v);
//        for (long i = 0; i < size; i++) {
////            cout << i << endl;
//            blitz::TinyVector<int,dim> position;
//            F value;
//
//            getline(cin, line);
//            stringstream ss(line);
//            for (int n=0; n<dim; n++)
//                ss >> position[n];
//            ss >> value;
//            A.setAt(position, value);
//        }
////        A.print();
////        cout << endl;
//        A.RST();
////        cout << ctr << endl;
////        cout << "Delta size: " << A.getDeltaSize() << endl;
//    }
////    cout << "\n\nNext line in buffer:" << endl;
////    getline(cin, line);
////    cout << line << endl;










//GENERALIZED TERNARY LEGENDRE FORM PAPER MULTI-DIM ARRAYS FOR WATERMARKING (SECTION IIB)
// Files are generated in the Sage program Generalized_Legendre

//    NTL::ZZ p(3);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const unsigned long dim = 2;
//
//
//    string line;
//    vector<long> v(dim);
//    long size = 1;
//    for (int i=0; i<dim; i++) {
//        cin >> v[i];
//        size = size * v[i];
//    }
//
////    cout << "Size " << size << endl;
//
////    Removing character in buffer after cin
//    getline(cin, line);
//
//    int ctr = 0;
//
//    while(getline(cin, line)) {
//        ctr++;
//        // Skipping following two lines
//        getline(cin, line); getline(cin, line);
//
//        MultiDimArray<F, dim> A(v);
//        for (long i = 0; i < size; i++) {
////            cout << i << endl;
//            blitz::TinyVector<int,dim> position;
//            F value;
//
//            getline(cin, line);
//            stringstream ss(line);
//            for (int n=0; n<dim; n++)
//                ss >> position[n];
//            ss >> value;
//            A.setAt(position, value);
//        }
////        A.print();
////        cout << endl;
//        A.RST();
//        cout << ctr << endl;
//        cout << "Delta size: " << A.getDeltaSize() << endl;
//    }
////    cout << "\n\nNext line in buffer:" << endl;
////    getline(cin, line);
////    cout << line << endl;






//GENERALIZED 3D BINARY LEGENDRE FORM PAPER MULTI-DIM ARRAYS FOR WATERMARKING (SECTION III-B)
// Files are generated in the Sage program Generalized_Legendre

//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const unsigned long dim = 3;
//
//
//    string line;
//    vector<long> v(dim);
//    long size = 1;
//    for (int i=0; i<dim; i++) {
//        cin >> v[i];
//        size = size * v[i];
//    }
//
////    cout << "Size " << size << endl;
//
////    Removing character in buffer after cin
//    getline(cin, line);
//
//    int ctr = 0;
//
//    while(getline(cin, line)) {
//        ctr++;
//        // Skipping following two lines
//        getline(cin, line); getline(cin, line);
//
//        MultiDimArray<F, dim> A(v);
//        for (long i = 0; i < size; i++) {
////            cout << i << endl;
//            blitz::TinyVector<int,dim> position;
//            F value;
//
//            getline(cin, line);
//            stringstream ss(line);
//            for (int n=0; n<dim; n++)
//                ss >> position[n];
//            ss >> value;
//            A.setAt(position, value);
//        }
////        A.print();
////        cout << endl;
//        A.RST();
////        cout << ctr << endl;
//        cout << "Delta size: " << A.getDeltaSize() << '\t' << A.getDeltaSize()/(1.0 * size) << endl;
//    }
////    cout << "\n\nNext line in buffer:" << endl;
////    getline(cin, line);
////    cout << line << endl;




//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
//    const unsigned long dim = 3;
//
//
//    string line;
//    vector<long> v(dim);
//    long size = 1;
//    for (int i=0; i<dim; i++) {
//        cin >> v[i];
//        size = size * v[i];
//    }
//
////    cout << "Size " << size << endl;
//
////    Removing character in buffer after cin
//    getline(cin, line);
//
//    int ctr = 0;
//
//    while(getline(cin, line)) {
//        ctr++;
//        // Skipping following two lines
//        getline(cin, line); getline(cin, line);
//
//        MultiDimArray<F, dim> A(v);
//        for (long i = 0; i < size; i++) {
////            cout << i << endl;
//            blitz::TinyVector<int,dim> position;
//            F value;
//
//            getline(cin, line);
//            stringstream ss(line);
//            for (int n=0; n<dim; n++)
//                ss >> position[n];
//            ss >> value;
//            A.setAt(position, value);
//        }
////        A.print();
////        cout << endl;
//        A.RST();
////        cout << ctr << endl;
//        cout << "Delta size: " << A.getDeltaSize() << '\t' << A.getDeltaSize()/(1.0 * size) << endl;
//    }
////    cout << "\n\nNext line in buffer:" << endl;
////    getline(cin, line);
////    cout << line << endl;
//
//    return 0;
}


// g++ -g -std=c++11 -pthread -march=native main.cpp  -lntl -lblitz -lgmp -lm



