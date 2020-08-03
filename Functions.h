
#ifndef RST_FUNCTIONS_H
#define RST_FUNCTIONS_H

#include <vector>
#include "NTL/ZZ_p.h"

using namespace std;

/******************************************************
*
*       CHECKS
*
*******************************************************/

bool isQuadratic(const long& n, const long& p);
bool isPrime(long& n);
bool isRoot(long& a, long& p);
bool isHere(vector<long> differences, long diff);
bool isCostas(vector<long>& shift);


/******************************************************
*
*       MISCELLANY
*
*******************************************************/

unsigned long readDim();
bool compare(vector<long>& v1, vector<long>& v2);
long powerMod(const long& x, const long& n, const long& mod);
long f(const long& x, const long& mod, long& a, long& b);
long LegendreComplexity(long& p);


/******************************************************
*
*       SEQUENCES
*
*******************************************************/


struct MSeq {
    vector<long> v;
    MSeq() {
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


struct ExpQuadraticSeq {
    long mod;
    long root;
    long a;
    explicit ExpQuadraticSeq(long p, long r, long A) {
        mod = p; root = r; a = A;
    }
    long operator () (const long& i) {
        return (a * powerMod(root*root, i, mod) + powerMod(root, i, mod)) % mod;
    }
};

struct LogQuadraticSeq {
    long mod;
    long root;
    vector<long> logTable;
    explicit LogQuadraticSeq(long p, long r) {
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

struct LosWelchSeq {
    long mod;
    long root;
    vector<long> log;
    explicit LosWelchSeq(long m, long r){
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


template <typename Func>
void printSeq(Func func1, long period){
    for (long i = 0; i < period; i++) {
        cout << func1(i) << " ";
    }
    cout << endl;
}

long LogRootSeq(long i);
long noShiftFunc(const long& i);

#endif //RST_FUNCTIONS_H
