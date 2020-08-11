
#ifndef RST_FUNCTIONS_H
#define RST_FUNCTIONS_H

#include <vector>
#include "NTL/ZZ_p.h"
#include <string>

using namespace std;

/******************************************************
*
*       CHECKS
*
*******************************************************/

bool isQuadratic(const int& n, const int& p);
bool isPrime(int& n);
bool isRoot(int& a, int& p);
bool isHere(vector<int> differences, int diff);
bool isCostas(vector<int>& shift);


/******************************************************
*
*       MISCELLANY
*
*******************************************************/

unsigned int readDim();
bool compare(vector<int>& v1, vector<int>& v2);
int powerMod(const int& x, const int& n, const int& mod);
int f(const int& x, const int& mod, int& a, int& b);
int LegendreComplexity(int& p);
int ConjectureComplexity(int& p, int& n1);
string satisfied_conjecture(int& d, int& conj);


/******************************************************
*
*       SEQUENCES
*
*******************************************************/


struct MSeq {
    vector<int> v;
    MSeq() {
        v = {0,0,1,0,1,1,1};
    }
    NTL::ZZ_p operator () (const int& i) const {
        return (NTL::ZZ_p) v[i];
    }
};

struct LegendreSeq {
    int mod;
    explicit LegendreSeq(int p) {
        mod = p;
    }
    NTL::ZZ_p operator () (const int& i) const {
        bool result;
        if ((i % mod) == 0)
            result = false;
        else
            result = isQuadratic(i,mod);
        return (NTL::ZZ_p) result;
    }
};

struct PolynomialSeq {
    int m;
    int alpha;
    int A;
    int B;
    explicit PolynomialSeq(int mod, int root, int a, int b) {
        m = mod; alpha = root; A = a; B = b;
    }
    int operator () (const int& i) {
        return f(powerMod(alpha, i, m), m, A, B);
    }
};


struct ExpQuadraticSeq {
    int mod;
    int root;
    int a;
    explicit ExpQuadraticSeq(int p, int r, int A) {
        mod = p; root = r; a = A;
    }
    int operator () (const int& i) {
        return (a * powerMod(root*root, i, mod) + powerMod(root, i, mod)) % mod;
    }
};

struct LogQuadraticSeq {
    int mod;
    int root;
    vector<int> logTable;
    explicit LogQuadraticSeq(int p, int r) {
        mod = p; root = r;
        logTable.resize(p);
        logTable[0] = -1;
        for (int i=0; i<mod-1; i++) {
            logTable[powerMod(root, i, mod)] = i;
        }
    }
    int operator () (const int& i) {
        int result = (4*powerMod(root*root, i, mod) + powerMod(root, i, mod)) % mod;
        return logTable[result];
    }
};


struct ExponentialSeq {
    int mod;
    int root;
    explicit ExponentialSeq(int m, int r){
        mod = m; root = r;
    }
    int operator () (const int& i) const {
        return powerMod(root, i, mod);
    }
};

struct OneSequence {
    int n;
    explicit OneSequence(int length) {
        n = length;
    }
    NTL::ZZ_p operator () (const int& i) const {
        if(i==0)    return (NTL::ZZ_p) 1;
        else        return (NTL::ZZ_p) 0;
    }
};

struct LosWelchSeq {
    int mod;
    int root;
    vector<int> log;
    explicit LosWelchSeq(int m, int r){
        mod = m;
        root = r;
        log.resize(m);
        for (int i=0; i<m-1; i++){
            log[powerMod(r,i,m)-1] = i;
        }
    }
    int operator () (const int& i) const {
        return log[i];
    }
};


template <typename Func>
void printSeq(Func func1, int period){
    for (int i = 0; i < period; i++) {
        cout << func1(i) << " ";
    }
    cout << endl;
}

int LogRootSeq(int i);
int noShiftFunc(const int& i);

#endif //RST_FUNCTIONS_H
