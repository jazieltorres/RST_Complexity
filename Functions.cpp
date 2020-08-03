//
// Created by Jaziel Torres Fuentes on 8/3/20.
//

#include "Functions.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>


using namespace std;


/******************************************************
*
*       CHECKS
*
*******************************************************/

bool isQuadratic(const long& n, const long& p){
    for (long i=0; i<p; i++){
        if ((i*i %p) == n) return true;
    }
    return false;
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

/******************************************************
*
*       MISCELLANY
*
*******************************************************/

unsigned long readDim() {
    unsigned long dim;
    cin >> dim;
    return dim;
}

bool compare(vector<long>& v1, vector<long>& v2){
    return lexicographical_compare(v1.begin(), v1.end(), v2.begin(), v2.end());
}

long powerMod(const long& x, const long& n, const long& mod){
    if (n == 0) return 1;
    long result = x;
    for (long i=1; i<n; i++){
        result = (result * x) % mod;
    }
    return result;
}

// This is the polynomial (degree 3) that generates the PolynomialSeq
long f(const long& x, const long& mod, long& a, long& b) {
    return( a * powerMod(x, 3, mod) +
            b * powerMod(x, 2, mod) +
            x )
          % mod;
}

long LegendreComplexity(long& p) {
    long criteria = p % 8;
    if(criteria == 1) return (p-1)/2;
    if(criteria == 3) return p;
    if(criteria == 5) return p-1;
    if(criteria == 7) return (p+1)/2;
    return 0;
}


long LogRootSeq(long i){
    return i;
}
long noShiftFunc(const long& i){
    return 0;
}