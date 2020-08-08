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

bool isQuadratic(const int& n, const int& p){
    for (int i=0; i<p; i++){
        if ((i*i %p) == n) return true;
    }
    return false;
}

bool isPrime(int& n){
    if (n<3) return false;
    int bound = floor(sqrt(n));
    for (int i=2; i<=bound; i++){
        if (n%i == 0) return false;
    }
    return true;
}

bool isRoot(int& a, int& p) {
    int num = a*a;
    int ord = 1;
    while (num % p != a) {
        num = num*a % p;
        ord++;
    }
    if (ord == p-1) return true;
    return false;
}

bool isHere(vector<int> differences, int diff){
    for (auto d : differences){
        if (d==diff) return true;
    }
    return false;
}

bool isCostas(vector<int>& shift) {
    for (int i=1; i<shift.size()-1; i++){
        vector<int> differences;
        int first_index(0);
        while (first_index + i < shift.size()) {
            int diff = shift[first_index + i] - shift[first_index];
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

unsigned int readDim() {
    unsigned int dim;
    cin >> dim;
    return dim;
}

bool compare(vector<int>& v1, vector<int>& v2){
    return lexicographical_compare(v1.begin(), v1.end(), v2.begin(), v2.end());
}

int powerMod(const int& x, const int& n, const int& mod){
    if (n == 0) return 1;
    int result = x;
    for (int i=1; i<n; i++){
        result = (result * x) % mod;
    }
    return result;
}

// This is the polynomial (degree 3) that generates the PolynomialSeq
int f(const int& x, const int& mod, int& a, int& b) {
    return( a * powerMod(x, 3, mod) +
            b * powerMod(x, 2, mod) +
            x )
          % mod;
}

int LegendreComplexity(int& p) {
    int criteria = p % 8;
    if(criteria == 1) return (p-1)/2;
    if(criteria == 3) return p;
    if(criteria == 5) return p-1;
    if(criteria == 7) return (p+1)/2;
    return 0;
}


int LogRootSeq(int i){
    return i;
}
int noShiftFunc(const int& i){
    return 0;
}