#include <iostream>
#include <vector>
#include "MExponent.h"

using namespace std;

int main() {
    vector<long> v1 = {0,3,4} ;
    vector<long> v2 = {1,2,0} ;
    MExponent m1(v1) ;
    MExponent m2(v2) ;
    MExponent m3 = m1+m2 ;
    if (m1.grlex_less(m2)) cout << m1 ;
    else cout << m2 ;
    // cout << m1+m2 << endl ; // no se puede
    return 0;
}
