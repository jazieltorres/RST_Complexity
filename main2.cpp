#include <iostream>
#include "MultiDimArray.cpp"

// using namespace blitz;
using namespace std;

int main() {

    vector<long> v1({0,3,4}) ;
    vector<long> v2({1,2,0}) ;
    MExponent<3> m1(v1) ;
    MExponent<3> m2(v2) ;
    MExponent<3> m3 = m1+m2 ;
    if (m1.grlex_less(m2)) cout << m1 << endl ;
    else cout << m2 << endl ;
    // cout << "Compiles!" << endl;
    
    return 0;
}


// g++ -std=c++11 -pthread -march=native main.cpp  -lntl -lblitz -lgmp -lm
