#include <iostream>
#include "MultiDimArray.hpp"

using namespace std;


int main() {
    blitz::TinyVector<long, 5> A(3);
    blitz::TinyVector<long, 5> B(4);
    cout << A << endl;
    return 0;
}


// g++ -std=c++11 -pthread -march=native main.cpp  -lntl -lblitz -lgmp -lm
