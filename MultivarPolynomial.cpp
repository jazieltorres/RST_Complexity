#ifndef RST_MULTIVARPOLYNOMIAL_H
#define RST_MULTIVARPOLYNOMIAL_H
#include "Monomial.cpp"
#include <vector>

template <typename F, int m>
class MultivarPolynomial {
private:
    vector< Monomial<m> > exponents;
    vector<F> coefficients;
    int length;

public:
    MultivarPolynomial();
    void add_term(Monomial<m>, F coefficient = 1);
    void print();
};

#endif //RST_MULTIVARPOLYNOMIAL_H


template <typename F, int m>
MultivarPolynomial<F,m>::MultivarPolynomial() {
    length = 0;
}

template <typename F, int m>
void MultivarPolynomial<F,m>::add_term(Monomial<m> exponent, F coefficient) {
    coefficients.push_back(coefficient);
    exponents.push_back(exponent);
    length++;
}

template <typename F, int m>
void MultivarPolynomial<F,m>::print() {
    if(m == 1) {
        for(int i=0; i<length-1; i++) {
            if (coefficients[i] != (F)1)
                cout << coefficients[i];
            cout << "X^{" << exponents[i] << "} + ";
        }
        if (coefficients[length-1] != (F)1)
            cout << coefficients[length-1];
        cout << "X^{" << exponents[length-1] << "} + ";
    }
    else if (m == 2) {
        for (int i=0; i<length-1; i++) {
            if (coefficients[i] != (F)1)
                cout << coefficients[i];
            if (exponents[i].exponent()[0] != 0)
                cout << "X^{" << exponents[i].exponent()[0] << "}";
            if (exponents[i].exponent()[1] != 0)
                cout << "Y^{" << exponents[i].exponent()[1] << "}";
            cout << " + ";
        }
        if (coefficients[length-1] != (F)1)
            cout << coefficients[length-1];
        if (exponents[length-1].exponent()[0] != 0)
            cout << "X^{" << exponents[length-1].exponent()[0] << "}";
        if (exponents[length-1].exponent()[1] != 0)
            cout << "Y^{" << exponents[length-1].exponent()[1] << "}";
        cout << endl;
    }
    else if (m == 3) {
        for (int i=0; i<length-1; i++) {
            if (coefficients[i] != (F)1)
                cout << coefficients[i];
            if (exponents[i].exponent()[0] != 0)
                cout << "X^{" << exponents[i].exponent()[0] << "}";
            if (exponents[i].exponent()[1] != 0)
                cout << "Y^{" << exponents[i].exponent()[1] << "}";
            if (exponents[i].exponent()[2] != 0)
                cout << "Z^{" << exponents[i].exponent()[2] << "}";
            cout << " + ";
        }
        if (coefficients[length-1] != (F)1)
            cout << coefficients[length-1];
        if (exponents[length-1].exponent()[0] != 0)
            cout << "X^{" << exponents[length-1].exponent()[0] << "}";
        if (exponents[length-1].exponent()[1] != 0)
            cout << "Y^{" << exponents[length-1].exponent()[1] << "}";
        if (exponents[length-1].exponent()[2] != 0)
            cout << "Z^{" << exponents[length-1].exponent()[2] << "}";
        cout << endl;
    }
    else {
        for (int i = 0; i < length; i++) {
            cout << coefficients[i] << "x" << exponents[i] << '\t';
        }
        cout << endl;
    }
}