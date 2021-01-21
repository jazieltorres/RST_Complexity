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
    for(int i=0; i<length; i++) {
        cout << coefficients[i] << "x" << exponents[i] << '\t';
    }
    cout << endl;
}