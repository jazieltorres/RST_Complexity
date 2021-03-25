/******************************************************************************************
*
*               CLASS FOR MULTIVARIATE POLYNOMIALS
*
* This class is used to store the polynomials in the Grobner basis of the ideal of
* polynomials satisfying the recurrence relation in the multidimensional array.
*
* The class assumes the monomials are added in ascending order, with respect to some
* monomial ordering. That is, coefficient[0] contains the constant term.
******************************************************************************************/

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
    int terms();
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
int MultivarPolynomial<F,m>::terms() {
    return length;
}

template <typename F, int m>
void MultivarPolynomial<F,m>::print() {
    if(length == 0) {
        cout << "Empty polynomial." << endl;
        return;
    }

    if (m<4) {
//        Verifying if polynomial has constant term
        bool constant_term = true;
        for (int i=0; i<m; i++){
            if (exponents[0].exponent()[i] != 0)    constant_term = false;
        }
//        Printing constant term (if there is one)
        if (constant_term)  {
            cout << coefficients[0];
            if (length > 1) cout << " + ";
        }


        char variables[] = {'X', 'Y', 'Z'};
//        if no constant term, start from 0, otherwise, start from 1
        for (int i = constant_term; i < length; i++) {
            if (coefficients[i] != (F) 1)   cout << coefficients[i];
            for (int j = 0; j < m; j++) {
                int power = exponents[i].exponent()[j];
                if (power != 0) {
                    cout << variables[j];
                    if (power != 1)
                        cout << "^{" << power << "}";
                }
            }
            if (i != length - 1)   cout << " + ";
        }
        cout << endl;
    }
    else {
        for (int i = 0; i < length; i++) {
            cout << coefficients[i] << "x" << exponents[i] << '\t';
        }
        cout << endl;
    }
}