#ifndef MONOMIAL_H
#define MONOMIAL_H


#include <vector>
#include <algorithm>
#include <iostream>
#include <blitz/array.h>


using namespace std;

template <int m>
class Monomial {
    private:
        blitz::TinyVector<int, m> exp ;
    public:
        // CONSTRUCTORS
        Monomial() = default;
        explicit Monomial(const vector<int>&) ;
        explicit Monomial(const blitz::TinyVector<int, m>&) ;

        // DATA ACCESS AND MANIPULATION
        const blitz::TinyVector<int, m>& exponent() const ; // read-only access of exponents
        blitz::TinyVector<int, m>& get_exponent(); // access to modify exponents

        // ASSIGNMENT AND ARITHMETIC
        void operator=(const Monomial&) ; // Overload assignment
        Monomial operator+(const Monomial&) const ; // add two exponents
        Monomial mod(const Monomial&) const ; // component wise mod
        bool equal(const Monomial&) const ;

        // MONOMIAL COMPARISONS
        bool lex_less(const Monomial&) const ; // compare by lex X1 > X2 > ... > Xn
        bool lex_less2(const Monomial&) const ; // compare by lex X1 < X2 < ... < Xn
        bool grlex_less(const Monomial&) const ; // compare by graded lex with lex_less
        bool grlex_less2(const Monomial&) const ; // compare by graded lex with lex_less2
        bool leq_d(const Monomial&) const ; // compare by divisibility

        // Overload << to print vector of exponents
        friend ostream& operator<<(ostream& out, const Monomial<m>& e) {
          out << "(" ;
          for(auto i = e.exponent().begin(); i != e.exponent().end()-1; i++)
              out << *i << ", " ;
          out << e.exponent()[e.exponent().length()-1] << ")";
          return out;
        }
} ;
#endif


/******************************************************
*
*           CONSTRUCTORS
*
*******************************************************/

template <int m>
Monomial<m>::Monomial(const vector<int>& e) {
    if (e.size() != m)
        cerr << "Vector of incorrect size. Has to be of size " << m << "." << endl;
    for (int i = 0; i < exp.length(); i++) exp[i] = e[i] ;
}

template <int m>
Monomial<m>::Monomial(const blitz::TinyVector<int, m>& e) {
    exp = e ;
}


/******************************************************
*
*           DATA ACCESS AND MANIPULATION
*
*******************************************************/

// read-only access
template <int m>
const blitz::TinyVector<int, m>& Monomial<m>::exponent() const {
    return exp ;
}

// access to modify
template <int m>
blitz::TinyVector<int, m>& Monomial<m>::get_exponent() {
    return exp ;
}

/******************************************************
*
*           ASSIGNMENT AND ARITHMETIC
*
*******************************************************/

// assignment
template <int m>
void Monomial<m>::operator=(const Monomial& e) {
    exp = e.exponent() ;
}

// addition
template <int m>
Monomial<m> Monomial<m>::operator+(const Monomial& e) const {
    return Monomial(exp + e.exponent()) ;
}

// mod component wise
template <int m>
Monomial<m> Monomial<m>::mod(const Monomial<m>& e) const {
    blitz::TinyVector<int, m> result ;
    for(int i = 0; i < exp.length(); i++)
        result[i] = exp[i] % e.exponent()[i] ;
    return Monomial(result) ;
}

// equal
template <int m>
bool Monomial<m>::equal(const Monomial<m>& e) const {
    for (int i=0; i<m; i++) {
        if (e.exponent()[i] != exp[i])
            return false;
    }
    return true;
}

/******************************************************
*
*           MONOMIAL COMPARISONS
*
*******************************************************/

// Lexicographic (lex) with X1 > X2 > ... > Xn
template <int m>
bool Monomial<m>::lex_less(const Monomial<m>& e) const {
    return lexicographical_compare(exp.begin(), exp.end(), e.exponent().begin(), e.exponent().end());
}

// Lexicographic (lex) with X1 < X2 < ... < Xn
template <int m>
bool Monomial<m>::lex_less2(const Monomial<m>& e) const {
    auto first1 = exp.end()-1;
    auto last1 = exp.begin()-1;
    auto first2 = e.exponent().end()-1;
    auto last2 = e.exponent().begin()-1;
    while (first1!=last1)
    {
        if (first2==last2 || *first2<*first1) return false;
        else if (*first1<*first2) return true;
        --first1; --first2;
    }
    return (first2!=last2);
}

// Graded lexicographic (grlex) with X1 > X2 > ... > Xn
template <int m>
bool Monomial<m>::grlex_less(const Monomial<m>& e) const {
    int sum1 = sum(exp) ;
    int sum2 = sum(e.exponent()) ;
    if (sum1 > sum2) return false;
    else if (sum1 < sum2) return true;
    else return lex_less(e);
}

// Graded lexicographic (grlex) with X1 < X2 < ... < Xn
template <int m>
bool Monomial<m>::grlex_less2(const Monomial<m>& e) const {
    int sum1 = sum(exp) ;
    int sum2 = sum(e.exponent()) ;
    if (sum1 > sum2) return false;
    else if (sum1 < sum2) return true;
    else return lex_less2(e);
}

// Divisibility partial ordering
template <int m>
bool Monomial<m>::leq_d(const Monomial<m>& e) const {
    for(int i = 0; i < exp.length(); i++) {
        if(exp[i] > e.exponent()[i]) return false;
    }
    return true ;
}


/******************************************************
*
*           PRINT (<< OVERLOAD)
*
*******************************************************/

// template <int m>
// ostream& operator<<(ostream& out, Monomial<m>& e) {
//     out << "(" ;
//     for(auto i = e.exponent().begin(); i != e.exponent().end()-1; i++)
//         out << *i << ", " ;
//     out << e.exponent()[e.exponent().length()-1] << ")";
//     return out;
// }

