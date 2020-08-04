#ifndef MEXPONENT_H
#define MEXPONENT_H


#include <vector>
#include <algorithm>
#include <iostream>
#include <blitz/array.h>


using namespace std;

template <int m>
class MExponent {
    private:
        blitz::TinyVector<int, m> exp ;
    public:
        // CONSTRUCTORS
        MExponent() = default;
        explicit MExponent(const vector<int>&) ;
        explicit MExponent(const blitz::TinyVector<int, m>&) ;

        // DATA ACCESS AND MANIPULATION
        const blitz::TinyVector<int, m>& getExp() const ; // read-only access of exponents
        blitz::TinyVector<int, m>& editExp(); // access to modify exponents

        // ASSIGNMENT AND ARITHMETIC
        void operator=(const MExponent&) ; // Overload assignment
        MExponent operator+(const MExponent&) const ; // add two exponents
        MExponent mod(const MExponent&) const ; // component wise mod

        // MONOMIAL COMPARISONS
        bool lex_less(const MExponent&) const ; // compare by lex Y<X
        bool lex_less2(const MExponent&) const ; // compare by lex X<Y
        bool grlex_less(const MExponent&) const ; // compare by graded lex
        bool leq_d(const MExponent&) const ; // compare by divisibility

        // Overload << to print vector of exponents
        friend ostream& operator<<(ostream& out, const MExponent<m>& e) {
          out << "(" ;
          for(auto i = e.getExp().begin(); i != e.getExp().end()-1; i++)
              out << *i << ", " ;
          out << e.getExp()[e.getExp().length()-1] << ")";
          return out;
        }
} ; // end MExponent
#endif // MEXPONENT_H


/******************************************************
*
*           CONSTRUCTORS
*
*******************************************************/

template <int m>
MExponent<m>::MExponent(const vector<int>& e) {
    for (int i = 0; i < exp.length(); i++) exp[i] = e[i] ;
}

template <int m>
MExponent<m>::MExponent(const blitz::TinyVector<int, m>& e) {
    exp = e ;
}


/******************************************************
*
*           DATA ACCESS AND MANIPULATION
*
*******************************************************/

// read-only access
template <int m>
const blitz::TinyVector<int, m>& MExponent<m>::getExp() const {
    return exp ;
}

// access to modify
template <int m>
blitz::TinyVector<int, m>& MExponent<m>::editExp() {
    return exp ;
}

/******************************************************
*
*           ASSIGNMENT AND ARITHMETIC
*
*******************************************************/

// assignment
template <int m>
void MExponent<m>::operator=(const MExponent& e) {
    exp = e.getExp() ;
}

// addition
template <int m>
MExponent<m> MExponent<m>::operator+(const MExponent& e) const {
    return MExponent(exp + e.getExp()) ;
}

// mod component wise
template <int m>
MExponent<m> MExponent<m>::mod(const MExponent<m>& e) const {
    blitz::TinyVector<int, m> result ;
    for(int i = 0; i < exp.length(); i++)
        result[i] = exp[i] % e.getExp()[i] ;
    return MExponent(result) ;
}

/******************************************************
*
*           MONOMIAL COMPARISONS
*
*******************************************************/

// lexicographic (lex) X<Y
template <int m>
bool MExponent<m>::lex_less(const MExponent<m>& e) const {
    return lexicographical_compare(exp.begin(), exp.end(), e.getExp().begin(), e.getExp().end());
}

//
template <int m>
bool MExponent<m>::lex_less2(const MExponent<m>& e) const {
    auto first1 = exp.end()-1;
    auto last1 = exp.begin()-1;
    auto first2 = e.getExp().end()-1;
    auto last2 = e.getExp().begin()-1;
    while (first1!=last1)
    {
        if (first2==last2 || *first2<*first1) return false;
        else if (*first1<*first2) return true;
        --first1; --first2;
    }
    return (first2!=last2);
}

// graded lexicographic (grlex)
template <int m>
bool MExponent<m>::grlex_less(const MExponent<m>& e) const {
    int sum1 = sum(exp) ;
    int sum2 = sum(e.getExp()) ;
    if (sum1 > sum2) return false;
    else if (sum1 < sum2) return true;
    else return lex_less(e);
}

// divisibility
template <int m>
bool MExponent<m>::leq_d(const MExponent<m>& e) const {
    for(int i = 0; i < exp.length(); i++) {
        if(exp[i] > e.getExp()[i]) return false;
    }
    return true ;
}


/******************************************************
*
*           PRINT (<< OVERLOAD)
*
*******************************************************/

// template <int m>
// ostream& operator<<(ostream& out, MExponent<m>& e) {
//     out << "(" ;
//     for(auto i = e.getExp().begin(); i != e.getExp().end()-1; i++)
//         out << *i << ", " ;
//     out << e.getExp()[e.getExp().length()-1] << ")";
//     return out;
// }

