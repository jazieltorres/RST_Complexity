#ifndef MEXPONENT2_H
#define MEXPONENT2_H


#include <iostream>
#include <numeric> // accumulate
#include <blitz/array.h>


using namespace std;



template <long m>
class MExponent {
    private:
        blitz::TinyVector<long, m> exp ;
    public:
        // CONSTRUCTORS
        MExponent() {}
        MExponent(const vector<long>&) ;
        MExponent(const blitz::TinyVector<long, m>&) ;

        // DATA ACCESS AND MANIPULATION
        const blitz::TinyVector<long, m>& getExp() const ; // read-only access of exponents
        blitz::TinyVector<long, m>& editExp(); // access to modify exponents

        // ASSIGNMENT AND ARITHMETIC
        void operator=(const MExponent&) ; // Overload assignment
        MExponent operator+(const MExponent&) const ; // add two exponents
        MExponent mod(const MExponent&) const ; // component wise mod

        // MONOMIAL COMPARISONS
        bool lex_less(const MExponent& m2) const ; // compare by lex
        bool grlex_less(const MExponent& m2) const ; // compare by graded lex
        bool leq_d(const MExponent& m2) const ; // compare by divisibility

        friend ostream& operator<<(ostream& out, const MExponent<m>&) ; // Overload << to print vector of exponents

} ; // end MExponent
#endif // MEXPONENT_H


/******************************************************
*
*           CONSTRUCTORS
*
*******************************************************/

template <long m>
MExponent<m>::MExponent(const vector<long>& e) {
    for (int i = 0; i<exp.length(); i++) exp[i] = e[i] ;
}

template <long m>
MExponent<m>::MExponent(const blitz::TinyVector<long, m>& e) {
    exp = e ;
}

/******************************************************
*
*           DATA ACCESS AND MANIPULATION
*
*******************************************************/

// read-only access
template <long m>
const blitz::TinyVector<long, m>& MExponent<m>::getExp() const {
    return exp ;
}

// access to modify
template <long m>
blitz::TinyVector<long, m>& MExponent<m>::editExp() {
    return exp ;
}

/******************************************************
*
*           ASSIGNMENT AND ARITHMETIC
*
*******************************************************/

// assignment
template <long m>
void MExponent<m>::operator=(const MExponent& e) {
    exp = e.getExp() ;
}

// addition
template <long m>
MExponent<m> MExponent<m>::operator+(const MExponent& e) const {
    return MExponent(exp + e.getExp()) ;
}

// mod component wise
template <long m>
MExponent<m> MExponent<m>::mod(const MExponent<m>& e) const {
    blitz::TinyVector<long, m> result ;
    for(int i = 0; i < exp.length(); i++)
        result[i] = exp[i] % e.getExp()[i] ;
    return MExponent(result) ;
}

/******************************************************
*
*           MONOMIAL COMPARISONS
*
*******************************************************/

// lexicographic (lex)
template <long m>
bool MExponent<m>::lex_less(const MExponent<m>& e) const {
    for (int i = 0; i < exp.length(); i++) {
        long x = exp[i] - e.getExp()[i] ;
        if (x < 0)  return true ;
        else if (x > 0) return false ;
    }
    return false ;
}

// graded lexicographic (grlex)
template <long m>
bool MExponent<m>::grlex_less(const MExponent<m>& e) const {
    int sum1 = sum(exp) ;
    int sum2 = sum(e.getExp()) ;
    if (sum1 < sum2 || (sum1 == sum2 && lex_less(e))) return true ;
    return false ;
}

// divisibility
template <long m>
bool MExponent<m>::leq_d(const MExponent<m>& e) const {
    for(int i = 0; i < exp.length(); i++) {
        if(exp[i] > e.getExp()[i]) return false ;
    }
    return true ;
}


/******************************************************
*
*           PRINT (<< OVERLOAD)
*
*******************************************************/

template <long m>
ostream& operator<<(ostream& out, MExponent<m>& e) {
    out << "(" ;
    for(auto i = e.getExp().begin(); i != e.getExp().end()-1; i++)
        out << *i << ", " ;
    out << e.getExp()[e.getExp().length()-1] << ")";
    return out;
}
//
//
// /******************************************************
// *           END OF MExponent IMPLEMENTATION
// *******************************************************/
