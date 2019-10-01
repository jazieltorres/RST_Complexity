#ifndef MEXPONENT2_H
#define MEXPONENT2_H


#include <vector>
#include <iostream>
#include <numeric> // accumulate
#include <blitz/array.h>


using namespace std;

template <unsigned long m>
class MExponent {
    private:
        blitz::TinyVector<unsigned long, m> exp ;
    public:
        // CONSTRUCTORS
        MExponent() {}
        MExponent(const vector<unsigned long>&) ;
        MExponent(const blitz::TinyVector<unsigned long, m>&) ;

        // DATA ACCESS AND MANIPULATION
        const blitz::TinyVector<unsigned long, m>& getExp() const ; // read-only access of exponents
        blitz::TinyVector<unsigned long, m>& editExp(); // access to modify exponents

        // ASSIGNMENT AND ARITHMETIC
        void operator=(const MExponent&) ; // Overload assignment
        MExponent operator+(const MExponent&) const ; // add two exponents
        MExponent mod(const MExponent&) const ; // component wise mod

        // MONOMIAL COMPARISONS
        bool lex_less(const MExponent& m2) const ; // compare by lex
        bool grlex_less(const MExponent& m2) const ; // compare by graded lex
        bool leq_d(const MExponent& m2) const ; // compare by divisibility

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

template <unsigned long m>
MExponent<m>::MExponent(const vector<unsigned long>& e) {
    for (unsigned long i = 0; i < exp.length(); i++) exp[i] = e[i] ;
}

template <unsigned long m>
MExponent<m>::MExponent(const blitz::TinyVector<unsigned long, m>& e) {
    exp = e ;
}

/******************************************************
*
*           DATA ACCESS AND MANIPULATION
*
*******************************************************/

// read-only access
template <unsigned long m>
const blitz::TinyVector<unsigned long, m>& MExponent<m>::getExp() const {
    return exp ;
}

// access to modify
template <unsigned long m>
blitz::TinyVector<unsigned long, m>& MExponent<m>::editExp() {
    return exp ;
}

/******************************************************
*
*           ASSIGNMENT AND ARITHMETIC
*
*******************************************************/

// assignment
template <unsigned long m>
void MExponent<m>::operator=(const MExponent& e) {
    exp = e.getExp() ;
}

// addition
template <unsigned long m>
MExponent<m> MExponent<m>::operator+(const MExponent& e) const {
    return MExponent(exp + e.getExp()) ;
}

// mod component wise
template <unsigned long m>
MExponent<m> MExponent<m>::mod(const MExponent<m>& e) const {
    blitz::TinyVector<unsigned long, m> result ;
    for(unsigned long i = 0; i < exp.length(); i++)
        result[i] = exp[i] % e.getExp()[i] ;
    return MExponent(result) ;
}

/******************************************************
*
*           MONOMIAL COMPARISONS
*
*******************************************************/

// lexicographic (lex)
template <unsigned long m>
bool MExponent<m>::lex_less(const MExponent<m>& e) const {
    for (unsigned i = 0; i < exp.length(); i++) {
        unsigned long x = exp[i] - e.getExp()[i] ;
        if (x < 0)  return true ;
        else if (x > 0) return false ;
    }
    return false ;
}

// graded lexicographic (grlex)
template <unsigned long m>
bool MExponent<m>::grlex_less(const MExponent<m>& e) const {
    int sum1 = sum(exp) ;
    int sum2 = sum(e.getExp()) ;
    if (sum1 < sum2 || (sum1 == sum2 && lex_less(e))) return true ;
    return false ;
}

// divisibility
template <unsigned long m>
bool MExponent<m>::leq_d(const MExponent<m>& e) const {
    for(unsigned long i = 0; i < exp.length(); i++) {
        if(exp[i] > e.getExp()[i]) return false;
    }
    return true ;
}


/******************************************************
*
*           PRINT (<< OVERLOAD)
*
*******************************************************/

// template <long m>
// ostream& operator<<(ostream& out, MExponent<m>& e) {
//     out << "(" ;
//     for(auto i = e.getExp().begin(); i != e.getExp().end()-1; i++)
//         out << *i << ", " ;
//     out << e.getExp()[e.getExp().length()-1] << ")";
//     return out;
// }
//
//
// /******************************************************
// *           END OF MExponent IMPLEMENTATION
// *******************************************************/
