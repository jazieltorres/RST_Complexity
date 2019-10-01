#ifndef MEXPONENT_H
#define MEXPONENT_H


#include <iostream>
#include <vector>
#include <numeric> // accumulate
using namespace std;

class MExponent
{
    public:
        // CONSTRUCTORS
        MExponent() ; // Default constructor
        MExponent(unsigned long& m);
        MExponent(const vector<unsigned long>&) ;

        // DATA ACCESS AND MANIPULATION
        const vector<long>& getExp() const ; // read-only access of exponents
        vector<long>& editExp(); // access to modify exponents
        void setExp(const vector<long>&) ; // set exponents

        // ASSIGNMENT AND ARITHMETIC
        void operator=(const MExponent&) ; // Overload assignment
        MExponent operator+(const MExponent&) const ; // add two exponents
        MExponent operator-(const MExponent&) const ; // subtract two exponents
        MExponent scalar_mult(const long&) const ; // scalar multiplication
        MExponent mod(const MExponent&) const ; // component wise mod

        // MONOMIAL COMPARISONS
        bool lex_less(const MExponent& m2) const ; // compare by lex
        bool grlex_less(const MExponent& m2) const ; // compare by graded lex
        bool leq_d(const MExponent& m2) const ; // compare by divisibility

        friend ostream& operator<<(ostream& out, const MExponent&) ; // Overload << to print vector of exponents


    private:
        vector<long> exp;

} ; // end MExponent
#endif // MEXPONENT_H



/******************************************************
*
*           CONSTRUCTORS
*
*******************************************************/

// Default
MExponent::MExponent() : exp(3) { }

MExponent::MExponent(unsigned long& m) {
    exp.resize(m);
}

MExponent::MExponent(const vector<unsigned long>& e) {
    exp = e ;
}

/******************************************************
*
*           DATA ACCESS AND MANIPULATION
*
*******************************************************/

// read-only access
const vector<long>& MExponent::getExp() const {
    return exp ;
}

// access to modify
vector<long>& MExponent::editExp() {
    return exp ;
}

// set exponents
void MExponent::setExp(const vector<long>& e) {
    exp = e ;
}

/******************************************************
*
*           ASSIGNMENT AND ARITHMETIC
*
*******************************************************/

// assignment
void MExponent::operator=(const MExponent& e) {
    exp = e.getExp() ;
}

// addition
MExponent MExponent::operator+(const MExponent& e) const {
    vector<long> result(exp.size()) ;
    for (int i = 0; i < exp.size(); i++) {
        result[i] = exp[i] + e.getExp()[i] ;
    }
    return MExponent(result) ;
}

// substraction
MExponent MExponent::operator-(const MExponent& e) const {
    vector<long> result(exp.size()) ;
    for (int i = 0; i < exp.size(); i++)
        result[i] = exp[i] - e.getExp()[i] ;
    return MExponent(result) ;
}

// scalar multiplication
MExponent MExponent::scalar_mult(const long& alpha) const {
    vector<long> result(exp.size()) ;
    for(int i = 0; i < exp.size(); i++)
        result[i] = exp[i] * alpha ;
    return MExponent(result) ;
}

// mod component wise
MExponent MExponent::mod(const MExponent& e) const {
    vector<long> result(exp.size()) ;
    for(int i = 0; i < exp.size(); i++)
        result[i] = exp[i] % e.getExp()[i] ;
    return MExponent(result) ;
}

/******************************************************
*
*           MONOMIAL COMPARISONS
*
*******************************************************/

// lexicographic (lex)
bool MExponent::lex_less(const MExponent& e) const {
    for (int i = 0; i < exp.size(); i++) {
        long x = exp[i] - e.getExp()[i] ;
        if (x < 0)  return true ;
        else if (x > 0) return false ;
    }
    return false ;
}

// graded lexicographic (grlex)
bool MExponent::grlex_less(const MExponent& e) const {
    int sum1 = accumulate(exp.begin(), exp.end(), 0) ;
    int sum2 = accumulate(e.getExp().begin(), e.getExp().end(), 0) ;
    if (sum1 < sum2 || (sum1 == sum2 && lex_less(e))) return true ;
    return false ;
}

// divisibility
bool MExponent::leq_d(const MExponent& e) const {
    for(int i = 0; i < exp.size(); i++) {
        if(exp[i] > e.getExp()[i]) return false ;
    }
    return true ;
}


/******************************************************
*
*           PRINT (<< OVERLOAD)
*
*******************************************************/

ostream& operator<<(ostream& out, MExponent& e) {
    out << "(" ;
    for(int i = 0; i < e.getExp().size()-1; i++)
        out << e.getExp()[i] << ", " ;
    out << e.getExp()[e.getExp().size()-1] << ")";
    return out;
}


/******************************************************
*           END OF MExponent IMPLEMENTATION
*******************************************************/
