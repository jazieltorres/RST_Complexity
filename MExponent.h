#include <iostream>
#include <vector>
#include <numeric>
using namespace std;

#ifndef MEXPONENT_H
#define MEXPONENT_H

// Class specification

class MExponent 
{
    public:
        void operator=(const MExponent&) ; // Overload assignment
        MExponent() ; // Default constructor
        MExponent operator+(const MExponent&) const ; // Overload + to add two exponents
        MExponent operator-(const MExponent&) const ; // Overload - to subtract two exponents
        void setExp(const vector<long>&) ; // Function to set the vector of exponents
        bool lex_less(const MExponent& m2) const ; // Function to compare by lex
        bool grlex_less(const MExponent& m2) const ; // Function to compare by graded lex
        bool leq_d(const MExponent& m2) const ; // Function to compare by division
        MExponent(const vector<long>&) ; // Constructor
        friend ostream& operator<<(ostream& out, const MExponent&) ; // Overload << to print vector of exponents
        const vector<long>& getExp() const ; // Function to retrieve vector of exponents

    private:
        vector<long> exp;

} ; // end MExponent

#endif // MEXPONENT_H

// Implementation

// Default constructor
MExponent::MExponent() : exp(3) { }

MExponent::MExponent(const vector<long>& e) {
    exp = e ;
 }

// Overload assignment of exponent
void MExponent::operator=(const MExponent& e) {
    exp = e.getExp() ;
}

// Overload + to add two exponents
MExponent MExponent::operator+(const MExponent& e) const {
    vector<long> result(exp.size()) ;
    for (int i = 0; i < exp.size(); i++) {
        result[i] = exp[i] + e.getExp()[i] ;
    }
    return MExponent(result) ;
}

 // Overload - to subtract two exponents
MExponent MExponent::operator-(const MExponent& e) const {
    vector<long> result(exp.size()) ;
    for (int i = 0; i < exp.size(); i++) {
        result[i] = exp[i] - e.getExp()[i] ;
    }
    return MExponent(result) ;
}

// Function to compare by lex
bool MExponent::lex_less(const MExponent& e) const {
    for (int i = 0; i < exp.size(); i++) {
        long x = exp[i] - e.getExp()[i] ;
        if (x < 0)  return true ;
        else if (x > 0) return false ;
    }
    return false ;
}

// Function to compare by graded lex
bool MExponent::grlex_less(const MExponent& e) const {
    int sum1 = accumulate(exp.begin(), exp.end(), 0) ;
    int sum2 = accumulate(e.getExp().begin(), e.getExp().end(), 0) ;
    if (sum1 < sum2 || (sum1 == sum2 && lex_less(e))) return true ;
    return false ;
}

// Function to compare less or equal by divisibility
bool MExponent::leq_d(const MExponent& e) const {
    for(int i = 0; i < exp.size(); i++) {
        if(exp[i] > e.getExp()[i]) return false ;
    }
    return true ;
}

// Overload << to print vector of exponents
ostream& operator<<(ostream& out, MExponent& e) {
    out << "(" ;
    for(int i = 0; i < e.getExp().size()-1; i++) out << e.getExp()[i] << ", " ;
    out << e.getExp()[e.getExp().size()-1] << ")";
    return out;
}

// Function to retrieve vector of exponents
const vector<long>& MExponent::getExp() const {
    return exp ;
}
