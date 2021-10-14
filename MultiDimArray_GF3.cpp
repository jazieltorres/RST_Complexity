#ifndef MULTIDIMARRAY_GF2_H
#define MULTIDIMARRAY_GF2_H

#include <iostream>
#include <functional>
#include <vector>
#include <forward_list>
#include <unordered_map>
#include "MultivarPolynomial.cpp"
#include <algorithm> //sort
#include "NTL/GF2.h"
#include "boost/dynamic_bitset.hpp"
#include <memory>
using namespace std;

typedef NTL::GF2 GF2;

template <int m>
class MultiDimArray_GF2 {
private:
    blitz::Array<GF2, m> A;
    blitz::TinyVector<int, m> period;
    vector< Monomial<m> > lead_monomials;
    vector< MultivarPolynomial<GF2,m> > groebner_basis;
    unsigned int size;
    int delta_size;
    int ordering_number;
public:
//  Constructor: receives a blitz array of dimension m, with entries in GF2.
    explicit MultiDimArray_GF2(blitz::Array<GF2,m>&);

//  Constructor: receives period vector and resizes array
    explicit MultiDimArray_GF2(const blitz::TinyVector<int,m>&);

//    Constructor: receives a sequence as vector (for one dimension only)
    explicit MultiDimArray_GF2(const vector<GF2>&);

//  Constructor: receives shift sequence, column sequence, and their periods, respectively.
    MultiDimArray_GF2(const function<int (int)>&, const function<GF2 (int)>&, int, int);

//  Constructor: same as above but shift sequence is a vector.
    MultiDimArray_GF2(const vector<int>&, GF2, const function<GF2 (int)>&, int);

//  Constructor: same as above but both sequences are vectors.
    MultiDimArray_GF2(const vector<int>&, const vector<GF2>&, const GF2);

//  Getters
    int dimension();
    int complexity();
    double normalized_complexity();
    int period_size();
    blitz::TinyVector<int, m> period_vector();
    string ordering_used();


    void set_at(const blitz::TinyVector<int,m>&, GF2&);
    void print_array();
    void print_basis();
    void draw_lead_monomials();

//  Rubio-Sweedler-Taylor algorithm for computing the linear complexity
    void RST(int ordering = 1);
    void RST_NEW(int ordering = 1);
    void RST_simple();
    void RST_simpleNEW();
    void RST_simpleNEW2();
} ;

//void insertBefore(vector< vector<GF2> >&, int&, const int&);

#endif




/******************************************************
*
*       CONSTRUCTORS
*
*******************************************************/


//  Constructor: receives a blitz array of dimension m, with entries in GF2.
template <int m>
MultiDimArray_GF2<m>::MultiDimArray_GF2(blitz::Array<GF2,m>& array){
    A.reference(array);
    size = 1;
    for (int i=0; i<m; i++){
        period[i] = A.extent(i);
        size = size * period[i];
    }
    delta_size = -1;
    ordering_number = 0;
}

//  Constructor: receives period vector and resizes array
template<int m>
MultiDimArray_GF2<m>::MultiDimArray_GF2(const blitz::TinyVector<int,m>& period_vector){
    period = period_vector;
    A.resize(period);
    size = 1;
    for (int i=0; i<m; i++){
        size = size * period[i];
    }
    delta_size = -1;
    ordering_number = 0;
}


template <int m>
MultiDimArray_GF2<m>::MultiDimArray_GF2(const vector<GF2>& seq){
    if (m != 1) {
        cout << "ERROR. Allowed dimension for this constructor: 1" << endl;
    }
    else {
        period = seq.size();
        size = seq.size();
        A.resize(period);
        delta_size = -1;
        ordering_number = 0;
        blitz::TinyVector<int, 1> index;
        for (int i = 0; i < size; i++) {
            index = i;
            A(index) = seq[i];
        }
    }
}


// MultiDimArray_GF2(Shift seq, column seq, horizontal period, column period)
template <int m>
MultiDimArray_GF2<m>::MultiDimArray_GF2(const function<int(int)>& func1, const function<GF2(int)>& func2,
                                        int n1, int n2) {
    if (m != 2) {
        cout << "ERROR. Allowed dimension for this constructor: 2" << endl;
    }
    else {
        period = n1, n2;
        A.resize(period);
        size = n1*n2;
        delta_size = -1;
        ordering_number = 0;
        blitz::TinyVector<int, 2> index;
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                index = i, j;
                if (func1(i) == -1) { // When the shift sequence has a star (-1)
                    if (n2 % 4 == 1)
                        A(index) = 1;
                    else
                        A(index) = 0;
                } else
                    A(index) = func2(((j - func1(i)) % n2 + n2) % n2);
            }
        }
    }
}


// Constructor for shift_seq as vector
template <int m>
MultiDimArray_GF2<m>::MultiDimArray_GF2(const vector<int>& shift_seq, GF2 constant,
                                        const function<GF2(int)>& column_seq, int vertical_period) {
    if (m == 2) {
        int horizontal_period = shift_seq.size();
        period = horizontal_period, vertical_period;
        A.resize(period);
        size = shift_seq.size() * vertical_period;
        delta_size = -1;
        ordering_number = 0;
        blitz::TinyVector<int, 2> index;
        for (int i = 0; i < horizontal_period; i++) {
            for (int j = 0; j < vertical_period; j++) {
                index = i, j;
                if (shift_seq[i] != -1) //shift not infinity
                    A(index) = column_seq(((j - shift_seq[i]) % vertical_period + vertical_period) % vertical_period);
                else
                    A(index) = constant;
            }
        }
    }
    else {
        cout << "ERROR: Allowed dimension for this constructor: 2" << endl;
    }
}


// Constructor for both sequences as vector
template <int m>
MultiDimArray_GF2<m>::MultiDimArray_GF2(const vector<int>& shift_seq, const vector<GF2>& column_seq, const GF2 constant) {
    if (m == 2) {
        int horizontal_period = shift_seq.size();
        int vertical_period = column_seq.size();
        period = horizontal_period, vertical_period;
        A.resize(period);
        size = horizontal_period * vertical_period;
        delta_size = -1;
        ordering_number = 0;
        blitz::TinyVector<int, 2> index;
        for (int i = 0; i < horizontal_period; i++) {
            for (int j = 0; j < vertical_period; j++) {
                index = i, j;
                if (shift_seq[i] != -1) //shift not infinity
                    A(index) = column_seq[((j - shift_seq[i]) % vertical_period + vertical_period) % vertical_period];
                else
                    A(index) = constant;
            }
        }
    }
    else {
        cout << "ERROR: Allowed dimension for this constructor: 2" << endl;
    }
}



/******************************************************
*
*       GETTERS
*
*******************************************************/
template <int m>
int MultiDimArray_GF2<m>::dimension() {
    return m;
}

template <int m>
int MultiDimArray_GF2<m>::complexity() {
    if (delta_size == -1) {
        this->RST();
        return delta_size;
    }
    else
        return delta_size;
}

template <int m>
double MultiDimArray_GF2<m>::normalized_complexity() {
    if (delta_size == -1) {
        this->RST();
        return static_cast<double>(delta_size)/size;
    }
    else
        return static_cast<double>(delta_size)/size;
}

template <int m>
int MultiDimArray_GF2<m>::period_size() {
    return size;
}

template <int m>
blitz::TinyVector<int, m> MultiDimArray_GF2<m>::period_vector() {
    return period;
}

template <int m>
string MultiDimArray_GF2<m>::ordering_used() {
    if (m == 1)
        return "The only monomial ordering in GF2[X]";
    switch (ordering_number) {
        case 1 :
            if (m == 2)       return "Graded lexicographic X > Y.";
            else if (m == 3)  return "Graded lexicographic X > Y > Z.";
            else              return "Graded lexicographic X1 > X2 > ... > Xn.";
        case 2 :
            if (m == 2)       return "Graded lexicographic X < Y.";
            else if (m == 3)  return "Graded lexicographic X < Y < Z.";
            else              return "Graded lexicographic X1 < X2 < ... < Xn.";
        case 3 :
            if (m == 2)       return "Lexicographic X > Y.";
            else if (m == 3)  return "Lexicographic X > Y > Z.";
            else              return "Lexicographic X1 > X2 > ... > Xn.";
        case 4:
            if (m == 2)       return "Lexicographic X < Y.";
            else if (m == 3)  return "Lexicographic X < Y < Z.";
            else              return "Lexicographic X1 < X2 < ... < Xn.";
        case 5 :
            if (m == 2)       return "Graded lexicographic X > Y.";
            else if (m == 3)  return "Graded lexicographic X > Y > Z.";
            else              return "Graded lexicographic X1 > X2 > ... > Xn.";
        case 6 :
            if (m == 2)       return "Graded lexicographic X < Y.";
            else if (m == 3)  return "Graded lexicographic X < Y < Z.";
            else              return "Graded lexicographic X1 < X2 < ... < Xn.";
        default:
            return "No ordering used yet.";
    }
}



/******************************************************
*
*       OTHER CLASS METHODS
*
*******************************************************/

template <int m>
void MultiDimArray_GF2<m>::set_at(const blitz::TinyVector<int, m>& coordinates, GF2& value) {
    A(coordinates) = value;
}

template <int m>
void MultiDimArray_GF2<m>::print_array() {
    blitz::TinyVector<int, m> index;
    for (int j=period[1]-1; j>=0; j--) {
        for (int i=0; i<period[0]; i++) {
            index = i,j;
            cout << A(index) << " ";
        }
        cout << endl;
    }
}

template <int m>
void MultiDimArray_GF2<m>::print_basis() {
    for (auto poly : groebner_basis)
        poly.print();
}

template <int m>
void MultiDimArray_GF2<m>::draw_lead_monomials() {
    for(int i=2*period[1]+1; i>0; i--){
        for(int j=0; j<period[0]; j++) {
            if (1 == i % 2) {
                cout << "--";
            }
            else {
                bool yes = false;
                for (auto e : lead_monomials){
                    if ((e.exponent()[0]== j) && (e.exponent()[1]==i/2-1)) yes = true;
                }
                if(yes) cout << "|*";
                else cout << "| ";
            }
        }
        if(1 == i%2) cout << " " << endl;
        else cout << "|" << endl;
    }
}


/******************************************************
*
*       MODULES FOR RST IMPLEMENTATION
*
*******************************************************/

template<int m>
void generateDivisors(int depth, int& index, Monomial<m>& e,
                      const Monomial<m>& bound, vector< Monomial<m> >& list){
    if (depth > 0){
        for(int i = 0; i <= bound.exponent()[depth-1]; i++) {
            e.get_exponent()[depth-1] = i;
            generateDivisors(depth-1, index, e, bound, list);
        }
    }
    else {
        list[index] = e;
        index++;
    }
}

template<int m>
unsigned long sizeOfDivisors(const Monomial<m> e) {
    unsigned long size(1);
    for (auto x : e.exponent()) {
        size = size * (x + 1);
    }
    return size;
}

// Row reduction with the identity matrix to obtain polynomial in the basis
static void reduce(vector< boost::dynamic_bitset<> >& M, vector< boost::dynamic_bitset<> >& Id) {
    if (M.empty())
        return;

    const unsigned int m = M.size() - 1;
    const unsigned int width = M[0].size();
    unsigned  int j = 0;

    boost::dynamic_bitset<> *mm = &(M[m]);
    boost::dynamic_bitset<> *idm = &(Id[m]);

    boost::dynamic_bitset<> *mi, *idi;

    mi = &(M[0]);
    idi = &(Id[0]);
    for (int i = 0; i < m; i++) {
        while ((*mi)[j] == 0 && j < width) {
            if ((*mm)[j] != 0) { // insert last row before row i
                for (int n = 0; n < (m - i); n++) {
                    M[m-n].swap(M[m-n-1]) ;
                }
                for (int n = 0; n < (m - i); n++) {
                    Id[m-n].swap(Id[m-n-1]) ;
                }
                return;
            }
            j++;
        }
        if (j == width) // implies matrix has row of zeros: ERROR
            return;

        if ((*mm)[j]) { //add row i to eliminate the one in M[m][j]
            (*mm) ^= (*mi);
            (*idm) ^= (*idi);
        }

        j++;
        mi++;
        idi++;
    }
}

// Row reduction just for computing the complexity
static void reduce(vector< boost::dynamic_bitset<> >& M) {
    if (M.empty())
        return;

    const unsigned int m = M.size() - 1;
    const unsigned int width = M[0].size();
    unsigned  int j = 0;

    boost::dynamic_bitset<> *mm = &(M[m]);
    boost::dynamic_bitset<> *mi;

    mi = &(M[0]);
    for (int i = 0; i < m; i++) {
        while ((*mi)[j] == 0 && j < width) {
            if ((*mm)[j] != 0) { // insert last row before row i
                for (int n = 0; n < (m - i); n++)
                    M[m-n].swap(M[m-n-1]) ;
                return;
            }
            j++;
        }
        if (j == width) // implies matrix has row of zeros: ERROR
            return;

        if ((*mm)[j]) { //add row i to eliminate the one in M[m][j]
            (*mm) ^= (*mi);
        }
        j++;
        mi++;
    }
}

static void printMatrix (const vector< boost::dynamic_bitset<> >& M) {
    if (M.empty()) return;
    int m = M.size();
    int n = M[0].size();

    if (n>15) n = 15;

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            cout << M[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

template<int m>
boost::dynamic_bitset<> rowPiAlpha(const Monomial<m>& alpha, const blitz::Array<GF2,m>& A,
                                   vector< Monomial<m> >& exponentsColumn, Monomial<m> period) {
    unsigned long n = exponentsColumn.size();
    GF2 one(1);
    boost::dynamic_bitset<> row(n);
    for (unsigned int i = 0; i < n; i++) {
        Monomial<m> index = exponentsColumn[i]+alpha;
        if (A(index.mod(period).exponent()) == one)
            row[i] = true;
    }
//    cout << "Row PI(" << alpha << "):" << endl;
//    for(int i=0; i<n; i++) {
//        cout << row[i] << " ";
//    } cout << endl << endl;
    return row;
}

template<int m>
boost::dynamic_bitset<>* rowPiAlphaNEW(const Monomial<m>& alpha, const blitz::Array<GF2,m>& A,
                                       vector< Monomial<m> >& exponentsColumn, Monomial<m> period) {
    unsigned long n = exponentsColumn.size();
    GF2 one(1);
    auto *ptr = new  boost::dynamic_bitset<>(n);
    for (unsigned int i = 0; i < n; i++) {
        Monomial<m> index = exponentsColumn[i]+alpha;
        if (A(index.mod(period).exponent()) == one)
            (*ptr)[i] = true;
    }
//    cout << "Row PI(" << alpha << "):" << endl;
//    for(int i=0; i<n; i++) {
//        cout << row[i] << " ";
//    } cout << endl << endl;
    return ptr;
}

// Functor used in exponentsRow.remove_if() for removing all the multiples of alpha
template<int m>
struct isMultiple {
    Monomial<m> alpha;
    explicit isMultiple(Monomial<m>& e) {alpha = e;}
    bool operator () (const Monomial<m>& e) const {
        return alpha.leq_d(e);
    }
};

template<int m>
MultivarPolynomial<GF2,m> get_polynomial(vector< boost::dynamic_bitset<> >& id_matrix, vector< Monomial<m> >& id_column) {
    MultivarPolynomial<GF2,m> poly;
    GF2 one(1);
    unsigned int n = id_matrix.size()-1;
    for (int i = 0; i < id_matrix[0].size(); i++) {
        if(id_matrix[n][i]) poly.add_term(id_column[i], one);
    }
    return poly;
}



/******************************************************
*
*               MONOMIAL ORDERINGS
*
*******************************************************/
// Graded lexicographic X1 > X2 > ... > Xn
template<int m>
bool ordering_1(const Monomial<m>& e1, const Monomial<m>& e2) {
    return (e1.grlex_less(e2));
}

// Graded lexicographic X1 < X2 < ... < Xn
template<int m>
bool ordering_2(const Monomial<m>& e1, const Monomial<m>& e2) {
    return (e1.grlex_less(e2));
}

// Lexicographic with X1 > X2 > ... > Xn
template<int m>
bool ordering_3(const Monomial<m>& e1, const Monomial<m>& e2) {
    return (e1.lex_less(e2));
}

// Lexicographic with X1 < X2 < ... < Xn
template<int m>
bool ordering_4(const Monomial<m>& e1, const Monomial<m>& e2) {
    return (e1.lex_less2(e2));
}


/****************************************************************************
*
*                       THE ONE AND ONLY: RST
*
* Calls RST with respect to a monomial ordering:
*   1. Graded lexicographic X1 > X2 > ... > Xn
*   2. Graded lexicographic X1 < X2 < ... < Xn
*   3. Lexicographic X1 > X2 > ... > Xn
*   4. Lexicographic X1 < X2 < ... < Xn
****************************************************************************/
template<int m>
void MultiDimArray_GF2<m>::RST(int ordering) {
    delta_size = 0;
    lead_monomials.empty();
    groebner_basis.empty();

    Monomial<m>  alpha; //dummy container
    Monomial<m>  e_period(period);
    Monomial<m>  bound;
    for (int i = 0; i < m; i++) {
        bound.get_exponent()[i] = 2 * period[i] - 1;
    }

    const unsigned long width = sizeOfDivisors(bound);
    const unsigned long id_width = sizeOfDivisors(e_period);

    vector< Monomial<m> >  exponentsColumn(width);
    vector< Monomial<m> >  id_column;
    forward_list< Monomial<m> >  exponentsRow;

    int index = 0;
    generateDivisors<m>(m, index, alpha, bound, exponentsColumn);

    switch (ordering) {
        case 1 :
            sort(exponentsColumn.begin(), exponentsColumn.end(), ordering_1<m>);
            break;
        case 2 :
            sort(exponentsColumn.begin(), exponentsColumn.end(), ordering_2<m>);
            break;
        case 3 :
            sort(exponentsColumn.begin(), exponentsColumn.end(), ordering_3<m>);
            break;
        case 4 :
            sort(exponentsColumn.begin(), exponentsColumn.end(), ordering_4<m>);
            break;
        default:
            cerr << "Invalid ordering argument in RST." << endl;
            delta_size = -1;
            return;
    }

    for (auto e = exponentsColumn.rbegin(); e != exponentsColumn.rend(); e++) {
        if (e->leq_d(e_period))
            exponentsRow.push_front(*e);
    }

//    cout << "exponent row" << endl;
//    for (auto e = exponentsRow.begin(); e != exponentsRow.end(); e++) {
//        cout << *e << " ";
//    } cout << endl;

    vector< boost::dynamic_bitset<> > matrix;
    vector< boost::dynamic_bitset<> > id_matrix;
    unsigned long ctr = 0;
    unsigned long height;
    while (!exponentsRow.empty()) {
        alpha = exponentsRow.front();
        matrix.push_back(rowPiAlpha<m>(alpha, A, exponentsColumn, e_period));

//        cout << "Row PI" << alpha << endl;
//        cout << rowPiAlpha<m>(alpha, A, exponentsColumn, e_period) << endl << endl;

//        Adding a row to the identity matrix
        boost::dynamic_bitset<> id_row(id_width);
        id_row[ctr] = true;
        id_matrix.push_back(id_row);

        id_column.push_back(alpha);

//        cout << "Before reduce" << endl;
//        printMatrix(matrix);

        reduce(matrix, id_matrix);


//        cout << "After reduce" << endl;
//        printMatrix(matrix);
//        cout << "\n\n" << endl;


        if (matrix[matrix.size()-1].none()) {
            lead_monomials.push_back(alpha);
            groebner_basis.push_back(get_polynomial(id_matrix, id_column));
            matrix.pop_back();
            id_matrix.pop_back();
            exponentsRow.remove_if(isMultiple<m>(alpha));
        }
        else {
            exponentsRow.pop_front();
            delta_size++;
        }
        ctr++;
    }
}


/******************************************************
*
*    RST THAT ONLY COMPUTES COMPLEXITY, NOT BASIS
*
*******************************************************/
template<int m>
void MultiDimArray_GF2<m>::RST_simple() {
    delta_size = 0;
    lead_monomials.empty();
    groebner_basis.empty();

    Monomial<m> alpha; //dummy container
    Monomial<m> e_period(period);
    Monomial<m> bound;
    for (int i = 0; i < m; i++) {
        bound.get_exponent()[i] = 2 * period[i] - 1;
    }

    const unsigned long width = sizeOfDivisors(bound);

    vector<Monomial<m> > exponentsColumn(width);
    vector<Monomial<m> > id_column;
    forward_list<Monomial<m> > exponentsRow;

    int index = 0;
    generateDivisors<m>(m, index, alpha, bound, exponentsColumn);
    sort(exponentsColumn.begin(), exponentsColumn.end(), ordering_1<m>);
    for (auto e = exponentsColumn.rbegin(); e != exponentsColumn.rend(); e++) {
        if (e->leq_d(e_period))
            exponentsRow.push_front(*e);
    }

//    cout << "exponent row" << endl;
//    for (auto e = exponentsRow.begin(); e != exponentsRow.end(); e++) {
//        cout << *e << " ";
//    } cout << endl;

    vector< boost::dynamic_bitset<> > matrix;

    while (!exponentsRow.empty()) {
        alpha = exponentsRow.front();
        matrix.push_back(rowPiAlpha<m>(alpha, A, exponentsColumn, e_period));

//        cout << "Row PI" << alpha << endl;
//        cout << rowPiAlpha<m>(alpha, A, exponentsColumn, e_period) << endl << endl;

//        cout << "Before reduce" << endl;
//        printMatrix(matrix);
        reduce(matrix);
//        cout << "After reduce" << endl;
//        printMatrix(matrix);
//        cout << "\n\n" << endl;

        if (matrix[matrix.size()-1].none()) {
            matrix.pop_back();
            exponentsRow.remove_if(isMultiple<m>(alpha));
        }
        else {
            exponentsRow.pop_front();
            delta_size++;
        }
    }
}


/******************************************************
*
*    RST SIMPLE EXPERIMENT
*
*******************************************************/
template<int m>
void MultiDimArray_GF2<m>::RST_simpleNEW() {
    delta_size = 0;
    lead_monomials.empty();
    groebner_basis.empty();

    Monomial<m> alpha; //dummy container
    Monomial<m> e_period(period);
    Monomial<m> bound;
    for (int i = 0; i < m; i++) {
        bound.get_exponent()[i] = 2 * period[i] - 1;
    }

    const unsigned long width = sizeOfDivisors(bound);

    vector<Monomial<m> > exponentsColumn(width);
    forward_list<Monomial<m> > exponentsRow;

    int index = 0;
    generateDivisors<m>(m, index, alpha, bound, exponentsColumn);
    sort(exponentsColumn.begin(), exponentsColumn.end(), ordering_1<m>);
    for (auto e = exponentsColumn.rbegin(); e != exponentsColumn.rend(); e++) {
        if (e->leq_d(e_period))
            exponentsRow.push_front(*e);
    }

    unordered_map<unsigned long, boost::dynamic_bitset<>*> pivot;
    boost::dynamic_bitset<>* last_row;

    while (!exponentsRow.empty()) {
        alpha = exponentsRow.front();
        last_row = rowPiAlphaNEW<m>(alpha, A, exponentsColumn, e_period);
        boost::dynamic_bitset<>::size_type j = last_row->find_first();
        bool new_pivot = false;
        while (!new_pivot and j != boost::dynamic_bitset<>::npos) {
            if (pivot.count(j)) {
                (*last_row) ^= *(pivot[j]);
                j = last_row->find_next(j);
            }
            else {
                pivot[j] = last_row;
                new_pivot = true;
            }
        }
        if (j == boost::dynamic_bitset<>::npos) { // implies last row is of all zeros
            delete last_row;
            exponentsRow.remove_if(isMultiple<m>(alpha));
        }
        else {
            exponentsRow.pop_front();
            delta_size++;
        }
    }

    for (auto &x : pivot) delete x.second;
}



template<int m>
void MultiDimArray_GF2<m>::RST_simpleNEW2() {
    delta_size = 0;
    lead_monomials.empty();
    groebner_basis.empty();

    Monomial<m> alpha; //dummy container
    Monomial<m> e_period(period);
    Monomial<m> bound;
    for (int i = 0; i < m; i++) {
        bound.get_exponent()[i] = 2 * period[i] - 1;
    }

    const unsigned long width = sizeOfDivisors(bound);

    vector<Monomial<m> > exponentsColumn(width);
    forward_list<Monomial<m> > exponentsRow;

    int index = 0;
    generateDivisors<m>(m, index, alpha, bound, exponentsColumn);
    sort(exponentsColumn.begin(), exponentsColumn.end(), ordering_1<m>);
    for (auto e = exponentsColumn.rbegin(); e != exponentsColumn.rend(); e++) {
        if (e->leq_d(e_period))
            exponentsRow.push_front(*e);
    }

    unordered_map<unsigned long, unique_ptr<boost::dynamic_bitset<> > > pivot;

    while (!exponentsRow.empty()) {
        alpha = exponentsRow.front();
        unique_ptr<boost::dynamic_bitset<> > last_row(rowPiAlphaNEW<m>(alpha, A, exponentsColumn, e_period));
        boost::dynamic_bitset<>::size_type j = last_row->find_first();
        bool new_pivot = false;
        while (!new_pivot and j != boost::dynamic_bitset<>::npos) {
            if (pivot.count(j)) {
                (*last_row) ^= *(pivot[j]);
                j = last_row->find_next(j);
            }
            else {
                pivot[j] = std::move(last_row);
                new_pivot = true;
            }
        }
        if (j == boost::dynamic_bitset<>::npos) { // implies last row is of all zeros
            exponentsRow.remove_if(isMultiple<m>(alpha));
        }
        else {
            exponentsRow.pop_front();
            delta_size++;
        }
    }
}
