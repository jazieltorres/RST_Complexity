#ifndef MULTIDIMARRAY_GF2_H
#define MULTIDIMARRAY_GF2_H

#include <iostream>
#include <functional>
#include <vector>
#include <string>
#include <forward_list>
#include "MultivarPolynomial.cpp"
#include <algorithm> //sort
#include "NTL/GF2.h"
#include "NTL/mat_GF2.h"
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
    void RST_general(int);
    void RST_optimized(int);


public:
//  Constructor: receives a blitz array of dimension m, with entries in GF2.
    explicit MultiDimArray_GF2(blitz::Array<GF2,m>&);

//  Constructor: receives period vector and resizes array
    explicit MultiDimArray_GF2(const blitz::TinyVector<int,m>&);

//  Constructor: receives shift sequence, column sequence, and their periods, respectively.
    MultiDimArray_GF2(const function<int (int)>&, const function<GF2 (int)>&, int, int);

//  Constructor: same as above but shift sequence is a vector.
    MultiDimArray_GF2(const vector<int>&, GF2, const function<GF2 (int)>&, int);

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
    void RST_simple();
} ;

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
        this->RST_optimized(1);
        return delta_size;
    }
    else
        return delta_size;
}

template <int m>
double MultiDimArray_GF2<m>::normalized_complexity() {
    if (delta_size == -1) {
        this->RST_optimized(1);
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
int sizeOfDivisors(const Monomial<m> e) {
    int size(1);
    for (auto x : e.exponent()) {
        size = size * (x + 1);
    }
    return size;
}

void insertBefore(vector< vector<GF2> >& M, int& index, const int& m) {
    for (int i = 0; i < (m - index); i++) {
        M[m-i].swap(M[m-i-1]) ;
    }
}

// Row reduction with the identity matrix to obtain polynomial in the basis
void reduce(vector< vector<GF2> >& M, vector< vector<GF2> >& Id, int column_bound) {
    const int m = M.size() - 1 ;
    GF2 lambda ;
    int j = 0 ;

    for (int i = 0; i < m; i++) {
        while (M[i][j] == 0 && j < column_bound) {
            if (M[m][j] != 0) {
                insertBefore(M, i, m);
                insertBefore(Id, i, m);
                return;
            }
            j++;
        }
        lambda = -M[m][j]/M[i][j];
        for (int k = j; k < column_bound; k++)
            M[m][k] = M[m][k] + lambda * M[i][k];
        for (int k = 0; k < Id[0].size(); k++)
            Id[m][k] = Id[m][k] + lambda * Id[i][k];
        j++;
    }
}

// Row reduction just for computing the complexity
void reduce(vector< vector<GF2> >& M, int column_bound) {
    const int m = M.size() - 1 ;
    GF2 lambda ;
    int j = 0 ;

    for (int i = 0; i < m; i++) {
        while (M[i][j] == 0 && j < column_bound) {
            if (M[m][j] != 0) {
                insertBefore(M, i, m);
                return;
            }
            j++;
        }
        lambda = -M[m][j]/M[i][j];
        for (int k = j; k < column_bound; k++)
            M[m][k] = M[m][k] + lambda * M[i][k];
        j++;
    }
}

void printMatrix (const vector< vector<GF2> >& M) {
    int m = M.size();
    int n = M[0].size();

    if (n>15) n = 15;

    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            cout << M[i][j] << " ";
        }
        cout << endl;
    }
}

template<int m>
vector<GF2> rowPiAlpha(const Monomial<m>& alpha, const blitz::Array<GF2,m>& A,
                     vector< Monomial<m> >& exponentsColumn, Monomial<m> period) {
    int n = exponentsColumn.size();
    vector<GF2> piAlpha;
    piAlpha.resize(n);
    for (int i = 0; i < n; i++) {
        Monomial<m> index = exponentsColumn[i]+alpha;
        piAlpha[i] = A(index.mod(period).exponent());
    }
    return piAlpha;
}

void addDimension(vector< vector<GF2> >& id_matrix) {
    int n, m;
    m = id_matrix.size();
    if (m == 0) {
        n = 0;
    }
    else {
        n = id_matrix[0].size();
    }

    vector<GF2> lastRow(n+1, (GF2)0);
    lastRow[n] = (GF2)1;
    for (int i = 0; i < m ; i++){
        id_matrix[i].push_back((GF2)0);
    }
    id_matrix.push_back(lastRow);
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
void printPoly(vector< vector<GF2> >& id_matrix, vector< Monomial<m> >& id_column) {
    int n = id_matrix.size()-1;
    for (int i = 0; i < id_matrix[0].size(); i++) {
        if(id_matrix[n][i] != (GF2)0) cout << id_matrix[n][i] << "x" << id_column[i] << "\t";
    }
    cout << endl;
}

template<int m>
MultivarPolynomial<GF2,m> get_polynomial(vector< vector<GF2> >& id_matrix, vector< Monomial<m> >& id_column) {
    MultivarPolynomial<GF2,m> poly;
    int n = id_matrix.size()-1;
    for (int i = 0; i < id_matrix[0].size(); i++) {
        if(id_matrix[n][i] != (GF2)0) poly.add_term(id_column[i], id_matrix[n][i]);
    }
    return poly;
}

bool isZeroRow (vector< vector<GF2> >& matrix, int column_bound) {
    int lastRow = matrix.size()-1;
    GF2 zero(0);
    for (int i=0; i<column_bound; i++) {
        if (matrix[lastRow][i] != zero)
            return false;
    }
    return true;
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

    vector<Monomial<m> > exponentsColumn(sizeOfDivisors(bound));
    vector<Monomial<m> > id_column;
    forward_list<Monomial<m> > exponentsRow;

    int index = 0;
    generateDivisors<m>(m, index, alpha, bound, exponentsColumn);

    sort(exponentsColumn.begin(), exponentsColumn.end(), ordering_1<m>);

//    After each iteration using grlex, we have to check up to (1 + the number of monomials that are in
//    exponentColumn but are jumped in exponentRow) less in row reduction.
//    Hence, we count the number of monomials that are jumped and, in the same loop, fill the list
//    exponentRow in sorted order.

//    skipped_monomials[i] will contain the number of monomials between exponentRow[i] and
//    exponentRow[i+1] (excluding both) that are in exponentColumn.
    vector<int> skipped_monomials(sizeOfDivisors(e_period), 0);

    auto it_period = exponentsColumn.rbegin();
    auto it_skipped = skipped_monomials.rbegin() - 1;
    while (!it_period->equal(e_period)) //Skip all monomials greater than e_period
        it_period++;
    for (auto e = it_period; e != exponentsColumn.rend(); e++) {
        if (e->leq_d(e_period)) {
            exponentsRow.push_front(*e);
            it_skipped++;
        } else
            *it_skipped += 1;
    }

    vector<vector<GF2> > matrix;
    int column_bound = exponentsColumn.size() + 1;
    int iteration = 0;
    while (!exponentsRow.empty()) {
        alpha = exponentsRow.front();
        matrix.push_back(rowPiAlpha<m>(alpha, A, exponentsColumn, e_period));
        id_column.push_back(alpha);
        column_bound = column_bound - 1 - skipped_monomials[iteration];
        reduce(matrix, column_bound);
        if (isZeroRow(matrix, column_bound)) {
            matrix.pop_back();
            exponentsRow.remove_if(isMultiple<m>(alpha));
        } else {
            exponentsRow.pop_front();
            delta_size++;
        }
        iteration++;
//        if ((delta_size%100)==0) cerr << "Going through... " << delta_size << endl;
    }
}

/******************************************************
*
*             THE ONE AND ONLY: RST
*
*******************************************************/
template<int m>
void MultiDimArray_GF2<m>::RST_general(int ordering) {
    delta_size = 0;
    lead_monomials.empty();
    groebner_basis.empty();

    Monomial<m>  alpha; //dummy container
    Monomial<m>  e_period(period);
    Monomial<m>  bound;
    for (int i = 0; i < m; i++) {
        bound.get_exponent()[i] = 2 * period[i] - 1;
    }

    vector< Monomial<m> >  exponentsColumn(sizeOfDivisors(bound));
    vector< Monomial<m> >  id_column;
    forward_list< Monomial<m> >  exponentsRow;

    int index = 0;
    generateDivisors<m>(m, index, alpha, bound, exponentsColumn);

    switch (ordering) {
        case 3 :
            sort(exponentsColumn.begin(), exponentsColumn.end(), ordering_3<m>);
            break;
        case 4 :
            sort(exponentsColumn.begin(), exponentsColumn.end(), ordering_4<m>);
            break;
        case 5 :
            sort(exponentsColumn.begin(), exponentsColumn.end(), ordering_1<m>);
            break;
        case 6 :
            sort(exponentsColumn.begin(), exponentsColumn.end(), ordering_2<m>);
            break;
        default:
            cerr << "Invalid ordering argument in RST." << endl;
            delta_size = -1;
            return;
    }

    for (auto e = exponentsColumn.rbegin(); e != exponentsColumn.rend(); e++) {
        if (e->leq_d(e_period)) {exponentsRow.push_front (*e);}
    }

    vector< vector<GF2> >  matrix;
    vector< vector<GF2> >  id_matrix;
    vector< vector <Monomial<m> > >  basis;
    while (!exponentsRow.empty()) {
        alpha = exponentsRow.front();
        matrix.push_back(rowPiAlpha<m>(alpha, A, exponentsColumn, e_period));
        addDimension(id_matrix);
        id_column.push_back(alpha);
        reduce(matrix, id_matrix, exponentsColumn.size());
        if (isZeroRow(matrix, exponentsColumn.size())) {
            lead_monomials.push_back(alpha);
//            printPoly(id_matrix, id_column);
            groebner_basis.push_back(get_polynomial(id_matrix, id_column));
            matrix.pop_back();
            id_matrix.pop_back();
            exponentsRow.remove_if(isMultiple<m>(alpha));
        }
        else {
            exponentsRow.pop_front();
            delta_size++;
        }
//        if ((delta_size%100)==0) cout << "Going through... " << delta_size << endl;
    }
}


/******************************************************
*
*           RST OPTIMIZED FOR GRLEX
*
*******************************************************/
template<int m>
void MultiDimArray_GF2<m>::RST_optimized(int ordering) {
    delta_size = 0;
    lead_monomials.empty();
    groebner_basis.empty();

    Monomial<m>  alpha; //dummy container
    Monomial<m>  e_period(period);
    Monomial<m>  bound;
    for (int i = 0; i < m; i++) {
        bound.get_exponent()[i] = 2 * period[i] - 1;
    }

    vector< Monomial<m> >  exponentsColumn(sizeOfDivisors(bound));
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
        default:
            cerr << "Invalid ordering argument in RST." << endl;
            delta_size = -1;
            return;
    }

//    After each iteration using grlex, we have to check up to (1 + the number of monomials that are in 
//    exponentColumn but are jumped in exponentRow) less in row reduction.
//    Hence, we count the number of monomials that are jumped and, in the same loop, fill the list
//    exponentRow in sorted order.

//    skipped_monomials[i] will contain the number of monomials between exponentRow[i] and 
//    exponentRow[i+1] (excluding both) that are in exponentColumn.
    vector<int> skipped_monomials(sizeOfDivisors(e_period), 0);

    auto it_period = exponentsColumn.rbegin();
    auto it_skipped = skipped_monomials.rbegin()-1;
    while (!it_period->equal(e_period)) //Skip all monomials greater than e_period
        it_period++;
    for (auto e = it_period; e != exponentsColumn.rend(); e++) {
        if (e->leq_d(e_period)) {
            exponentsRow.push_front(*e);
            it_skipped++;
        }
        else
            *it_skipped += 1;
    }

    vector< vector<GF2> >  matrix;
    vector< vector<GF2> >  id_matrix;
    vector< vector <Monomial<m> > >  basis;
    int column_bound = exponentsColumn.size()+1;
    int iteration = 0;
    while (!exponentsRow.empty()) {
        alpha = exponentsRow.front();
        matrix.push_back(rowPiAlpha<m>(alpha, A, exponentsColumn, e_period));
        addDimension(id_matrix);
        id_column.push_back(alpha);
        column_bound = column_bound - 1 - skipped_monomials[iteration];
        reduce(matrix, id_matrix, column_bound);
        if (isZeroRow(matrix, column_bound)) {
            lead_monomials.push_back(alpha);
//            printPoly(id_matrix, id_column);
            groebner_basis.push_back(get_polynomial(id_matrix, id_column));
            matrix.pop_back();
            id_matrix.pop_back();
            exponentsRow.remove_if(isMultiple<m>(alpha));
        }
        else {
            exponentsRow.pop_front();
            delta_size++;
        }
        iteration++;
//        if ((delta_size%100)==0) cout << "Going through... " << delta_size << endl;
    }
//
// Printing leading monomials in the array, whom closes the monomials in the delta set
//    cout << "Leading Monomials " << endl;
//    printStars(lead_monomials, period);
//
//    cout << "\n\n\n groebner basis" << endl;
//    for (auto poly : basis) {
//        printStars(poly, period);
//        cout << endl;
//    }
}


/****************************************************************************
*
*           RST GATEWAY FUNCTION
*
* Calls RST with respect to a monomial ordering:
*   1. Graded lexicographic X1 > X2 > ... > Xn with ladder optimization
*   2. Graded lexicographic X1 < X2 < ... < Xn with ladder optimization
*   3. Lexicographic X1 > X2 > ... > Xn
*   4. Lexicographic X1 < X2 < ... < Xn
*   5. Graded lexicographic X1 > X2 > ... > Xn
*   6. Graded lexicographic X1 < X2 < ... < Xn
****************************************************************************/
template<int m>
void MultiDimArray_GF2<m>::RST(int ordering) {
    if (ordering < 1 or ordering > 6) {
        cerr << "Invalid argument." << endl;
        return;
    }
    ordering_number = ordering;
    if (ordering <= 2)
        RST_optimized(ordering);
    else if (ordering <= 6)
        RST_general(ordering);

}