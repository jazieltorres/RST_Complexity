#ifndef MULTIDIMARRAY_H
#define MULTIDIMARRAY_H

#include <iostream>
#include <vector>
#include <string>
#include <forward_list>
#include "MExponent.cpp"
#include <algorithm> //sort
using namespace std;

template <typename F, int m>
class MultiDimArray {
    private:
        blitz::Array<F, m> A;
        int delta_size;
        blitz::TinyVector<int, m> period;
    public:
//      Constructor: receives a blitz array of dimension m, with entries in F.
        explicit MultiDimArray(blitz::Array<F,m>&);

//      Constructor: receives period vector and resizes array
        explicit MultiDimArray(const blitz::TinyVector<F,m>&);

//      Constructor: receives shift sequence, column sequence, and their periods, respectively.
        MultiDimArray(const function<int (int)>&, const function<F (int)>&, int, int);

//      Constructor: same as above but shift sequence is a vector.
        MultiDimArray(const vector<int>&, const function<F (int)>&, int, int);

        static const int dimension = m;
        void RST();
//        Complexity
        int getDeltaSize(){return delta_size;}
        void setAt(const blitz::TinyVector<int,m>&, F&);
        void print();
        void printPeriod();
} ;

#endif



/******************************************************
*
*       CONSTRUCTORS
*
*******************************************************/


//  Constructor: receives a blitz array of dimension m, with entries in F.
template <typename F, int m>
MultiDimArray<F,m>::MultiDimArray(blitz::Array<F,m>& array){
    A.reference(array);
    for (int i=0; i<m; i++){
        period[i] = A.extent(i);
    }
    delta_size = 0;
}

//  Constructor: receives period vector and resizes array
template<typename F, int m>
MultiDimArray<F,m>::MultiDimArray(const blitz::TinyVector<F,m>& period_vector){
    delta_size = 0;
    period = period_vector;
    A.resize(period);
}


// MultiDimArray(Shift seq, column seq, horizontal period, column period)
template <typename F, int m>
MultiDimArray<F,m>::MultiDimArray(const function<int(int)>& func1, const function<F(int)>& func2,
        int n1, int n2) {
    if (m != 2) {
        cout << "ERROR: Allowed dimension for this constructor: 2" << endl;
    }
    else {
        period = n1, n2;
        A.resize(period);
        delta_size = 0;
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


// Constructor for permutations
template <typename F, int m>
MultiDimArray<F,m>::MultiDimArray(const vector<int>& func1, const function<F(int)>& func2,
                                  int n1, int n2) {
    if (m != 2) {
        period = n1, n2;
        A.resize(period);
        delta_size = 0;
        blitz::TinyVector<int, 2> index;
        for (int i = 0; i < n1; i++) {
            for (int j = 0; j < n2; j++) {
                index = i, j;
                A(index) = func2(((j - func1[i]) % n2 + n2) % n2);
            }
        }
    }
    else {
        cout << "ERROR: Allowed dimension for this constructor: 2" << endl;
    }
}


// Constructor for logWelch with column of zeros at the end
//template <typename F, int m>
//MultiDimArray<F,m>::MultiDimArray(const function<int(int)>& func1, const function<F(int)>& func2,
//                                  int n1, int n2) {
//
//    period = {(int)n1+1, (int)n2};
//    cout << "period " << period[0] << " " << period[1] << endl;
//    delta_size = 0;
//    A.resize(n1+1, n2);
//    blitz::TinyVector<int, 2> index;
//    for (int i = 0; i < n1; i++) {
//        for (int j = 0; j < n2; j++) {
//            index = i,j;
//            A(index) = func2(((j-func1(i)) % n2 + n2) % n2);
//        }
//        for (int j=0; j<n2; j++){
//            index = n1, j;
//            A(index) = 0;
//        }
//    }
//}


template <typename F, int m>
void MultiDimArray<F,m>::setAt(const blitz::TinyVector<int, m>& coordinates, F& value) {
    A(coordinates) = value;
}

template <typename F, int m>
void MultiDimArray<F, m>::print() {
    blitz::TinyVector<int, m> index;
    for (int j=period[1]-1; j>=0; j--) {
        for (int i=0; i<period[0]; i++) {
            index = i,j;
            cout << A(index) << " ";
        }
        cout << endl;
    }
}

template <typename F, int m>
void MultiDimArray<F,m>::printPeriod() {
    cout << period << endl;
}


/******************************************************
*
*       OTHER MODULES FOR RST IMPLEMENTATION
*
*******************************************************/

template<int m>
void generateDivisors(int depth, int& index, MExponent<m>& e,
        const MExponent<m>& bound, vector< MExponent<m> >& list){
  if (depth > 0){
    for(int i = 0; i <= bound.getExp()[depth-1]; i++) {
        e.editExp()[depth-1] = i;
        generateDivisors(depth-1, index, e, bound, list);
    }
  }
  else {
    list[index] = e;
    index++;
  }
}

template<int m>
int sizeOfDivisors(const MExponent<m> e) {
    int size(1);
    for (auto x : e.getExp()) {
        size = size * (x + 1);
    }
    return size;
}

template<int m>
bool comp(const MExponent<m>& e1, const MExponent<m>& e2) {
    return (e1.lex_less(e2));
}

template<typename F>
void insertBefore(vector< vector<F> >& M, int& index, const int& m) {
    for (int i = 0; i < (m - index); i++) {
        M[m-i].swap(M[m-i-1]) ;
    }
}

template<typename F>
void reduce(vector< vector<F> >& M, vector< vector<F> >& Id) {
    const int m = M.size() - 1 ;
    const int n = M[0].size() ;
    F lambda ;
    int j = 0 ;

    for (int i = 0; i < m; i++) {
        while (M[i][j] == 0 && j < n) {
            if (M[m][j] != 0) {
                insertBefore(M, i, m);
                insertBefore(Id, i, m);
                return;
            }
            j++;
        }
        lambda = -M[m][j]/M[i][j];
        for (int k = j; k < n; k++)
            M[m][k] = M[m][k] + lambda * M[i][k];
        for (int k = 0; k < Id[0].size(); k++)
            Id[m][k] = Id[m][k] + lambda * Id[i][k];
        j++;
    }
} // end reduce

template<typename F>
void printMatrix (const vector< vector<F> >& M) {
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

template<typename F, int m>
vector<F> rowPiAlpha(const MExponent<m>& alpha, const blitz::Array<F,m>& A,
        vector< MExponent<m> >& exponentsColumn, MExponent<m> period) {
    int n = exponentsColumn.size();
    vector<F> piAlpha;
    piAlpha.resize(n);
    for (int i = 0; i < n; i++) {
        MExponent<m> index = exponentsColumn[i]+alpha;
        piAlpha[i] = A(index.mod(period).getExp());
    }
    return piAlpha;
}

template<typename F>
void addDimension(vector< vector<F> >& idMatrix) {
    int n, m;
    m = idMatrix.size();
    if (m == 0) {
        n = 0;
    }
    else {
        n = idMatrix[0].size();
    }

    vector<F> lastRow(n+1, (F)0);
    lastRow[n] = (F)1;
    for (int i = 0; i < m ; i++){
        idMatrix[i].push_back((F)0);
    }
    idMatrix.push_back(lastRow);
}

// Functor used in exponentsRow.remove_if for removing all the multiples of alpha
template<int m>
class isMultiple {
private:
    MExponent<m> alpha;
public:
    explicit isMultiple(MExponent<m>& e) {alpha = e;}
    bool operator () (const MExponent<m>& e) const {
        return alpha.leq_d(e);
    }
};

template<typename F, int m>
void printPoly(vector< vector<F> >& idMatrix, vector< MExponent<m> >& idColumn) {
    int n = idMatrix.size()-1;
    for (int i = 0; i < idMatrix[0].size(); i++) {
        if(idMatrix[n][i] != (F)0) cout << idMatrix[n][i] << "x" << idColumn[i] << "\t";
    }
    cout << endl;
}

template<typename F, int m>
vector< MExponent<m> > getPoly(vector< vector<F> >& idMatrix, vector< MExponent<m> >& idColumn) {
    vector< MExponent<m> > poly;
    int n = idMatrix.size()-1;
    for (int i = 0; i < idMatrix[0].size(); i++) {
        if(idMatrix[n][i] != (F)0) poly.push_back(idColumn[i]);
    }
    return poly;
}

template <int m>
void printStars(vector< MExponent<m> >& v, vector<int>& period) {
    for(int i=2*period[1]+1; i>0; i--){
        for(int j=0; j<period[0]; j++) {
            if (1 == i % 2) {
                cout << "--";
            }
            else {
                bool yes = false;
                for (auto e : v){
                    if ((e.getExp()[0]== j) && (e.getExp()[1]==i/2-1)) yes = true;
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
*       THE ONE AND ONLY: RST
*
*******************************************************/
template<typename F, int m>
void MultiDimArray<F,m>::RST() {
    MExponent<m> bound;
    for (int i = 0; i < m; i++) {
        bound.editExp()[i] = 2 * period[i] - 1;
    }

    vector< MExponent<m> >  exponentsColumn(sizeOfDivisors(bound)),
                            idColumn;
    forward_list< MExponent<m> > exponentsRow;


    int index = 0;
    MExponent<m> exp;
    MExponent<m> e_period(period);
    generateDivisors<m>(m, index, exp, bound, exponentsColumn);

    sort(exponentsColumn.begin(), exponentsColumn.end(), comp<m>);

    for (auto e = exponentsColumn.rbegin(); e != exponentsColumn.rend(); e++) {
        if (e->leq_d(e_period)) {exponentsRow.push_front (*e);}
    }

    vector< vector<F> > matrix, idMatrix;
    vector< MExponent<m> > leadingMonomials;
    vector< vector <MExponent<m> > > grobnerBasis;


    MExponent<m> alpha;
    while (!exponentsRow.empty()) {
        alpha = exponentsRow.front();

//        cout << "ALPHA:\t" << alpha << endl;

        matrix.push_back(rowPiAlpha<F,m>(alpha, A, exponentsColumn, e_period));
        addDimension(idMatrix);
        idColumn.push_back(alpha);




//        cout << "BEFORE REDUCE:" << endl;
//        printMatrix(matrix);
//        cout << endl;
//        printMatrix(idMatrix);

        reduce(matrix, idMatrix);

//        cout << "\nAFTER REDUCE:" << endl;
//        printMatrix(matrix);
//        cout << endl;
//        printMatrix(idMatrix);
//        cout << "\n" << endl;



        vector<F> zeroRow(matrix[0].size(), (F)0);

        if (matrix[matrix.size()-1] == zeroRow) {
            leadingMonomials.push_back(alpha);
//            printPoly(idMatrix, idColumn);
            grobnerBasis.push_back(getPoly(idMatrix, idColumn));
            matrix.pop_back();
            idMatrix.pop_back();

//            cout << "BEFORE REMOVE" << endl;
//            for (auto i=exponentsRow.begin(); i!=exponentsRow.end(); i++) cout << *i << " ";
//            cout << endl;

            exponentsRow.remove_if(isMultiple<m>(alpha));

//            cout << "AFTER REMOVE" << endl;
//            for (auto i=exponentsRow.begin(); i!=exponentsRow.end(); i++) cout << *i << " ";
//            cout << endl;
        }
        else {
            exponentsRow.pop_front();
            delta_size++;
        }
//        if ((delta_size%100)==0) cout << "Going through... " << delta_size << endl;
    }
//    cout << "Delta size: " <<  delta_size << endl;
//
// Printing leading monomials in the array, whom closes the monomials in the delta set
//    cout << "Leading Monomials " << endl;
//    printStars(leadingMonomials, period);
//
//    cout << "\n\n\n Grobner basis" << endl;
//    for (auto poly : grobnerBasis) {
//        printStars(poly, period);
//        cout << endl;
//    }

}
