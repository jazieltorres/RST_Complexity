#ifndef MULTIDIMARRAY_H
#define MULTIDIMARRAY_H

#include <iostream>
#include <vector>
#include <forward_list>
#include "MExponent2.cpp"
#include <algorithm>
using namespace std;

template <typename F, unsigned long m>
class MultiDimArray {
    private:
        blitz::Array<F, m> A;
        vector<unsigned long> period;
        unsigned long delta_size;
    public:
        explicit MultiDimArray(blitz::Array<F,m>&);
        static const unsigned long dimension = m;
        void RST();
} ;

#endif

template <typename F, unsigned long m>
MultiDimArray<F,m>::MultiDimArray(blitz::Array<F,m>& array){
    A.reference(array);
    period.resize(m);
    for (unsigned long i=0; i<m; i++){
        period[i] = A.extent(i);
    }
    delta_size = 0;
}


/******************************************************
*
*       OTHER MODULES FOR RST IMPLEMENTATION
*
*******************************************************/

template<unsigned long m>
void generateDivisors(unsigned long depth, unsigned long& index, MExponent<m>& e, const MExponent<m>& bound, vector< MExponent<m> >& list){
  if (depth > 0){
    for(unsigned long i = 0; i <= bound.getExp()[depth-1]; i++) {
        e.editExp()[depth-1] = i;
        generateDivisors(depth-1, index, e, bound, list);
    }
  }
  else {
    list[index] = e;
    index++;
  }
}

template<unsigned long m>
unsigned long sizeOfDivisors(const MExponent<m> e) {
    unsigned long size(1);
    for (auto x : e.getExp()) {
        size = size * (x + 1);
    }
    return size;
}

template<unsigned long m>
bool comp(const MExponent<m>& e1, const MExponent<m>& e2) {
    return (e1.grlex_less(e2));
}

template<typename F>
void insertBefore(vector< vector<F> >& M, unsigned long& index, const unsigned long& m) {
    for (unsigned long i = 0; i < (m - index); i++) {
        M[m-i].swap(M[m-i-1]) ;
    }
}

template<typename F>
void reduce(vector< vector<F> >& M, vector< vector<F> >& Id) {
    const unsigned long m = M.size() - 1 ;
    const unsigned long n = M[0].size() ;
    F lambda ;
    long j = 0 ;

    for (unsigned long i = 0; i < m; i++) {
        while (M[i][j] == 0 && j < n) {
            if (M[m][j] != 0) {
                insertBefore(M, i, m);
                insertBefore(Id, i, m);
                return;
            }
            j++;
        }
        lambda = -M[m][j]/M[i][j];
        for (unsigned long k = j; k < n; k++)
            M[m][k] = M[m][k] + lambda * M[i][k];
        for (unsigned long k = 0; k < Id[0].size(); k++)
            Id[m][k] = Id[m][k] + lambda * Id[i][k];
        j++;
    }
} // end reduce

template<typename F>
void printMatrix (const vector< vector<F> >& M) {
    unsigned long m = M.size();
    unsigned long n = M[0].size();

    for (unsigned long i = 0; i < m; i++){
        for (unsigned long j = 0; j < n; j++){
            cout << M[i][j] << " ";
        }
        cout << endl;
    }
}

template<typename F, unsigned long m>
vector<F> rowPiAlpha(const MExponent<m>& alpha, const blitz::Array<F,m>& A, vector< MExponent<m> >& exponentsColumn, MExponent<m> e_period) {
    unsigned long n = exponentsColumn.size();
    vector<F> piAlpha;
    piAlpha.resize(n);
    for (unsigned long i = 0; i < n; i++) {
        MExponent<m> index = exponentsColumn[i]+alpha;
        piAlpha[i] = A(index.mod(e_period).getExp());
    }
    return piAlpha;
}

template<typename F>
void addDimension(vector< vector<F> >& idMatrix) {
    unsigned long n;
    if(idMatrix.empty())    n = idMatrix.size();
    else                    n = idMatrix[0].size();

    vector<F> lastRow(n+1, (F)0);
    lastRow[n] = (F)1;
    for (unsigned long i = 0; i < n; i++){
        idMatrix[i].push_back((F)0);
    }
    idMatrix.push_back(lastRow);
}

// Functor used in exponentsRow.remove_if for removing all the multiples of alpha
template<unsigned long m>
class isMultiple {
private:
    MExponent<m> alpha;
public:
    explicit isMultiple(MExponent<m>& e) {alpha = e;}
    bool operator () (const MExponent<m>& e) const {
        return alpha.leq_d(e);
    }
};

template<typename F, unsigned long m>
void printPoly(vector< vector<F> >& idMatrix, vector< MExponent<m> >& exponentsColumn) {
    unsigned long n = idMatrix.size()-1;
    for (unsigned long i = 0; i < idMatrix[0].size(); i++){
        if(idMatrix[n][i] != (F)0) cout << idMatrix[n][i] << "x" << exponentsColumn[i] << "\t";
    }
    cout << endl;
}




/******************************************************
*
*       THE ONE AND ONLY: RST
*
*******************************************************/
template<typename F, unsigned long m>
void MultiDimArray<F,m>::RST() {
    MExponent<m> bound;
    for (unsigned long i = 0; i < m; i++) {
        bound.editExp()[i] = 2 * period[i] - 1;
    }

    vector< MExponent<m> > exponentsColumn(sizeOfDivisors(bound));
    forward_list< MExponent<m> > exponentsRow;

    unsigned long index = 0;
    MExponent<m> exp;
    generateDivisors<m>(m, index, exp, bound, exponentsColumn);
    sort(exponentsColumn.begin(), exponentsColumn.end(), comp<m>);

    MExponent<m> e_period(period);
    for (auto e = exponentsColumn.rbegin(); e != exponentsColumn.rend(); e++) {
        if (e->leq_d(e_period)) {exponentsRow.push_front (*e);}
    }

    vector< vector<F> > matrix;
    vector< vector<F> > idMatrix;

    while (!exponentsRow.empty()) {
        MExponent<m> alpha = exponentsRow.front();
        matrix.push_back(rowPiAlpha<F,m>(alpha, A, exponentsColumn, e_period));
        addDimension(idMatrix);


        // // // TEST
        // exponentsRow.pop_front();
        // alpha = exponentsRow.front();
        // matrix.push_back(rowPiAlpha<F,m>(alpha, A, exponentsColumn, e_period));
        // addDimension(idMatrix);
        // reduce(matrix, idMatrix);
        //
        // exponentsRow.pop_front();
        // alpha = exponentsRow.front();
        // matrix.push_back(rowPiAlpha<F,m>(alpha, A, exponentsColumn, e_period));
        // addDimension(idMatrix);
        //
        //
        //
        // cout << "Matrix before:\n";
        // printMatrix(matrix);
        // cout << "\n\n\n";
        // cout << "ID Matrix before:\n";
        // printMatrix(idMatrix);
        // cout << "\n\n\n";
        //
        // reduce(matrix, idMatrix);
        //
        // cout << "Matrix after:\n";
        // printMatrix(matrix);
        // cout << "\n\n\n";
        // cout << "ID Matrix after:\n";
        // printMatrix(idMatrix);
        // cout << "\n\n\n";
        // // END TEST

        reduce(matrix, idMatrix);
        vector<F> zeroRow(matrix[0].size(), (F)0);
        if (matrix[matrix.size()-1] == zeroRow) {
            printPoly(idMatrix, exponentsColumn);
            matrix.pop_back();
            idMatrix.pop_back();
            exponentsRow.remove_if(isMultiple<m>(alpha));
        }
        else {
            exponentsRow.pop_front();
            delta_size++;
        }

    }
}
