#ifndef MULTIDIMARRAY_H
#define MULTIDIMARRAY_H

#include <iostream>
#include <vector>
#include <forward_list>
#include <NTL/ZZ_p.h>
#include "MExponent2.cpp"
using namespace std;


// class MultiVarMonomial {
//     private:
//         NTL::ZZ_p coeff;
//         MExponent exp;
//     public:
//         MultiVarMonomial () {}
//         MultiVarMonomial (NTL::ZZ_p& c) {coeff = c;}
//         MultiVarMonomial (MExponent& e) {exp = e;}
//         MultiVarMonomial (NTL::ZZ_p& c, MExponent& e) {coeff = c; exp = e;}
//         void display(ostream& out) {out << coeff << "*" << exp;}
// };
//
// ostream& operator<<(ostream& out, const MultiVarMonomial&);

template <typename F, long m>
class MultiDimArray {
    private:
        blitz::Array<F, m> A;
        vector<long> period;
        long delta_size;
        // vector<MultiVarMonomial> grobner_basis;
    public:
        MultiDimArray() {}
        MultiDimArray(blitz::Array<F, m>&);
        static const long dimension = m;
        void RST();
} ;

#endif

template <typename F, long m>
MultiDimArray<F,m>::MultiDimArray(blitz::Array<F,m>& array){
    A.reference(array);
    period.resize(m);
    for (long i=0; i<m; i++){
        period[i] = A.extent(i);
    }
}


/******************************************************
*
*       OTHER MODULES FOR RST IMPLEMENTATION
*
*******************************************************/

template<long m>
void generateDivisors(long depth, long index, MExponent<m>& e, const vector<long>& period, vector< MExponent<m> >& list){
  if (depth > 0){
    for(int i=0; i <= period[depth-1]; i++){
        e.editExp()[depth-1] = i;
        generateDivisors(depth-1, index, e, period, list);
    }
  }
  else {
    list[index] = e;
    cout << &list[index] << endl ;
    index++;
  }
}

template<long m>
long sizeOfDivisors(const MExponent<m>& e) {
    long size(1);
    for (auto it = e.getExp().begin(); it != e.getExp().end(); it++) {
        size = size * (*it + 1);
    }
    return size;
}

template<long m>
bool comp(const MExponent<m>& e1, const MExponent<m>& e2) {
    return (e1.grlex_less(e2));
}

template<typename F>
void insertBefore(vector< vector<F> >& M, long& index, const long& m) {
    for (long i = 0; i < (m - index); i++) {
        M[m-i].swap(M[m-i-1]) ;
    }
}

template<typename F>
void reduce(vector< vector<F> >& M, vector< vector<F> >& Id) {
    const long m = M.size() - 1;
    const long n = M[0].size() - 1;
    const long n_id = Id[0].size() - 1;
    vector<F> lambda(m, (F)0);
    long j = 0 ;

    for (long i = 0; i < m; i++) {
        while (M[i][j] == 0 && j <= n) {
            for (long k = 0; k < i; k++) {
                M[m][j] += lambda[k] * M[k][j];
                if (j <= n_id) {Id[m][j] += lambda[k] * Id[k][j];}
            }
            if (M[m][j] != 0) {
                for (long l = j+1; l <= n; l++) {
                    for (long k = 0; k < i; k++) {
                        M[m][l] += lambda[k] * M[k][l];
                        if (l <= n_id) {Id[m][l] += lambda[k] * Id[k][l];}
                    }
                }
                insertBefore(M, i, m);
                insertBefore(Id, i, m);
                return;
            }
            j++;
        }
        lambda[i] = M[m][j];
        for (long k = 0; k < i; k++) {
            lambda[i] += lambda[k] * M[k][j];
        }
        lambda[i] = -lambda[i]/M[i][j];
        M[m][j] = 0;
        j++;
    }

    for (long l = j; l <= n; l++) {
        for (long k = 0; k < m; k++) {
            M[m][l] += lambda[k] * M[k][l];
            if (l <= n_id) {Id[m][l] += lambda[k] * Id[k][l];}
        }
    }

    return;
} // end reduce

template<typename F>
void printMatrix (const vector< vector<F> >& M) {
    long m = M.size();
    long n = M[0].size();

    for (long i = 0; i<m; i++){
        for (long j=0; j<n; j++){
            cout << M[i][j] << " ";
        }
        cout << endl;
    }
}

template<typename F, long m>
vector<F> rowPiAlpha(const MExponent<m>& alpha, const blitz::Array<F,m>& A, vector< MExponent<m> >& exponentsColumn, MExponent<m> e_period) {
    long n = exponentsColumn.size();
    vector<F> piAlpha;
    piAlpha.resize(n);
    for (long i=0; i<n; i++) {
        MExponent<m> index = exponentsColumn[i]+alpha;
        // cout << alpha << endl;
        // cout << exponentsColumn[i] << " + " << alpha << " = " << index << endl;
        // piAlpha[i] = A(index.mod(e_period).getExp());
    }
    return piAlpha;
}


/******************************************************
*
*       THE ONE AND ONLY: RST
*
*******************************************************/
template<typename F, long m>
void MultiDimArray<F,m>::RST() {
    MExponent<m> e_period(period);
    MExponent<m> e_bound(period);
    for (auto i = e_bound.editExp().begin(); i != e_bound.editExp().end(); i++) {
        *i = (2 * *i) - 1;
    }

    vector< MExponent<m> > exponentsColumn(sizeOfDivisors(e_bound));
    forward_list< MExponent<m> > exponentsRow;

    generateDivisors<m>(m, 0,e_bound,period,exponentsColumn);

    cout << "begin exponentsColumn" << endl ;
    for (auto e = exponentsColumn.begin(); e != exponentsColumn.end(); e++){
        cout << *e << endl ;
    }

    sort(exponentsColumn.begin(), exponentsColumn.end(), comp<m>);


    for (auto e = exponentsColumn.rbegin(); e != exponentsColumn.rend(); e++) {
        // cout << *e << endl ;
        if (e->leq_d(e_period)) {exponentsRow.push_front (*e);}
    }

    cout << "end exponentsColumn" << endl ;

    vector< vector<F> > matrix;
    vector< vector<F> > id_matrix;
    long delta_size(0);

    // while (!exponentsRow.empty()) {
        MExponent<m> alpha = exponentsRow.front();
        matrix.push_back(rowPiAlpha<F,m>(alpha, A, exponentsColumn, e_period));
        // exponentsRow.pop_front();
        // alpha = exponentsRow.front();
        // matrix.push_back(rowPiAlpha<F,m>(alpha, A, exponentsColumn, e_period));
        // printMatrix(matrix);
    // }
}
