#include <iostream>
#include <vector>
#include <forward_list>
#include "MExponent.h"
using namespace std;

#ifndef MULTIDIMARRAY_H
#define MULTIDIMARRAY_H

class MultiVarMonomial {
    double coeff;
    MExponent exp;
    MultiVarMonomial () {}
    MultiVarMonomial (double& c) {coeff = c;}
    MultiVarMonomial (MExponent& e) {exp = e;}
    MultiVarMonomial (double& c, MExponent& e) {coeff = c; exp = e;}
    void display(ostream& out) {out << coeff << "*" << exp;}
};

ostream& operator<<(ostream& out, const MultiVarMonomial&);

class MultiDimArray
{
    private:
        vector< vector<double> > A;
        long delta_size;
        vector<MultiVarMonomial> grobner_basis;
    public:
        void RST(const vector< vector<double> >&, const long&, const vector<long>&);
} ;
#endif

/******************************************************
*
*       OTHER MODULES FOR RST IMPLEMENTATION
*
*******************************************************/

void generateDivisors(long depth, long& index, MExponent& e, const vector<long>& period, vector<MExponent>& list){
  if (depth > 0){
    for(int i=0; i <= period[depth-1]; i++){
        e.editExp()[depth-1] = i;
        generateDivisors(depth-1, index, e, period, list);
    }
  }
  else {
    list[index] = e;
    index++;
  }
}


long sizeOfDivisors(const MExponent& e) {
    long size(1);
    for (auto it = e.getExp().begin(); it != e.getExp().end(); it++) {
        size = size * (*it + 1);
    }
    return size;
}

bool comp(const MExponent& e1, const MExponent& e2) {
    return (e1.grlex_less(e2));
}

void printM (const vector< vector<double> >& M) {
    for (long i = 0; i<M.size(); i++){
        for (long j=0; j<M[0].size(); j++){
            cout << M[i][j] << " ";
        }
        cout << endl;
    }
}

void insertBefore(vector< vector<double> >& M, long& index, const long& m) {
    for (long i = 0; i < (m - index); i++) {
        M[m-i].swap(M[m-i-1]) ;
    }
}


void reduce(vector< vector<double> >& M, vector< vector<double> >& Id) {
    const long m = M.size() - 1;
    const long n = M[0].size() - 1;
    const long n_id = Id[0].size() - 1;
    vector<double> lambda(m, 0);
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

vector<double> rowPiAlpha(const MExponent& alpha, const vector< vector<double> >& A, vector<MExponent>& exponentsColumn) {
    vector<double> piAlpha;

    return piAlpha;
}


/******************************************************
*
*       THE ONE AND ONLY: RST
*
*******************************************************/

void MultiDimArray::RST(const vector< vector<double> >& A, const long& m, const vector<long>& period) {
    MExponent e_period(period);
    MExponent e_bound(period);

    for (auto i = e_bound.editExp().begin(); i != e_bound.editExp().end(); i++) {
        *i = (2 * *i) - 1;
    }

    vector<MExponent> exponentsColumn(sizeOfDivisors(e_bound));
    forward_list<MExponent> exponentsRow;

    sort(exponentsColumn.begin(), exponentsColumn.end(), comp);

    for (auto e = exponentsColumn.rbegin(); e != exponentsColumn.rend(); e++) {
        if (e->leq_d(e_period)) {exponentsRow.push_front (*e);}
    }

    vector< vector<double> > matrix;// Later this will be ZZ_p from NTL
    vector< vector<double> > id_matrix;
    long delta_size(0);

    while (!exponentsRow.empty()) {
        MExponent alpha = exponentsRow.front();
        matrix.push_back(rowPiAlpha(alpha, A, exponentsColumn));
    }


}
