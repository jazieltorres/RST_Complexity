#ifndef MULTIDIMARRAY_H
#define MULTIDIMARRAY_H

#include <iostream>
#include <vector>
#include <forward_list>
#include <NTL/ZZ_p.h>
#include "MExponent2.hpp"
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

template <typename F, int m>
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

template <typename F, int dim>
MultiDimArray<F, dim>::MultiDimArray(blitz::Array<F,dim>& array){
    A.reference(array);
}



//
// /******************************************************
// *
// *       OTHER MODULES FOR RST IMPLEMENTATION
// *
// *******************************************************/
//
// void generateDivisors(long depth, long& index, MExponent& e, const vector<long>& period, vector<MExponent>& list){
//   if (depth > 0){
//     for(int i=0; i <= period[depth-1]; i++){
//         e.editExp()[depth-1] = i;
//         generateDivisors(depth-1, index, e, period, list);
//     }
//   }
//   else {
//     list[index] = e;
//     index++;
//   }
// }
//
//
// long sizeOfDivisors(const MExponent& e) {
//     long size(1);
//     for (auto it = e.getExp().begin(); it != e.getExp().end(); it++) {
//         size = size * (*it + 1);
//     }
//     return size;
// }
//
// bool comp(const MExponent& e1, const MExponent& e2) {
//     return (e1.grlex_less(e2));
// }
//
// void printM (const vector< vector<NTL::ZZ_p> >& M) {
//     for (long i = 0; i<M.size(); i++){
//         for (long j=0; j<M[0].size(); j++){
//             cout << M[i][j] << " ";
//         }
//         cout << endl;
//     }
// }
//
// void insertBefore(vector< vector<NTL::ZZ_p> >& M, long& index, const long& m) {
//     for (long i = 0; i < (m - index); i++) {
//         M[m-i].swap(M[m-i-1]) ;
//     }
// }
//
//
// void reduce(vector< vector<NTL::ZZ_p> >& M, vector< vector<NTL::ZZ_p> >& Id) {
//     const long m = M.size() - 1;
//     const long n = M[0].size() - 1;
//     const long n_id = Id[0].size() - 1;
//     vector<NTL::ZZ_p> lambda(m, (NTL::ZZ_p)0);
//     long j = 0 ;
//
//     for (long i = 0; i < m; i++) {
//         while (M[i][j] == 0 && j <= n) {
//             for (long k = 0; k < i; k++) {
//                 M[m][j] += lambda[k] * M[k][j];
//                 if (j <= n_id) {Id[m][j] += lambda[k] * Id[k][j];}
//             }
//             if (M[m][j] != 0) {
//                 for (long l = j+1; l <= n; l++) {
//                     for (long k = 0; k < i; k++) {
//                         M[m][l] += lambda[k] * M[k][l];
//                         if (l <= n_id) {Id[m][l] += lambda[k] * Id[k][l];}
//                     }
//                 }
//                 insertBefore(M, i, m);
//                 insertBefore(Id, i, m);
//                 return;
//             }
//             j++;
//         }
//         lambda[i] = M[m][j];
//         for (long k = 0; k < i; k++) {
//             lambda[i] += lambda[k] * M[k][j];
//         }
//         lambda[i] = -lambda[i]/M[i][j];
//         M[m][j] = 0;
//         j++;
//     }
//
//     for (long l = j; l <= n; l++) {
//         for (long k = 0; k < m; k++) {
//             M[m][l] += lambda[k] * M[k][l];
//             if (l <= n_id) {Id[m][l] += lambda[k] * Id[k][l];}
//         }
//     }
//
//     return;
// } // end reduce
//
//
//
//
//
// vector<F> rowPiAlpha(const MExponent& alpha, const blitz::Array<F, dim>& A, vector<MExponent>& exponentsColumn) {
//     vector<F> piAlpha;
//     return piAlpha;
// }
//
//
// /******************************************************
// *
// *       THE ONE AND ONLY: RST
// *
// *******************************************************/
//
// void MultiDimArray::RST() {
//     MExponent e_period(period);
//     MExponent e_bound(period);
//
//     for (auto i = e_bound.editExp().begin(); i != e_bound.editExp().end(); i++) {
//         *i = (2 * *i) - 1;
//     }
//
//     vector<MExponent> exponentsColumn(sizeOfDivisors(e_bound));
//     forward_list<MExponent> exponentsRow;
//
//     sort(exponentsColumn.begin(), exponentsColumn.end(), comp);
//
//     for (auto e = exponentsColumn.rbegin(); e != exponentsColumn.rend(); e++) {
//         if (e->leq_d(e_period)) {exponentsRow.push_front (*e);}
//     }
//
//     vector< vector<F> > matrix;
//     vector< vector<F> > id_matrix;
//     long delta_size(0);
//
//     while (!exponentsRow.empty()) {
//         MExponent alpha = exponentsRow.front();
//         matrix.push_back(rowPiAlpha(alpha, A, exponentsColumn));
//     }
// }
