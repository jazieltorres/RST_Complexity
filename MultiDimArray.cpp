#ifndef MULTIDIMARRAY_H
#define MULTIDIMARRAY_H

#include <iostream>
#include <vector>
#include <forward_list>
#include "MExponent.cpp"
#include <algorithm> //sort
using namespace std;

template <typename F, unsigned long m>
class MultiDimArray {
    private:
        blitz::Array<F, m> A;
        unsigned long delta_size;
    public:
        vector<unsigned long> period;
        explicit MultiDimArray(blitz::Array<F,m>&);
        MultiDimArray(const function<long (long)>&, const function<F (long)>&, long, long);
        MultiDimArray(const vector<long>&, const function<F (long)>&, long, long);
        static const unsigned long dimension = m;
        void RST();
        unsigned long getDeltaSize(){return delta_size;}
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

// MultiDimArray(Costas, Legendre, Costas period, Legendre period)
template <typename F, unsigned long m>
MultiDimArray<F,m>::MultiDimArray(const function<long(long)>& func1, const function<F(long)>& func2,
        long n1, long n2) {
    period = {(unsigned long)n1, (unsigned long)n2};
    delta_size = 0;
    A.resize(n1, n2);
    blitz::TinyVector<unsigned long, 2> index;
    for (long i = 0; i < n1; i++) {
        for (long j = 0; j < n2; j++) {
            index = i,j;
            if (func1(i) == -1)
                A(index) = 0;
            else
                A(index) = func2(((j-func1(i)) % n2 + n2) % n2);
        }
    }

// FOR PRINTING THE ARRAY
//    for (long j=n2-1; j>=0; j--) {
//        for (long i=0; i<n1; i++) {
//            index = i,j;
//            cout << A(index) << " ";
//        }
//        cout << endl;
//    }
}


// Constructor for permutations
template <typename F, unsigned long m>
MultiDimArray<F,m>::MultiDimArray(const vector<long>& func1, const function<F(long)>& func2,
                                  long n1, long n2) {
    period = {(unsigned long)n1, (unsigned long)n2};
    delta_size = 0;
    A.resize(n1, n2);
    blitz::TinyVector<unsigned long, 2> index;
    for (long i = 0; i < n1; i++) {
        for (long j = 0; j < n2; j++) {
            index = i,j;
            A(index) = func2(((j-func1[i]) % n2 + n2) % n2);
        }
    }

// FOR PRINTING THE ARRAY
//    for (long j=n2-1; j>=0; j--) {
//        for (long i=0; i<n1; i++) {
//            index = i,j;
//            cout << A(index) << " ";
//        }
//        cout << endl;
//    }
}

// Constructor for logWelch with column of zeros at the end
//template <typename F, unsigned long m>
//MultiDimArray<F,m>::MultiDimArray(const function<long(long)>& func1, const function<F(long)>& func2,
//                                  long n1, long n2) {
//
//    period = {(unsigned long)n1+1, (unsigned long)n2};
//    cout << "period " << period[0] << " " << period[1] << endl;
//    delta_size = 0;
//    A.resize(n1+1, n2);
//    blitz::TinyVector<unsigned long, 2> index;
//    for (long i = 0; i < n1; i++) {
//        for (long j = 0; j < n2; j++) {
//            index = i,j;
//            A(index) = func2(((j-func1(i)) % n2 + n2) % n2);
//        }
//        for (long j=0; j<n2; j++){
//            index = n1, j;
//            A(index) = 0;
//        }
//    }
//
//// FOR PRINTING THE ARRAY
//    for (long j=n2-1; j>=0; j--) {
//        for (long i=0; i<n1+1; i++) {
//            index = i,j;
//            cout << A(index) << " ";
//        }
//        cout << endl;
//    }
//}



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
    return (e1.lex_less(e2));
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

    if (n>15) n = 15;

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
    unsigned long n, m;
    m = idMatrix.size();
    if (m == 0) {
        n = 0;
    }
    else {
        n = idMatrix[0].size();
    }

    vector<F> lastRow(n+1, (F)0);
    lastRow[n] = (F)1;
    for (unsigned long i = 0; i < m ; i++){
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
void printPoly(vector< vector<F> >& idMatrix, vector< MExponent<m> >& idColumn) {
    unsigned long n = idMatrix.size()-1;
    for (unsigned long i = 0; i < idMatrix[0].size(); i++) {
        if(idMatrix[n][i] != (F)0) cout << idMatrix[n][i] << "x" << idColumn[i] << "\t";
    }
    cout << endl;
}

template<typename F, unsigned long m>
vector< MExponent<m> > getPoly(vector< vector<F> >& idMatrix, vector< MExponent<m> >& idColumn) {
    vector< MExponent<m> > poly;
    unsigned long n = idMatrix.size()-1;
    for (unsigned long i = 0; i < idMatrix[0].size(); i++) {
        if(idMatrix[n][i] != (F)0) poly.push_back(idColumn[i]);
    }
    return poly;
}

template <unsigned long m>
void printStars(vector< MExponent<m> >& v, vector<unsigned long>& period) {
    for(long i=2*period[1]+1; i>0; i--){
        for(long j=0; j<period[0]; j++) {
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
template<typename F, unsigned long m>
void MultiDimArray<F,m>::RST() {
    MExponent<m> bound;
    for (unsigned long i = 0; i < m; i++) {
        bound.editExp()[i] = 2 * period[i] - 1;
    }

    vector< MExponent<m> >  exponentsColumn(sizeOfDivisors(bound)),
                            idColumn;
    forward_list< MExponent<m> > exponentsRow;


    unsigned long index = 0;
    MExponent<m> exp;
    generateDivisors<m>(m, index, exp, bound, exponentsColumn);

    sort(exponentsColumn.begin(), exponentsColumn.end(), comp<m>);

    MExponent<m> e_period(period);
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
//            leadingMonomials.push_back(alpha);
//            printPoly(idMatrix, idColumn);
//            grobnerBasis.push_back(getPoly(idMatrix, idColumn));
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
