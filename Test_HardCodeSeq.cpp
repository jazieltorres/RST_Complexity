/******************************************************************************************
*
*               HARD-CODED SEQUENCES
*
* This code helps to compute linear complexity of arrays constructed using composition
* method, specifying (hard-coded) every parameter for each sequence.
******************************************************************************************/

#include "MultiDimArray.cpp"
#include "Sequences.h"
//#include "NTL/ZZ_p.h"
#include "NTL/GF2.h"
#include <chrono>
#include <iostream>

using namespace std;

int main() {
/******************************************************
* Array with entries in the field F
*******************************************************/
    typedef NTL::GF2 F;

//    NTL::ZZ p(2);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;


/******************************************************
* Dimension of the array
*******************************************************/
    const unsigned int dim = 2;


/******************************************************
* Parameters
*******************************************************/
    int p_legendre = 3;
    int p_shift = 11;
    int root = 7;

//    while(!isRoot(root, p_shift)) root++;

    int coeff = 4; // Coefficient for ExpQuadraticSeq


/******************************************************
*                  COLUMN SEQUENCES
*
* Uncomment desired sequence
*******************************************************/
    LegendreSeq column_seq(p_legendre);


/******************************************************
*                  SHIFT SEQUENCES
*
* Uncomment desired sequence
*******************************************************/
//    PolynomialSeq shift_seq(p_shift, root, 2, 3);
    ExpQuadraticSeq shift_seq(p_shift, root, coeff);
//    ExponentialSeq shift_seq(p_shift, root);
//    vector<int> shift_seq = {1,6,3,5,4,2};


/******************************************************
* Printing the sequences
*******************************************************/
    cout << "Column sequence:\n\t";
    for (int i = 0; i < p_legendre; i++) cout << column_seq(i) << " ";
    cout << endl;

    cout << "Shift sequence:\n\t";
    for (int i = 0; i < p_shift - 1; i++) cout << shift_seq(i) << " ";
    cout << endl << endl;


/******************************************************
* Creating the array
*******************************************************/
    MultiDimArray<F, dim> A(shift_seq, column_seq, p_shift - 1, p_legendre);


/******************************************************
* Computing complexity
*******************************************************/
    cout << "Computing the complexity" << endl;
    auto start = chrono::high_resolution_clock::now();
    A.RST(3);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
    cout << "Computation time:\t" << duration.count() << " ms" << endl << endl;

    cout << "Complexity:\t" << A.complexity() << endl;
    cout << "Normalized:\t" << A.normalized_complexity() << endl << endl;


/******************************************************
* Checking the conjecture:
*******************************************************/
    unsigned int d = A.complexity();
    int n1 = A.period_vector()[0];
    int criteria = p_legendre % 8;
    int expected;
    bool satisfied = true;

    if (criteria == 1 && ((p_legendre - 1) / 2 * n1 != d)) {
        satisfied = false;
        expected = (p_legendre - 1) / 2 * n1;
    }
    if (criteria == 3 && (n1 * (p_legendre - 1) + 1 != d)) {
        satisfied = false;
        expected = n1 * (p_legendre - 1) + 1;
    }
    if (criteria == 5 && (n1 * (p_legendre - 1) != d)) {
        satisfied = false;
        expected = n1 * (p_legendre - 1);
    }
    if (criteria == 7 && ((p_legendre - 1) / 2 * n1 + 1 != d)) {
        satisfied = false;
        expected = (p_legendre - 1) / 2 * n1 + 1;
    }

    if (!satisfied) {
        cout << "CONJECTURE FAILED." << endl;
        cout << "Expected complexity:\t" << expected << endl;
        cout << "Expected normalized:\t" << static_cast<double>(expected)/A.period_size() << endl;
    }
    else {
        expected = A.complexity();
        cout << "Conjecture satisfied" << endl << endl;
    }
    cout << "Grobner basis:" << endl;
    A.print_basis();

/******************************************************
* This line is for the google sheet:
*******************************************************/
//    cout << p_legendre << ',' << p_legendre << ',' << root << ',' << p_legendre << ',' << p_legendre-1 << ','
//         << A.period_size() << ',' << A.complexity() << ',' << expected << endl;

    return 0;
}