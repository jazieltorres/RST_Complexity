/********************************************************************
 * WARNING: This file is used by the researchers to run tests, so it
 * does not perform data validation.
 *
 * Description:
 * A main file used for computing the linear complexity of arrays in
 * two dimensions produce by composition. The shift sequence is
 * passed to the program in a text file.
 * This program receives a command line argument: the constant to
 * be replaced in the column with infinity shift.
 *
 *  Note #1: The base field must be hard-coded. It is not identified
 *           using the input file.
 *  Note #2: The number of tests in the file is not required. This
 *           program iterates using a while loop.
 *
 * The format is as follows:
 *      Fist line:
 *          A positive integer indicating the length of the shift
 *          sequence
 *      Following lines:
 *          A list of positive integers, separated by spaces,
 *          indicating the shift. The value -1 represents infinity.
 ********************************************************************/

#include "MultiDimArray.cpp"
#include "Sequences.h"
#include "NTL/GF2.h"
//#include "NTL/ZZ_p.h"
#include <iostream>
#include <string>
#include <chrono>

using namespace std;

int main(int argc, char *argv[]) {
/** Change desired field **/
//    NTL::ZZ p(3);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;

    typedef NTL::GF2 F;

/** Fixed dimension for this file **/
    const unsigned int dim = 2;

/************************************************
 *  Reading argument, if passed:
 *      constant for infinity shift column
 ************************************************/
    F constant = F(0);
    if (argc > 1) {
        istringstream ss(argv[1]);
        int i;
        if (ss >> i) {
            constant = F(i);
        }
        else {
            cerr << "Invalid argument: " << argv[1] << endl
                 << "Terminating program" << endl;
            return 0;
        }
    }

/************************************************
 *  Reading first line:
 *      Length of shift vector
 ************************************************/
    int horizontal_period, vertical_period;
    cin >> horizontal_period;
    cin >> vertical_period;
    cin.ignore();

/************************************************
 *  Processing each test:
 *      coordinates and value of an array
 ************************************************/
    int ctr = 0;
    string line;
    LegendreSeq column_seq(vertical_period);
    while (getline(cin, line)) {
        ctr++;
//        cout << "Test " << ctr << " -----------------------------" << endl;
//        cout << "Prime " << vertical_period << endl;
        vector<int> shift_seq(horizontal_period);
        istringstream ss(line);
        for (int i=0; i<horizontal_period; i++) {
            ss >> shift_seq[i];
        }

        MultiDimArray<F, dim> A(shift_seq, constant, column_seq, vertical_period);
        auto start = chrono::high_resolution_clock::now();
        A.RST_simple();
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
//        cout << "Duration:  " << duration.count() << " ms" << endl;
//        cout << "Complexity:  " << A.complexity() << endl;
//        cout << "Period size: " << A.period_size() << endl;
//        cout << "Normalized:  " << A.normalized_complexity() << endl;
//        cout << "Groebner Basis:" << endl;
//        A.print_basis();
//        cout << endl;

        cout << A.complexity() << endl;
    }

    return 0;
}