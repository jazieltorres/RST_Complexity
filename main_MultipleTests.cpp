/********************************************************************
 * WARNING: This file is used by the researchers to run tests, so it
 * does not perform data validation.
 *
 * Description:
 * A main file used for computing the linear complexity of
 * multidimensional periodic arrays. This main must receive a text
 * file with the necessary information, similar to the file main.cpp
 *
 *      The file "main.cpp" receives a text file with the data of just
 *      one array.
 *
 *      This "main_MultipleTests.cpp" receives a file with similar
 *      format but with the information of multiple arrays.
 *
 *  Note #1: The base field must be hard-coded. It is not identified
 *           using the input file.
 *  Note #2: The number of tests in the file is not required. This
 *           program iterates using a while loop.
 *
 * The format is as follows (assuming array with dimension n):
 *      Fist line:
 *          n positive integers indicating the length of the array in
 *          each dimension, each integer separated by a blank space.
 *      Then, for each array, a block with the following lines:
 *          One line with the test number.
 *          Two more lines of information:
 *              The first three lines of each block are not used for
 *              computing the complexity, but these lines should be
 *              used to identify each test. An example:
 *                      Test #1
 *                      Primitive polynomial: x^2 + x + 1
 *                      Coefficient: 3
 *          Following lines of the block:
 *              Each line containing n numbers indicating a coordinate,
 *              and the entry value at that coordinate, all separated
 *              by a blank space.
 ********************************************************************/

//#include "MultiDimArray.cpp"
#include "MultiDimArray_GF2.cpp"
#include "NTL/GF2.h"
//#include "NTL/ZZ_p.h"
#include <iostream>
#include <string>
#include <chrono>

using namespace std;

int main() {
/** Change desired field **/
//    NTL::ZZ p(3);
//    NTL::ZZ_p::init(p);
//    typedef NTL::ZZ_p F;
    typedef NTL::GF2 F;

/** Change desired dimension **/
    const unsigned int dim = 3;

/************************************************
 *  Processing first line:
 *      Length of array in each dimension
 ************************************************/
    string line;
    blitz::TinyVector<int, dim> v;
    int size = 1;
    for (int i = 0; i < dim; i++) {
        cin >> v[i];
        size = size * v[i];
    }

//    Removing character in buffer after cin
    getline(cin, line);

    int ctr = 0;

/************************************************
 *  Processing each test:
 *      coordinates and value of an array
 ************************************************/
    while (getline(cin, line)) { // Skipping line with counter number
        ctr++;
//        cout << "Test " << ctr << " -----------------------------" << endl;
        // Printing following two lines
        getline(cin, line); // First line of information
//        cout << line << endl;
        getline(cin, line); // Second line of information
//        cout << line << endl;

        MultiDimArray_GF2<dim> A(v);
        for (int i = 0; i < size; i++) {
            blitz::TinyVector<int, dim> position;
            F value;

            getline(cin, line);
            stringstream ss(line);
            for (int n = 0; n < dim; n++)
                ss >> position[n];
            ss >> value;
            A.set_at(position, value);
        }

        auto start = chrono::high_resolution_clock::now();
        A.RST_simple();
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "Duration:  " << duration.count() << " ms" << endl;
        cout << "Complexity:  " << A.complexity() << endl;
        cout << "Period size: " << A.period_size() << endl;
        cout << "Normalized:  " << A.normalized_complexity() << endl;
        cout << "Groebner Basis:" << endl;
        A.print_basis();
//        cout << endl;
//        A.print_array();
//        cout << endl << endl;


/**  Print just the complexity **/
//        cout << A.complexity() <<  endl;
    }

    return 0;
}