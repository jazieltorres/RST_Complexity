/********************************************************************
 * WARNING: This file is used by the researchers to run tests, so it
 * does not perform data validation.
 *
 * Description:
 * A main file used for computing the linear complexity of arrays in
 * two dimensions produce by composition of Exponential Quadratic and
 * Legendre column (commensurate).
 * This program receives a command line argument: the prime.
 ********************************************************************/

#include "MultiDimArray_GF2.cpp"
#include "Sequences.h"
#include "NTL/GF2.h"
#include "NTL/ZZ.h"
#include <iostream>
#include <string>
#include <chrono>

using namespace std;

int main(int argc, char *argv[]) {
    typedef NTL::GF2 F;

/** Fixed dimension for this file **/
    const unsigned int dim = 2;


/************************************************
 *  Reading arguments
 ************************************************/

    long prime=0, low_bound, up_bound;
    if (argc == 1) {
        cerr << "Not argument given: prime number or two numbers for an interval." << endl
             << "Terminating program." << endl;
        return 0;
    }
    else if (argc == 2) {
        istringstream ss(argv[1]);
        if (!(ss >> prime)) {
            cerr << "Invalid argument: " << argv[1] << endl
                 << "Terminating program." << endl;
            return 0;
        }
        else if (!NTL::ProbPrime(prime)) {
            cerr << "Invalid argument: " << argv[1] << " must be prime" << endl
                     << "Terminating program." << endl;
            return 0;
        }
    }
    else {
        istringstream ss1(argv[1]);
        istringstream ss2(argv[2]);
        if (!(ss1 >> low_bound) or !(ss2 >> up_bound)) {
            cerr << "Invalid argument: " << argv[1] << " or " << argv[2] << endl
                 << "Terminating program." << endl;
            return 0;
        }
        else if (low_bound >= up_bound) {
            cerr << "Invalid arguments: " << argv[1] << " must be less than " << argv[2] << endl
                     << "Terminating program." << endl;
            return 0;
        }
    }

    if (prime != 0) {
        low_bound = prime;
        up_bound = prime;
    }

    prime = NTL::NextPrime(low_bound);

    while (prime <= up_bound) {
        LegendreSeq column_seq(prime);

        long root = 2;
        while (!isRoot(root, prime))
            root++;

        ExpQuadraticSeq shift_seq(prime, root, 1);

        MultiDimArray_GF2<dim> A(shift_seq, column_seq, prime-1, prime);

    /******************************************************
    * Computing complexity
    *******************************************************/
        cout << "Prime " << prime << " ------------------------------" << endl;
        auto start = chrono::high_resolution_clock::now();
        A.RST_simple();
        auto stop = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
        cout << "Computation time:\t" << duration.count() << " ms" << endl;

        cout << "Complexity:\t" << A.complexity() << endl;
        cout << "Normalized:\t" << A.normalized_complexity() << endl;


    /******************************************************
    * Checking the conjecture:
    *******************************************************/
        unsigned int d = A.complexity();
        int n1 = A.period_vector()[0];
        int criteria = prime % 8;
        int expected;
        bool satisfied = true;

        if (criteria == 1 && ((prime - 1) / 2 * n1 != d)) {
            satisfied = false;
            expected = (prime - 1) / 2 * n1;
        }
        if (criteria == 3 && (n1 * (prime - 1) + 1 != d)) {
            satisfied = false;
            expected = n1 * (prime - 1) + 1;
        }
        if (criteria == 5 && (n1 * (prime - 1) != d)) {
            satisfied = false;
            expected = n1 * (prime - 1);
        }
        if (criteria == 7 && ((prime - 1) / 2 * n1 + 1 != d)) {
            satisfied = false;
            expected = (prime - 1) / 2 * n1 + 1;
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

//        A.print_basis();
//        cout << endl << endl;

        cerr << prime << ',' << prime << ',' << root << ',' << 1 << ','
             << A.complexity() << ',' << expected << endl;
        prime = NTL::NextPrime(prime+1);
    }

    return 0;
}