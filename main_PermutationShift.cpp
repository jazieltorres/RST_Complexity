#include "MultiDimArray_GF2.cpp"
#include "Sequences.h"
#include "NTL/GF2.h"
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char *argv[]) {
    typedef NTL::GF2 F;
    const unsigned int dim = 2;

    int p_legendre = 13;
    int shift_len = p_legendre-1;

    LegendreSeq column_seq(p_legendre);
    vector<int> shift(shift_len);
    for (int i = 0; i < shift_len; i++) {
        shift[i] = i;
    }

// Fixing the second entry in the permutation as the integer n, entered as command-line argument
    int n;
    if (argc == 2) {
        istringstream ss(argv[1]);
        if (!(ss >> n)) {
            cerr << "Invalid argument: " << argv[1] << endl
                 << "Terminating program." << endl;
            return 0;
        }
        // Putting n in the second position, (by swap)
        for (int i = n; i>2; i--) {
            shift[i] = shift[i-1];
        }
        shift[1] = n;
        shift[2] = 1;
    }

    cout << "Original shift:" << endl;
    for (auto x : shift) cout << x << " ";
    cout << endl << endl;


    bool found = false;


    do {
//        for (auto x : shift) cout << x << " ";
//        cout << endl;

        MultiDimArray_GF2<dim> A(shift, F(0), column_seq, p_legendre);
        A.RST_simple();

        unsigned int d = A.complexity();
        int n1 = shift_len;
        int criteria = p_legendre % 8;
        bool satisfied = true;

        if (criteria == 1 && ((p_legendre - 1) / 2 * n1 != d)) {
            satisfied = false;
        }
        if (criteria == 3 && (n1 * (p_legendre - 1) + 1 != d)) {
            satisfied = false;
        }
        if (criteria == 5 && (n1 * (p_legendre - 1) != d)) {
            satisfied = false;
        }
        if (criteria == 7 && ((p_legendre - 1) / 2 * n1 + 1 != d)) {
            satisfied = false;
        }

        if (!satisfied) {
            for (auto x : shift) cout << x << " ";
            cout << endl;
            found = true;
        }
    } while (next_permutation(shift.begin()+2, shift.end()) && !found);
}