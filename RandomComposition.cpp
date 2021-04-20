#include "MultiDimArray_GF2.cpp"
#include "NTL/GF2.h"
#include "NTL/ZZ_p.h"
#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <vector>

// Program receives command line argument: size of the array (square array)

using namespace std;

//int main(int argc, char *argv[]) {
//    int numTest = 5;
//    int len;
//
//    NTL::GF2 zero(0);
//    NTL::GF2 one(1);
//    NTL::GF2 constant;
//
//    if (argc < 2) {
//        cerr << "Not enough arguments." << endl;
//        return 0;
//    }
//
//    istringstream ss(argv[1]);
//    if (!(ss >> len)) {
//        cerr << "Invalid argument " << argv[1] << endl;
//        return 0;
//    }
//
//    if (argc == 3) {
//        istringstream ss1(argv[2]);
//        if(!(ss1 >> constant)) {
//            cerr << "Not viable conversion from" << argv[2] << " to GF2." << endl;
//            return 0;
//        }
//    }
//    else {
//        constant = zero;
//    }
//
//
//    blitz::TinyVector<int, 2> v(len, len);
//    srand (time(NULL));
//    for (int k=0; k<numTest; k++) {
//        vector<int> shift_seq(len);
//        for (int i=0; i<len; i++) {
//            int t = rand() % len+1;
//            if (t == len)
//                shift_seq[i] = -1;
//            else
//                shift_seq[i] = t;
//        }
//
////        cout << "Shift sequence:\t\t";
////        for (auto &x : shift_seq) {
////            cout << x << " ";
////        } cout << endl;
//
//        for (int i=0; i<numTest; i++) {
//            vector<GF2> column_seq(len);
//            for (int j=0; j<len; j++) {
//                bool value = rand() % 2;
//                if (value)
//                    column_seq[j] = one;
//                else
//                    column_seq[j] = zero;
//            }
//
////            cout << "Column sequence:\t";
////            for (auto &y : column_seq) {
////                cout << y << " ";
////            } cout << endl;
//
//            MultiDimArray_GF2<1> seq(column_seq);
//            seq.RST_simple();
//
//            MultiDimArray_GF2<2> A(shift_seq, column_seq, constant);
//            A.RST_simple();
////            cout << "Column complexity:     " << seq.complexity() << endl;
////            cout << "Array complexity:      " << A.complexity() << endl;
////            cout << "Normalized difference: " << seq.normalized_complexity() - A.normalized_complexity() << endl;
////            cout << "Array:" << endl;
////            A.print_array();
////            cout << endl;
//
//
//            cout << seq.normalized_complexity() - A.normalized_complexity() << endl;
//        }
//    }
//
//
//
//    return 0;
//}


int main(int argc, char *argv[]) {
    int numTest = 20;
    int len;

    NTL::GF2 zero(0);
    NTL::GF2 one(1);
    NTL::GF2 constant;

    if (argc < 2) {
        cerr << "Not enough arguments." << endl;
        return 0;
    }

    istringstream ss(argv[1]);
    if (!(ss >> len)) {
        cerr << "Invalid argument " << argv[1] << endl;
        return 0;
    }

    if (argc == 3) {
        istringstream ss1(argv[2]);
        if(!(ss1 >> constant)) {
            cerr << "Not viable conversion from" << argv[2] << " to GF2." << endl;
            return 0;
        }
    }
    else {
        constant = zero;
    }


    blitz::TinyVector<int, 2> v(len, len);
    srand (time(NULL));
    for (int k=0; k<numTest; k++) {
        vector<int> shift_seq(len);
        for (int i=0; i<len; i++) {
            int t = rand() % len+1;
            if (t == len)
                shift_seq[i] = -1;
            else
                shift_seq[i] = t;
        }

        vector<GF2> column_seq(len);
        for (int j=0; j<len; j++) {
            bool value = rand() % 2;
            if (value)
                column_seq[j] = one;
            else
                column_seq[j] = zero;
        }

        MultiDimArray_GF2<1> seq(column_seq);
        seq.RST_simple();

        MultiDimArray_GF2<2> A(shift_seq, column_seq, constant);
        A.RST_simple();

        cout << len << "," << seq.normalized_complexity() - A.normalized_complexity() << endl;
    }
    return 0;
}