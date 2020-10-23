#include "MultiDimArray.cpp"
#include "Functions.h"
#include "NTL/ZZ_p.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>


/*******************************************************************
 * This is the main file for complexity analysis of multidimensional
 * periodic arrays. This main receives a txt file with the data.
 * The data in the txt file must be formatted as follow
 * (Assuming array with dimension n):
 *      Fist line:
 *          A prime number (the array will have entries in ZZ_p)
 *      Second line:
 *          n numbers indicating the length of the array in
 *              each dimension, separated by blank space.
 *      Following lines:
 *          n numbers indicating a coordinate, and the entry value
 *          at that coordinate, all separated by a blank space.
 *******************************************************************/


using namespace std;

int main() {
    const unsigned int dim = 1;
    int prime;
    string line;
    istringstream ss;
    blitz::TinyVector<int, dim> length(dim);
    int size = 1;

/************************************************
 *  Processing first line:
 *      A prime number
 ************************************************/
    getline(cin, line);
    ss.str(line);
    if (ss >> prime) {
        if (prime != 2) {
            if (!isPrime(prime)) {
                cerr << "ERROR: First number is not prime. \n\nProcess finished." << endl;
                return 0;
            }
        }
    }
    else {
        cerr << "ERROR: First line is empty. \n\nProcess finished." << endl;
        return 0;
    }
    if (ss >> prime) {
        cerr << "ERROR: More than one number in first line. \n\nProcess finished." << endl;
        return 0;
    }

    NTL::ZZ p(prime);
    NTL::ZZ_p::init(p);
    typedef NTL::ZZ_p F;



/************************************************
 *  Processing second line:
 *      Length of array in each dimension
 ************************************************/
    getline(cin, line);
    ss.clear();
    ss.str(line);

    for (int i = 0; i < dim; i++) {
        if (ss >> length[i])
            size = size * length[i];
        else {
            cerr << "ERROR: Second line with less parameters than selected dimension (" << dim << ")." <<
                 "\n\nProcess finished." << endl;
            return 0;
        }
    }
    if (ss >> length[0]) {
        cerr << "ERROR: Second line with more parameters than selected dimension (" << dim << ")." <<
             "\n\nProcess finished." << endl;
        return 0;
    }



/************************************************
 *  Processing remaining lines:
 *      coordinates and value
 ************************************************/
    MultiDimArray<F, dim> Array(length);
    for (int i = 0; i < size; i++) {
        if (getline(cin, line)) {
            ss.clear();
            ss.str(line);
            blitz::TinyVector<int, dim> coordinate;
            F value;
            for (int n = 0; n < dim; n++) {
                if (ss >> coordinate[n]) {
                    if (coordinate[n] > length[n]) {
                        cerr << "ERROR: Coordinate out of range. \n\nProcess finished." << endl;
                        return 0;
                    }
                } else {
                    cerr << "ERROR: Less than " << dim << " coordinates. \n\nProcess finished." << endl;
                    return 0;
                }
            }
            if (ss >> value) {
                Array.set_at(coordinate, value);
            } else {
                cerr << "ERROR: Missing value at coordinate " << coordinate << "." <<
                     "\n\nProcess finished." << endl;
                return 0;
            }
            if (ss >> value) {
                cerr << "ERROR: More parameters than expected: " << dim << "coordinates + 1 value." <<
                     "\n\nProcess finished." << endl;
                return 0;
            }
        } else {
            cerr << "ERROR: Entries missing values: " << i <<
                 " lines of coordinate-value instead of " << size << ". \n\nProcess finished." << endl;
            return 0;
        }
    }
    Array.RST();
    cout << "Complexity: " << Array.complexity() << endl;
    cout << "Normalized Complexity" << Array.normalized_complexity() << endl;
    return 0;
}



