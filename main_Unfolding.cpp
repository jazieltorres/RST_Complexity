/******************************************************
 *
 *  FILE TO STUDY COMPLEXITY OF THE SEQUENCE OBTAINED
 *  BY UNFOLDING THE ARRAY WITH RESPECT TO A MONOMIAL
 *  ORDERING.
 *
 ******************************************************/

#include "MultiDimArray.cpp"
#include "Sequences.h"
#include "NTL/ZZ_p.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main() {
    int prime = 2;

    NTL::ZZ p(prime);
    NTL::ZZ_p::init(p);
    typedef NTL::ZZ_p F;
    const unsigned int dim = 2;


    string line;
    blitz::TinyVector<int, dim> v, top_corner;
    int size = 1;
    int sizePlus = 1;
    for (int i = 0; i < dim; i++) {
        cin >> v[i];
        top_corner[i] = v[i] - 1;
        size = size * v[i];
        sizePlus = sizePlus * (v[i] + 1);
    }

    Monomial<dim> alpha, e_period(top_corner);

    vector< Monomial<dim> > coordinates(size);
    int x = 0;
    generateDivisors<dim>(dim, x, alpha, e_period, coordinates);
    sort(coordinates.begin(), coordinates.end(), ordering_1<dim>);


//    Removing character in buffer after cin
    getline(cin, line);

    int ctr = 0;

    while (getline(cin, line)) {
        ctr++;
        // Skipping following two lines
        getline(cin, line);
        getline(cin, line);

        blitz::Array<F, dim> A;
        A.resize(v);
        for (int i = 0; i < size; i++) {
            blitz::TinyVector<int, dim> position;
            F value;

            getline(cin, line);
            stringstream ss(line);
            for (int n = 0; n < dim; n++)
                ss >> position[n];
            ss >> value;
            A(position) = value;
        }

        blitz::TinyVector<int,1> index = size;
        blitz::TinyVector<int,dim> at;
        F value;
        MultiDimArray<F,1> Unfolded_array(index);
        for (int i = 0; i < size; i++) {
            index[0] = i;
            at = coordinates[i].exponent();
            value = A(at);
            Unfolded_array.set_at(index, value);
        }

        Unfolded_array.RST();
        cout << "Test " << ctr << " -----------------------------" << endl;
        cout << "Complexity:\t\t" << Unfolded_array.complexity() << endl;
        cout << "Period size:\t" << Unfolded_array.period_size() << endl;
        cout << "Normalized:\t\t" << Unfolded_array.normalized_complexity() << endl;
//        cout << "Grobner Basis:" << endl;
        Unfolded_array.print_basis();
        cout << endl << endl;

    }
    return 0;
}