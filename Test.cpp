#include <NTL/ZZ_p.h>
#include <vector>

using namespace std;


void printM (const vector<vector<NTL::ZZ_p>>& M) {
    long m = M.size();
    long n = M[0].size();

    for (long i = 0; i<m; i++){
        for (long j=0; j<n; j++){
            cout << M[i][j] << " ";
        }
        cout << endl;
    }
}

void instertBefore(vector<vector<NTL::ZZ_p>>& M, long& index, long& m) {
    for (long i = 0; i < (m - index); i++) {
        M[m-i].swap(M[m-i-1]) ;
    }
}

void reduce(vector<vector<NTL::ZZ_p>>& M) {
    long m = M.size() - 1;
    long n = M[0].size()-1;
    const NTL::ZZ_p zero((long)0);
    vector<NTL::ZZ_p> lambda(m, zero);
    long j = 0 ;

    for (long i = 0; i < m; i++) {
        while (M[i][j] == 0 && j <= n) {
            for (long k = 0; k < i; k++) {
                M[m][j] += lambda[k] * M[k][j];
            }
            if (M[m][j] != 0) {
                for (long l = j+1; l <= n; l++) {
                    for (long k = 0; k < i; k++) {
                        M[m][l] += lambda[k] * M[k][l];
                    }
                }
                instertBefore(M, i, m);
                return;
            }
            j++;
        }
        lambda[i] = M[m][j];
        for (long k = 0; k < i; k++) {
            lambda[i] += lambda[k] * M[k][j];
        }
        lambda[i] = -lambda[i]/M[i][j];
        M[m][j] = 0;
        j++;
    }

    for (long l = j; l <= n; l++) {
        for (long k = 0; k < m; k++) {
            M[m][l] += lambda[k] * M[k][l];
        }
    }
    cout << "Lambda \t";
    for (long k = 0; k<lambda.size(); k++){
        cout << lambda[k] << " ";
    }
    cout << endl;
    return;
}



int main()
{
   NTL::ZZ p = (NTL::ZZ) 7;
   NTL::ZZ_p::init(p);

   vector<NTL::ZZ_p> M1({(NTL::ZZ_p)2, (NTL::ZZ_p)0, (NTL::ZZ_p)3});
   vector<NTL::ZZ_p> M2({(NTL::ZZ_p)0, (NTL::ZZ_p)0, (NTL::ZZ_p)1});
   vector<NTL::ZZ_p> M3({(NTL::ZZ_p)3, (NTL::ZZ_p)1, (NTL::ZZ_p)3});

   vector<vector<NTL::ZZ_p>> M({M1, M2, M3});

   reduce(M);
   printM(M);

   cout << "Corre..." << endl;

   return 0;
}



// g++ -g -O2 -std=c++11 -pthread -march=native Test.cpp -o Test -lntl -lgmp -lm
