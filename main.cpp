//
// Created by yuchen on 2/21/22.
//

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <tuple>
#include <map>

#include "linearalgebra.h"
#include "keyFunction.h"
#include "inout.h"
#include<random>
#include<ctime>
#include <time.h>

using namespace std;

int main (int argc, char *argv[]) {
    if(argc != 4 ){
        cout<<"Input error! Tensor Dimension, Rank and data file needed."<<endl;
        cout<<"Format: I1,I2,I3 r1,r2,r3 'DataFilePath'";
        return -1;
    }
    int* J = argv_split(argv[1]);
    int* K = argv_split(argv[2]);
    string filename = argv[3];

    // J and K
//   int J[3] = {100000,100000,100000}; // synthetic data;
//   int J[3] = {46952,46951,1592}; //Facebook: 46952 x 46951 x 1592 with 738079 nonzeros (2.1e-7)
//   int J[3] = {610,49961,8215}; //1M_MovieLen: 610, 49961, 8215 with 84159 nonzeros (3.36e-7)
//   int J[3] = {2100, 18744, 12647}; //Last: 2100 x 18744 x 12647 with 186479 nonzeros (3.7e-7)
//   int J[3] = {108035, 107253, 52955}; //Delicious: 108035 x 107253 x 52955 with 437593 nonzeros (7.1e-10)
//   int J[3] = {162541, 49994, 9083}; //10M_MovieLen: 162541, 49994, 9083 with 20503478 nonzeros (2.7e-7)

//    int K[3] = {10,10,10};
//    int K[3] = {20,20,10};
//    int K[3] = {15,25,15};
//    int K[3] = {20,20,20};
//    int K[3] = {20,20,10};
//    int K[3] = {20,20,10};

//    string filename = "/home/yuchen/Desktop/e5.txt";
//    string filename = "/home/yuchen/Desktop/data/Large_MovieLen.txt";
//    string filename = "/home/yuchen/Desktop/data/facebook.txt"; // /Users/yc
//    string filename = "/home/yuchen/Desktop/data/MovieLen.txt";
//    string filename = "/home/yuchen/Desktop/data/Last.txt";
//    string filename = "/home/yuchen/Desktop/data/Delicious.txt";


    if(J[0]<K[0]||J[1]<K[1]||J[2]<K[2]){
        cout<<"Rank couldn't be large than original size!";
        return 0;
    }
    // read X from tensor file
    map<tuple<int,int,int>,double> mytensor;
    fstream tensorfile;

    tensorfile.open(filename,ios::in);

    if(tensorfile.is_open()){
        cout<<"Open file: "<<filename<<endl;
        string line;
        while(getline(tensorfile,line)){
            istringstream ss(line);
            int index1;
            ss>>index1;
            int index2;
            ss>>index2;
            int index3;
            ss>>index3;
            tuple<int, int, int> index(index1,index2,index3);
            double value;
            ss>>value;
            mytensor[index] = value;
        }
    }
    else{
        cout<<"No file"<<endl;
        return 0;
    }
//    print_sptensor("Input tensor",mytensor);



    double **U = new double*[J[0]];
    for(int index = 0; index < J[0]; index++) {
        U[index] = new double[K[0]];
    }
    double **V = new double*[J[1]];
    for(int index = 0; index < J[1]; index++) {
        V[index] = new double[K[1]];
    }
    double **W = new double*[J[2]];
    for(int index = 0; index < J[2]; index++) {
        W[index] = new double[K[2]];
    }

    default_random_engine e(time(0));
    uniform_real_distribution<double> u(0,1);
    for(int j=0; j<J[0]; j++){
        for(int k=0; k<K[0]; k++) {
            U[j][k] = u(e);
        }
    }
    for(int j=0; j<J[1]; j++){
        for(int k=0; k<K[1]; k++) {
            V[j][k] = u(e);
        }
    }
    for(int j=0; j<J[2]; j++){
        for(int k=0; k<K[2]; k++) {
            W[j][k] = u(e);
        }
    }

    // update iteration
    int iteration = 10;
    double *** G = ttm_update(mytensor, U, V, W, K);
    double newLoss = norm_tensor(G,K), oldLoss(0), lossChange(0);

    cout<<"Initial norm of G: " << newLoss <<endl;
    cout<<"Itr\t\tTime\t\tLoss Change\t\tnorm(core)"<<endl;

//    ofstream outfile;
//    outfile.open("/home/yuchen/Desktop/hoqri_rank.txt", ios::app);
    clock_t start = clock();
    for(int itr=1; itr <= iteration; itr++){
        //        cout<< " " << (double)(clock() - start) / (double)CLOCKS_PER_SEC<<endl;

        // Calculate A1
        double **A1 = TTMcTC_update(mytensor,G,U,V,W,J,K,1);
        // Calculate A2
        double **A2 = TTMcTC_update(mytensor,G,U,V,W,J,K,2);
        // Calculate A3
        double **A3 = TTMcTC_update(mytensor,G,U,V,W,J,K,3);

        /* execute gramSchmidt to compute QR factorization */
        for(int j=0; j<J[0]; j++){
            for(int k=0; k<K[0]; k++) {
                U[j][k] = A1[j][k];
            }
        }
        for(int j=0; j<J[1]; j++){
            for(int k=0; k<K[1]; k++) {
                V[j][k] = A2[j][k];
            }
        }
        for(int j=0; j<J[2]; j++){
            for(int k=0; k<K[2]; k++) {
                W[j][k] = A3[j][k];
            }
        }
        QR(U,J[0], K[0]);
        QR(V,J[1], K[1]);
        QR(W,J[2], K[2]);

        /* execute gramSchmidt to compute QR factorization */
//        U = myQR(A1, J[0], K[0]);
//        V = myQR(A2, J[1], K[1]);
//        W = myQR(A3, J[2], K[2]);

        oldLoss = newLoss;
        G = ttm_update(mytensor, U, V, W, K);
        newLoss = norm_tensor(G,K);
        lossChange = abs(newLoss - oldLoss);

        double norm_G = norm_tensor(G, K);
        std::cout << itr<< "\t\t" << (double)(clock() - start) / (double)CLOCKS_PER_SEC<< "\t\t" << lossChange<< "\t\t" << norm_G << endl;
//        outfile << itr<< "\t\t" << (double)(clock() - start) / (double)CLOCKS_PER_SEC<< "\t\t" << lossChange<< "\t\t" << norm_G << endl;

        if(lossChange<1e-5)
            break;

        for(int i = 0; i < K[0]; i++) {
            delete[] A1[i];
        }
        for(int i = 0; i < K[1]; i++) {
            delete[] A2[i];
        }
        for(int i = 0; i < K[2]; i++) {
            delete[] A3[i];
        }
        delete[] A1;
        delete[] A2;
        delete[] A3;

    }
//    outfile.close();


    /* free memory */
    for(int i = 0; i < K[0]; i++) {
        delete[] U[i];
    }
    for(int i = 0; i < K[1]; i++) {
        delete[] V[i];
    }
    for(int i = 0; i < K[2]; i++) {
        delete[] W[i];
    }

    delete[] U;
    delete[] V;
    delete[] W;

    return 0;       // exit main
}
