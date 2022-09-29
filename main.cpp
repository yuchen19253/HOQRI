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
#include<random>
#include<ctime>
#include <time.h>

using namespace std;

int main () {
    bool full = 0; // use thin QR

    // J and K
//   int J[3] = {162541, 49994, 9083}; //10M_MovieLen: 162541, 49994, 9083 with 20503478 nonzeros
   int J[3] = {46952,46951,1592}; //facebook: 46952 x 46951 x 1592 with 738079 nonzeros
//   int J[3] = {108035, 107253, 52955}; //Delicious: 108035 x 107253 x 52955 with 437593 nonzeros
//   int J[3] = {2100, 18744, 12647}; //Last: 2100 x 18744 x 12647 with 186479 nonzeros
//   int J[3] = {610,49961,8215}; //1M_MovieLen: 610, 49961, 8215 with 84159 nonzeros
//   int J[3] = {5,5,5};

   int K[3] = {10,10,10};
//   int K[3] = {4,3,2};
    if(J[0]<K[0]||J[1]<K[1]||J[2]<K[2]){
        cout<<"Rank couldn't be large than original size!";
        return 0;
    }
    // read X from tensor file
    map<tuple<int,int,int>,double> mytensor;
    fstream tensorfile;
//    string filename = "/home/yuchen/Desktop/mytensor.txt";
//    string filename = "/home/yuchen/Desktop/data/Large_MovieLen.txt";
    string filename = "/home/yuchen/Desktop/data/facebook.txt"; // /Users/yc
//    string filename = "/home/yuchen/Desktop/MovieLen.txt";
//    string filename = "/home/yuchen/Desktop/Last.txt";
//    string filename = "/home/yuchen/Desktop/Delicious.txt";
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
    print_sptensor("Input tensor",mytensor);


    // initialize U, V and W
//    double U[J[0]][K[0]];
//    double V[J[1]][K[1]];
//    double W[J[2]][K[2]];
//    fill(U[0], U[0] + J[0] * K[0], 0);
//    fill(V[0], V[0] + J[1] * K[1], 0);
//    fill(W[0], W[0] + J[2] * K[2], 0);


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
    for(int i=0; i<J[0]; i++){
        for(int j=0; j<K[0]; j++) {
            U[i][j] = 0;
        }
    }
    for(int i=0; i<J[1]; i++){
        for(int j=0; j<K[1]; j++) {
            V[i][j] = 0;
        }
    }
    for(int i=0; i<J[2]; i++){
        for(int j=0; j<K[2]; j++) {
            W[i][j] = 0;
        }
    }
    int cnt = J[0]/K[0];
    for(int j=0; j<K[0]; j++){
        for(int i=0; i<cnt; i++) {
            U[j+i*K[0]][j]=1/ sqrt(cnt);
        }
    }
    cnt = J[1]/K[1];
    for(int j=0; j<K[1]; j++){
        for(int i=0; i<cnt; i++) {
            V[j+i*K[1]][j]=1/ sqrt(cnt);
        }
    }
    cnt = J[2]/K[2];
    for(int j=0; j<K[2]; j++){
        for(int i=0; i<cnt; i++) {
            W[j+i*K[2]][j]=1/ sqrt(cnt);
        }
    }
//    print_matrix("U", U, J[0], K[0]);
//    print_matrix("V", V, J[1], K[1]);
//    print_matrix("W", W, J[2], K[2]);
//
//    default_random_engine e(time(0));
//    uniform_real_distribution<double> u(0,1);
//    for(int j=0; j<J[0]; j++){
//        for(int k=0; k<K[0]; k++) {
//            U[j][k] = u(e);
//        }
//    }
//    for(int j=0; j<J[1]; j++){
//        for(int k=0; k<K[1]; k++) {
//            V[j][k] = u(e);
//        }
//    }
//    for(int j=0; j<J[2]; j++){
//        for(int k=0; k<K[2]; k++) {
//            W[j][k] = u(e);
//        }
//    }
//
//    print_matrix("U", U, J[0], K[0]);
//    print_matrix("V", V, J[1], K[1]);
//    print_matrix("W", W, J[2], K[2]);

    /* allocate memory for the matrices Rs */
//    double ** r1 = new double*[K[0]];
//    double ** r2 = new double*[K[1]];
//    double ** r3 = new double*[K[2]];
//    for(int index = 0; index < K[0]; index++) {
//        r1[index] = new double[K[0]];
//    }
//    for(int index = 0; index < K[1]; index++) {
//        r2[index] = new double[K[1]];
//    }
//    for(int index = 0; index < K[2]; index++) {
//        r3[index] = new double[K[2]];
//    }
//    qr(U, J[0], K[0]);
//    qr(V, J[1], K[1]);
//    qr(W, J[2], K[2]);

//    print_matrix("r1", r1, K[0],K[0]);
//    print_matrix("r2", r2, K[1],K[1]);
//    print_matrix("r3", r3, K[2],K[2]);

    // update iteration
    int iteration = 50;
    double *** G = ttm_update(mytensor, U, V, W, K);
    double newLoss = norm_tensor(G,K), oldLoss(0), lossChange(0);

    cout<<"Initial norm of G: " << newLoss <<endl;
    cout<<"Itr\t\tTime\t\tLoss Change\t\tnorm(core)"<<endl;

    clock_t start = clock();
    for(int itr=1; itr <= iteration; itr++){
//        cout<<"----------------------------------------- Iteration "<<itr+1<<" -----------------------------------------"<<endl;
//        print_sptensor("X"+ to_string(itr), mytensor);

        // Calculate G
//        cout<< "calculate G(mode1) " << (double)(clock() - start) / (double)CLOCKS_PER_SEC<<endl;
//        map<tuple<int,int,int>,double> G1=spttm_mode(mytensor,U,J,K[0],1);
//        cout<< "calculate G(mode2) "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
//        map<tuple<int,int,int>,double> G2=spttm_mode(G1,V,J,K[1],2);
//        cout<< "calculate G(mode3) " << (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
//        map<tuple<int,int,int>,double> G3=spttm_mode(G2,W,J,K[2],3);
//        print_sptensor("sparse G", G3);
//
//        cout<< "sparse to dense G "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
//        double *** G = sp2dense(G3,K);
//        print_densetensor("G",G,K);
//
//        map<tuple<int,int,int>,double> ttmG = spttm(mytensor, U, V, W, K);
//        cout<< "calculate sparse G "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
////        print_sptensor("ttm G", ttmG);
//        double *** G2 = sp2dense(ttmG,K);
//        print_densetensor("G",G2,K);

//        double *** G = ttm_update(mytensor, U, V, W, K);
//        print_densetensor("G",G,K);

//        cout<< "calculate dense G "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;

//        cout<< "calculate norm(G) "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;


        // Calculate A1
        double **A1 = TTMcTC_update(mytensor,G,U,V,W,J,K,1);
//        cout<< "calculate A1 "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
        // Calculate new U via QR
        for(int j=0; j<J[0]; j++){
            for(int k=0; k<K[0]; k++) {
                U[j][k] = A1[j][k];
            }
        }

        qr(U, J[0], K[0]);
//        print_matrix("newU", U, J[0],K[0]);
//        cout<< "calculate U "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;

        // Calculate A2
        double **A2 = TTMcTC_update(mytensor,G,U,V,W,J,K,2);
//        cout<< "calculate A2 "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
        // Calculate new V via QR
        for(int j=0; j<J[1]; j++){
            for(int k=0; k<K[1]; k++) {
                V[j][k] = A2[j][k];
            }
        }
        /* execute gramSchmidt to compute QR factorization */
        qr(V, J[1], K[1]);
//        cout<< "calculate V "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
        /* print the matrix Q resulting from gramSchmidt */
//        print_matrix("newV", V, J[1],K[1]);


        // Calculate A3
        double **A3 = TTMcTC_update(mytensor,G,U,V,W,J,K,3);
//        cout<< "calculate A3 "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
        // Calculate new W via QR
        for(int j=0; j<J[2]; j++){
            for(int k=0; k<K[2]; k++) {
                W[j][k] = A3[j][k];
            }
        }
        /* execute gramSchmidt to compute QR factorization */
        qr(W, J[2], K[2]);
//        cout<< "calculate W "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
        /* print the matrix Q resulting from gramSchmidt */
//        print_matrix("newW", W, J[2],K[2]);

        oldLoss = newLoss;
        G = ttm_update(mytensor, U, V, W, K);
        newLoss = norm_tensor(G,K);
        lossChange = abs(newLoss - oldLoss);

        double norm_G = norm_tensor(G, K);
        std::cout << itr<< "\t\t" << (double)(clock() - start) / (double)CLOCKS_PER_SEC<< "\t\t" << lossChange<< "\t\t" << norm_G << endl;


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



//    /* print the matrix R resulting from gramSchmidt */
//    print_matrix("R", r, K[0],K[0]);
//    std::cout << "R = " << std::endl;
//    for(int i = 0; i < K[0]; i++) {
//        for(int j = 0; j < K[0]; j++) {
//            printf("%9.6lg ", r[j][i]);
//        }
//        std::cout << std::endl;
//    }
//    std::cout << std::endl;

//
//    /* print numerical evidence that columns of Q are orthonormal */
//    printf("Numerical verification that {q_1, ..., q_%i} is an "
//           "orthonormal set:\n", K[0]);
//    for(int i = 0; i < K[0]; i++) {
//        for(int j = i; j < K[0]; j++) {
//            x = dot_product(a[i], a[j], J[0]);
//            printf("q_%i * q_%i = %lg\n", i + 1, j + 1, x);
//        }
//    }


//    double ***X;
//    X = new double ** [J[0]];
//    for (int i = 0;i < J[0];i++)
//    {
//        X[i] = new double*[J[1]];
//        for (int j = 0;j < J[1];j++)
//        {
//            X[i][j] = new double [J[2]];
//        }
//    }
//
//    for(int i=0; i<J[0]; i++){
//        for(int j=0; j<J[1]; j++){
//            for(int k=0; k<J[2]; k++) {
//                X[i][j][k] = u(e);
//            }
//        }
//    }
//    std::cout << "X = " << std::endl;
//    for(int i=0; i<J[2]; i++){
//        for(int j=0; j<J[0]; j++){
//            for(int k=0; k<J[1]; k++) {
//                printf("%9.6lg ", X[j][k][i]);
//            }
//            cout<<endl;
//        }
//        cout<< endl;
//    }
//    std::cout << std::endl;

//    for(int i = 0; i < K[0]; i++) {
//        delete[] r1[i];
//    }
//    for(int i = 0; i < K[1]; i++) {
//        delete[] r2[i];
//    }
//    for(int i = 0; i < K[2]; i++) {
//        delete[] r3[i];
//    }
//    delete[] r1;
//    delete[] r2;
//    delete[] r3;

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

//    for(int i=0; i<J[0]; i++){
//        for(int j=0; j<J[1]; j++){
//            delete[] X[i][j];
//        }
//    }
//    for(int i=0; i<J[0]; i++){
//        delete[] X[i];
//    }
//    delete[] X;

    return 0;       // exit main
}
