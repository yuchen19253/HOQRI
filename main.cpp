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
    int J[3] = {46952,46951,1592}; //facebook: 46952 x 46951 x 1592 with 738079 nonzeros
//    int J[3] = {610,49961,8215}; //MovieLen
//    int J[3] = {100,100,100};
    int K[3] = {10,10,10};
    // read X from tensor file
    map<tuple<int,int,int>,double> mytensor;
    fstream tensorfile;
    string filename = "/Users/yc/Desktop/facebook.txt"; // /home/yuchen
//    string filename = "/Users/yc/Desktop/MovieLen.txt";
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
    /* allocate memory for the matrices Rs */
    double ** r1 = new double*[K[0]];
    double ** r2 = new double*[K[1]];
    double ** r3 = new double*[K[2]];
    for(int index = 0; index < K[0]; index++) {
        r1[index] = new double[K[0]];
    }
    for(int index = 0; index < K[1]; index++) {
        r2[index] = new double[K[1]];
    }
    for(int index = 0; index < K[2]; index++) {
        r3[index] = new double[K[2]];
    }
//    print_matrix("U", U, J[0], K[0]);
//    print_matrix("V", V, J[1], K[1]);
//    print_matrix("W", W, J[2], K[2]);

    qr(U, J[0], K[0]);
    qr(V, J[1], K[1]);
    qr(W, J[2], K[2]);

//    print_matrix("U", U, J[0], K[0]);
//    print_matrix("V", V, J[1], K[1]);
//    print_matrix("W", W, J[2], K[2]);
//    print_matrix("r1", r1, K[0],K[0]);
//    print_matrix("r2", r2, K[1],K[1]);
//    print_matrix("r3", r3, K[2],K[2]);

    // update iteration
    int iteration = 1;
    double newLoss(0), oldLoss(0), lossChange(0);
    cout<<"Itr\t\tTime\t\tLoss Change\t\tnorm(core)"<<endl;
    for(int itr=0; itr <= iteration; itr++){
        clock_t start = clock();
        oldLoss = newLoss;
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

        map<tuple<int,int,int>,double> ttmG = spttm(mytensor, U, V, W, K);
        cout<< "calculate sparse G "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
//        print_sptensor("ttm G", ttmG);
        double *** G = sp2dense(ttmG,K);
        cout<< "calculate dense G "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
        newLoss = norm_tensor(G,K);
        cout<< "calculate norm(G) "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
        lossChange = abs(newLoss - oldLoss);



        // Calculate A1
        double **A1 = TTMcTC(mytensor,G,U,V,W,J,K,1);
        // Calculate new U via QR
        for(int j=0; j<J[0]; j++){
            for(int k=0; k<K[0]; k++) {
                U[j][k] = A1[j][k];
            }
        }
        cout<< "calculate A1 "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;

        qr(U, J[0], K[0]);
//        print_matrix("newU", U, J[0],K[0]);
        cout<< "calculate U "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;

        // Calculate A2
        double **A2 = TTMcTC(mytensor,G,U,V,W,J,K,2);
        // Calculate new V via QR
        for(int j=0; j<J[1]; j++){
            for(int k=0; k<K[1]; k++) {
                V[j][k] = A2[j][k];
            }
        }
        cout<< "calculate A2 "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
        /* execute gramSchmidt to compute QR factorization */
        qr(V, J[1], K[1]);
        cout<< "calculate V "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
        /* print the matrix Q resulting from gramSchmidt */
//        print_matrix("newV", V, J[1],K[1]);


        // Calculate A3
        double **A3 = TTMcTC(mytensor,G,U,V,W,J,K,3);
        // Calculate new W via QR
        for(int j=0; j<J[2]; j++){
            for(int k=0; k<K[2]; k++) {
                W[j][k] = A3[j][k];
            }
        }
        cout<< "calculate A3 "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
        /* execute gramSchmidt to compute QR factorization */
        qr(W, J[2], K[2]);
        cout<< "calculate W "<< (double)(clock() - start) / (double)CLOCKS_PER_SEC << endl;
        /* print the matrix Q resulting from gramSchmidt */
//        print_matrix("newW", W, J[2],K[2]);



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

        double norm_G = norm_tensor(G, K);
        std::cout << itr<< "\t\t" << (double)(clock() - start) / (double)CLOCKS_PER_SEC<< "\t\t" << lossChange<< "\t\t" << norm_G << endl;

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

    for(int i = 0; i < K[0]; i++) {
        delete[] r1[i];
    }
    for(int i = 0; i < K[1]; i++) {
        delete[] r2[i];
    }
    for(int i = 0; i < K[2]; i++) {
        delete[] r3[i];
    }
    delete[] r1;
    delete[] r2;
    delete[] r3;

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