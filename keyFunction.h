#ifndef KEYFUNCTION_H
#define KEYFUNCTION_H
#include <iostream>
#include <cmath>
#include <map>
#include <tuple>
using namespace std;

// QR factorization
void QR(double ** a, int m, int n);

//Tensor times matrix using sparse tensor through all three modes.
double*** ttm_update(map<tuple<int,int,int>,double> spX, double ** U, double ** V, double ** W, int K[]);

// TTMcTC function described in paper
double** TTMcTC_update(map<tuple<int,int,int>,double> X, double*** G, double **U, double **V, double **W, int J[], int K[], int mode);


#endif //KEYFUNCTION_H
