#include "linearalgebra.h"

double norm (double * x, int length) {
    double sum = 0;
    for(int i=0; i < length; i++) {
        sum += x[i] * x[i];
    }

    return sqrt(sum);
}

double norm_tensor(double *** X, int J[]){
    double sum=0;
    for(int i=0; i<J[0]; i++){
        for(int j=0; j<J[1]; j++){
            for(int k=0; k<J[2]; k++) {
                sum += X[i][j][k] * X[i][j][k];
            }
        }
    }
    return sqrt(sum);
}


void scalar_div(double * x, double r, int length, double * y) {
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i] = x[i]/r;
    }
    for(; i < length; i += 5) {
        y[i] = x[i]/r;
        y[i + 1] = x[i + 1]/r;
        y[i + 2] = x[i + 2]/r;
        y[i + 3] = x[i + 3]/r;
        y[i + 4] = x[i + 4]/r;
    }
}


void scalar_sub (double * x, double r, int length, double * y) {
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i] -= r * x[i];
    }
    for(; i < length; i += 5) {
        y[i] -= r * x[i];
        y[i + 1] -= r * x[i + 1];
        y[i + 2] -= r * x[i + 2];
        y[i + 3] -= r * x[i + 3];
        y[i + 4] -= r * x[i + 4];
    }
}

double dot_product (double * x, double * y, int length) {
    double sum = 0;
    for(int i = 0; i < length; i++) {
        sum += x[i] * y[i];
    }
    return sum;
}
