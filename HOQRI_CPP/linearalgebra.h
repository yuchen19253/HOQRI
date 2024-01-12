#ifndef HOQRI_LINEARALGEBRA_H
#define HOQRI_LINEARALGEBRA_H
#include <cmath>
using namespace std;

/* ----------------------- norm ----------------------- */
/*  l2-norm.
        first input: pointer to vector/tensor.
        second input: size of the first input.          */

double norm (double * x, int length);

double norm_tensor(double *** X, int J[]);

/* y_i = x_i / r (element-wise)
 * length: number of entries in x and in y. */
void scalar_div(double * x, double r, int length, double * y);


/* y_i = y_i - x_i * r (element-wise)
 * length: number of entries in x and in y.  */
void scalar_sub (double * x, double r, int length, double * y);


/*  Inner production of vector x and y.
 * sum = x^Ty
 * length: number of entries in x and in y.          */
double dot_product (double * x, double * y, int length);

#endif //HOQRI_LINEARALGEBRA_H
