//
// Created by yuchen on 3/14/22.
//

#ifndef HOQRI_LINEARALGEBRA_H
#define HOQRI_LINEARALGEBRA_H

#include <cmath>
using namespace std;

//print map of sparse tensor
void print_sptensor(string name, map<tuple<int,int,int>,double> &mapOfTuple)
{
    cout << name << " is: " << endl;
    if(mapOfTuple.empty()){
        cout<<"It's empty!"<<endl;
    } else {
        cout << "Index             " <<
             "Value\n";

        // Iterating over map using range-based loop
        for (auto pr: mapOfTuple)
            // pr points to current pair of mapOfTuple
            cout << get<0>(pr.first) <<
                 " " << get<1>(pr.first) <<
                 " " << get<2>(pr.first) <<
                 " " << pr.second << "\n";
    }
}

void print_densetensor(string name, double ***X, int J[]){
    std::cout << "Dense tensor "<< name << " = " << std::endl;
    for(int i=0; i<J[2]; i++){
        for(int j=0; j<J[0]; j++){
            for(int k=0; k<J[1]; k++) {
                printf("%9.6lg ", X[j][k][i]);
            }
            cout<<endl;
        }
        cout<< endl;
    }
    std::cout << std::endl;
}

//print matrix
void print_matrix(string name, double ** X, int I, int J){
    cout << name << " is: " << endl;
    for(int i=0; i<I; i++){
        for(int j=0; j<J; j++){
            printf("%9.6lg ", X[i][j]);
        }
        cout<< endl;
    }
    std::cout << std::endl;
}

/* ----------------------- norm ----------------------- */
/*  l2-norm.

    Input variables:
        x     : pointer to array.
        length: number of entries in x.                 */

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

/* ----------------------- scalar_div ----------------------- */
/*
    Input variables:
        x     : pointer to array whose components are to be
                 divided by r and stored in second array, y.
        r     : scalar used in division.
        length: number of entries in x and in y.
        y     : pointer to array in which the components
                 of x are to be stored.                        */

void scalar_div (double * x, double r, int length, double * y) {
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


/* ----------------------- scalar_sub ----------------------- */
/*
    Input variables:
        x     : pointer to array whose components are to be
                 multiplied by r then subtracted from the
                 components of the second array, y.
        r     : scalar used in multiplication.
        length: number of entries in x and in y.
        y     : pointer to array in which the components
                 of x are to be stored.                       */

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



/* --------------------- dot_product --------------------- */
/*  Inner production.

    Input variables:
        x     : pointer to first array.
        y     : pointer to second array.
        length: number of entries in x and in y.          */

double dot_product (double * x, double * y, int length) {
    double sum = 0;
    for(int i = 0; i < length; i++) {
        sum += x[i] * y[i];
    }
    return sum;
}



/* ----------------------- gramSchmidt ----------------------- */
/*  Given a matrix A of dimension m by n, this algorithm
    computes a QR decomposition of A, where Q is a unitary
    m by n matrix and R is a n by n upper triangular matrix
    and A = QR.

    Input variables:
        a   : pointer to array of arrays, the ith array of
                which should correspond to the ith column of the
                matrix A. During the algorithm, the columns of Q
                will replace the columns of A.
        r   : pointer to array of arrays in which the ith
                column of the upper triangular matrix R will be
                stored in the ith subarray of r.
        m   : number of columns in A.
        n   : number of rows in A.
        thin: TRUE  => thin QR factorization computed
              FALSE => full QR factorization computed

    Features: This implementation has time complexity O(m n^2)
    and requires O(1) additional memory.

    Remarks: Due to the nature of the problem, if A is nearly
    rank-deficient then the resulting columns of Q may not
    exhibit the orthogonality property.                        */
void gramSchmidt (double ** a, double ** r, int m, int n, bool full) {
    int i, j;
    double anorm, tol = 10e-7;

    for(i = 0; i < n; i++) {
        r[i][i] = norm(a[i], m);                  // r_ii = ||a_i||

        if(r[i][i] > tol) {
            scalar_div(a[i], r[i][i], m, a[i]);   // a_i = a_i/r_ii
        }
        else if(i == 0) { // set a[0] = [1 0 0 ... 0]^T
            a[i][0] = 1;
            for(j = 1; j < m; j++) {
                a[i][j] = 0;
            }
        }
        else{ // need to choose a_i orthogonal to < a_1, ... a_{i-1} >
            for(j = 0; j < m; j++) {
                a[i][j] = -a[0][i] * a[0][j];
            }
            a[i][i] += 1;

            for(j = 1; j < i; j++) {
                scalar_sub(a[j], a[j][i], m, a[i]);
            }

            anorm = norm(a[i], m);
            scalar_div(a[i], anorm, m, a[i]);
        }

        for(j = i+1; j < n; j++) {
            r[j][i] = dot_product(a[i], a[j], m); // r_ij = a_i*a_j
            scalar_sub(a[i], r[j][i], m, a[j]);   // a_j -= r_ij a_i
        }
    }

    /* if full QR factorization requested, we choose remaining
       columns of Q so that the m columns of Q form an
       orthonormal set                                          */
    if(full) {
        for(; i < m; i++) {
            for(j = 0; j < m; j++) {
                a[i][j] = -a[0][i] * a[0][j];
            }
            a[i][i] += 1;

            for(j = 1; j < i; j++) {
                scalar_sub(a[j], a[j][i], m, a[i]);
            }

            anorm = norm(a[i], m);
            scalar_div(a[i], anorm, m, a[i]);
        }
    }
}

void qr(double ** Q, int m, int n){
    double ** a = new double*[n];
    double ** r = new double*[n];
    for(int index = 0; index < n; index++) {
        a[index] = new double[m];
        r[index] = new double[n];
    }
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            a[j][i] = Q[i][j];
        }
    }
    gramSchmidt(a,r,m,n,0);

    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            Q[i][j] = a[j][i];
        }
    }
//
//    std::cout << "func Q = " << std::endl;
//    for(int i = 0; i < m; i++) {
//        for(int j = 0; j < n; j++) {
//            printf("%9.6lg ", Q[i][j]);
//        }
//        std::cout << std::endl;
//    }
//    std::cout << std::endl;
//
//    std::cout << "func R = " << std::endl;
//    for(int i = 0; i < n; i++) {
//        for(int j = 0; j < n; j++) {
//            printf("%9.6lg ", r[j][i]);
//        }
//        std::cout << std::endl;
//    }
//    std::cout << std::endl;
}



/* ----------------------- vec_copy ----------------------- */
/*  Given two arrays of the same length and their length,
    this function stores the values from the first array
    in the second array.

    Input variables:
        x     : pointer to array whose entries are to be
                 copied.
        y     : pointer to array in which the components
                 of x are to be stored.
        length: number of entries in x and in y.

    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.          */

void vec_copy (double * x, double * y, int length) {
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i] = x[i];
    }
    for(; i < length; i += 5) {
        y[i] = x[i];
        y[i + 1] = x[i + 1];
        y[i + 2] = x[i + 2];
        y[i + 3] = x[i + 3];
        y[i + 4] = x[i + 4];
    }
}


/* ------------------- partialvec_copy ------------------- */
/*  Given two arrays, the length of the second array, and
    an index this function stores the values from the
    subarray x[index : index + length] in the array
    y[0 : length].

    Input variables:
        x     : pointer to array whose entries are to be
                 copied.
        y     : pointer to array in which the components
                 of x are to be stored.
        length: number of entries in y.
        index : starting index of subarray of x to be
                copied to y.

    Example: Suppose x is a pointer to the array
    {1, 2, 3, 4, 5}, y is a pointer to the array {0, 0, 0},
    length = 3, and index = 2. Then after executing
    partialvec_copy(x, y, 3, 2), the array pointed to by
    y is now {3, 4, 5}.

    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.         */

void partialvec_copy (double * x, double * y, int length, int index) {
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i] = x[i + index];
    }
    for(; i < length; i += 5) {
        y[i] = x[i + index];
        y[i + 1] = x[i + index + 1];
        y[i + 2] = x[i + index + 2];
        y[i + 3] = x[i + index + 3];
        y[i + 4] = x[i + index + 4];
    }
}




/* --------------------- partialscalar_sub --------------------- */
/*  Given two arrays, the length of the second array, a scalar
    value, and an index, this function multiplies the values
    starting at the given index from the first array by the
    scalar value and then subtracts the computed components from
    the components the second array.

    Input variables:
        x     : pointer to array whose components are to be
                 multiplied by r then subtracted from the
                 components of the second array, y.
        r     : scalar used in multiplication.
        length: number of entries in y.
        index :
        y     : pointer to array in which the components
                 of x are to be stored.

    Example: Suppose x is a pointer to the array
    {1, 2, 3, 4, 5}, y is a pointer to the array {0, 0, 0},
    length = 3, r = -1, and index = 2. Then after executing
    partialscalar_sub(x, -1, 3, 2, y), the array pointed to
    by y is now {-3, -4, -5}.

    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.               */

void partialscalar_sub (double * x, double r, int length,
                        int index, double * y)
{
    int i, length5;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        y[i + index] -= r * x[i];
    }
    for(; i < length; i += 5) {
        y[i + index] -= r * x[i];
        y[i + index + 1] -= r * x[i + 1];
        y[i + index + 2] -= r * x[i + 2];
        y[i + index + 3] -= r * x[i + 3];
        y[i + index + 4] -= r * x[i + 4];
    }
}



/* ------------------ partialdot_product ------------------ */
/*  Given two arrays of the same length, their length, and
    an index this function returns the dot product of the
    two subarrays x[index : length] and y[index : length].

    Input variables:
        x     : pointer to first array.
        y     : pointer to second array.
        length: number of entries in x and in y.
        index : starting index for subarrays.

    Example: Suppose x is a pointer to the array
    {1, 2, 3, 4}, y is a pointer to the array {5, 6, 7, 8},
    length = 4, and index = 2. Then the value returned by
    executing partialdot_product(x, y, 4, 2) is 53, which
    is computed by
        x[2] * y[2] + x[3] * y[3] = 3 * 7 + 4 * 8
                                  = 21 + 32
                                  = 53.

    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.          */

double partialdot_product (double * x, double * y, int length, int index) {
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for(i = index; i < length5; i++) {
        sum += x[i] * y[i];
    }
    for(; i < length; i += 5) {
        sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
               + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
    }

    return sum;
}


/* -------------------- subdot_product -------------------- */
/*  Given two arrays, the length of the second array, and
    an index this function returns the dot product of the
    two subarrays x[index : index + length] and
    y[0 : length]. It is necessary that index + length is
    at most the length of the first array.

    Input variables:
        x     : pointer to first array.
        y     : pointer to second array.
        length: number of entries in y.
        index : starting index for subarray of x.

    Example: Suppose x is a pointer to the array
    {1, 2, 3, 4, 5}, y is a pointer to the array
    {-1, -2, -3}, length = 3, and index = 2. Then the value
    returned by executing subdot_product(x, y, 3, 2) is 53,
    which is computed by
            x[2] * y[0] + x[3] * y[1] + x[4] * y[2]

          =  3   *  -1  +  4   *  -2  +  5   *  -3

          = -    3      -      8      -      15

          = -26.

    Features: This implementation has time complexity
    O(length) and requires O(1) additional memory.          */

double subdot_product (double * x, double * y, int length, int index) {
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for(i = 0; i < length5; i++) {
        sum += x[i + index] * y[i];
    }
    for(; i < length; i += 5) {
        sum += x[i + index] * y[i] + x[i + index + 1] * y[i + 1]
               + x[i + index + 2] * y[i + 2]
               + x[i + index + 3] * y[i + 3]
               + x[i + index + 4] * y[i + 4];
    }

    return sum;
}

/* -------------------- find a vector in tensor by index -------------------- */

void vec_in_tensor(double *vec, double ***tensor, int dim[], int length){
    if(dim[0]==-1){
        for(int i=0; i<length;i++){
            vec[i] = tensor[i][dim[1]][dim[2]];
        }
    } else if(dim[1]==-1){
        for(int i=0; i<length;i++){
            vec[i] = tensor[dim[0]][i][dim[2]];
        }
    } else if(dim[2]==-1){
        for(int i=0; i<length;i++){
            vec[i] = tensor[dim[0]][dim[1]][i];
        }
    } else{
        cout<<"Wrong dimension! Expected -1 for mode-n vector."<<endl;
    }
    for(int i=0; i<length;i++){
        cout<< vec[i] << " ";
    }
}

/* -------------------- find a vector in matrix by index -------------------- */

void vec_in_matrix(double *vec, double **matrix, int index, int length){
    for(int i=0; i<length;i++){
        vec[i] = matrix[i][index];
    }
    for(int i=0; i<length;i++){
        cout<< vec[i] << " ";
    }
}


double*** ttm2(map<tuple<int,int,int>,double> spX, double ** U, double ** V, double ** W, int K[]){
    double*** G;
    G = new double ** [K[0]];
    for (int i = 0;i < K[0];i++)
    {
        G[i] = new double*[K[1]];
        for (int j = 0;j < K[1];j++)
        {
            G[i][j] = new double [K[2]];
        }
    }

    for(int i=0; i<K[0]; i++){
        for(int j=0; j<K[1]; j++){
            for(int k=0; k<K[2]; k++) {
                G[i][j][k] = 0;
            }
        }
    }

    map<tuple<int,int,int>,double>::iterator iter;
    for(int i = 0; i<K[0]; i++){
        for(int j = 0; j<K[1]; j++){
            for(int k = 0; k<K[2]; k++){
                double tmp = 0;
                for(iter=spX.begin(); iter!=spX.end(); iter++){
                    int index_i=get<0>(iter->first);
                    int index_j=get<1>(iter->first);
                    int index_k=get<2>(iter->first);
                    tmp += (iter->second) * U[index_i-1][i] * V[index_j-1][j] * W[index_k-1][k];
                }
                G[i][j][k]=tmp;
            }
        }
    }
    return G;
}

double*** ttm(map<tuple<int,int,int>,double> spX, double ** U, double ** V, double ** W, int K[]){
    double***  G;
    G = new double ** [K[0]];
    for (int i = 0;i < K[0];i++)
    {
        G[i] = new double*[K[1]];
        for (int j = 0;j < K[1];j++)
        {
            G[i][j] = new double [K[2]];
        }
    }

    for(int i=0; i<K[0]; i++){
        for(int j=0; j<K[1]; j++){
            for(int k=0; k<K[2]; k++) {
                G[i][j][k] = 0;
            }
        }
    }
    map<tuple<int,int,int>,double>::iterator iter;
    for(iter=spX.begin(); iter!=spX.end(); iter++){
        int index_i=get<0>(iter->first);
        int index_j=get<1>(iter->first);
        int index_k=get<2>(iter->first);
        for(int i = 0; i<K[0]; i++){
            for(int j = 0; j<K[1]; j++) {
                for (int k = 0; k < K[2]; k++) {
                    G[i][j][k]+= (iter->second) * U[index_i-1][i] * V[index_j-1][j] * W[index_k-1][k];
                }
            }
        }
    }

    return G;
}

double*** ttm_update(map<tuple<int,int,int>,double> spX, double ** U, double ** V, double ** W, int K[]){
    double***  G;
    G = new double ** [K[0]];
    for (int i = 0;i < K[0];i++)
    {
        G[i] = new double*[K[1]];
        for (int j = 0;j < K[1];j++)
        {
            G[i][j] = new double [K[2]];
        }
    }

    for(int i=0; i<K[0]; i++){
        for(int j=0; j<K[1]; j++){
            for(int k=0; k<K[2]; k++) {
                G[i][j][k] = 0;
            }
        }
    }
    map<tuple<int,int,int>,double>::iterator iter;
    for(iter=spX.begin(); iter!=spX.end(); iter++){
        int index_i=get<0>(iter->first);
        int index_j=get<1>(iter->first);
        int index_k=get<2>(iter->first);
        double tmp = (iter->second) ;
        for(int i = 0; i<K[0]; i++){
            double tmp_U = tmp * U[index_i - 1][i];
            for(int j = 0; j<K[1]; j++) {
                double tmp_UV = tmp_U * V[index_j - 1][j];
                for (int k = 0; k < K[2]; k++) {
                    G[i][j][k]+=tmp_UV*W[index_k - 1][k];
                }
            }
        }
    }

    return G;
}



map<tuple<int,int,int>,double> spttm(map<tuple<int,int,int>,double> spX, double ** U, double ** V, double ** W, int K[]){
    map<tuple<int,int,int>,double> G;
    map<tuple<int,int,int>,double>::iterator iter;
    for(int i = 0; i<K[0]; i++){
        for(int j = 0; j<K[1]; j++){
            for(int k = 0; k<K[2]; k++){
                tuple<int, int, int> index(i+1,j+1,k+1);
                double tmp = 0;
                for(iter=spX.begin(); iter!=spX.end(); iter++){
                    int index_i=get<0>(iter->first);
                    int index_j=get<1>(iter->first);
                    int index_k=get<2>(iter->first);
                    tmp += (iter->second) * U[index_i-1][i] * V[index_j-1][j] * W[index_k-1][k];
                }
                G[index]=tmp;
            }
        }
    }
    return G;
}

map<tuple<int,int,int>,double> spttm2(map<tuple<int,int,int>,double> spX, double ** U, double ** V, double ** W, int K[]){
    map<tuple<int,int,int>,double> G;
    map<tuple<int,int,int>,double>::iterator iter;
    for(iter=spX.begin(); iter!=spX.end(); iter++){
        int index_i=get<0>(iter->first);
        int index_j=get<1>(iter->first);
        int index_k=get<2>(iter->first);
        for(int i = 0; i<K[0]; i++){
            for(int j = 0; j<K[1]; j++) {
                for (int k = 0; k < K[2]; k++) {
                    double tmp = 0;
                    tuple<int, int, int> index(i+1,j+1,k+1);
                    if(!G.count(index))
                        G[index] = 0;
                    tmp = (iter->second) * U[index_i - 1][i] * V[index_j - 1][j] * W[index_k - 1][k];
                    G[index]+=tmp;
                }
            }
        }
    }

    return G;
}

map<tuple<int,int,int>,double> spttm_mode(map<tuple<int,int,int>,double> spT, double ** matrix, int J[], int change_size, int mode){
    map<tuple<int,int,int>,double> G;
    map<tuple<int,int,int>,double> X(spT);
//    tuple<int, int, int> index(index1,index2,index3);
//    double value;
//    G[index] = value;

    if(mode==1){
        while (!X.empty())
        {
            //extract vector which the first item of current X in
            tuple<int,int,int> t = X.begin()->first;
            map<int,double> vector;
            int index_i=get<0>(t);
            int index_j=get<1>(t);
            int index_k=get<2>(t);
//            cout << index_i<<' '<<index_j<<' '<<index_k<<endl;

            vector[index_i]=X.begin()->second;
//            cout<<"First add "<<X.begin()->second<<endl;

            //find same j and k
            map<tuple<int,int,int>,double>::iterator iter;
            for(iter=X.begin(); iter!=X.end(); iter++)
            {
                if(get<1>(iter -> first)==index_j && get<2>(iter -> first)==index_k){
                    vector[get<0>(iter->first)]=iter->second;
//                    cout<<"add "<<iter->second<<endl;
                }
            }

//            cout<<"----END extract vector----"<<vector.begin()->first<<", "<< vector.begin()->second<<endl;

            map<int,double>::iterator vector_iter;
            for(int new_i=1; new_i<=change_size;new_i++){
//                cout<<"in update for loop"<<endl;
                tuple<int, int, int> index(new_i,index_j,index_k);
                double tmp=0;
                for(vector_iter=vector.begin(); vector_iter!=vector.end(); vector_iter++){
//                    cout<<"in inner product for loop"<<endl;
                    X.erase(tuple<int,int,int>(vector_iter->first,index_j,index_k));
//                    cout<<"print current X"<<endl;
//                    print(X);
                    tmp += (vector_iter->second) * matrix[(vector_iter->first)-1][new_i-1];
                }
                G[index]=tmp;
//                cout<<"ADD "<< get<0>(index)<<", "<< get<1>(index)<<", "<< get<2>(index)<<" to G: "<<tmp<<endl;
            }
//            cout<<"--------G--------"<<endl;
//            print(G);
//            cout<<endl;
        }

    } else if(mode==2){
        while (!X.empty())
        {
            //extract vector which the first item of current X in
            tuple<int,int,int> t = X.begin()->first;
            map<int,double> vector;
            int index_i=get<0>(t);
            int index_j=get<1>(t);
            int index_k=get<2>(t);

            vector[index_j]=X.begin()->second;

            //find same j and k
            map<tuple<int,int,int>,double>::iterator iter;
            for(iter=X.begin(); iter!=X.end(); iter++)
            {
                if(get<0>(iter -> first)==index_i && get<2>(iter -> first)==index_k){
                    vector[get<1>(iter->first)]=iter->second;
//                    cout<<"add "<<iter->second<<endl;
                }
            }

//            cout<<"----END extract vector----"<<vector.begin()->first<<", "<< vector.begin()->second<<endl;

            map<int,double>::iterator vector_iter;
            for(int new_j=1; new_j<=change_size;new_j++){
//                cout<<"in update for loop"<<endl;
                tuple<int, int, int> index(index_i,new_j,index_k);
                double tmp=0;
                for(vector_iter=vector.begin(); vector_iter!=vector.end(); vector_iter++){
//                    cout<<"in inner product for loop"<<endl;
                    X.erase(tuple<int,int,int>(index_i,vector_iter->first,index_k));
//                    cout<<"print current X"<<endl;
//                    print(X);
                    tmp += (vector_iter->second) * matrix[(vector_iter->first)-1][new_j-1];
                }
                G[index]=tmp;
//                cout<<"ADD "<< get<0>(index)<<", "<< get<1>(index)<<", "<< get<2>(index)<<" to G: "<<tmp<<endl;
            }
//            cout<<"--------G--------"<<endl;
//            print(G);
//            cout<<endl;
        }

    } else if(mode==3){
        while (!X.empty())
        {
            //extract vector which the first item of current X in
            tuple<int,int,int> t = X.begin()->first;
            map<int,double> vector;
            int index_i=get<0>(t);
            int index_j=get<1>(t);
            int index_k=get<2>(t);

            vector[index_k]=X.begin()->second;

            //find same j and k
            map<tuple<int,int,int>,double>::iterator iter;
            for(iter=X.begin(); iter!=X.end(); iter++)
            {
                if(get<0>(iter -> first)==index_i && get<1>(iter -> first)==index_j){
                    vector[get<2>(iter->first)]=iter->second;
//                    cout<<"add "<<iter->second<<endl;
                }
            }

//            cout<<"----END extract vector----"<<vector.begin()->first<<", "<< vector.begin()->second<<endl;

            map<int,double>::iterator vector_iter;
            for(int new_k=1; new_k<=change_size;new_k++){
//                cout<<"in update for loop"<<endl;
                tuple<int, int, int> index(index_i,index_j,new_k);
                double tmp=0;
                for(vector_iter=vector.begin(); vector_iter!=vector.end(); vector_iter++){
//                    cout<<"in inner product for loop"<<endl;
                    X.erase(tuple<int,int,int>(index_i,index_j,vector_iter->first));
//                    cout<<"print current X"<<endl;
//                    print(X);
                    tmp += (vector_iter->second) * matrix[(vector_iter->first)-1][new_k-1];
                }
                G[index]=tmp;
//                cout<<"ADD "<< get<0>(index)<<", "<< get<1>(index)<<", "<< get<2>(index)<<" to G: "<<tmp<<endl;
            }
//            cout<<"--------G--------"<<endl;
//            print(G);
//            cout<<endl;
        }

    } else{
        cout<<"Wrong mode! Expected 1, 2 or 3."<<endl;
    }
    return G;
}

double *** sp2dense(map<tuple<int,int,int>,double> spX, int J[]){
    double ***X;
    X = new double ** [J[0]];
    for (int i = 0;i < J[0];i++)
    {
        X[i] = new double*[J[1]];
        for (int j = 0;j < J[1];j++)
        {
            X[i][j] = new double [J[2]];
        }
    }
    for(int i=0; i<J[0]; i++){
        for(int j=0; j<J[1]; j++){
            for(int k=0; k<J[2]; k++) {
                X[i][j][k] = 0;
            }
        }
    }
    map<tuple<int,int,int>,double>::iterator iter;
    for(iter=spX.begin(); iter!=spX.end(); iter++){
        int i = get<0>(iter->first);
        int j = get<1>(iter->first);
        int k = get<2>(iter->first);
        X[i-1][j-1][k-1]=iter->second;
//        cout<<"insert "<<i<<", "<<j<<", "<<k<<"=>"<<iter->second<<endl;
    }

    return X;
}

double** TTMcTC(map<tuple<int,int,int>,double> X, double*** G, double **U, double **V, double **W, int J[], int K[], int mode){
    double ** A = new double*[J[mode-1]];
    for(int index = 0; index < J[mode-1]; index++) {
        A[index] = new double[K[mode-1]];
    }
    for(int j=0; j<J[mode-1]; j++){
        for(int k=0; k<K[mode-1]; k++) {
            A[j][k] = 0;
        }
    }
    map<tuple<int,int,int>,double>::iterator iter;
    for(iter=X.begin(); iter!=X.end(); iter++){
        int i1 = get<0>(iter->first)-1;
        int i2 = get<1>(iter->first)-1;
        int i3 = get<2>(iter->first)-1;
        if(mode==1){
            for(int k1=0; k1<K[0]; k1++) {
                for (int k2 = 0; k2 < K[1]; k2++)
                    for (int k3 = 0; k3 < K[2]; k3++)
                        A[i1][k1] += (iter->second) * G[k1][k2][k3] * V[i2][k2] * W[i3][k3];
            }

        } else if(mode==2){
            for(int k1=0; k1<K[0]; k1++) {
                for (int k2 = 0; k2 < K[1]; k2++)
                    for (int k3 = 0; k3 < K[2]; k3++)
                        A[i2][k2] += (iter->second) * G[k1][k2][k3] * U[i1][k1] * W[i3][k3];
            }

        } else if(mode==3){
            for(int k1=0; k1<K[0]; k1++) {
                for (int k2 = 0; k2 < K[1]; k2++)
                    for (int k3 = 0; k3 < K[2]; k3++)
                        A[i3][k3] += (iter->second) * G[k1][k2][k3] * U[i1][k1] * V[i2][k2];
            }

        } else{
            cout<<"Wrong mode! Expected 1, 2 or 3."<<endl;
        }
    }

    return A;
}

double** TTMcTC_update(map<tuple<int,int,int>,double> X, double*** G, double **U, double **V, double **W, int J[], int K[], int mode){
    double ** A = new double*[J[mode-1]];
    for(int index = 0; index < J[mode-1]; index++) {
        A[index] = new double[K[mode-1]];
    }
    for(int j=0; j<J[mode-1]; j++){
        for(int k=0; k<K[mode-1]; k++) {
            A[j][k] = 0;
        }
    }
    map<tuple<int,int,int>,double>::iterator iter;
    for(iter=X.begin(); iter!=X.end(); iter++){
        int i1 = get<0>(iter->first)-1;
        int i2 = get<1>(iter->first)-1;
        int i3 = get<2>(iter->first)-1;
        double tmp = (iter->second);
        if(mode==1){
            for (int k3 = 0; k3 < K[2]; k3++){
                double tmp_W = tmp * W[i3][k3];
                for (int k2 = 0; k2 < K[1]; k2++) {
                    double tmp_VW = tmp_W * V[i2][k2];
                     for(int k1=0; k1<K[0]; k1++) {
                        A[i1][k1] += G[k1][k2][k3] * tmp_VW ;
                    }
                }
            }

        } else if(mode==2){
            for(int k1=0; k1<K[0]; k1++) {
                double tmp_U = tmp * U[i1][k1];
                for (int k3 = 0; k3 < K[2]; k3++)  {
                    double tmp_UW = tmp_U* W[i3][k3];
                    for (int k2 = 0; k2 < K[1]; k2++){
                        A[i2][k2] += G[k1][k2][k3] * tmp_UW;
                    }
                }
            }

        } else if(mode==3){
            for(int k1=0; k1<K[0]; k1++) {
                double tmp_U = tmp * U[i1][k1];
                for (int k2 = 0; k2 < K[1]; k2++) {
                    double tmp_UV = tmp_U * V[i2][k2];
                    for (int k3 = 0; k3 < K[2]; k3++) {
                        A[i3][k3] += G[k1][k2][k3] * tmp_UV;
                    }
                }
            }

        } else{
            cout<<"Wrong mode! Expected 1, 2 or 3."<<endl;
        }
    }

    return A;
}

double** TTMcTC2(map<tuple<int,int,int>,double> X, double*** G, double **U, double **V, double **W, int J[], int K[], int mode){
    double ** A = new double*[J[mode-1]];
    for(int index = 0; index < J[mode-1]; index++) {
        A[index] = new double[K[mode-1]];
    }
    for(int j=0; j<J[mode-1]; j++){
        for(int k=0; k<K[mode-1]; k++) {
            A[j][k] = 0;
        }
    }
    map<tuple<int,int,int>,double>::iterator iter;
    for(int k1=0; k1<K[0]; k1++) {
        for (int k2 = 0; k2 < K[1]; k2++) {
            for (int k3 = 0; k3 < K[2]; k3++) {
                for(iter=X.begin(); iter!=X.end(); iter++){
                    int i1 = get<0>(iter->first)-1;
                    int i2 = get<1>(iter->first)-1;
                    int i3 = get<2>(iter->first)-1;
                    if(mode==1){
                        A[i1][k1] += (iter->second) * G[k1][k2][k3] * V[i2][k2] * W[i3][k3];

                    } else if(mode==2){
                        A[i2][k2] += (iter->second) * G[k1][k2][k3] * U[i1][k1] * W[i3][k3];

                    } else if(mode==3){
                        A[i3][k3] += (iter->second) * G[k1][k2][k3] * U[i1][k1] * V[i2][k2];

                    } else{
                        cout<<"Wrong mode! Expected 1, 2 or 3."<<endl;
                    }
                }

            }
        }
    }
    return A;
}


#endif //HOQRI_LINEARALGEBRA_H