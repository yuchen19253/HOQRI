#include "keyFunction.h"

void QR(double ** a, int m, int n){
    int i, row;
    double anorm, tol = 10e-7;

    double ** r = new double*[n];
    for(int index = 0; index < n; index++) {
        r[index] = new double[n];
    }

    for(i = 0; i < n; i++) {
        double nm = 0;
        for(row=0; row<m; row++)
            nm += a[row][i]*a[row][i];
        r[i][i] = sqrt(nm);                  // r_ii = ||a_i||

        if(r[i][i] > tol) {
            for(row=0; row<m; row++)
                a[row][i] /= r[i][i];
        }
        else if(i == 0) { // set a[0] = [1 0 0 ... 0]^T
            a[0][i] = 1;
            for(row = 1; row < m; row++) {
                a[row][i] = 0;
            }
        }
        else{ // need to choose a_i orthogonal to < a_1, ... a_{i-1} >
            for(row=0; row<m%5; row++)
                a[row][i] = -a[i][0] * a[row][0];
            for(; row < m; row+=5) {
                a[row][i] = -a[i][0] * a[row][0];
                a[row+1][i] = -a[i][0] * a[row+1][0];
                a[row+2][i] = -a[i][0] * a[row+2][0];
                a[row+3][i] = -a[i][0] * a[row+3][0];
                a[row+4][i] = -a[i][0] * a[row+4][0];
            }
            a[i][i] += 1;

            for(int col = 1; col < i; col++) {
                for(row=0; row<m%5; row++)
                    a[row][i] -= a[i][col] * a[row][col];
                for(; row < m; row+=5) {
                    a[row][i] -= a[i][col] * a[row][col];
                    a[row+1][i] -= a[i][col] * a[row+1][col];
                    a[row+2][i] -= a[i][col] * a[row+2][col];
                    a[row+3][i] -= a[i][col] * a[row+3][col];
                    a[row+4][i] -= a[i][col] * a[row+4][col];
                }
            }

            nm = 0;
            for(row=0; row<m; row++)
                nm += a[row][i]*a[row][i];
            anorm = sqrt(nm);

            for(row=0; row<m; row++)
                a[row][i] /= anorm;
        }

        for(int col = i+1; col < n; col++) {
            double innerProd = 0;
            for(row=0; row<m; row++)
                innerProd += a[row][i] * a[row][col];
            r[col][i] = innerProd; // r_ij = a_i*a_j
            for(row=0; row<m%5; row++)
                a[row][col] -= r[col][i] * a[row][i];
            for(; row < m; row+=5) {
                a[row][col] -= r[col][i] * a[row][i];       // a_j -= r_ij*a_i
                a[row+1][col] -= r[col][i] * a[row+1][i];
                a[row+2][col] -= r[col][i] * a[row+2][i];
                a[row+3][col] -= r[col][i] * a[row+3][i];
                a[row+4][col] -= r[col][i] * a[row+4][i];
            }
        }
    }
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

