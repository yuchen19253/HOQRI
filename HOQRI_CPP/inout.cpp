#include "inout.h"
#include <string.h>

int* argv_split(char* argv){
    int i,j,k;
    char* result = NULL;
    result = strtok(argv, ",");
    std::stringstream s1(result);
    s1 >> i;
    result = strtok(NULL, ",");
    std::stringstream s2(result);
    s2 >> j;
    result = strtok(NULL, ",");
    std::stringstream s3(result);
    s3 >> k;
    return new int[3]{i,j,k};
}

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

