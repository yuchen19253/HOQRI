#ifndef INOUT_H
#define INOUT_H
#include<iostream>
#include<sstream>
#include<string>
#include<map>
using namespace std;

// split command line arguments into dimension or rank
int* argv_split(char* argv);

//print map of sparse tensor
void print_sptensor(string name, map<tuple<int,int,int>,double> &mapOfTuple);

//print the densetensor
void print_densetensor(string name, double ***X, int J[]);

//print matrix
void print_matrix(string name, double ** X, int I, int J);



#endif //INOUT_H
