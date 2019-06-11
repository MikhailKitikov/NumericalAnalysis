#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <iomanip>
#include <vector>
#include <cassert>
#include <ctime>
#include <stdlib.h>

using namespace std;

int Lab6_Solution();

void GenerateMatrix(vector<vector<float> >& matr, int sz);

bool Danilevski(int sz, vector<vector<float> >& matr, vector<vector<float> >& M);

void PrintMatrix(int sz, vector<vector<float> >& matr, ostream& os = cout);

float Trace(int n, vector<vector<float> >& matr);
