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

using namespace std;

int Lab5_Solution();

void GenerateMatrix(vector<vector<float> >& matr, int sz);

int Jacobi(int sz, vector<vector<float> >& matr, vector<float>& B, vector<float>& x, float eps, int k_max, float w);

int Relaxation(int sz, vector<vector<float> >& matr, vector<float>& B, vector<float>& x, float eps, int k_max, float w);

void Multiply(int sz, vector<vector<float> >& matr, vector<float>&X, vector<float>& B);

void PrintMatrix(int sz, vector<vector<float> >& matr, vector<float>& B);
