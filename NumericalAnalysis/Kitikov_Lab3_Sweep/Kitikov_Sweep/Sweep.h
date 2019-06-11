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

int Solve();

void GenerateMatrix(vector<vector<float> >& matr, int sz, int k, int m);

bool RightSweep(int sz, vector<vector<float> >& matr, vector<float>& b, vector<float>& alpha, vector<float>& beta);

float Norm(vector<float>& vec);

void CalculateAnswer(int sz, vector<float>& alpha, vector<float>& beta, vector<float>& ans);
