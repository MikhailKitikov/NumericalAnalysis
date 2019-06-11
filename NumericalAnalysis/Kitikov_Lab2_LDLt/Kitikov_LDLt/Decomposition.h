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

void GenerateMatrix(vector<vector<float> >& matr, int sz, int k);

bool Decompose(int sz, vector<vector<float> >& matr);

float Norm(vector<float>& vec);

void CalculateAnswer(int sz, vector<vector<float> >& matr, vector<float>& X, vector<float>& B, vector<float>& ans);
