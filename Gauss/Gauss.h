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

int Task1();

int Task2();

void GenerateMatrix(vector<vector<float> >& matr, int sz, int k);

void GaussWithoutSelection(int sz, vector<vector<float> >& matr, vector<float>& B);

void GaussWithSelection(int sz, vector<vector<float> >& matr, vector<float>& B);

float Norm(vector<float>& vec);

void Multiply(int sz, vector<vector<float> >& matr, vector<float>&X, vector<float>& B);

void CalculateAnswer(int sz, vector<vector<float> >& matr, vector<float>& X, vector<float>& B, vector<float>& ans);
