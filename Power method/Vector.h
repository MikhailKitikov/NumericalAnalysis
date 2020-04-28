#pragma once
#include "Libs.h"
using namespace std;

class Vector {

private:

	vector<float> vec;

public:

	/// конструктор по умолчанию
	Vector() {}

	/// конструктор по размеру
	Vector(int n) {
		assert(n > 0);
		vec.assign(n, 0);
	}

	/// конструктор по размеру + заполнение
	Vector(int n, float val) {
		assert(n > 0);
		vec.assign(n, val);
	}

	/// конструктор копирования
	Vector(Vector& arg) {
		vec = arg.vec;
	}

	/// изменение размера
	void resize(int n) {
		assert(n > 0);
		vec.resize(n);
	}

	/// размер вектора
	int size() { return vec.size(); }

	/// присваивание
	Vector& operator = (Vector& arg) {
		vec = arg.vec;
		return *this;
	}

	/// доступ к координате вектора
	float& operator [](int i) {
		assert(i >= 0 && i < this->size());
		return vec[i];
	}

	/// сложение векторов
	Vector operator + (Vector& arg) {
		Vector temp(*this);
		for (int i = 0; i < this->size(); ++i)
			temp[i] += arg[i];
		return temp;
	}

	/// вычитание векторов
	Vector operator - (Vector& arg) {
		Vector temp(*this);
		for (int i = 0; i < this->size(); ++i)
			temp[i] -= arg[i];
		return temp;
	}

	/// кубическая норма
	float cubicNorm() {
		float res = fabs(vec[0]);
		for (int i = 1; i < vec.size(); ++i) {
			res = max(res, fabs(vec[i]));
		}
		return res;
	}

	/// умножение на скаляр
	Vector operator * (float val) {
		Vector temp(*this);
		for (int i = 0; i < this->size(); ++i)
			temp[i] *= val;
		return temp;
	}

	/// деление на скаляр
	Vector operator / (float val) {
		assert(val != 0);
		Vector temp(*this);
		for (int i = 0; i < this->size(); ++i)
			temp[i] /= val;
		return temp;
	}

	/// индекс max координаты
	int maxIndex() {
		float vnorm = this->cubicNorm();
		return find_if(vec.begin(), vec.end(), [vnorm](float a) { return fabs(a) == vnorm; }) - vec.begin();
	}

};

/// печать вектора
ostream& operator << (ostream& os, Vector& arg) {
	cout << "(";
	for (int i = 0; i < arg.size(); ++i) {
		cout << setprecision(5) << arg[i] << ", ";
	}
	cout << ")";
	return os;
}

/// скалярное произведение
float ScalarProduct(Vector& a, Vector& b) {
	assert(a.size() == b.size());
	float res = 0;
	for (int i = 0; i < a.size(); ++i) {
		res += a[i] * b[i];
	}
	return res;
}