#pragma once

#include "iostream"
using namespace std;

class Rational
{
private:
	int a;
	unsigned int b;

	int NOD(int v1, int v2);
	void reduct();

public:
	Rational();
	Rational(int _a);
	Rational(int _a, unsigned int _b);
	~Rational();

	void print();
	bool isEven();
};

