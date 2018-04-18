#include "Rational.h"

int Rational::NOD(int v1, int v2)
{
	if (v1 > v2)
	{
		swap(v1, v2);
	}
	return 1;
}

void Rational::reduct()
{
	int nod = NOD(a, b);
	if (nod != 1)
	{
		a /= nod;
		b /= nod;
	}
}

Rational::Rational()
{
	a = 0;
	b = 1;
}

Rational::Rational(int _a)
{
	a = _a;
	b = 1;
}

Rational::Rational(int _a, unsigned int _b)
{
	a = _a;
	b = _b;
}

Rational::~Rational()
{}


void Rational::print()
{
	cout << a << "/" << b;
}

bool Rational::isEven()
{
	return !(a % 2);
}



