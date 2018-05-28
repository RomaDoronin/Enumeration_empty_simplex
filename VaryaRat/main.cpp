#define _CRT_SECURE_NO_WARNINGS // ��� ����������� ������ [c4996 deprecate]

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <ilcplex/ilocplex.h> // IBM ILOG CPLEX


/*----------------------------------------------------*/
// ���������
// ������� SMP - ��� ������ ��������� ��������� � ���������, MIP - � �������� ��� CPLEX

//#define SMP_DEBUG
//#define SMP_OLD_CODE
#define SMP_INIT_MODE_0 0
#define SMP_INIT_MODE_1 1

#define TEST_MOD

// ���� ������� �������
enum Target {
	MIP_MAX,
	MIP_MIN
};

// �������������� ����������� ��� MIP
enum AdditRestr {
	MIP_NONE_REST,
	MIP_XG0, // X greater than zero
	MIP_XL0 // X less than zero
};

/*----------------------------------------------------*/
// �������

#define SMP_NEW(mas, n, type) mas = new type*[n]; \
for (int newCount = 0; newCount < n; newCount++) \
	mas[newCount] = new type[n]

#define SMP_DELETE(mas, n) for (int deleteCount = 0; deleteCount < n; deleteCount++) \
delete[] mas[deleteCount]

#define SMP_COMPARE(val1, val2) if (val1 == val2) std::cout << "[OK] \n"; \
else std::cout << "[FAILED] \n"

/*----------------------------------------------------*/
// ��������� ������� � �������
class Simplex;
bool SimplexIsEmpty(Simplex *S);
void GetAllCovers(Simplex *S, int N);
int CalcSimplexWidth(Simplex *S);
int Determinant(int **mas, int m);
double** MatIntToDouble(int **A, int n);
double* VecIntToDouble(int *A, int n);
double* Gauss(double **a, double *y, int n);

/*----------------------------------------------------*/

// ����� ��������� [ Ax <= b ] -> [ AQx <= b ] -> [ Hx <= b ] (����� ������)
class Simplex {
public:
	int **H;         // ������� H
	int *b;          // ������ b

	int Multip;      // ���������
	int n;           // �����������

					 // ������ � (cTx <= c0)
	int *cover;
	int c0;

public:
	// �����������
	Simplex(int _n, int _M) : n(_n), Multip(_M)
	{
		b = new int[n];
		cover = new int[n];
		SMP_NEW(H, n, int);

		InitDef();
	}

	// ����������
	~Simplex()
	{
		delete[] b;
		delete[] cover;
		SMP_DELETE(H, n);
	}

	// ���������� ������� � ���������� �������
	int** ReplaceMatrixString(int rI)
	{
		int **Tm;
		SMP_NEW(Tm, n, int);

		for (int i = 0; i < n; i++)
			if (i == rI)
				for (int j = 0; j < n; j++)
					Tm[i][j] = cover[j];
			else
				for (int j = 0; j < n; j++)
					Tm[i][j] = H[j][i];

		return Tm;
	}

	// 
	int* ReplaceVectorString(int rI)
	{
		int *Tv = new int[n];

		for (int i = 0; i < n; i++)
			if (i == rI)
				Tv[i] = c0;
			else
				Tv[i] = b[i];

		return Tv;
	}

	/*----------------------------------------------------*/
	// ������ ��� ������ � �������

#ifdef SMP_OLD_CODE
	// ����������� ������� ����� ����� � n-������ ���������������, ������������ �� ����� ������� H
	void EnumerationIntPointsInParal(int m, int *C)
	{
		m++; // �������� ���������� �� ��������� ��������, �� ���� ������� ����� �� �� �� ����� �� n 
			 // ������� ���� ����� �����
		if (m != n)
			for (int i = 0; i <= ColomnSum(m); i++)
			{
				C[m] = i;
				EnumerationIntPointsInParal(m, C);
			}
		else
		{ // ���������� ����� �� �������������� �������������� �� ����� H
			double *t = new double[n];
			t = Gauss(MatIntToDouble(Transposition(), n), VecIntToDouble(C, n), n);

			// ��������� ��� t[i] ����������� (0,1]
			for (int i = 0; i < n; i++)
				if (t[i] <= 0 || t[i] > 1) {
					delete[] t;
					return;
				}

			delete[] t;

			// ���� ������, ������ ��������� 1-�� �������

			for (int i = 0; i < n; i++)
				cover[i] = -C[i];

			if (!CheckDeterminants()) // �������� 2-��� �������
				return;

			FindC0(); // ���� �������� �0

					  //!
			std::cout << std::endl << "A suitable cover: ( ";
			for (int i = 0; i < n; i++)
				std::cout << cover[i] << " ";
			std::cout << ") | (" << c0 << ")" << std::endl;
		}


	}

	// ���������� ����� ��������� ������ �������
	int ColomnSum(int j)
	{
		int sum = 0;
		for (int i = 0; i < n; i++)
			sum += H[i][j];

		return sum;
	}

	// ������������ ����� ����� !!! ���� ��� ������������� �����
	void EnumerationIntDot()
	{
		int det = Determinant(H, n);
		std::cout << std::endl << "Det(H) = " << det;
		std::cout << std::endl << "Number of integer dots: " << det << std::endl;

		int *C = new int[n];
		for (int i = 0; i < n; i++)
			C[i] = 0;

		// ������ �������� ����� ����� � ��������������� 
		for (int i = 0; i <= ColomnSum(0); i++) // ColomnSum(0) - ������ ��������������� �� ������ ���
		{
			C[0] = i; // ����������� ����� �� �������������� � ���������������
			EnumerationIntPointsInParal(0, C);
		}

		delete[] C;
	}
#endif

	// ������� ��������� ����� t[i]*Ht[m] ������ �� Ht[m][m]
	double Tau(int **Ht, double *t, int m)
	{
		double sum = 0;
		for (int i = m + 1; i < n; i++)
			sum += Ht[m][i] * t[i];

		return sum;
	}

	// ������������ ����� ����� ���������� ����������. ������ ��� ������������� ������������� H
	// t - ������ �������� ���������� �� ������� ����������� x
	// m - ������� ����������, ���������� ��� ��������� ��������
	// x - ������� ������ ������
	void EnumIntDot(int **Ht, double *t, int m, int *x)
	{
		double tauVal = Tau(Ht, t, m);
		int res = ceil(tauVal);
		if (ceil(tauVal) == tauVal)
			res++;

		t[m] = (res - tauVal) / Ht[m][m];

		while (t[m] <= 1)
		{
			x[m] = res;
			if (m != 0)
				EnumIntDot(Ht, t, m - 1, x);
			else
			{
				for (int i = 0; i < n; i++)
					cover[i] = -x[i];

				if (!CheckDeterminants()) // �������� 2-��� �������
				{
					std::cout << std::endl << "FAILED: CheckDeterminants()";
					return;
				}

				FindC0(); // ���� �������� �0
#ifdef TEST_MOD
				std::cout /*<< std::endl*/ << "A suitable cover: ( ";
				for (int i = 0; i < n; i++)
					std::cout << cover[i] << " ";
				std::cout << ") | (" << c0 << ")" << std::endl;
#endif

				// �������� �� �������
#ifdef TEST_MOD
				if (SimplexIsEmpty(this))
					std::cout << "Void check : Empty" << std::endl;
				else
					std::cout << "Void check : Not Empty" << std::endl;
#else
				SimplexIsEmpty(this);
#endif

				// ���������� ������
#ifdef TEST_MOD
				//std::cout << "Width: " << CalcSimplexWidth(this) << std::endl;
#else
				//CalcSimplexWidth(this);
#endif
			}

			res++;
			t[m] = (res - tauVal) / Ht[m][m];
		}
	}

	// �������� ������������� � ���������� ������� | 2-�� �������
	bool CheckDeterminants()
	{
		int **Tm;
		SMP_NEW(Tm, n, int);

		for (int l = 0; l < n; l++) // ���� ������ ������� ��������
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					if (i == l)
						Tm[i][j] = cover[j];
					else
						Tm[i][j] = H[i][j];

			if (abs(Determinant(Tm, n)) > Multip - 1)
			{
				SMP_DELETE(Tm, n);
				return false;
			}
		}

		SMP_DELETE(Tm, n);
		return true;
	}

	// cT ������� �� ����� ����� | 3-� �������

	// ������� ������ c0
	void FindC0()
	{
		double *V = Gauss(MatIntToDouble(H, n), VecIntToDouble(b, n), n); // ������� ������
		double delt = 0;

		for (int i = 0; i < n; i++)
			delt += V[i] * cover[i];

		if (ceil(delt) != delt)
			c0 = ceil(delt);
		else
			c0 = ceil(delt) + 1;
	}

	/*----------------------------------------------------*/

	// ������������� ��������� � ������
	void InitDef(int Mode = SMP_INIT_MODE_0) // Mode 0 - � �������������� ������� ��������� � ������� b, Mode 1 ���
	{
		for (int i = 0; i < n; i++)
		{
			if (Mode == SMP_INIT_MODE_0)
			{
				b[i] = 0;
				cover[i] = 0;
			}

			for (int j = 0; j < n; j++)
				if (i == j)
				{
					if (Mode == SMP_INIT_MODE_0) H[i][j] = 1;
				}
				else
					H[i][j] = 0;
		}

		c0 = 0;
	}

	// ������������� ��������� � ������, ���� ��� �������� ����� ������ �������� ������ (������� ���������)
	int IncSimp(int j) // j - ������� ���� ������� ����������� �� 1
	{
		if (H[j][j] + 1 < Multip)
		{
			H[j][j]++;
			return H[j][j];
		}

		if (j != 0)
			H[j][j] = IncSimp(j - 1);
		else
		{
			H[j][j] = 1;
			return H[j][j];
		}
	}

	// ������������� ������� b
	int IncB(int j)
	{
		if (b[j] + 1 < H[j][j])
		{
			b[j]++;
			return 0; //b[j];
		}

		if (j != 0)
			b[j] = IncB(j - 1);
		else
		{
			b[j] = 0;
			return b[j];
		}
	}

	// �������� �� ��� �(main diag) ������ �������������� ������
	bool CheckDelta(int &delta)
	{
		return (MultMainDiag() <= delta);
	}

	// ������ ���������
	void printSimplex()
	{
		std::cout << "> " << std::endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
				std::cout << H[i][j] << " ";
			std::cout << " | " << b[i] << std::endl;
		}
		std::cout << "< " << std::endl;
	}

	// ���������� ������������ ��������� ������� ���������
	int MultMainDiag()
	{
		int mult = 1;
		for (int i = 0; i < n; i++)
			mult *= H[i][i];

		return mult;
	}

	// ������� ������� ������������
	void EnumerationLowTr(int _i, int _j, int &_NumOfCone)
	{
		if (_j + 1 != _i) // �� �������?
			EnumerationLowTr(_i, _j + 1, _NumOfCone);
		else
			if (_i + 1 != n) // �� ������ ������?
				EnumerationLowTr(_i + 1, 0, _NumOfCone);

		while (H[_i][_j] + 1 < H[_i][_i])
		{
			H[_i][_j]++;
			_NumOfCone++;
#ifdef TEST_MOD
			std::cout << std::endl << std::endl << "Number of cone: " << _NumOfCone;
			std::cout << std::endl << "Determinate: " << MultMainDiag();
			printSimplex();
#else
			std::cout << "d";
#endif
			GetAllCovers(this, this->n);

			if (_j + 1 != _i) // �� �������?
				EnumerationLowTr(_i, _j + 1, _NumOfCone);
			else
				if (_i + 1 != n) // �� ������ ������?
					EnumerationLowTr(_i + 1, 0, _NumOfCone);
		}

		H[_i][_j] = 0;
	}
};

/*----------------------------------------------------*/
// �������������� �������

// ��������� ������� ��� i-� ������ � j-�� ������� | ��� ������� Determinant()
void GetMatr(int **mas, int **p, int i, int j, int m) {
	int ki, kj, di, dj;
	di = 0;
	for (ki = 0; ki < m - 1; ki++) { // �������� ������� ������
		if (ki == i) di = 1;
		dj = 0;
		for (kj = 0; kj < m - 1; kj++) { // �������� ������� �������
			if (kj == j) dj = 1;
			p[ki][kj] = mas[ki + di][kj + dj];
		}
	}
}

// ����������� ���������� ������������
int Determinant(int **mas, int m) {
	int i, j, d, k, n;
	int **p;
	SMP_NEW(p, m, int);
	j = 0; d = 0;
	k = 1; //(-1) � ������� i
	n = m - 1;
	if (m < 1) std::cout << "The determinant cannot be computed!";

	if (m == 1) {
		d = mas[0][0];
		SMP_DELETE(p, m);
		return(d);
	}
	if (m == 2) {
		d = mas[0][0] * mas[1][1] - (mas[1][0] * mas[0][1]);
		SMP_DELETE(p, m);
		return(d);
	}
	if (m > 2) {
		for (i = 0; i<m; i++) {
			GetMatr(mas, p, i, 0, m);
			//cout << mas[i][j] << endl;
			//PrintMatr(p, n);
			d = d + k * mas[i][0] * Determinant(p, n);
			k = -k;
		}
	}
	SMP_DELETE(p, m);
	return(d);
}

// ���������� ����������������� ������� H
int** Transposition(int **H, int n)
{
	int **Tm;
	SMP_NEW(Tm, n, int);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			Tm[i][j] = H[j][i];

	return Tm;
}

// ������� ������������ ������� ��� i-�� ������ � j-��� �������
int** StrikeOut(int **Mat, int n, int sI, int sJ)
{
	int **Tm;
	SMP_NEW(Tm, n - 1, int);

	int tI = 0, tJ = 0;
	for (int i = 0; i < n; i++)
		if (i != sI) {
			tJ = 0;
			for (int j = 0; j < n; j++)
				if (j != sJ) {
					Tm[tI][tJ] = Mat[i][j];
					tJ++;
				}
			tI++;
		}

	return Tm;
}

// ���������� ������������� �������
int** AttachedMat(int **H, int n)
{
	int **Res, **Tm;
	SMP_NEW(Res, n, int);

	for(int i = 0; i < n; i++)
		for (int j = 0; j < n; j++) {
			Tm = StrikeOut(H, n, i, j);
			Res[i][j] = Determinant(Tm, n - 1);
		}

	//SMP_DELETE(Tm, n - 1);
	return Res;
}

// ������� ������� �� INT � DOUBLE
double** MatIntToDouble(int **A, int n)
{
	double **Tm;
	SMP_NEW(Tm, n, double);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			Tm[i][j] = A[i][j];

	return Tm;
}

// ������� ������� �� INT � DOUBLE
double* VecIntToDouble(int *A, int n)
{
	double *Tm = new double[n];

	for (int i = 0; i < n; i++)
		Tm[i] = (double)A[i];

	return Tm;
}

// ����� ������ ��� ������� ���
double* Gauss(double **a, double *y, int n)
{
	double *x, max;
	int k, index;
	const double eps = 0.00001;  // ��������
	x = new double[n];
	k = 0;
	while (k < n)
	{
		// ����� ������ � ������������ a[i][k]
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				index = i;
			}
		}
		// ������������ �����
		if (max < eps)
		{
			// ��� ��������� ������������ ���������
			std::cout << "������� �������� ���������� ��-�� �������� ������� ";
			std::cout << index << " ������� A" << std::endl;
			return 0;
		}
		for (int j = 0; j < n; j++)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;
		// ������������ ���������
		for (int i = k; i < n; i++)
		{
			double temp = a[i][k];
			if (abs(temp) < eps) continue; // ��� �������� ������������ ����������
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue; // ��������� �� �������� ���� �� ����
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}
	// �������� �����������
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
	}
	return x;
}

// CPLEX - MIP [ObjFunc->Trg, ExprMat<=b]
int* CPLEX_MIP(int **ExprMat, int *b, int Size, int exprSize, int *ObjFunc, Target Trg, AdditRestr Ar = MIP_NONE_REST)
{
	int *result = new int[Size];

	IloEnv env; // �����
	std::stringstream logfile;

	try {
		// lowerBound � upperBound ����������� �� X
		int upperBound = INT_MAX, lowerBound = -INT_MAX;
		if (Ar == MIP_XG0)
			lowerBound = 0;
		else if (Ar == MIP_XL0)
			upperBound = 0;

		IloNumVarArray x(env, Size, lowerBound, upperBound, ILOINT);
		IloModel model(env); // ������
		IloCplex cplex(env);
		cplex.setOut(env.getNullStream());

		cplex.setParam(IloCplex::Threads, 1);
		//cplex.setParam(IloCplex::WorkMem, 1024);
		cplex.setParam(IloCplex::TreLim, 2048);
		cplex.setParam(IloCplex::Param::Parallel, 1);

		cplex.extract(model);

		IloNumVarArray startVar(env); // ����������
		IloNumArray startVal(env); // �������� ����������

		for (int i = 0; i < Size; i++) {
			startVar.add(x[i]);
			startVal.add(0);
		}

		// ������� ��������
		IloExpr obj(env);
		for (int i = 0; i < Size; i++)
			obj += ObjFunc[i] * startVar[i];
		if (Trg == MIP_MAX)
			model.add(IloMaximize(env, obj));
		else if (Trg == MIP_MIN)
			model.add(IloMinimize(env, obj));
		obj.end();

		// �����������
		for (int i = 0; i < exprSize; i++) {
			IloExpr expr(env);
			for (int j = 0; j < Size; j++)
				expr += ExprMat[i][j] * startVar[j];
			model.add(expr <= b[i]);
			expr.end();
		}

		// ������ MIP mode
		cplex.addMIPStart(startVar, startVal);


		// Solve
		cplex.solve();
		if (cplex.getStatus() == IloAlgorithm::Optimal) {
			cplex.getValues(startVal, startVar);
#ifdef TEST_MOD
			env.out() << "Values = " << startVal << std::endl;
#endif
		}

		// ��������� ��������� �� IloNumVarArray � Int
		const IloArray<double> &a = (IloArray<double> &)startVal;
		for (int i = 0; i < Size; i++)
			result[i] = a[i];

		startVal.end();
		startVar.end();
	}
	catch (IloException& e) {
		std::cerr << "C-Exp: " << e << std::endl;
	}
	catch (...) {
		std::cerr << "Unknown Exception" << std::endl;
	}

	//const std::string str = logfile.str();
	//std::cout << str << std::endl;
	env.end();

	return result;
}

#define NOT_WORKING -1

// ������� ��������� ������ ���������
int CalcSimplexWidth(Simplex *S)
{
	// ������� ������
	double **Vertex;
	SMP_NEW(Vertex, S->n + 1, double);

	for (int i = 0; i < S->n + 1; i++)
	{
		int ** rep_i;
		SMP_NEW(rep_i, S->n, int);
		rep_i = S->ReplaceMatrixString(i);

		double **rep_d;
		SMP_NEW(rep_d, S->n, double);
		rep_d = MatIntToDouble(rep_i, S->n);

		Vertex[i] = Gauss(rep_d, VecIntToDouble(S->ReplaceVectorString(i), S->n),
			S->n);
	}

	// S->ReplaceMatrixString(iV) -> -1 -> T -> *(-1) | b = {0} | n
	// S->ReplaceMatrixString(jV) -> -1 -> T | b = {0} | n
	// w = SUMM(Bi) | b = -1 | 1

	for(int iV = 0; iV < S->n; iV++)
		for (int jV = iV + 1; jV < S->n + 1; jV++)
		{
			int **ExprMat = new int*[2 * S->n + 1];
			for (int newCount = 0; newCount < 2 * S->n + 1; newCount++)
				ExprMat[newCount] = new int[S->n];


			//CPLEX_MIP();
		}

	SMP_DELETE(Vertex, S->n + 1);

	return NOT_WORKING;
}

/*----------------------------------------------------*/
// ����� ��� CPLEX MIP
#define TEST_COMPIRE(val1, val2)            \
if (val1 == val2) {                         \
	std::cout << "[ OK ]" << std::endl;     \
}                                           \
else {                                      \
	std::cout << "[ FAILED ]" << std::endl; \
}

void RunAllTests() {

	int Delta = 4;
	int N = 2;
	Simplex S(N, Delta + 1);
	int *res;

	/*Case 00*/
	S.H[0][0] = 2; S.H[0][1] = 0; S.b[0] = 1 - 1;
	S.H[1][0] = 1; S.H[1][1] = 2; S.b[1] = 0 - 1;

	S.cover[0] = -1;
	S.cover[1] = -1;

	res = CPLEX_MIP(S.H, S.b, S.n, S.n, S.cover, MIP_MIN);

	std::cout << "Case 00: "; TEST_COMPIRE(res[0] * S.cover[0] + res[1] * S.cover[1], 1); std::cout << std::endl;

	/*Case 01*/
	S.H[0][0] = 1; S.H[0][1] = 0; S.b[0] = 1;
	S.H[1][0] = -3; S.H[1][1] = 2; S.b[1] = 3;

	S.cover[0] = -3;
	S.cover[1] = -3;

	res = CPLEX_MIP(S.H, S.b, S.n, S.n, S.cover, MIP_MIN);

	std::cout << "Case 01: "; TEST_COMPIRE(res[0] * S.cover[0] + res[1] * S.cover[1], -12); std::cout << std::endl;

	/*Case 02*/
	S.H[0][0] = 1; S.H[0][1] = 0; S.b[0] = 1;
	S.H[1][0] = -3; S.H[1][1] = 2; S.b[1] = 3;

	S.cover[0] = 1;
	S.cover[1] = 1;

	res = CPLEX_MIP(S.H, S.b, S.n, S.n, S.cover, MIP_MAX);

	std::cout << "Case 02: "; TEST_COMPIRE(res[0] * S.cover[0] + res[1] * S.cover[1], 4); std::cout << std::endl;

	/*Case 03*/
	S.H[0][0] = 2; S.H[0][1] = 3; S.b[0] = 6;
	S.H[1][0] = 2; S.H[1][1] = -3; S.b[1] = 3;

	S.cover[0] = 3;
	S.cover[1] = 1;

	res = CPLEX_MIP(S.H, S.b, S.n, S.n, S.cover, MIP_MAX);

	std::cout << "Case 03: "; TEST_COMPIRE(res[0] * S.cover[0] + res[1] * S.cover[1], 4); std::cout << std::endl;

	/*Case 04*/
	Delta = 14;
	N = 3;
	Simplex S3(N, Delta + 1);

	S3.H[0][0] = 3; S3.H[0][1] = 2; S3.H[0][2] = 8; S3.b[0] = 11;
	S3.H[1][0] = 2; S3.H[1][1] = 0; S3.H[1][2] = 1; S3.b[1] = 5;
	S3.H[2][0] = 3; S3.H[2][1] = 3; S3.H[2][2] = 1; S3.b[2] = 13;

	S3.cover[0] = 11;
	S3.cover[1] = 5;
	S3.cover[2] = 4;

	res = CPLEX_MIP(S3.H, S3.b, S3.n, S3.n, S3.cover, MIP_MAX, MIP_XG0);

	std::cout << "Case 04: "; TEST_COMPIRE(res[0] * S3.cover[0] + res[1] * S3.cover[1] + res[2] * S3.cover[2], 32); std::cout << std::endl;
}

/*----------------------------------------------------*/
// �������� �������

void GetAllCovers(Simplex *S, int N)
{
	double *t = new double[N];
	int *x = new int[N];
	for (int i = 0; i < N; i++)
	{
		t[i] = 0;
		x[i] = 0;
	}

	// �� ���� ���� ����������������� ������� H, ��� ������ ������� t � x, �������� m - ����� ���������� 
	S->EnumIntDot(Transposition(S->H, S->n), t, N - 1, x);

	delete[] t;
	delete[] x;
}

// ������� ����������� �������� �� ������� (���������� ������ ��������� ����� �����)
bool SimplexIsEmpty(Simplex *S)
{
	/*S.H[0][0] = 1; S.H[0][1] = 0; S.b[0] = 0;
	S.H[1][0] = 3; S.H[1][1] = 4; S.b[1] = 2;

	S.cover[0] = -3; S.cover[1] = -3; S.c0 = -1;*/
	
	// ����������� b �� �������
	int *b_inc = new int[S->n];
	for (int i = 0; i < S->n; i++)
		b_inc[i] = S->b[i] - 1;


	int *res = CPLEX_MIP(S->H, b_inc, S->n, S->n, S->cover, MIP_MIN);

	int res_d = 0;

	for (int i = 0; i < S->n; i++) {
		res_d += res[i] * S->cover[i];
	}

#ifdef TEST_MOD
	std::cout << "Obj func: " << res_d << std::endl;
#endif

	return (S->c0 <= res_d);
}

void GetAllSimplex(int &n, int &delta)
{
#ifdef TEST_MOD
	std::cout << "Simplex" << std::endl;
	std::cout << "N: " << n << std::endl;
	std::cout << "Delta: " << delta << std::endl;
#endif

	Simplex s(n, delta + 1);

	int NumOfCone  = 0; // ���������� �������
	int DeltaCount = 0; // ������� ��� ��������� �� ������

	while(true)
	{
		if (s.CheckDelta(delta))
		{
			for (int i = 0; i < s.MultMainDiag(); i++)
			{
#ifdef TEST_MOD
				std::cout << std::endl << std::endl << "Number of cone: " << NumOfCone;
				std::cout << std::endl << "Determinate: " << s.MultMainDiag();
				s.printSimplex();
#else
				std::cout << "d";
#endif
				GetAllCovers(&s, s.n);

				s.EnumerationLowTr(1, 0, NumOfCone);
				s.InitDef(SMP_INIT_MODE_1);

				NumOfCone++;
				s.IncB(n - 1);
			}

			DeltaCount = 0;
		}
		else
			DeltaCount++;

		if ((s.IncSimp(n - 1) == 1) || (DeltaCount == delta))
			break;
	}

	std::cout << "Number of all cones: " << NumOfCone << std::endl;
}

/*----------------------------------------------------*/

//#define SOURCE_CODE
//#define TEST_SIMPLEX_IS_EMPTY
//#define EXAMPLE_2
#ifdef SOURCE_CODE
#define SMP_TWO
#endif

int main(int argc, char** argv)
{
#ifdef SOURCE_CODE

	int MaxDimensionSize = 15; // ������������ �����������, ��� GetAllSimplex(N, Delta)

#ifdef SMP_TWO
	// ������ ���������� ������
	int Delta = 4;
	int N = 2;
	Simplex S(N, Delta + 1);

	S.H[1][0] = 3; S.H[1][1] = 4;
	S.b[1] = 2;
#else
	// ������ ���������� ������
	int Delta = 12;
	int N = 3;
	Simplex S(N, Delta + 1);

	S.H[1][0] = 2; S.H[1][1] = 3;
	S.H[2][0] = 1; S.H[2][1] = 2; S.H[2][2] = 4;
	S.b[0] = 0; S.b[1] = 2; S.b[2] = 3;
#endif

	S.printSimplex();

	unsigned int start_time;
	unsigned int end_time;

#ifdef SMP_OLD_CODE
	std::cout << std::endl << "#OLD_CODE" << std::endl;
	start_time = clock();
	S.EnumerationIntDot();
	end_time = clock();
	std::cout << std::endl << "Time: " << (double)(end_time - start_time) / 1000 << std::endl;
	std::cout << std::endl << "#OLD_CODE" << std::endl;
#endif

	start_time = clock();
	GetAllCovers(&S, N);
	end_time = clock();
	std::cout << std::endl << "Time: " << (double)(end_time - start_time) / 1000 << std::endl;

#endif

	//RunAllTests();

	Simplex *s = new Simplex(2, 3);
	s->H[0][0] = -1; s->H[0][1] = 0; s->b = 0;
	s->H[1][0] = 0; s->H[1][1] = -1; s->b = 0;
	s->cover[0] = 1; s->cover[1] = 1; s->c0 = 1;

	//std::cout << "Width: " << CalcSimplexWidth(s);

	// ------------------------------------------------------------------------------------------------

	int sN = 8, eN = 8;
	int sDelta = 3, eDelta = 3;

	double start_time, end_time;

	for (int delta = sDelta; delta <=eDelta; delta++)
		for (int n = sN; n <= eN; n++) {
			std::cout << std::endl << "-----------------------------------------------------------";
			std::cout << std::endl << "Delta: " << delta;
			std::cout << std::endl << "N: " << n;
			std::cout << std::endl << "-----------------------------";
			start_time = clock();
			GetAllSimplex(n, delta);
			end_time = clock();
			std::cout << std::endl << "-----------------------------";
			std::cout << std::endl << "Time: " << (end_time - start_time) / 1000 << std::endl;
			std::cout << std::endl << "-----------------------------------------------------------" << std::endl;
		}

	return 0;
}

#ifdef TEST_SIMPLEX_IS_EMPTY
int main()
{
	int Delta = 4;
	int N = 2;
	Simplex S(N, Delta + 1);

	if (SimplexIsEmpty(S))
		std::cout << std::endl << "IT'S WORK!";
	else
		std::cout << std::endl << "Nope";

	std::cin >> N;
	return 0;
}
#endif

#ifdef EXAMPLE_2

int main()
{
	int nAgents = 2;
	int nTasks = 2;
	int nLevels = 2;
	double seconds = 30.0;

	IloEnv env;
	IloArray<IloArray<IloNumVarArray> > x(env, nAgents);

	for (int i = 0; i < nAgents; i++) {
		x[i] = IloArray<IloNumVarArray>(env, nTasks);
		for (int j = 0; j < nTasks; j++) {
			x[i][j] = IloNumVarArray(env, nLevels, 0, 1, ILOINT);
		}
	}

	IloModel model(env);
	IloExpr obj(env);

	for (int i = 0; i < nAgents; i++)
		for (int j = 0; j < nTasks; j++)
			for (int k = 0; k < nLevels; k++)
				obj += (i + j + k + 1) * x[i][j][k];

	model.add(IloMinimize(env, obj));
	obj.end();

	IloCplex cplex(env);
	cplex.setOut(env.getNullStream());

	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::WorkMem, 1024);
	cplex.setParam(IloCplex::TreLim, 2048);
	cplex.setParam(IloCplex::Param::Parallel, 1); /* Deterministic */

	cplex.extract(model);

	IloNumVarArray startVar(env);
	IloNumArray startVal(env);

	for (int i = 0; i < nAgents; i++)
		for (int j = 0; j < nTasks; j++)
			for (int k = 0; k < nLevels; k++) {
				startVar.add(x[i][j][k]);
				startVal.add(0);
			}

	cplex.addMIPStart(startVar, startVal);

	// Solve
	cplex.solve();
	if (cplex.getStatus() == IloAlgorithm::Optimal) {
		cplex.getValues(startVal, startVar);
		env.out() << "Values = " << startVal << std::endl;
	}

	startVar.end();
	startVal.end();

	return 0;
}

#endif