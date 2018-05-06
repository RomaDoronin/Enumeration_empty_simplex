#include <iostream>
#include <ctime>

/*----------------------------------------------------*/
// Директивы

//#define DEBUG
#define OLD_CODE
#define SMP_INIT_MODE_0 0
#define SMP_INIT_MODE_1 1

/*----------------------------------------------------*/
// Макросы

#define SMP_NEW(mas, n, type) mas = new type*[n]; \
for (int newCount = 0; newCount < n; newCount++) \
	mas[newCount] = new type[n]

#define SMP_DELETE(mas, n) for (int deleteCount = 0; deleteCount < n; deleteCount++) \
delete[] mas[deleteCount]

/*----------------------------------------------------*/
// Дополнительные функции

// Получение матрицы без i-й строки и j-го столбца | Для функции Determinant()
void GetMatr(int **mas, int **p, int i, int j, int m) {
	int ki, kj, di, dj;
	di = 0;
	for (ki = 0; ki < m - 1; ki++) { // Проверка индекса строки
		if (ki == i) di = 1;
		dj = 0;
		for (kj = 0; kj < m - 1; kj++) { // Проверка индекса столбца
			if (kj == j) dj = 1;
			p[ki][kj] = mas[ki + di][kj + dj];
		}
	}
}

// Рекурсивное вычисление определителя
int Determinant(int **mas, int m) {
	int i, j, d, k, n;
	int **p;
	SMP_NEW(p, m, int);
	j = 0; d = 0;
	k = 1; //(-1) в степени i
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

// Перевод матрицы из INT в DOUBLE
double** MatIntToDouble(int **A, int n)
{
	double **Tm;
	SMP_NEW(Tm, n, double);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			Tm[i][j] = A[i][j];

	return Tm;
}

// Перевод вектора из INT в DOUBLE
double* VecIntToDouble(int *A, int n)
{
	double *Tm = new double[n];

	for (int i = 0; i < n; i++)
		Tm[i] = (double)A[i];

	return Tm;
}

// Метод Гаусса для решения СЛУ
double* Gauss(double **a, double *y, int n)
{
	double *x, max;
	int k, index;
	const double eps = 0.00001;  // точность
	x = new double[n];
	k = 0;
	while (k < n)
	{
		// Поиск строки с максимальным a[i][k]
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
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			std::cout << "Решение получить невозможно из-за нулевого столбца ";
			std::cout << index << " матрицы A" << std::endl;
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
		// Нормализация уравнений
		for (int i = k; i < n; i++)
		{
			double temp = a[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
	}
	return x;
}

/*----------------------------------------------------*/
// Класс симплекса [ Ax <= b ] -> [ AQx <= b ] -> [ Hx <= b ] (Форма Эрмита)
class Simplex {
public:
	int **H;         // Матрица H
	int *b;          // Вектор b

	int Multip;      // Кратность
	int n;           // Размерность

	// Крышка с (cTx <= c0)
	int *cover;
	int c0;

public:
	// Конструктор
	Simplex(int _n, int _M) : n(_n), Multip(_M)
	{
		b = new int[n];
		cover = new int[n];
		SMP_NEW(H, n, int);

		InitDef();
	}

	// Деструктор
	~Simplex()
	{
		delete[] b;
		delete[] cover;
		SMP_DELETE(H, n);
	}

	// Возвращает транспонированную матрицу H
	int** Transposition()
	{
		int **Tm;
		SMP_NEW(Tm, n, int);

		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				Tm[i][j] = H[j][i];

		return Tm;
	}

	/*----------------------------------------------------*/
	// Методы для работы с крышкой

#ifdef OLD_CODE
	// Рекурсивный перебор целых точек в n-мерном параллелепипиде, составленный из строк матрицы H
	void EnumerationIntPointsInParal(int m, int *C)
	{
		m++; // Параметр отвечающий за остановку рекурсии, по сути считает чтобы мы мы не вышли за n 
		// Перебор всех целых точек
		if (m != n)
			for (int i = 0; i <= ColomnSum(m); i++)
			{
				C[m] = i;
				EnumerationIntPointsInParal(m, C);
			}
		else
		{ // Провеверка точек на пренадлежность паралелограмму из строк H
			double *t = new double[n];
			t = Gauss(MatIntToDouble(Transposition(), n), VecIntToDouble(C, n), n);

			// Проверяем что t[i] принадлежат (0,1]
			for (int i = 0; i < n; i++)
				if (t[i] <= 0 || t[i] > 1) {
					delete[] t;
					return;
				}

			delete[] t;
			
			// Если прошло, значит выполнено 1-ое УСЛОВИЕ

			for (int i = 0; i < n; i++)
				cover[i] = -C[i];

			if (!CheckDeterminants()) // Проверка 2-ого УСЛОВИЯ
				return;

			FindC0(); // Поск элемента с0

			//!
			std::cout << std::endl << "A suitable cover: ( ";
			for (int i = 0; i < n; i++)
				std::cout << cover[i] << " ";
			std::cout << ") | (" << c0 << ")" << std::endl;
		}


	}

	// Возвращает сумму элементов одного столбца
	int ColomnSum(int j)
	{
		int sum = 0;
		for (int i = 0; i < n; i++)
			sum += H[i][j];

		return sum;
	}

	// Перечисление целых точек !!! Пока без отрицательных чисел
	void EnumerationIntDot()
	{
		int det = Determinant(H, n);
		std::cout << std::endl << "Det(H) = " << det;
		std::cout << std::endl << "Number of integer dots: " << det << std::endl;

		int *C = new int[n];
		for (int i = 0; i < n; i++)
			C[i] = 0;

		// Начало перебора целых точек в параллелепипеде 
		for (int i = 0; i <= ColomnSum(0); i++) // ColomnSum(0) - Высота параллелепипеда по первой оси
		{
			C[0] = i; // Проверяемая точка на принадлежность к параллелепипеду
			EnumerationIntPointsInParal(0, C);
		}

		delete[] C;
	}
#endif

	// Функция считающая сумму t[i]*Ht[m] вправа от Ht[m][m]
	double Tau(int **Ht, double *t, int m)
	{
		double sum = 0;
		for (int i = m + 1; i < n; i++)
			sum += Ht[m][i] * t[i];

		return sum;
	}

	// Перечисление целых точек ПРАВИЛЬНАЯ реализация. Только для положительных коэффициентов H
	// t - вектор условных переменных по которым проверяются x
	// m - степень погружения, переменная для остановки рекурсии
	// x - искомый вектор крышек
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

				if (!CheckDeterminants()) // Проверка 2-ого УСЛОВИЯ
				{
					std::cout << std::endl << "FAILED: CheckDeterminants()";
					return;
				}

				FindC0(); // Поск элемента с0

				std::cout << std::endl << "A suitable cover: ( ";
				for (int i = 0; i < n; i++)
					std::cout << cover[i] << " ";
				std::cout << ") | (" << c0 << ")" << std::endl;
			}

			res++;
			t[m] = (res - tauVal) / Ht[m][m];
		}
	}

	// Проверка определителей с замененной строкой | 2-ое УСЛОВИЕ
	bool CheckDeterminants()
	{
		int **Tm;
		SMP_NEW(Tm, n, int);

		for (int l = 0; l < n; l++) // Цикл строки которую заменяем
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

	// cT состоит из целых чисел | 3-е УСЛОВИЕ

	// Функция поиска c0
	void FindC0()
	{
		double *V = Gauss(MatIntToDouble(H, n), VecIntToDouble(b, n), n); // Вершина конуса
		double delt = 0;

		for (int i = 0; i < n; i++)
			delt += V[i] * cover[i];	

		if (ceil(delt) != delt)
			c0 = ceil(delt);
		else
			c0 = ceil(delt) + 1;
	}

	/*----------------------------------------------------*/

	// Инициализация еденицами и нулями
	void InitDef(int Mode = SMP_INIT_MODE_0) // Mode 0 - с инициализацией главной диагонали и вектора b, Mode 1 без
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

	// Инкрементация симплекса с учетом, того что значения слева больше значений справа (Главная диагональ)
	int IncSimp(int j) // j - позиция бита который инкрементим на 1
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

	// Инкрементация вектора b
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

	// Проверка то что П(main diag) меньше установленного дельта
	bool CheckDelta(int delta)
	{
		return (MultMainDiag() <= delta);
	}

	// Печать симплекса
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

	// Возвращает произведение элементов главной диагонали
	int MultMainDiag()
	{
		int mult = 1;
		for (int i = 0; i < n; i++)
			mult *= H[i][i];

		return mult;
	}

	// Перебор нижнего треугольника
	void EnumerationLowTr(int _i, int _j, int &_NumOfCone)
	{
		if (_j + 1 != _i) // Не крайняя?
			EnumerationLowTr(_i, _j + 1, _NumOfCone);
		else
			if (_i + 1 != n) // Не нижняя строка?
				EnumerationLowTr(_i + 1, 0, _NumOfCone);

		while (H[_i][_j] + 1 < H[_i][_i])
		{
			H[_i][_j]++;
			printSimplex();
			_NumOfCone++;

			if (_j + 1 != _i) // Не крайняя?
				EnumerationLowTr(_i, _j + 1, _NumOfCone);
			else
				if (_i + 1 != n) // Не нижняя строка?
					EnumerationLowTr(_i + 1, 0, _NumOfCone);
		}

		H[_i][_j] = 0;
	}
};

/*----------------------------------------------------*/
// Тестовая функции

void GetAllSimplex(int n, int delta)
{
	std::cout << "Simplex" << std::endl;
	std::cout << "N: " << n << std::endl;
	std::cout << "Delta: " << delta << std::endl;

	Simplex s(n, delta + 1);

	int NumOfCone  = 0; // Количество конусов
	int DeltaCount = 0; // Счетчик для остановки по дельте
#ifdef DEBUG
	int Complexity = 0; // Сложность
#endif

	while(true)
	{
#ifdef DEBUG
		Complexity++;
#endif
		if (s.CheckDelta(delta))
		{
			for (int i = 0; i < s.MultMainDiag(); i++)
			{
				s.printSimplex();

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

#ifdef DEBUG
	std::cout << "Complexity: " << Complexity << std::endl;
#endif
	std::cout << "Number of cone: " << NumOfCone << std::endl;
}

void GetAllCovers(Simplex *S, int N)
{
	double *t = new double[N];
	int *x = new int[N];
	for (int i = 0; i < N; i++)
	{
		t[i] = 0;
		x[i] = 0;
	}

	// На вход идут транспонированная матрица H, два пустых вектора t и x, значение m - точка погружения 
	S->EnumIntDot(S->Transposition(), t, N - 1, x);

	delete[] t;
	delete[] x;
}

/*----------------------------------------------------*/

//#define TWO

int main(int argc, char** argv)
{
	//int MaxDimensionSize = 15; // Максивальная размерность

#ifdef TWO
	// Пример двухмерной задачи
	int Delta = 4;
	int N = 2;
	Simplex S(N, Delta + 1);

	S.H[1][0] = 3; S.H[1][1] = 4;
	S.b[1] = 2;
#else
	// Пример трехмерной задачи
	int Delta = 9;
	int N = 9;
	Simplex S(N, Delta + 1);

	S.H[1][0] = 1; S.H[1][1] = 2;
	S.H[2][0] = 3; S.H[2][1] = 3; S.H[2][2] = 4;
	S.b[1] = 1; S.b[2] = 1;
#endif

	S.printSimplex();

#ifdef OLD_CODE
	std::cout << std::endl << "#OLD_CODE" << std::endl;
	unsigned int start_time = clock();
	S.EnumerationIntDot();
	unsigned int end_time = clock();
	std::cout << std::endl << "Time: " << (double)(end_time - start_time) / 1000 << std::endl;
	std::cout << std::endl << "#OLD_CODE" << std::endl;
#endif

	start_time = clock();
	GetAllCovers(&S, N);
	end_time = clock();
	std::cout << std::endl << "Time: " << (double)(end_time - start_time) / 1000 << std::endl;

	std::cout << std::endl;
	//GetAllSimplex(N, Delta);

	return 0;
}