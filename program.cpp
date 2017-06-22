#define _CRT_SECURE_NO_WARNINGS
#define PRECYZJA 1e-9
#include<iostream>
#include<stdio.h>
#include<fstream>
#include<windows.h>
#include<cmath>
using namespace std;
ofstream wyjscie("wynik.csv");

double* obliczResiduum(double** tablica,double* X,int N,double* wyrazyWolne)
{
	double* residuum = new double[N];
	for (int i = 0;i < N;i++)
	{
		double sum = 0;
		for (int j = 0;j < N;j++)
			sum += tablica[i][j] * X[j];
		residuum[i] = sum-wyrazyWolne[i];
	}
	return residuum;
}
double normaEuklidesowa(double* tablica,int N)
{
	double max = 0;
	for (int i = 0;i < N;i++)
		max += tablica[i] * tablica[i];
	max = sqrt(max);
	return max;
}
void metodaJacobiego(double **tablica,int N,double* wyrazWolne)
{
	cout.precision(20);
	double *wyniki = new double[N];
	double *wyniki1 = new double[N];
	double *residuum = new double[N];
	double norma;
	for (int i = 0;i < N;i++)
		wyniki[i] = 0;

	int iloscIteracji = 0;
	int czasPocz = GetTickCount();
	do
	{
		for (int i = 0;i < N;i++)
			residuum[i] = 0;
		for (int i = 0;i < N;i++)
		{
			double sum = 0;
			for (int j = 0;j < N;j++)
				if(i!=j)
					sum += tablica[i][j] * wyniki[j];
			wyniki1[i]=(-sum + wyrazWolne[i]) / tablica[i][i];
		}

		iloscIteracji++;
		residuum = obliczResiduum(tablica, wyniki1, N, wyrazWolne);
		norma = normaEuklidesowa(residuum,N);
		for (int i = 0;i < N;i++)
			wyniki[i] = wyniki1[i];
		
	} while ((norma)>PRECYZJA);
	int czasKoniec = GetTickCount();

	//wyjscie << iloscIteracji << ";";
	wyjscie << czasKoniec - czasPocz <<endl;
	delete[] wyniki;
	delete[] wyniki1;
	delete[] residuum;
}
void metodaGaussaSeidela(double **tablica, int N, double* wyrazWolne)
{
	cout.precision(20);
	double *wyniki = new double[N];
	double *wyniki1 = new double[N];
	double *residuum = new double[N];
	double norma;
	for (int i = 0;i < N;i++)
		wyniki[i] = 0;
	int iloscIteracji = 0;
	int czasPocz = GetTickCount();
	do
	{
		for (int i = 0;i < N;i++)
			residuum[i] = 0;
		for (int i = 0;i < N;i++)
		{
			double sum = 0;
			for (int j = 0;j <= i-1;j++)
				sum += tablica[i][j] * wyniki1[j];
			for (int j = i+1 ;j < N;j++)
				sum += tablica[i][j] * wyniki[j];
			wyniki1[i] = (-sum + wyrazWolne[i]) / tablica[i][i];
		}

		iloscIteracji++;
		residuum = obliczResiduum(tablica, wyniki1, N, wyrazWolne);
		norma = normaEuklidesowa(residuum, N);
		for (int i = 0;i < N;i++)
			wyniki[i] = wyniki1[i];
		

	} while ((norma)>PRECYZJA);
	int czasKoniec = GetTickCount();

	
	//wyjscie << iloscIteracji << ";";
	wyjscie << czasKoniec - czasPocz<< ";";
	delete [] wyniki;
	delete [] wyniki1 ;
	delete[] residuum;

}

void metodaBezpodrednialepsza(double **tablica, int N, double* wyrazyWolne)
{
	cout.precision(20);
	double *poczWolne = new double[N];
	for (int i = 0;i < N;i++)
		poczWolne[i] = wyrazyWolne[i];
	double **poczTab = new double*[N];
	for (int i = 0;i < N;i++)
		poczTab[i] = new double[N];
	for (int i = 0;i < N;i++)
		for (int j = 0;j < N;j++)
			poczTab[i][j] = tablica[i][j];

	double *wyrazyWolne1 = new double[N];
	wyrazyWolne1[0] = wyrazyWolne[0];
	double *wyniki = new double[N];
	double **tablica1 = new double*[N];
	for (int i = 0;i < N;i++)
		tablica1[i] = new double[N];
	for (int i = 0;i < N;i++)
		tablica1[0][i] = tablica[0][i];
	int czasPocz = GetTickCount();
	//1 etap
	for (int k = 1;k < N;k++)
	{
		for (int i = k;i < N;i++)
		{
			for (int j = k;j < N;j++)
			{	
					tablica1[i][j] = tablica[i][j] - ((tablica[i][k - 1] / tablica[k - 1][k - 1]) * tablica[k - 1][j]);
			}
			wyrazyWolne1[i] = wyrazyWolne[i] - ((tablica[i][k - 1] / tablica[k - 1][k - 1])*wyrazyWolne[k - 1]);
		}


		for (int i = 0;i < N;i++)
		{
			for (int j = 0;j < N;j++)
				tablica[i][j] = tablica1[i][j];
		}
		for (int i = 0;i < N;i++)
			wyrazyWolne[i] = wyrazyWolne1[i];
		
	}
	// 2 etap


	wyniki[N - 1] = wyrazyWolne1[N - 1] / tablica1[N - 1][N - 1];
	for (int i = N - 1;i > 0;i--)
	{
		double pom = 0;
		for (int j = i;j < N;j++)
		{
			pom += tablica1[i - 1][j] * wyniki[j];
		}
		wyniki[i - 1] = (wyrazyWolne1[i - 1] - pom) / tablica1[i - 1][i - 1];
	}
	double *residuum = new double[N];
	double norma;
	//residuum=obliczResiduum(poczTab, wyniki, N, poczWolne);
	//norma = normaEuklidesowa(residuum, N);
	//wyjscie << norma << ";";
	int czasKoniec = GetTickCount();
	wyjscie << czasKoniec - czasPocz << endl;
}
double czy_zbiezny(double** tablica,int N)
{
	double** D = new double*[N];
	double** LU = new double*[N];
	double** wynik = new double*[N];
	for (int i = 0;i < N;i++)
	{
		D[i] = new double[N];
		LU[i] = new double[N];
		wynik[i] = new double[N];
	}
	for (int i = 0;i < N;i++)
	{
		for (int j = 0; j < N;j++)
		{
			if (i == j)
			{
				D[i][j] = 1 / tablica[i][j];
				LU[i][j] = 0;
			}
			else
			{
				D[i][j] = 0;
				LU[i][j] = -tablica[i][j];
			}	
			wynik[i][j] = 0;
		}
	}
	for (int i = 0;i < N;i++)
	{
		for (int j = 0;j < N;j++)
		{
			for (int k = 0;k < N;k++)
				wynik[i][j] += D[i][k] * LU[k][j];
		}
	}
	double  pom=0,max=0;
	for (int i = 0;i < N;i++)
	{
		for (int j = 0;j < N;j++)
		{
			max += wynik[i][j] * wynik[i][j];		
		}
	}
	max = sqrt(max);
	for (int i = 0; i < N; i++)
	{
		delete[] wynik[i];
		delete[] D[i];
		delete[] LU[i];
	}

	delete[] wynik;
	delete[] D;
	delete[] LU;

	
	return max;
}
int main()
{
	int N = 999;
	double a1=6, a2=-1, a3=-1,e=0;
	//wyjscie << "a1;zbieznoscCf" << endl;
	for (double i = 10;i >1 ;i -= 0.125)
	{
		a1 = i;
		double *wyrazyWolne = new double[N];
		double **tablica = new double*[N];
		for (int i = 0;i < N;i++)
			tablica[i] = new double[N];
		wyjscie << a1 << ";";
		for (int i = 0;i < N;i++)
			for (int j = 0;j < N;j++)
				switch (abs(i - j))
				{
				case 0:
					tablica[i][j] = a1;
					break;
				case 1:
					tablica[i][j] = a2;
					break;
				case 2:
					tablica[i][j] = a3;
					break;
				default:
					tablica[i][j] = 0;
				}

		for (int i = 0;i < N;i++)
			wyrazyWolne[i] = sin(i*(e + 1) / 50);

		//metodaGaussaSeidela(tablica, N, wyrazyWolne);
		//metodaJacobiego(tablica, N, wyrazyWolne);
		//metodaBezpodrednialepsza(tablica, N, wyrazyWolne);
		wyjscie << czy_zbiezny(tablica, N);
		for (int i = 0; i < N; i++)
		{
			delete[] tablica[i];
		}

		delete[] tablica;
		delete[] wyrazyWolne;
	}
	
			//wyjscie << a1 << ";";
			//metodaGaussaSeidela(tablica, N, wyrazyWolne);
			//	wyjscie << a1 << ";";
			//	metodaJacobiego(tablica, N, wyrazyWolne);
	
	//	metodaBezpodrednialepsza(tablica, N, wyrazyWolne);
	
		wyjscie.close();
	return 0;
}
