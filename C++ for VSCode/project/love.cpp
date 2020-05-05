#include <iostream>
#include <Eigen/dense>
using namespace Eigen;
using namespace std;

int main()
{
	Matrix2f P, Q, G;
	P(0, 1) = 1;
	P(1, 0) = 1;
	P(0, 0) = 0;
	P(1, 1) = 0;

	Q = P;

	G = Q * P;

	cout << P << endl;
}