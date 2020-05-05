#ifndef MATRIX_H_
#define MATRIX_H_
#include <iostream>
#define NUM 1000

namespace MATRIX
{
    class Matrix
    {
        private:
        double M[NUM][NUM]={0};
        int row;
        int column;

        public:
            Matrix();
            Matrix(const Matrix &A);
            Matrix(const double **N, int r, int c);
            double get(int r, int c) { return M[r][c]; };
            void set_position(int r, int c);
            void set_Num(int r, int c, double n);

            void operator=(const Matrix &A);
            Matrix operator+(const Matrix &A) const;
            Matrix operator-(const Matrix &A) const;
            Matrix operator*(const Matrix &A) const;
            Matrix operator*(const double n) const;

            void swag(int r1,int r2,char ch='c');
            void Primary_Matrix(int r, int c, double k,char ch='c');
            friend Matrix Matrix_multiply(const double **M, const double **N, int r, int c, int c1);
            friend void swag(double &a,double &b);

    };

} // namespace MATRIX
#endif