#ifndef __S21MATRIX_H__
#define __S21MATRIX_H__

#include <cmath>
#include <iostream>

class S21Matrix {
  friend double CalcDeterminant(const S21Matrix &);
  friend S21Matrix operator*(const double num, const S21Matrix &other);

 private:
  int rows_, cols_;
  double **matrix_;

  void Separate(const S21Matrix &m, int x, int y) noexcept;
  void CopyMatrix(double **matr) noexcept;
  void NewMatrix(int rows, int columns) noexcept;
  void FreeMatrix() noexcept;

 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix &other);
  S21Matrix(S21Matrix &&other);
  ~S21Matrix();

  S21Matrix operator+(const S21Matrix &other) const;
  S21Matrix operator-(const S21Matrix &other) const;
  S21Matrix operator*(const S21Matrix &other) const;
  S21Matrix operator*(const double num) const;
  S21Matrix &operator=(const S21Matrix &other);
  S21Matrix &operator+=(const S21Matrix &other);
  S21Matrix &operator-=(const S21Matrix &other);
  S21Matrix &operator*=(const S21Matrix &other);
  S21Matrix &operator*=(const double num);
  bool operator==(const S21Matrix &other) const noexcept;
  double &operator()(int row, int col);
  double &operator()(int row, int col) const;

  void SumMatrix(const S21Matrix &other) const;
  bool EqMatrix(const S21Matrix &other) const noexcept;
  void SubMatrix(const S21Matrix &other) const;
  void MulNumber(const double num) const noexcept;
  void MulMatrix(const S21Matrix &other);
  double Determinant() const;
  S21Matrix Transpose() const;
  S21Matrix CalcComplements() const;
  S21Matrix InverseMatrix() const;

  int GetRows() const noexcept;
  int GetCols() const noexcept;
  void SetRows(int rows);
  void SetCols(int cols);

  // void PrintMatrix();
  void InitMatrix() noexcept;
};

#endif