#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() {  // Базовый конструктор матрицы 3х3
  rows_ = 0;
  cols_ = 0;
  matrix_ = nullptr;
  // this->NewMatrix(rows_, cols_);
}

S21Matrix::S21Matrix(int rows, int cols)
    : rows_(rows), cols_(cols) {  // Параметризированный конструктор с
                                  // количеством строк и столбцов
  if (rows < 0 || cols < 0)
    throw std::length_error("Matrix parameters are incorrectly specified");
  this->NewMatrix(rows, cols);
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : S21Matrix(other.rows_, other.cols_) {  // Конструктор копирования
  this->CopyMatrix(other.matrix_);
}

S21Matrix::S21Matrix(S21Matrix &&other)
    : rows_(other.rows_),
      cols_(other.cols_),
      matrix_(other.matrix_) {  // Конструктор переноса
  other.matrix_ = nullptr;
  other.rows_ = 0;
  other.cols_ = 0;
}

S21Matrix::~S21Matrix() {  // Деструктор
  this->FreeMatrix();
}

bool S21Matrix::EqMatrix(
    const S21Matrix &other) const noexcept {  // Сравнение матриц
  bool res = 1;
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    // throw std::length_error("Different matrix dimensions");
    res = 0;
  } else {
    for (int i = 0; i < rows_ && res; i++) {
      for (int j = 0; j < cols_ && res; j++) {
        res = fabs(matrix_[i][j] - other.matrix_[i][j]) < 1e-7 ? res : 0;
      }
    }
  }
  return res;
}

void S21Matrix::SumMatrix(const S21Matrix &other) const {  // Cложение матриц
  if (rows_ != other.rows_ || cols_ != other.cols_ || rows_ == 0 ||
      other.rows_ == 0 || cols_ == 0 || other.cols_ == 0) {
    throw std::length_error("Different matrix dimensions");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) const {  // Вычетание матриц
  if (rows_ != other.rows_ || cols_ != other.cols_ || rows_ == 0 ||
      other.rows_ == 0 || cols_ == 0 || other.cols_ == 0) {
    throw std::length_error("Different matrix dimensions");
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(
    const double num) const noexcept {  // Умножение матрицы на число
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {  // Перемножение матриц
  if (cols_ != other.rows_) {
    throw std::length_error(
        "The dimensionality of matrix is not suitable for this operation");
  }
  S21Matrix tmp(rows_, other.cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      for (int k = 0; k < cols_; k++) {
        tmp.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  cols_ = other.cols_;
  std::swap(matrix_, tmp.matrix_);
}

S21Matrix S21Matrix::Transpose() const {  // Транспонирование матриц
  if (rows_ == 0 || cols_ == 0) {
    throw std::length_error(
        "The dimensionality of matrix is not suitable for this operation");
  }
  S21Matrix res(rows_, cols_);
  res.rows_ = rows_;
  res.cols_ = cols_;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      res.matrix_[j][i] = matrix_[i][j];
    }
  }
  return res;
}

double S21Matrix::Determinant() const {  // Вычисление определителя
  if (rows_ != cols_ || rows_ == 0 || cols_ == 0) {
    throw std::length_error("This matrix is not square");
  }
  double res = CalcDeterminant(*this);
  return res;
}

S21Matrix S21Matrix::CalcComplements() const {  // Вычисление матрицы
  // алгебраических дополнений
  if (rows_ != cols_ || rows_ == 0 || cols_ == 0) {
    throw std::length_error("This matrix is not square");
  }
  S21Matrix res(rows_, cols_);
  if (rows_ == 1) {
    res.matrix_[0][0] = matrix_[0][0];
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        S21Matrix tmp(rows_ - 1, cols_ - 1);
        tmp.Separate(*this, i, j);
        res.matrix_[i][j] = pow(-1, i + j) * tmp.Determinant();
      }
    }
  }
  return res;
}

S21Matrix S21Matrix::InverseMatrix() const {  // Вычисление обратной матрицы
  if (cols_ != rows_ || rows_ == 0 || cols_ == 0) {
    throw std::length_error("This matrix is not square");
  }
  S21Matrix res(rows_, cols_);
  if (rows_ == 1) {
    res.matrix_[0][0] = 1.0 / matrix_[0][0];
  } else {
    double x = this->Determinant();
    if (std::fabs(x) < 1e-7) {
      throw std::length_error("Determinant is less than zero");
    }
    res = this->Transpose().CalcComplements();
    res.MulNumber(1 / x);
  }
  return res;
}

void S21Matrix::Separate(const S21Matrix &m, int x, int y) noexcept {
  for (int i = 0, r = 0; r < m.rows_; r++) {
    if (r == x) continue;
    for (int j = 0, c = 0; c < m.cols_; c++) {
      if (c == y) continue;
      matrix_[i][j] = m.matrix_[r][c];
      j++;
    }
    i++;
  }
}

double CalcDeterminant(const S21Matrix &m) {
  double res = 0;
  if (m.rows_ == 1) {
    res = m.matrix_[0][0];
  } else {
    S21Matrix tmp(m.rows_ - 1, m.cols_ - 1);
    for (int i = 0; i < m.rows_; i++) {
      tmp.Separate(m, 0, i);
      res += pow(-1, i) * m.matrix_[0][i] * CalcDeterminant(tmp);
    }
  }
  return res;
}

S21Matrix operator*(const double num, const S21Matrix &other) {
  return other * num;
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) const {
  S21Matrix res(other);
  res.SumMatrix(*this);
  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) const {
  S21Matrix res(*this);
  res.SubMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) const {
  S21Matrix res(*this);
  res.MulMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const double num) const {
  S21Matrix res(*this);
  res.MulNumber(num);
  return res;
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  this->FreeMatrix();
  rows_ = other.rows_;
  cols_ = other.cols_;
  this->NewMatrix(rows_, cols_);
  this->CopyMatrix(other.matrix_);
  return *this;
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  this->SumMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  this->SubMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  this->MulMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const double num) {
  this->MulNumber(num);
  return *this;
}

bool S21Matrix::operator==(const S21Matrix &other) const noexcept {
  return this->EqMatrix(other);
}

double &S21Matrix::operator()(int row, int col) const {
  if (row < 0 || col < 0 || row > rows_ - 1 || col > cols_ - 1) {
    throw std::length_error("There is no such matrix element");
  }
  return matrix_[row][col];
}

double &S21Matrix::operator()(int row, int col) {
  if (row < 0 || col < 0 || row > rows_ - 1 || col > cols_ - 1) {
    throw std::length_error("There is no such matrix element");
  }
  return matrix_[row][col];
}

int S21Matrix::GetRows() const noexcept { return rows_; }

int S21Matrix::GetCols() const noexcept { return cols_; }

void S21Matrix::SetRows(int rows) {
  if (rows <= 0) {
    throw std::length_error("Are You Normal?");
  }
  if (rows_ != rows) {
    S21Matrix tmp(rows, cols_);
    for (int i = 0; i < rows_ && i < rows; i++) {
      for (int j = 0; j < cols_; j++) {
        tmp.matrix_[i][j] = matrix_[i][j];
      }
    }
    *this = tmp;
  }
}

void S21Matrix::SetCols(int cols) {
  if (cols <= 0) {
    throw std::length_error("Are You Normal?");
  }
  if (cols_ != cols) {
    S21Matrix tmp(rows_, cols);
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_ && j < cols; j++) {
        tmp.matrix_[i][j] = matrix_[i][j];
      }
    }
    *this = tmp;
  }
}

void S21Matrix::CopyMatrix(double **matr) noexcept {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = matr[i][j];
    }
  }
}

void S21Matrix::NewMatrix(int rows, int cols) noexcept {
  matrix_ = new double *[rows]();
  for (int i = 0; i < rows; i++) {
    matrix_[i] = new double[cols]();
  }
}

void S21Matrix::FreeMatrix() noexcept {
  if (matrix_) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete matrix_;
  }
  matrix_ = nullptr;
  rows_ = 0;
  cols_ = 0;
}

void S21Matrix::InitMatrix() noexcept {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = i + j;
    }
  }
  matrix_[0][0] = 1;
  // matrix_[0][1] = 1;
  // matrix_[0][2] = 2;
  // matrix_[1][0] = 3;
  // matrix_[1][1] = 4;
  // matrix_[1][2] = 5;
  // matrix_[2][0] = 6;
  // matrix_[2][1] = 7;
  // matrix_[2][2] = 8;
}

// void S21Matrix::PrintMatrix() {
//   for (int i = 0; i < rows_; i++) {
//     for (int j = 0; j < cols_; j++)
//       std::cout << matrix_[i][j] << " ";
//     std::cout << std::endl;
//   }
//   std::cout << std::endl;
// }