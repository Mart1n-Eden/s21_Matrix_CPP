#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

TEST(Test, Constructor) {
  S21Matrix mat;
  ASSERT_EQ(mat.GetRows(), 0);
  ASSERT_EQ(mat.GetCols(), 0);
}

TEST(Test, ConstructorError) {
  EXPECT_THROW(S21Matrix mat(-1, -1), std::length_error);
}

TEST(Test, MoveConstructor) {
  S21Matrix A(7, 7);
  S21Matrix B(std::move(A));
  EXPECT_EQ(7, B.GetRows());
  EXPECT_EQ(7, B.GetCols());
  EXPECT_EQ(0, A.GetRows());
  EXPECT_EQ(0, A.GetCols());
}

TEST(Test, CopyConstructor) {
  S21Matrix A(4, 3);
  S21Matrix B(A);
  EXPECT_EQ(4, B.GetRows());
  EXPECT_EQ(3, B.GetCols());
  EXPECT_EQ(1, A == B);
}

TEST(Test, SumMatrixTest) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  S21Matrix res(2, 2);

  a.InitMatrix();
  b.InitMatrix();
  res(0, 0) = 2;
  res(0, 1) = 2;
  res(1, 0) = 2;
  res(1, 1) = 4;

  a.SumMatrix(b);
  EXPECT_TRUE(a == res);
}

TEST(Test, SumMatrixError) {
  S21Matrix a;
  S21Matrix b;
  EXPECT_THROW(a.SumMatrix(b), std::length_error);
}

TEST(Test, SubMatrixTest) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  S21Matrix res(2, 2);

  a.InitMatrix();
  b.InitMatrix();

  a.SubMatrix(b);
  EXPECT_TRUE(a == res);
}

TEST(Test, SubMatrixError) {
  S21Matrix a;
  S21Matrix b;
  EXPECT_THROW(a.SubMatrix(b), std::length_error);
}

TEST(Test, MulMatrixTest) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  S21Matrix res(2, 2);

  a.InitMatrix();
  b.InitMatrix();
  res(0, 0) = 2;
  res(0, 1) = 3;
  res(1, 0) = 3;
  res(1, 1) = 5;

  a.MulMatrix(b);

  EXPECT_TRUE(a == res);
}

TEST(Test, MulMatrixError) {
  S21Matrix a(1, 2);
  S21Matrix b(3, 4);
  EXPECT_THROW(a.MulMatrix(b), std::length_error);
}

TEST(Test, MulNumberTest) {
  S21Matrix a(2, 2);
  S21Matrix res(2, 2);

  a.InitMatrix();
  a.MulNumber(2);

  res(0, 0) = 2;
  res(0, 1) = 2;
  res(1, 0) = 2;
  res(1, 1) = 4;

  EXPECT_TRUE(a == res);
}

TEST(Test, EqualTest) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);

  a.InitMatrix();
  b.InitMatrix();

  EXPECT_TRUE(a.EqMatrix(b));
}

TEST(Test, EQmatrix) {
  S21Matrix a(3, 3);
  S21Matrix b(2, 2);
  EXPECT_FALSE(a.EqMatrix(b));
}

TEST(Test, OperatorPlus) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);
  b.InitMatrix();

  a = a + b;

  EXPECT_EQ(1, a == b);
}

TEST(Test, OperatorMinus) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);

  a.InitMatrix();
  b.InitMatrix();
  a.MulNumber(2);
  a = a - b;

  EXPECT_EQ(1, a == b);
}

TEST(Test, OperatorMulMatrix) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  S21Matrix res(2, 2);

  a.InitMatrix();
  b.InitMatrix();
  res(0, 0) = 2;
  res(0, 1) = 3;
  res(1, 0) = 3;
  res(1, 1) = 5;

  a = a * b;
  EXPECT_EQ(1, a == res);
}

TEST(Test, OperatorMulNumber) {
  S21Matrix a(2, 2);
  S21Matrix res(2, 2);

  a.InitMatrix();
  a = a * 2;

  res(0, 0) = 2;
  res(0, 1) = 2;
  res(1, 0) = 2;
  res(1, 1) = 4;

  EXPECT_TRUE(a == res);
}

TEST(Test, OperatorPlusEq) {
  S21Matrix a(3, 3);
  S21Matrix b(3, 3);

  b.InitMatrix();
  a += b;

  EXPECT_EQ(1, a == b);
}

TEST(Test, OperatorMinusEq) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);

  a.InitMatrix();
  b.InitMatrix();
  a.MulNumber(2);
  a -= b;

  EXPECT_EQ(1, a == b);
}

TEST(Test, OperatorMulNumberEq) {
  S21Matrix a(2, 2);
  S21Matrix res(2, 2);

  a.InitMatrix();
  a *= 2;

  res(0, 0) = 2;
  res(0, 1) = 2;
  res(1, 0) = 2;
  res(1, 1) = 4;

  EXPECT_TRUE(a == res);
}

TEST(Test, OperatorNumberMulEq) {
  S21Matrix a(2, 2), b(2, 2);
  S21Matrix res(2, 2);

  a.InitMatrix();
  b = 2 * a;

  res(0, 0) = 2;
  res(0, 1) = 2;
  res(1, 0) = 2;
  res(1, 1) = 4;

  EXPECT_TRUE(b == res);
}

TEST(Test, OperatorMulMatrixEq) {
  S21Matrix a(2, 2);
  S21Matrix b(2, 2);
  S21Matrix res(2, 2);

  a.InitMatrix();
  b.InitMatrix();
  res(0, 0) = 2;
  res(0, 1) = 3;
  res(1, 0) = 3;
  res(1, 1) = 5;

  a *= b;
  EXPECT_EQ(1, a == res);
}

TEST(Test, InverseMatrixTest) {
  S21Matrix a(3, 3), b(3, 3), res(3, 3);

  a.InitMatrix();
  b = a.InverseMatrix();
  res(0, 0) = 1;
  res(0, 1) = -2;
  res(0, 2) = 1;
  res(1, 0) = -2;
  res(1, 1) = 0;
  res(1, 2) = 1;
  res(2, 0) = 1;
  res(2, 1) = 1;
  res(2, 2) = -1;

  EXPECT_EQ(1, res == b);
}

TEST(Test, InverseMatrixError) {
  S21Matrix a(1, 2);
  EXPECT_THROW(a.InverseMatrix(), std::length_error);
}

TEST(Test, DeterminantTest) {
  S21Matrix a(2, 2);
  double x = 0;

  a.InitMatrix();
  x = a.Determinant();

  EXPECT_EQ(1, x);
}

TEST(Test, DeterminantError) {
  S21Matrix a(1, 2);
  EXPECT_THROW(a.Determinant(), std::length_error);
}

TEST(Test, TransposeTest) {
  S21Matrix a(2, 2), b(2, 2), res(2, 2);

  a.InitMatrix();
  a(0, 0) = 3;
  a(0, 1) = 3;
  b = a.Transpose();
  res(0, 0) = 3;
  res(0, 1) = 1;
  res(1, 0) = 3;
  res(1, 1) = 2;

  EXPECT_EQ(1, res == b);
}

TEST(Test, TransposeError) {
  S21Matrix a;
  EXPECT_THROW(a.Transpose(), std::length_error);
}

TEST(Test, CalcComplementsTest) {
  S21Matrix a(3, 3), b(3, 3), res(3, 3);

  a.InitMatrix();
  b = a.CalcComplements();
  res(0, 0) = -1;
  res(0, 1) = 2;
  res(0, 2) = -1;
  res(1, 0) = 2;
  res(1, 1) = 0;
  res(1, 2) = -1;
  res(2, 0) = -1;
  res(2, 1) = -1;
  res(2, 2) = 1;

  EXPECT_EQ(1, res == b);
}

TEST(Test, CalcComplementsError) {
  S21Matrix a;
  EXPECT_THROW(a.CalcComplements(), std::length_error);
}

TEST(Test, SetCols) {
  S21Matrix a(1, 1);
  S21Matrix b(1, 2);
  a.SetCols(2);
  EXPECT_EQ(a.GetCols(), b.GetCols());
}

TEST(Test, SetColumnsError) {
  S21Matrix a;
  EXPECT_THROW(a.SetCols(0), std::length_error);
}

TEST(Test, SetRows) {
  S21Matrix a(1, 1);
  S21Matrix b(2, 1);
  a.SetRows(2);
  EXPECT_EQ(a.GetRows(), b.GetRows());
}

TEST(Test, SetRowsError) {
  S21Matrix a;
  EXPECT_THROW(a.SetRows(0), std::length_error);
}

TEST(Test, OperatorBracketsError) {
  S21Matrix b(2, 2);
  EXPECT_THROW((b(0, 2) = 2), std::length_error);
}

TEST(Test, OperatorBracketsConst) {
  const S21Matrix a(1, 1);
  a(0, 0) = 1;
  EXPECT_EQ(a(0, 0), 1);
}

TEST(Test, OperatorBracketsConstError) {
  const S21Matrix a(1, 1);
  EXPECT_THROW((a(0, 2) = 1), std::length_error);
}

int main(int argc, char *argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}