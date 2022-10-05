#pragma once

#include <cmath>
#include <cstring>
#include <exception>
#include <iostream>

class MatrixException : public std::exception {
 public:
  explicit MatrixException(const std::string& strExcep) : str(strExcep) {}
  ~MatrixException() throw() {}
  const char* what() const throw() { return str.c_str(); }

 private:
  std::string str;
};

template <class T>
class S21Matrix {
 public:
  //  CONSTRUCTORS
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix<T>& other);
  S21Matrix(S21Matrix<T>&& other);
  //  DESTRUCTOR
  ~S21Matrix();
  //  FUNCTIONS
  bool eq_matrix(const S21Matrix<T>& other);
  void sum_matrix(const S21Matrix<T>& other);
  void sub_matrix(const S21Matrix<T>& other);
  void mul_number(const T& val);
  void mul_matrix(const S21Matrix<T>& other);
  S21Matrix<T> transpose();
  void fillMatrix(const T& val);

  // OPERATOR
  S21Matrix& operator=(const S21Matrix<T>& other);
  S21Matrix& operator=(S21Matrix<T>&& other);
  S21Matrix<T>& operator+=(const S21Matrix& other);
  S21Matrix<T>& operator-=(const S21Matrix<T>& other);
  S21Matrix<T>& operator*=(const S21Matrix<T>& other);
  S21Matrix<T>& operator*=(const T& val);
  S21Matrix<T> operator+(const S21Matrix<T>& other);
  S21Matrix<T> operator-(const S21Matrix<T>& other);
  S21Matrix<T> operator*(const T& val);
  S21Matrix<T> operator*(const S21Matrix<T>& other);
  bool operator==(const S21Matrix<T>& other);
  T& operator()(int row, int col);
  T operator()(int row, int col) const;

  // GET SET
  int get_rows();
  int get_rows() const;
  void set_rows(int rows);
  int get_cols();
  int get_cols() const;
  void set_cols(int cols);

 private:
  // HELPER FOO
  void memoryAllocation();
  void memoryClean();
  void sumOrSub(const S21Matrix<T>& other, int sign);

  int _rows = 0, _cols = 0;
  T** _matrix = nullptr;
};

#include "s21_matrix_oop.inl"
