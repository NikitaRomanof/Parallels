#pragma once

#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <thread>
#include <mutex>

#include "../s21_matrix_oop.h"

class Gauss {
 public:
  Gauss();
  explicit Gauss(S21Matrix<double> mat);
  Gauss(S21Matrix<double> mat, S21Matrix<double> b);
  Gauss(const Gauss &other) = delete;
  Gauss(Gauss &&other) = delete;
  ~Gauss();
  Gauss &operator=(const Gauss &other) = delete;
  void GetLinearSystemFromFile(const std::string &filename);
  void GenerateLinearSystem(int size);
  S21Matrix<double> *SolveLinearSystemSingleThread();
  S21Matrix<double> *SolveLinearSystemMultiThread();
  void Clear();
  void PrintMatrix(const std::string &label);
  int GetSize();
  bool IsEmpty();

 private:
  int size_;
  bool multi_thread_;
  S21Matrix<double> *matrix_;
  S21Matrix<double> *b_matrix_;
  S21Matrix<double> *l_matrix_;  // матрица коэффициентов домножения
  S21Matrix<double> *result_;  // матрица (вектор) результатов

  // methods
  // single thread
  void FindLMatrix();
  void RecalculateBMatrix();
  void DivideAllRowsByFirst();
  void FindResult();
  // multi thread
  void FindLMatrixMultithread();
  void RecalculateBMatrixMultithread();
  void DivideAllRowsByFirstMultithread();
  void FindResultMultithread();
  const double EPS = 1e-7;
  int values_per_thread;
};
