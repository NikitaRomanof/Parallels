#pragma once

#include <map>
#include <random>
#include <vector>

#include "../s21_matrix_oop.h"

class Winograd {
 public:
  S21Matrix<double> WinogradSinglethread(S21Matrix<double> a,
                                                 S21Matrix<double> b);
  S21Matrix<double> WinogradMultithread(S21Matrix<double> a,
                                              S21Matrix<double> b);
  S21Matrix<double> WinogradConveyor(S21Matrix<double> a, S21Matrix<double> b);
  int RandomVal(int val);

 private:
  void CalculateRowsFactor(std::vector<double> *row_factor,
                           const S21Matrix<double> &matrixA, const int &d);
  void CalculateColsFactor(std::vector<double> *col_factor,
                           const S21Matrix<double> &matrixB, const int &d);
  S21Matrix<double> CalculateResult(const std::vector<double> &row_factor,
                                    const std::vector<double> &col_factor,
                                    const S21Matrix<double> &matrixA,
                                    const S21Matrix<double> &matrixB,
                                    const int &d);
  void CalculateRowsFactorConveyor(const S21Matrix<double> &matrixA,
                                   const S21Matrix<double> &matrixB,
                                   std::vector<double> **row_factor);
  void CalculateColsFactorConveyor(const S21Matrix<double> &matrixA,
                                   const S21Matrix<double> &matrixB,
                                   std::vector<double> **col_factor);
  void Calcfirst_matrix(const S21Matrix<double> &matrixA,
                       const S21Matrix<double> &matrixB,
                       S21Matrix<double> *first_matrix);
  void CalcResultMatrix(const std::vector<double> &row_factor,
                        const std::vector<double> &col_factor,
                        const S21Matrix<double> &first_matrix,
                        S21Matrix<double> *rez_matrix);
};
