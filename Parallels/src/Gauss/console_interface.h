#pragma once

#include <iostream>

#include "../s21_matrix_oop.h"
#include "gauss.h"

class ConsoleInterfaceGauss {
 public:
  void Start();

 private:
  void Read();
  void Menu();

  void PrintRoots(S21Matrix<double> roots);
  void GenerateRandomMatrix();
  void SolveGauss(bool multithreading, int N);
  int AskN();

  S21Matrix<double> matrix_;
  Gauss gauss;
  int size_;
  bool matrixGenerateOrLOad_;

  std::string end = "\u001b[0m";
  std::string end1 = "\u001b[0m\n";
  std::string style1 = "\u001b[1;35;5;117m";
  std::string style2 = "\u001b[1;33;5;117m";
  std::string style3 = "\u001b[1;32;5;117m";
  std::string style4 = "\u001b[1;31;5;117m";
};
