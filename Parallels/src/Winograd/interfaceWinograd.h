#pragma once

#include <iostream>

#include "../loader.h"
#include "../s21_matrix_oop.h"
#include "winograd.h"

class InterfaceWinograd {
 public:
  void Start();

 private:
  void Read();
  void Menu();

  template <typename T>
  void PrintMatrix(S21Matrix<T> res);
  void GenerateRandomMatrix();
  void WinogradSinglethread();
  void WinogradMultithread();
  void WinogradConv();
  void CompareWinograd();

  S21Matrix<double> matrixA_;
  S21Matrix<double> matrixB_;
  Winograd winograd_;
  Loader loader;
  bool matrix_generate_or_load_;

  std::string end = "\u001b[0m";
  std::string end1 = "\u001b[0m\n";
  std::string style1 = "\u001b[1;35;5;117m";
  std::string style2 = "\u001b[1;33;5;117m";
  std::string style3 = "\u001b[1;32;5;117m";
  std::string style4 = "\u001b[1;31;5;117m";
};
