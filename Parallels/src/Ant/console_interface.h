#pragma once

#include <iostream>

#include "../loader.h"
#include "s21_graph_algorithms.h"
#include "../s21_matrix_oop.h"

class ConsoleInterface {
 public:
  void Start();

 private:
  void Read();
  void Menu();
  template <typename T>
  void PrintVector(std::vector<T> path);
  template <typename T>
  void PrintMatrix(S21Matrix<T> res);
  void GenerateRandomMatrix();
  void Ant(bool multithreading, int count, int *n = 0);

  S21Matrix<double> matrix_;
  GraphAlgorithms algo_;
  Loader loader;
  bool matrix_generate_or_load_;

  std::string end = "\u001b[0m";
  std::string end1 = "\u001b[0m\n";
  std::string style1 = "\u001b[1;35;5;117m";
  std::string style2 = "\u001b[1;33;5;117m";
  std::string style3 = "\u001b[1;32;5;117m";
  std::string style4 = "\u001b[1;31;5;117m";
};
