#pragma once

#include <fstream>
#include <iostream>

#include "s21_matrix_oop.h"

class Loader {
 public:
  S21Matrix<double> loadGraphFromFile(const std::string& filename);
};
