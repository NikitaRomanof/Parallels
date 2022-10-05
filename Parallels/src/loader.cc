#include "loader.h"

#include <random>

S21Matrix<double> Loader::loadGraphFromFile(const std::string& filename) {
  std::ifstream fs(filename);
  if (!fs.is_open()) {
    throw std::runtime_error("Can`t open file");
  }

  if (fs.eof() == true) {
    throw std::runtime_error("Empty file");
  }

  int row = 0;
  fs >> row;
  int col = 0;
  fs >> col;

  if (fs.fail() == true) {
    throw std::runtime_error("Bad data on file");
  }

  S21Matrix<double> rez(row, col);

  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      fs >> rez(i, j);
      if (fs.fail() == true) {
        throw std::runtime_error("Bad data on file");
      }
    }
  }

  fs.close();
  return rez;
}
