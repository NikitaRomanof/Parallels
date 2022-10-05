#include "winograd.h"

#include <thread>

S21Matrix<double> Winograd::WinogradSinglethread(S21Matrix<double> a,
                                                         S21Matrix<double> b) {
  if (a.get_cols() != b.get_rows()) {
    throw std::invalid_argument(
        "Error on winogradWithoutMultithreadin. Columns are not equal to "
        "rows.");
  }
  if (a.get_cols() == 0 || a.get_rows() == 0 || b.get_cols() == 0 ||
      b.get_rows() == 0) {
    throw std::invalid_argument(
        "Error on winogradWithoutMultithreadin. Invalid matrix)");
  }

  std::vector<double> row_factor(a.get_rows());
  std::vector<double> column_factor(b.get_cols());
  int d = a.get_cols() / 2;

  CalculateRowsFactor(&row_factor, a, d);
  CalculateColsFactor(&column_factor, b, d);
  return CalculateResult(row_factor, column_factor, a, b, d);
}

void Winograd::CalculateRowsFactor(std::vector<double> *row_factor,
                                   const S21Matrix<double> &matrixA,
                                   const int &d) {
  for (int i = 0; i < matrixA.get_rows(); ++i) {
    row_factor->operator[](i) = 0.0;
    for (int j = 0; j < d; ++j) {
      row_factor->operator[](i) += matrixA(i, 2 * j) * matrixA(i, 2 * j + 1);
    }
  }
}

void Winograd::CalculateColsFactor(std::vector<double> *col_factor,
                                   const S21Matrix<double> &matrixB,
                                   const int &d) {
  for (int i = 0; i < matrixB.get_cols(); ++i) {
    col_factor->operator[](i) = 0.0;
    for (int j = 0; j < d; ++j) {
      col_factor->operator[](i) += matrixB(2 * j, i) * matrixB(2 * j + 1, i);
    }
  }
}

S21Matrix<double> Winograd::CalculateResult(
    const std::vector<double> &row_factor, const std::vector<double> &col_factor,
    const S21Matrix<double> &matrixA, const S21Matrix<double> &matrixB,
    const int &d) {
  S21Matrix<double> rez(matrixA.get_rows(), matrixB.get_cols());
  for (int i = 0; i < matrixA.get_rows(); ++i) {
    for (int j = 0; j < matrixB.get_cols(); ++j) {
      rez(i, j) = -row_factor[i] - col_factor[j];
      for (int k = 0; k < d; k++) {
        rez(i, j) += (matrixA(i, 2 * k) + matrixB(2 * k + 1, j)) *
                     (matrixA(i, 2 * k + 1) + matrixB(2 * k, j));
      }
      if (matrixA.get_cols() % 2 != 0) {
        rez(i, j) += matrixA(i, matrixA.get_cols() - 1) *
                     matrixB(matrixA.get_cols() - 1, j);
      }
    }
  }
  return rez;
}

S21Matrix<double> Winograd::WinogradMultithread(S21Matrix<double> a,
                                                      S21Matrix<double> b) {
  if (a.get_cols() != b.get_rows()) {
    throw std::invalid_argument(
        "Error on winogradWithoutMultithreadin. Columns are not equal to "
        "rows.");
  }
  if (a.get_cols() == 0 || a.get_rows() == 0 || b.get_cols() == 0 ||
      b.get_rows() == 0) {
    throw std::invalid_argument(
        "Error on winogradWithoutMultithreadin. Invalid matrix)");
  }
  S21Matrix<double> rez(a.get_rows(), b.get_cols());
  rez.fillMatrix(0.0);
  std::vector<double> row_factor(a.get_cols());
  std::vector<double> column_factor(b.get_rows());
  int d = a.get_cols() / 2;
  std::thread t1([&] { CalculateRowsFactor(&row_factor, a, d); });
  std::thread t2([&] { CalculateColsFactor(&column_factor, b, d); });
  t1.join();
  t2.join();
  return CalculateResult(row_factor, column_factor, a, b, d);
}

S21Matrix<double> Winograd::WinogradConveyor(S21Matrix<double> a,
                                             S21Matrix<double> b) {
  if (a.get_cols() != b.get_rows()) {
    throw std::invalid_argument(
        "Error on winogradWithoutMultithreadin. Columns are not equal to "
        "rows.");
  }
  if (a.get_cols() == 0 || a.get_rows() == 0 || b.get_cols() == 0 ||
      b.get_rows() == 0) {
    throw std::invalid_argument(
        "Error on winogradWithoutMultithreadin. Invalid matrix)");
  }

  std::map<int, double> rowF;
  std::vector<double> row_factor(a.get_cols(), INFINITY);
  std::vector<double> column_factor(b.get_rows(), INFINITY);
  S21Matrix<double> first(a.get_rows(), b.get_cols());
  first.fillMatrix(INFINITY);
  S21Matrix<double> rez(a.get_rows(), b.get_cols());

  std::thread ta([&] { CalculateRowsFactor(&row_factor, a, a.get_cols() / 2); });
  std::thread tb(
      [&] { CalculateColsFactor(&column_factor, b, a.get_cols() / 2); });
  std::thread tf([&] { Calcfirst_matrix(a, b, &first); });
  std::thread tr(
      [&] { CalcResultMatrix(row_factor, column_factor, first, &rez); });
  ta.join();
  tb.join();
  tf.join();
  tr.join();
  return rez;
}

void Winograd::Calcfirst_matrix(const S21Matrix<double> &matrixA,
                               const S21Matrix<double> &matrixB,
                               S21Matrix<double> *first_matrix) {
  for (int i = 0; i < matrixA.get_rows(); ++i) {
    for (int j = 0; j < matrixB.get_cols(); ++j) {
      first_matrix->operator()(i, j) = 0.0;
      for (int k = 0; k < matrixA.get_cols() / 2; k++) {
        first_matrix->operator()(i, j) +=
            (matrixA(i, 2 * k) + matrixB(2 * k + 1, j)) *
            (matrixA(i, 2 * k + 1) + matrixB(2 * k, j));
      }
      if (matrixA.get_cols() % 2 != 0) {
        first_matrix->operator()(i, j) += matrixA(i, matrixA.get_cols() - 1) *
                                         matrixB(matrixA.get_cols() - 1, j);
      }
    }
  }
}

void Winograd::CalcResultMatrix(const std::vector<double> &row_factor,
                                const std::vector<double> &col_factor,
                                const S21Matrix<double> &first_matrix,
                                S21Matrix<double> *rez_matrix) {
  for (int i = 0; i < first_matrix.get_rows(); ++i) {
    for (int j = 0; j < first_matrix.get_cols();) {
      if (std::isinf(first_matrix(i, j)) == false &&
          std::isinf(row_factor[i]) == false &&
          std::isinf(col_factor[j]) == false) {
        rez_matrix->operator()(i, j) =
            (-row_factor[i] - col_factor[j]) + first_matrix(i, j);
        ++j;
      }
    }
  }
}

int Winograd::RandomVal(int val) {
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_int_distribution<int> distr(0, val);
  return distr(eng);
}
