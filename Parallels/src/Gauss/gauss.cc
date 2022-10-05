#include "gauss.h"

Gauss::Gauss()
    : size_(0),
      multi_thread_(false),
      matrix_(nullptr),
      b_matrix_(nullptr),
      l_matrix_(nullptr),
      result_(nullptr),
      values_per_thread(1) {}

Gauss::Gauss(S21Matrix<double> mat)
    : size_(mat.get_rows()),
      multi_thread_(false),
      matrix_(nullptr),
      b_matrix_(nullptr),
      l_matrix_(nullptr),
      result_(nullptr),
      values_per_thread(1) {
  matrix_ = new S21Matrix<double>(mat);
}

Gauss::Gauss(S21Matrix<double> mat, S21Matrix<double> b)
    : size_(mat.get_rows()),
      multi_thread_(false),
      matrix_(nullptr),
      b_matrix_(nullptr),
      l_matrix_(nullptr),
      result_(nullptr),
      values_per_thread(1) {
  matrix_ = new S21Matrix<double>(mat);
  b_matrix_ = new S21Matrix<double>(b);
}

Gauss::~Gauss() { Clear(); }

void Gauss::GetLinearSystemFromFile(const std::string &filename) {
  std::ifstream fs(filename);
  if (!fs.is_open()) {
    throw std::runtime_error("Can`t open file");
  }

  if (fs.eof() == true) {
    throw std::runtime_error("Empty file");
  }

  int fileSize = 0;
  fs >> fileSize;

  if (fs.fail() == true) {
    throw std::runtime_error("Bad data on file");
  }

  size_ = fileSize;
  matrix_ = new S21Matrix<double>(fileSize, fileSize);
  b_matrix_ = new S21Matrix<double>(1, fileSize);

  for (int i = 0; i < size_; i++) {
    for (int j = 0; j < size_ + 1; j++) {
      if (j < size_)
        fs >> matrix_->operator()(i, j);
      else
        fs >> b_matrix_->operator()(0, j - size_ + i);
      if (fs.fail() == true) {
        throw std::runtime_error("Bad data on file");
      }
    }
  }
  fs.close();
}
void Gauss::GenerateLinearSystem(int size) {
  // cгенерируем матрицу
  size_ = size;
  std::random_device dev;
  std::mt19937 gen(dev());
  std::uniform_real_distribution<double> distr(-10, 10);
  if (matrix_ == nullptr) matrix_ = new S21Matrix<double>(size, size);
  if (b_matrix_ == nullptr) b_matrix_ = new S21Matrix<double>(1, size);
  for (int i = 0; i < size_; ++i) {
    for (int j = 0; j < size_; ++j) {
      while (matrix_->operator()(i, j) == 0) {
        matrix_->operator()(i, j) = (double)(int)distr(gen);
      }
    }
    b_matrix_->operator()(0, i) = (double)(int)distr(gen);
  }
}

S21Matrix<double> *Gauss::SolveLinearSystemSingleThread() {
  values_per_thread = (int)std::ceil(double(size_) / 4.0);
  if (0 == values_per_thread) values_per_thread = 1;
  if (IsEmpty() == false) {
    FindLMatrix();
    RecalculateBMatrix();
    DivideAllRowsByFirst();
    FindResult();
  } else {
    if (result_ != nullptr) {
      delete result_;
      result_ = nullptr;
    }
    result_ = new S21Matrix<double>(size_, size_);
  }
  return result_;
}

S21Matrix<double> *Gauss::SolveLinearSystemMultiThread() {
  values_per_thread = (int)std::ceil(double(size_) / 4.0);
  if (0 == values_per_thread) values_per_thread = 1;
  if (IsEmpty() == false) {
    FindLMatrixMultithread();
    RecalculateBMatrixMultithread();
    DivideAllRowsByFirstMultithread();
    FindResultMultithread();
  } else {
    if (result_ != nullptr) {
      delete result_;
      result_ = nullptr;
    }
    result_ = new S21Matrix<double>(size_, size_);
  }
  return result_;
}

void Gauss::Clear() {
  if (matrix_ != nullptr) {
    delete matrix_;
    matrix_ = nullptr;
  }
  if (b_matrix_ != nullptr) {
    delete b_matrix_;
    b_matrix_ = nullptr;
  }
  if (l_matrix_ != nullptr) {
    delete l_matrix_;
    l_matrix_ = nullptr;
  }
  if (result_ != nullptr) {
    delete result_;
    result_ = nullptr;
  }
}

bool Gauss::IsEmpty() {
  if (matrix_ == nullptr) return true;
  return false;
}
void Gauss::FindLMatrix() {
  if (l_matrix_ == nullptr) l_matrix_ = new S21Matrix<double>(size_, size_);
  for (int k = 0; k < size_; ++k) {
    for (int i = k + 1; i < size_; ++i) {
      l_matrix_->operator()(i, k) =
          matrix_->operator()(i, k) / matrix_->operator()(k, k);
    }
    for (int j = 0; j < size_; ++j) {
      for (int i = k + 1; i < size_; ++i) {
        matrix_->operator()(i, j) =
            matrix_->operator()(i, j) -
            l_matrix_->operator()(i, k) * matrix_->operator()(k, j);
      }
    }
  }
  for (int i = 0; i < size_; ++i) {
    for (int j = 0; j < size_; ++j) {
      if (std::fabs(matrix_->operator()(i, j)) < EPS)
        matrix_->operator()(i, j) = 0;  // можно заменить 1e-7 на константу
    }
  }
}

void Gauss::RecalculateBMatrix() {
  for (int k = 0; k < size_; ++k) {
    for (int i = k + 1; i < size_; ++i) {
      if (l_matrix_->operator()(i, k) != 0)
        b_matrix_->operator()(0, i) =
            b_matrix_->operator()(0, i) -
            b_matrix_->operator()(0, k) * l_matrix_->operator()(i, k);
    }
  }
}

void Gauss::DivideAllRowsByFirst() {
  for (int i = 0; i < size_; ++i) {
    double divider = 1;
    for (int j = 0; j < size_; ++j) {
      if (i == j) divider = matrix_->operator()(i, j);
      matrix_->operator()(i, j) = matrix_->operator()(i, j) / divider;
    }
    b_matrix_->operator()(0, i) = b_matrix_->operator()(0, i) / divider;
  }
}
void Gauss::FindResult() {
  // подстановка значений в уравнение, начиная с последнего
  if (result_ == nullptr) result_ = new S21Matrix<double>(1, size_);
  for (int i = size_ - 1; i >= 0; --i) {
    double subtractor = 0;
    for (int j = size_ - 1; j >= 0; --j) {
      if (i != j)
        subtractor += matrix_->operator()(i, j) * result_->operator()(0, j);
    }
    result_->operator()(0, i) = b_matrix_->operator()(0, i) - subtractor;
  }
}

void Gauss::FindLMatrixMultithread() {
  if (l_matrix_ == nullptr) l_matrix_ = new S21Matrix<double>(size_, size_);
  for (int k = 0; k < size_; ++k) {
    int start1 = k + 1;
    int start2 = start1 + values_per_thread;
    int start3 = start2 + values_per_thread;
    int start4 = start3 + values_per_thread;
    int end = start4 + values_per_thread;
    auto f = [&](int start, int end) {
      if (end > size_) end = size_;
      if (start < size_) {
        for (int i = start; i < end; ++i) {
          l_matrix_->operator()(i, k) =
              matrix_->operator()(i, k) / matrix_->operator()(k, k);
        }
      }
    };
    std::thread ta(f, start1, start2);
    std::thread tb(f, start2, start3);
    std::thread tc(f, start3, start4);
    std::thread td(f, start4, end);
    ta.join();
    tb.join();
    tc.join();
    td.join();
    for (int j = 0; j < size_; ++j) {
      for (int i = k + 1; i < size_; ++i) {
        matrix_->operator()(i, j) =
            matrix_->operator()(i, j) -
            l_matrix_->operator()(i, k) * matrix_->operator()(k, j);
      }
    }
  }
  for (int i = 0; i < size_; ++i) {
    for (int j = 0; j < size_; ++j) {
      if (std::fabs(matrix_->operator()(i, j)) < 1e-7)
        matrix_->operator()(i, j) = 0;  // можно заменить 1e-7 на константу
    }
  }
}

void Gauss::RecalculateBMatrixMultithread() {
  int start1 = 0;
  int start2 = start1 + values_per_thread;
  int start3 = start2 + values_per_thread;
  int start4 = start3 + values_per_thread;
  int end1 = start1 + values_per_thread;
  int end2 = start2 + values_per_thread;
  int end3 = start3 + values_per_thread;
  int end4 = start4 + values_per_thread;
  auto f = [&](int start, int end) {
    if (end > size_) end = size_;
    if (start < size_) {
      for (int k = start; k < end; ++k) {
        for (int i = k + 1; i < size_; ++i) {
          if (l_matrix_->operator()(i, k) != 0)
            b_matrix_->operator()(0, i) =
                b_matrix_->operator()(0, i) -
                b_matrix_->operator()(0, k) * l_matrix_->operator()(i, k);
        }
      }
    }
  };
  std::thread ta(f, start1, end1);
  std::thread tb(f, start2, end2);
  std::thread tc(f, start3, end3);
  std::thread td(f, start4, end4);
  ta.join();
  tb.join();
  tc.join();
  td.join();
}

void Gauss::DivideAllRowsByFirstMultithread() {
  int start1 = 0;
  int start2 = start1 + values_per_thread;
  int start3 = start2 + values_per_thread;
  int start4 = start3 + values_per_thread;
  int end = start4 + values_per_thread;
  auto f = [&](int start, int end) {
    if (end > size_) end = size_;
    if (start < size_) {
      for (int i = start; i < end; ++i) {
        double divider = 1;
        for (int j = 0; j < size_; ++j) {
          if (i == j) divider = matrix_->operator()(i, j);
          matrix_->operator()(i, j) = matrix_->operator()(i, j) / divider;
        }
        b_matrix_->operator()(0, i) = b_matrix_->operator()(0, i) / divider;
      }
    }
  };
  std::thread t1(f, start1, start2);
  std::thread t2(f, start2, start3);
  std::thread t3(f, start3, start4);
  std::thread t4(f, start4, end);
  t1.join();
  t2.join();
  t3.join();
  t4.join();
}

void Gauss::FindResultMultithread() {
  if (result_ == nullptr) result_ = new S21Matrix<double>(1, size_);
  int start1 = size_ - 1;
  int start2 = start1 - values_per_thread;
  int start3 = start2 - values_per_thread;
  int start4 = start3 - values_per_thread;
  int end1 = start1 - values_per_thread;
  int end2 = start2 - values_per_thread;
  int end3 = start3 - values_per_thread;
  int end4 = 0;
  auto f = [&](int start, int end) {
    if (start < 0) start = 0;
    if (end >= 0) {
      for (int i = start; i >= end; --i) {
        double subtractor = 0;
        for (int j = size_ - 1; j >= 0; --j) {
          if (i != j)
            subtractor += matrix_->operator()(i, j) * result_->operator()(0, j);
        }
        result_->operator()(0, i) = b_matrix_->operator()(0, i) - subtractor;
      }
    }
  };
  std::thread ta(f, start1, end1);
  std::thread tb(f, start2, end2);
  std::thread tc(f, start3, end3);
  std::thread td(f, start4, end4);
  ta.join();
  tb.join();
  tc.join();
  td.join();
}

void Gauss::PrintMatrix(const std::string &label) {
  std::cout << "\n" << label << "\n";
  for (int m = 0; m < matrix_->get_cols(); ++m) {
    for (int n = 0; n < matrix_->get_rows(); ++n) {
      std::cout.precision(3);
      std::cout.width(6);
      std::cout << matrix_->operator()(m, n) << " ";
    }
    std::cout << "= " << b_matrix_->operator()(0, m) << "\n";
  }
}

int Gauss::GetSize() { return size_; }
