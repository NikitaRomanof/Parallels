
// Constructor
template <typename T>
S21Matrix<T>::S21Matrix() {
  _rows = 2;
  _cols = 2;
  _matrix = nullptr;
  memoryAllocation();
}

template <typename T>
S21Matrix<T>::S21Matrix(int rows, int cols) {
  if (rows <= 0 || cols <= 0) {
    MatrixException nonExist(
        "The number of rows/columns of the matrix is less than 1");
    throw nonExist;
  }
  _rows = rows;
  _cols = cols;
  _matrix = nullptr;
  memoryAllocation();
}

template <typename T>
S21Matrix<T>::S21Matrix(const S21Matrix<T>& other) {
  _matrix = nullptr;
  _rows = 0;
  _cols = 0;
  *this = other;
}

template <typename T>
S21Matrix<T>::S21Matrix(S21Matrix<T>&& other) {
  _matrix = nullptr;
  _rows = 0;
  _cols = 0;
  *this = std::move(other);
}

//  Destructor
template <typename T>
S21Matrix<T>::~S21Matrix() {
  memoryClean();
}

// Function
template <typename T>
void S21Matrix<T>::memoryAllocation() {
  _matrix = new T*[_rows];
  for (int i = 0; i < _rows; i++) {
    _matrix[i] = new T[_cols]{};
  }
}

template <typename T>
void S21Matrix<T>::memoryClean() {
  if (_matrix != nullptr) {
    for (int i = 0; i < _rows; i++) {
      delete[] _matrix[i];
    }
    delete[] _matrix;
    _matrix = nullptr;
    _rows = 0;
    _cols = 0;
  }
}

template <typename T>
void S21Matrix<T>::sumOrSub(const S21Matrix<T>& other, int sign) {
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = _matrix[i][j] + (sign * other._matrix[i][j]);
    }
  }
}

template <typename T>
bool S21Matrix<T>::eq_matrix(const S21Matrix<T>& other) {
  bool rez = false;
  if (_rows == other._rows && _cols == other._cols) {
    rez = true;
    for (int i = 0; i < _rows && rez == true; i++) {
      for (int j = 0; j < _cols && rez == true; j++) {
        if (fabs(_matrix[i][j] - other._matrix[i][j]) >= 1e-7) {
          rez = false;
        }
      }
    }
  }
  return rez;
}

template <typename T>
void S21Matrix<T>::sum_matrix(const S21Matrix<T>& other) {
  if (_rows != other._rows || _cols != other._cols) {
    MatrixException inconsistency("Rows/columns of matrices are not equal");
    throw inconsistency;
  }
  sumOrSub(other, 1);
}

template <typename T>
void S21Matrix<T>::sub_matrix(const S21Matrix<T>& other) {
  if (_rows != other._rows || _cols != other._cols) {
    MatrixException inconsistency("Rows/columns of matrices are not equal");
    throw inconsistency;
  }
  sumOrSub(other, -1);
}

template <typename T>
void S21Matrix<T>::mul_number(const T& val) {
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = _matrix[i][j] * val;
    }
  }
}

template <typename T>
void S21Matrix<T>::mul_matrix(const S21Matrix<T>& other) {
  if (_cols != other._rows) {
    MatrixException nonExist("Columns are not equal to rows");
    throw nonExist;
  }
  S21Matrix<T> buf(_rows, other._cols);
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < other._cols; j++) {
      for (int n = 0; n < other._rows; n++) {
        buf._matrix[i][j] += _matrix[i][n] * other._matrix[n][j];
      }
    }
  }
  *this = buf;
}

template <typename T>
S21Matrix<T> S21Matrix<T>::transpose() {
  S21Matrix<T> buf(_cols, _rows);
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      buf._matrix[j][i] = _matrix[i][j];
    }
  }
  return buf;
}

// Operator
template <typename T>
S21Matrix<T>& S21Matrix<T>::operator=(const S21Matrix<T>& other) {
  if (this == &other) {
    MatrixException nonExist("Cannot assign an object's value to itself");
    throw nonExist;
  }
  memoryClean();
  _rows = other._rows;
  _cols = other._cols;
  memoryAllocation();

  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = other._matrix[i][j];
    }
  }
  return *this;
}

template <typename T>
S21Matrix<T>& S21Matrix<T>::operator=(S21Matrix<T>&& other) {
  if (this == &other) {
    MatrixException nonExist("Cannot assign an object's value to itself");
    throw nonExist;
  }
  memoryClean();
  std::swap(_rows, other._rows);
  std::swap(_cols, other._cols);
  std::swap(_matrix, other._matrix);
  return *this;
}

template <typename T>
S21Matrix<T>& S21Matrix<T>::operator+=(const S21Matrix<T>& other) {
  sum_matrix(other);
  return *this;
}

template <typename T>
S21Matrix<T>& S21Matrix<T>::operator-=(const S21Matrix<T>& other) {
  sub_matrix(other);
  return *this;
}

template <typename T>
S21Matrix<T>& S21Matrix<T>::operator*=(const T& val) {
  mul_number(val);
  return *this;
}

template <typename T>
S21Matrix<T>& S21Matrix<T>::operator*=(const S21Matrix<T>& other) {
  mul_matrix(other);
  return *this;
}

template <typename T>
bool S21Matrix<T>::operator==(const S21Matrix<T>& other) {
  return eq_matrix(other);
}

template <typename T>
S21Matrix<T> S21Matrix<T>::operator+(const S21Matrix<T>& other) {
  S21Matrix<T> buf(*this);
  buf.sum_matrix(other);
  return buf;
}

template <typename T>
S21Matrix<T> S21Matrix<T>::operator-(const S21Matrix<T>& other) {
  S21Matrix<T> buf(*this);
  buf.sub_matrix(other);
  return buf;
}

template <typename T>
S21Matrix<T> S21Matrix<T>::operator*(const T& val) {
  S21Matrix<T> buf(*this);
  buf.mul_number(val);
  return buf;
}

template <typename T>
S21Matrix<T> S21Matrix<T>::operator*(const S21Matrix<T>& other) {
  S21Matrix<T> buf(*this);
  buf.mul_matrix(other);
  return buf;
}

template <typename T>
T& S21Matrix<T>::operator()(int row, int col) {
  if (row < 0 || col < 0 || row >= _rows || col >= _cols) {
    MatrixException non("Incorrect row/col");
    throw non;
  }
  return _matrix[row][col];
}

template <typename T>
T S21Matrix<T>::operator()(int row, int col) const {
  if (row < 0 || col < 0 || row >= _rows || col >= _cols) {
    MatrixException non("Incorrect row/col");
    throw non;
  }
  return _matrix[row][col];
}

// GET SET
template <typename T>
int S21Matrix<T>::get_rows() {
  return _rows;
}

template <typename T>
int S21Matrix<T>::get_rows() const {
  return _rows;
}

template <typename T>
int S21Matrix<T>::get_cols() {
  return _cols;
}

template <typename T>
int S21Matrix<T>::get_cols() const {
  return _cols;
}

template <typename T>
void S21Matrix<T>::set_rows(int rows) {
  S21Matrix<T> buf(rows, _cols);
  int minRow = _rows > rows ? rows : _rows;
  for (int i = 0; i < minRow; i++) {
    for (int j = 0; j < _cols; j++) {
      buf._matrix[i][j] = _matrix[i][j];
    }
  }
  *this = buf;
}

template <typename T>
void S21Matrix<T>::set_cols(int cols) {
  S21Matrix<T> buf(_rows, cols);
  int minCol = _cols > cols ? cols : _cols;
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < minCol; j++) {
      buf._matrix[i][j] = _matrix[i][j];
    }
  }
  *this = buf;
}

template <typename T>
void S21Matrix<T>::fillMatrix(const T& val) {
  for (int i = 0; i < _rows; i++) {
    for (int j = 0; j < _cols; j++) {
      _matrix[i][j] = val;
    }
  }
}