#include "interfaceWinograd.h"

#include <chrono>

void InterfaceWinograd::Start() {
  std::cout << style1
            << "\n*******************************************************"
            << end << std::endl;
  std::cout << style2
            << "                Welcome to parallels                     "
            << end << std::endl;
  std::cout << style1
            << "*******************************************************" << end
            << std::endl;

  matrix_generate_or_load_ = false;
  Menu();
}

void InterfaceWinograd::Read() {
  while (true) {
    std::cout << style3
              << "Enter a filename to load the first matrix           " << end
              << std::endl;
    std::cout << std::endl;
    std::string str1;
    std::cin >> str1;
    if (std::cin.fail()) return;
    std::cout << style3
              << "Enter a filename to load the second matrix           " << end
              << std::endl;
    std::cout << std::endl;
    std::string str2;
    std::cin >> str2;
    if (std::cin.fail()) return;
    try {
      matrixA_ = loader.loadGraphFromFile("../datasets/" + str1);
      matrixB_ = loader.loadGraphFromFile("../datasets/" + str2);
      matrix_generate_or_load_ = true;
      break;
    } catch (std::exception &e) {
      std::cout << style4 << "\nError: " << e.what() << end << "\n\n";
      std::cout.flush();
      std::cout.clear();
      return;
    }
  }
  Menu();
}

void InterfaceWinograd::GenerateRandomMatrix() {
  std::cout << style3
            << "\nEnter the size of row of the Matrix A from 2 to 100 " << end1
            << std::endl;
  int sizeRowA;
  std::cin >> sizeRowA;
  std::cout << style3
            << "\nEnter the size of column of the Matrix A from 2 to 100 "
            << end1 << std::endl;
  int sizeColA;
  std::cin >> sizeColA;
  std::cout << style3
            << "\nEnter the size of row of the Matrix B from 2 to 100 " << end1
            << std::endl;
  int sizeRowB;
  std::cin >> sizeRowB;
  std::cout << style3
            << "\nEnter the size of column of the Matrix B from 2 to 100 "
            << end1 << std::endl;
  int sizeColB;
  std::cin >> sizeColB;
  if (std::cin.fail()) return;

  if (sizeRowA < 2 || sizeRowA > 100 || sizeRowB < 2 || sizeRowB > 100 ||
      sizeColA < 2 || sizeColA > 100 || sizeColB < 2 || sizeColB > 100 ||
      sizeColA != sizeRowB || std::cin.fail()) {
    std::cout << style4 << "\nError: Incorrect size matrix \n";
    std::cout.flush();
    std::cout.clear();
  } else {
    S21Matrix<double> buf1(sizeRowA, sizeColA);
    S21Matrix<double> buf2(sizeRowB, sizeColB);
    for (int i = 0; i < sizeRowA; ++i) {
      for (int j = 0; j < sizeColA; ++j) {
        if (buf1(i, j) == 0.0) {
          buf1(i, j) = winograd_.RandomVal(sizeRowA + sizeColA) + 1;
        }
      }
    }

    for (int i = 0; i < sizeRowB; ++i) {
      for (int j = 0; j < sizeColB; ++j) {
        if (buf2(i, j) == 0.0) {
          buf2(i, j) = winograd_.RandomVal(sizeRowB + sizeColB) + 1;
        }
      }
    }
    matrixA_ = std::move(buf1);
    matrixB_ = std::move(buf2);
    matrix_generate_or_load_ = true;
    std::cout << style1 << "Random " << matrixA_.get_rows() << "x"
              << matrixA_.get_cols() << " matrix A generated " << end1
              << std::endl;
    PrintMatrix(matrixA_);
    std::cout << style1 << "Random " << matrixB_.get_rows() << "x"
              << matrixB_.get_cols() << " matrix B generated " << end1
              << std::endl;
    PrintMatrix(matrixB_);
    Menu();
  }
}

void InterfaceWinograd::Menu() {
  std::cout << style1
            << "*******************************************************" << end
            << std::endl;
  std::cout << style2
            << "                           MENU                        " << end
            << std::endl;
  std::cout << style1
            << "*******************************************************" << end1
            << std::endl;
  std::cout << style3
            << " 1. ----> Generate A matrix and B matrix random                "
            << end << std::endl;
  std::cout << style3
            << " 2. ----> Load A matrix and B matrix from file                 "
            << end << std::endl;
  std::cout << style3
            << " 3. ----> Run Winograd algorithm without multithreading        "
            << end << std::endl;
  std::cout << style3
            << " 4. ----> Run Winograd algorithm with multithreading classic   "
            << end << std::endl;
  std::cout << style3
            << " 5. ----> Run Winograd algorithm with multithreading conveyor  "
            << end << std::endl;
  std::cout << style3
            << " 6. ----> Сompare Winograd algorithm with and without "
               "multithreading "
            << end << std::endl;
  std::cout << style3
            << " 0. ----> Exit from the program                              "
               "        "
            << end1 << std::endl;
  int choice;
  std::cin >> choice;
  if (choice == 1) {
    GenerateRandomMatrix();
  } else if (choice == 2) {
    Read();
  } else if (choice == 3) {
    WinogradSinglethread();
  } else if (choice == 4) {
    WinogradMultithread();
  } else if (choice == 5) {
    WinogradConv();
  } else if (choice == 6) {
    CompareWinograd();
  } else if (choice == 0) {
    //
  } else {
    std::cout << style4 << "\nError: incorrect input, try again " << end
              << std::endl;
    std::cout.flush();
    std::cout.clear();
  }
}

template <typename T>
void InterfaceWinograd::PrintMatrix(S21Matrix<T> res) {
  for (int i = 0; i < res.get_rows(); ++i) {
    for (int j = 0; j < res.get_cols(); ++j) {
      std::cout << style4 << res(i, j);
      if (res(i, j) < 10)
        std::cout << "   ";
      else if (res(i, j) >= 10 && res(i, j) < 100)
        std::cout << "  ";
      else
        std::cout << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void InterfaceWinograd::WinogradSinglethread() {
  if (matrix_generate_or_load_ == true) {
    int N;
    std::cout << style4
              << "\nenter how many times to keep track of the time (1 - 1000)"
              << end1 << std::endl;
    std::cin >> N;
    if (N > 0 && N < 1001) {
      S21Matrix<double> c2;
      std::chrono::system_clock::time_point ta1 =
          std::chrono::high_resolution_clock::now();
      for (int i = 0; i < N; ++i) {
        c2 = winograd_.WinogradSinglethread(matrixA_, matrixB_);
      }
      std::chrono::high_resolution_clock::time_point ta2 =
          std::chrono::high_resolution_clock::now();
      auto time =
          std::chrono::duration_cast<std::chrono::nanoseconds>(ta2 - ta1)
              .count();
      std::cout << style1 << "TIME: " << ((double)time / 1000000000.0) << end1
                << std::endl;
      std::cout << style1 << "MATRIX: " << std::endl;
      PrintMatrix(c2);
      Menu();
    } else {
      std::cout << style4 << "\nINCORRECT COUNT OF TIMES!" << end1 << std::endl;
      std::cout.flush();
      std::cout.clear();
      Menu();
    }

  } else {
    std::cout << style4 << "\nError: incorrect input, try again " << end
              << std::endl;
    std::cout.flush();
    std::cout.clear();
    Menu();
  }
}

void InterfaceWinograd::WinogradMultithread() {
  if (matrix_generate_or_load_ == true) {
    int N;
    std::cout << style4
              << "\nenter how many times to keep track of the time (1 - 1000)"
              << end1 << std::endl;
    std::cin >> N;
    if (N > 0 && N < 1001) {
      S21Matrix<double> c2;
      std::chrono::system_clock::time_point ta1 =
          std::chrono::high_resolution_clock::now();
      for (int i = 0; i < N; ++i) {
        c2 = winograd_.WinogradSinglethread(matrixA_, matrixB_);
      }
      std::chrono::high_resolution_clock::time_point ta2 =
          std::chrono::high_resolution_clock::now();
      auto time =
          std::chrono::duration_cast<std::chrono::nanoseconds>(ta2 - ta1)
              .count();
      std::cout << style1 << "TIME: " << ((double)time / 1000000000.0) << end1
                << std::endl;
      std::cout << style1 << "MATRIX: " << std::endl;
      PrintMatrix(c2);
      Menu();
    } else {
      std::cout << style4 << "\nINCORRECT COUNT OF TIMES!" << end1 << std::endl;
      std::cout.flush();
      std::cout.clear();
      Menu();
    }

  } else {
    std::cout << style4 << "\nError: incorrect input, try again " << end
              << std::endl;
    std::cout.flush();
    std::cout.clear();
    Menu();
  }
}

void InterfaceWinograd::WinogradConv() {
  if (matrix_generate_or_load_ == true) {
    int N;
    std::cout << style4
              << "\nenter how many times to keep track of the time (1 - 1000)"
              << end1 << std::endl;
    std::cin >> N;
    if (N > 0 && N < 1001) {
      S21Matrix<double> c2;
      std::chrono::system_clock::time_point ta1 =
          std::chrono::high_resolution_clock::now();
      for (int i = 0; i < N; ++i) {
        c2 = winograd_.WinogradConveyor(matrixA_, matrixB_);
      }
      std::chrono::high_resolution_clock::time_point ta2 =
          std::chrono::high_resolution_clock::now();
      auto time =
          std::chrono::duration_cast<std::chrono::nanoseconds>(ta2 - ta1)
              .count();
      std::cout << style1 << "TIME: " << ((double)time / 1000000000.0) << end1
                << std::endl;
      std::cout << style1 << "MATRIX: " << std::endl;
      PrintMatrix(c2);
      Menu();
    } else {
      std::cout << style4 << "\nINCORRECT COUNT OF TIMES!" << end1 << std::endl;
      std::cout.flush();
      std::cout.clear();
      Menu();
    }

  } else {
    std::cout << style4 << "\nError: incorrect input, try again " << end
              << std::endl;
    std::cout.flush();
    std::cout.clear();
    Menu();
  }
}

void InterfaceWinograd::CompareWinograd() {
  if (matrix_generate_or_load_ == true) {
    int N;
    std::cout << style4
              << "\nenter how many times to keep track of the time (1 - 1000)"
              << end1 << std::endl;
    std::cin >> N;
    if (N > 0 && N < 1001) {
      S21Matrix<double> c2;
      std::chrono::system_clock::time_point ta1 =
          std::chrono::high_resolution_clock::now();

      for (int i = 0; i < N; ++i) {
        c2 = winograd_.WinogradSinglethread(matrixA_, matrixB_);
      }

      std::chrono::high_resolution_clock::time_point ta2 =
          std::chrono::high_resolution_clock::now();

      auto timeA =
          std::chrono::duration_cast<std::chrono::nanoseconds>(ta2 - ta1)
              .count();
      S21Matrix<double> c3;
      ta1 = std::chrono::high_resolution_clock::now();

      for (int i = 0; i < N; ++i) {
        c3 = winograd_.WinogradMultithread(matrixA_, matrixB_);
      }

      ta2 = std::chrono::high_resolution_clock::now();

      auto timeB =
          std::chrono::duration_cast<std::chrono::nanoseconds>(ta2 - ta1)
              .count();
      S21Matrix<double> c4;
      ta1 = std::chrono::high_resolution_clock::now();

      for (int i = 0; i < N; ++i) {
        c4 = winograd_.WinogradConveyor(matrixA_, matrixB_);
      }

      ta2 = std::chrono::high_resolution_clock::now();

      auto timeC =
          std::chrono::duration_cast<std::chrono::nanoseconds>(ta2 - ta1)
              .count();
      std::cout << style1 << "MATRIX single: " << std::endl;
      PrintMatrix(c2);
      std::cout << style1 << "MATRIX with multithreading: " << std::endl;
      PrintMatrix(c3);
      std::cout << style1 << "MATRIX conveyor: " << std::endl;
      PrintMatrix(c4);
      std::cout << style1 << "Time single : " << ((double)timeA / 1000000000.0)
                << end1 << std::endl;
      std::cout << style1 << "Time with multithreading : "
                << ((double)timeB / 1000000000.0) << end1 << std::endl;
      std::cout << style1
                << "Time conveyor : " << ((double)timeC / 1000000000.0) << end1
                << std::endl;
      Menu();
    } else {
      std::cout << style4 << "\nINCORRECT COUNT OF TIMES!" << end1 << std::endl;
      std::cout.flush();
      std::cout.clear();
      Menu();
    }

  } else {
    std::cout << style4 << "\nError: incorrect input, try again " << end
              << std::endl;
    std::cout.flush();
    std::cout.clear();
    Menu();
  }
}
