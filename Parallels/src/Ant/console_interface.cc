#include "console_interface.h"

#include <chrono>

void ConsoleInterface::Start() {
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

void ConsoleInterface::Read() {
  while (true) {
    std::cout << style3 << " Enter file name           " << end << std::endl;
    std::cout << std::endl;
    std::string str;
    std::cin >> str;
    if (std::cin.fail() == false) {
      try {
        matrix_ = loader.loadGraphFromFile("../datasets/" + str);
        matrix_generate_or_load_ = true;
        break;
      } catch (std::exception &e) {
        std::cout << style4 << "\nError: " << e.what() << end << "\n\n";
        std::cout.flush();
        std::cout.clear();
        break;
      }
    } else {
      break;
    }
  }
  Menu();
}

void ConsoleInterface::GenerateRandomMatrix() {
  std::cout << style3 << "\nEnter the size of the square matrix from 2 to 100 "
            << end1 << std::endl;
  int size;
  std::cin >> size;

  if (size < 2 || size > 100 || std::cin.fail()) {
    std::cout << style4 << "\nError: Incorrect size matrix \n";
    std::cout.flush();
    std::cout.clear();
  } else {
    S21Matrix<double> buf(size, size);
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) {
        if (buf(i, j) == 0.0) {
          buf(i, j) = algo_.RandomVal(size) + 1;
          buf(j, i) = buf(i, j);
        }
      }
    }
    matrix_ = std::move(buf);
    matrix_generate_or_load_ = true;
    std::cout << style1 << "Random " << size << "x" << size
              << " matrix generated " << end1 << std::endl;
    PrintMatrix(matrix_);
    Menu();
  }
}

void ConsoleInterface::Menu() {
  std::cout << style1
            << "*******************************************************" << end
            << std::endl;
  std::cout << style2
            << "                           MENU                        " << end
            << std::endl;
  std::cout << style1
            << "*******************************************************" << end1
            << std::endl;
  std::cout
      << style3
      << " 1. ----> Generate random matrix                                 "
      << end << std::endl;
  std::cout
      << style3
      << " 2. ----> Load matrix from file                                  "
      << end << std::endl;
  std::cout
      << style3
      << " 3. ----> Run Ant algorithm without multithreading               "
      << end << std::endl;
  std::cout
      << style3
      << " 4. ----> Run Ant algorithm with multithreading                  "
      << end << std::endl;
  std::cout
      << style3
      << " 5. ----> Сompare ant algorithm with and without multithreading  "
      << end << std::endl;
  std::cout
      << style3
      << " 0. ----> Exit from the program                                  "
      << end1 << std::endl;
  int choice;
  std::cin >> choice;
  if (choice == 1) {
    GenerateRandomMatrix();
  } else if (choice == 2) {
    Read();
  } else if (choice == 3) {
    Ant(false, 3);
  } else if (choice == 4) {
    Ant(true, 4);
  } else if (choice == 5) {
    int n = 0;
    Ant(false, 5, &n);
    Ant(true, 6, &n);
  } else if (choice == 0) {
    //
  } else {
    std::cout << style4 << "\nError: incorrect input, try again " << end
              << std::endl;
    std::cout.flush();
    std::cout.clear();
  }
}

void ConsoleInterface::Ant(bool multithreading, int count, int *n) {
  if (matrix_generate_or_load_ == true) {
    try {
      int N;
      if (count != 6) {
        std::cout
            << style4
            << "\nenter how many times to keep track of the time (1 - 1000)"
            << end1 << std::endl;
        std::cin >> N;
        if (count == 5) *n = N;
      } else {
        N = *n;
      }

      if (N > 0 && N < 1001) {
        TsmResult rez;
        std::chrono::system_clock::time_point ta1 =
            std::chrono::high_resolution_clock::now();

        for (int i = 0; i < N; ++i) {
          if (i == 0) {
            rez = algo_.SolveTravelingSalesmanProblem(matrix_, multithreading);
          } else {
            TsmResult tmp =
                algo_.SolveTravelingSalesmanProblem(matrix_, multithreading);
            if (tmp.distance < rez.distance) {
              rez = tmp;
            }
          }
        }
        std::chrono::high_resolution_clock::time_point ta2 =
            std::chrono::high_resolution_clock::now();
        auto timeAnt =
            std::chrono::duration_cast<std::chrono::microseconds>(ta2 - ta1)
                .count();
        if (multithreading == false) {
          std::cout << style1 << "Run Ant algorithm without multithreading "
                    << end1 << std::endl;
        } else {
          std::cout << style1 << "Run Ant algorithm with multithreading "
                    << end1 << std::endl;
        }

        std::cout << style1 << "TIME: " << ((double)timeAnt / 1000000.0) << end1
                  << std::endl;
        std::cout << style1 << "PATH: ";
        PrintVector(rez.visit);
        std::cout << style1 << "DISTANCE: " << rez.distance << end1
                  << std::endl;
        std::cout << std::endl;
        if (count != 5) Menu();
      } else {
        std::cout << style4 << "\nINCORRECT COUNT OF REPS!" << end1
                  << std::endl;
        std::cout.flush();
        std::cout.clear();
      }
    } catch (std::exception &e) {
      if (count != 5) {
        std::cout << style4 << "\nError in Ant Foo: " << e.what() << std::endl;
        std::cout.flush();
        std::cout.clear();
        Menu();
      }
    }
  } else {
    if (count != 5) {
      std::cout << style4 << "\nError in Ant Foo: Matrix not loaded, try again "
                << end << std::endl;
      std::cout.flush();
      std::cout.clear();
      Menu();
    }
  }
}

template <typename T>
void ConsoleInterface::PrintVector(std::vector<T> path) {
  for (auto it : path) {
    std::cout << style3 << it << " ";
  }
  std::cout << end1 << std::endl;
}

template <typename T>
void ConsoleInterface::PrintMatrix(S21Matrix<T> res) {
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
