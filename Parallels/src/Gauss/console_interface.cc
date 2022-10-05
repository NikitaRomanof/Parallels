#include "console_interface.h"

#include <chrono>

void ConsoleInterfaceGauss::Start() {
  std::cout << style1
            << "\n*******************************************************"
            << end << std::endl;
  std::cout << style2
            << "                Welcome to parallels                     "
            << end << std::endl;
  std::cout << style1
            << "*******************************************************" << end
            << std::endl;

  matrixGenerateOrLOad_ = false;
  Menu();
}

void ConsoleInterfaceGauss::Read() {
  while (true) {
    std::cout << style3 << " Enter file name           " << end << std::endl;
    std::cout << std::endl;
    std::string str;
    std::cin >> str;
    try {
      gauss.GetLinearSystemFromFile("Gauss/" + str);
      matrixGenerateOrLOad_ = true;
      break;
    } catch (std::exception &e) {
      std::cout << style4 << "\nError: " << e.what() << end << "\n\n";
      std::cout.flush();
      std::cout.clear();
      return;
    }
  }
  size_ = gauss.GetSize();
  gauss.PrintMatrix("Linear system from file");
}

void ConsoleInterfaceGauss::GenerateRandomMatrix() {
  std::cout << "Which size of linear system?" << std::endl;
  std::cin >> size_;
  gauss.Clear();
  try {
    gauss.GenerateLinearSystem(size_);
    gauss.PrintMatrix("Generated linear system");
  } catch (MatrixException &e) {
    std::cout << e.what() << std::endl;
  }
}

void ConsoleInterfaceGauss::Menu() {
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
      << " 1. ----> Generate random linear system                          "
      << end << std::endl;
  std::cout
      << style3
      << " 2. ----> Load linear system from file                           "
      << end << std::endl;
  std::cout
      << style3
      << " 3. ----> Run Gauss algorithm without multithreading             "
      << end << std::endl;
  std::cout
      << style3
      << " 4. ----> Run Gauss algorithm with multithreading                "
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
  int N;
  if (choice == 1) {
    GenerateRandomMatrix();
    Menu();
  } else if (choice == 2) {
    Read();
    Menu();
  } else if (choice == 3) {
    if (gauss.IsEmpty() == false) {
      N = AskN();
      SolveGauss(false, N);
    } else {
      std::cout << "Linear system is unavailable\n";
    }
    Menu();
  } else if (choice == 4) {
    if (gauss.IsEmpty() == false) {
      N = AskN();
      SolveGauss(true, N);
    } else {
      std::cout << "Linear system is unavailable\n";
    }
    Menu();
  } else if (choice == 5) {
    if (gauss.IsEmpty() == false) {
      N = AskN();
      SolveGauss(false, N);
      SolveGauss(true, N);
    } else {
      std::cout << "Linear system is unavailable\n";
    }
    Menu();
  } else if (choice == 0) {
    // return;
  } else {
    std::cout << style4 << "\nError: incorrect input, try again " << end
              << std::endl;
    std::cout.flush();
    std::cout.clear();
  }
}

int ConsoleInterfaceGauss::AskN() {
  int N;
  std::cout << style4
            << "\nenter how many times to keep track of the time (1 - 1000)"
            << end1 << std::endl;
  std::cin >> N;
  return N;
}

void ConsoleInterfaceGauss::SolveGauss(bool multithreading, int N) {
  if (N > 0 && N < 1001) {
    std::chrono::system_clock::time_point ta1 =
        std::chrono::high_resolution_clock::now();
    S21Matrix<double> *res;
    if (multithreading == false) {
      for (int i = 0; i < N; ++i) {
        res = gauss.SolveLinearSystemSingleThread();
      }
    } else {
      for (int i = 0; i < N; ++i) {
        res = gauss.SolveLinearSystemMultiThread();
      }
    }
    std::chrono::high_resolution_clock::time_point ta2 =
        std::chrono::high_resolution_clock::now();
    auto time_gauss =
        std::chrono::duration_cast<std::chrono::microseconds>(ta2 - ta1)
            .count();
    if (multithreading == true) {
      std::cout << style1 << "Gauss algorithm time with multithreading " << end1
                << std::endl;
    } else {
      std::cout << style1 << "Gauss algorithm time without multithreading "
                << end1 << std::endl;
    }

    std::cout << style1 << "TIME: " << ((double)time_gauss / 1000000.0) << end1
              << std::endl;
    std::cout << std::endl;
    std::cout << style1 << "ROOTS: " << end1 << std::endl;
    PrintRoots(*res);
    std::cout << std::endl;
  }
}

void ConsoleInterfaceGauss::PrintRoots(S21Matrix<double> roots) {
  for (int i = 0; i < size_; ++i) {
    std::cout << style3 << roots(0, i) << " ";
  }
  std::cout << end1 << std::endl;
}
