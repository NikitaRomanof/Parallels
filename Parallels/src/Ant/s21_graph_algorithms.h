#pragma once

#include <random>
#include <set>
#include <vector>

#include "../s21_matrix_oop.h"

struct TsmResult {
  std::vector<int> visit;
  double distance;
};

class Ant {
 public:
  explicit Ant(int start_vertex) { tabu.insert(start_vertex); }
  std::set<int> tabu;

  double rand() {
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<double> distr(0, 1);
    return distr(eng);
  }
};

class GraphAlgorithms {
 public:
  //  Functions

  TsmResult SolveTravelingSalesmanProblem(const S21Matrix<double> &matrix,
                                          bool multithreading = true);
  int RandomVal(int val);

 private:
  S21Matrix<double> OneAntPath(const S21Matrix<double> &matrix,
                               const S21Matrix<double> &count_feromone,
                               int cur_vertex);
  std::vector<double> CalculateProbability(const S21Matrix<double> &matrix,
                                           int cur_vertex, double *sum_wish,
                                           const S21Matrix<double> &count_feromone,
                                           Ant ant);
  double CheckDistance(const S21Matrix<double> &matrix,
                       const double &cur_distance,
                       const S21Matrix<double> &count_feromone,
                       std::vector<int> *path);
};
