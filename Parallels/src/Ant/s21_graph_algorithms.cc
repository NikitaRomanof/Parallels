#include "s21_graph_algorithms.h"

#include <cmath>
#include <future>
#include <thread>

int GraphAlgorithms::RandomVal(int val) {
  std::random_device rd;
  std::default_random_engine eng(rd());
  std::uniform_int_distribution<int> distr(0, val);
  return distr(eng);
}

TsmResult GraphAlgorithms::SolveTravelingSalesmanProblem(
    const S21Matrix<double> &matrix, bool multithreading) {
  if (matrix.get_rows() < 2) {
    throw std::invalid_argument(
        "Error on solveTravelingSalesmanProblem. Graph < 2");
  }
  S21Matrix<double> count_feromone(matrix);
  count_feromone.fillMatrix(0.2);
  S21Matrix<double> deltaTau(matrix);
  deltaTau.fillMatrix(0.0);
  std::vector<int> path;
  double cur_distance = 0.0;
  const int CYCLE_STEPS = 40, ANT_COUNTS = 200;
  const double VAPE = 0.5;

  for (int count = 0; count <= CYCLE_STEPS; ++count) {
    if (multithreading == true) {
      for (int stepsAnt = 0; stepsAnt < ANT_COUNTS;) {
        std::vector<std::thread> threadPool;
        for (unsigned int tr = 0; tr < std::thread::hardware_concurrency();
             ++tr) {
          std::thread t1([&] {
            deltaTau += OneAntPath(matrix, count_feromone,
                                   RandomVal(matrix.get_rows() - 1));
          });
          threadPool.push_back(std::move(t1));
        }
        for (unsigned int tr = 0; tr < std::thread::hardware_concurrency();
             ++tr) {
          threadPool[tr].join();
        }
        stepsAnt += std::thread::hardware_concurrency();
      }
    } else {
      for (int stepsAnt = 0; stepsAnt < ANT_COUNTS; ++stepsAnt) {
        deltaTau += OneAntPath(matrix, count_feromone,
                               (RandomVal((matrix.get_rows() - 1))));
      }
    }
    count_feromone *= VAPE;
    count_feromone += deltaTau;
    deltaTau.fillMatrix(0.0);
    if (count == CYCLE_STEPS) {
      double newDistance =
          CheckDistance(matrix, cur_distance, count_feromone, &path);
      if (newDistance < 0.001) {
        count = 0;
      } else {
        cur_distance = newDistance;
      }
    }
  }
  TsmResult rez{path, cur_distance};
  return rez;
}

S21Matrix<double> GraphAlgorithms::OneAntPath(
    const S21Matrix<double> &matrix, const S21Matrix<double> &count_feromone,
    int cur_vertex) {
  const double Q = 10.0;
  S21Matrix<double> deltaTau(matrix);
  deltaTau.fillMatrix(0.0);
  Ant ant(cur_vertex);
  while ((int)ant.tabu.size() < matrix.get_rows()) {
    ant.tabu.insert(cur_vertex);
    double sum_wish = 0.0;
    std::vector<double> probability =
        CalculateProbability(matrix, cur_vertex, &sum_wish, count_feromone, ant);
    int tmpVertex = cur_vertex;
    double step = ant.rand();
    for (int n = 0; n < (int)probability.size(); ++n) {
      if (ant.tabu.count(n) == 0) {
        step -= probability[n];
        if (step < 0) {
          cur_vertex = n;
          break;
        }
      }
    }
    if (cur_vertex == tmpVertex && (int)ant.tabu.size() < matrix.get_rows()) {
      throw std::invalid_argument("Error on OneAntPath. Graph is incomplete");
    }
    if (matrix(tmpVertex, cur_vertex) > 0) {
      deltaTau(tmpVertex, cur_vertex) += Q / matrix(tmpVertex, cur_vertex);
    }
  }
  return deltaTau;
}

std::vector<double> GraphAlgorithms::CalculateProbability(
    const S21Matrix<double> &matrix, int cur_vertex, double *sum_wish,
    const S21Matrix<double> &count_feromone, Ant ant) {
  const int ALPHA = 1, BETTA = 1;
  std::vector<double> accesableWish;
  for (int j = 0; j < matrix.get_rows(); ++j) {
    double n = 0.0;
    double tau = 0.0;
    if (matrix(cur_vertex, j) > 0 && ant.tabu.count(j) == 0) {
      n = 1.0 / matrix(cur_vertex, j);
      tau = count_feromone(cur_vertex, j);
      *sum_wish += (pow(tau, ALPHA) * pow(n, BETTA));
    }
    accesableWish.push_back((pow(tau, ALPHA) * pow(n, BETTA)));
  }
  std::vector<double> probability;
  for (int k = 0; k < matrix.get_rows(); ++k) {
    if (matrix(cur_vertex, k) > 0.0 && ant.tabu.count(k) == 0) {
      probability.push_back((accesableWish[k] / *sum_wish));
    } else {
      probability.push_back(0.0);
    }
  }

  return probability;
}

double GraphAlgorithms::CheckDistance(const S21Matrix<double> &matrix,
                                      const double &cur_distance,
                                      const S21Matrix<double> &count_feromone,
                                      std::vector<int> *path) {
  int iter = 0;
  std::set<int> check;
  double newDistance = 0.0;
  path->push_back(iter);
  check.insert(iter);
  while ((int)path->size() < matrix.get_rows()) {
    int tmpIter = 0;
    double max = 0.0;
    for (int z = 0; z < matrix.get_rows(); ++z) {
      if (max < count_feromone(iter, z) && check.count(z) == 0) {
        max = count_feromone(iter, z);
        tmpIter = z;
      }
    }
    newDistance += matrix(iter, tmpIter);
    iter = tmpIter;
    path->push_back(iter);
    check.insert(iter);
  }
  if (fabs(newDistance - cur_distance) < 0.9) {
    newDistance = 0.00000;
    path->clear();
  } else {
    int tmp = path->operator[](path->size() - 1);
    newDistance += matrix(tmp, 0);
    path->push_back(0);
  }
  return newDistance;
}
