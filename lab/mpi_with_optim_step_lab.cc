#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <iostream>
#include <iomanip>
#include <vector>
#include <climits>
#include <cfloat>
#include <cmath>
#include <sstream>

#include "include/command_line_parser.h"
#include "include/table_printer.h"
#include "include/dirichlet_task.h"

enum Task { TEST, MAIN };

static const double kLeftBorder = 0.0;
static const double kRightBorder = 1.0;
static const double kBottomBorder = 0.0;
static const double kTopBorder = 1.0;

double GetExternalHeat(double x, double y, Task task);

double GetBorderCondition(Border border, double coordinate, Task task);

// eps - max|x[s+1]-x[s]| element wise.
void Solve(DirichletTask* dirichlet_task, std::vector<double>& result,
           int& n_processed_iters, double& achieved_eps,
           int max_n_iters = 1, double target_eps = 0);

void Print(int n_intervals_by_x, int n_intervals_by_y,
           std::vector<double>& robust_values,
           std::vector<double>& result, Task task,
           int n_processed_iters, double achieved_eps, bool print_tables);

DirichletTask* BuildDirichletTask(Task task, int n, int m);

void PrintAbout();

int main(int argc, char** argv) {
  CommandLineParser parser(argc, argv);
  if (parser.Exists("h") || argc == 1) {
    PrintAbout();
    return 0;
  }
  if (!parser.Exists("n")) {
    std::cout << "Set number of intervals by x [-n], [-h] for help"
              << std::endl;
    return 0;
  }
  if (!parser.Exists("m")) {
    std::cout << "Set number of intervals by y [-m], [-h] for help"
              << std::endl;
    return 0;
  }
  if (!parser.Exists("iters") && !parser.Exists("eps")) {
    std::cout << "Set one of stopping criterion, [-h] for help" << std::endl;
    return 0;
  }

  const int n = parser.Get<int>("n");
  const int m = parser.Get<int>("m");
  const int n_iters = parser.Get<int>("iters", INT_MAX);
  const double eps = parser.Get<double>("eps", 0);
  const Task task = (parser.Exists("main") ? MAIN : TEST);
  const bool print_tables = !parser.Exists("s");

  DirichletTask* dirichlet_task = BuildDirichletTask(task, n, m);

  std::vector<double> robust_values;
  std::vector<double> result;
  int n_processed_iters;
  double achieved_eps;

  Solve(dirichlet_task, result, n_processed_iters, achieved_eps, n_iters, eps);
  delete dirichlet_task; dirichlet_task = 0;

  const double h = (kRightBorder - kLeftBorder) / n;
  const double k = (kTopBorder - kBottomBorder) / m;
  if (task == TEST) {
    for (int j = 0; j <= m; ++j) {
      const double y = j * k;
      for (int i = 0; i <= n; ++i) {
        robust_values.push_back(exp(pow(sin(M_PI * i * h * y), 2)));
      }
    }
    Print(n, m, robust_values, result, TEST, n_processed_iters, achieved_eps,
          print_tables);
  } else {
    std::vector<double> result_on_dense_net;
    int n_processed_iters_dense_net;
    double achieved_eps_dense_net;

    dirichlet_task = BuildDirichletTask(MAIN, 2 * n, 2 * m);
    Solve(dirichlet_task, result_on_dense_net, n_processed_iters,
          achieved_eps_dense_net, n_iters, eps);
    delete dirichlet_task; dirichlet_task = 0;

    std::vector<double> extracted_values((n + 1) * (m + 1), 0);
    for (int j = 1; j < m; ++j) {
      const double y = kBottomBorder + j * k;
      extracted_values[j * (n + 1)] = GetBorderCondition(LEFT, y, MAIN);
      extracted_values[j * (n + 1) + n] = GetBorderCondition(RIGHT, y, MAIN);
      for (int i = 1; i < n; ++i) {
        const int idx = (2 * j - 1) * (2 * n - 1) + 2 * i - 1;
        extracted_values[j * (n + 1) + i] = result_on_dense_net[idx];
      }
    }
    for (int i = 0; i <= n; ++i) {
      const double x = kLeftBorder + i * h;
      extracted_values[i] = GetBorderCondition(BOTTOM, x, MAIN);
      extracted_values[(n + 1) * m + i] = GetBorderCondition(TOP, x, MAIN);
    }
    Print(n, m, extracted_values, result, MAIN, n_processed_iters,
          achieved_eps, print_tables);
  }
}

DirichletTask* BuildDirichletTask(Task task, int n, int m) {
  // Fill external heat for task.
  const double h = (kRightBorder - kLeftBorder) / n;
  const double k = (kTopBorder - kBottomBorder) / m;
  double* external_heat = new double[(m - 1) * (n - 1)];
  for (int j = 1; j < m; ++j) {
    const int offset = (j - 1) * (n - 1);
    const double y = kBottomBorder + j * k;
    for (int i = 1; i < n; ++i) {
      const double x = kLeftBorder + i * h;
      external_heat[offset + i - 1] = GetExternalHeat(x, y, task);
    }
  }

  DirichletTask* dirichlet_task = new DirichletTask(kLeftBorder, kRightBorder,
                                                    kTopBorder, kBottomBorder,
                                                    n, m, external_heat);
  delete[] external_heat; external_heat = 0;

  // Set borders.
  double* first_border = new double[n - 1];
  double* second_border = new double[n - 1];
  for (int i = 1; i < n; ++i) {
    const double x = kLeftBorder + i * h;
    first_border[i - 1] = GetBorderCondition(BOTTOM, x, task);
    second_border[i - 1] = GetBorderCondition(TOP, x, task);
  }
  dirichlet_task->UpdateBorder(BOTTOM, first_border);
  dirichlet_task->UpdateBorder(TOP, second_border);
  delete[] first_border;
  delete[] second_border;

  first_border = new double[m - 1];
  second_border = new double[m - 1];
  for (int j = 1; j < m; ++j) {
    const double y = kBottomBorder + j * k;
    first_border[j - 1] = GetBorderCondition(LEFT, y, task);
    second_border[j - 1] = GetBorderCondition(RIGHT, y, task);
  }
  dirichlet_task->UpdateBorder(LEFT, first_border);
  dirichlet_task->UpdateBorder(RIGHT, second_border);
  delete[] first_border; first_border = 0;
  delete[] second_border; second_border = 0;

  return dirichlet_task;
}

void Solve(DirichletTask* dirichlet_task, std::vector<double>& result,
           int& n_processed_iters, double& achieved_eps,
           int max_n_iters, double target_eps) {
  achieved_eps = DBL_MAX;
  for (n_processed_iters = 0;
       n_processed_iters < max_n_iters && achieved_eps > target_eps;
       ++n_processed_iters) {
    dirichlet_task->MPIIteration(achieved_eps);

    printf("\rProcessed iterations: %d, max|x[s+1]-x[s]| = %e",
           n_processed_iters + 1, achieved_eps);
    fflush(stdout);
  }
  std::cout << std::endl;

  dirichlet_task->GetState(result);
}

double GetBorderCondition(Border border, double coordinate, Task task) {
  switch (task) {
    case TEST:
      switch (border) {
        case LEFT: case BOTTOM: return 1;
        case RIGHT: case TOP: return exp(pow(sin(M_PI * coordinate), 2));
        default: return 0;
      }
    case MAIN:
      switch (border) {
        case LEFT: case RIGHT: return sin(M_PI * coordinate);
        case BOTTOM: case TOP: return coordinate * (1 - coordinate);
        default: return 0;
      }
    default: return 0;
  }
}

double GetExternalHeat(double x, double y, Task task) {
  switch (task) {
    case TEST: {
      double term = cos(2 * M_PI * x * y);
      return M_PI * M_PI * (x * x + y * y) * exp(pow(sin(M_PI * x * y), 2)) *
             (1 - pow(term, 2) + 2 * term);
    }
    case MAIN: return -pow(sin(M_PI * x * y), 2);
    default: return 0;
  }
}

void PrintAbout() {
  std::cout << "Iterative method for solving equation:\n"
             "div(grad( u(x, y) )) = -f(x, y)\n"
             "x in [0, 1]\n"
             "y in [0, 1]\n\n"

             "Test task [--test]:\n"
             "f(x, y) = -U*pi*pi*(xx+yy)*(sin^2(2pi*xy)+2*cos(2pi*xy))\n"
             "U(0, y) = 1;                  U(x, 0) = 1\n"
             "U(1, y) = exp( sin^2(pi*y) ); U(x, 1) = exp( sin^2(pi*x) )\n"
             "Equation solve is u(x, y) = exp( sin^2(pi*xy) )\n\n"

             "Main task [--main]:\n"
             "f(x, y) = sin^2(pi*xy)\n"
             "U(0, y) = sin(pi*y); U(x, 0) = x-x^2\n"
             "U(1, y) = sin(pi*y); U(x, 1) = x-x^2\n\n"

             "Method parameters:\n"
             "[-n] - number of intervals by x axis\n"
             "[-m] - number of intervals by y axis\n"
             "[-iters] - stopping criterion by number of iterations\n"
             "[-eps] - stopping criterion by accuracy\n"
             "[--s] - silence mode without tables printing"
          << std::endl;
}

void Print(int n_intervals_by_x, int n_intervals_by_y,
           std::vector<double>& robust_values,
           std::vector<double>& result, Task task,
           int n_processed_iters, double achieved_eps, bool print_tables) {
  const int n = n_intervals_by_x;
  const int m = n_intervals_by_y;
  const double h = (kRightBorder - kLeftBorder) / n;
  const double k = (kTopBorder - kBottomBorder) / m;

  std::cout << "\nNet step by x: " << h << std::endl;
  std::cout << "Net step by y: " << k << std::endl;
  
  // Extract accuracy.
  double max_diff = 0;
  int argmax_i = 0;
  int argmax_j = 0;
  int idx = 0;
  for (int j = 1; j < m; ++j) {
    for (int i = 1; i < n; ++i) {
      double diff = fabs(result[idx] - robust_values[j * (n + 1) + i]);
      ++idx;
      if (diff > max_diff) {
        max_diff = diff;
        argmax_i = i;
        argmax_j = j;
      }
    }
  }

  printf("Number of processed iterations: %d\n", n_processed_iters);
  printf("max|x[s+1]-x[s]| = %e\n", achieved_eps);
  std::cout << "max|V-" << (task == TEST ? "U" : "V2") << "| = " << std::flush;
  printf("%e (at x[%d]=%f, y[%d]=%f)\n", max_diff, argmax_i, h * argmax_i,
         argmax_j, k * argmax_j);

  if (!print_tables) return;

  // Robust values.
  std::vector<std::vector<std::string> > data(m + 2);
  for (int i = 0; i < m + 2; ++i) {
    data[i].resize(n + 2);
  }
  data[0][0] = (task == TEST? "Robust values" : "Results on V2");
  for (int i = 0; i <= n; ++i) {
    std::ostringstream ss;
    ss << i * h;
    data[0][i + 1] = ss.str();

    for (int j = 0; j <= m; ++j) {
      std::ostringstream ss;
      ss << robust_values[(m - j) * (n + 1) + i];
      data[j + 1][i + 1] = ss.str();
    }
  }
  for (int j = 0; j <= m; ++j) {
    std::ostringstream ss;
    ss << (m - j) * k;
    data[j + 1][0] = ss.str();
  }
  TablePrinter::Print(data);

  // Results.
  data[0][0] = "Results";
  data.erase(data.end());
  data.erase(data.begin() + 1);
  for (int i = 0; i < data.size(); ++i) {
    data[i].erase(data[i].end());
    data[i].erase(data[i].begin() + 1);
  }
  for (int i = 0; i < n - 1; ++i) {
    for (int j = 0; j < m - 1; ++j) {
      std::ostringstream ss;
      ss << result[(m - 2 - j) * (n - 1) + i];
      data[j + 1][i + 1] = ss.str();
    }
  }
  TablePrinter::Print(data);

  // Accuracy.
  data[0][0] = "Differences";
  for (int i = 0; i < n - 1; ++i) {
    for (int j = 0; j < m - 1; ++j) {
      std::ostringstream ss;
      ss << std::scientific
         << fabs(result[(m - 2 - j) * (n - 1) + i] -
                 robust_values[(m - j - 1) * (n + 1) + i + 1]);
      data[j + 1][i + 1] = ss.str();
    }
  }
  TablePrinter::Print(data);
}
