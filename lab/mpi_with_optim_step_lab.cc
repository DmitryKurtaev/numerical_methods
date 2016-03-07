#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <vector>
#include <climits>
#include <cfloat>
#include <cmath>
#include <sstream>

#include "include/command_line_parser.h"
#include "include/table_printer.h"

enum Border { TOP, RIGHT, BOTTOM, LEFT };
enum Task { TEST, MAIN };

static const double kLeftBorder = 0.0;
static const double kRightBorder = 1.0;
static const double kBottomBorder = 0.0;
static const double kTopBorder = 1.0;

double GetExternalHeat(double x, double y, Task task);

double GetBorderCondition(Border border, double coordinate, Task task);

void Solve(int n_intervals_by_x, int n_intervals_by_y,
           const std::vector<double>& robust_values,
           std::vector<double>& result, Task task,
           int n_iters = INT_MAX, double eps = DBL_MAX);

void Print(int n_intervals_by_x, int n_intervals_by_y,
           std::vector<double>& robust_values,
           std::vector<double>& result, Task task);

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

  const int n_intervals_by_x = parser.Get<int>("n");
  const int n_intervals_by_y = parser.Get<int>("m");
  Task task = (parser.Exists("main") ? MAIN : TEST);

  std::vector<double> robust_values;
  if (task == TEST) {
    const int n = n_intervals_by_y;
    const int m = n_intervals_by_y;
    const double h = (kRightBorder - kLeftBorder) / n;
    const double k = (kTopBorder - kBottomBorder) / m;
    for (int j = 1; j < n_intervals_by_y - 1; ++j) {
      const double y = j * k;
      for (int i = 1; i < n_intervals_by_x - 1; ++i) {
        robust_values.push_back(exp(pow(sin(M_PI * i * h * y), 2)));
      }
    }
    Print(n_intervals_by_x, n_intervals_by_y, robust_values, robust_values,
          TEST);
  }
  // Solve(n_intervals_by_x, n_intervals_by_y, )
}

void Solve(int n_intervals_by_x, int n_intervals_by_y,
           const std::vector<double>& robust_values,
           std::vector<double>& result, Task task,
           int n_iters, double eps) {
  const int n = n_intervals_by_y;
  const int m = n_intervals_by_y;
  const double h = (kRightBorder - kLeftBorder) / n;
  const double k = (kTopBorder - kBottomBorder) / m;
  const double inv_h_quad = 1.0 / h * h;
  const double inv_k_quad = 1.0 / k * k;
  const int dim = (n - 1) * (m - 1);
  const double step = 1.0 / ((2 + sin(M_PI / n) - cos(M_PI / n)) * inv_h_quad +
                             (2 + sin(M_PI / m) - cos(M_PI / m)) * inv_k_quad);

  double* x = new double[dim];
  memset(x, 4 * dim, 0);

  // Setup right part of equations system.
  double* b = new double[dim];
  memset(b, 4 * dim, 0);
  for (int j = 0; j < m - 1; ++j) {
    const double y = j * k;
    const int offset = j * (n - 1);
    b[offset] -= GetBorderCondition(LEFT, y, task) * inv_h_quad;
    b[offset + n - 2] -= GetBorderCondition(RIGHT, y, task) * inv_h_quad;
    for (int i = 0; i < n - 1; ++i) {
      b[offset + i] += GetExternalHeat(i * h, y, task);
    }
  }
  const int offset = (n - 1) * (m - 2); 
  for (int i = 0; i < n - 1; ++i) {
    const double x = i * h;
    b[i] -= GetBorderCondition(BOTTOM, x, task) * inv_k_quad;
    b[offset + i] -= GetBorderCondition(TOP, x, task) * inv_k_quad;
  }

  // Collect robust values.

  double accuracy = DBL_MAX;
  for (int iter = 0; iter < n_iters && (robust_values.size() != dim ||
                                        accuracy >= eps); ++iter) {
    // Do iteration.

    // |-- Bottom border--------------------------------------------------------
    //     |-- Left bottom point.
    double term = inv_k_quad * (x[n - 1] - 2 * x[0]) +
                  inv_h_quad * (x[1] - 2 * x[0]);
    x[0] += step * (b[0] - term);

    //     |-- Bottom line, center points.
    for (int i = 1; i < n - 2; ++i) {
      term = inv_k_quad * (x[n - 1 + i] - 2 * x[i]) + 
             inv_h_quad * (x[i + 1] + x[i - 1] - 2 * x[i]);
      x[i] += step * (b[i] - term);
    }

    //     |-- Right bottom point.
    term = inv_k_quad * (x[2 * n - 3] - 2 * x[n - 2]) + 
           inv_h_quad * (x[n - 3] - 2 * x[n - 2]);
    x[n - 2] += step * (b[n - 2] - term);

    // |-- Center lines---------------------------------------------------------
    for (int j = 1; j < m - 2; ++j) {
      const int offset = j * (n - 1);

      //   |-- Left border.
      double term = inv_k_quad * (x[offset - n + 1] + x[offset + n - 1]) +
                    inv_h_quad * x[offset + 1] +
                    2 * (inv_h_quad + inv_k_quad) * x[offset];
      x[offset] += step * (b[offset] - term);

      //   |-- Centers.
      for (int i = 1; i < n - 2; ++i) {
        term = inv_k_quad * x[offset - n + 1 + i] + 
               inv_k_quad * x[offset + n - 1 + i] + 
               inv_h_quad * x[offset + i + 1] +
               inv_h_quad * x[offset + i - 1] -
               2 * (inv_h_quad + inv_k_quad) * x[offset + i];
        x[offset + i] += step * (b[offset + i] - term);
      }

      //   |-- Right border.
      term = inv_k_quad * (x[offset - 1] + x[offset + 2 * n - 3]) +
             inv_h_quad * x[offset + n - 3] -
             2 * (inv_h_quad + inv_k_quad) * x[offset + n - 2];
      x[offset + n - 2] += step * (b[offset + n - 2] - term);
    }
    
    // |-- Top border-----------------------------------------------------------
    int offset = (m - 2) * (n - 1);

    //   |-- Left top point.
    term = inv_k_quad * x[offset - n + 1] + 
           inv_h_quad * x[offset + 1] +
           2 * (inv_h_quad + inv_k_quad) * x[offset];
    x[offset] += step * (b[offset] - term);

    //   |-- Centers.
    for (int i = 1; i < n - 2; ++i) {
      double term = inv_k_quad * x[offset - n + 1 + i] + 
                    inv_h_quad * x[offset + i + 1] +
                    inv_h_quad * x[offset + i - 1] -
                    2 * (inv_h_quad + inv_k_quad) * x[offset + i];
      x[offset + i] += step * (b[offset + i] - term);
    }

    //   |-- Right top point.
    offset = (n - 1) * (m - 1);
      term = inv_k_quad * x[offset - n] +
             inv_h_quad * x[offset - 2] -
             2 * (inv_h_quad + inv_k_quad) * x[offset - 1];
    x[offset - 1] += step * (b[offset - 1] - term);

    // Compute accuracy.
    if (robust_values.size() == dim) {
      accuracy = 0;
      for (int i = 0; i < dim; ++i) {
        accuracy = std::max(accuracy, fabs(x[i] - robust_values[i]));
      }
    }
  }

  delete[] b;
  delete[] x;
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
      double sin_quad = 1 - pow(term, 2);
      return M_PI * M_PI * (x * x + y * y) * exp(sin_quad) *
             (sin_quad + 2 * term);
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
             "[-eps] - stopping criterion by accuracy"
          << std::endl;
}

void Print(int n_intervals_by_x, int n_intervals_by_y,
           std::vector<double>& robust_values,
           std::vector<double>& result, Task task) {
  const int n = n_intervals_by_x;
  const int m = n_intervals_by_y;
  const double h = (kRightBorder - kLeftBorder) / n;
  const double k = (kTopBorder - kBottomBorder) / m;

  std::cout << "Net step by x: " << h << std::endl;
  std::cout << "Net step by y: " << k << std::endl;
  
  // Extract accuracy.
  double max_diff = 0;
  int argmax = 0;
  for (int i = 0; i < (n - 1) * (m - 1); ++i) {
    if (fabs(robust_values[i] - result[i]) > max_diff) {
      max_diff = fabs(robust_values[i] - result[i]);
      argmax = i;
    }
  }

  std::vector<std::vector<std::string> > data(m);
  for (int i = 0; i < m; ++i) {
    data[i].resize(n);
  }
  data[0][0] = "Robust values";
  for (int i = 0; i < n - 1; ++i) {
    std::ostringstream ss;
    ss << (i + 1) * h;
    data[0][i + 1] = ss.str();

    for (int j = 0; j < m - 1; ++j) {
      std::ostringstream ss;
      ss << robust_values[(m - 2 - j) * (n - 1) + i];
      data[j + 1][i + 1] = ss.str();
    }
  }
  for (int j = 0; j < m - 1; ++j) {
    std::ostringstream ss;
    ss << (m - 1 - j) * k;
    data[j + 1][0] = ss.str();
  }
  TablePrinter::Print(data);


    // std::cout << "max|V-" << (task == TEST ? "U" : "V2") << "| = "
    //         << max_diff << " (at x[" << argmax << "]=" h * argmax << ")"
    //         << std::endl;
}