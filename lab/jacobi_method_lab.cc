#include <math.h>
#include <stdio.h>

#include <cmath>
#include <climits>
#include <iostream>

#include "include/command_line_parser.h"

enum Border { TOP, RIGHT, BOTTOM, LEFT };
enum Filler { ZEROS };

static const double kLeftBorder = 0;
static const double kRightBorder = 1;
static const double kBottomBorder = 0;
static const double kTopBorder = 1;

inline double GetRobustValue(double x, double y);

inline double GetExternalHeat(double x, double y);

inline double GetBorderCondition(Border border, double coordinate);

void PrintAbout();

void Fill(double** numerical_func, double** external_heat,
          int n_intervals_by_x, int n_intervals_by_y,
          Filler numerical_func_fill_policy);

double NextIteration(double** numerical_func, double** external_heat,
                     int n_intervals_by_x, int n_intervals_by_y);

double SolveAccuracy(double** numerical_func, int n_intervals_by_x,
                     int n_intervals_by_y);

int main(int argc, char** argv) {
  CommandLineParser parser(argc, argv);
  if (parser.Exists("h") || argc == 1) {
    PrintAbout();
    return 0;
  }
  if (!parser.Exists("n")) {
    std::cout << "Set number of intervals [-n], [-h] for help" << std::endl;
    return 0;
  }
  if (!parser.Exists("m")) {
    std::cout << "Set number of intervals [-m], [-h] for help" << std::endl;
    return 0;
  }
  if (!parser.Exists("iters") && !parser.Exists("eps")) {
    std::cout << "Set one of stopping criterion, [-h] for help" << std::endl;
    return 0;
  }
  const int n_intervals_by_x = parser.Get<int>("n");
  const int n_intervals_by_y = parser.Get<int>("m");
  double** numerical_func = new double*[n_intervals_by_x + 1];
  double** external_heat = new double*[n_intervals_by_x + 1];
  for (int x = 0; x <= n_intervals_by_x; ++x) {
    numerical_func[x] = new double[n_intervals_by_y + 1];
    external_heat[x] = new double[n_intervals_by_y + 1];
  }

  Fill(numerical_func, external_heat,
       n_intervals_by_x, n_intervals_by_y, ZEROS);

  int max_iters = INT_MAX;
  double target_method_accuracy = 0;
  if (parser.Exists("iters")) {
    max_iters = parser.Get<int>("iters");
  }
  if (parser.Exists("eps")) {
    target_method_accuracy = parser.Get<double>("eps");
  }

  int iter;
  double method_accuracy;
  for (iter = 0; iter < max_iters; ++iter) {
    method_accuracy = NextIteration(numerical_func, external_heat,
                                    n_intervals_by_x, n_intervals_by_y);
    if (method_accuracy <= target_method_accuracy) {
      break;
    }
  }
  printf("Iterations elapsed: %d\nMethod accuracy: %e\nmax |u - v|: %e\n\n",
         iter, method_accuracy, SolveAccuracy(numerical_func,
                                              n_intervals_by_x,
                                              n_intervals_by_y));

  for (int x = 0; x <= n_intervals_by_x; ++x) {
    delete[] numerical_func[x];
    delete[] external_heat[x];
  }
  delete[] numerical_func;
  delete[] external_heat;
}

double NextIteration(double** numerical_func, double** external_heat,
                     int n_intervals_by_x, int n_intervals_by_y) {
  double old_values[n_intervals_by_x + 1][n_intervals_by_y + 1];
  for (int x = 0; x <= n_intervals_by_x; ++x) {
    for (int y = 0; y <= n_intervals_by_y; ++y) {
      old_values[x][y] = numerical_func[x][y];
    }
  }

  double quad_step_x = pow((kRightBorder - kLeftBorder) / n_intervals_by_x, 2);
  double quad_step_y = pow((kTopBorder - kBottomBorder) / n_intervals_by_y, 2);
  double denominator = 2 * (1.0 / quad_step_x + 1.0 / quad_step_y);
  double distance_from_prev_iter = 0;
  for (int x = 1; x < n_intervals_by_x; ++x) {
    for (int y = 1; y < n_intervals_by_y; ++y) {
      numerical_func[x][y] =
          ((old_values[x - 1][y] + old_values[x + 1][y]) / quad_step_x +
          (old_values[x][y - 1] + old_values[x][y + 1]) / quad_step_y -
          external_heat[x][y]) / denominator;
      distance_from_prev_iter += pow(numerical_func[x][y] -
                                     old_values[x][y], 2);
    }
  }
  return sqrt(distance_from_prev_iter);
}

void Fill(double** numerical_func, double** external_heat,
          int n_intervals_by_x, int n_intervals_by_y,
          Filler numerical_func_fill_policy) {
  double step_x = (kRightBorder - kLeftBorder) / n_intervals_by_x;
  double step_y = (kTopBorder - kBottomBorder) / n_intervals_by_y;

  // Fill borders condition.
  for (int x = 0; x <= n_intervals_by_x; ++x) {
    numerical_func[x][0] = GetBorderCondition(BOTTOM, x * step_x);
    numerical_func[x][n_intervals_by_y] =
        GetBorderCondition(TOP, x * step_x);
  }
  for (int y = 0; y <= n_intervals_by_y; ++y) {
    numerical_func[0][y] = GetBorderCondition(LEFT, y * step_y);
    numerical_func[n_intervals_by_x][y] =
        GetBorderCondition(RIGHT, y * step_y);
  }

  switch (numerical_func_fill_policy) {
    case ZEROS:
      for (int x = 1; x < n_intervals_by_x; ++x) {
        for (int y = 1; y < n_intervals_by_y; ++y) {
          numerical_func[x][y] = 0;
        }
      }
      break;
    default: break;
  }

  // Fill external heat values.
  for (int x = 0; x <= n_intervals_by_x; ++x) {
    for (int y = 0; y <= n_intervals_by_y; ++y) {
      external_heat[x][y] = GetExternalHeat(x * step_x, y * step_y);
    }
  }
}

double SolveAccuracy(double** numerical_func, int n_intervals_by_x,
                     int n_intervals_by_y) {
  double max_diff = 0;
  double step_x = (kRightBorder - kLeftBorder) / n_intervals_by_x;
  double step_y = (kTopBorder - kBottomBorder) / n_intervals_by_y;
  for (int x = 1; x < n_intervals_by_x; ++x) {
    for (int y = 1; y < n_intervals_by_y; ++y) {
      double diff = std::abs(numerical_func[x][y] - GetRobustValue(x * step_x,
                                                                   y * step_y));
      if (diff > max_diff) {
        max_diff = diff;
      }
    }
  }
  return max_diff;
}

inline double GetRobustValue(double x, double y) {
  return exp(pow(sin(M_PI * x * y), 2));
}

inline double GetExternalHeat(double x, double y) {
  return M_PI * M_PI * (x * x + y * y) * GetRobustValue(x, y) *
         (pow(sin(2 * M_PI * x * y), 2) + 2 * cos(2 * M_PI * x * y));
}

inline double GetBorderCondition(Border border, double coordinate) {
  switch (border) {
    case TOP: case RIGHT: return exp(pow(sin(M_PI * coordinate), 2));
    case BOTTOM: case LEFT: return 1;
    default: return 0;
  }
}

void PrintAbout() {
  std::cout << "Heat equation:\n"
               "div(grad( u(x, y) )) = -f(x, y)\n"
               "x in (0, 1)\n"
               "y in (0, 1)\n"
               "Equation solve is u(x, y) = exp( sin^2(pi * x * y) )\n\n"
               "Experiment №1: Jacobi method for solving "
               "finite difference scheme\n"
               "Method parameters:\n"
               "[-n] - number of intervals by x axis\n"
               "[-m] - number of intervals by y axis\n"
               "[-iters] - stopping criterion by number of iterations\n"
               "[-eps] - stopping criterion by accuracy\n"
            << std::endl;
}
