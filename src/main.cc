#include "include/plot.h"
#include "include/cubic_spline.h"
#include <cmath>
#include <stdio.h>
#include <stdlib.h>

Plot plot;

enum Border { LEFT, RIGHT };

// Functions:
// 0. f(x) = x^3 + 3x^2, if x in [-1, 0]
//           -x^3 + 3x^2, if x in [0, 1];
// 1. f(x) = sqrt(exp(x) - 1)
// 2. f(x) = sqrt(exp(x) - 1) + cos10x
// 3. f(x) = sqrt(exp(x) - 1) + cos100x
//
// Borders:
// 0. [-1, 1]
// 1. [1, 3]
//
// Borders conditions:
// 0. [0, 0] identity
// 1. [f''(lbound), f''(rbound)] robust

double Function(int func_id, double x);

double FunctionDerivate(int func_id, double x);

double FunctionSecondDerivate(int func_id, double x);

double GetBorderCondition(int condition_id,
                          int func_id,
                          double border);

double GetBorder(int borders_id, Border border);

void GenNodes(int number,
              double lbound,
              double rbound,
              int func_id,
              std::vector<double>* x,
              std::vector<double>* y);

void CheckAccuracy(CubicSpline& spline,
                   int func_id,
                   int borders_id);

void BuildSpline(int func_id,
                 int borders_id,
                 int borders_condition_id,
                 int number_nodes,
                 CubicSpline* spline);

void Visualize(CubicSpline& spline,
               int number_nodes,
               int func_id,
               int borders_id);

int main(int argc, char** argv) {
  double tmp = 4.942921e+01 / 9.854518e+00;
  tmp = 9.854518e+00 / 2.432854e+00;
  tmp = 2.432854e+00 / 6.077543e-01;
  const int func_id = atoi(argv[1]);
  const int borders_id = atoi(argv[2]);
  const int borders_condition_id = atoi(argv[3]);
  const int number_nodes = atoi(argv[4]);

  CubicSpline spline;
  BuildSpline(func_id,
              borders_id,
              borders_condition_id,
              number_nodes,
              &spline);
  CheckAccuracy(spline,
                func_id,
                borders_id);
  Visualize(spline,
            number_nodes,
            func_id,
            borders_id);

  return 0;
}

double Function(int func_id, double x) {
  switch (func_id) {
    case 0: {
      if (x >= -1 - kZeroLimit && x <= kZeroLimit) {
        return x * x * x + 3 * x * x;
      } else {
        if (x >= -kZeroLimit && x <= 1 + kZeroLimit) {
          return -x * x * x + 3 * x * x;
        } else {
          return 0;
        }
      }
    }
    case 1: return sqrt(exp(x) - 1);
    case 2: return Function(1, x) + cos(10 * x);
    case 3: return Function(1, x) + cos(100 * x);
    default: return 0;
  }
}

double FunctionDerivate(int func_id, double x) {
  switch (func_id) {
    case 0: {
      if (x >= -1 - kZeroLimit && x <= kZeroLimit) {
        return 3 * x * x + 6 * x;
      } else {
        if (x >= -kZeroLimit && x <= 1 + kZeroLimit) {
          return -3 * x * x + 6 * x;
        } else {
          return 0;
        }
      }
    }
    case 1: return 0.5 * exp(x) / sqrt(exp(x) - 1);
    case 2: return FunctionDerivate(1, x) - 10 * sin(10 * x);
    case 3: return FunctionDerivate(1, x) - 100 * sin(100 * x);
    default: return 0;
  }
}

double FunctionSecondDerivate(int func_id, double x) {
  switch (func_id) {
    case 0: {
      if (x >= -1 - kZeroLimit && x <= kZeroLimit) {
        return 6 * x + 6;
      } else {
        if (x >= -kZeroLimit && x <= 1 + kZeroLimit) {
          return -6 * x + 6;
        } else {
          return 0;
        }
      }
    }
    case 1: return (exp(2 * x) - 2 * exp(x)) / (4 * pow(exp(x) - 1, 1.5));
    case 2: return FunctionSecondDerivate(1, x) - 100 * cos(10 * x);
    case 3: return FunctionSecondDerivate(1, x) - 1e+4 * cos(100 * x);
    default: return 0;
  }
}

double GetBorderCondition(int condition_id,
                          int func_id,
                          double border) {
  switch (condition_id) {
    case 0: return 0;
    case 1: return FunctionSecondDerivate(func_id, border);
    default: return 0;
  }
}

double GetBorder(int borders_id, Border border) {
  switch (borders_id) {
    case 0: {
      if (border == LEFT) return -1.0;
      if (border == RIGHT) return 1.0;
      return 0.0;
    }
    case 1: {
      if (border == LEFT) return 1.0;
      if (border == RIGHT) return 3.0;
      return 0.0;
    }
    default: return 0.0;
  }
}

void GenNodes(int number,
              double lbound,
              double rbound,
              int func_id,
              std::vector<double>* x,
              std::vector<double>* y) {
  double step = (rbound - lbound) / number;
  x->resize(number + 1);
  y->resize(number + 1);
  double current_x = lbound;
  for (int i = 0; i < number + 1; ++i) {
    (*x)[i] = current_x;
    (*y)[i] = Function(func_id, current_x);
    current_x += step;
  }
}

void CheckAccuracy(CubicSpline& spline,
                   int func_id,
                   int borders_id) {
  const double kAccuracyStep = 0.0001;

  double lbound = GetBorder(borders_id, LEFT);
  double rbound = GetBorder(borders_id, RIGHT);
  double max_diff_original = 0;
  double max_diff_derivate = 0;
  for (double x = lbound; x <= rbound; x += kAccuracyStep) {
    double diff_original = std::abs(spline.GetValue(x) -
                                    Function(func_id, x));
    if (diff_original > max_diff_original) {
      max_diff_original = diff_original;
    }

    double diff_derivate = std::abs(spline.GetDerivate(x) -
                                    FunctionDerivate(func_id, x));
    if (diff_derivate > max_diff_derivate) {
      max_diff_derivate = diff_derivate;
    }
  }
  printf("max|f(x) - S(x)| = %e\n", max_diff_original);
  printf("max|f'(x) - S'(x)| = %e\n", max_diff_derivate);
  fflush(stdout);
}

void BuildSpline(int func_id,
                 int borders_id,
                 int borders_condition_id,
                 int number_nodes,
                 CubicSpline* spline) {
  double lbound = GetBorder(borders_id, LEFT);
  double rbound = GetBorder(borders_id, RIGHT);

  double lbound_xx = GetBorderCondition(borders_condition_id,
                                        func_id, lbound);
  double rbound_xx = GetBorderCondition(borders_condition_id,
                                        func_id, rbound);
  std::vector<double> nodes_x;
  std::vector<double> nodes_y;
  GenNodes(number_nodes,
           lbound,
           rbound,
           func_id,
           &nodes_x, &nodes_y);
  spline->Build(nodes_x, nodes_y,
                lbound_xx,
                rbound_xx);

  printf("Function id: %d\n", func_id);
  printf("x in [%lf; %lf]\n", lbound, rbound);
  printf("Conditions:\nS''(%lf) = %lf\nS''(%lf) = %lf\n",
         lbound, lbound_xx, rbound, rbound_xx);
  printf("Number of nodes: %d\n", number_nodes);
  printf("Net step: %lf\n", nodes_x[1] - nodes_x[0]);
  fflush(stdout);
}

void Visualize(CubicSpline& spline,
               int number_nodes,
               int func_id,
               int borders_id) {
  const double kDrawStep = 0.001;

  double lbound = GetBorder(borders_id, LEFT);
  double rbound = GetBorder(borders_id, RIGHT);
  std::vector<double> xs;
  std::vector<double> spline_y;
  std::vector<double> spline_derivates_y;
  std::vector<double> spline_second_derivates_y;
  std::vector<double> func_y;
  std::vector<double> func_derivates_y;
  std::vector<double> func_second_derivates_y;
  std::vector<double> nodes_x;
  std::vector<double> nodes_y;
  GenNodes(number_nodes,
           lbound,
           rbound,
           func_id,
           &nodes_x, &nodes_y);
  for (double x = lbound; x <= rbound; x += kDrawStep) {
    xs.push_back(x);
    spline_y.push_back(spline.GetValue(x));
    spline_derivates_y.push_back(spline.GetDerivate(x));
    spline_second_derivates_y.push_back(spline.GetSecondDerivate(x));
    func_y.push_back(Function(func_id, x));
    func_derivates_y.push_back(FunctionDerivate(func_id, x));
    func_second_derivates_y.push_back(FunctionSecondDerivate(func_id, x));
  }
  Plot plot;

  // Function and spline
  plot.Add(xs, spline_y, 0, 0.0, 0.0, 0.0, true);
  plot.Add(xs, func_y, 0, 0.0, 0.0, 1.0, true);
  plot.Add(nodes_x, nodes_y, 2, 0.0, 0.0, 0.5, false);
  plot.Show("Cubic spline. Originals");

  // Function and spline derivates
  plot.Clear();
  plot.Add(xs, spline_derivates_y, 0, 0.0, 0.0, 0.0, true);
  plot.Add(xs, func_derivates_y, 0, 0.0, 0.0, 1.0, true);
  plot.Show("Cubic spline. Derivates");

  // Function and spline second derivates
  plot.Clear();
  plot.Add(xs, spline_second_derivates_y, 0, 0.0, 0.0, 0.0, true);
  plot.Add(xs, func_second_derivates_y, 0, 0.0, 0.0, 1.0, true);
  plot.Show("Cubic spline. Second derivates");
}

