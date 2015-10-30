// Copyright 2015 Dmitry Kurtaev

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "include/plot.h"
#include "include/cubic_spline.h"
#include "include/command_line_parser.h"

Plot plot;

enum Border { LEFT, RIGHT };

double Function(int func_id, double x);

double FunctionDerivate(int func_id, double x);

double FunctionSecondDerivate(int func_id, double x);

double GetBorderCondition(int condition_id,
                          int func_id,
                          double border);

double GetBorder(int borders_id, Border border);

void GenNodes(int number_intervals,
              double lbound,
              double rbound,
              int func_id,
              std::vector<double>* x,
              std::vector<double>* y);

void CheckAccuracy(const CubicSpline& spline,
                   int func_id,
                   int borders_id);

void BuildSpline(int func_id,
                 int borders_id,
                 int borders_condition_id,
                 int number_intervals,
                 CubicSpline* spline);

void Visualize(const CubicSpline& spline,
               int number_intervals,
               int func_id,
               int borders_id);

int main(int argc, char** argv) {
  CommandLineParser parser(argc, argv);
  if (parser.Exists("h") || argc == 1) {
    std::cout << "Functions [-f]:\n"
                 "0. f(x) = x^3 + 3x^2, if x in [-1, 0]\n"
                 "          -x^3 + 3x^2, if x in [0, 1];\n"
                 "1. f(x) = sqrt(exp(x) - 1)\n"
                 "2. f(x) = sqrt(exp(x) - 1) + cos10x\n"
                 "3. f(x) = sqrt(exp(x) - 1) + cos100x\n"
                 "\n"
                 "Borders [-b]:\n"
                 "0. [-1, 1]\n"
                 "1. [1, 3]\n"
                 "\n"
                 "Borders conditions [-bc]:\n"
                 "0. [0, 0] identity\n"
                 "1. [f''(lbound), f''(rbound)] robust\n" << std::endl;
    return 0;
  }
  if (!parser.Exists("f")) {
    std::cout << "Set function [-f]" << std::endl;
    return 0;
  }
  if (!parser.Exists("b")) {
    std::cout << "Set borders [-b], [-h] for help" << std::endl;
    return 0;
  }
  if (!parser.Exists("bc")) {
    std::cout << "Set borders condition [-bc], [-h] for help" << std::endl;
    return 0;
  }
  if (!parser.Exists("n")) {
    std::cout << "Set number of intervals [-n], [-h] for help" << std::endl;
    return 0;
  }
  const int func_id = parser.Get<int>("f");
  const int borders_id = parser.Get<int>("b");
  const int borders_condition_id = parser.Get<int>("bc");
  const int number_intervals = parser.Get<int>("n");

  CubicSpline spline;
  BuildSpline(func_id,
              borders_id,
              borders_condition_id,
              number_intervals,
              &spline);
  CheckAccuracy(spline,
                func_id,
                borders_id);
  spline.PrintCoeffs();
  Visualize(spline,
            number_intervals,
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

void GenNodes(int number_intervals,
              double lbound,
              double rbound,
              int func_id,
              std::vector<double>* x,
              std::vector<double>* y) {
  double step = (rbound - lbound) / number_intervals;
  x->resize(number_intervals + 1);
  y->resize(number_intervals + 1);
  double current_x = lbound;
  for (int i = 0; i < number_intervals + 1; ++i) {
    (*x)[i] = current_x;
    (*y)[i] = Function(func_id, current_x);
    current_x += step;
  }
}

void CheckAccuracy(const CubicSpline& spline,
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
                 int number_intervals,
                 CubicSpline* spline) {
  double lbound = GetBorder(borders_id, LEFT);
  double rbound = GetBorder(borders_id, RIGHT);

  double lbound_xx = GetBorderCondition(borders_condition_id,
                                        func_id, lbound);
  double rbound_xx = GetBorderCondition(borders_condition_id,
                                        func_id, rbound);
  std::vector<double> nodes_x;
  std::vector<double> nodes_y;
  GenNodes(number_intervals,
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
  printf("Number of intervals: %d\n", number_intervals);
  printf("Net step: %lf\n", nodes_x[1] - nodes_x[0]);
  fflush(stdout);
}

void Visualize(const CubicSpline& spline,
               int number_intervals,
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
  GenNodes(number_intervals,
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

