#include <stdio.h>

#include <cmath>
#include <vector>
#include <iostream>

#include "include/plot.h"
#include "include/command_line_parser.h"

enum Limit { UPPER, LOWER };

enum Method { MIDDLE_RECTS, TRAPEZES, SIMPSON };

double GetFunction(unsigned id, double x);

double GetLimit(unsigned id, Limit limit);

double Integrate(unsigned id, Method method, int n_intervals);

double GetRobustIntegral(unsigned id);

void ShowResults(unsigned id);

void PrintAbout();

int main(int argc, char** argv) {
  // Parse command line.
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
  if (!parser.Exists("id")) {
    std::cout << "Set id of integral [-id], [-h] for help" << std::endl;
    return 0;
  }

  const int integral_id = parser.Get<int>("id");
  const int n_intervals = parser.Get<int>("n");

  std::vector<double> numerical_integrals(3);
  for (int i = 0; i < 3; ++i) {
    numerical_integrals[i] = Integrate(integral_id, (Method)i, n_intervals);
  }

  double robust_value = GetRobustIntegral(integral_id);
  printf("Robust value of integral: %e\n", robust_value);
  printf("Method of middle rectangles gives: %e", numerical_integrals[0]);
  printf("(diff %e)\n", robust_value - numerical_integrals[0]);
  printf("Method of trapezes gives:          %e", numerical_integrals[1]);
  printf("(diff %e)\n", robust_value - numerical_integrals[1]);
  printf("Simpson's method gives:            %e", numerical_integrals[2]);
  printf("(diff %e)\n", robust_value - numerical_integrals[2]);

  ShowResults(integral_id);

  return 0;
}

double Integrate(unsigned id, Method method, int n_intervals) {
  const double lower_limit = GetLimit(id, LOWER);
  const double step = (GetLimit(id, UPPER) - lower_limit) / n_intervals;

  // Compute integral.
  double res = 0;
  switch (method) {

    case MIDDLE_RECTS: {
      for (int i = 0; i < n_intervals; ++i) {
        res += step * GetFunction(id, lower_limit + step * (i + 0.5));
      }
      break;
    }

    case TRAPEZES: {
      double left_point = GetFunction(id, lower_limit);
      for (int i = 1; i <= n_intervals; ++i) {
        double right_point = GetFunction(id, lower_limit + step * i);
        res += step * (right_point + left_point) / 2;
        left_point = right_point;
      }
      break;
    }

    case SIMPSON: {
      double values[3];
      values[0] = GetFunction(id, lower_limit);
      values[1] = GetFunction(id, lower_limit + step);
      for (int i = 1; i < n_intervals; ++i) {
        values[2] = GetFunction(id, lower_limit + step * (i + 1));
        res += step * (values[0] + 4 * values[1] + values[2]) / 3;
        values[0] = values[1];
        values[1] = values[2];
      }
      break;
    }

    default: return 0;
  }
  return res;
}

double GetFunction(unsigned id, double x) {
  // 0: sin2x * cosx
  // 1: sin2x * cosx + cos10x
  // 2: sin2x * cosx + cos100x
  if (id < 3) {
    double res = sin(x);
    res = 2 * res * (1 - res * res);
    if (id != 0) {
      res += cos((90 * id - 80) * x);
    }
    return res;
  }
  return 0;
}

double GetLimit(unsigned id, Limit limit) {
  switch (id) {
    case 0: case 1: case 2:
      switch (limit) {
        case UPPER: return M_PI; break;
        case LOWER: return 1; break;
        default: return 0;
      }
      break;
    default: return 0;
  }
}

double GetRobustIntegral(unsigned id) {
  double term = 2.0 / 3 * (1 + pow(cos(1), 3));
  switch (id) {
    case 0: return term;
    case 1: return term - 0.1 * sin(10);
    case 2: return term - 0.01 * sin(100);
    default: return 0;
  }
}

void ShowResults(unsigned id) {
  static const double kDrawStep = 1e-3;

  std::vector<double> func_values;
  std::vector<double> xs;
  const double upper_limit = GetLimit(id, UPPER);
  for (double x = GetLimit(id, LOWER); x <= upper_limit; x += kDrawStep) {
    func_values.push_back(GetFunction(id, x));
    xs.push_back(x);
  }

  Plot plot;
  plot.Add(xs, func_values, 1, 0, 0, 0.9, true);
  plot.Show("Integrated function");
}

void PrintAbout() {
  std::cout << "Computes numerical integral of specific functions [-id]:\n"
               "0: sin2x * cosx,           x in [1, pi]\n"
               "1: sin2x * cosx + cos10x,  x in [1, pi]\n"
               "2: sin2x * cosx + cos100x, x in [1, pi]\n\n"

               "Method parameters:\n"
               "[-n] - number of intervals for computing local integrals"
            << std::endl;
}