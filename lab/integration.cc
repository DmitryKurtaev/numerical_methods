#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#include <cmath>
#include <vector>
#include <iostream>

#include "include/plot.h"
#include "include/command_line_parser.h"

enum Limit { UPPER, LOWER };

enum Method { MIDDLE_RECTS, TRAPEZES, SIMPSON };

double GetFunction(unsigned id, double x);

double GetWaveFunction(double x, double t);

double GetLimit(unsigned id, Limit limit);

double Integrate(unsigned id, Method method, int n_intervals);

// Computes using Simpson's formula.
double IntegrateWaveFunction(double x, double lower_limit, double upper_limit,
                             int n_intervals);

double IntegrateWaveFunction(double x, int n_intervals, double eps);

double GetRobustIntegral(unsigned id);

void ShowFunction(unsigned id);

void ShowWaveFunction(int n_intervals, double eps);

void PrintAbout();

void PrintWaiter();

static const unsigned kWaveFunctionId = 3;

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

  if (integral_id != kWaveFunctionId) {
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

    ShowFunction(integral_id);
  } else {
    ShowWaveFunction(n_intervals, 1e-3);
  }

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
      double left_point = GetFunction(id, lower_limit);
      for (int i = 0; i < n_intervals; i += 2) {
        double middle_point = GetFunction(id, lower_limit + step * (i + 1));
        double right_point = GetFunction(id, lower_limit + step * (i + 2));
        res += step * (left_point + 4 * middle_point + right_point) / 3;
        left_point = right_point;
      }
      break;
    }

    default: return 0;
  }
  return res;
}

double IntegrateWaveFunction(double x, double lower_limit, double upper_limit,
                             int n_intervals) {
  const double step = (upper_limit - lower_limit) / n_intervals;
  double left_point = GetWaveFunction(x, lower_limit);
  double res = 0;
  for (int i = 0; i < n_intervals; i += 2) {
    double middle_point = GetWaveFunction(x, lower_limit + step * (i + 1));
    double right_point = GetWaveFunction(x, lower_limit + step * (i + 2));
    res += step * (left_point + 4 * middle_point + right_point) / 3;
    left_point = right_point;
  }
  return res;
}

double IntegrateWaveFunction(double x, int n_intervals, double eps) {
  double lower_limit = GetLimit(kWaveFunctionId, LOWER);
  double upper_limit = GetLimit(kWaveFunctionId, UPPER);
  std::vector<double> lower_limits(1, lower_limit);
  std::vector<double> upper_limits(1, upper_limit);
  std::vector<double> target_accuracies(1, eps);
  std::vector<double> global_integrals(1);

  global_integrals[0] = IntegrateWaveFunction(x, lower_limit, upper_limit,
                                              n_intervals);
  double res = 0;
  do {
    lower_limit = lower_limits.back();
    upper_limit = upper_limits.back();
    eps = target_accuracies.back();
    const double global_integral = global_integrals.back();

    lower_limits.pop_back();
    upper_limits.pop_back();
    target_accuracies.pop_back();
    global_integrals.pop_back();

    const double middle = (lower_limit + upper_limit) / 2;
    const double left_integral = IntegrateWaveFunction(x, lower_limit, middle,
                                                       n_intervals);
    const double right_integral = IntegrateWaveFunction(x, middle, upper_limit,
                                                        n_intervals);

    if (fabs(left_integral + right_integral - global_integral) >= eps) {
      eps /= 2;
      global_integrals.push_back(left_integral);
      lower_limits.push_back(lower_limit);
      upper_limits.push_back(middle);
      target_accuracies.push_back(eps);

      global_integrals.push_back(right_integral);
      lower_limits.push_back(middle);
      upper_limits.push_back(upper_limit);
      target_accuracies.push_back(eps);
    } else {
      res += left_integral + right_integral;
    }
  } while (!global_integrals.empty());
  return res;
}

double GetFunction(unsigned id, double x) {
  // 0: sin2x * cosx
  // 1: sin2x * cosx + cos10x
  // 2: sin2x * cosx + cos100x
  if (id != kWaveFunctionId) {
    double res = sin(x);
    res = 2 * res * (1 - res * res);
    if (id != 0) {
      res += cos((90 * id - 80) * x);
    }
    return res;
  }
  return 0;
}

double GetWaveFunction(double x, double t) {
  static std::vector<double> params_a;
  static std::vector<double> params_b;
  static double alpha;

  if (params_a.empty()) {
    srand(time(0));
    params_a.resize(14);
    params_b.resize(14);
    for (int i = 0; i < 14; ++i) {
      params_a[i] = (static_cast<double>(rand()) / RAND_MAX) * 2 - 1;
      params_b[i] = (static_cast<double>(rand()) / RAND_MAX) * 2 - 1;
    }
    alpha = static_cast<double>(rand()) / RAND_MAX;
    std::cout << "Generated parameter alpha = " << alpha << std::endl;
  }

  double res = 0;
  for (int i = 0; i < 14; ++i) {
    const double arg = 2 * M_PI * i * (alpha - x) * t;
    res += params_a[i] * sin(arg) + params_b[i] * cos(arg);
  }
  return res;
}

double GetLimit(unsigned id, Limit limit) {
  switch (id) {
    case 0: case 1: case 2:
      switch (limit) {
        case UPPER: return M_PI;
        case LOWER: return 1;
        default: return 0;
      }
      break;
    case kWaveFunctionId:
      switch (limit) {
        case UPPER: return M_PI / 2;
        case LOWER: return -M_PI / 2;
        default: return 0;
      }
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

void ShowWaveFunction(int n_intervals, double eps) {
  static const unsigned kDrawPoints = 1025;
  static const double kLeftValue = 0.0;
  static const double kRightValue = 1.0;

  const double draw_step = (kRightValue - kLeftValue) / (kDrawPoints - 1);

  std::vector<double> func_values(kDrawPoints);
  std::vector<double> xs(kDrawPoints);

  for (int i = 0; i < kDrawPoints; ++i) {
    const double x = kLeftValue + i * draw_step;
    func_values[i] =  IntegrateWaveFunction(x, n_intervals, eps);
    xs[i] = x;

    std::cout << "\r";
    PrintWaiter();
    std::cout << " Computing " << i << '/' << (kDrawPoints - 1) << std::flush;
  }
  std::cout << std::endl;

  Plot plot;
  plot.Add(xs, func_values, 1, 0, 0, 0.9, true);
  plot.Show("Integrated function", "x", "g(x)");
}

void PrintWaiter() {
  static const unsigned kAnimPeriod = 1000;

  static int symbol_id = 7;
  static char symbols[8] = {'|', '/', '-', '\\', '|', '/', '-', '\\'};
  static timeval last_animation;

  timeval new_time;
  gettimeofday(&new_time, 0);
  if ((new_time.tv_sec - last_animation.tv_sec) * 1e+3 +
      (new_time.tv_usec - last_animation.tv_usec) / 1e+3 >= kAnimPeriod) {
    symbol_id = (symbol_id != 7 ? symbol_id + 1 : 0);
    last_animation = new_time;
  }
  std::cout << symbols[symbol_id] << std::flush;
}

void ShowFunction(unsigned id) {
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
