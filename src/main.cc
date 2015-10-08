#include "include/plot.h"
#include "include/cubic_spline.h"
#include <math.h>

Plot plot;

enum Border { LEFT, RIGHT };

// Test functions:
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
// 0. [0, 12]

double TestFunction(int func_id, double x);

double TestFunctionDerivate(int func_id, double x);

double TestFunctionSecondDerivate(int func_id, double x);

double GetBorderCondition(int condition_id, Border border);

double GetBorder(int borders_id, Border border);

void GenNodes(int number,
              double lbound,
              double rbound,
              int func_id,
              std::vector<double>* x,
              std::vector<double>* y);

int main(int argc, char** argv) {
  int func_id = 0;
  int borders_id = 0;
  int borders_condition_id = 0;
  int number_nodes = 10;
  std::vector<double> x;
  std::vector<double> y;
  GenNodes(number_nodes,
           GetBorder(borders_id, LEFT),
           GetBorder(borders_id, RIGHT),
           func_id,
           &x, &y);
  CubicSpline spline;
  spline.Build(x, y,
               GetBorderCondition(borders_condition_id, LEFT),
               GetBorderCondition(borders_condition_id, RIGHT));
  spline.Show();

  return 0;
}

double TestFunction(int func_id, double x) {
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
    case 2: return TestFunction(1, x) + cos(10 * x);
    case 3: return TestFunction(1, x) + cos(100 * x);
    default: return 0;
  }
}

double TestFunctionDerivate(int func_id, double x) {
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
    case 2: return TestFunctionDerivate(1, x) - 10 * sin(10 * x);
    case 3: return TestFunctionDerivate(1, x) - 100 * sin(100 * x);
    default: return 0;
  }
}

double TestFunctionSecondDerivate(int func_id, double x) {
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
    case 2: return TestFunctionSecondDerivate(1, x) - 100 * cos(10 * x);
    case 3: return TestFunctionSecondDerivate(1, x) - 1e+5 * cos(100 * x);
    default: return 0;
  }
}

double GetBorderCondition(int condition_id, Border border) {
  switch (condition_id) {
    case 0: {
      if (border == LEFT) return 0.0;
      if (border == RIGHT) return 12.0;
      return 0;
    }
    case 1: {
      return 0.0;
    }
    case 2: {
      return 0.0;
    }
    default: return 0.0;
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
  double step = (rbound - lbound) / (number + 1);
  x->resize(number + 2);
  y->resize(number + 2);
  double current_x = lbound;
  for (int i = 0; i < number + 2; ++i) {
    (*x)[i] = current_x;
    (*y)[i] = TestFunction(func_id, current_x);
    current_x += step;
  }
}

