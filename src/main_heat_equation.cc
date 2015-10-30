#include <math.h>
#include <string>
#include <vector>
#include <stdio.h>
#include <cmath>
#include "include/command_line_parser.h"
#include "include/system_solver.h"

std::string about = "Heat equation\n"
                    "(k(x)u'(x))' - q(x)u(x) = -f(x)\n"
                    "x in (0, 1); gap point x = 0.5\n"
                    "k(x) - heat conductivity\n"
                    "q(x) - heat transfer\n"
                    "f(x) - external heat\n"
                    "\n"
                    "Test task №1 [-t 0]:\n"
                    "-1/3 u''(x) = -2, x in (0, 1)\n"
                    "u(0) = 4, u(1) = 9\n"
                    "details:\n"
                    "k(x) = -1/3, x in (0, 1)\n"
                    "q(x) = 0,    x in (0, 1)\n"
                    "f(x) = 2,    x in (0, 1)\n"
                    "\n"
                    "Test task №2 [-t 1]:\n"
                    "2.25 u''(x) - u(x) = 0, x in (0, 0.5]\n"
                    "1.25 u''(x) - u(x) = -1, x in (0.5, 1)\n"
                    "u(0) = 0, u(1) = 1\n"
                    "details:\n"
                    "k(x) = 2.25, x in (0, 0.5]\n"
                    "     = 1.25, x in (0.5, 1)\n"
                    "q(x) = 1, x in (0, 1)\n"
                    "f(x) = 0, x in (0, 0.5]\n"
                    "     = 1, x in (0.5, 1)";

enum Border { LEFT, RIGHT };
enum IntegralType { BY_EXTERNAL_HEAT, BY_HEAT_TRANSFER, BY_HEAT_CONDUCTIVITY };

static const double kGapPoint = 0.5;
static const double kZeroLimit = 1e-10;

double GetHeatConductivity(int task_id, double x);

double GetHeatTransfer(int task_id, double x);

double GetExternalHeat(int task_id, double x);

double GetBorder(Border border);

double GetBorderHeat(int task_id, Border border);

void GenNodes(int number_intervals,
              std::vector<double>* nodes);

double GetIntegral(int task_id,
                   IntegralType type, int node_idx,
                   const std::vector<double>& nodes);

void Solve(int task_id, int number_intervals);

void PrintResults(int task_id,
                  const std::vector<double>& numerical_func,
                  const std::vector<double>& nodes);

double GetRobustValue(int task_id, double x);

int main(int argc, char** argv) {
  CommandLineParser parser(argc, argv);
  if (parser.Exists("h") || argc == 1) {
    std::cout << about << std::endl;
    return 0;
  }
  if (!parser.Exists("t")) {
    std::cout << "Select task [-t], [-h] for help" << std::endl;
    return 0;
  }
  if (!parser.Exists("n")) {
    std::cout << "Set number of intervals [-n], [-h] for help" << std::endl;
    return 0;
  }
  const int task_id = parser.Get<int>("f");
  const int number_intervals = parser.Get<int>("n");

  Solve(task_id, number_intervals);
}

void Solve(int task_id, int number_intervals) {
  int n = number_intervals;
  std::vector<double> nodes;
  GenNodes(n, &nodes);

  std::vector<double> right_part(n - 1);
  for (int i = 0; i < right_part.size(); ++i) {
    right_part[i] = -GetIntegral(task_id, BY_EXTERNAL_HEAT, i + 1, nodes);
  }
  right_part[0] -= GetIntegral(task_id, BY_HEAT_CONDUCTIVITY, 1, nodes) *
                   GetBorderHeat(task_id, LEFT);
  right_part[n - 2] -= GetIntegral(task_id, BY_HEAT_CONDUCTIVITY, n, nodes) *
                       GetBorderHeat(task_id, RIGHT);

  std::vector<double> sub_diag(n - 2);
  std::vector<double> main_diag(n - 1);
  for (int i = 0; i < main_diag.size(); ++i) {
    main_diag[i] = -GetIntegral(task_id, BY_HEAT_CONDUCTIVITY, i + 1, nodes)
                   -GetIntegral(task_id, BY_HEAT_CONDUCTIVITY, i + 2, nodes)
                   -GetIntegral(task_id, BY_HEAT_TRANSFER, i + 1, nodes);
  }
  for (int i = 0; i < sub_diag.size(); ++i) {
    sub_diag[i] = GetIntegral(task_id, BY_HEAT_CONDUCTIVITY, i + 2, nodes);
  }

  std::vector<double> numerical_func;
  SystemSolver system_solver;
  system_solver.Solve(sub_diag, main_diag, sub_diag,
                      right_part, &numerical_func);

  numerical_func.insert(numerical_func.begin(), GetBorderHeat(task_id, LEFT));
  numerical_func.push_back(GetBorderHeat(task_id, RIGHT));

  for (int i = 0; i < numerical_func.size(); ++i) {
    if (std::abs(numerical_func[i]) < kZeroLimit) {
      numerical_func[i] = 0;
    }
  }

  PrintResults(task_id, numerical_func, nodes);
}

void GenNodes(int number_intervals,
              std::vector<double>* nodes) {
  double lbound = GetBorder(LEFT);
  double rbound = GetBorder(RIGHT);
  nodes->resize(number_intervals + 1);
  double step = (rbound - lbound) / number_intervals;
  for (int i = 0; i < nodes->size(); ++i) {
    nodes->at(i)= lbound + i * step;
  }
}

double GetBorderHeat(int task_id, Border border) {
  switch (task_id) {
    case 0: {
      if (border == LEFT) return 4;
      if (border == RIGHT) return 9;
      return 0;
    }
    case 1: {
      if (border == LEFT) return 0;
      if (border == RIGHT) return 1;
      return 0;
    }
    default: return 0;
  }
}

double GetIntegral(int task_id,
                   IntegralType type, int node_idx,
                   const std::vector<double>& nodes) {
  double h = nodes[1] - nodes[0];
  switch (type) {
    case BY_EXTERNAL_HEAT: {
      if (nodes[node_idx] + h / 2 <= kGapPoint ||
          nodes[node_idx] - h / 2 >= kGapPoint) {
        return h * GetExternalHeat(task_id, nodes[node_idx]);
      } else {
        return (kGapPoint - nodes[node_idx] + h / 2) *
               GetExternalHeat(task_id,
                               (kGapPoint + nodes[node_idx] - h / 2) / 2) +
               (nodes[node_idx] + h / 2 - kGapPoint) *
               GetExternalHeat(task_id,
                               (nodes[node_idx] + h / 2 + kGapPoint) / 2);
      }
    }
    case BY_HEAT_TRANSFER: {
      if (nodes[node_idx] + h / 2 <= kGapPoint ||
          nodes[node_idx] - h / 2 >= kGapPoint) {
        return h * GetHeatTransfer(task_id, nodes[node_idx]);
      } else {
        return (kGapPoint - nodes[node_idx] + h / 2) *
               GetHeatTransfer(task_id,
                               (kGapPoint + nodes[node_idx] - h / 2) / 2) +
               (nodes[node_idx] + h / 2 - kGapPoint) *
               GetHeatTransfer(task_id,
                               (nodes[node_idx] + h / 2 + kGapPoint) / 2);
      }
    }
    case BY_HEAT_CONDUCTIVITY: {
      if (nodes[node_idx] <= kGapPoint ||
          nodes[node_idx - 1] >= kGapPoint) {
        return GetHeatConductivity(task_id, nodes[node_idx] - h / 2) / h;
      } else {
        return GetHeatConductivity(task_id,
                                   (kGapPoint + nodes[node_idx - 1]) / 2) /
               (kGapPoint - nodes[node_idx - 1]) +
               GetHeatConductivity(task_id,
                                   (kGapPoint + nodes[node_idx]) / 2) /
               (nodes[node_idx] - kGapPoint);
      }
    }
    default: return 0;
  }
}

double GetBorder(Border border) {
  switch (border) {
    case LEFT: return 0;
    case RIGHT: return 1;
  }
}

double GetHeatConductivity(int task_id, double x) {
  switch (task_id) {
    case 0: return -1.0 / 3;
    case 1: {
      if (x <= kGapPoint) {
        return 2.25;
      } else {
        return 1.25;
      }
    }
    default: return 0;
  }

/*  if (x <= kGapPoint) {
    return pow(x + 1, 2);
  } else {
    return pow(x, 2) + 1;
  }*/
}

double GetHeatTransfer(int task_id, double x) {
  switch (task_id) {
    case 0: return 0;
    case 1: return 1;
    default: return 0;
  }

/*  if (x <= kGapPoint) {
    return exp(-x + 0.5);
  } else {
    return exp(x - 0.5);
  }*/
}

double GetExternalHeat(int task_id, double x) {
  switch(task_id) {
    case 0: return 2;
    case 1: {
      if (x <= kGapPoint) {
        return 0;
      } else {
        return 1;
      }
    }
  }
 /* if (x <= kGapPoint) {
    return cos(M_PI * x);
  } else {
    return 1.0;
  }*/
}

void PrintResults(int task_id,
                  const std::vector<double>& numerical_func,
                  const std::vector<double>& nodes) {
  const int kMaxOutLines = 20;
  const int kNumberCells = 5;
  const std::string kFmtStr("|% 13s");
  std::string cells[] = {"i", "x[i]", "u(x[i])", "v(x[i])",
                         "u - v"};

  // Header.
  for (int i = 0; i < kNumberCells; ++i) {
    printf(kFmtStr.c_str(), cells[i].c_str());
  }
  printf("|\n");

  for (int i = 0; i < kNumberCells; ++i) {
    printf("|");
    for (int j = 0; j < 13; ++j) {
       printf("-");
    }
  }
  printf("|\n");

  // Lines.
  if (nodes.size() <= kMaxOutLines) {
    for (int i = 0; i < nodes.size(); ++i) {
      double robust = GetRobustValue(task_id, nodes[i]);
      printf("|% 13d|% 13e|% 13e|% 13e|% 13e|\n",
             i,
             nodes[i],
             robust,
             numerical_func[i],
             robust - numerical_func[i]);
    }
  } else {
    for (int i = 0; i < kMaxOutLines / 2; ++i) {
      double robust = GetRobustValue(task_id, nodes[i]);
      printf("|% 13d|% 13e|% 13e|% 13e|% 13e|\n",
             i,
             nodes[i],
             robust,
             numerical_func[i],
             robust - numerical_func[i]);
    }
    for (int i = 0; i < kNumberCells; ++i) {
      printf("|");
      for (int j = 0; j < 13; ++j) {
         printf(".");
      }
    }
    printf("|\n");
    for (int i = nodes.size() - 1 - kMaxOutLines / 2;
         i < nodes.size() - 1; ++i) {
      double robust = GetRobustValue(task_id, nodes[i]);
      printf("|% 13d|% 13e|% 13e|% 13e|% 13e|\n",
             i,
             nodes[i],
             robust,
             numerical_func[i],
             robust - numerical_func[i]);
    }
  }
  fflush(stdout);
}

double GetRobustValue(int task_id, double x) {
  switch (task_id) {
    case 0:
      return 3 * x * x + 2 * x + 4;
    case 1:
      if (x <= kGapPoint)
        return 0.53529341 * exp(2.0 / 3 * x) - 0.53529341 * exp(-2.0 / 3 * x);
      else
        return 0.28146359 * exp(2.0 / sqrt(5) * x) -
               1.68388261 * exp(-2.0 / sqrt(5) * x);
    default: return 0;
  }
}
