#include <math.h>
#include <string>
#include <vector>
#include <stdio.h>
#include <cmath>
#include "include/command_line_parser.h"
#include "include/system_solver.h"
#include "include/plot.h"

enum Border { LEFT, RIGHT };
enum Coefficient { EXTERNAL_HEAT, HEAT_TRANSFER, HEAT_CONDUCTIVITY };

static const double kGapPoint = 0.5;
static const double kLeftBorder = 0;
static const double kRightBorder = 1;
static const double kZeroLimit = 1e-10;

double GetCoefficient(Coefficient coeff_type, int task_id, double x);

double GetBorderHeat(int task_id, Border border);

void GenNodes(int number_intervals,
              std::vector<double>* nodes);

double GetIntegral(int task_id,
                   Coefficient coeff_type, int node_idx,
                   const std::vector<double>& nodes);

void Solve(int task_id, int number_intervals,
           bool visualize, bool keep_quiet,
           std::vector<double>* numerical_func,
           std::vector<double>* nodes);

void PrintResults(int task_id,
                  const std::vector<double>& numerical_func,
                  const std::vector<double>& robust_func,
                  const std::vector<double>& nodes);

double GetRobustValue(int task_id, double x);

void Visualize(int task_id,
               const std::vector<double>& numerical_func,
               const std::vector<double>& robust_func,
               const std::vector<double>& nodes);

double CheckAccuracy(int task_id,
                     const std::vector<double>& numerical_func,
                     const std::vector<double>& robust_func,
                     const std::vector<double>& nodes,
                     bool keep_quiet);

void Experiment(int experiment_id, int step,
                int task_id, int max_number_intervals, double mult);

void PrintAbout();

int main(int argc, char** argv) {
  CommandLineParser parser(argc, argv);
  if (parser.Exists("h") || argc == 1) {
    PrintAbout();
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

  const int task_id = parser.Get<int>("t");
  const int number_intervals = parser.Get<int>("n");
  const bool keep_quiet = parser.Exists("q");
  if (!parser.Exists("e")) {
    const bool visualize = parser.Exists("vis");

    printf("Task id: %d\n", task_id);
    printf("x in [%lf, %lf]\n", kLeftBorder, kRightBorder);
    printf("Number of nodes: %d\n", number_intervals + 1);
    printf("Net step: %lf\n", (kRightBorder - kLeftBorder) / number_intervals);

    std::vector<double> numerical_func;
    std::vector<double> nodes;
    std::vector<double> robust_func;
    Solve(task_id, number_intervals,
          visualize, keep_quiet,
          &numerical_func, &nodes);
    if (task_id != 2) {
      for (int i = 0; i < nodes.size(); ++i) {
        robust_func.push_back(GetRobustValue(task_id, nodes[i]));
      }
    } else {
      std::vector<double> tmp_nodes;
      std::vector<double> tmp_robust;
      Solve(task_id, 2 * number_intervals,
            false, true,
            &tmp_robust, &tmp_nodes);
      for (int i = 0; i < nodes.size(); ++i) {
        robust_func.push_back(tmp_robust[i * 2]);
      }
    }
    CheckAccuracy(task_id, numerical_func, robust_func, nodes, keep_quiet);
    if (!keep_quiet) {
      PrintResults(task_id, numerical_func, robust_func, nodes);
      if (visualize) {
        Visualize(task_id, numerical_func, robust_func, nodes);
      }
    }
  } else {
    const int step = parser.Get<int>("s");
    const double mult = parser.Get<double>("m");
    const int experiment_id = parser.Get<int>("e");
    Experiment(experiment_id, step,
               task_id, number_intervals, mult);
  }
}

void Experiment(int experiment_id, int step,
                int task_id, int number_intervals, double mult) {
  switch (experiment_id) {
    case 0: {
      std::vector<double> numbers;
      std::vector<double> accuracies;
      for (int n = 1000; n <= number_intervals; n += step) {
        numbers.push_back(n);
        std::vector<double> numerical_func;
        std::vector<double> nodes;
        Solve(task_id, n,
              false, true,
              &numerical_func, &nodes);
        std::vector<double> robust_func;
        if (task_id != 2) {
          for (int i = 0; i < nodes.size(); ++i) {
            robust_func.push_back(GetRobustValue(task_id, nodes[i]));
          }
        } else {
          std::vector<double> tmp_nodes;
          std::vector<double> tmp_robust;
          Solve(task_id, 2 * n,
                false, true,
                &tmp_robust, &tmp_nodes);
          for (int i = 0; i < nodes.size(); ++i) {
            robust_func.push_back(tmp_robust[i * 2]);
          }
        }
        double diff = CheckAccuracy(task_id, numerical_func,
                                    robust_func, nodes, true);
        accuracies.push_back(diff);
      }
      Plot plot;
      plot.Add(numbers, accuracies, 1, 0.5, 0, 0.9, true);
      plot.Show("Dependences between number of intervals and accuracy");
      break;
    }
    case 1: {
      double n = number_intervals;
      double last_accuracy = 1;
      while(true) {
        if (static_cast<int>(n) % 2 != 0) {
          ++n;
        }
        std::vector<double> numerical_func;
        std::vector<double> nodes;
        Solve(task_id, n,
              false, true,
              &numerical_func, &nodes);
        std::vector<double> robust_func;
        if (task_id != 2) {
          for (int i = 0; i < nodes.size(); ++i) {
            robust_func.push_back(GetRobustValue(task_id, nodes[i]));
          }
        } else {
          std::vector<double> tmp_nodes;
          std::vector<double> tmp_robust;
          Solve(task_id, 2 * n,
                false, true,
                &tmp_robust, &tmp_nodes);
          for (int i = 0; i < nodes.size(); ++i) {
            robust_func.push_back(tmp_robust[i * 2]);
          }
        }
        double accuracy = CheckAccuracy(task_id, numerical_func,
                                        robust_func, nodes, true);
        printf("% 10d % 10e (% 10lf)\n",
               static_cast<int>(n), accuracy, last_accuracy / accuracy);
        fflush(stdout);
        last_accuracy = accuracy;
        n *= mult;
      }
      break;
    }
    default: break;
  }
}

void Solve(int task_id, int number_intervals,
           bool visualize, bool keep_quiet,
           std::vector<double>* numerical_func,
           std::vector<double>* nodes) {
  const int n = number_intervals;
  nodes->clear();
  GenNodes(n, nodes);

  std::vector<double> right_part(n - 1);
  for (int i = 0; i < right_part.size(); ++i) {
    right_part[i] = -GetIntegral(task_id, EXTERNAL_HEAT, i + 1, *nodes);
  }
  right_part[0] -= GetIntegral(task_id, HEAT_CONDUCTIVITY, 1, *nodes) *
                   GetBorderHeat(task_id, LEFT);
  right_part[n - 2] -= GetIntegral(task_id, HEAT_CONDUCTIVITY, n, *nodes) *
                       GetBorderHeat(task_id, RIGHT);

  std::vector<double> sub_diag(n - 2);
  std::vector<double> main_diag(n - 1);
  for (int i = 0; i < main_diag.size(); ++i) {
    main_diag[i] = -GetIntegral(task_id, HEAT_CONDUCTIVITY, i + 1, *nodes)
                   -GetIntegral(task_id, HEAT_CONDUCTIVITY, i + 2, *nodes)
                   -GetIntegral(task_id, HEAT_TRANSFER, i + 1, *nodes);
  }
  for (int i = 0; i < sub_diag.size(); ++i) {
    sub_diag[i] = GetIntegral(task_id, HEAT_CONDUCTIVITY, i + 2, *nodes);
  }

  numerical_func->clear();
  SystemSolver system_solver;
  system_solver.Solve(sub_diag, main_diag, sub_diag,
                      right_part, numerical_func);

  numerical_func->insert(numerical_func->begin(), GetBorderHeat(task_id, LEFT));
  numerical_func->push_back(GetBorderHeat(task_id, RIGHT));

  for (int i = 0; i < numerical_func->size(); ++i) {
    if (std::abs(numerical_func->operator [](i)) < kZeroLimit) {
      numerical_func->operator [](i) = 0;
    }
  }
}

double CheckAccuracy(int task_id,
                     const std::vector<double>& numerical_func,
                     const std::vector<double>& robust_func,
                     const std::vector<double>& nodes,
                     bool keep_quiet) {
  double max_diff = 0;
  int max_diff_idx = 0;
  for (int i = 0; i < nodes.size(); ++i) {
    double diff = std::abs(robust_func[i] - numerical_func[i]);
    if (diff > max_diff) {
      max_diff = diff;
      max_diff_idx = i;
    }
  }
  if (!keep_quiet) {
    std::string fmt_str = (task_id != 2?
                          "max|u(x[i]) - v(x[i])| = %e (at x[%d] = %e)\n" :
                          "max|v(x[i]) - v2(x[i])| = %e (at x[%d] = %e)\n");
    printf(fmt_str.c_str(),
           max_diff,
           max_diff_idx,
           nodes[max_diff_idx]);
    fflush(stdout);
  }
  return max_diff;
}

void Visualize(int task_id,
               const std::vector<double>& numerical_func,
               const std::vector<double>& robust_func,
               const std::vector<double>& nodes) {
  const double kDrawStep = 0.001;

  std::vector<double> nodes_for_draw;

  GenNodes((kRightBorder - kLeftBorder) / kDrawStep, &nodes_for_draw);

  Plot plot;

  // Robust function and numerical result.
  if (task_id != 2) {
    std::vector<double> robust_func_for_draw(nodes_for_draw.size());
    for (int i = 0; i < nodes_for_draw.size(); ++i) {
      robust_func_for_draw[i] = GetRobustValue(task_id, nodes_for_draw[i]);
    }
    plot.Add(nodes_for_draw, robust_func_for_draw, 1, 0, 0, 1, true);
  } else {
    plot.Add(nodes, robust_func, 1, 0, 0, 1, true);
  }
  plot.Add(nodes, numerical_func, 1, 0.5, 0, 0.5, true);
  plot.Add(nodes, numerical_func, 2, 0.5, 0, 0.5, false);
  plot.Show("Blue - robust function. Purple - numerical result");

  // Difference between robust function and numerical result.
  std::vector<double> difference_in_nodes(nodes.size());
  for (int i = 0; i < nodes.size(); ++i) {
    difference_in_nodes[i] = robust_func[i] - numerical_func[i];
  }
  plot.Clear();
  plot.Add(nodes, difference_in_nodes, 1, 1, 0, 0, true);
  plot.Add(nodes, difference_in_nodes, 2, 1, 0, 0, false);
  plot.Show("Difference");
}

void GenNodes(int number_intervals,
              std::vector<double>* nodes) {
  nodes->resize(number_intervals + 1);
  double step = (kRightBorder - kLeftBorder) / number_intervals;
  for (int i = 0; i < nodes->size(); ++i) {
    nodes->at(i)= kLeftBorder + i * step;
  }
}

double GetBorderHeat(int task_id, Border border) {
  switch (task_id) {
    case 0: {
      if (border == LEFT) return 4;
      if (border == RIGHT) return 9;
      return 0;
    }
    case 1: case 2: {
      if (border == LEFT) return 0;
      if (border == RIGHT) return 1;
      return 0;
    }
    default: return 0;
  }
}

double GetIntegral(int task_id,
                   Coefficient coeff_type, int node_idx,
                   const std::vector<double>& nodes) {
  double h = nodes[1] - nodes[0];
  switch (coeff_type) {
    case EXTERNAL_HEAT: {
      if (nodes[node_idx] + h / 2 <= kGapPoint ||
          nodes[node_idx] - h / 2 >= kGapPoint) {
        return h * GetCoefficient(EXTERNAL_HEAT, task_id, nodes[node_idx]);
      } else {
        return (kGapPoint - nodes[node_idx] + h / 2) *
            GetCoefficient(EXTERNAL_HEAT, task_id,
                           (kGapPoint + nodes[node_idx] - h / 2) / 2) +
            (nodes[node_idx] + h / 2 - kGapPoint) *
            GetCoefficient(EXTERNAL_HEAT, task_id,
                           (nodes[node_idx] + h / 2 + kGapPoint) / 2);
      }
    }
    case HEAT_TRANSFER: {
      if (nodes[node_idx] + h / 2 <= kGapPoint ||
          nodes[node_idx] - h / 2 >= kGapPoint) {
        return h * GetCoefficient(HEAT_TRANSFER, task_id, nodes[node_idx]);
      } else {
        return (kGapPoint - nodes[node_idx] + h / 2) *
            GetCoefficient(HEAT_TRANSFER, task_id,
                           (kGapPoint + nodes[node_idx] - h / 2) / 2) +
            (nodes[node_idx] + h / 2 - kGapPoint) *
            GetCoefficient(HEAT_TRANSFER, task_id,
                           (nodes[node_idx] + h / 2 + kGapPoint) / 2);
      }
    }
    case HEAT_CONDUCTIVITY: {
      if (nodes[node_idx] <= kGapPoint ||
          nodes[node_idx - 1] >= kGapPoint) {
        return GetCoefficient(HEAT_CONDUCTIVITY, task_id,
                              nodes[node_idx] - h / 2) / h;
      } else {
        double first_point = (kGapPoint + nodes[node_idx - 1]) / 2;
        double second_point = (kGapPoint + nodes[node_idx]) / 2;
        double first_term = (kGapPoint - nodes[node_idx - 1]) /
            GetCoefficient(HEAT_CONDUCTIVITY, task_id, first_point);
        double second_term = (nodes[node_idx] - kGapPoint) /
            GetCoefficient(HEAT_CONDUCTIVITY, task_id, second_point);
        return 1.0 / (first_term + second_term);
      }
    }
    default: return 0;
  }
}

double GetCoefficient(Coefficient coeff_type, int task_id, double x) {
  switch (coeff_type) {
    case HEAT_CONDUCTIVITY: {
      switch (task_id) {
        case 0: return -1.0 / 3;
        case 1:
          if (x <= kGapPoint) return 2.25;
          else return 1.25;
        case 2:
          if (x <= kGapPoint) return pow(x + 1, 2);
          else return x * x + 1;
        default: return 0;
      }
    }
    case HEAT_TRANSFER: {
      switch (task_id) {
        case 0: return 0;
        case 1: return 1;
        case 2:
          if (x <= kGapPoint) return exp(0.5 - x);
          else return exp(x - 0.5);
        default: return 0;
      }
    }
    case EXTERNAL_HEAT: {
      switch(task_id) {
        case 0: return 2;
        case 1:
          if (x <= kGapPoint) return 0;
          else return 1;
        case 2:
          if (x <= kGapPoint) return cos(M_PI * x);
          else return 1.0;
        default: return 0;
      }
    }
    default: return 0;
  }
}

void PrintResults(int task_id,
                  const std::vector<double>& numerical_func,
                  const std::vector<double>& robust_func,
                  const std::vector<double>& nodes) {
  const int kMaxOutLines = 20;
  const int kNumberCells = 5;
  const std::string kFmtStr("|% 13s");
  std::string cells[] = {"i", "x[i]", "u(x[i])", "v(x[i])",
                         "u - v"};
  if (task_id == 2) {
    cells[2] = "v(x[i])";
    cells[3] = "v2(x[i])";
    cells[4] = "v - v2";
  }

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
      printf("|% 13d|% 13e|% 13e|% 13e|% 13e|\n",
             i,
             nodes[i],
             robust_func[i],
             numerical_func[i],
             robust_func[i] - numerical_func[i]);
    }
  } else {
    for (int i = 0; i < kMaxOutLines / 2; ++i) {
      printf("|% 13d|% 13e|% 13e|% 13e|% 13e|\n",
             i,
             nodes[i],
             robust_func[i],
             numerical_func[i],
             robust_func[i] - numerical_func[i]);
    }
    for (int i = 0; i < kNumberCells; ++i) {
      printf("|");
      for (int j = 0; j < 13; ++j) {
        printf(".");
      }
    }
    printf("|\n");
    for (int i = nodes.size() - 1 - kMaxOutLines / 2; i < nodes.size(); ++i) {
      printf("|% 13d|% 13e|% 13e|% 13e|% 13e|\n",
             i,
             nodes[i],
             robust_func[i],
             numerical_func[i],
             robust_func[i] - numerical_func[i]);
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
            1.68388261 * exp(-2.0 / sqrt(5) * x) + 1;
    default: return 0;
  }
}

void PrintAbout() {
  std::cout << "Heat equation\n"
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
               "     = 1, x in (0.5, 1)\n"
               "\n"
               "Main task [-t 2]:\n"
               "(k(x)u'(x))' - q(x)u(x) = -f(x)\n"
               "u(0) = 0, u(1) = 1\n"
               "details:\n"
               "k(x) = (x + 1)^2, x in (0, 0.5]\n"
               "     = x^2 + 1, x in (0.5, 1)\n"
               "q(x) = exp(0.5 - x), x in (0, 1)\n"
               "     = exp(x - 0.5), x in (0.5, 1)\n"
               "f(x) = cos(pi * x), x in (0, 0.5]\n"
               "     = 1, x in (0.5, 1)\n"
               "\n"
               "Experiment №1 [-e 0]\n"
               "Visualize, how number of nodes affects on accuracy\n"
               "[-t] - task id \n"
               "[-n] - maximal number of intervals\n"
               "[-s] - step"
               "\n"
               "Experiment №2 [-e 1]\n"
               "When accuracy is wrong\n"
               "[-t] - task id \n"
               "[-m] - multiplication factor "
               "(start number of intervals 2, next multiplies)" << std::endl;
}
