// Implements Dirichlet's task on combined area.
// 1-----------+
// |  block_3  |
// |           |
// 0.625-------+--------+
// |           |        |
// |           |        |
// |  block_1  | block_2|
// |           |        |
// 0-----------0.75-----1
// Better if borders are decomposing to degrees of 2:
// 0.75 = 1/2 + 1/4
// 0.625 = 1/2 + 1/8

#include <stdio.h>

#include <cmath>
#include <climits>
#include <cfloat>
#include <sstream>

#include "include/command_line_parser.h"
#include "include/dirichlet_task.h"
#include "include/table_printer.h"

enum Task { TEST, MAIN };

enum Block { FIRST, SECOND, THIRD };

static const double kGlobalBorders[4] = { 1.0, 1.0, 0.0, 0.0 };
static double kBlockBorders[3][4] = { { 0.625, 0.75, 0.0, 0.0 },
                                      { 0.625, 1.0, 0.0, 0.75 },
                                      { 1.0, 0.75, 0.625, 0.0 } };
double GetExternalHeat(double x, double y, Task task);

double GetBorderCondition(Border border, Block block, double coord, Task task);

DirichletTask* BuildDirichletTask(Task task, Block block, int n, int m);

bool BuildCombinedTask(std::vector<DirichletTask*>& dirichlet_tasks, Task task,
                       int n, int m);

void Print(int n_intervals_by_x, int n_intervals_by_y,
           std::vector<double>& robust_values,
           std::vector<double>& result, Task task,
           int n_processed_iters, double achieved_eps, bool print_tables);

void ShortPrint(int n, int m, Task task, int n_processed_iters,
                double achieved_eps, double max_diff);

void Solve(std::vector<DirichletTask*>& dirichlet_tasks, int& n_processed_iters,
           double& achieved_eps, int max_n_iters = 1, double target_eps = 0);

void Merge(std::vector<double>* blocks,
           std::vector<double>& result,
           int* blocks_n, int* blocks_m, int n, int m);

int main(int argc, char** argv) {
  CommandLineParser parser(argc, argv);
  if (parser.Exists("h") || argc == 1) {
    // PrintAbout();
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
  const double h = (kGlobalBorders[RIGHT] - kGlobalBorders[LEFT]) / n;
  const double k = (kGlobalBorders[TOP] - kGlobalBorders[BOTTOM]) / m;

  std::vector<DirichletTask*> dirichlet_tasks;
  bool success = BuildCombinedTask(dirichlet_tasks, task, n, m);
  if (!success) return 1;

  int n_processed_iters;
  double achieved_eps;
  Solve(dirichlet_tasks, n_processed_iters, achieved_eps, n_iters, eps);

  // Extract states of tasks.
  std::vector<double> local_results[3];
  int blocks_m[3];
  int blocks_n[3];
  for (int i = 0; i < 3; ++i) {
    dirichlet_tasks[i]->GetState(local_results[i]);
    dirichlet_tasks[i]->GetDimensions(blocks_n[i], blocks_m[i]);
    delete dirichlet_tasks[i];
  }

  if (task == TEST) {
    std::vector<double> robust_values((n + 1) * (m + 1));
    for (int j = 0; j <= m; ++j) {
      const double y = j * k;
      const int offset = j * (n + 1);
      for (int i = 0; i <= n; ++i) {
        robust_values[offset + i] = exp(pow(sin(M_PI * i * h * y), 2));
      }
    }

    // Remove values outside computation area.
    for (int j = blocks_m[SECOND] + 1; j <= m; ++j) {
      const int offset = j * (n + 1);
      for (int i = blocks_n[FIRST] + 1; i <= n; ++i) {
        robust_values[offset + i] = nan("");
      }
    }

    std::vector<double> result;
    result.reserve((n - 1) * (m - 1));
    std::vector<double>::iterator first_it = local_results[FIRST].begin();
    std::vector<double>::iterator second_it = local_results[SECOND].begin();
    for (int i = 0; i < blocks_m[FIRST] - 1; ++i) {
      result.insert(result.end(), first_it, first_it + blocks_n[FIRST] - 1);
      result.insert(result.end(), second_it, second_it + blocks_n[SECOND] - 1);
      first_it += blocks_n[FIRST] - 1;
      second_it += blocks_n[SECOND] - 1;
    }

    std::vector<double>::iterator third_it = local_results[THIRD].begin();
    for (int i = 0; i < blocks_m[THIRD]; ++i) {
      result.insert(result.end(), third_it, third_it + blocks_n[THIRD] - 1);

      int offset = (i + blocks_m[FIRST]) * (n + 1) + blocks_n[THIRD];
      result.insert(result.end(), robust_values.begin() + offset,
                    robust_values.begin() + offset + blocks_n[SECOND] - 1);
      third_it += blocks_n[THIRD] - 1;
    }

    Print(n, m, robust_values, result, TEST, n_processed_iters, achieved_eps,
          print_tables);
  } else {
    std::vector<DirichletTask*> dense_dirichlet_tasks;
    BuildCombinedTask(dense_dirichlet_tasks, MAIN, 2 * n, 2 * m);
    Solve(dense_dirichlet_tasks, n_processed_iters, achieved_eps, n_iters, eps);

    std::vector<double> dense_results[3];
    
    int dense_blocks_m[3];
    int dense_blocks_n[3];
    for (int k = 0; k < 3; ++k) {
      dense_dirichlet_tasks[k]->GetState(dense_results[k]);
      dense_dirichlet_tasks[k]->GetDimensions(dense_blocks_n[k],
                                              dense_blocks_m[k]);
      delete dense_dirichlet_tasks[k];
    }

    // Collect values.
    std::vector<double> dense_result;
    Merge(dense_results, dense_result, dense_blocks_n, dense_blocks_m, 2 * n,
          2 * m);

    // Remove each 2th row and column (starts from 1).
    int row = 0;
    for (int j = 0; j < 2 * m - 1; ++j) {
      if (j % 2 == 0) {
        for (int i = 0; i < 2 * n - 1; ++i) {
          dense_result.erase(dense_result.begin() + row * (n - 1));
        }
      } else {
        for (int i = 0; i < n; ++i) {
          dense_result.erase(dense_result.begin() + row * (n - 1) + i);
        }
        ++row;
      }
    }
    std::vector<double> origin_result;
    Merge(local_results, origin_result, blocks_n, blocks_m, n, m);

    double max_diff = 0;
    const int dim = origin_result.size();
    for (int i = 0; i < dim; ++i) {
      double diff = fabs(origin_result[i] - dense_result[i]);
      if (diff > max_diff) max_diff = diff;
    }
    ShortPrint(n, m, MAIN, n_processed_iters, achieved_eps, max_diff);
  }
  return 0;
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

double GetBorderCondition(Border border, Block block, double coord, Task task) {
  switch (task) {
    case TEST:
      switch (block) {
        case FIRST: {
          switch (border) {
            case LEFT: case BOTTOM: return 1;
            default: return 0;
          }
        }
        case SECOND: {
          switch (border) {
            case TOP: return exp(pow(sin(M_PI * 0.625 * coord), 2));
            case RIGHT: return exp(pow(sin(M_PI * coord), 2));
            case BOTTOM: return 1;
            default: return 0;
          }
        }
        case THIRD: {
          switch (border) {
            case TOP: return exp(pow(sin(M_PI * coord), 2));
            case LEFT: return 1;
            case RIGHT: return exp(pow(sin(M_PI * 0.75 * coord), 2));
            default: return 0;
          }
        }
      }
      break;
    case MAIN:
      switch (block) {
        case FIRST: {
          switch (border) {
            case LEFT: return sin(M_PI * coord);
            case BOTTOM: return coord * (1 - coord);
            default: return 0;
          }
        }
        case SECOND: {
          switch (border) {
            case TOP: return coord * sin(0.625 * M_PI);
            case RIGHT: return sin(M_PI * coord);
            case BOTTOM: return coord * (1 - coord);
            default: return 0;
          }
        }
        case THIRD: {
          switch (border) {
            case TOP: return coord * (1 - coord);
            case LEFT: return sin(M_PI * coord);
            case RIGHT: return 2 * sin(0.625 * M_PI) * (1 - coord) + 
                               0.5 * coord - 0.3125;
            default: return 0;
          }
        }
      }
      break;
    default: return 0;
  }
}

DirichletTask* BuildDirichletTask(Task task, Block block, int n, int m) {
  // Fill external heat for task.
  const double top = kBlockBorders[block][TOP];
  const double left = kBlockBorders[block][LEFT];
  const double right = kBlockBorders[block][RIGHT];
  const double bottom = kBlockBorders[block][BOTTOM];
  const double h = (right - left) / n;
  const double k = (top - bottom) / m;
  double* external_heat = new double[(m - 1) * (n - 1)];
  for (int j = 1; j < m; ++j) {
    const int offset = (j - 1) * (n - 1);
    const double y = bottom + j * k;
    for (int i = 1; i < n; ++i) {
      const double x = left + i * h;
      external_heat[offset + i - 1] = GetExternalHeat(x, y, task);
    }
  }

  DirichletTask* dirichlet_task = new DirichletTask(left, right, top, bottom,
                                                    n, m, external_heat);
  delete[] external_heat; external_heat = 0;

  // Set borders.
  double* first_border = new double[n - 1];
  double* second_border = new double[n - 1];
  for (int i = 1; i < n; ++i) {
    const double x = left + i * h;
    first_border[i - 1] = GetBorderCondition(BOTTOM, block, x, task);
    second_border[i - 1] = GetBorderCondition(TOP, block, x, task);
  }
  dirichlet_task->UpdateBorder(BOTTOM, first_border);
  dirichlet_task->UpdateBorder(TOP, second_border);
  delete[] first_border;
  delete[] second_border;

  first_border = new double[m - 1];
  second_border = new double[m - 1];
  for (int j = 1; j < m; ++j) {
    const double y = bottom + j * k;
    first_border[j - 1] = GetBorderCondition(LEFT, block, y, task);
    second_border[j - 1] = GetBorderCondition(RIGHT, block, y, task);
  }
  dirichlet_task->UpdateBorder(LEFT, first_border);
  dirichlet_task->UpdateBorder(RIGHT, second_border);
  delete[] first_border;
  delete[] second_border;

  return dirichlet_task;
}

bool BuildCombinedTask(std::vector<DirichletTask*>& dirichlet_tasks, Task task,
                       int n, int m) {
  dirichlet_tasks.resize(3);

  // Check that blocks border under net's node.
  const double h = (kGlobalBorders[RIGHT] - kGlobalBorders[LEFT]) / n;
  const double k = (kGlobalBorders[TOP] - kGlobalBorders[BOTTOM]) / m;

  for (int border = 0; border < 4; ++border) {
    const double step = (border % 2 ? h : k);
    for (int block = 0; block < 3; ++block) {
      const int idx = kBlockBorders[block][border] / step;
      if (idx * step != kBlockBorders[block][border]) {
        std::cout << "Block #" << block << " has border "
                  << kBlockBorders[block][border] << " which not on net with "
                  << "step " << step << ". (" << kBlockBorders[block][border]
                  << " in [" << idx * step << ", " << (idx + 1) * step <<
                  "])" << std::endl;
        return false;
      }
    }
  }

  // Split nodes to blocks. Let block_2 & block_3 has own borders.
  int blocks_n[3];
  blocks_n[FIRST] = kBlockBorders[FIRST][RIGHT] / h;
  blocks_n[THIRD] = blocks_n[FIRST];
  blocks_n[SECOND] = n - blocks_n[FIRST] + 1;

  int blocks_m[3];
  blocks_m[FIRST] = kBlockBorders[FIRST][TOP] / k;
  blocks_m[SECOND] = blocks_m[FIRST];
  blocks_m[THIRD] = m - blocks_m[FIRST] + 1;

  kBlockBorders[THIRD][BOTTOM] -= k;
  kBlockBorders[SECOND][LEFT] -= h;

  for (int i = 0; i < 3; ++i) {
    dirichlet_tasks[i] = BuildDirichletTask(task, Block(i), blocks_n[i],
                                            blocks_m[i]);
  }
  kBlockBorders[THIRD][BOTTOM] += k;
  kBlockBorders[SECOND][LEFT] += h;

  dirichlet_tasks[FIRST]->UpdateBorder(RIGHT, *dirichlet_tasks[SECOND]);
  dirichlet_tasks[SECOND]->UpdateBorder(LEFT, *dirichlet_tasks[FIRST]);
  dirichlet_tasks[FIRST]->UpdateBorder(TOP, *dirichlet_tasks[THIRD]);
  dirichlet_tasks[THIRD]->UpdateBorder(BOTTOM, *dirichlet_tasks[FIRST]);
  return true;
}

void Solve(std::vector<DirichletTask*>& dirichlet_tasks, int& n_processed_iters,
           double& achieved_eps, int max_n_iters, double target_eps) {
  achieved_eps = DBL_MAX;
  for (n_processed_iters = 0;
       n_processed_iters < max_n_iters && achieved_eps > target_eps;
       ++n_processed_iters) {

    dirichlet_tasks[0]->MPIIteration(achieved_eps);
    for (int i = 1; i < 3; ++i) {
      double current_eps;
      dirichlet_tasks[i]->MPIIteration(current_eps);
      achieved_eps = std::max(achieved_eps, current_eps);
    }

    dirichlet_tasks[FIRST]->UpdateBorder(RIGHT, *dirichlet_tasks[SECOND]);
    dirichlet_tasks[SECOND]->UpdateBorder(LEFT, *dirichlet_tasks[FIRST]);
    dirichlet_tasks[FIRST]->UpdateBorder(TOP, *dirichlet_tasks[THIRD]);
    dirichlet_tasks[THIRD]->UpdateBorder(BOTTOM, *dirichlet_tasks[FIRST]);
    
    printf("\rProcessed iterations: %d, max|x[s+1]-x[s]| = %e",
           n_processed_iters + 1, achieved_eps);
    fflush(stdout);
  }
  std::cout << std::endl;
}

void Print(int n_intervals_by_x, int n_intervals_by_y,
           std::vector<double>& robust_values,
           std::vector<double>& result, Task task,
           int n_processed_iters, double achieved_eps, bool print_tables) {
  const int n = n_intervals_by_x;
  const int m = n_intervals_by_y;
  const double h = (kGlobalBorders[RIGHT] - kGlobalBorders[LEFT]) / n;
  const double k = (kGlobalBorders[TOP] - kGlobalBorders[BOTTOM]) / m;

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
  data[0][0] = (task == TEST? "U" : "V2");
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
  TablePrinter::Print(data, 20, 11, 10);

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
  TablePrinter::Print(data, 20, 8, 12);

  // Accuracy.
  data[0][0] = "Diff.";
  for (int i = 0; i < n - 1; ++i) {
    for (int j = 0; j < m - 1; ++j) {
      std::ostringstream ss;
      ss << std::scientific
         << fabs(result[(m - 2 - j) * (n - 1) + i] -
                 robust_values[(m - j - 1) * (n + 1) + i + 1]);
      data[j + 1][i + 1] = ss.str();
    }
  }
  TablePrinter::Print(data, 20, 8, 12);
}

void ShortPrint(int n, int m, Task task, int n_processed_iters,
                double achieved_eps, double max_diff) {
  const double h = (kGlobalBorders[RIGHT] - kGlobalBorders[LEFT]) / n;
  const double k = (kGlobalBorders[TOP] - kGlobalBorders[BOTTOM]) / m;

  std::cout << "\nNet step by x: " << h << std::endl;
  std::cout << "Net step by y: " << k << std::endl;
  
  printf("Number of processed iterations: %d\n", n_processed_iters);
  printf("max|x[s+1]-x[s]| = %e\n", achieved_eps);
  std::cout << "max|V-" << (task == TEST ? "U" : "V2") << "| = " << std::flush;
  printf("%e\n", max_diff);
}

void Merge(std::vector<double>* blocks, std::vector<double>& result,
           int* blocks_n, int* blocks_m, int n, int m) {
  result.clear();
  result.reserve((n - 1) * (m - 1));
  std::vector<double>::iterator first_it = blocks[FIRST].begin();
  std::vector<double>::iterator second_it = blocks[SECOND].begin();
  for (int i = 0; i < blocks_m[FIRST] - 1; ++i) {
    result.insert(result.end(), first_it, first_it + blocks_n[FIRST] - 1);
    result.insert(result.end(), second_it, second_it + blocks_n[SECOND] - 1);
    first_it += blocks_n[FIRST] - 1;
    second_it += blocks_n[SECOND] - 1;
  }

  std::vector<double>::iterator third_it = blocks[THIRD].begin();
  std::vector<double> zeros(blocks_n[SECOND] - 1, 0);
  for (int i = 0; i < blocks_m[THIRD] - 1; ++i) {
    result.insert(result.end(), third_it, third_it + blocks_n[THIRD] - 1);
    result.insert(result.end(), zeros.begin(), zeros.end());
    third_it += blocks_n[THIRD] - 1;
  }
}
