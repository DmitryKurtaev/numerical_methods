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

#include <cmath>
#include <climits>

#include "include/command_line_parser.h"
#include "include/dirichlet_task.h"

enum Task { TEST, MAIN };

enum Block { FIRST, SECOND, THIRD };

static const double kGlobalBorders[4] = { 1.0, 1.0, 0.0, 0.0 };
static double kBlockBorders[3][4] = { { 0.625, 0.75, 0.0, 0.0 },
                                      { 0.625, 1.0, 0.0, 0.75 },
                                      { 1.0, 0.75, 0.625, 0.0 } };
double GetExternalHeat(double x, double y, Task task);

double GetBorderCondition(Border border, Block block, double coord, Task task);

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
        return 1;
      }
    }
  }

  // Split nodes to blocks. Let block_1 has own borders.
  int blocks_n[3];
  blocks_n[FIRST] = 1 + kBlockBorders[FIRST][RIGHT] / h;
  blocks_n[THIRD] = blocks_n[FIRST] - 1;
  blocks_n[SECOND] = n - blocks_n[THIRD];

  int blocks_m[3];
  blocks_m[FIRST] = 1 + kBlockBorders[FIRST][TOP] / k;
  blocks_m[SECOND] = blocks_m[FIRST] - 1;
  blocks_m[THIRD] = m - blocks_m[SECOND];

  kBlockBorders[FIRST][RIGHT] += h;
  kBlockBorders[FIRST][TOP] += k;

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
    default: return 0;
  }
}