#include <vector>

#include "opencv2/opencv.hpp"

#include "include/runge_kutta_solver.h"
#include "include/plot.h"

enum Task { TEST, FIRST_MAIN, SECOND_MAIN };

static const double kLeftBorder = 0.0;

const char* kCmdParams =
    "{ h | help | false | Print this message }"
    "{ t | task | 0 | Id of task (0-test, 1-first main, 2-second main) }"
    "{ n | n_iters | 16 | Number of iterations }"
    "{ ec | err_control | false | Use local error control }"
    "{ b | right_point | 1.0 | Right border of interval }"
    "{ is | init_state | 0.0 | Initial state (at left point) }"
    "{ id | init_deriv | 0.0 | Initial derivation (at left point }";

void GetTestRightPart(const std::vector<double>& point,
                      const std::vector<double>& state,
                      std::vector<double>* derivation);

void GetFirstMainRightPart(const std::vector<double>& point,
                           const std::vector<double>& state,
                           std::vector<double>* derivation);

void GetSecondMainRightPart(const std::vector<double>& point,
                            const std::vector<double>& state,
                            std::vector<double>* derivation);

void GetTestRobustValues(const std::vector<double>& points,
                         double initial_state,
                         std::vector<double>* states);

int main(int argc, char** argv) {
  cv::CommandLineParser parser(argc, argv, kCmdParams);
  if (argc == 1 || parser.get<bool>("help")) {
    parser.printParams();
    return 0;
  }

  const unsigned n_iters = parser.get<unsigned>("n_iters");
  const bool use_error_control = parser.get<bool>("err_control");
  const double right_border = parser.get<double>("right_point");
  const double initial_state = parser.get<double>("init_state");
  const double initial_derivation = parser.get<double>("init_state");
  const Task task = static_cast<Task>(parser.get<unsigned>("task"));
  const double step = (right_border - kLeftBorder) / n_iters;

  AbstractSolver* solver = 0;
  switch (task) {
    case TEST:
      solver = new RungeKuttaSolver(GetTestRightPart, step, 1, 1);
      break;
    case FIRST_MAIN:
      solver = new RungeKuttaSolver(GetTestRightPart, step, 1, 1);
      break;
    default:
      std::cout << "Unknows task" << std::endl;
      return 0;
  }

  std::vector<double> points(1, kLeftBorder);  // Grid.
  std::vector<double> states(1, initial_state);  // Results of method's work.
  std::vector<double> next_point(1, kLeftBorder);
  std::vector<double> next_state(1, initial_state);
  for (unsigned i = 0; i < n_iters; ++i) {
    solver->Step(std::vector<double>(next_point),
                 std::vector<double>(next_state),
                 &next_state,
                 &next_point);
    points.push_back(next_point[0]);
    states.push_back(next_state[0]);
  }

  Plot plot;
  plot.Add(points, states, 2, 0.5, 0.8, 0.4, true);
  plot.Show("Test", "x", "U(x)");

  delete solver;

  return 0;
}

void GetTestRightPart(const std::vector<double>& point,
                      const std::vector<double>& state,
                      std::vector<double>* derivation) {
  derivation->resize(1);
  derivation->operator [](0) = -5.5 * state[0];
}

void GetFirstMainRightPart(const std::vector<double>& point,
                           const std::vector<double>& state,
                           std::vector<double>* derivation) {
  derivation->resize(1);
  const double x = point[0];
  const double u = state[0];
  const double term = (1.0 + pow(x, 3)) / (1.0 + pow(x, 5));
  derivation->operator [](0) = u * (1.0 + u * (term - u * sin(10 * x)));
}

void GetSecondMainRightPart(const std::vector<double>& point,
                            const std::vector<double>& state,
                            std::vector<double>* derivation) {

}

void GetTestRobustValues(const std::vector<double>& points,
                         double initial_state,
                         std::vector<double>* states) {
  const unsigned size = points.size();
  states->resize(size);
  for (unsigned i = 0; i < size; ++i) {
    states->operator [](i) = initial_state * exp(-5.5 * points[i]);
  }
}
