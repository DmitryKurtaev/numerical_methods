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
    case TEST: solver = new RungeKuttaSolver(GetTestRightPart, step, 1, 1);
  }

  std::vector<double> points(n_iters + 1);
  std::vector<double> states(n_iters + 1);
  std::vector<double> current_point(1, kLeftBorder);
  std::vector<double> current_state(1, initial_state);
  std::vector<double> next_point;
  std::vector<double> next_state;
  points[0] = kLeftBorder;
  states[0] = initial_state;
  for (unsigned i = 0; i < n_iters; ++i) {
    solver->Step(current_point, current_state, &next_state, &next_point);
    points[i + 1] = current_point[0] = next_point[0];
    states[i + 1] = current_state[0] = next_state[0];
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
