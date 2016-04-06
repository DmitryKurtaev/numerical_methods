#include <vector>

#include "opencv2/opencv.hpp"

#include "include/runge_kutta_solver.h"
#include "include/plot.h"
#include "include/table_printer.h"

enum Task { TEST, FIRST_MAIN, SECOND_MAIN };

static const double kLeftBorder = 0.0;

const char* kCmdParams =
    "{ h | help | false | Print this message }"
    "{ t | task | 0 | Id of task (0-test, 1-first main, 2-second main) }"
    "{ n | n_iters | 16 | Number of iterations }"
    "{ eps | epsilon | 0.01 | Parameter for local error control method }"
    "{ ec | err_control | false | Use local error control }"
    "{ b | right_point | 1.0 | Right border of interval }"
    "{ is | init_state | 0.0 | Initial state (at left point) }"
    "{ id | init_deriv | 0.0 | Initial derivation (at left point }";

void GetTestRightPart(double point, const std::vector<double>& state,
                      std::vector<double>* derivation);

void GetFirstMainRightPart(double point, const std::vector<double>& state,
                           std::vector<double>* derivation);

void GetSecondMainRightPart(double point, const std::vector<double>& state,
                            std::vector<double>* derivation);

void GetTestRobustValues(const std::vector<double>& points,
                         const std::vector<double>& initial_state,
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
  const double eps = parser.get<double>("epsilon");
  if (task == FIRST_MAIN && initial_state > 0.6) {
    std::cout << "Reccomends U(0) is lower or equals 0.6" << std::endl;
  }

  AbstractSolver* solver = 0;
  switch (task) {
    case TEST:
      solver = new RungeKuttaSolver(GetTestRightPart, step, 1,
                                    GetTestRobustValues);
      break;
    case FIRST_MAIN:
      solver = new RungeKuttaSolver(GetFirstMainRightPart, step, 1);
      break;
    default:
      std::cout << "Unknows task" << std::endl;
      return 0;
  }

  // Method.
  std::vector<double> init_state(1, initial_state);
  if (use_error_control) {
//    unsigned step_increasing_counter = 0;
//    unsigned step_decreasing_counter = 0;
//    double local_error;
//    for (unsigned i = 0; i < n_iters && point <= right_border; ++i) {
//      solver->StepWithLocalErrorControl(
//            point, std::vector<double>(next_state), &next_state, eps,
//            &local_error, &step_increasing_counter, &step_decreasing_counter,
//            &point);
//      points.push_back(point);
//      states.push_back(next_state[0]);
//      local_errors.push_back(local_error);
//      step_increasing_counters.push_back(step_increasing_counter);
//      step_decreasing_counters.push_back(step_decreasing_counter);
//    }
  } else {
    solver->Solve(init_state, kLeftBorder, right_border, n_iters);
  }
  solver->ShowResults();

  delete solver;

  return 0;
}

void GetTestRightPart(double point, const std::vector<double>& state,
                      std::vector<double>* derivation) {
  derivation->resize(1);
  derivation->operator [](0) = -5.5 * state[0];
}

void GetFirstMainRightPart(double point, const std::vector<double>& state,
                           std::vector<double>* derivation) {
  derivation->resize(1);
  const double u = state[0];
  const double term = (1.0 + pow(point, 3)) / (1.0 + pow(point, 5));
  derivation->operator [](0) = u * (1.0 + u * (term - u * sin(10 * point)));
}

void GetSecondMainRightPart(double point, const std::vector<double>& state,
                            std::vector<double>* derivation) {

}

void GetTestRobustValues(const std::vector<double>& points,
                         const std::vector<double>& initial_state,
                         std::vector<double>* states) {
  const unsigned size = points.size();
  states->resize(size);
  for (unsigned i = 0; i < size; ++i) {
    states->operator [](i) = initial_state[0] * exp(-5.5 * points[i]);
  }
}
