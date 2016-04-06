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
      solver = new RungeKuttaSolver(GetTestRightPart, step, 1);
      break;
    case FIRST_MAIN:
      solver = new RungeKuttaSolver(GetFirstMainRightPart, step, 1);
      break;
    default:
      std::cout << "Unknows task" << std::endl;
      return 0;
  }

  // Method.
  std::vector<double> points(1, kLeftBorder);  // Grid.
  std::vector<double> states(1, initial_state);  // Results of method's work.
  std::vector<double> next_state(1, initial_state);
  for (unsigned i = 0; i < n_iters; ++i) {
    double point = points[points.size() - 1];
    solver->Step(point, std::vector<double>(next_state), &next_state, &point);
    points.push_back(point);
    states.push_back(next_state[0]);
  }

  // Output procedures.
  Plot plot;
  std::vector<std::vector<std::string> > table(points.size() + 1);
  std::vector<std::string> row;
  std::ostringstream ss;
  switch (task) {

    case TEST: {
      row.resize(6);
      row[0] = "i"; row[1] = "x[i]"; row[2] = "V[i]"; row[3] = "h[i]";
      row[4] = "U[i]"; row[5] = "|U[i] - V[i]|";
      table[0] = row;

      std::vector<double> robust_values;
      GetTestRobustValues(points, initial_state, &robust_values);
      double max_diff = 0;
      unsigned argmax = 0;
      for (int i = 0; i < points.size(); ++i) {
        double diff = fabs(robust_values[i] - states[i]);
        if (diff > max_diff) {
          max_diff = diff;
          argmax = i;
        }
        ss << i;                row[0] = ss.str(); ss.str("");
        ss << points[i];        row[1] = ss.str(); ss.str("");
        ss << states[i];        row[2] = ss.str(); ss.str("");
        ss << step;             row[3] = ss.str(); ss.str("");
        ss << robust_values[i]; row[4] = ss.str(); ss.str("");
        ss << diff;             row[5] = ss.str(); ss.str("");
        table[i + 1] = row;
      }
      TablePrinter::Print(table);
      std::cout << "Number of iterations: " << n_iters << std::endl;
//      std::cout << "b - x[n] = " << right_border - points.back() << std::endl;
      std::cout << "max|U[i] - V[i]| = " << max_diff << " at x[" << argmax
                << "] = " << points[argmax] << std::endl;

      plot.Add(points, states, 2, 0.8, 0.5, 0.4, true);
      plot.Add(points, robust_values, 2, 0.5, 0.8, 0.4, true);
      plot.Show("Test", "x", "U(x)");
      break;
    }

    case FIRST_MAIN: {
      plot.Add(points, states, 2, 0.8, 0.5, 0.4, true);
      plot.Show("First main task", "x", "U(x)");
      break;
    }

    default: break;
  }

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
                         double initial_state,
                         std::vector<double>* states) {
  const unsigned size = points.size();
  states->resize(size);
  for (unsigned i = 0; i < size; ++i) {
    states->operator [](i) = initial_state * exp(-5.5 * points[i]);
  }
}
