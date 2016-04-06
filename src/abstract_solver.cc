#include "include/abstract_solver.h"

#include <math.h>

#include <string>
#include <sstream>
#include <iostream>

#include "include/plot.h"
#include "include/table_printer.h"

AbstractSolver::AbstractSolver(
    void (*GetRightPart)(double point,
                         const std::vector<double>& state,
                         std::vector<double>* derivation),
    double step, unsigned state_dim,
    void (*GetRobustValues)(const std::vector<double>& points,
                            const std::vector<double>& init_state,
                            std::vector<double>* states))
  : GetRightPart_(GetRightPart), GetRobustValues_(GetRobustValues), step_(step),
    state_dim_(state_dim) {
}

void AbstractSolver::Solve(const std::vector<double>& init_state,
                           double init_point, double right_border,
                           unsigned max_n_iters) {
  Reset();
  points_.push_back(init_point);
  states_.push_back(init_state);

  std::vector<double> next_state;
  for (unsigned i = 0; i < max_n_iters && init_point <= right_border; ++i) {
    Step(init_point, states_.back(), &next_state, &init_point);
    points_.push_back(init_point);
    states_.push_back(next_state);
  }
}

void AbstractSolver::StepWithLocalErrorControl(double eps) {
  Reset();
//  std::vector<double> temp_state_full_step;
//  Step(point, state, temp_state_full_step);

//  std::vector<double> temp_state_half_step;
//  std::vector<double> temp_state_double_half_step;
//  step_ /= 2;
//  Step(point, state, temp_state_half_step, &point);
//  Step(point, temp_state_half_step, temp_state_double_half_step, next_point);
//  step_ *= 2;

//  *local_error = 0;
//  for (unsigned i = 0; i < state_dim_; ++i) {
//    double diff = temp_state_double_half_step[i] - temp_state_full_step[i];
//    *local_error += diff * diff;
//  }
//  *local_error = sqrt(*local_error);
//  const unsigned term = 1 << GetOrder();
//  *local_error /= (term - 1);

//  if (*local_error > eps) {
//    step_ /= 2;
//    *step_decreasing_counter += 1;
//    StepWithLocalErrorControl(point, state, next_state, eps, local_error,
//                              step_increasing_counter, step_decreasing_counter,
//                              next_point);
//    return;
//  } else {
//    next_state->resize(state_dim_);
//    for (unsigned i = 0; i < state_dim_; ++i) {
//      next_state->operator [](i) = temp_state_full_step[i];
//    }
//    if (*local_error < eps / (term * 2)) {
//      step_ *= 2;
//      *step_increasing_counter += 1;
//    }
//  }
}

void AbstractSolver::ShowResults() {
  const unsigned n_points = points_.size();
  std::vector<std::vector<std::string> > table;
  std::vector<std::string> row;
  std::vector<double> robust_values;

  // Header.
  row.push_back("i");
  row.push_back("x[i]");
  row.push_back("V[i]");
  if (!states_double_half_step_.empty()) {
    row.push_back("V2[i]");
    row.push_back("V[i] - V2[i]");
  }
  if (!local_errors_.empty()) {
    row.push_back("LE est.");
  }
  row.push_back("h[i]");
  if (!local_errors_.empty()) {
    row.push_back("Num step inc");
    row.push_back("Num step dec");
  }
  if (GetRobustValues_ != 0) {
    row.push_back("U[i]");
    row.push_back("|U[i] - V[i]|");
     GetRobustValues_(points_, states_.front(), &robust_values);
  }
  table.push_back(row);

  std::ostringstream ss;
  unsigned number_step_increasing = 0;
  unsigned number_step_decreasing = 0;
  for (unsigned i = 0; i < n_points; ++i) {
    row.clear();
    ss << i;             row.push_back(ss.str()); ss.str("");
    ss << points_[i];    row.push_back(ss.str()); ss.str("");
    ss << states_[i][0]; row.push_back(ss.str()); ss.str("");

    if (!states_double_half_step_.empty()) {
      ss << states_double_half_step_[i][0]; row.push_back(ss.str()); ss.str("");
      double diff = states_[i][0] - states_double_half_step_[i][0];
      ss << diff; row.push_back(ss.str()); ss.str("");
    }
    if (!local_errors_.empty()) {
      ss << local_errors_[i]; row.push_back(ss.str()); ss.str("");
    }

    double step = (i != 0 ? points_[i] - points_[i - 1] :
                            points_[1] - points_[0]);
    ss << step; row.push_back(ss.str()); ss.str("");

    if (!local_errors_.empty()) {
      double previous_step = (i > 1 ? points_[i - 1] - points_[i - 2] :
                                      points_[1] - points_[0]);
      if (step > previous_step) {
        ss << ++number_step_increasing; row.push_back(ss.str()); ss.str("");
      } else if (step < previous_step) {
        ss << ++number_step_decreasing; row.push_back(ss.str()); ss.str("");
      }
    }

    if (GetRobustValues_ != 0) {
      ss << robust_values[i]; row.push_back(ss.str()); ss.str("");
      ss << fabs(robust_values[i] - states_[i][0]);
      row.push_back(ss.str()); ss.str("");
    }

    table.push_back(row);
  }
  TablePrinter::Print(table);
  std::cout << "Number of iterations: " << n_points - 1 << std::endl;

  Plot plot;
  if (state_dim_ == 1) {
    std::vector<double> single_states(n_points);
    for (unsigned i = 0; i < n_points; ++i) {
      single_states[i] = states_[i][0];
    }
    if (GetRobustValues_ != 0) {
      plot.Add(points_, robust_values, 2, 0.5, 0.8, 0.4, true);
    }
    plot.Add(points_, single_states, 2, 0.8, 0.5, 0.4, true);
    plot.Show("", "x", "U(x)");
  }
}

void AbstractSolver::Reset() {
  points_.clear();
  states_.clear();
  states_double_half_step_.clear();
  local_errors_.clear();
}
