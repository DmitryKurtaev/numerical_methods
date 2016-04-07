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
    state_dim_(state_dim), init_step_(step) {
}

void AbstractSolver::Solve(const std::vector<double>& init_state,
                           double init_point, double right_border,
                           unsigned max_n_iters) {
  Reset();
  points_.push_back(init_point);
  states_.push_back(init_state);

  std::vector<double> next_state;
  for (unsigned i = 0; i < max_n_iters && init_point < right_border; ++i) {
    Step(init_point, states_.back(), &next_state, step_, &init_point);
    points_.push_back(init_point);
    states_.push_back(next_state);
  }
}

void AbstractSolver::SolveWithLocalErrorControl(
    const std::vector<double>& init_state, double init_point,
    double right_border, unsigned max_n_iters, double eps) {
  Reset();
  points_.push_back(init_point);
  states_.push_back(init_state);
  step_inc_history_.push_back(0);
  step_dec_history_.push_back(0);
  for (unsigned i = 0; i < max_n_iters && points_.back() < right_border; ++i) {
    step_inc_history_.push_back(step_inc_history_.back());
    step_dec_history_.push_back(step_inc_history_.back());
    StepWithLocalErrorControl(eps, &step_inc_history_.back(),
                              &step_dec_history_.back());
  }
}

void AbstractSolver::StepWithLocalErrorControl(double eps, unsigned* n_step_inc,
                                               unsigned* n_step_dec) {
  std::vector<double> state_full_step;
  Step(points_.back(), states_.back(), &state_full_step, step_);

  std::vector<double> state_half_step;
  std::vector<double> state_double_half_step;
  double point;
  Step(points_.back(), states_.back(), &state_half_step, step_ / 2, &point);
  Step(point, state_half_step, &state_double_half_step, step_ / 2, &point);

  double local_error = 0;
  for (unsigned i = 0; i < state_dim_; ++i) {
    double diff = state_double_half_step[i] - state_full_step[i];
    local_error += diff * diff;
  }
  local_error = sqrt(local_error);
  const unsigned term = 1 << GetOrder();
  local_error /= (term - 1);

  if (local_error > eps) {
    step_ /= 2;
    *n_step_dec += 1;
    StepWithLocalErrorControl(eps, n_step_inc, n_step_dec);
  } else {
    states_.push_back(state_full_step);
    states_double_half_step_.push_back(state_double_half_step);
    points_.push_back(point);
    local_errors_.push_back(local_error);
    if (local_error < eps / (term * 2)) {
      step_ *= 2;
      *n_step_inc += 1;
    }
  }
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
    states_double_half_step_.insert(states_double_half_step_.begin(),
                                    std::vector<double>(state_dim_, 0));
  }
  if (!local_errors_.empty()) {
    row.push_back("LE est.");
    local_errors_.insert(local_errors_.begin(), 0);
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
  std::ostringstream sc_ss;
  sc_ss << std::scientific;
  for (unsigned i = 0; i < n_points; ++i) {
    row.clear();
    ss << i;             row.push_back(ss.str()); ss.str("");
    ss << points_[i];    row.push_back(ss.str()); ss.str("");
    ss << states_[i][0]; row.push_back(ss.str()); ss.str("");

    if (!states_double_half_step_.empty()) {
      ss << states_double_half_step_[i][0];
      row.push_back(ss.str()); ss.str("");
      double diff = states_[i][0] - states_double_half_step_[i][0];
      sc_ss << diff; row.push_back(sc_ss.str()); sc_ss.str("");
    }

    if (!local_errors_.empty()) {
      sc_ss << local_errors_[i]; row.push_back(sc_ss.str()); sc_ss.str("");
    }

    double step = (i != 0 ? points_[i] - points_[i - 1] : init_step_);
    sc_ss << step; row.push_back(sc_ss.str()); sc_ss.str("");

    if (!local_errors_.empty()) {
      ss << step_inc_history_[i]; row.push_back(ss.str()); ss.str("");
      ss << step_dec_history_[i]; row.push_back(ss.str()); ss.str("");
    }

    if (GetRobustValues_ != 0) {
      ss << robust_values[i]; row.push_back(ss.str()); ss.str("");
      ss << fabs(robust_values[i] - states_[i][0]);
      row.push_back(ss.str()); ss.str("");
    }

    table.push_back(row);
  }
  TablePrinter::Print(table, 20, 11);
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
    plot.Show("Solution", "x", "U(x)");
  }
  if (!local_errors_.empty()) {
    plot.Clear();
    plot.Add(points_, local_errors_, 3, 0.6, 0.3, 0.6, true);
    plot.Show("Local error", "node", "local error");
  }
}

void AbstractSolver::Reset() {
  points_.clear();
  states_.clear();
  states_double_half_step_.clear();
  local_errors_.clear();
  step_inc_history_.clear();
  step_dec_history_.clear();
}
