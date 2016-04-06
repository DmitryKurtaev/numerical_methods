#include "include/abstract_solver.h"

AbstractSolver::AbstractSolver(
    void (*GetRightPart)(double point,
                         const std::vector<double>& state,
                         std::vector<double>* derivation),
    double step, unsigned state_dim)
  : GetRightPart_(GetRightPart), step_(step), state_dim_(state_dim) {
}

void AbstractSolver::StepWithLocalErrorControl(
    double point,
    const std::vector<double>& state,
    std::vector<double>* next_state,
    double eps,
    double* local_error,
    unsigned* step_increasing_counter,
    unsigned* step_decreasing_counter,
    double* next_point) {
  std::vector<double> temp_state_full_step;
  Step(point, state, temp_state_full_step);

  std::vector<double> temp_state_half_step;
  std::vector<double> temp_state_double_half_step;
  step_ /= 2;
  Step(point, state, temp_state_half_step, &point);
  Step(point, temp_state_half_step, temp_state_double_half_step, next_point);
  step_ *= 2;

  *local_error = 0;
  for (unsigned i = 0; i < state_dim_; ++i) {
    double diff = temp_state_double_half_step[i] - temp_state_full_step[i];
    *local_error += diff * diff;
  }
  *local_error = sqrt(*local_error);
  const unsigned term = 1 << GetOrder();
  *local_error /= (term - 1);

  if (*local_error > eps) {
    step_ /= 2;
    *step_decreasing_counter += 1;
    StepWithLocalErrorControl(point, state, next_state, eps, local_error,
                              step_increasing_counter, step_decreasing_counter,
                              next_point);
    return;
  } else {
    next_state->resize(state_dim_);
    for (unsigned i = 0; i < state_dim_; ++i) {
      next_state->operator [](i) = temp_state_full_step[i];
    }
    if (*local_error < eps / (term * 2)) {
      step_ *= 2;
      *step_increasing_counter += 1;
    }
  }
}
