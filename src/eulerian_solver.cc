#include "include/eulerian_solver.h"

#include <iostream>

EulerianSolver::EulerianSolver(
    void (*GetRightPart)(double point,
                         const std::vector<double>& state,
                         std::vector<double>* derivation),
    double step, unsigned state_dim)
  : AbstractSolver(GetRightPart, step, state_dim) {
}

void EulerianSolver::Step(double point,
                          const std::vector<double>& state,
                          std::vector<double>* next_state,
                          double* next_point) {
  if (state.size() != state_dim_)
    std::cout << "Passed state has unwanted dimension" << std::endl;

  next_state->resize(state_dim_);
  GetRightPart_(point, state, next_state);
  for (unsigned i = 0; i < state_dim_; ++i) {
    next_state->operator [](i) = state[i] + step_ * next_state->operator [](i);
  }
  if (next_point) *next_point = point + step_;
}
