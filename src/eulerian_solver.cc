#include "include/eulerian_solver.h"

#include <iostream>

EulerianSolver::EulerianSolver(
    void (*GetRightPart)(const std::vector<double>& point,
                         const std::vector<double>& state,
                         std::vector<double>* derivation),
    double step, unsigned point_dim, unsigned state_dim)
  : AbstractSolver(GetRightPart, step, point_dim, state_dim) {
}

void EulerianSolver::Step(const std::vector<double> &point,
                          const std::vector<double> &state,
                          std::vector<double>* next_state,
                          std::vector<double>* next_point) {
  if (point.size() != point_dim_)
    std::cout << "Passed point has unwanted dimension" << std::endl;
  if (state.size() != state_dim_)
    std::cout << "Passed state has unwanted dimension" << std::endl;

  next_state->resize(state_dim_);
  GetRightPart_(point, state, next_state);
  for (unsigned i = 0; i < state_dim_; ++i) {
    next_state->operator [](i) += state[i] + step_ * next_state->operator [](i);
  }

  if (next_point) {
    next_point->resize(point_dim_);
    for (unsigned i = 0; i < point_dim_; ++i) {
      next_point->operator [](i) = point[i] + step_;
    }
  }
}
