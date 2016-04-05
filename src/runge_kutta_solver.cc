#include "include/runge_kutta_solver.h"

#include <iostream>

RungeKuttaSolver::RungeKuttaSolver(
    void (*GetRightPart)(const std::vector<double>&,
                         const std::vector<double>&,
                         std::vector<double>*),
    double step, unsigned point_dim, unsigned state_dim)
  : AbstractSolver(GetRightPart, step, point_dim, state_dim),
    eulerian_solver_(GetRightPart, step, point_dim, state_dim) {
}

void RungeKuttaSolver::Step(const std::vector<double>& point,
                            const std::vector<double>& state,
                            std::vector<double>* next_state,
                            std::vector<double>* next_point) {
  if (point.size() != point_dim_)
    std::cout << "Passed point has unwanted dimension" << std::endl;
  if (state.size() != state_dim_)
    std::cout << "Passed state has unwanted dimension" << std::endl;

  next_state->resize(state_dim_);
  next_point->resize(point_dim_);
  GetRightPart_(point, state, next_state);

  std::vector<double> eulerian_solution;
  eulerian_solver_.Step(point, state, &eulerian_solution, next_point);

  std::vector<double> term;
  GetRightPart_(*next_point, eulerian_solution, &term);

  for (unsigned i = 0; i < state_dim_; ++i) {
    next_state->operator [](i) = state[i] + 0.5 * step_ *
                                 (next_state->operator [](i) + term[i]);
  }
}
