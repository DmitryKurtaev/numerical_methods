#include "include/runge_kutta_solver.h"

#include <iostream>

RungeKuttaSolver::RungeKuttaSolver(
    void (*GetRightPart)(double point,
                         const std::vector<double>& state,
                         std::vector<double>* derivation),
    double step, unsigned state_dim,
    void (*GetRobustValues)(const std::vector<double>& points,
                            const std::vector<double>& init_state,
                            std::vector<double>* states))
  : AbstractSolver(GetRightPart, step, state_dim, GetRobustValues),
    eulerian_solver_(GetRightPart, step, state_dim, GetRobustValues) {
}

void RungeKuttaSolver::Step(double point,
                            const std::vector<double>& state,
                            std::vector<double>* next_state,
                            double* next_point) {
  if (state.size() != state_dim_)
    std::cout << "Passed state has unwanted dimension" << std::endl;

  next_state->resize(state_dim_);
  GetRightPart_(point, state, next_state);

  std::vector<double> eulerian_solution;
  eulerian_solver_.Step(point, state, &eulerian_solution, &point);

  std::vector<double> term;
  GetRightPart_(point, eulerian_solution, &term);

  for (unsigned i = 0; i < state_dim_; ++i) {
    next_state->operator [](i) = state[i] + 0.5 * step_ *
                                 (next_state->operator [](i) + term[i]);
  }
  if (next_point) *next_point = point;
}
