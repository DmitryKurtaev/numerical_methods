#include "include/abstract_solver.h"

AbstractSolver::AbstractSolver(
    void (*GetRightPart)(double point,
                         const std::vector<double>& state,
                         std::vector<double>* derivation),
    double step, unsigned state_dim)
  : GetRightPart_(GetRightPart), step_(step), state_dim_(state_dim) {
}
