#include "include/abstract_solver.h"

AbstractSolver::AbstractSolver(
    void (*GetRightPart)(const std::vector<double>& point,
                         const std::vector<double>& state,
                         std::vector<double>* derivation),
    double step, unsigned point_dim, unsigned state_dim)
  : GetRightPart_(GetRightPart), step_(step), point_dim_(point_dim),
    state_dim_(state_dim) {
}
