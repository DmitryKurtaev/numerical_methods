#ifndef INCLUDE_EULERIAN_SOLVER_H_
#define INCLUDE_EULERIAN_SOLVER_H_

#include <vector>

#include "include/abstract_solver.h"

// Interface for differential equation solving.
class EulerianSolver : public AbstractSolver {
 public:
  EulerianSolver(void (*GetRightPart)(const std::vector<double>& point,
                                      const std::vector<double>& state,
                                      std::vector<double>* derivation),
                 double step, unsigned point_dim, unsigned state_dim);

  void Step(const std::vector<double>& point,
            const std::vector<double>& state,
            std::vector<double>* next_state,
            std::vector<double>* next_point = 0);
};

#endif  // INCLUDE_EULERIAN_SOLVER_H_
