#ifndef INCLUDE_EULERIAN_SOLVER_H_
#define INCLUDE_EULERIAN_SOLVER_H_

#include <vector>

#include "include/abstract_solver.h"

// Interface for differential equation solving.
class EulerianSolver : public AbstractSolver {
 public:
  EulerianSolver(void (*GetRightPart)(double point,
                                      const std::vector<double>& state,
                                      std::vector<double>* derivation),
                 double step, unsigned state_dim);

  void Step(double point,
            const std::vector<double>& state,
            std::vector<double>* next_state,
            double* next_point = 0);
};

#endif  // INCLUDE_EULERIAN_SOLVER_H_
