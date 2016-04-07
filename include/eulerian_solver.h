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
                 double step, unsigned state_dim,
                 void (*GetRobustValues)(const std::vector<double>& points,
                                         const std::vector<double>& init_state,
                                         std::vector<double>* states) = 0);

  void Step(double point,
            const std::vector<double>& state,
            std::vector<double>* next_state,
            double step,
            double* next_point = 0);

  unsigned GetOrder() { return 1; }
};

#endif  // INCLUDE_EULERIAN_SOLVER_H_
