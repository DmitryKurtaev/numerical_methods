#ifndef INCLUDE_RUNGE_KUTTA_SOLVER_H_
#define INCLUDE_RUNGE_KUTTA_SOLVER_H_

#include <vector>

#include "include/abstract_solver.h"
#include "include/eulerian_solver.h"

class RungeKuttaSolver : public AbstractSolver {
 public:
  RungeKuttaSolver(void (*GetRightPart)(double point,
                                        const std::vector<double>& state,
                                        std::vector<double>* derivation),
                 double step, unsigned state_dim);

  void Step(double point,
            const std::vector<double>& state,
            std::vector<double>* next_state,
            double* next_point = 0);

  unsigned GetOrder() { return 2; }

 private:
  EulerianSolver eulerian_solver_;
};

#endif
