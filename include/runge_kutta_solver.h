#ifndef INCLUDE_RUNGE_KUTTA_SOLVER_H_
#define INCLUDE_RUNGE_KUTTA_SOLVER_H_

#include <vector>

#include "include/abstract_solver.h"
#include "include/eulerian_solver.h"

class RungeKuttaSolver : public AbstractSolver {
 public:
  RungeKuttaSolver(void (*GetRightPart)(const std::vector<double>& point,
                                        const std::vector<double>& state,
                                        std::vector<double>* derivation),
                 double step, unsigned point_dim, unsigned state_dim);

  void Step(const std::vector<double>& point,
            const std::vector<double>& state,
            std::vector<double>* next_state,
            std::vector<double>* next_point = 0);

 private:
  EulerianSolver eulerian_solver_;
};

#endif
