#ifndef INCLUDE_ABSTRACT_SOLVER_H_
#define INCLUDE_ABSTRACT_SOLVER_H_

#include <vector>

// Interface for differential equation solving.
class AbstractSolver {
 public:
  AbstractSolver(void (*GetRightPart)(const std::vector<double>& point,
                                      const std::vector<double>& state,
                                      std::vector<double>* derivation),
                 double step, unsigned point_dim, unsigned state_dim);

  virtual void Step(const std::vector<double>& point,
                    const std::vector<double>& state,
                    std::vector<double>* next_state,
                    std::vector<double>* next_point = 0) = 0;

 protected:
  // Function called for filling derivations vector.
  void (*GetRightPart_)(const std::vector<double>& point,
                        const std::vector<double>& state,
                        std::vector<double>* derivation);
  double step_;
  unsigned point_dim_;
  unsigned state_dim_;
};

#endif  // INCLUDE_ABSTRACT_SOLVER_H_
