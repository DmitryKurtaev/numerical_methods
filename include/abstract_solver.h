#ifndef INCLUDE_ABSTRACT_SOLVER_H_
#define INCLUDE_ABSTRACT_SOLVER_H_

#include <vector>

// Interface for differential equation solving.
class AbstractSolver {
 public:
  AbstractSolver(void (*GetRightPart)(double point,
                                      const std::vector<double>& state,
                                      std::vector<double>* derivation),
                 double step, unsigned state_dim);

  virtual void Step(double point,
                    const std::vector<double>& state,
                    std::vector<double>* next_state,
                    double* next_point = 0) = 0;

  virtual unsigned GetOrder() = 0;

  void StepWithLocalErrorControl(double point,
                                 const std::vector<double>& state,
                                 std::vector<double>* next_state,
                                 double eps,
                                 double* local_error,
                                 unsigned* step_increasing_counter,
                                 unsigned* step_decreasing_counter,
                                 double* next_point = 0);

 protected:
  // Function called for filling derivations vector.
  void (*GetRightPart_)(double point,
                        const std::vector<double>& state,
                        std::vector<double>* derivation);
  double step_;
  unsigned state_dim_;
};

#endif  // INCLUDE_ABSTRACT_SOLVER_H_
