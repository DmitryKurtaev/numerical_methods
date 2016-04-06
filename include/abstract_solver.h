#ifndef INCLUDE_ABSTRACT_SOLVER_H_
#define INCLUDE_ABSTRACT_SOLVER_H_

#include <vector>

// Interface for differential equation solving.
class AbstractSolver {
 public:
  AbstractSolver(void (*GetRightPart)(double point,
                                      const std::vector<double>& state,
                                      std::vector<double>* derivation),
                 double step, unsigned state_dim,
                 void (*GetRobustValues)(const std::vector<double>& points,
                                         const std::vector<double>& init_state,
                                         std::vector<double>* states) = 0);

  virtual unsigned GetOrder() = 0;

  void Solve(const std::vector<double>& init_state, double init_point,
             double right_border, unsigned max_n_iters);

  void SolveWithLocalErrorControl(double eps);

  void ShowResults();

 protected:
  virtual void Step(double point,
                    const std::vector<double>& state,
                    std::vector<double>* next_state,
                    double* next_point = 0) = 0;

  // Function called for filling derivations vector.
  void (*GetRightPart_)(double point,
                        const std::vector<double>& state,
                        std::vector<double>* derivation);

  void (*GetRobustValues_)(const std::vector<double>& points,
                           const std::vector<double>& init_state,
                           std::vector<double>* states);
  double step_;
  unsigned state_dim_;
  std::vector<double> points_;
  std::vector<std::vector<double> > states_;

 private:
  void StepWithLocalErrorControl(double eps);

  void Reset();

  std::vector<std::vector<double> > states_double_half_step_;
  std::vector<double> local_errors_;
};

#endif  // INCLUDE_ABSTRACT_SOLVER_H_
