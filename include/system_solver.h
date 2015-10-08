#include <vector>

// This class solves linear system
// with three-diagonal coefficients matrix
class SystemSolver {
public:
  // Return empty vector if method have error.
  void Solve(const std::vector<double>& upper_diag,
             const std::vector<double>& main_diag,
             const std::vector<double>& lower_diag,
             const std::vector<double>& right_part,
             std::vector<double>* y);
private:
  void ForwardMove();

  void BackwardMove(std::vector<double>* y);

  int n_;
  double kappa_1_;
  double kappa_2_;
  double mu_1_;
  double mu_2_;
  std::vector<double> phi_;
  std::vector<double> A_;
  std::vector<double> B_;
  std::vector<double> C_;
  std::vector<double> alpha_;
  std::vector<double> beta_;
};
