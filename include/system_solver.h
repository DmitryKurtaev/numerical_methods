#include <vector>

// This class solves linear system
// with three-diagonal coefficients matrix
class SystemSolver {
public:
  // Return empty vector if method have error.
  void Solve(const std::vector<float>& upper_diag,
             const std::vector<float>& main_diag,
             const std::vector<float>& lower_diag,
             const std::vector<float>& right_part,
             std::vector<float>* y);
private:
  void ForwardMove();

  void BackwardMove(std::vector<float>* y);

  int n_;
  float kappa_1_;
  float kappa_2_;
  float mu_1_;
  float mu_2_;
  std::vector<float> phi_;
  std::vector<float> A_;
  std::vector<float> B_;
  std::vector<float> C_;
  std::vector<float> alpha_;
  std::vector<float> beta_;
};
