#include <vector>
#include "include/system_solver.h"

class CubicSpline {
public:
  void Build(const std::vector<float>& x,
             const std::vector<float>& y,
             float lbound_xx,
             float rbound_xx);

  void Show();
private:
  void SolveSystem(const std::vector<float>& x,
                   const std::vector<float>& y,
                   float lbound_xx,
                   float rbound_xx);

  void ComputeCoefficients(const std::vector<float>& y,
                           float lbound_xx);

  SystemSolver solver_;
  std::vector<float> a_;
  std::vector<float> b_;
  std::vector<float> c_;
  std::vector<float> d_;
  std::vector<float> h_;
  std::vector<float> x_;
  std::vector<float> y_;
};
