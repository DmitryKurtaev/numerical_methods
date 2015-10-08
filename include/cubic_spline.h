#include <vector>
#include <iostream>
#include "include/system_solver.h"

static const double kZeroLimit = 1e-10;

class CubicSpline {
public:
  void Build(const std::vector<double>& x,
             const std::vector<double>& y,
             double lbound_xx,
             double rbound_xx);

  void Show();

  void PrintCoeffs(std::ostream* s);

private:
  void SolveSystem(const std::vector<double>& x,
                   const std::vector<double>& y,
                   double lbound_xx,
                   double rbound_xx);

  void ComputeCoefficients(const std::vector<double>& y,
                           double lbound_xx);

  SystemSolver solver_;
  std::vector<double> a_;
  std::vector<double> b_;
  std::vector<double> c_;
  std::vector<double> d_;
  std::vector<double> h_;
  std::vector<double> x_;
  std::vector<double> y_;
};
