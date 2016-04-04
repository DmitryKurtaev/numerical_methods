#ifndef INCLUDE_CUBIC_SPLINE_H_
#define INCLUDE_CUBIC_SPLINE_H_

#include <vector>

#include "include/system_solver.h"

static const double kZeroLimit = 1e-10;

class CubicSpline {
 public:
  void Build(const std::vector<double>& x,
             const std::vector<double>& y,
             double lbound_xx,
             double rbound_xx);

  void Show() const;

  void PrintCoeffs() const;

  double GetValue(double x) const;

  double GetDerivate(double x) const;

  double GetSecondDerivate(double x) const;

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

#endif  // INCLUDE_CUBIC_SPLINE_H_
