#include "include/system_solver.h"
#include <stdlib.h>

void SystemSolver::Solve(const std::vector<double>& upper_diag,
                         const std::vector<double>& main_diag,
                         const std::vector<double>& lower_diag,
                         const std::vector<double>& right_part,
                         std::vector<double>* y) {
  y->clear();
  if (main_diag.size() < 2 ||
      upper_diag.size() != main_diag.size()- 1 ||
      lower_diag.size() != main_diag.size() - 1 ||
      right_part.size() != main_diag.size())
    return;
  if (main_diag[0] == 0 ||
      main_diag[main_diag.size() - 1] == 0)
    return;

  kappa_1_ = -upper_diag[0] / main_diag[0];
  kappa_2_ = -lower_diag[lower_diag.size() - 1] /
             main_diag[main_diag.size()  - 1];
  if (abs(kappa_1_) > 1 ||
      abs(kappa_2_) >= 1)
    return;

  mu_1_ = right_part[0] / main_diag[0];
  mu_2_ = right_part[right_part.size() - 1] / main_diag[main_diag.size() - 1];

  n_ = main_diag.size() - 1;
  A_.resize(n_ - 1);
  B_.resize(n_ - 1);
  C_.resize(n_ - 1);
  phi_.resize(n_ - 1);
  for (int i = 0; i < n_ - 1; ++i) {
    A_[i] = lower_diag[i];
    B_[i] = upper_diag[i + 1];
    C_[i] = -main_diag[i + 1];
    phi_[i] = -right_part[i + 1];
    if (A_[i] == 0 ||
        B_[i] == 0 ||
        abs(C_[i]) < abs(A_[i]) + abs(B_[i]))
      return;
  }

  ForwardMove();
  BackwardMove(y);

  phi_.clear();
  A_.clear();
  B_.clear();
  C_.clear();
  alpha_.clear();
  beta_.clear();
}

void SystemSolver::ForwardMove() {
  alpha_.resize(n_);
  beta_.resize(n_);
  alpha_[0] = kappa_1_;
  beta_[0] = mu_1_;
  for (int i = 1; i < n_; ++i) {
    double denominator = C_[i - 1] - alpha_[i - 1] * A_[i - 1];
    alpha_[i] = B_[i - 1] / denominator;
    beta_[i] = (phi_[i - 1] + beta_[i - 1] * A_[i - 1]) / denominator;
  }
}

void SystemSolver::BackwardMove(std::vector<double>* y) {
  y->resize(n_ + 1);
  (*y)[n_] = (kappa_2_ * beta_[n_ - 1] + mu_2_) /
      (1 - kappa_2_ * alpha_[n_ - 1]);
  for (int i = n_ - 1; i >= 0; --i) {
    (*y)[i] = alpha_[i] * (*y)[i + 1] + beta_[i];
  }
}
