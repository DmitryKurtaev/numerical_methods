#include "include/cubic_spline.h"
#include "include/plot.h"
#include <math.h>

void CubicSpline::Build(const std::vector<float>& x,
                        const std::vector<float>& y,
                        float lbound_xx,
                        float rbound_xx) {
  SolveSystem(x, y, lbound_xx, rbound_xx);
  c_.push_back(rbound_xx);
  ComputeCoefficients(y, lbound_xx);
  x_ = x;
  y_ = y;
}

void CubicSpline::SolveSystem(const std::vector<float>& x,
                              const std::vector<float>& y,
                              float lbound_xx,
                              float rbound_xx) {
  int n = x.size() - 1;
  h_.resize(n);
  for (int i = 0; i < n; ++i) {
    h_[i] = x[i + 1] - x[i];
  }

  std::vector<float> sub_diag(n - 2);
  std::vector<float> main_diag(n - 1);
  std::vector<float> right_part(n - 1);
  for (int i = 0; i < n - 2; ++i) {
    sub_diag[i] = h_[i + 1];
    main_diag[i] = 2 * (h_[i] + h_[i + 1]);
    right_part[i] = 6 * ((y[i] - y[i + 1]) / h_[i] +
                    (y[i + 2] - y[i + 1]) / h_[i + 1]);
  }
  main_diag[n - 2] = 2 * (h_[n - 2] + h_[n - 1]);
  right_part[n - 2] = 6 * ((y[n - 2] - y[n - 1]) / h_[n - 2] +
                      (y[n] - y[n - 1]) / h_[n - 1]) - h_[n - 1] * rbound_xx;
  right_part[0] -= h_[0] * lbound_xx;

  solver_.Solve(sub_diag,
                main_diag,
                sub_diag,
                right_part,
                &c_);
}

void CubicSpline::ComputeCoefficients(const std::vector<float>& y,
                                      float lbound_xx) {
  a_.resize(c_.size());
  b_.resize(c_.size());
  d_.resize(c_.size());
  a_[0] = y[1];
  b_[0] = (y[1] - y[0] +
          1.0 / 6 * pow(h_[0], 2) * (2 * c_[0] + lbound_xx)) / h_[0];
  d_[0] = (c_[0] - lbound_xx) / h_[0];
  for (int i = 1; i < y.size() - 1; ++i) {
    a_[i] = y[i + 1];
    b_[i] = (y[i + 1] - y[i] +
            1.0 / 6 * pow(h_[i], 2) * (2 * c_[i] + c_[i - 1])) / h_[i];
    d_[i] = (c_[i] - c_[i - 1]) / h_[i];
  }
}

void CubicSpline::Show() {
  const float spline_step = 0.001f;

  Plot plot;
  std::vector<float> spline_x;
  std::vector<float> spline_y;
  for (int i = 0; i < x_.size() - 1; ++i) {
    for (float x = x_[i]; x < x_[i + 1]; x += spline_step) {
      float y  = a_[i] +
                 b_[i] * (x - x_[i + 1]) +
                 0.5 * c_[i] * pow(x - x_[i + 1], 2) +
                 1.0 / 6 * d_[i] * pow(x - x_[i + 1], 3);
      spline_x.push_back(x);
      spline_y.push_back(y);
    }
  }
  plot.Add(spline_x, spline_y, 1, 0.0, 1.0, 0.0);
  plot.Add(x_, y_, 5, 1.0, 0.0, 0.0);
  plot.Show();
}
