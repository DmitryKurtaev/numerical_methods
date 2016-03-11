#include "include/dirichlet_task.h"

DirichletTast::DirichletTask(double left, double right, double top,
                             double bottom, int n_intervals_by_x,
                             int n_intervals_by_y, double* external_heat) 
: n_(n_intervals_by_x), m_(n_intervals_by_y), ex
  h_((right - left) / n_), k_((top - bottom) / m_) {

    borders_[TOP] = top;
    borders_[RIGHT] = right;
    borders_[BOTTOM] = bottom; 
    borders_[LEFT] = left;

    borders_condition_[TOP] = new double[n_ - 1];
    borders_condition_[RIGHT] = new double[m_ - 1];
    borders_condition_[BOTTOM] = new double[n_ - 1];
    borders_condition_[LEFT] = new double[m_ - 1];

    const int dim = (n_ - 1) * (m_ - 1);
    x_ = new double[dim];
    b_ = new double[dim];
    memset(x_, 0, sizeof(double) * dim);  // Initial solution.
}

DirichletTast::~DirichletTast() {
  delete[] b_;
  delete[] x_;
  for (int i = 0; i < 4; ++i) {
    delete borders_condition_[i];
  }
}

void DirichletTast::UpdateBorder(Border border, const DirichletTask& src) {
  double* mem = 0;
  switch (border) {
    case TOP: case BOTTOM: {
      mem = new double[n_ - 1];
      const int offset = (border == BOTTOM ? 0 : m_ * (n_ - 1));
      for (int i = 0; i < n_ - 1; ++i) {
        mem[i] = src.x_[offset + i];
      }
      break;
    }
    case LEFT: case RIGHT: {
      mem = new double[m_ - 1];
      const int offset = (border == LEFT ? 0 : n_ - 2);
      for (int i = 0; i < m_ - 1; ++i) {
        mem[i] = src.x_[i * (n_ - 1) + offset];
      }
      break;
    }
    default: return;
  }
  UpdateBorder(border, mem);
  delete[] mem;
}

void DirichletTast::UpdateBorder(Border border, double* src) {
  const double inv_k_quad = 1.0 / (k_ * k_);
  const double inv_h_quad = 1.0 / (h_ * h_);
  switch (border) {
    case TOP: case BOTTOM: {
      const offset = (border == BOTTOM ? 0 : (m_ - 2) * (n_ - 1));
      for (int i = 0; i < n_ - 1; ++i) {
        b[offset + i] = external_heat_[offset + i] - inv_k_quad * src[i];
      }
      const row = (border == BOTTOM ? 0 : m_ - 2);
      b[offset] -= inv_h_quad * borders_condition_[LEFT][row];
      b[offset + n_ - 2] -= inv_h_quad * borders_condition_[RIGHT][row];
      memcpy(borders_condition_[border], src, sizeof(double) * (n_ - 1));
      break;
    }
    case LEFT: case RIGHT: {
      const int offset = (border == LEFT ? 0 : n_ - 2);
      for (int i = 0; i < m_ - 1; ++i) {
        b[i * (n_ - 1) + offset] = external_heat_[i * (n_ - 1) + offset] -
                                   inv_h_quad * src[i];
      }
      b[offset] -= inv_k_quad * borders_condition_[BOTTOM][offset];
      b[m_ * (n_ - 1) + offset] -= inv_k_quad * borders_condition_[TOP][offset];
      memcpy(borders_condition_[border], src, sizeof(double) * (m_ - 1));
      break;
    }
    default: break;
  }
}