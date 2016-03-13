#include "include/dirichlet_task.h"

#include <stdlib.h>
#include <string.h>

#include <iostream>

#include "include/mpi_solver.h"

DirichletTask::DirichletTask(double left, double right, double top,
                             double bottom, int n_intervals_by_x,
                             int n_intervals_by_y, double* external_heat) {
    n_ = n_intervals_by_x;
    m_ = n_intervals_by_y;
    h_ = (right - left) / n_;
    k_ = (top - bottom) / m_;

    borders_[TOP] = top;
    borders_[LEFT] = left;
    borders_[RIGHT] = right;
    borders_[BOTTOM] = bottom;

    borders_condition_[TOP] = new double[n_ - 1];
    borders_condition_[LEFT] = new double[m_ - 1];
    borders_condition_[RIGHT] = new double[m_ - 1];
    borders_condition_[BOTTOM] = new double[n_ - 1];
    memset(borders_condition_[TOP], 0, sizeof(double) * (n_ - 1));
    memset(borders_condition_[LEFT], 0, sizeof(double) * (m_ - 1));
    memset(borders_condition_[RIGHT], 0, sizeof(double) * (m_ - 1));
    memset(borders_condition_[BOTTOM], 0, sizeof(double) * (n_ - 1));

    const int dim = (n_ - 1) * (m_ - 1);
    x_ = new double[dim];
    b_ = new double[dim];
    memset(x_, 0, sizeof(double) * dim);  // Initial solution.
    
    external_heat_ = new double[dim];
    memcpy(external_heat_, external_heat, sizeof(double) * dim);
    memcpy(b_, external_heat, sizeof(double) * dim);
}

DirichletTask::~DirichletTask() {
  delete[] external_heat_;
  delete[] b_;
  delete[] x_;
  for (int i = 0; i < 4; ++i) {
    delete[] borders_condition_[i];
  }
}

void DirichletTask::UpdateBorder(Border border, const DirichletTask& src) {
  double* mem = 0;
  switch (border) {
    case TOP: case BOTTOM: {
      if (n_ != src.n_) {
        std::cout << "Warning: UpdateBorder (" << n_ << " != " << src.n_ << ")"
                  << std::endl;
      }
      mem = new double[n_ - 1];
      const int offset = (border == BOTTOM ? (src.m_ - 2) * (src.n_ - 1) : 0);
      for (int i = 0; i < n_ - 1; ++i) {
        mem[i] = src.x_[offset + i];
      }
      break;
    }
    case LEFT: case RIGHT: {
      if (m_ != src.m_) {
        std::cout << "Warning: UpdateBorder (" << m_ << " != " << src.m_ << ")"
                  << std::endl;
      }
      mem = new double[m_ - 1];
      const int offset = (border == LEFT ? src.n_ - 2 : 0);
      for (int i = 0; i < m_ - 1; ++i) {
        mem[i] = src.x_[i * (src.n_ - 1) + offset];
      }
      break;
    }
    default: return;
  }
  UpdateBorder(border, mem);
  delete[] mem;
}

void DirichletTask::UpdateBorder(Border border, double* src) {
  const double inv_k_quad = 1.0 / (k_ * k_);
  const double inv_h_quad = 1.0 / (h_ * h_);
  switch (border) {
    case TOP: case BOTTOM: {
      const int row = (border == BOTTOM ? 0 : m_ - 2);
      const int offset = row * (n_ - 1);
      for (int i = 0; i < n_ - 1; ++i) {
        b_[offset + i] = external_heat_[offset + i] - inv_k_quad * src[i];
      }
      b_[offset] -= inv_h_quad * borders_condition_[LEFT][row];
      b_[offset + n_ - 2] -= inv_h_quad * borders_condition_[RIGHT][row];
      memcpy(borders_condition_[border], src, sizeof(double) * (n_ - 1));
      break;
    }
    case LEFT: case RIGHT: {
      const int offset = (border == LEFT ? 0 : n_ - 2);
      for (int i = 0; i < m_ - 1; ++i) {
        b_[i * (n_ - 1) + offset] = external_heat_[i * (n_ - 1) + offset] -
                                    inv_h_quad * src[i];
      }
      b_[offset] -= inv_k_quad * borders_condition_[BOTTOM][offset];
      b_[(m_ - 2) * (n_ - 1) + offset] -= inv_k_quad *
                                          borders_condition_[TOP][offset];
      memcpy(borders_condition_[border], src, sizeof(double) * (m_ - 1));
      break;
    }
    default: break;
  }
}

void DirichletTask::MPIIteration(double& achieved_eps) {
  MPISolver::Iteration(x_, b_, n_, m_, h_, k_, achieved_eps);
}

void DirichletTask::GetState(std::vector<double>& state) {
  const int dim = (n_ - 1) * (m_ - 1);
  state.resize(dim);
  for (int i = 0; i < dim; ++i) {
    state[i] = x_[i];
  }
}
