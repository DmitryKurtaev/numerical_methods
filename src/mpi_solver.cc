#include "include/mpi_solver.h"

#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <cmath>

void MPISolver::Iteration(double* x, double* b, int n, int m, double h,
                          double k, double& achieved_eps) {
  const int dim = (n - 1) * (m - 1);
  const double inv_h_quad = 1.0 / (h * h);
  const double inv_k_quad = 1.0 / (k * k);
  const double step = -0.5 / (inv_h_quad + inv_k_quad);

  double* new_x = new double[dim];
  memcpy(new_x, x, sizeof(double) * dim);

  // Do iteration.
  // |-- Bottom border--------------------------------------------------------
  //     |-- Left bottom point.
  int idx = 0;
  double term = inv_k_quad * x[idx + n - 1] +  // Top neighbor.
                inv_h_quad * x[idx + 1] -  // Right neighbor.
                2 * x[idx] * (inv_k_quad + inv_h_quad);
  new_x[idx] += step * (b[idx] - term);

  //     |-- Bottom line, center points.
  for (int i = 1; i < n - 2; ++i) {
    idx = i;
    term = inv_k_quad * x[idx + n - 1] +  // Top neighbor.
           inv_h_quad * x[idx + 1] +  // Right neighbor.
           inv_h_quad * x[idx - 1] -  // Left neighbor.
           2 * x[idx] * (inv_k_quad + inv_h_quad);
    new_x[idx] += step * (b[idx] - term);
  }

  //     |-- Right bottom point.
  idx = n - 2;
  term = inv_k_quad * x[idx + n - 1] +  // Top neighbor.
         inv_h_quad * x[idx - 1] -  // Left neighbor.
         2 * x[idx] * (inv_k_quad + inv_h_quad);
  new_x[idx] += step * (b[idx] - term);

  // |-- Center lines---------------------------------------------------------
  for (int j = 1; j < m - 2; ++j) {
    //   |-- Left border.
    idx = j * (n - 1);
    term = inv_k_quad * x[idx + n - 1] +  // Top neighbor.
           inv_k_quad * x[idx - n + 1] +  // Bottom neighbor.
           inv_h_quad * x[idx + 1] -  // Right neighbor.
           2 * (inv_h_quad + inv_k_quad) * x[idx];
    new_x[idx] += step * (b[idx] - term);

    //   |-- Centers.
    for (int i = 1; i < n - 2; ++i) {
      idx = j * (n - 1) + i;
      term = inv_k_quad * x[idx + n - 1] +  // Top neighbor.
             inv_k_quad * x[idx - n + 1] +  // Bottom neighbor.
             inv_h_quad * x[idx + 1] +  // Right neighbor.
             inv_h_quad * x[idx - 1] -  // Left neighbor.
             2 * (inv_h_quad + inv_k_quad) * x[idx];
      new_x[idx] += step * (b[idx] - term);
    }

    //   |-- Right border.
    idx = j * (n - 1) + n - 2;
    term = inv_k_quad * x[idx + n - 1] +  // Top neighbor.
           inv_k_quad * x[idx - n + 1] +  // Bottom neighbor.
           inv_h_quad * x[idx - 1] -  // Left neighbor.
           2 * (inv_h_quad + inv_k_quad) * x[idx];
    new_x[idx] += step * (b[idx] - term);
  }
  
  // |-- Top border-----------------------------------------------------------
  //     |-- Left top point.
  idx = (m - 2) * (n - 1);
  term = inv_k_quad * x[idx - n + 1] +  // Bottom neighbor.
         inv_h_quad * x[idx + 1] -  // Right neighbor.
         2 * (inv_h_quad + inv_k_quad) * x[idx];
  new_x[idx] += step * (b[idx] - term);

  //   |-- Centers.
  for (int i = 1; i < n - 2; ++i) {
    idx = (m - 2) * (n - 1) + i;
    term = inv_k_quad * x[idx - n + 1] +  // Bottom neighbor.
           inv_h_quad * x[idx + 1] +  // Right neighbor.
           inv_h_quad * x[idx - 1] -  // Left neighbor.
           2 * (inv_h_quad + inv_k_quad) * x[idx];
    new_x[idx] += step * (b[idx] - term);
  }

  //   |-- Right top point.
  idx = (m - 1) * (n - 1) - 1;
  term = inv_k_quad * x[idx - n + 1] +  // Bottom neighbor.
         inv_h_quad * x[idx - 1] -  // Left neighbor.
         2 * (inv_h_quad + inv_k_quad) * x[idx];
  new_x[idx] += step * (b[idx] - term);

  // Compute accuracy.
  achieved_eps = 0;
  for (int i = 0; i < dim; ++i) {
    achieved_eps = std::max(achieved_eps, fabs(x[i] - new_x[i]));
  }

  memcpy(x, new_x, sizeof(double) * dim);
  delete[] new_x;
}