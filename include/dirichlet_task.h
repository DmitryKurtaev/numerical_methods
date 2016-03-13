#ifndef INCLUDE_DIRICHLET_TASK_H_
#define INCLUDE_DIRICHLET_TASK_H_

#include <vector>

enum Border { TOP, RIGHT, BOTTOM, LEFT };

class DirichletTask {
 public:
  DirichletTask(double left, double right, double top, double bottom,
                int n_intervals_by_x, int n_intervals_by_y,
                double* external_heat);

  ~DirichletTask();

  void UpdateBorder(Border border, const DirichletTask& src);

  void UpdateBorder(Border border, double* src);

  void MPIIteration(double& achieved_eps);

  void GetState(std::vector<double>& state);

  void GetDimensions(int& n, int& m);

 private:
  double* x_;
  double* b_;
  double* external_heat_;
  double* borders_condition_[4];
  double borders_[4];
  double h_;
  double k_;
  int n_;
  int m_;
};

#endif  // INCLUDE_DIRICHLET_TASK_H_
