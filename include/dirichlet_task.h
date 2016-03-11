#ifndef INCLUDE_DIRICHLET_TASK_H_
#define INCLUDE_DIRICHLET_TASK_H_

enum Border { TOP, RIGHT, BOTTOM, LEFT };

class DirichletTask {
 public:
  DirichletTask(double left, double right, double top, double bottom,
                int n_intervals_by_x, int n_intervals_by_y,
                double* external_heat);

  ~DirichletTask();

  void UpdateBorder(Border border, const DirichletTask& src);

  void UpdateBorder(Border border, double* src);

 private:
  double* x_;
  double* b_;
  double* external_heat_;
  double* borders_condition_[4];
  double borders_[4];
  const double h_;
  const double k_;
  const int n_;
  const int m_;
};

#endif  // INCLUDE_DIRICHLET_TASK_H_