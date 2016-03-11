#ifndef INCLUDE_MPI_SOLVER_H_
#define INCLUDE_MPI_SOLVER_H_

class MPISolver {
 public:
  static void Iteration(double* x, double* b, int n, int m, double h, double k,
                        double& achieved_eps);

 private:
  
};

#endif