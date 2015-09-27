#include "include/plot.h"
#include "include/cubic_spline.h"

Plot plot;

int main(int argc, char** argv) {
//  float x[] = {-1, -0.5, 0, 0.5, 1};
//  float y[] = {1, 0.375, 0, 0.625, 3};
  float x[] = {-1, 2, 3, 4.5, 7};
  float y[] = {1, 0.375, 0.55, 2.625, 3};
  std::vector<float> xs(x, x + 5);
  std::vector<float> ys(y, y + 5);
  CubicSpline cubic_spline;
  cubic_spline.Build(xs, ys, -2, 10);
  cubic_spline.Show();
  plot.Add(xs, ys, 5, 1.0, 0, 0);
  plot.Show();
  return 0;
}


