// Copyright 2015 Dmitry Kurtaev

#ifndef INCLUDE_PLOT_H_
#define INCLUDE_PLOT_H_

#include <vector>
#include <string>

class Plot;
static Plot* current_plot;

// This class draws points sets
class Plot {
 public:
  Plot();

  void Clear();

  void Add(const std::vector<double>& x,
           const std::vector<double>& y,
           int points_size,
           float color_red,
           float color_green,
           float color_blue,
           bool uniform);

  void Show(const std::string& title);

  ~Plot();

 private:
  static const int kMarkersLength = 4;
  static const int kCharsShifts = 10;
  static const int kTopIndent = 10;
  static const int kBottomIndent = 50;
  static const int kLeftIndent = 120;
  static const int kRightIndent = 10;
  static const int kArrowsWidth = 5;
  static const int kArrowsHeight = 10;
  static const float kAxisesColor[];

  struct Set {
    std::vector<double> x;
    std::vector<double> y;
    int points_size;
    float color[3];
    bool uniform;
  };

  void DrawAxises();

  void DrawMarkers();

  void DrawPoints();

  void DrawString(std::string str, int x, int y);

  void DrawCircle(int x, int y, int radius);

  static void Display();

  static void Reshape(int width, int height);

  static void KeyPressed(unsigned char key, int x, int y);

  int view_width_;
  int view_height_;
  double min_x_;
  double max_x_;
  double min_y_;
  double max_y_;
  double y_ratio_;
  double x_ratio_;
  int first_marker_x_;
  int first_marker_y_;
  std::vector<Set> points_sets_;
  int window_handle_;
};

#endif  // INCLUDE_PLOT_H_
