#include <vector>
#include <string>

class Plot;
static Plot* current_plot;

// This class draws points sets
class Plot {
public:
  Plot();

  void Add(const std::vector<float>& x,
           const std::vector<float>& y,
           int points_size,
           float color_red,
           float color_green,
           float color_blue);

  void Show();

  ~Plot();

private:
  static const int kMarkersLength = 4;
  static const int kCharsShifts = 10;
  static const int kTopIndent = 10;
  static const int kBottomIndent = 50;
  static const int kLeftIndent = 80;
  static const int kRightIndent = 10;
  static const int kArrowsWidth = 5;
  static const int kArrowsHeight = 10;
  static const float kAxisesColor[];

  struct Set {
    std::vector<float> x;
    std::vector<float> y;
    int points_size;
    float color[3];
  };

  void DrawAxises();

  void DrawMarkers();

  void DrawPoints();

  void DrawString(std::string str, int x, int y);

  void DrawCircle(int x, int y, int radius);

  static void Display();

  static void Reshape(int width, int height);

  int view_width_;
  int view_height_;
  float min_x_;
  float max_x_;
  float min_y_;
  float max_y_;
  float y_ratio_;
  float x_ratio_;
  int first_marker_x_;
  int first_marker_y_;
  std::vector<Set> points_sets_;
};
