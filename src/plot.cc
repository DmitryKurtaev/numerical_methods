#include "include/plot.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include <string>
#include <sstream>
#include <iomanip>

#include <GL/freeglut.h>
#include <opencv2/opencv.hpp>

const float Plot::kAxisesColor[] = {0.0f, 0.0f, 0.0f};

Plot::Plot()
  : min_x_(DBL_MAX),
    max_x_(DBL_MIN),
    min_y_(DBL_MAX),
    max_y_(DBL_MIN),
    view_width_(500),
    view_height_(500) {
}

void Plot::Add(const std::vector<double>& x, const std::vector<double>& y,
               int points_size, float color_red, float color_green,
               float color_blue, bool uniform) {
  Set new_set;
  new_set.x = x;
  new_set.y = y;
  new_set.points_size = points_size;
  new_set.color[0] = color_red;
  new_set.color[1] = color_green;
  new_set.color[2] = color_blue;
  new_set.uniform = uniform;
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] < min_x_)
      min_x_ = x[i];
    if (x[i] > max_x_)
      max_x_ = x[i];
    if (y[i] < min_y_)
      min_y_ = y[i];
    if (y[i] > max_y_)
      max_y_ = y[i];
  }
  points_sets_.push_back(new_set);
}

void Plot::Show(const std::string& title, const std::string& xtitle,
                const std::string& ytitle, const std::string& saving_path) {
  xtitle_ = xtitle;
  ytitle_ = ytitle;
  saving_path_ = saving_path;

  // Init window.
  int argc = 0;
  glutInit(&argc, 0);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
  glutInitWindowSize(view_width_, view_height_);
  glutInitWindowPosition(0, 0);
  window_handle_ = glutCreateWindow(title.c_str());
  glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);

  glutDisplayFunc(Display);
  glutReshapeFunc(Reshape);
  glutKeyboardFunc(KeyPressed);

  // Init GL.
  glClearColor(1.0, 1.0, 1.0, 1.0);

  current_plot = this;

  if (saving_path_ != "") {
    glutHideWindow();
  }
  glutMainLoop();
}

void Plot::Display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  glOrtho(0, current_plot->view_width_, 0, current_plot->view_height_, 0, 1);
  current_plot->DrawAxises();
  current_plot->DrawPoints();
  glutSwapBuffers();

  if (current_plot->saving_path_ != "") {
    const int width = current_plot->view_width_;
    const int height = current_plot->view_height_;
    unsigned char* data = new unsigned char[width * height * 3];

    glReadPixels(0, 0, width, height, GL_BGR, GL_UNSIGNED_BYTE, data);
    cv::Mat mat(height, width, CV_8UC3, data);
    cv::flip(mat, mat, 0);
    cv::imwrite(current_plot->saving_path_, mat);
    current_plot->saving_path_ = "";
    KeyPressed(27, 0, 0);
  }
}

void Plot::Reshape(int width, int height) {
  current_plot->view_width_ = width;
  current_plot->view_height_ = height;
  glViewport(0, 0, width, height);
}

void Plot::DrawAxises() {
  glBegin(GL_LINES);

  glColor3fv(kAxisesColor);
  // Verticel axis.
  glVertex2i(kLeftIndent, kBottomIndent);
  glVertex2i(kLeftIndent, view_height_ - 1 - kTopIndent);

  // Verticel axis arrow.
  glVertex2i(kLeftIndent, view_height_ - 1 - kTopIndent);
  glVertex2i(kLeftIndent - kArrowsWidth,
             view_height_ - 1 - kTopIndent - kArrowsHeight);

  glVertex2i(kLeftIndent, view_height_ - 1 - kTopIndent);
  glVertex2i(kLeftIndent + kArrowsWidth,
             view_height_ - 1 - kTopIndent - kArrowsHeight);

  // Horizontal axis.
  glVertex2i(kLeftIndent, kBottomIndent);
  glVertex2i(view_width_ - 1 - kRightIndent, kBottomIndent);

  // Horizontal axis arrow.
  glVertex2i(view_width_ - 1 - kRightIndent, kBottomIndent);
  glVertex2i(view_width_ - 1 - kRightIndent - kArrowsHeight,
             kBottomIndent - kArrowsWidth);

  glVertex2i(view_width_ - 1 - kRightIndent, kBottomIndent);
  glVertex2i(view_width_ - 1 - kRightIndent - kArrowsHeight,
             kBottomIndent + kArrowsWidth);

  glEnd();
  DrawMarkers();
}

void Plot::DrawMarkers() {
  const int kMarkersSize = 5;
  const int kVerticalAxisMarkers = 5;
  const int kHorizontalAxisMarkers = 5;
  const int kGridSize = 3;
  const int kGridStep = 3;
  const float kGridColor[] = {0.5f, 0.5f, 0.5f};

  // Markers of vertical axis
  int markers_step = (view_height_ - 1 - kTopIndent - kArrowsHeight -
                      kBottomIndent) / kVerticalAxisMarkers;
  double values_step = (max_y_ - min_y_) / (kVerticalAxisMarkers - 1);
  y_ratio_ = values_step / markers_step;
  first_marker_y_ = kBottomIndent + markers_step / 2;
  for (int i = 0; i < kVerticalAxisMarkers; ++i) {
    int y = first_marker_y_ + i * markers_step;

    glBegin(GL_LINES);
    glColor3fv(kAxisesColor);
    glVertex2i(kLeftIndent - kMarkersSize, y);
    glVertex2i(kLeftIndent + kMarkersSize, y);
    glColor3fv(kGridColor);
    for (int x = kLeftIndent + kMarkersSize + kGridStep;
         x < view_width_ - kRightIndent;
         x += kGridSize + kGridStep) {
      glVertex2i(x, y);
      glVertex2i(x + kGridSize, y);
    }
    glEnd();

    glColor3fv(kAxisesColor);
    std::ostringstream oss;
    oss << std::setprecision(kMarkersLength) << min_y_ + values_step * i;
    std::string str(oss.str());
    int x = kLeftIndent - kMarkersSize - kCharsShifts * (str.length() + 1);
    DrawString(str, x, y);
  }

  // Markers of horizontal axis
  markers_step = (view_width_ - 1 - kRightIndent - kArrowsHeight -
                  kLeftIndent) / kVerticalAxisMarkers;
  values_step = (max_x_ - min_x_) / (kVerticalAxisMarkers - 1);
  x_ratio_ = values_step / markers_step;
  first_marker_x_ = kLeftIndent + markers_step / 2;
  for (int i = 0; i < kHorizontalAxisMarkers; ++i) {
    int x = first_marker_x_ + i * markers_step;

    glBegin(GL_LINES);
    glColor3fv(kAxisesColor);
    glVertex2i(x, kBottomIndent - kMarkersSize);
    glVertex2i(x, kBottomIndent + kMarkersSize);
    glColor3fv(kGridColor);
    for (int y = kBottomIndent + kMarkersSize + kGridStep;
         y < view_height_ - kTopIndent;
         y += kGridSize + kGridStep) {
      glVertex2i(x, y);
      glVertex2i(x, y + kGridSize);
    }
    glEnd();

    glColor3fv(kAxisesColor);
    std::ostringstream oss;
    oss << std::setprecision(kMarkersLength) << min_x_ + values_step * i;
    std::string str(oss.str());
    DrawString(str,
               x - 0.5 * kCharsShifts * (str.length() + 1), kBottomIndent / 2);
  }

  // Title of axises.
  glColor3fv(kAxisesColor);
  DrawString(xtitle_, 0.5 * (view_width_ - kCharsShifts * (xtitle_.length() + 1)),
             kBottomIndent / 4);
  DrawString(ytitle_, kLeftIndent / 4,
             0.5 * (view_height_ - kCharsShifts * (ytitle_.length() + 1)),
             true);
}

void Plot::DrawPoints() {
  for (int i = 0; i < points_sets_.size(); ++i) {
    glColor3fv(points_sets_[i].color);
    int prev_x;
    int prev_y;
    for (int j = 0; j < points_sets_[i].x.size(); ++j) {
      int x = first_marker_x_ + (points_sets_[i].x[j] - min_x_) / x_ratio_;
      int y = first_marker_y_ + (points_sets_[i].y[j] - min_y_) / y_ratio_;
      if (points_sets_[i].uniform) {
        if (j != 0) {
          glBegin(GL_LINES);
          glVertex2i(prev_x, prev_y);
          glVertex2i(x, y);
          glEnd();
        }
      } else {
        DrawCircle(x, y, points_sets_[i].points_size);
      }
      prev_x = x;
      prev_y = y;
    }
  }
}

void Plot::DrawString(std::string str, int x, int y, bool vertical) {
  const double kFontScaling = 0.1;

  glPushMatrix();
  for (int i = 0; i < str.length(); ++i) {
    glPushMatrix();
    if (vertical) {
      glTranslatef(x, y + i * kCharsShifts, 0);
      glRotatef(90, 0, 0, 1);
    } else {
      glTranslatef(x + i * kCharsShifts, y, 0);
    }
    glScalef(kFontScaling, kFontScaling, 1.0);
    glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, str[i]);
    glPopMatrix();
  }
  glPopMatrix();
}

Plot::~Plot() {
  Clear();
}

void Plot::Clear() {
  min_x_ = DBL_MAX;
  max_x_ = DBL_MIN;
  min_y_ = DBL_MAX;
  max_y_ = DBL_MIN;
  view_width_ = 500;
  view_height_ = 500;
  for (int i = 0; i < points_sets_.size(); ++i) {
    points_sets_[i].x.clear();
    points_sets_[i].y.clear();
  }
  points_sets_.clear();
}

void Plot::DrawCircle(int x, int y, int radius) {
  glBegin(GL_POINTS);
  int cur_x = x + radius;
  int cur_y = y;
  while (cur_x != x) {
    for (int i = 2 * x - cur_x; i <= cur_x; ++i) {
      glVertex2i(i, cur_y);
      glVertex2i(i, 2 * y - cur_y);
    }
    int weight_top = abs(pow(cur_x - x, 2) +
                         pow(cur_y - y + 1, 2) - pow(radius, 2));
    int weight_diag = abs(pow(cur_x - x - 1, 2) +
                          pow(cur_y - y + 1, 2) - pow(radius, 2));
    int weight_left = abs(pow(cur_x - x - 1, 2) +
                          pow(cur_y - y, 2) - pow(radius, 2));
    if (weight_top < weight_diag &&
        weight_top < weight_left) {
      ++cur_y;
    } else {
      if (weight_diag < weight_top &&
          weight_diag < weight_left) {
        ++cur_y;
        --cur_x;
      } else {
        --cur_x;
      }
    }
  }
  glVertex2i(x, y + radius);
  glVertex2i(x, y - radius);
  glEnd();
}

void Plot::KeyPressed(unsigned char key, int x, int y) {
  if (key == 27) {
    glutDestroyWindow(current_plot->window_handle_);
  }
}
