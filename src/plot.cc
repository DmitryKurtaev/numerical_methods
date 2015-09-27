#include "include/plot.h"
#include <float.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

const float Plot::kAxisesColor[] = {0.0f, 1.0f, 0.0f};

Plot::Plot()
  : min_x_(FLT_MAX),
    max_x_(FLT_MIN),
    min_y_(FLT_MAX),
    max_y_(FLT_MIN),
    view_width_(500),
    view_height_(500) {
}

void Plot::Add(const std::vector<float>& x,
               const std::vector<float>& y,
               int points_size,
               float color_red,
               float color_green,
               float color_blue) {
  Set new_set;
  new_set.x = x;
  new_set.y = y;
  new_set.points_size = points_size;
  new_set.color[0] = color_red;
  new_set.color[1] = color_green;
  new_set.color[2] = color_blue;
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

void Plot::Show() {
  // Init window.
  int argc = 0;
  glutInit(&argc, 0);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
  glutInitWindowSize(view_width_, view_height_);
  glutInitWindowPosition(0, 0);
  glutCreateWindow("Plot");

  glutDisplayFunc(Display);
  glutReshapeFunc(Reshape);

  // Init GL.
  glClearColor(0.0, 0.0, 0.0, 1.0);

  current_plot = this;

  glutMainLoop();
}

void Plot::Display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  glOrtho(0, current_plot->view_width_, 0, current_plot->view_height_, 0, 1);
  current_plot->DrawAxises();
  current_plot->DrawPoints();
  glutSwapBuffers();
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
  float values_step = (max_y_ - min_y_) / (kVerticalAxisMarkers - 1);
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
}

void Plot::DrawPoints() {
  for (int i = 0; i < points_sets_.size(); ++i) {
    for (int j = 0; j < points_sets_[i].x.size(); ++j) {
      glColor3fv(points_sets_[i].color);
      int x = first_marker_x_ + (points_sets_[i].x[j] - min_x_) / x_ratio_;
      int y = first_marker_y_ + (points_sets_[i].y[j] - min_y_) / y_ratio_;
      DrawCircle(x, y, points_sets_[i].points_size);
    }
  }
}

void Plot::DrawString(std::string str, int x, int y) {
  const float kFontScaling = 0.1;

  glPushMatrix();
  for (int i = 0; i < str.length(); ++i) {
    glPushMatrix();
    glTranslatef(x + i * kCharsShifts, y, 0);
    glScalef(kFontScaling, kFontScaling, 1.0);
    glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, str[i]);
    glPopMatrix();
  }
  glPopMatrix();
}

Plot::~Plot() {
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
