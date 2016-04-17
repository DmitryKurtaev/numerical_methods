cmake_minimum_required(VERSION 2.8)

set(LIBRARY numerical_methods)
project(${LIBRARY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(sources
  src/abstract_solver.cc
  src/command_line_parser.cc
  src/cubic_spline.cc
  src/dirichlet_task.cc
  src/eulerian_solver.cc
  src/mpi_solver.cc
  src/plot.cc
  src/runge_kutta_solver.cc
  src/system_solver.cc
  src/table_printer.cc
)

set(headers	
  include/abstract_solver.h
  include/command_line_parser.h
  include/cubic_spline.h
  include/dirichlet_task.h
  include/eulerian_solver.h
  include/mpi_solver.h
  include/plot.h
  include/runge_kutta_solver.h
  include/system_solver.h
  include/table_printer.h
)

find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})

find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

include_directories(${CMAKE_SOURCE_DIR})

add_subdirectory(lab)

add_library(${CMAKE_PROJECT_NAME} STATIC ${sources} ${headers})
target_link_libraries(${CMAKE_PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${OpenCV_LIBRARIES})

