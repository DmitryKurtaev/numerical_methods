cmake_minimum_required(VERSION 2.8)

set(LIBRARY numerical_methods)
project(${LIBRARY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(sources
  src/command_line_parser.cc
  src/cubic_spline.cc
  src/plot.cc
  src/system_solver.cc
)

set(headers	
  include/command_line_parser.h
  include/cubic_spline.h
  include/plot.h
  include/system_solver.h
)

find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})

find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})

include_directories(${CMAKE_SOURCE_DIR})

add_subdirectory(lab)

add_library(${CMAKE_PROJECT_NAME} STATIC ${sources} ${headers})
target_link_libraries(${CMAKE_PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

