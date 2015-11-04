cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(PROJECT_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include/")
set(PROJECT_SOURCE_DIR "${CMAKE_SOURCE_DIR}/src/")

# Spline.
project(cubic_spline)

set(sources
  src/command_line_parser.cc
  src/cubic_spline.cc
  src/main_cubic_spline.cc
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

add_executable(${CMAKE_PROJECT_NAME} ${sources} ${headers})
target_link_libraries(${CMAKE_PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

# Heat equation.
project(heat_equation)

set(sources
  src/command_line_parser.cc
  src/main_heat_equation.cc
  src/plot.cc
  src/system_solver.cc
)

set(headers
  include/command_line_parser.h
  include/plot.h
  include/system_solver.h
)

find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})

include_directories(${CMAKE_SOURCE_DIR})

add_executable(${CMAKE_PROJECT_NAME} ${sources} ${headers})
target_link_libraries(${CMAKE_PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})
