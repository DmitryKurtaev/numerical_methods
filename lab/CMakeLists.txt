include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../3rdparty")

set(labs
  cubic_spline_lab.cc
  heat_equation_lab.cc
  jacobi_method_lab.cc
  mpi_with_optim_step_lab.cc
  mpi_combined_areas.cc
)

foreach(lab ${labs})
  get_filename_component(target ${lab} NAME_WE)
  add_executable(${target} ${lab})
  target_link_libraries(${target} ${LIBRARY})
endforeach(lab)

