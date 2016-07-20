cmake_minimum_required(VERSION 2.8)

PROJECT(raytracer)
SET(VTK_DIR /usr/include/vtk-6.0)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(raytracer MACOSX_BUNDLE raytracer)

SET(CMAKE_CXX_LINK_FLAGS)
if(VTK_LIBRARIES)
  target_link_libraries(raytracer ${VTK_LIBRARIES})
else()
  target_link_libraries(raytracer vtkHybrid)
endif()


