cmake_minimum_required(VERSION 2.8)

PROJECT(Isolines)
SET(VTK_DIR C:/VTKsrc)
SET(CMAKE_VERBOSE_MAKEFILE ON)
find_package(VTK REQUIRED)
include(${VTK_USE_FILE} ${VTK_DIR}/Rendering)

add_executable(Isolines Isolines)


target_link_libraries(Isolines glu32)
target_link_libraries(Isolines opengl32)
if(VTK_LIBRARIES)
target_link_libraries(Isolines ${VTK_LIBRARIES})
else()
target_link_libraries(Isolines vtkHybrid)
endif()