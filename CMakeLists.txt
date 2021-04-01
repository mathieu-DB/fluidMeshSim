cmake_minimum_required(VERSION 3.1)
project(fluidMeshSim)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    ON)

find_package(LIBIGL REQUIRED QUIET)

# Add your project files
set(
        SRCFILES
        fluidMeshSim.cpp
        fluidMeshSim.h
        trimesh.h
        trimesh.cpp
        trimesh_types.h
        fluid.h
        fluid.cpp 
        "Source.h"
        "tools.h")

add_executable(${PROJECT_NAME} ${SRCFILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl_glfw igl::triangle igl::opengl_glfw_imgui)
