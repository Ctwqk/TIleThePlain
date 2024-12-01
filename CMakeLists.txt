cmake_minimum_required(VERSION  2.8)
project(hw3)
set(SOURCES src/main.cpp ${CMAKE_SOURCE_DIR}/src/HelperFunctions.cpp)
include_directories(
        ${CMAKE_SOURCE_DIR}/src
)


find_package(OpenCV REQUIRED)
include_directories( ${OpenCV_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIRS})

find_package(glfw3 3.3 REQUIRED)
include_directories(${GLFW_INCLUDE_DIRS})

set(IMGUI_DIR ${CMAKE_SOURCE_DIR}/src/imgui)
include_directories(${IMGUI_DIR})
file(GLOB IMGUI_SOURCES ${IMGUI_DIR}/*.cpp ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp)

add_executable(main ${SOURCES} ${IMGUI_SOURCES})

target_link_Libraries(main ${OpenCV_LIBS} Eigen3::Eigen ${OPENGL_LIBRARIES} glfw dl pthread)


#add_executable(vidDis vidDis.cpp threshold.cpp denoise.cpp region.cpp features.cpp csv_util.cpp NetWorkThings/dnnEmbedding.cpp dist.cpp)
#target_link_Libraries(vidDis ${OpenCV_LIBS})
