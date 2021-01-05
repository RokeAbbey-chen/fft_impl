cmake_minimum_required(VERSION 3.17)
project(fft_impl)

set(CMAKE_CXX_STANDARD 14)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
add_executable(fft_impl main.cpp base/complex.cpp base/complex.hpp fft/fft.cpp fft/fft.hpp base/image_util.cpp base/image_util.hpp boot/boot.cpp boot/boot.hpp fft/fft2.cpp)
target_link_libraries(fft_impl ${OpenCV_LIBS})
