cmake_minimum_required(VERSION 3.13)
project(test)

set(CMAKE_CXX_STANDARD 14)
file(GLOB sourcefiles "*.h" "*.cpp")
add_executable(out ${sourcefiles})