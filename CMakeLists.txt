cmake_minimum_required(VERSION 3.19)
project(try)

set(CMAKE_CXX_STANDARD 14)

add_executable(try main.cpp linearalgebra.h inout.cpp inout.h keyFunction.h keyFunction.cpp linearalgebra.cpp)