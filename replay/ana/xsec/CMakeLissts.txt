cmake_minimum_required(VERSION 2.8)
add_definitions("-std=c++11 -Wall -O3 -march=native")
include_directories("/usr/local/include/eigen3") # eigen33??????
add_executable(main test.cc)
