cmake_minimum_required(VERSION 3.14)
project(RST)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(RST
        main.cpp
        MExponent.h)
