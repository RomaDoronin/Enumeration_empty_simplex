cmake_minimum_required(VERSION 3.0)

project(Enumeration_empty_simplex CXX)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)

set(SRC
Src/main.cpp)

add_executable(Enumeration_empty_simplex  ${SRC})