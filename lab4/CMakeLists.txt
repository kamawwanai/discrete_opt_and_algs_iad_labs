cmake_minimum_required(VERSION 3.20)
project(week-1-pr LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(lab-4-pr
                src/main.cpp)

# target_include_directories(mst_prim_lab PRIVATE include)