cmake_minimum_required(VERSION 3.9)
project(syssoft5)

set(CMAKE_C_STANDARD 11)
include_directories(include)

set(CMAKE_BUILD_TYPE Release)
add_executable(buddhaBrot
        include/qdbmp.h
        src/main.c
        src/qdbmp.c)

target_link_libraries (buddhaBrot pthread)
