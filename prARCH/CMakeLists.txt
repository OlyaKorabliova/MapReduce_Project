cmake_minimum_required(VERSION 3.6)
project(threads_project)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES callable_example.cpp)


add_executable(threads_project ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(threads_project pthread)




