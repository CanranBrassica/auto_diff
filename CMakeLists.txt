cmake_minimum_required(VERSION 3.10)
project(kyadet)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_COMPILER /usr/bin/g++)

enable_testing()

include_directories(include)

add_subdirectory(test)
