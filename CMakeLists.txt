cmake_minimum_required(VERSION 3.6)
project(RKF4_5__dormandprince)

set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES main.c)
add_executable(RKF4_5__dormandprince ${SOURCE_FILES})

target_link_libraries(RKF4_5__dormandprince m)
