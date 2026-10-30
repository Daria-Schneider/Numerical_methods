cmake_minimum_required(VERSION 3.10)

project(Numerical_Methods_MAI LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra)
endif()

# все exe складываем в build/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(Lab1)
add_subdirectory(Lab2)
add_subdirectory(Lab3)
add_subdirectory(Lab4)
