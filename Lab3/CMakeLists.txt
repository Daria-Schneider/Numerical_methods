add_executable(lab3
    main.cpp
    src/utils.cpp
    src/interpolation.cpp
    src/spline.cpp
    src/least_squares.cpp
    src/differentiation.cpp
    src/integration.cpp
)

target_include_directories(lab3 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# путь к папке с входными данными
target_compile_definitions(lab3 PRIVATE LAB3_DATA_DIR="${PROJECT_SOURCE_DIR}/data/lab3/")
