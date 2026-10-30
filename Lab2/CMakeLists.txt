add_executable(lab2
    main.cpp
    src/utils.cpp
    src/nonlinear_equation.cpp
    src/nonlinear_system.cpp
)

target_include_directories(lab2 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# путь к папке с входными данными
target_compile_definitions(lab2 PRIVATE LAB2_DATA_DIR="${PROJECT_SOURCE_DIR}/data/lab2/")
