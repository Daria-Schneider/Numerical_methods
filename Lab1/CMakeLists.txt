add_executable(lab1
    main.cpp
    src/matrix.cpp
    src/file_io.cpp
    src/lu_decomposition.cpp
    src/tridiagonal_sweep.cpp
    src/iterative_methods.cpp
    src/jacobi_eigen.cpp
    src/qr_decomposition.cpp
)

target_include_directories(lab1 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# путь к папке с входными данными
target_compile_definitions(lab1 PRIVATE LAB1_DATA_DIR="${PROJECT_SOURCE_DIR}/data/lab1/")
