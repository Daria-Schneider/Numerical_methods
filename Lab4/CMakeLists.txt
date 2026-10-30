add_executable(lab4
    main.cpp
    src/utils.cpp
    src/cauchy.cpp
    src/boundary.cpp
)

target_include_directories(lab4 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# путь к папке с входными данными
target_compile_definitions(lab4 PRIVATE LAB4_DATA_DIR="${PROJECT_SOURCE_DIR}/data/lab4/")
