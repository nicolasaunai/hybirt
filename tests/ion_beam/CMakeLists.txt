cmake_minimum_required(VERSION 3.20.1)
project(test-ion_beam)
set(SOURCES test_ionbeam.cpp
    ${CMAKE_SOURCE_DIR}/src/pusher.hpp
    ${CMAKE_SOURCE_DIR}/src/vecfield.hpp
    ${CMAKE_SOURCE_DIR}/src/field.hpp
    ${CMAKE_SOURCE_DIR}/src/gridlayout.hpp
)
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} PRIVATE HighFive)
message(${PROJECT_NAME} " target: ${SOURCES}")
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src/)
