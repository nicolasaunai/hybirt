cmake_minimum_required (VERSION 3.20.1) 

project(HYBIRT VERSION 0.1 LANGUAGES CXX C)

# Release mode default
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()


set(HYBIRT_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


set(SOURCE_INC
   src/ampere.hpp
   src/boris.hpp
   src/faraday.hpp
   src/field.hpp
   src/ohm.hpp
   src/pusher.hpp
   src/utils.hpp
   src/vecfield.hpp
)

set (SOURCE_CPP
    src/hybirt.cpp
)

add_executable(hybirt ${SOURCE_INC} ${SOURCE_CPP})




