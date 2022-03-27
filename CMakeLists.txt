cmake_minimum_required(VERSION 3.14)

project(kanaval
    VERSION 1.0.0
    DESCRIPTION "C++ library for validating kana format files"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(kanaval INTERFACE)

target_include_directories(kanaval INTERFACE include/)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif()
endif()
