cmake_minimum_required(VERSION 3.14)

project(kanaval
    VERSION 1.0.0
    DESCRIPTION "C++ library for validating kana format files"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

add_library(kanaval INTERFACE)

add_subdirectory(extern)

target_include_directories(kanaval INTERFACE include/)

target_link_libraries(kanaval INTERFACE millijson)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    include(CTest)
    if(BUILD_TESTING)
        add_subdirectory(tests)
    endif()
endif()
