cmake_minimum_required(VERSION 3.18)

project(count_alkanes
    VERSION 1.0
    DESCRIPTION "Count alkyls and alkanes without chiral considered."
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# if(NOT CMAKE_BUILD_TYPE)
#     set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Build type" FORCE)
# endif()


file(GLOB SOURCES
    src/*.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE 
    ${PROJECT_SOURCE_DIR}/include
)


set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

message(STATUS "")
message(STATUS "=== Configuration Summary ===")
message(STATUS "Project: ${PROJECT_NAME}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Extra optimizations: ${ENABLE_EXTRA_OPTIMIZATIONS}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "=================================")
