cmake_minimum_required(VERSION 3.21)

# Configure language settings
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) 
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -funroll-loops -DNDEBUG")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -march=native -funroll-loops -DNDEBUG")

add_compile_options(
  -Wall              # basic warnings
  -Wextra            # extra warnings
  # -Wpedantic         # strict ISO C++ compliance
  # -Werror            # treat all warnings as errors
  # -Wconversion       # implicit conversions
  # -Wsign-conversion  # sign conversions
  -Wshadow           # variable shadowing
  -Wnull-dereference # null pointer dereference
  # -Wold-style-cast   # C-style casts
  -Wnon-virtual-dtor # non-virtual destructors in polymorphic classes
  -Woverloaded-virtual
  # -Wdouble-promotion
  -Wformat=2
  -Wundef            # use of undefined macros
)

set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH ON)
project(MultiDemo)

# Configure requirements for particular programs
find_package(PkgConfig REQUIRED)
if(NOT (PKGCONFIG_FOUND))
  message(FATAL_ERROR "Please Install PPkgConfig: CMake will Exit")
endif()
set(ENV{PKG_CONFIG_PATH})

# Link and build everything in 'main' directory
set(SOURCE_DIR "main")
set(INCLUDE_DIR "include")
file(GLOB SOURCE_FILES "${SOURCE_DIR}/*.cpp")
set(LINK_LIBRARIES "-pthread")

foreach(source_file ${SOURCE_FILES})
  get_filename_component(executable_name ${source_file} NAME_WE)
  add_executable(${executable_name} ${source_file})
  target_include_directories(${executable_name} PUBLIC ${INCLUDE_DIR})
  target_link_libraries(${executable_name} ${LINK_LIBRARIES})
endforeach()

# Configure and build tests
enable_testing()

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

project(RunTests)

add_executable(${PROJECT_NAME} "tests/MainTest.cpp")
target_link_libraries(${PROJECT_NAME} ${GTEST_LIBRARIES} pthread)
target_include_directories(${PROJECT_NAME} PUBLIC "include")

add_test(NAME GeneralTest COMMAND ${PROJECT_NAME})
