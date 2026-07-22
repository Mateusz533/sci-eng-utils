cmake_minimum_required(VERSION 3.21)

project(MultiDemo)

include(CTest)

# Configure language settings
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(COMMON_COMPILE_OPTIONS
    -Wall # basic warnings
    -Wextra # extra warnings
    # -Wpedantic         # strict ISO C++ compliance
    -Werror # treat all warnings as errors
    # -Wconversion       # implicit conversions
    # -Wsign-conversion  # sign conversions
    -Wshadow # variable shadowing
    -Wnull-dereference # null pointer dereference
    -Wold-style-cast # C-style casts
    -Wnon-virtual-dtor # non-virtual destructors in polymorphic classes
    -Woverloaded-virtual
    # -Wdouble-promotion
    -Wformat=2
    -Wundef # use of undefined macros
)

# Link and build everything in 'main' directory
set(INCLUDE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/main")

file(
    GLOB
    SOURCE_FILES
    CONFIGURE_DEPENDS
    "${SOURCE_DIRECTORY}/*.cpp"
)

foreach(source_file IN LISTS SOURCE_FILES)
    get_filename_component(executable_name ${source_file} NAME_WE)
    add_executable(${executable_name} ${source_file})
    target_include_directories(${executable_name} PRIVATE ${INCLUDE_DIRECTORY})
    target_compile_options(${executable_name} PRIVATE ${COMMON_COMPILE_OPTIONS})
endforeach()

# Configure and build tests
if(BUILD_TESTING)
    find_package(GTest REQUIRED)

    add_executable(RunTests "${CMAKE_CURRENT_SOURCE_DIR}/tests/MainTest.cpp")
    target_link_libraries(RunTests PRIVATE GTest::gtest)
    target_include_directories(RunTests PRIVATE ${INCLUDE_DIRECTORY} "${CMAKE_CURRENT_SOURCE_DIR}/tests/include")
    target_compile_options(RunTests PRIVATE ${COMMON_COMPILE_OPTIONS})

    include(GoogleTest)
    gtest_discover_tests(RunTests)
endif()
