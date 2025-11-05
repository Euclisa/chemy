cmake_minimum_required(VERSION 3.16)

project(ChemGraph VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

file(GLOB SOURCES CMAKE_CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/server/src/*.cpp")

find_package(PkgConfig REQUIRED)
pkg_check_modules(LIBPQXX REQUIRED libpqxx)
find_package(nlohmann_json REQUIRED)
find_package(fmt REQUIRED)

# Common include directories and libraries
set(COMMON_INCLUDES
    ${LIBPQXX_INCLUDE_DIRS}
    "${CMAKE_SOURCE_DIR}/server/include/"
)
set(COMMON_LIBS ${LIBPQXX_LIBRARIES})


add_executable(chem_debug ${SOURCES})
target_include_directories(chem_debug PRIVATE ${COMMON_INCLUDES})
target_link_libraries(chem_debug PRIVATE ${COMMON_LIBS} fmt::fmt)
target_compile_options(chem_debug PRIVATE -Wall -Wextra -Wpedantic -O0 -g)


add_executable(chem_release ${SOURCES})
target_include_directories(chem_release PRIVATE ${COMMON_INCLUDES})
target_link_libraries(chem_release PRIVATE ${COMMON_LIBS} fmt::fmt)
target_compile_options(chem_release PRIVATE -Wall -Wextra -Wpedantic -O3 -march=native -DNDEBUG)