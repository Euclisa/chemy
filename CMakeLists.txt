cmake_minimum_required(VERSION 3.16)

project(ChemGraph VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

file(GLOB SOURCES CMAKE_CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/server/src/*.cpp")

find_package(PkgConfig REQUIRED)
pkg_check_modules(LIBPQXX REQUIRED libpqxx)
find_package(nlohmann_json REQUIRED)
find_package(fmt REQUIRED)


include(FetchContent)

FetchContent_Declare(
    croaring
    GIT_REPOSITORY https://github.com/RoaringBitmap/CRoaring.git
    GIT_TAG v4.4.2
)

FetchContent_MakeAvailable(croaring)


# Common include directories and libraries
set(COMMON_INCLUDES
    ${LIBPQXX_INCLUDE_DIRS}
    "${CMAKE_SOURCE_DIR}/server/include/"
    "/usr/include/rdkit"
    "${croaring_SOURCE_DIR}/include"
    "${croaring_SOURCE_DIR}/cpp"
)
set(COMMON_LIBS
    ${LIBPQXX_LIBRARIES}
    fmt::fmt
    RDKitSmilesParse
    RDKitGraphMol
    RDKitDepictor
    RDKitMolDraw2D
    roaring
)


add_executable(chem_debug ${SOURCES})
target_include_directories(chem_debug PRIVATE ${COMMON_INCLUDES})
target_link_libraries(chem_debug PRIVATE ${COMMON_LIBS})
target_compile_options(chem_debug PRIVATE -Wall -Wextra -Wpedantic -O0 -g)
target_compile_definitions(chem_debug
    PRIVATE PROJECT_DIR="${CMAKE_SOURCE_DIR}"
)


add_executable(chem_release ${SOURCES})
target_include_directories(chem_release PRIVATE ${COMMON_INCLUDES})
target_link_libraries(chem_release PRIVATE ${COMMON_LIBS})
target_compile_options(chem_release PRIVATE -Wall -Wextra -Wpedantic -O3 -march=native -DNDEBUG)
target_compile_definitions(chem_release
    PRIVATE PROJECT_DIR="${CMAKE_SOURCE_DIR}"
)
