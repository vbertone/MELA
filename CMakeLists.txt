# BASIC DEFINITIONS ========================================================================

# define minimum version of cmake
cmake_minimum_required (VERSION 3.5)

# define project name and its language
project(mela CXX Fortran)

# define c++ standard and issue all the warning demanded by this standard
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -fPIC")
set(CMAKE_SHARED_LINKER_FLAGS -w)

if (NOT DEFINED CMAKE_MACOSX_RPATH)
   set(CMAKE_MACOSX_RPATH 0)
endif()

# VERSION ==============================================================================

set(mela_VERSION 2.0.2)

# CONFIG SCRIPT =========================================================================

# Configuration script
set(prefix ${CMAKE_INSTALL_PREFIX})
set(libdir "${prefix}/lib")
configure_file("${PROJECT_SOURCE_DIR}/bin/mela-config.in" "${PROJECT_SOURCE_DIR}/bin/mela-config")

# FINALIZE ==================================================================================

# generate list of source files
file(GLOB_RECURSE  source_files src/*f src/*cc)

# build test codes
enable_testing()
add_subdirectory(run)

configure_file(${PROJECT_SOURCE_DIR}/run/Reference.ini ${CMAKE_CURRENT_BINARY_DIR}/run/Reference.ini COPYONLY)

# define target library
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

add_library(mela SHARED ${source_files})
# install
install(TARGETS mela
           RUNTIME DESTINATION bin
           LIBRARY DESTINATION lib)
install(FILES ${PROJECT_SOURCE_DIR}/bin/mela-config DESTINATION bin
PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
