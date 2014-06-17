cmake_minimum_required(VERSION 2.8)
project(panache CXX)
enable_language(C) # Required for find_package(lapack)
enable_language(Fortran) # Test program


###################################################
# A simple macro that sets a boolean
#   based on whether or not a given path string
#   is empty or not. If it is not empty, it adds
#   to the list of include directories
###################################################
macro(PATH_TO_BOOL PATHVAR BOOLVAR)
  if(${PATHVAR} STREQUAL "")
    set(${BOOLVAR} False)
  else(${PATHVAR} STREQUAL "")
    set(${BOOLVAR} True)
    set(LIBINT_INCLUDEDIR ${${PATHVAR}}/include)
    set(LIBINT_LIBRARYDIR ${${PATHVAR}}/lib)
  endif(${PATHVAR} STREQUAL "")
endmacro(PATH_TO_BOOL PATHVAR BOOLVAR)


####################################################
# Build types - I only really use the standard ones
####################################################
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type - release, debug, etc" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug" "RelWithDebInfo" "MinSizeRel" "None")
endif(NOT CMAKE_BUILD_TYPE)


#####################################
# Compiler flags
#####################################
# PANACHE_CXX_FLAGS - flags for main library
# RUNTEST_CXX_FLAGS - flags for the C++ test program
# FRUNTEST_F90_FLAGS - flags for the Fortran90 test program

# PANACHE requires a few C++11 features, namely
#         std::shared_ptr and a few instances of 'auto'
# -DFC_SYMBOL - what fortran symbols to use in libqt
list(APPEND PANACHE_CXX_FLAGS "--std=c++11")
list(APPEND PANACHE_CXX_FLAGS "-Wall")
list(APPEND PANACHE_CXX_FLAGS "-DFC_SYMBOL=2")

# Ditto for test programs
list(APPEND RUNTEST_CXX_FLAGS "--std=c++11 -Wall")
list(APPEND FRUNTEST_F90_FLAGS "-warn all")



#####################################
# OpenMP
#####################################
set(PANACHE_OPENMP TRUE CACHE BOOL "Enable OpenMP")
if(PANACHE_OPENMP)
  message(STATUS "Enabling OpenMP for Panache")
  find_package(OpenMP)
  if(OPENMP_FOUND)
    message(STATUS "Found OpenMP. Remember to link to OpenMP libraries when building your final project")
    list(APPEND PANACHE_CXX_FLAGS "${OpenMP_CXX_FLAGS}") 
    list(APPEND RUNTEST_CXX_FLAGS "${OpenMP_CXX_FLAGS}") 
    # BUG - cmake doesn't support Fortran in FindOpenMP, so this is a big assumption
    list(APPEND FRUNTEST_F90_FLAGS "${OpenMP_CXX_FLAGS}") 
  else(OPENMP_FOUND)
    message(WARNING "OpenMP was not found! Disabling...")
  endif(OPENMP_FOUND)
endif(PANACHE_OPENMP)



#####################################
# Timing of panache stuff
#####################################
set(PANACHE_TIMING FALSE CACHE BOOL "Enable timing of some panache functionality")
if(PANACHE_TIMING)
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_TIMING")
  message(STATUS "Enabling timing of panache functionality")
endif(PANACHE_TIMING)



############################
# Find lapack/blas/mkl, etc
############################
find_package(LAPACK)
if(NOT LAPACK_FOUND)
  message(FATAL_ERROR "LAPACK MODULE NOT FOUND")
else(NOT LAPACK_FOUND)
  message(STATUS "LAPACK Libraries: ${LAPACK_LIBRARIES}")
endif(NOT LAPACK_FOUND)



#####################################
# Use the 64-bit or 32-bit integer interface to LAPACK/BLAS/MKL/etc
# If unset, 32-bit interface is used
#####################################
set(USE_64LAPACK FALSE CACHE BOOL "Use the 64-bit interface to Lapack/BLAS. Must be set if you link to 64-bit lapack/BLAS/MKL, etc")

if(USE_64LAPACK)
  list(APPEND PANACHE_CXX_FLAGS "-DUSE_64LAPACK")
  message(STATUS "Using 64-bit lapack interface")
else(USE_64LAPACK)
  message(STATUS "Using 32-bit lapack interface")
endif(USE_64LAPACK)



#####################################
# Prefetching from disk
#####################################
set(PANACHE_DISKPREFETCH FALSE CACHE BOOL "Prefetch next batch from disk in separate thread")

if(PANACHE_DISKPREFETCH)
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_DISKPREFETCH")
  message(WARNING "Enabling disk prefetching. THIS IS EXPERIMENTAL AND MAY REQUIRE A VERY NEW COMPILER")
endif(PANACHE_DISKPREFETCH)



#####################################
# Whether the C and Fortran interfaces should be 64-bit
#####################################
set(USE_64PANACHE FALSE CACHE BOOL "Build 64-bit C/Fortran interface")

if(USE_64PANACHE)
  list(APPEND PANACHE_CXX_FLAGS "-DUSE_64PANACHE")
  message(STATUS "Building 64-bit C/Fortran interface")
else(USE_64PANACHE)
  message(STATUS "Building 32-bit C/Fortran interface")
endif(USE_64PANACHE)



#####################################
# PANACHE uses the 'restrict' keyword for optimization
# This will test for the proper variation of the keyword
#####################################
include(cmake/Restrict.cmake)
test_restrict(RESTRICT)
list(APPEND PANACHE_CXX_FLAGS "-Drestrict=${RESTRICT}")




################################################################
# Path to various libints
# These should be set by the user using the command line or cmake gui
################################################################
set(LIBINT_PATH "" CACHE FILEPATH "Path to libint directory")
set(LIBINT2_PATH "" CACHE FILEPATH "Path to libint2 directory")

# Parse these using the macro
PATH_TO_BOOL(LIBINT_PATH USE_LIBINT)
PATH_TO_BOOL(LIBINT2_PATH USE_LIBINT2)

# You can't use both versions of libint...
if(USE_LIBINT AND USE_LIBINT2)
  message(FATAL_ERROR "Cannot use both Libint and Libint2!")
endif(USE_LIBINT AND USE_LIBINT2)


# Warn if not using something faster
if(NOT(USE_LIBINT OR USE_LIBINT2))
  message(WARNING "Using very slow ERI code. You probably want to use something else...")
  set(USE_SLOWERI True)
endif(NOT(USE_LIBINT OR USE_LIBINT2))



################################################################
# Print some messages and set the defines for use in the code
# Also add .cc files to be compiled based on the libint version
# as well as libint library filename
################################################################
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

if(USE_LIBINT)
  set(LIBINT_LIBRARY "${LIBINT_LIBRARYDIR}/libint.a")
  set(LIBINT_FILES "${LIBINT_FILES}" LibintTwoElectronInt.cc LibintERI.cc)
  message(STATUS "Using Libint headers from ${LIBINT_INCLUDEDIR}")
  message(STATUS "Using Libint library from ${LIBINT_LIBRARYDIR}")
  message(STATUS "Using Libint library ${LIBINT_LIBRARY}")
  list(APPEND PANACHE_CXX_FLAGS "-DUSE_LIBINT")
endif(USE_LIBINT)

if(USE_LIBINT2)
  # Note - must add directory with header files. libint2.h includes
  # other files via #include <somefile.h>. Therefore, just adding
  # .../include and then using #include <libint2/libint2.h> will not work!
  set(LIBINT_INCLUDEDIR ${LIBINT_INCLUDEDIR}/libint2)
  set(LIBINT_LIBRARY "${LIBINT_LIBRARYDIR}/libint2.a")
  set(LIBINT_FILES "${LIBINT_FILES}" Libint2TwoElectronInt.cc Libint2ERI.cc)
  message(STATUS "Using Libint2 headers from ${LIBINT_INCLUDEDIR}")
  message(STATUS "Using Libint2 library from ${LIBINT_LIBRARYDIR}")
  message(STATUS "Using Libint2 library ${LIBINT_LIBRARY}")
  list(APPEND PANACHE_CXX_FLAGS "-DUSE_LIBINT2")
endif(USE_LIBINT2)

if(USE_SLOWERI)
  message(STATUS "Using internal SlowERI")
  set(LIBINT_FILES "${LIBINT_FILES}" SlowTwoElectronInt.cc SlowERI.cc SlowERIBase.cc)
  list(APPEND PANACHE_CXX_FLAGS "-DUSE_SLOWERI")
endif(USE_SLOWERI)




######################
# Add subdirectories
######################
add_subdirectory(panache)
add_subdirectory(runtest)



