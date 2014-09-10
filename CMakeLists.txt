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
# LIBERD_F90_FLAGS - flags for libERD
# RUNTEST_CXX_FLAGS - flags for the C++ test program
# FRUNTEST_F90_FLAGS - flags for the Fortran90 test program

# PANACHE requires a few C++11 features, namely
#         std::shared_ptr and a few instances of 'auto'
# -DFC_SYMBOL - what fortran symbols to use in libqt
list(APPEND PANACHE_CXX_FLAGS "--std=c++11")
list(APPEND PANACHE_CXX_FLAGS "-Wall")
list(APPEND PANACHE_CXX_FLAGS "-DFC_SYMBOL=2")
list(APPEND PANACHE_CXX_INCLUDES "${CMAKE_SOURCE_DIR}")


# Ditto for test programs
list(APPEND RUNTEST_CXX_FLAGS "--std=c++11 -Wall")
list(APPEND RUNTEST_CXX_INCLUDES "${CMAKE_SOURCE_DIR}")

# Todo - other compiler types?
if(CMAKE_COMPILER_IS_GNUG77)
  list(APPEND FRUNTEST_F90_FLAGS "-Wall")
  list(APPEND LIBERD_F90_FLAGS "-Wall")
else(CMAKE_COMPILER_IS_GNUG77)
  list(APPEND FRUNTEST_F90_FLAGS "-warn all")
  list(APPEND LIBERD_F90_FLAGS "-warn all")
endif(CMAKE_COMPILER_IS_GNUG77)

# Standard link libraries
list(APPEND RUNTEST_LINK_LIBRARIES "stdc++" "m")
list(APPEND FRUNTEST_LINK_LIBRARIES "stdc++" "m")


#####################################
# OpenMP
#####################################
set(PANACHE_OPENMP TRUE CACHE BOOL "Enable OpenMP")
if(PANACHE_OPENMP)
  message(STATUS "Enabling OpenMP for Panache")
  find_package(OpenMP)
  if(OPENMP_FOUND)
    message(STATUS "Found OpenMP. Remember to link to OpenMP libraries when building your final project (if needed)")
    list(APPEND PANACHE_CXX_FLAGS "${OpenMP_CXX_FLAGS}") 
    list(APPEND RUNTEST_CXX_LINK_FLAGS "${OpenMP_CXX_FLAGS}") 

    # BUG - cmake doesn't support Fortran in FindOpenMP, so this is a big assumption
    list(APPEND FRUNTEST_F90_FLAGS "${OpenMP_CXX_FLAGS}") 
    list(APPEND LIBERD_F90_FLAGS "${OpenMP_CXX_FLAGS}") 
  else(OPENMP_FOUND)
    message(WARNING "OpenMP was not found! Disabling...")
  endif(OPENMP_FOUND)
endif(PANACHE_OPENMP)


#####################################
# Cyclops & MPI
#####################################
set(CYCLOPS_PATH "" CACHE FILEPATH "Path to Cyclops top directory")
PATH_TO_BOOL(CYCLOPS_PATH PANACHE_CYCLOPS)
if(PANACHE_CYCLOPS)
  message(STATUS "Enabling Cyclops/MPI for Panache")
  find_package(MPI)
  if(MPI_CXX_FOUND)
    message(STATUS "Found MPI. Remember to link to MPI libraries when building your final project")
    list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_CYCLOPS")
    list(APPEND PANACHE_CXX_FLAGS "${MPI_CXX_COMPILE_FLAGS}") 
    list(APPEND RUNTEST_CXX_INCLUDES "${CYCLOPS_PATH}/include") 
    #list(APPEND FRUNTEST_F90_FLAGS "${MPI_CXX_FLAGS}") 
    #list(APPEND LIBERD_F90_FLAGS "${MPI_CXX_FLAGS}")
    list(APPEND RUNTEST_LINK_FLAGS ${MPI_CXX_LINK_FLAGS})
    list(APPEND RUNTEST_LINK_LIBRARIES ${MPI_CXX_LIBRARIES})
    list(APPEND RUNTEST_LINK_LIBRARIES "${CYCLOPS_PATH}/lib/libctf.a")
  else(MPI_CXX_FOUND)
    message(ERROR "MPI was not found! This is required for Cyclops")
  endif(MPI_CXX_FOUND)
endif(PANACHE_CYCLOPS)


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
  list(APPEND RUNTEST_LINK_LIBRARIES ${LAPACK_LIBRARIES})
endif(NOT LAPACK_FOUND)



#####################################
# Use the 64-bit or 32-bit integer interface to LAPACK/BLAS/MKL/etc
# If unset, 32-bit interface is used
#####################################
set(PANACHE_LAPACK64 FALSE CACHE BOOL "Use the 64-bit interface to Lapack/BLAS. Must be set if you link to 64-bit lapack/BLAS/MKL, etc")

if(PANACHE_LAPACK64)
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_LAPACK64")
  list(APPEND RUNTEST_CXX_FLAGS "-DPANACHE_LAPACK64") 
  message(STATUS "Using 64-bit lapack interface")
else(PANACHE_LAPACK64)
  message(STATUS "Using 32-bit lapack interface")
endif(PANACHE_LAPACK64)



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
set(PANACHE_INTERFACE64 FALSE CACHE BOOL "Build 64-bit C/Fortran interface")

if(PANACHE_INTERFACE64)
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_INTERFACE64")
  list(APPEND RUNTEST_CXX_FLAGS "-DPANACHE_INTERFACE64")
  message(STATUS "Building 64-bit C/Fortran interface")
else(PANACHE_INTERFACE64)
  message(STATUS "Building 32-bit C/Fortran interface")
endif(PANACHE_INTERFACE64)



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
PATH_TO_BOOL(LIBINT_PATH PANACHE_USE_LIBINT)
PATH_TO_BOOL(LIBINT2_PATH PANACHE_USE_LIBINT2)

# You can't use both versions of libint...
if(PANACHE_USE_LIBINT AND PANACHE_USE_LIBINT2)
  message(FATAL_ERROR "Cannot use both Libint and Libint2!")
endif(PANACHE_USE_LIBINT AND PANACHE_USE_LIBINT2)

# Or libint + liberd
if(PANACHE_USE_LIBINT AND PANACHE_USE_LIBERD)
  message(FATAL_ERROR "Cannot use both Libint and LibERD!")
endif(PANACHE_USE_LIBINT AND PANACHE_USE_LIBERD)

if(PANACHE_USE_LIBINT2 AND PANACHE_USE_LIBERD)
  message(FATAL_ERROR "Cannot use both Libint2 and LibERD!")
endif(PANACHE_USE_LIBINT2 AND PANACHE_USE_LIBERD)


# Warn if not using something faster
if(NOT(PANACHE_USE_LIBINT OR PANACHE_USE_LIBINT2 OR PANACHE_USE_LIBERD))
  message(WARNING "Using very slow ERI code. You probably want to use something else...")
  set(PANACHE_USE_SLOWERI True)
endif(NOT(PANACHE_USE_LIBINT OR PANACHE_USE_LIBINT2 OR PANACHE_USE_LIBERD))



################################################################
# Print some messages and set the defines for use in the code
# Also add .cc files to be compiled based on the libint version
# as well as libint library filename
################################################################
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")


if(PANACHE_USE_LIBINT)
  set(LIBINT_LIBRARY "${LIBINT_LIBRARYDIR}/libint.a")
  set(LIBINT_FILES "${LIBINT_FILES}" LibintTwoElectronInt.cc LibintERI.cc)
  message(STATUS "Using Libint headers from ${LIBINT_INCLUDEDIR}")
  message(STATUS "Using Libint library from ${LIBINT_LIBRARYDIR}")
  message(STATUS "Using Libint library ${LIBINT_LIBRARY}")
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_USE_LIBINT")
  list(APPEND PANACHE_CXX_INCLUDES "${LIBINT_INCLUDEDIR}")
  list(APPEND RUNTEST_LINK_LIBRARIES ${LIBINT_LIBRARY})
endif(PANACHE_USE_LIBINT)


if(PANACHE_USE_LIBINT2)
  # Note - must add directory with header files. libint2.h includes
  # other files via #include <somefile.h>. Therefore, just adding
  # .../include and then using #include <libint2/libint2.h> will not work!
  set(LIBINT_INCLUDEDIR ${LIBINT_INCLUDEDIR}/libint2)
  set(LIBINT_LIBRARY "${LIBINT_LIBRARYDIR}/libint2.a")
  set(LIBINT_FILES "${LIBINT_FILES}" Libint2TwoElectronInt.cc Libint2ERI.cc)
  message(STATUS "Using Libint2 headers from ${LIBINT_INCLUDEDIR}")
  message(STATUS "Using Libint2 library from ${LIBINT_LIBRARYDIR}")
  message(STATUS "Using Libint2 library ${LIBINT_LIBRARY}")
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_USE_LIBINT2")
  list(APPEND PANACHE_CXX_INCLUDES "${LIBINT_INCLUDEDIR}")
  list(APPEND RUNTEST_LINK_LIBRARIES ${LIBINT_LIBRARY})
endif(PANACHE_USE_LIBINT2)


if(PANACHE_USE_LIBERD)
  message(STATUS "Using internal LibERD")
  set(LIBINT_FILES "${LIBINT_FILES}" ERDTwoElectronInt.cc ERDERI.cc)
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_USE_LIBERD")
  list(APPEND RUNTEST_LINK_LIBRARIES panache_erd)
  list(APPEND FRUNTEST_LINK_LIBRARIES panache_erd)
endif(PANACHE_USE_LIBERD)


if(PANACHE_USE_SLOWERI)
  message(STATUS "Using internal SlowERI")
  set(LIBINT_FILES "${LIBINT_FILES}" SlowTwoElectronInt.cc SlowERI.cc SlowERIBase.cc)
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_USE_SLOWERI")
endif(PANACHE_USE_SLOWERI)




######################
# Add subdirectories
######################
add_subdirectory(panache)

if(PANACHE_USE_LIBERD)
add_subdirectory(liberd)
endif(PANACHE_USE_LIBERD)

add_subdirectory(runtest)



