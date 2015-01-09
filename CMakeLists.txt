cmake_minimum_required(VERSION 2.8)
project(panache CXX)
enable_language(C) # Required for find_package(lapack)


###################################################
# A simple macro that sets a boolean
#   based on whether or not a given path string
#   is empty or not. If it is not empty, it adds
#   to the list of include directories
###################################################
macro(PATH_TO_BOOL PATHVAR BOOLVAR)
  if(${PATHVAR} STREQUAL "")
    set(${BOOLVAR} False)
  else()
    set(${BOOLVAR} True)
    set(LIBINT_INCLUDEDIR ${${PATHVAR}}/include)
    set(LIBINT_LIBRARYDIR ${${PATHVAR}}/lib)
  endif()
endmacro()



#####################################
# Various options
#####################################
set(PANACHE_F03_INTERFACE FALSE CACHE BOOL "Build interface to fortran")
set(PANACHE_OPENMP TRUE CACHE BOOL "Enable OpenMP")
set(PANACHE_TIMING FALSE CACHE BOOL "Enable timing of some panache functionality")
set(PANACHE_USE_LIBERD FALSE CACHE BOOL "Use libERD for integrals")
set(PANACHE_USE_SLOWERI FALSE CACHE BOOL "Use slow ERI calculation (for testing)")
set(PANACHE_PROFILE FALSE CACHE BOOL "Enable some profiling code")
set(PANACHE_LAPACK64 FALSE CACHE BOOL "Use the 64-bit interface to Lapack/BLAS. Must be set if you link to 64-bit lapack/BLAS/MKL, etc")
set(PANACHE_INTERFACE64 FALSE CACHE BOOL "Build 64-bit C/Fortran interface")



################################################################
# Path to various optional components
# These should be set by the user using the command line or cmake gui
################################################################
set(LIBINT_PATH "" CACHE FILEPATH "Path to libint directory")
set(LIBINT2_PATH "" CACHE FILEPATH "Path to libint2 directory")
set(CYCLOPS_PATH "" CACHE FILEPATH "Path to Cyclops top directory")

# Convert some paths to boolean variables
PATH_TO_BOOL(LIBINT_PATH PANACHE_USE_LIBINT)
PATH_TO_BOOL(LIBINT2_PATH PANACHE_USE_LIBINT2)
PATH_TO_BOOL(CYCLOPS_PATH PANACHE_CYCLOPS)

# Internally used
set(PANACHE_FORTRAN FALSE) # Build needs a fortran compiler in some way


#########################################
# Check for conflicting options 
# and select LibERD as default if needed
#########################################
set(OPTION_SUM 0)
if(PANACHE_USE_LIBINT)
  math(EXPR OPTION_SUM "${OPTION_SUM}+1")
endif()
if(PANACHE_USE_LIBINT2)
  math(EXPR OPTION_SUM "${OPTION_SUM}+1")
endif()
if(PANACHE_USE_LIBERD)
  math(EXPR OPTION_SUM "${OPTION_SUM}+1")
endif()
if(PANACHE_USE_SLOWERI)
  math(EXPR OPTION_SUM "${OPTION_SUM}+1")
endif()

if(OPTION_SUM GREATER 1)
  message(FATAL_ERROR "Cannot select more than one integral backend! You've selected ${OPTION_SUM}")
endif()

if(OPTION_SUM LESS 1)
  message(WARNING "No integral backed selected. Using libERD")
  set(PANACHE_USE_LIBERD TRUE)
endif()




################################################################
# Handle if we need a Fortran compiler or not
################################################################
if(PANACHE_F03_INTERFACE OR PANACHE_USE_LIBERD)
  message(STATUS "Need a fortran compiler...")
  set(PANACHE_FORTRAN TRUE)
  enable_language(Fortran)
endif()



####################################################
# Build types - I only really use the standard ones
####################################################
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type - release, debug, etc" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug" "RelWithDebInfo" "MinSizeRel" "None")
endif()



#####################################
# Compiler flags
#####################################
# PANACHE_CXX_FLAGS - flags for main library
# PANACHE_F03_INTERFACE_FLAGS - flags for F2003 interface
# LIBERD_F90_FLAGS - flags for libERD
# RUNTEST_CXX_FLAGS - flags for the C++ test program

# PANACHE requires a few C++11 features, namely
#         std::shared_ptr and a few instances of 'auto'
# -DFC_SYMBOL - what fortran symbols to use in libqt
#
# Note that paths should probably be relative to the panache/ subdirectory
# (since that is where add_library() is
# No harm in defining the fortran ones, even if PANACHE_F03_INTERFACE is not set
list(APPEND PANACHE_CXX_FLAGS "-std=c++11")
list(APPEND PANACHE_CXX_FLAGS "-Wall")
list(APPEND PANACHE_CXX_FLAGS "-DFC_SYMBOL=2")
list(APPEND PANACHE_CXX_INCLUDES "${CMAKE_SOURCE_DIR}")


# Ditto for test programs
list(APPEND RUNTEST_CXX_FLAGS "--std=c++11 -Wall")
list(APPEND RUNTEST_CXX_INCLUDES "${CMAKE_SOURCE_DIR}")

# Todo - other compiler types?
message(STATUS "PANACHE FORT: ${PANACHE_FORTRAN}")
if(PANACHE_FORTRAN)
  if("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    list(APPEND PANACHE_F03_INTERFACE_FLAGS "-Wall")
    list(APPEND LIBERD_F90_FLAGS "-Wall")
  elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    list(APPEND PANACHE_F03_INTERFACE_FLAGS "-warn all")
    list(APPEND LIBERD_F90_FLAGS "-warn all")
  else()
    message(WARNING "Unknown Fortran compiler ID: ${CMAKE_Fortran_COMPILER_ID}")
    list(APPEND PANACHE_F03_INTERFACE_FLAGS "-Wall")
    list(APPEND LIBERD_F90_FLAGS "-Wall")
  endif()
endif()


# Standard link libraries
list(APPEND RUNTEST_LINK_LIBRARIES "stdc++" "m")


#####################################
# Fortran interface
#####################################
if(PANACHE_F03_INTERFACE)
  # enable_language is done later
  message(STATUS "Enabling Fortran interface")
endif()



#####################################
# Cyclops & MPI
#####################################
if(PANACHE_CYCLOPS)
  message(STATUS "Enabling Cyclops/MPI for Panache")
  find_package(MPI)
  if(MPI_CXX_FOUND)
    message(STATUS "Found MPI. Remember to link to MPI libraries when building your final project")
    list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_CYCLOPS")
    list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_MPI")
    list(APPEND PANACHE_CXX_FLAGS "${MPI_CXX_COMPILE_FLAGS}") 
    list(APPEND PANACHE_CXX_INCLUDES "${CYCLOPS_PATH}/include") 
    list(APPEND RUNTEST_CXX_FLAGS "-DPANACHE_CYCLOPS")
    list(APPEND RUNTEST_CXX_FLAGS "-DPANACHE_MPI")
    list(APPEND RUNTEST_CXX_FLAGS "${MPI_CXX_COMPILE_FLAGS}") 
    list(APPEND RUNTEST_CXX_INCLUDES "${CYCLOPS_PATH}/include") 
    #list(APPEND LIBERD_F90_FLAGS "${MPI_CXX_FLAGS}")
    #list(APPEND PANACHE_F03_INTERFACE_FLAGS "${MPI_CXX_FLAGS}")
    list(APPEND RUNTEST_LINK_FLAGS ${MPI_CXX_LINK_FLAGS})
    list(APPEND RUNTEST_LINK_LIBRARIES ${MPI_CXX_LIBRARIES})
    list(APPEND RUNTEST_LINK_LIBRARIES "${CYCLOPS_PATH}/lib/libctf.a")
  else()
    message(FATAL_ERROR "MPI was not found! This is required for Cyclops")
  endif()
endif()



#####################################
# OpenMP
#####################################
if(PANACHE_CYCLOPS)
  message(STATUS "Enabling OpenMP because Cyclops is enabled")
  set(PANACHE_OPENMP TRUE)
endif()

if(PANACHE_OPENMP)
  message(STATUS "Enabling OpenMP for Panache")
  find_package(OpenMP)
  if(OPENMP_FOUND)
    message(STATUS "Found OpenMP. Remember to link to OpenMP libraries when building your final project (if needed)")
    list(APPEND PANACHE_CXX_FLAGS "${OpenMP_CXX_FLAGS}") 
    list(APPEND RUNTEST_CXX_LINK_FLAGS "${OpenMP_CXX_FLAGS}") 

    # BUG - cmake doesn't support Fortran in FindOpenMP, so this is a big assumption
    list(APPEND LIBERD_F90_FLAGS "${OpenMP_CXX_FLAGS}") 
  else()
    message(WARNING "OpenMP was not found! Disabling...")
  endif()
endif()

if(PANACHE_CYCLOPS AND NOT PANACHE_OPENMP)
  message(FATAL_ERROR "OpenMP was not found, but is required for Cyclops")
endif()


#####################################
# Timing of panache stuff
#####################################
if(PANACHE_TIMING)
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_TIMING")
  message(STATUS "Enabling timing of panache functionality")
endif()



############################
# Find lapack/blas/mkl, etc
############################
find_package(LAPACK)
if(NOT LAPACK_FOUND)
  message(FATAL_ERROR "LAPACK MODULE NOT FOUND")
else()
  message(STATUS "LAPACK Libraries: ${LAPACK_LIBRARIES}")
  list(APPEND RUNTEST_LINK_LIBRARIES ${LAPACK_LIBRARIES})
endif()



#####################################
# Use the 64-bit or 32-bit integer interface to LAPACK/BLAS/MKL/etc
# If unset, 32-bit interface is used
#####################################
if(PANACHE_LAPACK64)
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_LAPACK64")
  list(APPEND RUNTEST_CXX_FLAGS "-DPANACHE_LAPACK64") 
  message(STATUS "Using 64-bit lapack interface")
else()
  message(STATUS "Using 32-bit lapack interface")
endif()



#####################################
# Whether the C and Fortran interfaces should be 64-bit
#####################################
if(PANACHE_INTERFACE64)
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_INTERFACE64")
  list(APPEND RUNTEST_CXX_FLAGS "-DPANACHE_INTERFACE64")
  message(STATUS "Building 64-bit C/Fortran interface")
else()
  message(STATUS "Building 32-bit C/Fortran interface")
endif()



#####################################
# PANACHE uses the 'restrict' keyword for optimization
# This will test for the proper variation of the keyword
#####################################
include(cmake/Restrict.cmake)
test_restrict(RESTRICT)
list(APPEND PANACHE_CXX_FLAGS "-Drestrict=${RESTRICT}")







################################################################
# Print some messages and set the defines for use in the code
# and add libint library stuff if necessary
################################################################
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")


if(PANACHE_USE_LIBINT)
  set(LIBINT_LIBRARY "${LIBINT_LIBRARYDIR}/libint.a")
  message(STATUS "Using Libint headers from ${LIBINT_INCLUDEDIR}")
  message(STATUS "Using Libint library from ${LIBINT_LIBRARYDIR}")
  message(STATUS "Using Libint library ${LIBINT_LIBRARY}")
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_USE_LIBINT")
  list(APPEND PANACHE_CXX_INCLUDES "${LIBINT_INCLUDEDIR}")
  list(APPEND RUNTEST_LINK_LIBRARIES ${LIBINT_LIBRARY})
endif()


if(PANACHE_USE_LIBINT2)
  # Note - must add directory with header files. libint2.h includes
  # other files via #include <somefile.h>. Therefore, just adding
  # .../include and then using #include <libint2/libint2.h> will not work!
  set(LIBINT_INCLUDEDIR ${LIBINT_INCLUDEDIR}/libint2)
  set(LIBINT_LIBRARY "${LIBINT_LIBRARYDIR}/libint2.a")
  message(STATUS "Using Libint2 headers from ${LIBINT_INCLUDEDIR}")
  message(STATUS "Using Libint2 library from ${LIBINT_LIBRARYDIR}")
  message(STATUS "Using Libint2 library ${LIBINT_LIBRARY}")
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_USE_LIBINT2")
  list(APPEND PANACHE_CXX_INCLUDES "${LIBINT_INCLUDEDIR}")
  list(APPEND RUNTEST_LINK_LIBRARIES ${LIBINT_LIBRARY})
endif()


if(PANACHE_USE_LIBERD)
  message(STATUS "Using internal LibERD")
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_USE_LIBERD")
  list(APPEND RUNTEST_LINK_LIBRARIES panache_erd)
endif()


if(PANACHE_USE_SLOWERI)
  message(STATUS "Using internal SlowERI")
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_USE_SLOWERI")
endif()



######################
# Add subdirectories
######################
add_subdirectory(panache)

if(PANACHE_USE_LIBERD)
  add_subdirectory(liberd)
endif()

add_subdirectory(runtest)

