cmake_minimum_required(VERSION 2.8)
project(libpanache CXX)
#enable_language(Fortran) # liberd

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
    list(APPEND PANACHE_CXX_INCLUDEDIRS ${${PATHVAR}})
  endif(${PATHVAR} STREQUAL "")
endmacro(PATH_TO_BOOL PATHVAR BOOLVAR)


####################################################
# Build types - I only really use the standard ones
####################################################
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type - release, debug, etc" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug" "RelWithDebInfo" "MinSizeRel" "None")
endif(NOT CMAKE_BUILD_TYPE)


#####################################################
# Include the header files from the source directory
#####################################################
list(APPEND PANACHE_CXX_INCLUDEDIRS ${CMAKE_SOURCE_DIR})


#####################################
# Compiler flags
#####################################
# PANACHE requires a few C++11 features, namely
#         std::shared_ptr and a few instances of 'auto'
# -DFC_SYMBOL - what fortran symbols to use in libqt
list(APPEND PANACHE_CXX_FLAGS "--std=c++11")
list(APPEND PANACHE_CXX_FLAGS "-DFC_SYMBOL=2")

# Release, debug, etc
#if(CMAKE_BUILD_TYPE STREQUAL "Debug")
#  list(APPEND PANACHE_CXX_FLAGS "-g" "-O0")
#endif(CMAKE_BUILD_TYPE STREQUAL "Debug")
#
#if(CMAKE_BUILD_TYPE STREQUAL "Release")
#  list(APPEND PANACHE_CXX_FLAGS "-O2")
#endif(CMAKE_BUILD_TYPE STREQUAL "Release")


# Timing of panache stuff
set(PANACHE_TIMING FALSE CACHE BOOL "Enable timing of some panache functionality")
if(PANACHE_TIMING)
  list(APPEND PANACHE_CXX_FLAGS "-DPANACHE_TIMING")
  message(STATUS "Enabling timing of panache functionality")
endif(PANACHE_TIMING)



# Use the 64-bit or 32-bit integer interface to LAPACK/BLAS/MKL/etc
# If unset, 32-bit interface is used
set(USE_64LAPACK FALSE CACHE BOOL "Use the 64-bit interface to Lapack/BLAS. Must be set if you link to 64-bit lapack/BLAS/MKL, etc")

if(USE_64LAPACK)
  list(APPEND PANACHE_CXX_FLAGS "-DUSE_64LAPACK")
  message(STATUS "Using 64-bit lapack interface")
else(USE_64LAPACK)
  message(STATUS "Using 32-bit lapack interface")
endif(USE_64LAPACK)


# Whether the C and Fortran interfaces should be 64-bit
set(USE_64PANACHE FALSE CACHE BOOL "Build 64-bit C/Fortran interface")

if(USE_64PANACHE)
  list(APPEND PANACHE_CXX_FLAGS "-DUSE_64PANACHE")
  message(STATUS "Building 64-bit C/Fortran interface")
else(USE_64PANACHE)
  message(STATUS "Building 32-bit C/Fortran interface")
endif(USE_64PANACHE)


# PANACHE uses the 'restrict' keyword for optimization
# This will test for the proper variation of the keyword
include(cmake/Restrict.cmake)
test_restrict(RESTRICT)
list(APPEND PANACHE_CXX_FLAGS "-Drestrict=${RESTRICT}")




################################################################
# Path to various libints
# These should be set by the user using the command line or cmake gui
################################################################
set(LIBINT_INCLUDE_PATH "" CACHE FILEPATH "Path to libint include directory")
set(LIBINT2_INCLUDE_PATH "" CACHE FILEPATH "Path to libint2 include directory")

# Parse these using the macro
PATH_TO_BOOL(LIBINT_INCLUDE_PATH USE_LIBINT)
PATH_TO_BOOL(LIBINT2_INCLUDE_PATH USE_LIBINT2)

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
################################################################
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")

if(USE_LIBINT)
  message(STATUS "Using Libint headers from ${LIBINT_INCLUDE_PATH}")
  set(LIBINT_FILES "${LIBINT_FILES}" LibintTwoElectronInt.cc LibintERI.cc)
  list(APPEND PANACHE_CXX_FLAGS "-DUSE_LIBINT")
endif(USE_LIBINT)

if(USE_LIBINT2)
  message(STATUS "Using Libint2 headers from ${LIBINT2_INCLUDE_PATH}")
  set(LIBINT_FILES "${LIBINT_FILES}" Libint2TwoElectronInt.cc Libint2ERI.cc)
  list(APPEND PANACHE_CXX_FLAGS "-DUSE_LIBINT2")
endif(USE_LIBINT2)

if(USE_SLOWERI)
  message(STATUS "Using Libint2 headers from ${LIBINT2_INCLUDE_PATH}")
  set(LIBINT_FILES "${LIBINT_FILES}" SlowTwoElectronInt.cc SlowERI.cc SlowERIBase.cc)
  list(APPEND PANACHE_CXX_FLAGS "-DUSE_SLOWERI")
endif(USE_SLOWERI)




######################
# Create the library
######################
add_library(panache 
            AOIntegralsIterator.cc
            AOShellCombinationsIterator.cc
            BasisSet.cc
            BasisSetParser.cc
            CartesianIter.cc
            CharacterTable.cc
            CoordEntry.cc
            DFTensor.cc
            Dimension.cc
            ERDERI.cc
            ERDTwoElectronInt.cc
            FittingMetric.cc
            Fjt.cc
            GaussianShell.cc
            IntegralParameters.cc
            IrreducibleRepresentation.cc
            Math.cc
            Matrix.cc
            Molecule.cc
            Output.cc
            PetiteList.cc
            PointGroup.cc
            ShellInfo.cc
            ShellRotation.cc
            SolidHarmonic.cc
            SphericalTransform.cc
            SymRep.cc
            SymmOps.cc
            SymmetryOperation.cc
            TwoBodyAOInt.cc
            Vector.cc
            block_matrix.cc
            c_interface.cc
            fortran_interface.cc
            c_convert.cc
            dot.cc
            init_array.cc
            int_array.cc
            invert.cc
            lubksb.cc
            ludcmp.cc
            schmidt.cc
            sq_rsp.cc
            sq_to_tri.cc
            zero.cc
            Testing.cc
            Lapack.cc
            ${LIBINT_FILES}
            )

##############################################
# Panache library compiler flags, etc
##############################################
# Compiler definitions list need to be split
string(REPLACE ";" " " PANACHE_CXX_FLAGS "${PANACHE_CXX_FLAGS}")
set_target_properties(panache PROPERTIES COMPILE_FLAGS ${PANACHE_CXX_FLAGS})
set_target_properties(panache PROPERTIES INCLUDE_DIRECTORIES "${PANACHE_CXX_INCLUDEDIRS}")


# FOR NOW
file(GLOB HFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${HFILES} DESTINATION include/panache)
install(TARGETS panache ARCHIVE DESTINATION lib)
