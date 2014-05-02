cmake_minimum_required(VERSION 2.8)

include_directories(${CMAKE_SOURCE_DIR})

# A simple macro that sets a boolean
#   based on whether or not a given path string
#   is empty or not. If it is not empty, it adds
#   to the list of include directories
MACRO(PATH_TO_BOOL PATHVAR BOOLVAR)
  if(${PATHVAR} STREQUAL "")
    set(${BOOLVAR} False)
  else(${PATHVAR} STREQUAL "")
    set(${BOOLVAR} True)
    include_directories(${${PATHVAR}})
  endif(${PATHVAR} STREQUAL "")
ENDMACRO(PATH_TO_BOOL PATHVAR BOOLVAR)


# PANACHE requires a few C++11 features, namely
#         std::shared_ptr and a few instances of 'auto'
add_definitions("--std=c++11 -DFC_SYMBOL=2 -DPANACHE_DEVELOPER_GENERATE -g -O0")


# PANACHE uses the 'restrict' keyword for optimization
# This will test for the proper variation of the keyword
include(cmake/Restrict.cmake)
test_restrict(RESTRICT)
add_definitions("-Drestrict=${RESTRICT}")


# Path to various libints
# These should be set by the user using the command line or gui
set(LIBINT_INCLUDE_PATH "" CACHE FILEPATH "Path to libint include directory")
set(LIBINT2_INCLUDE_PATH "" CACHE FILEPATH "Path to libint2 include directory")


# Parse these using the macro
PATH_TO_BOOL(LIBINT_INCLUDE_PATH USE_LIBINT)
PATH_TO_BOOL(LIBINT2_INCLUDE_PATH USE_LIBINT2)

# You can't use both versions of libint...
if(USE_LIBINT AND USE_LIBINT2)
  MESSAGE(FATAL_ERROR "Cannot use both Libint and Libint2!")
endif(USE_LIBINT AND USE_LIBINT2)


# ... and you have to use at least one
if(NOT(USE_LIBINT OR USE_LIBINT2))
  MESSAGE(FATAL_ERROR "You must specify either the path to Libint or Libint2 header files!")
endif(NOT(USE_LIBINT OR USE_LIBINT2))


# Print some messages and set the defines for use in the code
# Also add .cc files to be compiled based on the libint version
if(USE_LIBINT)
  MESSAGE(STATUS "Using Libint headers from ${LIBINT_INCLUDE_PATH}")
  SET(LIBINT_FILES "${LIBINT_FILES}" TwoElectronInt.cc ERI.cc)
  add_definitions("-DUSE_LIBINT")
endif(USE_LIBINT)

if(USE_LIBINT2)
  MESSAGE(STATUS "Using Libint2 headers from ${LIBINT2_INCLUDE_PATH}")
  SET(LIBINT_FILES "${LIBINT_FILES}" TwoElectronInt2.cc ERI2.cc)
  add_definitions("-DUSE_LIBINT2")
endif(USE_LIBINT2)


# This is needed for the testing executable (for now)
find_package(LAPACK)
if(NOT LAPACK_FOUND)
  MESSAGE(FATAL_ERROR "LAPACK MODULE NOT FOUND")
else(NOT LAPACK_FOUND)
  MESSAGE(STATUS "LAPACK Libraries: ${LAPACK_LIBRARIES}")
endif(NOT LAPACK_FOUND)


# Create the library
add_library(panache 
            AOIntegralsIterator.cc
            AOShellCombinationsIterator.cc
            BasisSet.cc
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
            SlowTwoElectronInt.cc
            SolidHarmonic.cc
            SphericalTransform.cc
            SymRep.cc
            SymmOps.cc
            SymmetryOperation.cc
            TwoBodyAOInt.cc
            Vector.cc
            blas_intfc.cc
            blas_intfc23.cc
            block_matrix.cc
            c_interface.cc
            c_convert.cc
            dot.cc
            init_array.cc
            int_array.cc
            invert.cc
            lapack_intfc.cc
            lubksb.cc
            ludcmp.cc
            schmidt.cc
            sq_rsp.cc
            sq_to_tri.cc
            zero.cc
            Testing.cc
            ${LIBINT_FILES}
            )

add_executable(runtest runtest.cc)
target_link_libraries(runtest panache ${LAPACK_LIBRARIES})


# FOR NOW
FILE(GLOB HFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
INSTALL(FILES ${HFILES} DESTINATION include/panache)
#INSTALL(FILES GLOB *.h
#              BasisSet.h 
#              Molecule.h
#              c_interface.h
#              Matrix.h
#              DFTensor.h
#              Dimension.h
#              Exception.h
#              DESTINATION include/panache)

INSTALL(TARGETS panache
        ARCHIVE DESTINATION lib)
