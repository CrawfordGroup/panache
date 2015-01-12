##############################################
# Add files
##############################################
set(PANACHE_CXX_FILES
            AOIntegralsIterator.cc
            AOShellCombinationsIterator.cc
            BasisSet.cc
            BasisSetParser.cc
            CartesianIter.cc
            ThreeIndexTensor.cc
            DFTensor.cc
            CHTensor.cc
            FittingMetric.cc
            Fjt.cc
            GaussianShell.cc
            IntegralParameters.cc
            Math.cc
            Output.cc
            ShellInfo.cc
            SolidHarmonic.cc
            SphericalTransform.cc
            TwoBodyAOInt.cc
            c_interface.cc
            c_convert.cc
            Lapack.cc
            Reorder.cc
            storedqtensor/StoredQTensor.cc
            storedqtensor/LocalQTensor.cc
            storedqtensor/MemoryQTensor.cc
            storedqtensor/DiskQTensor.cc
            storedqtensor/StoredQTensorFactory.cc
)

##################################################
# Handle which files to include in the compilation
##################################################
if(PANACHE_CYCLOPS)
  list(APPEND PANACHE_CXX_FILES MPI.cc storedqtensor/CyclopsQTensor.cc)
else()
  list(APPEND PANACHE_CXX_FILES Serial.cc)
endif()

if(PANACHE_USE_LIBINT)
  list(APPEND PANACHE_CXX_FILES LibintTwoElectronInt.cc LibintERI.cc)
endif()

if(PANACHE_USE_LIBINT2)
  list(APPEND PANACHE_CXX_FILES Libint2TwoElectronInt.cc Libint2ERI.cc)
endif()

if(PANACHE_USE_LIBERD)
  list(APPEND PANACHE_CXX_FILES ERDTwoElectronInt.cc ERDERI.cc)
endif()

if(PANACHE_USE_SLOWERI)
  list(APPEND PANACHE_CXX_FILES SlowTwoElectronInt.cc SlowERI.cc SlowERIBase.cc)
endif()

if(PANACHE_F03_INTERFACE)
  set(PANACHE_F03_INTERFACE_FILES fortran/fortran_interface.f90)
endif()


##################################################
# Actual library
##################################################
add_library(panache 
            ${PANACHE_CXX_FILES}
            ${PANACHE_F03_INTERFACE_FILES}
            )

##############################################
# Panache library compiler flags, etc
##############################################
# Compiler definitions list need to be split
string(REPLACE ";" " " PANACHE_CXX_FLAGS "${PANACHE_CXX_FLAGS}")
set_source_files_properties(${PANACHE_CXX_FILES} PROPERTIES COMPILE_FLAGS "${PANACHE_CXX_FLAGS}")

# This will set the include directories for the fortran stuff as well, but it's
# all we can do. set_source_files_properties doesn't work with INCLUDE_DIRECTORIES
set_target_properties(panache PROPERTIES INCLUDE_DIRECTORIES "${PANACHE_CXX_INCLUDES}")
if(PANACHE_F03_INTERFACE)
  string(REPLACE ";" " " PANACHE_F03_INTERFACE_FLAGS "${PANACHE_F03_INTERFACE_FLAGS}")
  set_source_files_properties(${PANACHE_F03_INTERFACE_FILES} PROPERTIES COMPILE_FLAGS ${PANACHE_F03_INTERFACE_FLAGS})
endif()


# For now, install all header files
file(GLOB HFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${HFILES} DESTINATION include/panache)
install(TARGETS panache ARCHIVE DESTINATION lib)
