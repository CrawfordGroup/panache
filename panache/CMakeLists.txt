add_library(panache 
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
            fortran/fortran_interface.f90
            c_convert.cc
            Lapack.cc
            Reorder.cc
            storedqtensor/StoredQTensor.cc
            storedqtensor/LocalQTensor.cc
            storedqtensor/MemoryQTensor.cc
            storedqtensor/DiskQTensor.cc
            storedqtensor/StoredQTensorFactory.cc
            ${LIBINT_FILES}
            ${PARALLEL_FILES}
            ${CYCLOPS_FILES}
            )

##############################################
# Panache library compiler flags, etc
##############################################
# Compiler definitions list need to be split
string(REPLACE ";" " " PANACHE_CXX_FLAGS "${PANACHE_CXX_FLAGS}")
set_target_properties(panache PROPERTIES COMPILE_FLAGS ${PANACHE_CXX_FLAGS})
set_target_properties(panache PROPERTIES INCLUDE_DIRECTORIES "${PANACHE_CXX_INCLUDES}")

# For now, install all header files
file(GLOB HFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${HFILES} DESTINATION include/panache)
install(TARGETS panache ARCHIVE DESTINATION lib)
