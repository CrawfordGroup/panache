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
            FittingMetric.cc
            Fjt.cc
            GaussianShell.cc
            IntegralParameters.cc
            IrreducibleRepresentation.cc
            Math.cc
            Matrix.cc
            Molecule.cc
            Output.cc
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
            Lapack.cc
            ${LIBINT_FILES}
            )

##############################################
# Panache library compiler flags, etc
##############################################
# Compiler definitions list need to be split
string(REPLACE ";" " " PANACHE_CXX_FLAGS "${PANACHE_CXX_FLAGS}")
set_target_properties(panache PROPERTIES COMPILE_FLAGS ${PANACHE_CXX_FLAGS})
set_target_properties(panache PROPERTIES INCLUDE_DIRECTORIES "${LIBINT_INCLUDEDIR}")


# For now, install all header files
file(GLOB HFILES "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
install(FILES ${HFILES} DESTINATION include/panache)
install(TARGETS panache ARCHIVE DESTINATION lib)
