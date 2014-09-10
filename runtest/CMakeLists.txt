#############################
# Runtest C/C++ target
#############################
add_executable(runtest runtest.cc)

if(RUNTEST_LINK_LIBRARIES)
    target_link_libraries(runtest panache ${RUNTEST_LINK_LIBRARIES})
endif(RUNTEST_LINK_LIBRARIES)

if(RUNTEST_CXX_FLAGS)
    string(REPLACE ";" " " RUNTEST_CXX_FLAGS "${RUNTEST_CXX_FLAGS}")
    set_target_properties(runtest PROPERTIES COMPILE_FLAGS ${RUNTEST_CXX_FLAGS})
endif(RUNTEST_CXX_FLAGS)

set_target_properties(runtest PROPERTIES INCLUDE_DIRECTORIES "${RUNTEST_CXX_INCLUDES}")

if(RUNTEST_CXX_LINK_FLAGS)
    string(REPLACE ";" " " RUNTEST_CXX_LINK_FLAGS "${RUNTEST_CXX_LINK_FLAGS}")
    set_target_properties(runtest PROPERTIES LINK_FLAGS ${RUNTEST_CXX_LINK_FLAGS})
endif(RUNTEST_CXX_LINK_FLAGS)


#############################
# FRuntest F90 target
#############################
#add_executable(fruntest fruntest.f90)
#
#set_target_properties(fruntest PROPERTIES LINKER_LANGUAGE Fortran)
#target_link_libraries(fruntest panache ${FRUNTEST_LINK_LIBRARIES} ${LIBINT_LIBRARY} ${LAPACK_LIBRARIES})
#
#string(REPLACE ";" " " FRUNTEST_F90_FLAGS "${FRUNTEST_F90_FLAGS}")
#set_target_properties(fruntest PROPERTIES COMPILE_FLAGS ${FRUNTEST_F90_FLAGS})
#
#if(OPENMP_FOUND)
#  # BUG - findopenmp doesn't support Fortran. This is a big assumption
#  set_target_properties(fruntest PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS})
#endif(OPENMP_FOUND)
#
#
#install(TARGETS runtest fruntest RUNTIME DESTINATION bin)
install(TARGETS runtest RUNTIME DESTINATION bin)

