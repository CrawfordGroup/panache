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


install(TARGETS runtest RUNTIME DESTINATION bin)

