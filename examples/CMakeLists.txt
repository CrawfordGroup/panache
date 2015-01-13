#############################
# Runtest C/C++ target
#############################
add_executable(example_cpp example_cpp.cc)

if(RUNTEST_LINK_LIBRARIES)
    target_link_libraries(example_cpp panache ${RUNTEST_LINK_LIBRARIES})
endif(RUNTEST_LINK_LIBRARIES)

if(RUNTEST_CXX_FLAGS)
    string(REPLACE ";" " " RUNTEST_CXX_FLAGS "${RUNTEST_CXX_FLAGS}")
    set_target_properties(example_cpp PROPERTIES COMPILE_FLAGS ${RUNTEST_CXX_FLAGS})
endif(RUNTEST_CXX_FLAGS)

set_target_properties(example_cpp PROPERTIES INCLUDE_DIRECTORIES "${RUNTEST_CXX_INCLUDES}")

if(RUNTEST_CXX_LINK_FLAGS)
    string(REPLACE ";" " " RUNTEST_CXX_LINK_FLAGS "${RUNTEST_CXX_LINK_FLAGS}")
    set_target_properties(example_cpp PROPERTIES LINK_FLAGS ${RUNTEST_CXX_LINK_FLAGS})
endif(RUNTEST_CXX_LINK_FLAGS)


install(TARGETS example_cpp RUNTIME DESTINATION bin)

