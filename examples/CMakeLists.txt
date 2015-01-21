#############################
# Runtest C/C++ target
#############################


# see http://www.cmake.org/Wiki/CMake/Tutorials/Object_Library
add_library(commoncpp OBJECT common_cpp.cc)

add_executable(example_cpp example_cpp.cc $<TARGET_OBJECTS:commoncpp>)
add_executable(example_cpp_iterators example_cpp_iterators.cc $<TARGET_OBJECTS:commoncpp>)

if(RUNTEST_LINK_LIBRARIES)
    target_link_libraries(example_cpp panache ${RUNTEST_LINK_LIBRARIES})
    target_link_libraries(example_cpp_iterators panache ${RUNTEST_LINK_LIBRARIES})
endif(RUNTEST_LINK_LIBRARIES)

if(RUNTEST_CXX_FLAGS)
    string(REPLACE ";" " " RUNTEST_CXX_FLAGS "${RUNTEST_CXX_FLAGS}")
    set_target_properties(commoncpp PROPERTIES COMPILE_FLAGS ${RUNTEST_CXX_FLAGS})
    set_target_properties(example_cpp PROPERTIES COMPILE_FLAGS ${RUNTEST_CXX_FLAGS})
    set_target_properties(example_cpp_iterators PROPERTIES COMPILE_FLAGS ${RUNTEST_CXX_FLAGS})
endif(RUNTEST_CXX_FLAGS)

set_target_properties(commoncpp PROPERTIES INCLUDE_DIRECTORIES "${RUNTEST_CXX_INCLUDES}")
set_target_properties(example_cpp PROPERTIES INCLUDE_DIRECTORIES "${RUNTEST_CXX_INCLUDES}")
set_target_properties(example_cpp_iterators PROPERTIES INCLUDE_DIRECTORIES "${RUNTEST_CXX_INCLUDES}")

if(RUNTEST_CXX_LINK_FLAGS)
    string(REPLACE ";" " " RUNTEST_CXX_LINK_FLAGS "${RUNTEST_CXX_LINK_FLAGS}")
    set_target_properties(example_cpp PROPERTIES LINK_FLAGS ${RUNTEST_CXX_LINK_FLAGS})
    set_target_properties(example_cpp_iterators PROPERTIES LINK_FLAGS ${RUNTEST_CXX_LINK_FLAGS})
endif(RUNTEST_CXX_LINK_FLAGS)


install(TARGETS example_cpp RUNTIME DESTINATION bin)

