#############################
# Runtest C/C++ target
#############################


# see http://www.cmake.org/Wiki/CMake/Tutorials/Object_Library
add_library(commoncpp OBJECT common_cpp.cc)

if(PANACHE_FORTRAN)
    add_library(commonf90 OBJECT common_f90.f90)
endif(PANACHE_FORTRAN)

add_executable(example_cpp example_cpp.cc $<TARGET_OBJECTS:commoncpp>)
add_executable(example_cpp_iterators example_cpp_iterators.cc $<TARGET_OBJECTS:commoncpp>)

if(PANACHE_FORTRAN)
  add_executable(example_f90 example_f90.f90 $<TARGET_OBJECTS:commonf90>)
endif(PANACHE_FORTRAN)


if(RUNTEST_LINK_LIBRARIES)
    target_link_libraries(example_cpp panache ${RUNTEST_LINK_LIBRARIES})
    target_link_libraries(example_cpp_iterators panache ${RUNTEST_LINK_LIBRARIES})
endif(RUNTEST_LINK_LIBRARIES)

if(PANACHE_FORTRAN AND FRUNTEST_LINK_LIBRARIES)
    target_link_libraries(example_f90 panache ${FRUNTEST_LINK_LIBRARIES})
endif(PANACHE_FORTRAN AND FRUNTEST_LINK_LIBRARIES)

if(RUNTEST_CXX_FLAGS)
    string(REPLACE ";" " " RUNTEST_CXX_FLAGS "${RUNTEST_CXX_FLAGS}")
    set_target_properties(commoncpp PROPERTIES COMPILE_FLAGS ${RUNTEST_CXX_FLAGS})
    set_target_properties(example_cpp PROPERTIES COMPILE_FLAGS ${RUNTEST_CXX_FLAGS})
    set_target_properties(example_cpp_iterators PROPERTIES COMPILE_FLAGS ${RUNTEST_CXX_FLAGS})
endif(RUNTEST_CXX_FLAGS)

if(PANACHE_FORTRAN AND FRUNTEST_F90_FLAGS)
  set_target_properties(example_f90 PROPERTIES COMPILE_FLAGS ${FRUNTEST_F90_FLAGS})
endif(PANACHE_FORTRAN AND FRUNTEST_F90_FLAGS)

set_target_properties(commoncpp PROPERTIES INCLUDE_DIRECTORIES "${RUNTEST_CXX_INCLUDES}")
set_target_properties(example_cpp PROPERTIES INCLUDE_DIRECTORIES "${RUNTEST_CXX_INCLUDES}")
set_target_properties(example_cpp_iterators PROPERTIES INCLUDE_DIRECTORIES "${RUNTEST_CXX_INCLUDES}")

if(RUNTEST_CXX_LINK_FLAGS)
    string(REPLACE ";" " " RUNTEST_CXX_LINK_FLAGS "${RUNTEST_CXX_LINK_FLAGS}")
    set_target_properties(example_cpp PROPERTIES LINK_FLAGS ${RUNTEST_CXX_LINK_FLAGS})
    set_target_properties(example_cpp_iterators PROPERTIES LINK_FLAGS ${RUNTEST_CXX_LINK_FLAGS})
endif(RUNTEST_CXX_LINK_FLAGS)

if(PANACHE_FORTRAN AND FRUNTEST_F90_LINK_FLAGS)
    set_target_properties(example_f90 PROPERTIES LINK_FLAGS ${FRUNTEST_F90_LINK_FLAGS})
endif(PANACHE_FORTRAN AND FRUNTEST_F90_LINK_FLAGS)


install(TARGETS example_cpp example_cpp_iterators DESTINATION bin)

if(PANACHE_FORTRAN)
  install(TARGETS example_f90 DESTINATION bin)
endif(PANACHE_FORTRAN)


