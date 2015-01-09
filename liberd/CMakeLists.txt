file(GLOB LIBERD_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/*.f")

add_library(panache_erd 
            ${LIBERD_SOURCE}
            )

##############################################
# Panache library compiler flags, etc
##############################################
# Compiler definitions list need to be split
string(REPLACE ";" " " LIBERD_F90_FLAGS "${LIBERD_F90_FLAGS}")
set_target_properties(panache_erd PROPERTIES COMPILE_FLAGS ${LIBERD_F90_FLAGS})


# For now, install all header files
install(TARGETS panache_erd ARCHIVE DESTINATION lib)
