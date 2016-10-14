#CMakeList for src dir


#Add source headers
set(SOURCE_FILES include/stkITKAdapter.hxx)
#Add headers
set(INCLUDE_FILES include/stkITKAdapter.h)
message( STATUS "Source files: ${SOURCE_FILES}" )
#Add library
add_library( stkITKAdapter SHARED ${SOURCE_FILES} ${INCLUDE_FILES} )
set(ITK_DIR /software/ITK/build/)
find_package(ITK REQUIRED)

set_target_properties(stkITKAdapter PROPERTIES LINKER_LANGUAGE CXX)

