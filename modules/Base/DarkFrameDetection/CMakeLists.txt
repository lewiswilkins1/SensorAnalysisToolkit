#CMakeList for src dir

#Add source headers
set(SOURCE_FILES include/stkDarkFrameDetection.hxx)
#Add headers
set(INCLUDE_FILES include/stkDarkFrameDetection.h)
message( STATUS "Source files: ${SOURCE_FILES}" )
#Add library
add_library( stkDarkFrameDetection SHARED ${SOURCE_FILES} ${INCLUDE_FILES} )

set_target_properties(stkDarkFrameDetection PROPERTIES LINKER_LANGUAGE CXX)
