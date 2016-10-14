#CMakeList for src dir

find_package(ROOT)

#Add source headers
set(SOURCE_FILES include/stkImageHistogram.hxx)
#Add headers
set(INCLUDE_FILES include/stkImageHistogram.h )
message( STATUS "Source files: ${SOURCE_FILES}" )

set(PROJECT_RAW_IO_SOURCE /Users/lewish/SensorAnalysisToolkit_Eclipse)
set(TEST_LIB_DIR ${PROJECT_RAW_IO_SOURCE}/lib)

set(EXT_PROJECTS_DIR $ENV{HOME}/Software)
set(ROOT_INCLUDE_DIRS ${EXT_PROJECTS_DIR}/root-6.0.4/include)
set(ROOT_LIB_DIR ${EXT_PROJECTS_DIR}/root-6.0.4/lib )
include_directories(${ROOT_INCLUDE_DIRS})

message( STATUS "TEST LIB: ${TEST_LIB_DIR}")
#Add library
add_library( stkImageHistogram SHARED ${SOURCE_FILES} ${INCLUDE_FILES} )
#target_link_libraries(stkImageHistogram ${TEST_LIB_DIR}/libstkRawIO.dylib ${TEST_LIB_DIR}/libstkImage.dylib )

set_target_properties(stkImageHistogram PROPERTIES LINKER_LANGUAGE CXX)


