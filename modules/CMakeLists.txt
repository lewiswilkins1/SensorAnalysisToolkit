#Top level CMakeList for project

#cmake_minimum_required (VERSION 2.8)

#set(ROOT_DIR /panfs/panasas01/phys/phrfp/Software/root-6.0.4/cmake)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} /panfs/panasas01/phys/phrfp/Software/root-6.0.4/cmake/modules)

#project (SensorAnalysisToolKit)

#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
#set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )

add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/IO/RawIO)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/IO/RawIOImageStack)

add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/Base/Image)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/Base/ImageStack)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/Base/ITKAdapter)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/Base/DarkFrameDetection)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/ImageCalculations/ImageSum)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/ImageCalculations/ImageDivision)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/ImageCalculations/ImageSubtract)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/ImageCalculations/ImageSquare)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/ImageCalculations/ImageVariance)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/ImageCalculations/ImageMask)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/ImageCalculations/ImageFlatFieldCorrection)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/ImageCalculations/ImageBadPixelAverage)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/ImageCalculations/ImageMinus)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/ImageCalculations/ImageResize)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/ImageCalculations/ITKTools)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/Plotting/Histograms)
add_subdirectory ($ENV{HOME}/SensorAnalysisToolkit/modules/Plotting/HistogramAnalysis)