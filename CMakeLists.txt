cmake_minimum_required(VERSION 3.10)

project(nataf_gsa_project)


# sy- add mpi
# find_package(OpenMP)
 find_package(MPI REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

set(SOURCES main.cpp
	    ERADist.cpp
	    exponentialDist.cpp
	    normalDist.cpp
	    ERANataf.cpp
	    gammaDist.cpp
	    runGSA.cpp
	    runForward.cpp
	    RVDist.cpp
	    gumbelDist.cpp
	    truncExponentialDist.cpp
	    betaDist.cpp
	    jsonInput.cpp
	    uniformDist.cpp
	    chiSquaredDist.cpp
	    lognormalDIst.cpp
	    weibullDist.cpp
	    discreteDist.cpp
	    writeErrors.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib_armadillo/armadillo-10.1.0/include)
include_directories(${MPI_INCLUDE_PATH})

add_executable(nataf_gsa ${SOURCES})
target_link_libraries(nataf_gsa
	CONAN_PKG::boost
	CONAN_PKG::jsonformoderncpp	
	CONAN_PKG::nlopt
	CONAN_PKG::eigen
	CONAN_PKG::mkl-static
	${MPI_LIBRARIES}
	)


if(MPI_COMPILE_FLAGS)
  set_target_properties(nataf_gsa PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(nataf_gsa PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

#target_link_libraries(nataf_gsa
#	CONAN_PKG::boost
#	CONAN_PKG::jsonformoderncpp	
#	CONAN_PKG::nlopt
#	CONAN_PKG::mkl-static
#	CONAN_PKG::eigen)

install(TARGETS nataf_gsa DESTINATION bin)

#add_definitions(-D_ITERATOR_DEBUG_LEVEL=0) 
add_definitions(-D_ITERATOR_DEBUG_LEVEL=2) 

