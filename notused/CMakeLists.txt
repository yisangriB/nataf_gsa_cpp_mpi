cmake_minimum_required(VERSION 3.10)

project(nataf_gsa_mpi_project)

if (APPLE)
endif()
# sy- add mpi
# find_package(OpenMP)
find_package(MPI)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0)
#   set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++17")
#   set(CMAKE_CXX17_EXTENSION_COMPILE_OPTION "-std=gnu++17")
# elseif (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.1)
#   set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++1z")
#   set(CMAKE_CXX17_EXTENSION_COMPILE_OPTION "-std=gnu++1z")
# endif()

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
	    lognormalDist.cpp
	    weibullDist.cpp
	    discreteDist.cpp
	    writeErrors.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib_armadillo/armadillo-10.1.0/include ${MPI_INCLUDE_PATH})
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib_armadillo/armadillo-10.1.0/include)




add_executable(nataf_gsa_mpi ${SOURCES})
target_link_libraries(nataf_gsa_mpi
	CONAN_PKG::boost
	CONAN_PKG::jsonformoderncpp	
	CONAN_PKG::nlopt
	CONAN_PKG::eigen
	${MPI_LIBRARIES}
	)

#if(MPI_CXX_COMPILE_FLAGS)
#	set_target_properties(nataf_gsa_mpi PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
#endif()
#if(MPI_CXX_LINK_FLAGS)
#	set_target_properties(nataf_gsa_mpi PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
#endif()
#if(MPI_COMPILE_DEFINITIONS)
#	set_target_properties(nataf_gsa_mpi PROPERTIES
#			COMPILE_FLAGS "${MPI_COMPILE_DEFINITIONS}")
#endif()
#
#if(MPI_LINK_FLAGS)
#	set_target_properties(nataf_gsa_mpi PROPERTIES
#			COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
#			LINK_FLAGS "${MPI_LINK_FLAGS}")
#endif()

#target_link_libraries(nataf_gsa_mpi
#	CONAN_PKG::boost
#	CONAN_PKG::jsonformoderncpp	
#	CONAN_PKG::nlopt
#	CONAN_PKG::mkl-static
#	CONAN_PKG::eigen)

install(TARGETS nataf_gsa_mpi DESTINATION bin)

#add_definitions(-D_ITERATOR_DEBUG_LEVEL=0) 
add_definitions(-D_ITERATOR_DEBUG_LEVEL=0)

