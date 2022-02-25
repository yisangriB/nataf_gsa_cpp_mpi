cmake_minimum_required(VERSION 3.10)

project(nataf_gsa_project)

#if(LINUX) 
	set ( PARALLEL_LIB "MPI")
#else()
#	set ( PARALLEL_LIB "OPENMP")
#endif()


# sy- add mpi
# Find the package
if(PARALLEL_LIB STREQUAL  "MPI") 
	find_package(MPI REQUIRED)
else()
	find_package(OpenMP REQUIRED) # Find the package
endif() 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(LINUX)
	# sy - to check if filesystem is working?
	try_compile(HAS_FS "${CMAKE_BINARY_DIR}/temp" 
	"${CMAKE_SOURCE_DIR}/tmp/has_filesystem.cc" 
				CMAKE_FLAGS -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_STANDARD_REQUIRED=ON
				LINK_LIBRARIES stdc++fs)
	if(HAS_FS)
		message(STATUS "Compiler has filesystem support")
	else()
		message(FATAL_ERROR "Compiler is missing filesystem capabilities")
	endif(HAS_FS)
endif()

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)
# for visual studio, do on commandline: conan install . -s build_type=Debug -if=C:/Users/SimCenter/Sangri/nataf_gsa_cpp/out/build/x64-Debug/
# for for CLION, do on build_release folder

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

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib_armadillo/armadillo-10.1.0/include)

if(PARALLEL_LIB STREQUAL  "MPI") 
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
		add_compile_definitions(MPI)

else()
	add_executable(nataf_gsa ${SOURCES})
	target_link_libraries(nataf_gsa
		CONAN_PKG::boost
		CONAN_PKG::jsonformoderncpp	
		CONAN_PKG::nlopt
		CONAN_PKG::eigen
		CONAN_PKG::mkl-static
		OpenMP::OpenMP_CXX)
	add_compile_definitions(OPENMP)
endif() 

install(TARGETS nataf_gsa DESTINATION bin)
