cmake_minimum_required(VERSION 3.10)

# Set the target architecture.
# All modern x86/x64 processors support AVX2.
# Older x86/x64 processors may support SSE2 but not AVX2.
# Very old x86/x64 processors, or non x86/x64 
# processors, do not support any of the two.
set(ENABLE_SSE2 True)
set(ENABLE_AVX2 True)

#####################

# use gcc and g++ instead of clang
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
set(CMAKE_CXX_COMPILER "g++-10")
endif()

# set the project name
project(mesh_generator)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set (SOURCES
    src/makePolyhedralMesh.cpp
    src/delaunay.cpp
    src/conforming_mesh.cpp
    src/extended_predicates.cpp
    src/BSP.cpp
    src/inOutPartition.cpp
    src/winMem.cpp
)

# add the executable
add_executable(${PROJECT_NAME} 
    src/main.cpp
    ${SOURCES}
)

# make a static library for use within external programs
add_library(${PROJECT_NAME}_lib 
    ${SOURCES}
)

set(ALL_TARGETS ${PROJECT_NAME} ${PROJECT_NAME}_lib)

foreach (TARGET ${ALL_TARGETS})
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
# grant IEEE 754 compliance
	target_compile_options(${TARGET} PUBLIC "/fp:strict")
# use intrinsic functions
	target_compile_options(${TARGET} PUBLIC "/Oi")
# enable AVX2 instruction set
	target_compile_options(${TARGET} PUBLIC "/arch:AVX2")
# turn off annoying warnings
	target_compile_options(${TARGET} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:AVX2")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "/arch:SSE2")
	endif()
else()
# set standard optimization level
	target_compile_options(${TARGET} PUBLIC -O2)
# grant IEEE 754 compliance
	target_compile_options(${TARGET} PUBLIC -frounding-math)
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${PROJECT_NAME} PUBLIC "-mavx2")
		target_compile_options(${PROJECT_NAME} PUBLIC "-mfma")
	elseif(ENABLE_SSE2)
		target_compile_options(${PROJECT_NAME} PUBLIC "-msse2")
	endif()
endif()

# Public include directory
target_include_directories(${TARGET} PUBLIC
	src
)
endforeach(TARGET)

