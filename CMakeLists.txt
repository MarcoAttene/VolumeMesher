cmake_minimum_required(VERSION 3.10)
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

# Set the target architecture.
# All modern x86/x64 processors support AVX2.
# Older x86/x64 processors may support SSE2 but not AVX2.
# Very old x86/x64 processors, or non x86/x64 
# processors, do not support any of the two.
# Modern ARM processors support NEON.

# --- Uncomment max one of the following lines ---

#							#Generic build
#set(ENABLE_SSE2 True)		#Build for x86-64 with SSE2
set(ENABLE_AVX2 True)		#Build for x86-64 with AVX2
#set(ENABLE_ARMNEON True)	#Build for ARM with NEON

# --------- DON NOT EDIT BELOW THIS LINE! --------

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# disable SSE and AVX on MacOS
if(APPLE)
set(ENABLE_SSE2 False)
set(ENABLE_AVX2 False)
endif()

# AVX includes SSE -> make that explicit
if(ENABLE_AVX2)
set(ENABLE_SSE2 True)
endif()

# Check that no incompatible architectures are slelected
if(ENABLE_ARMNEON AND ENABLE_SSE2)
message(FATAL_ERROR "ARMNEON and SSE2/AVX2 cannot be both set!")
endif()

# set the project name
project(mesh_generator)

set (SOURCES
    src/makePolyhedralMesh.cpp
    src/conforming_mesh.cpp
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

# Dependencies
include(FetchContent)

if(EXISTS ../.local.cmake)
	set(nfg_SOURCE_DIR "../NFG/")
	set(indirect_predicates_SOURCE_DIR "../Indirect_Predicates/")
	set(delaunay3d_SOURCE_DIR "../Delaunay3D/")
else()
	FetchContent_Declare(nfg
		GIT_REPOSITORY https://github.com/MarcoAttene/nfg.git
		GIT_TAG main
		SOURCE_SUBDIR do-not-configure
	)

	FetchContent_Declare(indirect_predicates
		GIT_REPOSITORY https://github.com/MarcoAttene/Indirect_Predicates.git
		SOURCE_SUBDIR do-not-configure
	)

	FetchContent_Declare(Delaunay3D
		GIT_REPOSITORY https://github.com/MarcoAttene/Delaunay3D.git
		GIT_TAG main
		SOURCE_SUBDIR do-not-configure
	)

	FetchContent_MakeAvailable(nfg indirect_predicates delaunay3d)
	message(STATUS "NFG is available in " ${nfg_SOURCE_DIR})
	message(STATUS "Indirect_Predicates is available in " ${indirect_predicates_SOURCE_DIR})
	message(STATUS "Delaunay3D is available in " ${delaunay3d_SOURCE_DIR})
endif()

if(ENABLE_ARMNEON)
	FetchContent_Declare(simde
		GIT_REPOSITORY https://github.com/simd-everywhere/simde.git
		GIT_TAG v0.8.2
	)
	FetchContent_MakeAvailable(simde)
	message(STATUS "simde is available in " ${simde_SOURCE_DIR})

	target_include_directories(${PROJECT_NAME} PUBLIC ${simde_SOURCE_DIR}/simde)
	target_compile_options(${PROJECT_NAME} PUBLIC -DSIMDE_ENABLE_NATIVE_ALIASES)
endif()

# Compiler-specific options

set(ALL_TARGETS ${PROJECT_NAME} ${PROJECT_NAME}_lib)

foreach (TARGET ${ALL_TARGETS})
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	# grant IEEE 754 compliance; remove security controls; turn off annoying warnings
	target_compile_options(${TARGET} PUBLIC "/fp:strict;/GS-;/D _CRT_SECURE_NO_WARNINGS")
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${TARGET} PUBLIC "/arch:AVX2")
	elseif(ENABLE_SSE2)
		target_compile_options(${TARGET} PUBLIC "/arch:SSE2")
	elseif(ENABLE_ARMNEON)
		target_compile_options(${TARGET} PUBLIC "/arch:ARMV8.0")
	endif()
else()
	# set max optimization
	target_compile_options(${TARGET} PUBLIC -O3)
	# grant IEEE 754 compliance
	target_compile_options(${TARGET} PUBLIC -frounding-math)
	# the exact predicates rely on error-free transformations (e.g. Dekker
	# two-product): forbid the compiler from contracting a*b+c into an FMA, which
	# rounds once instead of twice and silently breaks those identities. This
	# mirrors the /fp:strict used on the MSVC path.
	target_compile_options(${TARGET} PUBLIC -ffp-contract=off)
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(${TARGET} PUBLIC "-mavx2")
		target_compile_options(${TARGET} PUBLIC "-mfma")
	elseif(ENABLE_SSE2)
		target_compile_options(${TARGET} PUBLIC "-msse2")
	elseif(ENABLE_ARMNEON)
		target_compile_options(${TARGET} PUBLIC "-march=armv8-a+simd")
	endif()
endif()

# Public include directory
target_include_directories(${TARGET} PUBLIC
	src
	${nfg_SOURCE_DIR}/include
	${indirect_predicates_SOURCE_DIR}/include
	${delaunay3d_SOURCE_DIR}/include
)
endforeach(TARGET)

