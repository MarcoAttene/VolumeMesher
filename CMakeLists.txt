cmake_minimum_required(VERSION 3.10)

# use gcc and g++ instead of clang
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
set(CMAKE_CXX_COMPILER "g++-10")
endif()

# set the project name
project(mesh_generator)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set (SOURCES
    src/makePolyhedralMesh.cpp
    src/delaunay.cpp
    src/conforming_mesh.cpp
    src/extended_predicates.cpp
    src/BSP.cpp
    src/inOutPartition.cpp
    Indirect_Predicates/implicit_point.cpp
    Indirect_Predicates/numerics.cpp
    Indirect_Predicates/predicates/hand_optimized_predicates.cpp
    Indirect_Predicates/predicates/indirect_predicates.cpp
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
# reserve enough stack size
	target_link_options(${TARGET} PUBLIC "/STACK:8421376")
# turn off annoying warnings
	target_compile_options(${TARGET} PUBLIC "/D _CRT_SECURE_NO_WARNINGS")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
# set standard optimization level
	target_compile_options(${TARGET} PUBLIC -O2)
# grant IEEE 754 compliance
	target_compile_options(${TARGET} PUBLIC -frounding-math)
# reserve enough stack size	
	target_compile_options(${TARGET} PUBLIC -Wl,-z,stacksize=8421376)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#Disable all optimizations
	target_compile_options(${TARGET} PUBLIC -O0)
# reserve enough stack size	
	target_compile_options(${TARGET} PUBLIC -Wl,-z,stacksize=8421376)
endif()

# Public include directory
target_include_directories(${TARGET} PUBLIC
	src
	Indirect_Predicates
	Indirect_Predicates/predicates
)
endforeach(TARGET)

