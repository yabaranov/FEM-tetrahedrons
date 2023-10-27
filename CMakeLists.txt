cmake_minimum_required(VERSION 3.27 FATAL_ERROR)

set(PROJECT_NAME FEM-tetrahedrons)
project(${PROJECT_NAME})

add_executable(${PROJECT_NAME} 
src/main.cpp
src/functions.h

src/grid/Grid.h
src/grid/Grid.cpp
src/grid/TransformCubeToHexagon.h
src/grid/TransformCubeToHexagon.cpp

src/vector/FiniteMatrixVector.h
src/vector/MatrixVector.h
src/vector/MatrixVector.cpp

src/SLAE/SparseSLAE.h
src/SLAE/SparseSLAE.cpp

src/FEM/FEM.h
src/FEM/FEM.cpp
src/FEM/LocalAssembly.h
src/FEM/LocalAssembly.cpp
)

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_17)

include_directories(external/rapidjson/include)