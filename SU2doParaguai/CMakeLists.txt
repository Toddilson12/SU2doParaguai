cmake_minimum_required (VERSION 3.8)

# 1. Direct the output to a clean Linux folder
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "$ENV{HOME}/CFD_Binaries")

# 2. Add executable and source files
add_executable (SU2DP 
    src/main.cpp
 "src/mesh/meshParser.cpp" "src/mesh/mesh.cpp" "src/models/FluidModel.cpp" "src/models/IdealPerfGas.cpp" "src/variables/Variables.cpp" "src/variables/EulerVars.cpp" "src/output/Output.cpp" "src/gradients/Gradient.cpp" "src/numerics/FluxNumerics.cpp")

# 3. Link the header folder
target_include_directories(SU2DP PUBLIC include)

# 4. Modern C++ standard
target_compile_features(SU2DP PUBLIC cxx_std_17)

# 5. Release Optimization Flags
if(CMAKE_BUILD_TYPE MATCHES "Release")
    target_compile_options(SU2DP PRIVATE -O3 -march=native -ffast-math -funroll-loops)
    set_property(TARGET SU2DP PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
