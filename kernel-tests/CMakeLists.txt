cmake_minimum_required(VERSION 3.20)
project(kernelTest LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
# set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_CUDA_ARCHITECTURES 89)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(USE_HPC_EXT "Use clang hpc-ext" OFF)

if (USE_HPC_EXT)
    set(CMAKE_CXX_COMPILER /home/zhouxing/Documents/Masterarbeit/llvm-hpc-ext/bin/clang++)
    add_compile_definitions(USE_HPC_EXT)
else()
    set(CMAKE_CXX_COMPILER /home/zhouxing/Documents/Masterarbeit/llvm-extended/bin/clang++)
    # set(CMAKE_CXX_COMPILER clang++)
    add_compile_definitions(USE_SOURCE_TO_SOURCE_TRANSFORM)
endif()

# set(CMAKE_CXX_COMPILER /home/zhouxing/Documents/Masterarbeit/llvm-extended/bin/clang++)
# set(CMAKE_CXX_COMPILER clang++)

set(DIMENSION {{DIMENSION}})

add_compile_definitions(GPUOffloadingOMPPacked)

aux_source_directory(. SRC)


#set(CMAKE_CXX_FLAGS "-std=c++17 -O3 -w -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -Xopenmp-target=nvptx64-nvidia-cuda -march=sm_89 -DDimensions=2 -DPeanoDebug=0 -L/home/zhouxing/Documents/Masterarbeit/Peano/src  -L/home/zhouxing/Documents/Masterarbeit/Peano/src/toolbox/loadbalancing  -L/home/zhouxing/Documents/Masterarbeit/Peano/src/toolbox/blockstructured  -L/home/zhouxing/Documents/Masterarbeit/Peano/src/exahype2 -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -Xopenmp-target=nvptx64-nvidia-cuda -march=sm_89 -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -L../../../..//src -lExaHyPE2Core2d  -lToolboxBlockstructured2d   -lToolboxLoadBalancing2d -lPeano4Core2d    -lTarch")
# set(CMAKE_CXX_FLAGS "-w -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -Xopenmp-target=nvptx64-nvidia-cuda -march=sm_89 -DDimensions=2 -DPeanoDebug=0")
set(CMAKE_CXX_FLAGS "-W -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -Xopenmp-target=nvptx64-nvidia-cuda -march=sm_89 -DDimensions=${DIMENSION} -DPeanoDebug=0")
# set(CMAKE_CXX_FLAGS "-W -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -Xopenmp-target=nvptx64-nvidia-cuda -march=sm_89 -DDimensions=${DIMENSION} -DPeanoDebug=0 -Wall -Wextra -Wpedantic -fsanitize=address")
# set(CMAKE_CXX_FLAGS "-g -W -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda -Xopenmp-target=nvptx64-nvidia-cuda -march=sm_89 -DDimensions=${DIMENSION} -DPeanoDebug=0 -Wall -Wextra -Wpedantic -fsanitize=memory,undefined -fsanitize-memory-track-origins")

add_executable(kernelTestMain ${SRC})
target_include_directories(kernelTestMain PUBLIC . .. ../../src)
target_link_directories(kernelTestMain PUBLIC ../../src ../../src/toolbox/loadbalancing ../../src/toolbox/blockstructured ../../src/exahype2)
if(${DIMENSION} EQUAL 2)
    target_link_libraries(kernelTestMain PUBLIC ExaHyPE2Core2d ToolboxBlockstructured2d ToolboxLoadBalancing2d Peano4Core2d Tarch)
elseif(${DIMENSION} EQUAL 3)
    target_link_libraries(kernelTestMain PUBLIC ExaHyPE2Core3d ToolboxBlockstructured3d ToolboxLoadBalancing3d Peano4Core3d Tarch)
endif()
