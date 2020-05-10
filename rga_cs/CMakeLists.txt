cmake_minimum_required(VERSION 3.5)
project(cross_section)


set(CMAKE_CXX_STANDARD 14)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net) #RooFit RooFitCore RooStats RootAuth
find_package(ROOT COMPONENTS MathMore)
include(${ROOT_USE_FILE})

include_directories(include)

file(GLOB cross_EXE
    "exe/*.C"
)

foreach(exe_file IN LISTS cross_EXE)
    get_filename_component(exe_name ${exe_file} NAME_WE)
    add_executable(${exe_name} ${exe_file})
    target_link_libraries(${exe_name} ${ROOT_LIBRARIES})
endforeach()

