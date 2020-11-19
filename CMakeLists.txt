cmake_minimum_required(VERSION 3.1)
project(final_project)

add_subdirectory(src/ext)

# include directories
include_directories(${PROJECT_SOURCE_DIR}/src
                    ${PROJECT_SOURCE_DIR}/src/ext
                    ${PROJECT_SOURCE_DIR}/src/ext/stb
)

set(ext_stb_srcs
    src/ext/stb/stb_image_write.h
    src/ext/stb/stb_image.h
)

add_definitions(-DDATA_DIR=R\"\(${PROJECT_SOURCE_DIR}/data\)\")

set(common_srcs
    src/array3D.h
    src/exceptions.h
    src/floatimage.cpp
    src/floatimage.h
    src/utils.h
    src/array.h
)

# Compile common library
add_library(common_lib STATIC
    ${common_srcs}
    ${ext_stb_srcs}
)
SOURCE_GROUP("src\\ext\\stb" FILES ${ext_stb_srcs})

set(final_project_srcs    src/final_project_main.cpp src/hdr.cpp src/filtering.cpp src/a2.cpp src/panoramicTrans.cpp)    # final_project
add_executable(final_project ${final_project_srcs})                                               # final_project
target_link_libraries(final_project common_lib)                                        # final_project
SOURCE_GROUP("src" FILES ${final_project_srcs})                                        # final_project
                                                                            # final_project

# Create the output data directory
add_custom_command(TARGET final_project POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory "${PROJECT_SOURCE_DIR}/data/output") # a5