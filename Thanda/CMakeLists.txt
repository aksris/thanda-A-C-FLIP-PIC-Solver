project(Thanda)
cmake_minimum_required(VERSION 2.8)
set (NUPARU "nuparu")
# add include and src directories to path
include_directories(
        ${NUPARU}/include
        ${NUPARU}/src
)
# You can add your own libraries here
# Add Nuparu library to path for OSX, linux and windows
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${NUPARU}/lib/osx)
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${NUPARU}/lib/linux /usr/lib64)
elseif(WIN32)
    set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${NUPARU}/lib/win)
endif()

#Add an include directory variable
set(GLFW_INCLUDE_DIR ${NUPARU}/include )
set(GLEW_INCLUDE_DIR ${NUPARU}/include )
set(GLFW_LIBRARY_DIR ${CMAKE_LIBRARY_PATH})
set(GLEW_LIBRARY_DIR ${CMAKE_LIBRARY_PATH})

find_package(OPENGL REQUIRED)
find_package(GLEW)
find_library(GLFW_LIBRARY "glfw3" HINTS ${GLFW_LIBRARY_DIR})
find_library(JSON_CPP "jsoncpp")

set(CORE_LIBS ${GLFW_LIBRARY} ${GLEW_LIBRARY} ${JSON_CPP} ${OPENGL_LIBRARY})

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(COREVIDEO "-framework CoreVideo")
    set(IOKIT "-framework IOKit")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -03 -m64 -msse2 -w")
    set(CORE_LIBS ${CORE_LIBS} ${COREVIDEO} ${IOKIT})
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lX11 -lXxf86vm -lXrandr -lpthread -lXi")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -03 -m64 -msse2 -w")
endif()

#if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
if(MSVC)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2")
    set(COMPILER_FLAGS CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE)
    foreach(COMPILER_FLAG ${COMPILER_FLAGS})
        string (REPLACE "/MD" "/MT" ${COMPILER_FLAG} "${${COMPILER_FLAG}}")
    endforeach()
endif()

set(CORE_SRC
    src/main.cpp

    src/particles.cpp

    src/json_parser.cpp

    src/scene/scene.cpp

    src/common/shader.cpp

    src/common/controls.cpp

    src/common/texture.cpp

    src/common/objloader.cpp

    src/common/vboindexer.cpp

)

add_executable(Thanda ${CORE_SRC})

target_link_libraries(Thanda ${CORE_LIBS})

