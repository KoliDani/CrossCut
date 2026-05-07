cmake_minimum_required(VERSION 3.15)

project(CrossCut
    VERSION 1.0.0
    LANGUAGES CXX
)

# Require C++17 for std::shared_ptr etc.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build shared library (DLL / SO / DYLIB)
add_library(CrossCut SHARED
    CrossCut.cpp
)

# Find JNI
find_package(JNI REQUIRED)
include_directories(${JNI_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Collect your source files
# Not implemented in this project
# set(SOURCES
#     jni/native/CrossCut.cpp
# )

# Export all symbols automatically on Windows (simple + safe)
if (WIN32)
    set_target_properties(CrossCut PROPERTIES
        WINDOWS_EXPORT_ALL_SYMBOLS ON
    )
endif()

# Public include paths
target_include_directories(CrossCut PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

# Require C++17
target_compile_features(CrossCut PUBLIC cxx_std_17)

# Optional: define symbol for dllexport if you use an API macro
target_compile_definitions(CrossCut PRIVATE POLYMESH_BUILD)
# List headers for IDE visibility
target_sources(CrossCut INTERFACE
    boolean/bool_engine.h
    boolean/bool_options.h
    data_structures/geo.h
    data_structures/mesh.h
    fem/mesh_engine.h
    fem/mesh_options.h
)