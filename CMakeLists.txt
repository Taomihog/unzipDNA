cmake_minimum_required(VERSION 3.10)

project(unzip)

# Set C++ standard (change to your preferred version)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add source files
file(GLOB SOURCES "src/*.cpp")
file(GLOB HEADERS "include/*.h")

# Create an executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Add include directories
target_include_directories(${PROJECT_NAME} PRIVATE "include")

# Optionally, link libraries if needed
# target_link_libraries(${PROJECT_NAME} PRIVATE YourLibraryName)

# Set additional compiler flags if needed
# target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra)