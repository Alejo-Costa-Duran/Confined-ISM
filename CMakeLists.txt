cmake_minimum_required(VERSION 3.0.0)
project(SwirlBirds VERSION 0.1.0)

include(CTest)
enable_testing()
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
set(SFML_STATIC_LIBRARIES TRUE)

include(CPack)
INCLUDE_DIRECTORIES(
    C:/dev/vcpkg/vcpkg/installed/x86-windows/include
)
LINK_DIRECTORIES(
    C:/dev/vcpkg/vcpkg/installed/x86-windows/lib
)

add_subdirectory(include)
find_package(SFML COMPONENTS system window graphics CONFIG REQUIRED)
add_executable(SwirlBirds main.cpp)
TARGET_LINK_LIBRARIES(SwirlBirds
  gsl
  sfml-window
  sfml-graphics
  Pvector
)