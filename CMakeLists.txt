cmake_minimum_required(VERSION 3.0.0)
project(SwirlBirds VERSION 0.1.0)

include(CTest)
enable_testing()
set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

include(CPack)


add_subdirectory(include)
find_package(SFML COMPONENTS System Window Graphics REQUIRED)

add_executable(SwirlBirds main.cpp)
TARGET_LINK_LIBRARIES(SwirlBirds
  gsl
  -lgslcblas
  sfml-window
  sfml-graphics
  sfml-system
  Pvector
)
