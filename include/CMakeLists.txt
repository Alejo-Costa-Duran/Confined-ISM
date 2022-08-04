add_library(
    Pvector
    Pvector.h
    ../src/PVector.cpp
    Flock.h
    ../src/Flock.cpp
    Bird.h
    ../src/Bird.cpp
    manager.h
    ../src/manager.cpp
)


TARGET_LINK_LIBRARIES(Pvector
    gsl
    sfml-window
    sfml-graphics
)
target_include_directories(Pvector PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")