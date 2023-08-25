add_library(
    Pvector
    PVector.h
    ../src/PVector.cpp
    Flock.h
    ../src/Flock.cpp
    Bird.h
    ../src/Bird.cpp
    Simulation.h
    ../src/Simulation.cpp
)


TARGET_LINK_LIBRARIES(Pvector
    gsl
    -lgslcblas
    sfml-window
    sfml-graphics
    sfml-system
)
target_include_directories(Pvector PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")