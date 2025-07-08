cmake_minimum_required(VERSION 3.29)
project(sutSpice2)

set(CMAKE_CXX_STANDARD 20)

add_executable(sutSpice2 main.cpp)
target_link_libraries(sutSpice2 -lmingw32 -lSDL2main -lSDL2 -lSDL2_image -lSDL2_ttf -lSDL2_mixer)