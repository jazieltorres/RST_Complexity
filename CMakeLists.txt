cmake_minimum_required(VERSION 3.14)
project(RST)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

include_directories(.)

add_executable(RST
        Test_1.cpp Functions.cpp Functions.h)

include_directories(/usr/local/Cellar/ntl/11.4.3/include)
link_directories(/usr/local/Cellar/ntl/11.4.3/lib)

include_directories(/usr/local/Cellar/blitz/1.0.2/include)
link_directories(/usr/local/Cellar/blitz/1.0.2/lib)


file(GLOB LIBRARIES "/usr/local/Cellar/ntl/11.4.3/lib/*.dylib")
message("LIBRARIES = ${LIBRARIES}")

TARGET_LINK_LIBRARIES(RST  ${LIBRARIES})

