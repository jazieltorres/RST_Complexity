cmake_minimum_required(VERSION 3.14)
project(RST)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

include_directories(.)

add_executable(RST
        main.cpp
        )

include_directories(/usr/local/Cellar/ntl/11.3.3/include)
link_directories(/usr/local/Cellar/ntl/11.3.3/lib)

include_directories(/usr/local/Cellar/blitz/HEAD-0e373f7/include)
link_directories(/usr/local/Cellar/blitz/HEAD-0e373f7/lib)


file(GLOB LIBRARIES "/usr/local/Cellar/ntl/11.3.3/lib/*.dylib")
message("LIBRARIES = ${LIBRARIES}")

TARGET_LINK_LIBRARIES(RST  ${LIBRARIES})