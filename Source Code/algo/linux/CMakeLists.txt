PROJECT(detector)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++0x -W -Wall -Wno-sign-compare -O2 -s -pipe -mmmx -msse -msse2 -msse3")

SET(SOURCES
  main.cpp
)

set(Boost_USE_STATIC_LIBS   ON)
find_package( Boost REQUIRED COMPONENTS serialization)

ADD_EXECUTABLE(detector ${SOURCES})
TARGET_LINK_LIBRARIES(detector ${Boost_LIBRARIES})
