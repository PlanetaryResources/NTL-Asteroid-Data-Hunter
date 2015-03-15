PROJECT(makedata)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++98 -W -Wall -Wno-sign-compare -O3 -s -pipe")

SET(SOURCES
  main.cpp
)

ADD_EXECUTABLE(makedata ${SOURCES})
TARGET_LINK_LIBRARIES(makedata)
