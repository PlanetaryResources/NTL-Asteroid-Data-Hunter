BOOST_PATH=/Users/wangjinggang/tc/AST-release/boost/boost_1_56_0/
CXXFLAGS=--std=c++0x -W -Wall -Wno-sign-compare -O2 -pipe -mmmx -msse -msse2 -msse3
all:
	g++ $(CXXFLAGS) -I"$(BOOST_PATH)" -L$(BOOST_PATH)/stage/lib -lboost_serialization -o detector main.cpp
