CXX = g++
CXXVER = -std=c++17
CXXWAR = -Wall -Wextra -Werror
INCLUDES = -I external


all: main.cpp
	$(CXX) $(CXXVER) $(CXXWAR) $(INCLUDES) -o main main.cpp




.PHONY: all