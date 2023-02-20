CXX = g++
#CXXVER = -std=c++17 #Handled by ROOT
CXXWAR = -Wall -Wextra -Werror


INCLUDES = -I external -I `root-config --incdir`/include


all: main.cpp
	echo $(INCLUDES)
	$(CXX) $(CXXVER) $(CXXWAR) $(INCLUDES) -o main main.cpp -O3 `root-config --glibs --cflags --libs`




.PHONY: all




