CXX = g++
#CXXVER = -std=c++17 #Handled by ROOT
CXXWAR = -Wall -Wextra -Werror


INCLUDES = -I external -I `root-config --incdir`/include


all: main.cpp
	$(CXX) $(CXXVER) $(CXXWAR) $(INCLUDES) -o main main.cpp -O3 `root-config --glibs --cflags --libs`

taskA1:: taskA1.cpp
	$(CXX) $(CXXVER) $(CXXWAR) $(INCLUDES) -o taskA1 taskA1.cpp -O3 `root-config --glibs --cflags --libs`

taskA2:: taskA2.cpp
	$(CXX) $(CXXVER) $(CXXWAR) $(INCLUDES) -o taskA2 taskA2.cpp -O3 `root-config --glibs --cflags --libs`

taskB1:: taskB1.cpp
	$(CXX) $(CXXVER) $(CXXWAR) $(INCLUDES) -o taskB1 taskB1.cpp -O3 `root-config --glibs --cflags --libs`



.PHONY: all




