CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -pedantic

all: nbody

nbody: main.cpp
	$(CXX) $(CXXFLAGS) -o nbody main.cpp

clean:
	rm -f nbody *.o
