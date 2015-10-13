
RM       = rm -f
CXX      ?= g++
CXXFLAGS ?= -Wall -Wextra -g -std=c++11

SRC = $(wildcard *.cpp)
OBJ = $(patsubst %.cpp,%.o,$(SRC))


default: xyz

clean:
	$(RM) xyz $(OBJ)

.PHONY: clean


xyz: $(OBJ)
	$(CXX) -o $@ $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

