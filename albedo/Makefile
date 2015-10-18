
RM       := rm -f
CXX      := g++
CXXFLAGS := -Wall -Wextra -std=c++11

SRC = $(wildcard *.cpp)
OBJ = $(patsubst %.cpp,%.o,$(SRC))
EXE = testsalbedo

ifdef EIGEN
 CXXFLAGS += -I$(EIGEN)
endif


default: $(EXE) 

clean:
	$(RM) $(EXE) $(OBJ)

.PHONY: clean

$(EXE): $(OBJ)
	$(CXX) -o $@ $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

