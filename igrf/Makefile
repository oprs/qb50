
OBJ = IGRF.o IGRF-mag.o

RM       := rm -f
CXX      := g++
CXXFLAGS := -Wall -Wextra -std=c++11

ifdef EIGEN
 CXXFLAGS += -I$(EIGEN)
endif

default: magtest

magtest: $(OBJ) magtest.o
	$(CXX) -o $@ $^

clean:
	$(RM) magtest magtest.o $(OBJ)

.PHONY: default clean


%.o: %.cpp
	$(QUIET) $(CXX) $(CXXFLAGS) -o $@ -c $<

