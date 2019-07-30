####################################################
programs := $(wildcard *.cpp)
objects:=$(patsubst %.cpp,%.o,$(programs))
execname=$(shell pwd | sed 's,^\(.*/\)\?\([^/]*\),\2,')
CXX=g++
OMPFLAGS = -fopenmp -O3

$(execname): $(objects)
	$(CXX) $(OMPFLAGS) -o $(execname) $(objects) 
$(objects): %.o: %.cpp
	$(CXX) -c $< -o $@
.PHONY : clean
clean:
	-rm $(execname) $(objects)  *.cpp~  *.h~  *.plt  
run:
	./$(execname)  
####################################################

