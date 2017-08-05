CXX:=g++
CXXFLAGS:=-std=c++11 -O3
CPPFLAGS:= -lboost_serialization
LIBS:=
INCLUDE:= -I../../include -I/usr/include/
DEPS:=$(wildcard ../../include/*.h)
SRC_FILES:=$(wildcard *.cpp)
TRG_FILES:=$(patsubst %.cpp,%,$(SRC_FILES))

all: $(TRG_FILES)

%: %.cpp $(DEPS)
	$(CXX) -o $@ $(CXXFLAGS) $(CPPFLAGS) $(INCLUDE) $< $(LIBS)

.PHONY: clean
clean:
	rm -rf $(TRG_FILES)
