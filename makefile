CXX:=g++
CXXFLAGS:=-std=c++11 -O3
CPPFLAGS:=
LIBS:= -lboost_serialization
INCLUDE:= -I./include
DEPS:=$(wildcard include/*.h)
SRC_FILES:=$(wildcard src/*.cpp)
TRG_FILES=src/triple-deck

all: $(TRG_FILES)

%: %.cpp $(DEPS)
	$(CXX) -o $@ $(CXXFLAGS) $(CPPFLAGS) $(INCLUDE) $< $(LIBS)

.PHONY: clean
clean:
	rm -rf $(TRG_FILES)
