.phony: all main clean

all: main

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
CXX = g++
else
CXX = g++-4.9
endif

CXX_FLAGS=-I. -std=c++0x -MMD -O2
LD_FLAGS =

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)

-include $(DEPS)

main: $(OBJS)
	$(CXX) -o $@ $^ $(LD_FLAGS)

%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -c -o $@ $<

clean:
	rm -f *.o *.d main
