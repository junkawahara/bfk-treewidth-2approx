CXX ?= g++
CXXFLAGS ?= -O2 -march=native -std=c++17 -Wall -Wextra

all: bin/tw2approx bin/verify_td bin/gen

bin/%: src/%.cpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $<

debug: CXXFLAGS = -O1 -g -std=c++17 -Wall -fsanitize=address,undefined
debug: clean all

clean:
	rm -rf bin

.PHONY: all clean debug
