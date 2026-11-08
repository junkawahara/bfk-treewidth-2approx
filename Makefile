CXX ?= g++
CXXFLAGS ?= -O2 -march=native -std=c++17 -Wall -Wextra

all: bin/tw2approx bin/verify_td bin/gen

bin/%: src/%.cpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $<

DEBUG_CXXFLAGS = -O1 -g -std=c++17 -Wall -fsanitize=address,undefined

# clean, then build: two sub-makes run one after the other, so "make -j debug" cannot
# run "clean" while (or after) "all" builds
debug:
	$(MAKE) clean
	$(MAKE) all CXXFLAGS="$(DEBUG_CXXFLAGS)"

clean:
	rm -rf bin

.PHONY: all clean debug
