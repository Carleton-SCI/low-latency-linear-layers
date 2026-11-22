CXX ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17

all: bin/dbs

bin/dbs: src/dbs.cpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $<

test: bin/dbs
	python3 -m unittest discover -s tests -v

verify:
	python3 paper/checks/verify_results.py

clean:
	rm -rf bin

.PHONY: all test verify clean
