# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -Wall -Wextra -std=c++11 

# Rule to compile main.cpp
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to compile test.cpp
test.o: test.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to link the executable for main.cpp
main: main.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule to link the executable for test.cpp
test: test.o
	$(CXX) $(CXXFLAGS) -o $@ $^

# Default rule
all: main test

# Rule to clean the directory
clean:
	rm -f main.o test.o main test