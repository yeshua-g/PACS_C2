# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -O3 -Wall -I. -Wno-conversion-null -Wno-deprecated-declarations -I../../pacs-examples/Examples/include -std=c++20

# Linker flags
LDFLAGS = -L../../pacs-examples/Examples/lib -Wl,-rpath=../../pacs-examples/Examples/lib

# Source and object files for main
MAIN_SRCS = main.cpp
MAIN_OBJS = $(MAIN_SRCS:.cpp=.o)

# Source and object files for test
TEST_SRCS = test.cpp
TEST_OBJS = $(TEST_SRCS:.cpp=.o)

# Rule to compile main.cpp
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to compile test.cpp
test.o: test.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to link the executable for main.cpp
main: main.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

# Rule to link the executable for test.cpp
test: test.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

# Rule to clean the object files
clean:
	rm -f $(MAIN_OBJS) $(TEST_OBJS)

# Rule to clean the executables
cleanall:
	rm -f main test
