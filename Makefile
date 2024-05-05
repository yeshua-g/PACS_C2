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

# Source and object files for test2
TEST2_SRCS = test2.cpp
TEST2_OBJS = $(TEST2_SRCS:.cpp=.o)

# Source and object files for extras
EXTRAS_SRCS = extras.cpp
EXTRAS_OBJS = $(EXTRAS_SRCS:.cpp=.o)

# Rule to compile main.cpp
main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to compile test.cpp
test.o: test.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to compile test2.cpp
test2.o: test2.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to compile extras.cpp
extras.o: extras.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to link the executable for main.cpp
main: main.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

# Rule to link the executable for test.cpp
test: test.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

# Rule to link the executable for test2.cpp
test2: test2.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

# Rule to link the executable for extras.cpp
extras: extras.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

# Rule to build all executables
all: main test test2 extras

# Rule to clean the object files and executables
clean: 
	rm -f $(MAIN_OBJS) $(TEST_OBJS) $(TEST2_OBJS) $(EXTRAS_OBJS) main test test2 extras


