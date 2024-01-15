CXX = g++ -O2
CXXFLAGS = -Wall -g

# Define the paths
SRCDIR = src
RUNDIR = run

TARGET = $(RUNDIR)/moldy
OBJS = $(RUNDIR)/main.o $(RUNDIR)/io.o $(RUNDIR)/sys_physics.o $(RUNDIR)/time_integration.o
HEADER = $(SRCDIR)/header.h

# Default target
all: $(TARGET)

# Phony target for moldy
.PHONY: moldy
moldy: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) -lm

# Compile .cpp files to .o files
$(RUNDIR)/main.o: $(SRCDIR)/main.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/main.cpp -o $(RUNDIR)/main.o

$(RUNDIR)/io.o: $(SRCDIR)/io.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/io.cpp -o $(RUNDIR)/io.o

$(RUNDIR)/sys_physics.o: $(SRCDIR)/sys_physics.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/sys_physics.cpp -o $(RUNDIR)/sys_physics.o

$(RUNDIR)/time_integration.o: $(SRCDIR)/time_integration.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) -c $(SRCDIR)/time_integration.cpp -o $(RUNDIR)/time_integration.o

clean:
	rm -f $(RUNDIR)/*.o $(TARGET)