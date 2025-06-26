# Define PETSc environment variables
PETSC_DIR ?= /home/jay/petsc
PETSC_ARCH ?= /home/jay/petsc/arch-linux-c-debug
# Compiler and flags
CXX = mpicxx
CXXFLAGS = -I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include
LDFLAGS = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc

# Source file and executable name
# uncomment the VischydroNode.cpp file to test the fillnode function
SRC = rough.cpp #VischydroNode.cpp
EXE = rough.exe 
DEL = *.exe

# Default target: compile the program
all: $(EXE)

# Rule to compile the executable
$(EXE): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(EXE) $(LDFLAGS)

# Rule to run the program with MPI
run: $(EXE)
	mpiexec -n 1 ./$(EXE)

# Clean up generated files
clean:
	rm $(DEL)
