# Makefile for 1D simualtion tools

# NB: I had to move all of the source code, including the libraries, into the src/
# directory. If I take them directly from lib then the USE statements in the modules
# try to find the .mod files in the directory they are in first. This means they get them
# from the library directory, preferentially compared to the build directory. This means
# things do not work if you compile with different precision compared to those .mod files
# that will exist in the library directory.

# Set the Fortran compiler
FC = gfortran

# Set the Fortran flags
FFLAGS = -Warray-bounds \
	-ffree-line-length-none \
	-fmax-errors=4 \
	-ffpe-trap=invalid,zero,overflow \
	-fimplicit-none \
	-fdefault-real-8 \
	-fdefault-double-8 \
	-std=gnu \
	-O3

# Build directory (to keep things neat)
BUILD = build

# Source code directory
SRC = src

# Executable
EXEC = exact1D

# Object files
_OBJS = constants.o \
	fix_polynomial.o \
	array_operations.o \
	random_numbers.o \
	numerology.o \
	file_info.o \
	sorting.o

# Append prefix of the build directory to all object files
OBJS = $(addprefix $(BUILD)/,$(_OBJS))

# Rule to create the executable
$(EXEC): $(OBJS) $(SRC)/exact1D.f90
	@echo
	@$(FC) --version
	$(FC) -o $@ $^ -J$(BUILD) $(LDFLAGS) $(FFLAGS)

# Rule to create the objects
$(BUILD)/%.o: $(SRC)/%.f90
	$(FC) -c -o $@ $< -J$(BUILD) $(LDFLAGS) $(FFLAGS)

# Clean command
clean:
	rm -f $(BUILD)/*.o
	rm -f $(BUILD)/*.mod
	rm -f $(EXEC)
