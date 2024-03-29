# Makefile for 1D simualtion tools

# NB: I had to move all of the source code, including the libraries, into the src/
# directory. If I take them directly from lib then the USE statements in the modules
# try to find the .mod files in the directory they are in first. This means they get them
# from the library directory, preferentially compared to the build directory. This means
# things do not work if you compile with different precision compared to those .mod files
# that will exist in the library directory.

# Set the Fortran compiler
FC = gfortran

# FFTW library locations
#FFTW_lib = /usr/local/lib
#FFTW_inc = /usr/local/include

# Set the Fortran flags
FFLAGS = \
	-Warray-bounds \
	-ffree-line-length-none \
	-fmax-errors=4 \
	-ffpe-trap=invalid,zero,overflow \
	-fimplicit-none \
	-fdefault-real-8 \
	-fdefault-double-8 \
	-std=gnu \
	-O3

#-L${FFTW_lib} \
#-I${FFTW_inc} \
#-lfftw3 \
#-lfftw3f \
#-lfftw3l

# Debugging flags
DEBUG_FLAGS = \
	-Wall \
	-fcheck=all \
	-fbounds-check \
	-fbacktrace \
	-Og

# Build directory (to keep things neat)
BUILD = build

# Source code directory
SRC = src

# Module directory
MOD = /Users/Mead/Physics/library/src

# Executable
EXEC = exact1D

# Object files
_OBJS = \
	constants.o \
	fix_polynomial.o \
	array_operations.o \
	logical_operations.o \
	random_numbers.o \
	file_info.o \
	sorting.o \
	table_integer.o

# Append prefix of the build directory to all object files
OBJS = $(addprefix $(BUILD)/,$(_OBJS))

# Rule to create the executable
$(EXEC): $(OBJS) $(SRC)/exact1D.f90
	@echo
	@$(FC) --version
	$(FC) -o $@ $^ -J$(BUILD) $(LDFLAGS) $(FFLAGS)

debug: FFLAGS += $(DEBUG_FLAGS)
debug: $(EXEC)

# Rule to create the objects
$(BUILD)/%.o: $(MOD)/%.f90
	$(FC) -c -o $@ $< -J$(BUILD) $(LDFLAGS) $(FFLAGS)

# Clean command
clean:
	rm -f $(BUILD)/*.o
	rm -f $(BUILD)/*.mod
	rm -f $(EXEC)
