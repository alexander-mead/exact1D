# Makefile for 1D simualtion tools

FC = gfortran

#FFLAGS = -Warray-bounds -ffree-line-length-none -fmax-errors=4 -ffpe-trap=invalid,zero,overflow -fimplicit-none -O3 -L/usr/local/lib -lfftw3

FFLAGS = -Warray-bounds -ffree-line-length-none -fmax-errors=4 -ffpe-trap=invalid,zero,overflow -fimplicit-none -O3

BUILD = build
SRC = src
LIB = /Users/Mead/Physics/library
EXEC = exact1D.e
OBJS = $(BUILD)/fix_polynomial.o $(BUILD)/array_operations.o $(BUILD)/random_numbers.o $(SRC)/exact1D.f90

$(EXEC): $(OBJS)
	$(FC) --version
	$(FC) $(FFLAGS) -I$(BUILD) -L$(BUILD) $^ -o $@

$(BUILD)/random_numbers.o: $(LIB)/random_numbers.f90
	$(FC) $(FFLAGS) -c $^
	mv random_numbers.o $(BUILD)/.
	mv random_numbers.mod $(BUILD)/.

$(BUILD)/array_operations.o: $(LIB)/array_operations.f90
	$(FC) $(FFLAGS) -c $^
	mv array_operations.o $(BUILD)/.
	mv array_operations.mod $(BUILD)/.

$(BUILD)/fix_polynomial.o: $(LIB)/fix_polynomial.f90
	$(FC) $(FFLAGS) -c $^
	mv fix_polynomial.o $(BUILD)/.
	mv fix_polynomial.mod $(BUILD)/.

clean:
	rm -f $(BUILD)/*.o
	rm -f $(BUILD)/*.mod
	rm -f $(EXEC)
