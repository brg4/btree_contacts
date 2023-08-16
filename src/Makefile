# Makefile for program contact_calc
F90	= ${FC}

OPTIONS	= -fdump-core -fcheck=all -Wunused -Wuninitialized -fbacktrace -fopenmp -lgomp -lgfortran -O3

OBJECTS = params.o contact_object.o contact_procedures.o read_write.o contact_calc.o
PROGRAM	= contact_calc 	#target- executable that is generated
ARCHIVE = contact_calc

help:
	@echo Try:
	@echo make $(PROGRAM) .... to build the program named $(PROGRAM)
	@echo make clean .... to clean up, removing object files and program $(PROGRAM)
	@echo make listing .... to make a printable listing ... you may want to adjust the page length
	@echo make archive .... to make an archive tar file you can transfer or submit


%.o:	%.f90
	$(F90) $(OPTIONS) -c $<

$(PROGRAM):	$(OBJECTS)	
	$(F90) $(OBJECTS) $(OPTIONS) -o $(PROGRAM)	

clean:
	rm $(OBJECTS) $(PROGRAM)
