AF_DIR  := $(HOME)/git/afivo
LIBDIRS := $(AF_DIR)/external_libraries/silo/lib \
$(AF_DIR)/external_libraries/hypre/lib $(HOME)/opt/lib64 $(AF_DIR)/lib_2d
INCDIRS := $(HOME)/opt/include/fortran_stdlib/GNU-13.2.1 $(AF_DIR)/lib_2d
LIBS	:= afivo silo HYPRE fortran_stdlib

include $(AF_DIR)/src/makerules.make

.PHONY: all clean lib

PROGS := poisson_solver test_m_solver

all:	$(PROGS)

clean:
	$(RM) $(PROGS) *.o *.mod

lib:	m_solver.o
	f2py -c -m poisson m_solver.f90 -I$(HOME)/opt/include/fortran_stdlib/GNU-13.2.1 \
	-I$(HOME)/git/afivo/lib_2d -L$(HOME)/git/afivo/lib_2d -lafivo m_solver_lib.o \
	-lHYPRE -lfortran_stdlib -lsilo -Llibs/lib -L$(HOME)/opt/lib64

$(PROGS): $(AF_DIR)/lib_2d/libafivo.a
m_solver.o: m_solver_lib.o
test_m_solver: m_solver.o m_solver_lib.o

