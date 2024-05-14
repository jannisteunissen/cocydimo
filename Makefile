.PHONY: all clean

all:	m_solver.o
	f2py -c -m poisson m_solver.f90 -I$(HOME)/opt/include/fortran_stdlib/GNU-13.2.1 \
	-I$(HOME)/git/afivo/lib_2d -L$(HOME)/git/afivo/lib_2d -lafivo m_solver_lib.o \
	-lHYPRE -lfortran_stdlib -lsilo -Llibs/lib -L$(HOME)/opt/lib64

clean:
	$(RM) *.so *.o *.mod

m_solver.o: m_solver_lib.o

