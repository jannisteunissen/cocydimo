AF_DIR=$(HOME)/git/afivo

LIBDIRS := $(AF_DIR)/lib_2d $(AF_DIR)/external_libraries/silo/lib	\
$(AF_DIR)/external_libraries/hypre/lib
INCDIRS := $(AF_DIR)/lib_2d
LIBS := afivo silo HYPRE gomp	# gomp required for f2py

include $(AF_DIR)/src/makerules.make

.PHONY: all clean

all:	m_solver.o
	f2py -c -m poisson m_solver.f90 m_solver_lib.o \
	$(addprefix -I,$(INCDIRS)) $(addprefix -L,$(LIBDIRS)) \
	$(addprefix -l,$(LIBS)) --f90flags="$(FFLAGS)"
clean:
	$(RM) *.so *.o *.mod

m_solver.o: m_solver_lib.o
