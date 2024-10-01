AF_DIR=$(HOME)/git/afivo

LIBDIRS_2D := $(AF_DIR)/lib_2d $(AF_DIR)/external_libraries/silo/lib	\
$(AF_DIR)/external_libraries/hypre/lib $(CURDIR)
LIBDIRS_3D := $(AF_DIR)/lib_3d $(AF_DIR)/external_libraries/silo/lib	\
$(AF_DIR)/external_libraries/hypre/lib $(CURDIR)
INCDIRS_2D := $(AF_DIR)/lib_2d $(CURDIR)
INCDIRS_3D := $(AF_DIR)/lib_3d $(CURDIR)
LIBS := afivo silo HYPRE gomp	# gomp required for f2py

include $(AF_DIR)/src/makerules.make

.PHONY: all clean

all:	m_solver.o m_solver_3d.o

m_solver.o: libsolver.a
	f2py -c -m poisson m_solver.f90 -lsolver \
	$(addprefix -I,$(INCDIRS_2D)) $(addprefix -L,$(LIBDIRS_2D)) \
	$(addprefix -l,$(LIBS)) --f90flags="$(FFLAGS)"

m_solver_3d.o: libsolver_3d.a
	f2py -c -m poisson_3d m_solver_3d.f90 -lsolver_3d \
	$(addprefix -I,$(INCDIRS_3D)) $(addprefix -L,$(LIBDIRS_3D)) \
	$(addprefix -l,$(LIBS)) --f90flags="$(FFLAGS)"

clean:
	$(RM) *.so *.o *.mod *.a

m_solver_lib.o: INCDIRS=$(INCDIRS_2D)
m_solver_lib.o: LIBDIRS=$(LIBDIRS_2D)
m_solver_lib_3d.o: INCDIRS=$(INCDIRS_3D)
m_solver_lib_3d.o: LIBDIRS=$(LIBDIRS_3D)

libsolver.a: m_solver_lib.o
	$(RM) $@
	$(AR) rcs $@ $^

libsolver_3d.a: m_solver_lib_3d.o
	$(RM) $@
	$(AR) rcs $@ $^
