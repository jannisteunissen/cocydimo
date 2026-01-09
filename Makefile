AF_DIR=${CURDIR}/afivo

FC := gfortran
FFLAGS := -O2 -g -std=f2008 -fopenmp -cpp -fPIC -Wall -Wno-unused-dummy-argument -Wl,-z,noexecstack
CFLAGS := -Wall -Winvalid-pch -O3 -fPIC -cpp

BUILD_DIR_2D := build_2d
BUILD_DIR_3D := build_3d

LIBDIRS_2D := $(AF_DIR)/lib_2d $(AF_DIR)/external_libraries/silo/lib	\
$(AF_DIR)/external_libraries/hypre/lib $(CURDIR)/$(BUILD_DIR_2D)
LIBDIRS_3D := $(AF_DIR)/lib_3d $(AF_DIR)/external_libraries/silo/lib	\
$(AF_DIR)/external_libraries/hypre/lib $(CURDIR)/$(BUILD_DIR_3D)
INCDIRS_2D := $(AF_DIR)/lib_2d
INCDIRS_3D := $(AF_DIR)/lib_3d
LIBS := afivo silo HYPRE gomp	# gomp required for f2py

AFIVO_LIB_2D := $(AF_DIR)/lib_2d/libafivo.a
AFIVO_LIB_3D := $(AF_DIR)/lib_3d/libafivo.a

.FORCE:
.PHONY: all lib2d lib3d clean .FORCE

all:	lib2d lib3d

lib2d: $(BUILD_DIR_2D)/libsolver.a
	FFLAGS="$(FFLAGS) -Dfndims=2" CFLAGS="$(CFLAGS) -Dfndims=2" f2py \
	-c -m poisson_2d m_solver.f90 cpp_macros.h -lsolver \
	$(addprefix -I,$(INCDIRS_2D)) $(addprefix -L,$(LIBDIRS_2D)) \
	$(addprefix -l,$(LIBS)) --build-dir $(BUILD_DIR_2D)

lib3d: $(BUILD_DIR_3D)/libsolver.a
	FFLAGS="$(FFLAGS) -Dfndims=3" CFLAGS="$(CFLAGS) -Dfndims=3" f2py \
	-c -m poisson_3d m_solver.f90 cpp_macros.h -lsolver \
	$(addprefix -I,$(INCDIRS_3D)) $(addprefix -L,$(LIBDIRS_3D)) \
	$(addprefix -l,$(LIBS)) --build-dir $(BUILD_DIR_3D)

clean:
	$(RM) *.so *.o *.mod *.a
	$(RM) -r $(BUILD_DIR_2D) $(BUILD_DIR_3D)

$(BUILD_DIR_2D)/m_solver_lib.o: INCDIRS=$(INCDIRS_2D)
$(BUILD_DIR_2D)/m_solver_lib.o: $(AFIVO_LIB_2D)
$(BUILD_DIR_3D)/m_solver_lib.o: INCDIRS=$(INCDIRS_3D)
$(BUILD_DIR_3D)/m_solver_lib.o: $(AFIVO_LIB_3D)

$(BUILD_DIR_2D)/libsolver.a: $(BUILD_DIR_2D)/m_solver_lib.o
	$(RM) $@
	$(AR) rcs $@ $^

$(BUILD_DIR_3D)/libsolver.a: $(BUILD_DIR_3D)/m_solver_lib.o
	$(RM) $@
	$(AR) rcs $@ $^

$(BUILD_DIR_2D)/%.o: %.f90
	mkdir -p $(BUILD_DIR_2D)
	$(FC) -c -o $@ $< $(FFLAGS) -J $(BUILD_DIR_2D) -Dfndims=2 $(addprefix -I,$(INCDIRS))

$(BUILD_DIR_3D)/%.o: %.f90
	mkdir -p $(BUILD_DIR_3D)
	$(FC) -c -o $@ $< $(FFLAGS) -J $(BUILD_DIR_3D) -Dfndims=3 $(addprefix -I,$(INCDIRS))

$(AFIVO_LIB_2D):
	$(MAKE) -C afivo lib_2d

$(AFIVO_LIB_3D):
	$(MAKE) -C afivo lib_3d
