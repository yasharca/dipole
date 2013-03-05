#######################################################################
#
#  This makefile creates the Fortran example programs for the
#  linear equation routines in SuperLU_DIST.
#
#  Creation date:   July 29, 2003   version 2.0
#  Modified:        Oct. 22, 2012   version 3.2
#
#######################################################################
.SUFFIXES: 
.SUFFIXES: .f90 .c .o
include make.inc
INCLUDEDIR = -I$(DSuperLUroot)/SRC

MKLROOT = ~sigurdkn/local/intel_mkl
LDFLAGS = -lhealpix -lhdf5_fortran -lhdf5 -lcfitsio -openmp -lgsl -lslatec -Wl,--start-group  $(MKLROOT)/mkl/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/mkl/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -openmp -lpthread
#FFLAGS   = -O3 -traceback
FFLAGS   = -O0 -g -C -traceback
#FFLAGS   = -vec_report0 -O0 -g -check bounds -check format -check pointers -check uninit -check output_conversion -assume byterecl -traceback -heap-arrays 16384
#FFLAGS  = -g -traceback -fp-stack-check -check all -fpe0 -warn all

#F90FLAGS	= $(FFLAGS) -qfree -qsuffix=f=f90  -qflag=w:w

F_MOD	= superlupara.o superlu_mod.o
C_DWRAP	= dcreate_dist_matrix.o superlu_c2f_dwrap.o
C_ZWRAP	= zcreate_dist_matrix.o superlu_c2f_zwrap.o

F_DEXM	= $(F_MOD) dhbcode1.o f_pddrive.o
F_ZEXM	= $(F_MOD) zhbcode1.o f_pzdrive.o
F_5x5 	= $(F_MOD) f_5x5.o sp_ienv.o

F_DPL   = $(F_MOD) dhbcode1.o f_dipolelike.o quiet_hdf_mod.o quiet_mapfile_mod.o
F_DPLMCMC   = $(F_MOD) dhbcode1.o f_dipolemcmc.o quiet_hdf_mod.o quiet_mapfile_mod.o mcmcmodule.o
F_TST   = $(F_MOD) dhbcode1.o f_test.o quiet_hdf_mod.o quiet_mapfile_mod.o
F_ZTST   = $(F_MOD) dhbcode1.o f_ztest.o quiet_hdf_mod.o quiet_mapfile_mod.o


all: f_dipolelike f_dipolemcmc f_test f_ztest

quiet_mapfile_mod.o: quiet_hdf_mod.o
f_dipolelike.o: quiet_hdf_mod.o quiet_mapfile_mod.o

f_dipolemcmc.o: quiet_hdf_mod.o quiet_mapfile_mod.o mcmcmodule.o

f_test.o: quiet_hdf_mod.o quiet_mapfile_mod.o

f_ztest.o: quiet_hdf_mod.o quiet_mapfile_mod.o

f_dipolelike: $(F_DPL) $(C_ZWRAP) $(DSUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(F_DPL) $(C_ZWRAP) $(LIBS) -o $@ $(LDFLAGS)

f_dipolemcmc: $(F_DPLMCMC) $(C_ZWRAP) $(DSUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(F_DPLMCMC) $(C_ZWRAP) $(LIBS) -o $@ $(LDFLAGS)

f_test: $(F_TST) $(C_DWRAP) $(DSUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(F_TST) $(C_DWRAP) $(LIBS) -o $@ $(LDFLAGS)

f_ztest: $(F_ZTST) $(C_ZWRAP) $(DSUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(F_ZTST) $(C_ZWRAP) $(LIBS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) $(BLASDEF) $(INCLUDEDIR) -c $< $(VERBOSE)

.f90.o:
	$(FORTRAN) $(F90FLAGS) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

clean:	
	rm -f *.o *.mod f_*drive f_5x5 f_dipolelike f_dipolemcmc f_test f_ztest


