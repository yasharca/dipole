############################################################################
#
#  Program:         SuperLU_DIST
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   February 4, 1999   version alpha
#
#  Modified:	    September 1, 1999  version 1.0
#                   March 15, 2003     version 2.0
#
# 		    January 18, 2006   Sam Adams
#                                      General Dynamics - Network Systems
#                       works for i386 Linux, with LAM-MPI 7.1.1 and GCC 4.
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT		= _linux


COMMANDER = ${HOME}/code/commander
COMMANDERLIB = -L$(COMMANDER)/src/commander -lcommander
COMMANDERINCLUDE = -I$(COMMANDER)/src/commander -I$(COMMANDER)/src/include

#
#  The name of the libraries to be created/linked to
#
DSuperLUroot 	= ${HOME}/sdss/SuperLU_DIST_3.2
DSUPERLULIB   	= $(DSuperLUroot)/lib/libsuperlu_dist_3.0.a
#
BLASDEF	     	= -DUSE_VENDOR_BLAS
#BLASLIB = ../lib/blas$(PLAT).a
MKLROOT = ~sigurdkn/local/intel_mkl
BLASLIB = -shared-intel -Wl,--start-group  $(MKLROOT)/mkl/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/mkl/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -openmp -lpthread

#BLASLIB      	= /usr/lib/libblas.so.3

METISLIB	= /mn/stornext/u2/yasharca/sdss/ParMetis-3.2.0
PARMETISLIB	= /mn/stornext/u2/yasharca/sdss/ParMetis-3.2.0

# Define the required Fortran libraries, if you use C compiler to link
FLIBS	 	= 

# Define all the libraries
LIBS            = $(DSUPERLULIB) $(BLASLIB) -L$(PARMETISLIB) -lparmetis -L$(METISLIB) -lmetis $(COMMANDERLIB)


INCLUDES = $(COMMANDERINCLUDE)

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         	= ar
ARCHFLAGS    	= cr
RANLIB       	= ranlib

############################################################################
# C compiler setup
CC           	= mpicc
# CFLAGS should be set to be the C flags that include optimization
CFLAGS          = -pipe -O2
#
# NOOPTS should be set to be the C flags that turn off any optimization
NOOPTS		= 
############################################################################
# FORTRAN compiler setup
FORTRAN         = mpif90
F90FLAGS	= 
############################################################################
LOADER	        = $(CC)
LOADOPTS	= 
############################################################################
#  C preprocessor defs for compilation (-DNoChange, -DAdd_, or -DUpCase)
#
#  Need follow the convention of how C calls a Fortran routine.
#
CDEFS        = -DAdd__ 
