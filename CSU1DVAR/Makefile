# makefile for batch processor
# a tab precedes all rules
#
# OBJS  - macro representing all dependent object modules
# 
CC = icc
FC= ifort

CFLAGS = -Zp8
FFLAGS = -g -fPIC -O2 -DLANGUAGE_FORTRAN

LIBES = -L/usr/lib64 -L/usr/local/lib -L/usr/local/hdf5/lib -lnetcdff -lnetcdf -lmfhdf -ldf -ljpeg -lhdf5_hl -lhdf5 -lz -limf -lhdf5_fortran #-lsz
CRTM_INC = /usr/local/CRTM/crtm_v2.2.3/include
CRTM_LIB = /usr/local/CRTM/crtm_v2.2.3/lib
INCFLAGS = -I/usr/local/include \
            -I$(CRTM_INC) \
	   -I/usr/local/hdf5/include -I$(EOS5_DIR)/include
LIBFLAGS = -L$(CRTM_LIB) -lcrtm 
FLAGS = $(INCFLAGS) $(LIBFLAGS) $(LIBES)

OUT = goe

OBJS=\
	define_intam2.o \
	GPM_arraydef_D.o \
	GPM_util_procedures.o \
	GPM_read_procedures.o \
	GPM_rain_procedures.o \
	GPROF_2017_V1.o \
	define_csu1dvar.o \
	miemod.o \
	crtm.o \
	radtran.o \
	csu1dvar_subr.o \
	csu1dvar.o \
	read_geos5_pp.o \
	sea_ice_nt_icecon.o \
	output_1dvar_nc.o \
	integrated_output.o \
	goe.o
#
# function depends on OBJS files and is the target file to be made
#
INTAM2:\
	$(OBJS)
	$(FC) $(CFLAGS) $(FFLAGS) $(OBJS) $(FLAGS) -o $(OUT)

	#rm -f *.o *.mod
#
#dependency of OBJS files
#
clean:
	#-rm -f *.o
	#-rm -f *.mod
# a tab precedes all rules
# default suffix dependency  .c files depend on .o
.suffixes: .f90 .f .c .o 
.c.o:
	$(CC) -c  $(CFLAGS) $<
sea_ice_nt_icecon.o: sea_ice_nt_icecon.c
	$(CC) -c  $(CFLAGS) $<
.f.o:
	$(FC) -c  $(FFLAGS) $(FLAGS) $<
%.o : %.f90
	$(FC) -c  $(FFLAGS) $(FLAGS) $<
