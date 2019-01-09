EXEC1 = geos5_pp
##############################################################################
MACHINE = LINUX
FC = ifort -traceback -warn nogeneral

# for L1C stuff:
CC = cc
CFLAGS = -g  -D$(MACHINE) -DX_WCHAR
TKLIB = -ltkcselect -ltkchdf4algs -ltkchdf4 -ltkchdf5algs -ltkchdf5 -ltkchelper -ltkctkTSDIS -ltkc
IOINC = -I$(TKIO)/inc/fortcode -I/$(TKIO)/inc/ccode -I/usr/include/libxml2
HDFINC = -I$(HDF5_INC)

FFLAGS =  -g -fPIC -D$(MACHINE) -O2 -DD -DLANGUAGE_FORTRAN -DH5_USE_16_API

LIBES = -lmfhdf -ldf -lhdf5 -lhdf5_hl -ljpeg -lz -lxml2 -lhdf5_fortran -lnetcdff -lnetcdf

###########################################################################
# HDF4_INC, HDF4_LIB, HDF5_INC, HDF5_LIB, TKIO environment variables should be 
# set in .cshrc
INC = -I/usr/local/include -I/usr/local/hdf5/include -I/usr/include/libxml2

LIBDIR = -L/usr/local/hdf5/lib -L/usr/local/lib -L/usr/lib64  -L$(TKIO)/lib
FLAGS = $(INC) $(IOINC) $(LIBDIR) $(LIBES) $(TKLIB)

##############################################################################
SOURCECODE= \
   $(EXEC1).F

OBJECTSC=\
	pp_definition.o \
	pp_commoncode.o \
	read_l1c.o \
	read_l1r.o \
	read_geos5.o \
	gpm.o \
	geos5_pp.o

############################################################################## 
all:            $(EXEC1)

$(EXEC1):       $(OBJECTSC)
	        $(FC) -o $(EXEC1) $(OBJECTSC)\
	        $(LIBDIR) $(TKLIB) $(LIBES)
	        rm -f *.o *.mod

clean:
	rm -f *.o *.mod

##############################################################################

.SUFFIXES:  .f90 .c .o
.F.o:
	$(FC) -c  $(IOINC) $(FFLAGS) $(HDFINC) $(FLAGS) $<
.f90.o:
	$(FC) -c  $(IOINC) $(FFLAGS) $(HDFINC) $(FLAGS) $<

