VPATH = .:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h .c


include ../Make.Inc/make.inc

# This is the hycom utility library
# NB: Path is relative to TMP dir
LIBS:= -L../../../lib/ -lhycnersc -lconfmap $(LIBS)

# This is the hycom utility includes (modules)
# NB: Path is relative to TMP dir
INCLUDE:= -I../../../include/ -I../ -I$(NETCDF_INCLUDE_OPTS)


CPPFLAGS:= $(CPPFLAGS) -UMULT_GRP_PER_FILE
#CPPFLAGS:= $(CPPFLAGS) 

# AIX (xlf) peculiarities. NB: Only necessary for Fortran files
ifneq (,$(findstring -DAIX,$(CPPFLAGS)))
   subs=-WF,-
   CPPFLAGS:=$(subst -,$(subs),$(CPPFLAGS))
endif

# Rules for running cpp and updating files in TMP directory
.F90.o:
	cd ./TMP ; $(CF90) -c $(CPPFLAGS) $(FFLAGS) $(F90FLG) $(INCLUDE) -o $*.o ../$<

.F.o:
	cd ./TMP ; $(CF77) -c $(CPPFLAGS) $(FFLAGS) $(F77FLG) $(INCLUDE) -o $*.o ../$<





TARGET1 = daily2proj 
TARGET2 = analy2proj 
TARGET3 = clim2proj

all: $(TARGET1) 
#all: $(TARGET1) $(TARGET2) $(TARGET3)

###############################################################################
# To generate daily2proj
OBJECTS1=\
	m_nearestpoint.o \
	mod_toproj.o\
	mod_rotate.o\
	m_strmf_eval.o \
	mod_levitus.o  \
	mod_read_netcdf.o  \
	m_read_mean_ssh.o \
	TP4_grid.o \
	p_daily_proj.o 

$(TARGET1): $(OBJECTS1)  
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET1)  $(OBJECTS1) $(LIBS) 

###############################################################################
# To generate analy2proj
OBJECTS2=\
	m_nearestpoint.o \
	mod_toproj.o\
	mod_rotate.o\
	m_strmf_eval.o \
	mod_levitus.o  \
	mod_read_netcdf.o  \
	m_read_mean_ssh.o \
	p_analy2proj.o 

$(TARGET2): $(OBJECTS2)  
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2)  $(OBJECTS2) $(LIBS) 
###############################################################################
# To generate woa2proj
OBJECTS3=\
	m_nearestpoint.o \
	mod_toproj.o\
	mod_rotate.o\
	m_strmf_eval.o \
	mod_levitus.o  \
	mod_read_netcdf.o  \
	m_read_mean_ssh.o \
        woa13_grid.o \
	p_woa2proj.o 

$(TARGET3): $(OBJECTS3)  
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET3)  $(OBJECTS3) $(LIBS) 
###############################################################################

clean:
	rm ./TMP/*.f* ./TMP/*.o TMP/*.h TMP/*.mod $(TARGET1) $(TARGET2) $(TARGET3) 

install : all
	mkdir -p ../../bin
	cp analy2proj ../../bin

