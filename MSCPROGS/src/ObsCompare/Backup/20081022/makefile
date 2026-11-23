VPATH = .:RCS:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h .c


MACHINE=$(shell uname -n | sed "s/\..*//")

include make.$(MACHINE)


# Rules for running cpp and updating files in TMP directory
.H.h:
	rm -f ./TMP/$*.h
	cat MODEL.CPP $*.H | $(CPP) $(CPPFLAGS)  > ./TMP/$*.h


.F90.o:
	rm -f ./TMP/$*.f90
	cat MODEL.CPP $*.F90 | $(CPP) $(CPPFLAGS)   > ./TMP/$*.f90
	cd ./TMP ; $(CF90) -c $(FFLAGS) $(F90FLG) -o $*.o $*.f90  

.F.o:
	rm -f ./TMP/$*.f
	cat MODEL.CPP $*.F | $(CPP) $(CPPFLAGS)  > ./TMP/$*.f
	cd ./TMP ; $(CF77) -c $(FFLAGS) $(F77FLG) -o $*.o $*.f  

# One call to c-routine  (retrieves file size)
.c.o:
	rm -f ./TMP/$*.c
	cat MODEL.CPP $*.c | $(CPP) $(CPPFLAGS)   >  ./TMP/$*.c
	cd ./TMP ; $(CC) -c $(CFLAGS)  -o $*.o $*.c




#TARGET1 = rawdaily2netcdf
TARGET  = argocmp
TARGET2 = argoreduce
TARGET3 = postprocess_mersea
TARGET4 = slacmp
TARGET5 = mldcmp

 
# Stay off daily2modstat until dbx bug is resolved
all: $(TARGET)  $(TARGET2)


###############################################################################
# To generate argocmp
OBJECTS=  mod_xc.o mod_za.o mod_common.o mod_year_info.o \
		  m_datetojulian.o mod_levitus.o \
		  mod_netcdf_pars.o m_handle_err.o m_def_netcdf_var.o \
		  mod_confmap.o m_initconfmap.o m_oldtonew.o m_spherdist.o m_pivotp.o \
		  m_bilincoeff.o m_nearestpoint.o m_calcconst.o m_calcu.o \
          mod_spline_calc.o m_read_coriolis.o m_dump_coriolis.o \
          m_pack_short.o m_year_day.o m_datetojulian.o \
          m_get_hyc_infiles.o m_get_depth_levels.o \
		  m_bigrid.o m_indxi.o m_indxj.o  \
		  m_get_grid.o mod_read_dailyab.o \
		  m_read_blkdat.o mod_read_rstab.o mod_read_weekab.o m_read_field_wrapper.o \
		  mod_readpak.o p_argocmp.o wtime.o m_get_grid.o \
		  machine.o machi_c.o \
		  wtime.o pakk.o

$(TARGET): $(OBJECTS) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET)  $(OBJECTS) $(LIBS) 
###############################################################################
# To generate argoreduce
OBJECTS2= p_argoreduce.o mod_sphere_tools.o

$(TARGET2): $(OBJECTS2) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2)  $(OBJECTS2) $(LIBS) 
###############################################################################
# To generate postprocess_mersea
OBJECTS3= p_postprocess_mersea.o 

$(TARGET3): $(OBJECTS3) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET3)  $(OBJECTS3) $(LIBS) 
###############################################################################
# To generate slacmp
OBJECTS4= mod_xc.o mod_za.o machine.o machi_c.o wtime.o m_getfiletype.o m_indxi.o m_indxj.o mod_common.o\
          m_year_day.o m_datetojulian.o m_read_blkdat.o m_bigrid.o m_get_grid.o mod_readpak.o mod_read_dailyab.o \
		  mod_year_info.o m_getfileinfo.o mod_sphere_tools.o \
          mod_confmap.o m_initconfmap.o m_oldtonew.o m_pivotp.o \
          mod_read_weekab.o mod_read_rstab.o mod_sphere_tools.o \
          m_read_field_wrapper.o m_duacsread.o p_slacompare.o pakk.o

$(TARGET4): $(OBJECTS4) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET4)  $(OBJECTS4)  $(LIBNCDAP)
###############################################################################
###############################################################################
OBJECTS5= mod_xc.o mod_za.o machine.o machi_c.o wtime.o m_getfiletype.o m_indxi.o m_indxj.o mod_common.o\
          m_year_day.o m_datetojulian.o m_read_blkdat.o m_bigrid.o m_get_grid.o mod_readpak.o mod_read_dailyab.o \
			 mod_netcdf_helpers.o m_parse_blkdat.o \
		    mod_year_info.o m_getfileinfo.o mod_sphere_tools.o \
          mod_confmap.o m_initconfmap.o m_oldtonew.o m_pivotp.o \
          mod_read_weekab.o mod_read_rstab.o mod_sphere_tools.o mod_regions.o \
          m_read_field_wrapper.o m_dBM_MLD_read.o p_mldcompare.o pakk.o

$(TARGET5): $(OBJECTS5) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET5)  $(OBJECTS5)  $(LIBS)
###############################################################################

clean:
	cd ./TMP ; rm *.f  *.o *.f90 *.h *.mod


new: depend

depend:
	./mkdepend.pl | sort -u > depends.file

include depends.file

