VPATH = .:RCS:TMP

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 .F .f .H .h .c


MACHINE=$(shell uname -n | sed "s/\..*//")

include make.$(MACHINE).ifc


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
TARGET2 = daily2regunc
#TARGET3 = daily2stationnc
#TARGET4 = daily2sectionnc
#TARGET5 = daily2obsstat
#TARGET6 = daily2modstat
#TARGET7 = daily2pak
#TARGET8 = levitus2obs
#TARGET9 = topazwatch
#TARGET10= daily2mosf
#TARGET11= dailymean
#TARGET12= dailyfixname
#TARGET13= dayofweek
#TARGET14= safotoobs
#TARGET15= coriolis2obs
#TARGET16= tide2obs
#TARGET17= innovstat

 
# Stay off daily2modstat until dbx bug is resolved
all: $(TARGET2) 


###############################################################################
# To generate daily2regunc
OBJECTS2= mod_xc.o mod_za.o mod_dimensions.o  mod_common.o mod_year_info.o \
		  mod_netcdf_pars.o m_handle_err.o m_def_netcdf_var.o \
		  mod_confmap.o m_initconfmap.o m_oldtonew.o m_spherdist.o m_pivotp.o \
          mod_toregugrid.o  mod_spline_calc.o m_rotate.o \
          m_strmf_eval.o m_mixlayer_depths.o m_pack_short.o m_year_day.o \
          m_mersea_prepare.o m_get_daily_infiles.o m_get_depth_levels.o \
		  m_get_regugrid.o m_datetojulian.o m_get_grid.o mod_read_dailyab.o \
		  p_daily2regunc.o wtime.o m_get_grid.o m_bilincoeff.o m_calcu.o \
		  m_calcconst.o m_bigrid.o m_indxi.o m_indxj.o 

$(TARGET2): $(OBJECTS2) 
	cd ./TMP ; $(LD) $(LINKFLAGS) -o ../$(TARGET2)  $(OBJECTS2) $(LIBS) 

###############################################################################
# To generate levitus2obs
OBJECTS8=mod_dimensions.o mod_measurement.o m_read_levitus.o m_read_PHC.o \
		 m_get_mod_grid.o mod_states.o m_initconfmap.o mod_confmap.o    \
		 m_oldtonew.o m_pivotp.o p_levitus2obs.o m_bilincoeff.o           \
		 mod_common.o mod_atmforc.o  mod_year_info.o m_tecfld.o m_insitu2theta.o
$(TARGET8): $(INC2)  $(OBJECTS8)  $(DEBOBJ77) $(DEBOBJ90)
	cd ./TMP ; $(LD) $(LINKFLAGS8) -o ../$(TARGET8)  $(OBJECTS8) $(LIBS) 
###############################################################################
$(TARGET9): $(INC2) $(OMOD) $(OBJECTS9)  $(DEBOBJ77) $(DEBOBJ90)
	cd ./TMP ; $(LD) $(LINKFLAGS9) -o ../$(TARGET9) $(OMOD) $(OBJECTS9) $(LIBS) 
###############################################################################
$(TARGET10): $(INC2) $(OMOD) $(OBJECTS10)  $(DEBOBJ77) $(DEBOBJ90)
	cd ./TMP ; $(LD) $(LINKFLAGS10) -o ../$(TARGET10) $(OMOD) $(OBJECTS10) $(LIBS) 
###############################################################################
OBJECTS11=mod_dimensions.o mod_common.o mod_daily_average.o mod_year_info.o \
		   m_blkdat.o m_gridsize.o m_geopar.o m_inicon.o m_dpthuv.o \
		   m_year_day.o m_get_daily_infiles.o m_interface_bug.o \
		   m_datetojulian.o m_indxi.o m_indxj.o p_dailymean.o m_spherdist.o \
		   m_mercrot.o m_blkini.o m_blkinr.o m_blkinvoid.o m_blkinl.o \
		   m_mercator.o m_bigrid.o m_sshcalc.o

$(TARGET11): $(OBJECTS11) $(INC2) 
	cd ./TMP ; $(LD) $(LINKFLAGS11) -o ../$(TARGET11) $(OBJECTS11) $(LIBS) 
###############################################################################
###############################################################################
OBJECTS12=mod_dimensions.o mod_common.o mod_daily_average.o mod_year_info.o \
		   m_blkdat.o m_gridsize.o m_geopar.o m_inicon.o m_dpthuv.o \
		   m_year_day.o m_get_daily_infiles.o m_interface_bug.o \
		   m_datetojulian.o m_indxi.o m_indxj.o p_dailyfixname.o m_spherdist.o \
		   m_mercrot.o m_blkini.o m_blkinr.o m_blkinvoid.o m_blkinl.o \
		   m_mercator.o m_bigrid.o m_sshcalc.o

$(TARGET12): $(OBJECTS12) $(INC2) 
	cd ./TMP ; $(LD) $(LINKFLAGS12) -o ../$(TARGET12) $(OBJECTS12) $(LIBS) 
###############################################################################
OBJECTS13=m_datetojulian.o p_dayofweek.o

$(TARGET13): $(OBJECTS13) $(INC2) 
	cd ./TMP ; $(LD) $(LINKFLAGS13) -o ../$(TARGET13) $(OBJECTS13) $(LIBS) 
###############################################################################
###############################################################################
# To generate safo2obs
OBJECTS14=mod_dimensions.o mod_measurement.o m_read_safo.o m_read_PHC.o \
		 m_get_mod_grid.o mod_states.o m_initconfmap.o mod_confmap.o    \
		 m_oldtonew.o m_pivotp.o p_safotoobs.o m_bilincoeff.o           \
		 mod_common.o mod_atmforc.o  mod_year_info.o m_tecfld.o m_insitu2theta.o
$(TARGET14): $(INC2)  $(OBJECTS14)  $(DEBOBJ77) $(DEBOBJ90)
	cd ./TMP ; $(LD) $(LINKFLAGS14) -o ../$(TARGET14)  $(OBJECTS14) $(LIBS) 
###############################################################################
# To generate coriolis2obs
OBJECTS15=mod_dimensions.o mod_measurement.o m_read_coriolis.o          \
		 m_get_mod_grid.o mod_states.o m_initconfmap.o mod_confmap.o    \
		 m_oldtonew.o m_pivotp.o p_coriolis2obs.o m_bilincoeff.o           \
		 mod_common.o mod_atmforc.o  mod_year_info.o m_tecfld.o         \
		 m_get_depth_interval.o m_datetojulian.o
$(TARGET15): $(INC2)  $(OBJECTS15)  
	cd ./TMP ; $(LD) $(LINKFLAGS15) -o ../$(TARGET15)  $(OBJECTS15) $(LIBS) 
###############################################################################
# To generate tide2obs
OBJECTS16=mod_dimensions.o mod_measurement.o m_read_tide.o          \
		 m_get_mod_grid.o mod_states.o m_initconfmap.o mod_confmap.o    \
		 m_oldtonew.o m_pivotp.o p_tide2obs.o m_bilincoeff.o           \
		 mod_common.o mod_atmforc.o  mod_year_info.o m_tecfld.o         \
		 m_get_depth_interval.o m_datetojulian.o rosame_filtre_demerliac_v00.o \
		 m_spherdist.o
$(TARGET16): $(INC2)  $(OBJECTS16)  
	cd ./TMP ; $(LD) $(LINKFLAGS16) -o ../$(TARGET16)  $(OBJECTS16) $(LIBS) 
###############################################################################
# To generate innovstat
OBJECTS17=mod_innovations_io.o m_datetojulian.o p_innovstat.o mod_sphere_tools.o \
		  m_get_regions.o
$(TARGET17): $(INC2)  $(OBJECTS17)  
	cd ./TMP ; $(LD) $(LINKFLAGS17) -o ../$(TARGET17)  $(OBJECTS17) $(LIBS) 
###############################################################################

clean:
	cd ./TMP ; rm *.f  *.o *.f90 *.h *.mod


new: depend

depend:
	mkdepend.pl | sort -u > depends.file

include depends.file

