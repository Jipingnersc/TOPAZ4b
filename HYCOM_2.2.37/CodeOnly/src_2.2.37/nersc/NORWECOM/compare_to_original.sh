#!/bin/bash

hycom_dir=/home/nersc/annettes/Progs/HYCOM_FELLES/
norwecom_dir=${hycom_dir}Eco/Norwecom/
mpi_norwecom_dir=/home/nersc/ceciha/Progs/HYCOM_2.1.03_MPI_NORWECOM/ECO-MPI/NOR05/
mpi_eco_dir=/home/nersc/ceciha/Progs/HYCOM_2.1.03_MPI_NORWECOM/ECO-MPI/

cmp biocom.h          ${norwecom_dir}biocom.h
cmp blkdat_eco.f      ${hycom_dir}blkdat_eco.F
cmp comdeck.h         ${norwecom_dir}comdeck.h
cmp globals.doc       ${norwecom_dir}globals.doc
cmp globals.h         ${norwecom_dir}globals.h
cmp iodeck.h          ${norwecom_dir}iodeck.h

cmp m_ecomodels_read_restart.F  ${norwecom_dir}m_ecomodels_read_restart.F
cmp m_ecomodels_save_restart.F  ${norwecom_dir}m_ecomodels_save_restart.F
cmp m_NOR05_affin.F             ${norwecom_dir}m_NOR05_affin.F 
cmp m_NOR05_area_pp.F           ${norwecom_dir}m_NOR05_area_pp.F
cmp m_NOR05_atmnit.F            ${mpi_norwecom_dir}m_NOR05_atmnit.F
cmp m_NOR05_biochm.F            ${mpi_norwecom_dir}m_NOR05_biochm.F
cmp m_NOR05_biosun.F            ${mpi_norwecom_dir}m_NOR05_biosun.F
cmp m_NOR05_botsource.F         ${mpi_norwecom_dir}m_NOR05_botsource.F
cmp m_NOR05_botstress.F         ${mpi_norwecom_dir}m_NOR05_botstress.F
cmp m_NOR05_checknan.F          ${mpi_norwecom_dir}m_NOR05_checknan.F
cmp m_NOR05_comprad.F           ${mpi_norwecom_dir}m_NOR05_comprad.F
cmp m_NOR05_daynumb.F           ${mpi_norwecom_dir}m_NOR05_daynumb.F
cmp m_NOR05_deadbot.F           ${norwecom_dir}m_NOR05_deadbot.F
cmp m_NOR05_declin.F            ${mpi_norwecom_dir}m_NOR05_declin.F
cmp m_NOR05_detritus.F          ${mpi_norwecom_dir}m_NOR05_detritus.F
cmp m_NOR05_fx.F                ${mpi_norwecom_dir}m_NOR05_fx.F
cmp m_NOR05_initialize_ecosys.F ${mpi_norwecom_dir}m_NOR05_initialize_ecosys.F
cmp m_NOR05_initialize_sed.F    ${mpi_norwecom_dir}m_NOR05_initialize_sed.F
cmp m_NOR05_irrsur.F            ${mpi_norwecom_dir}m_NOR05_irrsur.F
cmp m_NOR05_micomsink.F         ${mpi_norwecom_dir}m_NOR05_micomsink.F
cmp m_NOR05_nconservation.F     ${mpi_norwecom_dir}m_NOR05_nconservation.F
cmp m_NOR05_oxysat.F            ${mpi_norwecom_dir}m_NOR05_oxysat.F
cmp m_NOR05_phyto_growth.F      ${mpi_norwecom_dir}m_NOR05_phyto_growth.F
cmp m_NOR05_phyto_sed.F         ${mpi_norwecom_dir}m_NOR05_phyto_sed.F
#cmp m_NOR05_phytosink.F         ${norwecom_dir}m_NOR05_phytosink.F
#cmp m_NOR05_primp.F             ${norwecom_dir}m_NOR05_primp.F
cmp m_NOR05_radiation.F         ${mpi_norwecom_dir}m_NOR05_radiation.F
cmp m_NOR05_relaxation.F        ${norwecom_dir}m_NOR05_relaxation.F
cmp m_NOR05_resusp.F            ${mpi_norwecom_dir}m_NOR05_resusp.F
cmp m_NOR05_riverloads.F        ${mpi_norwecom_dir}m_NOR05_riverloads.F
cmp m_NOR05_sed_restart.F       ${norwecom_dir}m_NOR05_sed_restart.F
cmp m_NOR05_slim.F              ${mpi_norwecom_dir}m_NOR05_slim.F
cmp m_NOR05_snell.F             ${mpi_norwecom_dir}m_NOR05_snell.F
cmp m_NOR05_source.F            ${mpi_norwecom_dir}m_NOR05_source.F
cmp m_NOR05_trans.F             ${mpi_norwecom_dir}m_NOR05_trans.F
cmp m_NOR05_wqual.F             ${mpi_norwecom_dir}m_NOR05_wqual.F
#cmp mod_checknan.F              ${norwecom_dir}mod_checknan.F
cmp mod_necessary_ecovars.F90   ${mpi_eco_dir}mod_necessary_ecovars.F90
cmp sedcom.h                    ${mpi_eco_dir}sedcom.h
#cmp state.h                     ${mpi_eco_dir}state.h
