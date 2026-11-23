      program main
      implicit none
      include 'netcdf.inc'
c
      integer status,ncid,timeid,fieldid,start(3),count(3),stride(3),i,j
C      integer*2 ifield(londim,latdim),i2_mval
      real*8 time
      real*4 offset,scale_factor
c
c --- open netcdf file
      status=nf_open('icec.sfc.gauss.2010.nc',nf_nowrite,ncid)
      if (status.ne.nf_noerr) then
        write(*,'(4a)') 'nf_open: ',filename,',',nf_strerror(status)
        call xchalt('(read_atm_field)')
               stop '(read_atm_field)'
      endif
      end
