      module m_duacsread
      contains
      subroutine duacsread(slaout,lonout,latout,nlonin,nlatin,
     &                     fname,vname)
      ! Input arrays are hardcoded for now - pointer is a viable option
      !
      ! Routine to compare modeled sla with observations 
      ! Program takes three arguments  - first is model
      ! field (only DODS link for now), second is SLA
      ! field - also DODS. Third is some form of time 
      ! unit (to be decided)
      !
      ! This program must be linked against nc-dap
      implicit none
      include 'nc-dap.inc'

      integer, intent(in) :: nlonin, nlatin
      real   , intent(out) :: slaout(nlatin,nlonin)
      real   , intent(out) :: lonout(nlonin)
      real   , intent(out) :: latout(nlatin)
      character(len=*), intent(in) :: fname,vname

      integer :: ncid, stat, ndims, nvars, ngatts,varid,dvarid
      integer :: unlimdimid
      character(len=NF_MAX_NAME) :: varname
      integer, dimension(NF_MAX_VAR_DIMS) :: vardimids, vardimlen,
     &   varcoordtype
      character(len=NF_MAX_NAME), dimension(nf_max_var_dims) :: 
     &   vardimname
      integer :: xtype,varndims,nvaratts
      integer :: i,dimid(1),j
      double precision,  allocatable :: lats(:),lons(:)
      real*4,  allocatable :: sla(:,:)
      integer,  allocatable :: times(:)
      integer :: nlons,nlats, ntimes
      integer, allocatable :: indexes(:)
      integer :: timeind,start(3),count(3)

      stat = NF_OPEN(trim(fname), NF_NOWRITE , ncid)
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if
      print *,'Reading SLA from ',trim(fname)

      ! Get variable with name  
      stat = NF_INQ_VARID (NCID, vname, varid)
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if

      ! Get info on this variable
      stat = NF_INQ_VAR   (NCID, VARID, varname, xtype,
     &                           varndims, vardimids,
     &                           nvaratts)
      if (stat /= nf_noerr) then
         print *,nf_strerror(stat)
         call exit(3)
      end if

      ! Search for coordinate variables - this is terse, but should be
      ! safe
      print *,'Coordinate variable search:'
      do i=1,varndims


         ! Get dimension name
         stat = NF_INQ_DIMNAME (NCID,varDIMIDs(i), vardimname(i))
         if (stat /= nf_noerr) then
            print *,nf_strerror(stat)
            call exit(3)
         end if
         print *,' dimname:',trim(vardimname(i))

         ! And its length
         stat = NF_INQ_DIMLEN  (NCID, varDIMIDs(i), vardimlen(i))
         if (stat /= nf_noerr) then
            print *,nf_strerror(stat)
            call exit(3)
         end if
         !print *,'stat, dimlen :',stat,vardimlen(i)

         ! get corresponding variable (coordinate variable)
         stat = NF_INQ_VARID    (NCID, vardimname(i)(:), dVARID)
         if (stat /= nf_noerr) then
            print *,nf_strerror(stat)
            call exit(3)
         end if
         !print *,'inq_varname stat, dimvarid :',stat,dvarid

         ! If we are here, we have coord variables for this dim

         ! get type
         stat = NF_INQ_VARTYPE  (NCID, dVARID, varcoordtype(i))
         if (stat /= nf_noerr) then
            print *,nf_strerror(stat)
            call exit(3)
         end if
         !print *,'Coord var type:',varcoordtype(i)

         ! Check that variable is 1D
         stat = NF_INQ_VARNDIMS (NCID, dvarid, ndims)
         if (stat /= nf_noerr) then
            print *,nf_strerror(stat)
            call exit(3)
         end if
         !print *,'Coord var ndim:',ndims

         ! Check its dim id
         stat = NF_INQ_VARDIMID (NCID, dvarid, dimid)
         if (stat /= nf_noerr) then
            print *,nf_strerror(stat)
            call exit(3)
         end if
         if (dimid(1)/=vardimids(i)) then
            print *,'dimension coordinate  mismatch'
            call exit(4)
         end if
         !print *,'Coord var dimid, current dimid:',dimid,vardimids(i)

         ! If we are here, we have appropriate coord variables for this dim

         ! Get coord variable 
         if (trim(vardimname(i))=="NbLatitudes") then

            if (varcoordtype(i)/=NF_DOUBLE)then
               print *,'LAtitude isassumed double'
               print *,nf_double,varcoordtype(i)
               call exit(1)
            end if

            nlats=vardimlen(i)
            allocate(lats    (nlats))
            stat =  NF_get_VAR_double(NCID, dVARID, lats)
            !print *,lats
            if (stat/=NF_NOERR)then
               print *,'Error reading Latitude'
               call exit(2)
            end if

         elseif (trim(vardimname(i))=="NbLongitudes") then
            if (varcoordtype(i)/=NF_DOUBLE)then
               print *,'Longitude isassumed double'
               call exit(1)
            end if
            nlons=vardimlen(i)
            allocate(lons    (nlons))
            stat =  NF_get_VAR_double(NCID, dVARID, lons)
            !print *,lons
            if (stat/=NF_NOERR)then
               print *,'Error reading longitude'
               call exit(2)
            end if

         elseif (trim(vardimname(i))=="time") then
            if (varcoordtype(i)/=NF_INT)then
               print *,'time is assumed integer'
               call exit(1)
            end if
            ntimes=vardimlen(i)
            allocate(times    (ntimes))
            stat =  NF_get_VAR_int(NCID, dVARID, times)
            !print *,times
            if (stat/=NF_NOERR)then
               print *,'Error reading time'
               call exit(2)
            end if
         end if


      end do

      ! Testing
      timeind=50
      start=(/1,1,timeind/)
      !count=(/nlats,nlons,1/)
      !allocate(sla(nlats,nlons))
      allocate(sla(vardimlen(1),vardimlon(2)))
      count=(/vardimlen(1),vardimlen(2),1/)
      stat = NF_GET_VARA_REAL(NCID, VARID, START, COUNT, sla)
      if (stat/=NF_NOERR)then
         print *,'Error reading sla'
         call exit(2)
      end if

      !print *,sla

      if (nlats==nlatin .and. nlons==nlonin) then
         lonout=lons
         latout=lats
         slaout=sla
      else
         print *,'array shapes does not conform'
         print *,nlonin,nlons
         print *,nlatin,nlats
      end if
      deallocate(times,lons,lats,sla)


      end subroutine
      end module

