module m_tecfld
contains
subroutine tectest(fname,ii,jj,fld)
   implicit none
   character(len=*), intent(in) :: fname
   integer, intent(in) :: ii,jj
   real, intent(in) :: fld(ii,jj,5)
   integer m
   integer i,j

   open(10,file='tmed'//fname//'.dat',status='unknown')
      write(10,*)'TITLE = "',fname,'"'
      write(10,*)'VARIABLES = "i-index" "j-index" "',fname,'"'
      write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',ii,', J=',jj,', K=1'


! Positions in lon-lat
      write(10,'(30I4)')((i,i=1,ii),j=1,jj)
      write(10,'(30I4)')((j,i=1,ii),j=1,jj)

      write(10,900)((fld(i,j,1),i=1,ii),j=1,jj)

      do m=2,5
         write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',ii,', J=',jj,', K=1'
         write(10,*)'D=(1,2)'
         write(10,900)((fld(i,j,m),i=1,ii),j=1,jj)
      enddo

   close(10)
 900 format(10(1x,e12.5))

end subroutine tectest


subroutine tecfld(fname,fld,ii,jj,depths)
   implicit none
   character(len=*), intent(in) :: fname
   integer, intent(in) :: ii,jj
   real, intent(in) :: fld(ii,jj)
   real, intent(in) :: depths(ii,jj)
   integer i,j

   character*4 tag

   open(10,file='tec'//fname//'.dat',status='unknown')
      write(10,*)'TITLE = "',fname,'"'
      write(10,*)'VARIABLES = "xpos" "ypos" "depths" "',fname,'"'
      write(10,'(a,i4.4,a,i4.4,a)')' ZONE  F=BLOCK, I=',ii,', J=',jj,', K=1'


! Positions in lon-lat
      write(10,'(30I4)')((i,i=1,ii),j=1,jj)
      write(10,'(30I4)')((j,i=1,ii),j=1,jj)

! Depths
      write(10,900)((depths(i,j),i=1,ii),j=1,jj)

! Field
      write(10,900)((fld(i,j),i=1,ii),j=1,jj)

   close(10)
 900 format(10(1x,e12.5))

end subroutine tecfld

subroutine tecfldv(fname,fld,idim,jdim,deep,uvel,vvel)
   use mod_common
   implicit none
   character(len=*), intent(in) :: fname
   integer, intent(in) :: idim,jdim
   real, intent(in) :: fld(idim,jdim)
   real, intent(in) :: deep(idim,jdim)
   real, intent(in) :: uvel(idim,jdim),vvel(idim,jdim)
   integer i,j

   character*4 tag

   open(10,file='tec'//fname//'.dat',status='unknown')
      write(10,*)'TITLE = "',fname,'"'
      write(10,*)'VARIABLES = "xpos" "ypos" "lon" "lat" "deep" "',fname,'" "u" "v"'
      write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',idim,', J=',jdim,', K=1'


! Positions in lon-lat
      write(10,'(30I4)')((i,i=1,idim),j=1,jdim)
      write(10,'(30I4)')((j,i=1,idim),j=1,jdim)
      write(10,900)((modlon(i,j),i=1,idim),j=1,jdim)
      write(10,900)((modlat(i,j),i=1,idim),j=1,jdim)

! Depths
      write(10,900)((deep(i,j),i=1,idim),j=1,jdim)

! Field
      write(10,900)((fld(i,j),i=1,idim),j=1,jdim)
      write(10,900)((uvel(i,j),i=1,idim),j=1,jdim)
      write(10,900)((vvel(i,j),i=1,idim),j=1,jdim)
   close(10)
 900 format(10(1x,e12.5))

end subroutine tecfldv

subroutine tecfld2(fname,idim,jdim,deep,field)
   use mod_common
   implicit none
   character(len=*), intent(in) :: fname
   integer, intent(in) :: idim,jdim
   real, intent(in) :: deep(idim,jdim)
   real, intent(in) :: field(idim,jdim)
   integer i,j

   character*4 tag

   open(10,file=fname//'.dat',status='unknown')
      write(10,*)'TITLE = "',fname,'"'
      write(10,*)'VARIABLES = "i"  "j"  "lon" "lat" "deep" "field"'
      write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',idim,', J=',jdim,', K=1'


! Positions in Grid index
      write(10,'(30I4)')((i,i=1,idim),j=1,jdim)
      write(10,'(30I4)')((j,i=1,idim),j=1,jdim)


! Positions in Lon Lat
      write(10,900)((modlon(i,j),i=1,idim),j=1,jdim)
      write(10,900)((modlat(i,j),i=1,idim),j=1,jdim)
      write(10,900)((deep(i,j),i=1,idim),j=1,jdim)
      write(10,900)((field(i,j),i=1,idim),j=1,jdim)

   close(10)
900  format(10(1x,e12.5))

end subroutine tecfld2

subroutine tecfld3(fname,fld,idim,jdim)
   use mod_common
   implicit none
   character(len=*), intent(in) :: fname
   integer, intent(in) :: idim,jdim
   real, intent(in) :: fld(idim,jdim)
   integer i,j

   character*4 tag

   open(10,file='tec'//fname//'.dat',status='unknown')
      write(10,*)'TITLE = "',fname,'"'
      write(10,*)'VARIABLES = "i-index" "j-index" "',fname,'"'
      write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',idim,', J=',jdim,', K=1'

      write(10,'(30I4)')((i,i=1,idim),j=1,jdim)
      write(10,'(30I4)')((j,i=1,idim),j=1,jdim)
      write(10,900)((fld(i,j),i=1,idim),j=1,jdim)

   close(10)
 900 format(10(1x,e12.5))

end subroutine tecfld3

subroutine tecforce(frc)
   use mod_common
   use mod_atmforc
   implicit none
   type(forcing_fields), intent(in) :: frc
   character*4 tag
   integer idim,jdim,i,j

   idim=nx
   jdim=ny
    
   open(10,file='tecfrc.dat',status='unknown')
      write(10,*)'TITLE = "Forcing"'
      write(10,*)'VARIABLES = "i-index" "j-index" "lon" "lat" "depths"&
      & "uwind" "vwind" "taux" "tauy" "airt" "wndspd" "relhum" "clouds" "precip"&
      & "sss" "sst"'
      write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',idim,', J=',jdim,', K=1'

! Positions in lon-lat
      write(10,'(30I4)')((i,i=1,idim),j=1,jdim)
      write(10,'(30I4)')((j,i=1,idim),j=1,jdim)
      write(10,900)((modlon(i,j),i=1,idim),j=1,jdim)
      write(10,900)((modlat(i,j),i=1,idim),j=1,jdim)

! Depths
      write(10,900)((depths(i,j),i=1,idim),j=1,jdim)

! Fields
      write(10,900)((frc%uwind(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%vwind(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%taux(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%tauy(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%airtmp(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%wndspd(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%relhum(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%clouds(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%precip(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%sss(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%sst(i,j),i=1,idim),j=1,jdim)
   close(10)
 900 format(10(1x,e12.5))

end subroutine tecforce

subroutine tecaveforce(rungen,rt,frc)
   use mod_common
   use mod_atmforc
   use mod_year_info
   implicit none
   type(forcing_fields), intent(in) :: frc
   type(year_info), intent(inout) :: rt
   character(len=3), intent(inout) :: rungen
   character*4 tag
   integer idim,jdim,i,j

   idim=nx
   jdim=ny
    
   open(10,file=fileforce(rungen,rt),status='unknown')
!   open(10,file='FORCE.dat',status='unknown')
      write(10,*)'TITLE = "Forcing"'
      write(10,*)'VARIABLES = "i-index" "j-index" "lon" "lat" "depths"&
      & "taux" "tauy" "airt" "wndspd" "relhum" "clouds" "precip"&
      & "sss" "sst"'
      write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',idim,', J=',jdim,', K=1'

! Positions in lon-lat
      write(10,'(30I4)')((i,i=1,idim),j=1,jdim)
      write(10,'(30I4)')((j,i=1,idim),j=1,jdim)
      write(10,900)((modlon(i,j),i=1,idim),j=1,jdim)
      write(10,900)((modlat(i,j),i=1,idim),j=1,jdim)

! Depths
      write(10,900)((depths(i,j),i=1,idim),j=1,jdim)

! Fields
      write(10,900)((frc%taux(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%tauy(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%airtmp(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%wndspd(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%relhum(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%clouds(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%precip(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%sss(i,j),i=1,idim),j=1,jdim)
      write(10,900)((frc%sst(i,j),i=1,idim),j=1,jdim)
   close(10)
 900 format(10(1x,e12.5))

contains
      function fileforce(rungen,rt)
         use mod_dimensions
         implicit none
         character(len=19) fileforce
         type(year_info),  intent(in) :: rt
         character(len=3), intent(in) :: rungen
         fileforce(1:3)=rungen(1:3)
         fileforce(4:9)='FORCE_'
         fileforce(10:12)=rt%cdd(1:3)
         fileforce(13:16)=rt%cyy(1:4)
         fileforce(17:19)='.uf'
      end function fileforce
end subroutine tecaveforce

subroutine tecfldalt(fname,lat,lon,ii,jj,depths,field)
   implicit none
  character(len=*), intent(in) :: fname
  integer, intent(in) :: ii,jj
  real, intent(in) :: depths(ii,jj)
  real, intent(in) :: field(ii,jj)
  real, intent(in) :: lat(ii,jj)
  real, intent(in) :: lon(ii,jj)
  integer i,j

  character*4 tag

 open(10,file=fname//'.dat',status='unknown')
 write(10,*)'TITLE = "',fname,'"'
 write(10,*)'VARIABLES = "i"  "j"  "lon" "lat" "depths" "field"'
 write(10,'(a,i3,a,i3,a)')' ZONE  F=BLOCK, I=',ii,', J=',jj,', K=1'


! Positions in Grid index
 write(10,'(30I4)')((i,i=1,ii),j=1,jj)
 write(10,'(30I4)')((j,i=1,ii),j=1,jj)


! Positions in Lon Lat

 write(10,900)((lon(i,j),i=1,ii),j=1,jj)
 write(10,900)((lat(i,j),i=1,ii),j=1,jj)
 write(10,900)((depths(i,j),i=1,ii),j=1,jj)
 write(10,900)((field(i,j),i=1,ii),j=1,jj)

close(10)
900  format(10(1x,e12.5))

end subroutine tecfldalt
end module m_tecfld
