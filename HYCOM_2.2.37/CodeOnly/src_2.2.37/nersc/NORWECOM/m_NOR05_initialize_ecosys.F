      module m_NOR05_initialize_ecosys

      contains

      subroutine initialize_ecosys()
      use mod_xc
      use mod_za
      use mod_necessary_ecovars
      use m_NOR05_riverloads

      implicit none
      include 'common_blocks.h'
      integer :: i,j,k,o
      integer lslot
      character preambl(5)*79
      integer lgth
      real tiny
      !KAL - Should use existing margin, points up to ii+margin should
      !KAL - still be valid at this point 
      !margin=0
      tiny=2.23e-16

!    initialize nitrate from climatology files
      lgth=len_trim(flnmfor)
      call zaiopf(flnmforw(1:lgth)//'relax.nitr.a','old',931)
      if (mnproc.eq.1) then   ! .b file from 1st tile only
        open (unit=2931,file=flnmforw(1:lgth)//'relax.nitr.b',
     &      status='old',action='read')
        read(2931,'(a79)') preambl
      endif  !1st tile
      call preambl_print(preambl)

      do k=1,kdm
        call rdmonthck(bio(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,k,1,init),
     &                 931,1)
      enddo   
      close(931)
!     from molar to mgr:
      bio(:,:,:,1,init)=bio(:,:,:,1,init)*14.0
      bio(:,:,:,2,init)=bio(:,:,:,1,init)


!    initialize phosphate from climatology files
      lgth=len_trim(flnmfor)
      call zaiopf(flnmforw(1:lgth)//'relax.phos.a','old',932)
      if (mnproc.eq.1) then   ! .b file from 1st tile only
        open (unit=2932,file=flnmforw(1:lgth)//'relax.phos.b',
     &      status='old',action='read')
        read(2932,'(a79)') preambl
      endif  !1st tile
      call preambl_print(preambl)

      do k=1,kdm
        call rdmonthck(bio(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,k,1,ipho),
     &                 932,1)
      enddo   
      close(932)
!     from molar to mgr:      
      bio(:,:,:,1,ipho)=bio(:,:,:,1,ipho)*31.0
      bio(:,:,:,2,ipho)=bio(:,:,:,1,ipho)
      

!    initialize silicate from climatology files
      lgth=len_trim(flnmfor)
      call zaiopf(flnmforw(1:lgth)//'relax.sili.a','old',933)
      if (mnproc.eq.1) then   ! .b file from 1st tile only
        open (unit=2933,file=flnmforw(1:lgth)//'relax.sili.b',
     &      status='old',action='read')
        read(2933,'(a79)') preambl
      endif  !1st tile
      call preambl_print(preambl)

      do k=1,kdm
        call rdmonthck(bio(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,k,1,isil),
     &                 933,1)
      enddo   
      close(933)
!     from molar to mgr:
      bio(:,:,:,1,isil)=bio(:,:,:,1,isil)*28.0
      bio(:,:,:,2,isil)=bio(:,:,:,1,isil)

!    initialize oxygen from climatology files
      lgth=len_trim(flnmfor)
      call zaiopf(flnmforw(1:lgth)//'relax.oxyg.a','old',934)
      if (mnproc.eq.1) then   ! .b file from 1st tile only
        open (unit=2934,file=flnmforw(1:lgth)//'relax.oxyg.b',
     &      status='old',action='read')
        read(2934,'(a79)') preambl
      endif  !1st tile
      call preambl_print(preambl)

      do k=1,kdm
        call rdmonthck(bio(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,k,1,ioxy),
     &                 934,1)
      enddo   
      close(934)
!     from g/l to mgr:
      bio(:,:,:,1,ioxy)=bio(:,:,:,1,ioxy)*1330.0
      bio(:,:,:,2,ioxy)=bio(:,:,:,1,ioxy)


      ! KAL - initialization, so do the whole array
      do i=1-nbdy,idm+nbdy
      do j=1-nbdy,jdm+nbdy
      do k=1,kdm
      do o=1,nbio
         if(bio(i,j,k,1,o) .lt. 0.) then
            bio(i,j,k,1,o)=tiny
         endif   
         if(bio(i,j,k,2,o) .lt. 0.) then
            bio(i,j,k,2,o)=tiny
         endif
!Oxygene becomes negative when it is small   
         if (o==ioxy) then
            if (bio(i,j,k,1,o).le.10.0) then
               bio(i,j,k,1,o)=10.0
               bio(i,j,k,2,o)=10.0
            endif
         endif
      enddo   
      enddo   
      enddo   
      enddo  

      bio(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm,1:2,idet)=0.1
      bio(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm,1:2,isis)=0.1
      bio(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm,1:2,ifla)=0.1
      bio(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm,1:2,idia)=0.1 
#if defined DETPHO
      bio(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm,1:2,idetp)=0.1 
#endif
#if defined (ZOOPL)
      bio(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm,1:2,imeso)=0.1 
      bio(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm,1:2,imicro)=0.1 
#endif
      
      biodt_sou=baclin!AS*nphys_sou

! initialize river loads:
      call NOR05_readnutrloads()

      if (mod(86400,int(biodt_sou)) .ne. 0) then
         print *, 'mod_necessary_ecovars:'
         print *, 'baclin*nphys_sou should add up to 60*60*24=86400 s'
         print *, 'please redefine nphys_sou to ensure this...'
         stop
      else
         if (mnproc==1) print *, 'biosteps ok...'
      endif

      end subroutine initialize_ecosys
      end module m_NOR05_initialize_ecosys
