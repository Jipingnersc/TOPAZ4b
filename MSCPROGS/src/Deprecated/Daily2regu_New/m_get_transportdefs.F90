module m_get_transportdefs
contains

   subroutine get_transportdefs(transport_secnum,transport_classtype, &
      transport_classlimits, transport_direction,transport_name, &
      ntransport,nsec)
      implicit none

      character(len=6) , pointer :: transport_name(:)
      character(len=11), pointer :: transport_classtype(:)
      character(len=1) , pointer :: transport_direction(:)
      real,              pointer :: transport_classlimits(:,:)
      integer,           pointer :: transport_secnum(:)
      integer, intent(out) :: ntransport
      integer, intent(in) :: nsec

      character(len=*), parameter :: infile='transportdefs.in'

      integer          :: t1
      character(len=11):: t2
      real             :: t3, t4
      character(len=1) :: t5
      character(len=6) :: t6

      logical :: ex
      integer :: ios,itr,ntr



      ! First, get number of transport definitions
100   format(i6,a11,f7.1,f7.1,3x,a1,2x,a)
      inquire(exist=ex,file=infile)
      if (.not.ex) then
         print *,'Can not find file '//infile
         stop '(m_get_transportdefs.F90)'
      end if

      open(10,file=infile)
      ios=0
      ntr=0
      do while (ios==0)
         read(10,100,iostat=ios) t1,t2,t3,t4,t5,t6
         if(ios==0) ntr=ntr+1
      end do

      if (ntr<=0) then
         print *,'No transport definitions in '//infile
         stop '(m_get_transportdefs.F90)'
      end if

      ntransport=ntr
      allocate(transport_name       (ntransport))
      allocate(transport_classtype  (ntransport))
      allocate(transport_direction  (ntransport))
      allocate(transport_classlimits(ntransport,2))
      allocate(transport_secnum     (ntransport))
      rewind(10)
      do itr=1,ntransport
         read(10,100,iostat=ios) t1,t2,t3,t4,t5,t6
         write(*,100) t1,adjustl(t2),t3,t4,t5,t6

         transport_secnum     (itr)  =t1
         transport_classtype  (itr)  =t2
         transport_classlimits(itr,1)=t3
         transport_classlimits(itr,2)=t4
         transport_direction  (itr)  =t5
         transport_name       (itr)  =t6

      end do


      ! Consistency check: 
      if (maxval(transport_secnum(:))>nsec) then
         print *,'Section number specified in '//infile//' is higher than '// &
                 'available sections'
         print *,nsec,(maxval(transport_secnum(:)))
         stop '(m_get_transportdefs.F90)'
      end if


   end subroutine get_transportdefs
end module m_get_transportdefs






         


