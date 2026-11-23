      program tst

      integer iday, iyear, ihour, yrflag
      real*8  dtime

      yrflag=3
      dtime=365.5d0
      call forday(dtime,yrflag,iyear,iday,ihour)
      print *,iyear, iday, ihour
      call dayfor(dtime,yrflag,iyear,iday,ihour)
      print *,dtime

      end program


