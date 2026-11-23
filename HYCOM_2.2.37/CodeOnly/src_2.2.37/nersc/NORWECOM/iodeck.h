! iodeck.h 
!
!  computer:   switch for computer dependant i/o methodes
!              ('DEC-Alpha' or something else)
!  reclunit:   unit in bytes for recl= in open statements
!
      character*9 computer
      integer     reclunit
!     parameter (computer='decalpha', reclunit=4)
!     parameter (computer='ibmrs6000', reclunit=1)
      parameter (computer='originSGI', reclunit=1)
