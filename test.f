      IMPLICIT REAL*8(A-H,O-Z)

      do 100 i=0,2000 
         write(6,*) -2.0d0+4.0d0/2000.0d0*i,
     $        DERF(-2.0d0+4.0d0/2000.0d0*i)
 100  continue
      end
