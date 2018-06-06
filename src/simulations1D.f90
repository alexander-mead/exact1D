MODULE simulations1D

  IMPLICIT NONE
  
CONTAINS

   SUBROUTINE read_particles(x,v,n,infile,verbose)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: x(n), v(n)
    INTEGER, INTENT(IN) :: n
    CHARACTER(len=256), INTENT(IN) :: infile
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i

    IF(verbose) THEN
       WRITE(*,*) 'READ_PARTICLES: Inputting particle data'
       WRITE(*,*) 'READ_PARTICLES: File ', TRIM(infile)
    END IF

    x=0.
    v=0.

    OPEN(7,file=infile)
    DO i=1,n
       READ(7,*) x(i), v(i)
    END DO
    CLOSE(7)

    IF(verbose) THEN
       WRITE(*,*) 'READ_PARTICLES: Done'
       WRITE(*,*) 
    END IF
      
  END SUBROUTINE read_particles

  SUBROUTINE write_particles(x,v,n,outfile,verbose)

    IMPLICIT NONE
    REAL, INTENT(IN) :: x(n), v(n)
    INTEGER, INTENT(IN) :: n
    CHARACTER(len=256), INTENT(IN) :: outfile
    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i

    IF(verbose .EQV. .TRUE.) THEN
       WRITE(*,*) 'WRITE_PARTICLES: Outputting particle data'
       WRITE(*,*) 'WRITE_PARTICLES: File ', TRIM(outfile)
    END IF

    OPEN(7,file=outfile)
    DO i=1,n
       WRITE(7,*) x(i), v(i)
    END DO
    CLOSE(7)

    IF(verbose .EQV. .TRUE.) THEN
       WRITE(*,*) 'WRITE_PARTICLES: Done'
       WRITE(*,*)
    END IF

  END SUBROUTINE write_particles

  SUBROUTINE make_randoms(x,n,L)

    USE random_numbers

    !Makes a uniform grid of n particles on a line of length L
    IMPLICIT NONE
    REAL, INTENT(OUT) :: x(n)
    REAL, INTENT(IN) :: L
    INTEGER, INTENT(IN) :: n
    INTEGER :: i

    DO i=1,n
       x(i)=random_uniform(0.,L)
    END DO
    
  END SUBROUTINE make_randoms

  SUBROUTINE make_grid(x,n,L)

    !Makes a uniform grid of n particles on a line of length L
    IMPLICIT NONE
    REAL, INTENT(OUT) :: x(n)
    REAL, INTENT(IN) :: L
    INTEGER, INTENT(IN) :: n
    INTEGER :: i

    DO i=1,n
       x(i)=L*(REAL(i)-0.5)/REAL(n)
    END DO
    
  END SUBROUTINE make_grid

  SUBROUTINE enforce_periodicity(x,n,L)

    !Enfoce the periodicity of the line
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: x(n)
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: L
    INTEGER :: i, m

    DO i=1,n
       
       IF(x(i)<0.) THEN
          m=CEILING(-x(i)/L)
          x(i)=x(i)+REAL(m)*L
       ELSE IF(x(i)>L) THEN
          m=CEILING(x(i)/L)-1
          x(i)=x(i)-REAL(m)*L
       END IF

       IF(x(i)==0.) x=L

    END DO
    
  END SUBROUTINE enforce_periodicity

  SUBROUTINE k_FFT1D(i,m,kx,k,L)

     USE constants
     
    !Compute the k vector associated with positions in the Fourier mesh
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, m
    REAL*8, INTENT(OUT) :: kx, k
    REAL*8, INTENT(IN) :: L 

    !I could change it to work for odd mesh...
    !Note that this would be slightly non-trivial because the middle mode is missing
    IF(MOD(m,2) .NE. 0) STOP 'K_FFT1D: Mesh size needs to be even'

    kx=float(i-1)
    IF(i>m/2+1) kx=-float(m-i+1)
    kx=kx*2.*pi/L

    k=ABS(kx)

  END SUBROUTINE k_FFT1D

END MODULE simulations1D
