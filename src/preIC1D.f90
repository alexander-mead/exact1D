PROGRAM preIC1D

  USE random_numbers
  USE simulations1D

  IMPLICIT NONE
  INTEGER :: i, n
  INTEGER :: imode, iseed
  REAL :: L
  REAL, ALLOCATABLE :: x(:), v(:)
  CHARACTER(len=256) :: outfile

  LOGICAL, PARAMETER :: verbose=.TRUE.

  !Starting whitespace
  WRITE(*,*)

  CALL get_command_argument(1,outfile)
  IF(outfile=='') STOP 'PREIC1D: Error, specify an output file'

  !Choose particle configuration
  WRITE(*,*) 'PREIC1D: Choose particle generation method'
  WRITE(*,*) '1 - Random'
  WRITE(*,*) '2 - Grid'
  READ(*,*) imode
  WRITE(*,*)
  
  !Number of particles to generate
  WRITE(*,*) 'PREIC1D: Number of particles to generate'
  READ(*,*) n
  WRITE(*,*)

  !Box length
  WRITE(*,*) 'PREIC1D: Box length'
  READ(*,*) L
  WRITE(*,*)

  !Allocate position and velocity arrays
  ALLOCATE(x(n),v(n))

  !Always set velocity to zero in preICs
  v=0.

  !Create the preIC positions
  IF(imode==1) THEN
     WRITE(*,*) 'PREIC: Choose random seed (0 sets via clock)'
     READ(*,*) iseed
     WRITE(*,*)
     CALL RNG_set(iseed)
     !DO i=1,n
     !   x(i)=random_uniform(0.,L)
     !END DO
     CALL make_randoms(x,n,L)
  ELSE IF(imode==2) THEN
     !DO i=1,n
     !   x(i)=L*(REAL(i)-0.5)/REAL(n)
     !END DO
     CALL make_grid(x,n,L)
  ELSE
     STOP 'PREIC1D: Error, imode specified incorrectly'
  END IF

  !Write particle data
  CALL write_particles(x,v,n,outfile,verbose)
  
CONTAINS

!!$  SUBROUTINE write_particles(x,v,n,outfile)
!!$
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: x(n), v(n)
!!$    INTEGER, INTENT(IN) :: n
!!$    CHARACTER(len=256), INTENT(IN) :: outfile
!!$    INTEGER :: i
!!$
!!$    WRITE(*,*) 'WRITE_PARTICLES: Outputting particle data'
!!$    WRITE(*,*) 'WRITE_PARTICLES: File ', TRIM(outfile)
!!$    
!!$    OPEN(7,file=outfile)
!!$    DO i=1,n
!!$       WRITE(7,*) x(i), v(i)
!!$    END DO
!!$    CLOSE(7)
!!$
!!$    WRITE(*,*) 'WRITE_PARTICLES: Done'
!!$    WRITE(*,*) 
!!$    
!!$  END SUBROUTINE write_particles

END PROGRAM preIC1D
