PROGRAM exact1D

  USE random_numbers
  USE numerology
  USE array_operations
  USE sorting
  
  IMPLICIT NONE
  REAL, ALLOCATABLE :: x0(:), v0(:), a0(:), t0(:), tc(:)
  REAL, ALLOCATABLE :: t(:)
  INTEGER, ALLOCATABLE :: id(:)
  REAL :: E, K, V, E1, E2, Vir
  REAL :: ts, dx, G
  REAL :: vbar
  INTEGER :: i, j, m
  CHARACTER(len=256) :: infile, outfile

  ! Parameters
  LOGICAL, PARAMETER :: use_restart=.FALSE. ! Use previous simulation end
  LOGICAL, PARAMETER :: ran_offsets=.TRUE. ! Do random offsetting
  REAL, PARAMETER :: alpha=1. ! Amount to offset in terms of spacing (1 is maximum)
  LOGICAL, PARAMETER :: ran_velocities=.FALSE. ! Do random velocties
  LOGICAL, PARAMETER :: verbose=.TRUE. ! Set the speaking level
  LOGICAL, PARAMETER :: extra_verbose=.FALSE. ! Set the extra speaking level
  INTEGER :: n=3 ! Number of sheets
  REAL, PARAMETER :: L=0.5 ! Size of inital group
  REAL, PARAMETER :: t1=0. ! Initial time (is this necessary?)
  REAL, PARAMETER :: t2=10. ! Final time
  INTEGER, PARAMETER :: s=500 ! Number of time outputs
  INTEGER, PARAMETER :: iseed=1 ! Random-number seed
  !REAL :: G=1./REAL(n) ! Gravity strength (such that total mass is independent of number of sheets)

  ! Potential speed ups:
  ! Do not need array for the acceleration, could just use position info

  ! Potential problems:
  ! Pairs of sheets that collide at exactly the same time
  ! More than two sheets collide at exactly the same position
  ! Quantised positions and velocities due to numerical precision
  
  WRITE(*,*)
  WRITE(*,*) 'EXACT1D: 1D particle calculation'

  IF(use_restart) THEN

     infile='data/end.dat'
     CALL read_restart(x0,v0,a0,id,n,infile)

     G=1./REAL(n)

  ELSE
     
     CALL RNG_set(iseed)

     ! Set the number of sheets and allocate arrays
     WRITE(*,*) 'EXACT1D: Number of sheets:', n
     ALLOCATE(x0(n),v0(n),a0(n),id(n))
     WRITE(*,*)

     ! Fill the initial position array
     CALL fill_array(-L,L,x0,n)

     ! Do some random offsets
     IF(ran_offsets) CALL random_offsets(x0,L,n)

     ! Check that the sheets are arranged correctly
     DO i=1,n-1
        IF(x0(i)>x0(i+1)) STOP 'EXACT1D: Error, sheets are not initially aligned in x'
     END DO

     ! Set the location array (sheets must initially be organised in position)
     DO i=1,n
        id(i)=i
     END DO

     ! Set the velocities
     IF(ran_velocities) THEN
        ! Assign some random velocities
        DO i=1,n
           v0(i)=random_uniform(-1.,1.)
        END DO
        vbar=SUM(v0)/REAL(n) ! Calculate the mean velocity
        v0=v0-vbar ! Ensure that the mean velocity is zero
     ELSE
        ! Set the initial velocities to zero
        v0=0.
     ENDIF

     ! Set the initial accelerations, which only depend on the particle position
     ! Acceleration array does not change EVER if we keep the array order to be the positional order
     G=1./REAL(n)
     CALL fill_array((n-1)*G,-(n-1)*G,a0,n)     

  END IF

  ! Set the 'initial times' for all sheets to be t1
  ALLOCATE(t0(n))
  t0=t1

  ! Write initial conditions
  outfile='data/start.dat'
  CALL write_restart(t1,x0,v0,a0,t0,n,outfile)

  ! Write particle information to the screen
  IF(extra_verbose) THEN
     WRITE(*,*) '======================================='
     WRITE(*,*) '        i         x         v         a'
     WRITE(*,*) '======================================='
     DO i=1,n
        WRITE(*,fmt='(I10,3F10.5)') i, x0(i), v0(i), a0(i)
     END DO
     WRITE(*,*) '======================================='
     WRITE(*,*)
  END IF

  ! Fill the time array
  CALL fill_array(t1,t2,t,s)
  WRITE(*,*) 'EXACT1D: Initial time:', t1
  WRITE(*,*) 'EXACT1D: Final time:', t2
  WRITE(*,*) 'EXACT1D: Number of output times:', s
  WRITE(*,*)

  ! Calculate the initial collision times between all neighbouring pairs (n-1)
  ALLOCATE(tc(n-1))
  IF(extra_verbose) THEN
     WRITE(*,*) '======================================='
     WRITE(*,*) '        i       i+1                  tc'
     WRITE(*,*) '======================================='
  END IF 
  DO j=1,n-1  
     tc(j)=collision_time(x0(j),v0(j),a0(j),t0(j),x0(j+1),v0(j+1),a0(j+1),t0(j+1),t1)
     IF(extra_verbose) WRITE(*,fmt='(2I10,ES20.10)') j, j+1, tc(j)
  END DO
  IF(extra_verbose) THEN
     WRITE(*,*) '======================================='
     WRITE(*,*)
  END IF

  ! Open files for data
  OPEN(7,file='data/output.dat')
  OPEN(8,file='data/energy.dat')

  ! Make a nice display grid
  IF(verbose) THEN
     WRITE(*,*) '======================================================'
     WRITE(*,*) '     step           time          E=K+T            T/V'
     WRITE(*,*) '======================================================'
  END IF

  ! Calculate the initial energy of the system
  E1=kinetic_energy(t1,v0,a0,t0,n)+potential_energy(t1,x0,v0,a0,t0,n)

  ! Write data to file/screen
  CALL write_data(t(1),x0,v0,a0,t0,id,n)
  CALL write_energy(1,t(1),x0,v0,a0,t0,n)  

  ! Loop over output times
  DO i=1,s-1

     ! Calculate the next collision time
     ! This might be slightly wasteful. It could probably be ordered somehow
     IF(i==1) ts=MINVAL(tc)

     ! If the next collision time is before the output time then update the system
     DO WHILE (ts<t(i+1))
        CALL collision_update(ts,tc,x0,v0,a0,t0,id,n)
        ts=MINVAL(tc)
     END DO

     ! Otherwise write data for file/screen when the output time is before the next collision
     CALL write_data(t(i+1),x0,v0,a0,t0,id,n)
     CALL write_energy(i+1,t(i+1),x0,v0,a0,t0,n)
     
  END DO

  ! Close files for data
  CLOSE(7)
  CLOSE(8)

  ! Write final conditions
  outfile='data/end.dat'
  CALL write_restart(t2,x0,v0,a0,t0,n,outfile)

  ! Calculate the final energy
  E2=kinetic_energy(t2,v0,a0,t0,n)+potential_energy(t2,x0,v0,a0,t0,n)

  ! Write a nice grid to screen
  IF(verbose) THEN
     WRITE(*,*) '======================================================'
     WRITE(*,*)
  END IF

  ! Check that energy conservation is okay
  WRITE(*,*) 'EXACT1D: Initial energy:', E1
  WRITE(*,*) 'EXACT1D: Final energy:', E2
  WRITE(*,*) 'EXACT1D: Ratio of final to initial energy:', E2/E1
  WRITE(*,*)
  
CONTAINS

  SUBROUTINE collision_update(t,tc,x0,v0,a0,t0,id,n)

    IMPLICIT NONE
    REAL, INTENT(INOUT) :: x0(n), v0(n), t0(n), tc(n-1)
    REAL, INTENT(IN) :: a0(n), t
    INTEGER, INTENT(INOUT) :: id(n)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, m
    
    LOGICAL, PARAMETER :: check_collisions=.FALSE.

    ! Find which sheets are colliding
    m=MINLOC(tc,1)

    ! Check to see if more than one sheets are colliding
    IF(check_collisions) THEN
       DO i=1,n-1
          IF(tc(i)==t) WRITE(*,*) i, m, t, tc(i)
       END DO
    END IF

    ! Update the positions of the colliding sheets to be the collision position
    x0(m)=position(t,x0(m),v0(m),a0(m),t0(m))
    x0(m+1)=x0(m) ! Ensure that these are exactly the same (numerical roundoff)

    ! Swap the sheets so that array position corresponds to physical position
    ! DO NOT swap accelerations, as these depend on position only
    v0(m)=velocity(t,v0(m),a0(m),t0(m))
    v0(m+1)=velocity(t,v0(m+1),a0(m+1),t0(m+1))
    CALL swap_real(v0(m),v0(m+1))

    ! Fix the new t0 for these sheets to be the collision time
    t0(m)=t
    t0(m+1)=t

    ! Swap the ID numbers so as to keep track of the ordering of sheets
    CALL swap_int(id(m),id(m+1))

    ! Update the collision times
    DO i=m-1,m+1
       IF(i>=1 .AND. i<=n-1) THEN
          tc(i)=collision_time(x0(i),v0(i),a0(i),t0(i),x0(i+1),v0(i+1),a0(i+1),t0(i+1),t)
       END IF
    END DO
    
  END SUBROUTINE collision_update

  SUBROUTINE random_offsets(x,L,n)

    ! Give sheets random offsets
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: x(n)
    REAL, INTENT(IN) :: L
    INTEGER, INTENT(IN) :: n
    REAL :: dx
    INTEGER :: i

    ! Set the offset to be the spacing between the sheets
    dx=alpha*L/REAL(n)

    ! Do the offsetting
    DO i=1,n
       x(i)=x(i)+random_uniform(-dx,dx)
    END DO

  END SUBROUTINE random_offsets

  FUNCTION position(t,x,v,a,t0)

    ! Provides the position of the sheet at time t
    IMPLICIT NONE
    REAL :: position
    REAL, INTENT(IN) :: x, v, a, t0, t
    
    position=0.5*a*(t-t0)**2+v*(t-t0)+x
    
  END FUNCTION position

  FUNCTION velocity(t,v,a,t0)

    ! Provides the velocity of the sheet at time t
    IMPLICIT NONE
    REAL :: velocity
    REAL, INTENT(IN) :: v, a, t0, t 

    velocity=a*(t-t0)+v
    
  END FUNCTION velocity

  FUNCTION kinetic_energy(t,v0,a0,t0,n)

    ! Calculates the total kinetic energy of the set of sheets
    IMPLICIT NONE
    REAL :: kinetic_energy
    REAL, INTENT(IN) :: v0(n), a0(n), t0(n)
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: t
    REAL :: K, v
    INTEGER :: i

    ! This variable is used for sum
    K=0.

    ! Loop over all sheets and calculate kinetic energy of each
    DO i=1,n
       v=velocity(t,v0(i),a0(i),t0(i))
       K=K+0.5*v**2
    END DO

    ! Set the result
    kinetic_energy=K
    
  END FUNCTION kinetic_energy

  FUNCTION potential_energy(t,x0,v0,a0,t0,n)

    ! Calculates the total potential energy of the set of sheets
    IMPLICIT NONE
    REAL :: potential_energy
    REAL, INTENT(IN) :: x0(n), v0(n), a0(n), t0(n)
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: t
    REAL :: P, x(n)
    INTEGER :: i, j

    ! This variable is used for sum
    P=0.

    ! Loop over sheets and calculate position of each
    DO i=1,n
       x(i)=position(t,x0(i),v0(i),a0(i),t0(i))
    END DO

    ! Do a double-loop to calculate the total potential energy
    DO j=1,n
       DO i=j+1,n
          P=P+G*ABS(x(j)-x(i))
       END DO
    END DO

    ! Set the result
    potential_energy=P
    
  END FUNCTION potential_energy

  FUNCTION collision_time(x1,v1,a1,t1,x2,v2,a2,t2,t)

    ! Calculates the collision time between pair labelled '1' and '2'
    IMPLICIT NONE
    REAL :: collision_time
    REAL, INTENT(IN) :: x1, v1, a1, t1
    REAL, INTENT(IN) :: x2, v2, a2, t2
    REAL, INTENT(IN) :: t
    REAL :: discrim, tc1, tc2, tc, ts
    REAL :: A, B, C
    
    LOGICAL, PARAMETER :: cases=.FALSE.
    REAL, PARAMETER :: small=1d-14

    ! A, B and C for the quadratic formula
    A=0.5*(a1-a2)
    B=(v1-v2)-(a1*t1-a2*t2)
    C=(x1-x2)-(v1*t1-v2*t2)+(0.5*a1*t1**2-0.5*a2*t2**2)

    ! Quadratic discriminant
    discrim=B**2-4.*A*C

    ! Make the time slightly larger to accomodate round-off errors
    ! Is this necessary?
    ts=t!*(1.d0+small)

    ! Now do things based on the value of the discriminant
    IF(discrim<0.) THEN
       ! No solution (does this even happen?)
       collision_time=HUGE(tc)
       IF(cases) WRITE(*,*) 'Case 1: Discriminant is negative, no solution'
       STOP 'COLLISION_TIME: Case 1: Discriminant is negative, no solution'
    ELSE IF(discrim==0. .AND. a1==a2) THEN
       ! Solution is 'at infinity'
       collision_time=HUGE(tc)
       IF(cases) WRITE(*,*) 'Case 2: Solution is at infinifty'
       STOP 'COLLISION_TIME: Case 2: Solution is at infinifty'
    ELSE IF(discrim==0.) THEN
       ! Single solution
       tc=-0.5*B/A
       IF(tc<=ts) THEN
          collision_time=HUGE(tc)
          IF(cases) WRITE(*,*) 'Case 3: ?'
          STOP 'COLLISION_TIME: Case 3: ?'
       ELSE
          collision_time=tc
          IF(cases) WRITE(*,*) 'Case 4: Single solution'
       END IF
    ELSE
       ! Two solutions, choose the minimum one (I think tc1 is always > tc2)
       tc1=-B+sqrt(discrim)
       tc2=-B-sqrt(discrim)
       tc1=0.5*tc1/A
       tc2=0.5*tc2/A
       IF(tc1<=ts .AND. tc2<=ts) THEN
          ! Both solutions in the past (should this ever happen?)
          collision_time=HUGE(tc) 
          IF(cases) WRITE(*,*) 'Case 5: Both solutions are in the past'
          STOP 'COLLISION_TIME: Case 5: Both solutions are in the past'
       ELSE IF(tc1>ts .AND. tc2>ts) THEN
          collision_time=tc1
          IF(cases) WRITE(*,*) 'Case 6: Two solutions, picked minimum'
       ELSE IF(tc1>ts .AND. tc2<=ts) THEN
          collision_time=tc1
          IF(cases) WRITE(*,*) 'Case 7: Two solutions, picked minimum'
       ELSE IF(tc2>ts .AND. tc1<=ts) THEN
          collision_time=tc2
          IF(cases) WRITE(*,*) 'Case 8: Two solutions, picked minimum'
       ELSE
          STOP 'COLLISION_TIME: Error, quadratic not solved correctly'
       END IF
    END IF

    IF(collision_time==HUGE(tc)) STOP 'COLLISION_TIME: Error, this should never happen'
    
  END FUNCTION collision_time

  SUBROUTINE read_restart(x0,v0,a0,id,n,infile)

    USE file_info

    ! Read in a restart file
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(OUT) :: x0(:), v0(:), a0(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: id(:)
    INTEGER, INTENT(OUT) :: n
    CHARACTER(len=*), INTENT(IN) :: infile
    INTEGER :: i
    REAL :: x, v, a

    n=file_length(infile)
    ALLOCATE(x0(n),v0(n),a0(n),id(n))

    OPEN(7,file=infile)
    DO i=1,n
       READ(7,*) x0(i), v0(i), a0(i)
       id(i)=i
    END DO
    CLOSE(7)

  END SUBROUTINE read_restart

  SUBROUTINE write_restart(t,x0,v0,a0,t0,n,outfile)

    ! Figure out where the sheets are and write to file
    IMPLICIT NONE
    REAL, INTENT(IN) :: t
    REAL, INTENT(IN) :: x0(n), v0(n), a0(n), t0(n)
    INTEGER, INTENT(IN) :: n
    CHARACTER(len=*), INTENT(IN) :: outfile
    INTEGER :: i
    REAL :: x, v, a

    OPEN(7,file=outfile)
    DO i=1,n
       x=position(t,x0(i),v0(i),a0(i),t0(i))
       v=velocity(t,v0(i),a0(i),t0(i))
       a=a0(i)
       WRITE(7,*) x, v, a
    END DO
    CLOSE(7)

  END SUBROUTINE write_restart

  SUBROUTINE write_data(t,x,v,a,t0,id,n)

    USE sorting

    ! Figure out where the sheets are in the arrays and write to file
    IMPLICIT NONE
    REAL, INTENT(IN) :: t
    REAL, INTENT(IN) :: x(n), v(n), a(n), t0(n)
    INTEGER, INTENT(IN) :: id(n)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j, idx(n)
    REAL :: xn(n), vn(n)

    LOGICAL, PARAMETER :: bubble=.FALSE.

    IF(bubble) THEN

       CALL bubble_index(id,idx,n)
       DO i=1,n
          j=idx(i)
          !j=i
          xn(i)=position(t,x(j),v(j),a(j),t0(j))
          vn(i)=velocity(t,v(j),a(j),t0(j))
       END DO
       WRITE(7,*) t, (xn(i), vn(i), i=1,n)
       
    ELSE

       ! This is definitely a stupid way to write out data
       DO i=1,n
          DO j=1,n
             IF(id(j)==i) THEN
                xn(i)=position(t,x(j),v(j),a(j),t0(j))
                vn(i)=velocity(t,v(j),a(j),t0(j))
                EXIT
             END IF
          END DO
       END DO

       ! Write to file
       WRITE(7,*) t, (xn(i), vn(i), i=1,n)

    END IF
         
  END SUBROUTINE write_data

  SUBROUTINE write_energy(i,t,x0,v0,a,t0,n)

    IMPLICIT NONE
    REAL, INTENT(IN) :: x0(n), v0(n), a(n), t0(n)
    INTEGER, INTENT(IN) :: i, n
    REAL, INTENT(IN) :: t
    REAL :: K, V, E, Vir

    ! Calculate the energies
    K=kinetic_energy(t,v0,a,t0,n)
    V=potential_energy(t,x0,v0,a,t0,n)
    E=K+V
    Vir=K/V

    ! Write the energies
    IF(verbose) WRITE(*,fmt='(I10,3F15.7)') i, t, E, Vir  
    WRITE(8,*) t, K, V, E, Vir
    
  END SUBROUTINE write_energy

END PROGRAM exact1D
