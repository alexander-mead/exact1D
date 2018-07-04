PROGRAM exact1D

  USE random_numbers
  USE array_operations
  
  IMPLICIT NONE
  REAL, ALLOCATABLE :: x0(:), v0(:), a(:), t0(:), tc(:)
  REAL, ALLOCATABLE :: t(:)
  INTEGER, ALLOCATABLE :: id(:)
  REAL :: g, E, K, V, E1, E2, Vir
  REAL :: ts, dx
  INTEGER :: i, j, m

  !Parameters
  LOGICAL, PARAMETER :: iran=.TRUE. !Set the random number generator
  LOGICAL, PARAMETER :: verbose=.TRUE. !Set the speaking level
  INTEGER, PARAMETER :: n=101 !Number of sheets
  REAL, PARAMETER :: L=0.5 !Size of inital group
  REAL, PARAMETER :: t1=0. !Initial time (is this necessary?)
  REAL, PARAMETER :: t2=10. !Final time
  INTEGER, PARAMETER :: s=1001 !Number of time outputs
  INTEGER, PARAMETER :: iseed=0 !Random-number seed

  !Potential speed ups:
  !Do not need array for the acceleration, can just use position info

  !Potential problems:
  !Pairs of sheets that collide at exactly the same time
  !More than two sheets collide at exactly the same position 
  
  WRITE(*,*)
  WRITE(*,*) 'EXACT1D: 1D particle calculation'

  !Set the number of sheets and allocate arrays
  WRITE(*,*) 'EXACT1D: Number of sheets:', n
  ALLOCATE(x0(n),v0(n),a(n),t0(n),tc(n-1),id(n))
  WRITE(*,*)

  !Fill the initial position array
  CALL fill_array(-L,L,x0,n)

  !Do some random offsets
  IF(iran) THEN
     CALL RNG_set(iseed)
     dx=L/REAL(n)
     DO i=1,n
        !WRITE(*,*) i, random_uniform(REAL(-dx),REAL(dx))
        x0(i)=x0(i)+random_uniform(-dx,dx)
     END DO
     !STOP 'There is a problem with the random-number generator when upgrading the real precision to 8 or 16'
  END IF

  !Check that the sheets are arranged correctly
  DO i=1,n-1
     IF(x0(i)>x0(i+1)) STOP 'EXACT1D: Error, sheets are not initially aligned in x'
  END DO

  !Set the location array (sheets must initially be organised in position)
  DO i=1,n
     id(i)=i
  END DO

  !Set the initial velocities to zero
  v0=0.

  !Set the initial accelerations, which only depend on the particle position
  !These do not change ever if we keep the array order to be the positional order
  g=1./REAL(n)
  CALL fill_array((n-1)*g,-(n-1)*g,a,n)

  !Write particle information to the screen
  !IF(verbose) THEN
  !   WRITE(*,*) '======================================='
  !   WRITE(*,*) '        i         x         v         a'
  !   WRITE(*,*) '======================================='
  !   DO i=1,n
  !      WRITE(*,fmt='(I10,3F10.5)') i, x0(i), v0(i), a(i)
  !   END DO
  !   WRITE(*,*) '======================================='
  !   WRITE(*,*)
  !END IF

  !Fill the time array
  CALL fill_array(t1,t2,t,s)
  WRITE(*,*) 'EXACT1D: Initial time:', t1
  WRITE(*,*) 'EXACT1D: Initial time:', t2
  WRITE(*,*) 'EXACT1D: Number of output times:', s
  WRITE(*,*)

  !Set the 'initial times' for all sheets to be t1
  t0=t1

  !Calculate the initial collision times between all neighbouring pairs (n-1)
  !IF(verbose) THEN
  !   WRITE(*,*) '======================================='
  !   WRITE(*,*) '        i       i+1                  tc'
  !   WRITE(*,*) '======================================='
  !END IF
  DO j=1,n-1  
     tc(j)=collision_time(x0(j),v0(j),a(j),t0(j),x0(j+1),v0(j+1),a(j+1),t0(j+1),t1)
     !IF(verbose) WRITE(*,fmt='(2I10,ES20.10)') j, j+1, tc(j)
  END DO
  !IF(verbose) THEN
  !   WRITE(*,*) '======================================='
  !   WRITE(*,*)
  !END IF
  
  OPEN(7,file='data/output.dat')
  OPEN(8,file='data/energy.dat')
  IF(verbose) THEN
     WRITE(*,*) '======================================================'
     WRITE(*,*) '     step           time          E=K+T            T/V'
     WRITE(*,*) '======================================================'
  END IF
  K=kinetic_energy(v0,a,t0,n,t(1))
  V=potential_energy(x0,v0,a,t0,n,g,t(1))
  E=K+V
  E1=E
  Vir=K/V
  IF(verbose) WRITE(*,fmt='(I10,3F15.7)') 1, t(1), E, Vir
  WRITE(8,*) t(1), K, V, E, Vir
  CALL write_data(t1,x0,v0,a,t0,id,n)
  DO i=1,s-1
212  ts=MINVAL(tc)
     IF(ts<=t(i+1)) THEN
        !There is a collision between sheet k and k+1
        m=MINLOC(tc,1)
        x0(m)=position(x0(m),v0(m),a(m),t0(m),ts)
        x0(m+1)=x0(m) !Ensure these are exactly the same (numerical roundoff)
        !Swap the sheets so that array position corresponds to physical position
        !DO NOT swap accelerations, as these depend on position only
        v0(m)=velocity(v0(m),a(m),t0(m),ts)
        v0(m+1)=velocity(v0(m+1),a(m+1),t0(m+1),ts)
        CALL swap_real(v0(m),v0(m+1))
        !Fix the new t0 to be the collision time
        t0(m)=ts
        t0(m+1)=ts
        !Swap the ID numbers so as to keep track of sheets
        CALL swap_int(id(m),id(m+1))
        DO j=m-1,m+1
           IF(j>=1 .AND. j<=n-1) THEN
              tc(j)=collision_time(x0(j),v0(j),a(j),t0(j),x0(j+1),v0(j+1),a(j+1),t0(j+1),ts)
           END IF
        END DO
        GOTO 212
     ELSE
        !Output the positions
        K=kinetic_energy(v0,a,t0,n,t(i+1))
        V=potential_energy(x0,v0,a,t0,n,g,t(i+1))
        E=K+V
        Vir=K/V
        IF(verbose) WRITE(*,fmt='(I10,3F15.7)') i+1, t(i+1), E, Vir
        WRITE(8,*) t(i+1), K, V, E, Vir
        CALL write_data(t(i+1),x0,v0,a,t0,id,n)
     END IF
  END DO
  CLOSE(7)
  CLOSE(8)
  E2=E
  IF(verbose) THEN
     WRITE(*,*) '======================================================'
     WRITE(*,*)
  END IF

  WRITE(*,*) 'EXACT1D: Initial energy:', E1
  WRITE(*,*) 'EXACT1D: Final energy:', E2
  WRITE(*,*) 'EXACT1D: Energy conservation:', -1.+E2/E1
  WRITE(*,*)
  
CONTAINS

  FUNCTION position(x,v,a,t0,t)

    !Provides the position of the sheet at time t
    IMPLICIT NONE
    REAL :: position
    REAL, INTENT(IN) :: x, v, a, t0, t
    
    position=0.5*a*(t-t0)**2+v*(t-t0)+x
    
  END FUNCTION position

  FUNCTION velocity(v,a,t0,t)

    !Provides the velocity of the sheet at time t
    IMPLICIT NONE
    REAL :: velocity
    REAL, INTENT(IN) :: v, a, t0, t 

    velocity=a*(t-t0)+v
    
  END FUNCTION velocity

  FUNCTION kinetic_energy(v0,a,t0,n,t)

    !Calculates the total kinetic energy of the set of sheets
    IMPLICIT NONE
    REAL :: kinetic_energy
    REAL, INTENT(IN) :: v0(n), a(n), t0(n)
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: t
    REAL :: K, v
    INTEGER :: i

    !These variables will be used for sums
    K=0.

    !Loop over all sheets and calculate kinetic energy of each
    DO i=1,n
       v=velocity(v0(i),a(i),t0(i),t)
       K=K+0.5*v**2
    END DO
    
    kinetic_energy=K
    
  END FUNCTION kinetic_energy

  FUNCTION potential_energy(x0,v0,a,t0,n,g,t)

    !Calculates the total potential energy of the set of sheets
    IMPLICIT NONE
    REAL :: potential_energy
    REAL, INTENT(IN) :: x0(n), v0(n), a(n), t0(n)
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: t, g
    REAL :: P, x(n)
    INTEGER :: i, j

    !These variables will be used for sums
    P=0.

    !Loop over all sheets and calculate energy of each
    DO i=1,n
       x(i)=position(x0(i),v0(i),a(i),t0(i),t)
    END DO

    DO j=1,n
       DO i=j+1,n
          P=P+g*ABS(x(j)-x(i))
       END DO
    END DO

    potential_energy=P
    
  END FUNCTION potential_energy

  SUBROUTINE write_data(t,x,v,a,t0,id,n)

    !Figure out where the sheets are in the arrays and write to file
    IMPLICIT NONE
    REAL, INTENT(IN) :: t
    REAL, INTENT(IN) :: x(n), v(n), a(n), t0(n)
    INTEGER, INTENT(IN) :: id(n)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i, j
    REAL :: xn(n)
    REAL :: vn(n)

    !This is definitely a stupid way to write out data
    DO i=1,n
       DO j=1,n
          IF(id(j)==i) THEN
             xn(i)=position(x(j),v(j),a(j),t0(j),t)
             vn(i)=velocity(v(j),a(j),t0(j),t)
             CYCLE
          END IF
       END DO
    END DO

    !Write to file
    WRITE(7,*) t, (xn(j), vn(j), j=1,n)
       
  END SUBROUTINE write_data

  SUBROUTINE swap_real(a,b)

    !Swaps the values of variables a and b
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: a, b
    REAL :: c

    c=a
    a=b
    b=c
    
  END SUBROUTINE swap_real

  SUBROUTINE swap_int(a,b)

    !Swaps the values of variables a and b
    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: a, b
    INTEGER :: c    

    c=a
    a=b
    b=c
    
  END SUBROUTINE swap_int

  FUNCTION collision_time(x1,v1,a1,t1,x2,v2,a2,t2,t)

    !Calculates the collision time between pair labelled '1' and '2'
    IMPLICIT NONE
    REAL :: collision_time
    REAL, INTENT(IN) :: x1, v1, a1, t1
    REAL, INTENT(IN) :: x2, v2, a2, t2
    REAL, INTENT(IN) :: t
    REAL :: discrim, tc1, tc2, tc, ts
    REAL :: A, B, C
    LOGICAL, PARAMETER :: cases=.FALSE.
    REAL, PARAMETER :: small=1d-14

    !A, B and C for the quadratic formula
    A=0.5*(a1-a2)
    B=(v1-v2)-(a1*t1-a2*t2)
    C=(x1-x2)-(v1*t1-v2*t2)+(0.5*a1*t1**2-0.5*a2*t2**2)

    !Quadratic discriminant
    discrim=B**2-4.*A*C

    !Make the time slightly larger to accomodate round-off errors
    !Is this necessary?
    ts=t!*(1.d0+small)

    !Now do things based on the value of the discriminant
    IF(discrim<0.) THEN
       !No solution (does this even happen?)
       collision_time=HUGE(tc)
       IF(cases .EQV. .TRUE.) WRITE(*,*) 'Case 1'
       STOP 'COLLISION_TIME: Case 1'
    ELSE IF(discrim==0. .AND. a1==a2) THEN
       !Solution is 'at infinity'
       collision_time=HUGE(tc)
       IF(cases .EQV. .TRUE.) WRITE(*,*) 'Case 2'
       STOP 'COLLISION_TIME: Case 2'
    ELSE IF(discrim==0.) THEN
       !Single solution
       tc=-0.5*B/A
       IF(tc<=ts) THEN
          collision_time=HUGE(tc)
          IF(cases .EQV. .TRUE.) WRITE(*,*) 'Case 3'
          STOP 'COLLISION_TIME: Case 3'
       ELSE
          collision_time=tc
          IF(cases .EQV. .TRUE.) WRITE(*,*) 'Case 4'
       END IF
    ELSE
       !Two solutions, choose the minimum one (I think tc1 is always > tc2)
       tc1=-B+sqrt(discrim)
       tc2=-B-sqrt(discrim)
       tc1=0.5*tc1/A
       tc2=0.5*tc2/A
       IF(tc1<=ts .AND. tc2<=ts) THEN
          !Both solutions in the past (should this ever happen?)
          collision_time=HUGE(tc) 
          IF(cases .EQV. .TRUE.) WRITE(*,*) 'Case 5'
          STOP 'COLLISION_TIME: Case 5'
       ELSE IF(tc1>ts .AND. tc2>ts) THEN
          collision_time=tc1
          IF(cases .EQV. .TRUE.) WRITE(*,*) 'Case 6'
       ELSE IF(tc1>ts .AND. tc2<=ts) THEN
          collision_time=tc1
          IF(cases .EQV. .TRUE.) WRITE(*,*) 'Case 7'
       ELSE IF(tc2>ts .AND. tc1<=ts) THEN
          collision_time=tc2
          IF(cases .EQV. .TRUE.) WRITE(*,*) 'Case 8'
       ELSE
          STOP 'COLLISION_TIME: Error, quadratic not solved correctly'
       END IF
    END IF

    IF(collision_time==HUGE(tc)) STOP 'COLLISION_TIME: Error, this should never happen'
    
  END FUNCTION collision_time

!!$  SUBROUTINE fill_array(min,max,arr,n)
!!$
!!$    !Fills array 'arr' in equally spaced intervals
!!$    IMPLICIT NONE
!!$    INTEGER :: i
!!$    DOUBLE PRECISION, INTENT(IN) :: min, max
!!$    DOUBLE PRECISION, ALLOCATABLE, INTENT(INOUT) :: arr(:)
!!$    INTEGER, INTENT(IN) :: n
!!$
!!$    !Allocate the array, and deallocate it if it is full
!!$    IF(ALLOCATED(arr)) DEALLOCATE(arr)
!!$    ALLOCATE(arr(n))
!!$    arr=0.
!!$
!!$    IF(n==1) THEN
!!$       arr(1)=min
!!$    ELSE IF(n>1) THEN
!!$       DO i=1,n
!!$          arr(i)=min+(max-min)*REAL(i-1)/REAL(n-1)
!!$       END DO
!!$    END IF
!!$
!!$  END SUBROUTINE fill_array
!!$
!!$  SUBROUTINE RNG_set(seed)
!!$
!!$    !Seeds the RNG
!!$    IMPLICIT NONE
!!$    INTEGER :: int, timearray(3)
!!$    DOUBLE PRECISION :: rand
!!$    INTEGER, INTENT(IN) :: seed
!!$
!!$    WRITE(*,*) 'RNG_SET: Initialising random number generator'
!!$    WRITE(*,*) 'RNG_SET: seed:', seed
!!$
!!$    IF(seed==0) THEN
!!$       !This fills the time array using the system clock!
!!$       !If called within the same second the numbers will be identical!
!!$       CALL itime(timeArray)
!!$       !This then initialises the generator!
!!$       int=FLOOR(rand(timeArray(1)+timeArray(2)+timeArray(3)))
!!$    ELSE
!!$       !In this case you can keep track of the seed
!!$       int=FLOOR(rand(seed))
!!$    END IF
!!$    WRITE(*,*) 'RNG_SET: done'
!!$    WRITE(*,*)
!!$
!!$  END SUBROUTINE RNG_set
!!$
!!$   FUNCTION uniform(x1,x2)
!!$
!!$    !Produces a uniform random number between x1 and x2
!!$    IMPLICIT NONE
!!$    DOUBLE PRECISION :: uniform
!!$    DOUBLE PRECISION, INTENT(IN) :: x1,x2
!!$
!!$    !Rand is some inbuilt function
!!$    uniform=x1+(x2-x1)*(rand(0))
!!$
!!$  END FUNCTION uniform

END PROGRAM exact1D
