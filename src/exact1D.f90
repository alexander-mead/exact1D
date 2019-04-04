PROGRAM exact1D

  USE random_numbers
  USE numerology
  USE array_operations
  USE sorting

  ! Solves *exactly* (i.e., to machine precision) the evoltion of 'n' interacting sheets
  IMPLICIT NONE
  REAL, ALLOCATABLE :: x0(:), v0(:), a0(:), t0(:), collision_times(:)
  REAL, ALLOCATABLE :: output_times(:)
  INTEGER, ALLOCATABLE :: id(:)
  REAL :: initial_energy, final_energy, G
  REAL :: next_collision_time, next_collision_position
  REAL :: vbar
  INTEGER :: i, j
  CHARACTER(len=256) :: infile, outfile

  INTEGER :: number_of_collisions, number_of_sheets_colliding
  INTEGER :: ii, jj, j1, j2
  REAL :: xc
  REAL, ALLOCATABLE :: collision_velocity(:)
  INTEGER, ALLOCATABLE :: collision_times_index(:), collision_velocity_index(:), collision_ids(:)

  REAL, ALLOCATABLE :: phase_shot(:,:)
  REAL, ALLOCATABLE :: phase_space(:,:)

  ! Parameters
  LOGICAL, PARAMETER :: verbose=.FALSE.        ! Set the speaking level
  LOGICAL, PARAMETER :: extra_verbose=.FALSE.  ! Set the extra speaking level
  INTEGER, PARAMETER :: solution_method=1      ! Method to solve equations
  INTEGER, PARAMETER :: fp1=7                  ! File-pointer 1
  INTEGER, PARAMETER :: fp2=8                  ! File-pointer 2
  REAL, PARAMETER :: tolerance=1e4*epsilon(1.) ! Tolerance for 'exactly the same time/positon'
  INTEGER, PARAMETER :: index_method=2         ! Indexing method
  LOGICAL, PARAMETER :: use_restart=.FALSE.    ! Use previous simulation end

  INTEGER :: n=50                             ! Number of sheets
  REAL, PARAMETER :: L=1.                     ! Size of inital group
  INTEGER, PARAMETER :: iseed=0               ! Random-number seed
  LOGICAL, PARAMETER :: ran_offsets=.TRUE.    ! Do random offsetting
  REAL, PARAMETER :: alpha=0.01               ! Amount to offset in terms of spacing (1 is maximum)
  LOGICAL, PARAMETER :: ran_velocities=.TRUE. ! Do random velocties  
  REAL, PARAMETER :: beta=0.01                ! Random velocity maximum

  REAL, PARAMETER :: initial_time=0. ! Initial time (is this necessary?)
  REAL, PARAMETER :: final_time=10.  ! Final time
  INTEGER, PARAMETER :: ntimes=500   ! Number of time outputs

  INTEGER, PARAMETER :: mesh_phase=256 ! Phase-space mesh size

  ! Potential improvements:
  ! Give each sheet an 'acceleration contribution' and then ditch G

  ! Potential speed ups:
  ! Do not need array for the acceleration, could just use position info
  ! Use sorting algorithm for collision times
  ! Method 1 is way faster, but does not work if there are simultaneous collisions
  ! There is probably a nice hybrid method to use

  ! Potential problems:
  ! Different pairs of sheets collide at exactly the same time, but at different positions
  ! Quantised positions and velocities due to numerical precision

  WRITE(*,*)
  WRITE(*,*) 'EXACT1D: 1D particle calculation'
  WRITE(*,*) 'EXACT1D: Solving system using method:', solution_method
  IF(solution_method==2) THEN
     WRITE(*,*) 'EXACT1D: Position and time error tolerance parameter:', tolerance
     WRITE(*,*) 'EXACT1D: Indexing method', index_method
  END IF
  WRITE(*,*)

  IF(use_restart) THEN

     infile='data/end.dat'
     WRITE(*,*) 'EXACT1D: Starting simulation from restart file: ', TRIM(infile)
     CALL read_restart(x0,v0,a0,id,n,infile)
     G=1./REAL(n-1)
     WRITE(*,*)

  ELSE

     IF(ran_offsets .OR. ran_velocities) CALL RNG_set(iseed)

     ! Set the number of sheets and allocate arrays
     WRITE(*,*) 'EXACT1D: Size of initial group', L
     WRITE(*,*) 'EXACT1D: Number of sheets:', n
     ALLOCATE(x0(n),v0(n),a0(n),id(n))    

     IF(ran_offsets) THEN
        WRITE(*,*) 'EXACT1D: Sheets are randomly offset from equally spaced'
        WRITE(*,*) 'EXACT1D: Offset, in terms of the inter-sheet spacing:', alpha
     ELSE
        WRITE(*,*) 'EXACT1D: Sheets start equally spaced' 
     END IF
     IF(ran_velocities) THEN
        WRITE(*,*) 'EXACT1D: Sheets start with random velocities'
        WRITE(*,*) 'EXACT1D: Sheet random velocity:', beta
     ELSE
        WRITE(*,*) 'EXACT1D: Sheets start stationary'
     END IF

     ! Fill the initial position array
     CALL fill_array(-L/2.,L/2.,x0,n)

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
           v0(i)=random_uniform(-beta,beta)
        END DO
        vbar=SUM(v0)/REAL(n) ! Calculate the mean velocity
        v0=v0-vbar ! Ensure that the mean velocity is zero
     ELSE
        ! Set the initial velocities to zero
        v0=0.
     ENDIF

     ! Set the initial accelerations, which only depend on the particle position
     ! Acceleration array does not change EVER if we keep the array order to be the positional order
     ! Total Acceleration should be G/(n-1)
     G=1./REAL(n-1)
     CALL fill_array(G*(n-1),-G*(n-1),a0,n)     
     !CALL fill_array((n-1)/REAL(n-1),-(n-1)/REAL(n-1),a0,n)     

  END IF

  ! Set the 'initial times' for all sheets to be t1
  ALLOCATE(t0(n))
  t0=initial_time

  ! Write initial conditions
  outfile='data/start.dat'
  CALL write_restart(initial_time,x0,v0,a0,t0,n,outfile)

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
  CALL fill_array(initial_time,final_time,output_times,ntimes)
  WRITE(*,*) 'EXACT1D: Initial time:', initial_time
  WRITE(*,*) 'EXACT1D: Final time:', final_time
  WRITE(*,*) 'EXACT1D: Number of output times:', ntimes 
  WRITE(*,*)

  ! Calculate the initial collision times between all neighbouring pairs (n-1)
  ALLOCATE(collision_times(n-1))
  IF(extra_verbose) THEN
     WRITE(*,*) '======================================='
     WRITE(*,*) '        i       i+1     collision times'
     WRITE(*,*) '======================================='
  END IF
  DO j=1,n-1  
     collision_times(j)=collision_time(initial_time,x0(j),v0(j),a0(j),t0(j),x0(j+1),v0(j+1),a0(j+1),t0(j+1))
     IF(extra_verbose) WRITE(*,fmt='(2I10,ES20.10)') j, j+1, collision_times(j)
  END DO
  IF(extra_verbose) THEN
     WRITE(*,*) '======================================='
     WRITE(*,*)
  END IF

  ! Open files for data
  OPEN(fp1,file='data/output.dat')
  OPEN(fp2,file='data/energy.dat')

  ! Make a nice display grid
  IF(verbose) THEN
     WRITE(*,*) '============================================'
     WRITE(*,*) '          time          E=K+T            T/V'
     WRITE(*,*) '============================================'
  END IF

  ! Calculate the initial energy of the system
  initial_energy=kinetic_energy(initial_time,v0,a0,t0,n)+potential_energy(initial_time,x0,v0,a0,t0,n)

  ! Initial write of data to file/screen
  CALL write_data(output_times(1),x0,v0,a0,t0,id,n,fp1)
  CALL write_energy(output_times(1),x0,v0,a0,t0,n,fp2)
  ALLOCATE(phase_shot(mesh_phase,mesh_phase),phase_space(mesh_phase,mesh_phase))
  phase_space=0.
  CALL make_phase_space(2.*L,4.*L,output_times(1),x0,v0,a0,t0,n,phase_shot,mesh_phase)
  phase_space=phase_space+phase_shot 

  ! Loop over output times
  DO i=1,ntimes-1

     IF(solution_method==1) THEN

        IF(i==1) THEN

           ! Calculate the next collision time
           ! This might be slightly wasteful. It could probably be ordered somehow
           next_collision_time=MINVAL(collision_times)

        END IF

        ! If the next collision time is before the output time then update the system
        DO WHILE (next_collision_time<output_times(i+1))

           ! Update the array of collision times and the x0 etc.
           CALL collision_update(next_collision_time,collision_times,x0,v0,a0,t0,id,n)

           ! Get the next collision time
           next_collision_time=MINVAL(collision_times)

        END DO

     ELSE IF(solution_method==2) THEN

        IF(i==1) THEN

            ! Index the collision times and store the first collision time and position
           ALLOCATE(collision_times_index(n-1))     
           CALL index(collision_times,collision_times_index,n-1,index_method)
           next_collision_time=collision_times(collision_times_index(1))
           
        END IF

        DO WHILE(next_collision_time<output_times(i+1))

           ! Count the total number of sheet pair collisions (usually this will be only one)
           number_of_collisions=1
           DO j=2,n-1
              IF(ABS(collision_times(collision_times_index(j))-next_collision_time)<=tolerance) THEN
                 number_of_collisions=number_of_collisions+1
              END IF
           END DO

           IF(number_of_collisions==1) THEN

              CALL collision_update(next_collision_time,collision_times,x0,v0,a0,t0,id,n)

           ELSE

              ! There is always one more sheet colliding than the total number of collisons (assuming all the same position)
              number_of_sheets_colliding=number_of_collisions+1

              ! Calculate the position of the collision
              ii=collision_times_index(1)
              next_collision_position=position(next_collision_time,x0(ii),v0(ii),a0(ii),t0(ii))

              ! Make arrays of the coordinates of the colliding sheets
              ALLOCATE(collision_velocity(number_of_sheets_colliding),collision_ids(number_of_sheets_colliding))
              DO j=1,number_of_sheets_colliding

                 ! Indices of the sheets that are colliding
                 jj=MINVAL(collision_times_index)+j-1

                 ! Store the IDs of the sheets that are colliding
                 collision_ids(j)=id(jj)

                 ! Positions of sheets as they collide (should be identical)
                 ! TODO: Account for collisions that occur at the same time but at different places
                 xc=position(next_collision_time,x0(jj),v0(jj),a0(jj),t0(jj)) 

                 ! Store the collision velocities for later ordering
                 collision_velocity(j)=velocity(next_collision_time,v0(jj),a0(jj),t0(jj))

                 ! Write info to screen for debugging
                 !WRITE(*,*) j, xc, collision_velocity(j), collision_ids(j)

                 ! Check that the collisions really are happening at the same location
                 IF(ABS(xc-next_collision_position)>=tolerance) THEN
                    STOP 'EXACT1D: Error, different pairs colliding at different locations'
                 END IF

              END DO

              ! Index velocties and update new coordinates in velocity order       
              ALLOCATE(collision_velocity_index(number_of_sheets_colliding))
              CALL index(collision_velocity,collision_velocity_index,number_of_sheets_colliding,index_method)

              ! Write to screen useful information
              DO j=1,number_of_sheets_colliding
                 jj=collision_velocity_index(j)
                 !WRITE(*,*) 'Velocity:', j, id(j), collision_velocity(j), jj, collision_velocity(jj)
              END DO

              ! Set the new x0 etc. values for the sheets that have collided
              DO j=1,number_of_sheets_colliding
                 jj=MINVAL(collision_times_index)+j-1 ! Indices of sheets that are colliding
                 x0(jj)=next_collision_position ! New position is the collision position
                 v0(jj)=collision_velocity(collision_velocity_index(j)) ! New velocities go from lowest to highest
                 t0(jj)=next_collision_time ! New t is the collision time
                 id(jj)=collision_ids(collision_velocity_index(j)) ! New indices go from left to right
              END DO
              DEALLOCATE(collision_velocity,collision_ids,collision_velocity_index)

              ! Finally, update with the new collision times
              DO j=1,number_of_collisions
                 j1=MINVAL(collision_times_index)+j-1 ! Indices of sheets that are colliding
                 j2=j1+1 ! Next sheet                 
                 collision_times(j1)=collision_time(next_collision_time,x0(j1),v0(j1),a0(j1),t0(j1),x0(j2),v0(j2),a0(j2),t0(j2))
              END DO

           END IF

           ! Re-index the collision time array and get the new next-collision time
           CALL index(collision_times,collision_times_index,n-1,index_method)
           next_collision_time=collision_times(collision_times_index(1))

        END DO

     END IF

     ! Otherwise write data for file/screen when the output time is before the next collision
     CALL write_data(output_times(i+1),x0,v0,a0,t0,id,n,fp1)
     CALL write_energy(output_times(i+1),x0,v0,a0,t0,n,fp2)
     CALL make_phase_space(2.*L,4.*L,output_times(i+1),x0,v0,a0,t0,n,phase_shot,mesh_phase)
     phase_space=phase_space+phase_shot

  END DO

  ! Divide the phase space by the number of timeshots that went into it
  phase_space=phase_space/REAL(ntimes)
  outfile='data/phase_density.dat'
  CALL write_phase_space(2.*L,4.*L,phase_space,mesh_phase,outfile)

  ! Close files for data
  CLOSE(fp1)
  CLOSE(fp2)

  ! Write final conditions
  outfile='data/end.dat'
  CALL write_restart(final_time,x0,v0,a0,t0,n,outfile)

  ! Calculate the final energy
  final_energy=kinetic_energy(final_time,v0,a0,t0,n)+potential_energy(final_time,x0,v0,a0,t0,n)

  ! Write a nice grid to screen
  IF(verbose) THEN
     WRITE(*,*) '============================================'
     WRITE(*,*)
  END IF

  ! Check that energy conservation is okay
  WRITE(*,*) 'EXACT1D: Initial energy:', initial_energy
  WRITE(*,*) 'EXACT1D: Final energy:', final_energy
  WRITE(*,*) 'EXACT1D: Ratio of final to initial energy:', final_energy/initial_energy
  WRITE(*,*)

CONTAINS

  SUBROUTINE make_phase_space(L,U,t,x0,v0,a0,t0,n,p,m)

    IMPLICIT NONE   
    REAL, INTENT(IN) :: L, U, t, x0(n), v0(n), a0(n), t0(n)
    INTEGER, INTENT(IN) :: n, m
    REAL, ALLOCATABLE, INTENT(OUT) :: p(:,:)
    REAL :: xv(2,n), w(n), LU(2)
    INTEGER :: m2(2), i

    DO i=1,n
       xv(1,i)=position(t,x0(i),v0(i),a0(i),t0(i))
       xv(2,i)=velocity(t,v0(i),a0(i),t0(i))
    END DO

    w=1.
    m2=m
    LU(1)=L
    LU(2)=U
    DO  i=1,2
       xv(i,:)=xv(i,:)+LU(i)/2.
    END DO
    ALLOCATE(p(m2(1),m2(2)))
    CALL CIC2D(xv,n,LU,w,p,m2)
    p=p/(REAL(n)/REAL(mesh_phase)**2)

  END SUBROUTINE make_phase_space

  SUBROUTINE write_phase_space(L,U,p,m,outfile)

    IMPLICIT NONE
    REAL, INTENT(IN) :: L, U, p(m,m)
    INTEGER, INTENT(IN) :: m
    CHARACTER(len=*), INTENT(IN) :: outfile
    REAL :: x, v
    INTEGER :: i, j

    OPEN(10,file=outfile,status='replace')
    DO j=1,m
       DO i=1,m
          x=cell_position(i,L,m)-L/2.
          v=cell_position(j,U,m)-U/2.
          WRITE(10,*) x, v, p(i,j)
       END DO
    END DO
    CLOSE(10)

!!$    OPEN(11,file='')
!!$    DO
!!$
!!$    END DO
!!$    CLOSE(11)
    
  END SUBROUTINE write_phase_space

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
          tc(i)=collision_time(t,x0(i),v0(i),a0(i),t0(i),x0(i+1),v0(i+1),a0(i+1),t0(i+1))
       END IF
    END DO

  END SUBROUTINE collision_update

  FUNCTION collision_time(t,x1,v1,a1,t1,x2,v2,a2,t2)

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
       ! TODO: Speed up if tc1<tc2
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
!!$       collision_time=MIN(tc1,tc2)
    END IF

    IF(collision_time==HUGE(tc)) STOP 'COLLISION_TIME: Error, this should never happen'

  END FUNCTION collision_time

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

    ! Do a double-loop over all sheet pairs to calculate the total potential energy
    DO j=1,n
       DO i=j+1,n
          P=P+G*ABS(x(j)-x(i))!/REAL(n-1)
          !P=P+ABS(x(j)-x(i))
       END DO
    END DO

    ! Set the result
    potential_energy=P

  END FUNCTION potential_energy

  SUBROUTINE read_restart(x0,v0,a0,id,n,infile)

    USE file_info

    ! Read in a restart file
    IMPLICIT NONE
    REAL, ALLOCATABLE, INTENT(OUT) :: x0(:), v0(:), a0(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: id(:)
    INTEGER, INTENT(OUT) :: n
    CHARACTER(len=*), INTENT(IN) :: infile
    INTEGER :: i
    LOGICAL, PARAMETER :: verbose=.TRUE.

    n=file_length(infile,verbose)
    ALLOCATE(x0(n),v0(n),a0(n),id(n))

    OPEN(10,file=infile)
    DO i=1,n
       READ(10,*) x0(i), v0(i), a0(i)
       id(i)=i
    END DO
    CLOSE(10)

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

    OPEN(10,file=outfile)
    DO i=1,n
       x=position(t,x0(i),v0(i),a0(i),t0(i))
       v=velocity(t,v0(i),a0(i),t0(i))
       a=a0(i)
       WRITE(10,*) x, v, a
    END DO
    CLOSE(10)

  END SUBROUTINE write_restart

  SUBROUTINE write_data(t,x0,v0,a0,t0,id,n,fp)

    ! Figure out where the sheets are in the arrays and write to file
    IMPLICIT NONE
    REAL, INTENT(IN) :: t, x0(n), v0(n), a0(n), t0(n)
    INTEGER, INTENT(IN) :: id(n), n, fp
    INTEGER :: i, j, idx(n)
    REAL :: xn(n), vn(n)

    LOGICAL, PARAMETER :: bubble=.TRUE.

    IF(bubble) THEN

       CALL bubble_index(id,idx,n)
       DO i=1,n
          j=idx(i)
          xn(i)=position(t,x0(j),v0(j),a0(j),t0(j))
          vn(i)=velocity(t,v0(j),a0(j),t0(j))
       END DO
       WRITE(fp,*) t, (xn(i), vn(i), i=1,n)

    ELSE

       ! This is definitely a stupid way to write out data
       DO i=1,n
          DO j=1,n
             IF(id(j)==i) THEN
                xn(i)=position(t,x0(j),v0(j),a0(j),t0(j))
                vn(i)=velocity(t,v0(j),a0(j),t0(j))
                EXIT
             END IF
          END DO
       END DO

       ! Write to file
       WRITE(fp,*) t, (xn(i), vn(i), i=1,n)

    END IF

  END SUBROUTINE write_data

  SUBROUTINE write_energy(t,x0,v0,a,t0,n,fp)

    IMPLICIT NONE
    REAL, INTENT(IN) :: t, x0(n), v0(n), a(n), t0(n)
    INTEGER, INTENT(IN) :: n, fp
    REAL :: K, V, E, Vir

    ! Calculate the energies
    K=kinetic_energy(t,v0,a,t0,n)
    V=potential_energy(t,x0,v0,a,t0,n)
    E=K+V
    Vir=K/V

    ! Write the energies
    IF(verbose) WRITE(*,fmt='(3F15.7)') t, E, Vir  
    WRITE(fp,*) t, K, V, E, Vir

  END SUBROUTINE write_energy

  INTEGER FUNCTION NGP_cell(x,L,m)

    ! Find the integer coordinates of the cell that coordinate x is in
    IMPLICIT NONE
    REAL, INTENT(IN) :: x ! Particle position
    REAL, INTENT(IN) :: L ! Box size
    INTEGER, INTENT(IN) :: m ! Number of mesh cells in grid

    IF(x==0.) THEN
       ! Catch this edge case
       NGP_cell=1
    ELSE
       NGP_cell=CEILING(x*REAL(m)/L)
    END IF

    IF(NGP_cell<1 .OR. NGP_cell>m) THEN
       WRITE(*,*) 'NGP_CELL: Particle position [Mpc/h]:', x
       WRITE(*,*) 'NGP_CELL: Box size [Mpc/h]:', L
       WRITE(*,*) 'NGP_CELL: Mesh size:', m 
       WRITE(*,*) 'NGP_CELL: Assigned cell:', NGP_cell
       STOP 'NGP_CELL: Error, the assigned cell position is outside the mesh'
    END IF

  END FUNCTION NGP_cell

  REAL FUNCTION cell_position(i,L,m)

    ! Gets the coordinates of cell centre i in box of length L with m cells
    IMPLICIT NONE
    REAL, INTENT(IN) :: L
    INTEGER, INTENT(IN) :: i, m

    cell_position=L*(i-0.5)/REAL(m)

  END FUNCTION cell_position

  SUBROUTINE CIC2D(x,n,L,w,d,m)

    ! This could probably be usefully combined with CIC somehow
    IMPLICIT NONE
    REAL, INTENT(IN) :: x(2,n), w(n), L(2)
    INTEGER, INTENT(IN) :: n, m(2)
    REAL, INTENT(OUT) :: d(m(1),m(2))
    INTEGER :: i, j, ix1(2), ix2(2)
    REAL :: dx(2)

    LOGICAL, PARAMETER :: verbose=.FALSE.

    IF(verbose) THEN
       WRITE(*,*) 'CIC2D: Binning particles and creating density field'
       WRITE(*,*) 'CIC2D: Cells:', m
    END IF

    ! Set array to zero explicitly because it is used for a sum
    d=0.

    ! Loop over particles
    DO i=1,n

       ! Loop over coordinates in x
       DO j=1,2

          ! Get the cell interger coordinates
          ix1(j)=NGP_cell(x(j,i),L(j),m(j))

          ! dx, dy in cell units, away from cell centre
          dx(j)=(x(j,i)/L(j))*REAL(m(j))-(REAL(ix1(j))-0.5)

          ! Find CIC weights and enforce periodicity
          IF(dx(j)>0.) THEN
             ix2(j)=ix1(j)+1
             IF(ix2(j)>m(j)) ix2(j)=1
          ELSE
             ix2(j)=ix1(j)-1
             dx(j)=-dx(j)  
             IF(ix2(j)<1) ix2(j)=m(j)   
          END IF

       END DO

       ! Carry out CIC binning
       d(ix1(1),ix1(2))=d(ix1(1),ix1(2))+(1.-dx(1))*(1.-dx(2))*w(i)
       d(ix1(1),ix2(2))=d(ix1(1),ix2(2))+(1.-dx(1))*dx(2)*w(i)
       d(ix2(1),ix1(2))=d(ix2(1),ix1(2))+dx(1)*(1.-dx(2))*w(i)
       d(ix2(1),ix2(2))=d(ix2(1),ix2(2))+dx(1)*dx(2)*w(i)

    END DO

    IF(verbose) THEN
       WRITE(*,*) 'CIC2D: Binning complete'
       WRITE(*,*)
    END IF

  END SUBROUTINE CIC2D

END PROGRAM exact1D
