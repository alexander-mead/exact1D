PROGRAM PM1D

  !Necessary for FFTW!
  USE constants
  USE file_info
  USE array_operations
  USE string_operations
  USE fft
  USE cosmology_functions
  USE simulations1D

  !Always necessary
  IMPLICIT NONE

  !Necessary include statements for FFTW
  !INCLUDE '/usr/include/fftw3.f' !Linux
  !INCLUDE '/usr/local/include/fftw3.f' !OS X

  !Allocate variables in main program
  REAL, ALLOCATABLE :: x(:), v(:)
  REAL, ALLOCATABLE :: a(:)
  REAL, ALLOCATABLE :: d(:), F(:)
  REAL, ALLOCATABLE :: d_visual(:), da_visual(:,:)
  REAL, ALLOCATABLE :: phase(:,:)
  INTEGER :: i,j, n
  CHARACTER(len=256) :: infile, outfile, fext
  CHARACTER(len=256) :: density_base, particle_base, power_base
  CHARACTER(len=256) :: phase_base, histogram_base
  REAL :: da
  INTEGER :: icosmo
  TYPE(cosmology) :: cosm
  
  !Simulation parameters  
  REAL, PARAMETER :: a1=1 !Starting scale factor
  REAL, PARAMETER :: a2=1d4 !Final scale factor
  REAL, PARAMETER :: Om_m=1. !Omega_m
  REAL, PARAMETER :: Om_v=0. !Omega_v
  REAL, PARAMETER :: L=1. !Box size [Mpc/h]
  REAL, PARAMETER :: dx=0. !Shift [Mpc/h]
  INTEGER, PARAMETER :: m=2**20 !Force mesh
  INTEGER, PARAMETER :: s=128 !Time steps
  LOGICAL, PARAMETER :: log_time=.TRUE. !Log or linear time evolution
  LOGICAL, PARAMETER :: smooth_force=.FALSE. !Smooth the graviational force
  INTEGER, PARAMETER :: ibin=2 !Interpolation strategy (1 - NGP, 2 - CIC)
  LOGICAL, PARAMETER :: leapfrog=.TRUE. !Use leapfrog integration or not
  LOGICAL, PARAMETER :: KDK=.FALSE. !KDK or DKD
  
  !Output parameters
  LOGICAL, PARAMETER :: verbose=.TRUE. !Verbose mode
  INTEGER, PARAMETER :: m_visual=4096 !Visual density mesh (4096 is nice)
  INTEGER, PARAMETER :: s_visual=128 !Visual time mesh (3072 is nice)
  REAL, PARAMETER :: L1=0. !The lower limit for the output pictures
  REAL, PARAMETER :: L2=L  !The upper limit for the output pictures
  LOGICAL, PARAMETER :: visual_periodic=.TRUE. !Is the visual output a periodic volume?
  INTEGER, PARAMETER :: file_zeroes=4 !Number of zeroes in file names
  LOGICAL, PARAMETER :: output_particles=.FALSE.
  LOGICAL, PARAMETER :: output_density=.TRUE.
  LOGICAL, PARAMETER :: output_power=.FALSE.
  LOGICAL, PARAMETER :: output_histogram=.FALSE.
  LOGICAL, PARAMETER :: output_phase=.TRUE.

  !=====================================================!
  !============= Start the main program ================!
  !=====================================================!

  CALL get_command_argument(1,infile)
  IF(infile=='') STOP 'PM1D: Error, specify an input file'

  IF(verbose) THEN
     WRITE(*,*)
     WRITE(*,*) 'PM1D: 1-Dimensional particle-mesh simuilations'
     WRITE(*,*) 'PM1D: Verbose mode'
     WRITE(*,*) 'PM1D: Leapfrog integration:', leapfrog
     IF(leapfrog) THEN
        IF(KDK)  WRITE(*,*) 'PM1D: Kick-drift-kick scheme'
        IF(KDK .EQV. .FALSE.) WRITE(*,*) 'PM1D: Drift-kick-drift scheme'
     END IF
     WRITE(*,*) 'PM1D: Particle interpolation scheme:', ibin
     WRITE(*,*) 'PM1D: Number of zeroes in file names:', file_zeroes
     WRITE(*,*)
  END IF

  !Assign the cosmological model
  icosmo=1
  CALL assign_cosmology(icosmo,cosm)

  !Read in particle data
  n=file_length(infile)
  ALLOCATE(x(n),v(n))
  CALL read_particles(x,v,n,infile,verbose)

  !Write out initial conditions as particle data
  outfile='start.dat'
  CALL write_particles(x,v,n,outfile,verbose)

  !Shift particle positions (should only matter for visualisation...)
  x=x+dx
  CALL enforce_periodicity(x,n,L)

  !Initial and final 'a' values and time stepping
  IF(log_time) THEN
     CALL fill_array(log(a1),log(a2),a,s+1)
     a=exp(a)
  ELSE IF(log_time .EQV. .FALSE.) THEN
      CALL fill_array(a1,a2,a,s+1)
   END IF

  !Allocate density and force meshes for calculation
  ALLOCATE(d(m),F(m))

  !Write simulation information
  IF(verbose) THEN
     WRITE(*,*) 'PM1D: Starting a:', a(1)
     WRITE(*,*) 'PM1D: Final a:', a(s+1)
     WRITE(*,*) 'PM1D: Number of time steps:', s
     WRITE(*,*) 'PM1D: Starting da:', a(2)-a(1)
     WRITE(*,*) 'PM1D: Starting position shift [Mpc/h]', dx
     WRITE(*,*) 'PM1D: Force mesh:', m
     WRITE(*,*)
  END IF

  !File extensions for writing
  particle_base=TRIM('particles/data_')
  power_base=TRIM('power/data_')
  histogram_base=TRIM('histogram/data_')
  phase_base=TRIM('phase/data_')
  fext=TRIM('.dat')

  !Write directory locations to screen
  IF(verbose) THEN
     IF(output_particles) WRITE(*,*) 'PM1D: Particle output:', TRIM(particle_base)
     IF(output_power) WRITE(*,*) 'PM1D: Power output:', TRIM(power_base)
     IF(output_histogram) WRITE(*,*) 'PM1D: Histogram output:', TRIM(histogram_base)
     IF(output_phase) WRITE(*,*) 'PM1D: Phase output:', TRIM(phase_base)
     WRITE(*,*)
  END IF

  !Allocate meshes for visualisation
  IF(verbose) THEN
     WRITE(*,*) 'PM1D: Visual mesh:', m_visual
     WRITE(*,*) 'PM1D: Visual time steps:', s_visual
     IF(s_visual>s) STOP 'PM1D: Error, visual mesh is larger than time mesh'
     IF(MOD(s,s_visual) .NE. 0) STOP 'PM1D: Error, visual mesh is not 2**n less than time mesh'
     WRITE(*,*)
  END IF
  ALLOCATE(d_visual(m_visual),da_visual(m_visual,s_visual+1))
  d_visual=0.
  da_visual=0.

  !Write initial particle data
  IF(output_particles) THEN
     outfile=number_file_zeroes(particle_base,0,file_zeroes,fext)
     CALL write_particles(x,v,n,outfile,verbose)
  END IF

  !Calculate the initial density field for visualisation
  IF((output_histogram) .OR. (output_density)) THEN
     CALL calculate_density(x,n,L1,L2,L,d_visual,m_visual,visual_periodic)
     IF(output_density) THEN
        da_visual(:,1)=d_visual
     END IF
     IF(output_histogram) THEN
        outfile=number_file_zeroes(histogram_base,0,file_zeroes,fext)
        CALL histogram(1.+d_visual,m_visual,outfile)
     END IF
  END IF

  !Do the first power spectrum calculation
  IF(output_power) THEN
     outfile=number_file_zeroes(power_base,0,file_zeroes,fext)
     CALL calculate_density(x,n,0.,L,L,d,m,.TRUE.)
     CALL Pk1D(d,m,L,outfile)
  END IF

  !Do the first phasespace calculation
  IF(output_phase) THEN
     outfile=number_file_zeroes(phase_base,0,file_zeroes,fext)
     CALL calculate_phasespace(x,v,n,outfile)
  END IF

  !Do the integration from a(i) to a(i+1) for timesteps 1 to s
  DO i=1,s

     !Set the current time step
     da=a(i+1)-a(i)

     !Write to screen
     IF(verbose) THEN
        WRITE(*,fmt='(A17,I6,3F15.3)') 'PM1D: Timestep:', i, a(i), a(i+1), da
     END IF

     !Do the actual integration
     IF(leapfrog .EQV. .FALSE.) THEN
        !Use the density to calculate the forces in Fourier space
        !Then move particles around
        CALL calculate_density(x,n,0.,L,L,d,m,.TRUE.)      
        CALL calculate_forces(d,F,m,L,a(i),cosm)
        CALL move_particles(x,v,n,L,F,m,a(i),da,cosm)
        !TEST
        !CALL move_kick(x,v,n,L,F,m,a(i),da/2.)
        !CALL move_drift(x,v,n,L,a(i),da/2.)
        !CALL move_kick(x,v,n,L,F,m,a(i)+da/2.,da/2.)
        !CALL move_drift(x,v,n,L,a(i)+da/2.,da/2.)
        !TEST
     ELSE IF(leapfrog) THEN
        IF(KDK) THEN
           !Leap frog KDK
           CALL calculate_density(x,n,0.,L,L,d,m,.TRUE.)
           CALL calculate_forces(d,F,m,L,a(i),cosm)
           CALL move_kick(x,v,n,L,F,m,a(i),da/2.,cosm)
           CALL move_drift(x,v,n,L,a(i)+da/2.,da,cosm)
           CALL calculate_density(x,n,0.,L,L,d,m,.TRUE.)
           CALL calculate_forces(d,F,m,L,a(i)+da,cosm)
           CALL move_kick(x,v,n,L,F,m,a(i),da/2.,cosm)
        ELSE
           !Leap frog - DKD
           CALL move_drift(x,v,n,L,a(i),da/2.,cosm)
           CALL calculate_density(x,n,0.,L,L,d,m,.TRUE.)
           CALL calculate_forces(d,F,m,L,a(i)+da/2.,cosm)
           CALL move_kick(x,v,n,L,F,m,a(i)+da/2.,da,cosm)
           CALL move_drift(x,v,n,L,a(i)+da/2.,da/2.,cosm)
        END IF
     END IF

     !Output the data on an output step
     IF(MOD(i,s/s_visual)==0) THEN

        j=i/(s/s_visual)

        !Do a phasespace calculation
        IF(output_phase) THEN
           outfile=number_file_zeroes(phase_base,j,file_zeroes,fext)
           CALL calculate_phasespace(x,v,n,outfile)
        END IF

        !Do a power spectrum calculation
        IF(output_power) THEN
           outfile=number_file_zeroes(power_base,j,file_zeroes,fext)
           CALL Pk1D(d,m,L,outfile)
        END IF

        !Output the particle data
        IF(output_particles) THEN
           !Write particle data
           outfile=number_file_zeroes(particle_base,j,file_zeroes,fext)
           CALL write_particles(x,v,n,outfile,verbose)
        END IF

        !Calculate the initial density field for visualisation
        IF((output_histogram) .OR. (output_density)) THEN
           CALL calculate_density(x,n,L1,L2,L,d_visual,m_visual,visual_periodic)
           IF(output_density) THEN
              da_visual(:,j+1)=d_visual
           END IF
           IF(output_histogram) THEN
              outfile=number_file_zeroes(histogram_base,j,file_zeroes,fext)
              CALL histogram(1.+d_visual,m_visual,outfile)
           END IF
        END IF

     END IF

  END DO
  IF(verbose) WRITE(*,*)

  !Output the large density vs. a plot
  IF(output_density) THEN
     outfile='density/bigdensity.dat'
     CALL write_bigdensity(da_visual,m_visual,s_visual+1,L,outfile)
  END IF

  !Write the final particle positions out
  outfile='end.dat'
  CALL write_particles(x,v,n,outfile,verbose)

CONTAINS

  SUBROUTINE move_particles(x,v,n,L,F,m,a,da,cosm)

    IMPLICIT NONE
    REAL, INTENT(INOUT) :: x(n), v(n)
    REAL, INTENT(IN) :: F(m), L
    REAL, INTENT(IN) :: a, da
    INTEGER, INTENT(IN) :: n, m
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: dt, Fv, dx, cs, xc
    INTEGER :: i, j

    !Calculate the time interval
    dt=da/(sqrt(Hubble2(a,cosm))*a)

    !Update the positions
    !Shittiest integration imagineable (Euler method)
    DO i=1,n
       x(i)=x(i)+v(i)*dt
    END DO

    !Enforce periodicity
    CALL enforce_periodicity(x,n,L)

    !Calculat the physical mesh size
    IF(ibin==2) cs=L/REAL(m)

    !Update velocities
    DO i=1,n
       j=CEILING(float(m)*x(i)/L)
       IF(j<1 .OR. j>m) THEN
          WRITE(*,*) 'MOVE_PARTICLES: i:', i
          WRITE(*,*) 'MOVE_PARTICLES: x [Mpc/h]', x(i)
          WRITE(*,*) 'MOVE_PARTICLES: L [Mpc/h]', L
          WRITE(*,*) 'MOVE_PARTICLES: j', j
          WRITE(*,*) 'MOVE_PARTICLES: m', m
          STOP 'MOVE_PARTICLES: Error, position in mesh assigned incorrectly'
       END IF
       IF(ibin==1) THEN
          !NGP interpolation from grid to particles
          Fv=F(j)
       ELSE IF(ibin==2) THEN
          !CIC interpolation from grid to particles
          xc=L*(REAL(j)-0.5)/(REAL(m))
          !cs=L/REAL(m)
          dx=x(i)-xc
          dx=dx/cs
          IF(dx>0.) THEN
             Fv=F(j)*(1.-dx)
             IF(j==m) THEN
                Fv=Fv+F(1)*dx
             ELSE
                Fv=Fv+F(j+1)*dx
             END IF
          ELSE
             Fv=F(j)*(1.+dx)
             IF(j==1) THEN
                Fv=Fv-F(m)*dx
             ELSE
                Fv=Fv-F(j-1)*dx
             END IF
          END IF
       END IF
       v(i)=v(i)+(Fv-2.*sqrt(Hubble2(a,cosm))*v(i))*dt
    END DO

  END SUBROUTINE move_particles

  SUBROUTINE move_drift(x,v,n,L,a,da,cosm)

    !Change the particle positions
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: x(n)
    REAL, INTENT(IN) :: v(n), L
    REAL, INTENT(IN) :: a, da
    INTEGER, INTENT(IN) :: n
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: dt
    INTEGER :: i

    !Calculate the time interval
    dt=da/(sqrt(Hubble2(a,cosm))*a)

    !Update the positions
    DO i=1,n
       x(i)=x(i)+v(i)*dt
    END DO

    !Enforce periodicity
    CALL enforce_periodicity(x,n,L)

  END SUBROUTINE move_drift

  SUBROUTINE move_kick(x,v,n,L,F,m,a,da,cosm)

    !Change the particle velocities
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: v(n)
    REAL, INTENT(IN) :: x(n), F(m), L
    REAL, INTENT(IN) :: a, da
    INTEGER, INTENT(IN) :: n, m
    TYPE(cosmology), INTENT(INOUT) :: cosm
    REAL :: dt, Fv, dx, cs, xc
    INTEGER :: i, j

    !Calculate the time interval
    dt=da/(sqrt(Hubble2(a,cosm))*a)

    !Calculat the physical mesh size
    IF(ibin==2) cs=L/REAL(m)

    !Update the velocities
    DO i=1,n
       j=CEILING(float(m)*x(i)/L)
       IF(j<1 .OR. j>m) THEN
          WRITE(*,*) 'MOVE_PARTICLES: i:', i
          WRITE(*,*) 'MOVE_PARTICLES: x [Mpc/h]', x(i)
          WRITE(*,*) 'MOVE_PARTICLES: L [Mpc/h]', L
          WRITE(*,*) 'MOVE_PARTICLES: j', j
          WRITE(*,*) 'MOVE_PARTICLES: m', m
          STOP 'MOVE_PARTICLES: Error, position in mesh assigned incorrectly'
       END IF
       IF(ibin==1) THEN
          Fv=F(j)
       ELSE IF(ibin==2) THEN
          xc=L*(REAL(j)-0.5)/(REAL(m))
          dx=x(i)-xc
          dx=dx/cs
          !WRITE(*,*) i, j, xc, cs, dx
          !IF(i==10) STOP
          IF(dx>=0.) THEN
             Fv=F(j)*(1.-dx)
             IF(j==m) THEN
                Fv=Fv+F(1)*dx
             ELSE
                Fv=Fv+F(j+1)*dx
             END IF
          ELSE
             Fv=F(j)*(1.+dx)
             IF(j==1) THEN
                Fv=Fv-F(m)*dx
             ELSE
                Fv=Fv-F(j-1)*dx
             END IF
          END IF
       END IF
       v(i)=v(i)+(Fv-2.*sqrt(Hubble2(a,cosm))*v(i))*dt
       !WRITE(*,*) i, Fv, 2.*Hubble(a)*v(i)
       !IF(i==100) STOP
    END DO

  END SUBROUTINE move_kick

  SUBROUTINE calculate_forces(d,F,m,L,a,cosm)!plan_forward,plan_backward)

    IMPLICIT NONE
    REAL, INTENT(IN) :: d(m), L
    REAL, INTENT(IN) :: a
    REAL, INTENT(OUT) :: F(m)
    INTEGER, INTENT(IN) :: m
    TYPE(cosmology), INTENT(INOUT) :: cosm
    INTEGER :: i
    DOUBLE COMPLEX :: dc(m), Fc(m), Phi
    REAL :: kx, k, sig
    !INTEGER*8 :: plan_forward, plan_backward

    dc=(0.,0.)
    dc=d

    !Calculate delta_k
    CALL FFT1(dc,dc,m,-1)
    !CALL execute_FFT1(plan_forward)

    !Smoothing gravitational force over the mesh cell size
    sig=L/REAL(m)

    !Convert delta_k to forces, via the potential
    DO i=1,m
       IF(i==1) THEN
          Fc(i)=0.
       ELSE
          CALL k_FFT1D(i,m,kx,k,L)
          Phi=-(3./2.)*(Hubble2(a,cosm)**2)*Omega_m(a,cosm)*dc(i)/k**2
          Fc(i)=-(0.,1.)*kx*Phi/(a**2)
          IF(smooth_force) Fc(i)=Fc(i)*exp(-0.5*(k*sig)**2)
       END IF
    END DO

    !Convert the newly-created force back to real space
    CALL FFT1(Fc,Fc,m,1)
    !CALL execute_FFT1(plan_backward)

    !Normalise post FFT
    F=REAL(Fc)/REAL(m)!**2

  END SUBROUTINE calculate_forces

  SUBROUTINE calculate_density(x,n,L1,L2,L,d,m,periodic)

    IMPLICIT NONE
    REAL, INTENT(IN) :: x(n), L1, L2, L
    REAL, INTENT(OUT) :: d(m)
    INTEGER, INTENT(IN) :: n, m
    LOGICAL, INTENT(IN) :: periodic
    INTEGER :: i, j, j1, j2
    REAL :: dx, cs, xc

    !Set counting variables to zero
    d=0.

    !Calculate physical mesh size for CIC binning
    IF(ibin==2) cs=(L2-L1)/REAL(m) 
    
    !Do NGP binning
    DO i=1,n
       !j=CEILING(REAL(m)*x(i)/L)
       j=CEILING(REAL(m)*(x(i)-L1)/(L2-L1))
       IF((periodic) .AND. (j<1 .OR. j>m)) THEN
          WRITE(*,*) 'CALCULATE_DENSITY: i:', i
          WRITE(*,*) 'CALCULATE_DENSITY: x [Mpc/h]', x(i)
          WRITE(*,*) 'CALCULATE_DENSITY: L [Mpc/h]', L
          WRITE(*,*) 'CALCULATE_DENSITY: j', j
          WRITE(*,*) 'CALCULATE_DENSITY: m', m
          STOP 'CALCULATE_DENSITY: Error, position in mesh assigned incorrectly'
       ELSE IF((periodic .EQV. .FALSE.) .AND. (j<1 .OR. j>m)) THEN
          CYCLE
       END IF
       IF(ibin==1) THEN
          d(j)=d(j)+1.
       ELSE IF(ibin==2) THEN
          xc=L1+(L2-L1)*(REAL(j)-0.5)/(REAL(m))
          dx=x(i)-xc
          dx=dx/cs
          j1=j
          IF(dx>=0.) THEN
             j2=j+1
             IF(periodic .EQV. .TRUE. .AND. j2==m+1) j2=1
          ELSE IF(dx<0.) THEN
             j2=j-1
             IF(periodic .EQV. .TRUE. .AND. j2==0) j2=m
          END IF
          d(j1)=d(j1)+(1.-ABS(dx))
          IF(j2>0 .AND. j2<=m) d(j2)=d(j2)+ABS(dx)
       END IF
    END DO

    !Convert to overdensity
    d=-1.+d/(((L2-L1)*REAL(n)/L)/REAL(m))

!!$          IF(dx>0.) THEN
!!$             d(j)=d(j)+(1.-dx)
!!$             IF(j==m) THEN
!!$                d(1)=d(1)+dx
!!$             ELSE
!!$                d(j+1)=d(j+1)+dx
!!$             END IF
!!$          ELSE
!!$             d(j)=d(j)+(1.+dx)
!!$             IF(j==1) THEN
!!$                d(1)=d(1)-dx
!!$             ELSE
!!$                d(j-1)=d(j-1)-dx
!!$             END IF
!!$          END IF
    
!!$    ELSE IF(ibin==2) THEN
!!$       !Do CIC binning
!!$       DO i=1,n
!!$          !j=CEILING(REAL(m)*x(i)/L)
!!$          IF(j<1 .OR. j>m) THEN
!!$             WRITE(*,*) 'CALCULATE_DENSITY: i:', i
!!$             WRITE(*,*) 'CALCULATE_DENSITY: x [Mpc/h]', x(i)
!!$             WRITE(*,*) 'CALCULATE_DENSITY: L [Mpc/h]', L
!!$             WRITE(*,*) 'CALCULATE_DENSITY: j', j
!!$             WRITE(*,*) 'CALCULATE_DENSITY: m', m
!!$             STOP 'CALCULATE_DENSITY: Error, position in mesh assigned incorrectly'
!!$          END IF
!!$          xc=L*(REAL(j)-0.5)/(REAL(m))
!!$          cs=L/REAL(m)
!!$          dx=x(i)-xc
!!$          dx=dx/cs
!!$          IF(dx>0.) THEN
!!$             d(j)=d(j)+(1.-dx)
!!$             IF(j==m) THEN
!!$                d(1)=d(1)+dx
!!$             ELSE
!!$                d(j+1)=d(j+1)+dx
!!$             END IF
!!$          ELSE
!!$             d(j)=d(j)+(1.+dx)
!!$             IF(j==1) THEN
!!$                d(1)=d(1)-dx
!!$             ELSE
!!$                d(j-1)=d(j-1)-dx
!!$             END IF
!!$          END IF
!!$       END DO
!!$    END IF

  END SUBROUTINE calculate_density

  SUBROUTINE calculate_phasespace(x,v,n,outfile)

    IMPLICIT NONE
    REAL, INTENT(IN) :: x(n), v(n)
    INTEGER, INTENT(IN) :: n
    CHARACTER(len=256), INTENT(IN) :: outfile
    REAL, ALLOCATABLE :: p(:,:)
    INTEGER :: i, ix, iv, nx, nv
    REAL :: xmin, xmax, vmin, vmax, xc, vc, vr

    !Size of x region for phase space
    xmin=0.
    xmax=L

    !Size of v region for phase space
    !Should really use v in box units
    !vr=-0.001
    vr=1e-5  
    !vr=-MAXVAL(ABS(v))
    vmin=-vr
    vmax=vr

    !Mesh size for phase space output, allocate array
    !Set phasespace counting variable to zero
    nx=512
    nv=394
    ALLOCATE(p(nx,nv))
    p=0.

    !Simple NGP binning
    DO i=1,n
       ix=CEILING(REAL(nx)*(x(i)-xmin)/(xmax-xmin))
       iv=CEILING(REAL(nv)*(v(i)-vmin)/(vmax-vmin))
       IF(ix>0 .AND. ix<=nx .AND. iv>0 .AND. iv<=nv) THEN    
          p(ix,iv)=p(ix,iv)+1.
       END IF
    END DO

    !Convert counts to a density
    p=p/(REAL(n)/REAL(nx*nv))

    !Output phase space to file
    !Could really output array, below wastes lots of space
    OPEN(7,file=outfile)
    DO iv=1,nv
       vc=vmin+(vmax-vmin)*(REAL(iv)-0.5)/(REAL(nv-1))
       DO ix=1,nx
          xc=xmin+(xmax-xmin)*(REAL(ix)-0.5)/(REAL(nx-1))
          WRITE(7,*) REAL(xc), REAL(vc), REAL(p(ix,iv))
       END DO
    END DO
    CLOSE(7)

    !Deallocate phasespace array
    DEALLOCATE(p)

  END SUBROUTINE calculate_phasespace

  SUBROUTINE write_density(d,m,L,outfile)

    !Writes the density array to a file
    IMPLICIT NONE
    REAL, INTENT(IN) :: L
    REAL, INTENT(OUT) :: d(m)
    INTEGER, INTENT(IN) :: m
    CHARACTER(len=256) :: outfile
    REAL :: x
    INTEGER :: i

    OPEN(7,file=outfile)
    DO i=1,m
       x=L*(REAL(i)-0.5)/REAL(m)
       WRITE(7,*) REAL(x), REAL(d(i))
    END DO
    CLOSE(7)

  END SUBROUTINE write_density

  SUBROUTINE write_bigdensity(d,m,s,L,outfile)

    !Writes the density-time array to a file
    IMPLICIT NONE
    REAL, INTENT(IN) :: L
    REAL, INTENT(OUT) :: d(m,s)
    INTEGER, INTENT(IN) :: m, s
    CHARACTER(len=256) :: outfile
    REAL :: x, step
    INTEGER :: i, j

    WRITE(*,*) 'WRITE_BIGDENSITY: ', TRIM(outfile)
    OPEN(7,file=outfile)
    DO j=1,s
       DO i=1,m
          x=L*(REAL(i)-0.5)/REAL(m)
          step=(REAL(j)-0.5)/REAL(s)
          WRITE(7,*) REAL(x), REAL(step), REAL(d(i,j))
       END DO
    END DO
    CLOSE(7)
    WRITE(*,*) 'WRITE_BIGDENSITY: Done'
    WRITE(*,*)

  END SUBROUTINE write_bigdensity

!!$  FUNCTION Omega_m(a)
!!$
!!$    IMPLICIT NONE
!!$    REAL :: Omega_m
!!$    REAL, INTENT(IN) :: a
!!$
!!$    Omega_m=Om_m*a**(-3)/(Hubble(a)/H0)**2
!!$
!!$  END FUNCTION Omega_m

!!$  FUNCTION Hubble(a)
!!$
!!$    IMPLICIT NONE
!!$    REAL :: Hubble
!!$    REAL, INTENT(IN) :: a
!!$
!!$    Hubble=H0*sqrt(Om_m*a**(-3)+Om_v+(1.-Om_m-Om_v)*a**(-2))
!!$
!!$  END FUNCTION Hubble

!!$  SUBROUTINE read_particles(x,v,n,infile)
!!$
!!$    IMPLICIT NONE
!!$    REAL, INTENT(OUT) :: x(n), v(n)
!!$    INTEGER, INTENT(IN) :: n
!!$    CHARACTER(len=256), INTENT(IN) :: infile
!!$    INTEGER :: i
!!$
!!$    WRITE(*,*) 'READ_PARTICLES: Inputting particle data'
!!$    WRITE(*,*) 'READ_PARTICLES: File ', TRIM(infile)
!!$
!!$    x=0.
!!$    v=0.
!!$
!!$    OPEN(7,file=infile)
!!$    DO i=1,n
!!$       READ(7,*) x(i), v(i)
!!$    END DO
!!$    CLOSE(7)
!!$
!!$    WRITE(*,*) 'READ_PARTICLES: Done'
!!$    WRITE(*,*) 
!!$
!!$  END SUBROUTINE read_particles

!!$  SUBROUTINE write_particles(x,v,n,outfile)
!!$
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: x(n), v(n)
!!$    INTEGER, INTENT(IN) :: n
!!$    CHARACTER(len=256), INTENT(IN) :: outfile
!!$    INTEGER :: i
!!$
!!$    IF(verbose) THEN
!!$       WRITE(*,*) 'WRITE_PARTICLES: Outputting particle data'
!!$       WRITE(*,*) 'WRITE_PARTICLES: File ', TRIM(outfile)
!!$    END IF
!!$
!!$    OPEN(7,file=outfile)
!!$    DO i=1,n
!!$       WRITE(7,*) x(i), v(i)
!!$    END DO
!!$    CLOSE(7)
!!$
!!$    IF(verbose) THEN
!!$       WRITE(*,*) 'WRITE_PARTICLES: Done'
!!$       WRITE(*,*)
!!$    END IF
!!$
!!$  END SUBROUTINE write_particles

!!$  FUNCTION file_length(file_name)
!!$
!!$    IMPLICIT NONE
!!$    CHARACTER(len=256) :: file_name
!!$    INTEGER :: n, file_length
!!$
!!$    WRITE(*,*) 'FILE_LENGTH: File: ', TRIM(file_name)
!!$    OPEN(7,file=file_name)
!!$
!!$    !Newer version that lacks 'data' seems okay but not well tested
!!$    n=0
!!$    DO
!!$       n=n+1
!!$       READ(7,*, end=301)
!!$    END DO
!!$
!!$    !301 is just the label to jump to when the end of the file is reached
!!$
!!$301 CLOSE(7)
!!$
!!$    file_length=n-1
!!$
!!$    WRITE(*,*) 'FILE_LENGTH: Length:', file_length
!!$    WRITE(*,*)
!!$
!!$  END FUNCTION file_length

!!$  SUBROUTINE FFT1(in,out,nx,ifb)
!!$
!!$    !Wraps the 1D FFT
!!$    IMPLICIT NONE
!!$    DOUBLE COMPLEX, INTENT(IN) :: in(nx)
!!$    DOUBLE COMPLEX, INTENT(OUT) :: out(nx)
!!$    INTEGER, INTENT(IN) :: ifb
!!$    INTEGER, INTENT(IN) :: nx
!!$    INTEGER*8 :: plan
!!$
!!$    IF(ifb .NE. 1 .AND. ifb .NE. -1) THEN
!!$       WRITE(*,*) 'FFT1: Error - need to specify forwards or backwards'
!!$    END IF
!!$
!!$    !WRITE(*,*) 'FFT1: Starting FFT - creating plan'
!!$    IF(ifb==-1) THEN
!!$       CALL dfftw_plan_dft_1d(plan,nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
!!$    ELSE IF(ifb==1) THEN
!!$       CALL dfftw_plan_dft_1d(plan,nx,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
!!$    END IF
!!$
!!$    !This computes the FFT!
!!$    !WRITE(*,*) 'FFT1: Executing FFTW'
!!$    CALL dfftw_execute(plan)
!!$    !WRITE(*,*) 'FFT1: FFTW complete'
!!$
!!$    !And this destroys the plan!
!!$    CALL dfftw_destroy_plan(plan)
!!$    !WRITE(*,*) 'FFT1: Plan destroyed'
!!$    !WRITE(*,*)
!!$
!!$  END SUBROUTINE FFT1

!!$  SUBROUTINE k_FFT1D(i,m,kx,k,L)
!!$
!!$    !Compute the k vector associated with positions in the Fourier mesh
!!$    IMPLICIT NONE
!!$    INTEGER, INTENT(IN) :: i, m
!!$    REAL, INTENT(OUT) :: kx, k
!!$    REAL, INTENT(IN) :: L 
!!$
!!$    !I could change it to work for odd mesh...
!!$    !Note that this would be slightly non-trivial because the middle mode is missing
!!$    IF(MOD(m,2) .NE. 0) STOP 'K_FFT1D: Mesh size needs to be even'
!!$
!!$    kx=float(i-1)
!!$    IF(i>m/2+1) kx=-float(m-i+1)
!!$    kx=kx*2.*pi/L
!!$
!!$    k=ABS(kx)
!!$
!!$  END SUBROUTINE k_FFT1D

!!$  FUNCTION number_file_zeroes(fbase,i,num,fext)
!!$
!!$    !Number a file with zero padding
!!$    !Num specifies the number of digits
!!$    IMPLICIT NONE
!!$    CHARACTER(len=256) :: number_file_zeroes
!!$    CHARACTER(len=256), INTENT(IN) :: fbase, fext
!!$    CHARACTER(len=4) :: num4
!!$    CHARACTER(len=3) :: num3
!!$    CHARACTER(len=2) :: num2
!!$    CHARACTER(len=1) :: num1
!!$    INTEGER, INTENT(IN) :: i
!!$    INTEGER :: maxnum
!!$    INTEGER, INTENT(IN) :: num
!!$
!!$    maxnum=4
!!$
!!$    IF(i<0) STOP 'NUMBER_FILE_ZEROES: Error: cannot write negative number file names'
!!$
!!$    IF(num>maxnum) STOP 'NUMBER_FILE_ZEROES: Error: need to add extra number capacity'
!!$
!!$    IF(num==1) THEN
!!$
!!$       IF(i>=10) STOP 'NUMBER_FILE_ZEROES: Error: i is too large a number for space (i>num)'
!!$       WRITE(num1,fmt='(I0.1)') i
!!$       number_file_zeroes=TRIM(fbase)//num1//TRIM(fext)
!!$
!!$    ELSE IF(num==2) THEN
!!$
!!$       IF(i>=100) STOP 'NUMBER_FILE_ZEROES: Error: i is too large a number for space (i>num)'
!!$       WRITE(num2,fmt='(I0.2)') i
!!$       number_file_zeroes=TRIM(fbase)//num2//TRIM(fext)
!!$
!!$    ELSE IF(num==3) THEN
!!$
!!$       IF(i>=1000) STOP 'NUMBER_FILE_ZEROES: Error: i is too large a number for space (i>num)'
!!$       WRITE(num3,fmt='(I0.3)') i
!!$       number_file_zeroes=TRIM(fbase)//num3//TRIM(fext)
!!$
!!$    ELSE IF(num==4) THEN
!!$
!!$       IF(i>=10000) STOP 'NUMBER_FILE_ZEROES: Error: i is too large a number for space (i>num)'
!!$       WRITE(num4,fmt='(I0.4)') i
!!$       number_file_zeroes=TRIM(fbase)//TRIM(num4)//TRIM(fext)
!!$
!!$    END IF
!!$
!!$  END FUNCTION number_file_zeroes

!!$  SUBROUTINE fill_table(min,max,arr,n)
!!$
!!$    !Fills array 'arr' in equally spaced intervals
!!$    !I'm not sure if inputting an array like this is okay
!!$    IMPLICIT NONE
!!$    INTEGER :: i
!!$    REAL, INTENT(IN) :: min, max
!!$    REAL, ALLOCATABLE, INTENT(INOUT) :: arr(:)
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
!!$  END SUBROUTINE fill_table

!!$  SUBROUTINE fill_table8(min,max,arr,n)
!!$
!!$    !Fills array 'arr' in equally spaced intervals
!!$    !I'm not sure if inputting an array like this is okay
!!$    IMPLICIT NONE
!!$    INTEGER :: i
!!$    REAL, INTENT(IN) :: min, max
!!$    REAL, ALLOCATABLE, INTENT(INOUT) :: arr(:)
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
!!$  END SUBROUTINE fill_table8

  SUBROUTINE Pk1D(d,m,L,outfile)

    IMPLICIT NONE
    REAL, INTENT(IN) :: d(m), L
    INTEGER, INTENT(IN) :: m
    CHARACTER(len=256), INTENT(IN) :: outfile
    DOUBLE COMPLEX :: dc(m)
    REAL :: kx, modk
    REAL :: kmin, kmax
    INTEGER :: nk
    REAL, ALLOCATABLE :: kbin(:)
    REAL, ALLOCATABLE :: k8(:), pk8(:)
    INTEGER*8, ALLOCATABLE :: nbin(:)
    REAL, ALLOCATABLE :: k(:), pk(:), Dk(:)
    INTEGER :: i, j

    !Number of bins for P(k)
    nk=100

    !Minimum and maximum k for binning
    kmin=2.*pi/L
    kmax=2.*pi*REAL(m)/L
    CALL fill_table(log(kmin),log(kmax),kbin,nk+1)
    kbin=exp(kbin)

    !Expand first and last bins to ensure all modes are collected
    kbin(1)=kbin(1)*0.9999
    kbin(nk)=kbin(nk)*1.0001

    !Allocate the counting variables for modes, k and P(k)
    ALLOCATE(k8(nk),Pk8(nk),nbin(nk))
    nbin=0
    k8=0.
    pk8=0.

    !Set up a complex delta
    dc=(0.,0.)
    dc=d

    !Calculate delta_k
    CALL FFT1(dc,dc,m,-1)
    !CALL execute_FFT1(plan_forward)

    !Loop over all k (miss out k=0) and bin
    DO i=2,m
       CALL k_FFT1D(i,m,kx,modk,L)
       DO j=1,nk
          IF(modk>=kbin(j) .AND. modk<kbin(j+1)) THEN
             k8(j)=k8(j)+modk
             pk8(j)=pk8(j)+dc(i)*CONJG(dc(i))
             nbin(j)=nbin(j)+1
             EXIT
          END IF
       END DO
    END DO

    !Allocate arrays for actual k and actual P(k) and actual D(k)
    ALLOCATE(k(nk),pk(nk),Dk(nk))

    OPEN(7,file=outfile)
    DO i=1,nk
       IF(nbin(i)==0) THEN
          !No modes, do not write in this case
          k(i)=sqrt(kbin(i)*kbin(i+1))
          pk(i)=0.
       ELSE
          !Modes, k is average of all contributing k
          !P(k) is average of all contributing P(k)
          k(i)=k8(i)/REAL(nbin(i),8)
          pk(i)=pk8(i)/REAL(nbin(i),8)
          Dk(i)=2.*k(i)*L*pk(i)/(2.*pi)  
          WRITE(7,*) k(i), pk(i), Dk(i), nbin(i)
       END IF
    END DO
    CLOSE(7)

    DEALLOCATE(k,pk,Dk,kbin,k8,pk8,nbin)

  END SUBROUTINE PK1D

  SUBROUTINE histogram(d,m,outfile)

    IMPLICIT NONE
    REAL, INTENT(IN) :: d(m)
    INTEGER, INTENT(IN) :: m
    CHARACTER(len=256) :: outfile
    REAL :: dmin, dmax, delta, width
    INTEGER :: i, j, nd
    REAL, ALLOCATABLE :: nbin(:), dbin(:)

    !Parameters
    INTEGER, PARAMETER :: ilog=0

    IF(ilog==0) THEN
       dmin=0.
       dmax=10.
       nd=100
       width=(dmax-dmin)/REAL(nd)
       CALL fill_table(dmin,dmax,dbin,nd+1)
    ELSE IF(ilog==1) THEN
       dmin=1e-3
       dmax=1e3
       nd=100
       width=log(dmax/dmin)/REAL(nd)
       CALL fill_table(log(dmin),log(dmax),dbin,nd+1)
       dbin=exp(dbin)
    END IF

    ALLOCATE(nbin(nd))
    nbin=0

    !WRITE(*,*) SUM(d)
    !STOP

    DO i=1,m
       DO j=1,nd
          IF(d(i)>=dbin(j) .AND. d(i)<dbin(j+1)) THEN
             nbin(j)=nbin(j)+1
             EXIT
          END IF
       END DO
       !IF(j<=0) THEN
       !   CYCLE
       !ELSE IF(j>nd) THEN
       !   CYCLE
       !ELSE
       !END IF
    END DO

    OPEN(7,file=outfile)
    DO j=1,nd
       IF(ilog==0) THEN
          delta=0.5*(dbin(j)+dbin(j+1))
       ELSE IF(ilog==1) THEN
          delta=sqrt(dbin(j)*dbin(j+1))
       END IF
       WRITE(7,*) REAL(delta), REAL(nbin(j)/(width*REAL(m)))
    END DO
    CLOSE(7)

  END SUBROUTINE histogram

END PROGRAM PM1D
