PROGRAM IC1D

  USE constants
  USE random_numbers
  USE fft
  USE cosmology_functions
  USE simulations1D

  !Always necessary
  IMPLICIT NONE
  
  CHARACTER(len=256) :: infile, outfile
  REAL, ALLOCATABLE :: x(:), v(:), q(:)
  REAL, ALLOCATABLE :: d(:), s(:)
  INTEGER :: i
  INTEGER :: icosmo
  REAL :: average, max, rms
  TYPE(cosmology) :: cosm

  !Parameters
  REAL, PARAMETER :: Om_m=1. !Omega_m
  REAL, PARAMETER :: Om_v=0. !Omega_v
  REAL, PARAMETER :: L=1. !Box size in Mpc
  REAL, PARAMETER :: a=1. !Starting scale factor
  INTEGER, PARAMETER :: n=2**20 !Number of particles
  INTEGER, PARAMETER :: m=n !Mesh size
  REAL, PARAMETER :: index=2. !Spectral index
  REAL, PARAMETER :: amplitude=1e-6 !Spectral amplitude
  LOGICAL, PARAMETER :: verbose=.TRUE. !Does the code speak?
  LOGICAL, PARAMETER :: mode_average=.TRUE. !Use the exact sigma for P(k) or Rayleigh

  !Mathematical constants
  !REAL, PARAMETER :: pi=3.141592654 !Mathematical pi
  !REAL, PARAMETER :: H0=100. !H0 in h km/s/Mpc

  !Get the output file from the command line
  CALL get_command_argument(1,outfile)
  IF(outfile=='') STOP 'IC1D: Error, specify an output file'

  !Introduce the code
  IF(verbose .EQV. .TRUE.) THEN
     WRITE(*,*)
     WRITE(*,*) 'IC1D: 1D IC generator'
     WRITE(*,*) '====================='
     WRITE(*,*)
  END IF

  icosmo=1
  CALL assign_cosmology(icosmo,cosm)

  !Allocate the initial position, final position and velocity arrays 
  ALLOCATE(q(n),x(n),v(n))

  !Make the arrays for the displacement field the same size as the particle number
  ALLOCATE(d(m),s(m))

  !Make the density and displacement fields
  CALL make_displacement(d,s,m,L,a)

  !Write out useful crap
  IF(verbose .EQV. .TRUE.) THEN
     WRITE(*,*) 'IC1D: Box size [Mpc/h]:', L
     WRITE(*,*) 'IC1D: Number of particles:', n
     WRITE(*,*) 'IC1D: Mesh size for ICs:', m
     WRITE(*,*) 'IC1D: Spectral amplitude:', amplitude
     WRITE(*,*) 'IC1D: Spectral index:', index
     WRITE(*,*)
     average=SUM(ABS(s))/REAL(n)
     rms=sqrt(SUM(s**2)/REAL(n))
     max=MAXVAL(ABS(s))
     WRITE(*,*) 'IC1D: Average modulus displacement [Mpc/h]:', average
     WRITE(*,*) 'IC1D: RMS displacement [Mpc/h]:', rms
     WRITE(*,*) 'IC1D: Max modulus displacement [Mpc/h]:', max
     WRITE(*,*)
     WRITE(*,*) 'IC1D: Average modulus displacement IPS:', average/(L/REAL(n))
     WRITE(*,*) 'IC1D: RMS displacement IPS:', rms/(L/REAL(n))
     WRITE(*,*) 'IC1D: Max modulus displacement IPS:', max/(L/REAL(n))
     WRITE(*,*)
  END IF

  !Make a uniform grid of particles
  CALL make_grid(q,n,L)

  !Displace the particles
  CALL displace(q,x,v,s,n,a,cosm)

  !DO i=1,10
  !   WRITE(*,*) q(i), x(i), v(i)
  !END DO
  !STOP
  
  !Check for crossings
  CALL check(x,n)

  !Enforce the line periodicity
  CALL enforce_periodicity(x,n,L)
  
  !Write the particle data to a file
  CALL write_particles(x,v,n,outfile,verbose)
  
CONTAINS

  SUBROUTINE displace(q,x,v,s,n,a,cosm)

    !Displace particles according to the dispalcement field
    IMPLICIT NONE
    REAL, INTENT(IN) :: q(n), s(n)
    REAL, INTENT(OUT) ::  x(n), v(n)
    INTEGER, INTENT(IN) :: n
    REAL, INTENT(IN) :: a
    TYPE(cosmology), INTENT(INOUT) :: cosm 
    INTEGER :: i

    DO i=1,n
       x(i)=q(i)+s(i)
       v(i)=s(i)*sqrt(Hubble2(a,cosm))*a
    END DO

  END SUBROUTINE displace

  SUBROUTINE check(x,n)

    !Check that no crossings have occured
    IMPLICIT NONE
    REAL, INTENT(IN) ::  x(n)
    INTEGER, INTENT(IN) :: n
    INTEGER :: i

    DO i=1,n-1
       IF(x(i+1)<x(i)) THEN
          WRITE(*,*) 'Particle:', i
          WRITE(*,*) 'Position [Mpc/h:', x(i)
          WRITE(*,*) 'Particle:', i+1
          WRITE(*,*) 'Position [Mpc/h]:', x(i+1)
          STOP 'CHECK: Error, particles have crossed'
       END IF
    END DO

  END SUBROUTINE check

  SUBROUTINE make_displacement(d,s,m,L,a)

    !Make a Gaussian-random displacement field
    IMPLICIT NONE
    REAL, INTENT(OUT) :: d(m), s(m)
    INTEGER, INTENT(IN) :: m
    REAL, INTENT(IN) :: L, a
    INTEGER :: i
    REAL :: kx, k
    REAL :: R, phase, x
    DOUBLE COMPLEX, ALLOCATABLE :: dc(:), sc(:)

    !Parameters
    INTEGER, PARAMETER :: seed=0 !Seed for RNG

    ALLOCATE(dc(m),sc(m))

    CALL RNG_set(seed)

    !Make the fields in Fourier space
    DO i=1,m
       IF(i==1) THEN
          !Set the zero mode to zero
          d(i)=0.
       ELSE
          !CALL k_FFT(2,1,1,m,kx,ky,kz,k,L)
          CALL k_FFT1D(2,m,kx,k,L)
          x=spectrum(k,L,a,cosm)
          IF(mode_average .EQV. .TRUE.) THEN
             R=sqrt(pi/2.)*sqrt(x) !I think the sqrt(pi/2) is correct here
          ELSE  
             R=random_Rayleigh(sqrt(x))
          END IF
          phase=random_uniform(0.,2.*pi)
          dc(i)=R*exp((0.d0,1.d0)*phase)
          sc(i)=-(0.d0,1.d0)*kx*dc(i)/k**2
       END IF
    END DO

    !Convert to real space
    CALL FFT1(sc,sc,m,1)

    !Normalise post FFT (is this correct/necessary?)
    s=REAL(sc)/REAL(m)
    
  END SUBROUTINE make_displacement

!!$  SUBROUTINE enforce_periodicity(x,n,L)
!!$
!!$    !Enfoce the periodicity of the line
!!$    IMPLICIT NONE
!!$    REAL, INTENT(INOUT) :: x(n)
!!$    INTEGER, INTENT(IN) :: n
!!$    REAL, INTENT(IN) :: L
!!$    INTEGER :: i, m
!!$
!!$    DO i=1,n
!!$       
!!$       IF(x(i)<0.) THEN
!!$          m=CEILING(-x(i)/L)
!!$          x(i)=x(i)+REAL(m)*L
!!$       ELSE IF(x(i)>L) THEN
!!$          m=CEILING(x(i)/L)-1
!!$          x(i)=x(i)-REAL(m)*L
!!$       END IF
!!$
!!$       IF(x(i)==0.) x=L
!!$
!!$    END DO
!!$    
!!$  END SUBROUTINE enforce_periodicity

  FUNCTION spectrum(k,L,a,cosm)

    IMPLICIT NONE
    REAL :: spectrum
    REAL, INTENT(IN) :: k, a, L
    TYPE(cosmology), INTENT(INOUT) :: cosm

    spectrum=0.5*amplitude*(grow(a,cosm)**2)*(k*L/(2.*pi))**index
    
  END FUNCTION spectrum

!!$  FUNCTION growth(a)
!!$
!!$    IMPLICIT NONE
!!$    REAL :: growth
!!$    REAL, INTENT(IN) :: a
!!$
!!$    growth=a
!!$    
!!$  END FUNCTION growth

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

!!$  FUNCTION Hubble(a)
!!$
!!$    IMPLICIT NONE
!!$    REAL :: Hubble
!!$    REAL, INTENT(IN) :: a
!!$
!!$    Hubble=H0*sqrt(Om_m*a**(-3)+Om_v+(1.-Om_m-Om_v)*a**(-2))
!!$
!!$  END FUNCTION Hubble

END PROGRAM IC1D
