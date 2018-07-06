MODULE sorting

  INTERFACE bubble_index

     MODULE PROCEDURE bubble_index_real
     MODULE PROCEDURE bubble_index_int
     
  END INTERFACE bubble_index

  INTERFACE stupid_index

     MODULE PROCEDURE stupid_index_real
     MODULE PROCEDURE stupid_index_int
     
  END INTERFACE stupid_index

CONTAINS

  SUBROUTINE bubble_sort(a,n)

    IMPLICIT NONE
    REAL, INTENT(INOUT) :: a(n)
    INTEGER, INTENT(IN) :: n
    REAL :: hold
    INTEGER :: i, isort

    !This is a dodgy module!
    !I think it overwrites some of the things
    !in the original array somehow!
    STOP 'BUBBLE_SORT: This is dodgy, see the source code'

    DO
       isort=0
       DO i=1,n-1
          IF(a(i)>a(i+1)) THEN
             hold=a(i+1)
             a(i+1)=a(i)
             a(i)=hold
             isort=1
          END IF
       END DO
       IF(isort==0) EXIT
    END DO

  END SUBROUTINE bubble_sort

  SUBROUTINE selection_sort(a,n)

    IMPLICIT NONE
    REAL, INTENT(INOUT) :: a(n)
    INTEGER, INTENT(IN) :: n    
    REAL :: hold, min
    INTEGER :: i, j, minl

    DO i=1,n-1
       min=a(i)
       minl=i
       DO j=i+1,n
          IF(a(j)<min) THEN
             min=a(j)
             minl=j
          END IF
       END DO
       hold=a(i)
       a(i)=min
       a(minl)=hold
    END DO

  END SUBROUTINE selection_sort

  SUBROUTINE bubble_index_real(a,ind,n)!,verbose)

    ! Create an index array for a(:) that indexes from smallest to largest value
    IMPLICIT NONE    
    REAL, INTENT(IN) :: a(n)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: ind(n)
!!$    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i, isort, hold

!!$    IF(verbose) WRITE(*,*) 'BUBBLE_INDEX: Starting bubble index'

    DO i=1,n
       ind(i)=i
    END DO

    DO
       isort=0
       DO i=1,n-1
          IF(a(ind(i))>a(ind(i+1))) THEN
             hold=ind(i+1)
             ind(i+1)=ind(i)
             ind(i)=hold
             isort=1
          END IF
       END DO
       IF(isort==0) EXIT
    END DO

!!$    IF(verbose) THEN
!!$       WRITE(*,*) 'BUBBLE_INDEX: Bubble index finished'
!!$       WRITE(*,*)
!!$    END IF
      
  END SUBROUTINE bubble_index_real

  SUBROUTINE bubble_index_int(a,ind,n)!,verbose)

    ! Create an index array for integer a(:) that indexes from smallest to largest value
    IMPLICIT NONE   
    INTEGER, INTENT(IN) :: a(n), n
    INTEGER, INTENT(OUT) :: ind(n)
!!$    LOGICAL, INTENT(IN) :: verbose
    INTEGER :: i, isort, hold

!!$    IF(verbose) WRITE(*,*) 'Starting bubble index'

    DO i=1,n
       ind(i)=i
    END DO

    DO
       isort=0
       DO i=1,n-1
          IF(a(ind(i))>a(ind(i+1))) THEN
             hold=ind(i+1)
             ind(i+1)=ind(i)
             ind(i)=hold
             isort=1
          END IF
       END DO
       IF(isort==0) EXIT
    END DO

!!$    IF(verbose) THEN
!!$       WRITE(*,*) 'Bubble index finished'
!!$       WRITE(*,*)
!!$    END IF

  END SUBROUTINE bubble_index_int

  SUBROUTINE stupid_index_real(a,ind,n)

    IMPLICIT NONE
    REAL, INTENT(IN) :: a(n)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: ind(n)
    INTEGER :: i, j
    REAL :: b(n)

    b=a

    ! This is probably stupid
    DO i=1,n
       j=MINLOC(b,1)
       ind(i)=j
       b(j)=HUGE(b)
    END DO

  END SUBROUTINE stupid_index_real

  SUBROUTINE stupid_index_int(a,ind,n)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: a(n), n
    INTEGER, INTENT(OUT) :: ind(n)
    INTEGER :: i, j, b(n)

    b=a

    ! This is probably stupid
    DO i=1,n
       j=MINLOC(b,1)
       ind(i)=j
       b(j)=HUGE(b)
    END DO

  END SUBROUTINE stupid_index_int
  
END MODULE sorting
