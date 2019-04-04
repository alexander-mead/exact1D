MODULE sorting

  INTERFACE index
     MODULE PROCEDURE index_real
     MODULE PROCEDURE index_int     
  END INTERFACE index
  
  INTERFACE bubble_index
     MODULE PROCEDURE bubble_index_real
     MODULE PROCEDURE bubble_index_int     
  END INTERFACE bubble_index

  INTERFACE stupid_index
     MODULE PROCEDURE stupid_index_real
     MODULE PROCEDURE stupid_index_int     
  END INTERFACE stupid_index

CONTAINS

  SUBROUTINE sort(a,n,imeth)

    IMPLICIT NONE
    REAL, INTENT(INOUT) :: a(n)
    INTEGER, INTENT(IN) :: n, imeth

    IF(imeth==1) THEN
       CALL stupid_sort(a,n)
    ELSE IF(imeth==2) THEN
       CALL bubble_sort(a,n)
    ELSE IF(imeth==3) THEN
       CALL QsortC(a)
    ELSE
       STOP 'SORT: Error, imeth not specified correctly'
    END IF
    
  END SUBROUTINE sort

  SUBROUTINE bubble_sort(a,n)

    ! Bubble sort array 'a' into lowest to highest value
    IMPLICIT NONE
    REAL, INTENT(INOUT) :: a(n)
    INTEGER, INTENT(IN) :: n
    REAL :: hold
    INTEGER :: i
    LOGICAL :: sorted

    DO
       sorted=.TRUE.
       DO i=1,n-1
          IF(a(i)>a(i+1)) THEN
             hold=a(i+1)
             a(i+1)=a(i)
             a(i)=hold
             sorted=.FALSE.
          END IF
       END DO
       IF(sorted) EXIT
    END DO

  END SUBROUTINE bubble_sort

  SUBROUTINE stupid_sort(a,n)

    ! I have no idea what this is
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

  END SUBROUTINE stupid_sort

   RECURSIVE SUBROUTINE QsortC(A)

    ! Stolen from http://www.fortran.com/qsort_c.f95
    IMPLICIT NONE
    REAL, INTENT(INOUT), DIMENSION(:) :: A
    INTEGER :: iq

    if(size(A) > 1) then
       call Partition(A, iq)
       call QsortC(A(:iq-1))
       call QsortC(A(iq:))
    endif
    
  END SUBROUTINE QsortC

  SUBROUTINE Partition(A, marker)

    ! Stolen from http://www.fortran.com/qsort_c.f95
    IMPLICIT NONE
    REAL, INTENT(in out), DIMENSION(:) :: A
    INTEGER, INTENT(out) :: marker
    INTEGER :: i, j
    REAL :: temp
    REAL :: x      ! pivot point
    
    x = A(1)
    i= 0
    j= size(A) + 1

    do
       j = j-1
       do
          if (A(j) <= x) exit
          j = j-1
       end do
       i = i+1
       do
          if (A(i) >= x) exit
          i = i+1
       end do
       if (i < j) then
          ! exchange A(i) and A(j)
          temp = A(i)
          A(i) = A(j)
          A(j) = temp
       elseif (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do

  END SUBROUTINE Partition

  SUBROUTINE index_real(a,ind,n,imeth)
    
    ! Index the array 'a' from lowest to highest value
    IMPLICIT NONE
    REAL, INTENT(IN) :: a(n)
    INTEGER, INTENT(IN) :: n, imeth
    INTEGER, INTENT(OUT) :: ind(n)

    IF(imeth==1) THEN
       CALL stupid_index_real(a,ind,n)
    ELSE IF(imeth==2) THEN
       CALL bubble_index_real(a,ind,n)
    ELSE
       STOP 'INDEX_REAL: Error, imeth specified incorrectly'
    END IF

  END SUBROUTINE index_real

  SUBROUTINE index_int(a,ind,n,imeth)
    
    ! Index the array 'a' from lowest to highest value
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n, imeth, a(n)
    INTEGER, INTENT(OUT) :: ind(n)

    IF(imeth==1) THEN
       CALL stupid_index_int(a,ind,n)
    ELSE IF(imeth==2) THEN
       CALL bubble_index_int(a,ind,n)
    ELSE
       STOP 'INDEX_REAL: Error, imeth specified incorrectly'
    END IF

  END SUBROUTINE index_int

  SUBROUTINE bubble_index_real(a,ind,n)

    ! Create an index array for a(:) that indexes from smallest to largest value
    IMPLICIT NONE    
    REAL, INTENT(IN) :: a(n)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(OUT) :: ind(n)
    INTEGER :: i, isort, hold

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
      
  END SUBROUTINE bubble_index_real

  SUBROUTINE bubble_index_int(a,ind,n)!,verbose)

    ! Create an index array for integer a(:) that indexes from smallest to largest value
    IMPLICIT NONE   
    INTEGER, INTENT(IN) :: a(n), n
    INTEGER, INTENT(OUT) :: ind(n)
    INTEGER :: i, isort, hold

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

  LOGICAL FUNCTION check_sorted(a,n)

    ! Checks if array 'a' is sorted from highest to lowest
    IMPLICIT NONE
    REAL, INTENT(IN) :: a(n) ! Input array to check
    INTEGER, INTENT(IN) :: n ! Number of entried in array
    INTEGER :: i

    check_sorted=.TRUE.
    
    DO i=1,n-1
       IF(a(i)>a(i+1)) THEN
          check_sorted=.FALSE.
          EXIT
       END IF
    END DO
    
  END FUNCTION check_sorted

  LOGICAL FUNCTION check_sorted_index(a,j,n)

    ! Checks if array indices for 'a' are sorted from highest to lowest
    IMPLICIT NONE
    REAL, INTENT(IN) :: a(n)    ! Input array to check
    INTEGER, INTENT(IN) :: j(n) ! Input array indices to check
    INTEGER, INTENT(IN) :: n    ! Number of entried in array
    INTEGER :: i

    check_sorted_index=.TRUE.
    
    DO i=1,n-1
       IF(a(j(i))>a(j(i+1))) THEN
          check_sorted_index=.FALSE.
          EXIT
       END IF
    END DO

  END FUNCTION check_sorted_index

  SUBROUTINE reindex(a,j,n)

    ! Reindex the array 'a' with the new indices 'j'
    REAL, INTENT(INOUT) :: a(n) ! Input array to check
    INTEGER, INTENT(IN) :: j(n) ! Input array indices to check
    INTEGER, INTENT(IN) :: n    ! Number of entried in array
    REAL :: b(n)
    INTEGER :: i
    
    b=a ! Store the input array
     
    a=0. ! Delete the input array

    ! Loop over values and reindex
    DO i=1,n
       a(i)=b(j(i))
    END DO
    
  END SUBROUTINE reindex
  
END MODULE sorting
