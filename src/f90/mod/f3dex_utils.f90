MODULE f3dex_utils

  USE OMP_LIB

  USE healpix_types
  USE healpix_modules
  USE alm_tools
  USE fitstools
  IMPLICIT NONE

CONTAINS

  ! ---------------------------------------------------------------------------------------!     

  SUBROUTINE bubblesort(Y, N) 

    INTEGER(kind=I4B)  :: N, i
    REAL(kind=DP), DIMENSION(1:N) :: Y
    REAL(DP), DIMENSION(:), ALLOCATABLE :: R

    ALLOCATE(R(1:N))

    R = Y

    CALL BSORT(Y, R, N)

    DO i=1,N
       Y(N-i+1) = R(i)
    ENDDO

  END SUBROUTINE bubblesort

  ! ---------------------------------------------------------------------------------------!     

  SUBROUTINE BSORT (X, IY, N)

    IMPLICIT NONE    

    INTEGER(kind=I4B) 		      :: N, I, J, JMAX
    REAL(kind=DP), DIMENSION(1:N)     :: X
    REAL(kind=DP), DIMENSION(1:N,1:3) :: IY
    REAL(kind=DP)                     :: TEMP
    REAL(kind=DP), DIMENSION(1:3)     ::ITEMP

    JMAX=N-1
    DO 200 I=1,N-1
       !PRINT*,"     SORT : ",i,"/",N
       TEMP=1.E38
       DO 100 J=1,JMAX
          IF(X(J).GT.X(J+1)) GO TO 100
          TEMP=X(J)
          X(J)=X(J+1)
          X(J+1)=TEMP
          ITEMP=IY(J,:)
          IY(J,:)=IY(J+1,:)
          IY(J+1,:)=ITEMP
100       CONTINUE
          IF(TEMP.EQ.1.E38) GO TO 300	       
          JMAX=JMAX-1
200       CONTINUE

300 RETURN

  END SUBROUTINE BSORT

  ! ---------------------------------------------------------------------------------------!     

  FUNCTION wtime( )

    implicit none

    INTEGER ( kind = 4 ) clock_max
    INTEGER ( kind = 4 ) clock_rate
    INTEGER ( kind = 4 ) clock_reading
    REAL    ( kind = 8 ) wtime

    CALL system_clock ( clock_reading, clock_rate, clock_max )

    wtime = REAL ( clock_reading, kind = 8 ) &
         / REAL ( clock_rate, kind = 8 )

    RETURN
  END FUNCTION wtime

  ! ---------------------------------------------------------------------------------------!    

  SUBROUTINE PRINT_spectrum(spectr, nllim, nnlim, nlmax, nnmax, txt)

    INTEGER(I4B) :: nllim, nnlim, l, n, nnmax, nlmax
    REAL(DP), intent(IN), DIMENSION(0:nlmax,1:nnmax) :: spectr
    CHARACTER(len=20) :: nm
    CHARACTER(len=4) :: txt

    PRINT*,">> ",txt,"(l,n)"
111 FORMAT (A,A,I3,A,I3,A)  
    DO n = 1, nnlim
       DO l = 0, nllim
          WRITE(nm,111) trim(txt),'(', l,",", n,")"
          PRINT*,nm, spectr(l,n)
       ENDDO
    ENDDO


  END SUBROUTINE PRINT_spectrum

  ! ---------------------------------------------------------------------------------------!     

  SUBROUTINE PRINT_almn(almn, nllim, nmlim, nnlim, nlmax, nmmax, nnmax, txt)

    INTEGER(I4B) :: nllim, nmlim, nnlim, l, m, n, nnmax, nlmax, nmmax
    COMPLEX(DPC), intent(IN), DIMENSION(1:nnmax,0:nlmax,0:nmmax) :: almn
    CHARACTER(len=20) :: nm
    CHARACTER(len=4) :: txt

    PRINT*,">> ",txt,"(n,l,m)"
110 FORMAT (A,A,I3,A,I3,A,I3,A)    	
    DO n = 1, nnlim
       DO l = 0, nllim
          DO m = 0, l
             WRITE(nm,110) trim(txt),'(', n,",",l,",", m,")"
             PRINT*,nm, almn(n,l,m)
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE PRINT_almn

  ! ---------------------------------------------------------------------------------------!     

END MODULE f3dex_utils
