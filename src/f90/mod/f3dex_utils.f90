MODULE f3dex_utils

  USE OMP_LIB

  USE healpix_types
  USE healpix_modules
  USE alm_tools
  USE fitstools
  IMPLICIT NONE

CONTAINS

  ! ---------------------------------------------------------------------------------------!     

  !> Print preformated messages
  !!@param[in] code : code message
  !!@param[in] (msg,i,msg2,i2) : str/int/str/int chain 
  !!@param[in] start : If present, starting code message
  !!@param[in] fin : If present, ending code message 
  SUBROUTINE message(code, msg, i, i2, msg2, start, fin)
  
     CHARACTER(LEN=*) :: code
     CHARACTER(LEN=100) :: realmsg, conv
     INTEGER(I4B), OPTIONAL :: i, i2
     CHARACTER(LEN=*), OPTIONAL :: msg, msg2
     LOGICAL, OPTIONAL :: start, fin
     
     IF( PRESENT(start) .AND. start .EQV. .TRUE. ) THEN
     	 realmsg = "   Starting routine " // trim(code)
     	 !PRINT*," "
     	 PRINT*,"================================================"
     	 PRINT*,trim(realmsg)
     	 IF( present(msg) )  PRINT*, ("   " // trim(msg))
     	 PRINT*,"------------------------------------------------"
     ELSEIF( PRESENT(fin) .AND. fin .EQV. .TRUE. ) THEN
     	 realmsg = "   " // trim(code) // " terminated"
     	 PRINT*,"------------------------------------------------"
     	 PRINT*,trim(realmsg)
     	 IF( present(msg) )  PRINT*, ("   " // trim(msg))
     	 PRINT*,"================================================"
     	 !PRINT*," "
     ELSE
     	 IF( present(msg) ) THEN
     	    realmsg = trim(code) // " > " // trim(msg)
     	    IF( present(i) .AND. present(msg2) ) THEN
     	       WRITE(conv,'(I5)') i
     	       realmsg = trim(realmsg) // trim(conv) // "   " // trim(msg2)
     	    ENDIF
     	    IF( present(i2) ) THEN
     	       WRITE(conv,'(I5)') i2
     	       realmsg = trim(realmsg) // trim(conv)
     	    ENDIF
     	 ELSE
     	    realmsg = trim(code) 
     	 ENDIF
     	 PRINT*, trim(realmsg)
     ENDIF

  
  END SUBROUTINE 
    
    
  ! ---------------------------------------------------------------------------------------!     

  !> Assert if two doubles are equal
  !!@param[in] (r1,r2) : two doubles
  !!@param[in] prec : precision digits
  LOGICAL FUNCTION assert_DP(r1,r2,prec)

    REAL(DP) :: r1, r2
    INTEGER(I4B) :: prec

    assert_DP = ( floor( abs(r1 - r2)*10**prec ) .EQ. 0 )

    RETURN 
    
  END FUNCTION assert_DP

 ! ---------------------------------------------------------------------------------------!     

  !> Assert if two arrays are equal
  !!@param[in] (r1,r2) : two arrays
  !!@param[in] len : array length
  !!@param[in] prec : precision digits
  LOGICAL FUNCTION assert_DPARR(r1,r2,len,prec)

    INTEGER(I4B) :: prec, i, len
    REAL(DP),DIMENSION(0:len-1)  :: r1, r2
    LOGICAL :: res

    DO i=0,len-1
       !print*,r1(i), r2(i), prec
       res = res .AND. assert_DP(r1(i),r2(i),prec)
    ENDDO

    assert_DPARR = res

    RETURN 
    
  END FUNCTION assert_DPARR

  ! ---------------------------------------------------------------------------------------!     

  !> Classical bubble sort
  !!@param[in] Y : Initial vector
  !!@param[in] N : Array length
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

  !> Raw function doing bubble sort
  !!@param[in] X : Initial vector
  !!@param[in] IY : Temp array
  !!@param[in] N : Array length
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

  !> Measures system time
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

  !> Print spectrum in a preformated way
  !!@param[in] spectr : Power spectrum (two-array)
  !!@param[in] (nllim, nnlim) : Print limit
  !!@param[in] (nlmax, nnmax) : Inner bounds
  !!@param[in] txt : Additional text
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

  !> Print spectrum in a preformated way
  !!@param[in] almn : Power spectrum (two-array)
  !!@param[in] (nllim, nmlim, nnlim) : Print limit
  !!@param[in] (nlmax, nmmax, nnmax) : Inner bounds
  !!@param[in] txt : Additional text
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
