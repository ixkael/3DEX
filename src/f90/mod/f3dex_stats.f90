MODULE f3dex_stats

  USE OMP_LIB

  USE healpix_types
  USE healpix_modules
  USE alm_tools
  USE fitstools
  IMPLICIT NONE

CONTAINS

  ! ---------------------------------------------------------------------------------------!    

  SUBROUTINE almn2cln_naive(nnmax, nlmax, nmmax, almn, cln)

    INTEGER(I4B),                      intent(in) :: nlmax, nmmax, nnmax
    COMPLEX(DPC), dimension(1:nnmax,0:nlmax,0:nmmax), intent(in) :: almn
    REAL(DP)    , dimension(0:nlmax, 1:nnmax),  intent(out):: cln
    ! 
    INTEGER(I4B) :: n, l, ncl, na, mm
    COMPLEX(DPC)     :: dc
    REAL(DP), PARAMETER :: two = 2.000000000000000000_dp
    REAL(DP), PARAMETER :: one = 1.000000000000000000_dp

    ncl = size(cln, 2)
    na = size(almn, 1)
    cln = 0.0_DP

    DO n = 1, nnmax
       DO l = 0, nlmax
          mm = min(l, nmmax)
          dc = sum(almn(n,l,1:mm)*CONJG(almn(n,l,1:mm)))
          dc = (dc + CONJG(dc)) + almn(n,l,0)*almn(n,l,0)
          !PRINT*,l,n,dc,two,one
          cln(l,n) = abs(REAL(dc, KIND=DP)) / abs(two*l + one)

       ENDDO
    ENDDO

    RETURN

    ! ---------------------------------------------------------------------------------------!    

  END SUBROUTINE almn2cln_naive

  FUNCTION dzrel( cplxnb, cplxnbref )

    complex(DP) :: cplxnb, cplxnbref, zero
    REAL(DP) :: nbmod, nbmodref, dzrel

    zero = 0.0_dp
    nbmod = dsqrt( (REAL( cplxnb ))**2.0_dp + &
         & (aimag( cplxnb ))**2.0_dp )
    nbmodref = dsqrt( (REAL( cplxnbref ))**2.0_dp + &
         & (aimag( cplxnbref ))**2.0_dp )

    !PRINT*, nbmod, nbmodref, ( nbmod - nbmodref ), ( nbmod - nbmodref ) / nbmodref
    IF( nbmodref == zero  ) THEN
       dzrel = 0.0_dp
    ENDIF

    IF( ISNAN(nbmodref) .OR. ISNAN(nbmod)  ) THEN
       dzrel = 0.0_dp
    ELSE
       dzrel = abs( nbmod - nbmodref ) / nbmodref
    ENDIF

    IF( ISNAN(dzrel) .AND. nbmodref <= 0.000000000000001 ) THEN
       dzrel = 0.0_dp
    ENDIF

    RETURN 

  END FUNCTION dzrel

  ! ---------------------------------------------------------------------------------------!    

  FUNCTION dzabs( cplxnb, cplxnbref )

    complex(DP) :: cplxnb, cplxnbref
    REAL(DP) :: nbmod, nbmodref, dzabs

    nbmod = dsqrt( (REAL( cplxnb ))**2.0_dp + &
         & (aimag( cplxnb ))**2.0_dp )
    nbmodref = dsqrt( (REAL( cplxnbref ))**2.0_dp + &
         & (aimag( cplxnbref ))**2.0_dp )

    !PRINT*, nbmod, nbmodref, ( nbmod - nbmodref ), ( nbmod - nbmodref ) / nbmodref
    IF( ISNAN(nbmodref) ) THEN
       dzabs = 0.0_dp
    ELSE
       dzabs = abs( nbmod - nbmodref )
    ENDIF

    IF( ISNAN(dzabs) ) THEN
       dzabs = 0.0_dp
    ENDIF

    RETURN 

  END FUNCTION dzabs

  ! ---------------------------------------------------------------------------------------!    


END MODULE f3dex_stats
