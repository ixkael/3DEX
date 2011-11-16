MODULE f3dex_stats

  USE OMP_LIB

  USE f3dex_utils
  USE healpix_types
  USE healpix_modules
  USE alm_tools
  USE fitstools

  IMPLICIT NONE

CONTAINS

  ! ---------------------------------------------------------------------------------------!    

  !> Computes the module of a complex number
  REAL(DP) FUNCTION CMODULE(a)

    COMPLEX(DPC) :: a

    CMODULE = (REAL(a))**2 + (AIMAG(a))**2

  END FUNCTION CMODULE

  ! ---------------------------------------------------------------------------------------!    

  !> Naive estimator for the almn's
  !!@param[in] alm(0:nnmax,0:nlmax,0:nmmax)) : input almn's
  !!@param[in] (nnmax, nlmax, nmmax) : bounds
  !!@param[out] cl(0:nlmax,0:nnmax) : output power spectrum 
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
          cln(l,n) = abs(REAL(dc, KIND=DP)) / abs(two*l + one)
       ENDDO
    ENDDO

    RETURN

  END SUBROUTINE almn2cln_naive

  ! ---------------------------------------------------------------------------------------!   

  !> Relative error betwee two complex numbers
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

  !> Computes legendre coefs of a double using recurrence
  !!@param[in] x(0:nc-1) : abs
  !!@param[out] pl(0:nc-1,0:nlmax) : legendre matrix
  !!@paramts[in] (nc,nlmax) : dimensions
  SUBROUTINE legendre_vec(x,pl,nlmax)

    REAL(DP), DIMENSION(0:nlmax) :: pl
    REAL(DP) :: x
    INTEGER(I4B) :: nlmax, l

    pl(0) = 1.0
    pl(1) = x

    DO l=2, nlmax
       pl(l) = ( (2.0*REAL(l)-1.0)*x*pl(l-1) - REAL(l-1)*pl(l-2) ) / (REAL(l))
    ENDDO

  END SUBROUTINE legendre_vec

  ! ---------------------------------------------------------------------------------------!    

  !> Computes legendre coefs of an array using recurrence
  !!@param[in] x(0:nc-1) : abs
  !!@param[out] pl(0:nc-1,0:nlmax) : legendre matrix
  !!@paramts[in] (nc,nlmax) : dimensions
  SUBROUTINE legendre_mat(x,pl,nc,nlmax)

    REAL(DP), DIMENSION(0:nc-1, 0:nlmax) :: pl
    REAL(DP), DIMENSION(0:nc-1) :: x
    INTEGER(I4B) :: nc, nlmax, p

    DO p=0, nc-1
       CALL legendre_vec( x(p), pl(p,0:nlmax), nlmax )
    ENDDO

  END SUBROUTINE legendre_mat

  ! ---------------------------------------------------------------------------------------!    

  !> Generates <xxt> correlation angles map
  !!@param[in] (npix, nside) : number of pixels and Healpix parameter
  !!@param[out] ang_map(0:npix-1,0:npix-1) : angular map
  !!@param[in] nested : logical for Healpix mode
  !!@param[in] mask : optional mask
  SUBROUTINE gen_angles_map(npix, ang_map, nside, nested, mask)

    INTEGER(I4B) :: status, p, t, nside, loc, npix
    LOGICAL, OPTIONAL :: nested
    REAL(DP), DIMENSION(0:npix-1,0:npix-1) :: ang_map
    LOGICAL, DIMENSION(0:12*nside**2-1), OPTIONAL :: mask
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: thetaphis
    CHARACTER(LEN=*), PARAMETER :: code = 'gen_angles_map'
    REAL(DP) :: phi, theta1, theta2, phi1, phi2

    ALLOCATE(thetaphis(0:npix-1,1:2),stat = status)
    CALL assert_alloc(status,code,'thetas')

    loc = 0
    IF( present(nested) .AND. nested .EQV. .TRUE.) THEN
       DO p = 0, 12*nside*nside-1
          IF( present(mask) ) THEN
             IF( mask(p) .EQV. .TRUE. ) THEN
                CALL pix2ang_nest(nside, p, thetaphis(loc,1),thetaphis(loc,2) )
                loc = loc + 1
             ENDIF
          ELSE
             CALL pix2ang_nest(nside, p, thetaphis(p,1),thetaphis(p,2) )
          ENDIF
       ENDDO
    ELSE
       DO p = 0, 12*nside*nside-1
          IF( present(mask) ) THEN
             IF( mask(p) .EQV. .TRUE. ) THEN
                CALL pix2ang_ring(nside, p, thetaphis(loc,1),thetaphis(loc,2) )
                loc = loc + 1
             ENDIF
          ELSE
             CALL pix2ang_ring(nside, p, thetaphis(p,1),thetaphis(p,2) )
          ENDIF
       ENDDO
    ENDIF

    DO p = 0, npix-1
       DO t = 0, npix-1
          theta1 = thetaphis(p,1)
          theta2 = thetaphis(t,1)
          phi1 = thetaphis(p,2)
          phi2 = thetaphis(t,2)
          ang_map(p,t) = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2)
       ENDDO
    ENDDO

    DEALLOCATE(thetaphis)

  END SUBROUTINE gen_angles_map

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

  !> Creates eye matrix
  !!@param[out] A(s:e,s:e) : input matrix
  !!@param[in] (s,e) : bounds
  !!@paramts[in] val : value to put on the diagonal
  SUBROUTINE EYEMAT(A,s,e,val)

    INTEGER :: s, e, i
    REAL(DP), DIMENSION(s:e,s:e) :: A
    REAL(DP), OPTIONAL :: val
    IF(.NOT. present(val))then
       val=1.0
    ENDIF

    A = 0.0
    DO i=s,e
       A(i,i) = val
    ENDDO

  END SUBROUTINE EYEMAT

  ! ---------------------------------------------------------------------------------------!     

  !> Inverts matrix
  !!@param[inout] A(n,n) : input/output matrix
  !!@param[in] n : bounds
  !!@paramts[out] res : output result logical	  
  SUBROUTINE INVERSE(A,n,res)

    IMPLICIT NONE
    LOGICAL res
    INTEGER :: n, prec, i, j, IPIV(n), INFO
    DOUBLE PRECISION A(n,n)
    DOUBLE PRECISION WORK(n*n)

    CALL DGETRF( n, n, A, n, IPIV, INFO )
    CALL DGETRI( n, A, n, IPIV, WORK, n*n, INFO )

    IF( INFO .EQ. 0 ) THEN
       res=.TRUE.
    ELSE
       res=.FALSE.
    ENDIF

    RETURN

  END SUBROUTINE INVERSE


  ! ---------------------------------------------------------------------------------------!    

  !> Computes QML estimator of a healpix (masked) map
  !!@param[in] map(0:12*nside*nside-1) : input healpix map
  !!@param[out] cl(0:nlmax) : output power spectrum estimates
  !!@param[in] nside : healpix parameter
  !!@param[in] nlmax : input nlmax for the theory signal covariance 
  !!@param[in] nlmax_recon : lmax for the reconstruction and estimation
  !!@param[in] covarnoise : covariance matrix of the noise
  !!@param[in] cltheo : theoretical power spectrum
  !!@param[in] logicalmask(0:12*nside*nside-1) : (optional) healpix logical mask
  !!@param[out] covarest : (optional) output covariance of the noise 
  SUBROUTINE map2cl_QML( map, cl, nside, nlmax, nlmax_recon, covarnoise, cltheo, logicalmask, covarest )

    INTEGER(I4B) :: nside, nlmax, l, ll, locl, locll, i, ii, npix
    INTEGER(I4B) :: status, nlmax_recon
    LOGICAL :: res
    REAL(DP), DIMENSION(0:(12*nside*nside-1)) :: map
    REAL(DP), DIMENSION(0:nlmax) :: cltheo
    REAL(DP), DIMENSION(0:nlmax_recon) :: cl
    REAL(DP), DIMENSION(0:nlmax_recon), OPTIONAL :: covarest
    LOGICAL, DIMENSION(0:(12*nside*nside-1)), OPTIONAL :: logicalmask
    INTEGER(I4B), DIMENSION(0:(12*nside*nside-1)) :: subscr
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: El, F, covar, ang_map
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Pl
    REAL(DP), DIMENSION(:), ALLOCATABLE :: x, y
    REAL(DP), DIMENSION(0:12*nside**2-1,0:12*nside**2-1) :: covarnoise
    CHARACTER(LEN=*), PARAMETER :: code = 'map2cl_QML'
    REAL(DP) :: Fll, mapvec, trnoise, val                    

    CALL message(code, start=.TRUE.)

    IF( present(logicalmask) ) THEN
       npix = COUNT(logicalmask)
    ELSE
       npix = 12*nside*nside
    ENDIF
    CALL message(code, msg="Initial map has",i=npix,msg2="pixels")

    ALLOCATE(x(0:npix-1),stat = status)
    CALL assert_alloc(status,code,'x')
    x = PACK(map, mask=logicalmask)

    CALL message(code, msg="Generating masked angular map")
    ALLOCATE(ang_map(0:npix-1,0:npix-1),stat = status)
    CALL assert_alloc(status,code,'ang_map')
    CALL gen_angles_map(npix, ang_map, nside, nested=.FALSE., mask=logicalmask) 

    ALLOCATE(Pl(0:npix-1,0:npix-1,0:nlmax),stat = status)
    CALL assert_alloc(status,code,'Pl')
    DO l=0, nlmax
       CALL message(code, msg="Legendre matrices :", i=l, msg2="/", i2=nlmax)
       IF( l .EQ. 0 )THEN
          Pl(0:npix-1,0:npix-1,0) = 1.0_dp
       ELSE IF( l .EQ. 1 ) THEN
          Pl(0:npix-1,0:npix-1,1) = ang_map
       ELSE
          Pl(0:npix-1,0:npix-1,l) = &
               & ( (2.0_dp*REAL(l)-1.0_dp) * ang_map(0:npix-1,0:npix-1) * &
               & Pl(0:npix-1,0:npix-1,l-1) / ( (2.0_dp*REAL(l-1)+1.0_dp)/FOURPI ) &
               & - REAL(l-1) * Pl(0:npix-1,0:npix-1,l-2) / ((2.0_dp*REAL(l-2)+1.0_dp)/FOURPI) ) / REAL(l)
       END IF
       Pl(0:npix-1,0:npix-1,l) = Pl(0:npix-1,0:npix-1,l) * (2.0_dp*REAL(l)+1.0_dp)/FOURPI
    ENDDO
    DEALLOCATE( ang_map )

    CALL message(code, msg="Creating signal covariance matrix")
    ALLOCATE(covar(0:npix-1,0:npix-1),stat = status)
    CALL assert_alloc(status,code,'covar')
    covar = 0.0
    DO l=0, nlmax
       covar = covar + Pl(0:npix-1,0:npix-1,l)*cltheo(l)
    ENDDO

    CALL message(code, msg="Adding noise to signal covariance matrix")
    covar(0:npix-1,0:npix-1) = covar(0:npix-1,0:npix-1) + covarnoise(0:npix-1,0:npix-1)

    CALL message(code, msg="Inverting signal+noise matrix")
    CALL INVERSE(covar(0:npix-1,0:npix-1),npix,res)
    IF( res .EQV. .FALSE. )THEN
       PRINT*,"PROBLEM WITH INVERSION"
    ENDIF

    ALLOCATE(F(0:nlmax_recon,0:nlmax_recon),stat = status)
    CALL assert_alloc(status,code,'F')
    ALLOCATE(y(0:nlmax_recon),stat = status)
    CALL assert_alloc(status,code,'y')

    !$OMP PARALLEL DEFAULT(none) &
    !$OMP SHARED(nlmax_recon, map, covar, cl, Pl, npix, F, y, x, covarest) & 
    !$OMP PRIVATE(l, ll, Fll, El, status)

    ALLOCATE(El(0:npix-1,0:npix-1),stat = status)
    CALL assert_alloc(status,code,'El')  

    !$OMP DO
    DO l=0, nlmax_recon
       CALL message(code, msg="Power spectrum estimation :", i=l, msg2="/", i2=nlmax_recon)
       El(0:npix-1,0:npix-1) = MATMUL( covar(0:npix-1,0:npix-1), &
            & MATMUL( Pl(0:npix-1,0:npix-1,l), covar(0:npix-1,0:npix-1) ) )
       DO ll=0, nlmax_recon
          Fll = TRACEPROD(El(0:npix-1,0:npix-1),Pl(0:npix-1,0:npix-1,ll),0,npix-1)
          F(ll,l) = Fll    
       ENDDO
       y(l) = DOT_PRODUCT( x, MATMUL( El(0:npix-1,0:npix-1), x ) )
    ENDDO
    !$OMP END DO

    DEALLOCATE(El)

    !$OMP END PARALLEL

    CALL message(code, msg="Inverting Fisher matrix")
    CALL INVERSE(F(0:nlmax_recon,0:nlmax_recon),nlmax_recon+1,res)
    IF( res .EQV. .FALSE. )THEN
       PRINT*,"PROBLEM WITH INVERSION"
    ENDIF

    CALL message(code, msg="Evaluating final estimates")
    cl = MATMUL( F(0:nlmax_recon,0:nlmax_recon), y(0:nlmax_recon) )  

    IF( present(covarest) ) THEN
       DO l=0, nlmax_recon
          covarest(l) = F(l,l)
       ENDDO
    ENDIF

    DEALLOCATE(F, covar, Pl, x, y)

    CALL message(code, fin=.TRUE.)

  END SUBROUTINE map2cl_QML

  ! ---------------------------------------------------------------------------------------!    

  !> Computes QML estimator of a healpix (masked) map
  !!@param[in] map(0:12*nside*nside-1) : input healpix map
  !!@param[out] cl(0:nlmax) : output power spectrum estimates
  !!@param[in] nside : healpix parameter
  !!@param[in] nlmax : input nlmax for the theory signal covariance 
  !!@param[in] nlmax_recon : lmax for the reconstruction and estimation
  !!@param[in] logicalmask(0:12*nside*nside-1) : (optional) healpix logical mask
  !!@param[out] covarest : (optional) output covariance of the noise 
  SUBROUTINE map2cl_PCL( map, cl, nside, nlmax, nlmax_recon, logicalmask, covarest )

    INTEGER(I4B) :: nside, nlmax, l, ll, locl, locll, i, ii, npix
    INTEGER(I4B) :: status, nlmax_recon
    LOGICAL :: res
    REAL(DP), DIMENSION(0:(12*nside*nside-1)) :: map
    REAL(DP), DIMENSION(0:nlmax_recon) :: cl
    REAL(DP), DIMENSION(0:nlmax_recon), OPTIONAL :: covarest
    LOGICAL, DIMENSION(0:(12*nside*nside-1)), OPTIONAL :: logicalmask
    INTEGER(I4B), DIMENSION(0:(12*nside*nside-1)) :: subscr
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: El, F, ang_map
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Pl
    REAL(DP), DIMENSION(:), ALLOCATABLE :: x, y
    CHARACTER(LEN=*), PARAMETER :: code = 'map2cl_PCL'
    REAL(DP) :: Fll, mapvec, trnoise, val                    

    CALL message(code, start=.TRUE.)

    IF( present(logicalmask) ) THEN
       npix = COUNT(logicalmask)
    ELSE
       npix = 12*nside*nside
    ENDIF
    CALL message(code, msg="Initial map has",i=npix,msg2="pixels")

    ALLOCATE(x(0:npix-1),stat = status)
    CALL assert_alloc(status,code,'x')
    x = PACK(map, mask=logicalmask)

    CALL message(code, msg="Generating masked angular map")
    ALLOCATE(ang_map(0:npix-1,0:npix-1),stat = status)
    CALL assert_alloc(status,code,'ang_map')
    CALL gen_angles_map(npix, ang_map, nside, nested=.FALSE., mask=logicalmask) 

    ALLOCATE(Pl(0:npix-1,0:npix-1,0:nlmax),stat = status)
    CALL assert_alloc(status,code,'Pl')
    DO l=0, nlmax
       CALL message(code, msg="Legendre matrices :", i=l, msg2="/", i2=nlmax)
       IF( l .EQ. 0 )THEN
          Pl(0:npix-1,0:npix-1,0) = 1.0_dp
       ELSE IF( l .EQ. 1 ) THEN
          Pl(0:npix-1,0:npix-1,1) = ang_map
       ELSE
          Pl(0:npix-1,0:npix-1,l) = &
               & ( (2.0_dp*REAL(l)-1.0_dp) * ang_map(0:npix-1,0:npix-1) * &
               & Pl(0:npix-1,0:npix-1,l-1) / ( (2.0_dp*REAL(l-1)+1.0_dp)/FOURPI ) &
               & - REAL(l-1) * Pl(0:npix-1,0:npix-1,l-2) / ((2.0_dp*REAL(l-2)+1.0_dp)/FOURPI) ) / REAL(l)
       END IF
       Pl(0:npix-1,0:npix-1,l) = Pl(0:npix-1,0:npix-1,l) * (2.0_dp*REAL(l)+1.0_dp)/FOURPI
    ENDDO
    DEALLOCATE( ang_map )

    ALLOCATE(F(0:nlmax_recon,0:nlmax_recon),stat = status)
    CALL assert_alloc(status,code,'F')
    ALLOCATE(y(0:nlmax_recon),stat = status)
    CALL assert_alloc(status,code,'y')

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(nlmax_recon, map, cl, Pl, npix, F, y, x, covarest) & 
    !$OMP PRIVATE(l, ll, Fll)
    DO l=0, nlmax_recon
       CALL message(code, msg="Power spectrum estimation :", i=l, msg2="/", i2=nlmax_recon)
       DO ll=0, nlmax_recon
          Fll = TRACEPROD(Pl(0:npix-1,0:npix-1,l),Pl(0:npix-1,0:npix-1,ll),0,npix-1)
          F(ll,l) = Fll
       ENDDO
       y(l) = DOT_PRODUCT( x, MATMUL( Pl(0:npix-1,0:npix-1,l), x ) )
    ENDDO
    !$OMP END PARALLEL DO

    CALL message(code, msg="Inverting Fisher matrix")
    CALL INVERSE(F(0:nlmax_recon,0:nlmax_recon),nlmax_recon+1,res)
    IF( res .EQV. .FALSE. )THEN
       PRINT*,"PROBLEM WITH INVERSION"
    ENDIF

    CALL message(code, msg="Evaluating final estimates")
    cl(0:nlmax_recon) = MATMUL( F(0:nlmax_recon,0:nlmax_recon), y(0:nlmax_recon) )    

    IF( present(covarest) ) THEN
       DO l=0, nlmax_recon
          covarest(l) = F(l,l)
       ENDDO
    ENDIF

    DEALLOCATE(F, Pl, x, y)

    CALL message(code, fin=.TRUE.)

  END SUBROUTINE map2cl_PCL


  ! ---------------------------------------------------------------------------------------!   

  !> Parallel bilinear form
  !!@param[in] ( x(s:e),y(s:e) ) : left and right vectors
  !!@param[in] A(s:e,s:e) : middle matrix
  !!@param[in] (s,e) : bounds
  REAL(DP) FUNCTION BILIN_FORM_PAR( x, A, y, s, e )     

    INTEGER(I4B) :: s, e, i, j, k
    REAL(DP), DIMENSION(s:e) :: x, y
    REAL(DP), DIMENSION(s:e,s:e) :: A
    REAL(DP) :: val

    val = 0.0_dp
    !$OMP PARALLEL DO REDUCTION(+:val) &
    !$OMP SHARED(x,A,y,s,e) & 
    !$OMP PRIVATE(i,j)	
    DO i=s,e
       val = val + x(i) * DOT_PRODUCT( A(i,s:e), y(s:e) )
    ENDDO
    !$OMP END PARALLEL DO

    BILIN_FORM_PAR = val

    RETURN

  END FUNCTION BILIN_FORM_PAR

  ! ---------------------------------------------------------------------------------------!   

  !> Parallel matrix multiplication
  !!@param[in] A(s:d1,s:d2) : first matrix
  !!@param[in] B(s:d2,s:d3) : second matrix
  !!@param[in] C(s:d1,s:d3) : output matrix
  !!@param[in] (s,d1,d2,d3) : bounds
  SUBROUTINE MATMUL_PAR(A, B, C, s, d1, d2, d3)

    INTEGER(I4B) :: s, d1, d2, d3, i, j, k
    REAL(DP), DIMENSION(s:d1,s:d2) :: A
    REAL(DP), DIMENSION(s:d2,s:d3) :: B
    REAL(DP), DIMENSION(s:d1,s:d3) :: C

    !$OMP PARALLEL DO &
    !$OMP SHARED(A,B,C,s,d1,d2,d3) &
    !$OMP PRIVATE(i,j)
    DO i=s,d1
       DO j=s,d3
          C(i,j) = DOT_PRODUCT( A(i,s:d2), B(s:d2,j) )
       ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    RETURN
  END SUBROUTINE MATMUL_PAR

  ! ---------------------------------------------------------------------------------------!   

  !> TRACE
  !!@param[in] mat(s:e,s:e) : input matrix
  !!@param[in] (s,e) : bounds
  REAL(DP) FUNCTION TRACE(mat, s, e)

    INTEGER(I4B) :: i, s, e
    REAL(DP), INTENT(IN),  DIMENSION(s:e,s:e) :: mat
    REAL(DP) :: val

    val = 0.0
    !$OMP PARALLEL DO REDUCTION(+:val) &
    !$OMP SHARED(mat,s,e) & 
    !$OMP PRIVATE(i)	
    DO i=s,e 
       val = val + mat(i,i)
    ENDDO
    !$OMP END PARALLEL DO

    TRACE = val

    RETURN
  END FUNCTION TRACE

  ! ---------------------------------------------------------------------------------------!   

  !> TRACE of a product of matrices
  !!@param[in] mat1(s:e,s:e) : input matrix 1
  !!@param[in] mat2(s:e,s:e) : input matrix 2
  !!@param[in] (s,e) : bounds
  REAL(DP) FUNCTION TRACEPROD(mat1, mat2, s, e)

    INTEGER(I4B) :: i, s, e
    REAL(DP), INTENT(IN),  DIMENSION(s:e,s:e) :: mat1, mat2
    REAL(DP) :: val

    val = 0.0
    !$OMP PARALLEL DO REDUCTION(+:val) &
    !$OMP SHARED(mat1,mat2,s,e) &
    !$OMP PRIVATE(i)	
    DO i=s,e 
       val = val + DOT_PRODUCT(mat1(i,s:e),mat2(s:e,i))
    ENDDO
    !$OMP END PARALLEL DO

    TRACEPROD = val

    RETURN
  END FUNCTION TRACEPROD

  ! ---------------------------------------------------------------------------------------!   
  ! ---------------------------------------------------------------------------------------!   

END MODULE f3dex_stats
