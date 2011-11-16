MODULE f3dex_utests

  USE OMP_LIB

  USE healpix_types
  USE healpix_modules

  USE f3dex_transforms
  USE f3dex_cosmotools
  USE f3dex_stats
  USE f3dex_fitstools
  USE f3dex_utils

  IMPLICIT NONE

CONTAINS
    
  ! ---------------------------------------------------------------------------------------!  
   
  !> Test inverse functionalty
  !!@param[in] n : dimension
  !!@param[in] prec : precision digits
  LOGICAL FUNCTION utest_inverse(n, prec)
  
     INTEGER(I4B) :: n, prec, i
     REAL(DP), DIMENSION(0:n,0:n) :: Aini, A, Ainv, Id
     LOGICAL :: res1, res2, res3, res4
     A = 0.0

     A(0,0) = 3.0
     Id(0,0) = 1.0
     DO i=1, n-1
     	A(i,i) = 2.0
     	A(i+1,i-1) = 1.0
     	A(i-1,i+1) = 4.0
     	Id(i,i) = 1.0
     ENDDO
     A(n,n) = 3.0
     Id(n,n) = 1.0
     
     Aini = A
     ! Compute inverse
     CALL INVERSE(A(0:n,0:n),n+1,res1)     
     PRINT*,"> Compute inverse :", res1
     
     Ainv = A
     ! Back to A
     CALL INVERSE(A(0:n,0:n),n+1,res2)
     PRINT*,"> Back to original :", res1
     
     ! Assert back to original
     res3 = assert_DPARR(Aini,A,n+1,prec)
     PRINT*,"> Is valid original :", res3
     
     ! Product should be identity
     A = MATMUL(A,Ainv)
     res4 = assert_DPARR(Id,A,n+1,prec)
     PRINT*,"> A*Ainv = Id(n) :", res4
     
     utest_inverse = res1 .AND. res2 .AND. res3 .AND. res4
     
     RETURN
     
  END FUNCTION utest_inverse
  
  ! ---------------------------------------------------------------------------------------!  
   
  !> Test trave functionalty
  !!@param[in] ndim : dimension
  !!@param[in] prec : precision digits
  LOGICAL FUNCTION utest_trace(ndim, prec)

    INTEGER(I4B) :: ndim, i, prec
    REAL(DP),DIMENSION(0:ndim-1,0:ndim-1) :: M
    REAL(DP) :: res

    M = 0.0
    DO i=0,ndim-1
       M(i,i) = 1.0
    ENDDO
    res = TRACE(M,0,ndim-1)
    utest_trace = assert_DP( res, REAL(ndim,DP), prec )

  END FUNCTION utest_trace
  
  ! ---------------------------------------------------------------------------------------!  
   
  !> Test trace_prod functionalty
  !!@param[in] ndim : dimension
  !!@param[in] prec : precision digits
  LOGICAL FUNCTION utest_traceprod(ndim, prec)

    INTEGER(I4B) :: ndim, i, prec
    REAL(DP),DIMENSION(0:ndim-1,0:ndim-1) :: M, Mb
    REAL(DP) :: res

    M = 1.0
    Mb = 1.0
    DO i=0,ndim-1
       M(i,i) = 1.0
    ENDDO
    res = TRACEPROD(M,Mb,0,ndim-1)
    utest_traceprod = assert_DP( res, REAL(ndim**2,DP), prec )

  END FUNCTION utest_traceprod

  ! ---------------------------------------------------------------------------------------!  
   
  !> Test legendre functionalty
  !!@param[in] nc : dimension
  !!@param[in] prec : precision digits
  LOGICAL FUNCTION utest_legendre(nc, prec)

    IMPLICIT NONE

    INTEGER(I4B) :: i, prec
    INTEGER(I4B), PARAMETER :: nlmax = 10
    INTEGER(I4B) :: nc
    LOGICAL :: tempres, res
    REAL(DP), DIMENSION(0:nc-1) :: x
    REAL(DP), DIMENSION(0:nc-1,0:nlmax) :: pl_mat
    REAL(DP), DIMENSION(0:nlmax) :: pl_vec, solution

    DO i=0, nc-1
       x(i) = COS(REAL(i)*2*PI/nc)
    ENDDO
    
    res = .TRUE.

    CALL legendre_mat(x,pl_mat,nc,nlmax)
    DO i=0, nc-1
       solution(0) = 1.0
       solution(1) = x(i)
       solution(2) = 0.5 * ( 3 * x(i)**2 - 1 )
       solution(3) = 0.5 * ( 5 * x(i)**3 - 3*x(i) )
       solution(4) = (1.0/8.0) * ( 35.0 * x(i)**4 - 30.0* x(i)**2 + 3.0 )
       solution(5) = (1.0/8.0) * ( 63.0 * x(i)**5 - 70.0* x(i)**3 + 15.0*x(i) )
       solution(6) = (1.0/16.0) * ( 231.0 * x(i)**6 - 315.0* x(i)**4 + 105.0*x(i)**2 - 5.0 )
       solution(7) = (1.0/16.0) * ( 429.0 * x(i)**7 - 693.0* x(i)**5 + 315.0*x(i)**3 - 35.0*x(i) )
       solution(8) = (1.0/128.0) * ( 6435.0 * x(i)**8 - 12012.0* x(i)**6 + 6930.0*x(i)**4 &
                   & - 1260.0*x(i)**2 + 35.0 )
       solution(9) = (1.0/128.0) * ( 12155.0 * x(i)**9 - 25740.0* x(i)**7 + 18018.0*x(i)**5 &
                   & - 4620.0*x(i)**3 + 315.0*x(i) )
       solution(10) = (1.0/256.0) * ( 46189.0 * x(i)**10 - 109395.0* x(i)**8 + 90090.0*x(i)**6 &
                   & - 30030.0*x(i)**4 + 3465.0*x(i)**2 - 63.0 )
       CALL legendre_vec(x(i),pl_vec(0:nlmax),nlmax)
       tempres = assert_DPARR( pl_vec(0:nlmax), pl_mat(i,0:nlmax), nlmax, prec )
       res = res .AND. tempres
       tempres = assert_DPARR( pl_vec(0:nlmax), solution(0:nlmax), nlmax, prec )
       res = res .AND. tempres
    ENDDO
    utest_legendre = res

  END FUNCTION utest_legendre

  ! ---------------------------------------------------------------------------------------!     

  !> Test legendre angmatrices functionalty
  !!@param[in] nside : healpix dimension
  !!@param[in] nlmax : nlmax parameter
  !!@param[in] prec : precision digits
  LOGICAL FUNCTION utest_legendre_angmatrices(nside, nlmax, prec)
  
  INTEGER(I4B) :: l, nlmax, nside, prec, p, npix, status
  CHARACTER(LEN=*), PARAMETER :: code = 'utest_legendre_angmatrices'
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Pl1, Pl2, Pltemp, ang_map
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Plcomp 
  LOGICAL :: res
  
  npix = 12*nside**2
  
  ALLOCATE(ang_map(0:npix-1,0:npix-1),stat = status)
  CALL assert_alloc(status,code,'ang_map')
  ALLOCATE(Pl1(0:npix-1,0:npix-1),stat = status)
  CALL assert_alloc(status,code,'Pl1')
  ALLOCATE(Pl2(0:npix-1,0:npix-1),stat = status)
  CALL assert_alloc(status,code,'Pl2')
  ALLOCATE(Pltemp(0:npix-1,0:npix-1),stat = status)
  CALL assert_alloc(status,code,'Pltemp')
  ALLOCATE(Plcomp(0:npix-1,0:npix-1,0:nlmax),stat = status)
  CALL assert_alloc(status,code,'Plcomp')
    
  CALL gen_angles_map(npix, ang_map, nside)
  	  
  DO p=0, npix-1
     CALL legendre_mat(ang_map(p,0:npix-1),Plcomp(p,0:npix-1,0:nlmax),npix,nlmax)
  ENDDO
  
  res=.TRUE.
  DO l=0, nlmax
     IF( l .EQ. 0 )THEN
        Pl2 = 1.0
     ELSE IF( l .EQ. 1 ) THEN
        Pl1 = 1.0
        Pl2 = ang_map
     ELSE
       	PLtemp = Pl2
        Pl2 = ( (2.0*REAL(l)-1.0) * ang_map*Pl2 - REAL(l-1) * Pl1 ) / REAL(l)
       	Pl1 = Pltemp
    END IF
    res = res .AND. assert_DPARR( Plcomp(0:npix-1,0:npix-1,l),Pl2(0:npix-1,0:npix-1), nlmax, prec)
 ENDDO
  
  DEALLOCATE(Pl1, Pl2, Pltemp, ang_map)
  utest_legendre_angmatrices = res
  RETURN
  
  END FUNCTION utest_legendre_angmatrices
	  
  ! ---------------------------------------------------------------------------------------!     
  
  !> Test bilinear form functionalty
  !!@param[in] ndim : dimension
  !!@param[in] prec : precision digits
  LOGICAL FUNCTION utest_bilini_form_par( ndim, prec )
  
    INTEGER(I4B) :: ndim, i, j, prec
    REAL(DP),DIMENSION(0:ndim-1,0:ndim-1) :: M
    REAL(DP),DIMENSION(0:ndim-1) :: x, y
    REAL(DP) :: res1, res2
    INTEGER,PARAMETER :: seed = 86456
    REAL*8 :: tas, tae, tbs, tbe
          
    CALL srand(seed)

    DO i=0,ndim-1
       x(i) = rand()
       y(i) = rand()
       DO j=0, ndim-1
          M(i,j) = rand()
       ENDDO
    ENDDO
    
    tas = wtime()
    res1 = BILIN_FORM_PAR( x(0:ndim-1), M(0:ndim-1,0:ndim-1), &
    			 & y(0:ndim-1), 0, ndim-1 )
    tae = wtime()
    
    tbs = wtime()
    res2 = DOT_PRODUCT( x, MATMUL( M, y ) )
    tbe = wtime()
    
    !PRINT*,"NATIVE METHOD   :",tbe-tbs,"ms"
    !PRINT*,"PARALLEL METHOD :",tae-tas,"ms"
    
    IF( assert_DP( res1, res2 , prec ) .EQV. .FALSE.)THEN
    	PRINT*,"ERROR : ", res1,"!=", res2   
    ENDIF
    	   
    utest_bilini_form_par = assert_DP( res1, res2, prec )
    RETURN
  
  END FUNCTION utest_bilini_form_par
      
  ! ---------------------------------------------------------------------------------------!     
  
  !> Test parallel matrix multiplication functionalty
  !!@param[in] (s,d1,d2,d3) : start and three end dimensions
  !!@param[in] prec : precision digits
  LOGICAL FUNCTION utest_matmul_par( s, d1, d2, d3, prec )
  
    INTEGER(I4B) :: s, d1, d2, d3, i, j, k, prec
    REAL(DP), DIMENSION(s:d1,s:d2) :: A
    REAL(DP), DIMENSION(s:d2,s:d3) :: B
    REAL(DP), DIMENSION(s:d1,s:d3) :: C, D
    INTEGER,PARAMETER :: seed = 86456
    LOGICAL :: res   
    REAL*8 :: tas, tae, tbs, tbe
    
    CALL srand(seed)

    DO i=s,d1
       DO j=s, d2
       	  DO k= s, d3
             A(i,j) = rand()
             B(j,k) = rand()
          ENDDO
       ENDDO
    ENDDO
    
    tas = wtime()
    CALL MATMUL_PAR( A(s:d1,s:d2), B(s:d2,s:d3), C(s:d1,s:d3), &
    		   & s, d1, d2, d3 )
    tae = wtime()
    
    tbs = wtime()
    D = MATMUL( A, B )	
    tbe = wtime()
    
    DO i=s, d1
    	DO k=s,d3
    	   IF( assert_DP( C(i,k) , D(i,k) , prec ) .EQV. .FALSE.)THEN
    	      PRINT*,"ERROR : ", C(i,k), "!=" , D(i,k)	    
    	   ENDIF
    	   res = res .AND.  assert_DP( C(i,k) , D(i,k) , prec )
    	ENDDO
    ENDDO
    
    !PRINT*,"> NATIVE METHOD   :",tbe-tbs,"s"
    !PRINT*,"> PARALLEL METHOD :",tae-tas,"s"
    
    utest_matmul_par = res
    RETURN
  
  END FUNCTION utest_matmul_par
      
  ! ---------------------------------------------------------------------------------------!     
  
  	  
END MODULE f3dex_utests
