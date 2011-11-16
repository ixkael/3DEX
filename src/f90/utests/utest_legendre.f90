PROGRAM utest_legendre

  USE OMP_LIB

  USE healpix_types
  USE healpix_modules

  USE f3dex_transforms
  USE f3dex_cosmotools
  USE f3dex_stats
  USE f3dex_fitstools
  USE f3dex_utils

  IMPLICIT NONE

  INTEGER(I4B), PARAMETER :: nc=200
  INTEGER(I4B), PARAMETER :: prec=10
  CHARACTER(LEN=*), PARAMETER :: code = 'utest_legendre'
  INTEGER(I4B) :: i
  INTEGER(I4B), PARAMETER :: nlmax = 10
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

  CALL assert(res)

END PROGRAM utest_legendre
