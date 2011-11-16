PROGRAM utest_trace

  USE OMP_LIB

  USE healpix_types
  USE healpix_modules

  USE f3dex_transforms
  USE f3dex_cosmotools
  USE f3dex_stats
  USE f3dex_fitstools
  USE f3dex_utils

  IMPLICIT NONE

  INTEGER(I4B), PARAMETER :: ndim=200
  INTEGER(I4B), PARAMETER :: prec=10
  INTEGER(I4B) :: i
  REAL(DP),DIMENSION(0:ndim-1,0:ndim-1) :: M
  CHARACTER(LEN=*), PARAMETER :: code = 'utest_trace'
  REAL(DP) :: res
  LOGICAL :: res2

  CALL message(code,msg="Computes the trace of a large matrix")
  M = 0.0
  DO i=0,ndim-1
     M(i,i) = 1.0
  ENDDO
  res = TRACE(M,0,ndim-1)
  res2 = assert_DP( res, REAL(ndim,DP), prec )

  CALL assert(res2)

END PROGRAM utest_trace
