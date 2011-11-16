PROGRAM utest_inverse

  USE OMP_LIB

  USE healpix_types
  USE healpix_modules

  USE f3dex_transforms
  USE f3dex_cosmotools
  USE f3dex_stats
  USE f3dex_fitstools
  USE f3dex_utils

  IMPLICIT NONE

  INTEGER(I4B), PARAMETER :: n=200
  INTEGER(I4B), PARAMETER :: prec=10
  INTEGER(I4B) :: i
  REAL(DP), DIMENSION(0:n,0:n) :: Aini, A, Ainv, Id
  LOGICAL :: res1, res2, res3, res4, res
  CHARACTER(LEN=*), PARAMETER :: code = 'utest_inverse'
  A = 0.0

  CALL message(code,msg="Generates large matrix")

  A(0,0) = 3.0
  Id(0,0) = 1.0
  DO i=1, n-1
     A(i,i) = 2.0
     A(i+1,i-1) = 3.0
     A(i-1,i+1) = 4.0
     Id(i,i) = 1.0
  ENDDO
  A(n,n) = 3.0
  Id(n,n) = 1.0

  Aini = A
  CALL message(code,msg="Computes the inverse")
  CALL INVERSE(A(0:n,0:n),n+1,res1)     
  PRINT*,"> Compute inverse :", res1

  Ainv = A
  CALL message(code,msg="Inverts it again")
  CALL INVERSE(A(0:n,0:n),n+1,res2)
  PRINT*,"> Back to original :", res2

  CALL message(code,msg="Assert back to original")
  !res3 = assert_DPARR(Aini,A,n+1,prec)
  !PRINT*,"> Is valid original :", res3

  CALL message(code,msg="Product should be identity")
  A = MATMUL(A,Ainv)
  res4 = assert_DPARR(Id,A,n+1,prec)
  PRINT*,"> A*Ainv = Id(n) :", res4

  res = res1 .AND. res2 .AND. res4

  CALL assert(res)

END PROGRAM utest_inverse
