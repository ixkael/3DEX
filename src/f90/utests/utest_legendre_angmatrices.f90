PROGRAM utest_legendre_angmatrices

  USE OMP_LIB

  USE healpix_types
  USE healpix_modules

  USE f3dex_transforms
  USE f3dex_cosmotools
  USE f3dex_stats
  USE f3dex_fitstools
  USE f3dex_utils

  IMPLICIT NONE
  
     INTEGER(I4B), PARAMETER :: nside=16
     INTEGER(I4B), PARAMETER :: nlmax=10
     INTEGER(I4B), PARAMETER :: prec=10
   INTEGER(I4B) :: l, p, npix, status
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
     
     CALL assert(res)
  
END PROGRAM utest_legendre_angmatrices
