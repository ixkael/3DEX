PROGRAM utests 

  USE f3dex_transforms
  USE f3dex_cosmotools
  USE f3dex_stats
  USE f3dex_fitstools
  USE f3dex_utils
  USE f3dex_utests

  USE healpix_types
  USE healpix_modules

  ! ---------------------------------------------------------------------------------------!     

  LOGICAL :: res
  INTEGER(I4B), PARAMETER :: prec = 10
  INTEGER(I4B) :: tot, succ
  PRINT*," "
  
  tot=0
  succ=0
  ! ---------------------------------------------------------------------------------------!     
  PRINT*," "
  PRINT*,"F3DEX - TRANSFORM MODULE ======================="
  PRINT*,"Test : legendre functions & recurrence formula"
  PRINT*,"------------------------------------------------"
  res = utest_legendre(50,prec)
  PRINT*,"> Tested 50 occurences in cos[ 0 , 2pi ]"
  PRINT*,"> Result compared with analytical formula"
  PRINT*,"> Orders tested: l in [ 0 , 10 ]"
  PRINT*,"> Precision :",prec,"digits"
  tot=tot+1
  IF( res .EQV. .TRUE.) THEN
     succ=succ+1
     PRINT*,"-T---------- SUCCESS -------------------------T-"
  ELSE
     PRINT*,"-F---------- FAILURE ------------------------<F>"
  ENDIF
  PRINT*," "
  ! ---------------------------------------------------------------------------------------!     


  ! ---------------------------------------------------------------------------------------!     
  PRINT*," "
  PRINT*,"F3DEX - TRANSFORM MODULE ======================="
  PRINT*,"Test : matrix recurrence for legendre functions"
  PRINT*,"------------------------------------------------"
  res = utest_legendre_angmatrices(4,10,prec)
  PRINT*,"> Tested map with nside = 8"
  PRINT*,"> Result compared with raw recurrence"
  PRINT*,"> Orders tested: l in [ 0 , 10 ]"
  PRINT*,"> Precision :",prec,"digits"
  tot=tot+1
  IF( res .EQV. .TRUE.) THEN
     succ=succ+1
     PRINT*,"-T---------- SUCCESS -------------------------T-"
  ELSE
     PRINT*,"-F---------- FAILURE ------------------------<F>"
  ENDIF
  PRINT*," "
  ! ---------------------------------------------------------------------------------------!     

  	  
  ! ---------------------------------------------------------------------------------------!     
  PRINT*," "
  PRINT*,"F3DEX - STAT MODULE ============================"
  PRINT*,"Test : TRACE operator"
  PRINT*,"------------------------------------------------"
  res = utest_trace(50,prec)
  PRINT*,"> Tested matrix of dimension 50"
  PRINT*,"> Precision :",prec,"digits"
  tot=tot+1
  IF( res .EQV. .TRUE.) THEN
     succ=succ+1
     PRINT*,"-T---------- SUCCESS -------------------------T-"
  ELSE
     PRINT*,"-F---------- FAILURE ------------------------<F>"
  ENDIF
  PRINT*," "
  ! ---------------------------------------------------------------------------------------!   

  ! ---------------------------------------------------------------------------------------!     
  PRINT*," "
  PRINT*,"F3DEX - STAT MODULE ============================"
  PRINT*,"Test : TRACE PROD operator"
  PRINT*,"------------------------------------------------"
  res = utest_traceprod(50,prec)
  PRINT*,"> Tested two matrices of dimension 50"
  PRINT*,"> Precision :",prec,"digits"
  tot=tot+1
  IF( res .EQV. .TRUE.) THEN
     succ=succ+1
     PRINT*,"-T---------- SUCCESS -------------------------T-"
  ELSE
     PRINT*,"-F---------- FAILURE ------------------------<F>"
  ENDIF
  PRINT*," "
  ! ---------------------------------------------------------------------------------------!   
  	  
  ! ---------------------------------------------------------------------------------------!     
  !PRINT*," "
  !PRINT*,"F3DEX - STAT MODULE ============================"
  !PRINT*,"Test : Quadratic QML estimator"
  !PRINT*,"------------------------------------------------"
  !res = utest_map2cl_QML( 4, 4 )
  !tot=tot+1
  !IF( res .EQV. .TRUE.) THEN
  !   succ=succ+1
  !   PRINT*,"-T---------- SUCCESS -------------------------T-"
  !ELSE
  !   PRINT*,"-F---------- FAILURE ------------------------<F>"
  !ENDIF
  !PRINT*," "
  ! ---------------------------------------------------------------------------------------!   

  ! ---------------------------------------------------------------------------------------!     
  PRINT*," "
  PRINT*,"F3DEX - STAT MODULE ============================"
  PRINT*,"Test : Matrix inversion using LAPACK"
  PRINT*,"------------------------------------------------"
  PRINT*,"> Testing matrix of dimension 50"
  PRINT*,"> Precision :",prec,"digits"
  res = utest_inverse( 50, prec )
  tot=tot+1
  IF( res .EQV. .TRUE.) THEN
     succ=succ+1
     PRINT*,"-T---------- SUCCESS -------------------------T-"
  ELSE
     PRINT*,"-F---------- FAILURE ------------------------<F>"
  ENDIF
  PRINT*," "
  ! ---------------------------------------------------------------------------------------!   
  
  ! ---------------------------------------------------------------------------------------!     
  PRINT*," "
  PRINT*,"F3DEX - STAT MODULE ============================"
  PRINT*,"Test : Parallel bilinear form xAy"
  PRINT*,"------------------------------------------------"
  PRINT*,"> Testing random matrix/vectors of dimension 1000"
  PRINT*,"> Precision :",prec,"digits"
  res = utest_bilini_form_par( 1000, prec )
  tot=tot+1
  IF( res .EQV. .TRUE.) THEN
     succ=succ+1
     PRINT*,"-T---------- SUCCESS -------------------------T-"
  ELSE
     PRINT*,"-F---------- FAILURE ------------------------<F>"
  ENDIF
  PRINT*," "
  ! ---------------------------------------------------------------------------------------!   
  
  CALL assert(.FALSE.)
  ! ---------------------------------------------------------------------------------------!     
  PRINT*," "
  PRINT*,"F3DEX - STAT MODULE ============================"
  PRINT*,"Test : Parallel matrix multiplication"
  PRINT*,"------------------------------------------------"
  PRINT*,"> Testing random matrix/vectors of dimension 100"
  PRINT*,"> Precision :",prec,"digits"
  res = utest_matmul_par( 0, 100, 100, 100, prec )
  tot=tot+1
  IF( res .EQV. .TRUE.) THEN
     succ=succ+1
     PRINT*,"-T---------- SUCCESS -------------------------T-"
  ELSE
     PRINT*,"-F---------- FAILURE ------------------------<F>"
  ENDIF
  PRINT*," "
  ! ---------------------------------------------------------------------------------------!   
  
  PRINT*," "    
  PRINT*," "
  PRINT*,"F3DEX - UTESTS ============================"
  PRINT*,"PASSED TESTS :",succ,"ON",tot
  PRINT*," "    
  PRINT*," "
END PROGRAM utests
