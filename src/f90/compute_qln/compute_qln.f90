PROGRAM compute_qln

  !COMMAND :    

  USE f3dex_transforms
  USE f3dex_cosmotools
  USE f3dex_stats
  USE f3dex_fitstools
  USE f3dex_utils

  USE healpix_types
  USE healpix_modules

  CHARACTER(len=FILENAMELEN) ::  arg, qlnfile
  CHARACTER(len=80), DIMENSION(1:120) :: header
  INTEGER(i4b):: rdstatus, status, n_args
  REAL(DP), DIMENSION(:,:), allocatable :: qln
  INTEGER(I4B)  :: nnmax, nlmax
  CHARACTER(len=*), PARAMETER :: code = "gen_qln"
  LOGICAL(LGT) ::  fileexists
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: tod
  
  PRINT*," "
  PRINT*,"================================================"
  PRINT*,"COMPUTE_QLN"
  PRINT*," "

  n_args = nArguments()

  CALL getArgument(1, arg)
  qlnfile = trim(adjustl(arg))
  INQUIRE( file=qlnfile, exist=fileexists)
  IF( fileexists .eqv. .TRUE. ) CALL fatal_error("qlnfile : this file already exists!")

  CALL getArgument(2, arg)
  READ(arg,'(I4)') nnmax

  CALL getArgument(3, arg)
  READ(arg,'(I4)') nlmax

  ALLOCATE(qln(0:nlmax,1:nnmax),stat = status)
  CALL assert_alloc(status,code,'qln')

  PRINT*,"Generating qln coefficients"
  CALL gen_qln(qln, nnmax, nlmax)
  PRINT*,"Done."
  
  PRINT*," "
  PRINT*,"Exporting coefficients"
  header(:) = ' '
  CALL bitab2fits( qlnfile, qln, 0, nlmax, 1, nnmax )
  PRINT*,"Done."
 
  PRINT*," "
  PRINT*,"COMPUTE_QLN : terminated."
  PRINT*,"================================================"
  PRINT*," "

END PROGRAM compute_qln
