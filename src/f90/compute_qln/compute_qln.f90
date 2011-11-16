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
  CHARACTER(len=*), PARAMETER :: code = "compute_qln"
  LOGICAL(LGT) ::  fileexists
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: tod
  
  CALL message(code,start=.TRUE.)

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

  CALL message(code,msg="Generating qln coefficients")
  CALL gen_qln(qln, nnmax, nlmax)
  
  CALL message(code,msg="Exporting coefficients")
  header(:) = ' '
  CALL bitab2fits( qlnfile, qln, 0, nlmax, 1, nnmax )
 
  CALL message(code,fin=.TRUE.)

END PROGRAM compute_qln
