PROGRAM compute_legen_plm

  !COMMAND :    

  USE f3dex_transforms
  USE f3dex_cosmotools
  USE f3dex_stats
  USE f3dex_fitstools
  USE f3dex_utils

  USE healpix_types
  USE healpix_modules

  CHARACTER(len=FILENAMELEN) ::  arg, plmfile
  CHARACTER(len=80), DIMENSION(1:120) :: header
  INTEGER(i4b):: rdstatus, status, n_args
  REAL(DP), DIMENSION(:,:), allocatable :: plm, plm2
  INTEGER(I4B)  :: nlmax, nsmax
  CHARACTER(len=*), PARAMETER :: code = "compute_legen_plm"
  LOGICAL(LGT) ::  fileexists
  LOGICAL :: anynull
  REAL(DP) :: nullval

  CALL message(code,start=.TRUE.)
  
  n_args = nArguments()

  CALL getArgument(1, arg)
  plmfile = trim(adjustl(arg))
  INQUIRE( file=plmfile, exist=fileexists)
  IF( fileexists .eqv. .TRUE. ) CALL fatal_error("plmfile : this file already exists!")

  CALL getArgument(2, arg)
  READ(arg,'(I4)') nsmax
  
  CALL getArgument(3, arg)
  READ(arg,'(I4)') nlmax

  n_plm = nsmax*(nlmax+1)*(2*nlmax-nlmax+2)
  
  ALLOCATE(plm(0:n_plm-1,1:1))
  CALL assert_alloc(status,code,'plm')
  ALLOCATE(plm2(0:n_plm-1,1:1))
  CALL assert_alloc(status,code,'plm2')

  CALL message(code,msg="Generating plm coefficients")
  CALL plm_gen(nsmax, nlmax, nlmax, plm)

  CALL message(code,msg="Exporting coefficients")
  header(:) = ' '
  CALL write_plm(plm, n_plm, 1, header, 120, plmfile, nsmax, nlmax)
  
  CALL message(code,fin=.TRUE.)

END PROGRAM compute_legen_plm
