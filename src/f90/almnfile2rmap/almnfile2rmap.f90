PROGRAM almnfile2rmap

  !COMMAND : almnfile2map almninfile mapoutfile rho nside nnmax nlmax rmax
  !TEST : bin/almnfile2rmap out/almn.fits map.fits 400.0 256 10 10 2000.0

  USE f3dex_transforms
  USE f3dex_cosmotools
  USE f3dex_stats
  USE f3dex_fitstools
  USE f3dex_utils

  USE healpix_types
  USE healpix_modules
  USE misc_utils
  USE fitstools
  USE alm_tools

  REAL(DP) :: rmax              
  REAL*8 :: time1, time2
  REAL(DP) :: rho
  INTEGER(I4B) :: n_args, nnmax, nlmax, nmmax, nside, status
  LOGICAL(LGT) :: renormalization, fileexists
  CHARACTER(len=80), DIMENSION(1:120):: header
  CHARACTER(len=FILENAMELEN) :: headersfile, arg, almninfile, mapoutfile
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: map, kln, cln, plm 
  COMPLEX(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: almn
  CHARACTER(len=*), PARAMETER :: code = "almn2rmap"

  PRINT*,"=================================================="
  PRINT*," 3DEX F90 LIBRARY"
  PRINT*," STARTING FIELD RECONSTRUCTION"
  PRINT*,"=================================================="
  PRINT*," "

  !---------------------------------------------------------------------------
  ! Parameters extraction

  n_args = nArguments()

  CALL getArgument(1, arg)
  almninfile = trim(adjustl(arg))
  INQUIRE( file=almninfile, exist=fileexists)
  IF( fileexists .eqv. .FALSE. ) CALL fatal_error("almnfile : this file doesn't exist!")

  CALL getArgument(2, arg)
  mapoutfile = trim(adjustl(arg))
  INQUIRE( file=mapoutfile, exist=fileexists)
  IF( fileexists .eqv. .TRUE. ) CALL fatal_error("mapfile : this file already exists!")

  CALL getArgument(3, arg)
  READ(arg,'(F7.9)') rho
  PRINT*,"rho :",rho

  CALL getArgument(4, arg)
  READ(arg,'(I4)') nside

  CALL getArgument(5, arg)
  READ(arg,'(I4)') nnmax

  CALL getArgument(6, arg)
  READ(arg,'(I4)') nlmax

  IF( n_args == 7 ) then
     CALL getArgument(7, arg)
     READ(arg,'(F7.9)') rmax
     renormalization = .TRUE.
     PRINT*,rmax
  ELSE 
     rmax = 0.0
     renormalization = .FALSE.
  ENDIF

  !---------------------------------------------------------------------------
  ! Preliminaries

  nmmax = nlmax

  time1 = wtime()

  !Allocating memory
  ALLOCATE(almn(1:nnmax,0:nlmax,0:nmmax), stat = status )
  call assert_alloc(status,code,"almn")
  ALLOCATE(kln(0:nlmax,1:nnmax),stat = status)
  CALL assert_alloc(status,code,"kln")
  ALLOCATE(cln(0:nlmax,1:nnmax),stat = status)
  CALL assert_alloc(status,code,"cln")

  !Extraction from fits file
  PRINT*,"Starting extraction..."
  CALL fits2almn(almninfile, nnmax, nlmax, nmmax, almn, kln, cln, header,120)
  PRINT*,"Done"

  !If necessary, print the extracted data
  !CALL print_spectrum(cln, nlmax, nnmax, nlmax, nnmax," cln")
  !CALL print_spectrum(kln,nlmax, nnmax, nlmax, nnmax," kln")
  !CALL print_almn(almn, 3, 3, 3, nlmax, nmmax, nnmax,"almn")

  !Change the cln's normalization?
  IF( renormalization .eqv. .TRUE. ) THEN
     CALL gen_cln(cln, kln, nnmax, nlmax, rmax)
  ELSE
     cln = 1.0_dp
  ENDIF

  !Allocating data and harmonics
  ALLOCATE( map(0:(12*nside**2),1:1) )
  map = 0.0_dp
  ALLOCATE( plm(0:(nside*(nmmax+1)*(2*nlmax-nmmax+2)),1:1) )
  CALL plm_gen( nside, nlmax, nmmax, plm )

  !---------------------------------------------------------------------------
  ! Map reconstruction

  !Computing reconstructed map
  PRINT*,"Starting reconstruction..."
  CALL almn2rmap(map, almn, rho, nside, nnmax, nlmax, nmmax, kln, cln, plm)
  PRINT*,"Done"

  time2 = wtime()

  !---------------------------------------------------------------------------
  ! End of the procedure

  PRINT*,"=================================================="
  PRINT*," "
  PRINT*,"RECONSTRUCTION (almnfil2rmap) - SUMMARY"
  PRINT*," "
  PRINT*,"rho / nside", rho, nside
  PRINT*,"(nmax,lmax,mmax) :", nnmax, nlmax, nmmax
  PRINT*,"Elapsed time :",time2-time1
  PRINT*," "
  PRINT*,"=================================================="
  PRINT*," "
  !Fits output
  PRINT*,"Starting map output..."
  CALL write_minimal_header(header, 'MAP', nside=nside, ordering='Ring')
  CALL write_bintab( map, (12*nside**2), 1, header, 120, mapoutfile )
  PRINT*,"Done"
  PRINT*," "
  PRINT*,"=================================================="

END PROGRAM almnfile2rmap
