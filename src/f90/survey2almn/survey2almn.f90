PROGRAM survey2almn

  !COMMAND : bin/survey2almn test/horizon_02_30k.txt out/almn_02_50_1024.fits out/cln_02_50_1024.fits 50 50 1024 2000.0
  !TEST : bin/survey2almn test/horizon_01_30k.txt out/almn.fits out/cln.fits 50 50 256 2000.0
  ! export OMP_NUM_THREADS=8    

  USE f3dex_transforms
  USE f3dex_cosmotools
  USE f3dex_stats
  USE f3dex_fitstools
  USE f3dex_utils

  USE healpix_types
  USE healpix_modules

  CHARACTER(len=FILENAMELEN) ::  arg, clnfile, almnfile, inputfile
  CHARACTER(len=80), DIMENSION(1:120) :: header
  INTEGER(i4b):: rdstatus, i, status, n_args, ntot
  REAL(DP), DIMENSION(:,:), allocatable :: kln, cln
  REAL(DP), DIMENSION(:,:), allocatable :: survey, surveytemp
  INTEGER(I4B)  :: nnmax, nlmax, nmmax, nside, nbpts, nlimit
  REAL(kind=DP),     DIMENSION(1:2)  :: zbounds
  CHARACTER(len=*), PARAMETER :: code = "survey2almn"
  COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: almn, almntemp
  REAL(kind=DP), DIMENSION(:,:), ALLOCATABLE :: spectr
  LOGICAL(LGT) ::  fileexists
  REAL(DP) :: rmax
  REAL*8 :: time1, time2

  print*," "
  print*,"================================================"
  PRINT*," 3DEX F90 LIBRARY"
  PRINT*," STARTING FOURIER-BESSEL DECOMPOSITION"
  PRINT*,"================================================"
  PRINT*," "

  !---------------------------------------------------------------------------
  ! Parameters extraction

  n_args = nArguments()

  CALL getArgument(1, arg)
  inputfile = trim(adjustl(arg))
  INQUIRE( file=inputfile, exist=fileexists)
  IF( fileexists .eqv. .FALSE. ) CALL fatal_error("inputfile : this file doesn't exist!")

  CALL getArgument(2, arg)
  almnfile = trim(adjustl(arg))
  INQUIRE( file=almnfile, exist=fileexists)
  IF( fileexists .eqv. .TRUE. ) CALL fatal_error("almnfile : this file already exists!")

  CALL getArgument(3, arg)
  clnfile = trim(adjustl(arg))
  INQUIRE( file=clnfile, exist=fileexists)
  IF( fileexists .eqv. .TRUE. ) CALL fatal_error("clnfile : this file already exists!")

  CALL getArgument(4, arg)
  READ(arg,'(I4)') nnmax

  CALL getArgument(5, arg)
  READ(arg,'(I4)') nlmax

  CALL getArgument(6, arg)
  READ(arg,'(I4)') nside

  CALL getArgument(7, arg)
  READ(arg,'(F7.9)') rmax

  !---------------------------------------------------------------------------
  ! Preliminaries 

  !CALL OMP_SET_NUM_THREADS(16)

  nlimit = 10000000
  PRINT*,"rmax :",rmax
  ! Additional parameters
  iter_order = 0
  nmmax = nlmax
  zbounds = (/ -1.0 , 1.0 /)

  time1 = wtime()

  !Computing kln's
  PRINT*," "
  ALLOCATE(kln(0:nlmax,1:nnmax),stat = status)
  !DO i = 0, nlmax
  !	call logrange( kln(i:i,:) , 0.001_DP , 0.5_DP , nnmax )
  !ENDDO
  CALL gen_kln(kln, nnmax, nlmax, rmax)

  !Computing cln's
  PRINT*," "
  ALLOCATE(cln(0:nlmax,1:nnmax),stat = status)
  CALL assert_alloc(status,code,"cln")
  CALL gen_cln(cln, kln, nnmax, nlmax, rmax)
  PRINT*," "	

  !Allocating memory
  ALLOCATE( almn(1:nnmax,0:nlmax,0:nmmax), stat = status )
  CALL assert_alloc(status,code,"almn")
  ALLOCATE( almntemp(1:nnmax,0:nlmax,0:nmmax), stat = status )
  CALL assert_alloc(status,code,"almntemp")
  almn = cmplx(0.0,0.0,kind=dp)

  !---------------------------------------------------------------------------
  ! Fourier-Bessel decomposition (the survey is cut if necessary)

  !Open survey file
  OPEN(2,file=inputfile,iostat=rdstatus)
  !PRINT*,"Opening file status: ",rdstatus

  !Read the whole survey 
  DO WHILE(rdstatus == 0)  
     !Using a temporary array, limited to nlimit points
     ALLOCATE( surveytemp(1:nlimit,1:3) )
     CALL assert_alloc(status,code,"surveytemp")

     nbpts = 0
     surveytemp = 0.0_dp
     !Read the points
     DO i=1,nlimit
        READ(2,"(1X,F18.14,1X,F18.14,1X,F18.12)",iostat = rdstatus,end=61) (surveytemp(i,k),k=1,3) 
        nbpts = nbpts + 1
     ENDDO
61   CONTINUE  

     !Allocating and copying definitive array
     ALLOCATE( survey(1:nbpts,1:3) )
     CALL assert_alloc(status,code,"survey")
     survey(1:nbpts,1:3) = surveytemp(1:nbpts,1:3)           
     DEALLOCATE( surveytemp )

     !Performing Fourier-Bessel decomposition
     almntemp = cmplx(0.0,0.0,kind=dp)
     CALL survey2almn_srs( nside, nnmax, nlmax, nmmax, rmax, nbpts, &
          & zbounds, survey, kln, almntemp)

     !Combining with previous results      
     almn = almn + almntemp
     ntot = ntot + nbpts

     DEALLOCATE(survey)
     !End of the procedure
  ENDDO
  print*," "
  PRINT*,ntot, "points in the sample"
  print*," "
  CLOSE(2)

  !---------------------------------------------------------------------------
  ! Spectrum

  header(:) = ' '
  nlheader = 120
  print*,"------------------------------------------------"
  print*,"   Constructing the spectrum..."
  !Computing the power spectrum from the almn's
  ALLOCATE( spectr(0:nlmax,1:nnmax), stat = status )
  CALL assert_alloc(status,code,"spectrnl")
  spectr = 0.0_dp
  CALL almn2cln_naive( nnmax, nlmax, nmmax, almn, spectr)
  print*,"   Done."
  print*,"------------------------------------------------"
  print*," "
  print*,"   DONE: DECOMPOSITION AND SPECTRUM"
  !If necessary : print the results
  print*," "
  print*,"------------------------------------------------"
  print*,"   Printing the first terms of the decomposition"
  print*,"------------------------------------------------"
  !CALL print_spectrum(cln, 3, 3, nlmax, nnmax," cln")
  !CALL print_spectrum(kln, 3, 3, nlmax, nnmax," kln")
  !CALL print_almn(almn, 10, 10, 10, nlmax, nmmax, nnmax,"almn")
  !CALL print_spectrum(spectr, 10, 10, nnmax, nlmax,"sptr")
  print*,"------------------------------------------------"

  time2 = wtime()

  !---------------------------------------------------------------------------
  ! End of the procedure	
  print*," "	
  PRINT*,"================================================"
  PRINT*," "
  PRINT*,"DECOMPOSITION (survey2almn) - SUMMARY"
  PRINT*," "
  PRINT*,"Survey size :",ntot
  PRINT*,"nside / rmax", nside, rmax
  PRINT*,"(nmax,lmax,mmax) :", nnmax, nlmax, nmmax
  PRINT*,"Elapsed time :",time2-time1
  PRINT*," "
  PRINT*,"================================================"
  PRINT*," "
  !Fits output
  PRINT*,"Starting almn and cln output..."
  CALL almn2fits( almnfile, almn, kln, cln, nlmax, nmmax, nnmax, header, nlheader )
  CALL cln2fits( clnfile, spectr, kln, nlmax, nnmax, header, 120 )
  PRINT*,"Done"
  PRINT*," "
  PRINT*,"================================================"
  print*," "
  print*," "

END PROGRAM survey2almn
