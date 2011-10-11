PROGRAM almnfile2rmap

!COMMAND : almnfile2map almninfile mapoutfile rho nside nnmax nlmax rmax
!TEST : bin/almnfile2rmap out/almn.fits map.fits 400.0 256 10 10 2000.0
    
    USE f3dex
    USE cosmotools
    USE healpix_types
    USE healpix_modules
    USE misc_utils
    USE fitstools
    USE alm_tools
    
    REAL*8 :: time1, time2
    INTEGER(I4B) :: n_args, nnmax, nlmax, nmmax, status
    LOGICAL(LGT) :: renormalization, fileexists
    CHARACTER(len=80), DIMENSION(1:120):: header
    CHARACTER(len=FILENAMELEN) :: headersfile, arg, almninfile, clnoutfile
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: kln, cln, plm
    REAL(kind=DP), DIMENSION(:,:), ALLOCATABLE :: spectr 
    COMPLEX(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: almn
    CHARACTER(len=*), PARAMETER :: code = "almn2cln"
    
    PRINT*,"=================================================="
    PRINT*," 3DEX F90 LIBRARY"
    PRINT*," STARTING almn -> cln"
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
    clnoutfile = trim(adjustl(arg))
    INQUIRE( file=clnoutfile, exist=fileexists)
    IF( fileexists .eqv. .TRUE. ) CALL fatal_error("clnoutfile : this file already exists!")
  
    CALL getArgument(3, arg)
    READ(arg,'(I4)') nnmax
 
    CALL getArgument(4, arg)
    READ(arg,'(I4)') nlmax

        
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
    
    !Allocation
    ALLOCATE( spectr(0:nlmax,1:nnmax), stat = status )
    CALL assert_alloc(status,code,"spectrnl")
	
    !---------------------------------------------------------------------------
    ! Spectrum construction
	
    !Computing Spectrum
    PRINT*,"Starting reconstruction..."
    spectr = 0.0_dp
    CALL almn2cln( nnmax, nlmax, nmmax, almn, spectr)
    PRINT*,"Done"
    
    time2 = wtime()

    !---------------------------------------------------------------------------
    ! End of the procedure

    DEALLOCATE(almn, cln)
	
    PRINT*,"=================================================="
    PRINT*," "
    PRINT*,"SPECTRUM COMPUTATION (almn2cln) - SUMMARY"
    PRINT*," "
    PRINT*,"(nmax,lmax,mmax) :", nnmax, nlmax, nmmax
    PRINT*,"Elapsed time :",time2-time1
    PRINT*," "
    PRINT*,"=================================================="
    PRINT*," "
    !Fits output
    PRINT*,"Starting map output..."
    CALL cln2fits( clnoutfile, spectr, kln, nlmax, nnmax, header, 120 )
    PRINT*,"Done"
    PRINT*," "
    PRINT*,"=================================================="
    
    DEALLOCATE(spectr, kln)

END PROGRAM
