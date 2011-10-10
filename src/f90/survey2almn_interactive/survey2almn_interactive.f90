PROGRAM survey2almn_interactive

    USE f3dex
    USE cosmotools
    USE healpix_types
    USE healpix_modules
    USE misc_utils
    USE fitstools
    
    ! ------------------------------------------------------------------
    !!		Declarations
    ! ------------------------------------------------------------------
    
    character(len=80), DIMENSION(1:120) 	   :: header
    character(len=FILENAMELEN)   		   :: outfile
	
    integer(i4b)                  		   :: n_args
    character(len=FILENAMELEN)    		   :: infile, arg, surveyfile
    character(len=FILENAMELEN)    		   :: almnoutfile, clnoutfile
    character(len=*), PARAMETER   		   :: code = "Main program"
    integer(kind=i4b), DIMENSION(1:3)  		   :: iwhere
        
    real(kind=DP),   DIMENSION(:,:),   ALLOCATABLE :: cln, plm, spectrnl
    real(kind=DP),   DIMENSION(:,:),   ALLOCATABLE :: map
    real(kind=DP),   DIMENSION(:,:)  , ALLOCATABLE :: survey 
    real(kind=DP),   DIMENSION(:,:)    , ALLOCATABLE :: kln  
    complex(kind=DP),DIMENSION(:,:,:), ALLOCATABLE :: almnHP
    
    complex(kind=DP) :: temp
    complex(kind=DP), dimension(15,15,30) :: rholmn
    
    real(DP),DIMENSION(:,:), ALLOCATABLE :: testmap
        
    integer(i4b)  :: iter_order, irho, inradian, status, n_plm, nstemp
    integer(i4b)  :: nsmax, multiscale, sorted, rdstatus
    integer(i4b)  :: nnmax, nlmax, nmmax, nr, nbpts
    integer(i4b)  :: i, j, k, p, l, m, n
    real(kind=DP),     DIMENSION(1:2)  :: zbounds
    real*8	  :: tempval,time1,time2,time3,time4,time5,time6
    real(kind=DP) :: rmax, convfc, rhomean
    real(kind=DP) :: rho
    real(kind=DP), DIMENSION(:), ALLOCATABLE :: tempsurveyconv
    real(DP) :: h, omega_m, omega_l, omega_b, wa, w0, rmaxinput
    
    complex(kind=SP),DIMENSION(:,:,:), ALLOCATABLE :: almn
        
    
    ! ------------------------------------------------------------------
    !!	 	Parameters extraction
    ! ------------------------------------------------------------------
    	   
    print*,"========================================"
    print*,"    3DEX -"
    print*,"    Backward Fourier-Bessel Expansion"
    print*,"========================================"
    	
    time1 = wtime()
       
    n_args = nArguments()
	
    if( n_args == 1 ) then
        call getArgument(1, arg)
        infile = trim(adjustl(arg))
        call getParameters_survey2almn( multiscale, zbounds, nr, nsmax, nnmax, nlmax, nmmax, iter_order, &
	     & nbpts, iwhere, inradian, sorted, convfc, nbcolsinfile, surveyfile, infile, &
	     & almnoutfile, clnoutfile, h, w0, wa, omega_l, omega_b, omega_m, rmaxinput )
    else
    	infile = ""
    	call getParameters_survey2almn( multiscale, zbounds, nr, nsmax, nnmax, nlmax, nmmax, iter_order, &
	     & nbpts, iwhere, inradian, sorted, convfc, nbcolsinfile, surveyfile, infile, &
	     & almnoutfile, clnoutfile, h, w0, wa, omega_l, omega_b, omega_m, rmaxinput )
    endif	

    
    ! ------------------------------------------------------------------
    !!		Data extraction
    ! ------------------------------------------------------------------
    	    
    print*," "
    print*,"Extracting data..."
    
    ALLOCATE(survey(1:nbpts, 1:3),stat = status)
    call assert_alloc(status,code,"survey")
    
    call extractFromFile( nbpts, nr, iwhere, inradian, sorted, convfc, nbcolsinfile, surveyfile, survey )
    print*,"Data extraction...done"


	if(sorted == 1) then
	
    	print*," "
    	print*,"Conversion to comoving space"
   
   		ALLOCATE(tempsurveyconv(1:nbpts),stat = status)
    	CALL cosmo_z2s( h, omega_m, omega_l, omega_b, wa, w0, survey(:,3) , nbpts, tempsurveyconv )
    	!print*,tempsurveyconv
    	survey(:,3) = tempsurveyconv
    	DEALLOCATE(tempsurveyconv)
   	 	print*,"Conversion...done"
   	 	
   	endif

    rmax = max( rmaxinput,maxval(survey(:,3)) )
    print*,"RMAX",rmax
    
    	
    ! ------------------------------------------------------------------
    !!		Knl's and Cln's computation
    ! ------------------------------------------------------------------

    !Computing kln's
    PRINT*,"Computing the kln's"
    ALLOCATE(kln(0:nlmax,1:nnmax),stat = status)
    !DO i = 0, nlmax
    !   call logrange( kln(i:i,:) , 0.001_DP , 0.5_DP , nnmax )
    !ENDDO
    CALL gen_kln(kln, nnmax, nlmax, rmax)
    PRINT*,"kln's...done"
    
    !Computing cln's
    PRINT*,"Computing the cln's…"
    ALLOCATE(cln(0:nlmax,1:nnmax),stat = status)
    CALL assert_alloc(status,code,"cln")
    CALL gen_cln(cln, kln, nnmax, nlmax, rmax)
    PRINT*,"cln's done" 
  
    ! ------------------------------------------------------------------
    !!		3DEX_reversed calculations
    ! ------------------------------------------------------------------
    	
    print*,""
    print*,"-----------------------------------"
    print*,"Starting multiresolution"
    print*," spherical Healpix analysis"
    print*,"-----------------------------------"
    print*,""
    print*,"Entering almnHP computation..."

    time2 = wtime() 
    
    ALLOCATE( almnHP(1:nnmax,0:nlmax,0:nmmax), stat = status )
    call assert_alloc(status,code,"almnHP") 
    almnHP = cmplx( 0.0, 0.0, kind=DP )
    

    CALL survey2almn_srs( nsmax, nnmax, nlmax, nmmax, rmax, nbpts, &
   			     & zbounds, survey, kln, almnHP)
   	      
    !almnHP = almnHP / real(nbpts)
        
    print*,"almnHP computation... done"			
    time3 = wtime() 
    	
    header(:) = ' '
    nlheader = 120
    
    if( .false. ) then
    	
    	if( .false. ) then
    	
    	nbpts = 0
    	open(1,file="out/rholmn.out", status='old', iostat=rdstatus)
    	do while(rdstatus == 0)  
	    read(1,*,iostat = rdstatus) 
	    if (rdstatus == 0) then
	        nbpts = nbpts + 1
	    endif
	enddo 
	print*,"Number of points : ",nbpts
	rewind(unit=1,iostat=rdstatus)
    	do i=1, nbpts
    	    read(1,*) l,m,n,temp
    	    if( m > -1 ) then
    	    	print*,"----------------------------------------"
    	    	print*,temp
    	    	print*,almnHP(n,l,m)
    	    	print*,"----------------------------------------"
    	    	rholmn(l,m,n) = temp 
    	    endif
    	enddo
    	print*,"Additional import : done"
    	    	   
        endif
    	    
	    nstemp = 128
	    allocate( map(0:(12*nstemp**2),1:1) )
	    allocate( plm(0:(nstemp*(nmmax+1)*(2*nlmax-nmmax+2)),1:1) )
	    CALL plm_gen( nstemp, nlmax, nmmax, plm )
	    
	    print*,nstemp
	    	    
	    rho =100.0
	    CALL almn2rmap(map, almnHP, rho, nstemp, nnmax, nlmax, nmmax, kln, cln, plm)
	    call write_minimal_header(header, 'MAP', nside=nstemp, ordering='Ring') 
	    CALL write_bintab(map, (12*nstemp**2), 1, header, 120,'out/map_100.fits')
	    	    
	    rho =140.0
	    CALL almn2rmap(map, almnHP, rho, nstemp, nnmax, nlmax, nmmax, kln, cln, plm)
	    call write_minimal_header(header, 'MAP', nside=nstemp, ordering='Ring') 
	    CALL write_bintab(map, (12*nstemp**2), 1, header, 120,'out/map_140.fits')
		
	    rho =180.0
	    CALL almn2rmap(map, almnHP, rho, nstemp, nnmax, nlmax, nmmax, kln, cln, plm)
	    call write_minimal_header(header, 'MAP', nside=nstemp, ordering='Ring') 
	    CALL write_bintab(map, (12*nstemp**2), 1, header, 120,'out/map_180.fits')
	   	    
	    rho =220.0
	    CALL almn2rmap(map, almnHP, rho, nstemp, nnmax, nlmax, nmmax, kln, cln, plm)
	    call write_minimal_header(header, 'MAP', nside=nstemp, ordering='Ring') 
	    CALL write_bintab(map, (12*nstemp**2), 1, header, 120,'out/map_220.fits')
	    	    
	    rho =260.0
	    CALL almn2rmap(map, almnHP, rho, nstemp, nnmax, nlmax, nmmax, kln, cln, plm)
	    call write_minimal_header(header, 'MAP', nside=nstemp, ordering='Ring') 
	    CALL write_bintab(map, (12*nstemp**2), 1, header, 120,'out/map_260.fits')
	    	    
	    rho =300.0
	    CALL almn2rmap(map, almnHP, rho, nstemp, nnmax, nlmax, nmmax, kln, cln, plm)
	    call write_minimal_header(header, 'MAP', nside=nstemp, ordering='Ring') 
	    CALL write_bintab(map, (12*nstemp**2), 1, header, 120,'out/map_300.fits')

	    DEALLOCATE( map, plm , cln )
    endif
    
    
   

    
    
    !print*, almnHP
    
    !print*,spectrnl

    
    ! ------------------------------------------------------------------
    !!	 	OUTPUTS
    ! ------------------------------------------------------------------
  	  
    	
    print*,"-------------------------------------------------"
    print*,"-------------------------------------------------"
    print*,"-------------------------------------------------"
    print*,almnHP(1,0,0)
    print*,"-------------------------------------------------"
    print*,almnHP(1,1,0:1)
    print*,"-------------------------------------------------"
    print*,almnHP(1,2,0:2)
    print*,"-------------------------------------------------"
    print*,almnHP(1,3,0:3)
    print*,"-------------------------------------------------"
    print*,almnHP(1,4,0:3)
    print*,"-------------------------------------------------"
    print*,almnHP(1,5,0:3)
    print*,"-------------------------------------------------"
    print*,almnHP(1,6,0:3)
    print*,"-------------------------------------------------"
    print*,"-------------------------------------------------"
    print*,"-------------------------------------------------"
    print*,almnHP(2,0,0)
    print*,"-------------------------------------------------"
    print*,almnHP(2,1,0:1)
    print*,"-------------------------------------------------"
    print*,almnHP(2,2,0:2)
    print*,"-------------------------------------------------"
    print*,almnHP(2,3,0:3)
    print*,"-------------------------------------------------"
    print*,almnHP(2,4,0:3)
    print*,"-------------------------------------------------"
    print*,almnHP(2,5,0:3)
    print*,"-------------------------------------------------"
    print*,almnHP(2,6,0:3)
    print*,"-------------------------------------------------"
    print*,"-------------------------------------------------"
    print*,"-------------------------------------------------"

    
    
    print*," "
    print*,"Entering spectrum c(l,n) computation..."
    ALLOCATE( spectrnl(0:nlmax,1:nnmax), stat = status )
    CALL assert_alloc(status,code,"spectrnl")
   
    ! Call the appropriate function to compute the power spectrum
    CALL almn2cln( nnmax, nlmax, nmmax, almnHP, spectrnl)
    print*,"Spectrum computation... done"
    
    print*,"-------------------------------------------------"
    print*,"-------------------------------------------------"
    print*,"-------------------------------------------------"
    print*,spectrnl(1,0:3)
    print*,"-------------------------------------------------"
    print*,spectrnl(1,4:7)
    print*,"-------------------------------------------------"
    print*,spectrnl(2,0:3)
    print*,"-------------------------------------------------"
    print*,spectrnl(2,4:7)
    print*,"-------------------------------------------------"
    print*,spectrnl(3,0:3)
    print*,"-------------------------------------------------"
    print*,spectrnl(3,4:7)
    print*,"-------------------------------------------------"
    print*,"-------------------------------------------------"
    print*,"-------------------------------------------------"
    
    	
    	
    print*," "
    print*,"-----------------------------------"
    print*,"DATA OUTPUT"
    print*,"-----------------------------------"
    print*," "
    
    
    time4 = wtime()
    
    print*,"..."

        
    ! Output almn's
    if( almnoutfile /= "no" ) then
        CALL almn2fits( almnoutfile, almnHP, kln, cln, nlmax, nmmax, nnmax, header, 120 )
    endif
        
    
    ! Output cln's
    if( clnoutfile /= "no" ) then
    	CALL cln2fits( clnoutfile, spectrnl, kln, nlmax, nnmax, header, 120 )
    endif
    
    
    
    print*,"DONE!"
    
    time5 = wtime()
    
    ! ------------------------------------------------------------------
    !!	 	FLAGS
    ! ------------------------------------------------------------------
  	   
    print*,""
    print*,"-----------------------------------"
    print*,"PROGRAM TERMINATED"
    print*,"-----------------------------------"
    print*,""
    
    print*,"===================================================================="
    print*," "
    print*,"ELAPSED TIME"
    print*," "
    print*,"Data extraction, transformation and allocation : ",time2-time1
    print*," "
    print*,"MultiHealpix calculations, computation of almn : ",time3-time2
    print*," "
    print*,"Exportation : almns and cnls under fits format : ",time5-time4
    print*," "
    print*,"===================================================================="
    
    
    ! ------------------------------------------------------------------
    !!	 	DEALLOCATION
    ! ------------------------------------------------------------------
    	    
    deallocate( survey, almnHP, spectrnl )
    
END PROGRAM
