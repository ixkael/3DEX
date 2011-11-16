PROGRAM compute_qml_pl

  ! gfortran -O3 -fno-second-underscore -I$HEALPIX/include -L$HEALPIX/lib -L/usr/local/lib -lhealpix -lcfitsio -I$F3DEX/include -L$F3DEX/lib -l3dex -llapack -lbla -c compute_pl.f90 -o compute_pl.o

  ! gfortran -O3 -fno-second-underscore -o compute_pl compute_pl.o -I$HEALPIX/include -L$HEALPIX/lib -L/usr/local/lib -lhealpix -lcfitsio -I$F3DEX/include -L$F3DEX/lib -llapack -lblas -l3dex


  USE f3dex_transforms
  USE f3dex_cosmotools
  USE f3dex_stats
  USE f3dex_fitstools
  USE f3dex_utils
  USE f3dex_utests

  USE healpix_types
  USE healpix_modules

  INTEGER(I4B) :: nbpts, status, n_args,  npix, nside, nlmax, l
  CHARACTER(len=FILENAMELEN) ::  filename, arg, nlmaxchar, nsidechar, hpxmode
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Pl
  CHARACTER(len=*), PARAMETER   :: code = "compute_qml_pl"
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: ang_map, Pl1, Pl2
  LOGICAL :: fileexists
  
  CALL message(code,start=.TRUE.)
  
  n_args = nArguments()
  
  CALL getArgument(1, arg)
  filename = trim(adjustl(arg))
  INQUIRE( file=filename, exist=fileexists)
  IF( fileexists .eqv. .TRUE. ) CALL fatal_error("filename : this file already exists!")

  call getArgument(2, arg)
  read(arg,'(I4)') nlmax

  call getArgument(3, arg)
  read(arg,'(I4)') nside
  
  CALL getArgument(4, arg)
  hpxmode = trim(adjustl(arg))
  
  npix=12*nside**2 
  
  ALLOCATE(ang_map(0:npix-1,0:npix-1),stat = status)
  CALL assert_alloc(status,code,'ang_map')
  
  IF( hpxmode == "nested" .OR. hpxmode == "NESTED" ) THEN
    CALL message(code,msg="Nested scheme")
     CALL gen_angles_map(npix, ang_map, nside, nested=.TRUE.)
  ELSE IF( hpxmode == "ring" .OR. hpxmode == "RING" ) THEN
     CALL message(code,msg="Ring scheme")
     CALL gen_angles_map(npix, ang_map, nside, nested=.FALSE.)  
  ELSE
     PRINT*,"ERROR : you must provide NESTED or RING as fourth argument"
     STOP
  ENDIF
  
  ALLOCATE(Pl(0:npix-1,0:npix-1,0:nlmax),stat = status)
  CALL assert_alloc(status,code,'Pl')
  
  CALL message(code,msg="Generating pl matrices using recurrence")
  DO l=0, nlmax
   CALL message(code,msg="l",i=l,msg2="/",i2=nlmax)
       IF( l .EQ. 0 )THEN
    	    Pl(0:npix-1,0:npix-1,0) = 1.0
       ELSE IF( l .EQ. 1 ) THEN
       	    Pl(0:npix-1,0:npix-1,1) = ang_map
       ELSE
       	   Pl(0:npix-1,0:npix-1,l) = &
       	    & ( (2.0*REAL(l)-1.0) * ang_map * Pl(0:npix-1,0:npix-1,l-1) &
       	    & - REAL(l-1) * Pl(0:npix-1,0:npix-1,l-2) ) / REAL(l)
       END IF
       !Pl(0:npix-1,0:npix-1,l) = Pl(0:npix-1,0:npix-1,l)
  ENDDO
  
  DO l = 0, nlmax
     Pl(0:npix-1,0:npix-1,l) = Pl(0:npix-1,0:npix-1,l) * (2.0_dp*REAL(l)+1.0_dp)/FOURPI
  ENDDO

  !WRITE(nlmaxchar,'(I2)') nlmax
  !WRITE(nsidechar,'(I2)') nside
  !filename = "qm_plmat_0_" // trim(nlmaxchar) // "_" &
  !	     & // trim(nsidechar) // ".fits"
  !print*,"FILE : ",trim(filename)
  
  PRINT*,"Exporting coefficients"
  CALL tritab2fits( filename, Pl, 0, npix-1, 0, npix-1, 0, nlmax )
  PRINT*,"Done."
 
  !CALL fits2tritab( filename, Pl, 0, npix-1, 0, npix-1, 0, nlmax )

  DEALLOCATE( Pl, ang_map )
  
  CALL message(code,fin=.TRUE.)
  
END PROGRAM compute_qml_pl
