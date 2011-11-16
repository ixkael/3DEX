MODULE f3dex_cosmotools

  USE f3dex_utils
  USE f3dex_stats

  USE healpix_types
  USE healpix_modules

  IMPLICIT NONE

  REAL(DP), PARAMETER :: c = 299792458_dp  	! Speed of light
  REAL(DP), PARAMETER :: tcmb = 2726000_dp 	! CMB temperature [micro Kelvin]

  TYPE recognizedParams
     character(len=FILENAMELEN) :: surveyfile, almnoutfile, clnoutfile
     INTEGER(I4B) :: nbcolsinfile, inradian, sorted, nside, nblayers, multiscale
     REAL(DP) :: theta_cut_deg, convfc, h, omega_l, omega_m, omega_b, wa, w0
     REAL(DP) :: rmaxinput
     INTEGER(I4B) :: nnmax, nlmax, nmmax,iter_order, nbptsasked
     INTEGER(I4B), DIMENSION(1:3) :: iwhere
  END TYPE recognizedParams

CONTAINS  

  ! ---------------------------------------------------------------------------------------!    

  !> Equatorial to galactic coordinates conversion (in degrees)
  !!@param[in] (delta, alpha) : Doubles (equatorial)
  !!@param[in] (b, l) : Doubles (galactic)
  SUBROUTINE equatorial2galactic_deg( delta, alpha, b, l )

    REAL(DP) :: delta, alpha, b, l

    b = asin( sin(delta)*cos(62.6) - cos(delta)*sin(alpha-282.25)*sin(62.6) )
    l = 33 + atan( (sin(delta)*sin(62.6)+cos(delta)*sin(alpha-282.25)*cos(62.6)) / &
         & (cos(delta)*cos(alpha-282.25)))

  END SUBROUTINE equatorial2galactic_deg

  ! ---------------------------------------------------------------------------------------!    

  !> Inquire file and delete it if necessary
  !!@param[in] filename : string
  SUBROUTINE inquireAndDelete( filename )

    CHARACTER(len=FILENAMELEN) :: filename, cmd
    LOGICAL :: filefound

    INQUIRE(file=filename, exist=filefound)
    IF(filefound) THEN
       WRITE(cmd, '("rm ", A)' ) trim(filename)
       CALL system(cmd)
       PRINT*,"> Found and deleted: ",trim(filename)
    ENDIF

  END SUBROUTINE inquireAndDelete

  ! ---------------------------------------------------------------------------------------!    

  !> Read list of coordinates from a file, using a larger array
  !!@param[in] file : input file
  !!@param[in] file : (nbptsmax, nbcolsinfile) : bounds of the temp array / of the survey
  !!@param[out] surveytemp : output array
  !!@param[out] nbpts : actual number of coordinates read
  SUBROUTINE read_survey(file, nbptsmax, nbcolsinfile, nbpts, surveytemp)

    CHARACTER(len=FILENAMELEN) :: file
    INTEGER(I4B) :: rdstatus, status, nbptsmax, nbpts, nbcolsinfile, nbcols, i, k
    CHARACTER(len=*), PARAMETER   :: code = "read_survey"
    REAL(DP), DIMENSION(1:nbptsmax, 1:nbcolsinfile) :: surveytemp

    nbpts = 0
    OPEN(2,file=file,iostat=rdstatus)
    
    DO WHILE( nbpts < nbptsmax .AND. rdstatus == 0  )  
       DO i=1,nbptsmax
          nbpts = i-1
          READ(2,*,iostat = rdstatus,END=61) (surveytemp(i,k),k=1,nbcolsinfile) 
          nbpts = i
       ENDDO
    ENDDO
61  CONTINUE

  END SUBROUTINE read_survey

  ! ---------------------------------------------------------------------------------------!    

  !> Newton cotes integration using QUADRULE package
  !!@param[in] x(1:n) : input abs
  !!@param[in] f(1:n) : input coords
  !!@param[in] n : dimension
  REAL(DP) FUNCTION int_NC(x, f, n)
    ! Newton-Cotes integration, for f(x) unIFormly discretized

    REAL(DP), DIMENSION(1:n) :: x, f, weights
    INTEGER(I4B) :: n
    REAL(DP) :: integral

    ! Uses Newton-Cotes integration
    ! From the QUADRULE package

    CALL nc_com( 5, x(1), x(n), x, weights )

    int_NC = sum( weights*f )

    RETURN
  END FUNCTION int_NC

  ! ---------------------------------------------------------------------------------------!    

  !> Converts series of redshift values into radial coordinates
  !!@param[in] (h_in, omega_m_in, omega_l_in, omega_b_in, wa_in, w0_in) : cosmological parameters
  !!@param[in] z(1:nbpts_in) : input array
  !!@param[in] nbpts_in : input dimension
  !!@param[out] sk(1:nbpts_in) : output array
  SUBROUTINE cosmo_z2s( h_in, omega_m_in, omega_l_in, omega_b_in, wa_in, w0_in, z, nbpts_in, sk )   ! need in cosmo : h, omega_b, omega_m, omega_l

    INTEGER(I4B) :: nbpts_in, nbpts,i
    REAL(DP) :: h, omega_m, omega_l, omega_b, omega_k, wa, w0, sqrtk
    REAL(DP) :: h_in, omega_m_in, omega_l_in, omega_b_in, wa_in, w0_in
    COMMON h, omega_m, omega_l, omega_b, omega_k, wa, w0, sqrtk, nbpts

    REAL(DP) :: omega_dm, omega, gamma, h0, rh, k, r0

    REAL(DP), DIMENSION(1:nbpts_in) :: z, chi, sk

    h = h_in
    omega_m = omega_m_in
    omega_l = omega_l_in
    omega_b = omega_b_in
    wa = wa_in
    w0 = w0_in
    nbpts = nbpts_in

    h0 = 100.0_dp * h
    rh = c/(1000*h0)

    omega_dm = omega_m - omega_b 	! Dark matter density
    omega = omega_m + omega_l      	! Total density
    omega_k = 1.0_dp - omega           		! Curvature
    gamma = omega_m * h * &				! Sugiyama (1995, APJS, 100, 281)
         & exp(-omega_b*(1.0_dp + sqrt(2.0_dp*h)/omega_m)) 

    IF( omega > 1.001 ) THEN 			! closed
       k = 1.0_dp                         
       sqrtk = sqrt(omega-1.0_dp)   	! useful factor
    ENDIF                       			! with  a=r/r0
    IF( omega > 0.999 .and. omega < 1.001 ) THEN
       k = 0.0_dp                   		! flat
       sqrtk = 1.0_dp
    ENDIF
    IF( omega < 0.999 ) THEN        	! open
       k = -1.0_dp  
       sqrtk = sqrt(1.0_dp - omega)       
    ENDIF

    r0 = rh / sqrtk  

    ! comoving distance in units of R_0  [1] 
    PRINT*,"Computing chi..."
    !CALL chi_z1(z,chi)
    DO i=1, nbpts
       CALL chi_z0(z(i),chi(i))
    ENDDO

    ! sk(chi): radial comoving distance in units of R_0  [1],
    IF( omega > 1.001 ) THEN 
       sk = chi
    ENDIF
    IF( omega > 0.999 .and. omega < 1.001 ) THEN 
       sk = sin(chi)
    ENDIF
    IF( omega < 0.999 ) THEN 
       sk = sinh(chi)
    ENDIF

    sk = abs(sk) * r0

    RETURN

  END SUBROUTINE cosmo_z2s


  ! ---------------------------------------------------------------------------------------!

  !> Raw chi to rad conversion, using QUADPACK integration
  !!@param[in] zin : input value
  !!@param[out] chi : output value
  SUBROUTINE chi_z0(zin, chi)

    REAL(DP) :: zin, chi
    REAL, parameter :: a = 0.0E+00
    REAL :: abserr, b
    REAL, parameter :: epsabs = 0.001E+00
    REAL, parameter :: epsrel = 0.001E+00
    INTEGER :: ier
    INTEGER, parameter :: key = 6
    INTEGER :: neval, nb, i
    REAL, dimension(1:30) :: points
    REAL :: res

    b = REAL(zin)
    nb = 30
    DO i=1,nb
       points(i) = i*b/30
    ENDDO

    CALL qagp ( fz, a, b, nb, points, epsabs, epsrel, res, abserr, neval, ier )
    chi = res

  END SUBROUTINE chi_z0

  ! ---------------------------------------------------------------------------------------!

  !> Raw chi to rad conversion, using QUADPACK integration
  !!@param[in] zin : input value
  !!@param[out] chi : output value
  SUBROUTINE chi_z1(zin, chi)

    INTEGER(I4B) :: nbpts
    REAL(DP) :: h, omega_m, omega_l, omega_b, omega_k, wa, w0, sqrtk
    COMMON h, omega_m, omega_l, omega_b, omega_k, wa, w0, sqrtk, nbpts

    INTEGER(I4B) :: i, k, nb
    REAL(DP), DIMENSION(1:nbpts) :: zin,chi
    REAL(DP), DIMENSION(1:20) :: bpoints
    REAL, parameter :: a = 0.0E+00
    REAL :: abserr, b0, b1, res
    REAL, parameter :: epsabs = 0.001E+00
    REAL, parameter :: epsrel = 0.001E+00
    INTEGER :: ier, neval
    INTEGER, parameter :: key = 6

    nb = 20

    DO i = 1, nbpts

       b0 = REAL(zin(i-1)) !REAL(INT(zin(i-1) * 10000.0 + 0.5))  / 10000.0
       b1 = REAL(zin(i)) !REAL(INT(zin(i) * 10000.0 + 0.5))  / 10000.0
       res = 0.0

       IF( b1 > b0 ) THEN
          CALL qag ( fz, b0, b1, epsabs, epsrel, key, res, abserr, neval, ier )
          chi(i) = chi(i-1) + res
       ELSE
          CALL qag ( fz, a, b1, epsabs, epsrel, key, res, abserr, neval, ier )
          chi(i) = res
       ENDIF

       IF( isnan(res) ) THEN
          DO k=1, 20
             bpoints(k) = k*b1/20
          ENDDO
          CALL qagp ( fz, a, b1, nb, bpoints, epsabs, epsrel, res, abserr, neval, ier )
       ENDIF

    ENDDO

  END SUBROUTINE chi_z1

  ! ---------------------------------------------------------------------------------------!

  !> Raw chi to rad conversion, using QUADRULE integration
  !!@param[in] zin : input value
  !!@param[out] chi : output value
  SUBROUTINE chi_z2(zin, chi)

    INTEGER(I4B) :: nbpts
    REAL(DP) :: h, omega_m, omega_l, omega_b, omega_k, wa, w0, sqrtk
    COMMON h, omega_m, omega_l, omega_b, omega_k, wa, w0, sqrtk, nbpts

    INTEGER(I4B) :: i,j, norder, nsub
    REAL(DP), DIMENSION(1:nbpts) :: zin, chi
    REAL(DP), DIMENSION(1:50) :: weights, xtab
    DOuble precision a, b, xlo, xhi, res

    a = 0.0D+00

    norder = 5
    nsub = 1

    xlo = -1.0D+00
    xhi = +1.0D+00

    DO i = 1, nbpts

       b = zin(i) !REAL(INT(zin(i) * 10000.0 + 0.5))   / 10000.0

       !DO j=1,1000
       !   xtab(j)=j*b/1000
       !ENDDO

       CALL ncc_set ( norder, xtab, weights )

       CALL sum_sub( fz, a, b, nsub, norder, xlo, xhi, xtab, weights, res )

       chi(i) = res
       PRINT*,b,chi(i)

    ENDDO

  END SUBROUTINE chi_z2

  ! ---------------------------------------------------------------------------------------!

  !> Cosmological function
  !!@param[in] z : redshift
  REAL FUNCTION fz(z)

    INTEGER(I4B) :: nbpts
    REAL(DP) :: h, omega_m, omega_l, omega_b, omega_k, wa, w0, sqrtk
    COMMON h, omega_m, omega_l, omega_b, omega_k, wa, w0, sqrtk, nbpts

    REAL(DP) :: z

    fz = ( omega_m*((1+z)**3) + omega_k*((1+z)**2) + omega_l )**(-0.5_dp) 

    RETURN

  END FUNCTION fz

  ! ---------------------------------------------------------------------------------------!

  !> Scope and scan for parameters
  !!@param[in] paramfile : input file
  !!@param[out] params : output values (special TYPE)
  SUBROUTINE txtfile2parameters(paramfile,params)

    IMPLICIT NONE
    character(len=FILENAMELEN) :: paramfile, description
    TYPE(recognizedParams)     :: params
    LOGICAL(LGT)       	    :: fileexists
    INTEGER(i4b)       	    :: nbparam, rdstatus, i, k, status
    character(len=12)  	    :: tmp1
    character(len=30)  	    :: tmp2
    character(len=6)    	    :: name
    character(len=30)  	    :: value
    character(len=100) 	    :: chline,chline1
    type(paramfile_handle)      :: handle


    params%nblayers = 0
    params%convfc = 1.0
    params%almnoutfile = ""
    params%clnoutfile = ""
    params%nbptsasked = -1
    params%multiscale = -1
    !params%nbptsinfile = -1
    params%surveyfile = ""
    params%nbcolsinfile = -1
    params%inradian = -1
    params%sorted = -1
    params%nside = -1
    params%nnmax = -1
    params%nlmax = -1
    params%nmmax = -1
    params%iter_order = -1
    !params%coordsys = 'N'
    params%iwhere = -1
    params%theta_cut_deg = -1.0
    params%h = -1
    params%w0 = -10 
    params%wa = -10
    params%omega_l = -1 
    params%omega_b = -1
    params%omega_m = -1
    params%rmaxinput = -1.0_dp


    PRINT*,""
    PRINT*,"-----------------------------------"
    PRINT*,"Reading parameters"
    PRINT*,"-----------------------------------"
    PRINT*,""
    handle = parse_init('')
    IF( paramfile /= "" ) THEN
       PRINT*,"	"
       PRINT*,"File to be read : ",paramfile(1:20)
       inquire( file=paramfile, exist=fileexists)
       IF( fileexists ) THEN
          PRINT*,"Checking file : ok existing."
       ELSE
          CALL fatal_error("Incorrect input file (parameters)")
       ENDIF
       open(1,file= paramfile, status='old', iostat=rdstatus)
       PRINT*,"Opening file status: ",rdstatus
       DO while(rdstatus == 0)  
          read(1,*,iostat = rdstatus) 
          IF (rdstatus == 0) THEN
             nbparam = nbparam + 1
          ENDIF
       END DO
       PRINT*,"Number of parameters provided in file : ", nbparam
       !
       rewind(unit=1,iostat=rdstatus)
       PRINT*,"Extraction of parameters..."
       DO i=1, nbparam
          READ (1,*,END=51) tmp1, tmp2
          PRINT*,tmp1,tmp2
          name = tmp1(1:1)//tmp1(3:3)//tmp1(5:5)&
               &//tmp1(7:7)//tmp1(9:9)//tmp1(11:11)
          value = tmp2(1:1)//tmp2(3:3)//tmp2(5:5)&
               &//tmp2(7:7)//tmp2(9:9)//tmp2(11:11)&
               &//tmp2(13:13)//tmp2(15:15)//tmp2(17:17)&
               &//tmp2(19:19)//tmp2(21:21)//tmp2(23:23)&
               &//tmp2(25:25)//tmp2(27:27)//tmp2(29:29)
          !    
          IF(trim(name)=="survey".or.trim(name)=="SURVEY") THEN
             params%surveyfile = trim(value)
             PRINT*,"* ", trim(name)," detected :   ", trim(value)
          ENDIF
          !    
          IF(trim(name)=="almn".or.trim(name)=="ALMN") THEN
             params%almnoutfile = trim(value)
             PRINT*,"*   ", trim(name)," detected :   ", trim(value)
          ENDIF
          !    
          IF(trim(name)=="cln".or.trim(name)=="CLN") THEN
             params%clnoutfile = trim(value)
             PRINT*,"*    ", trim(name)," detected :   ", trim(value)
          ENDIF
          !    
          IF(trim(name)=="nbcols".or.trim(name)=="NBCOLS") THEN
             read(value,'(I2)') params%nbcolsinfile
             PRINT*,"* ", trim(name)," detected : ", params%nbcolsinfile
          ENDIF
          !    
          IF(trim(name)=="nlayrs".or.trim(name)=="NLAYRS") THEN
             read(value,'(I4)') params%nblayers
             PRINT*,"* ", trim(name)," detected : ", params%nblayers
          ENDIF
          !    
          IF(trim(name)=="phi".or.trim(name)=="PHI") THEN
             read(value,'(I2)') params%iwhere(1)
             PRINT*,"* ", trim(name)," detected : ", params%iwhere(1)
          ENDIF
          !    
          IF(trim(name)=="theta".or.trim(name)=="THETA") THEN
             read(value,'(I2)') params%iwhere(2)
             PRINT*,"* ", trim(name)," detected : ", params%iwhere(2)
          ENDIF
          !    
          IF(trim(name)=="z".or.trim(name)=="Z") THEN
             read(value,'(I2)') params%iwhere(3)
             PRINT*,"*     ", trim(name)," detected : ", params%iwhere(3)
          ENDIF
          !    
          IF(trim(name)=="radian".or.trim(name)=="RADIAN") THEN
             read(value,'(I2)') params%inradian
             PRINT*,"* ", trim(name)," detected : ", params%inradian 
          ENDIF
          !
          IF(trim(name)=="sorted".or.trim(name)=="SORTED") THEN
             read(value,'(I2)') params%sorted
             PRINT*,"* ", trim(name)," detected : ", params%sorted  
          ENDIF
          !
          IF(trim(name)=="equcut".or.trim(name)=="EQUCUT") THEN
             IF( trim(value) == "no" ) THEN
                params%theta_cut_deg = 0.0
             ELSE
                read(value,'(F5.3)') params%theta_cut_deg	
             ENDIF
             PRINT*,"* ", trim(name)," detected :           ", trim(value)  
          ENDIF
          !
          IF(trim(name)=="nside".or.trim(name)=="NSIDE") THEN
             read(value,'(I4)') params%nside
             PRINT*,"*  ", trim(name)," detected : ", params%nside 
          ENDIF
          !
          IF(trim(name)=="h".or.trim(name)=="H") THEN
             read(value,'(F5.3)') params%h
             PRINT*,"*      ", trim(name)," detected : ", params%h 
          ENDIF
          !
          IF(trim(name)=="omg_l".or.trim(name)=="OMG_L") THEN
             read(value,'(F5.3)') params%omega_l
             PRINT*,"*  ", trim(name)," detected : ", params%omega_l 
          ENDIF
          !
          IF(trim(name)=="omg_m".or.trim(name)=="OMG_M") THEN
             read(value,'(F10.3)') params%omega_m
             PRINT*,"*  ", trim(name)," detected : ", params%omega_m 
          ENDIF
          !
          IF(trim(name)=="omg_b".or.trim(name)=="OMG_B") THEN
             read(value,'(F5.3)') params%omega_b
             PRINT*,"*  ", trim(name)," detected : ", params%omega_b 
          ENDIF
          !
          IF(trim(name)=="wa".or.trim(name)=="WA") THEN
             read(value,'(F5.3)') params%wa
             PRINT*,"*     ", trim(name)," detected : ", params%wa 
          ENDIF
          !
          IF(trim(name)=="convfc".or.trim(name)=="CONVFC") THEN
             read(value,'(F11.4)') params%convfc
             PRINT*,"*   ", trim(name)," detected : ", params%convfc 
          ENDIF
          !
          IF(trim(name)=="w0".or.trim(name)=="W0") THEN
             read(value,'(F5.3)') params%w0
             PRINT*,"*     ", trim(name)," detected : ", params%w0 
          ENDIF
          !
          IF(trim(name)=="nmax".or.trim(name)=="NMAX") THEN
             read(value,'(I4)') params%nnmax
             PRINT*,"*   ", trim(name)," detected : ", params%nnmax 
          ENDIF
          !
          IF(trim(name)=="lmax".or.trim(name)=="LMAX") THEN
             read(value,'(I4)') params%nlmax
             PRINT*,"*   ", trim(name)," detected : ", params%nlmax   
          ENDIF
          !
          IF(trim(name)=="mmax".or.trim(name)=="MMAX") THEN
             read(value,'(I4)') params%nmmax
             PRINT*,"*   ", trim(name)," detected : ", params%nmmax   
          ENDIF
          !
          IF(trim(name)=="itrord".or.trim(name)=="ITRORD") THEN
             read(value,'(I2)')  params%iter_order
             PRINT*,"* ", trim(name)," detected : ", params%iter_order
          ENDIF
          !
          IF(trim(name)=="multsc".or.trim(name)=="MULTSC") THEN
             read(value,'(I4)')  params%multiscale
             PRINT*,"* ", trim(name)," detected : ", params%multiscale
          ENDIF
          !
          IF(trim(name)=="nbrpts".or.trim(name)=="NBRPTS") THEN
             read(value,'(I7)')  params%nbptsasked
             PRINT*,"* ", trim(name)," detected : ", params%nbptsasked
          ENDIF
          !
          IF(trim(name)=="rmax".or.trim(name)=="RMAX") THEN
             read(value,'(F11.4)')  params%rmaxinput
             PRINT*,"*   ", trim(name)," detected : ", params%rmaxinput
          ENDIF
          !
       ENDDO
       PRINT*,"-- DOne."
       PRINT*,"	"
51     continue
       close(1)
    ELSE
       description = concatnl("Enter input file (galaxy survey) name (eg: data.out) : ")
       params%surveyfile = parse_string(handle, 'surveyfile', default='', descr=description)
       inquire( file=params%surveyfile, exist=fileexists)
       IF( fileexists .eqv. .FALSE. ) CALL fatal_error("This file DOesn't exist!")
       IF (trim(params%surveyfile)=='') CALL fatal_error('Error : no input file provided')    
       PRINT*,"	"
       PRINT*,"	"
    ENDIF


  END SUBROUTINE txtfile2parameters

  ! ---------------------------------------------------------------------------------------!    


  !> Extracts parameters
  SUBROUTINE getParameters_survey2almn(multiscale, zbounds, nr, nside, nnmax, nlmax, nmmax, iter_order, &
       & nbpts, iwhere, inradian, sorted, convfc, nbcolsinfile, surveyfile, paramfile, &
       & almnoutfile, clnoutfile, h, w0, wa, omega_l, omega_b, omega_m, rmaxinput )

    ! ------------------------------------------------------------------------------
    ! Get all the parameters of the method
    ! ------------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER(i4b)       :: nbparam, rdstatus, i, k, status
    LOGICAL(LGT)       :: fileexists
    character(len=12)  :: tmp1
    character(len=30)  :: tmp2
    character(len=6)   :: name
    character(len=30)  :: value
    character(len=100) :: chline,chline1
    TYPE(recognizedParams)   :: params


    character(len=*), PARAMETER        :: code = "Parameters extraction"
    character(len=FILENAMELEN)         :: paramfile
    type(paramfile_handle)  	       :: handle
    character(len=FILENAMELEN)         :: surveyfile, description, almnoutfile, clnoutfile
    INTEGER(i4b)      		    	:: nbptsasked, nbptsinfile, nbcolsinfile
    INTEGER(i4b)  		      	    :: inradian, sorted, iter_order, nbpts, nr
    INTEGER(i4b)  		            :: nside, nnmax, nlmax, nmmax, multiscale
    character    		            :: coordsys
    INTEGER(kind=i4b), DIMENSION(1:3)  :: iwhere
    REAL(kind=DP)		            :: theta_cut_deg, cos_theta_cut, fsky, convfc
    REAL(kind=DP)		    :: rmaxinput, h, w0, wa, omega_l, omega_b, omega_m
    REAL(kind=DP),     DIMENSION(1:2)  :: zbounds

    ! ------------------------------------------------------------------
    !!		Initializations
    ! ------------------------------------------------------------------

    CALL txtfile2parameters(paramfile,params)

    handle = parse_init('')
    convfc = params%convfc
    almnoutfile = params%almnoutfile
    clnoutfile = params%clnoutfile
    nbptsasked = params%nbptsasked
    multiscale = params%multiscale
    nbptsinfile = nbptsasked
    surveyfile = params%surveyfile
    nr = params%nblayers
    nbcolsinfile = params%nbcolsinfile
    inradian = params%inradian
    sorted = params%sorted
    nside = params%nside
    nnmax = params%nnmax
    nlmax = params%nlmax
    nmmax = params%nmmax
    rmaxinput = params%rmaxinput
    iter_order = params%iter_order
    iwhere = params%iwhere
    theta_cut_deg = params%theta_cut_deg
    h = params%h
    w0 = params%w0
    wa = params%wa
    omega_l = params%omega_l
    omega_b = params%omega_b
    omega_m = params%omega_m


    ! ------------------------------------------------------------------
    !!		Ask for missing parameters
    ! ------------------------------------------------------------------

    inquire( file=surveyfile, exist=fileexists)
    IF( fileexists .eqv. .FALSE. ) THEN
       description = concatnl("Enter input file name (eg: data.out) : ")
       surveyfile = parse_string(handle, 'surveyfile', default='', descr=description)
       inquire( file=surveyfile, exist=fileexists)
       IF( fileexists .eqv. .FALSE. ) CALL fatal_error("This file DOesn't exist!")
       IF (trim(surveyfile)=='') CALL fatal_error('Error : no input file provided')
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( almnoutfile /= "no" ) THEN
       almnoutfile = "out/"//almnoutfile
       inquire( file=almnoutfile, exist=fileexists)
       IF( fileexists ) THEN
          PRINT*,"Error : almn outputfile already exists ; please provide another name"
          description = concatnl("Enter input file name (eg: almn.fits) : ")
          almnoutfile = parse_string(handle, 'almnoutfile', default='almn.fits', descr=description)
          almnoutfile = "out/"//almnoutfile
          inquire( file=almnoutfile, exist=fileexists)
          IF( fileexists .eqv. .TRUE. ) CALL fatal_error("This file already exists!")
          PRINT*,"	"
       ENDIF
    ENDIF
    !
    IF( clnoutfile /= "no" ) THEN
       clnoutfile = "out/"//clnoutfile
       inquire( file=clnoutfile, exist=fileexists)
       IF( fileexists ) THEN
          PRINT*,"Error : cln outputfile already exists ; please provide another name"
          description = concatnl("Enter input file name (eg: cln.fits) : ")
          clnoutfile = parse_string(handle, 'clnoutfile', default='cln.fits', descr=description)
          clnoutfile = "out/"//clnoutfile
          inquire( file=clnoutfile, exist=fileexists)
          IF( fileexists .eqv. .TRUE. ) CALL fatal_error("This file already exists!")
          PRINT*,"	"
       ENDIF
    ENDIF
    !
    IF( nbcolsinfile < 0 ) THEN
       description = concatnl("How many columns DOes this file contain? (Necessary for parsing)")
       nbcolsinfile = parse_int(handle, 'nbcolsinfile', default=-1, descr=description)
       IF( nbcolsinfile < 0 .or. nbcolsinfile > 20 ) CALL fatal_error("Error : Bad number.")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( iwhere(1)<1 .or. iwhere(1)>nbcolsinfile ) THEN
       WRITE(chline,"(a)") "In which column is Coord 1 ?"   
       description = concatnl( chline , "" , chline1 )
       iwhere(1) = parse_int(handle, 'Coord 1', default=-1, descr=description)
       IF( iwhere(1)<1 .or. iwhere(1)>nbcolsinfile ) CALL fatal_error("Incorrect column number")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( iwhere(2)<1 .or. iwhere(2)>nbcolsinfile ) THEN
       WRITE(chline,"(a)") "In which column is Coord 2 ?"
       description = concatnl( chline , "" , chline1 )
       iwhere(2) = parse_int(handle, 'Coord 2', default=-1, descr=description)
       IF( iwhere(2)<1 .or. iwhere(2)>nbcolsinfile ) CALL fatal_error("Incorrect column number")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( iwhere(3)<1 .or. iwhere(3)>nbcolsinfile ) THEN
       WRITE(chline,"(a)") "In which column is Coord 3 ?"   
       description = concatnl( chline , "" , chline1 )
       iwhere(3) = parse_int(handle, 'Coord 3', default=-1, descr=description)
       IF( iwhere(3)<1 .or. iwhere(3)>nbcolsinfile ) CALL fatal_error("Incorrect column number")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( inradian < 0 .or. inradian > 1 ) THEN
       description = concatnl("Are your data in radian?" ,&
	    &"","(0) No, in degree    (1) Yes, in radian")
       inradian = parse_int(handle, 'inradian', default=0, descr=description)
       IF( inradian < 0 .or. inradian > 1 ) CALL fatal_error("Error : Bad number.")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( sorted < 0 .or. sorted > 1 ) THEN
       description = concatnl("Is the array sorted by radial/redshIFt values? ",&
	    &"(0) No, no sorted    (1) Yes, sorted")
       sorted = parse_int(handle, 'sorted', default=0, descr=description)
       IF( sorted < 0 .or. sorted > 1 ) CALL fatal_error("Error : Bad number.")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( nside < 0 .or. nside > 4096 ) THEN
       WRITE(chline,"(a)") "We recommEND: (256 <= Nside <= 2048) a power of 2"
       description = concatnl(&
	    & "Enter the Nside parameter (nsmax) for the healpix analysis. ", &
	    & chline )
       nside = parse_int(handle, 'nside', default=64, descr=description)
       IF( nside<0 .or. nside>1024 ) CALL fatal_error("Incorrect boundaries (admitted : 4,8,16,64,128,256,512)")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( nnmax < 0 .or. nnmax >= 250 ) THEN
       WRITE(chline,"(a)") "We recommEND: (1 <= n <= n_max <= 250)"
       description = concatnl(&
	    & "Enter the maximum n range (n_max) for the analysis (Bessel FUNCTIONs). ", &
	    & chline )
       nnmax = parse_int(handle, 'nnmax', default=9, descr=description)
       IF( nnmax<0 .or. nnmax>500 ) CALL fatal_error("Incorrect boundaries (admitted : [1,250])")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( nlmax < 0 .or. nlmax > 256 ) THEN
       WRITE(chline1,"(a,i5)") "The map has Nside = ",nside
       WRITE(chline,"(a,i5,a)") "We recommEND: (0 <= l <= l_max <= ",200,")"
       description = concatnl(&
	    & chline1, &
	    & "Enter the maximum l range (l_max) for the analysis. ", &
	    & chline )
       nlmax = parse_int(handle, 'nlmax', default=2*nside, descr=description)
       IF( nlmax<0 .or. nlmax>4*nside ) CALL fatal_error("Incorrect boundaries (admitted : [0,4*nside])")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( nmmax < 0 .or. nmmax > 256 ) THEN
       WRITE(chline,"(a,i5,a)") "We recommEND: (0 <= m <= m_max <= ",nlmax,")"
       description = concatnl(&
	    & " Enter the maximum m range (m_max) for the analysis. ", &
	    & chline )
       nmmax = parse_int(handle, 'nmmax', default=nlmax, descr=description)
       IF( nmmax<0 .or. nmmax>nlmax ) CALL fatal_error("Incorrect boundaries (admitted : [1,nlmax])")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( nr < 4 ) THEN
       description = concatnl(&
	    & " Enter the number of layers (radial discretization) " )
       nr = parse_int(handle, 'nr', default=8, descr=description)
       IF( nr<4 ) CALL fatal_error("Incorrect boundaries (admitted : [4,nlmax])")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( h < 0 ) THEN
       description = concatnl(&
	    & " Cosmogoly - parameter : h " )
       h = parse_DOuble(handle, 'h', default=0.7_dp, descr=description)
       IF( h<0 ) CALL fatal_error("Error")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( omega_b < 0 ) THEN
       description = concatnl(&
	    & " Cosmogoly - parameter : omega_b " )
       omega_b = parse_DOuble(handle, 'omega_b', default=0.04_dp, descr=description)
       IF( omega_b<0 ) CALL fatal_error("Error")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( omega_m < 0 ) THEN
       description = concatnl(&
	    & " Cosmogoly - parameter : omega_m " )
       omega_m = parse_DOuble(handle, 'omega_m', default=0.3_dp, descr=description)
       IF( omega_m<0 ) CALL fatal_error("Error")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( omega_l < 0 ) THEN
       description = concatnl(&
	    & " Cosmogoly - parameter : omega_l " )
       omega_l = parse_DOuble(handle, 'omega_l', default=0.7_dp, descr=description)
       IF( omega_l<0 ) CALL fatal_error("Error")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( wa < -5 ) THEN
       description = concatnl(&
	    & " Cosmogoly - parameter : wa " )
       wa = parse_DOuble(handle, 'wa', default=0.0_dp, descr=description)
       IF( wa<-5 ) CALL fatal_error("Error")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( w0 < -5 ) THEN
       description = concatnl(&
	    & " Cosmogoly - parameter : w0 " )
       w0 = parse_DOuble(handle, 'w0', default=-1.0_dp, descr=description)
       IF( w0<-5 ) CALL fatal_error("Error")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( multiscale < 0 .or. multiscale > 4 ) THEN
       description = concatnl(&
	    & " Enter the number of layers (radial discretization) " )
       multiscale = parse_int(handle, 'multiscalenr', default=1, descr=description)
       IF( multiscale<0 .or. multiscale > 4 ) CALL fatal_error("Incorrect boundaries (admitted : [1,3])")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( rmaxinput < 0 ) THEN
       description = concatnl(&
	    & " Maximal radial value - rmax " )
       rmaxinput = parse_DOuble(handle, 'rmaxinput', default=-1.0_dp, descr=description)
       IF( rmaxinput<0 ) CALL fatal_error("Error")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( iter_order < 0 .or. iter_order > 10 ) THEN
       description = concatnl(&
	    & " Do you want : ", &
	    & " 0) a standard analysis", &
	    & " 1,2,3,4....) an iterative analysis", &
	    & " (enter order of iteration, 3rd order is usually optimal)")
       iter_order=parse_int(handle, 'iter_order', vmin=0, default=0, descr=description)
       select case (iter_order)
       case (0)
          PRINT*," Standard analysis"
       case default
          PRINT*," Iterative analysis"
       END select
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( theta_cut_deg < 0.0 .or. theta_cut_deg > 1.0 ) THEN
       description = concatnl(&
	    & " Enter the symmetric cut around the equator in DEGREES : ", &
	    & " (One ignores data within |b| < b_cut)     0 <= b_cut = ")
       theta_cut_deg = parse_DOuble(handle, 'theta_cut_deg', &
	    &                       vmin=0.0_dp, default=0.0_dp, descr=description)
       cos_theta_cut = SIN(theta_cut_deg/180.d0*PI) !counted from equator instead of from pole
       zbounds = (/ cos_theta_cut , -cos_theta_cut /)
       IF (theta_cut_deg<1e-4) zbounds = (/ -1.0_dp, 1.0_dp /) !keep all sphere IF cut not set
       fsky = (zbounds(2)-zbounds(1))/2.0_dp
       IF (fsky <= 0.0_dp) fsky = 1.0_dp + fsky
       write(*,"(a,f6.1,a)") "One keeps ",100.*fsky," % of the original map(s)"   
       PRINT*,"	"
       PRINT*,"	"
    ENDIF
    !
    IF( nbptsasked <= 0 ) THEN
       description = concatnl("How many galaxies DO you want to take into account?")
       nbptsasked = parse_int(handle, 'nbptsasked', default=-1, descr=description)
       IF( nbptsasked <= 0 ) CALL fatal_error("Error : Bad number.")
       PRINT*,"	"
       PRINT*,"	"
    ENDIF

    PRINT*,"-----------------------------------"
    PRINT*,"Starting extraction from file"
    PRINT*,"-----------------------------------"
    PRINT*,""
    PRINT*,"File to be read : ",surveyfile(1:10)
    inquire( file=surveyfile, exist=fileexists)
    IF( fileexists ) THEN
       PRINT*,"Checking file ... ok existing."
    ELSE
       CALL fatal_error("Error: unknown file")
    ENDIF
    open(2,file= surveyfile, status='old', iostat=rdstatus, form="formatted")
    PRINT*,"Opening file status: ",rdstatus
    nbptsinfile = 0
    DO while(rdstatus == 0)  
       read(2,*,iostat = rdstatus) 
       IF (rdstatus == 0) THEN
          nbptsinfile = nbptsinfile + 1
       ENDIF
    END DO
    PRINT*,"Number of points read : ", nbptsinfile
    nbpts = min(nbptsasked,nbptsinfile)
    PRINT*,"Final number of galaxies to be analyzed : ",nbpts
    close(2)

  END SUBROUTINE getParameters_survey2almn

  ! ---------------------------------------------------------------------------------------!    


  !> Extracts parameters from file
  SUBROUTINE extractFromFile( nbpts, nr, iwhere, inradian, sorted, convfc, nbcolsinfile, surveyfile, survey )

    IMPLICIT NONE

    ! ------------------------------------------------------------------------------
    ! Extract data from input file and convert it into valid healpix 3D maps
    ! ------------------------------------------------------------------------------

    character(len=*), PARAMETER    :: code = "Survey extraction"
    character(len=FILENAMELEN)     :: surveyfile
    INTEGER(I4B)		    :: nbcolsinfile, status, nbpts, k, sorted
    INTEGER(I4B)		    :: i, inradian, rdstatus, nsmax, nr
    REAL(kind=DP), DIMENSION(1:nbpts) 	     :: mapcop
    REAL(kind=DP), DIMENSION(:), allocatable:: temp
    REAL(kind=DP), DIMENSION(1:nbpts,1:3)   :: survey, surveybis
    INTEGER(i4b),  DIMENSION(1:3)           :: iwhere
    REAL(DP) :: convfc

    ALLOCATE(temp(1:nbcolsinfile),stat = status)
    CALL assert_alloc(status,code,"temp")

    survey = 0.0
    surveybis = 0.0
    open(2,file= trim(surveyfile), status='old', iostat=rdstatus, form="formatted")
    rewind(unit=2,iostat=rdstatus)
    PRINT*,""
    PRINT*,"Extraction..."
    DO while(rdstatus == 0)
       DO i=1, nbpts
          READ (2,*,iostat = rdstatus, END=61) (temp(k),k=1,nbcolsinfile) 
          surveybis(i,1) = temp(iwhere(1))
          surveybis(i,2) = temp(iwhere(2))
          surveybis(i,3) = temp(iwhere(3))
       ENDDO
    ENDDO
61  continue
    PRINT*,"-- DOne"
    close(2)
    deALLOCATE(temp)

    ! Convertion from degree to radian IF necessary
    IF( inradian == 0 ) THEN
       PRINT*,""
       PRINT*,"Convertion from degree to radian..."
       DO i=1,nbpts
          surveybis(i,1) = surveybis(i,1)*PI/180
          surveybis(i,2) = surveybis(i,2)*PI/180 
          surveybis(i,2) = HALFPI - surveybis(i,2)            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ENDDO
       PRINT*,"-- DOne"
    ENDIF

    ! ------------------------------------------------------------------
    !!		Sorting the survey
    ! ------------------------------------------------------------------ 

    IF( sorted == 0 ) THEN
       PRINT*,"Sorting values."
       mapcop(1:nbpts) = surveybis(1:nbpts,3)
       CALL BSORT (mapcop, surveybis, nbpts)
       DO i = 1, nbpts
          survey(i,1:3) = surveybis(nbpts-i+1,1:3)
          survey(i,3) = survey(i,3) / convfc
       ENDDO
    ELSE
       survey = surveybis
       survey(:,3) = survey(:,3) / convfc
    ENDIF

    RETURN

  END SUBROUTINE extractFromFile

  ! ---------------------------------------------------------------------------------------!    

END MODULE f3dex_cosmotools
