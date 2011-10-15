MODULE f3dex

  USE OMP_LIB

  USE healpix_types
  USE healpix_modules
  USE alm_tools
  USE fitstools
  IMPLICIT NONE

  real(kind=SP),     private, parameter :: s_bad_value = HPX_SBADVAL
  real(kind=DP),     private, parameter :: d_bad_value = HPX_DBADVAL

  TYPE recognizedParams
     character(len=FILENAMELEN) :: surveyfile, almnoutfile, clnoutfile
     INTEGER(I4B) :: nbcolsinfile, inradian, sorted, nside, nblayers, multiscale
     REAL(DP) :: theta_cut_deg, convfc, h, omega_l, omega_m, omega_b, wa, w0
     REAL(DP) :: rmaxinput
     INTEGER(I4B) :: nnmax, nlmax, nmmax,iter_order, nbptsasked
     INTEGER(I4B), DIMENSION(1:3) :: iwhere
  END TYPE recognizedParams


  !-------
CONTAINS
  !-------


  !===============================================================================
  !===============================================================================


  SUBROUTINE txtfile2parameters(paramfile,params)

    IMPLICIT NONE
    character(len=FILENAMELEN) :: paramfile, description
    TYPE(recognizedParams)     :: params
    logical(LGT)       	    :: fileexists
    integer(i4b)       	    :: nbparam, rdstatus, i, k, status
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


    print*,""
    print*,"-----------------------------------"
    print*,"Reading parameters"
    print*,"-----------------------------------"
    print*,""
    handle = parse_init('')
    if( paramfile /= "" ) then
       print*,"	"
       print*,"File to be read : ",paramfile(1:20)
       inquire( file=paramfile, exist=fileexists)
       if( fileexists ) then
          print*,"Checking file : ok existing."
       else
          call fatal_error("Incorrect input file (parameters)")
       endif
       open(1,file= paramfile, status='old', iostat=rdstatus)
       print*,"Opening file status: ",rdstatus
       do while(rdstatus == 0)  
          read(1,*,iostat = rdstatus) 
          if (rdstatus == 0) then
             nbparam = nbparam + 1
          endif
       end do
       print*,"Number of parameters provided in file : ", nbparam
       !
       rewind(unit=1,iostat=rdstatus)
       print*,"Extraction of parameters..."
       do i=1, nbparam
          READ (1,*,end=51) tmp1, tmp2
          print*,tmp1,tmp2
          name = tmp1(1:1)//tmp1(3:3)//tmp1(5:5)&
               &//tmp1(7:7)//tmp1(9:9)//tmp1(11:11)
          value = tmp2(1:1)//tmp2(3:3)//tmp2(5:5)&
               &//tmp2(7:7)//tmp2(9:9)//tmp2(11:11)&
               &//tmp2(13:13)//tmp2(15:15)//tmp2(17:17)&
               &//tmp2(19:19)//tmp2(21:21)//tmp2(23:23)&
               &//tmp2(25:25)//tmp2(27:27)//tmp2(29:29)
          !    
          if(trim(name)=="survey".or.trim(name)=="SURVEY") then
             params%surveyfile = trim(value)
             print*,"* ", trim(name)," detected :   ", trim(value)
          endif
          !    
          if(trim(name)=="almn".or.trim(name)=="ALMN") then
             params%almnoutfile = trim(value)
             print*,"*   ", trim(name)," detected :   ", trim(value)
          endif
          !    
          if(trim(name)=="cln".or.trim(name)=="CLN") then
             params%clnoutfile = trim(value)
             print*,"*    ", trim(name)," detected :   ", trim(value)
          endif
          !    
          if(trim(name)=="nbcols".or.trim(name)=="NBCOLS") then
             read(value,'(I2)') params%nbcolsinfile
             print*,"* ", trim(name)," detected : ", params%nbcolsinfile
          endif
          !    
          if(trim(name)=="nlayrs".or.trim(name)=="NLAYRS") then
             read(value,'(I4)') params%nblayers
             print*,"* ", trim(name)," detected : ", params%nblayers
          endif
          !    
          if(trim(name)=="phi".or.trim(name)=="PHI") then
             read(value,'(I2)') params%iwhere(1)
             print*,"* ", trim(name)," detected : ", params%iwhere(1)
          endif
          !    
          if(trim(name)=="theta".or.trim(name)=="THETA") then
             read(value,'(I2)') params%iwhere(2)
             print*,"* ", trim(name)," detected : ", params%iwhere(2)
          endif
          !    
          if(trim(name)=="z".or.trim(name)=="Z") then
             read(value,'(I2)') params%iwhere(3)
             print*,"*     ", trim(name)," detected : ", params%iwhere(3)
          endif
          !    
          if(trim(name)=="radian".or.trim(name)=="RADIAN") then
             read(value,'(I2)') params%inradian
             print*,"* ", trim(name)," detected : ", params%inradian 
          endif
          !
          if(trim(name)=="sorted".or.trim(name)=="SORTED") then
             read(value,'(I2)') params%sorted
             print*,"* ", trim(name)," detected : ", params%sorted  
          endif
          !
          if(trim(name)=="equcut".or.trim(name)=="EQUCUT") then
             if( trim(value) == "no" ) then
                params%theta_cut_deg = 0.0
             else
                read(value,'(F5.3)') params%theta_cut_deg	
             endif
             print*,"* ", trim(name)," detected :           ", trim(value)  
          endif
          !
          if(trim(name)=="nside".or.trim(name)=="NSIDE") then
             read(value,'(I4)') params%nside
             print*,"*  ", trim(name)," detected : ", params%nside 
          endif
          !
          if(trim(name)=="h".or.trim(name)=="H") then
             read(value,'(F5.3)') params%h
             print*,"*      ", trim(name)," detected : ", params%h 
          endif
          !
          if(trim(name)=="omg_l".or.trim(name)=="OMG_L") then
             read(value,'(F5.3)') params%omega_l
             print*,"*  ", trim(name)," detected : ", params%omega_l 
          endif
          !
          if(trim(name)=="omg_m".or.trim(name)=="OMG_M") then
             read(value,'(F10.3)') params%omega_m
             print*,"*  ", trim(name)," detected : ", params%omega_m 
          endif
          !
          if(trim(name)=="omg_b".or.trim(name)=="OMG_B") then
             read(value,'(F5.3)') params%omega_b
             print*,"*  ", trim(name)," detected : ", params%omega_b 
          endif
          !
          if(trim(name)=="wa".or.trim(name)=="WA") then
             read(value,'(F5.3)') params%wa
             print*,"*     ", trim(name)," detected : ", params%wa 
          endif
          !
          if(trim(name)=="convfc".or.trim(name)=="CONVFC") then
             read(value,'(F11.4)') params%convfc
             print*,"*   ", trim(name)," detected : ", params%convfc 
          endif
          !
          if(trim(name)=="w0".or.trim(name)=="W0") then
             read(value,'(F5.3)') params%w0
             print*,"*     ", trim(name)," detected : ", params%w0 
          endif
          !
          if(trim(name)=="nmax".or.trim(name)=="NMAX") then
             read(value,'(I4)') params%nnmax
             print*,"*   ", trim(name)," detected : ", params%nnmax 
          endif
          !
          if(trim(name)=="lmax".or.trim(name)=="LMAX") then
             read(value,'(I4)') params%nlmax
             print*,"*   ", trim(name)," detected : ", params%nlmax   
          endif
          !
          if(trim(name)=="mmax".or.trim(name)=="MMAX") then
             read(value,'(I4)') params%nmmax
             print*,"*   ", trim(name)," detected : ", params%nmmax   
          endif
          !
          if(trim(name)=="itrord".or.trim(name)=="ITRORD") then
             read(value,'(I2)')  params%iter_order
             print*,"* ", trim(name)," detected : ", params%iter_order
          endif
          !
          if(trim(name)=="multsc".or.trim(name)=="MULTSC") then
             read(value,'(I4)')  params%multiscale
             print*,"* ", trim(name)," detected : ", params%multiscale
          endif
          !
          if(trim(name)=="nbrpts".or.trim(name)=="NBRPTS") then
             read(value,'(I7)')  params%nbptsasked
             print*,"* ", trim(name)," detected : ", params%nbptsasked
          endif
          !
          if(trim(name)=="rmax".or.trim(name)=="RMAX") then
             read(value,'(F11.4)')  params%rmaxinput
             print*,"*   ", trim(name)," detected : ", params%rmaxinput
          endif
          !
       enddo
       print*,"-- done."
       print*,"	"
51     continue
       close(1)
    else
       description = concatnl("Enter input file (galaxy survey) name (eg: data.out) : ")
       params%surveyfile = parse_string(handle, 'surveyfile', default='', descr=description)
       inquire( file=params%surveyfile, exist=fileexists)
       if( fileexists .eqv. .false. ) call fatal_error("This file doesn't exist!")
       if (trim(params%surveyfile)=='') call fatal_error('Error : no input file provided')    
       print*,"	"
       print*,"	"
    endif


  END SUBROUTINE txtfile2parameters


  !===============================================================================
  !===============================================================================


  subroutine f90ftpcld(unit, colnum, frow, felem, np, data, status)
    integer(I4B), intent(in)  :: unit, colnum, frow, felem, np
    integer(I4B), intent(out) :: status
    real(DP),     intent(in), dimension(0:)  :: data
    call ftpcld(unit, colnum, frow, felem, np, data, status)
    return
  end subroutine f90ftpcld

  subroutine f90ftgcvd(unit, colnum, frow, felem, np, nullval, data, anynull, status)
    integer(I4B), intent(in)  :: unit, colnum, frow, felem, np
    integer(I4B), intent(out) :: status
    logical(LGT), intent(out) :: anynull
    real(DP), intent(out), dimension(0:) :: data
    real(DP), intent(in) :: nullval
    call ftgcvd(unit, colnum, frow, felem, np, nullval, data, anynull, status)
    return
  end subroutine f90ftgcvd

  subroutine f90ftgkyd(unit, keyword, value, comment, status)
    integer(I4B),     intent(in)  :: unit
    character(len=*), intent(in)  :: keyword
    integer(I4B),     intent(out) :: status
    character(len=*), intent(out) :: comment
    real(DP),         intent(out) :: value
    call ftgkyd(unit, keyword, value, comment, status)
    return
  end subroutine f90ftgkyd


  !===============================================================================
  !===============================================================================


  subroutine get_clean_header(unit, header, filename, error, xalso, xonly)

    INTEGER(I4B),                    intent(IN)           :: unit
    CHARACTER(LEN=*), DIMENSION(1:), INTENT(IN OUT)       :: header
    CHARACTER(LEN=*),                INTENT(IN)           :: filename
    INTEGER(I4B),                    intent(OUT)          :: error
    character(len=8), dimension(1:), intent(IN), optional :: xalso
    character(len=8), dimension(1:), intent(IN), optional :: xonly

    INTEGER(I4B) :: nlheader, status, i, n_excl
    CHARACTER(LEN=80) :: record
    CHARACTER(len=8), dimension(:), allocatable :: to_excl

    CHARACTER(len=8), dimension(1:21) :: def_excl

    def_excl=(/&
         & "SIMPLE  ","BITPIX  ","NAXIS   ",&
         & "NAXIS#  ","PCOUNT  ","GCOUNT  ",&
         & "EXTEND  ","ORIGIN  ","DATE*   ",&
         & "TFIELDS ","TFORM#  ",           & 
         & "TBCOL#  ","EXTNAME ","CTYPE#  ",&
         & "CRVAL#  ","CRPIX#  ","CDELT#  ",&
         & "XTENSION","INSTRUME","TELESCOP",&
         & "PDMTYPE "/)

    error = 0

    if (present(xonly)) then 
       n_excl = size(xonly)
       allocate(to_excl(1:n_excl))
       to_excl = xonly

    else if (present(xalso)) then
       n_excl = size(xalso) + size(def_excl)
       allocate(to_excl(1:n_excl))
       to_excl(1:size(def_excl)) = def_excl
       to_excl(size(def_excl)+1:n_excl) = xalso

    else
       n_excl = size(def_excl)
       allocate(to_excl(1:n_excl))
       to_excl = def_excl
    endif

    nlheader=size(header)
    ! go to end of fortran header
    do i = 1, nlheader
       if (trim(header(i)) == "") exit
    enddo
    ! go to top of fits file header
    status=0
    call ftgrec(unit,0_i4b,record,status)
    ! read in all header lines except those excluded
    do
       call ftgnxk(unit,'*',1_i4b,to_excl,n_excl,record,status)
       if (status > 0) exit ! end of header
       if (i > nlheader) then
          write(unit=*,fmt="(a,i5,a)") &
               & " WARNING : The header in "//  &
               &    trim(filename)//" has more than ", &
               &  nlheader," lines."
          print*," It will be truncated."
          error = 1
          exit
       endif
       header(i)=record
       i=i+1
    enddo
    status=0

    return
  end subroutine get_clean_header


  !===============================================================================
  !===============================================================================


  SUBROUTINE cln2fits( filename, cln, kln, nlmax, nnmax, header, nlheader )
    !
    INTEGER(I4B), INTENT(IN)      ::  nlmax, nnmax, nlheader
    REAL(DP), INTENT(IN)          ::  cln(0:nlmax,1:nnmax)
    REAL(DP), INTENT(IN)          ::  kln(0:nlmax,1:nnmax)
    CHARACTER(LEN=80), INTENT(IN), DIMENSION(1:nlheader) :: header
    CHARACTER(LEN=*),  INTENT(IN)           :: filename

    REAL(DP), DIMENSION(:), allocatable :: lindex,nindex, klntemp, clntemp
    INTEGER(I4B) ::  status,unit,blocksize,bitpix,naxis,naxes(1)
    INTEGER(I4B) ::  i,hdutype, lmax, mmax, lmin
    LOGICAL(LGT) ::  simple,extend
    CHARACTER(LEN=80) :: comment

    INTEGER(I4B), PARAMETER :: maxdim = 40 !number of columns in the extension
    INTEGER(I4B) :: nrows, npix, tfields, varidat, repeat
    INTEGER(I4B) :: frow,  felem, colnum, stride, istart, iend, k
    CHARACTER(LEN=20) :: ttype(maxdim), tform(maxdim), tunit(maxdim), extname
    CHARACTER(LEN=10) ::  card
    CHARACTER(LEN=2) :: stn
    INTEGER(I4B) :: itn, cnt, ncl
    integer(I4B) :: l, n
    character(len=filenamelen) sfilename
    character(len=1) :: pform

    status = 0
    unit = 140

    call ftinit(unit,filename,blocksize,status)

    simple=.true.
    bitpix=32     ! integer*4
    naxis=0       ! no image
    naxes(1)=0
    extend=.true. ! there is an extension

    ncl = 4

    ! primary header
    call ftphpr(unit,simple,bitpix,naxis,naxes,0_i4b,1_i4b,extend,status)

    !     creates an extension
    call ftcrhd(unit, status)
    !     writes required keywords
    nrows    = nnmax * nlmax  ! naxis1
    tfields  = ncl
    repeat   = 1
    do i=1, ncl
       ttype(i) = ''
       tform(i) = ''
       tform(i)(1:2) = '1D'
    enddo

    ttype(1)(1:9) = 'l index  '
    ttype(2)(1:9) = 'n index  '
    ttype(3)(1:9) = 'kln value'
    ttype(4)(1:9) = 'cln(l,n) '
    tunit(1:ncl) = ''      ! optional, will not appear
    extname  = ''      ! optional, will not appear
    varidat  = 0
    call ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
         &     extname, varidat, status)
    !     write the header literally, putting TFORMi at the desired place
    do i=1,nlheader
       card = header(i)
       if (card(1:5) == 'TTYPE') then ! if TTYPEi is explicitely given
          stn = card(6:7)
          read(stn,'(i2)') itn
          ! discard at their original location:
          call ftdkey(unit,'TTYPE'//stn,status)  ! old TTYPEi and  ! remove
          status = 0
          call ftdkey(unit,'TFORM'//stn,status)  !     TFORMi
          status = 0
          if (itn <= tfields) then ! only put relevant information 2008-08-27
             call putrec(unit,header(i), status)           ! write new TTYPE1
             status = 0
             comment = ''
             if (itn==1) then
                comment = 'data format of field: 4-byte INTEGER'
             else
                if (DP == SP) comment = 'data format of field: 4-byte REAL'
                if (DP == DP) comment = 'data format of field: 8-byte REAL'
             endif
             call ftpkys(unit,'TFORM'//stn,tform(itn),comment,status) ! and write new TFORM1 right after
          endif
       elseif (header(i)/=' ') then
          call putrec(unit,header(i), status)
       endif
       status = 0
    enddo
    call ftukyj(unit, 'MAX-LPOL', nlmax, 'Maximum L multipole order',  status)
    call ftukyj(unit, 'MAX-NPOL', nnmax, 'Maximum N multipole order', status)
    !     write the extension by blocks of rows ! EH, Dec 2004
    felem  = 1  ! starting position (element)
    call ftgrsz(unit, stride, status) ! find optimal stride in rows
    stride = max( stride, 1)
    allocate(lindex(0:stride-1))
    allocate(nindex(0:stride-1))
    allocate(klntemp(0:stride-1))
    allocate(clntemp(0:stride-1))
    lindex = 0.0
    nindex = 0.0
    klntemp = 0.0
    clntemp = 0.0

    cnt = 0
    istart = 1
    do l = 0, nlmax
       do n = 1, nnmax
          if( cnt < stride ) then
             lindex(cnt) = real(l)
             nindex(cnt) = real(n)
             klntemp(cnt) = kln(l,n)
             clntemp(cnt) = cln(l,n)
          endif
          if( cnt >= stride .or. cnt >= nlmax*nnmax ) then
             !print*,lindex,nindex
             call f90ftpcld(unit, 1_i4b, istart, felem, cnt, lindex, status)
             call f90ftpcld(unit, 2_i4b, istart, felem, cnt, nindex, status)
             call f90ftpcld(unit, 3_i4b, istart, felem, cnt, klntemp, status)
             call f90ftpcld(unit, 4_i4b, istart, felem, cnt, clntemp, status)
             istart = istart + cnt - 1
             cnt = 1
             lindex = 0.0
             nindex = 0.0
             klntemp = 0.0
             clntemp = 0.0
          endif
          cnt = cnt + 1
       enddo
    enddo

    deallocate(lindex,nindex,klntemp,clntemp)

    !     close the file and free the unit number
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    RETURN

  END SUBROUTINE cln2fits


  !===============================================================================
  !===============================================================================


  SUBROUTINE almn2fits( filename, almn, kln, cln, nlmax, nmmax, nnmax, header, nlheader )
    !
    INTEGER(I4B), INTENT(IN)      ::  nlmax, nmmax, nnmax, nlheader
    COMPLEX(DP), INTENT(IN)          ::  almn(1:nnmax,0:nlmax,0:nmmax)
    REAL(DP), INTENT(IN)          ::  kln(0:nlmax,1:nnmax)
    REAL(DP), INTENT(IN)          ::  cln(0:nlmax,1:nnmax)
    CHARACTER(LEN=80), INTENT(IN), DIMENSION(1:nlheader) :: header
    CHARACTER(LEN=*),  INTENT(IN)           :: filename

    REAL(DP), DIMENSION(:), allocatable :: lindex,mindex,nindex, klntemp, clntemp
    REAL(DP), DIMENSION(:), allocatable :: almntemp1, almntemp2
    INTEGER(I4B) ::  status,unit,blocksize,bitpix,naxis,naxes(1)
    INTEGER(I4B) ::  i,hdutype, lmax, mmax, lmin
    LOGICAL(LGT) ::  simple,extend
    CHARACTER(LEN=80) :: comment

    INTEGER(I4B), PARAMETER :: maxdim = 40 !number of columns in the extension
    INTEGER(I4B) :: nrows, npix, tfields, varidat, repeat, frow
    INTEGER(I4B) :: felem, colnum, stride, istart, iend, k
    CHARACTER(LEN=20) :: ttype(maxdim), tform(maxdim), tunit(maxdim), extname
    CHARACTER(LEN=10) ::  card
    CHARACTER(LEN=2) :: stn
    INTEGER(I4B) :: cnt, ncl
    integer(I4B) :: l, n, m, itn
    character(len=filenamelen) sfilename
    character(len=1) :: pform

    status = 0
    unit = 139

    call ftinit(unit,filename,blocksize,status)

    !     -----------------------------------------------------
    !     initialize parameters about the FITS image
    simple=.true.
    bitpix=32     ! integer*4
    naxis=0       ! no image
    naxes(1)=0
    extend=.true. ! there is an extension

    ncl = 7

    ! primary header
    call ftphpr(unit,simple,bitpix,naxis,naxes,0_i4b,1_i4b,extend,status)

    !     creates an extension
    call ftcrhd(unit, status)
    !     writes required keywords
    nrows    = nnmax * (nlmax+1) * (1+nmmax/2)  ! naxis1
    tfields  = ncl
    repeat   = 1
    do i=1, ncl
       ttype(i) = ''
       tform(i) = ''
       tform(i)(1:2) = '1D'
    enddo
    ttype(1)(1:9) = 'l index  '
    ttype(2)(1:9) = 'm index  '
    ttype(3)(1:9) = 'n index  '
    ttype(4)(1:9) = 'kln value'
    ttype(5)(1:9) = 'cln value'
    ttype(6)(1:15) = 'Real(almn(l,n))'
    ttype(7)(1:15) = 'Img(almn(l,n)) '
    tunit(1:ncl) = ''      ! optional, will not appear
    extname  = ''      ! optional, will not appear
    varidat  = 0
    call ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
	 &     extname, varidat, status)
    !     write the header literally, putting TFORMi at the desired place
    do i=1,nlheader
       card = header(i)
       if (card(1:5) == 'TTYPE') then ! if TTYPEi is explicitely given
	  stn = card(6:7)
	  read(stn,'(i2)') itn
   ! discard at their original location:
	  call ftdkey(unit,'TTYPE'//stn,status)  ! old TTYPEi and  ! remove
	  status = 0
	  call ftdkey(unit,'TFORM'//stn,status)  !     TFORMi
	  status = 0
	  if (itn <= tfields) then ! only put relevant information 2008-08-27
	     call putrec(unit,header(i), status)           ! write new TTYPE1
	     status = 0
	     comment = ''
	     if (itn==1) then
		comment = 'data format of field: 4-byte INTEGER'
	     else
		if (DP == SP) comment = 'data format of field: 4-byte REAL'
		if (DP == DP) comment = 'data format of field: 8-byte REAL'
	     endif
	     call ftpkys(unit,'TFORM'//stn,tform(itn),comment,status) ! and write new TFORM1 right after
	  endif
       elseif (header(i)/=' ') then
	  call putrec(unit,header(i), status)
       endif
       status = 0
    enddo
    call ftukyj(unit, 'MAX-LPOL', nlmax, 'Maximum L multipole order', status)
    call ftukyj(unit, 'MAX-MPOL', nmmax, 'Maximum M multipole order', status)
    call ftukyj(unit, 'MAX-NPOL', nnmax, 'Maximum N multipole order', status)
    !     write the extension by blocks of rows ! EH, Dec 2004
    felem  = 1  ! starting position (element)
    call ftgrsz(unit, stride, status) ! find optimal stride in rows
    stride = max( stride, 1)

    allocate(lindex(0:stride-1))
    allocate(mindex(0:stride-1))
    allocate(nindex(0:stride-1))
    allocate(klntemp(0:stride-1))
    allocate(clntemp(0:stride-1))
    allocate(almntemp1(0:stride-1))
    allocate(almntemp2(0:stride-1))

    lindex = 0.0
    mindex = 0.0
    nindex = 0.0
    klntemp = 0.0
    clntemp = 0.0
    almntemp1 = 0.0
    almntemp2 = 0.0

    cnt = 0
    istart = 1
    do l = 0, nlmax
       do m = 0, l
          do n = 1, nnmax
             lindex(cnt) = real(l,kind=DP)
             mindex(cnt) = real(m,kind=DP)
             nindex(cnt) = real(n,kind=DP)
             klntemp(cnt) = real(kln(l,n),kind=DP)
             clntemp(cnt) = real(cln(l,n),kind=DP)
             almntemp1(cnt) = real(almn(n,l,m))
             almntemp2(cnt) = aimag(almn(n,l,m))
             cnt = cnt + 1
             !print*,l,m,n,almn(n,l,m)
             if( cnt >= stride .or. cnt == nrows ) then
                call f90ftpcld(unit, 1_i4b, istart, felem, cnt, lindex, status)
                call f90ftpcld(unit, 2_i4b, istart, felem, cnt, mindex, status)
                call f90ftpcld(unit, 3_i4b, istart, felem, cnt, nindex, status)
                call f90ftpcld(unit, 4_i4b, istart, felem, cnt, klntemp, status)
                call f90ftpcld(unit, 5_i4b, istart, felem, cnt, clntemp, status)
                call f90ftpcld(unit, 6_i4b, istart, felem, cnt, almntemp1, status)
                call f90ftpcld(unit, 7_i4b, istart, felem, cnt, almntemp2, status)
                istart = istart + cnt - 1
                cnt = 0
                lindex = 0.0
                mindex = 0.0
                nindex = 0.0
                klntemp = 0.0
                clntemp = 0.0
                almntemp1 = 0.0
                almntemp2 = 0.0
             endif
          enddo
       enddo
    enddo

    deallocate(lindex,mindex,nindex,klntemp,clntemp,almntemp1,almntemp2)

    !     ----------------------
    !     close and exit
    !     ----------------------

    !     close the file and free the unit number
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    RETURN

  END SUBROUTINE almn2fits


  !===============================================================================
  !===============================================================================


  SUBROUTINE fits2almn(filename, nnmax, nlmax, nmmax, almn, kln, cln, header, nlheader ) 

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER(I4B),  INTENT(IN) :: nnmax,nlmax,nmmax,nlheader
    REAL(DP), INTENT(OUT), DIMENSION(0:nlmax,1:nnmax)  ::  kln
    REAL(DP), INTENT(OUT), DIMENSION(0:nlmax,1:nnmax) ::  cln
    COMPLEX(DP), DIMENSION(1:nnmax,0:nlmax,0:nmmax), INTENT(OUT) :: almn
    CHARACTER(LEN=80), INTENT(OUT), DIMENSION(1:nlheader) :: header
    REAL(DP)                                        :: nullval
    LOGICAL(LGT)                                    ::  anynull

    REAL(DP), DIMENSION(:), allocatable :: lindex,mindex,nindex, klntemp, clntemp
    REAL(DP), DIMENSION(:), allocatable :: almntemp1, almntemp2

    REAL(DP), DIMENSION(:,:), allocatable :: tab

    INTEGER(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis, extno
    INTEGER(I4B) :: npix, ncl
    CHARACTER(LEN=80) :: comment ! , record
    LOGICAL(LGT) :: extend
    INTEGER(I4B) :: nmove, hdutype ! , nkeys , nspace
    INTEGER(I4B) :: frow, imap
    INTEGER(I4B) :: datacode, repeat, width
    integer(I4B) :: i, l, m, n
    integer(i4b) :: nrow2read, nelem
    integer(i4b) :: i0, i1

    INTEGER(I4B), PARAMETER :: maxdim = 40 !number of columns in the extension
    INTEGER(I4B) :: nrows, tfields, varidat
    CHARACTER(LEN=20) :: ttype(maxdim), tform(maxdim), tunit(maxdim), extname
    CHARACTER(LEN=20) :: nmtemp

    !-----------------------------------------------------------------------
    status=0
    ncl = 7
    header=''
    extno = 1
    unit = 145
    naxes(1) = 1
    naxes(2) = 1
    nfound = -1
    anynull = .false.
    almn=0.
    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
    if (status > 0) then 
       call printerror(status)
       call fatal_error("Aborting.")
    endif
    !     -----------------------------------------

    !     determines the presence of image
    call ftgkyj(unit,'NAXIS', naxis, comment, status)
    if (status > 0) call printerror(status)

    !     determines the presence of an extension
    call ftgkyl(unit,'EXTEND', extend, comment, status)
    if (status > 0) status = 0 ! no extension : 
    !     to be compatible with first version of the code

    call assert (extend, 'No extension!')
    nmove = +extno
    call ftmrhd(unit, nmove, hdutype, status)
    !cc         write(*,*) hdutype

    call assert(hdutype==2, 'this is not a binary table')

    header = ""
    call get_clean_header( unit, header, filename, status)

    !        reads all the keywords
    call ftghbn(unit, maxdim, &
	 &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
	 &        status)

    if (tfields<ncl) then
       print *,'found ',tfields,' columns in the file'
       print *,'expected ',ncl
       call fatal_error
    endif
    call f90ftgkyd(unit, 'BAD_DATA', nullval, comment, status)
    if (status == 202) then ! bad_data not found
       if (DP == SP) nullval = s_bad_value ! default value
       if (DP == DP) nullval = d_bad_value ! default value
       status = 0
    endif
    !parse TFORM keyword to find out the length of the column vector
    call ftbnfm(tform(1), datacode, repeat, width, status)
    npix = nrows * repeat
    !    if (npix /= nalms) then
    !if (npix > nalms) then
    !print *,'found ',npix,' alms'
    !!       print *,'expected ',nalms
    !print *,'expected ',nalms,' or less'
    !call fatal_error
    !endif

    call ftgrsz(unit, nrow2read, status)
    nrow2read = max(nrow2read, 1)
    nelem = nrow2read * repeat
    i0 = 0_i8b

    allocate(tab(0:nelem-1,1:ncl))

    allocate(lindex(0:nelem-1))
    allocate(mindex(0:nelem-1))
    allocate(nindex(0:nelem-1))
    allocate(klntemp(0:nelem-1))
    allocate(clntemp(0:nelem-1))
    allocate(almntemp1(0:nelem-1))
    allocate(almntemp2(0:nelem-1))

    lindex = 0.0
    mindex = 0.0
    nindex = 0.0
    klntemp = 0.0
    clntemp = 0.0
    almntemp1 = 0.0
    almntemp2 = 0.0

    do frow = 1, nrows, nrow2read
       i1 = min(i0 + nrow2read * repeat, int(npix,i8b)) - 1_i8b
       nelem = i1 - i0 + 1
       !do imap = 1, ncl
       !   call f90ftgcvd(unit, imap, frow, 1_i4b, nelem, nullval, &
       !   &        tab(0:nelem-1,imap), anynull, status)
       !   call assert (.not. anynull, 'There are undefined values in the table!')
       !enddo

       call f90ftgcvd(unit, 1_i4b, frow, 1_i4b, nelem, nullval, &
            &        lindex(0:nelem-1), anynull, status)
       status = 0
       call f90ftgcvd(unit, 2_i4b, frow, 1_i4b, nelem, nullval, &
            &        mindex(0:nelem-1), anynull, status)
       status = 0
       call f90ftgcvd(unit, 3_i4b, frow, 1_i4b, nelem, nullval, &
            &        nindex(0:nelem-1), anynull, status)
       status = 0
       call f90ftgcvd(unit, 4_i4b, frow, 1_i4b, nelem, nullval, &
            &        klntemp(0:nelem-1), anynull, status)
       status = 0
       call f90ftgcvd(unit, 5_i4b, frow, 1_i4b, nelem, nullval, &
            &        clntemp(0:nelem-1), anynull, status)
       status = 0
       call f90ftgcvd(unit, 6_i4b, frow, 1_i4b, nelem, nullval, &
            &        almntemp1(0:nelem-1), anynull, status)
       status = 0
       call f90ftgcvd(unit, 7_i4b, frow, 1_i4b, nelem, nullval, &
            &        almntemp2(0:nelem-1), anynull, status) 
       status = 0
       call assert (.not. anynull, 'There are undefined values in the table!')

       do i = 0, nelem-1
          !print*,"entering writing"
          l = int(lindex(i),i4b)
          m = int(mindex(i),i4b)
          n = int(nindex(i),i4b)
          if( l >= 0 .and. l <= nlmax &
               & .and. m >= 0 .and. m <= l &
               & .and. n >= 1 .and. n <= nnmax ) then
             kln(l,n) = klntemp(i)
             cln(l,n) = clntemp(i)
             almn(n,l,m) = CMPLX( almntemp1(i), almntemp2(i), KIND=DP)
110          FORMAT(A,I3,A,I3,A,I3,A) 
             !WRITE(nmtemp,110),"(",l,",",m,",",n,")"
             !print*,nmtemp,almn(n,l,m)
          endif
       enddo

       i0 = i1 + 1_i8b
       !print*,i0,npix
    enddo

    deallocate(lindex,mindex,nindex,klntemp,clntemp,almntemp1,almntemp2)
    deallocate(tab)

    ! sanity check
    if (i0 /= npix) then
       print*,'something wrong during piece-wise reading'
       call fatal_error
    endif

    !     close the file
    call ftclos(unit, status)


    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)
    return

  END SUBROUTINE fits2almn


  !===============================================================================
  !===============================================================================


  SUBROUTINE print_almn(almn, nllim, nmlim, nnlim, nlmax, nmmax, nnmax, txt)

    INTEGER(I4B) :: nllim, nmlim, nnlim, l, m, n, nnmax, nlmax, nmmax
    COMPLEX(DPC), intent(IN), DIMENSION(1:nnmax,0:nlmax,0:nmmax) :: almn
    CHARACTER(len=20) :: nm
    CHARACTER(len=4) :: txt

    PRINT*,">> ",txt,"(n,l,m)"
110 FORMAT (A,A,I3,A,I3,A,I3,A)    	
    DO n = 1, nnlim
       DO l = 0, nllim
          DO m = 0, l
             WRITE(nm,110) trim(txt),'(', n,",",l,",", m,")"
             print*,nm, almn(n,l,m)
          ENDDO
       ENDDO
    ENDDO


  END SUBROUTINE print_almn

  !===============================================================================
  !===============================================================================

  SUBROUTINE print_spectrum(spectr, nllim, nnlim, nlmax, nnmax, txt)

    INTEGER(I4B) :: nllim, nnlim, l, n, nnmax, nlmax
    REAL(DP), intent(IN), DIMENSION(0:nlmax,1:nnmax) :: spectr
    CHARACTER(len=20) :: nm
    CHARACTER(len=4) :: txt

    PRINT*,">> ",txt,"(l,n)"
111 FORMAT (A,A,I3,A,I3,A)  
    DO n = 1, nnlim
       DO l = 0, nllim
          WRITE(nm,111) trim(txt),'(', l,",", n,")"
          PRINT*,nm, spectr(l,n)
       ENDDO
    ENDDO


  END SUBROUTINE print_spectrum

  !===============================================================================
  !===============================================================================



  SUBROUTINE almn2cln(nnmax, nlmax, nmmax, almn, cln)

    INTEGER(I4B),                      intent(in) :: nlmax, nmmax, nnmax
    COMPLEX(DPC), dimension(1:nnmax,0:nlmax,0:nmmax), intent(in) :: almn
    REAL(DP)    , dimension(0:nlmax, 1:nnmax),  intent(out):: cln
    ! 
    INTEGER(I4B) :: n, l, ncl, na, mm
    COMPLEX(DPC)     :: dc
    REAL(DP), PARAMETER :: two = 2.000000000000000000_dp
    REAL(DP), PARAMETER :: one = 1.000000000000000000_dp

    ncl = size(cln, 2)
    na = size(almn, 1)
    cln = 0.0_DP

    DO n = 1, nnmax
       DO l = 0, nlmax
          mm = min(l, nmmax)
          dc = sum(almn(n,l,1:mm)*CONJG(almn(n,l,1:mm)))
          dc = (dc + CONJG(dc)) + almn(n,l,0)*almn(n,l,0)
          !print*,l,n,dc,two,one
          cln(l,n) = abs(REAL(dc, KIND=DP)) / abs(two*l + one)

       ENDDO
    ENDDO

    RETURN

  END SUBROUTINE almn2cln


  !===============================================================================
  !===============================================================================

  SUBROUTINE getParameters_survey2almn(multiscale, zbounds, nr, nside, nnmax, nlmax, nmmax, iter_order, &
       & nbpts, iwhere, inradian, sorted, convfc, nbcolsinfile, surveyfile, paramfile, &
       & almnoutfile, clnoutfile, h, w0, wa, omega_l, omega_b, omega_m, rmaxinput )

    ! ------------------------------------------------------------------------------
    ! Get all the parameters of the method
    ! ------------------------------------------------------------------------------

    IMPLICIT NONE

    integer(i4b)       :: nbparam, rdstatus, i, k, status
    logical(LGT)       :: fileexists
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
    integer(i4b)      		    	:: nbptsasked, nbptsinfile, nbcolsinfile
    integer(i4b)  		      	    :: inradian, sorted, iter_order, nbpts, nr
    integer(i4b)  		            :: nside, nnmax, nlmax, nmmax, multiscale
    character    		            :: coordsys
    integer(kind=i4b), DIMENSION(1:3)  :: iwhere
    real(kind=DP)		            :: theta_cut_deg, cos_theta_cut, fsky, convfc
    real(kind=DP)		    :: rmaxinput, h, w0, wa, omega_l, omega_b, omega_m
    real(kind=DP),     DIMENSION(1:2)  :: zbounds

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
    if( fileexists .eqv. .false. ) then
       description = concatnl("Enter input file name (eg: data.out) : ")
       surveyfile = parse_string(handle, 'surveyfile', default='', descr=description)
       inquire( file=surveyfile, exist=fileexists)
       if( fileexists .eqv. .false. ) call fatal_error("This file doesn't exist!")
       if (trim(surveyfile)=='') call fatal_error('Error : no input file provided')
       print*,"	"
       print*,"	"
    endif
    !
    if( almnoutfile /= "no" ) then
       almnoutfile = "out/"//almnoutfile
       inquire( file=almnoutfile, exist=fileexists)
       if( fileexists ) then
          print*,"Error : almn outputfile already exists ; please provide another name"
          description = concatnl("Enter input file name (eg: almn.fits) : ")
          almnoutfile = parse_string(handle, 'almnoutfile', default='almn.fits', descr=description)
          almnoutfile = "out/"//almnoutfile
          inquire( file=almnoutfile, exist=fileexists)
          if( fileexists .eqv. .true. ) call fatal_error("This file already exists!")
          print*,"	"
       endif
    endif
    !
    if( clnoutfile /= "no" ) then
       clnoutfile = "out/"//clnoutfile
       inquire( file=clnoutfile, exist=fileexists)
       if( fileexists ) then
          print*,"Error : cln outputfile already exists ; please provide another name"
          description = concatnl("Enter input file name (eg: cln.fits) : ")
          clnoutfile = parse_string(handle, 'clnoutfile', default='cln.fits', descr=description)
          clnoutfile = "out/"//clnoutfile
          inquire( file=clnoutfile, exist=fileexists)
          if( fileexists .eqv. .true. ) call fatal_error("This file already exists!")
          print*,"	"
       endif
    endif
    !
    if( nbcolsinfile < 0 ) then
       description = concatnl("How many columns does this file contain? (Necessary for parsing)")
       nbcolsinfile = parse_int(handle, 'nbcolsinfile', default=-1, descr=description)
       if( nbcolsinfile < 0 .or. nbcolsinfile > 20 ) call fatal_error("Error : Bad number.")
       print*,"	"
       print*,"	"
    endif
    !
    if( iwhere(1)<1 .or. iwhere(1)>nbcolsinfile ) then
       WRITE(chline,"(a)") "In which column is Coord 1 ?"   
       description = concatnl( chline , "" , chline1 )
       iwhere(1) = parse_int(handle, 'Coord 1', default=-1, descr=description)
       if( iwhere(1)<1 .or. iwhere(1)>nbcolsinfile ) call fatal_error("Incorrect column number")
       print*,"	"
       print*,"	"
    endif
    !
    if( iwhere(2)<1 .or. iwhere(2)>nbcolsinfile ) then
       WRITE(chline,"(a)") "In which column is Coord 2 ?"
       description = concatnl( chline , "" , chline1 )
       iwhere(2) = parse_int(handle, 'Coord 2', default=-1, descr=description)
       if( iwhere(2)<1 .or. iwhere(2)>nbcolsinfile ) call fatal_error("Incorrect column number")
       print*,"	"
       print*,"	"
    endif
    !
    if( iwhere(3)<1 .or. iwhere(3)>nbcolsinfile ) then
       WRITE(chline,"(a)") "In which column is Coord 3 ?"   
       description = concatnl( chline , "" , chline1 )
       iwhere(3) = parse_int(handle, 'Coord 3', default=-1, descr=description)
       if( iwhere(3)<1 .or. iwhere(3)>nbcolsinfile ) call fatal_error("Incorrect column number")
       print*,"	"
       print*,"	"
    endif
    !
    if( inradian < 0 .or. inradian > 1 ) then
       description = concatnl("Are your data in radian?" ,&
	    &"","(0) No, in degree    (1) Yes, in radian")
       inradian = parse_int(handle, 'inradian', default=0, descr=description)
       if( inradian < 0 .or. inradian > 1 ) call fatal_error("Error : Bad number.")
       print*,"	"
       print*,"	"
    endif
    !
    if( sorted < 0 .or. sorted > 1 ) then
       description = concatnl("Is the array sorted by radial/redshift values? ",&
	    &"(0) No, no sorted    (1) Yes, sorted")
       sorted = parse_int(handle, 'sorted', default=0, descr=description)
       if( sorted < 0 .or. sorted > 1 ) call fatal_error("Error : Bad number.")
       print*,"	"
       print*,"	"
    endif
    !
    if( nside < 0 .or. nside > 4096 ) then
       WRITE(chline,"(a)") "We recommend: (256 <= Nside <= 2048) a power of 2"
       description = concatnl(&
	    & "Enter the Nside parameter (nsmax) for the healpix analysis. ", &
	    & chline )
       nside = parse_int(handle, 'nside', default=64, descr=description)
       if( nside<0 .or. nside>1024 ) call fatal_error("Incorrect boundaries (admitted : 4,8,16,64,128,256,512)")
       print*,"	"
       print*,"	"
    endif
    !
    if( nnmax < 0 .or. nnmax >= 250 ) then
       WRITE(chline,"(a)") "We recommend: (1 <= n <= n_max <= 250)"
       description = concatnl(&
	    & "Enter the maximum n range (n_max) for the analysis (Bessel functions). ", &
	    & chline )
       nnmax = parse_int(handle, 'nnmax', default=9, descr=description)
       if( nnmax<0 .or. nnmax>500 ) call fatal_error("Incorrect boundaries (admitted : [1,250])")
       print*,"	"
       print*,"	"
    endif
    !
    if( nlmax < 0 .or. nlmax > 256 ) then
       WRITE(chline1,"(a,i5)") "The map has Nside = ",nside
       WRITE(chline,"(a,i5,a)") "We recommend: (0 <= l <= l_max <= ",200,")"
       description = concatnl(&
	    & chline1, &
	    & "Enter the maximum l range (l_max) for the analysis. ", &
	    & chline )
       nlmax = parse_int(handle, 'nlmax', default=2*nside, descr=description)
       if( nlmax<0 .or. nlmax>4*nside ) call fatal_error("Incorrect boundaries (admitted : [0,4*nside])")
       print*,"	"
       print*,"	"
    endif
    !
    if( nmmax < 0 .or. nmmax > 256 ) then
       WRITE(chline,"(a,i5,a)") "We recommend: (0 <= m <= m_max <= ",nlmax,")"
       description = concatnl(&
	    & " Enter the maximum m range (m_max) for the analysis. ", &
	    & chline )
       nmmax = parse_int(handle, 'nmmax', default=nlmax, descr=description)
       if( nmmax<0 .or. nmmax>nlmax ) call fatal_error("Incorrect boundaries (admitted : [1,nlmax])")
       print*,"	"
       print*,"	"
    endif
    !
    if( nr < 4 ) then
       description = concatnl(&
	    & " Enter the number of layers (radial discretization) " )
       nr = parse_int(handle, 'nr', default=8, descr=description)
       if( nr<4 ) call fatal_error("Incorrect boundaries (admitted : [4,nlmax])")
       print*,"	"
       print*,"	"
    endif
    !
    if( h < 0 ) then
       description = concatnl(&
	    & " Cosmogoly - parameter : h " )
       h = parse_double(handle, 'h', default=0.7_dp, descr=description)
       if( h<0 ) call fatal_error("Error")
       print*,"	"
       print*,"	"
    endif
    !
    if( omega_b < 0 ) then
       description = concatnl(&
	    & " Cosmogoly - parameter : omega_b " )
       omega_b = parse_double(handle, 'omega_b', default=0.04_dp, descr=description)
       if( omega_b<0 ) call fatal_error("Error")
       print*,"	"
       print*,"	"
    endif
    !
    if( omega_m < 0 ) then
       description = concatnl(&
	    & " Cosmogoly - parameter : omega_m " )
       omega_m = parse_double(handle, 'omega_m', default=0.3_dp, descr=description)
       if( omega_m<0 ) call fatal_error("Error")
       print*,"	"
       print*,"	"
    endif
    !
    if( omega_l < 0 ) then
       description = concatnl(&
	    & " Cosmogoly - parameter : omega_l " )
       omega_l = parse_double(handle, 'omega_l', default=0.7_dp, descr=description)
       if( omega_l<0 ) call fatal_error("Error")
       print*,"	"
       print*,"	"
    endif
    !
    if( wa < -5 ) then
       description = concatnl(&
	    & " Cosmogoly - parameter : wa " )
       wa = parse_double(handle, 'wa', default=0.0_dp, descr=description)
       if( wa<-5 ) call fatal_error("Error")
       print*,"	"
       print*,"	"
    endif
    !
    if( w0 < -5 ) then
       description = concatnl(&
	    & " Cosmogoly - parameter : w0 " )
       w0 = parse_double(handle, 'w0', default=-1.0_dp, descr=description)
       if( w0<-5 ) call fatal_error("Error")
       print*,"	"
       print*,"	"
    endif
    !
    if( multiscale < 0 .or. multiscale > 4 ) then
       description = concatnl(&
	    & " Enter the number of layers (radial discretization) " )
       multiscale = parse_int(handle, 'multiscalenr', default=1, descr=description)
       if( multiscale<0 .or. multiscale > 4 ) call fatal_error("Incorrect boundaries (admitted : [1,3])")
       print*,"	"
       print*,"	"
    endif
    !
    if( rmaxinput < 0 ) then
       description = concatnl(&
	    & " Maximal radial value - rmax " )
       rmaxinput = parse_double(handle, 'rmaxinput', default=-1.0_dp, descr=description)
       if( rmaxinput<0 ) call fatal_error("Error")
       print*,"	"
       print*,"	"
    endif
    !
    if( iter_order < 0 .or. iter_order > 10 ) then
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
       end select
       print*,"	"
       print*,"	"
    endif
    !
    if( theta_cut_deg < 0.0 .or. theta_cut_deg > 1.0 ) then
       description = concatnl(&
	    & " Enter the symmetric cut around the equator in DEGREES : ", &
	    & " (One ignores data within |b| < b_cut)     0 <= b_cut = ")
       theta_cut_deg = parse_double(handle, 'theta_cut_deg', &
	    &                       vmin=0.0_dp, default=0.0_dp, descr=description)
       cos_theta_cut = SIN(theta_cut_deg/180.d0*PI) !counted from equator instead of from pole
       zbounds = (/ cos_theta_cut , -cos_theta_cut /)
       if (theta_cut_deg<1e-4) zbounds = (/ -1.0_dp, 1.0_dp /) !keep all sphere if cut not set
       fsky = (zbounds(2)-zbounds(1))/2.0_dp
       if (fsky <= 0.0_dp) fsky = 1.0_dp + fsky
       write(*,"(a,f6.1,a)") "One keeps ",100.*fsky," % of the original map(s)"   
       print*,"	"
       print*,"	"
    endif
    !
    if( nbptsasked <= 0 ) then
       description = concatnl("How many galaxies do you want to take into account?")
       nbptsasked = parse_int(handle, 'nbptsasked', default=-1, descr=description)
       if( nbptsasked <= 0 ) call fatal_error("Error : Bad number.")
       print*,"	"
       print*,"	"
    endif

    print*,"-----------------------------------"
    print*,"Starting extraction from file"
    print*,"-----------------------------------"
    print*,""
    print*,"File to be read : ",surveyfile(1:10)
    inquire( file=surveyfile, exist=fileexists)
    if( fileexists ) then
       print*,"Checking file ... ok existing."
    else
       call fatal_error("Error: unknown file")
    endif
    open(2,file= surveyfile, status='old', iostat=rdstatus, form="formatted")
    print*,"Opening file status: ",rdstatus
    nbptsinfile = 0
    do while(rdstatus == 0)  
       read(2,*,iostat = rdstatus) 
       if (rdstatus == 0) then
          nbptsinfile = nbptsinfile + 1
       endif
    end do
    print*,"Number of points read : ", nbptsinfile
    nbpts = min(nbptsasked,nbptsinfile)
    print*,"Final number of galaxies to be analyzed : ",nbpts
    close(2)

  END SUBROUTINE getParameters_survey2almn


  !===============================================================================
  !===============================================================================


  SUBROUTINE extractFromFile( nbpts, nr, iwhere, inradian, sorted, convfc, nbcolsinfile, surveyfile, survey )

    IMPLICIT NONE

    ! ------------------------------------------------------------------------------
    ! Extract data from input file and convert it into valid healpix 3D maps
    ! ------------------------------------------------------------------------------

    character(len=*), PARAMETER    :: code = "Survey extraction"
    character(len=FILENAMELEN)     :: surveyfile
    integer(I4B)		    :: nbcolsinfile, status, nbpts, k, sorted
    integer(I4B)		    :: i, inradian, rdstatus, nsmax, nr
    real(kind=DP), DIMENSION(1:nbpts) 	     :: mapcop
    real(kind=DP), DIMENSION(:), allocatable:: temp
    real(kind=DP), DIMENSION(1:nbpts,1:3)   :: survey, surveybis
    integer(i4b),  DIMENSION(1:3)           :: iwhere
    real(DP) :: convfc

    ALLOCATE(temp(1:nbcolsinfile),stat = status)
    call assert_alloc(status,code,"temp")

    survey = 0.0
    surveybis = 0.0
    open(2,file= trim(surveyfile), status='old', iostat=rdstatus, form="formatted")
    rewind(unit=2,iostat=rdstatus)
    print*,""
    print*,"Extraction..."
    do while(rdstatus == 0)
       do i=1, nbpts
          READ (2,*,iostat = rdstatus, end=61) (temp(k),k=1,nbcolsinfile) 
          surveybis(i,1) = temp(iwhere(1))
          surveybis(i,2) = temp(iwhere(2))
          surveybis(i,3) = temp(iwhere(3))
       enddo
    enddo
61  continue
    print*,"-- done"
    close(2)
    deallocate(temp)

    ! Convertion from degree to radian if necessary
    if( inradian == 0 ) then
       print*,""
       print*,"Convertion from degree to radian..."
       do i=1,nbpts
          surveybis(i,1) = surveybis(i,1)*PI/180
          surveybis(i,2) = surveybis(i,2)*PI/180 
          surveybis(i,2) = HALFPI - surveybis(i,2)            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       enddo
       print*,"-- done"
    endif

    ! ------------------------------------------------------------------
    !!		Sorting the survey
    ! ------------------------------------------------------------------ 

    if( sorted == 0 ) then
       print*,"Sorting values."
       mapcop(1:nbpts) = surveybis(1:nbpts,3)
       call BSORT (mapcop, surveybis, nbpts)
       do i = 1, nbpts
          survey(i,1:3) = surveybis(nbpts-i+1,1:3)
          survey(i,3) = survey(i,3) / convfc
       enddo
    else
       survey = surveybis
       survey(:,3) = survey(:,3) / convfc
    endif

    RETURN

  END SUBROUTINE extractFromFile


!!$    !===============================================================================
!!$    !===============================================================================
!!$    
!!$!> Multiresolution implementation of Fourier-Bessel decomposition. 
!!$!! DO NOT USE, USE survey2almn_srs instead.
!!$    SUBROUTINE survey2almn_mrs( multiscale, nsmax, nnmax, nlmax, nmmax, rmax, nbpts, &
!!$   			     & zbounds, survey, kln, almn)
!!$  
!!$    integer(I4B)  :: status, nbpts, ipring, nrr, ipring2, ipring4
!!$    integer(I4B)  :: nsmax, nnmax, nlmax, nmmax, n, p, multiscale, m, l, nsmax2, nsmax4
!!$    real(DP)	  :: rmax, tempval, order
!!$    !	
!!$    real(DP),    dimension(:,:), allocatable	     :: map,map2
!!$    real(DP), dimension(1:nbpts,1:3)         	     :: survey
!!$    real(DP),    dimension(1:2)       	             :: zbounds
!!$    real(DP),    dimension(0:nlmax,1:nnmax)         :: kln
!!$    real(DP),    dimension(:), allocatable        :: cl
!!$    !real(DP),    dimension(:,:,:) , allocatable     :: jln
!!$    real(DP)                                        :: jln
!!$    real(DP),    dimension(:,:)   , allocatable     ::  plm
!!$    complex(DP), dimension(1:nnmax,0:nlmax,0:nmmax) :: almn
!!$    complex(kind=DP), dimension(:,:,:), allocatable :: almn2
!!$    character(len=*), PARAMETER                     :: code = "survey2almn_srs"
!!$    !
!!$    almn = cmplx( 0.0, 0.0, kind=DP )
!!$    !
!!$    ALLOCATE(cl(0:nlmax),stat = status)
!!$    call assert_alloc(status,code,"cln")
!!$    !
!!$    cl = 1.0_dp
!!$    !
!!$    ALLOCATE(map(0:nlmax,0:(12*nsmax**2-1)),stat = status)
!!$    CALL assert_alloc(status,code,"map")
!!$    !
!!$    IF( .FALSE. ) THEN ! multiscale == 1 : TODO
!!$       !
!!$       CALL survey2almn_srs( nsmax, nnmax, nlmax, nmmax, rmax, nbpts, &
!!$   			     & zbounds, survey, kln, almn)
!!$       !
!!$    ELSE	! TODO !!!!!
!!$       !
!!$       print*,"Entering multiscale procedure..."
!!$       !
!!$       nsmax2 = 2*nsmax
!!$       nsmax4 = 2*nsmax2
!!$       !
!!$       ALLOCATE(map2(0:nlmax,0:(12*(nsmax2)**2-1)),stat = status)
!!$       CALL assert_alloc(status,code,"map2")
!!$       ALLOCATE(almn2(1:nnmax,0:nlmax,0:nmmax),stat = status)
!!$       CALL assert_alloc(status,code,"almn2")
!!$       almn2 = cmplx( 0.0, 0.0, kind=DP )
!!$       !
!!$       DO n = 1, nnmax
!!$          !
!!$          DO p = 1, nbpts
!!$             !
!!$             CALL ang2pix_ring(nsmax, survey(p,2), survey(p,1), ipring)
!!$             CALL ang2pix_ring(nsmax2, survey(p,2), survey(p,1), ipring2)
!!$             CALL ang2pix_ring(nsmax4, survey(p,2), survey(p,1), ipring4)
!!$             ! 	
!!$             DO l = 0, nlmax
!!$                !
!!$                jln = 0.0_dp
!!$                CALL BJL( l , kln(l,n)*survey(p,3) , jln ) 
!!$                !
!!$                IF(         ipring4 == 16*ipring+3 &
!!$                     & .OR. ipring4 == 16*ipring+6 &
!!$                     & .OR. ipring4 == 16*ipring+9 &
!!$                     & .OR. ipring4 == 16*ipring+12 ) THEN
!!$                   map(l,ipring) = map(l,ipring) + kln(l,n) * jln * cl(l)
!!$                ELSE
!!$                   map2(l,ipring2) = map(l,ipring2) + kln(l,n) * jln * cl(l)
!!$                ENDIF
!!$                !
!!$             ENDDO
!!$             !
!!$          ENDDO
!!$          !
!!$          CALL alnspring2almn( nsmax, nlmax, nmmax, map, &
!!$               & almn(n:n, 0:nlmax,0:nmmax), zbounds )
!!$          !
!!$          CALL alnspring2almn( nsmax2, nlmax, nmmax, map2, &
!!$               & almn2(n:n, 0:nlmax,0:nmmax), zbounds )
!!$          !
!!$       ENDDO
!!$       print*,"----------------------------------------"
!!$       !
!!$       almn = almn + almn2 !almn*12.0*(real(nsmax)**2.0) + almn2*12.0*(real(nsmax2)**2.0)
!!$       !
!!$       DEALLOCATE( map, map2, almn2 )
!!$       !
!!$    ENDIF
!!$    !
!!$    RETURN
!!$    !
!!$    END subroutine survey2almn_mrs
!!$    
!!$    !===============================================================================
!!$    !===============================================================================
!!$
!!$    
!!$    SUBROUTINE almn2alnspring_pre( nsmax, nlmax, nmmax, map, almn, plm )
!!$
!!$      integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
!!$      complex(DP), intent(IN),  dimension(1:1,0:nlmax,0:nmmax) :: almn
!!$      real(DP),   intent(OUT), dimension(0:nlmax,0:(12_i8b*nsmax)*nsmax-1) :: map
!!$      real(DP),     intent(IN),  dimension(0:)                  :: plm
!!$
!!$      integer(I4B) :: l, m, ith
!!$      integer(I8B) :: istart_south, istart_north, npix
!!$      integer(I4B) :: nrings, nphmx
!!$      integer(i4b), parameter :: SMAXCHK = 30
!!$
!!$      integer(I8B) :: n_lm, n_plm, i_mm
!!$      complex(DPC), dimension(-1:1)             :: b_ns
!!$      real(DP),     dimension(:,:), allocatable :: dalm
!!$      integer(i4b)                              :: ll, l_min, l_start
!!$      real(DP)                                  :: cth
!!$      complex(DPC), dimension(:,:), allocatable :: b_north, b_south
!!$      real(DP),     dimension(:),   allocatable :: ring
!!$      integer(i4b)             :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
!!$      integer(I8B), dimension(0:SMAXCHK-1) :: startpix
!!$      integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
!!$      real(DP),     dimension(0:SMAXCHK-1) :: sth
!!$
!!$      character(LEN=*), parameter :: code = 'almn2alnspring_pre'
!!$      integer(I4B) :: status
!!$      !=======================================================================
!!$
!!$      ! Healpix definitions
!!$      nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
!!$      npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
!!$      nphmx  = 4*nsmax           ! maximum number of pixels/ring
!!$      n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]
!!$      n_plm  = n_lm * nrings
!!$
!!$      !     --- allocates space for arrays ---
!!$      nchunks   = nrings/SMAXCHK + 1  ! number of chunks
!!$      chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK
!!$
!!$      allocate(b_north(0:nmmax, 0:chunksize-1),stat = status)
!!$      call assert_alloc(status,code,'b_north')
!!$
!!$      allocate(b_south(0:nmmax, 0:chunksize-1),stat = status)
!!$      call assert_alloc(status,code,'b_south')
!!$
!!$      if (.not. do_openmp()) then
!!$         allocate(dalm(0:1,0:nmmax), stat = status)
!!$         call assert_alloc(status,code,'dalm')
!!$         allocate(ring(0:nphmx-1),stat = status)
!!$         call assert_alloc(status,code,'ring')
!!$      endif
!!$      !     ------------ initiate variables and arrays ----------------
!!$
!!$      map = 0.0 ! set the whole map to zero
!!$
!!$      ! loop on chunks
!!$      do ichunk = 0, nchunks-1
!!$         lchk = ichunk * chunksize + 1
!!$         uchk = min(lchk+chunksize - 1, nrings)
!!$
!!$         do ith = lchk, uchk
!!$            ithl = ith - lchk !local index
!!$            ! get pixel location information
!!$            call get_pixel_layout(nsmax, ith, cth, sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
!!$         enddo
!!$
!!$         b_north(:,:) = 0_dpc ! pad with zeros
!!$         b_south(:,:) = 0_dpc
!!$
!!$         !$OMP parallel default(none) &
!!$         !$OMP shared(nsmax, nlmax, npix, nrings, nphmx, nmmax, sth, lchk, uchk, &
!!$         !$OMP plm, almn, n_lm, nph, startpix, kphi0, map) &
!!$         !$OMP private(dalm, b_north, b_south, b_ns, m, ith, ithl, nphl, l_min, l_start, &
!!$         !$OMP l, i_mm, istart_north, istart_south, ring, status)
!!$
!!$         if (do_openmp()) then
!!$            allocate(dalm(0:1,0:nmmax), stat = status)
!!$            call assert_alloc(status,code,'dalm')
!!$         endif
!!$
!!$         !$OMP do schedule(dynamic,1)
!!$
!!$         do l = 0, nlmax
!!$
!!$            do m = 0, l
!!$               dalm(0,m) =  real(almn(1,l,m),kind=dp)
!!$               dalm(1,m) = aimag(almn(1,l,m))
!!$            enddo ! loop on m 
!!$
!!$            do ithl = 0, uchk - lchk
!!$
!!$               do m = 0, l
!!$
!!$                  !l_min = l_min_ylm(m, sth(ithl))
!!$                  !if (nlmax >= l_min) then ! skip calculations when Ylm too small
!!$
!!$                  ith = ithl + lchk
!!$                  i_mm = n_lm * (ith-1) + ((2_I8B*nlmax + 3 - m)*m)/2 
!!$                  ! location of Ym,m for ring ith
!!$
!!$                  b_ns = 0.0_dpc
!!$
!!$                  ! odd values of (l+m)
!!$                  if (mod(m+l,2) == 0) then
!!$                     b_ns(-1) = cmplx(plm(i_mm+l-m) * dalm(0,m), plm(i_mm+l-m)*dalm(1,m), kind=DP)
!!$                  endif
!!$
!!$                  ! even values of (l+m)
!!$                  if (mod(m+l,2) == 1) then
!!$                     b_ns(1)  =  cmplx(plm(i_mm+l-m) * dalm(0,m), plm(i_mm+l-m)*dalm(1,m),kind=DP) 
!!$                  endif
!!$
!!$                  b_north(m,ithl) = b_ns(1) + b_ns(-1)
!!$                  b_south(m,ithl) = b_ns(1) - b_ns(-1)
!!$
!!$                  !endif ! test on nlmax
!!$
!!$               enddo ! loop on m
!!$
!!$               if (do_openmp()) then
!!$                  deallocate (dalm)
!!$               endif
!!$
!!$            enddo ! loop on ithl
!!$
!!$            if (do_openmp()) then
!!$               allocate(ring(0:nphmx-1),stat = status)
!!$               call assert_alloc(status,code,'ring')
!!$            endif
!!$
!!$            do ithl = 0, uchk - lchk
!!$
!!$               nphl = nph(ithl)
!!$               istart_north = startpix(ithl)
!!$               istart_south = npix-istart_north-nphl
!!$               ith  = ithl + lchk
!!$
!!$               call ring_synthesis(nsmax,nlmax,nmmax,b_north(0:nmmax,ithl),nphl,ring,kphi0(ithl))   ! north hemisph. + equator
!!$               map(l, istart_north:istart_north+nphl-1) = ring(0:nphl-1) 
!!$
!!$               if (ith < nrings) then
!!$                  call ring_synthesis(nsmax,nlmax,nmmax,b_south(0:nmmax,ithl),nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
!!$                  map(l, istart_south:istart_south+nphl-1) = ring(0:nphl-1) 
!!$               endif
!!$
!!$            enddo ! loop on rings (ithl)
!!$
!!$            if (do_openmp()) then
!!$               deallocate(ring)
!!$            endif
!!$
!!$         enddo ! loop on l
!!$
!!$         !$OMP end do
!!$
!!$         !$OMP end parallel
!!$
!!$      enddo    ! loop on chunks
!!$
!!$      !     --------------------
!!$      !     free memory and exit
!!$      !     --------------------
!!$      if (.not.do_openmp()) then
!!$         deallocate (ring,dalm)
!!$      endif
!!$      deallocate(b_north, b_south)
!!$      return
!!$
!!$    END SUBROUTINE almn2alnspring_pre
!!$
!!$
!!$    !===============================================================================
!!$    !===============================================================================
!!$
!!$    
!!$    SUBROUTINE alnspring2almn_iterative( nsmax, nlmax, nmmax, alnmap, &
!!$    	                & almn, zbounds, nb_iter )  
!!$
!!$    integer(I4B), intent(IN)   :: nsmax, nlmax, nmmax, nb_iter
!!$
!!$    real(DP), intent(IN), dimension(0:nlmax,0:(12_i8b*nsmax)*nsmax-1) :: alnmap
!!$
!!$    integer(I4B) :: n_plm, iter, l
!!$    real(DP), dimension(0:nlmax,0:(12_i8b*nsmax)*nsmax-1) :: alnmap_rec
!!$    
!!$    complex(DPC), dimension(1:1,0:nlmax,0:nmmax) :: almn
!!$
!!$    real(DP), intent(IN),  dimension(1:2), optional :: zbounds
!!$    real(DP), dimension(:,:), allocatable :: plm 
!!$
!!$    n_plm = nsmax*(nmmax+1)*(2*nlmax-nmmax+2) 
!!$    ALLOCATE(plm(0:n_plm-1,1:1))  
!!$    CALL plm_gen(nsmax, nlmax, nmmax, plm) 
!!$
!!$    alnmap_rec = alnmap
!!$    almn = 0.0
!!$    
!!$    DO iter=1,nb_iter
!!$       print*
!!$       print*
!!$       PRINT*, "Iteration : ", iter
!!$       print*
!!$        ! Improve almn : old_almn + d_almn ! almn = CALL alnspring2almn( alnmap_rec )
!!$        !CALL alnspring2almn_pre( nsmax, nlmax, nmmax, alnmap, almn, zbounds, plm )
!!$        CALL alnspring2almn( nsmax, nlmax, nmmax, alnmap_rec, almn, zbounds )
!!$        ! New reconstruction ! alnmap_rec = CALL almn2alnspring( almn )
!!$        CALL almn2alnspring_pre( nsmax, nlmax, nmmax, alnmap_rec, almn, plm(:,1) )
!!$        DO l=0,nlmax
!!$           print*, "Original : ", alnmap(l,0:6)
!!$           print*, "Reconstr : ", alnmap_rec(l,0:6)
!!$           print*
!!$        ENDDO
!!$        alnmap_rec = alnmap - alnmap_rec
!!$    ENDDO
!!$
!!$    END SUBROUTINE alnspring2almn_iterative
!!$
!!$
!!$    !===============================================================================
!!$    !===============================================================================


  SUBROUTINE alnspring2almn( nsmax, nlmax, nmmax, map, &
       & almn, zbounds )
    !=======================================================================
    !     computes the a(l,m) from a Temperature map for the HEALPIX pixelisation
    !        all from scratch
    !=======================================================================
    integer(I4B), intent(IN)                    :: nsmax, nlmax, nmmax
    real(DP),   intent(IN),  dimension(0:nlmax,0:(12_i8b*nsmax)*nsmax-1) :: map
    complex(DPC), dimension(1:1,0:nlmax,0:nmmax) :: almn
    real(DP),     intent(IN),  dimension(1:2),         optional :: zbounds

    real(DP), dimension(1:2)         :: zbounds_in
    real(DP), dimension(1:2*nsmax,1) :: w8ring_in
    integer(I4B) :: s, l, m, ith, scalem, scalel   ! alm related
    integer(I8B) :: istart_south, istart_north, npix  ! map related
    integer(I4B) :: nrings, nphmx
    real(DP)     :: omega_pix
    integer(i4b), parameter :: SMAXCHK = 30
    integer(kind=i4b), parameter :: RSMAX = 20, RSMIN = -20
    real(dp),          dimension(RSMIN:RSMAX) :: rescale_tab

    integer(I4B)                              :: par_lm
    real(DP)              :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)                      :: OVFLOW, UNFLOW
    real(DP),     dimension(-1:2)     :: phas_sd
    real(DP),     dimension(:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: mfac
    real(DP),     dimension(:,:), allocatable :: recfac

    integer(I4B)                              :: l_min
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:), allocatable :: phas_n, phas_s
    complex(DPC), dimension(:,:,:), allocatable :: phasl_n, phasl_s
    real(DP),     dimension(:),   allocatable :: ring
    integer(I4B)                   :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK-1) :: cth, sth
    logical(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it

    integer(I4B),     parameter :: LOG2LG   = 100
    real(KIND=DP),    parameter :: FL_LARGE = 2.0_dp **   LOG2LG
    real(KIND=DP),    parameter :: FL_SMALL = 2.0_dp ** (-LOG2LG)
    real(DP) :: logOVFLOW
    integer(i4b) :: smax

    character(LEN=*), PARAMETER :: code = 'MAP2ALMSPRING'
    integer(I4B) :: status
    !=======================================================================

    zbounds_in = (/-1.d0 , 1.d0/)
    if (present(zbounds)) zbounds_in = zbounds
    w8ring_in  = 1.d0

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    omega_pix = FOURPI / real(npix, kind=DP)  ! pixel area (identical for all pixels)

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    allocate(phasl_n(0:nlmax,0:nmmax,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phasl_n')

    allocate(phasl_s(0:nlmax,0:nmmax,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phasl_s')

    if (.not.do_openmp()) then
       allocate(ring(0:nphmx-1),stat = status)
       call assert_alloc(status,code,'ring')
       allocate(recfac(0:1,0:nlmax),stat = status)
       call assert_alloc(status,code,'recfac')
       allocate(dalm(1:2,0:nlmax),stat = status)
       call assert_alloc(status,code,'dalm')
       allocate(phas_n(0:nmmax,0:chunksize-1),stat = status)
       call assert_alloc(status,code,'phas_n')
       allocate(phas_s(0:nmmax,0:chunksize-1),stat = status)
       call assert_alloc(status,code,'phas_s')
    endif

    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax,mfac)

    logOVFLOW=log(FL_LARGE)
    smax = INT( log(MAX_DP) / logOVFLOW )

    if (smax > (RSMAX-1)) then
       print*,'Array rescale_tab too small in '//code
       print*,smax ,'>', RSMAX
       stop
    endif

    rescale_tab(RSMIN:RSMAX) = 0.0_dp
    do s = -smax, smax
       rescale_tab(s) = FL_LARGE ** s
    enddo
    rescale_tab(0) = 1.0_dp

    OVFLOW = rescale_tab(1)
    UNFLOW = rescale_tab(-1)
    almn(1:1,0:nlmax,0:nmmax) = 0.0 ! set the whole alm array to zero

    ! loop on chunks
    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          ! find out which rings are to be analysed
          call select_rings(cth(ithl), zbounds_in, keep_north(ithl), keep_south(ithl), keep_it(ithl))
       enddo

       phasl_n = 0_dpc
       phasl_s = 0_dpc

       !$OMP parallel default(none) &
       !$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
       !$OMP      lchk, uchk, nph, startpix, kphi0, w8ring_in, &
       !$OMP      phasl_n, phasl_s, keep_north, keep_south, chunksize,map) &
       !$OMP  private(ithl, nphl, istart_north, istart_south, l, ith, ring, status,  phas_n, phas_s)

       if (do_openmp()) then
          allocate(ring(0:nphmx-1),stat = status)
          call assert_alloc(status,code,'ring')
          allocate(phas_n(0:nmmax,0:chunksize-1),stat = status)
          call assert_alloc(status,code,'phas_n')
          allocate(phas_s(0:nmmax,0:chunksize-1),stat = status)
          call assert_alloc(status,code,'phas_s')
       endif

       phas_n = 0_dpc
       phas_s = 0_dpc

       !$OMP do schedule(dynamic,1)
       do l = 0, nlmax
          do ith = lchk, uchk
             ithl = ith - lchk !local index
             nphl = nph(ithl)
             istart_north = startpix(ithl)
             istart_south = npix-istart_north-nphl
             ! do Fourier Transform on rings
             if (keep_north(ithl)) then
                !
                ring(0:nphl-1) = map(l,istart_north:istart_north+nphl-1)! * w8ring_in(ith,1)
                call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_n(0,ithl), kphi0(ithl))
                phasl_n(l,0:nmmax,ithl) = phas_n(0:nmmax,ithl) + phasl_n(l,0:nmmax,ithl)
                !
             endif

             if (ith < nrings .and. keep_south(ithl)) then
                !
                ring(0:nphl-1) = map(l,istart_south:istart_south+nphl-1)! * w8ring_in(ith,1)
                call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_s(0,ithl), kphi0(ithl))
                phasl_s(l,0:nmmax,ithl) = phas_s(0:nmmax,ithl) + phasl_s(l,0:nmmax,ithl)
                !
             endif
          enddo ! loop on ring
       enddo
       !$OMP end do


       if (do_openmp()) then
          deallocate(ring)
          deallocate(phas_s,phas_n)
       endif
       !$OMP end parallel

       call init_rescale()
       OVFLOW = rescale_tab(1)
       UNFLOW = rescale_tab(-1)

       !$OMP parallel default(none) NUM_THREADS(8) &
       !$OMP shared(nlmax, nmmax, lchk, uchk, rescale_tab, ovflow, unflow, &
       !$OMP    cth, sth, mfac, almn, phasl_n, phasl_s, keep_it, omega_pix) &
       !$OMP private(recfac, dalm, phas_sd, status, m, ithl, l_min, &
       !$OMP   scalem, scalel, corfac, par_lm, lam_mm, lam_lm, lam_0, lam_1, lam_2, &
       !$OMP   cth_ring, l, php, phm)

       if (do_openmp()) then
          allocate(recfac(0:1,0:nlmax),stat = status)
          call assert_alloc(status,code,'recfac')
          allocate(dalm(1:2,0:nlmax),stat = status)
          call assert_alloc(status,code,'dalm')
       endif

       !$OMP do schedule(dynamic,1)
       do m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          call gen_recfac(nlmax, m, recfac)

          ! introduce double precision vector to perform summation over ith for each l
          dalm(1:2, m:nlmax ) = 0.0_dp

          do ithl = 0, uchk - lchk
             l_min = l_min_ylm(m, sth(ithl))
             if (keep_it(ithl) .and. nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001)
                ! determine lam_mm
                call compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)

                !           ---------- l = m ----------
                par_lm = 1

               	php = phasl_n(m,m,ithl) + phasl_s(m,m,ithl) ! sum  (if (l+m) even)
               	phm = phasl_n(m,m,ithl) - phasl_s(m,m,ithl) ! diff (if (l+m) odd)
               	phas_sd(-1:0) =  (/ real(phm, kind=dp), aimag(phm) /)
               	phas_sd(1:2) =  (/ real(php, kind=dp), aimag(php) /)
                if (m >= l_min) then
                   lam_lm = lam_mm * corfac !Actual lam_mm 
                   dalm(1:2, m) = dalm(1:2, m) + phas_sd(par_lm:par_lm+1) *lam_lm
                endif


                !           ---------- l > m ----------
                lam_0 = 0.0_dp
                lam_1 = 1.0_dp
                scalel=0
                cth_ring = cth(ithl)
                lam_2 = cth_ring * lam_1 * recfac(0,m)

                do l = m+1, nlmax

                   php = phasl_n(l,m,ithl) + phasl_s(l,m,ithl) ! sum  (if (l+m) even)
                   phm = phasl_n(l,m,ithl) - phasl_s(l,m,ithl) ! diff (if (l+m) odd)
                   phas_sd(-1:0) =  (/ real(phm, kind=dp), aimag(phm) /)
                   phas_sd(1:2) =  (/ real(php, kind=dp), aimag(php) /)

                   par_lm = - par_lm  ! = (-1)^(l+m)

                   if (l >= l_min) then
                      lam_lm = lam_2 * corfac * lam_mm
                      dalm(1:2, l) = dalm(1:2, l) &
                           &       +  phas_sd(par_lm:par_lm+1) *lam_lm
                   endif

                   lam_0 = lam_1 * recfac(1,l-1)
                   lam_1 = lam_2
                   lam_2 = (cth_ring * lam_1 - lam_0) * recfac(0,l)
                   if (abs(lam_2) > OVFLOW) then
                      lam_1 = lam_1*UNFLOW
                      lam_2 = lam_2*UNFLOW
                      scalel= scalel + 1
                      corfac= rescale_tab(max(scalem+scalel,RSMIN))
                   elseif (abs(lam_2) < UNFLOW) then
                      lam_1 = lam_1*OVFLOW
                      lam_2 = lam_2*OVFLOW
                      scalel= scalel - 1
                      corfac= rescale_tab(max(scalem+scalel,RSMIN))
                   endif

                enddo ! loop on l
             endif ! test on cut sky and nlmax
          enddo ! loop on ithl
          do l = m, nlmax
             almn(1, l, m) = almn(1, l, m) + cmplx(dalm(1, l), dalm(2, l), kind=DP)
          enddo
       enddo ! loop on m
       !$OMP end do
       if (do_openmp()) then
          deallocate (recfac,dalm)
       endif
       !$OMP end parallel
    enddo ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------

    deallocate(mfac)
    if (.not.do_openmp()) then
       deallocate (ring,recfac,dalm)
       deallocate(phas_n,phas_s)
    endif
    RETURN
    !
  END SUBROUTINE alnspring2almn


  !===============================================================================
  !===============================================================================

  SUBROUTINE survey2almn_opt( nsmax, nnmax, nlmax, nmmax, survey, nbpts, almn, kln, zbounds )
    !=======================================================================
    !     computes the a(n,l,m) decomposition from scratch
    !=======================================================================
    integer(I4B), intent(IN)                        :: nsmax, nnmax, nlmax, nmmax, nbpts
    real(DP),   intent(IN), dimension(1:nbpts,1:3)  :: survey
    complex(DPC), dimension(1:nlmax,0:nlmax,0:nmmax)    :: almn
    real(DP),   intent(IN), dimension(0:nlmax,1:nnmax) :: kln
    real(DP), intent(IN), dimension(1:2),  optional :: zbounds
    real(DP),     dimension(:,:),   allocatable     :: map_north, map_south

    real(DP), dimension(1:2)         :: zbounds_in
    real(DP), dimension(1:3)         :: surveytemp 
    real(DP), dimension(1:2*nsmax,1) :: w8ring_in
    integer(I4B) :: s, l, m, n, p, ipring, ith, scalem, scalel  
    integer(I8B) :: istart_south, istart_north, npix  
    integer(I8B) :: itotstart_north, itotstart_south, itotend_north, itotend_south
    integer(I4B) :: nrings, nphmx, sizenorth, sizesouth
    real(DP)     :: omega_pix, jln
    integer(i4b), parameter :: SMAXCHK = 30
    integer(kind=i4b), parameter :: RSMAX = 20, RSMIN = -20
    real(dp),          dimension(RSMIN:RSMAX) :: rescale_tab

    integer(I4B)                              :: par_lm
    real(DP)              :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    real(DP)                      :: OVFLOW, UNFLOW
    real(DP),     dimension(-1:2)     :: phas_sd
    real(DP),     dimension(:,:), allocatable :: dalm
    real(DP),     dimension(:),   allocatable :: mfac
    real(DP),     dimension(:,:), allocatable :: recfac
    real(DP), dimension(:), allocatable :: tempaccjln

    integer(I4B)                              :: l_min
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:), allocatable :: phas_n, phas_s
    complex(DPC), dimension(:,:,:), allocatable :: phasl_n, phasl_s
    real(DP),     dimension(:),   allocatable :: ring
    integer(I4B)                   :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK-1) :: cth, sth
    logical(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it

    integer(I4B),     parameter :: LOG2LG   = 100
    real(KIND=DP),    parameter :: FL_LARGE = 2.0_dp **   LOG2LG
    real(KIND=DP),    parameter :: FL_SMALL = 2.0_dp ** (-LOG2LG)
    real(DP) :: logOVFLOW
    integer(i4b) :: smax

    character(LEN=*), PARAMETER :: code = 'survey2almn_opt'
    integer(I4B) :: status

    !=======================================================================

    zbounds_in = (/-1.d0 , 1.d0/)
    if (present(zbounds)) zbounds_in = zbounds
    w8ring_in  = 1.d0

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    omega_pix = FOURPI / real(npix, kind=DP)  ! pixel area (identical for all pixels)

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(mfac(0:nmmax),stat = status)
    call assert_alloc(status,code,'mfac')

    allocate(phasl_n(0:nlmax,0:nmmax,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phasl_n')

    allocate(phasl_s(0:nlmax,0:nmmax,0:chunksize-1),stat = status)
    call assert_alloc(status,code,'phasl_s')

    if (.not.do_openmp()) then
       allocate(ring(0:nphmx-1),stat = status)
       call assert_alloc(status,code,'ring')
       allocate(recfac(0:1,0:nlmax),stat = status)
       call assert_alloc(status,code,'recfac')
       allocate(dalm(1:2,0:nlmax),stat = status)
       call assert_alloc(status,code,'dalm')
       allocate(phas_n(0:nmmax,0:chunksize-1),stat = status)
       call assert_alloc(status,code,'phas_n')
       allocate(phas_s(0:nmmax,0:chunksize-1),stat = status)
       call assert_alloc(status,code,'phas_s')
    endif

    !     ------------ initiate variables and arrays ----------------

    call gen_mfac(nmmax,mfac)

    logOVFLOW=log(FL_LARGE)
    smax = INT( log(MAX_DP) / logOVFLOW )

    if (smax > (RSMAX-1)) then
       print*,'Array rescale_tab too small in '//code
       print*,smax ,'>', RSMAX
       stop
    endif

    rescale_tab(RSMIN:RSMAX) = 0.0_dp
    do s = -smax, smax
       rescale_tab(s) = FL_LARGE ** s
    enddo
    rescale_tab(0) = 1.0_dp

    OVFLOW = rescale_tab(1)
    UNFLOW = rescale_tab(-1)
    almn(1:nnmax,0:nlmax,0:nmmax) = 0.0 ! set the whole alm array to zero

    print*,"================================================"
    print*,"   Computing 3D cones and partial maps"
    print*,"   > Modified FFT 3D algorithm"
    print*,"------------------------------------------------"

    ! loop on chunks
    do ichunk = 0, nchunks-1 ! for each chunk

       print*,"> Block : ",ichunk+1," on ", nchunks

       lchk = ichunk * chunksize + 1  ! first ring of the chunk
       uchk = min(lchk+chunksize - 1, nrings)  ! last ring of the chunk

       ! locate pixels
       do ith = lchk, uchk ! for each ring
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          ! find out which rings are to be analysed
          call select_rings(cth(ithl), zbounds_in, keep_north(ithl), keep_south(ithl), keep_it(ithl))
       enddo

       itotstart_north = startpix(0)
       itotend_north = startpix(uchk-lchk)+nph(uchk-lchk)-1
       itotstart_south = npix-startpix(uchk-lchk)-nph(uchk-lchk)
       itotend_south = npix-startpix(0)-1
       sizenorth = (itotend_north-itotstart_north)
       sizesouth = (itotend_south-itotstart_south)

       !print*,"North hemisphere:",itotstart_north,itotend_north
       !print*,"South hemisphere:",itotstart_south,itotend_south

       allocate(map_north(0:nlmax,0:sizenorth),stat = status)
       call assert_alloc(status,code,'map_north')
       allocate(map_south(0:nlmax,0:sizesouth),stat = status)
       call assert_alloc(status,code,'map_south')

       ! loop on n parameter 
       do n=1, nnmax

          !print*,"> N = ",n 
          map_north=0.0_dp
          map_south=0.0_dp

          !$OMP parallel default(none) &
          !$OMP private(p, l, ipring, jln, tempaccjln, status, surveytemp ) &
          !$OMP shared(n, survey, nsmax, nbpts, nlmax, kln, map_north, map_south, &
          !$OMP        itotend_north, itotstart_south, itotend_south, itotstart_north )

          allocate(tempaccjln(0:nlmax),stat = status)
          call assert_alloc(status,code,'tempaccjln')

          !$OMP do reduction(+:map_north, map_south) schedule(dynamic,1024) 
          do p=1, nbpts
             surveytemp(1:3) = survey(p,1:3)
             CALL ang2pix_ring(nsmax, surveytemp(2), surveytemp(1), ipring) 
             do l=0, nlmax
                CALL BJL( l , kln(l,n)*surveytemp(3), jln ) 
                tempaccjln(l) =  kln(l,n)*jln
             enddo
             if ( ipring .GE. itotstart_north .AND. ipring .LE. itotend_north ) then ! if north
                map_north(:,ipring-itotstart_north) =  map_north(:,ipring-itotstart_north) + tempaccjln(:)
             endif
             if ( ipring .GE. itotstart_south .AND. ipring .LE. itotend_south ) then ! if south
                map_south(:,ipring-itotstart_south) =  map_south(:,ipring-itotstart_south) + tempaccjln(:)
             endif
          enddo
          !$OMP end do

          deallocate(tempaccjln)

          !$OMP end parallel

          phasl_n = 0_dpc
          phasl_s = 0_dpc

          !$OMP parallel default(none) &
          !$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, n, &
          !$OMP      itotstart_north, itotend_north, itotstart_south, itotend_south, &
          !$OMP      lchk, uchk, nph, startpix, kphi0, w8ring_in, &
          !$OMP      phasl_n, phasl_s, keep_north, keep_south, chunksize, map_north, map_south ) &
          !$OMP  private(ithl, nphl, istart_north, istart_south, l, ith, ring, status,  phas_n, phas_s)

          if (do_openmp()) then
             allocate(ring(0:nphmx-1),stat = status)
             call assert_alloc(status,code,'ring')
             allocate(phas_n(0:nmmax,0:chunksize-1),stat = status)
             call assert_alloc(status,code,'phas_n')
             allocate(phas_s(0:nmmax,0:chunksize-1),stat = status)
             call assert_alloc(status,code,'phas_s')
          endif

          phas_n = 0_dpc
          phas_s = 0_dpc

          !$OMP do schedule(dynamic,1)
          do l = 0, nlmax

             do ith = lchk, uchk
                ithl = ith - lchk !local index
                nphl = nph(ithl)
                istart_north = startpix(ithl)
                istart_south = npix-istart_north-nphl
                ! do Fourier Transform on rings
                if (keep_north(ithl)) then
                   !
                   ring(0:nphl-1) = map_north(l,(istart_north-itotstart_north):(istart_north+nphl-1-itotstart_north))! * w8ring_in(ith,1)
                   call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_n(0,ithl), kphi0(ithl))
                   phasl_n(l,0:nmmax,ithl) = phas_n(0:nmmax,ithl) + phasl_n(l,0:nmmax,ithl)
                   !
                endif
                if (ith < nrings .and. keep_south(ithl) ) then
                   !
                   ring(0:nphl-1) = map_south(l,(istart_south-itotstart_south):&
                        &(istart_south+nphl-1-itotstart_south))! * w8ring_in(ith,1)
                   call ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_s(0,ithl), kphi0(ithl))
                   phasl_s(l,0:nmmax,ithl) = phas_s(0:nmmax,ithl) + phasl_s(l,0:nmmax,ithl)
                   !
                endif
             enddo ! loop on ring
          enddo ! loop on l
          !$OMP end do

          if (do_openmp()) then
             deallocate(ring)
             deallocate(phas_s,phas_n)
          endif
          !$OMP end parallel

          call init_rescale()
          OVFLOW = rescale_tab(1)
          UNFLOW = rescale_tab(-1)

          !$OMP parallel default(none) &
          !$OMP shared(nlmax, nmmax, lchk, uchk, rescale_tab, n, ovflow, unflow, &
          !$OMP    cth, sth, mfac, almn, phasl_n, phasl_s, keep_it, omega_pix) &
          !$OMP private(recfac, dalm, phas_sd, status, m, ithl, l_min, &
          !$OMP   scalem, scalel, corfac, par_lm, lam_mm, lam_lm, lam_0, lam_1, lam_2, &
          !$OMP   cth_ring, l, php, phm)

          if (do_openmp()) then
             allocate(recfac(0:1,0:nlmax),stat = status)
             call assert_alloc(status,code,'recfac')
             allocate(dalm(1:2,0:nlmax),stat = status)
             call assert_alloc(status,code,'dalm')
          endif

          !$OMP do reduction(+:almn) schedule(dynamic,1)
          do m = 0, nmmax
             ! generate recursion factors (recfac) for Ylm of degree m
             call gen_recfac(nlmax, m, recfac)

             ! introduce double precision vector to perform summation over ith for each l
             dalm(1:2, m:nlmax ) = 0.0_dp

             do ithl = 0, uchk - lchk
                l_min = 0 !l_min_ylm(m, sth(ithl))
                if (keep_it(ithl) .and. nlmax >= l_min) then ! avoid un-necessary calculations (EH, 09-2001)
                   ! determine lam_mm
                   call compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)

                   !           ---------- l = m ----------
                   par_lm = 1

                   php = phasl_n(m,m,ithl) + phasl_s(m,m,ithl) ! sum  (if (l+m) even)
                   phm = phasl_n(m,m,ithl) - phasl_s(m,m,ithl) ! diff (if (l+m) odd)
                   phas_sd(-1:0) =  (/ real(phm, kind=dp), aimag(phm) /)
                   phas_sd(1:2) =  (/ real(php, kind=dp), aimag(php) /)
                   if (m >= l_min) then
                      lam_lm = lam_mm * corfac !Actual lam_mm 
                      dalm(1:2, m) = dalm(1:2, m) + phas_sd(par_lm:par_lm+1) *lam_lm
                   endif


                   !           ---------- l > m ----------
                   lam_0 = 0.0_dp
                   lam_1 = 1.0_dp
                   scalel=0
                   cth_ring = cth(ithl)
                   lam_2 = cth_ring * lam_1 * recfac(0,m)

                   do l = m+1, nlmax

                      php = phasl_n(l,m,ithl) + phasl_s(l,m,ithl) ! sum  (if (l+m) even)
                      phm = phasl_n(l,m,ithl) - phasl_s(l,m,ithl) ! diff (if (l+m) odd)
                      phas_sd(-1:0) =  (/ real(phm, kind=dp), aimag(phm) /)
                      phas_sd(1:2) =  (/ real(php, kind=dp), aimag(php) /)

                      par_lm = - par_lm  ! = (-1)^(l+m)

                      if (l >= l_min) then
                         lam_lm = lam_2 * corfac * lam_mm
                         dalm(1:2, l) = dalm(1:2, l) &
                              &       +  phas_sd(par_lm:par_lm+1) *lam_lm
                      endif

                      lam_0 = lam_1 * recfac(1,l-1)
                      lam_1 = lam_2
                      lam_2 = (cth_ring * lam_1 - lam_0) * recfac(0,l)
                      if (abs(lam_2) > OVFLOW) then
                         lam_1 = lam_1*UNFLOW
                         lam_2 = lam_2*UNFLOW
                         scalel= scalel + 1
                         corfac= rescale_tab(max(scalem+scalel,RSMIN))
                      elseif (abs(lam_2) < UNFLOW) then
                         lam_1 = lam_1*OVFLOW
                         lam_2 = lam_2*OVFLOW
                         scalel= scalel - 1
                         corfac= rescale_tab(max(scalem+scalel,RSMIN))
                      endif

                   enddo ! loop on l
                endif ! test on cut sky and nlmax
             enddo ! loop on ithl
             do l = m, nlmax
                almn(n, l, m) = almn(n, l, m) + cmplx(dalm(1, l), dalm(2, l), kind=DP)
             enddo
          enddo ! loop on m
          !$OMP end do
          if (do_openmp()) then
             deallocate (recfac,dalm)
          endif
          !$OMP end parallel

       enddo ! loop on n

       deallocate(map_north, map_south)

    enddo ! loop on chunks

    print*,"------------------------------------------------"
    print*,"   Finished"
    print*,"================================================"

    !     --------------------
    !     free memory and exit
    !     --------------------

    deallocate(mfac)
    if (.not.do_openmp()) then
       deallocate (ring,recfac,dalm)
       deallocate(phas_n,phas_s)
    endif
    RETURN
    !
  END SUBROUTINE survey2almn_opt


  !===============================================================================
  !===============================================================================


  !> Performs the Fourier-Bessel decomposition (backward algorithm) of a discrete survey.
  !!@param[in] (nnmax, nlmax, nmmax) : decomposition triplet
  !!@param[in] nsmax : Healpix nside parameter 
  !!@param[in] (rmax, zbounds) : radial and angular cut 
  !!@param[in] kln : Hankel (scaling) coefficients
  !!@param[in] (survey, nbpts) : survey array and its length
  !!@param[out] almn : almn decomposition coefficients
  SUBROUTINE survey2almn_srs( nsmax, nnmax, nlmax, nmmax, rmax, nbpts, &
       & zbounds, survey, kln, almn)
    !
    integer(I4B)  :: status, nbpts, ipring, nrr, nb_iter
    integer(I4B)  :: nsmax, nnmax, nlmax, nmmax, n, p, k, m, l
    real(DP)	  :: rmax, tempval, order
    !	
    real(DP),    dimension(:,:), allocatable	     :: map
    real(DP), dimension(1:nbpts,1:3)         	     :: survey
    real(DP), dimension(1:3)                         :: surveytemp
    real(DP),    dimension(1:2)       	             :: zbounds
    real(DP),    dimension(0:nlmax,1:nnmax)         :: kln
    !real(DP),    dimension(:,:,:) , allocatable     :: jln
    real(DP)                                        :: jln
    real(DP),    dimension(:,:)   , allocatable     ::  plm
    complex(DP), dimension(1:nnmax,0:nlmax,0:nmmax) :: almn
    character(len=*), PARAMETER                     :: code = "survey2almn_srs"
    !

    almn = cmplx( 0.0, 0.0, kind=DP )

120 FORMAT(A,I3)

    ALLOCATE(map(0:nlmax,0:(12*nsmax**2-1)),stat = status)
    CALL assert_alloc(status,code,"map")

    print*,"================================================"
    print*,"   Computing multidimensional maps" 
    print*,"   > Modified FFT 3D algorithm"
    print*,"------------------------------------------------"

    DO n = 1, nnmax

       print*,"> Order n = ",n!,OMP_GET_THREAD_NUM()
       map = 0.0_dp

       !$OMP PARALLEL DO REDUCTION(+:map)  SCHEDULE(DYNAMIC,512) &
       !$OMP SHARED(nnmax,nsmax,nlmax,nmmax,almn,n,survey,kln,zbounds) &
       !$OMP PRIVATE(status,p,ipring,l,jln)

       DO p = 1, nbpts
          !
          surveytemp(1:3) = survey(p,1:3)
          CALL ang2pix_ring(nsmax, surveytemp(2), surveytemp(1), ipring)
          ! 	
          DO l = 0, nlmax
             !
             jln = 0.0_dp
             CALL BJL( l , kln(l,n)*surveytemp(3) , jln )
             !$OMP ATOMIC
             map(l,ipring) = map(l,ipring) + kln(l,n) * jln
             !
          ENDDO
       ENDDO

       !$OMP END PARALLEL DO

       CALL alnspring2almn( nsmax, nlmax, nmmax, map, &
            & almn(n:n, 0:nlmax,0:nmmax), zbounds )

    ENDDO

    print*,"------------------------------------------------"
    print*,"   Finished"
    print*,"================================================"


    DEALLOCATE( map )    

    RETURN
    !
  END subroutine survey2almn_srs


  !===============================================================================
  !===============================================================================

  SUBROUTINE fieldalmn2overalmn( almn, kln, rhomean, rmax, nnmax, nlmax, nmmax )

    IMPLICIT NONE

    INTEGER(I4B) :: nnmax, nlmax, nmmax, n
    REAL(DP) :: rmax, rhomean
    COMPLEX(DP),DIMENSION(1:nnmax,0:nlmax,0:nmmax) :: almn
    REAL(DP),DIMENSION(0:nlmax,1:nnmax) :: kln
    REAL(DP),DIMENSION(1:nnmax,1:1,1:1) :: jln

    print*, "rhomean", rhomean

    almn(:,:,:) = almn(:,:,:) / rhomean

    DO n = 1, nnmax
       CALL BJL( 1 , kln(1,n)*rmax , jln(n,1,1) )
       print*,sqrt(FOURPI)*(rmax**2.0)*jln(n,1,1)
       almn(n,0,0) = almn(n,0,0) - sqrt(FOURPI)*(rmax**2.0)*jln(n,1,1)
    ENDDO

    RETURN

  END SUBROUTINE fieldalmn2overalmn


  !===============================================================================
  !===============================================================================


  function wtime( )

    implicit none

    integer ( kind = 4 ) clock_max
    integer ( kind = 4 ) clock_rate
    integer ( kind = 4 ) clock_reading
    real    ( kind = 8 ) wtime

    call system_clock ( clock_reading, clock_rate, clock_max )

    wtime = real ( clock_reading, kind = 8 ) &
         / real ( clock_rate, kind = 8 )

    return
  end function wtime


  !===============================================================================
  !===============================================================================


  SUBROUTINE almn2rmap(map, almn, rho, nsmax, nnmax, nlmax, nmmax, kln, cln, plm)

    IMPLICIT NONE

    real(DP)             :: rho
    integer(I4B), intent(IN)                     :: nsmax, nnmax, nlmax, nmmax
    complex(DPC), intent(IN),  dimension(1:nnmax,0:nlmax,0:nmmax)      :: almn
    real(DP),     intent(OUT), dimension(0:(12_i8b*nsmax)*nsmax-1,1:1) :: map
    real(DP),     intent(IN),  dimension(0:( nsmax*(nmmax+1)*(2*nlmax-nmmax+2) ) )   :: plm

    integer(I4B) :: l, m, ith, n
    integer(I8B) :: istart_south, istart_north, npix
    integer(I4B) :: nrings, nphmx

    integer(I8B) :: n_lm, n_plm, i_mm
    real(DP),     dimension(0:nlmax,1:nnmax)  :: cln, kln,jlns
    complex(DPC), dimension(-1:1)             :: b_ns
    real(DP),     dimension(:,:), allocatable :: dalm
    integer(i4b)                              :: ll, l_min, l_start
    real(DP)                                  :: cth
    complex(DPC), dimension(:,:), allocatable :: b_north, b_south
    real(DP),     dimension(:),   allocatable :: ring
    integer(i4b)             :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    integer(i4b), parameter :: SMAXCHK = 50
    integer(I8B), dimension(0:SMAXCHK-1) :: startpix
    integer(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    real(DP),     dimension(0:SMAXCHK-1) :: sth

    character(LEN=*), parameter :: code = 'ALM2MAP'
    integer(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]
    n_plm  = n_lm * nrings

    !     --- allocates space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    allocate(b_north(0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_north')

    allocate(b_south(0:nmmax, 0:chunksize-1),stat = status)
    call assert_alloc(status,code,'b_south')

    if (.not. do_openmp()) then
       allocate(dalm(0:1,0:nlmax), stat = status)
       call assert_alloc(status,code,'dalm')
       allocate(ring(0:nphmx-1),stat = status)
       call assert_alloc(status,code,'ring')
    endif


    !     ------------ Compute required bessel functions ----------------

    CALL gen_jln(jlns, kln, rho, nnmax, nlmax)


    !     ------------ initiate variables and arrays ----------------

    map = 0.0 ! set the whole map to zero

    ! loop on chunks
    do ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       do ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          call get_pixel_layout(nsmax, ith, cth, sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
       enddo

       b_north(:,:) = 0_dpc ! pad with zeros
       b_south(:,:) = 0_dpc

       !$OMP parallel default(none) NUM_THREADS(8) &
       !$OMP shared(nlmax, nmmax, sth, lchk, uchk, plm, b_north, b_south, n_lm, jlns, cln, almn, nnmax) &
       !$OMP private(dalm, b_ns, status, m, ith, ithl, l_min, l_start, l, ll, i_mm)

       if (do_openmp()) then
          allocate(dalm(0:1,0:nlmax), stat = status)
          call assert_alloc(status,code,'dalm')
       endif


       !$OMP do schedule(dynamic,1)
       do m = 0, nmmax
    	  dalm = 0.0_dp
       ! extract needed alm under memory and CPU efficient form
          do ll = m, nlmax
             do n = 0, nnmax
                dalm(0,ll) = dalm(0,ll) + real(almn(n,ll,m),kind=dp) * jlns(ll,n) * cln(ll,n)
                dalm(1,ll) = dalm(1,ll) + aimag(almn(n,ll,m)) * jlns(ll,n) * cln(ll,n)
             enddo
          enddo
          do ithl = 0, uchk - lchk
             l_min = l_min_ylm(m, sth(ithl))
             if (nlmax >= l_min) then ! skip calculations when Ylm too small
                ith = ithl + lchk
                i_mm = n_lm * (ith-1) + ((2_I8B*nlmax + 3 - m)*m)/2 ! location of Ym,m for ring ith

                !           ---------- l = m ----------

                if (m >= l_min) then
                   b_ns( 1) = cmplx(plm(i_mm) * dalm(0,m), plm(i_mm) * dalm(1,m), kind=DP)
                   b_ns(-1) = 0.0_dpc
                else
                   b_ns = 0.0_dpc
                endif

                !           ---------- l > m ----------
                l_start = max(m+1, l_min) ! odd values of (l+m)
                if (mod(m+l_start,2) == 0) l_start = l_start+1
                do l = l_start, nlmax, 2
                   b_ns(-1) = b_ns(-1) + &
                        & cmplx(plm(i_mm+l-m) * dalm(0,l), plm(i_mm+l-m)*dalm(1,l), kind=DP)
                enddo ! loop on l

                l_start = max(m+2, l_min) ! even values of (l+m)
                if (mod(m+l_start,2) == 1) l_start = l_start+1
                do l = l_start, nlmax, 2
                   b_ns(1)  = b_ns(1) + &
                        & cmplx(plm(i_mm+l-m) * dalm(0,l), plm(i_mm+l-m)*dalm(1,l),kind=DP) 
                enddo ! loop on l

                b_north(m,ithl) = b_ns(1) + b_ns(-1)
                b_south(m,ithl) = b_ns(1) - b_ns(-1)
             endif ! test on nlmax
          enddo ! loop on rings (ithl)
       enddo ! loop on m
       !$OMP end do

       if (do_openmp()) then
          deallocate (dalm)
       endif
       !$OMP end parallel

       !$OMP parallel default(none) NUM_THREADS(8) &
       !$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
       !$OMP      lchk, uchk, b_north, b_south, nph, startpix, kphi0, map) &
       !$OMP  private(ithl, nphl, istart_north, istart_south, &
       !$OMP      ith, ring, status)
       if (do_openmp()) then
          allocate(ring(0:nphmx-1),stat = status)
          call assert_alloc(status,code,'ring')
       endif
       !$OMP do schedule(dynamic,1)
       do ithl = 0, uchk - lchk
          nphl = nph(ithl)
          istart_north = startpix(ithl)
          istart_south = npix-istart_north-nphl
          ith  = ithl + lchk

          call ring_synthesis(nsmax,nlmax,nmmax,b_north(0,ithl),nphl,ring,kphi0(ithl))   ! north hemisph. + equator
          map(istart_north:istart_north+nphl-1,1) = map(istart_north:istart_north+nphl-1,1) + ring(0:nphl-1)

          if (ith < nrings) then
             call ring_synthesis(nsmax,nlmax,nmmax,b_south(0,ithl),nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
             map(istart_south:istart_south+nphl-1,1) = map(istart_south:istart_south+nphl-1,1) + ring(0:nphl-1)
          endif
       enddo ! loop on ithl
       !$OMP end do
       if (do_openmp()) then
          deallocate(ring)
       endif
       !$OMP end parallel
    enddo    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    if (.not.do_openmp()) then
       deallocate (ring,dalm)
    endif
    deallocate(b_north, b_south)
    return

  END SUBROUTINE almn2rmap


  !===============================================================================
  !===============================================================================


  SUBROUTINE BSORT (X, IY, N)

    IMPLICIT NONE    

    integer(kind=I4B) 		       :: N, I, J, JMAX
    real(kind=DP), DIMENSION(1:N)     :: X
    real(kind=DP), DIMENSION(1:N,1:3) :: IY
    real(kind=DP)		       :: TEMP
    real(kind=DP), dimension(1:3)     ::ITEMP

    JMAX=N-1
    DO 200 I=1,N-1
       !print*,"     SORT : ",i,"/",N
       TEMP=1.E38
       DO 100 J=1,JMAX
          IF(X(J).GT.X(J+1)) GO TO 100
          TEMP=X(J)
          X(J)=X(J+1)
          X(J+1)=TEMP
          ITEMP=IY(J,:)
          IY(J,:)=IY(J+1,:)
          IY(J+1,:)=ITEMP
100       CONTINUE
          IF(TEMP.EQ.1.E38) GO TO 300	       
          JMAX=JMAX-1
200       CONTINUE	

300       RETURN

        END SUBROUTINE BSORT



        !===============================================================================
        !===============================================================================

        SUBROUTINE gen_qln(qln, nnmax, nlmax)

          IMPLICIT NONE

          INTEGER(I4B) :: nnmax, nlmax, l, nrr, maxrt, lastRoot, err
          PARAMETER(maxrt = 1000)
          REAL(DP) :: rmax, order,A,B, A_before
          REAl(DP), DIMENSION(0:nlmax,1:nnmax) :: qln
          !	REAL*8, DIMENSION(1:maxrt) :: roots
          REAL*8, DIMENSION(:), allocatable :: roots
          character(LEN=*), parameter :: code = 'BESSEL'
          integer(I4B) :: status

          print*,"------------------------------------------------"
          print*,"   Computing scaling coefficients..."
          print*,"-----------------------"


          !$OMP parallel &
          !$OMP shared(qln,rmax,nnmax,nlmax,status) &
          !$OMP private(roots,order,A,B,A_before,l,nrr,lastRoot,err)

          allocate(roots(1:maxrt),stat = status)
          call assert_alloc(status,code,'roots')

          !$OMP do schedule(dynamic,1)
          DO l = 0, nlmax

             order = real(l) + 0.5_dp
             lastRoot = 1

             PRINT*, "Order l =",l

             A=0.0
             B=2000.0

             DO WHILE ( lastRoot .LT. nnmax)

                roots = 0.0_sp
                A_before = A
                CALL ROOTBESSJ( order, A, B, maxrt, nrr, err, roots )
                !		print*,"number of roots found",nrr,err
                if(nrr .EQ. 0) then
                   B = B + 1000.0
                else if(err .NE. 1) then
                   print*,"err code :",err
                   A = A_before
                   B = (B - A)/2.0 + A
                else

                   if(lastRoot + nrr .LE. nnmax) then
                      qln(l,lastRoot:lastRoot + nrr - 1) = roots(1:nrr)
                      lastRoot = lastRoot + nrr - 1
                   else
                      qln(l,lastRoot:nnmax) = roots(1: nnmax - lastRoot + 1)
                      lastRoot = nnmax
                   endif
                   B = A + 2048.0
                endif

             ENDDO

          ENDDO
          !$OMP end do
          deallocate(roots)
          !$OMP end parallel

          print*,"-----------------------"
          print*,"   Done."
          print*,"------------------------------------------------"


        END SUBROUTINE gen_qln

        !===============================================================================
        !===============================================================================	


        SUBROUTINE gen_kln(kln, nnmax, nlmax, rmax)

          IMPLICIT NONE

          INTEGER(I4B) :: nnmax, nlmax, l, nrr
          REAL(DP) :: rmax, order
          REAl(DP), DIMENSION(0:nlmax,1:nnmax) :: kln
          !	REAL*8, DIMENSION(1:1000) :: roots

          !	DO l = 0, nlmax
          !
          !	    order = real(l) + 0.5_dp
          !
          !	    roots = 0.0_sp
          !
          !	    CALL ROOTBESSJ( order, nrr, roots )
          !
          !	    !print*,order,nrr,roots
          !	    !if (nrr <= nnmax) then
          !	    !   print*, "ERROR during roots computation"
          !	    !   STOP
          !	    !else
          !
          !	    kln(l,1:nnmax) = roots(1:(nnmax)) / rmax
          !
          !	    !endif
          !
          !	ENDDO

          call gen_qln(kln,nnmax,nlmax)

          kln = kln / rmax


        END SUBROUTINE gen_kln


        !===============================================================================
        !===============================================================================


        SUBROUTINE gen_plm(nsmax, nlmax, nmmax, plm)
          IMPLICIT NONE

          integer(I4B), intent(IN)                     :: nsmax, nlmax, nmmax
          real(DP),     intent(OUT),  dimension(0:( nsmax*(nmmax+1)*(2*nlmax-nmmax+2) ),1:1 )   :: plm

          call plm_gen(nsmax,nlmax,nmmax,plm)

        END SUBROUTINE gen_plm


        !===============================================================================
        !===============================================================================


        SUBROUTINE logrange( x , mn , mx , npts )

          IMPLICIT NONE

          real(DP), dimension(1:npts) :: x
          real(DP) :: mn, mx, rg
          integer :: npts, i

          rg = (log10(mx)-log10(mn))/real(npts)

          do i = 1, npts
             x(i) = (10)**(real(i)*rg + log10(mn))
          enddo

          RETURN
        END SUBROUTINE logrange


        !===============================================================================
        !===============================================================================	


        SUBROUTINE gen_cln(cln, kln, nnmax, nlmax, rmax)

          IMPLICIT NONE

          INTEGER(I4B) :: nnmax, nlmax, l, n, k
          REAL(DP) :: rmax, tempval
          REAl(DP), DIMENSION(0:nlmax,1:nnmax) :: kln
          REAl(DP), DIMENSION(0:nlmax,1:nnmax) :: cln

          print*,"------------------------------------------------"
          print*,"   Computing normalization coefficients..."

          !$OMP PARALLEL &
          !$OMP SHARED(cln,rmax,kln,nnmax,nlmax) &
          !$OMP PRIVATE(l,n,k,tempval)	  

          !$OMP DO SCHEDULE(DYNAMIC,1)	 
          DO l = 0, nlmax 
             DO n = 1, nnmax   
         	k = n+1
        	call BJL( k , kln(l,n)*rmax, tempval )
        	cln(l,n) = ( sqrt(2.0)*rmax**(-1.5) ) / (kln(l,n)*tempval)
             ENDDO
          ENDDO
          !$OMP END DO
          !$OMP END PARALLEL   	

          print*,"   Done."
          print*,"------------------------------------------------"

        END SUBROUTINE gen_cln


        !===============================================================================
        !===============================================================================	


        SUBROUTINE gen_jln(jlns, kln, rho, nnmax, nlmax)

          IMPLICIT NONE

          INTEGER(I4B) :: nnmax,nlmax,n,l
          REAL(DP) :: rho, jln
          REAL(DP), dimension(0:nlmax,1:nnmax) :: kln
          REAL(DP), dimension(0:nlmax,1:nnmax) :: jlns      

          !$OMP PARALLEL &
          !$OMP SHARED(jlns,kln,rho,nnmax,nlmax) &
          !$OMP PRIVATE(jln,l,n)	  

          !$OMP DO SCHEDULE(RUNTIME)
          do n=1,nnmax   
             do l=0,nlmax
	    	!print*,n,l,OMP_GET_THREAD_NUM()
		CALL BJL( l , kln(l,n)*rho , jln )
		jlns(l,n) = kln(l,n) * jln
             enddo
          enddo
          !$OMP END DO
          !$OMP END PARALLEL

          return

        END SUBROUTINE gen_jln


        !===============================================================================
        !===============================================================================	


        FUNCTION dzrel( cplxnb, cplxnbref )

          complex(DP) :: cplxnb, cplxnbref, zero
          real(DP) :: nbmod, nbmodref, dzrel

          zero = 0.0_dp
          nbmod = dsqrt( (real( cplxnb ))**2.0_dp + &
               & (aimag( cplxnb ))**2.0_dp )
          nbmodref = dsqrt( (real( cplxnbref ))**2.0_dp + &
               & (aimag( cplxnbref ))**2.0_dp )

          !print*, nbmod, nbmodref, ( nbmod - nbmodref ), ( nbmod - nbmodref ) / nbmodref
          IF( nbmodref == zero  ) THEN
    	     dzrel = 0.0_dp
          ENDIF

          IF( ISNAN(nbmodref) .OR. ISNAN(nbmod)  ) THEN
    	     dzrel = 0.0_dp
          ELSE
    	     dzrel = abs( nbmod - nbmodref ) / nbmodref
          ENDIF

          IF( ISNAN(dzrel) .AND. nbmodref <= 0.000000000000001 ) THEN
    	     dzrel = 0.0_dp
          ENDIF

          RETURN 

        END FUNCTION dzrel


        !===============================================================================
        !===============================================================================	


        FUNCTION dzabs( cplxnb, cplxnbref )

          complex(DP) :: cplxnb, cplxnbref
          real(DP) :: nbmod, nbmodref, dzabs

          nbmod = dsqrt( (real( cplxnb ))**2.0_dp + &
               & (aimag( cplxnb ))**2.0_dp )
          nbmodref = dsqrt( (real( cplxnbref ))**2.0_dp + &
               & (aimag( cplxnbref ))**2.0_dp )

          !print*, nbmod, nbmodref, ( nbmod - nbmodref ), ( nbmod - nbmodref ) / nbmodref
          IF( ISNAN(nbmodref) ) THEN
    	     dzabs = 0.0_dp
          ELSE
    	     dzabs = abs( nbmod - nbmodref )
          ENDIF

          IF( ISNAN(dzabs) ) THEN
    	     dzabs = 0.0_dp
          ENDIF

          RETURN 

        END FUNCTION dzabs


        !===============================================================================
        !===============================================================================	

        !> Compute the value of the l-th order spherical bessel function at x
        !!@param[in] L : order 
        !!@param[in] X : abscissa
        !!@param[out] JL : result \f$ j_l(x) \f$
        SUBROUTINE BJL(L,X,JL)


          !!== MODIFIED SUBROUTINE FOR SPHERICAL BESSEL FUNCTIONS.                       ==!!
          !!== CORRECTED THE SMALL BUGS IN PACKAGE CMBFAST&CAMB(for l=4,5, x~0.001-0.002)==!! 
          !!== CORRECTED THE SIGN OF J_L(X) FOR X<0 CASE                                 ==!!
          !!== WORKS FASTER AND MORE ACCURATE FOR LOW L, X<<L, AND L<<X cases            ==!! 
          !!== zqhuang@astro.utoronto.ca                                                 ==!!
          IMPLICIT NONE
          INTEGER L
          real*8 X,JL
          real*8 AX,AX2
          real,PARAMETER::LN2=0.6931471805599453094D0
          real,PARAMETER::ONEMLN2=0.30685281944005469058277D0
          real,PARAMETER::PID2=1.5707963267948966192313217D0
          real,PARAMETER::PID4=0.78539816339744830961566084582D0
          real,parameter::ROOTPI12 = 21.269446210866192327578D0
          real,parameter::GAMMA1 =   2.6789385347077476336556D0 !/* Gamma function of 1/3 */
          real,parameter::GAMMA2 =   1.3541179394264004169452D0 !/* Gamma function of 2/3 */
          real,PARAMETER::PI=3.141592653589793238463D0
          real*8 NU,NU2,BETA,BETA2,COSB
          real*8 sx,sx2
          real*8 cotb,cot3b,cot6b,secb,sec2b
          real*8 trigarg,expterm,L3

          IF(L.LT.0)THEN
             write(*,*) 'Can not evaluate Spherical Bessel Function with index l<0'
             STOP
          ENDIF
          AX=DABS(X)
          AX2=AX**2
          IF(L.LT.7)THEN
             IF(L.EQ.0)THEN
                IF(AX.LT.1.D-1)THEN
                   JL=1.D0-AX2/6.D0*(1.D0-AX2/20.D0)
                ELSE
                   JL=DSIN(AX)/AX
                ENDIF

             ELSEIF(L.EQ.1)THEN
                IF(AX.LT.2.D-1)THEN
                   JL=AX/3.D0*(1.D0-AX2/10.D0*(1.D0-AX2/28.D0))
                ELSE
                   JL=(DSIN(AX)/AX-DCOS(AX))/AX
                ENDIF
             ELSEIF(L.EQ.2)THEN
                IF(AX.LT.3.D-1)THEN
                   JL=AX2/15.D0*(1.D0-AX2/14.D0*(1.D0-AX2/36.D0))
                ELSE
                   JL=(-3.0D0*DCOS(AX)/AX-DSIN(AX)*(1.D0-3.D0/AX2))/AX
                ENDIF
             ELSEIF(L.EQ.3)THEN
                IF(AX.LT.4.D-1)THEN
                   JL=AX*AX2/105.D0*(1.D0-AX2/18.D0*(1.D0-AX2/44.D0))
                ELSE
                   JL=(DCOS(AX)*(1.D0-15.D0/AX2)-DSIN(AX)*(6.D0-15.D0/AX2)/AX)/AX
                ENDIF
             ELSEIF(L.EQ.4)THEN
                IF(AX.LT.6.D-1)THEN
                   JL=AX2**2/945.D0*(1.D0-AX2/22.D0*(1.D0-AX2/52.D0))
                ELSE
                   JL=(DSIN(AX)*(1.D0-(45.D0-105.D0/AX2)/AX2)+DCOS(AX)*(10.D0-105.D0/AX2)/AX)/AX
                ENDIF
             ELSEIF(L.EQ.5)THEN
                IF(AX.LT.1.D0)THEN
                   JL=AX2**2*AX/10395.D0*(1.D0-AX2/26.D0*(1.D0-AX2/60.D0))
                ELSE
                   JL=(DSIN(AX)*(15.D0-(420.D0-945.D0/AX2)/AX2)/AX-DCOS(AX)*(1.D0-(105.D0-945.0d0/AX2)/AX2))/AX
                ENDIF
             ELSE
                IF(AX.LT.1.D0)THEN
                   JL=AX2**3/135135.D0*(1.D0-AX2/30.D0*(1.D0-AX2/68.D0))
                ELSE
                   JL=(DSIN(AX)*(-1.D0+(210.D0-(4725.D0-10395.D0/AX2)/AX2)/AX2)+ &
                        DCOS(AX)*(-21.D0+(1260.D0-10395.D0/AX2)/AX2)/AX)/AX
                ENDIF
             ENDIF
          ELSE
             NU=0.5D0+L
             NU2=NU**2
             IF(AX.LT.1.D-40)THEN
                JL=0.D0
             ELSEIF((AX2/L).LT.5.D-1)THEN
                JL=DEXP(L*DLOG(AX/NU)-LN2+NU*ONEMLN2-(1.D0-(1.D0-3.5D0/NU2)/NU2/30.D0)/12.D0/NU) &
                     /NU*(1.D0-AX2/(4.D0*NU+4.D0)*(1.D0-AX2/(8.D0*NU+16.D0)*(1.D0-AX2/(12.D0*NU+36.D0))))
             ELSEIF((real(L)**2/AX).LT.5.D-1)THEN
                BETA=AX-PID2*(L+1)
                JL=(DCOS(BETA)*(1.D0-(NU2-0.25D0)*(NU2-2.25D0)/8.D0/AX2*(1.D0-(NU2-6.25)*(NU2-12.25D0)/48.D0/AX2)) &
                     -DSIN(BETA)*(NU2-0.25D0)/2.D0/AX* (1.D0-(NU2-2.25D0)*(NU2-6.25D0)/24.D0/AX2*(1.D0-(NU2-12.25)* &
                     (NU2-20.25)/80.D0/AX2)) )/AX   
             ELSE
                L3=NU**0.325
                IF(AX .LT. NU-1.31*L3) then
                   COSB=NU/AX
                   SX = DSQRT(NU2-AX2)
                   COTB=NU/SX
                   SECB=AX/NU
                   BETA=DLOG(COSB+SX/AX)
                   COT3B=COTB**3
                   COT6B=COT3B**2
                   SEC2B=SECB**2
                   EXPTERM=( (2.D0+3.D0*SEC2B)*COT3B/24.D0 &
                        - ( (4.D0+SEC2B)*SEC2B*COT6B/16.D0 &
                        + ((16.D0-(1512.D0+(3654.D0+375.D0*SEC2B)*SEC2B)*SEC2B)*COT3B/5760.D0 &
                        + (32.D0+(288.D0+(232.D0+13.D0*SEC2B)*SEC2B)*SEC2B)*SEC2B*COT6B/128.D0/NU)*COT6B/NU) &
                        /NU)/NU
                   JL=DSQRT(COTB*COSB)/(2.D0*NU)*DEXP(-NU*BETA+NU/COTB-EXPTERM)

                   !          /**************** Region 2: x >> l ****************/

                ELSEIF (AX .GT. NU+1.48*L3) then
                   COSB=NU/AX
                   SX=DSQRT(AX2-NU2)
                   COTB=NU/SX
                   SECB=AX/NU
                   BETA=DACOS(COSB)
                   COT3B=COTB**3
                   COT6B=COT3B**2
                   SEC2B=SECB**2
                   TRIGARG=NU/COTB-NU*BETA-PID4 &
                        -((2.0+3.0*SEC2B)*COT3B/24.D0  &
                        +(16.D0-(1512.D0+(3654.D0+375.D0*SEC2B)*SEC2B)*SEC2B)*COT3B*COT6B/5760.D0/NU2)/NU
                   EXPTERM=( (4.D0+sec2b)*sec2b*cot6b/16.D0 &
                        -(32.D0+(288.D0+(232.D0+13.D0*SEC2B)*SEC2B)*SEC2B)*SEC2B*COT6B**2/128.D0/NU2)/NU2
                   JL=DSQRT(COTB*COSB)/NU*DEXP(-EXPTERM)*DCOS(TRIGARG)

                   !          /***************** Region 3: x near l ****************/

                ELSE
                   BETA=AX-NU
                   BETA2=BETA**2
                   SX=6.D0/AX
                   SX2=SX**2
                   SECB=SX**0.3333333333333333d0
                   SEC2B=SECB**2
                   JL=( GAMMA1*SECB + BETA*GAMMA2*SEC2B &
                        -(BETA2/18.D0-1.D0/45.D0)*BETA*SX*SECB*GAMMA1 &
                        -((BETA2-1.D0)*BETA2/36.D0+1.D0/420.D0)*SX*SEC2B*GAMMA2   &
                        +(((BETA2/1620.D0-7.D0/3240.D0)*BETA2+1.D0/648.D0)*BETA2-1.D0/8100.D0)*SX2*SECB*GAMMA1 &
                        +(((BETA2/4536.D0-1.D0/810.D0)*BETA2+19.D0/11340.D0)*BETA2-13.D0/28350.D0)*BETA*SX2*SEC2B*GAMMA2 &
                        -((((BETA2/349920.D0-1.D0/29160.D0)*BETA2+71.D0/583200.D0)*BETA2-121.D0/874800.D0)* &
                        BETA2+7939.D0/224532000.D0)*BETA*SX2*SX*SECB*GAMMA1)*DSQRT(SX)/ROOTPI12
                ENDIF
             ENDIF
          ENDIF
          IF(X.LT.0.AND.MOD(L,2).NE.0)JL=-JL


        END SUBROUTINE BJL


        !===============================================================================	
        !===============================================================================


      END MODULE f3DEX
