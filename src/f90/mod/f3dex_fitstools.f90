MODULE f3dex_fitstools

  USE OMP_LIB

  USE healpix_types
  USE healpix_modules
  USE alm_tools
  USE fitstools
  IMPLICIT NONE
  
  REAL(kind=SP),     private, parameter :: s_bad_value = HPX_SBADVAL
  REAL(kind=DP),     private, parameter :: d_bad_value = HPX_DBADVAL

CONTAINS

  ! ---------------------------------------------------------------------------------------!     


  !> Write power spectrum to file
  !!@param[in] filename : fits file for output
  !!@param[in] cln : power spectrum
  !!@param[in] kln : k spectrum
  !!@param[in] (nlmax,nnmax) : n-l limits
  !!@param[in] (header,nlheader) : headers
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
    LOGICAL(LGT) ::  simple,extEND
    CHARACTER(LEN=80) :: comment

    INTEGER(I4B), PARAMETER :: maxdim = 40 !number of columns in the extension
    INTEGER(I4B) :: nrows, npix, tfields, varidat, repeat
    INTEGER(I4B) :: frow,  felem, colnum, stride, istart, iEND, k
    CHARACTER(LEN=20) :: ttype(maxdim), tform(maxdim), tunit(maxdim), extname
    CHARACTER(LEN=10) ::  card
    CHARACTER(LEN=2) :: stn
    INTEGER(I4B) :: itn, cnt, ncl
    INTEGER(I4B) :: l, n
    character(len=filenamelen) sfilename
    character(len=1) :: pform

    status = 0
    unit = 140

    CALL ftinit(unit,filename,blocksize,status)

    simple=.TRUE.
    bitpix=32     ! INTEGER*4
    naxis=0       ! no image
    naxes(1)=0
    extEND=.TRUE. ! there is an extension

    ncl = 4

    ! primary header
    CALL ftphpr(unit,simple,bitpix,naxis,naxes,0_i4b,1_i4b,extEND,status)

    !     creates an extension
    CALL ftcrhd(unit, status)
    !     writes required keywords
    nrows    = nnmax * nlmax  ! naxis1
    tfields  = ncl
    repeat   = 1
    DO i=1, ncl
       ttype(i) = ''
       tform(i) = ''
       tform(i)(1:2) = '1D'
    ENDDO

    ttype(1)(1:9) = 'l index  '
    ttype(2)(1:9) = 'n index  '
    ttype(3)(1:9) = 'kln value'
    ttype(4)(1:9) = 'cln(l,n) '
    tunit(1:ncl) = ''      ! optional, will not appear
    extname  = ''      ! optional, will not appear
    varidat  = 0
    CALL ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
         &     extname, varidat, status)
    !     write the header literally, putting TFORMi at the desired place
    DO i=1,nlheader
       card = header(i)
       IF (card(1:5) == 'TTYPE') THEN ! IF TTYPEi is explicitely given
          stn = card(6:7)
          read(stn,'(i2)') itn
          ! discard at their original location:
          CALL ftdkey(unit,'TTYPE'//stn,status)  ! old TTYPEi and  ! remove
          status = 0
          CALL ftdkey(unit,'TFORM'//stn,status)  !     TFORMi
          status = 0
          IF (itn <= tfields) THEN ! only put relevant information 2008-08-27
             CALL putrec(unit,header(i), status)           ! write new TTYPE1
             status = 0
             comment = ''
             IF (itn==1) THEN
                comment = 'data format of field: 4-byte INTEGER'
             ELSE
                IF (DP == SP) comment = 'data format of field: 4-byte REAL'
                IF (DP == DP) comment = 'data format of field: 8-byte REAL'
             ENDIF
             CALL ftpkys(unit,'TFORM'//stn,tform(itn),comment,status) ! and write new TFORM1 right after
          ENDIF
       ELSEIF (header(i)/=' ') THEN
          CALL putrec(unit,header(i), status)
       ENDIF
       status = 0
    ENDDO
    CALL ftukyj(unit, 'MAX-LPOL', nlmax, 'Maximum L multipole order',  status)
    CALL ftukyj(unit, 'MAX-NPOL', nnmax, 'Maximum N multipole order', status)
    !     write the extension by blocks of rows ! EH, Dec 2004
    felem  = 1  ! starting position (element)
    CALL ftgrsz(unit, stride, status) ! find optimal stride in rows
    stride = max( stride, 1)
    ALLOCATE(lindex(0:stride-1))
    ALLOCATE(nindex(0:stride-1))
    ALLOCATE(klntemp(0:stride-1))
    ALLOCATE(clntemp(0:stride-1))
    lindex = 0.0
    nindex = 0.0
    klntemp = 0.0
    clntemp = 0.0

    cnt = 0
    istart = 1
    DO l = 0, nlmax
       DO n = 1, nnmax
          IF( cnt < stride ) THEN
             lindex(cnt) = REAL(l)
             nindex(cnt) = REAL(n)
             klntemp(cnt) = kln(l,n)
             clntemp(cnt) = cln(l,n)
          ENDIF
          IF( cnt >= stride .or. cnt >= nlmax*nnmax ) THEN
             !PRINT*,lindex,nindex
             CALL f90ftpcld(unit, 1_i4b, istart, felem, cnt, lindex, status)
             CALL f90ftpcld(unit, 2_i4b, istart, felem, cnt, nindex, status)
             CALL f90ftpcld(unit, 3_i4b, istart, felem, cnt, klntemp, status)
             CALL f90ftpcld(unit, 4_i4b, istart, felem, cnt, clntemp, status)
             istart = istart + cnt - 1
             cnt = 1
             lindex = 0.0
             nindex = 0.0
             klntemp = 0.0
             clntemp = 0.0
          ENDIF
          cnt = cnt + 1
       ENDDO
    ENDDO

    deALLOCATE(lindex,nindex,klntemp,clntemp)

    !     close the file and free the unit number
    CALL ftclos(unit, status)

    !     check for any error, and IF so PRINT out error messages
    IF (status > 0) CALL PRINTerror(status)

    RETURN

  END SUBROUTINE cln2fits

  ! ---------------------------------------------------------------------------------------!     

  !> Write a_lmn's to file
  !!@param[in] filename : fits file for output
  !!@param[in] almn : Fourier-bessel decomposition
  !!@param[in] kln : k spectrum
  !!@param[in] cln : normalization
  !!@param[in] (nlmax,nmmax,nnmax) : n-m-l limits
  !!@param[in] (header,nlheader) : headers
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
    LOGICAL(LGT) ::  simple,extEND
    CHARACTER(LEN=80) :: comment

    INTEGER(I4B), PARAMETER :: maxdim = 40 !number of columns in the extension
    INTEGER(I4B) :: nrows, npix, tfields, varidat, repeat, frow
    INTEGER(I4B) :: felem, colnum, stride, istart, iEND, k
    CHARACTER(LEN=20) :: ttype(maxdim), tform(maxdim), tunit(maxdim), extname
    CHARACTER(LEN=10) ::  card
    CHARACTER(LEN=2) :: stn
    INTEGER(I4B) :: cnt, ncl
    INTEGER(I4B) :: l, n, m, itn
    character(len=filenamelen) sfilename
    character(len=1) :: pform

    status = 0
    unit = 139

    CALL ftinit(unit,filename,blocksize,status)

    !     -----------------------------------------------------
    !     initialize parameters about the FITS image
    simple=.TRUE.
    bitpix=32     ! INTEGER*4
    naxis=0       ! no image
    naxes(1)=0
    extEND=.TRUE. ! there is an extension

    ncl = 7

    ! primary header
    CALL ftphpr(unit,simple,bitpix,naxis,naxes,0_i4b,1_i4b,extEND,status)

    !     creates an extension
    CALL ftcrhd(unit, status)
    !     writes required keywords
    nrows    = nnmax * (nlmax+1) * (1+nmmax/2)  ! naxis1
    tfields  = ncl
    repeat   = 1
    DO i=1, ncl
       ttype(i) = ''
       tform(i) = ''
       tform(i)(1:2) = '1D'
    ENDDO
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
    CALL ftphbn(unit, nrows, tfields, ttype, tform, tunit, &
	 &     extname, varidat, status)
    !     write the header literally, putting TFORMi at the desired place
    DO i=1,nlheader
       card = header(i)
       IF (card(1:5) == 'TTYPE') THEN ! IF TTYPEi is explicitely given
	  stn = card(6:7)
	  read(stn,'(i2)') itn
   ! discard at their original location:
	  CALL ftdkey(unit,'TTYPE'//stn,status)  ! old TTYPEi and  ! remove
	  status = 0
	  CALL ftdkey(unit,'TFORM'//stn,status)  !     TFORMi
	  status = 0
	  IF (itn <= tfields) THEN ! only put relevant information 2008-08-27
	     CALL putrec(unit,header(i), status)           ! write new TTYPE1
	     status = 0
	     comment = ''
	     IF (itn==1) THEN
		comment = 'data format of field: 4-byte INTEGER'
	     ELSE
		IF (DP == SP) comment = 'data format of field: 4-byte REAL'
		IF (DP == DP) comment = 'data format of field: 8-byte REAL'
	     ENDIF
	     CALL ftpkys(unit,'TFORM'//stn,tform(itn),comment,status) ! and write new TFORM1 right after
	  ENDIF
       ELSEIF (header(i)/=' ') THEN
	  CALL putrec(unit,header(i), status)
       ENDIF
       status = 0
    ENDDO
    CALL ftukyj(unit, 'MAX-LPOL', nlmax, 'Maximum L multipole order', status)
    CALL ftukyj(unit, 'MAX-MPOL', nmmax, 'Maximum M multipole order', status)
    CALL ftukyj(unit, 'MAX-NPOL', nnmax, 'Maximum N multipole order', status)
    !     write the extension by blocks of rows ! EH, Dec 2004
    felem  = 1  ! starting position (element)
    CALL ftgrsz(unit, stride, status) ! find optimal stride in rows
    stride = max( stride, 1)

    ALLOCATE(lindex(0:stride-1))
    ALLOCATE(mindex(0:stride-1))
    ALLOCATE(nindex(0:stride-1))
    ALLOCATE(klntemp(0:stride-1))
    ALLOCATE(clntemp(0:stride-1))
    ALLOCATE(almntemp1(0:stride-1))
    ALLOCATE(almntemp2(0:stride-1))

    lindex = 0.0
    mindex = 0.0
    nindex = 0.0
    klntemp = 0.0
    clntemp = 0.0
    almntemp1 = 0.0
    almntemp2 = 0.0

    cnt = 0
    istart = 1
    DO l = 0, nlmax
       DO m = 0, l
          DO n = 1, nnmax
             lindex(cnt) = REAL(l,kind=DP)
             mindex(cnt) = REAL(m,kind=DP)
             nindex(cnt) = REAL(n,kind=DP)
             klntemp(cnt) = REAL(kln(l,n),kind=DP)
             clntemp(cnt) = REAL(cln(l,n),kind=DP)
             almntemp1(cnt) = REAL(almn(n,l,m))
             almntemp2(cnt) = aimag(almn(n,l,m))
             cnt = cnt + 1
             !PRINT*,l,m,n,almn(n,l,m)
             IF( cnt >= stride .or. cnt == nrows ) THEN
                CALL f90ftpcld(unit, 1_i4b, istart, felem, cnt, lindex, status)
                CALL f90ftpcld(unit, 2_i4b, istart, felem, cnt, mindex, status)
                CALL f90ftpcld(unit, 3_i4b, istart, felem, cnt, nindex, status)
                CALL f90ftpcld(unit, 4_i4b, istart, felem, cnt, klntemp, status)
                CALL f90ftpcld(unit, 5_i4b, istart, felem, cnt, clntemp, status)
                CALL f90ftpcld(unit, 6_i4b, istart, felem, cnt, almntemp1, status)
                CALL f90ftpcld(unit, 7_i4b, istart, felem, cnt, almntemp2, status)
                istart = istart + cnt - 1
                cnt = 0
                lindex = 0.0
                mindex = 0.0
                nindex = 0.0
                klntemp = 0.0
                clntemp = 0.0
                almntemp1 = 0.0
                almntemp2 = 0.0
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    deALLOCATE(lindex,mindex,nindex,klntemp,clntemp,almntemp1,almntemp2)

    !     ----------------------
    !     close and exit
    !     ----------------------

    !     close the file and free the unit number
    CALL ftclos(unit, status)

    !     check for any error, and IF so PRINT out error messages
    IF (status > 0) CALL PRINTerror(status)

    RETURN

  END SUBROUTINE almn2fits

  ! ---------------------------------------------------------------------------------------!     

  !> Extracts almn's from fits file
  !!@param[in] filename : fits file
  !!@param[out] almn : Fourier-bessel decomposition
  !!@param[out] kln : k spectrum
  !!@param[out] cln : normalization
  !!@param[in] (nlmax,nmmax,nnmax) : n-m-l limits
  !!@param[in] (header,nlheader) : headers
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
    LOGICAL(LGT) :: extEND
    INTEGER(I4B) :: nmove, hdutype ! , nkeys , nspace
    INTEGER(I4B) :: frow, imap
    INTEGER(I4B) :: datacode, repeat, width
    INTEGER(I4B) :: i, l, m, n
    INTEGER(i4b) :: nrow2read, nelem
    INTEGER(i4b) :: i0, i1

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
    anynull = .FALSE.
    almn=0.
    readwrite=0
    CALL ftopen(unit,filename,readwrite,blocksize,status)
    IF (status > 0) THEN 
       CALL PRINTerror(status)
       CALL fatal_error("Aborting.")
    ENDIF
    !     -----------------------------------------

    !     determines the presence of image
    CALL ftgkyj(unit,'NAXIS', naxis, comment, status)
    IF (status > 0) CALL PRINTerror(status)

    !     determines the presence of an extension
    CALL ftgkyl(unit,'EXTEND', extEND, comment, status)
    IF (status > 0) status = 0 ! no extension : 
    !     to be compatible with first version of the code

    CALL assert (extEND, 'No extension!')
    nmove = +extno
    CALL ftmrhd(unit, nmove, hdutype, status)
    !cc         write(*,*) hdutype

    CALL assert(hdutype==2, 'this is not a binary table')

    header = ""
    CALL get_clean_header( unit, header, filename, status)

    !        reads all the keywords
    CALL ftghbn(unit, maxdim, &
	 &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
	 &        status)

    IF (tfields<ncl) THEN
       PRINT *,'found ',tfields,' columns in the file'
       PRINT *,'expected ',ncl
       CALL fatal_error
    ENDIF
    CALL f90ftgkyd(unit, 'BAD_DATA', nullval, comment, status)
    IF (status == 202) THEN ! bad_data not found
       IF (DP == SP) nullval = s_bad_value ! default value
       IF (DP == DP) nullval = d_bad_value ! default value
       status = 0
    ENDIF
    !parse TFORM keyword to find out the length of the column vector
    CALL ftbnfm(tform(1), datacode, repeat, width, status)
    npix = nrows * repeat
    !    IF (npix /= nalms) THEN
    !IF (npix > nalms) THEN
    !PRINT *,'found ',npix,' alms'
    !!       PRINT *,'expected ',nalms
    !PRINT *,'expected ',nalms,' or less'
    !CALL fatal_error
    !ENDIF

    CALL ftgrsz(unit, nrow2read, status)
    nrow2read = max(nrow2read, 1)
    nelem = nrow2read * repeat
    i0 = 0_i8b

    ALLOCATE(tab(0:nelem-1,1:ncl))

    ALLOCATE(lindex(0:nelem-1))
    ALLOCATE(mindex(0:nelem-1))
    ALLOCATE(nindex(0:nelem-1))
    ALLOCATE(klntemp(0:nelem-1))
    ALLOCATE(clntemp(0:nelem-1))
    ALLOCATE(almntemp1(0:nelem-1))
    ALLOCATE(almntemp2(0:nelem-1))

    lindex = 0.0
    mindex = 0.0
    nindex = 0.0
    klntemp = 0.0
    clntemp = 0.0
    almntemp1 = 0.0
    almntemp2 = 0.0

    DO frow = 1, nrows, nrow2read
       i1 = min(i0 + nrow2read * repeat, int(npix,i8b)) - 1_i8b
       nelem = i1 - i0 + 1
       !DO imap = 1, ncl
       !   CALL f90ftgcvd(unit, imap, frow, 1_i4b, nelem, nullval, &
       !   &        tab(0:nelem-1,imap), anynull, status)
       !   CALL assert (.not. anynull, 'There are undefined values in the table!')
       !ENDDO

       CALL f90ftgcvd(unit, 1_i4b, frow, 1_i4b, nelem, nullval, &
            &        lindex(0:nelem-1), anynull, status)
       status = 0
       CALL f90ftgcvd(unit, 2_i4b, frow, 1_i4b, nelem, nullval, &
            &        mindex(0:nelem-1), anynull, status)
       status = 0
       CALL f90ftgcvd(unit, 3_i4b, frow, 1_i4b, nelem, nullval, &
            &        nindex(0:nelem-1), anynull, status)
       status = 0
       CALL f90ftgcvd(unit, 4_i4b, frow, 1_i4b, nelem, nullval, &
            &        klntemp(0:nelem-1), anynull, status)
       status = 0
       CALL f90ftgcvd(unit, 5_i4b, frow, 1_i4b, nelem, nullval, &
            &        clntemp(0:nelem-1), anynull, status)
       status = 0
       CALL f90ftgcvd(unit, 6_i4b, frow, 1_i4b, nelem, nullval, &
            &        almntemp1(0:nelem-1), anynull, status)
       status = 0
       CALL f90ftgcvd(unit, 7_i4b, frow, 1_i4b, nelem, nullval, &
            &        almntemp2(0:nelem-1), anynull, status) 
       status = 0
       CALL assert (.not. anynull, 'There are undefined values in the table!')

       DO i = 0, nelem-1
          !PRINT*,"entering writing"
          l = int(lindex(i),i4b)
          m = int(mindex(i),i4b)
          n = int(nindex(i),i4b)
          IF( l >= 0 .and. l <= nlmax &
               & .and. m >= 0 .and. m <= l &
               & .and. n >= 1 .and. n <= nnmax ) THEN
             kln(l,n) = klntemp(i)
             cln(l,n) = clntemp(i)
             almn(n,l,m) = CMPLX( almntemp1(i), almntemp2(i), KIND=DP)
110          FORMAT(A,I3,A,I3,A,I3,A) 
             !WRITE(nmtemp,110),"(",l,",",m,",",n,")"
             !PRINT*,nmtemp,almn(n,l,m)
          ENDIF
       ENDDO

       i0 = i1 + 1_i8b
       !PRINT*,i0,npix
    ENDDO

    deALLOCATE(lindex,mindex,nindex,klntemp,clntemp,almntemp1,almntemp2)
    deALLOCATE(tab)

    ! sanity check
    IF (i0 /= npix) THEN
       PRINT*,'something wrong during piece-wise reading'
       CALL fatal_error
    ENDIF

    !     close the file
    CALL ftclos(unit, status)


    !     check for any error, and IF so PRINT out error messages
    IF (status > 0) CALL PRINTerror(status)
    RETURN

  END SUBROUTINE fits2almn
  
        ! ---------------------------------------------------------------------------------------!    

        SUBROUTINE f90ftpcld(unit, colnum, frow, felem, np, data, status)
          INTEGER(I4B), intent(in)  :: unit, colnum, frow, felem, np
          INTEGER(I4B), intent(out) :: status
          REAL(DP),     intent(in), dimension(0:)  :: data
          CALL ftpcld(unit, colnum, frow, felem, np, data, status)
          RETURN
        END SUBROUTINE f90ftpcld

        ! ---------------------------------------------------------------------------------------!    

        SUBROUTINE f90ftgcvd(unit, colnum, frow, felem, np, nullval, data, anynull, status)
          INTEGER(I4B), intent(in)  :: unit, colnum, frow, felem, np
          INTEGER(I4B), intent(out) :: status
          LOGICAL(LGT), intent(out) :: anynull
          REAL(DP), intent(out), dimension(0:) :: data
          REAL(DP), intent(in) :: nullval
          CALL ftgcvd(unit, colnum, frow, felem, np, nullval, data, anynull, status)
          RETURN
        END SUBROUTINE f90ftgcvd

        ! ---------------------------------------------------------------------------------------!    

        SUBROUTINE f90ftgkyd(unit, keyword, value, comment, status)
          INTEGER(I4B),     intent(in)  :: unit
          character(len=*), intent(in)  :: keyword
          INTEGER(I4B),     intent(out) :: status
          character(len=*), intent(out) :: comment
          REAL(DP),         intent(out) :: value
          CALL ftgkyd(unit, keyword, value, comment, status)
          RETURN
        END SUBROUTINE f90ftgkyd

        ! ---------------------------------------------------------------------------------------!    

        SUBROUTINE get_clean_header(unit, header, filename, error, xalso, xonly)

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

          IF (present(xonly)) THEN 
             n_excl = size(xonly)
             ALLOCATE(to_excl(1:n_excl))
             to_excl = xonly

          ELSE IF (present(xalso)) THEN
             n_excl = size(xalso) + size(def_excl)
             ALLOCATE(to_excl(1:n_excl))
             to_excl(1:size(def_excl)) = def_excl
             to_excl(size(def_excl)+1:n_excl) = xalso

          ELSE
             n_excl = size(def_excl)
             ALLOCATE(to_excl(1:n_excl))
             to_excl = def_excl
          ENDIF

          nlheader=size(header)
          ! go to END of fortran header
          DO i = 1, nlheader
             IF (trim(header(i)) == "") exit
          ENDDO
          ! go to top of fits file header
          status=0
          CALL ftgrec(unit,0_i4b,record,status)
          ! read in all header lines except those excluded
          DO
             CALL ftgnxk(unit,'*',1_i4b,to_excl,n_excl,record,status)
             IF (status > 0) exit ! END of header
             IF (i > nlheader) THEN
                write(unit=*,fmt="(a,i5,a)") &
                     & " WARNING : The header in "//  &
                     &    trim(filename)//" has more than ", &
                     &  nlheader," lines."
                PRINT*," It will be truncated."
                error = 1
                exit
             ENDIF
             header(i)=record
             i=i+1
          ENDDO
          status=0

          RETURN
        END SUBROUTINE get_clean_header

  ! ---------------------------------------------------------------------------------------!       

  !> Write bi-tab to file
  !!@param[in] file : fits file for output
  !!@param[in] tab : bi-array to write
  !!@param[in] (dim1start, dim1en) : first dim bounds
  !!@param[in] (dim2start, dim2en) : first dim bounds
   SUBROUTINE bitab2fits( file, tab, dim1start, dim1end, dim2start, dim2end )
   
   CHARACTER(len=FILENAMELEN) ::  file
   CHARACTER(len=80), DIMENSION(1:120) :: header
   INTEGER(I4B) :: status,  dim1start, dim1end, dim2start, dim2end
   INTEGER(I8B) :: d1, d2, npix, loc, pix
   REAL(DP), DIMENSION(dim1start:dim1end,dim2start:dim2end) :: tab
   REAL(DP), DIMENSION(:,:), ALLOCATABLE :: temptab
   CHARACTER(len=*), PARAMETER :: code = "bitab2fits"
  
   npix = (dim1end-dim1start+1)*(dim2end-dim2start+1)
   header(:) = ' '
   loc = 0
   pix = 0
   
   ALLOCATE( temptab(0:npix-1,1:3),stat = status )
   CALL assert_alloc(status,code,'temptab')
   
   DO d1=dim1start, dim1end
      DO d2=dim2start, dim2end
         temptab(pix,1) = REAL(d1)
         temptab(pix,2) = REAL(d2)
         temptab(pix,3) = tab(d1,d2)
         pix = pix + 1
      ENDDO
   ENDDO
   CALL write_bintabh(temptab, npix, 3, header, 120, file, firstpix=loc)
   
   DEALLOCATE(temptab)

   END SUBROUTINE bitab2fits
      
  ! ---------------------------------------------------------------------------------------!       

  !> Read bi-tab from file
  !!@param[in] file : input fits file
  !!@param[in] tab : bi-array to write
  !!@param[in] (dim1start, dim1en) : first dim bounds
  !!@param[in] (dim2start, dim2en) : first dim bounds
   SUBROUTINE fits2bitab( file, tab, dim1start, dim1end, dim2start, dim2end )
   
   CHARACTER(len=FILENAMELEN) ::  file
   CHARACTER(len=80), DIMENSION(1:120) :: header
   INTEGER(I8B) :: d1, d2, loc, npix, c1, c2, pix
   INTEGER(I4B) :: status,  dim1start, dim1end, dim2start, dim2end
   REAL(DP), DIMENSION(dim1start:dim1end,dim2start:dim2end) :: tab
   REAL(DP), DIMENSION(:,:), ALLOCATABLE :: temptab
   REAL(DP) :: c3
   CHARACTER(len=*), PARAMETER :: code = "fits2bitab"
   
   npix = (dim1end-dim1start+1)*(dim2end-dim2start+1)
   header(:) = ' '
   loc = 0
   tab = 0.0
   pix = 0
   
   ALLOCATE( temptab(0:npix-1,1:3),stat = status )
   CALL assert_alloc(status,code,'temptab')
   
   CALL input_tod(file, temptab, npix, 3, firstpix=loc)
   DO d1=dim1start, dim1end
      DO d2=dim2start, dim2end
         c1 = INT(temptab(pix,1))
         c2 = INT(temptab(pix,2))
         c3 = temptab(pix,3)
         IF( (c1 .GE. 0.0) .AND. (c2 .GE. 0.0) .AND. &
           & (c1 .GE. dim1start) .AND. (c1 .LE. dim1end) .AND. & 
           & (c2 .GE. dim2start) .AND. (c2 .LE. dim2end) ) THEN
             tab(c1,c2) = c3
         ENDIF
         pix = pix + 1
      ENDDO
   ENDDO
   
   DEALLOCATE(temptab)
   
   END SUBROUTINE fits2bitab
   
  ! ---------------------------------------------------------------------------------------!       

  !> Write tri-tab to file
  !!@param[in] file : fits file for output
  !!@param[in] tab : tri-array to write
  !!@param[in] (dim1start, dim1en) : first dim bounds
  !!@param[in] (dim2start, dim2en) : first dim bounds
  !!@param[in] (dim3start, dim3en) : first dim bounds
   SUBROUTINE tritab2fits( file, tab, dim1start, dim1end, dim2start, dim2end, dim3start, dim3end )
   
   CHARACTER(len=FILENAMELEN) ::  file
   CHARACTER(len=80), DIMENSION(1:2) :: header
   INTEGER(I4B) :: status,  dim1start, dim1end, dim2start, dim2end, dim3start, dim3end
   INTEGER(I8B) :: d1, d2, d3, npix, loc, pix
   REAL(DP), DIMENSION(dim1start:dim1end,dim2start:dim2end,dim3start:dim3end) :: tab
   REAL(DP), DIMENSION(:,:), ALLOCATABLE :: temptab
   CHARACTER(len=*), PARAMETER :: code = "tritab2fits"
  
   npix = (dim2end-dim2start+1)*(dim3end-dim3start+1)
   header(:) = ' '
   loc = 0
   
   ALLOCATE( temptab(0:npix-1,1:4),stat = status )
   CALL assert_alloc(status,code,'temptab')
   
   DO d1=dim1start, dim1end
      temptab=0.0
      pix = 0
      DO d2=dim2start, dim2end
      	 DO d3=dim3start, dim3end
            temptab(pix,1) = REAL(d1)
            temptab(pix,2) = REAL(d2)
            temptab(pix,3) = REAL(d3)
            temptab(pix,4) = tab(d1,d2,d3)
            pix = pix + 1
         ENDDO
      ENDDO
      CALL write_bintabh(temptab, npix, 4, header, 2, file, firstpix=loc)
      loc = loc + npix
   ENDDO
   
   DEALLOCATE(temptab)

   END SUBROUTINE tritab2fits
      
  ! ---------------------------------------------------------------------------------------!       

  !> Read tri-tab from file
  !!@param[in] file : input fits file
  !!@param[in] tab : tri-array to write
  !!@param[in] (dim1start, dim1en) : first dim bounds
  !!@param[in] (dim2start, dim2en) : first dim bounds
  !!@param[in] (dim3start, dim3en) : first dim bounds
   SUBROUTINE fits2tritab( file, tab, dim1start, dim1end, dim2start, dim2end, dim3start, dim3end )
   
   CHARACTER(len=FILENAMELEN) ::  file
   CHARACTER(len=80), DIMENSION(1:2) :: header
   INTEGER(I8B) :: d1, d2, d3, loc, npix, c1, c2, c3, pix
   INTEGER(I4B) :: status,  dim1start, dim1end, dim2start, dim2end, dim3start, dim3end
   REAL(DP), DIMENSION(dim1start:dim1end,dim2start:dim2end,dim3start:dim3end) :: tab
   REAL(DP), DIMENSION(:,:), ALLOCATABLE :: temptab
   REAL(DP) :: c4
   CHARACTER(len=*), PARAMETER :: code = "fits2tritab"
   
   npix = (dim2end-dim2start+1)*(dim3end-dim3start+1)
   header(:) = ' '
   loc = 0
   tab = 0.0
   
   ALLOCATE( temptab(0:npix-1,1:4),stat = status )
   CALL assert_alloc(status,code,'temptab')
   
   DO d1=dim1start, dim1end
      temptab = -1.0
      CALL input_tod(file, temptab, npix, 4, firstpix=loc)
      pix = 0
      DO d2=dim2start, dim2end
         DO d3=dim3start, dim3end
            c1 = INT(temptab(pix,1))
            c2 = INT(temptab(pix,2))
            c3 = INT(temptab(pix,3))
            c4 = REAL(temptab(pix,4))
            IF( (c1 .GE. 0.0) .AND. (c2 .GE. 0.0) .AND. &
              & (c1 .GE. dim1start) .AND. (c1 .LE. dim1end) .AND. & 
              & (c2 .GE. dim2start) .AND. (c2 .LE. dim2end) .AND. & 
              & (c3 .GE. dim3start) .AND. (c3 .LE. dim3end) ) THEN
                tab(c1,c2,c3) = c4
            ENDIF
            pix = pix + 1
         ENDDO
      ENDDO
      loc = loc + npix
   ENDDO
   
   DEALLOCATE(temptab)
   
   END SUBROUTINE fits2tritab 
      
  ! ---------------------------------------------------------------------------------------!       
      
END MODULE f3dex_fitstools
