MODULE f3dex_transforms

  USE OMP_LIB

  USE f3dex_utils
  USE f3dex_stats

  USE healpix_types
  USE healpix_modules
  USE alm_tools
  USE fitstools
  IMPLICIT NONE

  REAL(kind=SP),     private, parameter :: s_bad_value = HPX_SBADVAL
  REAL(kind=DP),     private, parameter :: d_bad_value = HPX_DBADVAL

CONTAINS

  ! ---------------------------------------------------------------------------------------!    

!!$    
!!$    SUBROUTINE almn2alnspring_pre( nsmax, nlmax, nmmax, map, almn, plm )
!!$
!!$      INTEGER(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
!!$      complex(DP), intent(IN),  dimension(1:1,0:nlmax,0:nmmax) :: almn
!!$      REAL(DP),   intent(OUT), dimension(0:nlmax,0:(12_i8b*nsmax)*nsmax-1) :: map
!!$      REAL(DP),     intent(IN),  dimension(0:)                  :: plm
!!$
!!$      INTEGER(I4B) :: l, m, ith
!!$      INTEGER(I8B) :: istart_south, istart_north, npix
!!$      INTEGER(I4B) :: nrings, nphmx
!!$      INTEGER(i4b), parameter :: SMAXCHK = 30
!!$
!!$      INTEGER(I8B) :: n_lm, n_plm, i_mm
!!$      complex(DPC), dimension(-1:1)             :: b_ns
!!$      REAL(DP),     dimension(:,:), allocatable :: dalm
!!$      INTEGER(i4b)                              :: ll, l_min, l_start
!!$      REAL(DP)                                  :: cth
!!$      complex(DPC), dimension(:,:), allocatable :: b_north, b_south
!!$      REAL(DP),     dimension(:),   allocatable :: ring
!!$      INTEGER(i4b)             :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
!!$      INTEGER(I8B), dimension(0:SMAXCHK-1) :: startpix
!!$      INTEGER(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
!!$      REAL(DP),     dimension(0:SMAXCHK-1) :: sth
!!$
!!$      character(LEN=*), parameter :: code = 'almn2alnspring_pre'
!!$      INTEGER(I4B) :: status
!!$      !=======================================================================
!!$
!!$      ! Healpix definitions
!!$      nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
!!$      npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
!!$      nphmx  = 4*nsmax           ! maximum number of pixels/ring
!!$      n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]
!!$      n_plm  = n_lm * nrings
!!$
!!$      !     --- ALLOCATEs space for arrays ---
!!$      nchunks   = nrings/SMAXCHK + 1  ! number of chunks
!!$      chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK
!!$
!!$      ALLOCATE(b_north(0:nmmax, 0:chunksize-1),stat = status)
!!$      CALL assert_alloc(status,code,'b_north')
!!$
!!$      ALLOCATE(b_south(0:nmmax, 0:chunksize-1),stat = status)
!!$      CALL assert_alloc(status,code,'b_south')
!!$
!!$      IF (.not. DO_openmp()) THEN
!!$         ALLOCATE(dalm(0:1,0:nmmax), stat = status)
!!$         CALL assert_alloc(status,code,'dalm')
!!$         ALLOCATE(ring(0:nphmx-1),stat = status)
!!$         CALL assert_alloc(status,code,'ring')
!!$      ENDIF
!!$      !     ------------ initiate variables and arrays ----------------
!!$
!!$      map = 0.0 ! set the whole map to zero
!!$
!!$      ! loop on chunks
!!$      DO ichunk = 0, nchunks-1
!!$         lchk = ichunk * chunksize + 1
!!$         uchk = min(lchk+chunksize - 1, nrings)
!!$
!!$         DO ith = lchk, uchk
!!$            ithl = ith - lchk !local index
!!$            ! get pixel location information
!!$            CALL get_pixel_layout(nsmax, ith, cth, sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
!!$         ENDDO
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
!!$         IF (DO_openmp()) THEN
!!$            ALLOCATE(dalm(0:1,0:nmmax), stat = status)
!!$            CALL assert_alloc(status,code,'dalm')
!!$         ENDIF
!!$
!!$         !$OMP DO schedule(dynamic,1)
!!$
!!$         DO l = 0, nlmax
!!$
!!$            DO m = 0, l
!!$               dalm(0,m) =  REAL(almn(1,l,m),kind=dp)
!!$               dalm(1,m) = aimag(almn(1,l,m))
!!$            ENDDO ! loop on m 
!!$
!!$            DO ithl = 0, uchk - lchk
!!$
!!$               DO m = 0, l
!!$
!!$                  !l_min = l_min_ylm(m, sth(ithl))
!!$                  !IF (nlmax >= l_min) THEN ! skip calculations when Ylm too small
!!$
!!$                  ith = ithl + lchk
!!$                  i_mm = n_lm * (ith-1) + ((2_I8B*nlmax + 3 - m)*m)/2 
!!$                  ! location of Ym,m for ring ith
!!$
!!$                  b_ns = 0.0_dpc
!!$
!!$                  ! odd values of (l+m)
!!$                  IF (mod(m+l,2) == 0) THEN
!!$                     b_ns(-1) = cmplx(plm(i_mm+l-m) * dalm(0,m), plm(i_mm+l-m)*dalm(1,m), kind=DP)
!!$                  ENDIF
!!$
!!$                  ! even values of (l+m)
!!$                  IF (mod(m+l,2) == 1) THEN
!!$                     b_ns(1)  =  cmplx(plm(i_mm+l-m) * dalm(0,m), plm(i_mm+l-m)*dalm(1,m),kind=DP) 
!!$                  ENDIF
!!$
!!$                  b_north(m,ithl) = b_ns(1) + b_ns(-1)
!!$                  b_south(m,ithl) = b_ns(1) - b_ns(-1)
!!$
!!$                  !ENDIF ! test on nlmax
!!$
!!$               ENDDO ! loop on m
!!$
!!$               IF (DO_openmp()) THEN
!!$                  deALLOCATE (dalm)
!!$               ENDIF
!!$
!!$            ENDDO ! loop on ithl
!!$
!!$            IF (DO_openmp()) THEN
!!$               ALLOCATE(ring(0:nphmx-1),stat = status)
!!$               CALL assert_alloc(status,code,'ring')
!!$            ENDIF
!!$
!!$            DO ithl = 0, uchk - lchk
!!$
!!$               nphl = nph(ithl)
!!$               istart_north = startpix(ithl)
!!$               istart_south = npix-istart_north-nphl
!!$               ith  = ithl + lchk
!!$
!!$               CALL ring_synthesis(nsmax,nlmax,nmmax,b_north(0:nmmax,ithl),nphl,ring,kphi0(ithl))   ! north hemisph. + equator
!!$               map(l, istart_north:istart_north+nphl-1) = ring(0:nphl-1) 
!!$
!!$               IF (ith < nrings) THEN
!!$                  CALL ring_synthesis(nsmax,nlmax,nmmax,b_south(0:nmmax,ithl),nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
!!$                  map(l, istart_south:istart_south+nphl-1) = ring(0:nphl-1) 
!!$               ENDIF
!!$
!!$            ENDDO ! loop on rings (ithl)
!!$
!!$            IF (DO_openmp()) THEN
!!$               deALLOCATE(ring)
!!$            ENDIF
!!$
!!$         ENDDO ! loop on l
!!$
!!$         !$OMP END DO
!!$
!!$         !$OMP END parallel
!!$
!!$      ENDDO    ! loop on chunks
!!$
!!$      !     --------------------
!!$      !     free memory and exit
!!$      !     --------------------
!!$      IF (.not.DO_openmp()) THEN
!!$         deALLOCATE (ring,dalm)
!!$      ENDIF
!!$      deALLOCATE(b_north, b_south)
!!$      RETURN
!!$
!!$    END SUBROUTINE almn2alnspring_pre
!!$
!!$
!!$! ---------------------------------------------------------------------------------------!    
!!$    
!!$    SUBROUTINE alnspring2almn_iterative( nsmax, nlmax, nmmax, alnmap, &
!!$    	                & almn, zbounds, nb_iter )  
!!$
!!$    INTEGER(I4B), intent(IN)   :: nsmax, nlmax, nmmax, nb_iter
!!$
!!$    REAL(DP), intent(IN), dimension(0:nlmax,0:(12_i8b*nsmax)*nsmax-1) :: alnmap
!!$
!!$    INTEGER(I4B) :: n_plm, iter, l
!!$    REAL(DP), dimension(0:nlmax,0:(12_i8b*nsmax)*nsmax-1) :: alnmap_rec
!!$    
!!$    complex(DPC), dimension(1:1,0:nlmax,0:nmmax) :: almn
!!$
!!$    REAL(DP), intent(IN),  dimension(1:2), optional :: zbounds
!!$    REAL(DP), dimension(:,:), allocatable :: plm 
!!$
!!$    n_plm = nsmax*(nmmax+1)*(2*nlmax-nmmax+2) 
!!$    ALLOCATE(plm(0:n_plm-1,1:1))  
!!$    CALL plm_gen(nsmax, nlmax, nmmax, plm) 
!!$
!!$    alnmap_rec = alnmap
!!$    almn = 0.0
!!$    
!!$    DO iter=1,nb_iter
!!$       PRINT*
!!$       PRINT*
!!$       PRINT*, "Iteration : ", iter
!!$       PRINT*
!!$        ! Improve almn : old_almn + d_almn ! almn = CALL alnspring2almn( alnmap_rec )
!!$        !CALL alnspring2almn_pre( nsmax, nlmax, nmmax, alnmap, almn, zbounds, plm )
!!$        CALL alnspring2almn( nsmax, nlmax, nmmax, alnmap_rec, almn, zbounds )
!!$        ! New reconstruction ! alnmap_rec = CALL almn2alnspring( almn )
!!$        CALL almn2alnspring_pre( nsmax, nlmax, nmmax, alnmap_rec, almn, plm(:,1) )
!!$        DO l=0,nlmax
!!$           PRINT*, "Original : ", alnmap(l,0:6)
!!$           PRINT*, "Reconstr : ", alnmap_rec(l,0:6)
!!$           PRINT*
!!$        ENDDO
!!$        alnmap_rec = alnmap - alnmap_rec
!!$    ENDDO
!!$
!!$    END SUBROUTINE alnspring2almn_iterative
!!$
!!$

  ! ---------------------------------------------------------------------------------------!    

  SUBROUTINE alnspring2almn( nsmax, nlmax, nmmax, map, &
       & almn, zbounds )
    !=======================================================================
    !     computes the a(l,m) from a Temperature map for the HEALPIX pixelisation
    !        all from scratch
    !=======================================================================
    INTEGER(I4B), intent(IN)                    :: nsmax, nlmax, nmmax
    REAL(DP),   intent(IN),  dimension(0:nlmax,0:(12_i8b*nsmax)*nsmax-1) :: map
    complex(DPC), dimension(1:1,0:nlmax,0:nmmax) :: almn
    REAL(DP),     intent(IN),  dimension(1:2),         optional :: zbounds

    REAL(DP), dimension(1:2)         :: zbounds_in
    REAL(DP), dimension(1:2*nsmax,1) :: w8ring_in
    INTEGER(I4B) :: s, l, m, ith, scalem, scalel   ! alm related
    INTEGER(I8B) :: istart_south, istart_north, npix  ! map related
    INTEGER(I4B) :: nrings, nphmx
    REAL(DP)     :: omega_pix
    INTEGER(i4b), parameter :: SMAXCHK = 30
    INTEGER(kind=i4b), parameter :: RSMAX = 20, RSMIN = -20
    REAL(dp),          dimension(RSMIN:RSMAX) :: rescale_tab

    INTEGER(I4B)                              :: par_lm
    REAL(DP)              :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    REAL(DP)                      :: OVFLOW, UNFLOW
    REAL(DP),     dimension(-1:2)     :: phas_sd
    REAL(DP),     dimension(:,:), allocatable :: dalm
    REAL(DP),     dimension(:),   allocatable :: mfac
    REAL(DP),     dimension(:,:), allocatable :: recfac

    INTEGER(I4B)                              :: l_min
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:), allocatable :: phas_n, phas_s
    complex(DPC), dimension(:,:,:), allocatable :: phasl_n, phasl_s
    REAL(DP),     dimension(:),   allocatable :: ring
    INTEGER(I4B)                   :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    INTEGER(I8B), dimension(0:SMAXCHK-1) :: startpix
    INTEGER(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    REAL(DP),     dimension(0:SMAXCHK-1) :: cth, sth
    LOGICAL(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it

    INTEGER(I4B),     parameter :: LOG2LG   = 100
    REAL(KIND=DP),    parameter :: FL_LARGE = 2.0_dp **   LOG2LG
    REAL(KIND=DP),    parameter :: FL_SMALL = 2.0_dp ** (-LOG2LG)
    REAL(DP) :: logOVFLOW
    INTEGER(i4b) :: smax

    character(LEN=*), PARAMETER :: code = 'MAP2ALMSPRING'
    INTEGER(I4B) :: status
    !=======================================================================

    zbounds_in = (/-1.d0 , 1.d0/)
    IF (present(zbounds)) zbounds_in = zbounds
    w8ring_in  = 1.d0

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    omega_pix = FOURPI / REAL(npix, kind=DP)  ! pixel area (identical for all pixels)

    !     --- ALLOCATEs space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    ALLOCATE(mfac(0:nmmax),stat = status)
    CALL assert_alloc(status,code,'mfac')

    ALLOCATE(phasl_n(0:nlmax,0:nmmax,0:chunksize-1),stat = status)
    CALL assert_alloc(status,code,'phasl_n')

    ALLOCATE(phasl_s(0:nlmax,0:nmmax,0:chunksize-1),stat = status)
    CALL assert_alloc(status,code,'phasl_s')

    IF (.not.DO_openmp()) THEN
       ALLOCATE(ring(0:nphmx-1),stat = status)
       CALL assert_alloc(status,code,'ring')
       ALLOCATE(recfac(0:1,0:nlmax),stat = status)
       CALL assert_alloc(status,code,'recfac')
       ALLOCATE(dalm(1:2,0:nlmax),stat = status)
       CALL assert_alloc(status,code,'dalm')
       ALLOCATE(phas_n(0:nmmax,0:chunksize-1),stat = status)
       CALL assert_alloc(status,code,'phas_n')
       ALLOCATE(phas_s(0:nmmax,0:chunksize-1),stat = status)
       CALL assert_alloc(status,code,'phas_s')
    ENDIF

    !     ------------ initiate variables and arrays ----------------

    CALL gen_mfac(nmmax,mfac)

    logOVFLOW=log(FL_LARGE)
    smax = INT( log(MAX_DP) / logOVFLOW )

    IF (smax > (RSMAX-1)) THEN
       PRINT*,'Array rescale_tab too small in '//code
       PRINT*,smax ,'>', RSMAX
       stop
    ENDIF

    rescale_tab(RSMIN:RSMAX) = 0.0_dp
    DO s = -smax, smax
       rescale_tab(s) = FL_LARGE ** s
    ENDDO
    rescale_tab(0) = 1.0_dp

    OVFLOW = rescale_tab(1)
    UNFLOW = rescale_tab(-1)
    almn(1:1,0:nlmax,0:nmmax) = 0.0 ! set the whole alm array to zero

    ! loop on chunks
    DO ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       DO ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          CALL get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          ! find out which rings are to be analysed
          CALL select_rings(cth(ithl), zbounds_in, keep_north(ithl), keep_south(ithl), keep_it(ithl))
       ENDDO

       phasl_n = 0_dpc
       phasl_s = 0_dpc

       !$OMP parallel default(none) &
       !$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
       !$OMP      lchk, uchk, nph, startpix, kphi0, w8ring_in, &
       !$OMP      phasl_n, phasl_s, keep_north, keep_south, chunksize,map) &
       !$OMP  private(ithl, nphl, istart_north, istart_south, l, ith, ring, status,  phas_n, phas_s)

       IF (DO_openmp()) THEN
          ALLOCATE(ring(0:nphmx-1),stat = status)
          CALL assert_alloc(status,code,'ring')
          ALLOCATE(phas_n(0:nmmax,0:chunksize-1),stat = status)
          CALL assert_alloc(status,code,'phas_n')
          ALLOCATE(phas_s(0:nmmax,0:chunksize-1),stat = status)
          CALL assert_alloc(status,code,'phas_s')
       ENDIF

       phas_n = 0_dpc
       phas_s = 0_dpc

       !$OMP DO schedule(dynamic,1)
       DO l = 0, nlmax
          DO ith = lchk, uchk
             ithl = ith - lchk !local index
             nphl = nph(ithl)
             istart_north = startpix(ithl)
             istart_south = npix-istart_north-nphl
             ! DO Fourier Transform on rings
             IF (keep_north(ithl)) THEN
                !
                ring(0:nphl-1) = map(l,istart_north:istart_north+nphl-1)! * w8ring_in(ith,1)
                CALL ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_n(0,ithl), kphi0(ithl))
                phasl_n(l,0:nmmax,ithl) = phas_n(0:nmmax,ithl) + phasl_n(l,0:nmmax,ithl)
                !
             ENDIF

             IF (ith < nrings .and. keep_south(ithl)) THEN
                !
                ring(0:nphl-1) = map(l,istart_south:istart_south+nphl-1)! * w8ring_in(ith,1)
                CALL ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_s(0,ithl), kphi0(ithl))
                phasl_s(l,0:nmmax,ithl) = phas_s(0:nmmax,ithl) + phasl_s(l,0:nmmax,ithl)
                !
             ENDIF
          ENDDO ! loop on ring
       ENDDO
       !$OMP END DO


       IF (DO_openmp()) THEN
          deALLOCATE(ring)
          deALLOCATE(phas_s,phas_n)
       ENDIF
       !$OMP END parallel

       CALL init_rescale()
       OVFLOW = rescale_tab(1)
       UNFLOW = rescale_tab(-1)

       !$OMP parallel default(none) NUM_THREADS(8) &
       !$OMP shared(nlmax, nmmax, lchk, uchk, rescale_tab, ovflow, unflow, &
       !$OMP    cth, sth, mfac, almn, phasl_n, phasl_s, keep_it, omega_pix) &
       !$OMP private(recfac, dalm, phas_sd, status, m, ithl, l_min, &
       !$OMP   scalem, scalel, corfac, par_lm, lam_mm, lam_lm, lam_0, lam_1, lam_2, &
       !$OMP   cth_ring, l, php, phm)

       IF (DO_openmp()) THEN
          ALLOCATE(recfac(0:1,0:nlmax),stat = status)
          CALL assert_alloc(status,code,'recfac')
          ALLOCATE(dalm(1:2,0:nlmax),stat = status)
          CALL assert_alloc(status,code,'dalm')
       ENDIF

       !$OMP DO schedule(dynamic,1)
       DO m = 0, nmmax
          ! generate recursion factors (recfac) for Ylm of degree m
          CALL gen_recfac(nlmax, m, recfac)

          ! introduce DOuble precision vector to perform summation over ith for each l
          dalm(1:2, m:nlmax ) = 0.0_dp

          DO ithl = 0, uchk - lchk
             l_min = l_min_ylm(m, sth(ithl))
             IF (keep_it(ithl) .and. nlmax >= l_min) THEN ! avoid un-necessary calculations (EH, 09-2001)
                ! determine lam_mm
                CALL compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)

                !           ---------- l = m ----------
                par_lm = 1

               	php = phasl_n(m,m,ithl) + phasl_s(m,m,ithl) ! sum  (IF (l+m) even)
               	phm = phasl_n(m,m,ithl) - phasl_s(m,m,ithl) ! dIFf (IF (l+m) odd)
               	phas_sd(-1:0) =  (/ REAL(phm, kind=dp), aimag(phm) /)
               	phas_sd(1:2) =  (/ REAL(php, kind=dp), aimag(php) /)
                IF (m >= l_min) THEN
                   lam_lm = lam_mm * corfac !Actual lam_mm 
                   dalm(1:2, m) = dalm(1:2, m) + phas_sd(par_lm:par_lm+1) *lam_lm
                ENDIF


                !           ---------- l > m ----------
                lam_0 = 0.0_dp
                lam_1 = 1.0_dp
                scalel=0
                cth_ring = cth(ithl)
                lam_2 = cth_ring * lam_1 * recfac(0,m)

                DO l = m+1, nlmax

                   php = phasl_n(l,m,ithl) + phasl_s(l,m,ithl) ! sum  (IF (l+m) even)
                   phm = phasl_n(l,m,ithl) - phasl_s(l,m,ithl) ! dIFf (IF (l+m) odd)
                   phas_sd(-1:0) =  (/ REAL(phm, kind=dp), aimag(phm) /)
                   phas_sd(1:2) =  (/ REAL(php, kind=dp), aimag(php) /)

                   par_lm = - par_lm  ! = (-1)^(l+m)

                   IF (l >= l_min) THEN
                      lam_lm = lam_2 * corfac * lam_mm
                      dalm(1:2, l) = dalm(1:2, l) &
                           &       +  phas_sd(par_lm:par_lm+1) *lam_lm
                   ENDIF

                   lam_0 = lam_1 * recfac(1,l-1)
                   lam_1 = lam_2
                   lam_2 = (cth_ring * lam_1 - lam_0) * recfac(0,l)
                   IF (abs(lam_2) > OVFLOW) THEN
                      lam_1 = lam_1*UNFLOW
                      lam_2 = lam_2*UNFLOW
                      scalel= scalel + 1
                      corfac= rescale_tab(max(scalem+scalel,RSMIN))
                   ELSEIF (abs(lam_2) < UNFLOW) THEN
                      lam_1 = lam_1*OVFLOW
                      lam_2 = lam_2*OVFLOW
                      scalel= scalel - 1
                      corfac= rescale_tab(max(scalem+scalel,RSMIN))
                   ENDIF

                ENDDO ! loop on l
             ENDIF ! test on cut sky and nlmax
          ENDDO ! loop on ithl
          DO l = m, nlmax
             almn(1, l, m) = almn(1, l, m) + cmplx(dalm(1, l), dalm(2, l), kind=DP)
          ENDDO
       ENDDO ! loop on m
       !$OMP END DO
       IF (DO_openmp()) THEN
          deALLOCATE (recfac,dalm)
       ENDIF
       !$OMP END parallel
    ENDDO ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------

    deALLOCATE(mfac)
    IF (.not.DO_openmp()) THEN
       deALLOCATE (ring,recfac,dalm)
       deALLOCATE(phas_n,phas_s)
    ENDIF
    RETURN
    !
  END SUBROUTINE alnspring2almn

  ! ---------------------------------------------------------------------------------------!    

  SUBROUTINE survey2almn_opt( nsmax, nnmax, nlmax, nmmax, survey, nbpts, almn, kln, zbounds )
    !=======================================================================
    !     computes the a(n,l,m) decomposition from scratch
    !=======================================================================
    INTEGER(I4B), intent(IN)                        :: nsmax, nnmax, nlmax, nmmax, nbpts
    REAL(DP),   intent(IN), dimension(1:nbpts,1:3)  :: survey
    complex(DPC), dimension(1:nlmax,0:nlmax,0:nmmax)    :: almn
    REAL(DP),   intent(IN), dimension(0:nlmax,1:nnmax) :: kln
    REAL(DP), intent(IN), dimension(1:2),  optional :: zbounds
    REAL(DP),     dimension(:,:),   allocatable     :: map_north, map_south

    REAL(DP), dimension(1:2)         :: zbounds_in
    REAL(DP), dimension(1:3)         :: surveytemp 
    REAL(DP), dimension(1:2*nsmax,1) :: w8ring_in
    INTEGER(I4B) :: s, l, m, n, p, ipring, ith, scalem, scalel  
    INTEGER(I8B) :: istart_south, istart_north, npix  
    INTEGER(I8B) :: itotstart_north, itotstart_south, itotEND_north, itotEND_south
    INTEGER(I4B) :: nrings, nphmx, sizenorth, sizesouth
    REAL(DP)     :: omega_pix, jln
    INTEGER(i4b), parameter :: SMAXCHK = 30
    INTEGER(kind=i4b), parameter :: RSMAX = 20, RSMIN = -20
    REAL(dp),          dimension(RSMIN:RSMAX) :: rescale_tab

    INTEGER(I4B)                              :: par_lm
    REAL(DP)              :: lam_mm, lam_lm, lam_0, lam_1, lam_2, corfac, cth_ring
    REAL(DP)                      :: OVFLOW, UNFLOW
    REAL(DP),     dimension(-1:2)     :: phas_sd
    REAL(DP),     dimension(:,:), allocatable :: dalm
    REAL(DP),     dimension(:),   allocatable :: mfac
    REAL(DP),     dimension(:,:), allocatable :: recfac
    REAL(DP), dimension(:), allocatable :: tempaccjln

    INTEGER(I4B)                              :: l_min
    complex(DPC)                              :: php, phm
    complex(DPC), dimension(:,:), allocatable :: phas_n, phas_s
    complex(DPC), dimension(:,:,:), allocatable :: phasl_n, phasl_s
    REAL(DP),     dimension(:),   allocatable :: ring
    INTEGER(I4B)                   :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    INTEGER(I8B), dimension(0:SMAXCHK-1) :: startpix
    INTEGER(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    REAL(DP),     dimension(0:SMAXCHK-1) :: cth, sth
    LOGICAL(LGT), dimension(0:SMAXCHK-1) :: keep_north, keep_south, keep_it

    INTEGER(I4B),     parameter :: LOG2LG   = 100
    REAL(KIND=DP),    parameter :: FL_LARGE = 2.0_dp **   LOG2LG
    REAL(KIND=DP),    parameter :: FL_SMALL = 2.0_dp ** (-LOG2LG)
    REAL(DP) :: logOVFLOW, r
    INTEGER(i4b) :: smax
    INTEGER(i4b), dimension(1:nbpts) :: iprings

    character(LEN=*), PARAMETER :: code = 'survey2almn_opt'
    INTEGER(I4B) :: status

    !=======================================================================

    zbounds_in = (/-1.d0 , 1.d0/)
    IF (present(zbounds)) zbounds_in = zbounds
    w8ring_in  = 1.d0

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    omega_pix = FOURPI / REAL(npix, kind=DP)  ! pixel area (identical for all pixels)

    !     --- ALLOCATEs space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    ALLOCATE(mfac(0:nmmax),stat = status)
    CALL assert_alloc(status,code,'mfac')

    ALLOCATE(phasl_n(0:nlmax,0:nmmax,0:chunksize-1),stat = status)
    CALL assert_alloc(status,code,'phasl_n')

    ALLOCATE(phasl_s(0:nlmax,0:nmmax,0:chunksize-1),stat = status)
    CALL assert_alloc(status,code,'phasl_s')

    IF (.not.DO_openmp()) THEN
       ALLOCATE(ring(0:nphmx-1),stat = status)
       CALL assert_alloc(status,code,'ring')
       ALLOCATE(recfac(0:1,0:nlmax),stat = status)
       CALL assert_alloc(status,code,'recfac')
       ALLOCATE(dalm(1:2,0:nlmax),stat = status)
       CALL assert_alloc(status,code,'dalm')
       ALLOCATE(phas_n(0:nmmax,0:chunksize-1),stat = status)
       CALL assert_alloc(status,code,'phas_n')
       ALLOCATE(phas_s(0:nmmax,0:chunksize-1),stat = status)
       CALL assert_alloc(status,code,'phas_s')
    ENDIF

    !     ------------ initiate variables and arrays ----------------

    CALL gen_mfac(nmmax,mfac)

    logOVFLOW=log(FL_LARGE)
    smax = INT( log(MAX_DP) / logOVFLOW )

    IF (smax > (RSMAX-1)) THEN
       PRINT*,'Array rescale_tab too small in '//code
       PRINT*,smax ,'>', RSMAX
       stop
    ENDIF

    rescale_tab(RSMIN:RSMAX) = 0.0_dp
    DO s = -smax, smax
       rescale_tab(s) = FL_LARGE ** s
    ENDDO
    rescale_tab(0) = 1.0_dp

    OVFLOW = rescale_tab(1)
    UNFLOW = rescale_tab(-1)
    almn(1:nnmax,0:nlmax,0:nmmax) = 0.0 ! set the whole alm array to zero

    !$OMP parallel DO schedule(dynamic,1024) &
    !$OMP shared(survey, nsmax, nbpts, iprings) &
    !$OMP private(p)
    DO p = 1, nbpts
       CALL ang2pix_ring(nsmax, survey(p,2), survey(p,1), iprings(p)) 
    ENDDO
    !$OMP END parallel DO

    CALL message(code, start=.TRUE., msg="ModIFied FFT 3D algorithm")

    ! loop on chunks
    DO ichunk = 0, nchunks-1 ! for each chunk

       CALL message(code, msg="Block :",i=ichunk+1,msg2=" on ", i2=nchunks)

       lchk = ichunk * chunksize + 1  ! first ring of the chunk
       uchk = min(lchk+chunksize - 1, nrings)  ! last ring of the chunk

       ! locate pixels
       DO ith = lchk, uchk ! for each ring
          ithl = ith - lchk !local index
          ! get pixel location information
          CALL get_pixel_layout(nsmax, ith, cth(ithl), sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
          ! find out which rings are to be analysed
          CALL select_rings(cth(ithl), zbounds_in, keep_north(ithl), keep_south(ithl), keep_it(ithl))
       ENDDO

       itotstart_north = startpix(0)
       itotEND_north = startpix(uchk-lchk)+nph(uchk-lchk)-1
       itotstart_south = npix-startpix(uchk-lchk)-nph(uchk-lchk)
       itotEND_south = npix-startpix(0)-1
       sizenorth = (itotEND_north-itotstart_north)
       sizesouth = (itotEND_south-itotstart_south)

       !PRINT*,"North hemisphere:",itotstart_north,itotEND_north
       !PRINT*,"South hemisphere:",itotstart_south,itotEND_south

       ALLOCATE(map_north(0:nlmax,0:sizenorth),stat = status)
       CALL assert_alloc(status,code,'map_north')
       ALLOCATE(map_south(0:nlmax,0:sizesouth),stat = status)
       CALL assert_alloc(status,code,'map_south')

       ! loop on n parameter 
       DO n=1, nnmax

          !PRINT*,"> N = ",n 
          map_north=0.0_dp
          map_south=0.0_dp

          !$OMP parallel default(private) &
          !$OMP private(p, l, ipring, jln, tempaccjln, status, r ) &
          !$OMP shared(n, survey, nsmax, nbpts, nlmax, kln, map_north, map_south, &
          !$OMP   iprings, itotEND_north, itotstart_south, itotEND_south, itotstart_north )

          ALLOCATE(tempaccjln(0:nlmax),stat = status)
          CALL assert_alloc(status,code,'tempaccjln')

          !$OMP DO reduction(+:map_north, map_south) schedule(dynamic,1024) 
          DO p=1, nbpts
             r = survey(p,3)
             ipring = iprings(p)
             DO l=0, nlmax
                CALL BJL( l , kln(l,n)*r, jln ) 
                tempaccjln(l) =  kln(l,n)*jln
             ENDDO

             IF ( ipring .GE. itotstart_north .AND. ipring .LE. itotEND_north ) THEN ! IF north
                map_north(:,ipring-itotstart_north) =  map_north(:,ipring-itotstart_north) + tempaccjln(:)
             ENDIF
             IF ( ipring .GE. itotstart_south .AND. ipring .LE. itotEND_south ) THEN ! IF south
                map_south(:,ipring-itotstart_south) =  map_south(:,ipring-itotstart_south) + tempaccjln(:)
             ENDIF
          ENDDO
          !$OMP END DO

          deALLOCATE(tempaccjln)

          !$OMP END parallel

          !PRINT*,map_north(6,0:6)

          phasl_n = 0_dpc
          phasl_s = 0_dpc

          !$OMP parallel default(none) &
          !$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, n, &
          !$OMP      itotstart_north, itotEND_north, itotstart_south, itotEND_south, &
          !$OMP      lchk, uchk, nph, startpix, kphi0, w8ring_in, &
          !$OMP      phasl_n, phasl_s, keep_north, keep_south, chunksize, map_north, map_south ) &
          !$OMP  private(ithl, nphl, istart_north, istart_south, l, ith, ring, status,  phas_n, phas_s)

          IF (DO_openmp()) THEN
             ALLOCATE(ring(0:nphmx-1),stat = status)
             CALL assert_alloc(status,code,'ring')
             ALLOCATE(phas_n(0:nmmax,0:chunksize-1),stat = status)
             CALL assert_alloc(status,code,'phas_n')
             ALLOCATE(phas_s(0:nmmax,0:chunksize-1),stat = status)
             CALL assert_alloc(status,code,'phas_s')
          ENDIF

          phas_n = 0_dpc
          phas_s = 0_dpc

          !$OMP DO schedule(dynamic,1)
          DO l = 0, nlmax

             DO ith = lchk, uchk
                ithl = ith - lchk !local index
                nphl = nph(ithl)
                istart_north = startpix(ithl)
                istart_south = npix-istart_north-nphl
                ! DO Fourier Transform on rings
                IF (keep_north(ithl)) THEN
                   !
                   ring(0:nphl-1) = map_north(l,(istart_north-itotstart_north):(istart_north+nphl-1-itotstart_north))! * w8ring_in(ith,1)
                   CALL ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_n(0,ithl), kphi0(ithl))
                   phasl_n(l,0:nmmax,ithl) = phas_n(0:nmmax,ithl) + phasl_n(l,0:nmmax,ithl)
                   !
                ENDIF
                IF (ith < nrings .and. keep_south(ithl) ) THEN
                   !
                   ring(0:nphl-1) = map_south(l,(istart_south-itotstart_south):&
                        &(istart_south+nphl-1-itotstart_south))! * w8ring_in(ith,1)
                   CALL ring_analysis(nsmax, nlmax, nmmax, ring, nphl, phas_s(0,ithl), kphi0(ithl))
                   phasl_s(l,0:nmmax,ithl) = phas_s(0:nmmax,ithl) + phasl_s(l,0:nmmax,ithl)
                   !
                ENDIF
             ENDDO ! loop on ring
          ENDDO ! loop on l
          !$OMP END DO

          IF (DO_openmp()) THEN
             deALLOCATE(ring)
             deALLOCATE(phas_s,phas_n)
          ENDIF
          !$OMP END parallel

          CALL init_rescale()
          OVFLOW = rescale_tab(1)
          UNFLOW = rescale_tab(-1)

          !$OMP parallel default(none) &
          !$OMP shared(nlmax, nmmax, lchk, uchk, rescale_tab, n, ovflow, unflow, &
          !$OMP    cth, sth, mfac, almn, phasl_n, phasl_s, keep_it, omega_pix) &
          !$OMP private(recfac, dalm, phas_sd, status, m, ithl, l_min, &
          !$OMP   scalem, scalel, corfac, par_lm, lam_mm, lam_lm, lam_0, lam_1, lam_2, &
          !$OMP   cth_ring, l, php, phm)

          IF (DO_openmp()) THEN
             ALLOCATE(recfac(0:1,0:nlmax),stat = status)
             CALL assert_alloc(status,code,'recfac')
             ALLOCATE(dalm(1:2,0:nlmax),stat = status)
             CALL assert_alloc(status,code,'dalm')
          ENDIF

          !$OMP DO reduction(+:almn) schedule(dynamic,1)
          DO m = 0, nmmax
             ! generate recursion factors (recfac) for Ylm of degree m
             CALL gen_recfac(nlmax, m, recfac)

             ! introduce DOuble precision vector to perform summation over ith for each l
             dalm(1:2, m:nlmax ) = 0.0_dp

             DO ithl = 0, uchk - lchk
                l_min = 0 !l_min_ylm(m, sth(ithl))
                IF (keep_it(ithl) .and. nlmax >= l_min) THEN ! avoid un-necessary calculations (EH, 09-2001)
                   ! determine lam_mm
                   CALL compute_lam_mm(mfac(m), m, sth(ithl), lam_mm, corfac, scalem)

                   !           ---------- l = m ----------
                   par_lm = 1

                   php = phasl_n(m,m,ithl) + phasl_s(m,m,ithl) ! sum  (IF (l+m) even)
                   phm = phasl_n(m,m,ithl) - phasl_s(m,m,ithl) ! dIFf (IF (l+m) odd)
                   phas_sd(-1:0) =  (/ REAL(phm, kind=dp), aimag(phm) /)
                   phas_sd(1:2) =  (/ REAL(php, kind=dp), aimag(php) /)
                   IF (m >= l_min) THEN
                      lam_lm = lam_mm * corfac !Actual lam_mm 
                      dalm(1:2, m) = dalm(1:2, m) + phas_sd(par_lm:par_lm+1) *lam_lm
                   ENDIF


                   !           ---------- l > m ----------
                   lam_0 = 0.0_dp
                   lam_1 = 1.0_dp
                   scalel=0
                   cth_ring = cth(ithl)
                   lam_2 = cth_ring * lam_1 * recfac(0,m)

                   DO l = m+1, nlmax

                      php = phasl_n(l,m,ithl) + phasl_s(l,m,ithl) ! sum  (IF (l+m) even)
                      phm = phasl_n(l,m,ithl) - phasl_s(l,m,ithl) ! dIFf (IF (l+m) odd)
                      phas_sd(-1:0) =  (/ REAL(phm, kind=dp), aimag(phm) /)
                      phas_sd(1:2) =  (/ REAL(php, kind=dp), aimag(php) /)

                      par_lm = - par_lm  ! = (-1)^(l+m)

                      IF (l >= l_min) THEN
                         lam_lm = lam_2 * corfac * lam_mm
                         dalm(1:2, l) = dalm(1:2, l) &
                              &       +  phas_sd(par_lm:par_lm+1) *lam_lm
                      ENDIF

                      lam_0 = lam_1 * recfac(1,l-1)
                      lam_1 = lam_2
                      lam_2 = (cth_ring * lam_1 - lam_0) * recfac(0,l)
                      IF (abs(lam_2) > OVFLOW) THEN
                         lam_1 = lam_1*UNFLOW
                         lam_2 = lam_2*UNFLOW
                         scalel= scalel + 1
                         corfac= rescale_tab(max(scalem+scalel,RSMIN))
                      ELSEIF (abs(lam_2) < UNFLOW) THEN
                         lam_1 = lam_1*OVFLOW
                         lam_2 = lam_2*OVFLOW
                         scalel= scalel - 1
                         corfac= rescale_tab(max(scalem+scalel,RSMIN))
                      ENDIF

                   ENDDO ! loop on l
                ENDIF ! test on cut sky and nlmax
             ENDDO ! loop on ithl
             DO l = m, nlmax
                almn(n, l, m) = almn(n, l, m) + cmplx(dalm(1, l), dalm(2, l), kind=DP)
             ENDDO
          ENDDO ! loop on m
          !$OMP END DO
          IF (DO_openmp()) THEN
             deALLOCATE (recfac,dalm)
          ENDIF
          !$OMP END parallel

       ENDDO ! loop on n

       deALLOCATE(map_north, map_south)

    ENDDO ! loop on chunks

    CALL message(code,fin=.TRUE.)

    !     --------------------
    !     free memory and exit
    !     --------------------

    deALLOCATE(mfac)
    IF (.not.DO_openmp()) THEN
       deALLOCATE (ring,recfac,dalm)
       deALLOCATE(phas_n,phas_s)
    ENDIF
    RETURN
    !
  END SUBROUTINE survey2almn_opt

  ! ---------------------------------------------------------------------------------------!    

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
    INTEGER(I4B)  :: status, nbpts, ipring, nrr, nb_iter
    INTEGER(I4B)  :: nsmax, nnmax, nlmax, nmmax, n, p, k, m, l
    REAL(DP)	  :: rmax, tempval, order
    !	
    REAL(DP),    dimension(:,:), allocatable	     :: map
    REAL(DP), dimension(1:nbpts,1:3)         	     :: survey
    REAL(DP), dimension(1:3)                         :: surveytemp
    REAL(DP),    dimension(1:2)       	             :: zbounds
    REAL(DP),    dimension(0:nlmax,1:nnmax)         :: kln
    !REAL(DP),    dimension(:,:,:) , allocatable     :: jln
    REAL(DP)                                        :: jln
    REAL(DP),    dimension(:,:)   , allocatable     ::  plm
    complex(DP), dimension(1:nnmax,0:nlmax,0:nmmax) :: almn
    character(len=*), PARAMETER                     :: code = "survey2almn_srs"
    !

    almn = cmplx( 0.0, 0.0, kind=DP )

120 FORMAT(A,I3)

    ALLOCATE(map(0:nlmax,0:(12*nsmax**2-1)),stat = status)
    CALL assert_alloc(status,code,"map")

    PRINT*,"================================================"
    PRINT*,"   Computing multidimensional maps" 
    PRINT*,"   > ModIFied FFT 3D algorithm"
    PRINT*,"------------------------------------------------"

    DO n = 1, nnmax

       PRINT*,"> Order n = ",n!,OMP_GET_THREAD_NUM()
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

       !PRINT*,map(6,0:6)

       CALL alnspring2almn( nsmax, nlmax, nmmax, map, &
            & almn(n:n, 0:nlmax,0:nmmax), zbounds )

    ENDDO

    PRINT*,"------------------------------------------------"
    PRINT*,"   Finished"
    PRINT*,"================================================"


    DEALLOCATE( map )    

    RETURN
    !
  END SUBROUTINE survey2almn_srs

  ! ---------------------------------------------------------------------------------------!    

  SUBROUTINE fieldalmn2overalmn( almn, kln, rhomean, rmax, nnmax, nlmax, nmmax )

    IMPLICIT NONE

    INTEGER(I4B) :: nnmax, nlmax, nmmax, n
    REAL(DP) :: rmax, rhomean
    COMPLEX(DP),DIMENSION(1:nnmax,0:nlmax,0:nmmax) :: almn
    REAL(DP),DIMENSION(0:nlmax,1:nnmax) :: kln
    REAL(DP),DIMENSION(1:nnmax,1:1,1:1) :: jln

    PRINT*, "rhomean", rhomean

    almn(:,:,:) = almn(:,:,:) / rhomean

    DO n = 1, nnmax
       CALL BJL( 1 , kln(1,n)*rmax , jln(n,1,1) )
       PRINT*,sqrt(FOURPI)*(rmax**2.0)*jln(n,1,1)
       almn(n,0,0) = almn(n,0,0) - sqrt(FOURPI)*(rmax**2.0)*jln(n,1,1)
    ENDDO

    RETURN

  END SUBROUTINE fieldalmn2overalmn

  ! ---------------------------------------------------------------------------------------!    

  SUBROUTINE almn2rmap(map, almn, rho, nsmax, nnmax, nlmax, nmmax, kln, cln, plm)

    IMPLICIT NONE

    REAL(DP)             :: rho
    INTEGER(I4B), intent(IN)                     :: nsmax, nnmax, nlmax, nmmax
    complex(DPC), intent(IN),  dimension(1:nnmax,0:nlmax,0:nmmax)      :: almn
    REAL(DP),     intent(OUT), dimension(0:(12_i8b*nsmax)*nsmax-1,1:1) :: map
    REAL(DP),     intent(IN),  dimension(0:( nsmax*(nmmax+1)*(2*nlmax-nmmax+2) ) )   :: plm

    INTEGER(I4B) :: l, m, ith, n
    INTEGER(I8B) :: istart_south, istart_north, npix
    INTEGER(I4B) :: nrings, nphmx

    INTEGER(I8B) :: n_lm, n_plm, i_mm
    REAL(DP),     dimension(0:nlmax,1:nnmax)  :: cln, kln,jlns
    complex(DPC), dimension(-1:1)             :: b_ns
    REAL(DP),     dimension(:,:), allocatable :: dalm
    INTEGER(i4b)                              :: ll, l_min, l_start
    REAL(DP)                                  :: cth
    complex(DPC), dimension(:,:), allocatable :: b_north, b_south
    REAL(DP),     dimension(:),   allocatable :: ring
    INTEGER(i4b)             :: nchunks, chunksize, ichunk, lchk, uchk, ithl, nphl
    INTEGER(i4b), parameter :: SMAXCHK = 50
    INTEGER(I8B), dimension(0:SMAXCHK-1) :: startpix
    INTEGER(I4B), dimension(0:SMAXCHK-1) :: nph, kphi0
    REAL(DP),     dimension(0:SMAXCHK-1) :: sth

    character(LEN=*), parameter :: code = 'ALM2MAP'
    INTEGER(I4B) :: status
    !=======================================================================

    ! Healpix definitions
    nrings = 2*nsmax           ! number of isolatitude rings on N. hemisphere + equat
    npix   = (12_I8B*nsmax)*nsmax    ! total number of pixels on the sphere
    nphmx  = 4*nsmax           ! maximum number of pixels/ring
    n_lm   = ((nmmax+1_I8B)*(2*nlmax-nmmax+2))/2 !number of (l,m) with m in[0,M] and l in [m,L]
    n_plm  = n_lm * nrings

    !     --- ALLOCATEs space for arrays ---
    nchunks   = nrings/SMAXCHK + 1  ! number of chunks
    chunksize = (nrings+nchunks-1)/nchunks ! <= SMAXCHK

    ALLOCATE(b_north(0:nmmax, 0:chunksize-1),stat = status)
    CALL assert_alloc(status,code,'b_north')

    ALLOCATE(b_south(0:nmmax, 0:chunksize-1),stat = status)
    CALL assert_alloc(status,code,'b_south')

    IF (.not. DO_openmp()) THEN
       ALLOCATE(dalm(0:1,0:nlmax), stat = status)
       CALL assert_alloc(status,code,'dalm')
       ALLOCATE(ring(0:nphmx-1),stat = status)
       CALL assert_alloc(status,code,'ring')
    ENDIF


    !     ------------ Compute required bessel FUNCTIONs ----------------

    CALL gen_jln(jlns, kln, rho, nnmax, nlmax)


    !     ------------ initiate variables and arrays ----------------

    map = 0.0 ! set the whole map to zero

    ! loop on chunks
    DO ichunk = 0, nchunks-1
       lchk = ichunk * chunksize + 1
       uchk = min(lchk+chunksize - 1, nrings)

       DO ith = lchk, uchk
          ithl = ith - lchk !local index
          ! get pixel location information
          CALL get_pixel_layout(nsmax, ith, cth, sth(ithl), nph(ithl), startpix(ithl), kphi0(ithl))
       ENDDO

       b_north(:,:) = 0_dpc ! pad with zeros
       b_south(:,:) = 0_dpc

       !$OMP parallel default(none) NUM_THREADS(8) &
       !$OMP shared(nlmax, nmmax, sth, lchk, uchk, plm, b_north, b_south, n_lm, jlns, cln, almn, nnmax) &
       !$OMP private(dalm, b_ns, status, m, ith, ithl, l_min, l_start, l, ll, i_mm)

       IF (DO_openmp()) THEN
          ALLOCATE(dalm(0:1,0:nlmax), stat = status)
          CALL assert_alloc(status,code,'dalm')
       ENDIF


       !$OMP DO schedule(dynamic,1)
       DO m = 0, nmmax
    	  dalm = 0.0_dp
       ! extract needed alm under memory and CPU efficient form
          DO ll = m, nlmax
             DO n = 0, nnmax
                dalm(0,ll) = dalm(0,ll) + REAL(almn(n,ll,m),kind=dp) * jlns(ll,n) * cln(ll,n)
                dalm(1,ll) = dalm(1,ll) + aimag(almn(n,ll,m)) * jlns(ll,n) * cln(ll,n)
             ENDDO
          ENDDO
          DO ithl = 0, uchk - lchk
             l_min = l_min_ylm(m, sth(ithl))
             IF (nlmax >= l_min) THEN ! skip calculations when Ylm too small
                ith = ithl + lchk
                i_mm = n_lm * (ith-1) + ((2_I8B*nlmax + 3 - m)*m)/2 ! location of Ym,m for ring ith

                !           ---------- l = m ----------

                IF (m >= l_min) THEN
                   b_ns( 1) = cmplx(plm(i_mm) * dalm(0,m), plm(i_mm) * dalm(1,m), kind=DP)
                   b_ns(-1) = 0.0_dpc
                ELSE
                   b_ns = 0.0_dpc
                ENDIF

                !           ---------- l > m ----------
                l_start = max(m+1, l_min) ! odd values of (l+m)
                IF (mod(m+l_start,2) == 0) l_start = l_start+1
                DO l = l_start, nlmax, 2
                   b_ns(-1) = b_ns(-1) + &
                        & cmplx(plm(i_mm+l-m) * dalm(0,l), plm(i_mm+l-m)*dalm(1,l), kind=DP)
                ENDDO ! loop on l

                l_start = max(m+2, l_min) ! even values of (l+m)
                IF (mod(m+l_start,2) == 1) l_start = l_start+1
                DO l = l_start, nlmax, 2
                   b_ns(1)  = b_ns(1) + &
                        & cmplx(plm(i_mm+l-m) * dalm(0,l), plm(i_mm+l-m)*dalm(1,l),kind=DP) 
                ENDDO ! loop on l

                b_north(m,ithl) = b_ns(1) + b_ns(-1)
                b_south(m,ithl) = b_ns(1) - b_ns(-1)
             ENDIF ! test on nlmax
          ENDDO ! loop on rings (ithl)
       ENDDO ! loop on m
       !$OMP END DO

       IF (DO_openmp()) THEN
          deALLOCATE (dalm)
       ENDIF
       !$OMP END parallel

       !$OMP parallel default(none) NUM_THREADS(8) &
       !$OMP  shared(nsmax, nlmax, nmmax, npix, nrings, nphmx, &
       !$OMP      lchk, uchk, b_north, b_south, nph, startpix, kphi0, map) &
       !$OMP  private(ithl, nphl, istart_north, istart_south, &
       !$OMP      ith, ring, status)
       IF (DO_openmp()) THEN
          ALLOCATE(ring(0:nphmx-1),stat = status)
          CALL assert_alloc(status,code,'ring')
       ENDIF
       !$OMP DO schedule(dynamic,1)
       DO ithl = 0, uchk - lchk
          nphl = nph(ithl)
          istart_north = startpix(ithl)
          istart_south = npix-istart_north-nphl
          ith  = ithl + lchk

          CALL ring_synthesis(nsmax,nlmax,nmmax,b_north(0,ithl),nphl,ring,kphi0(ithl))   ! north hemisph. + equator
          map(istart_north:istart_north+nphl-1,1) = map(istart_north:istart_north+nphl-1,1) + ring(0:nphl-1)

          IF (ith < nrings) THEN
             CALL ring_synthesis(nsmax,nlmax,nmmax,b_south(0,ithl),nphl,ring,kphi0(ithl)) ! south hemisph. w/o equat
             map(istart_south:istart_south+nphl-1,1) = map(istart_south:istart_south+nphl-1,1) + ring(0:nphl-1)
          ENDIF
       ENDDO ! loop on ithl
       !$OMP END DO
       IF (DO_openmp()) THEN
          deALLOCATE(ring)
       ENDIF
       !$OMP END parallel
    ENDDO    ! loop on chunks

    !     --------------------
    !     free memory and exit
    !     --------------------
    IF (.not.DO_openmp()) THEN
       deALLOCATE (ring,dalm)
    ENDIF
    deALLOCATE(b_north, b_south)
    RETURN

  END SUBROUTINE almn2rmap

  ! ---------------------------------------------------------------------------------------!    

  !> Computes series of roots of Bessel functions
  !!@param[out] qln(0:nlmax,1:nnmax) : roots
  !!@param[in] (nnmax,nlmax) : bounds
  SUBROUTINE gen_qln(qln, nnmax, nlmax)

    IMPLICIT NONE

    INTEGER(I4B) :: nnmax, nlmax, l, nrr, maxrt, lastRoot, err
    PARAMETER(maxrt = 1000)
    REAL(DP) :: rmax, order,A,B, A_before
    REAl(DP), DIMENSION(0:nlmax,1:nnmax) :: qln
    !	REAL*8, DIMENSION(1:maxrt) :: roots
    REAL*8, DIMENSION(:), allocatable :: roots
    character(LEN=*), parameter :: code = 'gen_qln'
    INTEGER(I4B) :: status

    CALL message(code,start=.TRUE.)

    !$OMP parallel &
    !$OMP shared(qln,rmax,nnmax,nlmax,status) &
    !$OMP private(roots,order,A,B,A_before,l,nrr,lastRoot,err)

    ALLOCATE(roots(1:maxrt),stat = status)
    CALL assert_alloc(status,code,'roots')

    !$OMP DO schedule(dynamic,1)
    DO l = 0, nlmax

       order = REAL(l) + 0.5_dp
       lastRoot = 1

       PRINT*, "Order l =",l

       A=0.0
       B=2000.0

       DO WHILE ( lastRoot .LT. nnmax)

          roots = 0.0_sp
          A_before = A
          CALL ROOTBESSJ( order, A, B, maxrt, nrr, err, roots )
          !		PRINT*,"number of roots found",nrr,err
          IF(nrr .EQ. 0) THEN
             B = B + 1000.0
          ELSE IF(err .NE. 1) THEN
             PRINT*,"err code :",err
             A = A_before
             B = (B - A)/2.0 + A
          ELSE

             IF(lastRoot + nrr .LE. nnmax) THEN
                qln(l,lastRoot:lastRoot + nrr - 1) = roots(1:nrr)
                lastRoot = lastRoot + nrr - 1
             ELSE
                qln(l,lastRoot:nnmax) = roots(1: nnmax - lastRoot + 1)
                lastRoot = nnmax
             ENDIF
             B = A + 2048.0
          ENDIF

       ENDDO

    ENDDO
    !$OMP END DO
    deALLOCATE(roots)
    !$OMP END parallel

    CALL message(code,fin=.TRUE.)

  END SUBROUTINE gen_qln

  ! ---------------------------------------------------------------------------------------!    

  !> Computes k spectrum
  !!@param[out] kln(0:nlmax,1:nnmax) : k spectrum
  !!@param[in] (nnmax,nlmax) : bounds
  !!@param[in] rmax : R max value
  SUBROUTINE gen_kln(kln, nnmax, nlmax, rmax)

    IMPLICIT NONE

    INTEGER(I4B) :: nnmax, nlmax, l, nrr
    REAL(DP) :: rmax, order
    REAl(DP), DIMENSION(0:nlmax,1:nnmax) :: kln
    !	REAL*8, DIMENSION(1:1000) :: roots

    !	DO l = 0, nlmax
    !
    !	    order = REAL(l) + 0.5_dp
    !
    !	    roots = 0.0_sp
    !
    !	    CALL ROOTBESSJ( order, nrr, roots )
    !
    !	    !PRINT*,order,nrr,roots
    !	    !IF (nrr <= nnmax) THEN
    !	    !   PRINT*, "ERROR during roots computation"
    !	    !   STOP
    !	    !ELSE
    !
    !	    kln(l,1:nnmax) = roots(1:(nnmax)) / rmax
    !
    !	    !ENDIF
    !
    !	ENDDO

    CALL gen_qln(kln,nnmax,nlmax)

    kln = kln / rmax


  END SUBROUTINE gen_kln

  ! ---------------------------------------------------------------------------------------!    

  !> Generates a logrange
  SUBROUTINE logrange( x , mn , mx , npts )

    IMPLICIT NONE

    REAL(DP), dimension(1:npts) :: x
    REAL(DP) :: mn, mx, rg
    INTEGER :: npts, i

    rg = (log10(mx)-log10(mn))/REAL(npts)

    DO i = 1, npts
       x(i) = (10)**(REAL(i)*rg + log10(mn))
    ENDDO

    RETURN
  END SUBROUTINE logrange

  ! ---------------------------------------------------------------------------------------!    

  !> Computes series of normalization coefficients
  !!@param[out] cln(0:nlmax,1:nnmax) : output normalization coefficients
  !!@param[in] kln(0:nlmax,1:nnmax) : k spectrum
  !!@param[in] (nnmax,nlmax) : bounds
  !!@param[in] rmax : R max value
  SUBROUTINE gen_cln(cln, kln, nnmax, nlmax, rmax)

    IMPLICIT NONE

    INTEGER(I4B) :: nnmax, nlmax, l, n, k
    REAL(DP) :: rmax, tempval
    REAl(DP), DIMENSION(0:nlmax,1:nnmax) :: kln
    REAl(DP), DIMENSION(0:nlmax,1:nnmax) :: cln
    CHARACTER(LEN=*), PARAMETER :: code = 'gen_cln'

    CALL message(code,start=.TRUE.)

    !$OMP PARALLEL &
    !$OMP SHARED(cln,rmax,kln,nnmax,nlmax) &
    !$OMP PRIVATE(l,n,k,tempval)	  

    !$OMP DO SCHEDULE(DYNAMIC,1)	 
    DO l = 0, nlmax 
       DO n = 1, nnmax   
          k = n+1
          CALL BJL( k , kln(l,n)*rmax, tempval )
          cln(l,n) = ( sqrt(2.0)*rmax**(-1.5) ) / (kln(l,n)*tempval)
       ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL   	

    CALL message(code,fin=.TRUE.)

  END SUBROUTINE gen_cln

  ! ---------------------------------------------------------------------------------------!    

  !> Computes series of jl(kln r)
  !!@param[out] jlns(0:nlmax,1:nnmax) : jlns 
  !!@param[in] kln(0:nlmax,1:nnmax) : k spectrum
  !!@param[in] rho : double radius value
  !!@param[in] (nnmax,nlmax) : bounds
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
    DO n=1,nnmax   
       DO l=0,nlmax
          !PRINT*,n,l,OMP_GET_THREAD_NUM()
          CALL BJL( l , kln(l,n)*rho , jln )
          jlns(l,n) = kln(l,n) * jln
       ENDDO
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL

    RETURN

  END SUBROUTINE gen_jln

  ! ---------------------------------------------------------------------------------------!    

  !> Compute the value of the l-th order spherical bessel FUNCTION at x
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
    REAL*8 X,JL
    REAL*8 AX,AX2
    REAL,PARAMETER::LN2=0.6931471805599453094D0
    REAL,PARAMETER::ONEMLN2=0.30685281944005469058277D0
    REAL,PARAMETER::PID2=1.5707963267948966192313217D0
    REAL,PARAMETER::PID4=0.78539816339744830961566084582D0
    REAL,parameter::ROOTPI12 = 21.269446210866192327578D0
    REAL,parameter::GAMMA1 =   2.6789385347077476336556D0 !/* Gamma FUNCTION of 1/3 */
    REAL,parameter::GAMMA2 =   1.3541179394264004169452D0 !/* Gamma FUNCTION of 2/3 */
    REAL,PARAMETER::PI=3.141592653589793238463D0
    REAL*8 NU,NU2,BETA,BETA2,COSB
    REAL*8 sx,sx2
    REAL*8 cotb,cot3b,cot6b,secb,sec2b
    REAL*8 trigarg,expterm,L3

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
       ELSEIF((REAL(L)**2/AX).LT.5.D-1)THEN
          BETA=AX-PID2*(L+1)
          JL=(DCOS(BETA)*(1.D0-(NU2-0.25D0)*(NU2-2.25D0)/8.D0/AX2*(1.D0-(NU2-6.25)*(NU2-12.25D0)/48.D0/AX2)) &
               -DSIN(BETA)*(NU2-0.25D0)/2.D0/AX* (1.D0-(NU2-2.25D0)*(NU2-6.25D0)/24.D0/AX2*(1.D0-(NU2-12.25)* &
               (NU2-20.25)/80.D0/AX2)) )/AX   
       ELSE
          L3=NU**0.325
          IF(AX .LT. NU-1.31*L3) THEN
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

          ELSEIF (AX .GT. NU+1.48*L3) THEN
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

  ! ---------------------------------------------------------------------------------------!    

END MODULE f3dex_transforms
