
!> @file
!! \brief The driver program to solve a linear system with default options.
!!
!! <pre>
!! -- Distributed SuperLU routine (version 3.2) --
!! Lawrence Berkeley National Lab, Univ. of California Berkeley.
!! October, 2012
!! </pre>
!
      program f_dipolelike
! 
! Purpose
! =======
!
! The driver program F_PDDRIVE.
!
! This example illustrates how to use F_PDGSSVX with the full
! (default) options to solve a linear system.
! 
! Seven basic steps are required:
!   1. Create C structures used in SuperLU_DIST
!   2. Initialize the MPI environment and the SuperLU process grid
!   3. Set up the input matrix and the right-hand side
!   4. Set the options argument
!   5. Call f_pdgssvx
!   6. Release the process grid and terminate the MPI environment
!   7. Release all structures
!
!
      use superlu_mod
      use healpix_types
      use fitstools
      use pix_tools
      use alm_tools
      use head_fits
      use quiet_mapfile_mod
      use quiet_hdf_mod
      implicit none
      include 'mpif.h'
      integer maxn, maxnz, maxnrhs
      !parameter ( maxn = 5, maxnz = 100, maxnrhs = 10 )
      integer, allocatable, dimension(:)      :: rowind, colptr
      double complex, allocatable, dimension(:)       :: values, b, berr, diagonals
      double complex                                  :: logdet_M2, logdet_Siso, logdet
      integer n_slu, m_slu, nnz, nprow, npcol, ldb, init
      integer*4 iam, info, ierr, ldb4, nrhs
      integer  numberrow, numbercol, nnz_local, firstrow

      integer(superlu_ptr) :: grid
      integer(superlu_ptr) :: options
      integer(superlu_ptr) :: ScalePermstruct
      integer(superlu_ptr) :: LUstruct
      integer(superlu_ptr) :: SOLVEstruct
      integer(superlu_ptr) :: A
      integer(superlu_ptr) :: stat

  integer(i4b) 				       :: i, j, n, nside, npix, L_max, NDIM, IER
  integer 				       :: l, m, lprime, mprime, firstelm
  real(dp)       			       :: nullval, alpha, p_theta, p_phi, dist
  real(dp)                                     :: M2MIN, M2MAX
  real(dp), allocatable, dimension(:)          :: p_vector, n_vector, THRCOF
  complex(dpc), allocatable, dimension(:)      :: s_obs_lu, S_iso_lu
  real(dp), allocatable, dimension(:,:)        :: map_obs, gamma_pix, dw8, der1, cl
  complex(dpc), allocatable, dimension(:,:)    :: Mod_lu
  complex(dpc)                                 :: temp, chi2, loglikelihood
  complex(dpc), allocatable, dimension(:,:,:)  :: gamma_LM, s_obs_LM
  real(dp), allocatable, dimension(:,:,:,:,:,:):: Y
  real(dp), dimension(2)                       :: z
  character(len=80), dimension(1:60) 	       :: header, comment, card
  logical                                      :: anynull

! Initialize MPI environment 
      call mpi_init(ierr)

! Check malloc
!      call f_check_malloc(iam)

! Create Fortran handles for the C structures used in SuperLU_DIST
      call f_create_gridinfo_handle(grid)
      call f_create_options_handle(options)
      call f_create_ScalePerm_handle(ScalePermstruct)
      call f_create_LUstruct_handle(LUstruct)
      call f_create_SOLVEstruct_handle(SOLVEstruct)
      call f_create_SuperMatrix_handle(A)
      call f_create_SuperLUStat_handle(stat)

! Initialize the SuperLU_DIST process grid
      nprow = 1
      npcol = 1
      call f_superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, grid)

! Bail out if I do not belong in the grid. 
      call get_GridInfo(grid, iam=iam)
      if ( iam >= nprow * npcol ) then 
         go to 100
      endif
      if ( iam == 0 ) then 
         write(*,*) ' Process grid ', nprow, ' X ', npcol
      endif

!if ( iam == 0 ) then 
  nside = 32
  npix  = 12*nside**2

  allocate(map_obs(0:npix-1,1))
  call read_bintab('test_modns32alpha075.fits',map_obs,npix,1,nullval,anynull)

  read(*,*) alpha
  !alpha = 2.d0
  p_theta = 107 * DEG2RAD
  p_phi = 226 * DEG2RAD

  allocate(p_vector(3))
  allocate(n_vector(3))
  allocate(gamma_pix(0:npix-1,1))
  call ang2vec(p_theta, p_phi, p_vector)
  do n = 0, npix-1
     call pix2vec_ring(nside, n, n_vector)
     call angdist(p_vector, n_vector, dist)
     gamma_pix(n,1) = 1+alpha*cos(dist)
  enddo
  deallocate(p_vector)
  deallocate(n_vector)

  !L_max = 100
  L_max =nside*2
  allocate(gamma_LM(1, 0:2*L_max, 0:2*L_max))
  call map2alm_iterative(nside, 2*L_max, 2*L_max, 2, gamma_pix, gamma_LM)

  allocate(cl(0:L_max,1))
  call fits2cl('cl.fits', cl, L_max, 1, header)

  allocate(s_obs_LM(1, 0:L_max, 0:L_max))
  call map2alm_iterative(nside, L_max, L_max, 2, map_obs, s_obs_LM)

  allocate(S_iso_lu(1:L_max*L_max+2*L_max+1))
  do l =0, L_max
     do m = -l, l
        S_iso_lu(l*(l+1)+m+1) = cl(l,1)
     enddo
  enddo

      maxn = L_max*L_max+2*L_max+1
      write(*,*)'maxn = ',maxn
      !maxn = 5
      maxnz = 200000
      maxnrhs = 10
  allocate(values(maxnz))
  allocate(b(maxn))
  allocate(berr(maxnrhs))
  allocate(rowind(maxnz))
  allocate(colptr(maxn))

  allocate(s_obs_lu(1:L_max*L_max+2*L_max+1))
  do l =0, L_max
     do m = 0, l
        s_obs_lu(l*(l+1)+m+1) = s_obs_LM(1,l,m)
     enddo
     do m = -l, -1
        s_obs_lu(l*(l+1)+m+1) = conjg(s_obs_LM(1,l,-m))*(-1)**m
     enddo
  enddo
  !write(*,*)s_obs_lu

  !maxnz = 0
  nnz = 0
  temp = 0.d0
  i = 1
  j = 1
  !allocate(Mod_lu(1:L_max*L_max+2*L_max+1, 1:L_max*L_max+2*L_max+1))

  do lprime = 0, L_max
     do mprime = -lprime, lprime
        firstelm = 0
        do l =0, L_max
           do m = -l, l
             !if (m-mprime .ge. 0) Mod_lu(l*(l+1)+m+1,lprime*(lprime+1)+mprime+1) = Mod(real(l,dp),real(m,dp),real(lprime,dp),real(mprime,dp),gamma_LM(1,1,m-mprime),L_max)
             if (m-mprime .ge. 0) temp = Mod(real(l,dp),real(m,dp),real(lprime,dp),real(mprime,dp),gamma_LM(1,1,m-mprime),L_max)
             !if (m-mprime .lt. 0) Mod_lu(l*(l+1)+m+1,lprime*(lprime+1)+mprime+1) = Mod(real(l,dp),real(m,dp),real(lprime,dp),real(mprime,dp),(-1)**(m-mprime)*conjg(gamma_LM(1,1,mprime-m)),L_max)
             if (m-mprime .lt. 0) temp = Mod(real(l,dp),real(m,dp),real(lprime,dp),real(mprime,dp),(-1)**(m-mprime)*conjg(gamma_LM(1,1,mprime-m)),L_max)
             !if (Mod_lu(l*(l+1)+m+1,lprime*(lprime+1)+mprime+1) .ne. 0.d0) maxnz = maxnz + 1
             if (temp .ne. 0.d0) then
                !write(*,*)'temp = ',temp
                nnz = nnz + 1
                values(i) = temp
                rowind(i) = l*(l+1)+m+1
                if (firstelm == 0) then
                   colptr(j) = i
                   j = j + 1
                endif
                i = i + 1
                firstelm = 1
             endif
             !if (temp == 0.d0) write(*,*)'temp = ',temp
          enddo
        enddo
     enddo
  enddo
  colptr(j) = nnz + 1
  !write(*,*)'j = ',j
  !write(*,*)Mod_lu(1,1), Mod_lu(99,99)
  !write(*,*)maxnz
  !write(*,*)nnz

  deallocate(cl)
  deallocate(map_obs)
  deallocate(gamma_pix)
  deallocate(gamma_LM)
  deallocate(s_obs_LM)

! Read Harwell-Boeing matrix, and adjust the pointers and indices
! to 0-based indexing, as required by C routines.

      m_slu = maxn
      n_slu = maxn
      !nnz = maxnz

      if ( iam == 0 ) then 
      !   open(file = "../EXAMPLE/cg20.cua", status = "old", unit = 5)
      !   call zhbcode1(m, n, nnz, values, rowind, colptr)
      !   close(unit = 5)

         do i = 1, n_slu+1
            colptr(i) = colptr(i) - 1
         enddo
         do i = 1, nnz
            rowind(i) = rowind(i) - 1
         enddo

      endif
!endif
      
! Distribute the matrix to the process gird
      call  f_zcreate_dist_matrix(A, m_slu, n_slu, nnz, values, rowind, colptr, grid)

! Setup the right hand side
      call get_CompRowLoc_Matrix(A, nrow=numberrow, ncol=numbercol, nnz_loc=nnz_local, nrow_loc=ldb, fst_row=firstrow)
      !do i = 1, ldb
      !   b(i) = 1.0
      !enddo
      do i = 1, ldb
         b(i)=s_obs_lu(firstrow+i)
      enddo
      nrhs = 1
      ldb4 = ldb
      write(*,*)'ldb = ',ldb
      !write(*,*)'b(100) = ',b(100)

      !chi2 = 0.d0
      !do i = 5, ldb
      !   chi2 = chi2 + b(i)*conjg(b(i))*S_iso_lu(i)**(-1)
         !write(*,*)'S_iso_lu = ',S_iso_lu(i)
      !enddo
      !write(*,*)'chi2 = ',chi2

! Set the default input options
      call f_set_default_options(options)

! Change one or more options
!      call set_superlu_options(options,Fact=FACTORED)
!      call set_superlu_options(options,ParSymbFact=YES)

! Modify one or more options
      call set_superlu_options(options,ColPerm=NATURAL)
      call set_superlu_options(options,RowPerm=NOROWPERM)

! Initialize ScalePermstruct and LUstruct
      call get_SuperMatrix(A, nrow=m_slu, ncol=n_slu)
      call f_ScalePermstructInit(m_slu, n_slu, ScalePermstruct)
      call f_LUstructInit(m_slu, n_slu, LUstruct)

! Initialize the statistics variables
      call f_PStatInit(stat)

! Call the linear equation solver
      call f_pzgssvx(options, A, ScalePermstruct, b, ldb4, nrhs, &
                     grid, LUstruct, SOLVEstruct, berr, stat, info)

      if (info == 0) then
         write (*,*) 'Backward error: ', (berr(i), i = 1, nrhs)
      else
         write(*,*) 'INFO from f_pdgssvx = ', info
      endif
      !write(*,*)b
      !write(*,*)'b(100) = ',b(100)
      
      logdet_Siso = 0.d0
      chi2 = 0.d0
      do i = 5, ldb
         chi2 = chi2 + b(i)*conjg(b(i))*S_iso_lu(i)**(-1)
         logdet_Siso = logdet_Siso+log(S_iso_lu(i))
         !write(*,*)'logS_iso_lu = ',S_iso_lu(i)
      enddo
      write(*,*)'chi2 = ',chi2
      write(*,*)'logdet_Siso = ',logdet_Siso

      allocate(diagonals(ldb))
      diagonals = 0.d0
      call f_zGetDiagU(ldb, LUstruct, grid, diagonals)
      !write(*,*)'diagonals = ',diagonals
      logdet_M2 = 0.d0
      do i = 1, ldb
         logdet_M2 = logdet_M2+log(diagonals(i)*conjg(diagonals(i)))
      enddo
      write(*,*)'logdet_M2 = ',logdet_M2
      !write(*,*)'ldb = ',ldb
      deallocate(diagonals)

      logdet = logdet_Siso+logdet_M2
      write(*,*)'logdet = ',logdet
      
      loglikelihood = -0.5*chi2-0.5*logdet
      write(*,*)'loglikelihood = ',loglikelihood

  !deallocate(Mod_lu)
  deallocate(s_obs_lu)
  deallocate(S_iso_lu)
  deallocate(values)
  deallocate(b)
  deallocate(berr)
  deallocate(rowind)
  deallocate(colptr)

! Deallocate the storage allocated by SuperLU_DIST
      call f_PStatFree(stat)
      call f_Destroy_CompRowLoc_Mat_dist(A)
      call f_ScalePermstructFree(ScalePermstruct)
      call f_Destroy_LU(n_slu, grid, LUstruct)
      call f_LUstructFree(LUstruct)
      call get_superlu_options(options, SolveInitialized=init)
      if (init == YES) then
         call f_zSolveFinalize(options, SOLVEstruct)
      endif

! Release the SuperLU process grid
100   call f_superlu_gridexit(grid)

! Deallocate the C structures pointed to by the Fortran handles
      call f_destroy_gridinfo_handle(grid)
      call f_destroy_options_handle(options)
      call f_destroy_ScalePerm_handle(ScalePermstruct)
      call f_destroy_LUstruct_handle(LUstruct)
      call f_destroy_SOLVEstruct_handle(SOLVEstruct)
      call f_destroy_SuperMatrix_handle(A)
      call f_destroy_SuperLUStat_handle(stat)

! Check malloc
!      call f_check_malloc(iam)

! Terminate the MPI execution environment
      call mpi_finalize(ierr)

      !stop

contains

FUNCTION Mod(l,m,lprime,mprime,gamma_LM,L_max)

IMPLICIT NONE
   complex(dpc)                                 :: Mod
   integer(i4b)                                 :: L_max, NDIM, IER
   complex(dpc)                                 :: gamma_LM
   real(dp), allocatable, dimension(:)          :: p_vector, n_vector, THRCOF
   real(dp)                                     :: l, m, lprime, mprime, M2MIN, M2MAX

   if ((m == mprime) .and. (l == lprime)) then
      Mod = 1.d0
   else if ((abs(l-lprime) == 1) .and. (abs(m-mprime) .le. abs(l-lprime))) then

      NDIM = 2*L_max
      allocate(THRCOF(1:NDIM))
      THRCOF = 0.d0
      CALL DRC3JM (l, lprime, 1.d0, 0.d0, M2MIN, M2MAX, THRCOF, NDIM, IER)
      Mod = (-1)**m*gamma_LM*sqrt((2*l+1)*(2*lprime+1)*3/4/pi)*THRCOF(1-M2MIN)
      CALL DRC3JM (l, lprime, 1.d0, -m, M2MIN, M2MAX, THRCOF, NDIM, IER)
      Mod = Mod*THRCOF(1-M2MIN+mprime)
      !write(*,*)'here'
      deallocate(THRCOF)
   else
      Mod = 0.d0
   endif

END FUNCTION Mod

      end program f_dipolelike
