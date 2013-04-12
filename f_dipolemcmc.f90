
!> @file
!! \brief The driver program to solve a linear system with default options.
!!
!! <pre>
!! -- Distributed SuperLU routine (version 3.2) --
!! Lawrence Berkeley National Lab, Univ. of California Berkeley.
!! October, 2012
!! </pre>
!
      program f_dipolemcmc
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
      use healpix_types
      use fitstools
      use pix_tools
      use alm_tools
      use head_fits
      use quiet_mapfile_mod
      use quiet_hdf_mod
      use MCMCmodule
      implicit none
      include 'mpif.h'

      integer(i4b) 				   :: i, j, k, n
      integer(i4b)                                 :: nside, npix, L_max
      !parameter (nside = 16, npix = 12*nside**2, L_max = nside*2)
      parameter (nside = 256, npix = 12*nside**2, L_max = 64)
      integer 				           :: l, m, absm, nummap
      real(dp)       			           :: nullval, start, finish
      complex(dpc)                                 :: s_obs_lu(1:L_max*L_max+2*L_max+1), S_iso_lu(1:L_max*L_max+2*L_max+1)
      real(dp), allocatable, dimension(:,:)        :: map_obs, cl
      real(dp)                                     :: cl_temp(1:5)
      complex(dpc), allocatable, dimension(:,:,:)  :: s_obs_LM
      character(len=80), dimension(1:60) 	   :: header, comment, card
      logical                                      :: anynull
      real(dp)                                     :: bestfitparams(nn), f, ff, x_test(nn)
      type (mcmcParams) mctrl
      integer                                      :: status, clock, seed, ierr
      character (len=19)                           :: fstr1, fstr2
      CHARACTER*100				   :: filename
      
      common  s_obs_lu, S_iso_lu

      external   ff

! Initialize MPI environment 
      call mpi_init(ierr)

  nummap = 1

do k=1, nummap

  allocate(map_obs(0:npix-1,1))

  !write (filename, 10) k
!10    format ('modmaps_2/test_ns16lmax32_notsynfast_mod_', I0, '.fits')
  !call read_bintab(filename,map_obs,npix,1,nullval,anynull)
  !call read_bintab('commander_wmaponly_cl_cmb_test1_mix.fits',map_obs,npix,1,nullval,anynull)
  !call read_bintab('maps_ns256_lmax200/test_ns256lmax200_notsynfast_1.fits',map_obs,npix,1,nullval,anynull)
  call read_bintab('sim.fits',map_obs,npix,1,nullval,anynull)

  allocate(s_obs_LM(1, 0:L_max, 0:L_max))
  call map2alm_iterative(nside, L_max, L_max, 2, map_obs, s_obs_LM)
  deallocate(map_obs)

  i = 1
  absm = 0
     do l = absm, L_max
        s_obs_lu(i) = s_obs_LM(1,l,absm)
        i = i + 1
     enddo

  do absm = 1, L_max
     do l = absm, L_max
        s_obs_lu(i) = s_obs_LM(1,l,absm)
        i = i + 1
        s_obs_lu(i) = conjg(s_obs_LM(1,l,absm))*(-1)**(-absm)
        i = i + 1
     enddo
  enddo

  deallocate(s_obs_LM)

  allocate(cl(0:L_max,1))
  !call fits2cl('cl.fits', cl, L_max, 1, header)

  open(29, file = "camb_61139232_scalcls.dat")
  do i =2, L_max
    read(29,*)l,cl_temp
    cl(i,1)=cl_temp(1)*2*pi/i/(i+1)
  enddo
  cl(0,1) = 0.d0
  cl(1,1) = 0.d0
  close(29)

  i = 1
  absm = 0
     do l = absm, L_max
        S_iso_lu(i) = cl(l,1)
        i = i + 1
     enddo

  do absm = 1, L_max
     do l = absm, L_max
        S_iso_lu(i) = cl(l,1)
        i = i + 1
        S_iso_lu(i) = cl(l,1)
        i = i + 1
     enddo
  enddo

  deallocate(cl)

  call random_seed

! Initialise MCMC
  call init_MCMC

! Now call MCMC
  call MCMC(ff,nn,mctrl,bestfitparams,f,status)

  write(*,*)'The best-fit likelihood is ',f

enddo

! Terminate the MPI execution environment
      call mpi_finalize(ierr)

end program f_dipolemcmc

function ff(n,x_in)

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

!     Input:
      integer n
      real*8 x_in(n)

      integer(superlu_ptr) :: grid
      integer(superlu_ptr) :: options
      integer(superlu_ptr) :: ScalePermstruct
      integer(superlu_ptr) :: LUstruct
      integer(superlu_ptr) :: SOLVEstruct
      integer(superlu_ptr) :: A
      integer(superlu_ptr) :: stat

      integer maxn, maxnz, maxnrhs
      !parameter ( maxn = 5, maxnz = 100, maxnrhs = 10 )
      integer, allocatable, dimension(:)                              :: rowind, colptr
      double complex, allocatable, dimension(:)                       :: values, b, c, berr, diagonals
      double complex                                                  :: logdet_M2, logdet_Siso, logdet
      integer n_slu, m_slu, nnz, nprow, npcol, ldb, init
      integer*4 iam, info, ierr, ldb4, nrhs
      integer  numberrow, numbercol, nnz_local, firstrow
      integer(i4b) 				                      :: i, j, k
      integer 				                              :: l, m, lprime, mprime, absm, absmprime, firstelm
      integer(i4b)                                                    :: nside, npix, L_max
      !parameter (nside = 16, npix = 12*nside**2, L_max = nside*2)
      parameter (nside = 256, npix = 12*nside**2, L_max = 64)
      complex(dpc), allocatable, dimension(:,:,:)                     :: gamma_LM
      complex(dpc)                                                    :: gamma_00, gamma_1M(-1:1)
      complex(dpc)                                                    :: s_obs_lu(1:L_max*L_max+2*L_max+1), S_iso_lu(1:L_max*L_max+2*L_max+1)
      real(dp)       			                              :: alpha, p_theta, p_phi, dist
      real(dp), allocatable, dimension(:)                             :: p_vector, n_vector
      real(dp), allocatable, dimension(:,:)                           :: gamma_pix
      complex(dpc)                                                    :: temp, chi2, loglikelihood
      real(dp)                                                        :: f
      real(dp)                                                        :: start, finish, start_tot, finish_tot

      common  s_obs_lu, S_iso_lu

!     Output:
      real*8 ff

 call cpu_time(start_tot)

! Initialize MPI environment 
      !call mpi_init(ierr)

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

      alpha = x_in(1)
      !alpha = x_in(1) * 0.5
      !alpha = 0.4
      p_theta = x_in(2) * DEG2RAD
      !p_theta = x_in(2) * 180 * DEG2RAD
      !p_theta = 107 * DEG2RAD
      p_phi = x_in(3) * DEG2RAD
      !p_phi = x_in(3) * 360 * DEG2RAD
      !p_phi = 226 * DEG2RAD

 !call cpu_time(start)

      gamma_00 = sqrt(4*pi)
      gamma_1M(-1) = sqrt(2*pi/3)*alpha*CMPLX(sin(p_theta)*cos(p_phi),sin(p_theta)*sin(p_phi))
      gamma_1M(0) = sqrt(4*pi/3)*alpha*cos(p_theta)
      gamma_1M(1) = -sqrt(2*pi/3)*alpha*CMPLX(sin(p_theta)*cos(p_phi),-sin(p_theta)*sin(p_phi))

 !call cpu_time(finish)
 !write(*,*)'finish-start (1) =',finish-start

      maxn = L_max*L_max+2*L_max+1
      !maxn = 5
      maxnz = 200000
      maxnrhs = 10

      allocate(values(maxnz))
      allocate(b(maxn))
      allocate(berr(maxnrhs))
      allocate(rowind(maxnz))
      allocate(colptr(maxn))

      nnz = 0
      temp = 0.d0
      i = 1
      j = 1

 !call cpu_time(start)

      do absmprime = 0, L_max
         do lprime = absmprime, L_max
           mprime = absmprime
           do while (1 .eq. 1)
            firstelm = 0
            k = 0
            do absm = 0, L_max
              do l = absm, L_max
               m = absm
               do while (1 .eq. 1)
               k = k+1
                   if (abs(m-mprime) .gt. 1) then
                     temp = 0.d0
                   else
                     temp = Mod(real(l,dp),real(m,dp),real(lprime,dp),real(mprime,dp),gamma_1M,L_max)
                   endif
                   if (temp .ne. 0.d0) then
                     nnz = nnz + 1
                     values(i) = temp
                     !rowind(i) = l*(l+1)+m+1
                     rowind(i) = k
                     if (firstelm == 0) then
                         colptr(j) = i
                         j = j + 1
                     endif
                     i = i + 1
                     firstelm = 1
                   endif
                   if ((m .eq. 0) .or. (m .eq. -absm)) exit
                   m = m-2*absm
                 enddo
                enddo
             enddo
           if ((mprime .eq. 0) .or. (mprime .eq. -absmprime)) exit
           mprime = mprime-2*absmprime
           enddo
         enddo
      enddo
      colptr(j) = nnz + 1

 !call cpu_time(finish)
 !write(*,*)'finish-start (2) =',finish-start

! Read Harwell-Boeing matrix, and adjust the pointers and indices
! to 0-based indexing, as required by C routines.

 !call cpu_time(start)

      m_slu = maxn
      n_slu = maxn

      if ( iam == 0 ) then
         do i = 1, n_slu+1
            colptr(i) = colptr(i) - 1
         enddo
         do i = 1, nnz
            rowind(i) = rowind(i) - 1
         enddo
      endif

 !call cpu_time(finish)
 !write(*,*)'finish-start (3) =',finish-start

 !call cpu_time(start)

! Distribute the matrix to the process gird
      call  f_zcreate_dist_matrix(A, m_slu, n_slu, nnz, values, rowind, colptr, grid)

! Setup the right hand side
      call get_CompRowLoc_Matrix(A, nrow=numberrow, ncol=numbercol, nnz_loc=nnz_local, nrow_loc=ldb, fst_row=firstrow)
      do i = 1, ldb
         b(i)=s_obs_lu(firstrow+i)
      enddo
      nrhs = 1
      ldb4 = ldb
      write(*,*)'ldb = ',ldb

 !call cpu_time(finish)
 !write(*,*)'finish-start (4) =',finish-start

! Set the default input options
      call f_set_default_options(options)

! Change one or more options
      !call set_superlu_options(options,Fact=DOFACT)
      !call set_superlu_options(options,ParSymbFact=YES)
      !call set_superlu_options(options,IterRefine=NO)
      !call set_superlu_options(options,Equil=NO)
      !call set_superlu_options(options,SymPattern=YES)
      call set_superlu_options(options,PrintStat=NO)
      !call set_superlu_options(options,num_lookaheads=5)

! Modify one or more options
      call set_superlu_options(options,ColPerm=MMD_AT_PLUS_A)
      call set_superlu_options(options,RowPerm=NOROWPERM)

! Initialize ScalePermstruct and LUstruct
      call get_SuperMatrix(A, nrow=m_slu, ncol=n_slu)
      call f_ScalePermstructInit(m_slu, n_slu, ScalePermstruct)
      call f_LUstructInit(m_slu, n_slu, LUstruct)

! Initialize the statistics variables
      call f_PStatInit(stat)

!      allocate(c(maxn))
!      c = b

 !call cpu_time(start)

! Call the linear equation solver
      call f_pzgssvx(options, A, ScalePermstruct, b, ldb4, nrhs, &
                     grid, LUstruct, SOLVEstruct, berr, stat, info)

 !call cpu_time(finish)
 !write(*,*)'finish-start (5) =',finish-start

      !call set_superlu_options(options,Fact=SamePattern)

 !call cpu_time(start)

 !     call f_pzgssvx(options, A, ScalePermstruct, c, ldb4, nrhs, &
 !                    grid, LUstruct, SOLVEstruct, berr, stat, info)

 !call cpu_time(finish)
 !write(*,*)'finish-start (5.5) =',finish-start

!write(*,*)'b(10) = ',b(10)
!write(*,*)'c(10) = ',c(10)

 !call cpu_time(start)

      if (info == 0) then
         write (*,*) 'Backward error: ', (berr(i), i = 1, nrhs)
      else
         write(*,*) 'INFO from f_pdgssvx = ', info
      endif

      logdet_Siso = 0.d0
      chi2 = 0.d0
      do i = 1, ldb
         if (S_iso_lu(i) .ne. 0.d0) then
            chi2 = chi2 + b(i)*conjg(b(i))*S_iso_lu(i)**(-1)
            logdet_Siso = logdet_Siso+log(S_iso_lu(i))
         endif
      enddo
      write(*,*)'chi2 = ',chi2
      write(*,*)'logdet_Siso = ',logdet_Siso

      allocate(diagonals(ldb))
      diagonals = 0.d0
      call f_zGetDiagU(ldb, LUstruct, grid, diagonals)
      logdet_M2 = 0.d0
      do i = 1, ldb
         logdet_M2 = logdet_M2+log(diagonals(i)*conjg(diagonals(i)))
      enddo
      write(*,*)'logdet_M2 = ',logdet_M2
      deallocate(diagonals)

      logdet = logdet_Siso+logdet_M2
      write(*,*)'logdet = ',logdet
      
      loglikelihood = -0.5*chi2-0.5*logdet
      write(*,*)'loglikelihood = ',loglikelihood

 !call cpu_time(finish)
 !write(*,*)'finish-start (6) =',finish-start

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
      !call mpi_finalize(ierr)

      ff = -loglikelihood
      !ff = chi2

 call cpu_time(finish_tot)
 write(*,*)'finish-start (ff) =',finish_tot-start_tot

contains

function Mod(l,m,lprime,mprime,gamma_1M,L_max) !Gaunt_integral

   IMPLICIT NONE
   complex(dpc)                                 :: Mod
   integer(i4b)                                 :: L_max, NDIM, IER
   complex(dpc)                                 :: gamma_1M(-1:1)
   real(dp), allocatable, dimension(:)          :: p_vector, n_vector, THRCOF
   real(dp)                                     :: l, m, lprime, mprime, M2MIN, M2MAX

   if ((m == mprime) .and. (l == lprime)) then
      Mod = 1.d0
   else if ((abs(l-lprime) == 1) .and. (abs(m-mprime) .le. abs(l-lprime))) then
      NDIM = 2*L_max
      allocate(THRCOF(1:NDIM))
      THRCOF = 0.d0
      CALL DRC3JM (l, lprime, 1.d0, 0.d0, M2MIN, M2MAX, THRCOF, NDIM, IER)
      Mod = (-1)**m*gamma_1M(m-mprime)*sqrt((2*l+1)*(2*lprime+1)*3/4/pi)*THRCOF(1-M2MIN)
      CALL DRC3JM (l, lprime, 1.d0, -m, M2MIN, M2MAX, THRCOF, NDIM, IER)
      Mod = Mod*THRCOF(1-M2MIN+mprime)
      deallocate(THRCOF)
   else
      Mod = 0.d0
   endif

end function Mod



end function ff
