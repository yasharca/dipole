
!> @file
!! \brief The driver program to solve a linear system with default options.
!!
!! <pre>
!! -- Distributed SuperLU routine (version 3.2) --
!! Lawrence Berkeley National Lab, Univ. of California Berkeley.
!! October, 2012
!! </pre>
!
      program f_test
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
      parameter ( maxn = 5, maxnz = 100, maxnrhs = 10 )
      integer rowind(maxnz), colptr(maxn)
      real*8  values(maxnz), b(maxn), berr(maxnrhs), B_rhs(maxn), X_sol(maxn), diagonals(maxn), det
      integer n, m, nnz, nprow, npcol, ldb, init
      integer*4 iam, info, i, ierr, ldb4, nrhs
      integer  numberrow, numbercol, nnz_local, firstrow
 
      integer(superlu_ptr) :: grid
      integer(superlu_ptr) :: options
      integer(superlu_ptr) :: ScalePermstruct
      integer(superlu_ptr) :: LUstruct
      integer(superlu_ptr) :: SOLVEstruct
      integer(superlu_ptr) :: A
      integer(superlu_ptr) :: stat

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

! Read Harwell-Boeing matrix, and adjust the pointers and indices
! to 0-based indexing, as required by C routines.

      m = 5
      n = 5
      nnz = 12
      B_rhs(1) = 1.0
      B_rhs(2) = 2.0
      B_rhs(3) = 3.0
      B_rhs(4) = 4.0
      B_rhs(5) = 5.0

      if ( iam == 0 ) then 
      !   open(file = "../EXAMPLE/g20.rua", status = "old", unit = 5)
      !   call dhbcode1(m, n, nnz, values, rowind, colptr)
      !   close(unit = 5)
!

      values(1) = 19.0
      values(2) = 12.0
      values(3) = 12.0
      values(4) = 21.0
      values(5) = 12.0
      values(6) = 12.0
      values(7) = 21.0
      values(8) = 16.0
      values(9) = 21.0
      values(10) = 5.0
      values(11) = 21.0
      values(12) = 18.0

      rowind(1) = 1
      rowind(2) = 2
      rowind(3) = 5
      rowind(4) = 2
      rowind(5) = 3
      rowind(6) = 5
      rowind(7) = 1
      rowind(8) = 3
      rowind(9) = 1
      rowind(10) = 4
      rowind(11) = 4
      rowind(12) = 5

      colptr(1) = 1
      colptr(2) = 4
      colptr(3) = 7
      colptr(4) = 9
      colptr(5) = 11
      colptr(6) = 13

         do i = 1, n+1
            colptr(i) = colptr(i) - 1
         enddo
         do i = 1, nnz
            rowind(i) = rowind(i) - 1
         enddo

      !write(*,*)'I am here!'
      endif

      
! Distribute the matrix to the process gird
      call  f_dcreate_dist_matrix(A, m, n, nnz, values, rowind, colptr, grid)

! Setup the right hand side
      call get_CompRowLoc_Matrix(A, nrow=numberrow, ncol=numbercol, nnz_loc=nnz_local, nrow_loc=ldb, fst_row=firstrow)
      !do i = 1, ldb
      !   b(i) = 1.0
      !enddo
      nrhs = 1
      ldb4 = ldb
      write(*,*)'ldb = ',ldb
      write(*,*)'nrow = ',numberrow
      write(*,*)'ncol = ',numbercol
      write(*,*)'nnz_loc = ',nnz_local
      write(*,*)'fst_row = ',firstrow
      do i = 1, ldb
         b(i)=B_rhs(firstrow+i)
      enddo
      write(*,*)'b(2) = ',b(2)

! Set the default input options
      call f_set_default_options(options)

! Change one or more options
!      call set_superlu_options(options,Fact=FACTORED)
!      call set_superlu_options(options,ParSymbFact=YES)

! Modify one or more options
      call set_superlu_options(options,ColPerm=NATURAL)
      call set_superlu_options(options,RowPerm=NOROWPERM)

! Initialize ScalePermstruct and LUstruct
      call get_SuperMatrix(A, nrow=m, ncol=n)
      call f_ScalePermstructInit(m, n, ScalePermstruct)
      call f_LUstructInit(m, n, LUstruct)

! Initialize the statistics variables
      call f_PStatInit(stat)

! Call the linear equation solver
      call f_pdgssvx(options, A, ScalePermstruct, b, ldb4, nrhs, &
                     grid, LUstruct, SOLVEstruct, berr, stat, info)

      if (info == 0) then
         write (*,*) 'Backward error: ', (berr(i), i = 1, nrhs)
      else
         write(*,*) 'INFO from f_pdgssvx = ', info
      endif
      !write(*,*)b
      do i = 1, ldb
         X_sol(firstrow+i)=b(i)
      enddo
      write(*,*)'b(2) = ',b(2)
      write(*,*)'n = ',n
      write(*,*)'LUstruct = ',LUstruct
      write(*,*)'diagonals = ',diagonals

      call f_GetDiagU(n, LUstruct, grid, diagonals)
      write(*,*)'diagonals = ',diagonals
      det = 1.d0
      do i = 1, n
         det = det*diagonals(i)
      enddo
      write(*,*)'det = ',det

! Deallocate the storage allocated by SuperLU_DIST
      call f_PStatFree(stat)
      call f_Destroy_CompRowLoc_Mat_dist(A)
      call f_ScalePermstructFree(ScalePermstruct)
      call f_Destroy_LU(n, grid, LUstruct)
      call f_LUstructFree(LUstruct)
      call get_superlu_options(options, SolveInitialized=init)
      if (init == YES) then
         call f_dSolveFinalize(options, SOLVEstruct)
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
      !write(*,*)X_sol

      stop
      end
