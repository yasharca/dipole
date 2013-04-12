
module MCMCmodule
      use healpix_types
      use fitstools
      use pix_tools
      use alm_tools
      use head_fits
      use quiet_mapfile_mod
      use quiet_hdf_mod

      implicit none

!     Number of parameters - remember to change the hardset parmin and parmax at ~line 150 if you change this
      integer, parameter :: nn = 3
!     Format specificiers for output files
      character (len=19), private :: fstr1, fstr2
!     Global debug flag
      logical, private, parameter :: debug = .false.

      type mcmcParams
        logical :: outputMpcty = .false.
        integer :: fileunit = 89
        double precision :: initialWidths(nn) = 0.05, crudeTolerance = 1d-2
        !character (len=200) :: chainame = 'mcmchains/mcmchain_20_wmap09_test1_lmax64_mix_euler_posalpha_aroudpeak_widesteps_3.out'
        !character (len=200) :: chainame = 'mcmchains/mcmchain_12_ns256lmax200_notsynfast_map1_orig_euler_posalpha_aroudpeak_widesteps_3.out'
        character (len=200) :: chainame = 'mcmchains/sim.out'
        !character (len=200) :: chainame = 'mcmchains/mcmchain_wmap09_test1_4_lmax80.out'
      end type

!     User-supplied uniform random number generator
      !double precision :: urand
      !external urand
!
!     Function urand should not take any arguments.  If the user wishes
!     to be able to initialize urand, so that the same sequence of
!     random numbers can be repeated, this capability could be imple-
!     mented with a separate subroutine, and called from the user's
!     driver program.  An example urand function (and initialization
!     subroutine) which uses the function ran0 (the "minimal standard"
!     random number generator of Park and Miller [Comm. ACM 31, 1192-
!     1201, Oct 1988; Comm. ACM 36 No. 7, 105-110, July 1993]) is
!     provided.

contains

subroutine init_MCMC

     write(fstr1,'("(I5",i4,"(X, F20.10))")')nn+1
     write(fstr2,'("(",i4,"(X, F20.10))")')nn+1

end subroutine init_MCMC


subroutine MCMC(ff,n,mctrl,x,f,status)
!=======================================================================
!     Optimization (maximization) of user-supplied "fitness" function
!     ff  over n-dimensional parameter space  x  using a basic 
!     Metropolis-Hastings Markov Chain Monte Carlo method.  ff must be
!     positive definite.
!
!     Uses a seperable, adaptive, n-dimensional Gaussian proposal function.
!
!
!     Input:
      integer, intent(IN) :: n
      type(mcmcParams), intent(IN) :: mctrl
      real*8     :: ff
      external  ff

!	INTERFACE
!    		!the likelihood function
!    		real*8 function ff(n,x)
!            integer n
!            real*8 x(n)
!		end function ff
!    	end INTERFACE	  

!
!      o Integer  n  is the parameter space dimension, i.e., the number
!        of adjustable parameters. 
!
!      o Structure mctrl is a set of control flags and parameters, to
!        control the behavior of the MCMC, and also printed
!        output.  A default value will be used for any control variable
!        which is not set following declaration.  The
!        elements of mctrl and their defaults are:
!
!     outputMpcty	logical		(default false)
!        Flag indicating whether to print repeated points to
!        output files (false) or to just print unique points with
!        multiplicities (true)
!
!     propRefinLen	integer		(default 200)
!        Length over which to calculate the variance of the
!        sampled fitness function.  This variance is used to 
!        set the width of the propodal function at every 
!        RevPropLen-th step.
!
!     revPropLen	integer		(default 10)
!        The number of steps to use the current proposal function
!        for before revising it.
!
!     fileunit		integer		(default 89)
!	 The unit number for the output file
!
!     initialWidths(n)	double prec(n)	(defaults 0.01)
!        Initial partial widths of the proposal function in each
!        of the n parameter directions.
!
!     chainame		string		(default mcchain.out)
!        The chain output path and filename 

!      o Function  ff  is a user-supplied scalar function of n vari-
!        ables, which must have the calling sequence f = ff(n,x), where
!        x is a real parameter array of length n.  This function must
!        be written so as to bound all parameters to the interval [0,1];
!        that is, the user must determine a priori bounds for the para-
!        meter space, and ff must use these bounds to perform the appro-
!        priate scalings to recover true parameter values in the
!        a priori ranges.
!
!        By convention, ff should return higher values for more optimal
!        parameter values (i.e., individuals which are more "fit").
!        For example, in fitting a function through data points, ff
!        could return the inverse of chi**2.
!
!        In most cases initialization code will have to be written
!        (either in a driver or in a separate subroutine) which loads
!        in data values and communicates with ff via one or more labeled
!        common blocks.
!
!     Output:
      double precision, intent(OUT) :: x(n)
      real*8, intent(OUT) :: f
      integer, intent(OUT) :: status
!
!      o Array  x(1:n)  is the "fittest" (optimal) solution found,
!         i.e., the solution which maximizes fitness function ff
!
!      o Scalar  f  is the value of the fitness function at x
!
!      o Integer  status  is an indicator of the success or failure
!         of the call to MCMC (0=success; non-zero=failure)
!
!     Working variables
      logical:: stepOn
      integer :: i, j, stepCount=0, mpcty, f_prev
      double precision :: proposalWidths(n)
      double precision :: parmin(n), parmax(n), x_current(n), x_temp(n)
      real*8 :: f_current, start_tot, finish_tot


	parmin(1) = 0.d0
	parmax(1) = 0.5

	parmin(2) = 0.d0
	parmax(2) = 180

	parmin(3) = 0.d0
	parmax(3) = 360

      !Set status to fail utill success
      status = 1

      !Open output file for writing
      open(unit=mctrl%fileunit, file = trim(mctrl%chainame), status = 'REPLACE')

      !Randomly populate the history and choose a random starting point 
        do j = 1, n
          x_current(j) = urand()
        enddo
        do i = 1,n
          x_temp(i) = x_current(i)*(parmax(i)-parmin(i))+parmin(i)
        enddo
        x_temp(1) = 0.07
        x_temp(2) = 112
        x_temp(3) = 224
      f_current = ff(n,x_temp)
      x = x_temp
      f = f_current
      !write(*,*)x_current
      !Initialise the proposal widths
      proposalWidths = mctrl%initialWidths
      !proposalWidths(1) = 0.05
      proposalWidths(2) = 25 * DEG2RAD
      proposalWidths(3) = 25 * DEG2RAD
       
      !do j = 1, 10000
      do j = 1, 50000
      !write(*,*)x_current
        !do i = 1,n
        !  x_temp(i) = x_current(i)*(parmax(i)-parmin(i))+parmin(i)
        !enddo

 call cpu_time(start_tot)

        stepOn = MetroHastyStep(ff,f_current,n,x_temp,proposalWidths)

        write(mctrl%fileunit,fstr2) (x_temp(i), i=1,n),f_current
        !Update the best fit
        if (f .gt. f_current) then
          f = f_current
          x = x_temp
        endif

 call cpu_time(finish_tot)
 write(*,*)'finish-start (MetroHastyStep) =',finish_tot-start_tot
          
      enddo

      close(mctrl%fileunit)
      write(*,*) 'Chain has converged (!)'

end subroutine MCMC


logical function MetroHastyStep(ff,f,n,x,proposalWidths)
!     Does a single Metropolis-Hastings MCMC step

      integer, intent(IN) :: n
      integer :: i
      double precision, intent(IN)    :: proposalWidths(n)
      double precision, intent(INOUT) :: x(n)
      double precision                :: xnew(n), proposedJump(n), alpha, theta, phi
      double precision, dimension(3,3):: M_0, M_delta
      double precision, dimension(3)  :: nv
      real*8, intent(INOUT) :: f
      real*8 :: ff, f_new
      external ff

      alpha = x(1)
      theta = x(2) * DEG2RAD
      phi = x(3) * DEG2RAD

      nv(1) = sin(theta)*cos(phi)
      nv(2) = sin(theta)*sin(phi)
      nv(3) = cos(theta)

      !write(*,*)'theta = ',theta
      !write(*,*)'phi = ',phi
      !call compute_euler_matrix_zyz(phi, theta, 0.d0, M_0)
      !nv = matmul(M_0,matmul(transpose(M_0),nv))
      !xnew(2) = acos(nv(3))
      !xnew(3) = atan2(nv(2),nv(1))
      !if (xnew(3) .lt. 0.d0) xnew(3) = xnew(3) + 2*pi
      !write(*,*)'theta 2 = ',xnew(2)
      !write(*,*)'phi 2 = ',xnew(3)

      do 
        !do i = 1, n
        !  proposedJump(i) = gaussdev()*proposalWidths(i)
        !enddo
          proposedJump(1) = gaussdev()*proposalWidths(1)
        xnew(1) = alpha + proposedJump(1)
        !if (all(xnew .le. 1.d0 .and. xnew .ge. 0.d0)) exit
        if ((xnew(1) .le. 1.d0) .and. (xnew(1) .ge. 0.d0)) exit
      enddo

        !do i = 2, n
        !  proposedJump(i) = gaussdev()*proposalWidths(i)
        !enddo

          proposedJump(2) = gaussdev()*proposalWidths(2)
          !proposedJump(3) = gaussdev()*proposalWidths(3)
          proposedJump(3) = urand()*2*pi

      call compute_euler_matrix_zyz(phi, theta, 0.d0, M_0)
      call compute_euler_matrix_zyz(proposedJump(3), proposedJump(2), 0.d0, M_delta)
      nv = matmul(M_0,matmul(M_delta,matmul(transpose(M_0),nv)))
      xnew(2) = acos(nv(3))
      xnew(3) = atan2(nv(2),nv(1))
      if (xnew(3) .lt. 0.d0) xnew(3) = xnew(3) + 2*pi
      xnew(2) = xnew(2) / DEG2RAD
      xnew(3) = xnew(3) / DEG2RAD

      f_new = ff(n,xnew)

      if (debug) then
        write(*,*)
        write(*,*) 'f_old, f_new:'
        write(*, *)f
        write(*, *)f_new
      endif

      if (f_new .le. f) then
        !Accept new point immediately if the new likelihood is greater than the old
        MetroHastyStep = .true.
        if (debug) write(*,*) 'New point accepted: f_new > f_old'
      else
        !If the new likelihood is less or equal to the old, use the likelihood ratio
        !and a random number to decide if the new point is accepted
        if (log(urand()) .le. (f-f_new)) then
          MetroHastyStep = .true.
          if (debug) write(*,*) 'New point accepted, but f_new < f_old'
        else
          MetroHastyStep = .false.
          if (debug) write(*,*) 'New point rejected, f_new < f_old'
        endif
      endif

      if (MetroHastyStep) then
        !Update the likelihood and parameter values if the step is accepted
        f = f_new
        x = xnew
      endif

end function MetroHastyStep

! Convention: First psi around z, then theta around y, then phi around z
subroutine compute_euler_matrix_zyz(phi, theta, psi, euler_matrix)
    implicit none

    real(dp),                 intent(in)  :: phi, theta, psi
    real(dp), dimension(3,3), intent(out) :: euler_matrix

    real(dp) :: sphi, cphi, sth, cth, spsi, cpsi

    sphi = sin(phi)
    cphi = cos(phi)

    sth  = sin(theta)
    cth  = cos(theta)

    spsi = sin(psi)
    cpsi = cos(psi)

    euler_matrix(1,1) = -sphi * spsi + cth * cphi * cpsi
    euler_matrix(1,2) = -sphi * cpsi - cth * cphi * spsi
    euler_matrix(1,3) =                sth * cphi
    euler_matrix(2,1) =  cphi * spsi + cth * sphi * cpsi
    euler_matrix(2,2) =  cphi * cpsi - cth * sphi * spsi
    euler_matrix(2,3) =                sth * sphi
    euler_matrix(3,1) =              - sth * cpsi
    euler_matrix(3,2) =                sth * spsi
    euler_matrix(3,3) =                cth

end subroutine compute_euler_matrix_zyz

!subroutine convert_euler_matrix_to_angles_zyz(M, phi, theta, psi)
!    implicit none

!    real(dp), dimension(3,3), intent(in)  :: M
!    real(dp),                 intent(out) :: phi, theta, psi

!    real(dp) :: sth, cth, spsi, sphi

!    ! Get theta
!    cth = M(3,3)
!    if (cth > 1.d0) then
!       theta = 0.d0
!    else if (cth < -1.d0) then
!       theta = pi
!    else
!       theta = acos(cth)
!    end if

!    if (abs(cth) < 0.9999998d0) then
!       phi = atan2(M(2,3),M(1,3))
!       psi = -atan2(-M(3,2),-M(3,1))
!    else
!       ! Psi and phi are degenerate; define phi = 0
!       phi = 0.d0
!       psi = atan2(-M(1,2),M(2,2))
!    end if

!end subroutine convert_euler_matrix_to_angles_zyz

double precision function gaussdev()
!     Returns gaussianly-distributed random numbers with unit variance
!     Based on Numerical Recipes routine gasdev

      logical, save :: gotSpare=.false.
      double precision, save :: spare
      double precision :: fac, rsq, v1, v2

      if (gotSpare) then
        gotSpare = .false.
        gaussdev = spare
      else
        do
          v1 = 2.d0*urand() - 1.d0			!Get a point in the unit square
          v2 = 2.d0*urand() - 1.d0
          rsq = v1*v1 + v2*v2
          if (rsq .lt. 1.d0 .and. rsq .ne. 0.d0) exit	!See if it is in the unit circle
        enddo
        fac = sqrt(-2.d0*log(rsq)/rsq)			!Do Box-Muller transformation
        spare = v1 * fac
        gotSpare = .true.
        gaussdev = v2 * fac
      endif

end function gaussdev


double precision pure function stdev(data,n)

! Calculates standard deviation of data
! Based on Numerical Recipes moment function

      integer, intent(IN) :: n
      integer :: j
      double precision, intent(IN) :: data(n)
      double precision :: ave,var,p,s,ep

      if (n.le.1) then 
        stdev = 0.d0
        return
      endif
      s=sum(data)
      ave=s/n
      var=0.
      ep=0.
      do j=1,n
        s=data(j)-ave
        ep=ep+s
        p=s*s
        var=var+p
      enddo
      var=(var-ep**2/n)/(n-1)
      stdev=sqrt(var)
      return
      
end function stdev

function ran2(seed)
!=====================================================================
      INTEGER seed,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (seed.le.0) then
        seed=max(-seed,1)
        idum2=seed
        do 11 j=NTAB+8,1,-1
          k=seed/IQ1
          seed=IA1*(seed-k*IQ1)-k*IR1
          if (seed.lt.0) seed=seed+IM1
          if (j.le.NTAB) iv(j)=seed
11      continue
        iy=iv(1)
      endif
      k=seed/IQ1
      seed=IA1*(seed-k*IQ1)-k*IR1
      if (seed.lt.0) seed=seed+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=seed
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)

end function ran2

function urand()
!=====================================================================
!     Return the next pseudo-random deviate from a sequence which is
!     uniformly distributed in the interval [0,1]
!
!     Uses the function ran0, the "minimal standard" random number
!     generator of Park and Miller (Comm. ACM 31, 1192-1201, Oct 1988;
!     Comm. ACM 36 No. 7, 105-110, July 1993).
!=====================================================================
      implicit none
!
!     Input - none
!
!     Output
      real*8     urand
!
!     Local
      integer  iseed, rndcount
      !real*8     ran2
      !external ran2
!
!     Common block to make iseed visible to rninit (and to save
!     it between calls)
      common /rnseed/ iseed, rndcount
      !write(*,*)'iseed = ',iseed
      !write(*,*)'rndcount = ',rndcount
      !if (rndcount == 0) then
        !CALL SYSTEM_CLOCK(COUNT=iseed)
        !iseed = -1
        !rndcount = 1
      !endif
!
      !write(*,*)'iseed = ',iseed
      !write(*,*)'rndcount = ',rndcount
      !urand = ran2(iseed)
       call random_number(urand)
      !write(*,*)'urand = ',urand
      
end function urand

end module MCMCmodule
