    module BayOp
    use CalcLike
    use bobyqa_module
    use ParamPointSet
    use HelperRoutines

    implicit none
    private

	Type, extends(THelperRoutines) :: TGaussianRegressor
        real(mcp) :: prior_amp = 1.0_mcp	! prior widht will be changed by bobyqa
	contains
	procedure :: Kernel => TGaussianRegressor_Kernel
	procedure :: GPR => TGaussianRegressor_GPR
	procedure :: Optimise => TGaussianRegressor_Optimise
	end Type TGaussianRegressor

    Type, extends(TGaussianRegressor) :: TBayOptimisor
 	contains 
	procedure :: Init => TBayOptimisor_Init
	procedure :: EI => TBayOptimisor_EI
	procedure :: propose => TBayOptimisor_propose
	procedure :: BayOp => TBayOptimisor_BayOp
	procedure :: Sampler => TBayOptimisor_Sampler
	end Type TBayOptimisor
	
	public  TGaussianRegressor, TBayOptimisor 
	contains

!*****************************************************************************************
!>	TGaussianRegressor
!	
!	Class with the routines for GPR.
!	A Gaussian kernel with n+1 parameters is used.
!	The hyperparameters are determined by finding the MLE using the bobyqa minimizer.
!	
!	System of linear equations are solved with MKL routines	
	
	function TGaussianRegressor_Kernel(this, X, Y, hypers) result(kernel)
	class(TGaussianRegressor) :: this
	real(mcp), intent(in), dimension (:) :: X, Y
	real(mcp), dimension(size(X), size(Y)) :: kernel
	real(mcp), intent(in) :: hypers
	integer :: i,j
	! Returns gaussian kernel without the signal variance in front. Only hyperparameter is correlation length.
	! Prior width is added later

	! kernel = exp(-0.5 ((x-y)/l)^2
	do j=1,size(Y)
			do i=1, size(X)
					kernel(i,j) = exp(-0.5_mcp * ( (X(i) - Y(j)) / hypers)**2 )
			end do
	end do

	end function TGaussianRegressor_Kernel

	subroutine TGaussianRegressor_GPR(this, hypers, XData, YData, XPred, alpha, mu, std, n_dim, YMax)
	class(TGaussianRegressor) :: this
	real(mcp), intent(in), dimension(:) :: hypers
	real(mcp), intent(in), dimension(:,:) :: XData	! x values of data
	real(mcp), intent(in), dimension(:) :: YData 	! likelihood at given XData
	real(mcp), intent(in), dimension(:,:) :: XPred 	! x values we make predictions for
	real(mcp), intent(inout) :: YMax
	integer, intent(in) :: n_dim ! dimensions of input
	real(mcp), dimension( size(YData),size(YData) ) :: K, L	! Kernel of data and Cholesky decompositon
	real(mcp), dimension(size(YData), size(Xpred,1)) :: Ks, v	! Kernel(data,pred)
	real(mcp), dimension(:) :: alpha, mu, std
	real(mcp) :: amplitude	! amplitude of kernel = std of data
	integer :: i 
	! Routine to make predictions via Gaussian Process Regression given the hyper parameters
	! Is called after Optimise when hypers have been determined
	
	amplitude = (this%prior_amp * hypers(n_dim + 1))**2 ! scaled with std of data to get hypers approx. 1
	
	! K = kernel(XData,XData)
	do i=1, n_dim
		if (i == 1) then
			K = amplitude * this%kernel(Reshape( XData(:,i), [size(YData)] ), Reshape( XData(:,i), [size(YData)] ), hypers(i))
		else
			K = K * this%kernel(Reshape( XData(:,i), [size(YData)] ), Reshape( XData(:,i), [size(YData)] ), hypers(i))
		end if
	end do
	! Add jitter on the diagonal to make Cholesky decomposition stable
	do i=1, size(YData)
		K(i,i) = K(i,i) + amplitude * 1e-10_mcp ! amplitude * 1e-10_mcp is only for stability. Change this value to add noise in function
	end do
	
	! Ks = kernel(XData,XPred)
	do i=1, n_dim
		if (i == 1) then
			Ks =amplitude * this%kernel(Reshape( XData(:,i), [size(YData)] ), Reshape( XPred(:,i), [size(XPred,1)] ), hypers(i))
		else
			Ks = Ks * this%kernel(Reshape( XData(:,i), [size(YData)] ), Reshape( XPred(:,i), [size(XPred,1)] ), hypers(i))
		end if
	end do
	
	! Choleksy decomposition of K
	L = this%Cholesky(K)

	! Calculate mu and std by solving linear systems
	alpha = this%solvvec(TRANSPOSE(L), this%solvvec(L,Ydata)) ! K^-1 Ydata 
	v = this%solvmat(L, Ks)! L^-1 * Ks
	mu = Reshape(MATMUL(TRANSPOSE(Ks), Reshape(alpha,[size(YData),1])),[size(XPred,1)]) !Ks * K^-1 * Ydata = mu
 	do i=1, size(XPred,1)
	        std(i)=SQRT( amplitude - DOT_PRODUCT( v(:,i),v(:,i) ) )
	end do   
	YMax = MAXVAL(YData) ! Corresponding maximum of data. Get here so that mu and YMax are consistent with same normalisation 
	end subroutine TGaussianRegressor_GPR

	subroutine TGaussianRegressor_Optimise(this, hypers, XData, YData, XPred,  mu, std, n_dim, bobyyes, YMax)
	class(TGaussianRegressor) :: this
	real(mcp), intent(inout), dimension(:) :: hypers
	real(mcp), intent(in) :: Xdata(:,:), Ydata(:), XPred(:,:) 
	real(mcp) :: Ydata_norm(size(Ydata)), mean_norm, std_norm ! Data value with mean 0
	integer, intent(in) :: n_dim
	logical, intent(inout) :: bobyyes ! if true run bobyqa, if false skip bobyqa and only do GPR
	real(mcp), intent(inout) :: YMax
	real(mcp), dimension(size(Ydata)) :: alpha ! K^-1 Ydata
	real(mcp), dimension(:) :: mu, std
	real(mcp) :: nll, Trace_L
	! Bobyqa parameter
	real(mcp), dimension(size(hypers)) :: xl, xu   !xl,xu are the upper and lower boundaries
	real(mcp),parameter :: bdl    = 1e-5_mcp ! lower boundary passed to xl
	real(mcp),parameter :: bdu    = 5e2_mcp ! upper boundary passed to xu
	integer,parameter  :: iprint = 1 ! 0 for no output, 1-3 for more output
	integer,parameter  :: maxfun = 1e5 ! maximum number of calls to function
	real(mcp),parameter :: rhobeg = 5e0_mcp ! precision of first run
	real(mcp),parameter :: rhoend = 1e-5_mcp ! final precision
	real(mcp) :: time0, time1
	integer :: i, j, n, npt
	! Calling Bobyqa to find hypers via MLE using BOBYQA solver
	! Afterwards GPR is called


	n = size(hypers)
	do  i = 1,n ! asigning values to boundaries
		xl (i) = bdl
		xu (i) = bdu
	end do
	! Rerun bobyqa if hyperparameter equals boundary
	do i=1, n
		if (hypers(i) == bdu .or. hypers(i) == bdl) then
			bobyyes = .true.
		end if
	end do
	
	npt = (n+1)*(n+2)/2 ! there are some conditions on this parameter, but no idea. Higher npt = more runs
!	npt =  2*n + 1 ! This is most commonly used, but often made bobyqa end in local,non-global maximum
	
	! rescaling the hyper for the error to be of order 1
	call this%Normalize_Y(Ydata, Ydata_norm, mean_norm, std_norm)
	this%prior_amp = std_norm 
	
	
	if (bobyyes) then ! devide whether to run bobyqa or not
		call CPU_time(time0)
		call BOBYQA(n, npt, hypers, xl, xu, rhobeg, rhoend, iprint, maxfun, calfun)
		! if hypers end at prior boundaries, re-run bobyqa with different initial values
		do j=1, n
			if (hypers(j) == bdu) then
				do i=1,n
					hypers(i) = bdl
				end do
				call BOBYQA(n, npt, hypers, xl, xu, rhobeg, rhoend, iprint, maxfun, calfun)
			else if (hypers(j) == bdl) then
	                        do i=1,n
                                        hypers(i) = bdl
                                end do 
				call BOBYQA(n, npt, hypers, xl, xu, rhobeg, rhoend, iprint, maxfun, calfun)
			end if
		end do
		call CPU_time(time1)
		if (Feedback > 0) write(*,*) 'Time for Bobyqa: ', time1-time0
	end if
	call CPU_time(time0)
	call this%GPR(hypers, XData, YData_norm, XPred, alpha, mu, std, n_dim, YMax)
	call CPU_time(time1)
	if (Feedback > 0) write(*,*) 'Time for GPR: ', time1-time0
	write(*,*) 'Best fit point is: ', MAXVAL(YData)
	
	contains 
	
	function NegLogLike(hypers, XData, YData, XPred, n_dim) result(nll)
	real(mcp) :: nll
	real(mcp), intent(in) :: XData(:,:), YData(:), XPred(:,:)
	real(mcp), intent(in) :: hypers(:)
	integer, intent(in) :: n_dim
	integer :: i
	real(mcp) :: amplitude, Trace_L
	real(mcp), dimension( size(YData),size(YData) ) :: K,L
	! Calculates the maximum likelihood estimation (MLE)
	amplitude = (this%prior_amp * hypers(n_dim+1))**2
	do i=1, n_dim
		if (i == 1) then
			K = amplitude * this%kernel(Reshape( XData(:,i), [size(YData)] ), Reshape( XData(:,i), [size(YData)] ), hypers(i))
		else
			K = K * this%kernel(Reshape( XData(:,i), [size(YData)] ), Reshape( XData(:,i), [size(YData)] ), hypers(i))
		end if
	end do
	do i=1, size(YData)
		K(i,i) = K(i,i) + amplitude*1e-10_mcp
	end do
	
	L = this%cholesky(K)
	Trace_L=0
	do i=1,size(Ydata)
		Trace_L = Trace_L + Log(L(i,i))
	end do
	
	alpha = this%solvvec(TRANSPOSE(L), this%solvvec(L,Ydata)) ! K^-1 Ydata
	
	nll = DOT_PRODUCT(Ydata,alpha)/2 + Trace_L 		
	end function NegLogLike
	
	subroutine calfun(n, theta, f)
	integer, intent(in)	:: n
	real(mcp) :: nll
	real(mcp), dimension(:), intent(in) :: theta
	real(mcp), intent(out) :: f
	! Bobyqa requires a certain shape of the function that is called
	f = NegLogLike(theta, XData, YData_norm, XPred, n_dim)
	
	
	end subroutine calfun
	end subroutine TGaussianRegressor_Optimise


!*****************************************************************************************
!>	TBayOptimisor
!
!	This class contains the routines for Bayesian Optimisation
!	Expected Improvement is used for the acqusition function
!	At the bottom, wrapper routine is called
!

    subroutine TBayOptimisor_Init(this, LikeCalculator)
    class(TBayOptimisor) :: this
    class(TLikeCalculator), target:: LikeCalculator
	! Initialise the LikeCalculator to used in calculation
    this%LikeCalculator => LikeCalculator
    end subroutine TBayOptimisor_Init	

	subroutine TBayOptimisor_EI(this, mu, YMax, std, xi, EI)
	class(TBayOptimisor) :: this
	real(mcp), intent(in), dimension(:) :: mu, std
	real(mcp), intent(in) :: YMax, xi
	real(mcp), dimension(size(mu)) :: Z, pdf, cdf
	real(mcp), intent(out), dimension(size(mu)) :: EI
	! Calculates Expected Improvement

	Z = (mu - YMax - xi)/std
	
	call this%normal_pdf(Z,pdf)
	call vdcdfnorm(size(Z),Z, cdf)

	EI = (mu - YMax - xi) * cdf + std * pdf
	if (Feedback > 0) Write(*,*) 'Maximal Value of EI is: ', MaxVal(EI)
	if (Feedback > 1) Write(*,*) 'Maximal Improvment of Mu is: ', MaxVal(mu)-YMax
	if (Feedback > 1) Write(*,*) 'Expected Improvement at this point is: ', EI(MAXLOC(mu,1))
	end subroutine TBayOptimisor_EI

	subroutine TBayOptimisor_propose(this, mu, YMax, std, xi, XPred ,next_X, EI)
	class(TBayOptimisor) :: this
	real(mcp), intent(in), dimension(:,:) :: XPred
	real(mcp), intent(in), dimension(:) :: mu, std
	real(mcp), intent(inout), dimension(:) :: EI
	real(mcp), intent(in) :: YMax, xi
	real(mcp), dimension(size(mu)) :: pdf, cdf
	real(mcp), intent(out), dimension(1,size(XPred,2)) :: next_X
	integer :: Location
	! Find maximumal expected improvement and returns next_X as position for the next measurement	

	call this%EI(mu, YMax, std, xi, EI)
	Location = MAXLOC(EI, 1)
	next_X(1,:) = XPred(  Location, :)
        if (Feedback > 1) Write(*,*) 'Highest Ei is at: ', XPred(MAXLOC(mu,1), :)
	end subroutine TBayOptimisor_propose
	
	
	subroutine TBayOptimisor_BayOp(this, Params, XData, XPred, hypers, which_param, XBest, xi, Cutoff_EI, n_iterations, n_random, YMax, OutputName, do_profile)
	class(TBayOptimisor) :: this
	Type(ParamSet) Params
	real(mcp), intent(inout), allocatable :: XData(:,:), XPred(:,:), XBest(:,:), hypers(:)
	real(mcp), intent(in) :: xi, Cutoff_EI
	real(mcp), intent(inout) :: YMax
	integer, intent(in) :: n_iterations, n_random, which_param(:)
	character(len=*), intent(inout) :: OutputName
	character(len=40) :: OutputProfile	
	real(mcp), allocatable, dimension(:) :: Ydata, YData_loop, mu, std, Ei
	real(mcp), allocatable :: XData_loop(:,:), Next_X(:,:), XPred_profile(:,:)
	logical :: bobyyes = .True.
	logical :: terminate
	logical, intent(in) :: do_profile
	integer :: MemGB, n_grid, n_dim, n_max, i, ii, j, ppos, jj 
	real(mcp) :: EI_previous, time0, time1, time2, time3, time4, meanY, stdY
	! Routine generates random samples and then iterates over GPR+EI for n_iterations or when all points have been removed from grid
			
	n_grid = size(XPred, 1) ! get grid size
	n_dim = size(XPred, 2) ! get dimension
	EI_previous = 0	! will be changed after every run
	allocate(XPred_profile(n_grid, n_dim))
	XPred_profile = XPred ! make copy for profile likelihood clean up
	allocate(YData(n_iterations))
	
	! Get random samples to initalise the GPR
	call this%RandomSample(Params, n_random, n_grid, XPred, XData, YData, which_param)
	
	allocate(next_X(1,n_dim))
	allocate(mu(n_grid-n_random), std(n_grid-n_random), EI(n_grid-n_random))

	! Starting the GPR+EI algorithm
	do i=n_random+1, n_iterations	
		call CPU_time(time0)
		
		write(*,*) 'Getting Sample Number:', i
		! Allocating the loop quantity and assigning value (MKL complained when using mostly empty matrix for DPOTRF)
		allocate(YData_loop(i-1))
		allocate(XData_loop(i-1, n_dim))
		YData_loop(1:i-1) = YData(1:i-1)
		do ii=1, n_dim
			XData_loop(1:i-1, ii) = XData(1:i-1, ii)
		end do
		! write results to output file.
		call this%CreateOutput(XData_loop, YData_loop, OutputName)
		call CPU_time(time1)
		jj=i-n_random

		! Deciding if rerun bobyqa or not
		if (i < 300) then
			bobyyes = .True.
		else if (i > 301 .and. i<501 .and. MOD(i,5)==0) then
			bobyyes = .True.
		else if (i>501 .and. i<1001 .and. MOD(i,50)==0) then
			bobyyes = .True.
		else if (i>1001 .and. i<2501 .and. MOD(i,250)==0) then
			bobyyes = .True.
                else if (i>2501 .and. MOD(i,500)==0) then
			bobyyes = .True.
		else
			bobyyes = .False.
		end if

		! Update hyperparameters(if bobyyes=T) + GPR
		! Note that after this%Optimise we are working with mu and Ymax normalised
		call this%Optimise(hypers, XData_loop, YData_loop, XPred, mu, std, n_dim, bobyyes, YMax)
		!call this%CreateOutput(XPred(:,1), XPred(:,2), mu, PredOutput)                
  		
		deallocate(YData_loop, XData_loop)

		call CPU_time(time2)
		! if (Feedback > 1) write(*,*) 'Time to propose with to bobyqa+GPR: ', time2-time0
		! Use EI to find next best point ! YMax and mu are normalised
		call this%propose(mu, YMax, std, xi, XPred ,next_X, Ei) 
		call CPU_time(time3)
		if (Feedback > 0) write(*,*) 'Time to propose with EI: ', time3-time2
		! Add next point to XData
		do ii=1, n_dim
			XData(i,ii) = next_X(1,ii)
			Params%P(which_param(ii)) = next_X(1,ii)
		end do
		write(*,*) 'Sampling the likelihood at:', next_X
		! Sample Data	
		YData(i) = this%LogLike(Params)
		! For some feature models CAMB fails and returns -1e31 for the likelihood when use HMCode (halofit_version >= 5) in
		! /batch3/params_CMB_default.ini 
		! Turn the value into a bad point one std below the mean (mu-std) instead of -1e31.
		! Otherwise GPR does not work with the -1e31 outlier
		if (YData(i) < -1e10) then
				meanY = sum(Ydata(1:i-1))/(i-1)
				stdY = SUM(Ydata(1:i-1)**2) / (i-2)
				stdY = SQRT(stdY)
				YData(i) = meanY - stdY
				if (YData(i)>0.0) then
						YData(i)=-5
				end if
		end if
		write(*,*) 'Likelihood is:', YData(i)
		! Remove Grid points with low EI	
		if (i>4*n_random .and. ABS(MaxVal(EI) - EI_previous) < 1.0_mcp) then	! Only remove points after having a decent amount of samples 
				call this%RemoveGrid(XPred, EI, n_dim, Cutoff_EI)	! If there is too much change in EI, don't remove points. Fit might be bad
		else if (i<2*n_random .and. ABS(MaxVal(EI) - EI_previous) < 1.0_mcp) then ! Remove very bad points immediately
				call this%RemoveGrid(XPred, EI, n_dim, Cutoff_EI * 1e-50_mcp)
		else if (i<3*n_random .and. ABS(MaxVal(EI) - EI_previous) < 1.0_mcp) then	! Gradually start removing points strictlier
				call this%RemoveGrid(XPred, EI, n_dim, Cutoff_EI * 1e-15_mcp)
		else if (i<4 *n_random .and. ABS(MaxVal(EI) - EI_previous) < 1.0_mcp) then
				call this%RemoveGrid(XPred, EI, n_dim, Cutoff_EI * 1e-5_mcp)
		end if		
		
		! Update the Max_EI for the next run
		EI_previous = MaxVal(EI)
		! Remove Sampled point
		n_max = MAXLOC(EI,1)
		call this%RemoveData(XPred, n_max, n_dim)
		
		! update the size of mu, std, EI to new grid size after points have been removed
		deallocate( mu, std, EI )
		allocate( mu(size(XPred,1)), std(size(XPred,1)), EI(size(XPred,1))) ! allocate quantities for next run		
		
		! Finish loop
		write(*,*) n_grid-size(mu), 'out of', n_grid, ' points have been removed from the grid'
		call CPU_time(time4)
		if (Feedback > 0) write(*,*) 'Time for one round of BayOp: ', time4-time0

		j=i	! If loop finishes early, there will be empty entries in arrays. Integer j counts how many non-empty entries there are.
		! End loop if there are no more points to sample
		if (size(XPred,1)==0) exit	
		if (MAXVAL(EI)< 1e-20_mcp .and. EI_previous<1e-20) then
			write(*,*) 'MAXVAL(EI)< 1e-20_mcp) ending BayOp before 500 samples'
			exit
		end if
		!call this%Memory(MemGB) ! If want to see memory use uncomment function in helper routine
	end do ! BayOp is done
	
	
	! Get best point
	n_max = MAXLOC(YData(1:j),1) ! We need to restrict values from 1 to j, because 0 values are larger than the negative likelihood values
	XBest(1,:) = XData(n_max, :) 
	YMax = MAXVAL(YData(1:j))
	
	! If do_profile is true, additional checks are made
	! The grid initial grid is restored and a few more runs are made
	! Sometimes updated GPR gives new promising candidates
	! Those have to be sampled to get a nice profile likelihood
	! Also returns mu and std for the whole grid
	if (do_profile) then
		bobyyes = .FALSE. ! We do not want to update the hypers just do a clean up
		if (allocated(mu)) deallocate(mu)
		if (allocated(std)) deallocate(std)
		if (allocated(EI)) deallocate(EI)
		allocate( mu(size(XPred_profile,1)), std(size(XPred_profile,1)), EI(size(XPred_profile,1))) ! allocate quantities for next run		
		j = j+1 ! This will be start value for next loop

		do i=j,j+n_random
			jj = i !count number of samples
			write(*,*) 'Cleanup: Getting Sample Number:', i
			allocate(YData_loop(i-1))
			allocate(XData_loop(i-1, n_dim))
			YData_loop(1:i-1) = YData(1:i-1)
			do ii=1, n_dim
				XData_loop(1:i-1, ii) = XData(1:i-1, ii)
			end do
			call this%Optimise(hypers, XData_loop, YData_loop, XPred_profile, mu, std, n_dim, bobyyes, YMax) ! only need GPR actually
	                call this%CreateOutput(XData_loop, YData_loop, OutputName)
			meanY = sum(Ydata_loop)/(i-1) ! Normally we only use normalised mu
			deallocate(YData_loop, XData_loop)
			! Get next data point
			call this%propose(mu, YMax, std, xi, XPred_profile ,next_X, Ei)
				
        	        ! Add next data point 
			do ii=1, n_dim
                	        XData(i,ii) = next_X(1,ii)
                        	Params%P(which_param(ii)) = next_X(1,ii) ! Change parameters for sampling
                	end do		
			write(*,*) 'Sampling the likelihood at:', next_X
			! Sample Data	
			YData(i) = this%LogLike(Params)

			! For some feature models CAMB fails and returns -1e31 for the likelihood. This turns the value into a bad point one std below the mean. Otherwise GPR does not work with the -1e31 outlier
			if (YData(i) < -1e10) then
				meanY = sum(Ydata(1:i-1))/(i-1)
				stdY = SUM(Ydata(1:i-1)**2) / (i-2)
				stdY = SQRT(stdY)
				YData(i) = meanY - stdY
				if (YData(i)>0.0) then
					YData(i)=-5
				end if
			end if
			write(*,*) 'Likelihood is:', YData(i)			
			if (MAXVAL(EI)< 1e-3_mcp) then
				write(*,*) 'Clean up done. Profilelikelihood will be printed'
				exit
			end if

		end do
		mu = mu + meanY ! Remove normalisation of mu for printing
		ppos = scan(trim(OutputName),".", BACK=.true.)					
		OutputProfile = OutputName(1:ppos-1) // '_profile.txt'	
		print*, 'New OutputName for profile:', OutputProfile
		call this%OutputProfile(XPred_profile, mu, std, OutputProfile)
		
		n_max = MAXLOC(YData(1:jj),1)
		XBest(1,:) = XData(n_max, :)
	        YMax = MAXVAL(YData(1:jj))

	end if
	end subroutine TBayOptimisor_BayOp
	

	! Master routine for this module
	subroutine TBayOptimisor_Sampler(this, Params, n_dim, n_random,do_profile, use_refine, OutputName, baseline, xi, Cutoff_EI)
	class(TBayOptimisor) :: this
	Type(ParamSet) Params
	real(mcp), allocatable, dimension(:) :: hypers
	real(mcp), allocatable, dimension(:,:) :: XData, XPred, XBest
	real(mcp), intent(in) :: xi, Cutoff_EI, baseline
	real(mcp) :: YMax, YMaxOld
	integer, parameter :: n_iterations = 3000	! Hard coded maximum of number of iterations of BayOp. Should finish before 1000. Can be increased if necessary.
	integer, allocatable, dimension(:) :: which_param, n_input
	integer, intent(in) :: n_dim 
	integer, intent(inout) ::n_random
	logical, intent(inout) :: use_refine, do_profile
	character(len=*), intent(inout) :: OutputName
	integer :: n_grid, ii, n_max, n_refine, ppos, i
	
	this%baseline = baseline ! Asign baseline value from .ini file to TBayOptimisor class variable
	print*, 'Baseline is: ', this%baseline
	! Calls Initialise in HelperRoutines. This will generate the grid and allocate the parameter. Also decides which parameters are varied in BayOp (which_param)
	call this%Initialise(Params, hypers, XData, XPred, XBest, YMax, n_iterations, which_param, n_input, n_dim, n_random, n_grid, ii, n_max, n_refine, i)
	
	print*, 'Cutoff_EI is:', Cutoff_EI
	! This calls the main routine to start the Bayesian Optimisation. Including random samples, GPR, Bobyqa and EI
	call this%BayOp(Params, XData, XPred, hypers, which_param, XBest, xi, Cutoff_EI, n_iterations, n_random, YMax, OutputName, do_profile)	
	write(*,*) "Best fit point at:", XBest
	write(*,*) "Likelihood at best fit:", YMax
	
	! Make refined search around maximum
	if (use_refine) then
		do_profile = .FALSE.
		n_random = n_random/3 ! less points are needed because function is rather convex in refine area
		YMaxOld = YMax	! Save old value to compare later
		n_refine = 5	! How much smaller should the new stepwidth be. Default is 1/5
		n_refine = n_refine*4+1	! Refine grid is 4 steps around the maximum (2 right, 2 left)
		deallocate(XPred)
		allocate(XPred(n_refine**n_dim, n_dim))	! Allocate new Grid again
		call this%Grid_refine(n_dim, n_refine, which_param, XBest, XPred)	! Create Grid around maximum
		write(*,*) 'Parameter ranges are:'
		do ii=1, n_dim
			write(*,*) ii, XPred(1,ii), XPred(n_refine**ii,ii)
		end do

		! Change output name to differentiate between before and refine
		ppos = scan(trim(OutputName),".", BACK=.true.)					
		if (ii > 0) OutputName = OutputName(1:ppos-1) // '_refine.txt'	
		print*, 'New OutputName is', OutputName
		call this%BayOp(Params, XData, XPred, hypers, which_param, XBest, xi, Cutoff_EI, n_iterations, n_random, YMax, OutputName, do_profile) ! Call BayOp with new parameters
		
		write(*,*) "Best fit point at:", XBest
		write(*,*) "Likelihood at best fit:", YMax
		write(*,*) "Improvement due to refinement:", YMax-YMaxOld
	end if
	! Do some checks if best fit value is at prior boundary. If prior direction is split. This can be ignored. Otherwise there might be better fit outside the prior range
	do i=1, n_dim
		if (XBest(1,i) == BaseParams%Pmin(which_param(i))) then
			write(*,*) 'Warning: Parameter', i, ' has best fit at prior boundary' 
		                write(*,*) 'Consider extending the range'
		else if (XBest(1,i) == BaseParams%PMax(which_param(i))) then
			write(*,*) 'Warning: Parameter', i, ' has best fit at prior boundary'
			write(*,*) 'Consider extending the range'
		end if
	end do
	end subroutine TBayOptimisor_Sampler
	

	end module BayOp 
