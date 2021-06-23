    module BayOp
    use CalcLike
    use bobyqa_module
    use ParamPointSet
    use HelperRoutines

    implicit none
    private

	Type, extends(THelperRoutines) :: TGaussianRegressor
        real(mcp) :: prior_amp = 1.0_mcp
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

	
! In 2D we have multiple kernels each with their own correlation length. The kernel only has 1 hyperparameter, the amplitude gets multiplied in the end
	function TGaussianRegressor_Kernel(this, X, Y, hypers) result(kernel)
	class(TGaussianRegressor) :: this
	real(mcp), intent(in), dimension (:) :: X, Y
	real(mcp), dimension(size(X), size(Y)) :: kernel
	real(mcp), intent(in) :: hypers
	integer :: i,j

	! kernel = exp(-0.5 ((x-y)/l)^2
	do j=1,size(Y)
			do i=1, size(X)
					kernel(i,j) = exp(-0.5_mcp * ( (X(i) - Y(j)) / hypers)**2 )
			end do
	end do

	end function TGaussianRegressor_Kernel

! full subroutine for 2D Gaussian
	subroutine TGaussianRegressor_GPR(this, hypers, XData, YData, XPred, alpha, mu, std, input_dim, YMax)
	class(TGaussianRegressor) :: this
	real(mcp), intent(in), dimension(:) :: hypers
	real(mcp), intent(in), dimension(:,:) :: XData	! x values of data
	real(mcp), intent(in), dimension(:) :: YData 	! likelihood at given XData
	real(mcp), intent(in), dimension(:,:) :: XPred 	! x values we make predictions for
	real(mcp), intent(inout) :: YMax
	integer, intent(in) :: input_dim ! dimensions of input
	real(mcp), dimension( size(YData),size(YData) ) :: K, L	! Kernel of data and Cholesky decompositon
	real(mcp), dimension(size(YData), size(Xpred,1)) :: Ks, v	! Kernel(data,pred)
	real(mcp), dimension(:) :: alpha, mu, std
	real(mcp) :: amplitude	! amplitude of kernel = std of data
	integer :: i 

        amplitude = (this%prior_amp * hypers(input_dim + 1))**2 
	! K = kernel(XData,XData)
	do i=1, input_dim
		if (i == 1) then
			K = amplitude * this%kernel(Reshape( XData(:,i), [size(YData)] ), Reshape( XData(:,i), [size(YData)] ), hypers(i))
		else
			K = K * this%kernel(Reshape( XData(:,i), [size(YData)] ), Reshape( XData(:,i), [size(YData)] ), hypers(i))
		end if
	end do
	! add jitter on the diagonal to make cholesky stable
	do i=1, size(YData)
		K(i,i) = K(i,i) + amplitude * 1e-10_mcp
	end do
	
	! Ks = kernel(XData,XPred)
	do i=1, input_dim
		if (i == 1) then
			Ks =amplitude * this%kernel(Reshape( XData(:,i), [size(YData)] ), Reshape( XPred(:,i), [size(XPred,1)] ), hypers(i))
		else
			Ks = Ks * this%kernel(Reshape( XData(:,i), [size(YData)] ), Reshape( XPred(:,i), [size(XPred,1)] ), hypers(i))
		end if
	end do
	! Choleksy decomposition of K
	L = this%Cholesky(K)

	! Calculate mu and std
	alpha = this%solvvec(TRANSPOSE(L), this%solvvec(L,Ydata)) ! K^-1 Ydata 
	v = this%solvmat(L, Ks)! L^-1 * Ks
	mu = Reshape(MATMUL(TRANSPOSE(Ks), Reshape(alpha,[size(YData),1])),[size(XPred,1)]) !Ks * K^-1 * Ydata = mu
 	do i=1, size(XPred,1)
	        std(i)=SQRT( amplitude - DOT_PRODUCT( v(:,i),v(:,i) ) )
	end do   
	YMax = MAXVAL(YData) ! Corresponding maximum of data. Get here so that mu and YMax are consistent with same normalisation 
	end subroutine TGaussianRegressor_GPR

! Calling Bobyqa, but this time with 3 parameter to be fixed
	subroutine TGaussianRegressor_Optimise(this, hypers, XData, YData, XPred,  mu, std, n_dim, bobyyes, YMax)
	class(TGaussianRegressor) :: this
	real(mcp), intent(inout), dimension(:) :: hypers
	real(mcp), intent(in) :: Xdata(:,:), Ydata(:), XPred(:,:) 
	real(mcp) :: Ydata_norm(size(Ydata)), mean_norm, std_norm ! Data value with mean 0
	integer, intent(in) :: n_dim
	logical, intent(inout) :: bobyyes
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
	integer :: i, j, n, npt

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
	
	npt = (n+1)*(n+2)/2
!	npt =  2*n + 1 ! there are some conditions on this parameter, but no idea. Just try. Higher npt = more runs
	call this%Normalize_Y(Ydata, Ydata_norm, mean_norm, std_norm)

	this%prior_amp = std_norm ! rescaling the hyper for the error to be of order 1
	if (bobyyes) then
		call BOBYQA(n, npt, hypers, xl, xu, rhobeg, rhoend, iprint, maxfun, calfun)

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
	end if
	call this%GPR(hypers, XData, YData_norm, XPred, alpha, mu, std, n_dim, YMax)
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
	
	amplitude = (this%prior_amp  * hypers(n_dim + 1))**2
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
	
	f = NegLogLike(theta, XData, YData_norm, XPred, n_dim)
	
	
	end subroutine calfun
	end subroutine TGaussianRegressor_Optimise

!*****************************************************************************************
!>	TBayOptimisor

! Initialise the LikeCalculator to used in calculation, but actually no idea what it does
    subroutine TBayOptimisor_Init(this, LikeCalculator)
    class(TBayOptimisor) :: this
    class(TLikeCalculator), target:: LikeCalculator
    this%LikeCalculator => LikeCalculator
    end subroutine TBayOptimisor_Init	

! calculates expected improvement
	subroutine TBayOptimisor_EI(this, mu, YMax, std, xi, EI)
	class(TBayOptimisor) :: this
	real(mcp), intent(in), dimension(:) :: mu, std
	real(mcp), intent(in) :: YMax, xi
	real(mcp), dimension(size(mu)) :: Z, pdf, cdf
	real(mcp), intent(out), dimension(size(mu)) :: EI

	Z = (mu - YMax - xi)/std
	
	call this%normal_pdf(Z,pdf)
	call vdcdfnorm(size(Z),Z, cdf)

	EI = (mu - YMax - xi) * cdf + std * pdf
	Write(*,*) 'Maximal Value of EI is: ', MaxVal(EI)
	Write(*,*) 'Maximal Improvment of Mu is: ', MaxVal(mu)-YMax
	Write(*,*) 'Expected Improvement at this point is: ', EI(MAXLOC(mu,1))
end subroutine TBayOptimisor_EI

! find maximumal expected improvement and returns next_X as position for the next measurement	
	subroutine TBayOptimisor_propose(this, mu, YMax, std, xi, XPred ,next_X, EI)
	class(TBayOptimisor) :: this
	real(mcp), intent(in), dimension(:,:) :: XPred
	real(mcp), intent(in), dimension(:) :: mu, std
	real(mcp), intent(inout), dimension(:) :: EI
	real(mcp), intent(in) :: YMax, xi
	real(mcp), dimension(size(mu)) :: pdf, cdf
	real(mcp), intent(out), dimension(1,size(XPred,2)) :: next_X
	integer :: Location
	
	call this%EI(mu, YMax, std, xi, EI)
	Location = MAXLOC(EI, 1)
	next_X(1,:) = XPred(  Location, :)
        Write(*,*) 'Expected Improvement at YMax is: ', XPred(MAXLOC(mu,1), :)
	end subroutine TBayOptimisor_propose
	
	
! Master routine for this module
	subroutine TBayOptimisor_BayOp(this, Params, XData, XPred, hypers, which_param, XBest, xi, Cutoff_EI, n_iterations, n_random, YMax, OutputName)
	class(TBayOptimisor) :: this
	Type(ParamSet) Params
	real(mcp), intent(inout), allocatable :: XData(:,:), XPred(:,:), XBest(:,:), hypers(:)
	real(mcp), intent(in) :: xi, Cutoff_EI
	real(mcp), intent(inout) :: YMax
	integer, intent(in) :: n_iterations, n_random, which_param(:)
	character(len=*), intent(inout) :: OutputName
	character(len=40) :: PredOutput
	real(mcp), allocatable, dimension(:) :: Ydata, YData_loop, mu, std, Ei
	real(mcp), allocatable :: XData_loop(:,:), Next_X(:,:)
	logical :: bobyyes = .True.
	integer :: MemGB, n_Grid, n_dim, n_max, i, ii, j, ppos, jj 
	real(mcp) :: EI_previous, time0, time1, time2, time3, time4, meanY, stdY
			
	n_Grid = size(XPred, 1)
	n_dim = size(XPred, 2)
	EI_previous = 0	! will be changed after every run
	allocate(YData(n_iterations))
	
	call this%RandomSample(Params, n_random, n_Grid, XPred, XData, YData, which_param)
	
	allocate(next_X(1,n_dim))
	allocate(mu(n_Grid-n_random), std(n_Grid-n_random), EI(n_Grid-n_random))
	
        !ppos = scan(trim(OutputName),".", BACK=.true.)                          
	!PredOutput = OutputName(1:ppos-1) // '_Pred.txt'                                                          
	!write(*,*) 'PredOutput is', PredOutput
	
	do i=n_random+1, n_iterations	
		call CPU_time(time0)
		
		write(*,*) 'Getting Sample Number:', i
		! Allocating the loop quantity and assigning value(using mostly empty arrays leads to problems when inverting)
		allocate(YData_loop(i-1))
		allocate(XData_loop(i-1, n_dim))
		YData_loop(1:i-1) = YData(1:i-1)
		do ii=1, n_dim
			XData_loop(1:i-1, ii) = XData(1:i-1, ii)
		end do
		! write results to output file
		if (n_dim == 3) then
			call this%CreateOutput(XData_loop(:,1), XData_loop(:,2), YData_loop, XData_loop(:,3), XData_loop(:,3), XData_loop(:,3), OutputName)
		else if (n_dim == 4) then
			call this%CreateOutput(XData_loop(:,1), XData_loop(:,2), YData_loop, XData_loop(:,3), XData_loop(:,4), XData_loop(:,4), OutputName)
		else if (n_dim == 5) then
			call this%CreateOutput(XData_loop(:,1), XData_loop(:,2), YData_loop, XData_loop(:,3), XData_loop(:,4), XData_loop(:,5), OutputName)
		end if
		call CPU_time(time1)
		! Run GPR
		jj=i-n_random
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
		call this%Optimise(hypers, XData_loop, YData_loop, XPred, mu, std, n_dim, bobyyes, YMax)
		!call this%CreateOutput(XPred(:,1), XPred(:,2), mu, PredOutput)                
  		
		deallocate(YData_loop, XData_loop)

		call CPU_time(time2)
		if (Feedback > 0) write(*,*) 'Time to propose with to bobyqa+GPR: ', time2-time0
		! Use EI to find next best point
		call this%propose(mu, YMax, std, xi, XPred ,next_X, Ei) 
		call CPU_time(time3)
		if (Feedback > 0) write(*,*) 'Time to propose with BayOp: ', time3-time2
		! Add next point to XData
		do ii=1, n_dim
			XData(i,ii) = next_X(1,ii)
			Params%P(which_param(ii)) = next_X(1,ii)
		end do
		write(*,*) 'Sampling the likelihood at:', next_X
		! Sample Data	
		YData(i) = this%LogLike(Params)
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
				call this%RemoveGrid(XPred, EI, n_dim, Cutoff_EI)	! If therre is too much change in EI, don't remove points. Fit might be bad
		else if (ABS(MaxVal(EI) - EI_previous) < 1.0_mcp) then			! Remove very bad points immediately
				call this%RemoveGrid(XPred, EI, n_dim, Cutoff_EI * 1e-50_mcp)
		else if (i>2*n_random .and. ABS(MaxVal(EI) - EI_previous) < 1.0_mcp) then	! Gradually start removing points strictlier
				call this%RemoveGrid(XPred, EI, n_dim, Cutoff_EI * 1e-15_mcp)
		else if (i>3*n_random .and. ABS(MaxVal(EI) - EI_previous) < 1.0_mcp) then
				call this%RemoveGrid(XPred, EI, n_dim, Cutoff_EI * 1e-5_mcp)
		end if		
! Update the Max_EI for the next run
		EI_previous = MaxVal(EI)
		! Remove Sampled point
		n_max = MAXLOC(EI,1)
		call this%RemoveData(XPred, n_max, n_dim)
				
		deallocate( mu, std, EI )
		allocate( mu(size(XPred,1)), std(size(XPred,1)), EI(size(XPred,1))) ! allocate quantities for next run		
		
		! Finish loop
		write(*,*) n_Grid-size(mu), 'out of', n_Grid, ' points have been removed from the grid'
		call CPU_time(time4)
		if (Feedback > 0) write(*,*) 'Time for one round of BayOp: ', time4-time0
		if (size(XPred,1)==0) exit	! End loop if there are no more points to sample
		if (MAXVAL(EI)< 1e-20_mcp .and. EI_previous<1e-20) then
			write(*,*) 'MAXVAL(EI)< 1e-20_mcp) ending BayOp before 500 samples'
			exit
		end if
		if (Feedback > 0) call this%Memory(MemGB)
		j=i	! If loop finishes early, there will be empty entries in arrays. Integer j counts how many non-empty entries there are.
	end do
	! Get best point
	n_max = MAXLOC(YData(1:j),1) ! We need to restrict values from 1 to j, because 0 values are larger than the negative likelihood values
	XBest(1,:) = XData(n_max, :) 
	YMax = MAXVAL(YData(1:j))
	end subroutine TBayOptimisor_BayOp
	
	! Master routine for this module
	subroutine TBayOptimisor_Sampler(this, Params, input_dim, n_random, use_refine, OutputName)
	class(TBayOptimisor) :: this
	Type(ParamSet) Params
	real(mcp), allocatable, dimension(:) :: hypers
	real(mcp), allocatable, dimension(:,:) :: XData, XPred, XBest
	real(mcp) :: xi, Cutoff_EI, YMax, YMaxOld
	integer, parameter :: n_iterations = 3000
	integer, allocatable, dimension(:) :: which_param, n_input ! Parameter varied, 17 = As
	integer, intent(in) :: input_dim 
	integer, intent(inout) ::n_random
	logical, intent(in) :: use_refine
	character(len=*), intent(inout) :: OutputName
	integer :: n_Grid, ii, n_max, n_refine, ppos, i
	
	allocate( which_param(input_dim))
	allocate( n_input(input_dim) )
		
	write(*,*) 'Number of parameters sampled over: ', input_dim
	if (input_dim == 2) then
		which_param = [25, 26]
	else if (input_dim == 3) then
		which_param = [25, 26, 27]
	else if (input_dim == 4) then
		which_param = [25, 26, 27, 28]
        else if (input_dim == 5) then
                which_param = [25, 26, 27, 28, 29] 
	else
		call MpiStop('Only BayOpDimension 2-5 is supported. Change in .ini file')
	end if

	print*, 'Parameters sampled over are:'
	if (input_dim == 2) then
		print*, 'Amplitude, Frequency'
	else if (input_dim == 3) then
		print*, 'Amplitude, Frequency, Phase'
	else if (input_dim == 4) then
		print*, 'Ampltidue, Frequency, Phase, NewP5'
	else if (input_dim == 5) then                                 
		print*, 'Ampltidue, Frequency, Phase, NewP4, NewP5' 
	end if
	do ii=1, input_dim
		print*, BaseParams%PWidth(which_param(ii))
	end do

	if (Params%P(30) == 1) then
		print*, "Using Linear Modulation to Power Spectrum"
	else if (Params%P(30) == 2) then                                                                                                    
		print*, "Using Logarithmic Modulation to Power Spectrum"
	else if (Params%P(30) == 3) then                                                                                                    
		print*, "Using Running Logarithmic Modulation to Power Spectrum"      
	else if (Params%P(30) == 4) then                                                                                                    
		print*, "Using Primordial Standard Clock to Power Spectrum for expansions"   
	else if (Params%P(30) == 6) then
		print*, "Using Primordial Standard Clock to Power Spectrum for contractions (Bounce) with w=p*w/k^(1/p)"
        else if (Params%P(30) == 5) then                                                                                                
		print*, "Using Primordial Standard Clock to Power Spectrum for contractions (Pyro) with w=p*w/k_p^(1/p)"                                                         
	else
		print*, "No idea what you are sampling, error incoming"
	end if

	! Derive some integers and parameter lengths needed to calculate the Grid
	do ii=1, input_dim
		n_input(ii) = ( BaseParams%PMax(which_param(ii)) - BaseParams%Pmin(which_param(ii)) )/BaseParams%PWidth(which_param(ii)) + 1
	end do 
	if (input_dim == 1) then
		n_max = n_input(1)
	else
		n_max = MAXVAL(n_input) 
	end if
	n_Grid = n_input(1)
	do ii=2, input_dim
		n_Grid = n_input(ii) * n_Grid
	end do
	
	! allocate input parameter
	allocate( Xdata(n_iterations, input_dim))
	allocate( XPred(n_Grid, input_dim) )
	allocate( hypers(input_dim + 1) )
	do ii=1, input_dim+1
		hypers(ii)=1
	end do
	allocate( XBest(1, input_dim) )

	call this%Grid(n_max, n_input, input_dim, which_param, XPred)
	write(*,*) 'Parameter ranges are:' 	
	i = n_input(1)
	do ii=1, input_dim
		if (ii==1) then
			write(*,*) ii, XPred(1,ii), XPred(i,ii)
		else
			i = i * n_input(ii)
			write(*,*) ii, XPred(1,ii), XPred(i,ii)
		end if
	end do
	xi=0.0_mcp
	Cutoff_EI=1e-5_mcp
	call this%BayOp(Params, XData, XPred, hypers, which_param, XBest, xi, Cutoff_EI, n_iterations, n_random, YMax, OutputName)
	write(*,*) "Best fit point at:", XBest
	write(*,*) "Likelihood at best fit:", YMax
	
	! Make refined search around maximum
	if (use_refine) then
		Cutoff_EI=1e-5_mcp
		YMaxOld = YMax
		n_refine = 21
		deallocate(XPred)
		allocate(XPred(n_refine**input_dim, input_dim))
		call this%Grid_refine(input_dim, n_refine, which_param, XBest, XPred)
	        write(*,*) 'Parameter ranges are:'
       		do ii=1, input_dim
			write(*,*) ii, XPred(1,ii), XPred(n_refine**ii,ii)
		end do
		ppos = scan(trim(OutputName),".", BACK=.true.)		
		if (ii > 0) OutputName = OutputName(1:ppos-1) // '_refine.txt'
		print*, 'New OutputName is', OutputName
		call this%BayOp(Params, XData, XPred, hypers, which_param, XBest, xi, Cutoff_EI, n_iterations, n_random, YMax, OutputName)
		
		write(*,*) "Best fit point at:", XBest
		write(*,*) "Likelihood at best fit:", YMax
		write(*,*) "Improvement due to refinement:", YMax-YMaxOld
	end if
	do i=1, input_dim
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
