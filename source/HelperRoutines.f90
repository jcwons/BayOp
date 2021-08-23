	module HelperRoutines
! In this module all the routine that are needed for running the code are in. This contain linear algebra functions using LAPACK,
! as well as simple functions to create the grid and some other statistic functions. There is no hard coding in here.
! Only if you want to change the routines used to solve the linear equations. Otherwise changes should all be done in BayOp.f90
! DPOTRF is used for the Cholesky decomposition
! DGESV is used to solve the linear equation. (Note that this is not the most efficient method, 
! but doesn't really matter and I was too lazy to change)

	use CalcLike
	use ParamPointSet
    	implicit none
    	private

	Type, extends(TLikelihoodUser) :: THelperRoutines
	real(mcp) :: const_pi = 3.1415926535897932384626433832795_mcp
!	real(mcp) :: baseline = -12272.602152412439_mcp	! Planck TTTEE, simall, low-EE
!        real(mcp) :: baseline = -12277.05138222037_mcp  ! TTTEE + lensing
!	real(mcp) :: baseline = -1184.788_mcp	! For the test likelihood
	real(mcp) :: baseline 
	contains
	procedure :: Initialise => THelperRoutines_Initialise
!	procedure :: Memory => THelperRoutines_Memory ! Memory function for testing
	procedure :: Linspace => THelperRoutines_Linspace
	procedure :: Linspace_refine => THelperRoutines_Linspace_refine
	procedure :: GridMaker => THelperRoutines_GridMaker
	procedure :: Grid => THelperRoutines_Grid
	procedure :: Grid_refine => THelperRoutines_Grid_Refine
	procedure :: RandomSample => THelperRoutines_RandomSample
	procedure :: PrintMatrix => THelperRoutines_PrintMatrix ! optional/only for testing	
	procedure :: Cholesky => THelperRoutines_Cholesky
	procedure :: LogLike => THelperRoutines_LogLike
	procedure :: SolvVec => THelperRoutines_SolvVec
	procedure :: SolvMat => THelperRoutines_SolvMat
	procedure :: Normalize_Y => THelperRoutines_Normalize_Y 
	procedure :: CreateOutput => THelperRoutines_CreateOutput
	procedure :: OutputProfile => THelperRoutines_OutputProfile
	procedure :: RemoveData => THelperRoutines_RemoveData
	procedure :: normal_pdf => THelperRoutines_normal_pdf
	procedure :: RemoveGrid => THelperRoutines_RemoveGrid
	end Type THelperRoutines
	
	public THelperRoutines 
	contains
!*****************************************************************************************
!>	THelperRoutines

! Initialise BayOp. Get dimensions for quanitites, get grid, print some information
	subroutine THelperRoutines_Initialise(this, Params, hypers, XData, XPred, XBest, YMax, n_iterations, which_param, n_input, input_dim, n_random, n_Grid, ii, n_max, n_refine, i)
	class(THelperRoutines) :: this
	Type(ParamSet) Params
	real(mcp), intent(inout), allocatable, dimension(:) :: hypers
	real(mcp), intent(inout), allocatable, dimension(:,:) :: XData, XPred, XBest
	real(mcp), intent(inout) ::  YMax
	integer, intent(in) :: n_iterations 
	integer, intent(inout), allocatable, dimension(:) :: which_param, n_input ! Parameter varied, 17 = As
	integer, intent(in) :: input_dim 
	integer, intent(inout) ::n_random
	integer, intent(inout) :: n_Grid, ii, n_max, n_refine, i


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
	
	end subroutine THelperRoutines_Initialise
! Subroutine to call current Memory usage
!	 subroutine THelperRoutines_Memory(this, valueRSS)
!	 use ifport !if on intel compiler
!	 class(THelperRoutines) :: this
!	 integer, intent(out) :: valueRSS
!	 character(len=200):: filename=' '
!	 character(len=80) :: line
!	 character(len=8)  :: pid_char=' '
!	 integer :: pid
!	 logical :: ifxst

!	 valueRSS=-1    ! return negative number if not found

	 !--- get process ID

!	 pid=getpid()
!	 write(pid_char,'(I8)') pid
!	 filename='/proc/'//trim(adjustl(pid_char))//'/status'

	! !--- read system file

!	 inquire (file=filename,exist=ifxst)
!	 if (.not.ifxst) then
!	   write (*,*) 'system file does not exist'
!	   return
!	 endif
!
!	 open(unit=100, file=filename, action='read')
!	 do
!	   read (100,'(a)',end=120) line
!	   if (line(1:6).eq.'VmRSS:') then
!		  read (line(7:),*) valueRSS
!		  exit
!	   endif
!	 enddo
!	 120 continue
!	 close(100)
!	 write(*,*) valueRSS , 'Used memory'
!	 return
!	 end subroutine THelperRoutines_Memory
	
! Subroutine creating the Grid in higher dimension
	subroutine THelperRoutines_Linspace(this, Array, input_dim, which_param, n_input)
	class(THelperRoutines) :: this
	real(mcp), intent(inout) :: Array(:,:) ! dimension will be (max(n_input),input_dim)
	integer, intent(in), dimension(:) :: which_param ! which parameters are sampled over
	integer, intent(in) :: input_dim ! how many parameter are sampled over
	integer, intent(in), dimension(:) :: n_input ! grid size of each sample parameter
	integer :: nn, i

	! The outer loop runs over the input dimension and selects different columns
	! The inner loop runs over the rows and fills up with equidistant steps
	! Note that different columns will have different lenght depending on the choices in test.ini
	do nn=1, input_dim
		do i=1, n_input(nn)-1 
			Array(i, nn)= BaseParams%Pmin(which_param(nn)) + (i-1)*BaseParams%PWidth(which_param(nn))
		end do
		Array(n_input(nn), nn) = BaseParams%PMax(which_param(nn)) ! Last value turns out weird, so assign it manually
	end do
	end subroutine THelperRoutines_Linspace

! Grid width of refine is defined in here	
	subroutine THelperRoutines_Linspace_refine(this, Array, XBest, input_dim, which_param, n_refine)
	class(THelperRoutines) :: this
	real(mcp), intent(inout) :: Array(:,:) ! dimension will be (max(n_input),input_dim)
	real(mcp), intent(in) :: XBest(:,:)
	integer, intent(in), dimension(:) :: which_param ! which parameters are sampled over
	integer, intent(in) :: input_dim, n_refine
	integer :: j, i

	do j=1, input_dim
		do i=1, n_refine
			! If loop assures that new grid is inside prior range
			if ( XBest(1,j)- 2 * BaseParams%PWidth(which_param(j)) .LE. BaseParams%PMin(which_param(j)) ) then
				Array(i, j) = BaseParams%PMin(which_param(j)) + 4 * (i-1) * BaseParams%PWidth(which_param(j))/(n_refine-1)
				Array(1, j) = BaseParams%PMin(which_param(j))
			else if ( XBest(1,j)+ 2 * BaseParams%PWidth(which_param(j)) .GE. BaseParams%PMax(which_param(j)) ) then
				Array(i, j) = BaseParams%PMax(which_param(j)) - 4 * (n_refine-i) * BaseParams%PWidth(which_param(j))/(n_refine-1)
				Array(n_refine, j) = BaseParams%PMax(which_param(j))
			else
				Array(i, j) = XBest(1,j) - 2 * BaseParams%PWidth(which_param(j)) + 4 *(i-1) * BaseParams%PWidth(which_param(j))/(n_refine-1)
			end if
		end do
	end do
	
	end subroutine THelperRoutines_Linspace_refine
	
! Grid making routine turns two arrays(:,:) into a grid. size(b,1) has to be 1
	subroutine THelperRoutines_GridMaker(this, a, b, c)
	class(THelperRoutines) :: this
	real(mcp), intent(in) :: a(:,:), b(:,:)	
	real(mcp), allocatable, intent(inout) :: c(:,:)
	integer :: i, j, m

	allocate(c( size(a,1) * size(b,1), size(a,2)+1))
	m=1
	do j=1, size(b,1)
		do i=1, size(a,1)
			c(m, 1:size(a,2)) = a(i,:)
			c(m, size(c,2)) = b(j,1)
			m = m + 1
		end do
	end do
	end subroutine THelperRoutines_GridMaker	
		
! Grid for input parameters using Max, Min and Width provided from .ini file
	subroutine THelperRoutines_Grid(this, n_max, n_input, input_dim, which_param, XPred)
	class(THelperRoutines) :: this
	real(mcp), intent(inout) :: XPred(:,:)
	integer, intent(in) :: input_dim, n_max, which_param(:)
	integer, intent(in), dimension(:) :: n_input
	real(mcp), allocatable :: Array(:,:), A(:,:), temp(:,:), temp2(:,:)
	integer :: i,j,k,m,o,p
	! Allocating the temporary array containing all coordinates (not all combinations)
	allocate(Array(n_max,input_dim))
	! Creating the temporary array
	call this%Linspace(Array, input_dim, which_param, n_input)
	allocate(A (n_input(1),1) )
	do i=1, n_input(1)
		A(i,1) = Array(i,1)
	end do 
	do i=2, input_dim
		If (allocated(temp)) then
			deallocate(temp)
		end if
		allocate(temp2(n_input(i),1))
		do j=1, n_input(i)
			temp2(j,1) = Array(j,i)
		end do
		call this%GridMaker(A, temp2, temp)
		deallocate(A, temp2)
		allocate(A(size(temp,1),size(temp,2)))
		A = temp
	end do
	Xpred = temp
	write(*,*)'Steps in each dimension:', n_input
	write(*,*)'Total number of grid points', shape(XPred)
	end subroutine THelperRoutines_Grid
	
! 	Grid for refined !!! Note: Should be able to reuse Grid, but I didn't write the subroutine general enough smarter to write new routine
	subroutine THelperRoutines_Grid_Refine(this, input_dim, n_refine, which_param, XBest, XPred)
	class(THelperRoutines) :: this
	real(mcp), intent(inout) :: XPred(:,:)
	real(mcp), intent(in) :: XBest(:,:)
	integer, intent(in) :: input_dim, n_refine, which_param(:)
	real(mcp), allocatable :: Array(:,:), A(:,:), temp(:,:), temp2(:,:)
	integer :: i,j,k
	
	allocate(Array(n_refine, input_dim))
	call this%Linspace_refine(Array,XBest, input_dim, which_param, n_refine)
	
	allocate(A(n_refine,1))
	do i=1, n_refine
		A(i,1) = Array(i,1)
	end do
	do i=2, input_dim
		if (allocated(temp)) then
			deallocate(temp)
		end if
		allocate(temp2(n_refine,1))
		do j=1, n_refine
			temp2(j,1)= Array(j,i)
		end do
		call this%GridMaker(A, temp2, temp)
		deallocate(A, temp2)
		allocate(A(size(temp,1),size(temp,2)))
		A = temp
	end do
	XPred = temp
	end subroutine THelperRoutines_Grid_Refine

! Create random point to sample from Grid
	subroutine THelperRoutines_RandomSample(this, Params, n_random, n_Grid, XPred, XData, YData, which_param)
	class(THelperRoutines) :: this
	class(ParamSet) :: Params
	real(mcp), intent(inout), allocatable :: XPred(:,:), XData(:,:), YData(:)
	real(mcp) :: rng, mean, std
	integer, intent(in) :: n_Grid, n_random, which_param(:) ! no need for n_Grid
	integer :: i, k , ii
	

	write(*,*) 'Taking', n_random,  'random samples from the Grid before starting the Baysian Optimisation'
	do i=1, n_random
		print*,i
	! Create random number between 1 and n_Grid
		call RANDOM_NUMBER(rng) ! gives random number between 1 and 0
		k = FLOOR(size(XPred,1)*rng) + 1
		XData(i,:) = XPred(k,:)
	! Calculate likelihood at random point k
		do ii=1, size(which_param)
			Params%P(which_param(ii))=XData(i,ii)
		end do
		YData(i) = this%LogLike(Params)
		call this%RemoveData(XPred, k, size(which_param))
		if (YData(i) < -1e10) then
			if (i==1) then
				mean = 0
				std = 0
			else if (i==2) then
				mean = sum(Ydata(1:i-1))/(i-1)
				std = 0
			else
				mean = sum(Ydata(1:i-1))/(i-1)
				std = SUM(Ydata(1:i-1)**2) / (i-2)
			end if
			std = SQRT(std)
			YData(i) = mean - std
			if (YData(i)>0.0) then
				YData(i)=-5
			end if
		end if 
                if (Feedback>0) then
                        write(*,*) 'Sampling Likelihood at:', XData(i,:)
                        write(*,*) 'Likelihood is:', Ydata(i)
                end if
	end do
	! change Parameters back to normal
	do i=1, size(which_param)
		Params%P(which_param(i)) = BaseParams%center(which_param(i))
	end do
	end subroutine THelperRoutines_RandomSample	

! Helper Routine to print matrix output nicer
	subroutine THelperRoutines_PrintMatrix(this, A)
	class(THelperRoutines) :: this
	real(mcp), intent(in), dimension(:,:)::A
	integer :: i,j,m,n
	n = size(A,1)
	m = size(A,2)
	do, i=1,n
			write(*,*) ( A(i,j), j=1,m)
	enddo
	end subroutine THelperRoutines_PrintMatrix
	
! Cholesky Decomposition into lower triangular matrix
	function THelperRoutines_Cholesky(this, M) result(A)
	class(THelperRoutines) :: this
	real(mcp), intent(in), dimension(:,:) :: M
	real(mcp), dimension(size(M,1),size(M,2)) :: A
	integer :: n, info,i,j
	character :: UPLO
	n = size(M,1)
	UPLO = 'L' ! change to 'U'for upper triangular matrix
	! A will be overwritten by output of dpotrf, so copy
	A = M
	call dpotrf(UPLO, N, A, N, info)
	!check for success
	if (UPLO == 'L') then
		do i=1, n
			do j=i+1,n
				A(i,j)=0
			end do
		end do
		else
		do i=1, n
			do j=i+1,n
				A(j,i)=0
			end do
		end do
	end if
	if (info /= 0) then
		write (0,*) "solve: dpotrf returned an error code (", info, ")"
		error stop
	end if
	end function THelperRoutines_Cholesky

! Returns the LogLike of given Params 
	function THelperRoutines_LogLike(this, Params) result(logLike)
    class(THelperRoutines) :: this
	Type(ParamSet) Params, Trial
    real(mcp) :: logLike
	Trial = Params
	logLike = -this%LikeCalculator%GetLogLike(Trial) - this%baseline
!	if (logLike<-100) then
!		logLike = -100
!	end if
	call Trial%Clear(Params)
	end function THelperRoutines_LogLike
	
! Using dgesv to solve linear equation M*x=b with x, b rank 1
	function THelperRoutines_SolvVec(this, M, b) result(x)
	class(THelperRoutines) :: this
	real(mcp), intent(in) :: M(:,:)
	real(mcp) :: b(:)
	real(mcp), dimension(size(M,1),size(M,2)) :: A
	real(mcp), dimension(size(b)) :: x
	integer :: N, NRHS, info
	integer :: ipiv(size(B,1))

	N = size(M, 1)
        NRHS = 1
	! A and x are overwritten on output of dgesv, so copy
	A = M
	x = b
	!call DTRTRS('L', 'N', 'N', N, NRHS, A, N, x, N, info) 
	call dgesv(N, 1, A, N, ipiv, x, N, info)
	! check for success
	if (info /= 0) then
		write (0,*) "solve: dgesv returned an error code (", info, ")"
		error stop
	end if
	end function THelperRoutines_SolvVec
	
! Using dgesv to solve linear equation M*X=B with X,B rank 2
	function THelperRoutines_SolvMat(this, M, B) result(X)
	class(THelperRoutines) :: this
	real(mcp), intent(in) :: M(:,:)
	real(mcp), intent(in) :: B(:,:)
	real(mcp), dimension(size(M,1),size(M,2)) :: A
	real(mcp), dimension(size(B,1),size(B,2)) :: X
	integer :: N, NRHS,  info
	integer :: ipiv(size(B,1))

	N = size(M, 1)
	NRHS = size(B,2)
	! A and x are overwritten on output of dgesv, so copy
	A = M
	X = B
	!call gesv(A, B)
	call dgesv(N, NRHS, A, N, ipiv, X, N, info)
	!call DTRTRS('L', 'N', 'N', N, NRHS, A, N, x, N, info)
	!check for success
	if (info /= 0) then
		write (0,*) "solve: dgesv returned an error code (", info, ")"
		error stop
	end if
	end function THelperRoutines_SolvMat

	subroutine THelperRoutines_Normalize_Y(this, Ydata, Ydata_norm, mean, std)
	class(THelperRoutines) :: this
	real(mcp), intent(in) :: Ydata(:)
	real(mcp), intent(inout) :: Ydata_norm(size(Ydata))
	real(mcp), intent(out) :: mean, std
	
	mean = sum(Ydata)/size(Ydata)
	Ydata_norm = Ydata - mean
        if (Feedback>1) write(*,*) 'mean of data is: ', mean
        std = SUM(Ydata_norm**2) / (SIZE(Ydata)-1) ! this is actually the variance
	std = SQRT(std) 
        if (Feedback>1) write(*,*) 'std of data is: ', SQRT(std)     

	end subroutine THelperRoutines_Normalize_Y

! Creates output file for 3-5 dimensional array1 and 1D array2
	subroutine THelperRoutines_CreateOutput(this, array1, array2, text)
	class(THelperRoutines) :: this
	integer :: i
	real(mcp), intent(in) :: array1(:,:), array2(:)
	character(len =*) :: text
	open(1, file=text)
	do i=1, size(array2)
		if (size(array1,2)==3) then
			write(1,'(6E15.5)') array1(i,1), array1(i,2), array1(i,3), array2(i)
		else if (size(array1,2)==4) then
			write(1,'(6E15.5)') array1(i,1), array1(i,2), array1(i,3), array1(4,i), array2(i)
		else if (size(array1,2)==5) then
			write(1,'(6E15.5)') array1(i,1), array1(i,2), array1(i,3), array1(i,4), array1(i,5), array2(i)
		end if
	end do
	close(1)
	end subroutine THelperRoutines_CreateOutput

! Creates output file for 3-5 dimensional array1 and 1D array2, array3
        subroutine THelperRoutines_OutputProfile(this, array1, array2, array3, text)
        class(THelperRoutines) :: this
        integer :: i
        real(mcp), intent(in) :: array1(:,:), array2(:), array3(:)
        character(len =*) :: text
        open(1, file=text)
        do i=1, size(array2)
                if (size(array1,2)==3) then
                        write(1,'(6E15.5)') array1(i,1), array1(i,2), array1(i,3), array2(i), array3(i)
                else if (size(array1,2)==4) then
                        write(1,'(6E15.5)') array1(i,1), array1(i,2), array1(i,3), array1(4,i), array2(i), array3(i)
                else if (size(array1,2)==5) then
                        write(1,'(6E15.5)') array1(i,1), array1(i,2), array1(i,3), array1(i,4), array1(i,5), array2(i), array3(i)
                end if
        end do
        close(1)
        end subroutine THelperRoutines_OutputProfile

! Remove Sampled point from Grid
	subroutine THelperRoutines_RemoveData(this, XPred, n, input_dim)
	class(THelperRoutines) :: this
	real(mcp), intent(inout), allocatable, dimension(:,:) :: XPred 
	integer , intent(in) :: n, input_dim
	integer :: nn
	real(mcp), allocatable, dimension(:,:) :: temp
	
	allocate(temp(size(XPred,1), size(XPred,2))) ! Copy XPred, because it will be overwritten
	temp = XPred
	deallocate(XPred)
	allocate( XPred( size(temp,1)-1,size(temp,2) ) ) ! reduce size of XPred by 1 as we remove 1 value
	do nn=1, input_dim ! loop over all input dimensions
		XPred(:,nn) = [temp(1:n-1,nn), temp(n+1:size(temp,1),nn)]
	end do
	end subroutine THelperRoutines_RemoveData

! gives probability density function of the normal Gauss
	subroutine THelperRoutines_normal_pdf(this, X, pdf)
	class(THelperRoutines) :: this
	real(mcp), dimension(:) :: pdf
	real(mcp), dimension(:) :: X
	
	pdf = exp( -0.5_mcp * X * X)/ Sqrt(2.0_mcp * this%const_pi)
	return
	end subroutine THelperRoutines_normal_pdf


! Removes points with low expected Improvement
	subroutine THelperRoutines_RemoveGrid(this, grid, EI, input_dim, Cutoff_EI)
	class(THelperRoutines) :: this
	real(mcp), allocatable, intent(inout), dimension(:,:) :: grid
	real(mcp), allocatable, intent(inout), dimension(:) :: EI
	real(mcp), allocatable, dimension(:,:) :: temp_grid
	real(mcp), allocatable, dimension(:) :: temp_array
	real(mcp), intent(in) :: Cutoff_EI
	integer, intent(in) :: input_dim
	integer :: ii
	
	allocate(temp_grid(size(grid,1), size(grid,2))) ! Copy XPred, because it will be overwritten
	temp_grid = grid
	deallocate(grid)
	do ii=1, input_dim ! loop over all input dimensions
		temp_array = PACK(temp_grid(:,ii), EI>Cutoff_EI) ! Remove all points for which EI is below cutoff
		if (.not. allocated(grid)) allocate(grid(size(temp_array),input_dim) ) ! allocate XPred during the first run
		grid(:,ii)=temp_array
		deallocate(temp_array)
	end do
	!reduce also the size of EI
	allocate(temp_array(size(EI)))
	temp_array = EI
	deallocate(EI)
	EI = PACK(temp_array, temp_array>Cutoff_EI)
	end subroutine THelperRoutines_RemoveGrid

end module HelperRoutines
