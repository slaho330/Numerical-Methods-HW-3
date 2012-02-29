program heatIM
  implicit none
  real, parameter :: l = 1
  integer, parameter :: nd = 100
  real, parameter :: alpha = 0.001
  real, parameter :: dt = 0.005
  integer, parameter :: nsteps = 1
  real            :: d, C, bc
  real, dimension(nd+1, nd+1, nd+1)  :: x

  print*, 'nsteps:', nsteps
  d = l/nd !distance between grid points is l/nd
  bc = 0 !the constant value for the boundary conditions
  C = alpha * dt/(d**2) 

  call initialize(d, nd, x, bc)
  call Jacobi(x, nd, C, nsteps)
  call initialize(d, nd, x, bc)
  call GaussSeidel(x, nd, C, nsteps)
  call initialize(d, nd, x, bc)
  call SORelax(x, nd, C, nsteps)
  print*, "done!"
end program heatIM

subroutine initialize(d, nd, x, bc)
  implicit none
  integer :: i, j, k, l, nd
  real    :: d, xx, yy, zz, bc
  real, dimension(nd+1, nd+1, nd+1) :: x

  !initialize
  do i = 2, nd
     do j = 2, nd
        do k = 2, nd
           xx = d * (i-1)
           yy = d * (j-1)
           zz = d * (k-1)
           x(i,j,k) = exp(-(5*xx - 2.5)**2) * exp(-(5*yy - 2.5)**2) * exp(-(5*zz - 2.5)**2)
        enddo
     enddo
  enddo

  !boundary conditions
  do i = 1, nd+1
     do j = 1, nd+1
        x(1, i, j) = bc
        x(nd+1, i, j) = bc
     enddo
  enddo  
  do i = 1, nd+1
     do j = 1, nd+1
        x(i, 1, j) = bc
        x(i, nd+1, j) = bc
     enddo
  enddo
  do i = 1, nd+1
     do j = 1, nd+1
        x(i, j, 1) = bc
        x(i, j, nd+1) = bc
     enddo
  enddo

end subroutine initialize

subroutine Jacobi(x, nd, C, nsteps)
  implicit none
  integer :: i,j,k,n,nsteps,m,mits,nd, clock_rate, clock_max, t1, t2
  real :: C, mean
  real, dimension(nd+1, nd+1, nd+1) :: xold, x, xnew
  
  xold = x        !set the initial x vector to xold
  mits = 1000     !maximum number of iterations at each nstep

  do n = 1, nsteps
     call system_clock(t1, clock_rate, clock_max)
     x = xold !initial guess
     print*, 'step:', n
     do m = 1, mits !iterate until m has reached the max iterations
        do i = 2, nd
           do j = 2, nd
              do k = 2, nd
                 !for Jacobi IM, the value of x^m+1 depends on x^m (m is iteration number)
                 xnew(i,j,k) = C/(6*C+1)*(x(i-1,j,k) + x(i+1,j,k) + x(i,j-1,k) + x(i,j+1,k)&
                      + x(i,j,k-1) + x(i,j,k+1)) + 1/(6*C+1) * xold(i,j,k)
              enddo
           enddo
        enddo

        !check the error
        mean = sum(abs(x-xnew))/((nd+1)**3) !mean of the errors
        if (mean.lt.(.000001)) then !if the mean is less than threshold, next nstep
           EXIT
        endif
        x = xnew  !if mean over threshold, xnew becomes the new guess for iterating
     enddo
     call system_clock(t2, clock_rate, clock_max)
     print*,'Jacobi timer:'
     write(*,*) real(t2-t1)/real(clock_rate)

     if (m.ge.mits) then !check if a solution for the timestep is not found in mits
        print*, "reached max iteration with error over threshhold"
        EXIT
     endif

     xold = x !x is set as xold, the x value for this time step.
  enddo
end subroutine Jacobi

subroutine GaussSeidel(x, nd, C, nsteps)
  implicit none
  integer :: i,j,k,n,nsteps,m,mits,nd,t1,t2,clock_rate,clock_max
  real :: C, mean
  real, dimension(nd+1, nd+1, nd+1) :: xold, x, xnew
  
  xold = x
  mits = 1000
   
  do n = 1, nsteps
     call system_clock(t1, clock_rate, clock_max)
     x = xold !initial guess 
     print*, 'step:', n
     do m = 1, mits !iterate until error is small between two guesses
        do i = 2, nd
           do j = 2, nd
              do k = 2, nd
                 !for Gauss Seidel, the xnew value at each point is calculated using only xnew, not lastx
                 xnew(i,j,k) = C/(6*C+1)*(xnew(i-1,j,k) + xnew(i+1,j,k) + xnew(i,j-1,k) + xnew(i,j+1,k)&
                      + xnew(i,j,k-1) + xnew(i,j,k+1)) + 1/(6*C+1) * xold(i,j,k)
              enddo
           enddo
        enddo
        !check the error (here x is last x, so the error is between this iteration and the last)
        mean = sum(abs(x-xnew))/((nd+1)**3) !mean of the error
        if (mean.lt.(.000001)) then !if error is acceptably small then next nstep
           EXIT
        endif
        x = xnew
     enddo
     call system_clock(t2, clock_rate, clock_max)
     print*,'Gauss-Seidel timer:'
     write(*,*) real(t2-t1)/real(clock_rate)

     if (m.ge.mits) then !check if a solution for the timestep is not found in mits
        print*, "reached max iteration with error over threshhold"
        EXIT
     endif

     xold = x !set the guess chosen by error calc as xold for next nstep
  enddo
end subroutine GaussSeidel

subroutine SORelax(x, nd, C, nsteps)
  implicit none
  integer :: i,j,k,n,nsteps,m,mits,nd,t1,t2,clock_rate,clock_max
  real :: C, mean, w
  real, dimension(nd+1, nd+1, nd+1) :: xold, x, xnew
  
  w = 1.65 !omega for scaling between the two components of P
  xold = x
  mits = 1000

  do n = 1, nsteps
     call system_clock(t1, clock_rate, clock_max)
     x = xold !initial guess
     print*, 'step:',n
     do m = 1, mits
        do i = 2, nd
           do j = 2, nd
              do k = 2, nd
                 !for this algorithm, x is x^m and xnew is x^(m+1) from lecture
                 !xnew at each point is a scaled sum of x at the point and xnew around the point
                 xnew(i,j,k) = (1-w)*x(i,j,k) + (w*C)/(6*C+1)*(xnew(i-1,j,k) + xnew(i+1,j,k) + xnew(i,j-1,k)&
                      + xnew(i,j+1,k) + xnew(i,j,k-1) + xnew(i,j,k+1)) + w/(6*C+1) * xold(i,j,k)
              enddo
           enddo
        enddo
        !check the error
        mean = sum(abs(x-xnew))/((nd+1)**3) !mean of the errors
        if (mean.lt.(.000001)) then !if error is acceptibly small
           EXIT
        endif
        x = xnew !set xnew as the next guess for x
     enddo
     call system_clock(t2, clock_rate, clock_max)
     print*,'S.O.Relaxation timer:'
     write(*,*) real(t2-t1)/real(clock_rate)   
     
     xold = x !x(x guess) becomes the value of x at this nstep (xold), next nstep 
  enddo
end subroutine SORelax
