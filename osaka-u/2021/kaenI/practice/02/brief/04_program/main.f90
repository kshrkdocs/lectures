!-----------------------------------------------------------------------
program main
!-----------------------------------------------------------------------
  implicit none
!
  !real(8), parameter :: dt = 0.001
  !integer, parameter :: nt = 5000 
  real(8), parameter :: dt = 0.00001
  integer, parameter :: nt = 500000 
!
  real(8) :: k(-2:2), l1, l2, alpha, beta
  real(8) :: x(0:nt), y(0:nt), z(0:nt)
  real(8) :: xana(0:nt)
  real(8) :: xnum(0:nt), ynum(0:nt), znum(0:nt)
  real(8) :: xwol(0:nt)

  real(8) :: c0, c1, c2, c3
  real(8) :: Q, t
  real(8) :: a0, a1, a2, a3, a4, a5, a6, a7
  real(8) :: b0, b1, b2, b3

  integer :: i, j
 

  k(-1) = 2.0
  k( 1) = 1.0
  k(-2) = 3.0
  k( 2) = 4.0

  x(0)  = 1.0
  y(0)  = 3.0
  z(0)  = 2.0
  c0    = x(0) + y(0) + z(0)

  alpha = k(1) + k(-1) + k(2) + k(-2)
  alpha = alpha * 0.5d0
  beta  = k(1) * k(2) + k(1) * k(-2) + k(-1) * k(-2)

  l1 = - alpha + sqrt(alpha**2 - beta)
  l2 = - alpha - sqrt(alpha**2 - beta)

  Q  = k(1)*x(0) - k(-1)*y(0) + l1 * x(0)

  c1 = k(-1) * k(-2) * c0 / beta

  c2 =  x(0) &
      - Q / (l1 - l2) &
      + k(-1)*k(-2)*c0 / (l2*(l1-l2)) &
      - k(-1)*k(-2)*c0 / beta

  c3 =   Q/(l1-l2) &
       - k(-1)*k(-2)*c0 / (l2*(l1-l2))


  ! ... analytical solution ...
  !
  do i = 0, nt
    t = i * dt
    xana(i) = c1 + c2 * exp(l1*t) + c3 * exp(l2*t)
!    write(6,'(2f15.7)') t, xana(i)
  end do

  !
  ! ... numerical solution ...
  !
  xnum(0) = x(0)
  ynum(0) = y(0)
  znum(0) = z(0)

  do i = 0, nt - 1
    xnum(i+1) = xnum(i) &
              - k( 1) * xnum(i) * dt &
              + k(-1) * ynum(i) * dt
    ynum(i+1) = ynum(i) &
              + k(  1)       * xnum(i)*dt &
              - (k(-1)+k(2)) * ynum(i)*dt &
              + k(-2)        * znum(i)*dt
    znum(i+1) = znum(i) &
              + k(2)  * ynum(i) * dt &
              - k(-2) * znum(i) * dt
  end do
!
!
!
  a0 = 1.0d0 / 78.0d0
  a1 = -(69.0d0 + 25.0d0 * sqrt(3.0d0))
  a2 = 216.0d0
  a3 = -69.0d0  + 25.0d0 * sqrt(3.0d0)

  b0 = -(5.0d0 + 2.0d0 * sqrt(3.0d0))
  b1 = 4.0d0 * sqrt(3.0d0)
  b2 =   5.0d0 + 2.0d0 * sqrt(3.0d0) 

  do i = 0, nt
    t = i * dt
    xwol(i) = a0 * exp(b0*t) * (a1*exp(b1*t) + a2*exp(b2*t) + a3)
  end do

  write(6,'(" alpha = ", f15.7)') alpha
  write(6,'(" c1 = ", f15.7,"   exact = ",f15.7)') c1, 216.0d0 / 78.0d0
  write(6,'(" l1 = ", f15.7,"   exact = ",f15.7)') l1, -(5.0+2.0d0*sqrt(3.0d0))+4.0d0*sqrt(3.0d0) 
  write(6,'(" l2 = ", f15.7,"   exact = ",f15.7)') l2, -(5.0+2.0d0*sqrt(3.0d0))

  open(10,file="result")
    do i = 0, nt
      t = i * dt
      if (mod(i,100) == 0) then
        write(10,'(4f20.10)') t, xana(i), xnum(i), xwol(i) 
      end if
    end do
  close(10)

!
end program main
!-----------------------------------------------------------------------
