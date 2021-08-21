!------------------------------------------------------------------------------
program main
!------------------------------------------------------------------------------
  implicit none

  integer :: n = 500
  integer :: ig, im
  real(8) :: x, dx

  dx = 6.0d0 / dble(n)

  do ig = 0, 500
    x = -3.0d0 + dx * ig
    write(6,'(5000f15.7)') x, fpsum(1,x), fpsum(2,x), fpsum(10,x), &
                           fpsum(20,x), fpsum(1000,x) 

  end do

  contains
!------------------------------------------------------------------------------
    function fpsum(M, x)
!------------------------------------------------------------------------------
      implicit none

      real(8), parameter  :: pi = acos(-1.0d0)

      integer, intent(in) :: M
      real(8), intent(in) :: x

      real(8) :: fpsum

      integer :: i, j, k
      real(8) :: val, fact

      fpsum = 0.0d0

      !fact =  -1.0d0

      fpsum = 0.5d0
      do i = 1, M
        fact = 2.0d0 * i - 1.0d0 
        !val   = fact / (2.0d0 * pi * i) * sin(2.0d0*pi*i*x) 
        val   = - 4.0d0/(pi*pi) * 1.0d0/(fact*fact) * cos(fact * pi * x) 
        fpsum = fpsum + val 
      end do
      !fpsum = fpsum * 4.0d0/pi


    end function fpsum
!------------------------------------------------------------------------------

end program main
!------------------------------------------------------------------------------
