use iso_c_binding
interface
  subroutine FloodDepressions(dem,width,height,no_data) bind(C, name="FloodDepressions")
    use iso_c_binding
    implicit none
    real (c_float), intent (out), dimension(*):: dem
    integer (c_int), VALUE:: width
    integer (c_int), VALUE:: height
    real (c_float), VALUE:: no_data
  end subroutine
end interface

REAL, DIMENSION(20, 20) :: A       ! A is a 5 x 10 x 5 array 

do j=1,20
  do i=1,20
    A(i,j) = j*20+i
  enddo
enddo

call random_number (A)

A(4,4)=0.1

call FloodDepressions(A,20,20,0.1)

end
