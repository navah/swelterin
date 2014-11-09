
module globals
save

integer, parameter :: testInt = 31, xn =60, yn = 60, altnum = 167, mstep = 10000, cell = 6
integer, parameter :: tn = 10000000, g_pri = 5, g_sec = 69, g_sol = 15, g_med = 7
real(8) :: x_min = 0.0D+00, x_max = 3000.0D+00, y_min = -1300.0D+00, y_max = 0.0D+00
real(8) :: t_min = 0.0D+00, t_max = 10000000000000.0 
real(8) :: dt, dx, dy, dt0 = 0.001
real(8) :: dPsi, psiLast(xn,yn)

real(8) :: ki = .56, ra = 100.0, viscosity = 4.3e-4, cp = 4186.0, alpha =4.0e-4, k
real(8) :: lambda = 2.6
real(8) :: g = 9.8
real(8) :: rho_fluid = 1000.0
real(8) :: one_array(xn,yn) = 1.0
real(8) :: pec = 200.0
real(8) :: dTdx(xn,yn), dTdy(xn,yn), dPdx(xn,yn), dPdy(xn,yn)



contains
  
  
! ----------------------------------------------------------------------------------%%
!
! CHECK SOMETHING
!
! ----------------------------------------------------------------------------------%%

subroutine check(status)
integer, intent ( in) :: status

	!if(status /= nf90_noerr) then 
	  !print *, trim(nf90_strerror(status))
	  !stop "Stopped"
	!end if
	
end subroutine check  
  
  
! ----------------------------------------------------------------------------------%%
!
! BAND
!
! SUMMARY: Gaussian elimination for the banded case, when you know the structure of
! 		   the band, which you usually do
!
! INPUTS: a(n,m) : original matrix, but just the band; n rows deep, m columns wide
!		  		   relation between matrix A and band a is A(i,j) = a(i,j-i+(m+1)/2)
!		  m : width of the band (MUST BE ODD! no 4-wide bands. makes sense)
!		  n : order of the matrix A or # of rows in a
!         
!
! RETURNS: band(n,m) : which contains upper triangular and lower triangular in the
!       			   columns that are not the central column
!
! ----------------------------------------------------------------------------------%%

  
function band(a,m,n)

implicit none
integer g,h,i,j,k,m,n,r
real(8) :: a(n,m), band(n,m)
real eps


band = a
r = (m+1)/2
eps = 1.0e-6

do 20 k = 1,n
if (abs(band(k,r)) .le. eps) goto 99
band(k,r) = 1.0/band(k,r)
	h = r-1
	i = k+1
10 	if (h .lt. 1 .or. i .gt. n) goto 20
	band(i,h) = band(i,h) * band(k,r)
		j = h+1
		g = r+1
30		if (g .gt. m .or. j .gt. (r+n-i)) goto 40
		band(i,j) = band(i,j) - band(i,h)*band(k,g)
		j = j+1
		g = g+1
		goto 30
		
40	continue
	i = i+1
	h = h-1
	goto 10
20 continue
return

99 m=0
return
end function band
  
  
function solve(a,b,m,n)

implicit none
integer i,j,k,m,n,r
real(8) :: a(n,m), b(n), solve(n)

solve = b
r = (m+1)/2

do 100 k=1,n-1
	i = k+1
	j = r-1
110	if (j .lt. 1 .or. i .gt. n) goto 100
	solve(i) = solve(i) - a(i,j)*solve(k)
	i = i+1
	j = j-1
	goto 110
100 continue

do 120 k=n,1,-1
	i = k+1
	j = r+1
130	if (j .gt. m .or. i .gt. n) goto 140
	solve(k) = solve(k) - a(k,j)*solve(i)
	i = i+1
	j=j+1
	goto 130
	
140 continue
solve(k) = solve(k)*a(k,r)

120 continue
return 
end function

! ----------------------------------------------------------------------------------%%
!
! LINSPACE
! 
! SUMMARY: Makes a vector of dimension n, starting with a_first, ending with a_last
!          with n-2 evenly spaced items along the way.
!
! INPUTS: n : number of elements
!         a_first : value of first element
!         a_last : value of last element
!
! RETURNS: linspace(a)
!
! ----------------------------------------------------------------------------------%%

function linspace ( n, a_first, a_last )

  implicit none
  integer, intent(in) :: n 
   real (8) :: linspace(n)
   integer :: i
   real (8) :: a(n)
   real (8) , intent(in) ::  a_first, a_last

  if ( n == 1 ) then
    a(1) = ( a_first + a_last ) / 2.0D+00
  else
    do i = 1, n
      a(i) = ( real ( n - i,     kind = 8 ) * a_first &
             + real (     i - 1, kind = 8 ) * a_last ) &
             / real ( n     - 1, kind = 8 )
    end do
  end if
  linspace = a
  return
  
end function linspace







! ----------------------------------------------------------------------------------%%
!
! FOUND GAUSSIAN ELIMINATION FUNCTION
!
! SUMMARY: Solve Ax = b style linear system using gaussian elimination
!
! INPUTS: a : A concatenated to b
!         row : n-by-n dimension
!
! RETURNS: Gsselm(row)
!
! ----------------------------------------------------------------------------------%%

function Gsselm(a,row)
	implicit none
	INTEGER, INTENT(IN) :: row
	REAL(8) , INTENT(IN OUT)  ::  a(:,:)   	!Assume shape (:)
	REAL(8) , DIMENSION(row) :: Gsselm

	INTEGER i,j,k
	INTEGER, DIMENSION(2) :: shap
	REAL(8) , ALLOCATABLE :: swap_ik(:)
	REAL(8)  :: tmp

	ALLOCATE (swap_ik(row+1))

! 	Initialise
	swap_ik(:) = 0.0d0            ! Whole vector initialized to zero
	tmp = 0.0d0

! Check dimensions of input matrix
	shap = SHAPE(a)
	if ( (shap(1) .NE. row) .OR.  (shap(2) .NE. row+1) ) then
 !  	call break()
	end if


!/*   Gaussian Elimination - Row Reduction of matrix  */
	do k=1, row-1                             ! total of row-1 operations

!/*  Pivotal strategy - SWAP rows to make pivotal element a[k][k] have the
!    greatest magnitude in its column. This prevents unnecessary division by
!    a small number.                           */
!   	do i = k+1, row
!      	if ( (dabs(a(i,k))-dabs(a(k,k))).gt.eps  ) then
!         	do j = k, row+1                     !/* If pivotal element is not */
!            	swap_ik(j) = a(k,j)              !/* the highest then  */
!          	   a(k,j) = a(i,j)                  !/* swap i'th and k'th rows */
!            	a(i,j) = swap_ik(j)
!         	end do 		!j-loop
!         end if
!   	end do 				!i-loop


!/*   If the Matrix is SINGULAR then EXIT program      */
!  		IF ( dabs(a(k,k)) < EPS ) then
!   		print *,'After swapping rows to make the pivotal element become the'
!      	print *,'highest magnitude element in its column, its magnitude is'
!	      print *,'still extremely small.'
!	      print *,'Hence this is a SINGULAR MATRIX - no unique solution or '
!	      print *,'check that the input dimensions are correct.'
!   	   call break()
!	   END if
   


!/*      Perform row-reduction with pivotal element a[k][k]     */
		do i = k+1, row
			do j = row+1, k, -1			!/* starting from end of column */
	     	 	a(i,j) = a(i,j) - a(k,j) / a(k,k) * a(i,k)
			end DO 							!/* end of j loop     */
		end do 	 							!/* end of 2nd i loop */

	end DO 									!/* end of k loop     */
!  At this point, the bottom triangle is Zero


!/*   Back Substitution - Solutions of equations   */
	Gsselm(row) = a(row,row+1) / a(row,row)
	do k = row-1, 1, -1
   	tmp = 0.0d0
   	do j = k+1, row
      	tmp = tmp + a(k,j)*Gsselm(j)
	   end do 							!j-loop
   	Gsselm(k) = ( a(k,row+1) - tmp ) / a(k,k)
	end do 								!k-loop

   deallocate (swap_ik)
	RETURN
END function Gsselm



! ----------------------------------------------------------------------------------%%
!
! TRIDIAG
!
! SUMMARY: Solves Ax = b style linear system where A is tridiagonal matrix
!
! INPUTS: aa : tridiagonal matrix
!         nn : nn-by-nn dimension
!
! RETURNS: tridiag(nn)
!
! ----------------------------------------------------------------------------------%%





function tridiag(a,b,c,d,nn)
integer :: nn, k, km1, i
real(8) :: aa(nn,nn), a(nn), b(nn), c(nn), d(nn), tridiag(nn), xm
!      dimension a(nn),b(nn),c(nn),d(nn)
      !tridiag = aa(:,nn+1)
      tridiag = d
      a(1) = 0.0
      !b(1) = aa(1,1)
      !c(1) = aa(1,2)
      !a(nn) = aa(nn,nn-1)
      !b(nn) = aa(nn,nn)
      c(nn) = 0.0
      !do k = 2,nn-1
      !	a(k) = aa(k,k-1)
      !	b(k) = aa(k,k)
      !	c(k) = aa(k,k+1)
      !end do

      
!      if(nn .eq. 1) then
!      tridiag(1)=d(1)/b(1)
!      write(*,*) "shit"
!      end if
      
      do k = 2,nn
      km1 = k - 1
      
      if(b(k-1) .eq. 0.0) then
      write(*,*) "what"
      end if
      
      xm  = a(k)/b(km1)
      b(k)  = b(k) - xm*c(km1)
      tridiag(k)  = tridiag(k) - xm*tridiag(km1)
      end do
      
      tridiag(nn)   = tridiag(nn)/b(nn)

      k = nn
      do i = 2,nn
      k = nn + 1 - i
      tridiag(k) = (tridiag(k) - c(k)*tridiag(k+1))/b(k)
      end do
      
!      format(/3x,'diagonal element .eq. 0 in tridag at k = ',i2/)
      return
end function tridiag











! ----------------------------------------------------------------------------------%%
!
! TIMESTAMP
!
! SUMMARY: Write a timestamp
!
! RETURNS: timestamp
!
! ----------------------------------------------------------------------------------%%

subroutine timestamp ( )

  implicit none
  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestamp


! ----------------------------------------------------------------------------------%%
!
! GIVE ME A NUMBER I AM NOT USING SO I CAN USE IT TO MAKE MYSELF SOME FILES
!
! ----------------------------------------------------------------------------------%%

function get_unit ( )

  implicit none
  integer :: i, ios, get_unit
  logical lopen
  get_unit = 0
  do i = 1, 99
    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then
      inquire ( unit = i, opened = lopen, iostat = ios )
      if ( ios == 0 ) then
        if ( .not. lopen ) then
          get_unit = i
          return
        end if
      end if
    end if
  end do
  return
end function get_unit



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Subroutine for smoothing of a time-series data .......            !!!
!!! Input to the subroutine is 1-d array and output is also the same  !!!
!!! Input:                                                            !!!
!!!             x       = one dimensional array; x(n)                 !!!
!!!             n       = no. of elements in input array              !!!
!!!             m       = no. of points smoothing (should be odd no.) !!!
!!! Output:                                                           !!!
!!!             y       = one dimensional array; y(n)                 !!!
!!! Version: March 24, 2010 (ATJ/DBS)                                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function moving_average (x1, n, m)

implicit none
real(8) :: x1(n), moving_average(n), sumx
integer :: m, i, j, i1, i2, k1, k2, iflag, nd, n

!        Dimension x1(n), moving_average(n)

	If (mod(m,2) == 0) m = m + 1

        Do i = 1, n

        k2 = (m/2)

 151    Continue

        iflag = 0

        i1 = i - k2 
        i2 = i + k2

        If ((i1 <= 0).or.(i2 > n)) iflag = 1 
        If ((i1 <= 0).or.(i2 > n)) k2 = k2 -1

        If (iflag == 1) Go to 151

        sumx = 0.0

        Do j = i1, i2
        sumx = sumx + x1(j)
        Enddo

        nd = (i2 - i1) + 1

        moving_average(i) = sumx/float(nd)

        Enddo

        Return

End function moving_average

end module globals
