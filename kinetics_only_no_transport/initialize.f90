module initialize

use globals

save

real(8) :: x(xn), y(yn), t(tn), ytemp(xn,yn)
real(8) :: permeability(xn,yn), permx(xn,yn), permy(xn,yn), permLong((xn-2)*(yn-2))
real(8) :: permeability0(xn,yn)
real(8) :: rhoLong((xn-2)*(yn-2))
real(8) :: permxLong((xn-2)*(yn-2)), permyLong((xn-2)*(yn-2))
real(8) :: rho0(xn,yn)
real(8) :: bcx0(xn,2), bcy0(2,yn), bcxPsi(xn,2), bcyPsi(2,yn), ic0(xn,yn)
real(8) :: kMat(xn,yn), lambdaMat(xn,yn), porosity(xn,yn), permeable(xn)
real(8) :: sedx
real(8) :: wt1, wt2

contains

! ----------------------------------------------------------------------------------%%
!
! SUBROUTINE TO INITIALIZE, JUST CALL IT
!
! ----------------------------------------------------------------------------------%%
  
subroutine init ()
use globals
integer :: m,n


! SET UP THINGS THAT CAN'T BE DONE IN THE MODULE FOR WHATEVER REASON
dx = ( x_max - x_min ) / real ( xn - 1, kind = 8 ) 
x = linspace ( xn, x_min, x_max )
dy = ( y_max - y_min ) / real ( yn - 1, kind = 8 ) 
y = linspace ( yn, y_min, y_max )
dt = ( t_max - t_min ) / real ( tn - 1, kind = 8 ) 
t = linspace ( tn, t_min, t_max)
  
  
! BOUNDARY CONDITIONS
ic0(:,:) = 280.0 ! IC
do j=1,yn
do i=1,xn
	ic0(i,j) = 280.0 + 5.0*cos(3.14*x(i)/750)
end do
end do

bcx0(:,1) = 280.0 ! bottom
do i=1,xn
	bcx0(i,1) = 280.0 + 5.0*cos(3.14*x(i)/750)
end do
bcx0(:,2) = 280.0 ! top

bcy0(1,:) = 280.0 ! left
bcy0(2,:) = 280.0 ! right

bcyPsi(1,1:yn) = 0.0 ! left
bcyPsi(2,1:yn) = 0.0 ! right
bcxPsi(1:xn,1) = 0.0 ! bottom
bcxPsi(1:xn,2) = 0.0 ! top


! PERMEABILITY
lambdaMat = 2.6

permeability = 1e-21
do i=1,yn
do j=1,xn
	
! 	if (y(i) .ge. -800.0) then
! 		permeability(j,i) = 1e-13
! 	end if
!
!	if ((y(i) .ge. -200.0) .and. (x(j) .ge. 500.0) .and. (x(j) .le. (x_max-500.0))) then
!		permeability(j,i) = 4e-15
!	end if

	if (y(i) .ge. y_min) then
	permeability(j,i) = (0.5+0.5*tanh((y(i)+((800.0)))/50.0))*1e-13 &
	&+ (1.0 - (0.5+0.5*tanh((y(i)+((800.0)))/50.0)))*1e-21
	end if

	sedx = 400.0-600.0*( ( (x(j)/(x_max/2.0)) - 1.0) **2.0) ! must be the parabola
	sedx = 200.0

	if ((y(i) .gt. -500.0) .and. (x(j) .gt. 500.0) .and. (x(j) .lt. (x_max-500.0))) then
	permeability(j,i) = (0.5+0.5*tanh((y(i)+sedx)/50.0))*4e-15 &
	&+ (1.0 - (0.5+0.5*tanh((y(i)+sedx)/50.0)))*1e-13
	end if

	if ((y(i) .gt. -150.0) .and. (x(j) .le. x_max/2)) then
	permeability(j,i) = (0.5+0.5*tanh((x(j)-500.0)/50.0))*4e-15 &
	&+ (1.0 - (0.5+0.5*tanh((x(j)-500.0)/50.0)))*1e-13
	end if

	if ((y(i) .gt. -150.0) .and. (x(j) .gt. (x_max/2)-3.0)) then
	permeability(j,i) = (0.5+0.5*tanh((x(j)-2500.0)/50.0))*1e-13 &
	&+ (1.0 - (0.5+0.5*tanh((x(j)-2500.0)/50.0)))*4e-15
	end if
	
	if (y(i) .ge. -200.0) then
		lambdaMat(:,i) = 2.6
	end if
	
end do
end do



! ! first pass
! wt1 = 0.8
! wt2 = 0.05
! permeability = log(permeability)
! permeability0 = permeability
! do i=2,xn-1
! do j=2,yn-1
! 	permeability(i,j) = permeability0(i,j)*wt1 + permeability0(i-1,j)*wt2 + permeability0(i+1,j)*wt2 &
! 	& + permeability0(i,j-1)*wt2 + permeability0(i,j+1)*wt2
! end do
! end do
!
! ! second pass
! permeability0 = permeability
! do i=2,xn-1
! do j=2,yn-1
! 	permeability(i,j) = permeability0(i,j)*wt1 + permeability0(i-1,j)*wt2 + permeability0(i+1,j)*wt2 &
! 	& + permeability0(i,j-1)*wt2 + permeability0(i,j+1)*wt2
! end do
! end do
!
! ! third pass
! permeability0 = permeability
! do i=2,xn-1
! do j=2,yn-1
! 	permeability(i,j) = permeability0(i,j)*wt1 + permeability0(i-1,j)*wt2 + permeability0(i+1,j)*wt2 &
! 	& + permeability0(i,j-1)*wt2 + permeability0(i,j+1)*wt2
! end do
! end do
!
! permeability = exp(permeability)




! this should fix it, i just realized it might not.
permeability(:,yn) = permeability(:,yn-1)



! HEAT TRANSFER PARAMETERS
kMat = 2.6/(1000.0*4186.0)
ki=2.6/(1000.0*4186.0)


return

end subroutine init

end module initialize