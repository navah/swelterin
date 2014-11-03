! littleChunk.f90 
! gfortran -O3 -g -I/usr/local/include -c alteration.f90

! gfortran -O3 -g -I/usr/local/include -c littleChunk.f90
! gfortran -O3 -g -I/usr/local/include -o littleChunk littleChunk.o alteration.o -L/usr/local/lib -liphreeqc
! ./littleChunk

PROGRAM main
	use alteration
	
	implicit none
	
	! inputs
	real(8) :: temp, timestep, primary(5), secondary(69), solute(15), medium(7)
	
	! other stuff
	integer :: i, j
	real(8) ::  alt0(1,167) 
	
	! initial conditions
	timestep = 20000000000.0
	temp = 5.0
	
	primary(1) = 12.96 ! feldspar
	primary(2) = 6.96 ! augite
	primary(3) = 1.26 ! pigeonite
	primary(4) = .4 ! magnetite
	primary(5) = 96.77 ! basaltic glass
	
	secondary = 0.0
	
	! ocean today (w. sources)
	solute(1) = 7.8 ! ph
	solute(2) = .0023 ! Alk 1.6e-3
	solute(3) = .03860 ! water mass
	solute(4) = .0023 !1.2e-2 ! H2CO3
	solute(5) = .0105 ! Ca
	solute(6) = .0533 ! Mg
	solute(7) = .468 ! Na
	solute(8) = .00997 ! K
	solute(9) = 0.0 !1.2e-6 ! Fe
	solute(10) = .0281 ! 1.0e-4 ! S(6)
	solute(11) = 0.0 !2.0e-4 ! Si
	solute(12) = .546 ! Cl
	solute(13) = 1.0e-10 ! Al
	solute(14) = .00234 ! HCO3-
	solute(15) = 0.0 ! CO3-2
	
	! properties of the medium
	medium(1) = 0.0 ! phi 
	medium(2) = 0.0 ! s_sp
	medium(3) = .3860 ! water_volume
	!medium(:,1:(xn/cell)*(5/13),3) = .03 ! water_volume
	medium(4) = 0.0 ! rho_s
	medium(5) = 1.0 ! rxn toggle
	!medium(:,yn/cell,5) = 0.0
	medium(6) = 0.0 ! x-coord
	medium(7) = 0.0 ! y-coord
	
	write(*,*) "doing something..."
	
	do i=1,500
		write(*,*) i
		alt0 = alter(temp,timestep,primary,secondary,solute,medium)
		write(*,*) alt0
		!PARSING
! 		solute = (/ alt0(1,2), alt0(1,3), alt0(1,4), alt0(1,5), alt0(1,6), &
! 		alt0(1,7), alt0(1,8), alt0(1,9), alt0(1,10), alt0(1,11), alt0(1,12), &
! 		alt0(1,13), alt0(1,14), alt0(1,15), 0.0D+00/)
!
!
!
! 		secondary = (/ alt0(1,16), alt0(1,18), alt0(1,20), &
! 		alt0(1,22), alt0(1,24), alt0(1,26), alt0(1,28), alt0(1,30), alt0(1,32), alt0(1,34), &
! 		alt0(1,36), alt0(1,38), alt0(1,40), alt0(1,42), alt0(1,44), alt0(1,46), alt0(1,48), &
! 		alt0(1,50), alt0(1,52), alt0(1,54), alt0(1,56), alt0(1,58), alt0(1,60), alt0(1,62), &
! 		alt0(1,64), alt0(1,66), alt0(1,68), alt0(1,70)/)
!
!
! 		primary = (/ alt0(1,72), alt0(1,74), alt0(1,76), alt0(1,78), alt0(1,80)/)
		
!		write(*,*) "calcite"
! 		write(*,*) secondary(16)
 		write(*,*) "pri 1"
 		write(*,*) primary(5)
	
	end do
	
END PROGRAM main