
! ----------------------------------------------------------------------------------%%
!
! MASSACR 
! 
! SUMMARY: main method runs fluid dynamic simulation coupled to geochemical
!          simulation and writes selected model output to file
! 
! TO RUN: make -f theMakeFile
!		  mpirun -np 4 ./massacr
!
! ----------------------------------------------------------------------------------%%

PROGRAM main

use globals
use initialize
use alteration
!use netcdf

implicit none

include 'mpif.h'

! functions within massacr.f90
interface
	
	! solves thermal energy equation
	function h_next (h, psi, rho_in, flux)
		use globals
		use initialize
		! integers
		integer :: i, j, n, ii, c1, c2, c3, c4, m=5
		! inputs 
		real(8) :: sx, sy, qx, qy, rho_in(xn,yn), flux(xn,2)
		! velocity stuff
		real(8) :: u(xn,yn), v(xn,yn), uLong((xn-2)*(yn-2)), vLong((xn-2)*(yn-2))
		real(8) ::  velocities0(xn,2*yn)
		! matrix stuff
		real(8) :: h(xn,yn), h_next(xn,yn), psi(xn,yn)
		real(8) :: aBand((xn-2)*(yn-2),5), bBand((xn-2)*(yn-2),5)
		real(8) :: h0(xn,yn), uVec((xn-2)*(yn-2)), h_nextRow((xn-2)*(yn-2))
		real(8) :: kMatLong((xn-2)*(yn-2))
	end function h_next

	! solves streamfunction vorticity equation
	function psi_next (h, rhs0, psi,top_in,rho_in)
		use globals
		use initialize
		! integers
		integer :: i, j, n, m=5
		! inputs
		real(8) :: rhs0(xn,yn), rhs1(xn,yn), rhsLong((xn-2)*(yn-2))
		real(8) :: h(xn,yn), psi(xn,yn), rho_in(xn,yn)
		! matrix stuff
		real(8) :: uVec((xn-2)*(yn-2)), psiLong((xn)*(yn)), psi_nextRow((xn-2)*(yn-2))
		real(8) :: psi_next(xn,yn)
		real(8) :: top_in(xn,1)
		! tridiag
		real(8) :: aa((xn-2)*(yn-2),(xn-2)*(yn-2)), a((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
		real(8) :: bb((xn-2)*(yn-2),(xn-2)*(yn-2)), b((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
	end function psi_next
	
	! transports solutes
	function solute_next(sol, uTransport, vTransport)
		use globals
		use initialize
		! integers
		integer :: i, j, ii, n, m
		! inputs
		real(8) :: sol(xn/cell,yn/cell), sol0(xn/cell,yn/cell)
		real(8) :: uTransport(xn/cell,yn/cell), vTransport(xn/cell,yn/cell)
		! solver stuff
		real(8) :: uLong((xn/cell-2)*(yn/cell-2)), vLong((xn/cell-2)*(yn/cell-2))
		real(8) :: aBand((xn/cell-2)*(yn/cell-2),5), bBand((xn/cell-2)*(yn/cell-2),5)
		real(8) :: qx, qy, solute_next(xn/cell,yn/cell), vec((xn/cell-2)*(yn/cell-2))
		real(8) :: sol_nextRow((xn/cell-2)*(yn/cell-2))
	end function solute_next
	
	! runs geochemical alteration model for a single (coarse) cell
	function alt_next (temp, timestep, primaryList, secondaryList, soluteList, mediumList)
		use globals
		use initialize
		use alteration
		! declare yo shit
		integer :: order
		real(8) :: temp, timestep
		real(8) :: alt_next(1,167)
		real(8) :: alter0(1,167)
		real(8) :: primaryList(g_pri), secondaryList(g_sec), soluteList(g_sol), mediumList(g_med)
	end function alt_next

	! calculates fluid density
	function rho_next (h_in)
		use globals
		use initialize
		integer :: i,j
		real(8) :: h_in(xn,yn), rho_next(xn,yn)
	end function rho_next

	! calculates velocities from streamfunction values
	function velocities(psi)
		use globals
		use initialize
		implicit none
		real(8) :: velocities(xn,2*yn), psi(xn,yn)
		real(8) :: u0(xn,yn), v0(xn,yn)
	end function velocities

        ! calculates velocities from COARSE streamfunction values
	function velocitiesCoarse(psiCoarse)
		use globals
		use initialize
		implicit none
		real(8) :: velocitiesCoarse(xn/cell,2*yn/cell), psiCoarse(xn/cell,yn/cell)
		real(8) :: u0(xn/cell,yn/cell), v0(xn/cell,yn/cell)
	end function velocitiesCoarse

	! calculates partial derivative of any 1D or 2D array
	function partial(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(8) :: array(rows,cols), d1, d2, d
		real(8) :: partial(rows,cols)
	end function partial

	! writes 2D array to file
	function write_matrix ( m, n, table, filename )
		use globals
		implicit none
		integer :: m, n, j, output_status, unit0
		character ( len = * ) filename
		character ( len = 30 ) string
		real(4)  :: table(m,n) , write_matrix
	end function write_matrix

	! writes 1D array to file
	function write_vec ( n, vector, filename )
		use globals
		implicit none
		integer :: n, j, output_status, unit0
		character ( len = * ) filename 
		real(4)  :: vector(n), write_vec
	end function write_vec

end interface

!--------------DECLARE EVERYTHING 

! dependent variable arrays
real(8) :: h(xn,yn), psi(xn,yn) ! xn rows deep & yn columns wide
real(8) :: hmat(xn,(yn*tn/mstep)), psimat(xn,(yn*tn/mstep))
real(8) :: velocities0(xn,2*yn)
real(8) :: umat(xn,(yn*tn/mstep)), vmat(xn,(yn*tn/mstep))
real(8) :: u(xn,yn), uLong(xn*yn), v(xn,yn), vLong(xn*yn)

! material properties
real(8) :: rho(xn,yn), flux(xn,2)
real(8) :: rhs0(xn,yn)
integer :: unit

! netCDF & output stuff
integer :: xInt, yInt, tInt, hInt, uInt, vInt
integer :: ncid
integer :: x_dimid, y_dimid, t_dimid, h_dimid, u_dimid, v_dimid
integer :: x_varid, y_varid, t_varid, h_varid, u_varid, v_varid
integer :: i, j, ii, m, n, jj
real(8) :: yep

! benchmark stuff
real(8) :: nusseltLocalv(xn,1), nuBar

! geochemical alteration stuff
real(8) :: alt0(1,altnum)
real(8) :: primary(xn/cell,yn/cell,g_pri), primaryMat(xn/cell,yn*tn/(cell*mstep),g_pri)
real(8) :: secondary(xn/cell,yn/cell,g_sec), secondaryMat(xn/cell,yn*tn/(cell*mstep),g_sec)
real(8) :: solute(xn/cell,yn/cell,g_sol), soluteMat(xn/cell,yn*tn/(cell*mstep),g_sol)
real(8) :: medium(xn/cell,yn/cell,g_med), mediumMat(xn/cell,yn*tn/(cell*mstep),g_med)

! solute transport stuff
real(8) :: uTransport(xn/cell,yn/cell), vTransport(xn/cell,yn/cell)
real(8) :: uCoarse(xn/cell,yn/cell), vCoarse(xn/cell,yn/cell)
real(8) :: solTemp(xn/cell,yn/cell), soluteOcean(g_sol)
real(8) :: psiCoarse(xn/cell,yn/cell), velocitiesCoarse0(xn/cell,2*yn/cell)

! message passing stuff
integer, parameter :: max_rows = 10000000
integer, parameter :: send_data_tag = 2001, return_data_tag = 2002
integer :: my_id, root_process, ierr, status(MPI_STATUS_SIZE)
integer :: num_procs, an_id, num_rows_to_receive
integer :: avg_rows_per_process, num_rows, num_rows_to_send
integer :: end_row, sender, start_row, num_rows_received
real(8) :: vector(max_rows), vector2(max_rows), partial_sum, sum
real(8) :: local_mean, global_mean
real(8) :: hLocal((xn/cell)*(yn/cell)), dt_local
integer :: order

! formatted message passing arrays
real(8) :: hLong((xn/cell)*(yn/cell))
real(8) :: priLong((xn/cell)*(yn/cell),g_pri), priLocal((xn/cell)*(yn/cell),g_pri)
real(8) :: secLong((xn/cell)*(yn/cell),g_sec), secLocal((xn/cell)*(yn/cell),g_sec)
real(8) :: solLong((xn/cell)*(yn/cell),g_sol), solLocal((xn/cell)*(yn/cell),g_sol)
real(8) :: medLong((xn/cell)*(yn/cell),g_med), medLocal((xn/cell)*(yn/cell),g_med)
real(8) :: priLongBit((xn/cell)*(yn/cell)), priLocalBit((xn/cell)*(yn/cell))
real(8) :: secLongBit((xn/cell)*(yn/cell)), secLocalBit((xn/cell)*(yn/cell))
real(8) :: solLongBit((xn/cell)*(yn/cell)), solLocalBit((xn/cell)*(yn/cell))
real(8) :: medLongBit((xn/cell)*(yn/cell)), medLocalBit((xn/cell)*(yn/cell))

!--------------GEOCHEMICAL INITIAL CONDITIONS

! primary minerals [mol]
primary(:,:,:) = 0.0
primary(:,:,1) = 1.29600 ! feldspar
primary(:,:,2) = .69600 ! augite
primary(:,:,3) = .12600 ! pigeonite
primary(:,:,4) = .04000 ! magnetite
primary(:,:,5) = 9.67700 ! basaltic glass

! secondary minerals [mol]
secondary(:,:,:) = 0.0

! !columbia river composition
! ! hydrothermal solute concentrations [mol/kgw]
! solute(:,:,1) = 7.8 ! ph
! solute(:,:,2) = 8.451 ! pe
! solute(:,:,3) = 2.3e-3 ! Alk 1.6e-3
! solute(:,:,4) = 2.200e-3 !1.2e-2 ! H2CO3
! solute(:,:,5) = 6.0e-3 ! Ca
! !solute(:,:,5) = 0.0e-3 ! Ca
! solute(:,:,6) = 2.0e-5 ! Mg
! solute(:,:,7) = 1.0e-3 ! Na
! solute(:,:,8) = 1.0e-4 ! K
! solute(:,:,9) = 1.2e-6 ! Fe
! solute(:,:,10) = 0.0 !1.0e-4 ! S(6)
! solute(:,:,11) = 2.0e-4 ! Si
! solute(:,:,12) = 3.0e-4 ! Cl
! solute(:,:,13) = 1.0e-6 ! Al
! solute(:,:,14) = 2.200e-3 ! HCO3-
! solute(:,:,15) = 0.0 ! CO3-2

! ! hydrothermal solute concentrations [mol/kgw]
! solute(:,:,1) = 7.8 ! ph
! solute(:,:,2) = 8.451 ! pe
! solute(:,:,3) = 2.3e-3 ! Alk 1.6e-3
! solute(:,:,4) = 2.200e-3 !1.2e-2 ! H2CO3
! solute(:,:,5) = .0103 ! Ca
! solute(:,:,6) = .0528 ! Mg
! solute(:,:,7) = .469 ! Na
! solute(:,:,8) = .0102 ! K
! solute(:,:,9) = 0.0 ! Fe
! solute(:,:,10) = 1e-6 ! 1.0e-4 ! S(6)
! solute(:,:,11) = 0.0 ! Si
! solute(:,:,12) = .546 ! Cl
! solute(:,:,13) = 0.0 ! Al
! solute(:,:,14) = 2.200e-3 ! HCO3-
! solute(:,:,15) = 0.0 ! CO3-2

! ocean today (w. sources)
solute(:,:,1) = 7.8 ! ph
solute(:,:,2) = .0023 ! Alk 1.6e-3
solute(:,:,3) = .03860 ! water mass
solute(:,:,4) = .0023 !1.2e-2 ! H2CO3
solute(:,:,5) = .0105 ! Ca
solute(:,:,6) = .0533 ! Mg
solute(:,:,7) = .468 ! Na
solute(:,:,8) = .00997 ! K
solute(:,:,9) = 0.0 !1.2e-6 ! Fe
solute(:,:,10) = .0281 ! 1.0e-4 ! S(6)
solute(:,:,11) = 0.0 !2.0e-4 ! Si
solute(:,:,12) = .546 ! Cl
solute(:,:,13) = 1.0e-10 ! Al
solute(:,:,14) = .00234 ! HCO3-
solute(:,:,15) = 0.0 ! CO3-2

! seawater solute concentrations [mol/kgw]
soluteOcean = (/ solute(1,1,1), solute(1,1,2), solute(1,1,3), solute(1,1,4), solute(1,1,5), & 
			  & solute(1,1,6), solute(1,1,7), solute(1,1,8), solute(1,1,9), solute(1,1,10), &
			  & solute(1,1,11), solute(1,1,12), solute(1,1,13), solute(1,1,14), solute(1,1,15) /)

! properties of the medium
medium(:,:,1) = 0.0 ! phi 
medium(:,:,2) = 0.0 ! s_sp
medium(:,:,3) = .03860 ! water_volume
!medium(:,1:(xn/cell)*(5/13),3) = .03 ! water_volume
medium(:,:,4) = 0.0 ! rho_s
medium(:,:,5) = 1.0 ! rxn toggle
medium(:,yn/cell,5) = 0.0
medium(:,:,6) = 0.0 ! x-coord
medium(:,:,7) = 0.0 ! y-coord

write(*,*) "testing..."

!--------------INITIALIZE ALL PROCESSORS

! process #0 is the root process
root_process = 0

! initialize a process
call MPI_INIT ( ierr )

! find out the process ID and how many processes were started so far
call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

! print out current processor id
write(*,*) "my_id:", my_id
write(*,*) " "

! what to do if you are the master processor
if (my_id .eq. root_process) then 

!--------------DO STUFF WITH THE MASTER PROCESSOR

! initialize domain geometry
call init()

! fill coordinate arrays
do i = 1,(xn/cell)
	do ii = 1,(yn/cell)
		medium(i,ii,6) = x(i*cell) 
		medium(i,ii,7) = y(ii*cell)
	end do
end do

! boundary & initial condtions for flow equations
psi=0.0
psi(1,1:yn) = bcyPsi(1,1:yn)
psi(xn,1:yn) = bcyPsi(2,1:yn)
psi(1:xn,1) = bcxPsi(1:xn,1)
psi(1:xn,yn) = bcxPsi(1:xn,2)

h(1:xn,1:yn) = ic0
h(1,1:yn) = bcy0(1,1:yn)
h(xn,1:yn) = bcy0(2,1:yn)
h(1:xn,1) = bcx0(1:xn,1)
h(1:xn,yn) = bcx0(1:xn,2)

! put initial values into array for file
hmat(1:xn,1:yn) = h
psimat(1:xn,1:yn) = psi
umat(1:xn,1:yn) = u
vmat(1:xn,1:yn) = v

uTransport = 0.0
vTransport = 0.0
uCoarse = 0.0
vCoarse = 0.0

! read steady state results from file
write(*,*) "reading..."
OPEN(UNIT=11, FILE="steady_h.txt")
DO i = 1,xn
	READ(11,*) (h(i,ii),ii=1,yn)
END DO
h = transpose(h)

write(*,*) "reading..."
OPEN(UNIT=12, FILE="steady_u.txt")
DO i = 1,xn
	READ(12,*) (u(i,ii),ii=1,yn)
END DO
u = transpose(u)

write(*,*) "reading..."
OPEN(UNIT=13, FILE="steady_v.txt")
DO i = 1,xn
	READ(13,*) (v(i,ii),ii=1,yn)
END DO
v = transpose(v)

write(*,*) "reading..."
OPEN(UNIT=14, FILE="steady_psi.txt")
DO i = 1,xn
	READ(14,*) (psi(i,ii),ii=1,yn)
END DO
psi = transpose(psi)

! this is the main loop that does all the solving for tn timesteps
do j = 2, tn
	
	! print current timestep
	if (mod(j,mstep/10) .eq. 0) then
		write(*,*) j
	end if

	! THERMAL BOUNDARY CONDITIONS
	
	! bottom
	do i = 1,xn
		flux(i,1) = h(i,2) +((.1))*dy/(2.6)
	end do
	
	! top
	flux(:,2) = 280.0
	do i = 1,xn/2
		flux(i,2) = 274.0
	end do
	
! 	! stop fluid dynamic model after 19k timesteps
! 	if (j .lt. 19000) then
! 		! solve thermal energy equation
! 		rho = rho_next(h)
! 		h = h_next(h, psi,rho, flux)
!
! 		! put thermal energy boundary conditions in for next timestep
! 		h(1,:) = (4.0/3.0)*h(2,:) - (1.0/3.0)*h(3,:) ! left
! 		h(xn,:) = (4.0/3.0)*h(xn-1,:) - (1.0/3.0)*h(xn-2,:) ! right
! 		h(:,1) = flux(:,1)
! 		h(:,yn) = flux(:,2)
!
! 		! solve streamfunction-vorticity equation
! 		rhs0 = (1.0/(viscosity))*g*rho_fluid*alpha*partial(h,xn,yn,dx,dy,1)
! 		psi = psi_next(h, rhs0, psi, permeable, rho)
!
! 		! put in streamfunction boundary conditions for next timestep
! 		psi(1,1:yn) = bcyPsi(1,1:yn) ! left
! 		psi(xn,1:yn) = bcyPsi(2,1:yn) ! right
! 		psi(1:xn,1) = bcxPsi(1:xn,1) ! bottom
! 		psi(1:xn,yn) = bcxPsi(1:xn,2) ! top
! 		psi(:,yn) = ((4.0/3.0)*psi(:,yn-1) - (1.0/3.0)*psi(:,yn-2))/1.0
! 		permeable = psi(:,yn)
!
! 		! get velocities from streamfunction
! 		velocities0 = velocities(psi)
! 		u = velocities0(1:xn,1:yn)
! 		v = velocities0(1:xn,yn+1:2*yn)
! 	end if
!
! 	!write sample steady-state simulation to file
! 	!only gotta do this once
! 	if (j .eq. 19000) then
! 		 yep = write_matrix(xn, yn, real(h,kind=4), 'steady_h.txt')
! 		 yep = write_matrix(xn, yn, real(u,kind=4), 'steady_u.txt')
! 		 yep = write_matrix(xn, yn, real(v,kind=4), 'steady_v.txt')
! 		 yep = write_matrix(xn, yn, real(psi,kind=4), 'steady_psi.txt')
! 		 yep = write_matrix(xn, yn, real(permeability,kind=4), 'steady_permeability.txt')
! 	end if
	
	! interpolate fine grid onto coarse grid
	do i = 1,xn/cell
	do ii = 1,yn/cell
        psiCoarse(i,ii) = psi(i*cell,ii*cell)
	end do
	end do

    velocitiesCoarse0 = velocitiesCoarse(psiCoarse)
    uCoarse = velocitiesCoarse0(1:xn/cell,1:yn/cell)
    vCoarse = velocitiesCoarse0(1:xn/cell,yn/cell+1:2*yn/cell)

	uTransport = uTransport + uCoarse
	vTransport = vTransport + vCoarse
	
	! things only done every mth timestep go here
	if (mod(j,mstep) .eq. 0) then
	
		yep = write_matrix( 2, 1, real((/j*1.0, tn/), kind = 4), 'upStep1.txt' )
	
		! make coarse grid average velocities
		uTransport = uTransport/mstep
		vTransport = vTransport/mstep
	
		
!!!!!!!!!! CHANGE TRANSPORT !!!!!!!!!!!!!!
			
! 			! convert pH, pe to concentrations
 			!do i=1,xn/cell
 			!	do ii=1,yn/cell
 			!		solute(i,ii,1) = 10**(-solute(i,ii,1))
 			!		!solute(i,ii,2) = 10**(-solute(i,ii,2))
 			!	end do
 			!end do

                        ! actual boundary conditions
		do n=1,g_sol
			!solute(:,yn/cell,n) = (soluteOcean(n)) ! top
			do i=1,(yn/cell)
				if (vTransport(i,yn/cell) .lt. 0.0) then
					solute(i,yn/cell,n) = (soluteOcean(n)) ! last
                                    else
                                       solute(i,yn/cell,n)=(4.0/3.0)*solute(i,yn/cell-1,n)-(1.0/3.0)*solute(i,yn/cell-2,n)

				end if
			end do
			!solute(:,1,n) = (soluteOcean(n)) ! last
			!solute(1,:,n) = (soluteOcean(n)) ! last
			!solute(yn/cell,:,n) = (soluteOcean(n)) ! last
			solute(:,1,n) = (4.0/3.0)*solute(:,2,n) - (1.0/3.0)*solute(:,3,n) ! bottom
			!solute(1,:,n) = (4.0/3.0)*solute(2,:,n) - &
			!					& (1.0/3.0)*solute(3,:,n)  ! left
			!solute(yn/cell,:,n) = (4.0/3.0)*solute(yn/cell-1,:,n) - &
			!					& (1.0/3.0)*solute(yn/cell-2,:,n)  ! right
		end do



			n=1 ! ph
	 		solTemp = solute(:,:,n)
	 		solute(:,:,n) = solute_next(solTemp,uTransport,vTransport)
! 			n=2 ! pe
! 	 		solTemp = solute(:,:,n)
! 	 		solute(:,:,n) = solute_next(solTemp,uTransport,vTransport)
			n=4 ! c
	 		solTemp = solute(:,:,n)
	 		solute(:,:,n) = solute_next(solTemp,uTransport,vTransport)
			n=5 ! ca
	 		solTemp = solute(:,:,n)
	 		solute(:,:,n) = solute_next(solTemp,uTransport,vTransport)
			n=6 ! mg
	 		solTemp = solute(:,:,n)
	 		solute(:,:,n) = solute_next(solTemp,uTransport,vTransport)
			n=7 ! na
	 		solTemp = solute(:,:,n)
	 		solute(:,:,n) = solute_next(solTemp,uTransport,vTransport)
			n=8 ! k
	 		solTemp = solute(:,:,n)
	 		solute(:,:,n) = solute_next(solTemp,uTransport,vTransport)
			n=9 ! fe
	 		solTemp = solute(:,:,n)
	 		solute(:,:,n) = solute_next(solTemp,uTransport,vTransport)
			n=10 ! s issues
	 		solTemp = solute(:,:,n)
	 		solute(:,:,n) = solute_next(solTemp,uTransport,vTransport)
			n=11 ! si
	 		solTemp = solute(:,:,n)
	 		solute(:,:,n) = solute_next(solTemp,uTransport,vTransport)
			n=12 ! cl
	 		solTemp = solute(:,:,n)
	 		solute(:,:,n) = solute_next(solTemp,uTransport,vTransport)
			n=13 ! al
	 		solTemp = solute(:,:,n)
	 		solute(:,:,n) = solute_next(solTemp,uTransport,vTransport)


! 			!! convert [H+], [e-] back to pH, pe
 			!do i=1,xn/cell
 			!	do ii=1,yn/cell
 			!		solute(i,ii,1) = -log10(solute(i,ii,1))
! 			!		!solute(i,ii,2) = -log10(solute(i,ii,2))
 			!	end do
 			!end do

!!!!!!!!!! CHANGE TRANSPORT !!!!!!!!!!!!!!

write(*,*) maxval(solute(:,:,4))

		!-TRANSPOSE 1
		! stretch everything out
		!hLong = reshape(h(1:xn-1:cell,1:yn-1:cell), (/(xn/cell)*(yn/cell)/)) ! for cell > 1
		hLong = reshape(h(1:xn:cell,1:yn:cell), (/(xn/cell)*(yn/cell)/)) ! for cell = 1
		priLong = reshape(primary, (/(xn/cell)*(yn/cell), g_pri/))
		secLong = reshape(secondary, (/(xn/cell)*(yn/cell), g_sec/))
		solLong = reshape(solute, (/(xn/cell)*(yn/cell), g_sol/))
		medLong = reshape(medium, (/(xn/cell)*(yn/cell), g_med/))
	
		!--------------MESSAGE DISTRIBUTING FROM MASTER TO SLAVES
		do an_id = 1, num_procs - 1
		
			! put number of rows in vector here for hLong
			num_rows = (xn/cell)*(yn/cell)
			avg_rows_per_process = num_rows / (num_procs-1)
	        start_row = ( (an_id-1) * avg_rows_per_process) + 1
	        end_row = start_row + avg_rows_per_process - 1
	        if (an_id .eq. (num_procs - 1)) end_row = num_rows
	        num_rows_to_send = (end_row - start_row + 1)
		
			! send size of temperature array chunk to processor an_id
	        call MPI_SEND( num_rows_to_send, 1, MPI_INTEGER, &
			an_id, send_data_tag, MPI_COMM_WORLD, ierr)
		
			! send timestep size to processor an_id
	        call MPI_SEND( dt, 1, MPI_DOUBLE_PRECISION, &
			an_id, send_data_tag, MPI_COMM_WORLD, ierr)
		
			! send temperature array chunk to processor an_id
	        call MPI_SEND( hLong(start_row), num_rows_to_send, MPI_DOUBLE_PRECISION, &
			an_id, send_data_tag, MPI_COMM_WORLD, ierr)
		
			! send primary array chunk to processor an_id
			do ii = 1,g_pri
				priLongBit = priLong(:,ii)
	        	call MPI_SEND( priLongBit(start_row), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, send_data_tag, MPI_COMM_WORLD, ierr)
			end do
		
			! send secondary array chunk to processor an_id
			do ii = 1,g_sec
				secLongBit = secLong(:,ii)
	        	call MPI_SEND( secLongBit(start_row), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, send_data_tag, MPI_COMM_WORLD, ierr)
			end do
		
			! send solute array chunk to processor an_id
			do ii = 1,g_sol
				solLongBit = solLong(:,ii)
	        	call MPI_SEND( solLongBit(start_row), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, send_data_tag, MPI_COMM_WORLD, ierr)
			end do
			
			! send medium array chunk to processor an_id
			do ii = 1,g_med
				medLongBit = medLong(:,ii)
	        	call MPI_SEND( medLongBit(start_row), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, send_data_tag, MPI_COMM_WORLD, ierr)
			end do
			
			write(*,*) "DONE SENDING TO PROCESSOR", an_id

		end do

		!--------------MESSAGE RECEIVING FROM SLAVE PROCESSORS
		do an_id = 1, num_procs - 1
		
			! get the size of each chunk again
			num_rows = (xn/cell)*(yn/cell)
			avg_rows_per_process = num_rows / (num_procs-1)
	        start_row = ( (an_id-1) * avg_rows_per_process) + 1
	        end_row = start_row + avg_rows_per_process - 1
	        if (an_id .eq. (num_procs - 1)) end_row = num_rows
	        num_rows_to_send = (end_row - start_row + 1)
		
			! receive primary chunk
			do ii = 1,g_pri
				! receive it
				call MPI_RECV( priLocal(:,ii), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
				! fill it
				priLong(start_row:end_row,ii) = priLocal(1:num_rows_to_send,ii)
			end do
		
			! receive secondary chunk
			do ii = 1,g_sec
				! receive it
				call MPI_RECV( secLocal(:,ii), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
				! fill it
				secLong(start_row:end_row,ii) = secLocal(1:num_rows_to_send,ii)
			end do
		
			! receive solute chunk
			do ii = 1,g_sol
				! receive it
				call MPI_RECV( solLocal(:,ii), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
				! fill it
				solLong(start_row:end_row,ii) = solLocal(1:num_rows_to_send,ii)
			end do
			
			! receive medium chunk
			do ii = 1,g_med
				! receive it
				call MPI_RECV( medLocal(:,ii), num_rows_to_send, MPI_DOUBLE_PRECISION, &
				an_id, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
				! fill it
				medLong(start_row:end_row,ii) = medLocal(1:num_rows_to_send,ii)
			end do

			write(*,*) "DONE RECEIVING FROM PROCESSOR", an_id
		
		end do
		

		! reset coarse grid velocities for next timestep
		uTransport = 0.0
		vTransport = 0.0
	
		!--------------MASTER PROCESSOR SAVES OUTPUT TO BE WRITTEN TO FILE
	
		! put stretched vectors back into 2d arrays
		primary = reshape(priLong,(/(xn/cell), (yn/cell), g_pri/))
		secondary = reshape(secLong,(/(xn/cell), (yn/cell), g_sec/))
		solute = reshape(solLong,(/(xn/cell), (yn/cell), g_sol/))
		medium = reshape(medLong,(/(xn/cell), (yn/cell), g_med/))
		!-TRANSPOSE 2

		! add timestep's output to output arrays
		 hmat(1:xn,1+yn*(j/mstep-1):1+yn*(j/mstep)) = h
		 psimat(1:xn,1+yn*(j/mstep-1):1+yn*(j/mstep)) = psi
		 umat(1:xn,1+yn*(j/mstep-1):1+yn*(j/mstep)) = u
		 vmat(1:xn,1+yn*(j/mstep-1):1+yn*(j/mstep)) = v
		 primaryMat(1:xn/cell,1+(yn/cell)*(j/mstep-1):1+(yn/cell)*(j/mstep),:) = primary
		 secondaryMat(1:xn/cell,1+(yn/cell)*(j/mstep-1):1+(yn/cell)*(j/mstep),:) = secondary
		 soluteMat(1:xn/cell,1+(yn/cell)*(j/mstep-1):1+(yn/cell)*(j/mstep),:) = solute
		 mediumMat(1:xn/cell,1+(yn/cell)*(j/mstep-1):1+(yn/cell)*(j/mstep),:) = medium
	 
	end if 
	! end mth timestep loop, finally

end do 
! end all timestep loop




!--------------WRITE EVERYTHING TO FILE

yep = write_vec ( xn, real(x,kind=4), 'x.txt' )
yep = write_vec ( yn, real(y,kind=4), 'y.txt' )
!yep = write_vec ( tn, real(t, kind=4), 't.txt' )
!yep = write_matrix ( xn, yn*tn/mstep, real(hmat, kind = 4), 'hMat.txt' )
!yep = write_matrix ( xn, yn*tn/mstep, real(psimat,kind=4), 'psiMat.txt' )
!yep = write_matrix ( xn, yn*tn/mstep, real(umat, kind = 4), 'uMat.txt' )
!yep = write_matrix ( xn, yn*tn/mstep, real(vmat,kind=4), 'vMat.txt' )
yep = write_matrix ( xn, yn, real(rho,kind=4), 'rho.txt' )
yep = write_matrix ( xn, yn,real(permeability,kind=4), 'permeability.txt' )

! secondary minerals
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,1),kind=4), 'sec_stilbite.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,2),kind=4), 'sec_sio2.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,3),kind=4), 'sec_kaolinite.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,4),kind=4), 'sec_albite.txt')
 yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,5),kind=4), 'sec_saponite.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,6),kind=4), 'sec_celadonite.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,7),kind=4), 'sec_clinoptilolite.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,8),kind=4), 'sec_pyrite.txt')
 yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,9),kind=4), 'sec_mont_na.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,10),kind=4), 'sec_goethite.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,11),kind=4), 'sec_dolomite.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,12),kind=4), 'sec_smectite.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,13),kind=4), 'sec_dawsonite.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,14),kind=4), 'sec_magnesite.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,15),kind=4), 'sec_siderite.txt')
 yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,16),kind=4), 'sec_calcite.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,17),kind=4), 'sec_quartz.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,18),kind=4), 'sec_kspar.txt')
 yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,19),kind=4), 'sec_saponite_na.txt')
 yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,20),kind=4), 'sec_nont_na.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,21),kind=4), 'sec_nont_mg.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,22),kind=4), 'sec_nont_k.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,23),kind=4), 'sec_nont_h.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,24),kind=4), 'sec_nont_ca.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,25),kind=4), 'sec_muscovite.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,26),kind=4), 'sec_mesolite.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,27),kind=4), 'sec_hematite.txt')
! yep = write_matrix(xn/cell, yn*tn/(cell*mstep), real(secondaryMat(:,:,28),kind=4), 'sec_diaspore.txt')

! solute concentrations
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,1),kind=4), 'sol_ph.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,3),kind=4), 'sol_w.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,2),kind=4), 'sol_alk.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,4),kind=4), 'sol_c.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,5),kind=4), 'sol_ca.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,6),kind=4), 'sol_mg.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,7),kind=4), 'sol_na.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,8),kind=4), 'sol_k.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,9),kind=4), 'sol_fe.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,10),kind=4), 'sol_s.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,11),kind=4), 'sol_si.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,12),kind=4), 'sol_cl.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,13),kind=4), 'sol_al.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,14),kind=4), 'sol_hco3.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(soluteMat(:,:,15),kind=4), 'sol_co3.txt' )

! primary minerals
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(primaryMat(:,:,1),kind=4), 'pri_feldspar.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(primaryMat(:,:,2),kind=4), 'pri_augite.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(primaryMat(:,:,3),kind=4), 'pri_pigeonite.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(primaryMat(:,:,4),kind=4), 'pri_magnetite.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(primaryMat(:,:,5),kind=4), 'pri_glass.txt' )

! medium properties
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(mediumMat(:,:,1),kind=4), 'med_phi.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(mediumMat(:,:,2),kind=4), 'med_s_sp.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(mediumMat(:,:,3),kind=4), 'med_v_water.txt' )
yep = write_matrix ( xn/cell, yn*tn/(cell*mstep), real(mediumMat(:,:,4),kind=4), 'med_rho_solid.txt' )


!--------------WRITE TO NETCDF FILES

! net cdf is on hold for the time being to accomodate the new geochemical output

! call check( nf90_create('thermalNRG.nc', NF90_CLOBBER, ncid) )
! call check( nf90_def_dim(ncid, "x", xn, x_dimid) )
! call check( nf90_def_dim(ncid, "y", yn, y_dimid) )

! call check(nf90_def_var(ncid, "h", NF90_FLOAT, (/ x_dimid, y_dimid /), h_varid) )
! call check(nf90_def_var(ncid, "u", NF90_FLOAT, (/ x_dimid, y_dimid /), u_varid) )
! call check(nf90_def_var(ncid, "v", NF90_FLOAT, (/ x_dimid, y_dimid /), v_varid) )
! call check(nf90_enddef(ncid))

! call check( nf90_put_var(ncid, h_varid, h) )
! call check( nf90_put_var(ncid, u_varid, u) )
! call check( nf90_put_var(ncid, v_varid, v) )

!  
! CALL check(nf90_put_att(ncid, x_varid, "units", "meter")) 
! CALL check(nf90_put_att(ncid, y_varid, "units", "meter")) 
! CALL check(nf90_put_att(ncid, h_varid, "units", "K"))
! CALL check(nf90_put_att(ncid, u_varid, "units", "m/s"))
! CALL check(nf90_put_att(ncid, v_varid, "units", "m/s"))
!  
! CALL check(nf90_put_att(ncid, x_varid, "axis", "X")) 
! CALL check(nf90_put_att(ncid, y_varid, "axis", "Y")) 


! Close the file. This frees up any internal netCDF resources
! associated with the file, and flushes any buffers.
!  call check( nf90_close(ncid) )

write(*,*) " "
write(*,*) "ALL DONE!"

! what to do if you are a slave processor
else
	
	!--------------SLAVE PROCESSOR RECEIVES MESSAGE
	
	! message receiving has to happen every mth step
	do jj = 1, tn/mstep
		! here is a slave process, each process must receive a chunk of the h array and 
		! take the local mean, print it, send it back.
		
		! receive size of temperature array chunk
		call MPI_RECV ( num_rows_to_receive, 1 , MPI_INTEGER, &
		root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
		
		! receive timestep size
		call MPI_RECV ( dt_local, 1 , MPI_DOUBLE_PRECISION, &
		root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
		
		! receive temperature array chunk, save in local hLocal
		call MPI_RECV ( hLocal, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
		root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
		num_rows_received = num_rows_to_receive

		! receive primary array chunk, save in local priLocal
		do ii = 1,g_pri
			call MPI_RECV ( priLocalBit, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
			root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			priLocal(:,ii) = priLocalBit
		end do
	
		! receive secondary array chunk, save in local secLocal
		do ii = 1,g_sec
			call MPI_RECV ( secLocalBit, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
			root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			secLocal(:,ii) = secLocalBit
		end do
	
		! receive solute chunk, save in local solLocal
		do ii = 1,g_sol
			call MPI_RECV ( solLocalBit, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
			root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			solLocal(:,ii) = solLocalBit
		end do
		
		! receive medium chunk, save in local solLocal
		do ii = 1,g_med
			call MPI_RECV ( medLocalBit, num_rows_to_receive, MPI_DOUBLE_PRECISION, &
			root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
			medLocal(:,ii) = medLocalBit
		end do
	
		!--------------SLAVE PROCESSOR RUNS GEOCHEMICAL MODEL

		! slave processor loops through each coarse cell
		do m=1,num_rows_to_receive
			
			yep = write_matrix ( g_pri, 1, real(priLocal(m,:), kind = 4), 'upPri0.txt' )
			yep = write_matrix ( g_sec, 1, real(secLocal(m,:), kind = 4), 'upSec0.txt' )
			yep = write_matrix ( g_sol, 1, real(solLocal(m,:), kind = 4), 'upSol0.txt' )
			yep = write_matrix ( g_med, 1, real(medLocal(m,:), kind = 4), 'upMed0.txt' )
			
			
			
			
			if (medLocal(m,5) .eq. 1.0) then
                           write(*,*) medLocal(m,6:7)
				! run the phreeqc alteration model
 				alt0 = alt_next(hLocal(m),dt_local*mstep,priLocal(m,:), &
 						    secLocal(m,:),solLocal(m,:),medLocal(m,:))
			end if


			! parse the phreeqc output
			! changed 5 -> 3
			solLocal(m,:) = (/ alt0(1,2), alt0(1,3), alt0(1,4), alt0(1,5), alt0(1,6), &
			alt0(1,7), alt0(1,8), alt0(1,9), alt0(1,10), alt0(1,11), alt0(1,12), &
			alt0(1,13), alt0(1,14), alt0(1,15), 0.0D+00/)

			secLocal(m,:) = (/ alt0(1,16), alt0(1,18), alt0(1,20), &
			alt0(1,22), alt0(1,24), alt0(1,26), alt0(1,28), alt0(1,30), alt0(1,32), alt0(1,34), &
			alt0(1,36), alt0(1,38), alt0(1,40), alt0(1,42), alt0(1,44), alt0(1,46), alt0(1,48), &
			alt0(1,50), alt0(1,52), alt0(1,54), alt0(1,56), alt0(1,58), alt0(1,60), alt0(1,62), &
			alt0(1,64), alt0(1,66), alt0(1,68), alt0(1,70), &
			alt0(1,72), alt0(1,74), alt0(1,76), alt0(1,78), alt0(1,80), alt0(1,82), alt0(1,84), &
			alt0(1,86), alt0(1,88), alt0(1,90), alt0(1,92), alt0(1,94), alt0(1,96), alt0(1,98), &
			alt0(1,100), alt0(1,102), alt0(1,104), alt0(1,106), alt0(1,108), alt0(1,110), alt0(1,112), &
			alt0(1,114), alt0(1,116), alt0(1,118), alt0(1,120), alt0(1,122), alt0(1,124), alt0(1,126), &
			alt0(1,128), alt0(1,130), alt0(1,132), &
			alt0(1,134), alt0(1,136), alt0(1,138), alt0(1,140), alt0(1,142), alt0(1,144), alt0(1,146), &
			alt0(1,148), alt0(1,150), alt0(1,152) /)

			!priLocal(m,:) = (/ alt0(1,72), alt0(1,74), alt0(1,76), alt0(1,78), alt0(1,80)/)
			priLocal(m,:) = (/ alt0(1,154), alt0(1,156), alt0(1,158), alt0(1,160), alt0(1,162)/)

			!medLocal(m,1:4) = (/ alt0(1,82), alt0(1,83), alt0(1,84), alt0(1,4)/)
			medLocal(m,1:4) = (/ alt0(1,164), alt0(1,165), alt0(1,4), alt0(1,166)/)

			! print something you want to look at
			!write(*,*) medLocal(m,3) ! water
			
			!write(*,*) alt0
			
		end do
		
		! write something interesting to a file to look while the model is running
		yep = write_matrix ( num_rows_to_receive, 5, real(priLocal, kind = 4), 'realTime.txt' )
	
		!--------------SLAVE PROCESSOR SENDS ALTERED MESSAGE TO MASTER PROCESSOR
	
		! send primary array chunk back to root process
		do ii = 1,g_pri
			call MPI_SEND( priLocal(:,ii), num_rows_received, MPI_DOUBLE_PRECISION, root_process, &
			return_data_tag, MPI_COMM_WORLD, ierr)
		end do
		
		! send secondary array chunk back to root process
		do ii = 1,g_sec
			call MPI_SEND( secLocal(:,ii), num_rows_received, MPI_DOUBLE_PRECISION, root_process, &
			return_data_tag, MPI_COMM_WORLD, ierr)
		end do
		
		! send solute array chunk back to root process
		do ii = 1,g_sol
			call MPI_SEND( solLocal(:,ii), num_rows_received, MPI_DOUBLE_PRECISION, root_process, &
			return_data_tag, MPI_COMM_WORLD, ierr)
		end do
		
		! send medium array chunk back to root process
		do ii = 1,g_med
			call MPI_SEND( medLocal(:,ii), num_rows_received, MPI_DOUBLE_PRECISION, root_process, &
			return_data_tag, MPI_COMM_WORLD, ierr)
		end do
		
		write(*,*) "SLAVE PROCESSOR IS DONE WITH WORK"
	
	! done with looping through coarse timesteps
	end do

! end loop through processors	
end if 


! close up shop
call MPI_FINALIZE ( ierr )


END PROGRAM main



! ----------------------------------------------------------------------------------%%
!
! H_NEXT
!
! SUMMARY: computes the 2D temperature profile for the current timestep
!
! INPUTS: h(xn,yn) : temperature profile of previous timestep
!         psi(xn,yn) : 2D streamfunction array
!         rho_in(xn,yn) : 2D density array
!         flux(xn,2) : top and bottom heat flux boundary conditions
!
! RETURNS: h_next(xn,yn) : temperature profile of current timestep
!
! ----------------------------------------------------------------------------------%%

function h_next (h, psi, rho_in, flux)
	
use globals
use initialize
implicit none

interface
	
	function velocities(psi)
		use globals
		implicit none
		integer :: i,ii,j,jj
		real(8) :: uu(xn+1,yn+1), vv(xn+1,yn+1), velocities(xn,2*yn), psi(xn,yn)
		real(8) :: u0(xn,yn), v0(xn,yn), dLeft, dRight, dTop, dBottom
	end function
	
end interface

! declare errthing

! integers
integer :: i, j, n, ii, m=3
! inputs 
real(8) :: sx, sy, qx, qy, rho_in(xn,yn), flux(xn,2)
! velocity stuff
real(8) :: u(xn,yn), v(xn,yn), uLong((xn-2)*(yn-2)), vLong((xn-2)*(yn-2))
real(8) ::  velocities0(xn,2*yn)
! matrix stuff
real(8) :: h(xn,yn), h_next(xn,yn), psi(xn,yn)
! real(8) :: aa((xn-2)*(yn-2),(xn-2)*(yn-2)), a((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
real(8) :: aBand((xn-2)*(yn-2),5), bBand((xn-2)*(yn-2),5)
! real(8) :: bb((xn-2)*(yn-2),(xn-2)*(yn-2)), b((xn-2)*(yn-2),(xn-2)*(yn-2)+1)
real(8) :: h0(xn,yn), uVec((xn-2)*(yn-2)), h_nextRow((xn-2)*(yn-2))
real(8) :: kMatLong((xn-2)*(yn-2))
real(8) :: mn(xn,yn)
real(8) :: sxMat(xn,yn), syMat(xn,yn), sxLong((xn-2)*(yn-2)), syLong((xn-2)*(yn-2))
  
  
mn = h

! calculate velocities from streamfunction values
velocities0 = velocities(psi)
u = velocities0(1:xn,1:yn)
v = velocities0(1:xn,yn+1:2*yn)
uLong = reshape(u(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
vLong = reshape(transpose(v(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))

! print maximum fluid velocities at each timestep
write(*,*) " "
write(*,*) "max u"
write(*,*) maxval(abs(u))
write(*,*) "max v"
write(*,*) maxval(abs(v))

! print stability conditions at each timestep
write(*,*) " "
write(*,*) "velocity check"
write(*,"(F10.5)") (dt*maxval(abs(u)))*rho_fluid/(dx)
write(*,"(F10.5)") (dt*maxval(abs(v)))*rho_fluid/(dy)
write(*,*) "conduction check"
write(*,"(F10.5)") (2.0*dt*maxval(lambdaMat))/(4186.0*dy*dy)
write(*,*) " "



h0 = h

qx = dt*rho_fluid/(dx)
qy = dt*rho_fluid/(dy)
sx = (2.0*dt*lambda)/(4186.0*dx*dx)
sy = (2.0*dt*lambda)/(4186.0*dy*dy)
sxMat = (2.0*dt*lambdaMat)/(4186.0*dx*dx)
syMat = (2.0*dt*lambdaMat)/(4186.0*dy*dy)

! vertical boundary conditions
h(2,:) = h(2,:) + h0(1,:)*sxMat(1,:)/2.0  ! left
h(xn-1,:) = h(xn-1,:) + h0(xn,:)*sxMat(xn,:)/2.0  ! right
 
uVec = reshape(h(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
sxLong = reshape(sxMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))
syLong = reshape(syMat(2:xn-1,2:yn-1), (/(xn-2)*(yn-2)/))

! make the band
aBand = 0.0
do i = 1,(xn-2)*(yn-2)
	  
	aBand(i,2) = 1.0+sxLong(i)
	if (i-1 .gt. 0) then
		aBand(i,1) = -sxLong(i)/2.0 - uLong(i)*qx/2.0
	end if
	if (i+1 .le. (xn-2)*(yn-2)) then
		aBand(i,3) = -sxLong(i)/2.0 + uLong(i)*qx/2.0
	end if

	! first edge
	if (any(mod((/i-1/),xn-2) .eq. 0.0)) then
		aBand(i,2) = 1.0 + sxLong(i) - uLong(i)*qx
		if (i .gt. 1) then
			aBand(i,1) =  0.0
		end if
		if (i .lt. (xn-2)*(yn-2)) then
			aBand(i,3) = -sxLong(i)/2.0 + uLong(i)*qx
		end if
	end if

	! last edge
	if (any(mod((/i/),xn-2) .eq. 0.0)) then
		aBand(i,2) = 1.0 + sxLong(i) - uLong(i)*qx
		if (i .gt. 1) then
			aBand(i,3) =  0.0
		end if
		if (i .le. (xn-2)*(yn-2)) then
			aBand(i,1) = -sxLong(i)/2.0 + uLong(i)*qx
		end if
	end if
  
end do
  
! make sure solver doesn't go out of bounds
do i=1,((xn-2)-1)
	ii = i*(xn-2)
	aBand(ii,3) = 0.0
	aBand(ii+1,1) = 0.0
end do
  
!!!!!!!!!!!! THIS !!!!!!!!!!!
h_nextRow = tridiag(aBand(:,1),aBand(:,2),aBand(:,3),uVec,(xn-2)*(yn-2))
h(2:xn-1,2:yn-1) = reshape(h_nextRow, (/xn-2, yn-2/))
sxMat(2:xn-1,2:yn-1) = reshape(sxLong, (/xn-2, yn-2/))
syMat(2:xn-1,2:yn-1) = reshape(syLong, (/xn-2, yn-2/))
sxLong = reshape(transpose(sxMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))
syLong = reshape(transpose(syMat(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))

! horizontal boundary conditions
h(:,2) = h(:,2) + flux(:,1)*syMat(:,1)/2.0 ! bottom
h(:,xn-1) = h(:,xn-1) + flux(:,2)*syMat(:,yn)/2.0 ! top

h_nextRow = reshape(transpose(h(2:xn-1,2:yn-1)), (/(xn-2)*(yn-2)/))

! make the band
bBand = 0.0
do i = 1,(xn-2)*(yn-2)
	bBand(i,2) = 1.0+syLong(i)
	if (i-1 .gt. 0) then
		bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qy/2.0
	end if
	if (i+1 .le. (xn-2)*(yn-2)) then
		bBand(i,3) = -syLong(i)/2.0 + vLong(i)*qy/2.0
	end if

	! first edge
	if (any(mod((/i-1/),xn-2) .eq. 0.0)) then
		bBand(i,2) = 1.0 + syLong(i) - vLong(i)*qy
		if (i .gt. 1) then
			bBand(i,1) =  0.0
		end if
		if (i .lt. (xn-2)*(yn-2)) then
			bBand(i,3) = -syLong(i)/2.0 + vLong(i)*qy
		end if
	end if

	! last edge
	if (any(mod((/i/),xn-2) .eq. 0.0)) then
		bBand(i,2) = 1.0 + syLong(i) + vLong(i)*qy
		if (i .gt. 1) then
			bBand(i,3) =  0.0
		end if
		if (i .le. (xn-2)*(yn-2)) then
			bBand(i,1) = -syLong(i)/2.0 - vLong(i)*qy
		end if
	end if
end do
  
! make sure solver doesn't go out of bounds
do i=1,((xn-2)-1)
	ii = i*(xn-2)
	bBand(ii,3) = 0.0
	bBand(ii+1,1) = 0.0
end do

h_nextRow = tridiag(bBand(:,1),bBand(:,2),bBand(:,3),h_nextRow,(xn-2)*(yn-2))
h_next(2:xn-1,2:yn-1) = transpose(reshape(h_nextRow, (/xn-2, yn-2/)))

! check out how this equation is converging to steady-state
write(*,*) "deltaT"
write(*,*) maxval(abs((mn(2:xn-1,2:yn-1)-h_next(2:xn-1,2:yn-1))/h_next(2:xn-1,2:yn-1)))

return

end function h_next

! ----------------------------------------------------------------------------------%%
!
! PSI_NEXT
!
! SUMMARY: computes the 2D streamfunction array of the current timestep
!
! INPUTS: h(xn,yn) : temperature profile
!         rhs0(xn,yn) : right hand side of streamfunction-vorticity equation
!         psi(xn,yn) : 2D streamfunction array of previous timestep
!         top_in(xn,1) : permeable upper boundary
!         rho_in(xn,yn) : 2D density array
!
! RETURNS: psi_next(xn,yn): 2D streamfunction array for current timestep
!
! ----------------------------------------------------------------------------------%%

function psi_next (h, rhs0, psi, top_in, rho_in)

use globals
use initialize
implicit none

interface
	
	function partial(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(8) :: array(rows,cols), d1, d2, d
		real(8) :: partial(rows,cols)
	end function partial

	function write_matrix ( m, n, table, filename )
		use globals
		implicit none
		integer :: m, n, j, output_status, unit0
		character ( len = * ) filename
		character ( len = 30 ) string
		real(4)  :: table(m,n) , write_matrix
	end function write_matrix
	
end interface

! declare errthing

! integers
integer :: i, j, ii, n, m
! inputs
real(8) :: rhs0(xn,yn), rhs1(xn,yn), rhsLong((xn-2)*(yn-2))
real(8) :: h(xn,yn), psi(xn,yn), rho_in(xn,yn)
! matrix stuff
real(8) :: uVec((xn-2)*(yn-2)), psiLong((xn)*(yn)), psi_nextRow((xn-2)*(yn-2))
real(8) :: psi_next(xn,yn)
real(8) :: top_in(xn,1)
real(8) :: mn(xn,yn)
! back to band
real(8) :: aBand0((xn-2)*(yn-2),2*(xn-2) + 1)

mn = psi

! calculate coefficients in solver matrix 
permx = partial((1/(permeability)),xn,yn,dx,dy,1)
permy = partial((1/(permeability)),xn,yn,dx,dy,2)

rhoLong = reshape(rho_in(2:xn-1,2:yn-1),(/(xn-2)*(yn-2)/))
permLong = reshape(permeability(2:xn-1,2:yn-1),(/(xn-2)*(yn-2)/))
permxLong = reshape(permx(2:xn-1,2:yn-1),(/(xn-2)*(yn-2)/))
permyLong = reshape(permy(2:xn-1,2:yn-1),(/(xn-2)*(yn-2)/))


rhs1 = rhs0

! account for boundary cells in central cells
rhs1(2,:) = rhs1(2,:) 
rhs1(xn-1,:) = rhs1(xn-1,:) 
rhs1(:,2) = rhs1(:,2) 
rhs1(:,yn-1) = rhs1(:,yn-1)
rhs1(:,yn-1) = rhs1(:,yn-1) +&
& top_in(:,1)*(1.0/(dy*dy*(permeability(:,yn))))

uVec = reshape(rhs1(2:xn-1,2:yn-1),(/(xn-2)*(yn-2)/))

psi_next = 0.0

! make the band
aBand0 = 0.0
m = 2*(xn-2) + 1
do i = 1,(xn-2)*(yn-2)
	
	! diagonal
	aBand0(i,(m+1)/2) = (2.0)/(permLong(i)*dx*dx) + (2.0)/(permLong(i)*dy*dy)
	
	! off-diagonals
	if (i .gt. 1) then
		aBand0(i,((m+1)/2)-1) = (-1.0)/(permLong(i)*dx*dx) + (permxLong(i))/(2.0*dx)
	end if
	if (i .lt. (xn-2)*(yn-2)) then
		aBand0(i,((m+1)/2)+1) = (-1.0)/(permLong(i)*dx*dx) - (permxLong(i))/(2.0*dx)
	end if
	
	! further off-diagonals
	if (i .le. (xn-2)*(yn-2)-(xn-2)) then
		aBand0(i,m) = (-1.0)/(permLong(i)*dy*dy) - (permyLong(i))/(2.0*dy)
	end if
	if (i .ge. (xn-2)) then
		aBand0(i,1) = (-1.0)/(permLong(i)*dy*dy) + (permyLong(i))/(2.0*dy)
	end if
	
end do
  
! make sure solver doesn't go out of bounds
do i=1,((xn-2)-1)
	ii = i*(xn-2)
	aBand0(ii,((m+1)/2)+1) = 0.0
	aBand0(ii+1,((m+1)/2)-1) = 0.0
end do
  
! use the banded solver here
aBand0 = band(aBand0,m,(xn-2)*(yn-2))
psi_nextRow = solve(aBand0,uVec,m,(xn-2)*(yn-2))
psi_next(2:xn-1,2:yn-1) = reshape(psi_nextRow, (/xn-2, yn-2/))

! here lies a jacobified version of the solver i used for a while
! i'm saving it because: 1) it can be parallelized so if it's up to snuff
! 						 that would be really useful and
!						 2) i was working on a multi-grid solver for
!						 three consecutive nights, late, late into the
! 						 night and it might be fun to revisit

!do n=1,800
!do i=2,xn-1
!do j=2,yn-1

! this is good.
!psi_next(i,j) = ((dx*dx*permeability(i,j)*rho_in(i,j))/4.0)&
!&*(psi_next(i+1,j)/(permeability(i,j)*rho_in(i,j)*dx*dx)&
!&+psi_next(i-1,j)/(permeability(i,j)*rho_in(i,j)*dx*dx)&
!&+(permy(i,j)*psi_next(i,j+1) - permy(i,j)*psi_next(i,j-1))/(dx)&
!&+psi_next(i,j+1)/(permeability(i,j)*rho_in(i,j)*dx*dx)&
!&+psi_next(i,j-1)/(permeability(i,j)*rho_in(i,j)*dx*dx)&
!&+rhs1(i,j))

!!!!!!!!! THIS BIT !!!!!!!!!!!!!
!psi_next(i,j) = (1/((4.0/(dx*dx*(permeability(i,j)*rho_in(i,j))))))&
!&*(psi_next(i+1,j)/(permeability(i,j)*rho_in(i,j)*dx*dx)&
!&+psi_next(i-1,j)/(permeability(i,j)*rho_in(i,j)*dx*dx)&
!&+(permy(i,j)*psi_next(i,j+1) - permy(i,j)*psi_next(i,j-1))/(dy)&
!&-(permx(i,j)*psi_next(i+1,j) + permx(i,j)*psi_next(i-1,j))/(dx)&
!&+psi_next(i,j+1)/(permeability(i,j)*rho_in(i,j)*dy*dy)&
!&+psi_next(i,j-1)/(permeability(i,j)*rho_in(i,j)*dy*dy)&
!&+rhs1(i,j))


!end do
!end do
!end do

write(*,*) "deltaPSI"
write(*,*) maxval(abs((mn(2:xn-1,2:yn-1)-psi_next(2:xn-1,2:yn-1))/psi_next(2:xn-1,2:yn-1)))


return

end function psi_next


! ----------------------------------------------------------------------------------%%
!
! SOLUTE_NEXT
!
! SUMMARY: transports solutes on coarse mesh
!
! INPUTS: sol(xn/cell,yn/cell) : 2D array of initial solute concentrations
!         uTransport(xn/cell,yn/cell) : lateral velocities (coarse mesh)
!         vTransport(xn/cell,yn/cell) : vertical velocities (coarse mesh)
!
! RETURNS: solute_next(xn/cell,yn/cell): 2D array of solute concentrations
!
! ----------------------------------------------------------------------------------%%

function solute_next (sol, uTransport, vTransport)
	
use globals
use initialize
implicit none

! declare errthing

! integers
integer :: i, j, ii, n, m
! inputs
real(8) :: sol(xn/cell,yn/cell), sol0(xn/cell,yn/cell)
real(8) :: uTransport(xn/cell,yn/cell), vTransport(xn/cell,yn/cell)
! solver stuff
! real(8) :: uLong((xn/cell-2)*(yn/cell-2)), vLong((xn/cell-2)*(yn/cell-2))
! real(8) :: aBand((xn/cell-2)*(yn/cell-2),5), bBand((xn/cell-2)*(yn/cell-2),5)
! real(8) :: qx, qy, solute_next(xn/cell,yn/cell), vec((xn/cell-2)*(yn/cell-2))
! real(8) :: sol_nextRow((xn/cell-2)*(yn/cell-2))
real(8) :: uLong((xn/cell)*(yn/cell)), vLong((xn/cell)*(yn/cell))
real(8) :: aBand((xn/cell)*(yn/cell),5), bBand((xn/cell)*(yn/cell),5)
real(8) :: qx, qy, solute_next(xn/cell,yn/cell), vec((xn/cell)*(yn/cell))
real(8) :: sol_nextRow((xn/cell)*(yn/cell))



sol0 = sol

qx = dt*mstep/(dx*cell)
qy = dt*mstep/(dy*cell)
!qx = 0.0
!qy = 0.0

write(*,*) "SOLUTE COURANT NUMBERS"
write(*,*) qx*maxval(abs(uTransport))
write(*,*) qy*maxval(abs(vTransport))

! uLong = reshape(uTransport(2:xn/cell-1,2:yn/cell-1), (/(xn/cell-2)*(yn/cell-2)/))
! vLong = reshape(transpose(vTransport(2:xn/cell-1,2:yn/cell-1)), (/(xn/cell-2)*(yn/cell-2)/))


uLong = reshape(uTransport(1:xn/cell,1:yn/cell), (/(xn/cell)*(yn/cell)/))
vLong = reshape(transpose(vTransport(1:xn/cell,1:yn/cell)), (/(xn/cell)*(yn/cell)/))

! ! kinda lame finite difference advection scheme (order dx**2)
! do i=2,(xn/cell)-1
! do ii=2,(yn/cell)-1
! 	solute_next(i,ii) = sol(i,ii) - uTransport(i,ii)*qx*(sol(i+1,ii)-sol(i-1,ii)) &
! 	& - vTransport(i,ii)*qy*(sol(i,ii+1)-sol(i,ii-1))
! end do
! end do

!!  fully 2d lax-wendroff advection scheme test (order dx**3?)
! do i=2,(xn/cell)-1
! 	do ii=2,(yn/cell)-1
! 		solute_next(i,ii) = sol(i,ii) &
! 		& - uTransport(i,ii) * qx * (sol(i+1,ii)-sol(i-1,ii)) &
! 		& - vTransport(i,ii) * qy * (sol(i,ii+1)-sol(i,ii-1)) &
! 		& - ((uTransport(i,ii)*qx)**2) * 0.5 * (sol(i+1,ii) - 2.0*sol(i,ii) + sol(i-1,ii)) &
! 		& - ((vTransport(i,ii)*qy)**2) * 0.5 * (sol(i,ii+1) - 2.0*sol(i,ii) + sol(i,ii-1)) &
! 		& + uTransport(i,ii) * vTransport(i,ii) * qx * qy * .25 * &
! 		& (sol(i+1,ii+1) - sol(i-1,ii+1) - sol(i+1,ii-1) + sol(i-1,ii-1))
!
! 		! get rid of out-of-bounds truncation errors
! 		if (solute_next(i,ii) .gt. maxval(sol0)) then
! 			solute_next(i,ii) = maxval(sol0)
! 		end if
! 		if (solute_next(i,ii) .le. 0.0) then
! 			solute_next(i,ii) = 1e-10
! 		end if
!
! 	end do
! end do
!
! do i=2,(xn/cell)-1
! 	! top edge
! 	if (vTransport(i,yn/cell) .ge. 0.0) then
! 		solute_next(i,yn/cell) = sol(i,yn/cell) &
! 		& - uTransport(i,yn/cell) * qx * (sol(i,yn/cell)-sol(i-1,yn/cell)) &
! 		& - vTransport(i,yn/cell) * qy * (sol(i,yn/cell)-sol(i,yn/cell-1))
!
! 		! get rid of out-of-bounds truncation errors
! 		if (solute_next(i,yn/cell) .gt. maxval(sol0)) then
! 			solute_next(i,yn/cell) = maxval(sol0)
! 		end if
! 		if (solute_next(i,yn/cell) .le. 0.0) then
! 			solute_next(i,yn/cell) = 1e-10
! 		end if
! 	end if
!
! 	! bottom edge
! 	solute_next(i,1) = sol(i,1) &
! 	& - uTransport(i,1) * qx * (sol(i+1,1)-sol(i,1)) &
! 	& - vTransport(i,1) * qy * (sol(i,1+1)-sol(i,1))
!
! 	! get rid of out-of-bounds truncation errors
! 	if (solute_next(i,1) .gt. maxval(sol0)) then
! 		solute_next(i,1) = maxval(sol0)
! 	end if
! 	if (solute_next(i,1) .le. 0.0) then
! 		solute_next(i,1) = 1e-10
! 	end if
!
! 	! left edge
! 	solute_next(1,i) = sol(1,i) &
! 	& - vTransport(1,i) * qy * (sol(1,i+1)-sol(1,i)) &
! 	& - uTransport(1,i) * qx * (sol(1+1,i)-sol(1,i))
!
! 	! get rid of out-of-bounds truncation errors
! 	if (solute_next(1,i) .gt. maxval(sol0)) then
! 		solute_next(1,i) = maxval(sol0)
! 	end if
! 	if (solute_next(1,i) .le. 0.0) then
! 		solute_next(1,i) = 1e-10
! 	end if
!
! 	! right edge
! 	solute_next(xn/cell,i) = sol(xn/cell,i) &
! 	& - vTransport(xn/cell,i) * qy * (sol(xn/cell,i)-sol(xn/cell,i-1)) &
! 	& - uTransport(xn/cell,i) * qx * (sol(xn/cell,i)-sol((xn/cell)-1,i))
!
! 	! get rid of out-of-bounds truncation errors
! 	if (solute_next(xn/cell,i) .gt. maxval(sol0)) then
! 		solute_next(xn/cell,i) = maxval(sol0)
! 	end if
! 	if (solute_next(xn/cell,i) .le. 0.0) then
! 		solute_next(xn/cell,i) = 1e-10
! 	end if
!
!
! end do

! here lies a fully implicit version of the same thing. it seems to work just as well, but
! in case i want to parallelize this part later, i might as well use the explicit method.
! also boundary conditions are easier to deal with in the explicit method

! do ii=2,(yn/cell)-1
! 	! right edge
! 	solute_next(xn/cell,ii) = sol(xn/cell,ii) &
! 	& - uTransport(xn/cell,ii) * qx * (sol(xn/cell,ii)-sol(xn/cell,ii-1)) &
! 	& - vTransport(xn/cell,ii) * qy * (sol(xn/cell,ii)-sol(xn/cell-1,ii))
!
! 	! left edge
! 	solute_next(1,ii) = sol(1,ii) &
! 	& - uTransport(1,ii) * qx * (sol(1,ii+1)-sol(1,ii)) &
! 	& - vTransport(1,ii) * qy * (sol(1+1,ii)-sol(1,ii))
! end do

! not using this implicit method
! VERTICAL BOUNDARY CONDITIONS
sol(2,:) = sol(2,:) ! left
sol(xn/cell-1,:) = sol(xn/cell-1,:) ! right

! vec = reshape(sol(2:xn/cell-1,2:yn/cell-1), (/(xn/cell-2)*(yn/cell-2)/))
vec = reshape(sol(1:xn/cell,1:yn/cell), (/(xn/cell)*(yn/cell)/))

! MAKE THE BAND
aBand = 0.0
do i = 1,(xn/cell)*(yn/cell)

	!aBand(i,2) = 1.0
	!if (i-1 .gt. 0) then
	!aBand(i,1) = - uLong(i)*qx/2.0
	!end if
	!if (i+1 .le. (xn/cell)*(yn/cell)) then
	!aBand(i,3) = uLong(i)*qx/2.0
	!end if
        ! switching solvers for a change
        !aBand(i,2) = 1.0 - uLong(i)*qx
	!if (i .gt. 1) then
	!aBand(i,1) =  0.0
	!end if
	!if (i .le. (xn/cell)*(yn/cell)) then
	!aBand(i,3) = uLong(i)*qx
	!end if

	! first edge
	if ((any(mod((/i-1/),xn/cell) .eq. 0.0)) .OR. (uLong(i) .le. 0.0)) then
		aBand(i,2) = 1.0 - uLong(i)*qx
		if (i .gt. 1) then
		aBand(i,1) =  0.0
		end if
		if (i .lt. (xn/cell)*(yn/cell)) then
		aBand(i,3) = uLong(i)*qx
		end if
	end if

	! last edge
	if ((any(mod((/i/),xn/cell) .eq. 0.0)) .OR. (uLong(i) .gt. 0.0)) then
		aBand(i,2) = 1.0 + uLong(i)*qx
		if (i .gt. 1) then
		aBand(i,3) =  0.0
		end if
		if (i .le. (xn/cell)*(yn/cell)) then
		aBand(i,1) = - uLong(i)*qx
		end if
	end if

	! first edge
	if (any(mod((/i-1/),xn/cell) .eq. 0.0)) then
		aBand(i,2) = 1.0 - uLong(i)*qx
		if (i .gt. 1) then
		aBand(i,1) =  0.0
		end if
		if (i .lt. (xn/cell)*(yn/cell)) then
		aBand(i,3) = uLong(i)*qx
		end if
	end if

	! last edge
	if (any(mod((/i/),xn/cell) .eq. 0.0)) then
		aBand(i,2) = 1.0 + uLong(i)*qx
		if (i .gt. 1) then
		aBand(i,3) =  0.0
		end if
		if (i .le. (xn/cell)*(yn/cell)) then
		aBand(i,1) = - uLong(i)*qx
		end if
	end if

end do

do i=1,((xn/cell)-1)
	ii = i*(xn/cell)
	aBand(ii,3) = 0.0
	aBand(ii+1,1) = 0.0
end do

!!!!!!!!!!!! THIS !!!!!!!!!!!
sol_nextRow = tridiag(aBand(:,1),aBand(:,2),aBand(:,3),vec,(xn/cell)*(yn/cell))
sol(1:xn/cell,1:yn/cell) = reshape(sol_nextRow, (/xn/cell, yn/cell/))

! HORIZONTAL BOUNDARY CONDITIONS
sol(:,2) = sol(:,2) ! bottom
sol(:,xn/cell-1) = sol(:,xn/cell-1) ! top

sol_nextRow = reshape(transpose(sol(1:xn/cell,1:yn/cell)), (/(xn/cell)*(yn/cell)/))



! MAKE THE BAND
bBand = 0.0
do i = 1,(xn/cell)*(yn/cell)
	!bBand(i,2) = 1.0
	!if (i-1 .gt. 0) then
	!bBand(i,1) = - vLong(i)*qy/2.0
	!end if
	!if (i+1 .le. (xn/cell)*(yn/cell)) then
	!bBand(i,3) = vLong(i)*qy/2.0
	!end if
        !switching solvers for a change
        !bBand(i,2) = 1.0 - vLong(i)*qy
	!if (i .gt. 1) then
	!bBand(i,1) =  0.0
	!end if
	!if (i .le. (xn/cell)*(yn/cell)) then
	!bBand(i,3) = vLong(i)*qy
	!end if

	! first edge x 2
	if ((any(mod((/i-1/),xn/cell) .eq. 0.0)) .OR. (vLong(i) .le. 0.0)) then
		bBand(i,2) = 1.0 - vLong(i)*qy
		if (i .gt. 1) then
		bBand(i,1) =  0.0
		end if
		if (i .lt. (xn/cell)*(yn/cell)) then
		bBand(i,3) = vLong(i)*qy
		end if
	end if

	! last edge x 2
	if ((any(mod((/i/),xn/cell) .eq. 0.0)) .OR. (vLong(i) .gt. 0.0)) then
		bBand(i,2) = 1.0 + vLong(i)*qy
		if (i .gt. 1) then
		bBand(i,3) =  0.0
		end if
		if (i .le. (xn/cell)*(yn/cell)) then
		bBand(i,1) = - vLong(i)*qy
		end if
	end if

	! first edge x 2
	if (any(mod((/i-1/),xn/cell) .eq. 0.0)) then
		bBand(i,2) = 1.0 - vLong(i)*qy
		if (i .gt. 1) then
		bBand(i,1) =  0.0
		end if
		if (i .lt. (xn/cell)*(yn/cell)) then
		bBand(i,3) = vLong(i)*qy
		end if
	end if

	! last edge x 2
	if (any(mod((/i/),xn/cell) .eq. 0.0)) then
		bBand(i,2) = 1.0 + vLong(i)*qy
		if (i .gt. 1) then
		bBand(i,3) =  0.0
		end if
		if (i .le. (xn/cell)*(yn/cell)) then
		bBand(i,1) = - vLong(i)*qy
		end if
	end if
	
end do

do i=1,(((xn/cell))-1)
	ii = i*((xn/cell))
	bBand(ii,3) = 0.0
	bBand(ii+1,1) = 0.0
end do

sol_nextRow = tridiag(bBand(:,1),bBand(:,2),bBand(:,3),sol_nextRow,((xn/cell))*((yn/cell)))
solute_next(1:(xn/cell),1:(yn/cell)) = transpose(reshape(sol_nextRow, (/(xn/cell), (yn/cell)/)))




!do i=1,xn/cell
!do ii=1,yn/cell
!	if (solute_next(i,ii) .gt. maxval(sol0)) then
!		solute_next(i,ii) = maxval(sol0)
!	end if
!	if (solute_next(i,ii) .lt. minval(sol0)) then
!		solute_next(i,ii) = minval(sol0)
!	end if
!end do
!end do


!do i=2,xn/cell-1
!do ii=2,yn/cell-1
!	if (solute_next(i,ii) .gt. maxval(sol0(i-1:i+1,ii-1:ii+1))) then
!		solute_next(i,ii) = sol0(i,ii)
!	end if
!	if (solute_next(i,ii) .lt. minval(sol0(i-1:i+1,ii-1:ii+1))) then
!		solute_next(i,ii) = sol0(i,ii)
!	end if
!end do
!end do

return

end function solute_next

! ----------------------------------------------------------------------------------%%
!
! ALT_NEXT
!
! SUMMARY: solves for equilibrium at a single grid cell using PHREEQC
!
! INPUTS: temp : temperature of grid cell
!         timestep : time elapsed
!         primaryList(5) : amounts of primary minerals
!         secondaryList(16) : amounts of secondary minerals
!         soluteList(11) : concentrations of solutes
!
! RETURNS: alt_next(1,altnum): returns everything from PHREEQC in a big pile
!          and it gets parsed in the main method's geochem loop
!
! ----------------------------------------------------------------------------------%%

function alt_next (temp, timestep, primaryList, secondaryList, soluteList, mediumList)
use globals
use initialize
use alteration
implicit none

interface

end interface

! declare errthing

integer :: order
real(8) :: temp, timestep
real(8) :: alt_next(1,167)
real(8) :: alter0(1,167)
real(8) :: primaryList(g_pri), secondaryList(g_sec), soluteList(g_sol), mediumList(g_med)

! use the alteration module
alter0 = alter(temp-272.9, timestep, primaryList, secondaryList, soluteList, mediumList)

! rename it for a reason that i now forget
alt_next = alter0

end function alt_next


! ----------------------------------------------------------------------------------%%
!
! RHO_NEXT
!
! SUMMARY : solves for density using linear thermally expansive equation of state
!
! INPUTS : h_in(xn,yn) : 2D temperature array of current timestep
!
! RETURNS : rho_next(xn,yn) : 2D density array of current timestep
!
! ----------------------------------------------------------------------------------%%

function rho_next(h_in)
	
use globals
use initialize
implicit none
  
! declare errthing
integer :: i,j
real(8) :: h_in(xn,yn), rho_next(xn,yn)


do i=1,xn
	do j = 1,yn
		rho_next(i,j) = rho_fluid*(1.0 - alpha*((h_in(i,j)-273.0)))
	end do 
end do

return

end function rho_next



! ----------------------------------------------------------------------------------%%
!
! VELOCITIES
!
! SUMMARY : computes the darcy velocity (specific discharge) from the streamfunction
!           using finite difference partial derivatives
!
! INPUTS : psi(xn,yn) : 2D streamfunction array of current timestep
!
! RETURNS : velocities(xn,2*yn) : both u and v velocities in one matrix
!
! ----------------------------------------------------------------------------------%%


function velocities(psi)
	
use globals
use initialize
implicit none

interface

	function partial(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(8) :: array(rows,cols), d1, d2, d
		real(8) :: partial(rows,cols)
	end function partial

end interface

! declare errthing
real(8) :: velocities(xn,2*yn), psi(xn,yn)
real(8) :: u0(xn,yn), v0(xn,yn)

u0 = partial(psi,xn,yn,dx,dy,2)
v0 = -partial(psi,xn,yn,dx,dy,1)

velocities(1:xn,1:yn) = u0
velocities(1:xn,yn+1:2*yn) = v0

return
end function velocities


! calculates velocities from COARSE streamfunction values
function velocitiesCoarse(psiCoarse)
	use globals
	use initialize
	implicit none
	real(8) :: velocitiesCoarse(xn/cell,2*yn/cell), psiCoarse(xn/cell,yn/cell)
	real(8) :: u0(xn/cell,yn/cell), v0(xn/cell,yn/cell)

interface

	function partial(array,rows,cols,d1,d2,dim)
		use globals
		use initialize
		implicit none
		integer :: rows, cols, dim, i, j, ii, jj
		real(8) :: array(rows,cols), d1, d2, d
		real(8) :: partial(rows,cols)
	end function partial

end interface

u0 = partial(psiCoarse,xn/cell,yn/cell,dx*cell,dy*cell,2)
v0 = -partial(psiCoarse,xn/cell,yn/cell,dx*cell,dy*cell,1)

velocitiesCoarse(1:xn/cell,1:yn/cell) = u0
velocitiesCoarse(1:xn/cell,yn/cell+1:2*yn/cell) = v0

return
end function velocitiesCoarse



! ----------------------------------------------------------------------------------%%
!
! PARTIAL
!
! SUMMARY : versatile function for solving for second-order accurate partial 
!           derivatives of 1D or 2D arrays with respect to specified dimension
!           
! INPUTS : array(rows,cols) : array to be partially differentiated
!          rows : number of rows
!          cols : number of columns
!          d1 : grid spacing in first dimension
!          d2 : grid spacing in second dimension
!          dim : dimension you differentiate w.r.t.
!
! ----------------------------------------------------------------------------------%%

function partial(array,rows,cols,d1,d2,dim)
	
use globals
use initialize
implicit none

! declare errthing
integer :: rows, cols, dim, i, j, ii, jj
real(8) :: array(rows,cols), d1, d2, d
real(8) :: partial(rows,cols)

! figure out which direction derivative goes (dx or dy)
if (dim .eq. 1) then
	ii = 1
	jj = 0
	d = d1
	
	! compute edges beforehand
	partial(1,:) = ( -3.0*array(1,:) + 4.0*array(2,:) -array(3,:)) / (2.0*d)
	partial(rows,:) = ( 3.0*array(rows,:) - 4.0*array(rows-1,:) + array(rows-2,:) ) / (2.0*d)
end if

if (dim .eq. 2) then
	ii = 0
	jj = 1
	d = d2
	
	! compute edges beforehand 
	partial(:,1) = ( -3.0*array(:,1) + 4.0*array(:,2) -array(:,3)) / (2.0*d)
	partial(:,cols) = ( 3.0*array(:,cols) - 4.0*array(:,cols-1) + array(:,cols-2) ) / (2.0*d)
end if

! use central difference method ignoring edges (already done)
do i = 2-jj,rows-1+jj
    do j = 2-ii,cols-1+ii
    	partial(i,j) = (array(i+ii,j+jj)-array(i-ii,j-jj))/(2.0*d)
	end do
end do

return
end function partial






! ----------------------------------------------------------------------------------%%
!
! WRITE_VEC
!
! SUMMARY: Write linspace-style vector to file
!
! INPUTS: n : dimension
!         vector : vector with data
!         filename : file name
!
! RETURNS: write_vec
!
! ----------------------------------------------------------------------------------%%

function write_vec ( n, vector, filename )
use globals
  implicit none
  integer :: n, j, output_status, unit0
  character ( len = * ) filename 
  real(4)  :: vector(n), write_vec



  unit0 = get_unit ()
  open ( unit = unit0, file = filename, status = 'replace', iostat = output_status )
  if ( output_status /= 0 ) then
    write ( *, '(a,i8)' ) 'COULD NOT OPEN OUTPUT FILE "' // &
      trim ( filename ) // '" USING UNIT ', unit0
    unit0 = -1
    stop
  end if
  

  if ( 0 < n ) then
    do j = 1, n
      write ( unit0, '(2x,g24.16)' ) vector(j)
    end do

  end if


  close ( unit = unit0 )
  write_vec = 1.0
  return
end function write_vec




! ----------------------------------------------------------------------------------%%
!
! WRITE_MATRIX
!
! SUMMARY: Write 2d array to file
!
! INPUTS: m,n : 2d dimensinons
!         table : 2d array with data
!         filename : file name
!
! RETURNS: write_matrix
!
! ----------------------------------------------------------------------------------%%

function write_matrix ( m, n, table, filename )
use globals
  implicit none
  integer :: m, n, j, output_status, unit0
  character ( len = * ) filename
  character ( len = 30 ) string
  real(4)  :: table(m,n) , write_matrix

  unit0 = get_unit ()
  open ( unit = unit0, file = filename, &
    status = 'replace', iostat = output_status )

  if ( output_status /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_WRITE - Fatal error!'
    write ( *, '(a,i8)' ) 'Could not open the output file "' // &
      trim ( filename ) // '" on unit ', unit0
    unit0 = -1
    stop
  end if


	write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) '(', m, 'g', 24, '.', 16, ')'

    do j = 1, n
      write ( unit0, string ) table(1:m,j)
    end do


  close ( unit = unit0 )
  write_matrix = 2.0
  return
end function write_matrix

