!cell.f90

! gfortran -c -I/usr/local/include -L/usr/local/lib -liphreeqc cell.f90
! gfortran -I/usr/local/include -L/usr/local/lib -liphreeqc cell.o
! ./a.out

! gfortran -I/usr/local/include -L/usr/local/lib -liphreeqc cell.o -o c10


program cell
!use globals
INCLUDE "IPhreeqc.f90.inc"

!implicit none
INTEGER(KIND=4) :: id, all=85, its=50, its0=1
INTEGER(KIND=4) :: i, j, jj
CHARACTER(LEN=72000) :: line
character(len=71200) :: inputz0
character(len=4) :: fake_in
real(8) :: alter(1,59), mix1= .9, mix2=0.1
real(8), allocatable :: outmat(:,:)

! REAL GRABS
real(8) :: glass ! primary
real(8) :: siderite ! secondary
real(8) :: temp, timestep, primary(5), secondary(28), solute(15), solute0(15) ! important information

! STRINGS
character(len=50) :: s_siderite, s_kaolinite, s_goethite, s_dolomite, s_celadonite ! secondary
character(len=50) :: s_sio2, s_albite, s_calcite, s_mont_na, s_smectite, s_saponite ! secondary
character(len=50) :: s_stilbite, s_dawsonite, s_magnesite, s_clinoptilolite, s_pyrite ! secondary
character(len=50) :: s_quartz, s_kspar, s_saponite_na, s_nont_na, s_nont_mg, s_nont_k ! secondary
character(len=50) :: s_nont_h, s_nont_ca, s_muscovite, s_mesolite, s_hematite, s_diaspore ! 
character(len=50) :: s_feldspar, s_pigeonite, s_augite, s_glass, s_magnetite ! primary
character(len=50) :: s_temp, s_timestep ! important information
character(len=50) :: s_ph, s_ca, s_mg, s_na, s_k, s_fe, s_s, s_si, s_cl, s_al, s_alk, s_co2 ! solutes
character(len=50) :: s_hco3, s_co3, s_pe
real(8) :: water
character(len=50) :: s_water

! command line arguments
character(len=100) :: intemp
character(len=100) :: infile
integer :: in

!in = iargc()
!call getarg(1,intemp)
!call getarg(2,infile)
intemp = "10.0"
infile = "cellOut.txt"

glass = 2.0
siderite = 0.0



! function alter ( temp, timestep, primary, secondary, solute )
! INITIALIZE CHEMISTRY, STUFF DONE THAT IS PASSED TO WHAT WAS ORIGINALLY A SUBROUTINE OR WHATEVER
primary(1) = 129.6 ! feldspar
primary(2) = 69.6 ! augite
primary(3) = 12.6 ! pigeonite
primary(4) = 4.0 ! magnetite
primary(5) = 967.7 ! basaltic glass


secondary(:) = 0.0

solute(2) = 6.0e-4 ! Ca
solute(3) = 2.0e-5 ! Mg
solute(4) = 1.0e-3 ! Na
solute(5) = 1.0e-4 ! K
solute(6) = 1.2e-6 ! Fe
solute(7) = 1.0e-4 ! S(6)
solute(8) = 2.0e-4 ! Si
solute(9) = 3.0e-4 ! Cl
solute(10) = 1.0e-6 ! Al
solute(11) = 2.3e-3 ! Alk 1.6e-3 
solute(12) = 2.200e-3 !1.2e-2 ! H2CO3
solute(13) = 0.0 ! HCO3-
solute(14) = 0.0 ! CO3-2
solute(1) = 7.8 ! ph
solute(15) = 8.451 ! pe

solute0 = solute



! flushing timestep
timestep = 3.14e9
temp = 10.0
water = 5.0

! SOLUTES TO STRINGS
write(s_ph,'(F25.10)') solute(1)
write(s_ca,'(F25.10)') solute(2)
write(s_mg,'(F25.10)') solute(3)
write(s_na,'(F25.10)') solute(4)
write(s_k,'(F25.10)') solute(5)
write(s_fe,'(F25.10)') solute(6)
write(s_s,'(F25.10)') solute(7)
write(s_si,'(F25.10)') solute(8)
write(s_cl,'(F25.10)') solute(9)
write(s_al,'(F25.10)') solute(10)
write(s_alk,'(F25.10)') solute(11)
write(s_co2,'(F25.10)') solute(12)
write(s_hco3,'(F25.10)') solute(13)
write(s_co3,'(F25.10)') solute(14)
write(s_pe,'(F25.10)') solute(15)
write(s_water,'(F25.10)') water

! PRIMARIES TO STRINGS
write(s_feldspar,'(F25.10)') primary(1)
write(s_augite,'(F25.10)') primary(2)
write(s_pigeonite,'(F25.10)') primary(3)
write(s_magnetite,'(F25.10)') primary(4)
write(s_glass,'(F25.10)') primary(5)

! SECONDARIES TO STRINGS
write(s_siderite,'(F25.10)') secondary(1)
write(s_kaolinite,'(F25.10)') secondary(2)
write(s_goethite,'(F25.10)') secondary(3)
write(s_dolomite,'(F25.10)') secondary(4)
write(s_celadonite,'(F25.10)') secondary(5)
write(s_sio2,'(F25.10)') secondary(6)
write(s_albite,'(F25.10)') secondary(7)
write(s_calcite,'(F25.10)') secondary(8)
write(s_mont_na,'(F25.10)') secondary(9)
write(s_smectite,'(F25.10)') secondary(10)
write(s_saponite,'(F25.10)') secondary(11)
write(s_stilbite,'(F25.10)') secondary(12)
write(s_dawsonite,'(F25.10)') secondary(13)
write(s_magnesite,'(F25.10)') secondary(14)
write(s_clinoptilolite,'(F25.10)') secondary(15)
write(s_pyrite,'(F25.10)') secondary(16)
write(s_quartz,'(F25.10)') secondary(17)
write(s_kspar,'(F25.10)') secondary(18)
write(s_saponite_na,'(F25.10)') secondary(19)
write(s_nont_na,'(F25.10)') secondary(20)
write(s_nont_mg,'(F25.10)') secondary(21)
write(s_nont_k,'(F25.10)') secondary(22)
write(s_nont_h,'(F25.10)') secondary(23)
write(s_nont_ca,'(F25.10)') secondary(24)
write(s_muscovite,'(F25.10)') secondary(25)
write(s_mesolite,'(F25.10)') secondary(26)
write(s_hematite,'(F25.10)') secondary(27)
write(s_diaspore,'(F25.10)') secondary(28)

! OTHER INFORMATION TO STRINGS
write(s_temp,'(F25.10)') temp
write(s_timestep,'(F25.10)') timestep

write(*,*) timestep




! NOW YOU HAVE TO LOOP THIS ENTIRE THING AND ADD EACH OUTPUT TO A MATRIX EACH TIME
allocate(outmat(its*its0+1,all))
do j = 1,its

! ----------------------------------%%
! INITIAL AQUEOUS PHASE CONSITUENTS
! ----------------------------------%%
 
inputz0 = "SOLUTION 1 " //CHAR(13)//CHAR(10)// &
&"    pH " // trim(s_pH) //CHAR(13)//CHAR(10)// &
&"    pe " // trim(s_pe) //CHAR(13)//CHAR(10)// &
&"    units   mol/kgw" //CHAR(13)//CHAR(10)// &
&"    temp " // intemp //CHAR(13)//CHAR(10)// &
&"    Ca " // trim(s_ca) //CHAR(13)//CHAR(10)// &
&"    Mg " // trim(s_mg) //CHAR(13)//CHAR(10)// &
&"    Na " // trim(s_na) //CHAR(13)//CHAR(10)// &
&"    K " // trim(s_k) //CHAR(13)//CHAR(10)// &
&"    Fe " // trim(s_fe) //CHAR(13)//CHAR(10)// &
&"    S(6) "// trim(s_s) // " as SO4" //CHAR(13)//CHAR(10)// &
&"    Si " // trim(s_si) //CHAR(13)//CHAR(10)// &
&"    Cl " // trim(s_cl) //CHAR(13)//CHAR(10)// &
&"    Al " // trim(s_al) //CHAR(13)//CHAR(10)// &
!&"    C " // trim(s_co2) //CHAR(13)//CHAR(10)// &
&"    Alkalinity " // trim(s_alk) // " as HCO3" //CHAR(13)//CHAR(10)// &
!&"    -water		5.0	# kg" //CHAR(13)//CHAR(10)// &
&"    -water "// trim(s_water) //CHAR(13)//CHAR(10)// &
&"END" //CHAR(13)//CHAR(10)// &

! ----------------------------------%%
! HYDROTHERMAL MINERAL CHOICES
! ----------------------------------%%
  
&"EQUILIBRIUM_PHASES 1" //CHAR(13)//CHAR(10)// &
!&"    CO2(g) -3.25 100" //CHAR(13)//CHAR(10)// &
&"    Kaolinite 0.0 " // trim(s_kaolinite) //CHAR(13)//CHAR(10)// &
&"    Goethite 0.0 " // trim(s_goethite) //CHAR(13)//CHAR(10)// &
&"    Celadonite 0.0 " // trim(s_celadonite) //CHAR(13)//CHAR(10)// &
&"    SiO2(am) 0.0 " // trim(s_sio2) //CHAR(13)//CHAR(10)// &
&"    Albite 0.0 " // trim(s_albite) //CHAR(13)//CHAR(10)// &
&"    Calcite 0.0 " // trim(s_calcite) //CHAR(13)//CHAR(10)// &
&"    Montmor-Na 0.0 " // trim(s_mont_na) //CHAR(13)//CHAR(10)// &
&"    Saponite-Mg 0.0 " // trim(s_saponite) //CHAR(13)//CHAR(10)// &
&"    Stilbite 0.0 " // trim(s_stilbite) //CHAR(13)//CHAR(10)// &
&"    Clinoptilolite-Ca 0.0 " // trim(s_clinoptilolite) //CHAR(13)//CHAR(10)// &
&"    Pyrite 0.0 " // trim(s_pyrite) //CHAR(13)//CHAR(10)// &
&"    Quartz 0.0 " // trim(s_quartz) //CHAR(13)//CHAR(10)// &
&"    K-Feldspar 0.0 " // trim(s_kspar) //CHAR(13)//CHAR(10)// &

! NEW MINS

! &"    Dolomite 0.0 " // trim(s_dolomite) //CHAR(13)//CHAR(10)// &
 &"    Saponite-Na 0.0 " // trim(s_saponite_na) //CHAR(13)//CHAR(10)// &
 &"    Nontronite-Na 0.0 " // trim(s_nont_na) //CHAR(13)//CHAR(10)// &
 &"    Nontronite-Mg 0.0 " // trim(s_nont_mg) //CHAR(13)//CHAR(10)// &
 &"    Nontronite-K 0.0 " // trim(s_nont_k) //CHAR(13)//CHAR(10)// &
 &"    Nontronite-H 0.0 " // trim(s_nont_h) //CHAR(13)//CHAR(10)// &
 &"    Nontronite-Ca 0.0 " // trim(s_nont_ca) //CHAR(13)//CHAR(10)// &
 &"    Muscovite 0.0 " // trim(s_muscovite) //CHAR(13)//CHAR(10)// &
 &"    Mesolite 0.0 " // trim(s_mesolite) //CHAR(13)//CHAR(10)// &
 &"    Hematite 0.0 " // trim(s_hematite) //CHAR(13)//CHAR(10)// &
 &"    Diaspore 0.0 " // trim(s_diaspore) //CHAR(13)//CHAR(10)// &


!  &"    Dawsonite 0.0 " // trim(s_dawsonite) //CHAR(13)//CHAR(10)// &
!  &"    Magnesite 0.0 " // trim(s_magnesite) //CHAR(13)//CHAR(10)// &
!  &"    Quartz 0.0 0.0" //CHAR(13)//CHAR(10)// &
!  &"    Smectite-high-Fe-Mg 0.0 " // trim(s_smectite) //CHAR(13)//CHAR(10)// &
!  &"    Dolomite 0.0 " // trim(s_dolomite) //CHAR(13)//CHAR(10)// &
!&"    Siderite 0.0 " // trim(s_siderite) //CHAR(13)//CHAR(10)// &

! ----------------------------------%%
! CALCULATE POROSITY AND STUFF
! ----------------------------------%%

&"CALCULATE_VALUES" //CHAR(13)//CHAR(10)// &

&"R(sum)" //CHAR(13)//CHAR(10)// &
&"-start" //CHAR(13)//CHAR(10)// &
&"10 sum = EQUI('Stilbite')*832.2/2.15 + EQUI('SiO2(am)')*60.0/2.62" //&
&"+ EQUI('Kaolinite')*258.2/2.6 + EQUI('Albite')*262.3/2.62" // &
&"+ EQUI('Saponite-Mg')*385.537/2.4 + EQUI('Celadonite')*396.8/3.0" // &
&"+ EQUI('Clinoptilolite-Ca')*1344.49/2.62 + EQUI('Pyrite')*120.0/4.84" // &
&"+ EQUI('Montmor-Na')*103.8/5.3 + EQUI('Goethite')*88.8/3.8" // &
&"+ EQUI('Dolomite')*184.3/2.84 + EQUI('Smectite-high-Fe-Mg')*425.7/2.7" // &
&"+ EQUI('Dawsonite')*144.0/2.42 + EQUI('Magnesite')*84.3/3.0" // &
&"+ EQUI('Siderite')*115.8/3.96 + EQUI('Calcite')*100.0/2.71" // &
&"+ EQUI('Quartz')*60.0/2.62 + EQUI('K-Feldspar')*193.0/2.56" // &
&"+ KIN('Plagioclase')*270.0/2.68 + KIN('Augite')*230.0/3.4" // &
&"+ KIN('Pigeonite')*239.6/3.38 + KIN('Magnetite')*231.0/5.15" // &
&"+ KIN('BGlass')*46.5/2.92" // &
&"" //CHAR(13)//CHAR(10)// &
&"100 SAVE sum" //CHAR(13)//CHAR(10)// &
&"-end" //CHAR(13)//CHAR(10)// &
  
&"R(phi)" //CHAR(13)//CHAR(10)// &
&"-start" //CHAR(13)//CHAR(10)// &
&"10 phi = 1.0-(CALC_VALUE('R(sum)')/(CALC_VALUE('R(sum)')+(TOT('water')*1000.0)))" //&
!&"10 phi = 0.1" //&
&"" //CHAR(13)//CHAR(10)// &
&"100 SAVE phi" //CHAR(13)//CHAR(10)// &
&"-end" //CHAR(13)//CHAR(10)// &

&"R(water_volume)" //CHAR(13)//CHAR(10)// &
&"-start" //CHAR(13)//CHAR(10)// &
&"10 water_volume = TOT('water')" //&
&"" //CHAR(13)//CHAR(10)// &
&"100 SAVE water_volume" //CHAR(13)//CHAR(10)// &
&"-end" //CHAR(13)//CHAR(10)// &
  

&"R(rho_s)" //CHAR(13)//CHAR(10)// &
&"-start" //CHAR(13)//CHAR(10)// &
!&"10 rho_s = CALC_VALUE('R(sum)')" //CHAR(13)//CHAR(10)// &
&"10 rho_s = EQUI('Stilbite')*2.15 + EQUI('SiO2(am)')*2.62" //&
&"+ EQUI('Kaolinite')*2.6 + EQUI('Albite')*2.62" // &
&"+ EQUI('Saponite-Mg')*2.4 + EQUI('Celadonite')*3.0" // &
&"+ EQUI('Clinoptilolite-Ca')*2.62 + EQUI('Pyrite')*4.84" // &
&"+ EQUI('Montmor-Na')*5.3 + EQUI('Goethite')*3.8" // &
&"+ EQUI('Dolomite')*2.84 + EQUI('Smectite-high-Fe-Mg')*2.7" // &
&"+ EQUI('Dawsonite')*2.42 + EQUI('Magnesite')*3.0" // &
&"+ EQUI('Siderite')*3.96 + EQUI('Calcite')*2.71" // &
&"+ EQUI('Quartz')*2.62 + EQUI('k-Feldspar')*2.56" // &
&"+ KIN('Plagioclase')*2.68 + KIN('Augite')*3.4" // &
&"+ KIN('Pigeonite')*3.38 + KIN('Magnetite')*5.15" // &
&"+ KIN('BGlass')*2.92" //CHAR(13)//CHAR(10)// &
&"20 rho_s = rho_s/ (EQUI('Stilbite') + EQUI('SiO2(am)')" //&
&"+ EQUI('Kaolinite') + EQUI('Albite')" // &
&"+ EQUI('Saponite-Mg') + EQUI('Celadonite')" // &
&"+ EQUI('Clinoptilolite-Ca') + EQUI('Pyrite')" // &
&"+ EQUI('Montmor-Na') + EQUI('Goethite')" // &
&"+ EQUI('Dolomite') + EQUI('Smectite-high-Fe-Mg')" // &
&"+ EQUI('Dawsonite') + EQUI('Magnesite')" // &
&"+ EQUI('Siderite') + EQUI('Calcite')" // &
&"+ EQUI('Quartz') + EQUI('K-Feldspar')" // &
&"+ KIN('Plagioclase') + KIN('Augite')" // &
&"+ KIN('Pigeonite') + KIN('Magnetite')" // &
&"+ KIN('BGlass'))" //CHAR(13)//CHAR(10)// &
&"30 rho_s = rho_s * 1000000.0" //CHAR(13)//CHAR(10)// &
&"100 SAVE rho_s" //CHAR(13)//CHAR(10)// &
&"-end" //CHAR(13)//CHAR(10)// &
  
  
&"R(s_sp)" //CHAR(13)//CHAR(10)// &
&"-start" //CHAR(13)//CHAR(10)// &
!&"10 s_sp = (CALC_VALUE('R(phi)')/(1.0-CALC_VALUE('R(phi)')))*400.0/CALC_VALUE('R(rho_s)')" //&
&"10 s_sp = 1.53e-5" //&
&"" //CHAR(13)//CHAR(10)// &
&"100 SAVE s_sp" //CHAR(13)//CHAR(10)// &
&"-end" //CHAR(13)//CHAR(10)// &

! ----------------------------------%%
! PRIMARY (KINETIC) CONSTITUENTS
! ----------------------------------%%

&"KINETICS" //CHAR(13)//CHAR(10)// &
&"Plagioclase" //CHAR(13)//CHAR(10)// &
&"-m0 " // trim(s_feldspar) //CHAR(13)//CHAR(10)// &
&"Augite" //CHAR(13)//CHAR(10)// &
&"-m0 " // trim(s_augite) //CHAR(13)//CHAR(10)// &
&"Pigeonite" //CHAR(13)//CHAR(10)// &
&"-m0 " // trim(s_pigeonite) //CHAR(13)//CHAR(10)// &
&"Magnetite" //CHAR(13)//CHAR(10)// &
&"-m0 " // trim(s_magnetite) //CHAR(13)//CHAR(10)// &
&"BGlass" //CHAR(13)//CHAR(10)// &
&"-f CaO 0.025 Fe2O3 0.0475 MgO 0.065 " //&
& "Na2O 0.0125 K2O 0.005 Al2O3 0.034999 SiO2 0.5" //CHAR(13)//CHAR(10)// &
&"-m0 " // trim(s_glass) //CHAR(13)//CHAR(10)// &
!&"-m0 0.0"  //CHAR(13)//CHAR(10)// &

! SO2 0.003

&"    -step " // trim(s_timestep) // " in 1" //CHAR(13)//CHAR(10)// &
!&"    -step 3.14e11 in 1" //CHAR(13)//CHAR(10)// &

&"INCREMENTAL_REACTIONS true" //CHAR(13)//CHAR(10)// &

&"Use solution 1" //CHAR(13)//CHAR(10)// &

    
! ----------------------------------%%
! KINETIC DISSOLUTION RATE LAWS
! ----------------------------------%%
	
&"RATES" //CHAR(13)//CHAR(10)// &

&"BGlass" //CHAR(13)//CHAR(10)// &
&"-start" //CHAR(13)//CHAR(10)// &
&"    10 rate0=M*46.5*CALC_VALUE('R(s_sp)')*0.1*(1e4)*(2.51189e-6)*exp(-25.5/(.008314*TK))" // &
&"*(((ACT('H+')^3)/(ACT('Al+3')))^.333)" //CHAR(13)//CHAR(10)// &
&"    20 save rate0 * time" //CHAR(13)//CHAR(10)// &
&"-end" //CHAR(13)//CHAR(10)// &

&"Plagioclase" //CHAR(13)//CHAR(10)// &
&"-start" //CHAR(13)//CHAR(10)// &
&"    10 rate = (1-SR('Plagioclase'))*M*270.0*CALC_VALUE('R(s_sp)')*0.1*(((1.58e-9)"//&
&"*exp(-53.5/(.008314*TK))*(ACT('H+')^0.541) +(3.39e-12)*exp(-57.4/(.008314*TK)) +"//&
&"(4.78e-15)*exp(-59.0/(.008314*TK))*(ACT('H+'))^-0.57))"//CHAR(13)//CHAR(10)//&
&"    20 save rate * time"//CHAR(13)//CHAR(10)//&
&"-end" //CHAR(13)//CHAR(10)// &

&"Augite" //CHAR(13)//CHAR(10)// &
&"-start" //CHAR(13)//CHAR(10)// &
&"    10 rate0 = (1-SR('Augite'))*M*230.0*CALC_VALUE('R(s_sp)')*0.1*(((1.58e-7)" // &
&"*exp(-78.0/(.008314*TK))*(ACT('H+')^0.7)+(1.07e-12)*exp(-78.0/(.008314*TK))))" //CHAR(13)//CHAR(10)// & 
&"    20 save rate0 * time" //CHAR(13)//CHAR(10)// &
&"-end" //CHAR(13)//CHAR(10)// &

&"Pigeonite" //CHAR(13)//CHAR(10)// &
&"-start" //CHAR(13)//CHAR(10)// &
&"    10 rate0 = (1-SR('Pigeonite'))*M*236.0*CALC_VALUE('R(s_sp)')*0.1*(((1.58e-7)" // &
&"*exp(-78.0/(.008314*TK))*(ACT('H+')^0.7)+(1.07e-12)*exp(-78.0/(.008314*TK))))"//CHAR(13)//CHAR(10)// &
&"    20 save rate0 * time" //CHAR(13)//CHAR(10)// &
&"-end" //CHAR(13)//CHAR(10)// &

&"Magnetite" //CHAR(13)//CHAR(10)// &
&"-start" //CHAR(13)//CHAR(10)// &
&"    10 rate0 = (1-SR('Magnetite'))*M*231.0*CALC_VALUE('R(s_sp)')*0.1*(((2.57e-9)" // &
&"*exp(-18.6/(.008314*TK))*(ACT('H+')^0.279)+(1.66e-11)*exp(-18.6/(.008314*TK))))" //CHAR(13)//CHAR(10)// &
&"    20 save rate0 * time" //CHAR(13)//CHAR(10)// &
&"-end" //CHAR(13)//CHAR(10)// &


! ----------------------------------%%
! DEFINE THE KIND OF OUTPUT
! ----------------------------------%%

&"DUMP" //CHAR(13)//CHAR(10)// &
&"    -solution 1" //CHAR(13)//CHAR(10)// &
&"    -equilibrium_phases" //CHAR(13)//CHAR(10)// &

  &"SELECTED_OUTPUT" //CHAR(13)//CHAR(10)// &
  &"    -reset false" //CHAR(13)//CHAR(10)// &
  &"    -high_precision true" //CHAR(13)//CHAR(10)// &
  &"    -k plagioclase augite pigeonite magnetite bglass" //CHAR(13)//CHAR(10)// &
  &"    -ph" //CHAR(13)//CHAR(10)// &
  &"    -pe" //CHAR(13)//CHAR(10)// &
  &"    -molalities Ca+2 Mg+2 Na+ K+ Fe+3 SO4-2 SiO2 Cl- Al+3 HCO3-" //CHAR(13)//CHAR(10)// &
  &"    -totals C" //CHAR(13)//CHAR(10)// &
  &"    -alkalinity" //CHAR(13)//CHAR(10)// &
!  &"    -molalities HCO3-" //CHAR(13)//CHAR(10)// &
  &"    -p stilbite sio2(am) kaolinite albite saponite-mg celadonite Clinoptilolite-Ca" //CHAR(13)//CHAR(10)// &
  &"    -p pyrite Montmor-Na goethite dolomite Smectite-high-Fe-Mg Dawsonite" //CHAR(13)//CHAR(10)// &
  &"    -p magnesite siderite calcite quartz k-feldspar" //CHAR(13)//CHAR(10)// &
  &"    -p saponite-na Nontronite-Na Nontronite-Mg Nontronite-K Nontronite-H " //CHAR(13)//CHAR(10)// &
  &"    -p Nontronite-Ca muscovite mesolite hematite diaspore" //CHAR(13)//CHAR(10)// &
  &"    -calculate_values R(phi) R(s_sp) R(water_volume) R(rho_s)" //CHAR(13)//CHAR(10)// &
  &"    -time" //CHAR(13)//CHAR(10)// &
  &"END"
  
  
! INITIALIZE STUFF
id = CreateIPhreeqc()
IF (id.LT.0) THEN
	STOP
END IF

IF (SetSelectedOutputStringOn(id, .TRUE.).NE.IPQ_OK) THEN
	CALL OutputErrorString(id)
	STOP
END IF

IF (SetOutputStringOn(id, .TRUE.).NE.IPQ_OK) THEN
	CALL OutputErrorString(id)
	STOP
END IF
  
IF (LoadDatabase(id, "llnl.dat").NE.0) THEN
	CALL OutputErrorString(id)
	STOP
END IF

! RUN INPUT
IF (RunString(id, inputz0).NE.0) THEN
	CALL OutputErrorString(id)
	STOP
END IF
  
! PRINT DUMP/OUTPUT
DO i=1,GetOutputStringLineCount(id)
	call GetOutputStringLine(id, i, line)
	!write(*,*) trim(line)
END DO

! WRITE AWAY
!allocate(outmat(GetSelectedOutputStringLineCount(id)+1,all))
DO i=1,GetSelectedOutputStringLineCount(id)
	call GetSelectedOutputStringLine(id, i, line)
	! HEADER BITS YOU MAY WANT
	if ((i .eq. 1) .and. (j .eq. 1)) then
 	   !write(12,*) trim(line)
	   write(*,*) trim(line) ! PRINT LABELS FOR EVERY FIELD (USEFUL)
	end if
	! MEAT
	if (i .gt. 1) then
		!write(*,*) trim(line)
		!write(*,*) "line:", i
		!write(*,*) "timestep:", its0*j+i-its0-1
		outmat(its0*j+i-its0-1,1) = its0*j+i-its0-1
		read(line,*) outmat(its0*j+i-its0-1,2:)
	end if
END DO
! DESTROY INSTANCE
IF (DestroyIPhreeqc(id).NE.IPQ_OK) THEN
	STOP
END IF

jj = its0*j

write(*,*) "TIMESTEP:", jj
!write(*,*) primary
!write(*,*) secondary
!write(*,*) solute

! PUT IN VALUES FOR THE NEXT TIMESTEP
primary(1) = outmat(jj,73) ! feldspar
primary(2) = outmat(jj,75) ! augite
primary(3) = outmat(jj,77) ! pigeonite
primary(4) = outmat(jj,79) ! magnetite
primary(5) = outmat(jj,81) ! basaltic glass
secondary(1) = outmat(jj,45)
secondary(2) = outmat(jj,21)
secondary(3) = outmat(jj,35)
secondary(4) = outmat(jj,37)
secondary(5) = outmat(jj,27)
secondary(6) = outmat(jj,19)
secondary(7) = outmat(jj,23)
secondary(8) = outmat(jj,47)
secondary(9) = outmat(jj,33)
secondary(10) = outmat(jj,39)
secondary(11) = outmat(jj,25)
secondary(12) = outmat(jj,17)
secondary(13) = outmat(jj,41)
secondary(14) = outmat(jj,43)
secondary(15) = outmat(jj,29)
secondary(16) = outmat(jj,31)
secondary(17) = outmat(jj,49)
secondary(18) = outmat(jj,51)
secondary(19) = outmat(jj,53)
secondary(20) = outmat(jj,55)
secondary(21) = outmat(jj,57)
secondary(22) = outmat(jj,59)
secondary(23) = outmat(jj,61)
secondary(24) = outmat(jj,63)
secondary(25) = outmat(jj,65)
secondary(26) = outmat(jj,67)
secondary(27) = outmat(jj,69)
secondary(28) = outmat(jj,71)
solute(1) = -log10(mix1*10.0**(-outmat(jj,3)) + mix2*10.0**(-solute0(1)))
!solute(1) = outmat(jj,3)
solute(2) = outmat(jj,7)*mix1 + solute0(2)*mix2
solute(3) = outmat(jj,8)*mix1 + solute0(3)*mix2
solute(4) = outmat(jj,9)*mix1 + solute0(4)*mix2
solute(5) = outmat(jj,10)*mix1 + solute0(5)*mix2
solute(6) = outmat(jj,11)*mix1 + solute0(6)*mix2
solute(7) = outmat(jj,12)*mix1 + solute0(7)*mix2
solute(8) = outmat(jj,13)*mix1 + solute0(8)*mix2
solute(9) = outmat(jj,14)*mix1 + solute0(9)*mix2
solute(10) = outmat(jj,15)*mix1 + solute0(10)*mix2
solute(11) = outmat(jj,5)*mix1 + solute0(11)*mix2
!solute(12) = (outmat(jj,6))*mix1 + solute0(12)*mix2
solute(15) = -log10(mix1*10.0**(-outmat(jj,4)) + mix2*10.0**(-solute0(15)))
!solute(15) = outmat(jj,4)
water = outmat(jj,85) 

! SOLUTES TO STRINGS
write(s_ph,'(F25.10)') solute(1)
write(s_ca,'(F25.10)') solute(2)
write(s_mg,'(F25.10)') solute(3)
write(s_na,'(F25.10)') solute(4)
write(s_k,'(F25.10)') solute(5)
write(s_fe,'(F25.10)') solute(6)
write(s_s,'(F25.10)') solute(7)
write(s_si,'(F25.10)') solute(8)
write(s_cl,'(F25.10)') solute(9)
write(s_al,'(F25.10)') solute(10)
write(s_alk,'(F25.10)') solute(11)
write(s_co2,'(F25.10)') solute(12)
write(s_hco3,'(F25.10)') solute(13)
write(s_co3,'(F25.10)') solute(14)
write(s_water,'(F25.10)') water


! PRIMARIES TO STRINGS
write(s_feldspar,'(F25.10)') primary(1)
write(s_augite,'(F25.10)') primary(2)
write(s_pigeonite,'(F25.10)') primary(3)
write(s_magnetite,'(F25.10)') primary(4)
write(s_glass,'(F25.10)') primary(5)

! SECONDARIES TO STRINGS
write(s_siderite,'(F25.10)') secondary(1)
write(s_kaolinite,'(F25.10)') secondary(2)
write(s_goethite,'(F25.10)') secondary(3)
write(s_dolomite,'(F25.10)') secondary(4)
write(s_celadonite,'(F25.10)') secondary(5)
write(s_sio2,'(F25.10)') secondary(6)
write(s_albite,'(F25.10)') secondary(7)
write(s_calcite,'(F25.10)') secondary(8)
write(s_mont_na,'(F25.10)') secondary(9)
write(s_smectite,'(F25.10)') secondary(10)
write(s_saponite,'(F25.10)') secondary(11)
write(s_stilbite,'(F25.10)') secondary(12)
write(s_dawsonite,'(F25.10)') secondary(13)
write(s_magnesite,'(F25.10)') secondary(14)
write(s_clinoptilolite,'(F25.10)') secondary(15)
write(s_pyrite,'(F25.10)') secondary(16)
write(s_quartz,'(F25.10)') secondary(17)
write(s_kspar,'(F25.10)') secondary(18)
write(s_saponite_na,'(F25.10)') secondary(19)
write(s_nont_na,'(F25.10)') secondary(20)
write(s_nont_mg,'(F25.10)') secondary(21)
write(s_nont_k,'(F25.10)') secondary(22)
write(s_nont_h,'(F25.10)') secondary(23)
write(s_nont_ca,'(F25.10)') secondary(24)
write(s_muscovite,'(F25.10)') secondary(25)
write(s_mesolite,'(F25.10)') secondary(26)
write(s_hematite,'(F25.10)') secondary(27)
write(s_diaspore,'(F25.10)') secondary(28)

!write(*,*) s_dolomite
write(*,*) s_calcite

! END BIG LOOP
end do






!write(*,*) outmat(:,43)

! WRITE TO FILE
OPEN(UNIT=12, FILE=infile, ACTION="write", STATUS="replace") 
do i=1,its*its0
	write(12,*) outmat(i,:)
end do

! ALL DONE!
write(*,*) "this is phreeqing me out"




end program cell