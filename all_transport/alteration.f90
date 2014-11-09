module alteration
!use globals
INCLUDE "IPhreeqc.f90.inc"
save


! declare things here

contains
	
	
! ----------------------------------------------------------------------------------%%
!
! ALTERATION
! 
! SUMMARY: this function... alters basalt
!
! INPUTS: solutes : solutes in the aqueous phase to be defined in SOLUTION 1 block
!         primaries : primary (kinetic) constituents available for dissolution
!                     but this is the only place they will change, so, maybe
!                     i can get around using them as an actual input
!         secondaries : secondary (hydrothermally generated) minerals
!                     but this is the only place they will change, so, maybe
!                     i can get around using them as an actual input
!         timestep : length of timestep in seconds
!         temperature : of grid cell
!         porosity : of the grid cell
!         
!
! RETURNS: alt0 array to be parsed by main method
!          
!          
!
! 1 time
! 2 pH
! 3 pe
! 4 Al(k)
! 5 C
! 6 m_Ca+2
! 7 m_Mg+2
! 8 m_Na
! 9 m_K+
! 10 m_Fe+3
! 11 m_SO4-2
! 12 m_SiO
! 13 m_Cl-
! 14 m_Al+3
! 15 m_HCO3-
! 16 stilbit
! 17 d_stilbite
! 18 sio2(am)
! 19 d_sio2(am)
! 20 kaolinit
! 21 d_kaolinite
! 22 albite
! 23 d_albite
! 24 saponite-m
! 25 d_saponite-mg
! 26 celadonite
! 27 d_celadonite
! 28 Clinoptilolite-C
! 29 d_Clinoptilolite-Ca
! 30 pyrite
! 31 d_pyrite
! 32 Montmor-N
! 33 d_Montmor-Na
! 34 goethite
! 35 d_goethite
! 36 dolomit
! 37 d_dolomite
! 38 Smectite-high-Fe-Mg
! 39 d_Smectite-high-Fe-Mg
! 40 Dawsonit
! 41 d_Dawsonite
! 42 magnesite
! 43 d_magnesite
! 44 siderit
! 45 d_siderite
! 46 calcite
! 47 d_calcite
! 48 quart
! 49 d_quartz
! 50 k-feldspar
! 51 d_k-feldspar
! 52 saponite-n
! 53 d_saponite-na
! 54 Nontronite-Na
! 55 d_Nontronite-Na
! 56 Nontronite-M
! 57 d_Nontronite-Mg
! 58 Nontronite-K
! 59 d_Nontronite-K
! 60 Nontronite-
! 61 d_Nontronite-H
! 62 Nontronite-Ca
! 63 d_Nontronite-Ca
! 64 muscovit
! 65 d_muscovite
! 66 mesolite
! 67 d_mesolite
! 68 hematit
! 69 d_hematite
! 70 diaspore
! 71 d_diaspore
! 72 k_plagioclas
! 73 dk_plagioclase
! 74 k_augite
! 75 dk_augite
! 76 k_pigeonit
! 77 dk_pigeonite
! 78 k_magnetite
! 79 dk_magnetite
! 80 k_bglas
! 81 dk_bglass
! 82 V_R(phi)
! 83 V_R(s_sp)
! 84 V_R(water_volume
! 85 V_R(rho_s)
!
! ----------------------------------------------------------------------------------%%

function alter ( temp, timestep, primary, secondary, solute, medium)

implicit none
INTEGER(KIND=4) :: id, all=167
INTEGER(KIND=4) :: i, order
CHARACTER(LEN=63200) :: line
character(len=63200) :: inputz0
character(len=4) :: fake_in
real(8) :: alter(1,167)
real(8), allocatable :: outmat(:,:)


!!!!! YOU'RE GONNA HAVE TO CONVERT THINGS FROM FLOATS TO STRINGS AND BACK AGAIN.
!!!!! HERE IS HOW I DID IT OVER THE SUMMER

!!!!! in the inputz0 input string, variables are concatenated as such: // trim(s_siderite) //

!!!!! in the MEAT of the output, you gotta do:
! MEAT
! if (i .gt. 1) then
!              read(line,*) temp, plag, d_plag, aug, d_aug, pig, d_pig, glass, d_glass, &
!              & mag, d_mag, si_plag, si_aug, si_pig, si_glass, si_mag
!              write(*,'(11F16.5)') temp, plag, d_plag, aug, d_aug, pig, d_pig, glass, d_glass, &
!              & mag, d_mag, si_plag, si_aug, si_pig, si_glass, si_mag
! end if
!!!!! instead of 
! MEAT
! if (i .gt. 1) then
! 	read(line,*) outmat(i,:)
! 	!write(12,"(4F12.5)") outmat(i,:)
! 	write(12,*) outmat(i,:)
! 	write(*,*) trim(line)
! end if

! REAL GRABS
real(8) :: glass ! primary
real(8) :: siderite ! secondary
real(8) :: temp, timestep, primary(5), secondary(69), solute(15), medium(7) ! important information

! STRINGS

character(len=50) :: s_akermanite, s_analcime, s_annite, s_clinozoisite, s_dicalcium_silicate
character(len=50) :: s_diopside, s_epidote, s_ettringite, s_ferrite_ca, s_foshagite
character(len=50) :: s_gismondine, s_gyrolite, s_hedenbergite, s_hillebrandite, s_larnite
character(len=50) :: s_laumontite, s_lawsonite, s_merwinite, s_monticellite, s_natrolite
character(len=50) :: s_okenite, s_phlogopite, s_prehnite, s_pseudowollastonite, s_rankinite 
character(len=50) :: s_scolecite, s_tobermorite_9a, s_tremolite, s_wollastonite, s_xonotlite 
character(len=50) :: s_zoisite
character(len=50) :: s_andradite, s_troilite, s_pyrrhotite, s_minnesotaite, s_fayalite, s_daphnite_7a
character(len=50) :: s_daphnite_14a, s_cronstedtite_7a, s_greenalite, s_aragonite
character(len=50) :: s_siderite, s_kaolinite, s_goethite, s_dolomite, s_celadonite ! secondary
character(len=50) :: s_sio2, s_albite, s_calcite, s_mont_na, s_smectite, s_saponite ! secondary
character(len=50) :: s_stilbite, s_dawsonite, s_magnesite, s_clinoptilolite, s_pyrite ! secondary
character(len=50) :: s_quartz, s_kspar, s_saponite_na, s_nont_na, s_nont_mg, s_nont_k ! secondary
character(len=50) :: s_nont_h, s_nont_ca, s_muscovite, s_mesolite, s_hematite, s_diaspore ! 
character(len=50) :: s_feldspar, s_pigeonite, s_augite, s_glass, s_magnetite ! primary
character(len=50) :: s_temp, s_timestep ! important information
character(len=50) :: s_ph, s_ca, s_mg, s_na, s_k, s_fe, s_s, s_si, s_cl, s_al, s_alk, s_co2 ! solutes
character(len=50) :: s_hco3, s_co3, s_pe
character(len=50) :: s_water, s_w ! medium
real(8) :: water

glass = 2.0
siderite = 0.0

! if (solute(6) .lt. 1e-10) then
! 	solute(6) = 1e-10
! end if

! SOLUTES TO STRINGS
write(s_ph,'(F25.10)') solute(1)
write(s_alk,'(F25.10)') solute(2)
write(s_water,'(F25.10)') solute(3)
write(s_co2,'(F25.10)') solute(4)
write(s_ca,'(F25.10)') solute(5)
write(s_mg,'(F25.10)') solute(6)
write(s_na,'(F25.10)') solute(7)
write(s_k,'(F25.10)') solute(8)
write(s_fe,'(F25.10)') solute(9)
write(s_s,'(F25.10)') solute(10)
write(s_si,'(F25.10)') solute(11)
write(s_cl,'(F25.10)') solute(12)
write(s_al,'(F25.10)') solute(13)
write(s_hco3,'(F25.10)') solute(14)
write(s_co3,'(F25.10)') solute(15)

! MEDIUM TO STRINGS
write(s_w,'(F25.10)') solute(3)

! PRIMARIES TO STRINGS
write(s_feldspar,'(F25.10)') primary(1)
write(s_augite,'(F25.10)') primary(2)
write(s_pigeonite,'(F25.10)') primary(3)
write(s_magnetite,'(F25.10)') primary(4)
write(s_glass,'(F25.10)') primary(5)

! SECONDARIES TO STRINGS
!secondary = 0.0
write(s_stilbite,'(F25.10)') secondary(1)
write(s_sio2,'(F25.10)') secondary(2)
write(s_kaolinite,'(F25.10)') secondary(3)
write(s_albite,'(F25.10)') secondary(4)
write(s_saponite,'(F25.10)') secondary(5)
write(s_celadonite,'(F25.10)') secondary(6)
write(s_clinoptilolite,'(F25.10)') secondary(7)
write(s_pyrite,'(F25.10)') secondary(8)
write(s_mont_na,'(F25.10)') secondary(9)
write(s_goethite,'(F25.10)') secondary(10)
write(s_dolomite,'(F25.10)') secondary(11)
write(s_smectite,'(F25.10)') secondary(12)
write(s_dawsonite,'(F25.10)') secondary(13)
write(s_magnesite,'(F25.10)') secondary(14)
write(s_siderite,'(F25.10)') secondary(15)
write(s_calcite,'(F25.10)') secondary(16)
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
! MINERALS ADDED 10/18/2104
write(s_akermanite,'(F25.10)') secondary(29)
write(s_analcime,'(F25.10)') secondary(30)
write(s_annite,'(F25.10)') secondary(31)
write(s_clinozoisite,'(F25.10)') secondary(32)
write(s_dicalcium_silicate,'(F25.10)') secondary(33)
write(s_diopside,'(F25.10)') secondary(34)
write(s_epidote,'(F25.10)') secondary(35)
write(s_ettringite,'(F25.10)') secondary(36)
write(s_ferrite_ca,'(F25.10)') secondary(37)
write(s_foshagite,'(F25.10)') secondary(38)
write(s_gismondine,'(F25.10)') secondary(39)
write(s_gyrolite,'(F25.10)') secondary(40)
write(s_hedenbergite,'(F25.10)') secondary(41)
write(s_hillebrandite,'(F25.10)') secondary(42)
write(s_larnite,'(F25.10)') secondary(43)
write(s_laumontite,'(F25.10)') secondary(44)
write(s_lawsonite,'(F25.10)') secondary(45)
write(s_merwinite,'(F25.10)') secondary(46)
write(s_monticellite,'(F25.10)') secondary(47)
write(s_natrolite,'(F25.10)') secondary(48)
write(s_okenite,'(F25.10)') secondary(49)
write(s_phlogopite,'(F25.10)') secondary(50)
write(s_prehnite,'(F25.10)') secondary(51)
write(s_pseudowollastonite,'(F25.10)') secondary(52)
write(s_rankinite,'(F25.10)') secondary(53)
write(s_scolecite,'(F25.10)') secondary(54)
write(s_tobermorite_9a,'(F25.10)') secondary(55)
write(s_tremolite,'(F25.10)') secondary(56)
write(s_wollastonite,'(F25.10)') secondary(57)
write(s_xonotlite,'(F25.10)') secondary(58)
write(s_zoisite,'(F25.10)') secondary(59)
!! NEXT ROUND
write(s_andradite,'(F25.10)') secondary(60)
write(s_troilite,'(F25.10)') secondary(61)
write(s_pyrrhotite,'(F25.10)') secondary(62)
write(s_minnesotaite,'(F25.10)') secondary(63)
write(s_fayalite,'(F25.10)') secondary(64)
write(s_daphnite_7a,'(F25.10)') secondary(65)
write(s_daphnite_14a,'(F25.10)') secondary(66)
write(s_cronstedtite_7a,'(F25.10)') secondary(67)
write(s_greenalite,'(F25.10)') secondary(68)
write(s_aragonite,'(F25.10)') secondary(69)


! OTHER INFORMATION TO STRINGS
write(s_temp,'(F25.10)') temp
write(s_timestep,'(F25.10)') timestep




! ----------------------------------%%
! INITIAL AQUEOUS PHASE CONSITUENTS
! ----------------------------------%%

inputz0 = "SOLUTION 1 " //NEW_LINE('')// &
!&"    pH " // trim(s_pH) //NEW_LINE('')// &
!&"    pe " // trim(s_pe) //NEW_LINE('')// &
&"    units   mol/kgw" //NEW_LINE('')// &
&"    temp" // trim(s_temp) //NEW_LINE('')// &
&"    Ca " // trim(s_ca) //NEW_LINE('')// &
&"    Mg " // trim(s_mg) //NEW_LINE('')// &
&"    Na " // trim(s_na) //NEW_LINE('')// &
&"    K " // trim(s_k) //NEW_LINE('')// &
&"    Fe " // trim(s_fe) //NEW_LINE('')// &
&"    S "// trim(s_s)  //NEW_LINE('')// &
&"    Si " // trim(s_si) //NEW_LINE('')// &
&"    Cl " // trim(s_cl) //NEW_LINE('')// &
&"    Al " // trim(s_al) //NEW_LINE('')// &
&"    C " // trim(s_co2) //NEW_LINE('')// &
!&"    C " // trim(s_hco3) // "as HCO3-" //NEW_LINE('')// &
&"    Alkalinity " // trim(s_alk) //NEW_LINE('')// &
!&"    -water		.0386 # kg" //NEW_LINE('')// &
&"    -water "// trim(s_w) // " # kg" //NEW_LINE('')// &


! inputz0 = "SOLUTION 1 " //NEW_LINE('')// &
! &"    pH 8.22" //NEW_LINE('')// &
! &"    units   ppm" //NEW_LINE('')// &
! &"    temp 10.0" //NEW_LINE('')// &
! &"    Ca 412.3" //NEW_LINE('')// &
! &"    Mg 1291.8" //NEW_LINE('')// &
! &"    Na 10768.0" //NEW_LINE('')// &
! &"    K 399.1" //NEW_LINE('')// &
! &"    Fe .002" //NEW_LINE('')// &
! &"    Mn .0002" //NEW_LINE('')// &
! &"    S(6) 2712.0 as SO4" //NEW_LINE('')// &
! &"    Si 4.28" //NEW_LINE('')// &
! &"    Cl 19353.0" //NEW_LINE('')// &
! &"    NO3 0.29" //NEW_LINE('')// &
! &"    NH4 0.03" //NEW_LINE('')// &
! &"    Alkalinity 141.682 as HCO3" //NEW_LINE('')// &
! &"    -water		5.0	# kg" //NEW_LINE('')// &
! &"END" //NEW_LINE('')// &

! ----------------------------------%%
! HYDROTHERMAL MINERAL CHOICES
! ----------------------------------%%

&"EQUILIBRIUM_PHASES 1" //NEW_LINE('')// &
!&"    CO2(g) -3.25 1000" //NEW_LINE('')// &
&"    Kaolinite 0.0 " // trim(s_kaolinite)  //NEW_LINE('')// &
&"    Goethite 0.0 " // trim(s_goethite) //NEW_LINE('')// &
&"    Celadonite 0.0 " // trim(s_celadonite) //NEW_LINE('')// &
&"    SiO2(am) 0.0 " // trim(s_sio2) //NEW_LINE('')// &
&"    Albite 0.0 " // trim(s_albite) //NEW_LINE('')// &
&"    Calcite 0.0 " // trim(s_calcite) //NEW_LINE('')// &
&"    Montmor-Na 0.0 " // trim(s_mont_na) //NEW_LINE('')// &
&"    Saponite-Mg 0.0 " // trim(s_saponite) //NEW_LINE('')// &
&"    Stilbite 0.0 " // trim(s_stilbite) //NEW_LINE('')// &
&"    Clinoptilolite-Ca 0.0 " // trim(s_clinoptilolite) //NEW_LINE('')// &
&"    Pyrite 0.0 " // trim(s_pyrite) //NEW_LINE('')// &
&"    Quartz 0.0 " // trim(s_quartz) //NEW_LINE('')// &
&"    K-Feldspar 0.0 " // trim(s_kspar) //NEW_LINE('')// &

! NEW MINS

! &"    Dolomite 0.0 " // trim(s_dolomite) //NEW_LINE('')// &
 &"    Saponite-Na 0.0 " // trim(s_saponite_na) //NEW_LINE('')// &
 &"    Nontronite-Na 0.0 " // trim(s_nont_na) //NEW_LINE('')// &
 &"    Nontronite-Mg 0.0 " // trim(s_nont_mg) //NEW_LINE('')// &
 &"    Nontronite-K 0.0 " // trim(s_nont_k) //NEW_LINE('')// &
 &"    Nontronite-H 0.0 " // trim(s_nont_h) //NEW_LINE('')// &
 &"    Nontronite-Ca 0.0 " // trim(s_nont_ca) //NEW_LINE('')// &
 &"    Muscovite 0.0 " // trim(s_muscovite) //NEW_LINE('')// &
! &"    Mesolite 0.0 " // trim(s_mesolite) //NEW_LINE('')// & !!!!!!!!
 &"    Hematite 0.0 " // trim(s_hematite) //NEW_LINE('')// &
 &"    Diaspore 0.0 " // trim(s_diaspore) //NEW_LINE('')// &

! &"    Dawsonite 0.0 " // trim(s_dawsonite) //NEW_LINE('')// &
! &"    Anhydrite 0.0 " // trim(s_magnesite) //NEW_LINE('')// & ! formerly magnesite


!!!! &"    Quartz 0.0 0.0" //NEW_LINE('')// &
 &"    Smectite-high-Fe-Mg 0.0 " // trim(s_smectite) //NEW_LINE('')// &
! &"    Dolomite 0.0 " // trim(s_dolomite) //NEW_LINE('')// &
! &"    Siderite 0.0 " // trim(s_siderite) //NEW_LINE('')// &

!! MINERALS ADDED 10/18/2014
 &"    Akermanite 0.0 " // trim(s_akermanite) //NEW_LINE('')// &
 &"    Analcime 0.0 " // trim(s_analcime) //NEW_LINE('')// &
 &"    Annite 0.0 " // trim(s_annite) //NEW_LINE('')// &
 &"    Clinozoisite 0.0 " // trim(s_clinozoisite) //NEW_LINE('')// &
 &"    Dicalcium_silicate 0.0 " // trim(s_dicalcium_silicate) //NEW_LINE('')// &
 &"    Diopside 0.0 " // trim(s_diopside) //NEW_LINE('')// &
 &"    Epidote  0.0 " // trim(s_epidote ) //NEW_LINE('')// &
 &"    Ettringite 0.0 " // trim(s_ettringite) //NEW_LINE('')// &
 &"    Ferrite-Ca 0.0 " // trim(s_ferrite_ca) //NEW_LINE('')// &
 &"    Foshagite 0.0 " // trim(s_foshagite) //NEW_LINE('')// &
 &"    Gismondine 0.0 " // trim(s_gismondine) //NEW_LINE('')// &
 &"    Gyrolite 0.0 " // trim(s_gyrolite) //NEW_LINE('')// &
 &"    Hedenbergite 0.0 " // trim(s_hedenbergite) //NEW_LINE('')// &
 &"    Hillebrandite 0.0 " // trim(s_hillebrandite) //NEW_LINE('')// &
 &"    Larnite 0.0 " // trim(s_Larnite) //NEW_LINE('')// &
 &"    Laumontite 0.0 " // trim(s_laumontite) //NEW_LINE('')// &
 &"    Lawsonite 0.0 " // trim(s_lawsonite) //NEW_LINE('')// &
 &"    Merwinite 0.0 " // trim(s_merwinite) //NEW_LINE('')// &
 &"    Monticellite 0.0 " // trim(s_monticellite) //NEW_LINE('')// &
 &"    Natrolite 0.0 " // trim(s_natrolite) //NEW_LINE('')// &
 &"    Okenite 0.0 " // trim(s_okenite) //NEW_LINE('')// &
 &"    Phlogopite 0.0 " // trim(s_phlogopite) //NEW_LINE('')// &
 &"    Prehnite 0.0 " // trim(s_prehnite) //NEW_LINE('')// &
 &"    Pseudowollastonite 0.0 " // trim(s_pseudowollastonite) //NEW_LINE('')// &
 &"    Rankinite 0.0 " // trim(s_rankinite) //NEW_LINE('')// &
 &"    Scolecite 0.0 " // trim(s_scolecite) //NEW_LINE('')// &
 &"    Tobermorite-9A 0.0 " // trim(s_tobermorite_9a) //NEW_LINE('')// &
 &"    Tremolite 0.0 " // trim(s_tremolite) //NEW_LINE('')// &
 &"    Wollastonite 0.0 " // trim(s_wollastonite) //NEW_LINE('')// &
 &"    Xonotlite 0.0 " // trim(s_xonotlite) //NEW_LINE('')// &
 &"    Zoisite 0.0 " // trim(s_zoisite) //NEW_LINE('')// &
 !! NEXT ROUND
 &"   Andradite 0.0 " // trim(s_andradite) //NEW_LINE('')// &
 &"   Troilite 0.0 " // trim(s_troilite) //NEW_LINE('')// &
 &"   Pyrrhotite 0.0 " // trim(s_pyrrhotite) //NEW_LINE('')// &
 &"   Minnesotaite 0.0 " // trim(s_minnesotaite) //NEW_LINE('')// &
 &"   Fayalite 0.0 " // trim(s_fayalite) //NEW_LINE('')// &
 &"   Daphnite-7a 0.0 " // trim(s_daphnite_7a) //NEW_LINE('')// &
 &"   Daphnite-14a 0.0 " // trim(s_daphnite_14a) //NEW_LINE('')// &
 &"   Cronstedtite-7a 0.0 " // trim(s_cronstedtite_7a) //NEW_LINE('')// &
! &"   Greenalite 0.0 " // trim(s_greenalite) //NEW_LINE('')// &
! &"   Aragonite 0.0 " // trim(s_aragonite) //NEW_LINE('')// &
! no dolomite, aragonite
!&"SAVE solution 1"  //NEW_LINE('')// &
!&"SAVE equilibrium_phases 1"  //NEW_LINE('')// &
!&"END"  //NEW_LINE('')// &


! ----------------------------------%%
! KINETIC DISSOLUTION RATE LAWS
! ----------------------------------%%
	
&"RATES" //NEW_LINE('')// &

&"BGlass" //NEW_LINE('')// &
&"-start" //NEW_LINE('')// &
!CALC_VALUE('R(s_sp)')
&"    10 rate0=M*46.5*1.53e-5*0.01*(1e4)*(2.51189e-6)*exp(-25.5/(.008314*TK))" // &
&"*(((ACT('H+')^3)/(ACT('Al+3')))^.333)" //NEW_LINE('')// &
&"    20 save rate0 * time" //NEW_LINE('')// &
&"-end" //NEW_LINE('')// &

&"Plagioclase" //NEW_LINE('')// &
&"-start" //NEW_LINE('')// &
!(1-SR('Plagioclase'))*
&"    10 rate = (1-SR('Plagioclase'))*M*270.0*(1.53e-5)*0.01*(((1.58e-9)"//&
&"*exp(-53.5/(.008314*TK))*(ACT('H+')^0.541) +(3.39e-12)*exp(-57.4/(.008314*TK)) +"//&
&"(4.78e-15)*exp(-59.0/(.008314*TK))*(ACT('H+'))^-0.57))"//NEW_LINE('')//&
&"    20 save rate * time"//NEW_LINE('')//&
&"-end" //NEW_LINE('')// &

&"Augite" //NEW_LINE('')// &
&"-start" //NEW_LINE('')// &
&"    10 rate0 = (1-SR('Augite'))*M*230.0*(1.53e-5)*0.01*(((1.58e-7)" // &
&"*exp(-78.0/(.008314*TK))*(ACT('H+')^0.7)+(1.07e-12)*exp(-78.0/(.008314*TK))))" //NEW_LINE('')// & 
&"    20 save rate0 * time" //NEW_LINE('')// &
&"-end" //NEW_LINE('')// &

&"Pigeonite" //NEW_LINE('')// &
&"-start" //NEW_LINE('')// &
&"    10 rate0 = (1-SR('Pigeonite'))*M*236.0*(1.53e-5)*0.01*(((1.58e-7)" // &
&"*exp(-78.0/(.008314*TK))*(ACT('H+')^0.7)+(1.07e-12)*exp(-78.0/(.008314*TK))))"//NEW_LINE('')// &
&"    20 save rate0 * time" //NEW_LINE('')// &
&"-end" //NEW_LINE('')// &

&"Magnetite" //NEW_LINE('')// &
&"-start" //NEW_LINE('')// &
&"    10 rate0 = (1-SR('Magnetite'))*M*231.0*(1.53e-5)*0.01*(((2.57e-9)" // &
&"*exp(-18.6/(.008314*TK))*(ACT('H+')^0.279)+(1.66e-11)*exp(-18.6/(.008314*TK))))" //NEW_LINE('')// &
&"    20 save rate0 * time" //NEW_LINE('')// &
&"-end" //NEW_LINE('')// &

! ----------------------------------%%
! PRIMARY (KINETIC) CONSTITUENTS
! ----------------------------------%%

!&"Use solution 1" //NEW_LINE('')// &
!&"Use equilibrium_phases 1" //NEW_LINE('')// &
&"KINETICS 1" //NEW_LINE('')// &
&"Plagioclase" //NEW_LINE('')// &
&"-m0 " // trim(s_feldspar) //NEW_LINE('')// &
&"Augite" //NEW_LINE('')// &
&"-m0 " // trim(s_augite) //NEW_LINE('')// &
&"Pigeonite" //NEW_LINE('')// &
&"-m0 " // trim(s_pigeonite) //NEW_LINE('')// &
&"Magnetite" //NEW_LINE('')// &
&"-m0 " // trim(s_magnetite) //NEW_LINE('')// &
&"BGlass" //NEW_LINE('')// &

! !pham phudge
! &"-f CaO 0.025 Fe2O3 0.0475 MgO 0.065 " //&
! & "Na2O 0.0125 K2O 0.005 Al2O3 0.0525 SiO2 0.5" //NEW_LINE('')// &

! ! pham
! &"-f Ca 0.015 Fe 0.0095 Mg 0.065 " //&
! & "Na 0.025 K 0.001 Al 0.105 Si 0.5 O 1.35" //NEW_LINE('')// &

! ! pham
! &"-f Ca 0.015 Fe 0.0095 Mg 0.065 " //&
! &"Na 0.025 K 0.001 Al 0.105 " //&
! & "S 0.003 Si 0.5 O 1.35" //NEW_LINE('')// &

! ! grove and kinzler 1992 (thomspon et al 1980)
&"-f CaO 0.182 SiO2 0.861 Al2O3 0.16 " //&
& "FeO .121 MgO 0.195 K2O 0.00265 " //&
& "Na2O 0.0573" //NEW_LINE('')// &
!
! Al2O3



! grove and kinzler 1992 (thomspon et al 1980)
! &"-f Ca 0.182 Si 0.861 Al 0.32 " //&
! & "Fe 0.121 Mg 0.195 K 0.0053 " //&
! & "Na 0.1146 O 2.88095" //NEW_LINE('')// &
!2.88095
!2.75995
!2.63895
&"-m0 " // trim(s_glass) //NEW_LINE('')// &

&"    -step " // trim(s_timestep) // " in 1" //NEW_LINE('')// &

&"INCREMENTAL_REACTIONS true" //NEW_LINE('')// &


    
!
! ! ----------------------------------%%
! ! CALCULATE POROSITY AND STUFF
! ! ----------------------------------%%
!
! &"CALCULATE_VALUES" //NEW_LINE('')// &
!
! &"R(sum)" //NEW_LINE('')// &
! &"-start" //NEW_LINE('')// &
! ! &"10 sum = EQUI('Stilbite')*832.2/2.15 + EQUI('SiO2(am)')*60.0/2.62" //&
! ! &"+ EQUI('Kaolinite')*258.2/2.6 + EQUI('Albite')*262.3/2.62" // &
! ! &"+ EQUI('Saponite-Mg')*385.537/2.4 + EQUI('Celadonite')*396.8/3.0" // &
! ! &"+ EQUI('Clinoptilolite-Ca')*1344.49/2.62 + EQUI('Pyrite')*120.0/4.84" // &
! ! &"+ EQUI('Montmor-Na')*103.8/5.3 + EQUI('Goethite')*88.8/3.8" // &
! ! &"+ EQUI('Dolomite')*184.3/2.84 + EQUI('Smectite-high-Fe-Mg')*425.7/2.7" // &
! ! &"+ EQUI('Dawsonite')*144.0/2.42 + EQUI('Anhydrite')*84.3/3.0" // &
! ! &"+ EQUI('Siderite')*115.8/3.96 + EQUI('Calcite')*100.0/2.71" // &
! ! &"+ EQUI('Quartz')*60.0/2.62 + EQUI('K-Feldspar')*193.0/2.56" // &
! ! &"+ KIN('Plagioclase')*270.0/2.68 + KIN('Augite')*230.0/3.4" // &
! ! &"+ KIN('Pigeonite')*239.6/3.38 + KIN('Magnetite')*231.0/5.15" // &
! ! &"+ KIN('BGlass')*46.5/2.92" // &
!  &"10 sum = 5.0" //&
! &"" //NEW_LINE('')// &
! &"100 SAVE sum" //NEW_LINE('')// &
! &"-end" //NEW_LINE('')// &
!
! &"R(phi)" //NEW_LINE('')// &
! &"-start" //NEW_LINE('')// &
! !&"10 phi = 1.0-(CALC_VALUE('R(sum)')/(CALC_VALUE('R(sum)')+(TOT('water')*1000.0)))" //&
! &"10 phi = 0.1" //&
! &"" //NEW_LINE('')// &
! &"100 SAVE phi" //NEW_LINE('')// &
! &"-end" //NEW_LINE('')// &
!
! &"R(water_volume)" //NEW_LINE('')// &
! &"-start" //NEW_LINE('')// &
! !&"10 water_volume = SOLN_VOL / RHO" //NEW_LINE('')// &
! &"10 water_volume = 0.3" //NEW_LINE('')// &
! &"100 SAVE water_volume" //NEW_LINE('')// &
! &"-end" //NEW_LINE('')// &
!
!
! &"R(rho_s)" //NEW_LINE('')// &
! &"-start" //NEW_LINE('')// &
! ! &"10 rho_s = EQUI('Stilbite')*2.15 + EQUI('SiO2(am)')*2.62" //&
! ! &"+ EQUI('Kaolinite')*2.6 + EQUI('Albite')*2.62" // &
! ! &"+ EQUI('Saponite-Mg')*2.4 + EQUI('Celadonite')*3.0" // &
! ! &"+ EQUI('Clinoptilolite-Ca')*2.62 + EQUI('Pyrite')*4.84" // &
! ! &"+ EQUI('Montmor-Na')*5.3 + EQUI('Goethite')*3.8" // &
! ! &"+ EQUI('Dolomite')*2.84 + EQUI('Smectite-high-Fe-Mg')*2.7" // &
! ! &"+ EQUI('Dawsonite')*2.42 + EQUI('Anhydrite')*3.0" // &
! ! &"+ EQUI('Siderite')*3.96 + EQUI('Calcite')*2.71" // &
! ! &"+ EQUI('Quartz')*2.62 + EQUI('k-Feldspar')*2.56" // &
! ! &"+ KIN('Plagioclase')*2.68 + KIN('Augite')*3.4" // &
! ! &"+ KIN('Pigeonite')*3.38 + KIN('Magnetite')*5.15" // &
! ! &"+ KIN('BGlass')*2.92" //NEW_LINE('')// &
! ! &"20 rho_s = rho_s/ (EQUI('Stilbite') + EQUI('SiO2(am)')" //&
! ! &"+ EQUI('Kaolinite') + EQUI('Albite')" // &
! ! &"+ EQUI('Saponite-Mg') + EQUI('Celadonite')" // &
! ! &"+ EQUI('Clinoptilolite-Ca') + EQUI('Pyrite')" // &
! ! &"+ EQUI('Montmor-Na') + EQUI('Goethite')" // &
! ! &"+ EQUI('Dolomite') + EQUI('Smectite-high-Fe-Mg')" // &
! ! &"+ EQUI('Dawsonite') + EQUI('Anhydrite')" // &
! ! &"+ EQUI('Siderite') + EQUI('Calcite')" // &
! ! &"+ EQUI('Quartz') + EQUI('K-Feldspar')" // &
! ! &"+ KIN('Plagioclase') + KIN('Augite')" // &
! ! &"+ KIN('Pigeonite') + KIN('Magnetite')" // &
! ! &"+ KIN('BGlass'))" //NEW_LINE('')// &
! !&"30 rho_s = rho_s * 1000000.0" //NEW_LINE('')// &
! &"10 rho_s = 2.93e6" //NEW_LINE('')// &
! &"100 SAVE rho_s" //NEW_LINE('')// &
! &"-end" //NEW_LINE('')// &
!
!
! &"R(s_sp)" //NEW_LINE('')// &
! &"-start" //NEW_LINE('')// &
! !&"10 s_sp = (CALC_VALUE('R(phi)')/(1.0-CALC_VALUE('R(phi)')))*400.0/CALC_VALUE('R(rho_s)')" //&
! &"10 s_sp = 1.53e-5" //&
! &"" //NEW_LINE('')// &
! &"100 SAVE s_sp" //NEW_LINE('')// &
! &"-end" //NEW_LINE('')// &
!
! !&"SAVE solution 1"  //NEW_LINE('')// &
! !&"SAVE equilibrium_phases 1"  //NEW_LINE('')// &
!
! !&"END" //NEW_LINE('')// &


! ----------------------------------%%
! KNOBS FOR SOLVERS AND STUFF
! ----------------------------------%%
	
	
!&"END" //NEW_LINE('')// &
! &"KNOBS" //NEW_LINE('')// &
! &"-diagonal_scale true" //NEW_LINE('')// &


! ----------------------------------%%
! DEFINE THE KIND OF OUTPUT
! ----------------------------------%%

!&"DUMP" //NEW_LINE('')// &
!&"    -solution 1" //NEW_LINE('')// &
!&"    -equilibrium_phases" //NEW_LINE('')// &

  &"SELECTED_OUTPUT" //NEW_LINE('')// &
  &"    -reset false" //NEW_LINE('')// &
!  &"    -high_precision true" //NEW_LINE('')// &
  &"    -k plagioclase augite pigeonite magnetite bglass" //NEW_LINE('')// &
  &"    -ph" //NEW_LINE('')// &
  &"    -pe false" //NEW_LINE('')// &
!  &"    -molalities Ca+2 Mg+2 Na+ K+ Fe+3 SO4-2 SiO2 Cl- Al+3 HCO3-" //NEW_LINE('')// &
  &"    -totals C" //NEW_LINE('')// &
  &"    -totals Ca Mg Na K Fe S Si Cl Al " //NEW_LINE('')// &
  &"    -molalities HCO3-" //NEW_LINE('')// &
  &"    -water true" //NEW_LINE('')// &
  &"    -alkalinity" //NEW_LINE('')// &
  &"    -p stilbite sio2(am) kaolinite albite saponite-mg celadonite Clinoptilolite-Ca" //NEW_LINE('')// &
  &"    -p pyrite Montmor-Na goethite dolomite Smectite-high-Fe-Mg Dawsonite" //NEW_LINE('')// &
  &"    -p Anhydrite siderite calcite quartz k-feldspar" //NEW_LINE('')// &
  &"    -p saponite-na Nontronite-Na Nontronite-Mg Nontronite-K Nontronite-H " //NEW_LINE('')// &
  &"    -p Nontronite-Ca muscovite mesolite hematite diaspore" //NEW_LINE('')// &
  ! NEW MINERALS ADDED 10/18/2014
  &"    -p Akermanite Analcime Annite Clinozoisite Dicalcium_silicate" //NEW_LINE('')// &
  &"    -p Diopside Epidote Ettringite Ferrite-Ca Foshagite Gismondine" //NEW_LINE('')// &
  &"    -p Gyrolite Hedenbergite Hillebrandite Larnite Laumontite" //NEW_LINE('')// &
  &"    -p Lawsonite Merwinite Monticellite Natrolite Okenite" //NEW_LINE('')// &
  &"    -p Phlogopite Prehnite Pseudowollastonite Rankinite Scolecite" //NEW_LINE('')// &
  &"    -p Tobermorite-9A Tremolite Wollastonite Xonotlite  Zoisite" //NEW_LINE('')// &
  ! NEXT ROUND
  &"    -p Andradite Troilite Pyrrhotite Minnesotaite Fayalite Daphnite-7A" //NEW_LINE('')// &
  &"    -p Daphnite-14A Cronstedtite-7A Greenalite Aragonite" //NEW_LINE('')// &
  &"    -p Siderite Magnesite" //NEW_LINE('')// &
  !&"    -calculate_values R(phi) R(s_sp) R(water_volume) R(rho_s)" //NEW_LINE('')// &
  &"    -time" //NEW_LINE('')// &
&"END"
  
! write(*,*) "primary"
! write(*,*) primary
! write(*,*) "secondary"
! write(*,*) secondary
! write(*,*) "solute"
! write(*,*) solute
! write(*,*) "medium"
! write(*,*) medium
! write(*,*) "temp"
! write(*,*) temp
  
  
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
	write(*,*) "primary"
	write(*,*) primary
	write(*,*) "secondary"
	write(*,*) secondary
	write(*,*) "solute"
	write(*,*) solute
	write(*,*) "medium"
	write(*,*) medium
	write(*,*) "temp"
	write(*,*) temp
	!STOP
END IF



! ! !if errror print dump. look for other phases.
! if ((temp .GT. 35.0) .and. (temp .LT. 40.0)) then
! 	DO i=1,GetOutputStringLineCount(id)
! 		call GetOutputStringLine(id, i, line)
! 		write(*,*) trim(line)
! 	END DO
! end if
  
! PRINT DUMP/OUTPUT
DO i=1,GetOutputStringLineCount(id)
	call GetOutputStringLine(id, i, line)
 	if ((medium(6) .gt. 2800.0) .and. (medium(7) .lt. -1100.00)) then
 		!write(*,*) trim(line)
 	end if
	!write(*,*) trim(line)
END DO
  
! NOW KINDA USELESS PRINT STATEMENTS FOR WRITING TO FILES
!WRITE(*,*) "WRITING TO 2D ARRAY AND OUTPUT FILES"
!WRITE(*,*) "NUMBER OF LINES:"
!WRITE(*,*) GetSelectedOutputStringLineCount(id)
  
  
! OPEN FILE (don't need no file) (USEFUL)
!OPEN(UNIT=12, FILE="testMat.txt", ACTION="write", STATUS="replace") 
  
! WRITE AWAY
allocate(outmat(GetSelectedOutputStringLineCount(id)+1,all))
!write(*,*) size(outmat,1)
!write(*,*) size(outmat,2)
DO i=1,GetSelectedOutputStringLineCount(id)
	call GetSelectedOutputStringLine(id, i, line)
	! HEADER BITS YOU MAY WANT
	if (i .eq. 1) then
 	   !write(12,*) trim(line)
	   !write(*,*) "cell"
	   !write(*,*) trim(line) ! PRINT LABELS FOR EVERY FIELD (USEFUL)
		
		if ((medium(6) .gt. 2800.0) .and. (medium(7) .lt. -1100.00)) then
			!write(*,*) trim(line)
		end if
		
	end if
	! MEAT
	if (i .gt. 1) then
		if ((medium(6) .gt. 2800.0) .and. (medium(7) .lt. -1100.00)) then
			!write(*,*) trim(line)
			!write(*,*) " "
			!write(*,*) "line!"
		end if
		
		read(line,*) outmat(i,:)
		
		
		!write(12,"(4F12.5)") outmat(i,:)
		!!!!!write(12,*) outmat(i,:) ! this writes to file, which i don't need (USEFUL)
		!write(*,*) trim(line) ! PRINT EVERY GOD DAMN LINE
	end if
END DO
  

! OUTPUT TO THE MAIN MASSACR METHOD
alter(1,:) = outmat(3,:)
!write(*,*) outmat
write(*,*) alter(1,4)

IF (RunString(id, inputz0).NE.0) THEN
	alter(1,:) = 0.0
END IF

! DESTROY INSTANCE
IF (DestroyIPhreeqc(id).NE.IPQ_OK) THEN
	STOP
END IF

! ALL DONE!
!write(*,*) "phreeqed out"

!write(*,*) outmat(2,54)


return
  
end function alter





end module alteration
