      module radtran

      USE define_csu1dvar
      USE define_intam2
      !USE Fastem
      USE miemod

      implicit  none

      contains

      subroutine run_eddington(tempav,mixrat,clwc,swp,rwp,tskin,
     + emu,edm,sfc_wind,tbout)
      
      ! This subroutine uses the Eddington07 forward radiative transfer model
      ! to calculate the upwelling TOA TBs for a given atmosphere using GMI
      ! frequencies and incidence angles
      ! Inputs:
      ! 	tempav(1:nz) - average of temperature of layer (K)
      !		mixrat(1:nz) - average water vapor mixing ratio of each layer (g/kg)
      !		clwc(1:nz) - cloud liquid water content (kg/m^2/layer)
      ! 	ciwc(1:nz) - cloud ice water content - zero for now!
      ! 	tskin - skin temperature (K)
      !		sfc_wind - surface wind speed (m/s)
      !		emis_list(1:nchan) - emissivities for each GMI channel - added 6/9
      !		ang(1:2) - low(1) and high(2) frequency view angles - added 6/13
      ! Output:
      !		tbout(1:nchan) - upwelling TOA brightness temps for each sensor channel
      
      ! include   'parameters.inc'
c    
c     Input atmospheric model variables
c 
      ! 16 layers, ocean only (for now)
      integer, parameter	::	sfc_type = 0	! land = sfc_type 1
      !integer, parameter :: nz=16 ! same as retrieval
      real	view_angle, azz
      integer	pol(nfreq)
      real      tskin, sfc_wind, emis, refl
      real      emu,edm!,enw
      real      hgt_lev(0:nz), press_lev(0:nz), temp_lev(0:nz)
      real      vapor_pressure(nz),rwp(nz), cloud_water(nz)
      real      rain_water(nz), cloud_ice(nz), snow(nz)
      real	tempav(nz), mixrat(nz), clwc(nz), swp(nz)
      real      freq(nfreq)
      real	lapse		! added for interpolating to temp_lev
      real 	fullemis(4), fullrefl(4) 	! Added for dealing with FM6 emissivities
      real	trans, optd		! Added for calculating diffuse reflection
		
c      
c     Internally assigned/computed cloud parameters
c
      real      T_avg, P_avg, atm_ext
      real      kext_clw, salb_clw, asym_clw
      real      kext_rain, salb_rain, asym_rain
      real      kext_ci, salb_ci, asym_ci
      real      kext_snow, salb_snow, asym_snow
      real      kexttot(nz), salbtot(nz), asymtot(nz)
      logical   lambert
      real      ebar
      real      Tb, Tb_out(nfreq,0:1)
      real	tbout(nch)
      real :: reff,tau_est(nz)
      real :: rainr

c
c     Loop counters 
c
      integer    nl, nf, npol, i, k, nc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      ! Set freqs and pols
      !if (satcode.eq.'TMI' .and. nch.eq.9) then
      !  freq = (/10.65,19.35,21.3,37.0,85.5/)
      !  pol = (/2,2,1,2,2/)
      !endif
      !if (satcode.eq.'GMI' .and. nch.eq.9) then
      !  freq = (/10.65,18.7,23.8,36.64,89.0/)
!	pol = (/2,2,1,2,2/)
!      endif
      if (satcode.eq.'GMI' .and. nch.eq.13) then
        freq = (/10.65,18.7,23.8,36.64,89.0,166.0,186.31,190.31/)
	pol = (/2,2,1,2,2,2,1,1/)
      endif

c     ! Set up atmosphere
      !pecz(:) = lz(:)
      !pplev(:)= lpress(:)
      !pplev(nz+1) = 997.2 ! SLP
      do k = 1, nz
         rain_water(k) = rwp(nz+1-k)*1000.0/dz(nz+1-k)	! convert to g/m^3
         cloud_ice(k) = 0
         snow(k) = swp(nz+1-k)*1000.0/dz(nz+1-k)	! convert to g/m^3
         !snow(k) = ciwc(k) !0
         press_lev(k) = pplev(nz+1-k)			! flip so level 0 is surface
         hgt_lev(k) = pecz(nz+1-k)/1000.0			! convert to km
         cloud_water(k) = clwc(nz+1-k) * 1000.0 / dz(nz+1-k)	! convert to g/m^3
      end do

      press_lev(0) = pplev(nz+1)
      hgt_lev(0) = pecz(nz+1)/1000.0				! km
      ! temps at bottom and top of each layer needed for eddington - 
      ! linear interpolation used for now - is there a better way to do this?
      lapse = (tempav(15)-tempav(16)) / (0.5*(dz(16)+dz(15)))	! K/m
      temp_lev(0) = tempav(16) - lapse*0.5*dz(16)
      do k = 1, nz-1
         lapse=(tempav(16-k)-tempav(nz+1-k))/(0.5*(dz(nz+1-k)+dz(16-k)))
         temp_lev(k) = tempav(nz+1-k) + lapse*0.5*dz(nz+1-k)
      end do
      temp_lev(nz) = tempav(1) + lapse*0.5*dz(1)
      
      
!      write(*,*) ' height ',hgt_lev
!      write(*,*) ' dz ',dz
!      write(*,*) ' pecz ',pecz
!      write(*,*) ' pressure ',press_lev
!      write(*,*) ' temperature ',temp_lev
!      write(*,*) ' rwc ',rain_water
!      write(*,*) ' cwc ',cloud_water

      nc = 1	! counter
      lambert = .false.
c
c     Loop over each frequency          
c
      saver = 0.0 !initialize
      do nf = 1, nfreq
        if (nf.gt.5) then
	  view_angle = eia(oepix,oelin,14) 	! 183 eia
	else if (satcode.eq.'GMI') then
	  view_angle = eia(oepix,oelin,1)	! 10V eia
	!else	! TMI has eia differences
	!  view_angle = sat_eia(oepix,oelin,3)	! 19V, default
	!  if (nc.eq.1) view_angle = sat_eia(oepix,oelin,1)
	!  if (nc.eq.2) view_angle = sat_eia(oepix,oelin,2)
	!  if (nc.eq.6) view_angle = sat_eia(oepix,oelin,7) !37V eia
	endif
	optd = 0.0	! Reset optical depth
        do nl = 1, nz  
          P_avg = (press_lev(nl) - press_lev(nl-1))/
     +                  log(press_lev(nl)/press_lev(nl-1))
	  !P_avg = (press_lev(nl) + press_lev(nl-1))/2.0
          !write(*,*),'Pavg: ',P_avg
          T_avg = tempav(nz+1-nl)
          
          vapor_pressure(nl) = (mixrat(nz+1-nl) * 0.001 * P_avg) /
     +					  (0.622 + mixrat(nz+1-nl) * 0.001)    

          call absorb_clr(freq(nf), T_avg, P_avg, vapor_pressure(nl), 
     +                    atm_ext)
c
          call mie_clw(freq(nf), T_avg, cloud_water(nl), 
     +                   kext_clw, salb_clw, asym_clw)
c

! emu,edm,enw should be included in new call
          call mie_rain(dble(freq(nf)), dble(T_avg),
     >     dble(rain_water(nl)),dble(edm),dble(emu),
     >    kext_rain, salb_rain, asym_rain,rainr)
        !print*,'mie rain ext: ',kext_rain,nl,nf

        if(nf==nfreq) then
        !print*,rainr,rain_water(nl),nl,saver
          if(pplev(nz+1)>=990) then
            if(nl<=3.and.rain_water(nl)>0.0) saver=saver+rainr !btm 3 lyrs
          endif
          if(pplev(nz+1)<990.and.pplev(nz+1)>960) then
            if(nl<=4.and.nl>1.and.rain_water(nl)>0.0) saver=saver+rainr
          endif
          if(pplev(nz+1)<=960) then
            if(nl<=5.and.nl>2.and.rain_water(nl)>0.0) saver=saver+rainr
          endif
        endif ! helps avoid weird RWC/dz problems if screen for lyr depth
        ! average last 3 vertical layers together -- no freq dependence!
c
          call mie_ci(freq(nf), T_avg, cloud_ice(nl), 
     +                   kext_ci, salb_ci, asym_ci)
c
          call mie_snow(freq(nf), T_avg, snow(nl), 
     +                   kext_snow, salb_snow, asym_snow)
c
c
c         Sum the scattering parameters for each species 
c
          kexttot(nl) = atm_ext + kext_clw + kext_rain + kext_ci + 
     +                      kext_snow 
c     
          if ( kexttot(nl) .gt. 0. ) then
             salbtot(nl) = (salb_clw*kext_clw + salb_rain*kext_rain +
     +                      salb_ci*kext_ci   + salb_snow*kext_snow )/ 
     +                      kexttot(nl)
          else
             salbtot(nl) = 0.
          endif
c              
          if ( salbtot(nl) .gt. 0. ) then
            asymtot(nl) = ( asym_clw*salb_clw*kext_clw +
     +                       asym_rain*salb_rain*kext_rain +
     +                       asym_ci*salb_ci*kext_ci +
     +                       asym_snow*salb_snow*kext_snow ) /
     +                       (salbtot(nl)*kexttot(nl))
          else
            asymtot(nl) = 0.
          endif 
	  ! Add optical depth of layer to the sum
	  optd = optd + kexttot(nl) * 
     +		(hgt_lev(nl)-hgt_lev(nl-1)) /
     +		cos(view_angle*0.0174533)
          tau_est(nl) = kexttot(nl)*
     +      (hgt_lev(nl)-hgt_lev(nl-1))!/cos(view_angle*0.0174533) !optd
        end do        ! end loop over nl
	
	
	! Find emissivities - Fastem6 model used as default
        psss = 35.0 ! not important now
        azz = 0.0 ! whatever, not important
!	call compute_fastem(6,freq(nf),view_angle,tskin,psss,
!     +  		    sfc_wind,fullemis,fullrefl,
!     +			    azz,trans)
!+			    sataz(oepix,oelin),trans)
     

        if ((pol(nf) .eq. 1) .or. (pol(nf) .eq. 2)) then  	! vertical polarization
          call emit(freq(nf), 1, tskin, sfc_wind, view_angle,
     +					emis, ebar)			! To get ebar - is there a better way?
          !emis = fullemis(1)
          emis = save_emis(nc)
	  !refl = fullrefl(1)
	  refl = 1-emis
        !print*,'emis ',emis
         !print*,'kext_ed ',nc,kexttot(:)
         !print*,'tau_edd ',nc,tau_est(:)
         !print*,'tau_crt ',nc,save_tau(nc,:)
        
     	  call eddington(nz, tskin, hgt_lev, temp_lev, kexttot,
     +			salbtot, asymtot, emis, ebar, lambert,
     +			view_angle, freq(nf), sfc_wind, 1, tb,
     +			refl)
c                
c          Verify that Tb are in bounds (including NaN)
c                                
          if ( (Tb .gt.  50.) .and. (Tb .lt. 350.) ) then
              tbout(nc) = Tb
          else
            write(*,*) 'Tb = ',Tb, oepix,oelin
            write(*,*) ' Tb is outside physical range'
            print*,kexttot,freq(nf)
            stop
          endif 
          nc = nc+1
        end if
          
        if ((pol(nf) .eq. 0) .or. (pol(nf) .eq. 2)) then  	! horizontal polarization
          call emit(freq(nf), 0, tskin, sfc_wind, view_angle,
     +					emis, ebar)			
          !emis = fullemis(2)
          emis = save_emis(nc)
	  !refl = fullrefl(2)
	  refl = 1-emis
        !print*,'emis ',emis
         !print*,'kext_ed ',nc,kexttot(:)
         !print*,'tau_edd ',nc,tau_est(:)
         !print*,'tau_crt ',nc,save_tau(nc,:)
     	  call eddington(nz, tskin, hgt_lev, temp_lev, kexttot,
     +			salbtot, asymtot, emis, ebar, lambert,
     +			view_angle, freq(nf), sfc_wind, 2, tb,
     +			refl)                    
          if ( (Tb .gt.  50.) .and. (Tb .lt. 350.) ) then
            tbout(nc) = Tb
          else
            write(*,*) 'Tb = ',Tb,oepix,oelin
            write(*,*) ' Tb is outside physical range'
            stop
          endif 
          nc = nc+1
        end if

      end do  ! end loop over nfreq
      !print*,'saver edd',saver,rainr
      !print*,'**',rainr,rain_water(nz)
      test_out = tbout

      end subroutine run_eddington
      
      ! The remainder of this file includes all of the subroutines (in alphabetical order)
      ! needed for run_eddington. They are taken as-is from the Eddington07 folder. It 
      ! would be good to separate them if I can figure out how to make them all link up
      ! and work within the context of the OE - Rick 
      
!------------------------------------------------------------------
      subroutine absorb_clr(freqy, temp, pres, vapor_pres, kabs_clear)

c     Phil Rosenkranz's absorption code provided to C. kummerow in 2003.
c     Subroutine to provide calling sequence to water vapor, oxygen
c     and air absorption as computed by Phil Rosenkrantz.  The individual
c     absorption codes have not been formally published but rather
c     constitute ongoing work. 
c
c     INPUT
c       freqy      :  frequency [GHz]
c       temp       :  layer avg. temperature [K]
c       pres       :  layer avg. pressure [mb]
c       vapor_pres :  layer average water vapor pressure [mb]
c     OUTPUT
c       kabs       :  absorption coefficient [km^-1]
c
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      

      real  kabs_H2O, kabs_O2, kabs_N2, kabs_clear
      real  freqy, rho, temp, pres, vapor_pres
      
c     Compute the vapor density, rho [g/m^3]      
      rho = vapor_pres*100*18/(8.314*temp)  

      call abs_H2O(temp, pres, rho, freqy, kabs_H2O)      
      call abs_O2 (temp, pres, rho, freqy, kabs_O2) 
      call abs_N2 (temp, pres, freqy, kabs_N2)
      
      kabs_clear = kabs_H2O + kabs_O2 + kabs_N2
      
      return
      end subroutine absorb_clr
      
      Subroutine abs_H2O(T,P,RHO,F,ABH2O)
C
C     C. Kummerow, 8/2003.  Changed function to subroutine     
C     Copyright (c) 2002 Massachusetts Institute of Technology
C
C     NAME- ABH2O    LANGUAGE- FORTRAN 77
C
C     PURPOSE- COMPUTE ABSORPTION COEF IN ATMOSPHERE DUE TO WATER VAPOR
C 
      IMPLICIT NONE
C     CALLING SEQUENCE PARAMETERS-
C     SPECIFICATIONS
      REAL T,P,RHO,F,ABH2O
C      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
C      T       KELVIN    I   TEMPERATURE
C      P       MILLIBAR  I   PRESSURE              .1 TO 1000
C      RHO     G/M**3    I   WATER VAPOR DENSITY
C      F       GHZ       I   FREQUENCY             0 TO 800
C      ABH2O   NEPERS/KM O   ABSORPTION COEFFICIENT
C
C   REFERENCES-
C   P.W. ROSENKRANZ, RADIO SCIENCE V.33, PP.919-928 (1998); V.34, P.1025 (1999).
C
C   LINE INTENSITIES SELECTION THRESHOLD=
C     HALF OF CONTINUUM ABSORPTION AT 1000 MB.
C   WIDTHS MEASURED AT 22,183,380 GHZ, OTHERS CALCULATED.
C     A.BAUER ET AL.ASA WORKSHOP (SEPT. 1989) (380GHz).
c     M. TRETYAKOV et al., J. MOLEC. SPEC. (2003)
C
C   REVISION HISTORY-
C    DATE- OCT.6, 1988  P.W.ROSENKRANZ - EQS AS PUBL. IN 1993.
C          OCT.4, 1995  PWR- USE CLOUGH'S DEFINITION OF LOCAL LINE
C                   CONTRIBUTION,  HITRAN INTENSITIES, ADD 7 LINES.
C          OCT. 24, 95  PWR -ADD 1 LINE.
C          JULY 7, 97   PWR -SEPARATE COEFF. FOR SELF-BROADENING, 
C                       REVISED CONTINUUM.
C        Aug. 28, 2002  PWR - CORRECTED LINE INTENSITIES
C        Mar. 2, 2003   PWR - LINE SHIFT
C
C   LOCAL VARIABLES:
      INTEGER NLINES,I,J
      PARAMETER (NLINES=15)
      REAL DF(2),S1(NLINES),B2(NLINES),W3(NLINES),FL(NLINES),X(NLINES),
     & WS(NLINES),XS(NLINES),SR(NLINES)
      REAL PVAP,PDA,DEN,TI,TI2,SUM,WIDTH,WSQ,S,BASE,RES,CON,SHIFT
C     LINE FREQUENCIES:
      DATA FL/22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508,
     & 443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360,
     & 620.7008, 752.0332, 916.1712/
C     LINE INTENSITIES AT 300K:
      DATA S1/ .1314E-13, .2279E-11, .8058E-13, .2701E-11, .2444E-10,
     & .2185E-11, .4637E-12, .2568E-10, .8392E-12, .3272E-11, .6676E-12,
     & .1535E-08, .1711E-10, .1014E-08, .4238E-10/
C     T COEFF. OF INTENSITIES:
      DATA B2/ 2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405,
     & 3.597, 2.379, 2.852, .159, 2.391, .396, 1.441/
C     AIR-BROADENED WIDTH PARAMETERS AT 300K:
      DATA W3/.00281, .00287, .0023, .00278, .00287, .0021, .00186,
     & .00263, .00215, .00236, .0026, .00321, .00244, .00306, .00267/
C     T-EXPONENT OF AIR-BROADENING:
      DATA X/.69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69,
     & .71, .68, .70/
C     SELF-BROADENED WIDTH PARAMETERS AT 300K:
      DATA WS/.01349, .01491, .0108, .0135, .01541, .0090, .00788,
     & .01275, .00983, .01095, .01313, .01320, .01140, .01253, .01275/
C     T-EXPONENT OF SELF-BROADENING:
      DATA XS/ .61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72,
     & 1.0, .68, .84, .78/
C     RATIO OF SHIFT TO WIDTH
      DATA SR/ 0., -.017, 13*0./
C
      IF(RHO.LE.0.) THEN
        ABH2O = 0.
        RETURN
      ENDIF
      PVAP = RHO*T/217.
      PDA = P -PVAP
      DEN = 3.335E16*RHO ! const includes isotopic abundance
      TI = 300./T
      TI2 = TI**2.5
C
C      CONTINUUM TERMS
      CON = (5.43E-10*PDA*TI**3 + 1.8E-8*PVAP*TI**7.5)*PVAP*F*F 
C
C      ADD RESONANCES
      SUM = 0.
      DO 30 I=1,NLINES
      WIDTH = W3(I)*PDA*TI**X(I) + WS(I)*PVAP*TI**XS(I)
      SHIFT = SR(I)*WIDTH  ! unknown temperature dependence
      WSQ = WIDTH*WIDTH
      S = S1(I)*TI2*EXP(B2(I)*(1.-TI))
      DF(1) = F - FL(I) - SHIFT
      DF(2) = F + FL(I) + SHIFT
C  USE CLOUGH'S DEFINITION OF LOCAL LINE CONTRIBUTION
      BASE = WIDTH/(562500. + WSQ)
C  DO FOR POSITIVE AND NEGATIVE RESONANCES
      RES = 0.
      DO 20 J=1,2
      IF(ABS(DF(J)).LT.750.) RES = RES + WIDTH/(DF(J)**2+WSQ) - BASE
20    CONTINUE
30    SUM = SUM + S*RES*(F/FL(I))**2
      ABH2O = .3183E-4*DEN*SUM + CON
      RETURN
      END subroutine abs_H2O
      
      subroutine abs_O2(TEMP,PRES,VAPDEN,FREQY,O2ABS)
C
C     C. Kummerow, 8/2003.  Changed function to subroutine      
C  Copyright (c) 2003 Massachusetts Institute of Technology
C
C     PURPOSE: RETURNS ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
C              IN NEPERS/KM
C
C      5/1/95  P. Rosenkranz 
C      11/5/97  P. Rosenkranz - 1- line modification.
c      12/16/98 pwr - updated submm freq's and intensities from HITRAN96
c      8/21/02  pwr - revised width at 425
c      3/20/03  pwr - 1- line mixing and width revised
C
c     IMPLICIT NONE
C
C     ARGUMENTS:
      REAL TEMP,PRES,VAPDEN,FREQY
C
C     NAME    UNITS    DESCRIPTION        VALID RANGE
C
C     TEMP    KELVIN   TEMPERATURE        UNCERTAIN, but believed to be
c                                          valid for atmosphere
C     PRES   MILLIBARS PRESSURE           3 TO 1000
C     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
C                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
C     FREQ    GHZ      FREQUENCY          0 TO 900
C
C     REFERENCES FOR EQUATIONS AND COEFFICIENTS:
C     P.W. Rosenkranz, CHAP. 2 and appendix, in ATMOSPHERIC REMOTE SENSING
C      BY MICROWAVE RADIOMETRY (M.A. Janssen, ed., 1993).
C     H.J. Liebe et al, JQSRT V.48, pp.629-643 (1992).
c     M.J. Schwartz, Ph.D. thesis, M.I.T. (1998).
c     A.F. Krupnov et al, J. Mol. Spect. v.215, pp.309-311 (2002).
C     M.Yu. Tretyakov et al, J. Mol. Spect. (2003 preprint).
C     SUBMILLIMETER LINE INTENSITIES FROM HITRAN96.
c
c     This version differs from Liebe's MPM92 in these significant respects:
c     1. The 1- line has the width and mixing coefficient measured by 
c      Tretyakov et al. 
c     2. It modifies the 1- line width temperature dependence to (1/T)**0.9
c     3. It uses the same temperature dependence (X) for submillimeter 
c      line widths as in the 60 GHz band: (1/T)**0.8 
c     4. The 425 GHz line width is from Krupnov et al.
C
c     Local variables:
      REAL TH,TH1,B,PRESWV,PRESDA,DEN,DENS,DFNR,SUM,STR,Y,SF1,SF2,FCEN
      INTEGER K
      REAL X,WB300,W300(40),F(40),Y300(40),S300(40),V(40),BE(40)
      REAL O2ABS, DF
      COMMON /O2COM/ X,WB300,W300,F,Y300,S300,V,BE
C      LINES ARE ARRANGED 1-,1+,3-,3+,ETC. IN SPIN-ROTATION SPECTRUM
      DATA F/118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,
     2  59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,
     3  56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,
     4  55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,
     5  53.5957, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,
     6  52.0214, 67.3696, 51.5034, 67.9009, 368.4984, 424.7632,
     7  487.2494, 715.3931, 773.8397, 834.1458/
        DATA S300/.2936E-14,.8079E-15, .2480E-14,.2228E-14,
     &  .3351E-14,.3292E-14, .3721E-14,.3891E-14,
     &  .3640E-14,.4005E-14, .3227E-14,.3715E-14,
     &  .2627E-14,.3156E-14, .1982E-14,.2477E-14,
     &  .1391E-14,.1808E-14, .9124E-15,.1230E-14,
     &  .5603E-15,.7842E-15, .3228E-15,.4689E-15,
     &  .1748E-15,.2632E-15, .8898E-16,.1389E-15,
     &  .4264E-16,.6899E-16, .1924E-16,.3229E-16,
     &  .8191E-17,.1423E-16, .6494E-15, .7083E-14, .3025E-14,
     &  .1835E-14, .1158E-13, .3993E-14/
      DATA BE/.009,.015, .083,.084, 2*.212, 2*.391, 2*.626,
     & 2*.915, 2*1.260, 1.660,1.665, 2.119,2.115, 2.624,2.625,
     & 2*3.194, 2*3.814, 2*4.484, 2*5.224, 2*6.004, 2*6.844,
     & 2*7.744, .048, .044, .049, .145, .141, .145/
C      WIDTHS IN MHZ/MB
      DATA WB300/.56/, X/.8/
      DATA W300/1.67, 1.646, 1.468, 1.449, 1.382, 1.360,
     & 1.319, 1.297, 1.266, 1.248, 1.221, 1.207, 1.181, 1.171,
     & 1.144, 1.139, 1.110, 1.108, 1.079, 1.078, 2*1.05,
     & 2*1.02,2*1.00,2*.97,2*.94,2*.92,2*.89, 3*1.64, 3*1.81/
      DATA Y300/  -0.036,  0.2408, -0.3486,  0.5227,
     & -0.5430,  0.5877, -0.3970,  0.3237, -0.1348,  0.0311,
     &  0.0725, -0.1663,  0.2832, -0.3629,  0.3970, -0.4599,
     &  0.4695, -0.5199,  0.5187, -0.5597,  0.5903, -0.6246,
     &  0.6656, -0.6942,  0.7086, -0.7325,  0.7348, -0.7546,
     &  0.7702, -0.7864,  0.8083, -0.8210,  0.8439, -0.8529, 6*0./
      DATA V/  0.0079, -0.0978,  0.0844, -0.1273,
     &  0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584,
     &  0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675,
     &  0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590,
     &  0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091,
     &  0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545, 6*0./
C
      TH = 300./TEMP
      TH1 = TH-1.
      B = TH**X
      PRESWV = VAPDEN*TEMP/217.
      PRESDA = PRES -PRESWV
      DEN = .001*(PRESDA*B + 1.1*PRESWV*TH)
      DENS = .001*(PRESDA*TH**.9 + 1.1*PRESWV*TH)
      DFNR = WB300*DEN
      SUM = 1.6E-17*FREQY*FREQY*DFNR/(TH*(FREQY*FREQY + DFNR*DFNR))
      DO 32 K=1,40
      IF(K.EQ.1) THEN !exception for 1- line
        DF = W300(1)*DENS
      ELSE
        DF = W300(K)*DEN
      ENDIF
      FCEN = F(K)
      Y = .001*PRES*B*(Y300(K)+V(K)*TH1)
      STR = S300(K)*EXP(-BE(K)*TH1)
      SF1 = (DF + (FREQY-FCEN)*Y)/((FREQY-FCEN)**2 + DF*DF)
      SF2 = (DF - (FREQY+FCEN)*Y)/((FREQY+FCEN)**2 + DF*DF)
32    SUM = SUM + STR*(SF1+SF2)*(FREQY/F(K))**2
      O2ABS = .5034E12*SUM*PRESDA*TH**3/3.14159
      O2ABS = AMAX1(O2ABS,0.)
      RETURN
      END subroutine abs_O2

      Subroutine ABS_N2(T,P,F,ABSN2)
      
      REAL  T, P, F, ABSN2, TH, FDEPEN, BF

C
C     C. Kummerow, 8/2003.  Changed function to subroutine      
C  Copyright (c) 2002 Massachusetts Institute of Technology
C     ABSN2 = COLLISION-INDUCED ABSORPTION COEFFICIENT (NEPER/KM)

C     IN AIR

C     T = TEMPERATURE (K)

C     P = PRESSURE (MB)

C     F = FREQUENCY (GHZ)(valid 0-1000 GHz)

C

c     5/22/02 P.Rosenkranz

c

C     Equations based on:

C      Borysow, A, and L. Frommhold, 

C      Astrophysical Journal, v.311, pp.1043-1057 (1986)

C     with modification of 1.29 to account for O2-O2 and O2-N2

c     collisions, as suggested by

C      J.R. Pardo, E.Serabyn, J.Cernicharo, J. Quant. Spectros.

c      Radiat. Trans. v.68, pp.419-433 (2001).

c

      TH = 300./T

      FDEPEN = .5 + .5/(1.+(F/450.)**2)

      BF = 6.5E-14*FDEPEN*P*P*F*F*TH**3.6

      ABSN2 = 1.29*BF

      RETURN

      END subroutine ABS_N2
      
! --------------------------------------------------------------------------------------

 	  subroutine eddington(nlyr, Tskin, Z, lyrtemp, kext, salb,
     +        asym, emis, ebar, lambert, ang, freq, ssws, pol, tb,
     +		reflect)

C
C     CHRIS KUMMEROW
C     INCLUDES ASYMPTOTIC EXPRESSIONS FOR TERM3, TERM4, AND TERM5 IN
C     THE LIMIT OF SMALL EFFECTIVE OPTICAL DEPTHS; BILL OLSON FEB, 1995.
C
C     FIXED BUG IN TERM3 OF UPWELLING RADIANCE CK (July 2007)
C
      logical lambert
      integer  nlyr
      real  lyrtemp(0:nlyr)
      real  Tskin
      real  emis, ebar
      real  ang	! added, view angle
      real  freq, ssws  ! added for getting scattering correction term
      integer  pol  ! added input for omega
      real  omega(2)	! added
      real  kext(nlyr), salb(nlyr), asym(nlyr)
      real  L(nlyr), H(nlyr), B0(nlyr), B1(nlyr)
      real  W(2*nlyr,2*nlyr), BB(2*nlyr), DP(nlyr), DM(nlyr)
      real  Z(0:nlyr), IOUT(0:nlyr), I_IN(nlyr+1,100), MU, NU
      real  TB, FISOT, RCOND
      real  XNU, XA, XB, XC, XD, YA, YB, YC, DDZ
      real  TERM1, TERM2, TERM3, TERM4, TERM5, XIUP
      integer J, I, NANG, NN
      real  reflect	! added for testing
c
      data FISOT / 2.7 / 
c
c     Plane parallel Eddington computes Tb for one vertical hydrometeor
c     column at a time.  The code therefore loops over nx & ny
c
                            ! Construct the 1-D vertical column
          DO 20  J = 1,NLYR
            B0(J) = lyrtemp(J-1)
            B1(J) = (lyrtemp(J) - lyrtemp(J-1))/(Z(J) - Z(J-1))
            L(J) = SQRT(3.*KEXT(J)*KEXT(J)*(1. - SALB(J))*
     $            (1. - SALB(J)*ASYM(J)))
            H(J) = 1.5*KEXT(J)*(1. - SALB(J)*ASYM(J))
  20      CONTINUE
C
C         FILL IN THE MATRIX ELEMENTS WHICH FORM THE BOUNDARY CONDITIONS
C         AT THE TOP, BOTTOM AND LAYER INTERFACES OF THE CLOUD.  THERE ARE
C         TWO QUANTITIES, D+ "DP", AND D- "DM" AT EACH BOUNDARY, SO THE 
C         MATRIX HAS DIMENSIONS  2*NLYR BY 2*NLYR.
C         ORDER OF D'S:  D+(1),D-(1),D+(2),D-(2), .... , D+(NLYR),D-(NLYR)

C         SET ALL MATRIX ELEMENTS TO ZERO	
          DO 45  I = 1,2*NLYR
            DO 45  J = 1,2*NLYR
              W(I,J) = 0.0
  45      CONTINUE	

C         FILL IN THE NON-ZERO MATRIX ELEMENTS
          W(1,1)   = ((EBAR - 2.)*L(1)/H(1)) + EBAR
          W(1,2)   = ((2. - EBAR)*L(1)/H(1)) + EBAR
          DO 50  I = 2,2*(NLYR-1),2
            W(I,I-1)   =  (1. - L(I/2)/H(I/2))*EXP(+L(I/2)*(Z(I/2)-
     $                    Z(I/2-1)))
            W(I,I  )   =  (1. + L(I/2)/H(I/2))*EXP(-L(I/2)*(Z(I/2)-
     $                    Z(I/2-1)))
            W(I,I+1)   = -(1. - L(I/2+1)/H(I/2+1))
            W(I,I+2)   = -(1. + L(I/2+1)/H(I/2+1))

            W(I+1,I-1) =  (1. + L(I/2)/H(I/2))*EXP(+L(I/2)*(Z(I/2)-
     $                    Z(I/2-1)))
            W(I+1,I)   =  (1. - L(I/2)/H(I/2))*EXP(-L(I/2)*(Z(I/2)-
     $                    Z(I/2-1)))
            W(I+1,I+1) = -(1. + L(I/2+1)/H(I/2+1))
            W(I+1,I+2) = -(1. - L(I/2+1)/H(I/2+1))
  50      CONTINUE
          W(2*NLYR,2*NLYR-1) =  (1. + L(NLYR)/H(NLYR))*EXP(+L(NLYR)*
     $                          (Z(NLYR)-Z(NLYR-1)))
          W(2*NLYR,2*NLYR)   =  (1. - L(NLYR)/H(NLYR))*EXP(-L(NLYR)*
     $                          (Z(NLYR)-Z(NLYR-1)))

C         FILL IN THE ROW OF CONSTANTS IN THE LINEAR EQUATIONS
          BB(1)    = EBAR*Tskin - EBAR*B0(1) - 
     $               (EBAR - 2.)*B1(1)/H(1)
          DO 55  I = 2,2*(NLYR-1),2
            BB(I)   =  + B1(I/2)/H(I/2) - B1(I/2+1)/H(I/2+1)
            BB(I+1) =  - B1(I/2)/H(I/2) + B1(I/2+1)/H(I/2+1)
  55      CONTINUE
          BB(2*NLYR)  =  FISOT - B0(NLYR) - B1(NLYR)*(Z(NLYR) - 
     $                   Z(NLYR-1) + 1/H(NLYR))
C
C         MATRIX INVERSION IN DONE IN SUBROUTINE LINPAK
          CALL LINPAK(NLYR, W, BB, RCOND)
C
          DO 60  I = 1,NLYR
            DP(I) = BB(2*I-1)
            DM(I) = BB(2*I)
  60      CONTINUE

C         AFTER D'S ARE KNOWN, CALCULATE SURFACE RADIANCE

          MU =  cos(ang*0.01745329)
          NU = -MU
          
C         FOR THE FOLLOWING CALCULATIONS, REFER TO APPENDIX B OF THESIS
C         *************************************************************
      
          IF ( LAMBERT ) THEN 
C           CALCULATE THE DOWNWELLING FLUX AT 81 ANGLES
            NANG = 81
            DO 997  NN = 1,NANG
              XNU = -(2.*NN - 1.)/(NANG*2.)
              I_IN(NLYR+1,NN) = FISOT
C             LOOP THROUGH THE REMAINING LAYERS
              DO 100  J = NLYR,1,-1
C               CALCULATE RADIANCE FROM TOP OF LAYER "J"
                XA = B0(J) - 1.5*SALB(J)*ASYM(J)*XNU*B1(J)/H(J)
                XB = B1(J)
                XC = SALB(J)*DP(J)*(1. - 1.5*ASYM(J)*XNU*L(J)/H(J))               
                XD = SALB(J)*DM(J)*(1. + 1.5*ASYM(J)*XNU*L(J)/H(J))
                YA = KEXT(J)/XNU
                YB = YA + L(J)
                YC = YA - L(J)
                DDZ = Z(J) - Z(J-1)

                TERM1 = I_IN(J+1,NN)*EXP(YA*DDZ)
                TERM2 = XA*(1. - EXP(YA*DDZ))
                IF (ABS(YA*DDZ) .GT. 1.E-05) THEN
                  TERM3 = XB/YA*(EXP(YA*DDZ) - YA*DDZ - 1.)
                  ! previous to 7/2007 fix was:
                  ! TERM3 = XB/YA*(EXP(YA*DZ)*(1. - YA*DZ) - 1.)
                ELSE
                  TERM3 = 0.
                  ! previous to 7/2007 was: TERM3=-XB*YA*DZ*DZ
                ENDIF
                IF (ABS(YB*DDZ) .GT. 1.E-05 ) THEN
                  TERM4 = XC*YA/YB*(1. - EXP(YB*DDZ))
                ELSE
                  TERM4 = -XC*YA*DDZ
                ENDIF
                IF (ABS(YC*DDZ) .GT. 1.E-05 ) THEN
                  TERM5 = XD*YA/YC*(1. - EXP(YC*DDZ))
                ELSE
                  TERM5 = -XD*YA*DDZ
                ENDIF
                I_IN(J,NN) = TERM1 + TERM2 + TERM3 + TERM4 + TERM5
  100         CONTINUE
  997       CONTINUE
C
C           CALCULATE THE TOTAL DOWNWELLING FLUX REACHING THE SURFACE
            XIUP = 0.
            DO 47  NN = 1,NANG
              XIUP = XIUP + I_IN(1,NN)*(1./NANG)*(2.*NN-1.)/(2.*NANG)
  47        CONTINUE
            XIUP = 2.*XIUP
           
          ELSE

C           CALCULATE THE DOWNWELLING FLUX AT ANGLE MU ONLY
            NN = 22            ! THIS IS A DUMMY INDEX FOR I_IN
            I_IN(NLYR+1,NN) = FISOT
C           LOOP THROUGH THE REMAINING LAYERS
            DO 110  J = NLYR,1,-1
C             CALCULATE RADIANCE FROM TOP OF LAYER "J"
              XA = B0(J) - 1.5*SALB(J)*ASYM(J)*NU*B1(J)/H(J)
              XB = B1(J)
              XC = SALB(J)*DP(J)*(1. - 1.5*ASYM(J)*NU*L(J)/H(J))               
              XD = SALB(J)*DM(J)*(1. + 1.5*ASYM(J)*NU*L(J)/H(J))
              YA = KEXT(J)/NU
              YB = YA + L(J)
              YC = YA - L(J)
              DDZ = Z(J) - Z(J-1)
	      
              TERM1 = I_IN(J+1,NN)*EXP(YA*DDZ)
              TERM2 = XA*(1. - EXP(YA*DDZ))
              IF (ABS(YA*DDZ) .GT. 1.E-05) THEN
                TERM3 = XB/YA*(EXP(YA*DDZ) - YA*DDZ - 1.)
                ! previous to 7/2007 fix was:
                ! TERM3 = XB/YA*(EXP(YA*DZ)*(1. - YA*DZ) - 1.)
              ELSE
                TERM3 = 0.
                ! previous to 7/2007 was: TERM3=-XB*YA*DZ*DZ
              ENDIF
              IF (ABS(YB*DDZ) .GT. 1.E-05 ) THEN
                TERM4 = XC*YA/YB*(1. - EXP(YB*DDZ))
              ELSE
                TERM4 = -XC*YA*DDZ
              ENDIF
              IF (ABS(YC*DDZ) .GT. 1.E-05 ) THEN
                TERM5 = XD*YA/YC*(1. - EXP(YC*DDZ))
              ELSE
                TERM5 = -XD*YA*DDZ
              ENDIF
              I_IN(J,NN) = TERM1 + TERM2 + TERM3 + TERM4 + TERM5
  110       CONTINUE
            XIUP = I_IN(1,22)
C
          ENDIF
C
	  IOUT(0) = EMIS*Tskin + reflect*XIUP
          DO 101  J = 1,NLYR
C           CALCULATE THE UPWELLING RADIANCES AT THE TOP OF EACH LAYER J
            XA = B0(J) - 1.5*SALB(J)*ASYM(J)*MU*B1(J)/H(J)
            XB = B1(J)
            XC = SALB(J)*DP(J)*(1. - 1.5*ASYM(J)*MU*L(J)/H(J))               
            XD = SALB(J)*DM(J)*(1. + 1.5*ASYM(J)*MU*L(J)/H(J))
            YA = KEXT(J)/MU
            YB = YA + L(J)
            YC = YA - L(J)
            DDZ = Z(J) - Z(J-1)
	    	    
            TERM1 = IOUT(J-1)*EXP(-YA*DDZ)
            TERM2 = XA*(1. - EXP(-YA*DDZ))
            IF ( ABS(YA*DDZ) .GT. 1.E-05 ) THEN
              TERM3 = XB/YA*(EXP(-YA*DDZ) + YA*DDZ - 1.)
            ELSE
              TERM3 = 0.
            ENDIF
            IF ( ABS(YB*DDZ) .GT. 1.E-05 ) THEN
              TERM4 = XC*YA/YB*(EXP( (YB-YA)*DDZ ) - EXP(-YA*DDZ) )
            ELSE            
              TERM4 = XC*YA*DDZ*EXP(-YA*DDZ)
            ENDIF
            IF (ABS(YC*DDZ) .GT. 1.E-05 ) THEN
              TERM5 = XD*YA/YC*EXP(-YA*DDZ)*(EXP(YC*DDZ) - 1.)
            ELSE
              TERM5 = XD*YA*DDZ*EXP(-YA*DDZ)
            ENDIF
            IOUT(J) = TERM1 + TERM2 + TERM3 + TERM4 + TERM5
  101     CONTINUE
C
          TB = IOUT(NLYR)

      return
      end subroutine eddington
      
! --------------------------------------------------------------------------------------

      SUBROUTINE  EMIT( F, NPOL, TS, W, ANGLE, EMIS, EBAR )

      INTEGER, PARAMETER ::   NANG = 21
      REAL  ANG(NANG), MU(NANG), ESUM(NANG)
      REAL  F, TS, W, ANGLE, EMIS, EBAR
      REAL  PI, S, EV, EH, EMISH, EMISV, ANGLES, SUM
      REAL  EAVG, DMU, AVMU
      INTEGER NPOL, I
      
      PI = 2.*ASIN(1.0)
      S = 35.                                    ! SALINITY IN PPM
      ! ANGLE = ACOS(UMU)*180./PI
C      
C     CALCULATE EMIS AT GIVEN ANGLE
      CALL EMISS (F,ANGLE,S,TS,W,EV,EH,EMISH,EMISV)
      IF ( NPOL .EQ. 0 ) EMIS = EH
      IF ( NPOL .EQ. 1 ) EMIS = EV
C
C     CALCULATE EMIS AT VARIOUS ANGLES
      DO 58  I = 1,NANG
       ANG(I) = 4.*( I - 1 )
       MU(I) = COS(ANG(I)*PI/180.)
       ANGLES = ANG(I)
       CALL EMISS (F,ANGLES,S,TS,W,EV,EH,EMISH,EMISV)
       ESUM(I) =  EV + EH 
  58  CONTINUE
C
C     CALCULATE EBAR
      SUM = 0.0
      DO 59  I = 1,NANG-1
       EAVG = 0.5*( ESUM(I) + ESUM(I+1) )
       DMU = MU(I) - MU(I+1)
       AVMU = 0.5*( MU(I) + MU(I+1) )
       SUM = SUM + EAVG*AVMU*DMU
  59  CONTINUE
      EBAR = SUM
      RETURN
      END SUBROUTINE EMIT
C
C
      SUBROUTINE  DIECON(S,T,FREQ,E1,E2)

      save sold, told, st, s2, t2, sst, stt, sstt, es, tau, sigma
      REAL  TWOPI, EINF, SOLD, TOLD
      REAL S, T, FREQ, E1, E2
      REAL ST, S2, T2, SST, STT, SSTT, ES, TAU, SIGMA, ZNU
      REAL OMEGA, DEN
      DATA TWOPI /6.283185307/ , EINF /4.9/ , SOLD /0.0/ , TOLD /-99./
      IF (S .EQ. SOLD .AND. T .EQ. TOLD) GO TO 10
      ST = S*T
      S2 = S*S
      T2 = T*T
      SST = S2*T
      STT = T2*S
      SSTT = S2*T2
      ES = 88.-4.339E-01*S+1.71E-03*S2-4.035E-01*T+8.065E-04*T2+6.170
     $  E-03 * ST-8.910E-05*SST-6.934E-05*STT+1.439E-06*SSTT
C
      TAU = (18.70-7.924E-02*S+6.35E-04*S2-5.489E-01*T+5.758E-03*T2+
     $1.889E-03*ST-7.209E-06*SST-5.299E-07*STT-2.101E-07*SSTT)*1.0E-12
C
      SIGMA = (7.788E-03*S-1.672E-06*S2-8.570E-15*T+2.996E-16*T2+4.059E
     $     -04 * ST-3.215E-06*SST-1.423E-06*STT+3.229E-08*SSTT)*1.0E11
C
   10 ZNU = FREQ*1.E09
      OMEGA = TWOPI*ZNU
      DEN = 1. + (OMEGA * TAU) ** 2
      E1 = (ES-EINF)/DEN+EINF
      E2 = (ES-EINF)*OMEGA*TAU/DEN+2.*SIGMA/ZNU
      SOLD = S
      TOLD = T
C
      RETURN
      END SUBROUTINE DIECON
C
C
      SUBROUTINE  EMISS(F,ANGLE,S,TS,W,EV,EH,EMISH,EMISV)
      REAL F, ANGLE, S, TS, W, EV, EH, EMISH, EMISV
      REAL DTR, T, THETA, CMHU, CMHU2, FACT1, B, B2, ARG
      REAL PSIZ, COSPZ, F11, PSIY, G, COSPY, F22
      REAL RSH, RSV, FOAM, GH, GV, A1, RFV, RFH, Y
      REAL SQRTF, CORRV, CORRH, RRV, RRH, RV, RH
      REAL E1, E2
      DATA DTR / 0.01745329252 /
      T = TS-273.16
      THETA = ANGLE*DTR
C
      CALL DIECON(S,T,F,E1,E2)
C
      CMHU = COS (THETA)
      CMHU2 = CMHU*CMHU
      FACT1 = E1+CMHU2-1.
      B = (FACT1*FACT1+E2*E2)**0.25
      B2 = B*B
      ARG = E2/FACT1
      PSIZ = 0.5*ATAN(ARG)
      COSPZ = COS(PSIZ)*CMHU*B
      F11 = CMHU2+B2+2.*COSPZ
      PSIY = ATAN(E2/E1)
      G = CMHU*SQRT(E1*E1+E2*E2)
      COSPY = COS(PSIY-PSIZ)*G*B
      F22 = G*G+B2+2.*COSPY
C
C     FOR SPECULAR SURFACES THE HORIZONTAL AND VERTICAL EMISSIVITY
C     CALCULATION
C
      EMISH = 4.*COSPZ/F11
      EMISV = 4.*COSPY/F22
C
C     FROM HERE THE EFFECT OF SURFACE ROUGHNESS AND FOAM ARE INCLUDED
C     BASED ON HOLLINGER MODEL AND FOAM MODEL OF STOGRYN
C
      RSH = 1.-EMISH
      RSV = 1.-EMISV
C
C     P0 = 1.707476E-2+8.560329E-4*F+1.120024E-5*F*F
C     P1 = -1.500792E-2+1.820672E-3*F-4.633806E-5*F*F
C     P2 = 2.442217E-4-2.282022E-6*F+4.194212E-7*F*F
C
C     FOAM = (P0+P1*W+P2*W*W)
C
      FOAM = 7.751E-06 * W ** 3.231
C
      GH = 1.-1.748E-3*ANGLE-7.336E-5*ANGLE**2+1.044E-7*ANGLE**3
      GV = 1.-9.946E-4*ANGLE+3.218E-5*ANGLE**2-1.187E-6*ANGLE**3
     $   +7.0E-20*ANGLE**10
C
      A1 = (208.0+1.29*F)/TS
C
C     RFV = 1.-A1*GV-0.005*F
C     RFH = 1.-A1*GH-0.005*F
C
      RFV = 1. - A1 * GV
      RFH = 1. - A1 * GH
C
      Y = 7.32E-02*ANGLE
C
C     TS SURFACE TEMP IS IN DEGREE KELVIN
C
      SQRTF = SQRT(F)
C
C     CORRV = (W*(1.17E-01-2.09E-03*EXP(Y))*SQRTF/TS)-0.00065*F
C     CORRH = (W*(1.15E-01+3.80E-05*ANGLE**2)*SQRTF/TS)-0.00065*F
C
      CORRV = (W*(1.17E-01-2.09E-03*EXP(Y))*SQRTF/TS)
      CORRH = (W*(1.15E-01+3.80E-05*ANGLE**2)*SQRTF/TS)
C
      RRV = RSV-CORRV
      RRH = RSH-CORRH
C
      RV = RRV*(1.-FOAM)+RFV*FOAM
      RH = RRH*(1.-FOAM)+RFH*FOAM
C
      EH = 1.-RH
      EV = 1.-RV
C
      RETURN
      END SUBROUTINE EMISS
      
! --------------------------------------------------------------------------------------

      SUBROUTINE ICEOPTIC (freqy, temp, epsreal, epsimag)

**    Hufford (1991), see Brussard and Watson (1995), p.297


C     Input & output variables; eps the dielectric constant, epsilon
      real     freqy, temp
      real     epsreal, epsimag
      
C     internal variables 
      real     t_ice, theta, A , B    
**
      epsreal = 3.15

      if (temp .gt. 273.16) then
        t_ice = 273.16
      else
        t_ice = temp
      endif

      theta  = 300.0 / t_ice
      A      = 1.0 E-04 * (50.4 + 62.0 * (theta - 1.0))
     +         * EXP (-22.1 * (theta - 1.0))
      B      = 1.0 E-04 * (0.633 / theta - 0.131)
     +         + (7.36 E-04 * theta / (theta - 0.9927))
     +         * (7.36 E-04 * theta / (theta - 0.9927))
      epsimag = A / freqy + B * freqy
**
      return
      end subroutine iceoptic

! -------------------------------------------------------------------------------------

      SUBROUTINE LINPAK(NLYR, A, B, RCOND)
      
      INTEGER  NLYR
      REAL RCOND
      REAL A(2*nlyr, 2*nlyr), B(2*nlyr), Z(2*nlyr)
      INTEGER LDA, N, IPVT(2*nlyr)
      
      LDA = 2*nlyr
      N   = 2*NLYR

      CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
      CALL SGESL(A,LDA,N,IPVT,B,0)
      RETURN
      END SUBROUTINE LINPAK
      
      SUBROUTINE SGECO(A,LDA,N,IPVT,RCOND,Z)
C***BEGIN PROLOGUE  SGECO
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D2A1
C***KEYWORDS  CONDITION,FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION AND ESTIMATES
C            THE CONDITION NUMBER OF THE MATRIX.
C***DESCRIPTION
C
C     SGECO FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW SGECO BY SGESL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW SGECO BY SGEDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW SGECO BY SGEDI.
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   REAL
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       REAL(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK SGEFA
C     BLAS SAXPY,SDOT,SSCAL,SASUM
C     FORTRAN ABS,AMAX1,SIGN
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  SASUM,SAXPY,SDOT,SGEFA,SSCAL
C***END PROLOGUE  SGECO
      INTEGER LDA,N,IPVT(1)
      REAL A(LDA,1),Z(1)
      REAL RCOND
C
      REAL EK,T,WK,WKM
      REAL ANORM,S,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
C
C     COMPUTE 1-NORM OF A
C
C***FIRST EXECUTABLE STATEMENT  SGECO
      ANORM = 0.0E0
      DO 10 J = 1, N
         ANORM = AMAX1(ANORM,SASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL SGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0E0
      DO 20 J = 1, N
         Z(J) = 0.0E0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0E0) EK = SIGN(EK,-Z(K))
         IF (ABS(EK-Z(K)) .LE. ABS(A(K,K))) GO TO 30
            S = ABS(A(K,K))/ABS(EK-Z(K))
            CALL SSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = ABS(WK)
         SM = ABS(WKM)
         IF (A(K,K) .EQ. 0.0E0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0E0
            WKM = 1.0E0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + ABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + ABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + SDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0E0) GO TO 110
            S = 1.0E0/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
C
      YNORM = 1.0E0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL SAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (ABS(Z(K)) .LE. 1.0E0) GO TO 130
            S = 1.0E0/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (ABS(Z(K)) .LE. ABS(A(K,K))) GO TO 150
            S = ABS(A(K,K))/ABS(Z(K))
            CALL SSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0E0) Z(K) = 1.0E0
         T = -Z(K)
         CALL SAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0E0/SASUM(N,Z,1)
      CALL SSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
      RETURN
      END SUBROUTINE SGECO
      
      SUBROUTINE SGESL(A,LDA,N,IPVT,B,JOB)
C***BEGIN PROLOGUE  SGESL
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D2A1
C***KEYWORDS  LINEAR ALGEBRA,LINPACK,MATRIX,SOLVE
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  SOLVES THE REAL SYSTEM A*X=B OR TRANS(A)*X=B
C            USING THE FACTORS OF SGECO OR SGEFA
C***DESCRIPTION
C
C     SGESL SOLVES THE REAL SYSTEM
C     A * X = B  OR  TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY SGECO OR SGEFA.
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE OUTPUT FROM SGECO OR SGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM SGECO OR SGEFA.
C
C        B       REAL(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
C        OR SGEFA HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY,SDOT
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  SAXPY,SDOT
C***END PROLOGUE  SGESL
      INTEGER LDA,N,IPVT(1),JOB
      REAL A(LDA,1),B(1)
C
      REAL T
      INTEGER K,KB,L,NM1
C***FIRST EXECUTABLE STATEMENT  SGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL SAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL SAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = SDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + SDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END SUBROUTINE SGESL
      
      FUNCTION SASUM(N,SX,INCX)
C
C***PURPOSE  SUM OF MAGNITUDES OF S.P VECTOR COMPONENTS
C
C     --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
C       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
C     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
C
C     --OUTPUT--
C    SASUM  SINGLE PRECISION RESULT (ZERO IF N .LE. 0)
C
C     RETURNS SUM OF MAGNITUDES OF SINGLE PRECISION SX.
C     SASUM = SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))
C
      INTEGER N, INCX, NS, M, MP1, I
      REAL SX(1), SASUM
C
      SASUM = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I=1,NS,INCX
          SASUM = SASUM + ABS(SX(I))
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 6.
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SASUM = SASUM + ABS(SX(I))
   30 CONTINUE
      IF( N .LT. 6 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        SASUM = SASUM + ABS(SX(I)) + ABS(SX(I + 1)) + ABS(SX(I + 2))
     1  + ABS(SX(I + 3)) + ABS(SX(I + 4)) + ABS(SX(I + 5))
   50 CONTINUE
      RETURN
      END FUNCTION SASUM
      
      SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
C
C***PURPOSE  S.P. COMPUTATION Y = A*X + Y
C
C     --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
C       SA  SINGLE PRECISION SCALAR MULTIPLIER
C       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
C     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
C       SY  SINGLE PRECISION VECTOR WITH N ELEMENTS
C     INCY  STORAGE SPACING BETWEEN ELEMENTS OF SY
C
C     --OUTPUT--
C       SY  SINGLE PRECISION RESULT (UNCHANGED IF N .LE. 0)
C
C     OVERWRITE SINGLE PRECISION SY WITH SINGLE PRECISION SA*SX +SY.
C     FOR I = 0 TO N-1, REPLACE  SY(LY+I*INCY) WITH SA*SX(LX+I*INCX) +
C       SY(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N
C       AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C
      INTEGER N, INCX, INCY
      REAL SX(1),SY(1),SA
      INTEGER IX, IY, I, M, MP1, NS
C
      IF(N.LE.0.OR.SA.EQ.0.E0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SY(IY) = SY(IY) + SA*SX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SY(I) = SY(I) + SA*SX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        SY(I) = SY(I) + SA*SX(I)
        SY(I + 1) = SY(I + 1) + SA*SX(I + 1)
        SY(I + 2) = SY(I + 2) + SA*SX(I + 2)
        SY(I + 3) = SY(I + 3) + SA*SX(I + 3)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          SY(I) = SA*SX(I) + SY(I)
   70     CONTINUE
      RETURN
      END SUBROUTINE SAXPY
      
      FUNCTION SDOT(N,SX,INCX,SY,INCY)
C
C***PURPOSE  S.P. INNER PRODUCT OF S.P. VECTORS
C
C     --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
C       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
C     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
C       SY  SINGLE PRECISION VECTOR WITH N ELEMENTS
C     INCY  STORAGE SPACING BETWEEN ELEMENTS OF SY
C
C     --OUTPUT--
C     SDOT  SINGLE PRECISION DOT PRODUCT (ZERO IF N .LE. 0)
C
C     RETURNS THE DOT PRODUCT OF SINGLE PRECISION SX AND SY.
C     SDOT = SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      INTEGER X, INCX, INCY, IX, IY, I, M, MP1, NS, N
      REAL SX(1),SY(1), SDOT
      
C
      SDOT = 0.0E0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1)5,20,60
    5 CONTINUE
C
C        CODE FOR UNEQUAL INCREMENTS OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        SDOT = SDOT + SX(IX)*SY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SDOT = SDOT + SX(I)*SY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SDOT = SDOT + SX(I)*SY(I) + SX(I + 1)*SY(I + 1) +
     1   SX(I + 2)*SY(I + 2) + SX(I + 3)*SY(I + 3) + SX(I + 4)*SY(I + 4)
   50 CONTINUE
      RETURN
C
C        CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS=N*INCX
      DO 70 I=1,NS,INCX
        SDOT = SDOT + SX(I)*SY(I)
   70   CONTINUE
      RETURN
      END FUNCTION SDOT
      
      SUBROUTINE SSCAL(N,SA,SX,INCX)
C
C***PURPOSE  S.P. VECTOR SCALE X = A*X
C
C     --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
C       SA  SINGLE PRECISION SCALE FACTOR
C       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
C     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
C
C     --OUTPUT--
C       SX  SINGLE PRECISION RESULT (UNCHANGED IF N .LE. 0)
C
C     REPLACE SINGLE PRECISION SX BY SINGLE PRECISION SA*SX.
C     FOR I = 0 TO N-1, REPLACE SX(1+I*INCX) WITH  SA * SX(1+I*INCX)
C
      INTEGER N, INCX, NS, I, M, MP1
      REAL SA,SX(1)
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      NS = N*INCX
          DO 10 I = 1,NS,INCX
          SX(I) = SA*SX(I)
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        SX(I) = SA*SX(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        SX(I) = SA*SX(I)
        SX(I + 1) = SA*SX(I + 1)
        SX(I + 2) = SA*SX(I + 2)
        SX(I + 3) = SA*SX(I + 3)
        SX(I + 4) = SA*SX(I + 4)
   50 CONTINUE
      RETURN
      END SUBROUTINE SSCAL
      
      FUNCTION ISAMAX(N,SX,INCX)
C
C***PURPOSE  FIND LARGEST COMPONENT OF S.P. VECTOR
C
C     --INPUT--
C        N  NUMBER OF ELEMENTS IN INPUT VECTOR(S)
C       SX  SINGLE PRECISION VECTOR WITH N ELEMENTS
C     INCX  STORAGE SPACING BETWEEN ELEMENTS OF SX
C
C     --OUTPUT--
C   ISAMAX  SMALLEST INDEX (ZERO IF N .LE. 0)
C
C     FIND SMALLEST INDEX OF MAXIMUM MAGNITUDE OF SINGLE PRECISION SX.
C     ISAMAX =  FIRST I, I = 1 TO N, TO MINIMIZE  ABS(SX(1-INCX+I*INCX)
C
      INTEGER N, INCX, II, I, NS, ISAMAX
      REAL SX(1),SMAX,XMAG
C
      ISAMAX = 0
      IF(N.LE.0) RETURN
      ISAMAX = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20
C
C        CODE FOR INCREMENTS NOT EQUAL TO 1.
C
      SMAX = ABS(SX(1))
      NS = N*INCX
      II = 1
          DO 10 I=1,NS,INCX
          XMAG = ABS(SX(I))
          IF(XMAG.LE.SMAX) GO TO 5
          ISAMAX = II
          SMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN
C
C        CODE FOR INCREMENTS EQUAL TO 1.
C
   20 SMAX = ABS(SX(1))
      DO 30 I = 2,N
         XMAG = ABS(SX(I))
         IF(XMAG.LE.SMAX) GO TO 30
         ISAMAX = I
         SMAX = XMAG
   30 CONTINUE
      RETURN
      END FUNCTION ISAMAX
      
      SUBROUTINE SGEFA(A,LDA,N,IPVT,INFO)
C***BEGIN PROLOGUE  SGEFA
C***DATE WRITTEN   780814   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  D2A1
C***KEYWORDS  FACTOR,LINEAR ALGEBRA,LINPACK,MATRIX
C***AUTHOR  MOLER, C. B., (U. OF NEW MEXICO)
C***PURPOSE  FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.
C***DESCRIPTION
C
C     SGEFA FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.
C
C     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SAXPY,SSCAL,ISAMAX
C***REFERENCES  DONGARRA J.J., BUNCH J.R., MOLER C.B., STEWART G.W.,
C                 *LINPACK USERS  GUIDE*, SIAM, 1979.
C***ROUTINES CALLED  ISAMAX,SAXPY,SSCAL
C***END PROLOGUE  SGEFA
      INTEGER LDA,N,IPVT(1),INFO
      REAL A(LDA,N)
C
      REAL T
      INTEGER J,K,KP1,L,NM1
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
C***FIRST EXECUTABLE STATEMENT  SGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = ISAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0E0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0E0/A(K,K)
            CALL SSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL SAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0E0) INFO = N
      RETURN
      END SUBROUTINE SGEFA
      
! --------------------------------------------------------------------------------------

      SUBROUTINE MG_ELLIPS (FINCL, EMATRIX, EINCL, EMG)

      IMPLICIT NONE
      SAVE

** Maxwell-Garnett formula for effective permittivity of 2-component media 
** (elliptical inclusions) P. Bauer 1996
**
** FINCL     volume fraction of inclusions
** EMATRIX   permittivity of matrix
** EINCL     permittivity of inclusions
**
** EMG       effective permittivity

      REAL    FINCL
      COMPLEX EMATRIX, EINCL, EMG, GAMMA, Q

      Q     = (EINCL / (EINCL - EMATRIX)) 
     +      * CLOG (EINCL / EMATRIX) - 1.0
      GAMMA = 2.0 * EMATRIX * Q / (EINCL - EMATRIX)

      EMG = ((1.0 - FINCL) * EMATRIX + FINCL * GAMMA * EINCL)
     +    / (1.0 - FINCL + FINCL * GAMMA) 

      RETURN
      END SUBROUTINE MG_ELLIPS
      
! --------------------------------------------------------------------------------------

      subroutine mie_ci(freqy, temp, lwc, ksca, asca, gsca)
      

c     Compute the extinction, absorption, asymmetry parameter and
c     backscatter for a given water content of cloud ice in [g/m^3].

c     Input:
c     freqy		frequency of radiation [GHz]
c     temp		temperature of particles [K]
c     lwc		water content of cloud ice distribution [g/m**3]

c     Output:
c     ksca		extinction coefficient [1/km]
c     asca		single-scatter albedo []
c     gsca		asymmetry factor []
c     pbck		backscatter phase function/(4*pi) []
c

      implicit none

      real     freqy, temp, lwc, reff_cloudice
      real     ksca, asca, gsca, pbck

      real     pi, wavel, x
      real     densice, density_cloudice
      real     rad, dropmass, num
      real     qext, qsca, asym ,qbsca
      real     bext, bsca, bsym, bq11
      
      real     epsreal, epsimag
      real     fincl      
      complex  eice, eair, emg, cref

      data     pi  /3.14159265/
      data     densice    /0.917e+3/            ! [kg/m^3]
      data     reff_cloudice  / 0.10 /          ! [mm]
      data     density_cloudice / 0.917E+3 /    ! [kg/m^3]

c      
c     Assign some useful constants
      wavel = 300./freqy
      
c
c     Begin by checking if hydrometeors of this species are present.
c     If not, set scattering parameters to zero and return.
c
      if(lwc .lt. 0.001) then
        ksca=0.
        asca=0.
        gsca=0.
        pbck=0.
        return
      endif
c
c     Initialize the scattering parameters
c
      bext=0.
      bsca=0.
      bsym=0.
      bq11=0.
c
c     Compute scattering properties
c
      rad = reff_cloudice
      if (density_cloudice .gt. densice) density_cloudice = densice
      dropmass = 4/3.*pi*density_cloudice*rad*rad*rad*1.E-09
      num = lwc/dropmass

c     Get complex refractive index of cloud ice
c
      call iceoptic(freqy,temp,epsreal,epsimag)
      eice = cmplx(epsreal,epsimag)
      eair = cmplx(1.0006,0.0)

c     calculate dielectric constant of cloud ice as
c     ice matrix with air inclusions, using Maxwell-Garnett 
c     mixing for 2-component media w. elliptical inclusions 
      fincl = 1. - density_cloudice/densice
      call mg_ellips(fincl, eice, eair, emg)
      cref = csqrt(emg)
c
c     call Mie program
c 
      x = 2.*pi*rad/wavel
      call mie_sphere(x,cref,qsca,qext,asym,qbsca)

      bext = num*qext*pi*rad*rad*1.e-6
      bsca = num*qsca*pi*rad*rad*1.e-6
      bsym = num*qsca*asym*pi*rad*rad*1.e-6
      bq11 = num*qbsca*pi*rad*rad*1.e-6
c
c     check for distribution with very small extinction;
c     set parameters to zero to avoid numerical problems
      if( bext .gt. 1.e-6) then

        ksca=bext
        asca=bsca/bext
        gsca=bsym/bsca
        pbck=bq11/bsca

      else

        ksca=0.
        asca=0.
        gsca=0.
        pbck=0.

      end if

      return
      end subroutine mie_ci

! --------------------------------------------------------------------------------------

      subroutine mie_clw(freqy, temp, lwc, ksca, asca, gsca)
      

c     Compute the extinction, absorption, asymmetry parameter and
c     backscatter for a given water content of cloud water in [g/m^3]. 
c     Code assumes that cloud water drops are mono-disperse with an 
c     effective radius reff supplied in the parameter file.

c     Input:
c     freqy		frequency of radiation [GHz]
c     temp		temperature of particles [K]
c     lwc		water content of cloud water distribution [g/m**3]
c     reff              effective radius of particles [mm]

c     Output:
c     ksca		extinction coefficient [1/km]
c     asca		single-scatter albedo []
c     gsca		asymmetry factor []
c     pbck		backscatter phase function/(4*pi) []
c

      implicit none

      real    freqy, temp, lwc, reff_cloudwater
      real    ksca, asca, gsca, pbck

      real    pi, wavel, x
      real    densliq, density
      real    rad, dropmass
      real    num
      real    qext, qsca, asym, qbsca
      real    bext, bsca, bsym, bq11
      
      real     epsreal, epsimag
      complex  ewat, cref
      real :: ss=0.0
      !real*8 :: ss=0.0

c      
c     Assign some useful constants
      data       pi /3.14159265/
      data       densliq /1.0e+3/                 ! kg/m^3
      data       reff_cloudwater  / 0.10 /        ! [mm]
  
      wavel = 300./freqy
c
c     Begin by checking if hydrometeors of this species are present.
c     If not, set scattering parameters to zero and return.
c
      if(lwc .lt. 0.001) then
        ksca=0.
        asca=0.
        gsca=0.
        pbck=0.
        return
      endif
c
c     Initialize the scattering parameters
c
      bext=0.
      bsca=0.
      bsym=0.
      bq11=0.
c
c     Compute scattering properties
c
      rad = reff_cloudwater
      density = densliq
      dropmass = 4./3.*pi*density*rad*rad*rad*1.E-09
      num = lwc/dropmass      
c
c     Get complex refractive index of liquid water
c
      call watoptic(freqy,temp,ss,epsreal,epsimag)
!      call watoptic(dble(freqy),dble(temp),ss,dble(epsreal),
!     >              dble(epsimag))
      ewat = cmplx(epsreal,epsimag)
      cref = csqrt(ewat)
c
c     call Mie program
c
      x = 2.*pi*rad/wavel 
      call mie_sphere(x,cref,qsca,qext,asym,qbsca)

      bext=num*qext*pi*rad*rad*1.e-6
      bsca=bsca+num*qsca*pi*rad*rad*1.e-6
      bsym=bsym+num*qsca*asym*pi*rad*1.e-6
      bq11=bq11+num*qbsca*pi*rad*rad*1.e-6
c
c     check for distribution with very small extinction;
c     set parameters to zero to avoid numerical problems
c
      if( bext .gt. 1.e-6) then
        ksca=bext
        asca=bsca/bext
        gsca=bsym/bsca
        pbck=bq11/bsca
      else
        ksca=0.
        asca=0.
        gsca=0.
        pbck=0.
      end if

      return
      end subroutine mie_clw
      
! --------------------------------------------------------------------------------------

      
! ---------------------------------------------------------------------------------------

       subroutine mie_snow(freqy, temp, swc, ksca, asca, gsca)
c      
c     Compute the extinction, absorption, asymmetry parameter and
c     backscatter for a given water content of snow in [g/m^3], and
c     a particle size distribution n(D) with intercept n0s.

c     Input:
c     freqy		frequency of radiation [GHz]
c     temp		temperature of particles [K]
c     iwc		ice water content of snow distribution [g/m**3]

c     Output:
c     ksca		extinction coefficient [1/km]
c     asca		single-scatter albedo []
c     gsca		asymmetry factor []
c     pbck		backscatter phase function/(4*pi) []
c

      implicit none

      integer  i      
      real     freqy, temp, swc
      real     ksca, asca, gsca, pbck
      real     wavel, pi, densice, diam, diam_increment, density, fincl
      real     n0_snow, density_snow, lam, num, x
      real     qext, qsca, asym, qbsca
      real     bext, bsca, bsym, bq11      
      real     epsreal, epsimag      
      complex  eice, eair, emg, cref
      
      data     pi /3.14159265/
      data     densice /0.917e+3/               ! [kg/m^3]  
      data     n0_snow / 1.6e+7 /               ! [1/m^4]
      data     density_snow / 0.1E+3 /          ! [kg/m^3]  
      !data     diam_increment / 0.10 /          ! mm    
      data     diam_increment / 0.20 /          ! mm    
c      
c     Assign some useful constants
      wavel = 300./freqy
      
c
c     Begin by checking if hydrometeors of this species are present.
c     If not, set scattering parameters to zero and return.
c
      if(swc .lt. 0.0025) then
        ksca=0.
        asca=0.
        gsca=0.
        pbck=0.
        return
      endif

c
c     If hydrometeors are present, initialize the scattering parameters
c
      bext=0.
      bsca=0.
      bsym=0.
      bq11=0.

c
c     Loop over particle sizes:

c     increments of particle diameter are 0.10 mm; the particle
c     size distribution is expressed as a particle number density,
c     num, per diameter increment; the original psd is
c     n(D) = n0s * exp(-lam * D) where n is the number density
c     per diameter increment, n0s is the distribution intercept,
c     lam is the slope of the distribution of ln(n(D)), and D is
c     the particle diameter. It is assumed here that n0_snow 
c     and swc are prescribed, and lam is therefore constrained 
c     to yield the prescribed water content:
c     lwc = integral {n(D) * pi * (D**3) * density(D) * dD/6}
c     therefore:
c     lam=(n0s*pi*density_snow/lwc)**(0.25)

c
      do i=0,100 !200
        diam = 0.050 + diam_increment*float(i)
        x =  pi*diam/wavel 
        lam = (N0_snow*pi*density_snow/(swc*(1.e-3)))**(0.25)
        num = N0_snow*exp(-lam*diam*(1.e-3))
c
c       complex refractive index of snow
c
        call iceoptic(freqy,temp,epsreal,epsimag)
        eice = cmplx(epsreal,epsimag)
        eair = cmplx(1.0006,0.0)

c       calculate dielectric constant of snow as an ice matrix  
c       with air inclusions, using Maxwell-Garnett mixing for  
c       2-component media w. elliptical inclusions 
        fincl = 1. - (density_snow/densice)
        call mg_ellips(fincl, eice, eair, emg)
        cref = csqrt(emg)
c
c       call Mie program
c 
        call mie_sphere(x,cref,qsca,qext,asym,qbsca)

        qext = qext 
        qsca = qsca 
        asym = asym 
        qbsca = qbsca 
        
c
c       integrate over particle size distribution;

        bext=bext+num*qext*pi*0.25*diam*diam*diam_increment*1.e-6
        bsca=bsca+num*qsca*pi*0.25*diam*diam*diam_increment*1.e-6
        bsym=bsym+num*qsca*asym*pi*0.25*diam*diam*diam_increment*1.e-6
        bq11=bq11+num*qbsca*pi*0.25*diam*diam*diam_increment*1.e-6

      end do

c
c     check for distribution with very small extinction;
c     set parameters to zero to avoid numerical problems
      if( bext .gt. 1.e-9) then

        ksca=bext
        asca=bsca/bext
        gsca=bsym/bsca
        pbck=bq11/bsca

      else

        ksca=0.
        asca=0.
        gsca=0.
        pbck=0.

      end if

      return
      end subroutine mie_snow
      
! --------------------------------------------------------------------------------------

      SUBROUTINE MIE_SPHERE (X, MIN, QSCAT, QEXTI, ASYM, QBSCAT)
**
      implicit    none
      SAVE
**
**    Mie Routine P. Bauer 
**
      integer    limitx
      PARAMETER (LIMITX = 1500)

**
      REAL        X
      REAL        MR, MI, N1, N2
      REAL        QSCAT, QEXTI, QABSO, ASYM, QBSCAT
**
      REAL        RFAC1, RFAC2
      REAL        RHELP1(2), RHELP2(2)
**
      COMPLEX     M, MX, MIN
      COMPLEX     CHELP1, CHELP2, CFAC1, CFAC2, CBSCAT
**
      COMPLEX     DN(0:LIMITX), WN(-1:LIMITX)
      COMPLEX     AN(LIMITX), BN(LIMITX)
**
      INTEGER     NEND
      INTEGER     I100, I101
**
      EQUIVALENCE (CHELP1, RHELP1 (1))
      EQUIVALENCE (CHELP2, RHELP2 (1))
**
************************************************************************
**
      M      = CONJG (MIN)
      CHELP1 = M
      MR     =        RHELP1 (1)
      MI     = -1.0 * RHELP1 (2)
**      
      MX   = M  * X
      N1   = MR * X
      N2   = MI * X
**
      IF (X .LE. 20000.0) NEND = X + 4.00 * X ** (1.0 / 3.0) + 2.0
      IF (X .LE.  4200.0) NEND = X + 4.05 * X ** (1.0 / 3.0) + 2.0
      IF (X .LE.     8.0) NEND = X + 4.00 * X ** (1.0 / 3.0) + 1.0
      IF (NEND .LE.    5) NEND = 5
      IF (NEND .GT. LIMITX) NEND = LIMITX
**
      RFAC1      = SIN  (N1) * SIN  (N1) + SINH (N2) * SINH (N2)
      RHELP1 (1) = SIN  (N1) * COS  (N1) / RFAC1
      RHELP1 (2) = SINH (N2) * COSH (N2) / RFAC1
**
      DN (0) = CHELP1
**
      RHELP1 (1) =             COS (X)
      RHELP1 (2) = -1.0 E+00 * SIN (X)
      RHELP2 (1) =             SIN (X)
      RHELP2 (2) =             COS (X)
**
      WN (-1) = CHELP1
      WN ( 0) = CHELP2
**
      QEXTI  = 0.0
      QSCAT  = 0.0
      QBSCAT = 0.0
      QABSO  = 0.0
      ASYM   = 0.0 
      CBSCAT = CMPLX (0.0,0.0)
**
      DO 100 I100 = 1, NEND
         DN (I100) = -1.0 * I100 / MX
     +             +  1.0 / (I100 / MX - DN (I100 - 1))
         WN (I100) = WN (I100 - 1) * (2.0 * I100 - 1.0) / X
     +             - WN (I100 - 2)
**
         CFAC1 = DN (I100) / M + I100 / X
         CFAC2 = M * DN (I100) + I100 / X
**
         CHELP1 = WN (I100)
         CHELP2 = WN (I100 - 1)
**
         AN (I100) = (CFAC1 * RHELP1 (1) - RHELP2 (1))
     +             / (CFAC1 * CHELP1     - CHELP2    )
         BN (I100) = (CFAC2 * RHELP1 (1) - RHELP2 (1))
     +             / (CFAC2 * CHELP1     - CHELP2    )
**
         CHELP1 = AN (I100)
         CHELP2 = BN (I100)
**
         RFAC1 = RHELP1 (1) + RHELP2 (1)
         RFAC2 = CABS (AN (I100)) * CABS (AN (I100))
     +         + CABS (BN (I100)) * CABS (BN (I100))
**
         QEXTI  = QEXTI  + (2.0 * I100 + 1.0) * RFAC1
         QSCAT  = QSCAT  + (2.0 * I100 + 1.0) * RFAC2
         CBSCAT = CBSCAT + (2.0 * I100 + 1.0) * (-1.0) ** I100
     +          * (AN (I100) - BN (I100))
**
         IF (I100 .EQ. 1) GO TO 100
**
         CHELP1 = AN (I100 - 1) * CONJG (AN (I100))
     +          + BN (I100 - 1) * CONJG (BN (I100))
         CHELP2 = AN (I100 - 1) * CONJG (BN (I100 - 1))
**
         I101 = I100 - 1
         RFAC1  = I101 * (I101 + 2) / (I101 + 1.0)
         RFAC2  = (2.0 * I101 + 1.0) / (I101 * (I101 + 1.0))
**
         ASYM = ASYM + RFAC1 * RHELP1 (1) + RFAC2 * RHELP2 (1)
100   CONTINUE
**
      QEXTI  = QEXTI * 2.0 / (X * X)
      QSCAT  = QSCAT * 2.0 / (X * X)
      ASYM   = ASYM  * 4.0 / (X * X * QSCAT)
      QBSCAT = CABS (CBSCAT) * CABS (CBSCAT) / (X * X)
      IF (QSCAT .GT. QEXTI) QSCAT = QEXTI
**
      RETURN
      END SUBROUTINE MIE_SPHERE
      
! ---------------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  end module radtran
