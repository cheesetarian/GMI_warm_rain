      module read_l1c
       use pp_definition
       use pp_commoncode

! taken pp_read_GMI_L1C.F from GPM folks' GMI pp, made some light
! modifications for it to work with geos5 pp, D Duncan 1/11/16
            
       contains
       
!-------------------------------------------------------------------------

      subroutine read_L1C_GMI
      implicit none

#include 'TKHEADERS.h'
#include 'TK_1CGMI.h'
#include 'tkiofortdeclare.h'

      record /TKINFO/ hdfin
      record /L1CGMI_SWATHS/ swath
      common /L1CGMI_SWATHS/ swath

!--- Declaration of internal subroutine variables.

      integer :: ic, ipix, chan, j, i, qualhist(maxchans),idate
      logical :: fexists,misschan
      integer :: lastchar,ret=0,lin=1

      character (len=100) ::   file

!---check if file exists -------------------
      
      inquire(file=input_file,exist = fexists)
      if(fexists .eq. .false.) call reprt_error(6)           !input file does not exist
      
!--- open and start reading HDF file

      ret = TKopen(input_file,"1CGMI",TKREAD,"HDF5",jobid,hdfin, 1)
      if(ret .ne. 0) call reprt_error(8)                     !error in TKopen             

!--- Read global attribute metadata fields

      ret = TKgetmetaint(hdfin,"S2_SwathHeader","NumberScansGranule",
     >                  nscans)     
      if(ret .ne. 0) call reprt_error(9)                     !error in Tkgetmetaint1
      if(nscans .eq. 0) call reprt_error(7)

      ret = TKgetmetaint(hdfin,"S2_SwathHeader","NumberPixels", 
     >                  npix)     
      if(ret .ne. 0) call reprt_error(10)                     !error in Tkgetmetaint1

!--- find the orbit number = granule,  It's in the global attributes :FileHeader
!--- SatelliteName=GPM, InstrumentName=GMI, GranuleNumber=000914
!--- FileName=1C.GPM.GMI.XCAL2012-N.20120702-S020330-E033606.000914.V00A.HDF5

      ret = TKgetmetaint(hdfin,"FileHeader","GranuleNumber",granule)     
      if(ret .ne. 0) call reprt_error(11)                     !error in Tkgetmetaint1
      
      ret = TKgetmetastring(hdfin,"FileHeader","SatelliteName",
     >                  csatellite)     
      if(ret .ne. 0) call reprt_error(12)                     !error in Tkgetmetaint1
      ret = TKgetmetastring(hdfin,"FileHeader","InstrumentName",
     >                  csensor)     
      if(ret .ne. 0) call reprt_error(12)                     !error in Tkgetmetaint1

      ret = TKgetmetastring(hdfin,"FileHeader","FileName",
     >                 radfile)     
      if(ret .ne. 0) call reprt_error(12)                     !error in Tkgetmetaint1

      write(*,*)'  granule        : ',granule
      write(*,*)'  nscans         : ',nscans
      write(*,*)'  npix           : ',npix
      write(*,*)'  nchans         : ',nchans
      write(*,*)'  maxchans       : ',maxchans

!--- Allocate data arrays
 
      allocate(sclat(nscans))
      allocate(sclon(nscans))
      allocate(scalt(nscans))
      allocate(scorient(nscans))
      allocate(stdtime(nscans,6))

      allocate(lat(npix, nscans))
      allocate(lon(npix, nscans))

      allocate(Tb(npix, nscans, maxchans))
      allocate(eia(npix, nscans, maxchans))
      allocate(sglinta(npix,nscans))
      allocate(qualflag(npix,nscans))

! other allocations done in l1r routine, to match:
      allocate(sunglint(npix,nscans),satAZ(npix,nscans))
      allocate(tai93time(nscans),lo_flag(npix,nscans,4))

      satAZ(:,:) = -999.9
      qualflag = 0
      
!--- trim the sensor name and satellite name

      j = len_trim(csensor)
      do i = 1,j
        if(iachar(csensor(i:i)) .eq. 0)   exit
      enddo
      sensor_name = csensor(1:i-1)
              
      j = len_trim(csatellite)
      do i = 1,j
        if(iachar(csatellite(i:i)) .eq. 0)   exit
      enddo
      sat_name = csatellite(1:i-1)
      
      radfile = trim(radfile)

      write(*,'(a,a)')'   sat_name       : ',sat_name
      write(*,'(a,a)')'   sensor_name    : ',sensor_name
      write(*,'(a,a)')'   radiometer file: ',radfile
      
!--- Read variables

      do ic=1,nscans
        ret = TKreadScan(hdfin,swath)
        if(ret .ne. 0) call reprt_error(13)                      !error in TkreadScan

        sclat(lin)   = swath%S2%SCstatus%SClatitude 
        sclon(lin)   = swath%S2%SCstatus%SClongitude
        scalt(lin)   = swath%S2%SCstatus%SCaltitude      !in KM
        scorient(lin)= swath%S1%SCstatus%SCorientation

        do ipix = 1, npix 
          lat(ipix,lin) = swath%S1%Latitude(ipix)
          lon(ipix,lin) = swath%S1%Longitude(ipix)
          sunglint(ipix,lin)= swath%S1%sunGlintAngle(1,ipix)
          sglinta(ipix,lin)= int(swath%S1%sunGlintAngle(1,ipix))

          if(swath%S1%Quality(ipix) .ne. 0  .or.             ! set quality variable to 
     >	     swath%S2%Quality(ipix) .ne. 0) then             ! the value in S1, if 0
             if(swath%S1%Quality(ipix) .ne. 0) then	         ! set then check S2 and set
	         qualflag(ipix,lin) = swath%S1%quality(ipix)
	      else
	         qualflag(ipix,lin) = swath%S2%quality(ipix)
             endif
	  endif
	  
c      write(6,*) lin,ipix,swath%S1%quality(ipix),swath%S2%quality(ipix),
c     >           quality(ipix,lin)	  
          Tb(ipix,lin,1)    = swath%S1%Tc(1,ipix)         !10V
          Tb(ipix,lin,2)    = swath%S1%Tc(2,ipix)         !10H
          Tb(ipix,lin,3)    = swath%S1%Tc(3,ipix)         !19V
          Tb(ipix,lin,4)    = swath%S1%Tc(4,ipix)         !19H
          Tb(ipix,lin,5)    = swath%S1%Tc(5,ipix)         !23V
          Tb(ipix,lin,6)    = miss_flt
          Tb(ipix,lin,7)    = swath%S1%Tc(6,ipix)         !37V
          Tb(ipix,lin,8)    = swath%S1%Tc(7,ipix)         !37H	  	  
          Tb(ipix,lin,9)    = swath%S1%Tc(8,ipix)         !89V
          Tb(ipix,lin,10)   = swath%S1%Tc(9,ipix)         !89H

          Tb(ipix,lin,11)   = swath%S2%Tc(1,ipix)         !166V
          Tb(ipix,lin,12)   = swath%S2%Tc(2,ipix)         !166H
          Tb(ipix,lin,13)   = miss_flt                    !no 183v +/-1
          Tb(ipix,lin,14)   = swath%S2%Tc(3,ipix)         !183v +/-3
          Tb(ipix,lin,15)   = swath%S2%Tc(4,ipix)         !183v +/-7

          eia(ipix,lin,1:5)  = swath%S1%incidenceAngle(1,ipix)   !10V,10H,19V,19H,23V
          eia(ipix,lin,6)    = miss_flt                          !23H
          eia(ipix,lin,6:10) = swath%S1%incidenceAngle(1,ipix)   !37V,37H,89V,89H
          eia(ipix,lin,11:12)= swath%S2%incidenceAngle(1,ipix)   !166V, 166H
          eia(ipix,lin,13)   = miss_flt                          !183+/-1
          eia(ipix,lin,14:15)= swath%S2%incidenceAngle(1,ipix)   !186,190
	  
        enddo

        stdtime(lin,1)    = swath%S1%ScanTime%Year       !YYYY
        stdtime(lin,2)    = swath%S1%ScanTime%Month      !MM
        stdtime(lin,3)    = swath%S1%ScanTime%DayOfMonth !DD
        stdtime(lin,4)    = swath%S1%ScanTime%Hour       !HH 
        stdtime(lin,5)    = swath%S1%ScanTime%Minute     !MM
        stdtime(lin,6)    = swath%S1%ScanTime%Second     !SS
        tai93time(lin)   = -999.9 ! missing from GMI (for now)

        lin = lin+1
      enddo
      write(*,'(a,i4,1x,5(i2,1x))')'   first time     : ', stdtime(1,:)


!--- quality control of the Tbs in case missing value .ne. miss_flt
!--- and bounds check on Lat/Lon

      do lin = 1,nscans
        do ipix = 1, npix
          do chan = 1,maxchans
            if(Tb(ipix,lin,chan) .lt. 0) Tb(ipix,lin,chan) = miss_flt
          enddo
	  if(lat(ipix,lin).gt.90  .or. lat(ipix,lin).lt.-90  .or.
     >       lon(ipix,lin).gt.180 .or. lon(ipix,lin).lt.-180) then
                lat(ipix,lin) = miss_flt
		lon(ipix,lin) = miss_flt
          endif	  
        enddo
      enddo


c      do i = 1,10
c        write(6,'(I3,1x,26F7.2)')i,lat(91:115,i)
c	write(6,'(I3,1x,26F7.2)')i,lon(91:115,i)	
c      enddo
c      stop

!--- check L1C quality flag for bad data.   Don't check flags (-4 and -5) for missing chans
!--- also check pixel's Tbs against chan_avail.  If chan should be there but isn't, then
!--- set all Tbs to missing.

      do lin=1,nscans
        do ipix = 1,npix

c	  if((qualflag(ipix,lin) .eq. -1 .or.       !Data missing from file
c     >        qualflag(ipix,lin) .eq. -2 .or.      !one or more Tbs out of range 
c     >        qualflag(ipix,lin) .eq. -3 .or.      !Geolocation error
cc     >       qualflag(ipix,lin) .eq. -4.or.       !One Channel missing
cc     >       qualflag(ipix,lin) .eq. -5.or.       !More than one Channel missing
c     >        qualflag(ipix,lin) .eq. -6 .or.      !Lat/lon error
c     >        qualflag(ipix,lin) .eq. -7) .and.    !'non-normal' status error
c     >        qualflag(ipix,lin) .ne. -10 ) then   ! L1CR product - high freq channels no match

          if(qualflag(ipix,lin) .le. -1 .and. 
     >       qualflag(ipix,lin) .gt. -10) then 
                  Tb(ipix,lin,:)  = miss_flt 
                  eia(ipix,lin,:) = miss_flt		 
c                  write(6,*) lin,ipix,qualflag(ipix,lin),Tb(ipix,lin,:)
          endif 

        enddo
      enddo

!--- special to set frequencies to missing when channels fail        

c      idate = stdtime(1,1)*10000 + stdtime(1,2)*100 + stdtime(1,3)
c      write(6,*)' idate = ', idate       
c      if(idate .ge. 20200101) then  !example: 89H failure
c          chan_freq(10)  = miss_flt
c	   chan_avail(10) = 0
c          write(6,*)' setting to missing for GMI - 89H GHZ channels'
c      endif       
      
!--- histogram the number of pixels per channel, if all pixels unexpectedly missing, then
!--- call error(4), which will stop the proprocessor and generate operator error.

      qualhist = 0
      do chan = 1,maxchans
        do lin = 1,nscans
          do ipix = 1, npix	  
            if(Tb(ipix,lin,chan) .ne. miss_flt) then
                qualhist(chan) = qualhist(chan) + 1
            endif
          enddo
        enddo	
      enddo 
      misschan = .false.
      do chan = 1,maxchans
        if(qualhist(chan) .eq. 0 .and. chan_avail(chan)) then	           	    
	    write(*,*)' All Tbs unexpectedly missing, channel : ',chan
	    misschan = .true.
	endif
      enddo
      if(misschan) call reprt_error(4)

!--- Close the input data file

      ret = TKclose(hdfin)
      if(nscans .eq. 0) call reprt_error(14)

      return
      end  subroutine read_L1C_GMI
      
!---------------------------------------------------------------------------------

      end module read_l1c
