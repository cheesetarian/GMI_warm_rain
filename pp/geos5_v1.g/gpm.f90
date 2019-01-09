      module gpm
      use pp_definition
      use pp_commoncode

      contains

!------------------------------------------------------------------------------

      subroutine add_ancillary
       implicit none

!--- this routine will add the other required ancillary data to the pixel
!--- this includes model parameters(3), and sfccode

       integer           :: ipix, iscan, istat, ilat,ilon
       integer           :: time(3),i,j,k, hr, ichr, icdy
       
       real,parameter    :: srinc    = .03125, sppd = 16.           !srinc= 1/2 spacing
       integer,parameter :: nlon     = 5760,  nlat = 2880           ! for geo area
       
       real     :: newlon, tlon       
       integer  :: nglatlon,ngLR, ngTwb, ngtcwv, ngskint, ngT2m, tpixs
       
       !integer,parameter :: ysizeemis= 359, xsizeemis=720           !emis = 1/2 degree grid
       !real              :: eppd = 2.0, erinc = 0.25                !erinc = 1/2 spacing
       !integer           :: elat, elon, nmon, cnt(4)=0              ! for geo area
       
       dir_ancil='../input/'
       !dir_ancil='/tdata1/dduncan/transfer/'
       !dir_ancil='/home/drandel/gpm/ppancillary/'
       
!--- add CSU surface Codes----------------------
       write(*,*)' starting surface classification map for sensor: ',&    !make the CSU surface
                trim(sensor_name)                                       !classification map
       
       call GPM_make_surfacemap(istat)
       if(istat .ne. 0) call reprt_error(24)
              
       write(*,*)'  adding CSU surface codes'
       do ipix = 1,npix
         do iscan = 1, nscans
	       newlon = lon(ipix,iscan)
	       if(newlon .gt. 180) newlon = newlon - 360.
           
	   if(lat(ipix,iscan) .gt. -91.0 .and. newlon .ge.-180.0) then !check for good lat/lon
              ilat= nint(sppd* ( 90.0 + srinc - lat(ipix,iscan)))    !lat-lons in
              ilon= nint(sppd* ( 180.0 + srinc + newlon))             !-180:180 space	       
              if(ilon .eq. nlon+1) ilon = 1  
              if(ilat .eq. nlat+1) ilat = nlat
        
              if(ilon .gt. nlon .or. ilat .gt. nlat) then
                 write(*,*)' Lat/Lon OB :', ilat, ilon,lat(ipix,iscan),lon(ipix,iscan)                 
                 call reprt_error(25)
              else       
                 sfccode(ipix,iscan) = sfcgrid(ilon,ilat)
              endif      
           else
              sfccode(ipix,iscan) = miss_byt
           endif

           ! add-ons from GEOS5 ---
           if(sic(ipix,iscan)<5 .and.lo_flag(ipix,iscan,rz)<=2.and.snowci(ipix,iscan)<1) then
              sfccode(ipix,iscan) = 1 ! ocean (0-2% land okay)
           endif
           if(lo_flag(ipix,iscan,rz)<=2.and.skint(ipix,iscan)<268.and.sic(ipix,iscan)<5)then
              sfccode(ipix,iscan) = 2 !for frozen inland lakes -- sea ice, say
           endif
           if(sic(ipix,iscan)>=5 .and.lo_flag(ipix,iscan,rz)<=50)then
              sfccode(ipix,iscan) = 2 !GPM sea-ice
           endif
           if(sic(ipix,iscan)>=5 .and.lo_flag(ipix,iscan,rz)>50)then
              sfccode(ipix,iscan) =  10 !GPM light snow cover
           endif
           if(snowci(ipix,iscan)>=1 .and. lo_flag(ipix,iscan,rz)>50) then !snow and land...
             sfccode(ipix,iscan) = 11-snowci(ipix,iscan) !set to snow class
           endif
           if((snowci(ipix,iscan)>=1.or.sic(ipix,iscan)>=5).and.sfccode(ipix,iscan)>=13) then  !any snow or sic over coast
             sfccode(ipix,iscan) = 11 !set to min snow
           endif
           if(sic(ipix,iscan)>=5 .and. (lo_flag(ipix,iscan,rz)<=50.and.lo_flag(ipix,iscan,rz)>2)) then
             sfccode(ipix,iscan) = 14 ! sea-ice boundary
           endif
           if(sfccode(ipix,iscan)>14) then
             sfccode(ipix,iscan)=miss_byt
           endif
        ! force greenland and antarctica (from gpm code)
              if(lat(ipix,iscan) <= -60) then              !setsnow cover to 11 over
                !if(sfccode(ipix,iscan)>= 8 .and.sfccode(ipix,iscan)<= 10) then ! snow over antarctica
                if(lo_flag(ipix,iscan,rz)>95) then ! snow over antarctica
                    sfccode(ipix,iscan) = 8
                endif
              endif
              if(lat(ipix,iscan) >= 60.and.lat(ipix,iscan)<=82 .and. & !set snow cover to 11 for surface Greenland
                  lon(ipix,iscan)>= -66 .and. lon(ipix,iscan) <= -23 ) then
                !if(sfccode(ipix,iscan) >= 8  .and. sfccode(ipix,iscan) <= 10) then
                if(lo_flag(ipix,iscan,rz)>95) then
                  sfccode(ipix,iscan) = 8
                endif
              endif

         enddo !scan
       enddo !pixel
             

      return
      end subroutine add_ancillary

!-----------------------------------------------------------------------------    

      subroutine GPM_read_emissclass_grds(istat)
       implicit none
       
!--- this routine reads in the emissivity class grids converted from the Filipe ascii file
!--- all 12 fields are put into the emissgrd array which is defined in definition
!--- and therefore common and not needed to be passed

       integer :: i,j,k,mon,ii,jj,istat, ios, time(6), lun,ifound
       character(len=96) :: infile
       character(len=2)  :: cmon
       integer(2)        :: emissfix(720,359)
       istat = 0
       
!--- define input filename from month       

       write(6,'(a,a)')'   reading emiss class files'
       do i = 1,12
         write(cmon,'(i2.2)') i
         infile= trim(dir_ancil) // 'emiss_class_' // cmon // '.dat'
!	 write(6,*)i, trim(infile)
       
!--- open emiss class file

         call gprof_lun(lun)
         open(unit=lun, file=infile, access='stream', status='old',iostat=ios)
         if(ios .ne. 0) then 
             write(6,*)' error opening emiss class file'
	     istat = -1
	     return
         endif

!--- read in the emissivity field for month i
      
         read(lun,iostat=ios) emissgrd(:,:,i)
         if(ios .ne. 0) then
           write(6,*)' error reading emiss class file'
	   istat = -2
         endif
         close(lun)
       enddo

!--- fill some missing ocean/water values with closest land emiss code,
!interate till more filled

         emissfix = emissgrd(:,:,mon)
         do k = 0,3
           do j = 2,358
             do i = 2,719
               if(emissgrd(i,j,mon) .eq. 0) then
                   ifound = 0
                   do ii = -1,1
                     do jj = -1,1
                       if(ifound .eq. 0) then
                         if(emissgrd(i+ii,j+jj,mon) .ne. 0) then
                           emissfix(i,j) = emissgrd(i+ii,j+jj,mon)
                           ifound = 1
                         endif
                       endif
                     enddo
                   enddo
               endif
             enddo
           enddo
           emissgrd(:,:,mon) = emissfix
         enddo


      return
      end subroutine GPM_read_emissclass_grds

!-------------------------------------------------------------------------
      
      subroutine GPM_read_landmask(lmskfile, lsmask, istat)
       implicit none         
       
       character(len=96) :: lmskfile
       integer(1) :: lsmask(5760,2880)         !max size for lmaskres is 16 ppd

       integer(2)       :: ixsize,iysize,irecl,lmaskres,test(5760,2880) 
       character(len=2) :: cres, sres
       integer          :: lun, ios, i,j, istat

       istat = 0
 
!--- get Mask grid resolution from filename

       read(lmskfile(12:13),'(i2)',iostat=ios) lmaskres
       if(ios .ne. 0) then
          write(6,*)' error decoding lmaskres '
	  istat = -1
	  return
       endif
       write(*,*)'   landmask resolution : ',lmaskres
       
!--- set size of landmask file       
       
       ixsize = lmaskres * 360 
       iysize = lmaskres * 180
       irecl = ixsize / 4
       
!       write(6,*)' ixsize,iysize = ', ixsize,iysize

!--- open landmask file

       write(*,*)'   reading landmask    : ',trim(dir_ancil) //trim(lmskfile)       
       call gprof_lun(lun)
       open(unit=lun,file= trim(dir_ancil)//trim(lmskfile),form='unformatted', status='old', recl= irecl,access='direct',iostat=ios)
       if(ios.ne.0) then
            write(6,*)' error opening landmask file'
	    istat = -2
	    return
       endif

!--- read in land-sea mask file

       do j = 1,iysize
         read(lun,rec=j,iostat=ios)(lsmask(i,j),i=1,ixsize)
         if(ios.ne.0) then
	      write(6,*)' error reading landmask file'
	      istat = -3
	      return
	 endif
       enddo
       close(lun)            

       return
      end subroutine GPM_read_landmask
!-----------------------------------------------------------------------	

      subroutine GPM_make_surfacemap(istat)
       implicit none

!---   this routine makes the CSU surface map classifications.  It's a combination of
!---   the MODIS surface landmask, autosnow land and sea-ice snow cover, and a monthly 
!---   climatology of surface emissivity.

!---   the AUTOSNOW data can be read from the separate hemispheric files using the call:
!---    GPM_read_autosnow   call,  OR

!---   the AUTOSNOW data can be read from the global grids using the call :
!---    GPM_read_autosnow_global


       character(len=10) ::sensor
       character(len=128):: landmaskfile

       integer    :: i,j,k,ios,jj,ii,off,ix,iy,ie,is,je,js,isr,nmon,mon
       real       :: lat,lon
       integer(1) :: lsmaskgrd(5760,2880), tsnowgrid(9000,4500)  !16 degree mask, temp file
       integer    :: istat,time(6),indexexp
       integer    :: icnt(10)

       real,parameter    :: emisres = 8.0, snowres= 0.64
       integer,parameter :: ysizeemis= 359, xsizeemis=720
       integer,parameter :: ysizesnow= 4500, xsizesnow=9000
       integer           :: modeli, modelj
       integer           :: emisi,emisj, snowi,snowj

       integer           :: modellat, modellon               ! for gridpoint location
                                                             ! in high res T2m array
       !sfcgrid = 0
       icnt = 0
       
!---  set the time to first scan time, and sensor

       time = stdtime(1,:)
       mon = time(2)
       sensor = trim(sensor_name)
        print*,'SENSOR NAME: ',sensor

!---  set landmask by sensor
       
       if(trim(sensor) .eq. 'AMSRE') landmaskfile='landmask51_16.bin'
       if(trim(sensor) .eq. 'AMSR2') landmaskfile='landmask42_16.bin'
       if(trim(sensor) .eq. 'GMI')   landmaskfile='landmask32_16.bin'
       if(trim(sensor) .eq. 'SSMI ') landmaskfile='landmask69_16.bin'
       if(trim(sensor) .eq. 'SSMIS') landmaskfile='landmask74_16.bin'
       if(trim(sensor) .eq. 'TMI')   then
           if(prfdbfile(1:5) .eq. 'TMIPR') then
	         landmaskfile='landmask60_16.bin'
           elseif(prfdbfile(1:5) .eq. 'TMIPO') then
	         landmaskfile='landmask68_16.bin'
	   endif
       endif
       
       if(trim(sensor) .eq. 'MHS1')  landmaskfile='landmask20_16.bin'
       if(trim(sensor) .eq. 'MHS2')  landmaskfile='landmask23_16.bin'
       if(trim(sensor) .eq. 'MHS3')  landmaskfile='landmask28_16.bin'
       if(trim(sensor) .eq. 'MHS4')  landmaskfile='landmask32_16.bin'
       if(trim(sensor) .eq. 'MHS5')  landmaskfile='landmask37_16.bin'        

!---  read in the landmask, which is static over the year

       write(6,*)'  calling GPM_read_landmask for: ',trim(landmaskfile)
       call GPM_read_landmask(landmaskfile,lsmaskgrd,istat)
       if(istat .ne. 0) call reprt_error(21)
              
!---  read in the emissclass grids, and the autosnow grid
       
       write(6,*)'  calling GPM_read_emissclass_grd'       
       call GPM_read_emissclass_grds(istat)
       if(istat .ne. 0) call reprt_error(22)

!       write(6,*)'  calling GPM_read_autosnow_global'       
!       call GPM_read_autosnow_global(time(1:3), istat)
       
!       write(6,*)'  calling GPM_read_autosnow for hemispheric files'      
!       call GPM_read_autosnow(time(1:3), istat)

!       if(istat .ne. 0) call reprt_error(23)

!--- Insert landmask into final array

!       do j = 1,nscans
!        do i = 1,npixl
!         if(lo_flag(i,j,rz).ge.0.and.lo_flag(i,j,rz).le.2) then
!           sfccode(i,j)=10 !(0-2% = ocean)
!           icnt(1)=icnt(1)+1
!         endif
!         if(lo_flag(i,j,rz).gt.2.and.lo_flag(i,j,rz).le.95) then
!           sfccode(i,j)=30 ! (3-95% = coast)
!           icnt(3)=icnt(3)+1
!         endif
!         if(lo_flag(i,j,rz).gt.95) then
!           sfccode(i,j)=20  ! (96-100% = land)
!           icnt(2)=icnt(2)+1
!         endif
!        enddo
!       enddo

        ! replaced by using lo_flag from l1r data
!       write(6,*)'  starting landmask insertion'
!       do j = 1,2880
!         do i = 1,5760 
!	  if(lsmaskgrd(i,j) .ge. 0 .and. &      !0-2% land = ocean exit
!             lsmaskgrd(i,j) .le. 2) then
!	           sfcgrid(i,j) = 10 
!		   icnt(1) = icnt(1) + 1
!	  elseif(lsmaskgrd(i,j) .gt. 2 .and. &  !2-95%     = mixed pixel (coast)
!                 lsmaskgrd(i,j) .le. 95) then
!                   sfcgrid(i,j) = 30 
!		   icnt(3) = icnt(3) + 1
!	  elseif(lsmaskgrd(i,j) .gt. 95) then   !95-100%   = land
!	           sfcgrid(i,j) = 20 
!		   icnt(2) = icnt(2) + 1
!	  endif
!         enddo
!       enddo

!--- Insert the emiss class and autosnow values into the sfcgrid

       write(6,*)'  starting autosnow and emiss class insertion'
       do j = 1,2880
         emisj = int(1+ float(j-1)/emisres)
         snowj = int(1+ float(j-1)/snowres)
	 
	 if(emisj .eq. 360) emisj = 359
	 if(emisj .gt. ysizeemis) then
	    write(6,*)' emisj OB ', emisj
	    cycle
	 endif
	 if(snowj .gt. ysizesnow) then
	    write(6,*)' snowj OB ', snowj
	    cycle
	 endif
	 
         do i = 1,5760	   
           emisi = int(1+ float(i-1)/emisres)	   
           snowi = int(1+ float(i-1)/snowres)   
              
	   if(emisi .gt. xsizeemis) then
	      write(6,*)' emisi OB ',emisi
	      cycle
	   endif
	   if(snowi .gt. xsizesnow) then
	      write(6,*)' snowi OB ',snowi
	      cycle
	   endif

	   if(sfcgrid(i,j) .eq. 10) sfcgrid(i,j) = 1                       !ocean=1

           !if(sfcgrid(i,j) .eq. 1 .and. snowgrid(snowi,snowj).eq.2) then   !ocean, snow over land
           if(sfcgrid(i,j) .eq. 1 .and. snowgrid(snowi,snowj).ge.1) then   !ocean, snow over land
	        sfcgrid(i,j) = 2                                           !set to sea-ice
           endif
     
           ! sea ice taken care of from geos5 data
           !if(snowgrid(snowi,snowj) .eq. 3) sfcgrid(i,j) = 2               !sea ice
           !if(snowgrid(snowi,snowj) .eq. 5) sfcgrid(i,j) = 14              !sea-ice edge boundary
	     
           if(sfcgrid(i,j).eq.20) then                                     !snowcover over land 

!              if(snowgrid(snowi,snowj) .eq. 2 .or. 
!     >           snowgrid(snowi,snowj) .eq. 3) then
              if(snowgrid(snowi,snowj) >=1) then
	           if(emissgrd(emisi,emisj,mon) .ge. 6 .and.  &             !snow classes 6-9
                      emissgrd(emisi,emisj,mon) .le. 9) then               !and snowgrid yes
                     sfcgrid(i,j) = emissgrd(emisi,emisj,mon) + 2 
                   endif
	      
	           if((emissgrd(emisi,emisj,mon) .ge.  1 .and. &            !no snow emis classes
                       emissgrd(emisi,emisj,mon) .le.  5) .or. &            !and snowgrid=yes
                       emissgrd(emisi,emisj,mon) .eq. 10) then
		           sfcgrid(i,j) = 10
                   endif
		   
		   if(emissgrd(emisi,emisj,mon) .eq. 0) then                !antarctica, greenland
		      sfcgrid(i,j) = 8                                      !to max snow
		   endif	      
              else 
	           if(emissgrd(emisi,emisj,mon) .ge. 1 .and.   &           !vegetation 1-5  & no snow
                      emissgrd(emisi,emisj,mon) .le. 5) then               !sets sfcgrid to 3-8
	               sfcgrid(i,j) = emissgrd(emisi,emisj,mon) + 2          
	             !print*,sfcgrid(i,j),emissgrd(emisi,emisj,mon),i,j
                   endif
		  
		   if(emissgrd(emisi,emisj,mon) .ge. 6 .and.  &            !snow classes 6-9 & no snow
                      emissgrd(emisi,emisj,mon) .le. 9) then               !find latest no snow class                      
		       
		       do k = 1,11
	                 nmon = time(2) - k
		         if(nmon .le. 0) nmon = nmon + 12		       
		         if(emissgrd(emisi,emisj,nmon) .lt. 6) then
		            sfcgrid(i,j) = emissgrd(emisi,emisj,nmon)+2
	                    exit
	                 endif
			 if(emissgrd(emisi,emisj,nmon) .eq. 10) then
		            sfcgrid(i,j) = 12
	                    exit
	                 endif 
                       enddo		       
		       if(sfcgrid(i,j) .eq. 20) then                      !no unfrozen class found
			    sfcgrid(i,j) = 11                             !in all months, set to
		       endif                                              !min snow		   
		   endif 	  
		                
		   if(emissgrd(emisi,emisj,mon) .eq. 10) then           !standing water emiss class 10
		      sfcgrid(i,j) = 12                                 !sets sfcgrid to 12
                   endif
		   
		   if(emissgrd(emisi,emisj,mon) .eq. 0) then            !land class missing emis
		      sfcgrid(i,j) = 1                                  !set to ocean
		   endif
		   
              endif
	      
	      if(sfcgrid(i,j) .eq. 20) then
	        write(6,*) emisi,emisj,emissgrd(emisi,emisj,mon)
		sfcgrid(i,j) = 12                                 !if any more missing set to
	      endif                                               !inland water
	      	      	      
	   endif
          
	   if(sfcgrid(i,j) .eq. 30) sfcgrid(i,j) = 13                   !coast           

!	   if(sfcgrid(i,j).eq.13 .and. snowgrid(snowi,snowj).eq.2)then  !coast with land snow
	   if(sfcgrid(i,j).eq.13 .and. snowgrid(snowi,snowj).ge.1)then  !coast with land snow
                sfcgrid(i,j) = 9
	   endif           
          
           !if(sfcgrid(i,j).eq.13 .and. snowgrid(snowi,snowj).eq.3) then  !coast with sea-ice 
	   !     sfcgrid(i,j) = 2                                
	   !endif

 	   	   
	 enddo  !lon
       enddo   !lat
 

        write(6,*)' num ocean points in landmask  : ',icnt(1)
        write(6,*)' num land  points in landmask  : ',icnt(2)
        write(6,*)' num coast points in landmask  : ',icnt(3)    
        write(6,*)' num points ocean/ice bndry    : ',icnt(4)


      return
      end subroutine GPM_make_surfacemap

!--------------------------------------------------------------------------
!----------------------------------------------------------------------------------

       end module gpm
       
