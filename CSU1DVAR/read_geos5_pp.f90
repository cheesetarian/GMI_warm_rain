      subroutine read_geos5_pp(inp_file,xtime) 

       ! new geos5 preprocessor input 
       use define_intam2
       use GPM_arraydef_D
       use GPM_util_procedures!_D
       use define_csu1dvar
       use csu1dvar_subr


       character(len=100) :: inp_file

       integer :: iscan,ipix,reccnt=0,i,j,ios,icnt(maxchans)=0,rlun
       integer :: ic, idate, xtime
       character(len=128) :: blank = ' '      
       real,parameter  :: Tbmin=50.0, Tbmax=325.0
       logical         :: igood
       
!---  input orbit header structure
       
       type :: OrbitHdr              
           character(len=12) :: psatellite   !orbit header 440bytes
           character(len=12) :: psensor
           character(len=12) :: ppp_version  
           character(len=128):: radfile
           character(len=128):: pdbfile
           character(len=128):: calfile
           integer           :: granule
           integer           :: nscans
           integer           :: npixels
           integer           :: nchans
           real              :: chan_freqs(maxchans)
           character(len=10) :: convolve
           character(len=30) :: comment
       end type OrbitHdr
          
!---  input time structure

       type :: Date         
           integer(2):: year
           integer(2):: month
           integer(2):: day
           integer(2):: hour
           integer(2):: minute
           integer(2):: second
       end type Date    

!---  input scan header structure
       
       type :: ScanHdr
           type(Date)         :: ScanDate
           real               :: Sclat
           real               :: Sclon
           real               :: Scalt
           real*8             :: Scorient
           real*8             :: TAI93 !time in s since 1993
       end type ScanHdr

!---  input pixel structure

       type :: DataRec     
           real     :: latitude
           real     :: longitude      
           real     :: Tbs(maxchans)       !Tb channels       
           real     :: eia(maxchans)       !earth incident angle
           real     :: Twb                 !Wet Bulb Temperature
           real     :: tcwv 
           real     :: skint               !skin temperature
           real     :: T2m                 !2 meter temperature index
           real     :: azimuth             !azimuth angle
           real     :: windmag             !magnitude of wind speed at 10m
           real     :: slp                 !sea level pressure
           integer*2          :: tprof(nz)        !t profile (K*100)
           integer*2          :: wvprof(nz)       !wv mixratio (g/kg*1000)
           integer*2          :: aight(nz)       !height (m)
           integer*2          :: qualflag         !quality flag from L1R
           integer*2          :: winddir          !direction, deg from N
           integer*1          :: lo_pct(4) ! land pct at each FOV size
           integer*1          :: sunglint_angle
           integer*1          :: surface_type_index
           integer*1          :: snow_cover_index
           integer*1          :: orolift_index  
       end type DataRec
                
!---  assign structure names
       
       type(OrbitHdr) :: orbhdr
       type(ScanHdr)  :: scnhdr
       type(DataRec)  :: px



!---  open input preprocessor file

       call gprof_lun(rlun)
       open(unit=rlun,file=trim(inp_file),access='stream',status='old',iostat=ios)
        if(ios .ne. 0) then
          write(*,*),' Error opening pp file'
          stop
        endif


!---  read orbit header

       read(rlun, iostat=ios)  orbhdr       !read out orbit header
       !if(ios .ne. 0) call GPM_reprt_err(11,rlun,blank)
!reading preproc orbhdr

!---  write orbit header information

       !write(log_lun,'(a,a)')' ORBITAL HEADER INFO'
       !write(*,*)' Radiometer file= ',trim(orbhdr%radfile)
       !write(log_lun,'(a,a)')' profile DB file= ',trim(orbhdr%pdbfile)
       !write(*,*)' Calibrtion file= ',trim(orbhdr%calfile)
       !write(log_lun,'(a,a)')' comment       = ',trim(orbhdr%comment)       
       !write(*,*)' PP version    = ',trim(orbhdr%pp_version)
       !write(*,*)' satellite     = ',trim(orbhdr%psatellite)
       !write(*,*)' sensor        = ',trim(orbhdr%psensor)
       !write(log_lun,'(a,i7)')' granule num  = ',orbhdr%granule
       !write(*,*)' number scans = ',orbhdr%nscans
       !write(*,*)' pixels/scan  = ',orbhdr%npixels       
       !write(*,*)' number chans = ',orbhdr%nchans
       !write(*,*) ' chan freqs = ',orbhdr%chan_freqs
!--- assign these for global variables

       !satellite   = orbhdr%psatellite
       sat_name    = orbhdr%psatellite
       !sensor      = orbhdr%psensor
       sensor_name = orbhdr%psensor
       if(sensor_name.eq.'AMSR2') then
         satcode='AM2'
         nretvar = 2 + npc
         if(retrievesst) nretvar=3+npc
       endif
       if(sensor_name.eq.'GMI') then
         satcode='GMI'
         nretvar = 2 + npc
        !print*,'nretvar: ',nretvar
         ! no sst retrieval override, for now
       endif
       !ppversion   = orbhdr%ppp_version
       convolve    = orbhdr%convolve !trim("None")
       !orig_file   = orbhdr%radfile
       radfile     = orbhdr%radfile
       l1fi = trim(radfile)
       pdbfile     = orbhdr%pdbfile
       cal_file    = orbhdr%calfile
       npix        = orbhdr%npixels          
       npixs       = orbhdr%npixels          
       nscans      = orbhdr%nscans
       granule     = orbhdr%granule
       nch         = orbhdr%nchans
       !nchannels   = orbhdr%nchans
        
       if(nch==13) nfreq=8 ! for eddington code
              
!--- define channel availability array
  
       chan_avail = 0
       do i = 1,maxchans
         if(orbhdr%chan_freqs(i) .gt. 0) chan_avail(i) = 1
       enddo 
       !if(nch .ne. sum(chan_avail)) then
       !  if(nch.eq.9) chan_avail=(/1,1,1,1,1,0,1,1,1,1,0,0,0,0,0/)
       !  if(nch.eq.7) chan_avail=(/0,0,1,1,1,0,1,1,1,1,0,0,0,0,0/)
       !  if(nch.eq.5) chan_avail=(/0,0,1,1,1,0,1,1,0,0,0,0,0,0,0/)
       !endif
       !write(log_lun,*)'  chan_avail= ',chan_avail(:)
       avail(:) = chan_avail(:)

       !allocate all necessary arrays
       !call alloc
! allocations (some unique to intam2)
      if(xtime.eq.1) then
       allocate (tai93time(nscans))
       allocate (sclat(nscans),sclon(nscans))
       allocate (scalt(nscans),stdtime(nscans,6))
       allocate (scorient(nscans))

       allocate (lat(npixs,nscans),lon(npixs,nscans))
       allocate (tbs(npixs,nscans,maxchans))
        tbs = miss_flt
       allocate (oetbs(npixs,nscans,maxchans),eia(npixs,nscans,maxchans))
       allocate (twb(npixs,nscans))
       !allocate (LR(npixs,nscans))
       allocate (ppskint(npixs,nscans),pptcwv(npixs,nscans))
       allocate (ppT2m(npixs,nscans),ppwind(npixs,nscans))
       allocate (L1Cqualflag(npixs,nscans))

       allocate (sataz(npixs,nscans),lo_flag(npixs,nscans,4)) ! new
       allocate (sst(npixs,nscans))

       allocate (sglinta(npixs,nscans),sfccode(npixs,nscans))
       allocate (CAPE(npixs,nscans)) ! new for 2017V1
       !allocate (orolifti(npixs,nscans),snowci(npixs,nscans))

       allocate (pixel_status(npixs,nscans))

       allocate (save_slp(npixs,nscans),save_wdir(npixs,nscans))
       allocate (save_tprof(npixs,nscans,nz),pp_wvprof(npixs,nscans,nz))
       allocate (pph(npixs,nscans,nz))

       ! was in prep_oe before:
      allocate(oe_output(npix,nscans,9),screen(npix,nscans)) !last dimension hard-coded!
      allocate(Tb_diff(npix,nscans,maxchans))
      allocate(save_sal(npix,nscans))
      allocate(save_iter(npix,nscans),poster(npix,nscans,8))
      allocate(toffsets(nch),s_toffsets(nbins,nch),s_sy(nch,nch,nbins))
      allocate(sy(nch,nch),sy_i(nch,nch))
      allocate(rs_sgm(nch,3),rc_sgm(nch,3))
      allocate(rc_syl(nch,nch),rc_sym(nch,nch),rc_syh(nch,nch)) !rain convective
      allocate(rs_syl(nch,nch),rs_sym(nch,nch),rs_syh(nch,nch)) !rain stratiform
      allocate(rs_r(nch,nch,3),rc_r(nch,nch,3))
      allocate(rsl_sgm(nch,3),rcl_sgm(nch,3))
      allocate(rcl_syl(nch,nch),rcl_sym(nch,nch),rcl_syh(nch,nch)) !rain convective
      allocate(rsl_syl(nch,nch),rsl_sym(nch,nch),rsl_syh(nch,nch)) !rain stratiform
      allocate(rsl_r(nch,nch,3),rcl_r(nch,nch,3))
      allocate(save_emis(nch),save_tau(nch,nz))
      allocate(Amat_out(npix,nscans,nvar),save_sw(npix,nscans))

        if(scend.gt.nscans) scend=nscans
        if(pxend.gt.npix)   pxend=npix

      endif

!--- loop over and read all scans

       igood = .false.
       do iscan = 1, nscans
       
         read(rlun, iostat=ios)  scnhdr !read in scan header
         if(ios .ne. 0) print*,'scnhdr read promblem'
        
         stdtime(iscan,1) = scnhdr%scandate%year
         stdtime(iscan,2) = scnhdr%scandate%month 
         stdtime(iscan,3) = scnhdr%scandate%day
         stdtime(iscan,4) = scnhdr%scandate%hour
         stdtime(iscan,5) = scnhdr%scandate%minute
         stdtime(iscan,6) = scnhdr%scandate%second

         first_good_date = 0
         if(first_good_date .eq. 0) then
            idate = stdtime(iscan,1)*10000 + stdtime(iscan,2)*100 + stdtime(iscan,3)
            if(idate .gt. 19870101 .and. idate.lt.20500101) then
               first_good_date = idate
            endif
         endif

         sclat(iscan)     = scnhdr%sclat
         sclon(iscan)     = scnhdr%sclon
         scalt(iscan)     = scnhdr%scalt    
         scorient(iscan)  = scnhdr%scorient !am2 is forward always
         tai93time(iscan) = scnhdr%TAI93

!--     read all pixel data in each scan

         do ipix = 1, npix
         
           read(rlun, iostat=ios)  px     !read in pixel header
           !if(ios .ne. 0) call GPM_reprt_err(14,ipix,blank)
          
                reccnt = reccnt + 1                     
                lat(ipix,iscan)       = px%latitude
                lon(ipix,iscan)       = px%longitude
                !Tbb(ipix,iscan,:)     = px%tbs
                if(satcode=='AM2')  tbs(ipix,iscan,1:13)=px%tbs(3:15) !due to 6GHz
                if(satcode.ne.'AM2')tbs(ipix,iscan,:)=px%tbs(:)
                oetbs(ipix,iscan,:)   = px%tbs
                eia(ipix,iscan,:)     = px%eia
                twb(ipix,iscan)       = px%twb
                pptcwv(ipix,iscan)    = px%tcwv
                ppskint(ipix,iscan)   = px%skint
                ppT2m(ipix,iscan)     = px%T2m
                sataz(ipix,iscan)     = px%azimuth !from AMSR L1R, not others
                ppwind(ipix,iscan)    = px%windmag
                save_slp(ipix,iscan)  = px%slp
                save_tprof(ipix,iscan,:) = 0.01*px%tprof(:)
                pp_wvprof(ipix,iscan,:)  = 0.001*px%wvprof(:)
                pph(ipix,iscan,:)     = px%aight(:)
                L1Cqualflag(ipix,iscan) = px%qualflag
                save_wdir(ipix,iscan) = px%winddir
                lo_flag(ipix,iscan,:) = px%lo_pct(:)
                sglinta(ipix,iscan)   = px%sunglint_angle
                sfccode(ipix,iscan)   = px%surface_type_index
                !orolifti(ipix,iscan)  = px%orolift_index
                !snowci(ipix,iscan)    = px%snow_cover_index

         enddo   !ipix
       enddo  !iscan
        
       !if(satcode.eq.'AM2'.or.satcode.eq.'AME') then
       !  if(convolve.eq.'res06') lo_flag(:,:) = land_fov(:,:,1) !6/7GHz FOV for AMSR
       !  if(convolve.eq.'res10') lo_flag(:,:) = land_fov(:,:,2) !10GHz FOV for AMSR
       !  if(convolve.eq.'res23') lo_flag(:,:) = land_fov(:,:,3) !19/23GHz FOV for AMSR
       !  if(convolve.eq.'res36') lo_flag(:,:) = land_fov(:,:,4) !36GHz FOV for AMSR
       !  if(convolve.ne.'res36'.and.convolve.ne.'res23'.and.convolve.ne.'res10'.and.convolve.ne.'res06') then
       !    write(6,*),'convolution incorrectly chosen in pp!'
       !    stop
       !  endif
       !else
       !  lo_flag(:,:) = land_fov(:,:,1) !should be same at all FOVs if from L1C
       !  convolve = 'None' ! currently no convolution done on GMI/TMI/etc.
       !endif
       !write(*,*)' number of sat pixels read  = ',reccnt
       !write(log_lun,*)' number of sat pixels read  = ',reccnt
       !if(reccnt .eq. 0) call GPM_reprt_err(15,reccnt,blank)
       if(reccnt .eq. 0) print*,'Oh, no!'
        
      return
      end subroutine read_geos5_pp
