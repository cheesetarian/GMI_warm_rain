      module output_1dvar_nc

!------ contains subroutines for reading in pre-processor files,
!------  reading and processing ancillary data, and outputting

      use define_csu1dvar
      use define_intam2
      use GPM_arraydef_D
      use netcdf

      implicit none

      contains

      subroutine output_nc
      !subroutine output_nc(outdebug,outrt)

      implicit none
      logical :: outrt = .false. , outdebug = .false.
      !logical, intent ( in) :: outdebug, outrt

      integer, parameter :: CMPRS = 1 !level of compression!
      integer :: ctime(8)
      integer :: i, ipix, iscan, ich, ichan,cchan, itime, ieia
      integer :: time_varid, sclat_varid, sclon_varid, scalt_varid, scorient_varid
      integer :: lat_varid, lon_varid, press_varid, eia_varid, tbsim_varid, tbobs_varid
      integer :: tbdif_varid, sal_varid, slp_varid, tai93_varid
      integer :: tpw_varid, wsp_varid, lwp_varid, wdir_varid, sst_varid, reysst_varid
      integer :: chi_varid, chan_varid, rr_varid, iwp_varid, rwp_varid!eof_varid
      integer :: sun_varid, land_varid, qual_varid, niter_varid,sw_varid,amat_varid
      integer :: wprof_varid, tprof_varid, lprof_varid, iprof_varid, post_varid
      integer :: ncid, pix_dims(2), time_dims(2), prof_dims(3)!, eof_dims(3)
      integer :: eia_dims(3), tb_dims(3), chan_dims(2), post_dims(3)
      integer :: pix_dimid, time_dimid, scan_dimid, str_dimid, chan_dimid, eia_dimid
      integer :: npc_dimid, var_dimid, post_dimid, amat_dims(3)
      integer :: nlev_dimid, nz_dimid
      integer :: pix_dim, time_dim, str_dim, scan_dim, chan_dim, eia_dim
      integer :: npc_dim, nz_dim, nlev_dim, post_dim, var_dim
      character(len=40)  :: sdate, sprior
      character(len=40)  :: sfreq
      character(len=10)  :: sgran,siter,smiss,sthresh,sthresh1, em_version
      character(len=6)   :: schan(maxchans)
      integer :: neia=1

      integer,   allocatable :: scantime(:,:)
      integer*1, allocatable :: scorientout(:)
      integer*1, allocatable :: qflag(:,:)
      integer*1, allocatable :: landout(:,:)
      integer*2, allocatable :: sunout(:,:)
      integer*1, allocatable :: iterout(:,:)
      integer*1, allocatable :: swout(:,:)

      real, allocatable :: latout(:,:)
      real, allocatable :: lonout(:,:)
      real, allocatable :: tpwout(:,:)
      real, allocatable :: wspout(:,:)
      real, allocatable :: lwpout(:,:)
      real, allocatable :: iwpout(:,:)
      real, allocatable :: rwpout(:,:)
      real, allocatable :: sstout(:,:)
      real, allocatable :: chiout(:,:)

      real, allocatable :: reyout(:,:)
      real, allocatable :: salout(:,:)
      real, allocatable :: slpout(:,:)
      integer*2, allocatable :: dirout(:,:)
      real, allocatable :: post(:,:,:)
      !real, allocatable :: eia(:,:,:)
      real, allocatable :: tbsim(:,:,:)
      real, allocatable :: tbobs(:,:,:)
      real, allocatable :: tbdif(:,:,:)
      real, allocatable :: amatout(:,:,:)
      real, allocatable :: rrout(:,:)
!      real, allocatable :: eofout(:,:,:)
!      real, allocatable :: wvp_prof(:,:,:)
!      real, allocatable :: tmp_prof(:,:,:)
!      real, allocatable :: lwp_prof(:,:,:)
      !real, allocatable :: iwp_prof(:,:,:)

      !---Allocate arrays---
      allocate (scantime(nscans,6))
      allocate (qflag(nscans,npix))
      allocate (landout(nscans,npix))
      allocate (iterout(nscans,npix))
      allocate (swout(nscans,npix))
      allocate (sunout(nscans,npix))
      allocate (scorientout(nscans))

      allocate (latout(nscans,npix))
      allocate (lonout(nscans,npix))
      allocate (tpwout(nscans,npix))
      allocate (wspout(nscans,npix))
      allocate (dirout(nscans,npix))
      allocate (lwpout(nscans,npix))
      allocate (iwpout(nscans,npix))
      allocate (rwpout(nscans,npix))
      allocate (reyout(nscans,npix))
      allocate (sstout(nscans,npix))
      allocate (chiout(nscans,npix))
      allocate (rrout(nscans,npix))

      allocate (salout(nscans,npix))
      allocate (slpout(nscans,npix))
      !allocate (eia(nscans,npix,neia))
      allocate (post(nscans,npix,4))
!      allocate (tbsim(nscans,npix,nch))
!      allocate (tbobs(nscans,npix,nch))
      allocate (tbdif(nscans,npix,nch))
      !allocate (eofout(nscans,npix,npc))
!      allocate (wvp_prof(nscans,npix,nz))
!      allocate (tmp_prof(nscans,npix,nz))
!      allocate (lwp_prof(nscans,npix,nz))
      !allocate (iwp_prof(nscans,npix,nz))
      allocate (amatout(nscans,npix,nvar))

      !---Copy data into output arrays---

      scantime(:,:) = stdtime(:,:)
      do iscan=1,nscans
        !do itime=1,6
        !  scantime(iscan,itime) = stdtime(iscan,itime)
        !enddo
        scorientout(iscan)       = 0.0 !int(scorient(iscan))
        do ipix=1,npix

          latout(iscan,ipix)     = lat(ipix,iscan)
          lonout(iscan,ipix)     = lon(ipix,iscan)
          tpwout(iscan,ipix)     = ftrunc(oe_output(ipix,iscan,1), 100.0)
          lwpout(iscan,ipix)     = ftrunc(oe_output(ipix,iscan,2), 100.0)
          wspout(iscan,ipix)     = ftrunc(oe_output(ipix,iscan,3), 100.0)
          rrout(iscan,ipix)      = ftrunc(oe_output(ipix,iscan,4), 100.0)
          chiout(iscan,ipix)     = ftrunc(oe_output(ipix,iscan,5), 100.0)
          sstout(iscan,ipix)     = ftrunc(oe_output(ipix,iscan,6), 100.0)
          iwpout(iscan,ipix)     = ftrunc(oe_output(ipix,iscan,8), 100.0)
          rwpout(iscan,ipix)     = ftrunc(oe_output(ipix,iscan,9), 100.0)
!          dirout(iscan,ipix)     = int(save_wdir(ipix,iscan))
          !do i=1,npc
          !  eofout(iscan,ipix,i) = ftrunc(oe_output(ipix,iscan,6+npc), 100.0)
          !enddo
          post(iscan,ipix,1)     = ftrunc(poster(ipix,iscan,8), 100.0)
          post(iscan,ipix,2)     = ftrunc(poster(ipix,iscan,4), 100.0)
          post(iscan,ipix,3)     = ftrunc(poster(ipix,iscan,5), 100.0)
          post(iscan,ipix,4)     = ftrunc(poster(ipix,iscan,7), 100.0)
          !do i=1,nvar+1
          !  post(iscan,ipix,i)   = ftrunc(poster(ipix,iscan,i), 100.0)
          !enddo
!          do i=1,nz
!            wvp_prof(iscan,ipix,i) = ftrunc(save_mrprof(ipix,iscan,i), 100.0)
!            tmp_prof(iscan,ipix,i) = ftrunc(save_tprof(ipix,iscan,i), 100.0)
!            lwp_prof(iscan,ipix,i) = ftrunc(save_clwc(ipix,iscan,i), 100.0)
!            !iwp_prof(iscan,ipix,i) = ftrunc(save_ciwc(ipix,iscan,i), 100.0)
!          enddo

!          salout(iscan,ipix)       = ftrunc(save_sal(ipix,iscan),100.0)
!          slpout(iscan,ipix)       = ftrunc(save_slp(ipix,iscan),100.0)
          sunout(iscan,ipix)       = sglinta(ipix,iscan)
          landout(iscan,ipix)      = lo_flag(ipix,iscan,1)
          iterout(iscan,ipix)      = save_iter(ipix,iscan)
          reyout(iscan,ipix)       = ftrunc(sst(ipix,iscan),100.0)
          swout(iscan,ipix)        = save_sw(ipix,iscan)

! quality flag codes: 0 = great, 1 = converged but not great, 
        !  2 = no convergence, 3 = TPW quality check screened, 
        !  4 = sun glint, 5 = not run
        ![, 6 = land/ice contamination likely ]
          qflag(iscan,ipix) = 5
          if (oe_output(ipix,iscan,5).le.chis1 .and. oe_output(ipix,iscan,1).gt.0.0) qflag(iscan,ipix)=0
          if (oe_output(ipix,iscan,5).gt.chis1 .and. oe_output(ipix,iscan,1).gt.0.0) qflag(iscan,ipix)=1
          if (oe_output(ipix,iscan,1) .eq. miss_flt) qflag(iscan,ipix)=2
          if (oe_output(ipix,iscan,1) .eq. -996) qflag(iscan,ipix)=3
          !if (Tb_diff(ipix,iscan,1).gt.3.and.Tb_diff(ipix,iscan,2).gt.3 &
          !   oe_output(ipix,iscan,3).lt.15) qflag(iscan,ipix)=6
          if ((sglinta(ipix,iscan).ge.0 .and. sglinta(ipix,iscan).lt.20)) then
             qflag(iscan,ipix)=4
             !tpwout(iscan,ipix) = miss_flt 
             !lwpout(iscan,ipix) = miss_flt 
             !wspout(iscan,ipix) = miss_flt 
             !!iwpout(iscan,ipix) = miss_flt 
             !chiout(iscan,ipix) = miss_flt 
             !sstout(iscan,ipix) = miss_flt 
          endif
          do ieia=1,nvar
            amatout(iscan,ipix,ieia) = ftrunc(Amat_out(ipix,iscan,ieia),100.0)
          enddo
!          do ieia=1,neia
!            eia(iscan,ipix,ieia) = ftrunc(eia_out(ipix,iscan,ieia), 100.0)
!          enddo
!          cchan=1
          do ichan=1,nch !maxchans !nch
            !if (chan_avail(ichan) .eq. 1) then
            !if (avail(ichan) .eq. 1) then
!              tbobs(iscan,ipix,cchan) = ftrunc(Tbb(ipix,iscan,ichan), 100.0)
!              tbsim(iscan,ipix,cchan) = ftrunc(Tb_sim(ipix,iscan,cchan), 100.0)
              tbdif(iscan,ipix,ichan) = ftrunc(Tb_diff(ipix,iscan,ichan), 100.0)
              !tbdif(iscan,ipix,cchan) = ftrunc(Tb_diff(ipix,iscan,cchan), 100.0)
              !cchan=cchan+1
            !else
            !  tbobs(iscan,ipix,ichan) = miss_flt
            !  tbsim(iscan,ipix,ichan) = miss_flt
            !  tbdif(iscan,ipix,ichan) = miss_flt
            !endif
          enddo
        enddo
      enddo 

      ich = 1
      do ichan=1,maxchans
        !if (chan_avail(ichan) .eq. 1) then
        if (avail(ichan) .eq. 1) then
          schan(ich) = '' !chan_freqs(ichan)
          ich = ich + 1
        endif
      enddo

      !---open output file---

      call check(nf90_create(csu1dvar_out, NF90_NETCDF4, ncid))

      !---specify diminsions---

      pix_dim   = npix
      scan_dim  = nscans
      chan_dim  = nch
      time_dim  = 6
      str_dim   = 7
      eia_dim   = nEIA
      var_dim   = nvar
      post_dim  = 4 !nvar+1
      !npc_dim   = npc
      nz_dim    = nz
      nlev_dim  = nz+1
      call check(nf90_def_dim(ncid, "npix",    pix_dim,  pix_dimid))
      call check(nf90_def_dim(ncid, "nscan",  scan_dim, scan_dimid))
      call check(nf90_def_dim(ncid, "nchan",  chan_dim, chan_dimid))
      call check(nf90_def_dim(ncid, "ntime",  time_dim, time_dimid))
      call check(nf90_def_dim(ncid, "strlen",  str_dim,  str_dimid))
      call check(nf90_def_dim(ncid, "neia",    eia_dim,  eia_dimid))
      call check(nf90_def_dim(ncid, "nvar",    var_dim,  var_dimid))
      call check(nf90_def_dim(ncid, "npost",  post_dim, post_dimid))
      !call check(nf90_def_dim(ncid, "neof",    npc_dim,  npc_dimid))
      call check(nf90_def_dim(ncid, "nlayer",   nz_dim,   nz_dimid))
      call check(nf90_def_dim(ncid, "nlevel", nlev_dim, nlev_dimid))
      time_dims = (/ scan_dimid, time_dimid /)
      pix_dims  = (/ scan_dimid, pix_dimid /)
      !eof_dims  = (/ scan_dimid, pix_dimid, npc_dimid /)
      prof_dims = (/ scan_dimid, pix_dimid, nz_dimid /)
      post_dims = (/ scan_dimid, pix_dimid, post_dimid /)
      amat_dims = (/ scan_dimid, pix_dimid, var_dimid /)
      eia_dims  = (/ scan_dimid, pix_dimid, eia_dimid /)
      tb_dims   = (/ scan_dimid, pix_dimid, chan_dimid /)
      chan_dims = (/ str_dimid,  chan_dimid /)

      !write(6,*) 'nch       = ',nch
      !write(6,*) 'chan_dim  = ',chan_dim
      !write(6,*) 'tb_dims   = ',tb_dims
      !write(6,*) 'chan_dims = ',chan_dims

      !---write global metadata---

      write(sgran,'(i6.6)') granule
      write(siter,'(f5.2)') avg_iter
      write(smiss,'(f7.1)') miss_flt
      write(sthresh,'(f6.1)') chout !chisq_out_thresh
      write(sthresh1,'(f4.1)') chis1
      call date_and_time(values=ctime)

      !write(6,'(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2,"UTC")') &
      !      ctime(1),ctime(2),ctime(3),ctime(5),ctime(6),ctime(7)
      write(sdate,'(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2)') &
            ctime(1),ctime(2),ctime(3),ctime(5),ctime(6),ctime(7)
      if(merra) sprior=trim("MERRA reanalysis")
      if(.not.merra.and.useanalysis) sprior=trim("ECMWF Interim")
      if(.not.useanalysis) sprior=trim("Climatology")
      em_version = "FASTEM6"
      !if(run_rss) em_version = "RSS_2012"

      call check(nf90_put_att(ncid, NF90_GLOBAL, "Satellite",     trim(sat_name)))!ellite)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Sensor",           trim(sensor_name)))!atcode)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Granule",               trim(sgran)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Version",        trim(code_version)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "RT_Model",        trim(rtm_version)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Emis_Model",       trim(em_version)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Preprocessor_Version",trim(pp_version)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "L1_File",           trim(radfile)))!orig_file)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Convolution_Resolution",trim('undefined'))) !trim(convolve)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Calibration_File",   trim(cal_file)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "SY_File",             trim(sy_file)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Missing_Value",         trim(smiss)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "ChiSquare_Threshold", trim(sthresh)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "ChiSquare_Threshold_High_Quality",trim(sthresh1)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Source_of_Prior",      trim(sprior)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Average_Iterations",    trim(siter)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Creation_Date",         trim(sdate)))

      !---declare variables---
! -- modified for ver2.3 to CF-compliant, D Duncan, 9/28/16
      call check(nf90_def_var(ncid, "sc_time",           NF90_SHORT, time_dims,     time_varid))
      call check(nf90_def_var(ncid, "sc_latitude",       NF90_FLOAT, scan_dimid,    sclat_varid))
      call check(nf90_def_var(ncid, "sc_longitude",      NF90_FLOAT, scan_dimid,    sclon_varid))
      call check(nf90_def_var(ncid, "sc_altitude",       NF90_FLOAT, scan_dimid,    scalt_varid))
      call check(nf90_def_var(ncid, "sc_orientation",    NF90_BYTE, scan_dimid, scorient_varid))
      call check(nf90_def_var(ncid, "time",              NF90_DOUBLE, scan_dimid,    tai93_varid))
      call check(nf90_def_var(ncid, "channel_name",      NF90_CHAR, chan_dims,     chan_varid))

      call check(nf90_def_var(ncid, "latitude",          NF90_FLOAT, pix_dims,      lat_varid))
      call check(nf90_def_var(ncid, "longitude",         NF90_FLOAT, pix_dims,      lon_varid))
      call check(nf90_def_var(ncid,"atmosphere_mass_content_of_water_vapor", NF90_FLOAT, pix_dims,tpw_varid))
      call check(nf90_def_var(ncid, "wind_speed",        NF90_FLOAT, pix_dims,      wsp_varid))
      call check(nf90_def_var(ncid,"atmosphere_mass_content_of_cloud_liquid_water", NF90_FLOAT, pix_dims, lwp_varid))
      call check(nf90_def_var(ncid, "iwp",            NF90_FLOAT, pix_dims,      iwp_varid))
      call check(nf90_def_var(ncid, "rwp",            NF90_FLOAT, pix_dims,      rwp_varid))
      call check(nf90_def_var(ncid, "chi_squared",       NF90_FLOAT, pix_dims,      chi_varid))
      call check(nf90_def_var(ncid, "rain_rate",         NF90_FLOAT, pix_dims,      rr_varid))
      call check(nf90_def_var(ncid, "sea_surface_subskin_temperature", NF90_FLOAT,    pix_dims,      sst_varid))
      call check(nf90_def_var(ncid, "iterations",        NF90_BYTE, pix_dims,    niter_varid))
      call check(nf90_def_var(ncid, "switchflag",        NF90_BYTE, pix_dims,       sw_varid))
      call check(nf90_def_var(ncid, "sunglint_angle",    NF90_FLOAT, pix_dims,      sun_varid))
      call check(nf90_def_var(ncid, "land_area_fraction",NF90_BYTE, pix_dims,     land_varid))
      call check(nf90_def_var(ncid, "quality_flag",      NF90_BYTE, pix_dims,     qual_varid))
      call check(nf90_def_var(ncid, "posterior_error",   NF90_FLOAT, post_dims,     post_varid))
      call check(nf90_def_var(ncid, "tb_difference",     NF90_FLOAT, tb_dims,    tbdif_varid))
      call check(nf90_def_var(ncid, "amatrix",           NF90_FLOAT, amat_dims,    amat_varid))
      call check(nf90_def_var(ncid, "reynolds_sst",      NF90_FLOAT, pix_dims,   reysst_varid))

!      if (outrt) then
!        call check(nf90_def_var(ncid, "eia",          NF90_FLOAT,    eia_dims,      eia_varid))
!        call check(nf90_def_var(ncid, "wdir",         NF90_SHORT,    pix_dims,     wdir_varid))
!        call check(nf90_def_var(ncid, "wvp_prof",     NF90_FLOAT,   prof_dims,    wprof_varid))
!        call check(nf90_def_var(ncid, "tmp_prof",     NF90_FLOAT,   prof_dims,    tprof_varid))
!        call check(nf90_def_var(ncid, "lwp_prof",     NF90_FLOAT,   prof_dims,    lprof_varid))
!        !call check(nf90_def_var(ncid, "iwp_prof",     NF90_FLOAT,   prof_dims,    iprof_varid))
!        call check(nf90_def_var(ncid, "tbobs",        NF90_FLOAT,     tb_dims,    tbobs_varid))
!        call check(nf90_def_var(ncid, "tbsim",        NF90_FLOAT,     tb_dims,    tbsim_varid))
!        !call check(nf90_def_var(ncid, "EOF_Coeffs",   NF90_FLOAT,    eof_dims,      eof_varid))
!        call check(nf90_def_var(ncid, "sssal",        NF90_FLOAT,    pix_dims,      sal_varid))
!        call check(nf90_def_var(ncid, "slp",          NF90_FLOAT,    pix_dims,      slp_varid))
!      endif

      !---compress variables---

      if (CMPRS .gt. 0) then
        call check(nf90_def_var_deflate(ncid,     time_varid, 0, 1, CMPRS))
        !call check(nf90_def_var_deflate(ncid,    press_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    sclat_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    sclon_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    scalt_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid, scorient_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    tai93_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      lat_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      lon_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      tpw_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      wsp_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      lwp_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      iwp_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      rwp_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      chi_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,       rr_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      sst_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    niter_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,       sw_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      sun_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,     land_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,     qual_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,     post_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    tbdif_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,     amat_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,   reysst_varid, 0, 1, CMPRS))
!        if (outrt) then
!          call check(nf90_def_var_deflate(ncid,    eia_varid, 0, 1, CMPRS))
!          call check(nf90_def_var_deflate(ncid,   wdir_varid, 0, 1, CMPRS))
!          call check(nf90_def_var_deflate(ncid,  wprof_varid, 0, 1, CMPRS))
!          call check(nf90_def_var_deflate(ncid,  tprof_varid, 0, 1, CMPRS))
!          call check(nf90_def_var_deflate(ncid,  lprof_varid, 0, 1, CMPRS))
!          !call check(nf90_def_var_deflate(ncid,  iprof_varid, 0, 1, CMPRS))
!          call check(nf90_def_var_deflate(ncid,  tbobs_varid, 0, 1, CMPRS))
!          call check(nf90_def_var_deflate(ncid,  tbsim_varid, 0, 1, CMPRS))
!          !call check(nf90_def_var_deflate(ncid,    eof_varid, 0, 1, CMPRS))
!          call check(nf90_def_var_deflate(ncid,    sal_varid, 0, 1, CMPRS))
!          call check(nf90_def_var_deflate(ncid,    slp_varid, 0, 1, CMPRS))
!        endif
      endif
      call check(nf90_enddef(ncid))

      !---specify units---

      !call check(nf90_put_att(ncid,   press_varid, "units",      "mb"))
      call check(nf90_put_att(ncid,   sclat_varid, "units", "degrees"))
      call check(nf90_put_att(ncid,   sclon_varid, "units", "degrees"))
      call check(nf90_put_att(ncid,   scalt_varid, "units",      "km"))
      call check(nf90_put_att(ncid,   tai93_varid, "units",     "sec"))
      call check(nf90_put_att(ncid,     lat_varid, "units", "degrees"))
      call check(nf90_put_att(ncid,     lon_varid, "units", "degrees"))
      call check(nf90_put_att(ncid,     tpw_varid, "units",      "mm"))
      call check(nf90_put_att(ncid,     wsp_varid, "units",     "m/s"))
      call check(nf90_put_att(ncid,     lwp_varid, "units",   "g/m^2"))
      call check(nf90_put_att(ncid,     iwp_varid, "units",   "g/m^2"))
      call check(nf90_put_att(ncid,     rwp_varid, "units",   "g/m^2"))
      call check(nf90_put_att(ncid,      rr_varid, "units",   "mm/hr"))
      call check(nf90_put_att(ncid,     sst_varid, "units",       "K"))
      call check(nf90_put_att(ncid,     sun_varid, "units", "degrees"))
      call check(nf90_put_att(ncid,    land_varid, "units", "percent"))
      call check(nf90_put_att(ncid,   tbdif_varid, "units",       "K"))
      !call check(nf90_put_att(ncid,    amat_varid, "units",       "-"))
      call check(nf90_put_att(ncid,  reysst_varid, "units",       "K"))
!      if (outrt) then
!        call check(nf90_put_att(ncid,   eia_varid, "units", "degrees"))
!        call check(nf90_put_att(ncid,  wdir_varid, "units", "degrees"))
!        call check(nf90_put_att(ncid, wprof_varid, "units",    "g/kg"))
!        call check(nf90_put_att(ncid, tprof_varid, "units",       "K"))
!        call check(nf90_put_att(ncid, lprof_varid, "units","kg/m^2/layer"))
!        !call check(nf90_put_att(ncid, iprof_varid, "units",    "g/kg"))
!        call check(nf90_put_att(ncid, tbobs_varid, "units",       "K"))
!        call check(nf90_put_att(ncid, tbsim_varid, "units",       "K"))
!        !call check(nf90_put_att(ncid,   eof_varid, "units",       ""))
!        call check(nf90_put_att(ncid, sal_varid, "units",       "psu"))
!        call check(nf90_put_att(ncid, slp_varid, "units",        "mb"))
!      endif

      !---specify long name---

      call check(nf90_put_att(ncid,     time_varid, "long_name", "Scan time (year, month, day, hour, minute, second)"))
      !call check(nf90_put_att(ncid,    press_varid, "long_name", "Level pressure"))
      call check(nf90_put_att(ncid,     chan_varid, "long_name", "Channel specification"))
      call check(nf90_put_att(ncid,    sclat_varid, "long_name", "Spacecraft latitude"))
      call check(nf90_put_att(ncid,    sclon_varid, "long_name", "Spacecraft longitude"))
      call check(nf90_put_att(ncid,    scalt_varid, "long_name", "Spacecraft altitude"))
      call check(nf90_put_att(ncid, scorient_varid, "long_name", "Spacecraft orientation (0=forward, 180=aft)"))
      call check(nf90_put_att(ncid,    tai93_varid, "long_name", "Seconds since 1993"))
      call check(nf90_put_att(ncid,      lat_varid, "long_name", "Latitude"))
      call check(nf90_put_att(ncid,      lon_varid, "long_name", "Longitude"))
      call check(nf90_put_att(ncid,      tpw_varid, "long_name", "Total precipitable water"))
      call check(nf90_put_att(ncid,      wsp_varid, "long_name", "Wind speed"))
      call check(nf90_put_att(ncid,      lwp_varid, "long_name", "Liquid water path"))
      call check(nf90_put_att(ncid,      iwp_varid, "long_name", "Ice water path"))
      call check(nf90_put_att(ncid,      rwp_varid, "long_name", "Rain water path"))
      call check(nf90_put_att(ncid,      chi_varid, "long_name", "Chi-squared"))
      call check(nf90_put_att(ncid,       rr_varid, "long_name", "Rain Rate"))
      call check(nf90_put_att(ncid,      sst_varid, "long_name", "Sea surface temperature"))
      call check(nf90_put_att(ncid,    niter_varid, "long_name", "Number of iterations"))
      call check(nf90_put_att(ncid,       sw_varid, "long_name", "NR Strat Conv"))
      call check(nf90_put_att(ncid,      sun_varid, "long_name", "Sun glint angle"))
      call check(nf90_put_att(ncid,     land_varid, "long_name", "Percentage of land coverage"))
      call check(nf90_put_att(ncid,     qual_varid, "long_name", "Quality flag"))
      call check(nf90_put_att(ncid,     post_varid, "long_name", "Posteriori errors (sigma)"))
      call check(nf90_put_att(ncid,   reysst_varid, "long_name", "Reynolds sea surface temperature"))

!      if (outrt) then
!        call check(nf90_put_att(ncid,    eia_varid, "long_name", "Earth incidence angle"))
!        call check(nf90_put_att(ncid,   wdir_varid, "long_name", "Wind direction"))
!        call check(nf90_put_att(ncid,  wprof_varid, "long_name", "Water vapor mixing ratio"))
!        call check(nf90_put_att(ncid,  tprof_varid, "long_name", "Mean layer temperature"))
!        call check(nf90_put_att(ncid,  lprof_varid, "long_name", "Cloud liquid water content"))
!        !call check(nf90_put_att(ncid,  wprof_varid, "long_name", "Cloud ice water content"))
!        call check(nf90_put_att(ncid,  tbobs_varid, "long_name", "Observed brightness temperatures"))
!        call check(nf90_put_att(ncid,  tbsim_varid, "long_name", "Simulated brightness temperatures"))
!        !call check(nf90_put_att(ncid,    eof_varid, "long_name", "Water vapor profile EOF coefficients"))
!        call check(nf90_put_att(ncid,    sal_varid, "long_name", "Sea Surface Salinity"))
!        call check(nf90_put_att(ncid,    slp_varid, "long_name", "Sea Level Pressure"))
!      endif

      !---write variables---

      call check(nf90_put_var(ncid,     time_varid,     scantime))
      !call check(nf90_put_var(ncid,    press_varid,       lpress))
      call check(nf90_put_var(ncid,     chan_varid, schan(1:nch)))
      call check(nf90_put_var(ncid,    sclat_varid,        sclat))
      call check(nf90_put_var(ncid,    sclon_varid,        sclon))
      call check(nf90_put_var(ncid,    scalt_varid,        scalt))
      call check(nf90_put_var(ncid, scorient_varid,  scorientout))
      call check(nf90_put_var(ncid,    tai93_varid,    tai93time))
      call check(nf90_put_var(ncid,      lat_varid,       latout))
      call check(nf90_put_var(ncid,      lon_varid,       lonout))
      call check(nf90_put_var(ncid,      tpw_varid,       tpwout))
      call check(nf90_put_var(ncid,      wsp_varid,       wspout))
      call check(nf90_put_var(ncid,      lwp_varid,       lwpout))
      call check(nf90_put_var(ncid,      iwp_varid,       iwpout))
      call check(nf90_put_var(ncid,      rwp_varid,       rwpout))
      call check(nf90_put_var(ncid,      chi_varid,       chiout))
      call check(nf90_put_var(ncid,       rr_varid,        rrout))
      call check(nf90_put_var(ncid,      sst_varid,       sstout))
      call check(nf90_put_var(ncid,    niter_varid,      iterout))
      call check(nf90_put_var(ncid,       sw_varid,        swout))
      call check(nf90_put_var(ncid,      sun_varid,       sunout))
      call check(nf90_put_var(ncid,     land_varid,      landout))
      call check(nf90_put_var(ncid,     qual_varid,        qflag))
      call check(nf90_put_var(ncid,     post_varid,         post))
      call check(nf90_put_var(ncid,    tbdif_varid,        tbdif))
      call check(nf90_put_var(ncid,     amat_varid,      amatout))
      call check(nf90_put_var(ncid,   reysst_varid,       reyout))

!      if (outrt) then
!        call check(nf90_put_var(ncid,    eia_varid,          eia))
!        call check(nf90_put_var(ncid,   wdir_varid,       dirout))
!        call check(nf90_put_var(ncid,  wprof_varid,     wvp_prof))
!        call check(nf90_put_var(ncid,  tprof_varid,     tmp_prof))
!        call check(nf90_put_var(ncid,  lprof_varid,     lwp_prof))
!        !call check(nf90_put_var(ncid,  iprof_varid,     iwp_prof))
!        call check(nf90_put_var(ncid,  tbobs_varid,        tbobs))
!        call check(nf90_put_var(ncid,  tbsim_varid,        tbsim))
!        !call check(nf90_put_var(ncid,    eof_varid,       eofout))
!        call check(nf90_put_var(ncid,    sal_varid,       salout))
!        call check(nf90_put_var(ncid,    slp_varid,       slpout))
!      endif

      !---close file---

      call check(nf90_close(ncid))

      !---Deallocate arrays---

        ! currently lots of vars not deallocated!
      deallocate (scantime)
      deallocate (qflag)
      deallocate (landout)
      deallocate (sunout)
      deallocate (scorientout)
      deallocate (latout)
      deallocate (lonout)
      deallocate (tpwout)
      deallocate (wspout)
      deallocate (dirout)
      deallocate (lwpout)
      deallocate (rrout)
      deallocate (iwpout)
      deallocate (rwpout)
      deallocate (reyout)
      deallocate (sstout)
      deallocate (chiout)
!      deallocate (salout)
!      deallocate (slpout)
!      deallocate (eia)
      deallocate (post)
!      deallocate (tbobs)
!      deallocate (tbsim)
      deallocate (tbdif)
      deallocate (amatout,swout)
      !deallocate (eofout)
!      deallocate (wvp_prof)
!      deallocate (tmp_prof)
!      deallocate (lwp_prof)
      !deallocate (iwp_prof)

      return
      end subroutine output_nc

      subroutine check(status)

      integer, intent ( in) :: status
    
      if (status /= nf90_noerr) then 
         print *, trim(nf90_strerror(status))
         stop "Stopped"
      end if

      end subroutine check  

      function ftrunc(rval, rscale)

      real, intent (in) :: rval
      real, intent (in) :: rscale
      real              :: ftrunc

      if (rval .lt. 0.0) then
        ftrunc = real( int( rval * rscale - 0.5)) / rscale
      else
        ftrunc = real( int( rval * rscale + 0.5)) / rscale
      endif
      
      end function ftrunc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module output_1dvar_nc
