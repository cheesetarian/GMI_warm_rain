      module integrated_output
        ! common netCDF output file for all algorithms, 
        !  D Duncan CSU, 9/29/16

      use define_csu1dvar
      use define_intam2
      use GPM_arraydef_D
      use netcdf

      implicit none

      contains

      subroutine output_int

      implicit none

      integer, parameter :: CMPRS = 1 !level of compression!
      integer :: ctime(8)
      integer :: i, ipix, iscan, itime
      integer :: lat_varid, lon_varid, tai93_varid, time_varid
      integer :: tpw_varid, wsp_varid, lwp_varid, sst_varid, chi_varid!,reysst_varid
      integer :: sun_varid, land_varid, qual_varid
      integer :: sic_varid, pop_varid, precip_varid, fprcp_varid, gcode_varid
      integer :: swe_varid, sdepth_varid, smoist_varid, stemp_varid
      integer :: ncid, pix_dims(2), time_dims(2)
      integer :: pix_dimid, time_dimid, scan_dimid
      integer :: pix_dim, time_dim, scan_dim
      character(len=40)  :: sdate!, sprior
      character(len=10)  :: sgran!,siter,smiss,sthresh,sthresh1, em_version

      integer*2, allocatable :: scantime(:,:), sunout(:,:)
      integer*1, allocatable :: qflag(:,:), landout(:,:,:)

      real, allocatable :: latout(:,:), lonout(:,:)
      real, allocatable :: tpwout(:,:), wspout(:,:), lwpout(:,:)
      real, allocatable :: sstout(:,:), chiout(:,:)!, reyout(:,:)

      real, allocatable :: sicout(:,:), sweout(:,:), sdepthout(:,:)
      real, allocatable :: precipout(:,:), fprcpout(:,:)
      real, allocatable :: stempout(:,:), smoistout(:,:)
      integer*1,allocatable :: popout(:,:),gcodeout(:,:)


      !---Allocate arrays---
      allocate (scantime(nscans,6),sunout(nscans,npix))
      allocate (qflag(nscans,npix),landout(nscans,npix,4))

      allocate (latout(nscans,npix),lonout(nscans,npix))
      allocate (tpwout(nscans,npix),wspout(nscans,npix))
      allocate (lwpout(nscans,npix))!,reyout(nscans,npix))
      allocate (sstout(nscans,npix),chiout(nscans,npix))

      allocate(sicout(nscans,npix),sweout(nscans,npix),sdepthout(nscans,npix))
      allocate(stempout(nscans,npix),smoistout(nscans,npix))
      allocate(precipout(nscans,npix),fprcpout(nscans,npix))
      allocate(gcodeout(nscans,npix),popout(nscans,npix))

      !---Copy data into output arrays---

      scantime(:,:) = stdtime(:,:)
      do iscan=1,nscans
        do ipix=1,npix

          latout(iscan,ipix)     = lat(ipix,iscan)
          lonout(iscan,ipix)     = lon(ipix,iscan)
          tpwout(iscan,ipix)     = ftrunc(oe_output(ipix,iscan,1), 100.0)
          lwpout(iscan,ipix)     = ftrunc(oe_output(ipix,iscan,2), 100.0)
          wspout(iscan,ipix)     = ftrunc(oe_output(ipix,iscan,3), 100.0)
          chiout(iscan,ipix)     = ftrunc(oe_output(ipix,iscan,5), 100.0)
          sstout(iscan,ipix)     = ftrunc(oe_output(ipix,iscan,6), 100.0)
          sunout(iscan,ipix)     = sglinta(ipix,iscan)
          landout(iscan,ipix,:)  = lo_flag(ipix,iscan,:)
          !reyout(iscan,ipix)     = ftrunc(sst(ipix,iscan),100.0)
          qflag(iscan,ipix)      = int(outarr(ipix,iscan,13)) ! output flag
          sicout(iscan,ipix)     = ftrunc(outarr(ipix,iscan,12), 100.0)
          !sweout(iscan,ipix)     = ftrunc(swe(ipix,iscan), 100.0)
          !sdepthout(iscan,ipix)  = ftrunc(sd(ipix,iscan), 100.0)
          smoistout(iscan,ipix)  = ftrunc(outarr(ipix,iscan,14), 100.0)
          stempout(iscan,ipix)   = ftrunc(outarr(ipix,iscan,15), 100.0)
          popout(iscan,ipix)     = int(pop(ipix,iscan)* 100.0)
          precipout(iscan,ipix)  = ftrunc(sfcprcp(ipix,iscan),1000.0)
          fprcpout(iscan,ipix)   = ftrunc(frzprcp(ipix,iscan), 100.0)
          gcodeout(iscan,ipix)   = sfccode(ipix,iscan) !modified, not original from pp!


! quality flag codes: 0 = great, 1 = converged but not great, 
        !  2 = no convergence, 3 = TPW quality check screened, 
        !  4 = sun glint, 5 = not run
        ![, 6 = land/ice contamination likely ]
!          qflag(iscan,ipix) = 5
!          if (oe_output(ipix,iscan,5).le.chis1 .and. oe_output(ipix,iscan,1).gt.0.0) qflag(iscan,ipix)=0
!          if (oe_output(ipix,iscan,5).gt.chis1 .and. oe_output(ipix,iscan,1).gt.0.0) qflag(iscan,ipix)=1
!          if (oe_output(ipix,iscan,1) .eq. miss_flt) qflag(iscan,ipix)=2
!          if (oe_output(ipix,iscan,1) .eq. -996) qflag(iscan,ipix)=3
!          !if (Tb_diff(ipix,iscan,1).gt.3.and.Tb_diff(ipix,iscan,2).gt.3 &
!          !   oe_output(ipix,iscan,3).lt.15) qflag(iscan,ipix)=6
!          if ((sglinta(ipix,iscan).ge.0 .and. sglinta(ipix,iscan).lt.20)) then
!             qflag(iscan,ipix)=4
!          endif
        enddo
      enddo 

      !---open output file---

      call check(nf90_create(out_file, NF90_NETCDF4, ncid))

      !---specify diminsions---

      pix_dim   = npix
      scan_dim  = nscans
      time_dim  = 6
      call check(nf90_def_dim(ncid, "npix",    pix_dim,  pix_dimid))
      call check(nf90_def_dim(ncid, "nscan",  scan_dim, scan_dimid))
      call check(nf90_def_dim(ncid, "ntime",  time_dim, time_dimid))
      time_dims = (/ scan_dimid, time_dimid /)
      pix_dims  = (/ scan_dimid, pix_dimid /)


      !---write global metadata---

      write(sgran,'(i6.6)') granule
      !write(siter,'(f5.2)') avg_iter
      !write(smiss,'(f7.1)') miss_flt
      !write(sthresh,'(f4.1)') chisq_out_thresh
      !write(sthresh1,'(f4.1)') chis1
      call date_and_time(values=ctime)

      !write(6,'(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2,"UTC")') &
      !      ctime(1),ctime(2),ctime(3),ctime(5),ctime(6),ctime(7)
      write(sdate,'(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2)') &
            ctime(1),ctime(2),ctime(3),ctime(5),ctime(6),ctime(7)
      !if(merra) sprior=trim("MERRA reanalysis")
      !if(.not.merra.and.useanalysis) sprior=trim("ECMWF Interim")
      !if(.not.useanalysis) sprior=trim("Climatology")
      !em_version = "FASTEM6"
      !if(run_rss) em_version = "RSS_2012"

      call check(nf90_put_att(ncid, NF90_GLOBAL, "Satellite",     trim(sat_name)))!ellite)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Sensor",           trim(sensor_name)))!atcode)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Granule",               trim(sgran)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Version",        trim(code_version)))
      !call check(nf90_put_att(ncid, NF90_GLOBAL, "RT_Model",        trim(rtm_version)))
      !call check(nf90_put_att(ncid, NF90_GLOBAL, "Emis_Model",       trim(em_version)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Preprocessor_Version",trim(pp_version)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "L1_File",           trim(radfile)))!orig_file)))
      !call check(nf90_put_att(ncid, NF90_GLOBAL, "Convolution_Resolution",trim('undefined'))) !trim(convolve)))
      !call check(nf90_put_att(ncid, NF90_GLOBAL, "Calibration_File",   trim(cal_file)))
      !call check(nf90_put_att(ncid, NF90_GLOBAL, "SY_File",             trim(sy_file)))
      !call check(nf90_put_att(ncid, NF90_GLOBAL, "Missing_Value",         trim(smiss)))
      !call check(nf90_put_att(ncid, NF90_GLOBAL, "ChiSquare_Threshold", trim(sthresh)))
      !call check(nf90_put_att(ncid, NF90_GLOBAL, "ChiSquare_Threshold_High_Quality",trim(sthresh1)))
      !call check(nf90_put_att(ncid, NF90_GLOBAL, "Source_of_Prior",      trim(sprior)))
      !call check(nf90_put_att(ncid, NF90_GLOBAL, "Average_Iterations",    trim(siter)))
      call check(nf90_put_att(ncid, NF90_GLOBAL, "Creation_Date",         trim(sdate)))

      !---declare variables---
! -- modified for ver2.3 to CF-compliant, D Duncan, 9/28/16
      call check(nf90_def_var(ncid, "sc_time",           NF90_SHORT, time_dims,     time_varid))
      call check(nf90_def_var(ncid, "time",              NF90_DOUBLE, scan_dimid,    tai93_varid))
      call check(nf90_def_var(ncid, "latitude",          NF90_FLOAT, pix_dims,      lat_varid))
      call check(nf90_def_var(ncid, "longitude",         NF90_FLOAT, pix_dims,      lon_varid))
      call check(nf90_def_var(ncid, "sunglint_angle",    NF90_SHORT, pix_dims,      sun_varid))
      call check(nf90_def_var(ncid, "land_area_fraction",NF90_BYTE, pix_dims,     land_varid))
      !ocean vars:
      call check(nf90_def_var(ncid,"atmosphere_mass_content_of_water_vapor", NF90_FLOAT, pix_dims,tpw_varid))
      call check(nf90_def_var(ncid, "wind_speed",        NF90_FLOAT, pix_dims,      wsp_varid))
      call check(nf90_def_var(ncid,"atmosphere_mass_content_of_cloud_liquid_water", NF90_FLOAT, pix_dims, lwp_varid))
      call check(nf90_def_var(ncid, "chi_squared",       NF90_FLOAT, pix_dims,      chi_varid))
      call check(nf90_def_var(ncid, "sea_surface_subskin_temperature", NF90_FLOAT,    pix_dims,      sst_varid))
      call check(nf90_def_var(ncid, "quality_flag",      NF90_BYTE, pix_dims,     qual_varid))
      !call check(nf90_def_var(ncid, "reynolds_sst",      NF90_FLOAT, pix_dims,   reysst_varid))
      !precip vars:
      call check(nf90_def_var(ncid, "surface_precip_rate", NF90_FLOAT, pix_dims, precip_varid))
      call check(nf90_def_var(ncid, "surface_frozen_precip_rate", NF90_FLOAT, pix_dims, fprcp_varid))
      call check(nf90_def_var(ncid, "precip_probability", NF90_FLOAT, pix_dims, pop_varid))
      call check(nf90_def_var(ncid, "gpm_sfc_code", NF90_BYTE, pix_dims, gcode_varid))
      !sic/swe/sm vars:
      call check(nf90_def_var(ncid, "sea_ice_concentration", NF90_FLOAT, pix_dims, sic_varid))
      call check(nf90_def_var(ncid, "snow_depth", NF90_FLOAT, pix_dims, sdepth_varid))
      call check(nf90_def_var(ncid, "snow_water_equivalent", NF90_FLOAT, pix_dims, swe_varid))
      call check(nf90_def_var(ncid, "soil_temperature", NF90_FLOAT, pix_dims, stemp_varid))
      call check(nf90_def_var(ncid, "soil_moisture", NF90_FLOAT, pix_dims, smoist_varid))

      !---compress variables---

      if (CMPRS .gt. 0) then
        call check(nf90_def_var_deflate(ncid,     time_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    tai93_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      lat_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      lon_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      tpw_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      wsp_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      lwp_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      chi_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      sst_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      sun_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,     land_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,     qual_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,   precip_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    fprcp_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      pop_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    gcode_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      sic_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,   sdepth_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,      swe_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,    stemp_varid, 0, 1, CMPRS))
        call check(nf90_def_var_deflate(ncid,   smoist_varid, 0, 1, CMPRS))
      endif
      call check(nf90_enddef(ncid))

      !---specify units---

      call check(nf90_put_att(ncid,   tai93_varid, "units",       "s"))
      call check(nf90_put_att(ncid,     lat_varid, "units","degree_north"))
      call check(nf90_put_att(ncid,     lon_varid, "units","degree_east"))
      call check(nf90_put_att(ncid,     tpw_varid, "units",  "kg m-2"))
      call check(nf90_put_att(ncid,     wsp_varid, "units",   "m s-1"))
      call check(nf90_put_att(ncid,     lwp_varid, "units",  "kg m-2"))
      call check(nf90_put_att(ncid,     sst_varid, "units",       "K"))
      call check(nf90_put_att(ncid,     sun_varid, "units",  "degree"))
      call check(nf90_put_att(ncid,    land_varid, "units",       "1"))
      call check(nf90_put_att(ncid,  precip_varid, "units", "mm hr-1"))
      call check(nf90_put_att(ncid,   fprcp_varid, "units", "mm hr-1"))
      call check(nf90_put_att(ncid,     pop_varid, "units", "percent"))
      call check(nf90_put_att(ncid,     sic_varid, "units", "percent"))
      call check(nf90_put_att(ncid,     swe_varid, "units",  "kg m-2"))
      call check(nf90_put_att(ncid,  sdepth_varid, "units",      "cm"))
      call check(nf90_put_att(ncid,   stemp_varid, "units",       "K"))
      call check(nf90_put_att(ncid,  smoist_varid, "units",  "kg m-2"))

      !---specify long name---

      call check(nf90_put_att(ncid,     time_varid, "long_name", "Scan time (year, month, day, hour, minute, second)"))
      call check(nf90_put_att(ncid,    tai93_varid, "long_name", "Seconds since 1993"))
      call check(nf90_put_att(ncid,      lat_varid, "long_name", "Latitude"))
      call check(nf90_put_att(ncid,      lon_varid, "long_name", "Longitude"))
      call check(nf90_put_att(ncid,      tpw_varid, "long_name", "Total precipitable water"))
      call check(nf90_put_att(ncid,      wsp_varid, "long_name", "Wind speed"))
      call check(nf90_put_att(ncid,      lwp_varid, "long_name", "Liquid water path"))
      call check(nf90_put_att(ncid,      chi_varid, "long_name", "Chi-squared"))
      call check(nf90_put_att(ncid,      sst_varid, "long_name", "Sea surface temperature"))
      call check(nf90_put_att(ncid,      sun_varid, "long_name", "Sun glint angle"))
      call check(nf90_put_att(ncid,     land_varid, "long_name", "Percentage of land coverage"))
      call check(nf90_put_att(ncid,     qual_varid, "long_name", "Quality flag"))
      call check(nf90_put_att(ncid,   precip_varid, "long_name", "Surface precipitation rate"))
      call check(nf90_put_att(ncid,    fprcp_varid, "long_name", "Frozen precipitation rate"))
      call check(nf90_put_att(ncid,      pop_varid, "long_name", "Probability of precipitation"))
      call check(nf90_put_att(ncid,    gcode_varid, "long_name", "GPM surface type"))
      call check(nf90_put_att(ncid,      sic_varid, "long_name", "Sea ice concentration"))
      call check(nf90_put_att(ncid,      swe_varid, "long_name", "Snow water equivalent"))
      call check(nf90_put_att(ncid,   sdepth_varid, "long_name", "Snow depth"))
      call check(nf90_put_att(ncid,    stemp_varid, "long_name", "Soil temperature"))
      call check(nf90_put_att(ncid,   smoist_varid, "long_name", "Soil moisture"))


      !---write variables---

      call check(nf90_put_var(ncid,     time_varid,     scantime))
      call check(nf90_put_var(ncid,    tai93_varid,    tai93time))
      call check(nf90_put_var(ncid,      lat_varid,       latout))
      call check(nf90_put_var(ncid,      lon_varid,       lonout))
      call check(nf90_put_var(ncid,      tpw_varid,       tpwout))
      call check(nf90_put_var(ncid,      wsp_varid,       wspout))
      call check(nf90_put_var(ncid,      lwp_varid,       lwpout))
      call check(nf90_put_var(ncid,      chi_varid,       chiout))
      call check(nf90_put_var(ncid,      sst_varid,       sstout))
      call check(nf90_put_var(ncid,      sun_varid,       sunout))
      call check(nf90_put_var(ncid,     land_varid,      landout))
      call check(nf90_put_var(ncid,     qual_varid,        qflag))
      call check(nf90_put_var(ncid,   precip_varid,    precipout))
      call check(nf90_put_var(ncid,    fprcp_varid,     fprcpout))
      call check(nf90_put_var(ncid,      pop_varid,       popout))
      call check(nf90_put_var(ncid,    gcode_varid,     gcodeout))
      call check(nf90_put_var(ncid,      sic_varid,       sicout))
      call check(nf90_put_var(ncid,      swe_varid,       sweout))
      call check(nf90_put_var(ncid,   sdepth_varid,    sdepthout))
      call check(nf90_put_var(ncid,    stemp_varid,     stempout))
      call check(nf90_put_var(ncid,   smoist_varid,    smoistout))

      !---close file---

      call check(nf90_close(ncid))

      !---Deallocate arrays---

        ! currently lots of vars not deallocated!
      deallocate (scantime)
      deallocate (qflag)
      deallocate (landout)
      deallocate (sunout)
      deallocate (latout)
      deallocate (lonout)
      deallocate (tpwout)
      deallocate (wspout)
      deallocate (lwpout)
      !deallocate (reyout)
      deallocate (sstout)
      deallocate (chiout)

      return
      end subroutine output_int

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

      end module integrated_output
