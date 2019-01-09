      program geos5_pp
!
!-------------------------------------------------------------------
!
!   AMSR2 pre processor for gprof2014, csu1dvar, and integrated
!       algorithms. 
!       adapted from work by D. Randel, CSU, 2013. 
!       modified to current form, D Duncan, CSU, Dec 2016.
!
!-------------------------------------------------------------------
     
      USE pp_definition
      USE pp_commoncode
      USE read_l1r
      USE read_l1c
      USE read_geos5
      USE gpm
      use netcdf
      
      implicit none
 
      integer :: n_arguments,iargc,ios,idate,i , ncid, stat

!--- Get command line inputs
    
      n_arguments = iargc() 
      write(6,*)' n_arguments = ', n_arguments
      if(n_arguments.ne.4 .and. n_arguments.ne.5 ) then
          write(*,*)' error in command line arguments '
          write(*,*)' <sens> <input file> <output_file> <res_opt>'&
                    '<cal_file (optional)>'
          write(*,*)'res_opt choices: res06,res10,res23,res36'
          stop 2
      endif
      call getarg(1, sens)           !sensor code, e.g. AM2,GMI,etc.
      call getarg(2, input_file)      !input filename
      call getarg(3, output_file)     !output filename
      call getarg(4, res_opt)         !res option 1-4, default is 1 (res06GHz)
      if(n_arguments .eq. 5) then
          call getarg(5, cal_file)        !calibration filename in ppancillary
      else
          cal_file = 'no calibration table used'
      endif

!--- echo startup information

      write(*,'(a,a)')' sensor code          : ', trim(sens)
      write(*,'(a,a)')' input  file name     : ', trim(input_file)
      write(*,'(a,a)')' output file name     : ', trim(output_file)
      write(*,'(a,a)')' resolution option    : ', trim(res_opt)
      write(*,'(a,a)')' calibration table    : ', trim(cal_file)

!--- Set database names for each product

      pp_version = '1701geos5'
      !prfdbfile  = 'GMI_ECMWF_V2.2_NXobs'    !CLIMATOLOGY database -- for now

      if(sens .eq. 'AM2') then
        sat_name = 'GCOM-W1'
        sensor_name = 'AMSR2'
        nchans = 11
        chan_freq = (/6.925,-9999,10.65,10.65,18.7,18.7,23.8,23.8,36.5,36.5,89.0,89.0,-9999,-9999,-9999/)
        prfdbfile = 'AMSR2_GANAL_V3' ! use OLD ganal dtb for now
      endif
      if(sens .eq. 'GMI') then
        sat_name = 'GPM'
        sensor_name = 'GMI'
        nchans = 13
        chan_freq = (/10.65,10.65,18.7,18.7,23.8,-9999,36.64,36.64,89.0,89.0,166.0,166.0,-9999,183.31-3,183.31-7/)
        prfdbfile = 'GMI_GANAL_V3' ! use OLD ganal dtb for now
      endif
      write(*,'(a,a)')' pp_version = ',trim(pp_version)
      write(*,'(a,a)')' prfdbfile  = ',trim(prfdbfile)

!--- set chan avail flags here, from chan_freq in pp_definition
      chan_avail = 1
      do i = 1,maxchans
        if(chan_freq(i) .eq. -9999.0) chan_avail(i) = 0
      enddo
      write(*,'(a,15I2)')'  chan_avail = ', chan_avail  

!--- read Tb data
      write(*,*)' starting read_l1 data'
      if(sens .eq. 'AM2') then
        write(*,*)' starting read l1r '
        call read_L1R_ppam2
        write(*,*)' fill tbs'
        call fill_low2high_qctbs
        !--- Use sun azimuth and elevation, and EIA, to calculate
        !       the sunglint angle (not given in L1R directly) 
        write(*,*)' calculating sunglint'
        call calc_sunglint
        write(*,*)' done sunglint'
      endif
      if(sens .eq. 'GMI') then
        call read_L1C_GMI
        call read_landmask    ! no land info in L1C
        call assign_sfctype   ! use landmask, assign to lo_flag
      endif

!--- call calibration routine
      if(n_arguments .eq. 5) then !careful here...
        write(*,*)' starting GPM_calibration'
        call GPM_calibration
      endif
       
!--- add ancillary datasets
      write(*,*)' reading GEOS5 FP IT data'
      call geos5 !( trim(gdir)//trim(fi1), trim(gdir)//trim(fi2))

      savegran = granule ! for some reason 'granule' gets corrupted during add_ancillary?
!--- get GPM surface codes
      write(*,*)' getting GPM sfc codes'
      call add_ancillary
      
!--- output all data to pp file
      write(*,*)' starting output_pp'
      call output_pp
         
!--- deallocate all variables
      call dealloc
!--- all done!
      write(*,*) 'successful completion!'

      end
