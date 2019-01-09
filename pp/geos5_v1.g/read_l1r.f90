      module read_l1r
       use pp_definition
       use pp_commoncode

       contains

      subroutine read_L1R_ppam2
      USE HDF5 ! This module contains all necessary modules
      USE ISO_C_BINDING

      implicit none

!-- Declaration of HDF5 functions and variables

        INTEGER(HID_T) :: file_id       ! File identifier
        INTEGER(HID_T) :: dset_id       ! Dataset identifier
        INTEGER(HID_T) :: attr_id       ! Attribute identifier
        INTEGER(HID_T) :: space_id      ! Dataspace identifier
        INTEGER :: hdferr
        INTEGER(HSIZE_T), DIMENSION(2)   :: dmaxdims  ! max dimensions for data sets
        INTEGER :: len

        TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, TARGET :: rdata     !Read buffer
        INTEGER(hsize_t), DIMENSION(1:1) ::  dims = (/1/)
        INTEGER(hsize_t), DIMENSION(2)   :: ddims  ! data sets are 2-D arrays
        INTEGER(hsize_t), DIMENSION(3) :: lodims ! l/o flags are in 3D array
        TYPE(C_PTR) :: f_ptr
        CHARACTER(len=256, kind=c_char), POINTER :: data ! A pointer to Ca Fortran string


!-- Scanheader variables
      integer(kind=knd2) :: stime(6)
      integer :: iyear, imonth, iday, ihour, imin, isec
      character*15 :: cdate, ctime

!--- Declaration of internal subroutine variables.

      integer :: ipix, lin, chan, i, j
      integer :: ret=0
      logical :: fexists
      
!----------------------------------------------------------
        npix = npixl 

!---check if file exists -------------------
      
      inquire(file=input_file,exist = fexists)
      if(fexists .eq. .false.) call reprt_error(10)           !input file does not exist
     
!
! Initialize FORTRAN interface.
!
      CALL h5open_f(hdferr)
      if(hdferr .eq. -1) call reprt_error(9)                   !error in open             
!
! Open file and attribute.
!
      CALL h5fopen_f(input_file, H5F_ACC_RDONLY_F, file_id, hdferr)
      if(hdferr .eq. -1) call reprt_error(9)                   !error in open             
!
!-----------------
!
! Read in NumberofScans a variable-length string attribute
!
      ALLOCATE(rdata(1:dims(1)))
      attrname = "NumberOfScans"
      CALL h5aopen_f(file_id, attrname, attr_id, hdferr)
      f_ptr = C_LOC(rdata(1))
      CALL h5aread_f(attr_id, H5T_STRING, f_ptr, hdferr)
!
      CALL C_F_POINTER(rdata(1), data)
      len = 0
      DO WHILE(DATA(len+1:len+1).NE.C_NULL_CHAR)
         len = len + 1
      ENDDO
      read (data(1:len),*) nscans
      CALL h5aclose_f(attr_id , hdferr)   ! close and release resources
!
!-----------------
!
! Read in Number of OverlapScans a variable-length string attribute 
!
      attrname = "OverlapScans"
      CALL h5aopen_f(file_id, attrname, attr_id, hdferr)
      f_ptr = C_LOC(rdata(1))
      CALL h5aread_f(attr_id, H5T_STRING, f_ptr, hdferr)
!
      CALL C_F_POINTER(rdata(1), data)
      len = 0
      DO WHILE(DATA(len+1:len+1).NE.C_NULL_CHAR)
         len = len + 1
      ENDDO
      read (data(1:len),*) oscans
      CALL h5aclose_f(attr_id , hdferr)   ! close and release resources
!
!-----------------
!
! Read in Start Orbit Number a variable-length string attribute 
!
      attrname = "StartOrbitNumber"
      CALL h5aopen_f(file_id, attrname, attr_id, hdferr)
      f_ptr = C_LOC(rdata(1))
      CALL h5aread_f(attr_id, H5T_STRING, f_ptr, hdferr)
!
      CALL C_F_POINTER(rdata(1), data)
      len = 0
      DO WHILE(DATA(len+1:len+1).NE.C_NULL_CHAR)
         len = len + 1
      ENDDO
      read (data(1:len),*) granule
      CALL h5aclose_f(attr_id , hdferr)   ! close and release resources
      write(*,*)'  Orbit number     :  ',granule
!
!-----------------
!
! Read in Observation Start DateTime a variable-length string attribute 
!   YYYY-MM-DDTHH:MM:SS.SSSZ      2013-01-01T00:19:49.420Z
! NOTE (DD, 9/1/15) -- 'ScanTime' is also in files, a double that gives
!  seconds since 0:00 1/1 1993 (TAI)
!
      attrname = "ObservationStartDateTime"
      CALL h5aopen_f(file_id, attrname, attr_id, hdferr)
      f_ptr = C_LOC(rdata(1))
      CALL h5aread_f(attr_id, H5T_STRING, f_ptr, hdferr)
!
      CALL C_F_POINTER(rdata(1), data)
      len = 0
      DO WHILE(DATA(len+1:len+1).NE.C_NULL_CHAR)
         len = len + 1
      ENDDO
      read (data(1:10),*)  cdate                            ! YYYY-MM-DD 
      read (data(12:19),*) ctime                            ! HH:MM:SS
      CALL h5aclose_f(attr_id , hdferr)   ! close and release resources
      write(*,*)'  Start Date        :  ',cdate
      write(*,*)'  Start Time        :  ',ctime
!
!-----------------
!
! Read in GranuleID, filename,  a variable-length string attribute 
!     GW1AM2_201301010019_185D_L1SGRTBR_1110110
!
      attrname = "GranuleID"
      CALL h5aopen_f(file_id, attrname, attr_id, hdferr)
      f_ptr = C_LOC(rdata(1))
      CALL h5aread_f(attr_id, H5T_STRING, f_ptr, hdferr)
!
      CALL C_F_POINTER(rdata(1), data)
      len = 0
      DO WHILE(DATA(len+1:len+1).NE.C_NULL_CHAR)
         len = len + 1
      ENDDO
      read (data(1:len),*) radfile
      CALL h5aclose_f(attr_id , hdferr)   ! close and release resources
      write(*,*)'  radfile        :  ',radfile
!
!-----------------
!
! Read in Number of Satellite Altitude a variable-length string
! attribute 
!
      attrname = "SatelliteAltitude"
      CALL h5aopen_f(file_id, attrname, attr_id, hdferr)
      f_ptr = C_LOC(rdata(1))
      CALL h5aread_f(attr_id, H5T_STRING, f_ptr, hdferr)
!
      CALL C_F_POINTER(rdata(1), data)
      len = 0
      DO WHILE(DATA(len+1:len+1).NE.C_NULL_CHAR)
         len = len + 1
      ENDDO
      read (data(1:5),*) alt
      CALL h5aclose_f(attr_id , hdferr)   ! close and release resources
      write(*,*)'  Altitude in KM :  ',alt
!
!-----------------

!--- Need the number of scans in the data set arrays
        dsetname = "Brightness Temperature (original,89GHz-A,H)"
        CALL h5dopen_f (file_id, dsetname, dset_id, hdferr)
        CALL h5dget_space_f(dset_id, space_id, hdferr)
        CALL h5sget_simple_extent_dims_f(space_id,ddims,dmaxdims,hdferr)
!---  ddims(1)=243 or 486 number of pixels in a scan
!---  ddims(2)=2017 or so number of scans saved
        CALL h5sclose_f(space_id, hdferr)
        CALL h5dclose_f(dset_id , hdferr)
        sscans=ddims(2)
      write(*,*)'  Num of unique scans       (nscans) : ',nscans
      write(*,*)'  Num of overlap scans      (oscans) : ',oscans
      write(*,*)'  Num of scans in inputfile (sscans) : ',sscans

!-----------------
 
!--- allocate all input data and pixel flag arrays
      allocate (Tb(npixl,nscans,maxchans))      
      allocate (eia(npixl,nscans,maxchans))
      allocate (lat(npixl, nscans))
      allocate (lon(npixl, nscans))

      allocate(sclat(nscans))
      allocate(sclon(nscans))
      allocate(scalt(nscans))
      allocate(scorient(nscans))
      allocate(stdtime(nscans,6))
      allocate(sglinta(npix,nscans))
      allocate(sunglint(npix,nscans))
      allocate(qualflag(npix,nscans))
        qualflag = 0 ! maybe get from l1r later?
      !allocate(qflg(npix,nscans,5)) 

      ! new for L1R version of pp (vars needed for sglint calc, esp)
      allocate (satAZ(npixl,nscans))
      allocate (lofl(npixl,sscans,4),lo_flag(npixl,nscans,4)) 
      allocate (sunAZ(npixl,nscans),sunEL(npixl,nscans))
      allocate(sun_az(npixl,sscans),sun_el(npixl,sscans))
      allocate(az_low(npixl,sscans))

!--- ScanTime in seconds from 1993 
      allocate(t93t(sscans),tai93time(nscans)) ! real*8

!-- low-res channels

      allocate(tb6v_low(npixl,sscans))
      allocate(tb6h_low(npixl,sscans))
      allocate(tb10v_low(npixl,sscans))
      allocate(tb10h_low(npixl,sscans))
      allocate(tb19v_low(npixl,sscans))
      allocate(tb19h_low(npixl,sscans))
      allocate(tb23v_low(npixl,sscans))
      allocate(tb23h_low(npixl,sscans))
      allocate(tb37v_low(npixl,sscans))
      allocate(tb37h_low(npixl,sscans))
      allocate(tb89v_low(npixl,sscans))
      allocate(tb89h_low(npixl,sscans))

      allocate(eia_low(npixl,sscans))

!--- high res arrays for low-res channels 
!--- for the final number of scans

      allocate(tb6v(npixl,nscans))
      allocate(tb6h(npixl,nscans))
      allocate(tb10v(npixl,nscans))
      allocate(tb10h(npixl,nscans))
      allocate(tb19v(npixl,nscans))
      allocate(tb19h(npixl,nscans))
      allocate(tb23v(npixl,nscans))
      allocate(tb23h(npixl,nscans))
      allocate(tb37v(npixl,nscans))
      allocate(tb37h(npixl,nscans))
      allocate(tb89v(npixl,nscans))
      allocate(tb89h(npixl,nscans))

!--- high res pixels for original 89A lat/lon

      !allocate(tb89vA(npixh, sscans),tb89hA(npixh, sscans))

      allocate(lat_i2(npixh, sscans))  ! lat/lon from 89A
      allocate(lon_i2(npixh, sscans))  ! lat/lon from 89A

      j = len_trim(radfile)
      do i = 1,j
        if(iachar(radfile(i:i)) .eq. 0)   exit
      enddo
      radfile = radfile(1:i-1)

      write(*,*) '  sat_name       : ',sat_name
      write(*,*) '  sensor_name    : ',sensor_name
      write(*,*) '  radiometer file: ',trim(radfile)
      write(*,*) '  altitude km    : ',alt


      ! read in Land/Ocean flag (% of land in channel FOV)
        dsetname = "Land_Ocean Flag 6 to 36"
        lodims(1) = npixl !npix!4
        lodims(2) = nscans !sscans
        lodims(3) = 4!npix 
        CALL h5dopen_f (file_id, dsetname, dset_id, hdferr)
        CALL h5dread_f(dset_id,H5T_NATIVE_INTEGER,lofl,ddims,hdferr)
        if(hdferr .eq. -1) call reprt_error(58)           ! error
        CALL h5dclose_f(dset_id, hdferr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These are the channels for AMSR2-L1R
! 6.925GHz 7.3GHz 10.65GHz 18.7GHz 23.8GHz 36.5GHz 89.0GHz-A 89.0GHz-B
! All are pixels = 243 except original 89 have pixels = 486
! Code it with 243 pixels, leaving out every other pix for Original 89
! and lat/lon
! Brightness Temperatures have SCALE FACTOR = 0.01 
! Going to ignore 89-B for now
! Going to leave out freqs 6.9 and 7.3 GHz
!   ----- NO LONGER! (for 6.9GHz)
!
! Four choices of resolution: 
!  1) res06  
!  2) res10  
!  3) res23  
!  4) res36  
!  x) original
!
      select case(res_opt)
      case("res06")                ! res06
        rz = 1
        dsetname6h =  "Brightness Temperature (res06,6.9GHz,H)"
        dsetname6v =  "Brightness Temperature (res06,6.9GHz,V)"
        dsetname10h = "Brightness Temperature (res06,10.7GHz,H)"
        dsetname10v = "Brightness Temperature (res06,10.7GHz,V)"
        dsetname19h = "Brightness Temperature (res06,18.7GHz,H)"
        dsetname19v = "Brightness Temperature (res06,18.7GHz,V)"
        dsetname23h = "Brightness Temperature (res06,23.8GHz,H)"
        dsetname23v = "Brightness Temperature (res06,23.8GHz,V)"
        dsetname37h = "Brightness Temperature (res06,36.5GHz,H)"
        dsetname37v = "Brightness Temperature (res06,36.5GHz,V)"
        dsetname89h = "Brightness Temperature (res06,89.0GHz,H)"
        dsetname89v = "Brightness Temperature (res06,89.0GHz,V)"
      case("res10")                ! res10
        rz = 2
        dsetname6h =  "Brightness Temperature (res06,6.9GHz,H)"
        dsetname6v =  "Brightness Temperature (res06,6.9GHz,V)"
        dsetname10h = "Brightness Temperature (res10,10.7GHz,H)"
        dsetname10v = "Brightness Temperature (res10,10.7GHz,V)"
        dsetname19h = "Brightness Temperature (res10,18.7GHz,H)"
        dsetname19v = "Brightness Temperature (res10,18.7GHz,V)"
        dsetname23h = "Brightness Temperature (res10,23.8GHz,H)"
        dsetname23v = "Brightness Temperature (res10,23.8GHz,V)"
        dsetname37h = "Brightness Temperature (res10,36.5GHz,H)"
        dsetname37v = "Brightness Temperature (res10,36.5GHz,V)"
        dsetname89h = "Brightness Temperature (res10,89.0GHz,H)"
        dsetname89v = "Brightness Temperature (res10,89.0GHz,V)"
      case("res23")                ! res23
        rz = 3
        dsetname6h =  "Brightness Temperature (res06,6.9GHz,H)"
        dsetname6v =  "Brightness Temperature (res06,6.9GHz,V)"
        dsetname10h = "Brightness Temperature (res10,10.7GHz,H)"
        dsetname10v = "Brightness Temperature (res10,10.7GHz,V)"
        dsetname19h = "Brightness Temperature (res23,18.7GHz,H)"
        dsetname19v = "Brightness Temperature (res23,18.7GHz,V)"
        dsetname23h = "Brightness Temperature (res23,23.8GHz,H)"
        dsetname23v = "Brightness Temperature (res23,23.8GHz,V)"
        dsetname37h = "Brightness Temperature (res23,36.5GHz,H)"
        dsetname37v = "Brightness Temperature (res23,36.5GHz,V)"
        dsetname89h = "Brightness Temperature (res23,89.0GHz,H)"
        dsetname89v = "Brightness Temperature (res23,89.0GHz,V)"
      case("res36")                ! res36
        rz = 4
        dsetname6h =  "Brightness Temperature (res06,6.9GHz,H)"
        dsetname6v =  "Brightness Temperature (res06,6.9GHz,V)"
        dsetname10h = "Brightness Temperature (res10,10.7GHz,H)"
        dsetname10v = "Brightness Temperature (res10,10.7GHz,V)"
        dsetname19h = "Brightness Temperature (res23,18.7GHz,H)"
        dsetname19v = "Brightness Temperature (res23,18.7GHz,V)"
        dsetname23h = "Brightness Temperature (res23,23.8GHz,H)"
        dsetname23v = "Brightness Temperature (res23,23.8GHz,V)"
        dsetname37h = "Brightness Temperature (res36,36.5GHz,H)"
        dsetname37v = "Brightness Temperature (res36,36.5GHz,V)"
        dsetname89h = "Brightness Temperature (res36,89.0GHz,H)"
        dsetname89v = "Brightness Temperature (res36,89.0GHz,V)"
!      case("original")             ! original
!        rz = 1 !ie if original res, use 6/7GHz FOV screening
!        dsetname6h =  "Brightness Temperature (res06,6.9GHz,H)"
!        dsetname6v =  "Brightness Temperature (res06,6.9GHz,V)"
!        dsetname10h = "Brightness Temperature (res10,10.7GHz,H)"
!        dsetname10v = "Brightness Temperature (res10,10.7GHz,V)"
!        dsetname19h = "Brightness Temperature (res23,18.7GHz,H)"
!        dsetname19v = "Brightness Temperature (res23,18.7GHz,V)"
!        dsetname23h = "Brightness Temperature (res23,23.8GHz,H)"
!        dsetname23v = "Brightness Temperature (res23,23.8GHz,V)"
!        dsetname37h = "Brightness Temperature (res36,36.5GHz,H)"
!        dsetname37v = "Brightness Temperature (res36,36.5GHz,V)"
!        dsetname89hA= "Brightness Temperature (original,89GHz-A,H)"
!        dsetname89vA= "Brightness Temperature (original,89GHz-A,V)"
      end select
      write(*,*)'  ResolutionOption ',rz,') ', res_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--- Read variables

        dsetname  = "Latitude of Observation Point for 89A"
        CALL h5dopen_f (file_id, dsetname, dset_id, hdferr)
        ddims(1) = npixh     ! 243 or 486 number of pixels in a scan
!---    ddims(2) = sscans  = 2017 number of scans saved in data file
        CALL h5dread_f(dset_id, H5T_NATIVE_REAL, lat_i2, ddims, hdferr)
      if(hdferr .eq. -1) call reprt_error(68)           ! error in latitude 
        CALL h5dclose_f(dset_id, hdferr)

        dsetname  = "Longitude of Observation Point for 89A"
        CALL h5dopen_f (file_id, dsetname, dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_REAL, lon_i2, ddims, hdferr)
      if(hdferr .eq. -1) call reprt_error(68)           ! error in longitude 
        CALL h5dclose_f(dset_id, hdferr)

      write(*,*)'   lat_i2(1,1)       :  ',lat_i2(1,1)
      write(*,*)'   lat_i2(2,1)       :  ',lat_i2(2,1)
      write(*,*)'   lat_i2(392,1)     :  ',lat_i2(392,1)
      write(*,*)'   lon_i2(1,1)       :  ',lon_i2(1,1)
      write(*,*)'   lon_i2(1,21)      :  ',lon_i2(1,21)
      write(*,*)'   lon_i2(2,1)       :  ',lon_i2(2,1)
      write(*,*)'   lon_i2(392,1)     :  ',lon_i2(392,1)

        ddims(1) = npixl     ! 243 or 486 number of pixels in a scan
        ddims(2) = sscans    ! 2017 number of scans saved in data file

! NO 6GHz in this version!

        CALL h5dopen_f (file_id, dsetname6h, dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb6h_low, ddims,hdferr)
      if(hdferr .eq. -1) call reprt_error(60)           ! error
        CALL h5dclose_f(dset_id, hdferr)

        CALL h5dopen_f (file_id, dsetname6v, dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb6v_low, ddims,hdferr)
      if(hdferr .eq. -1) call reprt_error(59)           ! error
        CALL h5dclose_f(dset_id, hdferr)

        CALL h5dopen_f (file_id, dsetname10h, dset_id, hdferr)
        CALL h5dread_f(dset_id,H5T_NATIVE_INTEGER,tb10h_low,ddims, hdferr)
      if(hdferr .eq. -1) call reprt_error(60)           ! error
        CALL h5dclose_f(dset_id, hdferr)

        CALL h5dopen_f (file_id, dsetname10v, dset_id, hdferr)
        CALL h5dread_f(dset_id,H5T_NATIVE_INTEGER,tb10v_low,ddims,hdferr)
      if(hdferr .eq. -1) call reprt_error(59)           ! error
        CALL h5dclose_f(dset_id, hdferr)

        CALL h5dopen_f (file_id, dsetname19h, dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb19h_low, ddims,hdferr)
      if(hdferr .eq. -1) call reprt_error(62)           ! error
        CALL h5dclose_f(dset_id, hdferr)

        CALL h5dopen_f (file_id, dsetname19v, dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb19v_low, ddims,hdferr)
      if(hdferr .eq. -1) call reprt_error(61)           ! error
        CALL h5dclose_f(dset_id, hdferr)

        CALL h5dopen_f (file_id, dsetname23h, dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb23h_low, ddims,hdferr)
      if(hdferr .eq. -1) call reprt_error(63)           ! error
        CALL h5dclose_f(dset_id, hdferr)

        CALL h5dopen_f (file_id, dsetname23v, dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb23v_low, ddims,hdferr)
      if(hdferr .eq. -1) call reprt_error(63)           ! error
        CALL h5dclose_f(dset_id, hdferr)

        CALL h5dopen_f (file_id, dsetname37h, dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb37h_low, ddims,hdferr)
      if(hdferr .eq. -1) call reprt_error(65)           ! error
        CALL h5dclose_f(dset_id, hdferr)

        CALL h5dopen_f (file_id, dsetname37v, dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb37v_low, ddims,hdferr)
      if(hdferr .eq. -1) call reprt_error(64)           ! error
        CALL h5dclose_f(dset_id, hdferr)

      if (rz .eq. 4) then                 ! then 89 has npix high
!        ddims(2)=npixh
!        CALL h5dopen_f (file_id, dsetname89hA, dset_id, hdferr)
!        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb89hA, ddims,
!     >                  hdferr)
!      if(hdferr .eq. -1) call reprt_error(70)           ! error
!        CALL h5dclose_f(dset_id, hdferr)
!
!        CALL h5dopen_f (file_id, dsetname89vA, dset_id, hdferr)
!        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb89vA, ddims,
!     >                  hdferr)
!      if(hdferr .eq. -1) call reprt_error(69)           ! error
!        CALL h5dclose_f(dset_id, hdferr)
      else
        ddims(2)=npixl
        CALL h5dopen_f (file_id, dsetname89h, dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb89h_low, ddims,hdferr)
      if(hdferr .eq. -1) call reprt_error(70)           ! error
        CALL h5dclose_f(dset_id, hdferr)

        CALL h5dopen_f (file_id, dsetname89v, dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, tb89v_low, ddims,hdferr)
      if(hdferr .eq. -1) call reprt_error(69)           ! error
        CALL h5dclose_f(dset_id, hdferr)
      endif

        dsetname = "Earth Incidence"
        ddims(2)=npixl
        CALL h5dopen_f (file_id, dsetname, dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, eia_low, ddims,hdferr)
      if(hdferr .eq. -1) call reprt_error(58)           ! error
        CALL h5dclose_f(dset_id, hdferr)

        dsetname = "Earth Azimuth"
        ddims(2)=npixl
        CALL h5dopen_f (file_id, dsetname, dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, az_low, ddims, hdferr)
      if(hdferr .eq. -1) call reprt_error(58)           ! error
        CALL h5dclose_f(dset_id, hdferr)

        dsetname = "Sun Azimuth"
        ddims(2)=npixl
        CALL h5dopen_f (file_id, dsetname, dset_id, hdferr)
        CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER,sun_az,ddims,hdferr)
      if(hdferr .eq. -1) call reprt_error(58)           ! error
        CALL h5dclose_f(dset_id, hdferr)

        dsetname = "Sun Elevation"
        ddims(2)=npixl
        CALL h5dopen_f (file_id, dsetname, dset_id, hdferr)
        CALL h5dread_f(dset_id,H5T_NATIVE_INTEGER,sun_el,ddims,hdferr)
      if(hdferr .eq. -1) call reprt_error(58)           ! error
        CALL h5dclose_f(dset_id, hdferr)

        dsetname = "Scan Time"
        dims=sscans
        CALL h5dopen_f (file_id, dsetname, dset_id, hdferr)
        CALL h5dread_f(dset_id,H5T_NATIVE_REAL_8, t93t ,dims,hdferr)
      if(hdferr .eq. -1) call reprt_error(55)           ! error
        CALL h5dclose_f(dset_id, hdferr)

      do lin=1,nscans
        tai93time(lin) = t93t(lin+oscans)

        if(tai93time(lin) .GT. 0.0) then
          call convtime_amsr(tai93time(lin), stime)
          do j = 1, 6
            stdtime(lin, j) = stime(j)
          enddo
!          itime = int(time(lin))
        else
          do j = 1, 6
            stdtime(lin, j) = MISS_INT
          enddo
        endif

      enddo    ! lin nscans

      write(*,'(a,f20.2)')' time sec,scan1   :',tai93time(1)
      write(*,'(a,f20.2)')' time sec,scan2   :',tai93time(2)
      write(*,'(a,f20.2)')' time sec,scan1001:',tai93time(1001)
      write(*,'(a,i4,1x,5(i2,1x))')' scan1 time     : ',stdtime(1,:)
      write(*,'(a,i4,1x,5(i2,1x))')' scan2 time     : ',stdtime(2,:)
      write(*,'(a,i4,1x,5(i2,1x))')' scan1001 time  : ',stdtime(1001,:)

!-- break dowwn the date & time char arrays into integers

!      read (cdate(1:4),*)  iyear
!      read (cdate(6:7),*)  imonth
!      read (cdate(9:10),*) iday
!      read (ctime(1:2),*)  ihour
!      read (ctime(4:5),*)  imin
!      read (ctime(7:8),*)  isec
!      write(*,*) 'iyear = ',iyear
!      write(*,*) 'imonth= ',imonth
!      write(*,*) 'iday  = ',iday
!      write(*,*) 'ihour = ',ihour
!      write(*,*) 'imin  = ',imin
!      write(*,*) 'isec  = ',isec
!
!      do lin=1,nscans
!        stdtime(lin,1)    = iyear       !YYYY
!        stdtime(lin,2)    = imonth      !MM
!        stdtime(lin,3)    = iday        !DD
!        stdtime(lin,4)    = ihour       !HH 
!        stdtime(lin,5)    = imin        !MM
!        stdtime(lin,6)    = isec        !SS
!      enddo

      write(*,'(a,i4,1x,5(i2,1x))')'   first time     : ', stdtime(1,:)
      
!--- check for good scans

      if(nscans .eq. 0) call reprt_error(19)         ! nscans = 0

!    Terminate access to the SD interface and close the file.

        CALL h5fclose_f(file_id , hdferr)
      if(hdferr .ne. 0) call reprt_error(8)                !error in close             

      return
      end  subroutine read_L1R_ppam2

!------------------------------------------------------------------------------

      subroutine fill_low2high_qctbs           
       implicit    none

       integer  ::  chan, lin, pix, i, j

!--- apply scale factors and offsets

        write(*,*)'starting fl2hc'
      do lin=1,nscans     ! nscans = 1977
        i=lin+oscans      ! oscans = 20  overlap scans
        sclat(lin)  = lat_i2(243,i)*(1.00)      !A-SCAN lat/lon, using mid point 
        sclon(lin)  = lon_i2(243,i)*(1.00)      !A-SCAN lat/lon, using mid point
        scalt(lin)  = alt                       !in KM
        scorient(lin)  = 0                      !AMSR2 always fwd!
        do pix = 1,npixl
          lat(pix,lin) = lat_i2(pix*2,i) *(1.00)  ! put lat and lon into arrays
          lon(pix,lin) = lon_i2(pix*2,i) *(1.00)
          !--- copy 89 if option 4, original resolution
          !if (rz .eq. 4) then                 ! then 89 has npix
          !high
          !  tb89v(pix,lin)    = tb89vA(pix,i)         !89V_Ascan at
          !  npix
          !  tb89h(pix,lin)    = tb89hA(pix,i)         !89H_Ascan at
          !  npix
          !endif
        enddo
        do j = 1, npixl!/2 
          pix=j!*2-1
          tb6v(pix,lin)     = tb6v_low(j,i)          !6V
          tb6h(pix,lin)     = tb6h_low(j,i)          !6H
          tb10v(pix,lin)    = tb10v_low(j,i)         !10V
          tb10h(pix,lin)    = tb10h_low(j,i)         !10H
          tb19v(pix,lin)    = tb19v_low(j,i)         !19V
          tb19h(pix,lin)    = tb19h_low(j,i)         !19H
          tb23v(pix,lin)    = tb23v_low(j,i)         !23V
          tb23h(pix,lin)    = tb23h_low(j,i)         !23H
          tb37v(pix,lin)    = tb37v_low(j,i)         !37V
          tb37h(pix,lin)    = tb37h_low(j,i)         !37H           

!-- Same Earth_Incidence angle for each channel

          eia(pix,lin,:)    = eia_low(j,i) *(0.01)
          satAZ(pix,lin)    = az_low(j,i)*(0.01)+180.0
          !lo_flag(pix,lin)  = lofl(j,i,1) !6GHz FOV hardcoded!
          lo_flag(pix,lin,1)  = lofl(j,i,1) !6GHz FOV
          lo_flag(pix,lin,2)  = lofl(j,i,2) !10GHz FOV
          lo_flag(pix,lin,3)  = lofl(j,i,3) !19/23GHz FOV
          lo_flag(pix,lin,4)  = lofl(j,i,4) !36GHz FOV
          !lo_flag(pix,lin)  = lofl(j,i,rz)

          sunAZ(pix,lin)    = sun_az(j,i)*0.01
          sunEL(pix,lin)    = sun_el(j,i)*0.01

!--- replicate each pixel across the low-res scan to make it high-res

!          tb10v(pix+1,lin)    = tb10v_low(j,i)         !10V
!          tb10h(pix+1,lin)    = tb10h_low(j,i)         !10H
!          tb19v(pix+1,lin)    = tb19v_low(j,i)         !19V
!          tb19h(pix+1,lin)    = tb19h_low(j,i)         !19H
!          tb23v(pix+1,lin)    = tb23v_low(j,i)         !23V
!          tb23h(pix+1,lin)    = tb23h_low(j,i)         !23H
!          tb37v(pix+1,lin)    = tb37v_low(j,i)         !37V
!         tb37h(pix+1,lin)    = tb37h_low(j,i)         !37H           


!--- replicate 89 if not option 4, original resolution

          if (rz .ne. 4) then                 ! then 89 has npix high
            tb89v(pix,lin)    = tb89v_low(j,i)         !89V
            tb89h(pix,lin)    = tb89h_low(j,i)         !89H           
            !tb89v(pix+1,lin)  = tb89v_low(j,i)         !89V
            !tb89h(pix+1,lin)  = tb89h_low(j,i)         !89H           
          endif
        enddo
      enddo

        write(*,*)'more starting fl2hc'

!--- replicate each pixel across the low-res scan to make it high-res

        do lin = 1, nscans
          do pix = 1, npixl
            Tb(pix,lin,:) = miss_flt ! initialize!
            Tb(pix,lin,1) = tb6v(pix,lin)  *(0.01)
            Tb(pix,lin,2) = tb6h(pix,lin)  *(0.01)
            Tb(pix,lin,3) = tb10v(pix,lin) *(0.01)
            Tb(pix,lin,4) = tb10h(pix,lin) *(0.01)
            Tb(pix,lin,5) = tb19v(pix,lin) *(0.01)
            Tb(pix,lin,6) = tb19h(pix,lin) *(0.01)
            Tb(pix,lin,7) = tb23v(pix,lin) *(0.01)
            Tb(pix,lin,8) = tb23h(pix,lin) *(0.01)
            Tb(pix,lin,9) = tb37v(pix,lin) *(0.01)
            Tb(pix,lin,10)= tb37h(pix,lin) *(0.01)
            Tb(pix,lin,11)= tb89v(pix,lin) *(0.01)
            Tb(pix,lin,12)= tb89h(pix,lin) *(0.01)
          enddo
        enddo
      write(*,*)'   lat(1,1)       :  ',lat(1,1)
      write(*,*)'   lon(1,1)       :  ',lon(1,1)
      write(*,*)'   azimuth(1,1)   :  ',satAZ(1,1)
      write(*,*)'   lon_i2(1,21)   :  ',lon_i2(1,21)
      write(*,*)'   tb10v(1,1)      :  ',Tb(1,1,1)
      write(*,*)'   tb10h(1,1)      :  ',Tb(1,1,2)
      write(*,*)'   tb89v(1,1)     :  ',Tb(1,1,9)
      write(*,*)'   tb89h(1,1)     :  ',Tb(1,1,10)
      write(*,*)'   lat(201,1)     :  ',lat(201,1)
      write(*,*)'   lon(201,1)     :  ',lon(201,1)
      write(*,*)'   lo_flag(201,1) :  ',lo_flag(201,1,1)
      write(*,*)'   azimuth(201,1) :  ',satAZ(201,1)
      write(*,*)'   tb10v(201,1)    :  ',Tb(201,1,1)
      write(*,*)'   tb10h(201,1)    :  ',Tb(201,1,2)
      write(*,*)'   tb89v(201,1)   :  ',Tb(201,1,9)
      write(*,*)'   tb89h(201,1)   :  ',Tb(201,1,10)

!--- quality control the Tbs

        do chan = 1,maxchans
          do lin = 1,nscans
            do pix = 1, npixl
              if(Tb(pix,lin,chan) .lt. 0) Tb(pix,lin,chan) = miss_flt
            enddo
          enddo
        enddo  

      return
      end subroutine fill_low2high_qctbs

!----------------------------------------------------------------------

      subroutine convtime_amsr(amsr_time, time)
      implicit none

      real(kind=knd8)    :: amsr_time
      integer            :: lin, i
      integer(kind=knd2) :: time(6)                 ! YYYY MM DD HH MM

      call amsr_to_std(amsr_time, time)

      return
      end subroutine convtime_amsr

!--------------------------------------------------------------------------------

      subroutine amsr_to_std(sec93, time)
      implicit none

      real(kind=knd8)   :: sec93
      integer(kind=knd2):: time(6)

      integer,parameter :: sec_per_day = 86400
      integer,parameter :: sec_per_hour= 3600
      integer,parameter :: sec_per_min = 60

      integer            ::   sec
      integer(kind=knd2) ::   year, month

      sec = sec93
      year = 1993

      do while ((sec - sec_in_year(year)) .GE. 0)
          sec  = sec - sec_in_year(year)
          year = year + 1
      enddo
      time(1) = year

      month = 1
      do while ((sec - sec_in_month(month,year)) .GE. 0)
          sec  = sec - sec_in_month(month,year)
          month= month + 1
      enddo
      time(2) = month

      ! days
      time(3) = sec / sec_per_day + 1
      sec = mod(sec, sec_per_day)

      ! hours
      time(4) = sec / sec_per_hour
      sec = mod(sec, sec_per_hour)

      ! minutes
      time(5) = sec / sec_per_min
      sec = mod(sec, sec_per_min)

      ! seconds
      time(6) = sec

      return
      end subroutine amsr_to_std

!----------------------------------------------------------------------------


      integer function sec_in_month(month, year)
      implicit none

      !logical :: isleapyr
      integer*2,intent(in) :: month, year
      integer,parameter ::   sec_per_day = 86400
      integer   ::  days(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)


      if (month .NE. 2) then
            ! All months except February
            sec_in_month = days(month) * sec_per_day
      else
            ! February
            if (isleapyr(year)) then
                  sec_in_month = 29 * sec_per_day
            else
                  sec_in_month = 28 * sec_per_day
            end if
      end if
      return
      end function sec_in_month

!---------------------------------------------------------------------------

      integer function sec_in_year(year)
      implicit none

      integer*2,intent(in) :: year
      !logical :: isleapyr

      if (isleapyr(year)) then
            sec_in_year = 31622400
      else
            sec_in_year = 31536000
      end if

      return
      end function sec_in_year
!_______________________-------------------____________________________
      logical function isleapyr(year)
      implicit none

      integer*2,intent(in) :: year

      if(((mod(year,4)==0).AND.(mod(year,100)/=0)).OR.(mod(year,400)==0))then
        isleapyr = .TRUE.
      else
        isleapyr = .FALSE.
      end if

      return
      end function isleapyr

!---------------------------------------------------------------------------
      end module read_l1r
