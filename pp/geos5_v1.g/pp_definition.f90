      MODULE pp_definition
       implicit none
       save
       
!--- set variable definitions for Integer*2, Integer*1 and Real*8
       integer, parameter :: knd1 = 1
       integer, parameter :: knd2 = 2
       real,    parameter :: knd8 = 8
       real, parameter    :: miss_flt = -9999.0    !missing real value
       real, parameter    :: miss_int = -9999      !missing integer value 
       real, parameter    :: miss_byt = -99        !missing byte value 

!--- Set sensor specific variables

       character(len=12)  :: pp_version
       character(len=128) :: prfdbfile 

       character(len=3)   :: sens ! sensor code, given as input arg
       character(len=12)  :: sat_name != 'GCOM-W1'         !satellite name
       character(len=12)  :: sensor_name != 'AMSR2'       !sensor_name
       
       integer            :: nchans      != 11           !number of retrieval sat channels
       integer,parameter  :: maxchans    = 15           !num total sat channels in output
       
! define chan freq and availability = if channels drop out, then make it missing here
       
       real :: chan_freq(maxchans) != (/6.925,-9999,10.65,10.65,18.7,18.7,23.8,23.8,36.5,36.5,89.0,89.0,-9999,-9999,-9999/)
!       real :: chan_freq(maxchans) = (/10.65,10.65,18.7,18.7,23.8,23.8,36.5, 36.5, 89.0, 89.0,-9999.0,-9999.0,-9999.0,-9999.0,-9999.0/)
     
       logical :: chan_avail(maxchans)

      character(len=40) :: dsetname6h,  dsetname6v
      character(len=40) :: dsetname10h, dsetname10v
      character(len=40) :: dsetname19h, dsetname19v
      character(len=40) :: dsetname23h, dsetname23v
      character(len=40) :: dsetname37h, dsetname37v
      character(len=40) :: dsetname89h, dsetname89v
      character(len=43) :: dsetname89hA,dsetname89vA

!--- Set input and outfile names, and jobid so are common
 
       character(len=256) :: input_file, output_file  !read in as arguments
       character(len=128) :: jobid, cal_file
       character(len=256) :: gdir, gf1,gf2,gf3,gf4

       character(len=256) :: attrname, dsetname ! attribute name and dataset name
       !integer            :: apply_cal  = 0 ! apply cal table (1=yes, 0=no)
       ! User can decide which resolution to use (4 choices):
       !  1) res10  2) res23  3) res36  4) original
       character(len=10)  :: res_opt
       integer            :: rz  = 1 ! default is res10!

!-- misc variables

       integer :: nscans                        !number of scans read from file
       integer :: granule,savegran              !orbit number
       integer :: npix                          !number of pixels in scan 
       integer :: warning = 0                   !warning index

       integer :: oscans                        !number of overlap scans read from file (20)
       integer :: sscans                        !number of scans saved in data file   (2017) 
       integer :: npixh     = 486               !number of pixels in scan from file (high res)
       integer :: npixl     = 243               !number of pixels in scan from file (low res)
       real    :: alt                           !altitude of sat
       
!-- directory specifications for ancillary data-----------------------12
       
       character(len=128):: dir_ancil
       character(len=128):: dir_ingest
       character(len=96) :: dir_autosnow  = '/autosnow/'

!-- array for the surface map, and autosnow grid, emissclass grids

       integer(1) :: sfcgrid(5760,2880)  
       integer(1) :: snowgrid(9000,4500)
       integer(2) :: emissgrd(720,359,12) 

!-- high resolution model output grids

       real,allocatable :: LRhi(:,:)
       real,allocatable :: Twbhi(:,:)
       real,allocatable :: tcwvhi(:,:)
       real,allocatable :: tskinhi(:,:)
       real,allocatable :: T2mhi(:,:)
                            
!-- allocate variables for reading L1C
   
       character(len=23),  allocatable :: datetime(:)             !1 time for each scan
      
       character(len=128) :: radfile
       character(len=100) :: cgranule
       character(len=100) :: csatellite
       character(len=100) :: csensor
                       
!--- allocatable array declarations ------------------------------------        

       integer, allocatable :: stdtime(:,:)    !time
       real, allocatable    :: sclat(:)        ! spacecraft latitude [-90, +90]
       real, allocatable    :: sclon(:)        ! spacecraft longitude [-180, +180]
       real, allocatable    :: scalt(:)        ! spacecraft altitude (km) 
       real*8,allocatable   :: scorient(:)     ! spacecraft orientation (fwd/bkwd)
       real*8,allocatable   :: tai93time(:),t93t(:) ! double values,scantime in seconds


       real, allocatable :: Tb(:,:,:)          ! Tbs read from file
       real, allocatable :: lat(:,:)           ! pixel latitude [-90, +90]
       real, allocatable :: lon(:,:)           ! pixel longitude [-180, +180]
       real,allocatable  :: eia(:,:,:)         ! earth incident angle for each channel  

       integer,allocatable  :: lofl(:,:,:),lo_flag(:,:,:)      ! land/ocean pct
       real,allocatable  :: satAZ(:,:)
       integer,allocatable  :: eia_low(:,:)       ! earth incident angle
       integer,allocatable :: az_low(:,:)      ! azimuth angle (needed for wind dir)
       integer,allocatable  :: tb6v_low(:,:),tb6h_low(:,:)
       integer,allocatable  :: tb10v_low(:,:),tb10h_low(:,:)
       integer,allocatable  :: tb19v_low(:,:),tb19h_low(:,:)
       integer,allocatable  :: tb23v_low(:,:),tb23h_low(:,:)
       integer,allocatable  :: tb37v_low(:,:),tb37h_low(:,:)
       integer,allocatable  :: tb89v_low(:,:),tb89h_low(:,:)
       integer,allocatable  :: tb6v(:,:),tb6h(:,:)
       integer(2),allocatable  :: tb10v(:,:),tb10h(:,:)
       integer(2),allocatable  :: tb19v(:,:),tb19h(:,:)
       integer(2),allocatable  :: tb23v(:,:),tb23h(:,:)
       integer(2),allocatable  :: tb37v(:,:),tb37h(:,:)
       integer(2),allocatable  :: tb89v(:,:),tb89h(:,:)
       !integer(2),allocatable  :: tb89vA(:,:),tb89hA(:,:)
       real,allocatable  :: lat_i2(:,:),lon_i2(:,:)

       real, allocatable :: sunAZ(:,:)     ! sun azimuth (degrees)
       real, allocatable :: sunEL(:,:)     ! sun elevation (degrees)
       real, allocatable :: sunglint(:,:)  ! sunglint angle (degrees)
       integer, allocatable :: sun_az(:,:),sun_el(:,:)
      
       integer(kind=knd1), allocatable :: sglinta(:,:)    !sunglint angle        
       integer(kind=knd1), allocatable :: sfccode(:,:)    !surface type code
              
       real              , allocatable :: T2m(:,:)        !2 meter temperature
       integer(kind=knd2), allocatable :: snowci(:,:)     !snow cover index
       integer(kind=knd2), allocatable :: orolifti(:,:)   !orographic lifting index
       real,               allocatable :: Twb(:,:)        !2 meter wet bulb
       real,               allocatable :: LR(:,:)         !lowest 500m Lapse Rate
       real,               allocatable :: tcwv(:,:)       !total column water vapor index
       real,               allocatable :: skint(:,:)      !surface skin temp             
       integer*1,          allocatable :: sic(:,:)        !sea ice frac (%)          
       integer,            allocatable :: qflg(:,:,:)     !L1C quality flag all chans             
       integer,            allocatable :: qualflag(:,:)   !L1C output quality flag

       ! new for geos5 version
       real,      allocatable :: wind(:,:) !m/s
       integer(2),allocatable :: winddir(:,:) !deg from N
       real,      allocatable :: tprof(:,:,:) !K
       real,      allocatable :: wvmr(:,:,:) !g/kg
       integer(2),allocatable :: height(:,:,:) ! new
       real,      allocatable :: slp(:,:) !hPa
       integer,parameter :: nz=16

       ! for GMI
       integer(1),allocatable :: lsmask(:,:)
       
!-----------------------------------------------------------------------
      END MODULE pp_definition
