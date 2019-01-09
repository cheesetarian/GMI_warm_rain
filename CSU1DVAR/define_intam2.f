      MODULE define_intam2

        use CRTM_MODULE

       implicit none
       save

!--- these are shared variables, used by more than one 
!---  subroutine and therefore commonly defined here
        
       character(12) :: code_version = 'goe_v10' ! version of code
       integer            :: log_lun     !log file logical unit number

       integer       :: ocut = 0 ! % land in fov cutoff for ocean definition
!       VARS TO BE MODIFIED FOR TESTING/RUNNING 
       !integer  :: scstart=2500, pxstart=200, scend=2800, pxend=201
       !integer  :: scstart=1, pxstart=97, scend=3200, pxend=123 !DPR swath
       !integer  :: scstart=1, pxstart=26, scend=3200, pxend=201
       integer  :: scstart=2945, pxstart=37, scend=2945, pxend=39
       real :: maxla=70,minla=-65,maxlo=180.0,minlo=-180.0
       !real :: maxla=48,minla=44,maxlo=-2,minlo=-8.0 ! france case!
       !real :: maxla=54,minla=46,maxlo=-14,minlo=-19 ! mid-atlantic/ireland case!
       !real :: maxla=60.5,minla=56.5,maxlo=-143.0,minlo=-149.0 !Alaska PAIH
!      Common definitions:
!______________________________________________________________________
       integer,parameter :: knd1 = 1
       integer,parameter :: knd2 = 2
       integer,parameter :: knd8 = 8
       integer,parameter :: maxchans = 15
       real                 :: miss_flt = -9999.9
       integer, parameter   :: miss_int = -9999    !missing integer value
       integer, parameter   :: miss_byt = -99      !missing byte value
       character(len=100) :: input_file,in_file1,in_file2,l1fi,out_file 
       character(len=256) :: csu1dvar_out, gprof_out
! --------
       integer :: npix, nscans, pix, lin, oscans !last one is new
       real, allocatable    :: lat(:,:), lon(:,:) ! (pix,scan)
       real, allocatable    :: sclat(:), sclon(:), scalt(:) !spacecraft lat/lon/altitude
       real,allocatable :: tbs(:,:,:)
       real,allocatable :: oetbs(:,:,:)
       real,allocatable :: eia(:,:,:)
       integer, allocatable :: stdtime(:,:) ! time
       real*8, allocatable :: tai93time(:) ! time in sec from 1993
       integer(1),allocatable :: lo_flag(:,:,:) !landmask array
       !integer(2),allocatable ::sfc_type(:,:),lo_flag(:,:) !landmask array
       !integer(1), allocatable :: sglint(:,:), qualflag(:,:) 
       integer(kind=knd1),allocatable :: sglinta(:,:)      !
       real,allocatable     :: sst(:,:),sstgrid(:,:) !ice and sst grids
       real,allocatable     :: seaice(:,:) ! sea ice saved
       real,allocatable     :: outarr(:,:,:) ! retrievals' output: RR,POP,SIC,CHISQ

       character(len=12) :: sensor_name

       integer              :: chan_avail(maxchans)
       integer              :: avail(maxchans) ! whether ch exists or not for 1dvar
       integer(4) :: granule

       logical,allocatable :: torun_ocean(:,:)
       logical,allocatable :: torun_land(:,:)
       logical,allocatable :: torun_gprof(:,:)

      END MODULE define_intam2
