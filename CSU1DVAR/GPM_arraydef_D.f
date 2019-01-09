      MODULE gpm_arraydef_D

!------ GPM array definition module, modified D Duncan to work with the
!AMSR2 team integrated code (mostly just commenting out variable
!declarations that are now done elsewhere)

       use CRTM_MODULE
       use define_intam2

       implicit none
       save
             
       !integer,parameter :: knd1 = 1    
       !integer,parameter :: knd2 = 2
       !integer,parameter :: knd8 = 8

!--- misc variables       

       !integer            :: log_lun         !log file logical unit number
       integer            :: nprfbin
       character(len=256) :: output_file,log_file
       !character(len=256) :: input_file,output_file,log_file
       character(len=256) :: dir_anc
       character(len=12)  :: alg_version
       
       !real,parameter     :: miss_flt = -9999.0
       !integer,parameter  :: miss_int = -9999
       integer,parameter  :: miss_int2= -999
       !integer,parameter  :: miss_byt = -99
       
       integer,parameter  :: minT2m    =220, maxT2m     = 320
       integer,parameter  :: mintcwv   =  0, maxtcwv    =  78   !78 for GPM V2
       integer,parameter  :: minsfccode=  1, maxsfccode =  14   !14 for GPM V2
       
       integer            :: histpp(minsfccode:maxsfccode)    !sfccode histogram
       integer            :: binnfcnt                     ! num bins not found
       logical            :: EIAoffset                    !apply eiaoffsets to dbase Tbs?
       logical            :: probcut                      !apply probability cutoffs to precips
       real               :: chansens(15,14)              !15 channel sensitivities for 14 sfccodes
       real               :: chansensfull(15,14,mintcwv:maxtcwv)   !15 freqs X 14 sfccodes X 0:78 tcwv       
       real               :: rainsnowtwb(3,131)           !rain snow twb array

!--- structure profile variables

       integer            :: profstructflag    !0=no profile,1=yes profile       
       integer,parameter  :: nspecies = 5
       integer,parameter  :: ntemps   = 12
       integer,parameter  :: nlyrs    = 28
       integer,parameter  :: nclusts  = 40
       integer,parameter  :: nprfs    = 80
       integer,parameter  :: maxclusters  = 800     
       !integer,parameter  :: maxchans   = 15       ! maximum channels
       
       real,allocatable   :: ProfileClusters(:,:,:,:)
       real,allocatable   :: ProfileClustersNR(:,:,:,:)
             
!---   define subset lat-lons variables

       real                :: latmin, latmax ,lonmin, lonmax      
	     
!---  Preprocessor data variables

       real               :: chan_freqs(maxchans)
       real               :: chan_errs(maxchans)
       character(len=12)  :: sat_name                !orbhdr
       !character(len=12)  :: sensor_name
       character(len=12)  :: pp_version
       character(len=128) :: radfile
       character(len=128) :: pdbfile
       character(len=128) :: calfile
       !integer            :: granule
       !integer            :: nscans
       integer            :: npixs
       integer            :: nchans
                 
       !integer,allocatable :: stdtime(:,:)            !scnhdr
       !real,allocatable    :: sclat(:)
       !real,allocatable    :: sclon(:)
       !real,allocatable    :: scalt(:)  
                                           	     
       !real,allocatable :: lat(:,:)
       !real,allocatable :: lon(:,:)
       !real,allocatable :: tbs(:,:,:)    
       !real,allocatable :: eia(:,:,:)
       real,allocatable :: Twb(:,:)
       !real,allocatable :: LR(:,:)
       real,allocatable :: ppskint(:,:)
       real,allocatable :: pptcwv(:,:)
       real,allocatable :: ppT2m(:,:)
       integer(kind=knd1),allocatable :: L1Cqualflag(:,:)
       
       !integer(kind=knd1),allocatable :: sglinta(:,:)      !pixel rec
       integer(kind=knd1),allocatable :: sfccode(:,:)
       integer(kind=knd2),allocatable :: CAPE(:,:) 
       
!---  database variables       
       
       real                           :: db_nominal_eia(15)       
       integer,allocatable            :: idx(:,:,:)         !pointer for skint/tcwv bins      
       integer,allocatable            :: db_nprofiles(:,:)
       
       integer,allocatable            :: db_occur(:)  
       integer,allocatable            :: db_raincount(:)     
       integer(kind=knd2),allocatable :: db_tbs(:,:)     
       integer(kind=knd2),allocatable :: db_tbdelta(:,:)
       real,allocatable               :: db_sfcprcp(:)   
       real,allocatable               :: db_cnvprcp(:)          
       real,allocatable               :: db_rwp(:)	
       real,allocatable               :: db_cwp(:)	
       real,allocatable               :: db_iwp(:)
       
       real,allocatable               :: db_tcwv(:)
       real,allocatable               :: db_T2m(:)       
       
       integer(kind=knd2),allocatable :: db_rwc(:,:)      
       integer(kind=knd2),allocatable :: db_cwc(:,:) 
       integer(kind=knd2),allocatable :: db_iwc(:,:)             
       integer(kind=knd2),allocatable :: db_swc(:,:)             
       integer(kind=knd2),allocatable :: db_gwc(:,:)      
              
!---  assigned and computed variables

       integer(kind=knd1),allocatable:: pixel_status(:,:)      !quality variables
       integer(kind=knd1),allocatable:: quality_flag(:,:)

       real,allocatable              :: sfcprcp(:,:)           !rain and rain diagnostics       
       real,allocatable              :: cnvprcp(:,:)
       real,allocatable              :: frzprcp(:,:)
       real,allocatable              :: pop(:,:)
       
       real,allocatable              :: mlprcp(:,:)
       real,allocatable              :: prcp1stT(:,:)
       real,allocatable              :: prcp2ndT(:,:)

       real,allocatable              :: rwp(:,:)             ! path variables
       real,allocatable              :: cwp(:,:)
       real,allocatable              :: iwp(:,:)

       real,allocatable              :: rwc(:,:,:)            ! profile species variables
       real,allocatable              :: cwc(:,:,:)
       real,allocatable              :: iwc(:,:,:)
       real,allocatable              :: swc(:,:,:)
       real,allocatable              :: gwc(:,:,:)

       integer(kind=knd2),allocatable:: prfT2mindex(:,:)
       integer(kind=knd2),allocatable:: prfnum(:,:,:)         !profile definitions
       real,allocatable              :: prfscale(:,:,:)
	
       !integer                       :: chan_avail(5)        !channel availability
       real,allocatable              :: tbdelta_eia(:,:,:)   !delta Tb offsets from nominal
       
!---  probability cutoff arrays

       real, allocatable ::  probcutoff(:,:,:)
       real, allocatable ::  missvolfrac(:,:,:)
       integer(kind=knd1), allocatable ::  probquality(:,:,:)       
     
     
      end module GPM_arraydef_D
