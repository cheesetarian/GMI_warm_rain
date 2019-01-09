      MODULE define_csu1dvar

        USE CRTM_MODULE

        implicit none
        save

!--- these are shared variables, used by more than one 
!---  subroutine and therefore commonly defined here
        
       character(10) :: csu1dvar_version = 'ver3.0' ! version of code
       character(10) :: rtm_version = 'CRTMv2.2.3' ! RT code. needs to match Makefile!


!______________________________________________________________________
       logical :: useanalysis    = .true.  ! if true, reads in analysis data
       !logical :: dontoutput     = .false. ! if true, code stops after OE step
       logical :: runseaice      = .false. ! if true, force SST=271 for pix w/ reynolds seaice
       logical :: includeicecloud= .false. ! should match IWP being retrieved or not!
       logical :: merra          = .false. ! merra first guess (not EC)
       logical :: retrievesst    = .false.  !switch for sst retrieval--overriden for GMI (for now)
       !logical :: outrt,run_rss,fastem4 ! set by script calling code (see oeorbit.f)

       real :: latdsd = 30.0 ! cutoff for low vs. high lat DSDs
       real :: dlat  ! global var to pass lat through opt_est/noscatter

!       VARS TO BE MODIFIED FOR TESTING/RUNNING 
!       integer  :: scstart=1, pxstart=1, scend=400, pxend=243
       ! if end values are too high, nscans/npix will be end values
       
       integer,parameter    :: nz = 16  !# layers in atmosphere -- not negotiable now
       integer,parameter    :: nbins = 33 , npc = 2 !CHANGE 
         ! SST bins, # PCs in EOF LUT files, SST bins for Sy matrices/offsets
       integer,parameter    :: nn = 8 ! max # of iterations to find convergence
       integer,parameter    :: nvar= 7 !# of parameters with potential to be retrieved
       integer              :: nretvar != 3+npc ! # params retrieved -overridden in read_geos5
       real                 :: conv_factor = 20.0 !value to divide nchan
                               !by to determine convergence (Rodgers Eq5.33)
       real                 :: varmult = 1.8 ! a priori variance multiplier (wsp only)
       real                 :: cld_effrad = 11 ! cld water droplet eff radius [um]
       real                 :: ice_effrad = 60.0 ! cld ice eff radius [um]
       real                 :: chis1 = 1.0 ! chisq threshold for 'High Quality'
       real :: chout = 4.0 !threshold above which output = missing
       !real :: chisq_out_thresh = 10.0 !threshold above which output = missing
! reminder of variable order: magn EOFs 1-3, WIND, log(LWP), IWP, SST.
!     Bounds on retrieved elements -- all that are capable of being retrieved!
      real :: x_min(nvar) = (/-4.0,-4.0,-4.0, 0.5,-5.0, 0.0, 270.0/)
      real :: x_max(nvar) = (/ 4.0, 4.0, 4.0,30.0,-0.3, 0.3, 307.0/) 
      !real :: x_max(nvar) = (/ 4.0, 4.0, 4.0,30.0,-0.2, 0.3, 307.0/) ! increased to 307K
        ! changed to [-5,-.2] LWP limits
      integer*2 :: lwmax(16,22) ! sst (2K) x tpw (3mm), max set at .3kg
        ! drizzle max/mins for coefficients come from curvefitting and 2%/98% values
      real :: drzmi1(11) = -6 !-3.0
      real :: drzma1(11) =  9!7.0
      real :: drzmi2(11) = -7 !-5.0
      real :: drzma2(11) =  11! 7.0
      real :: drzmi3(11) = -8 !-4.0
      real :: drzma3(11) =  10!7.0
      !real :: drzmi1(11) = (/-1.8,-1.8,-1.9,-1.9,-1.9,-2.0,-2.1,-2.1,-2.2,-2.3,-2.2/) !RWC1
      !real :: drzma1(11) = (/ 7.8, 7.4, 7.8, 6.7, 6.7, 6.1, 6.3, 6.0, 5.3, 6.4, 6.5/)
      !real :: drzmi2(11) = (/-3.8,-4.6, -16,-8.9,-9.9, -11,-9.2,-6.0,-8.0,-7.0,-2.7/) !RWC2
      !real :: drzma2(11) = (/ 4.8, 5.0,  10,  10, 9.0, 7.9, 5.8, 2.3, 6.3, 8.0, 3.7/)
      !real :: drzmi3(11) = (/-2.8,-2.8,-2.9,-2.9,-2.9,-2.8,-2.8,-2.7,-2.7,-2.3,-2.3/) !IWC1
      !real :: drzma3(11) = (/ 5.6, 6.1, 6.8, 6.8, 7.3, 7.8, 8.5, 9.3, 9.7, 8.3, 7.9/)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        real :: final_tpw, final_rwp, final_iwp, final_lwp,lower_pw
        !integer,parameter :: knd1 = 1
        !integer,parameter :: knd2 = 2
        !integer,parameter :: knd8 = 8
        !integer,parameter :: maxchans = 15
        !real              :: c_sigma(maxchans)    !channel errors ! --------
       real,allocatable   :: oe_tbs(:), test_out(:),test_out2(:),test_out3(:) ! just for testing sim vs obs tbs
       !real               :: oe_tbs(9), test_out(9)
       !character(len=256) :: orig_file
       !character(len=100) :: input_file, output_file, spec_file 
       !character(len=100) :: cal_file, lmskfile, sy_file
       character(len=100) :: cal_file, sy_file
       !character(len=10)  :: satellite, sensor, nchannels, nfreqs
       !character(len=20)  :: ppversion
       character(len=3)   :: satcode != 'AM2'
       character(len=10)  :: convolve
       !character(len=1)   :: ascdesc ! A or D, (N if neither) from filename
       integer            :: nch,nfreq !=11 !HARD CODED!! !
       !integer(2)         :: npix, nscans
       integer            :: oepix, oelin
       integer            :: first_good_date  !long integer string of YYYYMMDD'      
       !integer            :: avail(maxchans) ! whether ch exists or not
       !character(6)       :: freqs(maxchans) ! whether ch exists or not

!--- allocatable arrays used by multiple subroutines
       real, allocatable    :: Tb_diff(:,:,:),save_slp(:,:),Amat_out(:,:,:)
       !real, allocatable    :: Tbb(:,:,:), Tb_diff(:,:,:), Tb_sim(:,:,:)! (pix,scan,nchan)
       real, allocatable    :: poster(:,:,:)  ! posteriori error, stddev
       !real, allocatable    :: sat_eia(:,:,:)  ! earth incidence angle (pix,scan,chan)
       !real, allocatable    :: eia_out(:,:,:)  ! earth incidence angle (pix,scan,2)
       real, allocatable    :: save_sal(:,:),save_wdir(:,:) ! SLP, Wind dir
       !real, allocatable    :: save_clwc(:,:,:),save_ciwc(:,:,:) ! CLW/CIW
       !real, allocatable    :: lat(:,:), lon(:,:)!, sfc_type(:,:) ! (pix,scan)
       !real, allocatable    :: sclat(:), sclon(:), scalt(:) !spacecraft lat/lon/altitude
       real*8, allocatable    :: scorient(:) ! spacecraft forward/backward
       integer, allocatable :: scad(:) !spacecraft orbital direction, asc/descending
       !integer, allocatable :: stdtime(:,:) ! time
       !integer(1),allocatable ::lsmask(:,:),sfc_type(:,:),lo_flag(:,:) !landmask array
       integer(1), allocatable :: save_iter(:,:),save_sw(:,:)
       !integer(1), allocatable :: sglint(:,:), qualflag(:,:) ! new for BCMX
       !real,allocatable     :: sst(:,:),sstgrid(:,:) !ice and sst grids
       !logical, allocatable :: icegrid(:,:)
       !real,allocatable     :: freeqs(:) !frequencies (only needed if using RSS model)
       !integer    :: nfreeqs
       real,allocatable :: save_emis(:),save_tau(:,:)

       !real                 :: miss_flt = -9999.9
       !integer, parameter   :: miss_int = -9999    !missing integer value
       !integer, parameter   :: miss_byt = -99      !missing byte value

!--- OE variables passed back to retrieval
       real                 :: xr(nvar)  ! retrieval variables in OE
       real                 :: chisq  ! Chi Squared
       real, allocatable    :: oe_output(:,:,:),screen(:,:)
       real                 :: last_oeout(nvar) ! used for a closer first guess 
       integer              :: loeop=-9,loeos=-9 ! save pix/scan for fg (initialize)
       integer              :: run_count=0
       real                 :: avg_iter ! averge # iterations required for retrieval


!--- variables for ERA-derived LUT table of means/sigmas
      character(len=90) erafile
      integer*2,allocatable :: era_wm(:,:),era_ws(:,:)
      integer*2 :: ssdex
      real, allocatable :: mprof(:,:), eofs(:,:,:), peofs(:,:)
      real :: mrmp(nz)
      integer,parameter :: gm=4 !LUT grid multiplier -- set in LUT creation!
      integer,parameter :: nlon=512*gm,nlat=256*gm!,nmo=12
      real :: losize=.703125, lasize=.701760 ! need to match LUT creation!
      ! gsize is exact for lons, lats start at 89.463 and go by ~.701760
      real,allocatable :: sy(:,:), sy_i(:,:)
      real,allocatable :: toffsets(:),s_toffsets(:,:),s_sy(:,:,:)
      ! RWP Sy regression vars:
      real*8,allocatable :: rc_syl(:,:),rc_sym(:,:),rc_syh(:,:) 
      real*8,allocatable :: rs_syl(:,:),rs_sym(:,:),rs_syh(:,:) 
      real*8,allocatable :: rcl_syl(:,:),rcl_sym(:,:),rcl_syh(:,:) 
      real*8,allocatable :: rsl_syl(:,:),rsl_sym(:,:),rsl_syh(:,:) 
      real*8,allocatable :: rc_sgm(:,:),rs_sgm(:,:)  ! nchx3
      real*8,allocatable :: rc_r(:,:,:),rs_r(:,:,:) !nchxnchx3
      real*8,allocatable :: rcl_sgm(:,:),rsl_sgm(:,:)  ! nchx3
      real*8,allocatable :: rcl_r(:,:,:),rsl_r(:,:,:) !nchxnchx3
      real :: rc_rlo,rc_rme,rc_rhi,rs_rlo,rs_rme,rs_rhi
      real :: rcl_rlo,rcl_rme,rcl_rhi,rsl_rlo,rsl_rme,rsl_rhi
      real :: saver !save rain rate

!--- CRTM-specific shared variables
        INTEGER :: n_channels  ! dimension 'L'
        INTEGER :: n_profiles  ! dimension 'M'
        INTEGER :: n_sensors   ! dimension 'N'
        integer :: n_layers
        TYPE(CRTM_ChannelInfo_type), ALLOCATABLE :: chinfo(:) !N
        TYPE(CRTM_Geometry_type),ALLOCATABLE :: geo(:) !M
        TYPE(CRTM_Options_type), ALLOCATABLE :: opt(:) !M
        ! Forward declarations
        TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm(:) !M
        TYPE(CRTM_Surface_type), ALLOCATABLE :: sfc(:) !M
        TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:) !LxM      
       
! -- P-level vars, in meters: deltaz between layers, layer t/b height
        real :: lpress(nz+1)= (/100,  200, 300,  400,  500,  550,  600,&
            650,  700,  750,  800,  850, 900, 925, 950, 975, 1000/)
        real :: dp(nz) =     (/ 100,  100, 100,  100,   50,   50,   50,&
             50,   50,   50,   50,   50,  25,  25,  25,  25/)
        !real :: dz(nz) =     (/2965, 2669, 1997, 1629,  717,  665, 621,&
        !    582,  549,  518,  491,  466, 225, 219, 215, 210/)
        !real :: lz(nz+1) =   (/14822,11858,9189, 7192, 5563, 4846,4181,&
        !   3560, 2978, 2429, 1911, 1420, 954, 729, 509, 295, 85/)
        real :: dz(nz) =     (/4042, 2509, 1933, 1584,  697, 648, 605, &
            566, 532,  503,  475,  447, 214, 210, 206, 201/)
        real :: lz(nz+1) =   (/15437,11395,8886, 6953, 5369, 4672,4024,&
           3419, 2853, 2321, 1818, 1343, 896, 682, 472, 266, 65/)
        ! altitudes in [m], derived from mean geopotential at P levels
        real :: pressave(nz)

        real :: eof_lw_cov(3,nbins),mr_sigmas(nbins,3)

        ! calculating satellite azimuthal angle
        real,allocatable :: sataz(:,:)
        ! variables read in from binary analysis file
        real :: winddir(512,256), pwinddir
        real :: windmag(512,256), pwindmag
        real :: eclwp(512,256), peclwp
        real :: eciwp(512,256), peciwp
        real :: ecslp(512,256), pecslp
        real :: ect(512,256,nz+1), pect(nz+1), apect(nz)
        real :: ecmr(512,256,nz+1), pecmr(nz+1)
        real :: pplev(nz+1),pecz(nz+1)
        real :: sss(720,360), psss=0 ! sea surface salinity annual mean grid

        real,allocatable :: ppwind(:,:)
        real,allocatable :: save_tprof(:,:,:),pp_wvprof(:,:,:)
        integer(2),allocatable :: pph(:,:,:)

        ! new for drizz
        real :: rwpm(11,14),piwpm(11,14),rwpe(11,14,2),piwpe(11,14,2)
        real :: drizsa(3,3,11), lwmaxy
        integer :: drdx ! index for drizzle params
        real :: rcoef(11,2) !a,b coeffs for rr calculation
        !real :: rwpdif


! added for using eddington/mie stuff:
      REAL, PARAMETER :: EPI = 3.141592653589793
!--- variables for the mie scattering lookup table access
!-----------------

      real,parameter      :: ndensityALL_list=10

      real,parameter      :: tmin_ALL=223.15, tmax_ALL=303.15
      real,parameter      :: tinc_ALL=0.1

      real,parameter      :: xmin_ALL=1.0e-3, xmax_ALL=2.0e2
      real                :: xinc_ALL
      real,parameter      :: nxlist_ALL=184
      real                :: xlist_ALL(nxlist_ALL)

      real,parameter      :: dielec_rmin_ALL=1.05   !1.085046
      real,parameter      :: dielec_rmax_ALL=8.90   !8.722480
      real                :: dielec_rinc_ALL
      real,parameter      :: ndielec_rlist_ALL= 108 !5%=44  1%=206, 2%=108 increments
      real                :: dielec_rlist_ALL(ndielec_rlist_ALL)


      real,parameter      :: dielec_imin_ALL=6.0e-6  !6.8959e-6
      real,parameter      :: dielec_imax_ALL=3.15e+0  !3.044035
      real                :: dielec_iinc_ALL
      real,parameter      :: ndielec_ilist_ALL= 666    !5%=270 1%=1324, 2%=666 increments
      real                :: dielec_ilist_ALL(ndielec_ilist_ALL)


      type                :: crefindex_structure_ALL
          real            :: real_cref
          real            :: imag_cref
          real            :: loc_cref      ! This is the sorting INDEX
      end type crefindex_structure_ALL
      type(crefindex_structure_ALL) :: crefindex_ALL(ndielec_rlist_ALL,ndielec_ilist_ALL)

      real,parameter      :: ncrefindex_ALL = 16496.  !5%=3387.0 2%=16496  increments
      type                :: mietable_ALL_structure
          real            :: real_cref
          real            :: imag_cref
          real            :: xsize
          real            :: qsca
          real            :: qext
          real            :: asym
          real            :: qbsca
      end type mietable_ALL_structure
      type(mietable_ALL_structure) :: mietable_ALL_db(ncrefindex_ALL,nxlist_ALL)

        
      END MODULE define_csu1dvar
