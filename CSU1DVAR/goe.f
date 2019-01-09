      program goe
!
!----------------------------------------------------------------
! integrated amsr2 algorithm -- intam2
! 
! version history
! begun 1/25/16 D Duncan, CSU
! last update 12/7/17 D Duncan, Chalmers
!
! master program that calls all routines
!
! name changed to 'goe' (gprof/oe) d duncan, Jan 2017
!
!----------------------------------------------------------------

       use define_intam2
       use define_csu1dvar
       use GPM_arraydef_D
       use csu1dvar
       use csu1dvar_subr
       use GPROF_2017_V1
       use output_1dvar_nc
       use integrated_output

      implicit none

      integer   :: n_arguments, iargc, ios, olun, c, d, stat,s,p
      integer   :: gdat, isfc
      integer(2):: sflag
      real :: la,lo
      real :: invalid, taz
      real :: sic
      real :: Tb(maxchans) ! for sic retr 
      integer :: q_1dvar
      character(len=256) :: gprof_log, gprof_anc
      integer :: nout = 15 ! number of orbit swath variables output

        INTERFACE
          SUBROUTINE icecon_nt2_(lat2,lon2,v12,h12,v22,v32,v82,h82,
     >                          sic2,gdat2,xsi2,ysi2,inv2,stat2) BIND(C)
            USE, INTRINSIC :: ISO_C_BINDING
            IMPLICIT NONE
            REAL(C_FLOAT) :: lat2, lon2
            REAL(C_FLOAT) :: v12,h12,v22,v32,v82,h82
            REAL(C_FLOAT) :: inv2,sic2
            INTEGER(C_INT) :: gdat2,xsi2,ysi2,stat2
          END SUBROUTINE icecon_nt2_
        END INTERFACE

        gprof_log = 'test_log'
        gprof_anc = 'binary/'
        !gprof_anc = '/cdata1/dduncan/intam2/'

        invalid = miss_flt
        gdat = 0

!--- Get command line inputs - <infile1 <infile2> <outfiles> 

      n_arguments = iargc()
      if(n_arguments .eq. 5) then
        call getarg(1, in_file1)      !input filename (pp file1)
        call getarg(2, in_file2)      !input filename (pp file2)
        call getarg(3, out_file)      !output filename (all)
        call getarg(4, gprof_out)     !output filename (precip)
        call getarg(5, csu1dvar_out)  !output filename (ocean)
      else    ! not correct number of arguments, n_arguments .ne. 3 or 2
         write(*,*)' error in command line arguments '
         write(*,*)' <inputfile1> <inputfile2> <out_file> 
     >               <gprof_out> <csu1dvar_out>'
         stop
      endif

!--- echo startup information

      !write(*,*)' input file name     : ', trim(in_file1)
      !write(*,*)' input file name (2) : ', trim(in_file2)
      !write(*,*)' output file name    : ', trim(out_file)
      !!!write(*,*)' resolution option   : ', trim(res_opt)
        
!--- read Tb data from L1R file
      !write(*,*)' starting read_pp'
      input_file = trim(in_file1)
      call read_geos5_pp(input_file,1) ! new pp call
      !write(*,*)' done read_pp'

! calculate azimuthal angle if non-AMSR
      if(satcode.eq.'GMI' .or. satcode.eq.'TMI') then ! AM2 has azimuth in L1 files!!
        allocate(scad(nscans))
        ! new section for winds read-in...
        scad(1) = 1 ! starts ascending, say
        scad(nscans) = -1 ! ends descending, say
        do s = 2, nscans-1
          if((sclat(s+1)-sclat(s-1)).gt.0) scad(s)=1
          if((sclat(s+1)-sclat(s-1)).lt.0) scad(s)=-1
          if((sclat(s+1)-sclat(s-1)).eq.0) scad(s)=1 ! ?
        enddo
        do s = 1, nscans
          do p = 1, npix
            call calc_az(sclat(s),p,scad(s),scorient(s),taz)
            sataz(p,s) = taz
          enddo
        enddo
        deallocate(scad)
      endif


       allocate(outarr(npix,nscans,nout))
       outarr(:,:,:)  = -99 ! initialize
       outarr(:,:,13) = 0
       outarr(:,:,14) = -1 ! Soil Moisture
       outarr(:,:,15) = miss_flt ! Soil Moisture
       !allocate(torun_ocean(npix,nscans) ! set to all true at first
       allocate(torun_gprof(npix,nscans)) ! set to all true at first
       torun_gprof(:,:) = .true.
       allocate(torun_land(nscans,npix)) ! set to all true at first
       torun_land(:,:) = .true.
       !torun_ocean(:,:) = .true.
       allocate(seaice(npix,nscans))
       seaice(:,:) = miss_flt

!--- call GPROF, run as normal
      !write(*,*)' calling GPROF2017v1'
      !call run_gprof(gprof_out,gprof_log,gprof_anc) ! outfile,logfile,anc_dir
      torun_gprof(:,:) = .false. ! now that all pixels run, set to false
      !write(*,*)' GPROF run'

      !write(*,*)' importing sst data'
      call read_reynolds_sst
      !write(*,*)' assigning sst data'
      call assign_sst
      !write(*,*)' prepping 1dvar'
      call prep_oe
      !write(*,*)' initialize crtm'
      call crtm_innit
              
      OPEN(UNIT=5,file='binary/binary_MIE_original_ALL_mis-99_2pct.tbl',
     >  access='stream',status='old')
      READ(5) dielec_rlist_ALL
      READ(5) dielec_ilist_ALL
      READ(5) crefindex_ALL
      READ(5) xlist_ALL
      READ(5) mietable_ALL_db
      CLOSE(5)
      !write(*,*)' calling 1dvar'
      call run_csu1dvar
        !print*,'outputting 1dvar'
       call output_nc ! output CSU1DVAR data to netcdf file
      stop

!!!!!!!!!! what follows is not part of 1DVAR warm rain retrieval, but
!part of integrated code


      !torun_ocean(:,:) = .false. ! now that all pixels run, set to false
      outarr(:,:,7)  = oe_output(:,:,5) ! chisq 
      outarr(:,:,8)  = oe_output(:,:,1) ! TPW
      outarr(:,:,9)  = oe_output(:,:,2) ! LWP
      outarr(:,:,10) = oe_output(:,:,3) ! WINd
      outarr(:,:,11) = sst(:,:) ! just Reynolds!! !oe_output(:,:,6) ! SST
      !write(*,*)' 1dvar complete'

      ! Land (SM) algorithm:
      !ancdir = '/cdata1/dduncan/AM2_sciteam/AMSR2_Land/anc/' !hard-coded!
      !l2bdir = 'out/' ! output directory, normally command line argument
      !pmcver = '4' ! hard-coded, 'product maturity code', normally command line argument

        !print*,'starting SM'
!      call sm_main
        !print*,'done SM'

      !write(*,*)' starting read_pp (again)'
        input_file = trim(in_file2)
      call read_geos5_pp(input_file,2) ! new pp call
      !write(*,*)' done read_pp (again)'
!--- start looping through pixels, applying Tb thresholds
        !print*,'npix,nscans: ',npix,nscans
        do lin = scstart, scend !1, 400 !nscans
          do pix = pxstart, pxend !1, npix
            
            ! assign cutoffs to lo_flag to define land/coast/ocean
            sflag = lo_flag(pix,lin,3) ! use 19GHz for SIC!

            Tb(:) = tbs(pix,lin,:)

            la = lat(pix,lin)
            lo = lon(pix,lin)

            ! flip soil moisture around:
!            outarr(pix,lin,14) = sca_vsm(lin+oscans,pix)
!            if(sflag > 50) ! simple for now, should screen more!
!     >       outarr(pix,lin,15) = stemp_sca(lin+oscans,pix)

            sic = 0 ! initialize every time
            if(sfcprcp(pix,lin).ge.2.0 .or. pop(pix,lin).gt.0.5)
     >          seaice(pix,lin) = -96 ! precip present
            if(sst(pix,lin).le.270 .or. sst(pix,lin).ge.278) 
     >          seaice(pix,lin) = -98 ! too warm, bad sst
            if(sflag.gt.50) seaice(pix,lin) = -97
            if(sflag.le.50 .and. 
     >          ((sst(pix,lin)<=278.0 .and. la >=  30).or.
     >           (sst(pix,lin)<=275.0 .and. la <= -40)) .and.
     >          sst(pix,lin).gt.270 .and. pop(pix,lin).le.0.5 .and. 
     >          sfcprcp(pix,lin).lt.2.0) then 
         call icecon_nt2_(la,lo,Tb(3),Tb(4),Tb(5),Tb(7),Tb(9),Tb(10),
     >                  sic,gdat,1,1,invalid,stat)
              !print*,'SIC: ',sic,lin,pix,sst(pix,lin)
              seaice(pix,lin) = sic 
            endif

        ! just for test runs, no rerunning other algs to save time!
       call output_nc ! output CSU1DVAR data to netcdf file
        !print*,'outputting gprof'
       call GPM_output ! output GPROF data structure to binary file
       !print*,'outputting integrated'
       call output_int
        stop

      ! rerun gprof for confirmed sea ice if it didn't have right sfc
      ! type before (and sea ice / coast GPM sfc type):
        if(sic.gt.0) then 
          if(sflag .gt. 5  .and. sfccode(pix,lin).ne.14) then
            sfccode(pix,lin) = 14 !GPROF land/ice code
            torun_gprof(pix,lin) = .true. ! turn rerun pix on
            outarr(pix,lin,13) = 1 ! flag to rerun gprof for ice/land
          endif
          if(sflag .le. 5 .and. sfccode(pix,lin).ne.2) then
            sfccode(pix,lin) = 2 !GPROF ice code
            torun_gprof(pix,lin) = .true. ! turn rerun pix on
            outarr(pix,lin,13) = 2 ! flag to rerun gprof for sea ice
          endif
        ! if 1dvar confirms open ocean, override sea ice
          if(oe_output(pix,lin,1) .gt. 0 .and. 
     >       oe_output(pix,lin,5) .lt. 4.0) then
            seaice(pix,lin) = 0
            outarr(pix,lin,13) = 3 ! simple flag for rerun/override (for now)
          endif
        endif
      ! rerun gprof for confirmed open ocean (i.e. 1dvar converged)
        if(oe_output(pix,lin,1).gt.0 .and. sfccode(pix,lin).ne.1 
     >     .and. oe_output(pix,lin,5).lt.2.0) then
          sfccode(pix,lin) = 1 !GPROF ocean code
          torun_gprof(pix,lin) = .true. ! turn rerun pix on
          outarr(pix,lin,13) = 4 ! simple flag for rerun/override (for now)
        endif
      ! provision for rerunning 1dvar? how and why?
        if(oe_output(pix,lin,1).gt.0 .and. pop(pix,lin).gt.0.5
     >     .and. sfcprcp(pix,lin).ge.0.5) then
          outarr(pix,lin,7:10) = miss_flt
          outarr(pix,lin,13) = 5 ! flag for override 1dvar from gprof
        endif
      ! 'grayzone' flag for pixels with no 1dvar convergence, no rain,
      ! no sea ice detected
        if( (oe_output(pix,lin,1).lt.0.or.oe_output(pix,lin,5).gt.4)
     >    .and. pop(pix,lin) .le. 0.5 .and. sfcprcp(pix,lin) .lt. 0.5
     >    .and. sic .eq. 0) then
          outarr(pix,lin,13) = 7 ! flag for ocean gray zone!
        endif


          enddo ! pixel loop
        enddo ! scan loop

        ! call GPROF again for problem pixels
        !print*,'calling GPROF again'
       !do isfc = minsfccode, maxsfccode
        do isfc = 1,2 !minsfccode, maxsfccode
         !print*,isfc
         call GPM_read_dbase(isfc)
         if(EIAoffset) call GPM_calc_EIA_offset
         call GPM_rain(isfc)
       enddo  !isfccode
       outarr(:,:,1) = sfcprcp(:,:)
       outarr(:,:,2) = frzprcp(:,:)
       outarr(:,:,3) = cwp(:,:)
       outarr(:,:,4) = iwp(:,:)
       outarr(:,:,5) = pop(:,:) * 100.0
       outarr(:,:,6) = real(sfccode(:,:))
       print*,'done gprof again'

        outarr(:,:,12) = seaice(:,:)

!--- call qual_czech afterward to verify the ocean product altogether
        !call qual_czech

!        call gprof_lun(olun)
!        open(unit=olun,file=trim(out_file),access='stream',
!     >       status='unknown',iostat=ios)
!        if(ios.ne.0)write(*,*)' problem outputting file'
!        
!        write(unit=olun,iostat=ios) real(npix)
!        write(unit=olun,iostat=ios) real(nscans)
!        write(unit=olun,iostat=ios) lat(:,:)
!        write(unit=olun,iostat=ios) lon(:,:)
!        write(unit=olun,iostat=ios) lo_flag(:,:,1)
!        do c = 1, nout
!          write(unit=olun,iostat=ios) outarr(:,:,c)
!        enddo
!        write(unit=olun,iostat=ios) save_slp(:,:)
!        write(unit=olun,iostat=ios) Tb_diff(:,:,1:nch)
!        !write(unit=olun,iostat=ios) Tb_diff(:,:,7:10)
!        close(unit=olun)

        ! deallocate arrays from memory
        !deallocate(sst,lat,lon,outarr)!,time,lat_i2,lon_i2) 
        !deallocate(lo_flag,lofl) 
        
        !print*,'outputting 1dvar'
       call output_nc ! output CSU1DVAR data to netcdf file
        !print*,'outputting gprof'
       call GPM_output ! output GPROF data structure to binary file
       !print*,'outputting integrated'
       call output_int

       !write(*,*)'COMPLETED'

      end
