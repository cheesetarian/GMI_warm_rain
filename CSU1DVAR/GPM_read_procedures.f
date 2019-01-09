      module GPM_read_procedures
        ! not modified other than different arraydef include below,
        !  and a tweak to get right chansens filename
        ! got rid of GPM_read_preprocessor subr (1/20/17 dduncan)
        ! updated to GPROF2017_V1 1/24/17 dduncan
       use GPM_arraydef_D
       use GPM_util_procedures
       use define_intam2

       contains

!-------------------------------------------------------------------------

      subroutine GPM_init_arrays
       implicit none
      
       integer :: i
       
!---  initialization procedures to allocate retrieval arrays

!--- allocate memory for the output parameters
      
       allocate (quality_flag(npixs,nscans))
       allocate (sfcprcp(npixs,nscans))
       allocate (cnvprcp(npixs,nscans))
       allocate (frzprcp(npixs,nscans))
       
       allocate (mlprcp(npixs,nscans))
       allocate (pop(npixs,nscans))
       allocate (prcp1stT(npixs,nscans))
       allocate (prcp2ndT(npixs,nscans))

       allocate (rwp(npixs,nscans))
       allocate (cwp(npixs,nscans))
       allocate (iwp(npixs,nscans))
       
       allocate (rwc(npixs,nscans,nlyrs))
       allocate (cwc(npixs,nscans,nlyrs))
       allocate (iwc(npixs,nscans,nlyrs))
       allocate (swc(npixs,nscans,nlyrs))
       allocate (gwc(npixs,nscans,nlyrs))
       
       allocate (prfT2mindex(npixs,nscans))
       allocate (prfnum(npixs,nscans,nspecies))
       allocate (prfscale(npixs,nscans,nspecies))
       
       quality_flag   = 0
       
       sfcprcp  = miss_flt
       cnvprcp  = miss_flt
       frzprcp  = miss_flt       
       rwp      = miss_flt
       cwp      = miss_flt
       iwp      = miss_flt       
       mlprcp   = miss_flt
       pop      = miss_flt
       prcp1stT = miss_flt
       prcp2ndT = miss_flt
       rwc      = miss_flt
       cwc      = miss_flt
       iwc      = miss_flt
       swc      = miss_flt
       gwc      = miss_flt
       
       prfnum         = 0
       prfscale       = 1.0
             
       write(log_lun,*)'  init_arrays : ', npixs,nscans

      end subroutine GPM_init_arrays

!---------------------------------------------------------------------12

      subroutine GPM_read_chansens
       implicit none

!--- routine to read in the sensor specific channel sensitivity ascii file

       integer :: tlun, ios,i,j, tchan
       character(len=256) :: sensfile
       character(len=100) :: chansens_header
       character(len=100) :: desc
       
!---  open chansen file for specific sensor and version number

       sensfile = trim(dir_anc)// trim(pdbfile) // '_chansens.txt'
       
       
       call GPM_lun(tlun)       
       open(tlun,file=sensfile,form='formatted',status='old',iostat=ios)
       if(ios .ne. 0) call GPM_reprt_err(27,tlun,sensfile)     !err opening chan sens file
      
!--- read in header lines      

       do i = 1,15
        read(tlun,'(a)')  chansens_header      
c        write(log_lun,'(a,i2.2,a)')'   chansens file header[',i,'] = ', 
c     >                       trim(chansens_header)
c        if(i .eq. 1) then
c	  write(log_lun,'(a,a100)')'   ',trim(chansens_header(1:40))
c	endif
       enddo

!--- read in the NeDt values (1st data line)

       read(tlun,'(i7,15F7.2,a)',iostat=ios)tchan,chan_errs,desc      !1st line is sensor err
c       write(log_lun,'(i7,15F7.2,a)') tchan,chan_errs(:), trim(desc) 
	        
!--- read in the 15 channels X 14 surface types
       
       do i = 1,14                                                     !14 surface classes V1-5
         read(tlun,'(i7,15F7.2,a)',iostat=ios) tchan,chansens(:,i),desc
c 	 write(log_lun,'(i7,15F7.2,a)') tchan,chansens(:,i),trim(desc)
       enddo         
       if(ios .ne. 0) call GPM_reprt_err(28,tlun,sensfile)     !err reading chan sens file
              
!--- fill entire 'full' table with chansens for each TCWV - just a replication process for now,
!--- later, the chansen will depend on tcwv, and a new full table will be read in

       do i = mintcwv, maxtcwv
         chansensfull(:,:,i) = chansens(:,:)    !chansensfull= (15 freqs,14 sfccode,0:78 tcwvs)
       enddo	      
	      	      
       close(tlun)
	     
      return
      end subroutine GPM_read_chansens

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      
      subroutine GPM_hist_preprocessor
       implicit none
     
!--- this routine will check the orbit preprocessor data and creates
!--- a sfccode histogram.  Used to know if all 14 surface codes need to be
!--- used in the GPM_rain subroutine.

       integer :: isfccode, iscan, ipix, tp,pw,sc
       character(len=128) :: blank = ' '
       logical :: chkhist
       
!--- loop over and hist all scans and pixels

       histpp = 0
       chkhist = .false.
              
       do iscan = 1, nscans
         do ipix = 1, npixs
	   if(pixel_status(ipix,iscan) .ne. 0) cycle        !ignore bad points	   	   
	   isfccode = sfccode(ipix,iscan)
	   
	   if(isfccode .lt. 1 .or. isfccode .gt. 14) then
	       write(log_lun,*)'bad sfccode: ',ipix,iscan,isfccode
	   endif	   
	   histpp(isfccode) = histpp(isfccode) + 1	   
	   chkhist = .true.
	 enddo   !ipix
       enddo  !iscan

D      write(log_lun,*)' histogram sfccode information'
D      do isfccode = minsfccode,maxsfccode
D        write(log_lun,*)' ', isfccode,'=',histpp(isfccode)
D      enddo

       if(.not. chkhist) call GPM_reprt_err(17,isfccode,blank)    !err no hist sfctypes

       return
      end subroutine GPM_hist_preprocessor

!-----------------------------------------------------------------------

      subroutine GPM_read_dbase(rsfccode)
       implicit none
       
       integer :: sfccode
       integer(kind=knd8) :: cntclusters
       integer(kind=knd2) :: trwc, tcwc, tiwc, tswc, tgwc
       
       integer  :: db_lun, ios, i, j, k, n, nf
       integer  :: maxbins, ibin, nclustbin
       character(len=256) :: dbfile
       character(len=5)   :: satcode_dtb, sensor_dtb
       character(len=2)   :: csfccode
       integer            :: rsfccode, isfccode, iT2m, itcwv       
       character(len=128) :: blank = ' '
       integer :: ntotprf, idxcnt = 0, pos=0, zero=0

!--- allocate memory for profile numbers in database

      allocate(db_nprofiles(minT2m:maxT2m,mintcwv:maxtcwv))

!--- define database filename 

      write(csfccode,'(i2.2)')  rsfccode
      dbfile = trim(dir_anc)// trim(pdbfile)// '_' // csfccode // '.dtb'

!--- Read Tb profile database

      call GPM_lun(db_lun)
D     write(log_lun,'(a,a)')'   reading profile database      : ',
D    >                 trim(dbfile)
      open ( unit=db_lun, file=dbfile, access='stream',
     >       status='old',iostat=ios)
      if(ios.ne.0) call GPM_reprt_err(40,db_lun,trim(dbfile))    !err reading prf dbase

!--- read total number of profiles from 1st line in Database
     
      read(db_lun,iostat=ios) ntotprf
      if(ios.ne.0) call GPM_reprt_err(41,ntotprf,blank)         !err reading totprfs prf dbase
D     write(log_lun,*)'   number clusters in DB header : ',ntotprf

!--- allocate memory for T2m/tcwv prof index file

      allocate(idx(minT2m:maxT2m,mintcwv:maxtcwv,maxclusters))
      
!---- allocate memory for each variable, for all profiles
      
      allocate(db_occur(ntotprf))
      allocate(db_raincount(ntotprf))
      allocate(db_tbs(ntotprf,maxchans))
      allocate(db_tbdelta(ntotprf,maxchans))
      allocate(db_sfcprcp(ntotprf))
      allocate(db_cnvprcp(ntotprf))
      allocate(db_rwp(ntotprf))
      allocate(db_cwp(ntotprf))
      allocate(db_iwp(ntotprf))
           
      allocate(db_tcwv(ntotprf))
      allocate(db_T2m(ntotprf)) 
       
      if(profstructflag .eq. 1) then
          allocate(db_rwc(ntotprf,nlyrs))
          allocate(db_cwc(ntotprf,nlyrs))
          allocate(db_iwc(ntotprf,nlyrs))	
          allocate(db_swc(ntotprf,nlyrs))	
          allocate(db_gwc(ntotprf,nlyrs))
      endif

!--- read in database
	 
      maxbins = (maxT2m-minT2m+1) * (maxtcwv-mintcwv+1)
D     write(log_lun,*)'   max num of bins to read in   : ', maxbins

      idxcnt = 0
      cntclusters = 0           
      do ibin = 1, maxbins

        read(db_lun,iostat=ios) satcode_dtb,sensor_dtb	      !sat code at top of each bin
c	write(log_lun,*) rsfccode,satcode_dtb, sensor_dtb 

	if(trim(satcode_dtb) .eq. 'AM2') satcode_dtb = 'AMSR2'
	if(trim(satcode_dtb) .eq. 'AME') satcode_dtb = 'AMSRE'
	if(trim(satcode_dtb) .eq. 'TMIPO') satcode_dtb = 'TMI'
	
c	if(trim(satcode_dtb) .ne. trim(sensor_name)) then
c	    write(log_lun,*)' satcodes=', satcode_dtb,' ',sensor_name
c	    call GPM_reprt_err(43,db_lun,blank)              !err mis-match in prf dbase sensor
c        endif

        read(db_lun,iostat=ios) db_nominal_eia	             !satellite sensors nominal eia
        if(ios.ne.0) call GPM_reprt_err(44,db_lun,blank)     !err reading db_nominal_eia
c	write(log_lun,'(15F8.2)') db_nominal_eia	
		
	read(db_lun,iostat=ios) nclustbin,isfccode,iT2m,itcwv    !num clusters in bin
        if(ios.ne.0) call GPM_reprt_err(45,isfccode,blank)       !err reading bin header info
c	write(log_lun,*) nclustbin,isfccode,iT2m,itcwv,maxclusters	
	
	if(nclustbin.gt.maxclusters) then      
          write(log_lun,*) '  nclustbin = ', nclustbin
          call GPM_reprt_err(46,nclustbin,blank)                !maxcluster check err
	endif
	
	db_nprofiles(iT2m,itcwv) = nclustbin	
        if(db_nprofiles(iT2m,itcwv) .gt. 0) then    !if there are profiles are in this bin
	   do n = 1, db_nprofiles(iT2m,itcwv)
	     idxcnt = idxcnt + 1	     
             idx(iT2m,itcwv,n) = idxcnt            !store in a linear array
	     
	     if(profstructflag .eq. 1) then
                read(db_lun,iostat=ios) 
     >           db_occur(idxcnt),
     >           db_raincount(idxcnt),
     >           db_tbs(idxcnt,:),
     >           db_tbdelta(idxcnt,:),
     >           db_sfcprcp(idxcnt),
     >           db_cnvprcp(idxcnt), 
     >           db_rwp(idxcnt),
     >           db_cwp(idxcnt),
     >           db_iwp(idxcnt),
     >           db_tcwv(idxcnt),
     >           db_T2m(idxcnt),
     >           db_rwc(idxcnt,:),
     >           db_cwc(idxcnt,:),
     >           db_iwc(idxcnt,:), 
     >           db_swc(idxcnt,:), 
     >           db_gwc(idxcnt,:)
             else
	        read(db_lun,iostat=ios) 
     >           db_occur(idxcnt),
     >           db_raincount(idxcnt),
     >           db_tbs(idxcnt,:),
     >           db_tbdelta(idxcnt,:),
     >           db_sfcprcp(idxcnt),
     >           db_cnvprcp(idxcnt), 
     >           db_rwp(idxcnt),
     >           db_cwp(idxcnt),
     >           db_iwp(idxcnt),
     >           db_tcwv(idxcnt),
     >           db_T2m(idxcnt),
     >           (trwc,k=1,nlyrs),             !if no layers, 'throw away'
     >           (tcwc,k=1,nlyrs),             !the layered database values
     >           (tiwc,k=1,nlyrs),
     >           (tswc,k=1,nlyrs),      
     >           (tgwc,k=1,nlyrs)                       
	     
	     endif    
             cntclusters = cntclusters + db_occur(idxcnt)
	     if(ios.ne.0) call GPM_reprt_err(47,idxcnt,blank)
	   enddo
	   
	endif
      enddo       
	    
D     write(log_lun,*)'   idxcnt cluster count         : ', idxcnt
D     write(log_lun,*)'   total num prfs in clusters   : ', cntclusters
       
      close(db_lun)
      return            
      end subroutine GPM_read_dbase

!------------------------------------------------------------------------------------

      subroutine  GPM_read_probability_cutoff
       implicit none
!
!---  this routine will read in the file which defines the probability of
!---  precip (pop) cutoff and rain rate fraction decrease that the cutoff 
!---  creates.  The POP is then applied to the surface precip, frozen precip
!---  and convective precip values in the GPM_apply_probability_cutoff routine.

!--- written by Dave Randel  Colorado State University
!--- October, 2016
       
       integer :: i,j,k, plun, ios
       integer :: pf_t2min,  pf_t2max,  pf_tcwvmin, pf_tcwvmax
       integer :: pf_sfcmin, pf_sfcmax, pf_probmin, pf_probmax
       character(len=256) :: probcutfile

       if(.not. probcut) then
           write(log_lun,*)' Probability Cutoffs NOT read'
	   return
       else
           write(log_lun,*)' Probability Cutoffs read'
       endif

!---  create filename from ancillary directory and sensor_name

       probcutfile = trim(dir_anc)// trim(sensor_name)// '_' // 
     >                  'probfrac.dat'
       write(log_lun,*)' probability cutoff file : ',trim(probcutfile)

!---  open probability cutoff file
       
       call GPM_lun(plun)
       open(unit=plun,file=trim(probcutfile),access='stream',
     >      status='old',iostat = ios, readonly)     
       if(ios .ne. 0) call GPM_reprt_err(50,plun,probcutfile)    !opening prob cut file

!---  read file header info

       read(plun,iostat=ios) pf_t2min,pf_t2max,pf_tcwvmin,pf_tcwvmax,
     >                       pf_sfcmin,pf_sfcmax,pf_probmin,pf_probmax
       if(ios .ne. 0) call GPM_reprt_err(51,plun,probcutfile)    !read prob cut file header

       write(log_lun,'(a,8i6)')'   probfile header= ',pf_t2min,pf_t2max,
     >                        pf_tcwvmin,pf_tcwvmax,pf_sfcmin,pf_sfcmax,
     >                        pf_probmin,pf_probmax
     
!---  allocate probability arrays

       allocate (probcutoff(pf_t2min:pf_t2max,pf_tcwvmin:pf_tcwvmax,
     >                      pf_sfcmin:pf_sfcmax))
       allocate (missvolfrac(pf_t2min:pf_t2max,pf_tcwvmin:pf_tcwvmax,
     >                      pf_sfcmin:pf_sfcmax))
       allocate (probquality(pf_t2min:pf_t2max,pf_tcwvmin:pf_tcwvmax,
     >                      pf_sfcmin:pf_sfcmax)) 

!--- read in the cutoff arrays

       read(plun,iostat=ios) probcutoff, missvolfrac, probquality
       if(ios .ne. 0) call GPM_reprt_err(52,plun,probcutfile)    !read prob cut file arrays       
 
       close(plun)

       return       
      end subroutine GPM_read_probability_cutoff      
       
!-----------------------------------------------------------------------

      subroutine GPM_read_profile_clusters
       implicit none

       integer            :: rlun, ispec,itemp,iclust, ios       
       character(len=128) :: fname 
       character(len=128) :: fspec(2) = (/'GPM_profile_clustersV3.dat',
     >                                 'GPM_profile_clustersNRV3.dat'/)
                                      
       allocate (ProfileClusters(nspecies,ntemps,nlyrs,nclusts))
       allocate (ProfileClustersNR(2,ntemps,nlyrs,nclusts))

!--- initialize cluster array

       profileclusters   = miss_flt
       profileclustersNR = miss_flt
       
       if(profstructflag .eq. 0) then
         write(log_lun,*)'   profile cluster info set to missing'
	 return
       endif
        
!--- find available lun

       call GPM_lun(rlun)
       
!--- read in raining clusters

       fname = trim(dir_anc) // trim(fspec(1))       
       open(rlun,file=fname,access='stream',status='old',iostat=ios)
       if(ios .ne. 0) call GPM_reprt_err(22,1,trim(fname))   !err open non-rain prf clusters
       do itemp = 1,ntemps
         do ispec = 1,nspecies
	   do iclust = 1,nclusts
	     read(rlun,iostat=ios)
     >            ProfileClusters(ispec,itemp,:,iclust)
             if(ios .ne. 0) call GPM_reprt_err(23,1,trim(fname))
	   enddo
	 enddo
       enddo       

!--- read in non-raining clusters

       fname = trim(dir_anc) // trim(fspec(2))       
       open(rlun,file=fname,access='stream',status='old',iostat=ios)
       if(ios .ne. 0) call GPM_reprt_err(22,1,trim(fname))   !err open non-rain prf clusters
       do itemp = 1,ntemps
         do ispec = 1,2                     !only has CWC, IWC
	   do iclust = 1,nclusts
	     read(rlun,iostat=ios)
     >            ProfileClustersNR(ispec,itemp,:,iclust)
             if(ios .ne. 0) call GPM_reprt_err(23,1,trim(fname))
	   enddo
	 enddo
       enddo

       return
      end subroutine GPM_read_profile_clusters	 	 

!------------------------------------------------------------------------------------

      subroutine GPM_read_rainsnowtwb
       implicit none
       character(len=128) :: fspec= 'rainsnowtwb.txt'
       character(len=256) :: infile
       character(len=20) :: headstr
       integer :: rlun,ios,i
                                      
!--- find available lun

       call GPM_lun(rlun)

!--- open rainsnowwtb.txt file

       infile = trim(dir_anc) // trim(fspec) 

       open(rlun,file=infile,form='formatted',status='old',readonly,
     >      iostat=ios)
       if(ios .ne. 0) call GPM_reprt_err(24,1,trim(infile))       

!--- read in data from file

       read(rlun,'(a20)',iostat=ios) headstr               !reads header line of ascii file
c       write(log_lun,*)' headstr= ',headstr
       do i = 1,131
         read(rlun,'(F4.1,2F8.2)',iostat=ios) rainsnowtwb(1:3,i)
         if(ios .ne. 0) call GPM_reprt_err(25,1,trim(infile))
	 
	 rainsnowtwb(1,i) = rainsnowtwb(1,i) + 273.15      !change table values to Kelvin
	 
c	 write(log_lun,*) rainsnowtwb(1:numrows,i)
       enddo

       return
       end subroutine GPM_read_rainsnowtwb

!------------------------------------------------------------------------------------

      end module GPM_read_procedures
