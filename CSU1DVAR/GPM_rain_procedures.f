      module GPM_rain_procedures
!       apart from use statement below, only a few modifications.
!       - custom loop start and stop points, commonly defined. (DD, 3/21/16)
!       version of GPROF is 2014_V2, copied over 3/21/16
!       - logical array determines whether gprof is run for given pixel, as
!       way to rerun pixels where sfccode has changed,[ also no deallocation of
!       database variables in here] (DD, 4/4/16)
!       updated to GPROF2017_V1 1/24/17
       use GPM_arraydef_D
       use GPM_util_procedures
       use define_intam2
       
       contains

!-----------------------------------------------------------------------

      subroutine GPM_rain(isfccode)
      implicit none

!---  retrieval routine for all surfaces

      integer :: ip, i, j, k, n, icnt(8), iexpand, totclustpix
      integer :: lin, pix, nch, itcwv, iskint, isfccode, s, t,ss,tt
      
      integer         :: T2m1, T2m2, tcwv1, tcwv2     
      real(kind=knd8) :: baywgt,invtoterr(15,mintcwv:maxtcwv),toterr
      real(kind=knd8) :: wgt,wgttot,wgtcond,wgttotcond,rms2,wgtscl,mvf
      real(kind=knd8) :: wgttest
      
      real(kind=knd8)    :: sumsfcprcp, sumppT2m=0
      real               :: cntsfcprcp, cntpixsfc, prof(nlyrs)
      integer(kind=knd8) :: numbersig
      real               :: avesfcprcp, numgoodchan
            	    
      real(kind=knd8) :: psfc, pcnv, prwp, pcwp, piwp, ptcwv, pprob
      real(kind=knd8) :: prwc(nlyrs),pcwc(nlyrs),piwc(nlyrs)
      real(kind=knd8) :: pswc(nlyrs), pgwc(nlyrs)          
      real(kind=knd8) :: pcsfc, pccnv
      
      real(kind=knd8) :: maxwgt
      real            :: missfactor
      integer         :: irainrate
      integer         :: maxexp 
      real            :: rrhist(0:2000), rrcdf
      
      logical         :: modsnowTbs = .true.
      integer         :: ilo, itbw
      real            :: snowprob
      
      sumsfcprcp = 0.0
      cntsfcprcp = 0.0
      icnt = 0

!---  don't use this channel in Tb matching RMS
c      write(log_lun,*)' setting 37V to missing'
c      chan_avail(7) = 0

!--- calculate satellite channel weights
  
      numgoodchan = 0.
      do j = mintcwv, maxtcwv
        do i = 1,15		     !loop over all  15 channels
          if(chan_avail(i)) then
              if(chansensfull(i,isfccode,j).ne.999.99 .and. j.eq.1) then
	      	numgoodchan = numgoodchan + 1.
	      endif		  
              toterr = chansensfull(i,isfccode,j)**2+chan_errs(i)**2
	      invtoterr(i,j) = 1./toterr			  !just for multiplying
          else						          ! in bayesian
              invtoterr(i,j) = 0.0                 ! (channel,tcwv)
          endif
        enddo
      enddo

      write(log_lun,'(15F8.2)') chan_errs(:)
      write(log_lun,'(15F8.2)') chansens(:,isfccode)
      write(log_lun,*)'   num of good channels = ',nint(numgoodchan)

!--- start the loop and rainfall retrieval over all pixels

!      do lin = 1, nscans       
!       do pix = 1, npixs
      do lin = scstart, scend
       do pix = pxstart, pxend
       
        maxwgt    = 0.0
	rrcdf     = 0.0
	rrhist    = 0

!---   QC the pixels + must be the correct sfccode

        if(.not.torun_gprof(pix,lin)) cycle  ! NEW for INTAM2 version
        if(pixel_status(pix,lin) .ne. 0) cycle        !don't do pixel if QC'd out	
	if(isfccode .ne. sfccode(pix,lin)) cycle      !retrieve only when matching dbase sfccode
          icnt(1) = icnt(1) + 1
	  
	if(histpp(isfccode) .eq. 0) then              !no profiles with this sfccode
	     sfcprcp(pix,lin) = miss_flt
	     cnvprcp(pix,lin) = miss_flt
	     frzprcp(pix,lin) = miss_flt
	     pixel_status(pix,lin) = 3                !surface code/ histogram mismatch 
	     cycle
	endif	  
        icnt(2) = icnt(2) + 1

	if(pptcwv(pix,lin).lt.0 .or. ppskint(pix,lin).le.0 .or.
     >       sfccode(pix,lin).le.0) then                  !if tcwv,skint,sfccode pixel missing
  	      icnt(3)  = icnt(3) + 1                         !so...
              sfcprcp(pix,lin) = miss_flt                 !set rain,frzprcp,cnvprcp to missing
	      cnvprcp(pix,lin) = miss_flt
	      frzprcp(pix,lin) = miss_flt
	      pixel_status(pix,lin) = 4                   !missing preprocess tcwv,skint,or sfccode
	      cycle                                       !go to next pixel                    
        endif
	
!---   limit observed tcwv, T2m to the max values

        if(pptcwv(pix,lin)  .gt. maxtcwv)  pptcwv(pix,lin) = maxtcwv
	if(ppT2m(pix,lin)   .gt. maxT2m)   ppT2m(pix,lin)  = maxT2m

!---   set limits for MRMS snow cover databases  (8-11)

        if(modsnowTbs) then
	    if(isfccode.eq.8  .and. ppT2m(pix,lin).lt.256) 
     >                                            ppT2m(pix,lin)=256
            if(isfccode.eq.9  .and. ppT2m(pix,lin).lt.254) 
     >                                            ppT2m(pix,lin)=254
	    if(isfccode.eq.10 .and. ppT2m(pix,lin).lt.252) 
     >                                            ppT2m(pix,lin)=252
	    if(isfccode.eq.11 .and. ppT2m(pix,lin).lt.250) 
     >                                            ppT2m(pix,lin)=250
	endif
	
!---   Zero out all the bayesian averages sums

        wgttot = 0.0
	wgttotcond = 0.0				    
	psfc=0.0; pcnv=0.0; pcsfc=0.0; pccnv=0.0
	prwp=0.0; pcwp=0.0; piwp=0.0; ptcwv=0.0; pprob=0	              	    
	prwc=0.0; pcwc=0.0; piwc=0.0; pswc=0.0; pgwc=0.0          !set summation values = zero
	cntpixsfc=0
  	  
!---   set the T2m/tcwv start/end to include bins around pixel's value
	  	
	T2m1 = nint(ppT2m(pix,lin) - 1)          ! use +/- 1 T2m
	T2m2 = nint(ppT2m(pix,lin) + 1)          
	tcwv1 = nint(pptcwv(pix,lin))            !NO expansion in TCWV - done in Dbase
	tcwv2 = nint(pptcwv(pix,lin))   

	if(T2m1 .lt. minT2m) T2m1 = minT2m
 	if(T2m2 .gt. maxT2m) T2m2 = maxT2m           
	if(tcwv1 .lt. mintcwv) tcwv1 = mintcwv
	if(tcwv2 .gt. maxtcwv) tcwv2 = maxtcwv
	  
!---   start bayesian averging over all the clustered profiles in the bin space

	do s = T2m1, T2m2            !loop over defined T2m bins
	  do t = tcwv1, tcwv2	     !loop over defined tcwv bins
	    
	    do ip = 1, db_nprofiles(s,t)         !all profiles in each bin 
              rms2 = 0.0
              if(EIAoffset) then                                       !eiaoffset flag set in main
		 do nch = 1,maxchans                                   !this is the main RMS calc
	           if(db_tbs(idx(s,t,ip),nch) .gt. 0.0) then           !If dbase chan tb good
                         rms2= rms2+ (invtoterr(nch,t) *               !  then include in RMS
     >	                    (tbs(pix,lin,nch) +
     >                      (tbdelta_eia(pix,lin,nch) *                  
     >                      (db_tbdelta(idx(s,t,ip),nch)/100.)) -        !tbdelta_eia calculated
     >                      (db_tbs(idx(s,t,ip),nch))/100.)**2)          !earlier
			
                   endif
		 enddo		      
              else		      
                 do nch = 1,maxchans		                       !main RMS loop if eiaoffset      
	           if(db_tbs(idx(s,t,ip),nch) .gt. 0.0) then            !is off. If Database chan Tb          
		      rms2 = rms2 + (invtoterr(nch,t)*(tbs(pix,lin,nch)-   !is good, include chan in RMS 
     >                              (db_tbs(idx(s,t,ip),nch)/100.))**2)
                   endif
		 enddo
	      endif

!---         start the summing and diagnostics calculation
		
	      baywgt = exp(-rms2 / 2.0) 			 !bayesian pixel weight
	     
	      if(baywgt .gt. maxwgt) then                        !most likely precip - 
                  maxwgt = baywgt                                ! is sfcprcp where the
		  mlprcp(pix,lin) = db_sfcprcp(idx(s,t,ip))      ! weight is the greatest
	      endif 	
		
	      if(db_sfcprcp(idx(s,t,ip)) .ge. 0.0) cntpixsfc=cntpixsfc+1
	     
	      wgt     = baywgt * db_occur(idx(s,t,ip))
	      wgtcond = baywgt * db_occur(idx(s,t,ip)) * 
     >                           db_raincount(idx(s,t,ip))
     
              psfc  = psfc  + (wgt     * db_sfcprcp(idx(s,t,ip)))
	      pcsfc = pcsfc + (wgtcond * db_sfcprcp(idx(s,t,ip)))
	      
              pcnv  = pcnv +  (wgt     * db_cnvprcp(idx(s,t,ip)))
	      pccnv = pccnv + (wgtcond * db_cnvprcp(idx(s,t,ip)))
	      
	      pprob = pprob + (wgt*((db_raincount(idx(s,t,ip))*1.0)/   !rain probability
     >       			    (db_occur(idx(s,t,ip))*1.0)))

	      irainrate = db_sfcprcp(idx(s,t,ip)) * 100 	 !create a rainrate histogram
	      if(irainrate .gt. 2000) irainrate = 2000  	 !from 0.0 to 20mm/hr
	      rrhist(irainrate) = rrhist(irainrate) + wgt	 !scaled by 100 to 0-2000

              prwp = prwp  + (wgt * db_rwp(idx(s,t,ip)))	 !path variables
              pcwp = pcwp  + (wgt * db_cwp(idx(s,t,ip)))		
	      piwp = piwp  + (wgt * db_iwp(idx(s,t,ip)))		

              if(profstructflag .eq. 1) then				   
             	  prwc(:)= prwc(:) + (wgt * db_rwc(idx(s,t,ip),:))
             	  pcwc(:)= pcwc(:) + (wgt * db_cwc(idx(s,t,ip),:))
                  piwc(:)= piwc(:) + (wgt * db_iwc(idx(s,t,ip),:))
             	  pswc(:)= pswc(:) + (wgt * db_swc(idx(s,t,ip),:))
             	  pgwc(:)= pgwc(:) + (wgt * db_gwc(idx(s,t,ip),:))
	      endif		 
           
	      wgttotcond = wgttotcond + wgtcond
	      wgttot = wgttot + wgt
	         
	    enddo   ! ip profile loop
	      
	  enddo   ! t = tcwv loop
	enddo   ! s = T2m loop
	
        wgttest = wgttot
	if(probcut) wgttest = wgttotcond
     
        if(wgttest .eq. 0) then              !no solution for pixel	       
	    sfcprcp(pix,lin) = miss_flt  
	    cnvprcp(pix,lin) = miss_flt
	    frzprcp(pix,lin) = miss_flt
	     
	    rwp(pix,lin)     = miss_flt   
	    cwp(pix,lin)     = miss_flt   
	    iwp(pix,lin)     = miss_flt
	    pop(pix,lin)     = miss_flt
	     
	    mlprcp(pix,lin)	 = miss_flt	 
	    prcp1stT(pix,lin)   = miss_flt
	    prcp2ndT(pix,lin)   = miss_flt
            
            if(profstructflag .eq. 1) then
          	 rwc(pix,lin,:) = miss_flt
          	 cwc(pix,lin,:) = miss_flt
          	 iwc(pix,lin,:) = miss_flt
          	 swc(pix,lin,:) = miss_flt
          	 gwc(pix,lin,:) = miss_flt
	    endif
	    pixel_status(pix,lin) = 5  	  !no bayesian solution for the pixel 
 	    icnt(4) = icnt(4) + 1
               
	else         	      	  

!---       surface and convective precipitation

            if(probcut) then
	        sfcprcp(pix,lin) = pcsfc  / wgttotcond        !apply conditional weights
	        cnvprcp(pix,lin) = pccnv  / wgttotcond
	    else
	        sfcprcp(pix,lin) = psfc  / wgttot
	        cnvprcp(pix,lin) = pcnv  / wgttot
	    endif

!---       Precip phase descrimination by Twb

            ilo = 2                              !array index to use land values from twb table
            if(isfccode .eq. 1) ilo = 3          !array index to use ocean values from twb table
	    snowprob = 0.0                       !assign frzprcp - from twb table, amount 
	    do i = 1,131                                      ! of sfcprcp that is frozen
              if(Twb(pix,lin) .le. rainsnowtwb(1,i)) then
	           snowprob = 1.0 - (rainsnowtwb(ilo,i)/100.)
	           itbw = i
		   exit
	      endif
            enddo
	    if(Twb(pix,lin) .ge. rainsnowtwb(1,131)) then
	        snowprob = 0.0
		itbw = 131
	    endif	     
	    frzprcp(pix,lin) = sfcprcp(pix,lin) * snowprob     !assign frzprcp value

!---       probability and threshold cutoffs if turned on
            
	    pop(pix,lin)     = pprob / wgttot                 !ave probability of precip
	       
	    if(probcut) then
	       sumppT2m = 0.0
	       t = tcwv1
	       do i = T2m1,T2m2            !ave the 3 T2ms
	         sumppT2m = sumppT2m + i
	       enddo
	       s = nint(sumppT2m / 3.)	    
	       if(pop(pix,lin)*100. .le. probcutoff(s,t,isfccode)) then 
		  sfcprcp(pix,lin) = 0.0                              !set precip to 0.0
  		  frzprcp(pix,lin) = 0.0
		  cnvprcp(pix,lin) = 0.0
	       else	           		      
                  if(missvolfrac(s,t,isfccode) .ne. -1.) then
		     missfactor = 1.0 / (1.0-missvolfrac(s,t,isfccode))   !missing mult factor
   		  else
   		     missfactor = 1.0
		  endif
		  		
		  if(probquality(s,t,isfccode) .eq. 1) then           !prob quality poor
		     quality_flag(pix,lin)=2 
	          endif

	          sfcprcp(pix,lin) = sfcprcp(pix,lin) *  missfactor    !add back rain
		  cnvprcp(pix,lin) = cnvprcp(pix,lin) *  missfactor    !removed due to
		  frzprcp(pix,lin) = frzprcp(pix,lin) *  missfactor    !probability cutoff
               endif
	    endif

!---       Precip distribution stats
	     
	    prcp1stT(pix,lin)   = miss_flt	  !find rainrate at 1st and 2nd Tertial of
	    prcp2ndT(pix,lin)   = miss_flt	  !the rainrate PDF 
	    if(cntpixsfc .gt. 0 .and. wgttot.gt.0) then 
	        do i = 0,2000  	       
	 	  rrcdf = rrcdf + (rrhist(i) / wgttot)
	          if(prcp1stT(pix,lin) .eq. miss_flt) then
	              if(rrcdf .ge. 0.3333) prcp1stT(pix,lin)= i/100.
	          endif
	          if(rrcdf .ge. 0.6666) then
	              prcp2ndT(pix,lin) = i / 100.
	              exit
	 	  endif		   
	        enddo  	       
	    endif 

!---       Path Variables	    
	    rwp(pix,lin)  = prwp  / wgttot
	    cwp(pix,lin)  = pcwp  / wgttot
	    iwp(pix,lin)  = piwp  / wgttot

!---       Profile selection if profiles turned on (profstrucflag)

	    if(profstructflag .eq. 1) then          
	 	rwc(pix,lin,:) = prwc(:) / wgttot / 1000.	 !g/m**3		
	        cwc(pix,lin,:) = pcwc(:) / wgttot / 1000.       		  
	        iwc(pix,lin,:) = piwc(:) / wgttot / 1000.       					  
	        swc(pix,lin,:) = pswc(:) / wgttot / 1000.       					  
	        gwc(pix,lin,:) = pgwc(:) / wgttot / 1000.       !missing in GPM V5      			
		call GPM_select_profile_clusters(pix,lin)       !gets profile cluster info
	    endif	     
	endif
	  
	if(sfcprcp(pix,lin) .ne. MISS_FLT) then
	     sumsfcprcp = sumsfcprcp + sfcprcp(pix,lin)	     	     
	     cntsfcprcp = cntsfcprcp + 1
        endif
	
       enddo   !NPIXS pixels		
      enddo   !LIN  scans

      if(cntsfcprcp .gt. 0) then
           avesfcprcp = sumsfcprcp / cntsfcprcp
      endif
 
      deallocate (idx,db_nprofiles, db_occur, db_raincount, db_tbs)
      deallocate (db_tbdelta,db_sfcprcp,db_cnvprcp)
      deallocate (db_rwp,db_cwp,db_iwp,db_tcwv, db_T2m)
      if(profstructflag) deallocate(db_rwc,db_cwc,db_iwc,db_swc,db_gwc)
      if(EIAoffset) deallocate (tbdelta_EIA)
	
      write(log_lun,'(a,F8.3,2x,i8)')'   average sfcprcp, num pixs : ',
     >                                 avesfcprcp,
     >                    nint(cntsfcprcp)
      write(log_lun,*)'   num passing sfccode       : ', icnt(1)
      write(log_lun,*)'   num passing histpp        : ', icnt(2)
      write(log_lun,*)'   num bad tcwv/T2m params   : ', icnt(3)
      write(log_lun,*)'   num no bayesian solution  : ', icnt(4)

      return
      end subroutine GPM_rain
      
!------------------------------------------------------------------      

      subroutine GPM_select_profile_clusters(pix,lin)
       implicit none

       real       :: rms, min_diff, sum, zsum
       integer    :: ispec, pix, lin, k, m,i
       integer    :: tempindex       
       real       :: nprof(nlyrs)
       integer    :: cluster_match
       
!--- check T2m for this pixel, if missing, no cluster retrieval
      
       if(ppT2m(pix,lin) .eq. miss_int) return       
       tempindex = nint((ppT2m(pix,lin)-268.) / 3.)       !set T2m to 1-12 temp index
       if(tempindex .le. 0)  tempindex = 1
       if(tempindex .gt. 12) tempindex = 12       
       prfT2mindex(pix,lin) = tempindex                                !1-12 temp index

!--- match hydrometer profiles
 
       do ispec = 1, nspecies
         cluster_match = 0; min_diff = 1000.0; sum = 0.0
         prfscale(pix,lin,ispec) = 0.0

         if(sfcprcp(pix,lin) .le. 0.01) then
            if(ispec.eq.1 .or. ispec.ge.4)then         ! if no sfcprcp, set RWC,
	         rwc(pix,lin,:) = 0.0                  ! SWC,and GWC profiles = 0
	         swc(pix,lin,:) = 0.0
		 gwc(pix,lin,:) = 0.0
	         cycle
            endif
         endif
        
         do k = 1,nlyrs
           select case(ispec)
            case(1)
              sum = sum + rwc(pix,lin,k)       ! sum up values in each profile 
            case(2)
              sum = sum + cwc(pix,lin,k)       
            case(3)
              sum = sum + iwc(pix,lin,k)       
            case(4)
              sum = sum + swc(pix,lin,k)       
            case(5)
              sum = sum + gwc(pix,lin,k)       
           end select
         enddo 
	
         if(sum .eq. 0) cycle                         !profile is all zeros, skip the rest
	 
	 prfscale(pix,lin,ispec) = sum
	 select case(ispec)                           !normalizes the clusters to 0-1
	   case(1)                                    !same as in cluster database files
	     nprof(:) = rwc(pix,lin,:) / sum
	   case(2)
	     nprof(:) = cwc(pix,lin,:) / sum
	   case(3)
	     nprof(:) = iwc(pix,lin,:) / sum
	   case(4)
	     nprof(:) = swc(pix,lin,:) / sum
	   case(5)
	     nprof(:) = gwc(pix,lin,:) / sum
	 end select
	   
	 do m = 1, nclusts
	   zsum = 0.0	     
           if(sfcprcp(pix,lin) .le. 0.01) then
              do k = 1,nlyrs
	        if(m .eq. 2 .or. m .eq. 3) then               !only cwc,iwc in NR(1:2)		   
                   zsum=zsum + profileclustersNR(ispec-1,tempindex,k,m)
		else
		   zsum = 0.0
		endif
	      enddo
           else
              do k = 1,nlyrs
                zsum = zsum + profileclusters(ispec,tempindex,k,m)
              enddo
	   endif
	   						
	   if(zsum .eq. 0.0) cycle		       ! no sum in this cluster
           
	   if(sfcprcp(pix,lin) .le. 0.01) then
	      rms = 0.0         
	      do k = 1,nlyrs
                rms =rms +(profileclustersNR(ispec-1,tempindex,k,m) - 
     >  		nprof(k))**2
              enddo
	   else
	      rms = 0.0          
	      do k = 1,nlyrs
                rms =rms +(profileclusters(ispec,tempindex,k,m) - 
     >  		nprof(k))**2
              enddo
           endif
	  
	   min_diff = min(min_diff, rms)		! identify cluster profile
           if(min_diff .eq. rms) cluster_match = m	! closest to observed
			      
         enddo   !cluster
	 	
	 if(sfcprcp(pix,lin) .le. 0.01) then	                    
	       prfnum(pix,lin,ispec) = cluster_match + 40      !cluster pointer for NR
	 else
	       prfnum(pix,lin,ispec) = cluster_match          !cluster pointer for Raining
	 endif

c	 write(log_lun,*)'best cluster match = ',prfnum(pix,lin,ispec),
c     >        prfscale(pix,lin,ispec), prfT2mindex(pix,lin)
      
       enddo    !ispec
       
      return
      end subroutine GPM_select_profile_clusters

!--------------------------------------------------------------------------------    

      end module GPM_rain_procedures
