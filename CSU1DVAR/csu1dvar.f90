      module csu1dvar

      
      USE define_csu1dvar
      USE define_intam2
      USE CRTM 
      USE csu1dvar_subr
      USE radtran
      USE GPM_arraydef_D

      implicit  none

!       non-raining + raining 1DVAR retrieval algorithm for 
!        microwave imagers
!      
!       Last updated 
!       D Duncan, Chalmers, Dec 2017
!       
!----------------------------------------------------------------------

      contains

      subroutine run_csu1dvar
       
! oe variables taken from pixel and scan numbers passed in the call  
      real         ::  oe_lat, oe_lon, oe_sst
      real         ::  worst
      
      integer   :: z,lodex,ladex!index for lon,lat
      logical     :: flag,mk(4)
      integer :: tdex,sdex 
      real :: init_tpw 

!     Arrays for Apriori Parameters and associated errors
      real  :: xa(nvar), sa(nvar,nvar), noretrieval 
      real  :: xmin(nretvar), xmax(nretvar)
      
!     Retrieval/Error diagnostics
      real         Amatrix(nvar), sigma(nvar)

!     Counters, etc.
      integer      nit, a, b, c, d, e, mc, i
      integer :: tot_iter=0,notrun_count=0 ! counters for oe
      real :: dmr(nz), tpwerr_lyr(nz,3), tpw_posterr
      real :: extra_ssd, extra_mr(nz), extra_eo(nz,6),extra_td,extra_sd
      real :: apecmr(nz)
      real :: offs(maxchans),fssd
      real :: rrr, chsq
      integer :: swic, n1,n2,n3, fl1,fl2,fl3
      real :: xr1(nvar),xr2(nvar),xr3(nvar), sig1(nvar),sig2(nvar),sig3(nvar)
      real :: A1(nvar),A2(nvar),A3(nvar), ch1,ch2,ch3
      real :: rwpdif,frp,rwp_pass,lwp_pass
      real :: saver2,saver3,saverwp2,saverwp3,saveiwp2,saveiwp3,savetp2,savetp3

      noretrieval = 0.00001
      mk = (/1,1,0,1/) !don't sum over 183+/-3 for high freq scatt calc

      allocate(oe_tbs(nch),test_out(nch))
      allocate(test_out2(nch),test_out3(nch))

      sa(:,:) = 0.0  ! initialize Sa matrix - set all to zero.
      
      ! read in lwmax and drizzle sa/eofs
      call read_lwmax
      call read_drizz

      mc = maxchans  !if(trim(sensor).eq.'AMSR2') mc = 12 ! change to 14 if using 7GHz too
        
        !initialize oe_output fields, tb diff to -997 (not run!)
            oe_output(:,:,:) = -997
            Tb_diff(:,:,:) = -997
            Amat_out(:,:,:) = -997
            save_iter(:,:) = 0 
            save_sw(:,:) = 9 

        ! save_sw codes: 
        !0 = non-raining converged, 
        !1 = stratiform rain, 
        !2 = convective rain, 
        !3 = regression-based RR, 
        !4 = too much scattering, 
        !5 = rain retrs failed, revert to NR solution (not used currently), 
        !6 = none of retrievals converged, 
        !9 = not run (land or not in chosen domain)


!!!!!!!! ---- BEGIN LOOP OVER PIX/SCANS ---- !!!!!!!

            do oelin = scstart, scend 
              !write(*,*),'--   scan = ',oelin
              do oepix = pxstart, pxend
               !print*,'---   pix = ',oepix

               ! define lat,sst,tbs for this pixel before calling OE
               oe_lat = lat(oepix,oelin) 
               dlat = oe_lat ! save to global var to pass more easily
                !print*,'dlat',dlat,abs(dlat),latdsd
               oe_lon = lon(oepix,oelin)

               ! new for goe9, lat/lon targeting (set define_intam2.f)
               if(oe_lat>maxla.or.oe_lat<minla.or.oe_lon>maxlo.or.oe_lon<minlo) then
                 oe_output(oepix,oelin,:) = -997 ! sea ice marker, ehh
                 notrun_count= notrun_count+1
                 cycle
               endif

               if(oe_lon.lt.0) oe_lon=lon(oepix,oelin)+360.0
               oe_sst = sst(oepix,oelin)
              ! fix SST, force retrieval for sea ice pixels
               if(oe_sst.lt.273.and.lo_flag(oepix,oelin,1).le.ocut) then !6GHz FOV
                 if(runseaice) then
                   oe_sst=272.5
                 else
                   oe_output(oepix,oelin,:) = -997 ! sea ice marker!
                   notrun_count= notrun_count+1
                   cycle 
                 endif
               endif ! assumes that sfc_type=1 is only other type!

               ! check for 100% ocean, sun glint, SST
               if(lo_flag(oepix,oelin,1)>ocut .or.&
                  oe_sst<271 .or. eia(oepix,oelin,1)<0 .or. &
                  (eia(oepix,oelin,14)<0.and.satcode=='GMI') .or. &
                  (minval(eia(oepix,oelin,:))<0.and.satcode=='TMI').or.&
              (oelin<=(nscans-2).and.(lo_flag(oepix,oelin+1,1)>ocut).or.&
              (oelin<=(nscans-1).and.&
               (lo_flag(oepix,oelin+1,1)>ocut)).or.&
                (oelin>=2 .and.(lo_flag(oepix,oelin-1,1)>ocut))   .or. &
                (oepix<=(npix-1).and.lo_flag(oepix+1,oelin,1)>ocut).or. &
                (oepix>=2 .and. lo_flag(oepix-1,oelin,1)>ocut) .or. &
                (oepix>=2.and.oelin>=2 .and. & ! +/- 1 diagonal
                 lo_flag(oepix-1,oelin-1,1)>ocut)     .or. &
                (oepix<=(npix-1).and.oelin>=2 &
                .and.lo_flag(oepix+1,oelin-1,1)>ocut) .or. &
                (oepix>=2.and.oelin<=(nscans-1) &
                .and.lo_flag(oepix-1,oelin+1,1)>ocut) .or. &
                (oepix<=(npix-1).and.oelin<=(nscans-1) &
                .and.lo_flag(oepix+1,oelin+1,1)>ocut) )) then
                  oe_output(oepix,oelin,:) = -998 ! land/bad pixel marker!
                  notrun_count= notrun_count+1
                  cycle ! high quality oceanic Tbs only!
               endif
               
               ssdex = floor((oe_sst-271.0)/1.0)+1
               if(ssdex.lt.1) ssdex = 1         ! first SST bin if <271K...
               if(ssdex.gt.nbins) ssdex = nbins ! max SST bin

               c = 1
               fssd = ((oe_sst-271.5)+1.0)-real(ssdex)
               do d = 1, mc !max # channels
                 if(avail(d).eq.1) then
                   if((fssd.lt.0.and.ssdex.eq.1) .or. fssd.eq.0 .or. &
                      (fssd.gt.0.and.ssdex.eq.nbins)) then
                     offs(c)=s_toffsets(ssdex,c)
                   endif
                   if(fssd.lt.0) offs(c)=s_toffsets(ssdex,c) + &
                   (s_toffsets(ssdex-1,c)-s_toffsets(ssdex,c))*abs(fssd)
                   if(fssd.gt.0) offs(c)=s_toffsets(ssdex,c) + &
                   (s_toffsets(ssdex+1,c)-s_toffsets(ssdex,c))*abs(fssd)

                   oe_tbs(c) = oetbs(oepix,oelin,d) + offs(c)  
                   c = c+1
                 endif
               enddo
               worst = minval(oe_tbs)
               if(worst .eq. miss_flt) then
                 oe_output(oepix,oelin,:) = -998
                 cycle
               endif
        
        ! Grid starts now at 0E, 90N !!
        lodex = nint(oe_lon/(losize/gm))+1 ! should be exact!!
        ladex = nint(abs(89.4628-oe_lat)/(lasize/gm))+1 ! more complicated...
        if(ladex.le.0 .or. ladex.ge.nlat+1) ladex=1 !just in case...

            ! assign pixel values from analysis (in pp file)
      pwindmag = ppwind(oepix,oelin)
      pwinddir = save_wdir(oepix,oelin)
      pecslp   = save_slp(oepix,oelin)
      pplev(1:nz) = lpress(1:nz)
      pplev(nz+1) = pecslp
      if(pecslp<lpress(nz)) pplev(nz)=pecslp-7.5
      if(pecslp<lpress(nz-1)) pplev(nz-1)=pecslp-15.0
      pecz(1:nz) = real(pph(oepix,oelin,:)) ! tops of layers!
      pecz(nz+1) = 0.0 ! set to sea sfc
      dp(nz) = pplev(nz+1)-pplev(nz)
      apect(:)  = save_tprof(oepix,oelin,:)
      apecmr(:) = pp_wvprof(oepix,oelin,:)
      do i = 1,nz ! get layer averages 
        if(pecz(i+1)>=pecz(i)) pecz(i+1)=pecz(i)-10.0 ! no neg/zero layers
        pressave(i) = (pplev(i+1)+pplev(i))/2.0
        dz(i) = pecz(i)-pecz(i+1)
        if(dz(i)<=0)dz(i)=0.01 ! make very thin
      enddo
      if(minval(apect).lt.190) then
        print*,'Tprof too low!',oepix,oelin,apect(:)
        stop
      endif

      if(sss(floor(oe_lon*2.0)+1,floor((90-oe_lat)*2.0)+1).lt.20) then
        if(psss.lt.20) psss= 35.0 ! default if clim grid point missing
        ! if psss already defined, just use that again (mis_val=-999.9)
      else
        psss = sss(floor(oe_lon*2.0)+1,floor((90-oe_lat)*2.0)+1) ! new,ver1.2
      endif

      save_slp(oepix,oelin)  = pplev(nz+1)
      !save_wdir(oepix,oelin) = pwinddir
      save_sal(oepix,oelin)  = psss
    
      ! Apriori error variances 
      ! Values of 0.00001 mean the parameters are fixed in the retrieval, ie not retrieved.
      ! Must not set any values to 0.0 since internal matrix 
      ! inversion will fail in the retrieval algorithm.

      ! reminder of variable order: EOF1-3 coeffs, WIND, LWP, IWP, SST.
      !sa(:,:) = 0.0
      ! only works if lwp is 5th state variable!
      do e = 1, npc ! take stddev values from EC analysis
        sa(e,e) = mr_sigmas(ssdex,e)
        sa(e,5) = eof_lw_cov(e,ssdex)
        sa(5,e) = eof_lw_cov(e,ssdex) ! OFFDIAGONAL LWP/EOF COVARIANCES
      enddo

        ! varmult for expanding variance if desired. set in definition file.
      sa(4,4) = (varmult*real(era_ws(lodex,ladex))*.01)**2 !10m WIND SIGMA SQUARED (m^2/s^2)
      sa(5,5) = 2.0**2 ! LOG10(LWP) SIGMA SQUARED
      sa(6,6) = noretrieval
      if(nretvar.eq.6) sa(7,7) = 0.7**2 ! -- SST SIGMA SQUARED (K^2)
      if(nretvar.eq.5) sa(7,7) = noretrieval ! -- SST SIGMA SQUARED (K^2)
      if(nretvar.eq.4) sa(7,7) = noretrieval ! -- SST SIGMA SQUARED (K^2)
      if(npc<3) then
        do e = npc+1, 3 
          sa(e,e) = noretrieval
        enddo
      endif

      d = 0
      do e = 1, nvar 
        if(sa(e,e).ne.noretrieval) then
          d=d+1
          xmin(d) = x_min(e)
          xmax(d) = x_max(e)
        endif
        if(e.eq.7.and.sa(e,e).ne.noretrieval) then
          xmin(d) = oe_sst - 2.0 ! set special bounds on SST retrieval
          xmax(d) = oe_sst + 2.0
        endif
      enddo
      if(d.ne.nretvar) then 
        print*,'nretvar doesnt match given variances!'
        stop 
      endif
        

!  Apriori state vector
! NOTE! -- Prior values for EOF coefficients can matter a lot!
!       -- Have to be non-zero, and it's better if 1 or 2 are
!       -- negative. Big values will bias the answer, and 
!       -- values too near zero cause more avg iterations needed!
       ! reminder of variable order: EOF1-3, WIND, LWP, IWP, SST.
      xa(1) =  0.02  ! EOF1 - something nominal, small, non-zero
      xa(2) = -0.01  ! EOF2 - something nominal, small, non-zero
      xa(3) =  0.01  ! EOF3 - something nominal, small, non-zero
      xa(4) = real(era_wm(lodex,ladex))*.01 !10 METER SFC WIND IN M/S
      xa(5) = -3.0 ! (0.001mm)         ![LOG10(LWP)] 
      xa(6) = 0.0001 ! whether retrieving or not, set here
      xa(7) = oe_sst        !SST IN K 
      mrmp(:) = mprof(ssdex,:) ! 
      peofs(:,:) = eofs(ssdex,:,:)

          peofs(:,:) = eofs(ssdex,:,:) ! initialize
          extra_eo(:,:) = 0.0
          ! smoothing of mrmp/eofs to lessen grid artifacts...
          extra_ssd = (oe_sst-271.5+1) - ssdex
          if(oe_sst .gt. 271.5 .and. oe_sst .lt. 303.5) then
            do i = 1, nz 
              if(extra_ssd.gt.0) then 
                extra_mr(i)=(mprof(ssdex+1,i)-mprof(ssdex,i))*extra_ssd
                do e = 1, npc
                  extra_eo(i,e) = (eofs(ssdex+1,i,e)-eofs(ssdex,i,e))*extra_ssd
                enddo
              endif
              if(extra_ssd.lt.0) then
                extra_mr(i)=(mprof(ssdex,i)-mprof(ssdex-1,i))*extra_ssd
                do e = 1, npc
                  extra_eo(i,e) = (eofs(ssdex,i,e)-eofs(ssdex-1,i,e))*extra_ssd
                enddo
              endif
              if(extra_ssd.eq.0) then
                extra_mr(i)   = 0.0
                extra_eo(i,:) = 0.0
              endif
            enddo
            mrmp(:) = mrmp(:) + extra_mr(:)
            peofs(:,:) = peofs(:,:) + extra_eo(:,:)
          endif
          if(oe_sst .gt.303.5) then ! special for extending high SST bin
            do i = 1, nz  
              extra_mr(i) = (mprof(ssdex,i)-mprof(ssdex-1,i))*extra_ssd
              do e = 1, npc
                extra_eo(i,e) = (eofs(ssdex,i,e)-eofs(ssdex-1,i,e))*extra_ssd
              enddo
            enddo
            mrmp(:) = mrmp(:) + extra_mr(:)
            peofs(:,:) = peofs(:,:) + extra_eo(:,:)
          endif

      if(useanalysis) then
        xa(4) = pwindmag !10 METER SFC WIND IN M/S
        mrmp(:) = apecmr(:)
      endif

      drdx = floor((oe_sst-271.0)/3.0)+1 !rain coeff index
      if(drdx>11)drdx=11
      if(drdx<1)drdx=1

! new part that determines LWP max value, DD 9/15/16
      sdex = floor((oe_sst-271.0)/2.0)+1 ! rounds down. 272.9 yields 1, etc.
      init_tpw=sum(1.0/(1000.0*9.81) * &
                     mrmp(:)*(dp(:))*100.0)
      tdex = floor((init_tpw)/3.0)+1 ! rounds down. 2.9 yields 1, etc.
      if(sdex<1)  sdex = 1         ! first bin if <1
      if(sdex>16) sdex = 16 ! max bin
      if(tdex<1)  tdex = 1         ! first bin if <1
      if(tdex>22) tdex = 22 ! max bin
        !print*,lwmax(sdex,tdex),sdex,tdex,oe_sst,init_tpw
        !xmax(5) = alog10(lwmax(sdex,tdex)*.001) ! no interpolation (for now)

      lwmaxy = lwmax(sdex,tdex)*.001 !into kg/m^2 
        ! interpolate the LWmax cutoff --new goe9w
        !! 16x22 (SST 271-272.9 up to 301+) (TPW 0-2.9 to 63+)
        extra_sd = ((oe_sst-271+1.0) - sdex*2)/2.0
        extra_td = ((init_tpw+1.5)   - tdex*3)/3.0
        if(sdex<16.and.extra_sd>0) lwmaxy = lwmaxy + &
                    extra_sd*(lwmax(sdex+1,tdex)-lwmax(sdex,tdex))*.001
        if(sdex>1.and.extra_sd<0)  lwmaxy = lwmaxy + &
                    extra_sd*(lwmax(sdex,tdex)-lwmax(sdex-1,tdex))*.001
        if(tdex<22.and.extra_td>0) lwmaxy = lwmaxy + &
                    extra_td*(lwmax(sdex,tdex+1)-lwmax(sdex,tdex))*.001
        if(tdex>1.and.extra_td<0)  lwmaxy = lwmaxy + &
                    extra_td*(lwmax(sdex,tdex)-lwmax(sdex,tdex-1))*.001
        !print*,'LWmax: ',lwmaxy,lwmax(sdex,tdex)*.001
        final_lwp = 0.0 !initialize path vars
        final_rwp = 0.0
        final_iwp = 0.0
        swic = 0 ! signifies NR retrieval
        rwpdif = 0.0 ! in kg/m^2
        nit = 0
        rrr = 0.0 ! init
        rwp_pass = 0.0
        lwp_pass = 0.0
        
        call opt_est(swic,nz,noretrieval, &
                      oe_tbs,xa,sa,xr,xmax,xmin, &
                      flag,Amatrix,chsq,sigma,nit)
        ch1    = chsq
        A1(:)  = Amatrix(:)
        xr1(:) = xr(:) ! save output to *1
        fl1    = flag
        n1     = nit
        sig1(:)= sigma(:)

      ! FIRST, assess whether NR retrieval successful or warm rain needs
      !  to be run
        if (flag .and. ch1<chout .and. ch1>0) then 
           oe_output(oepix,oelin,1) = final_tpw ! TPW
           oe_output(oepix,oelin,2) = 10**xr1(5) *1000 !CLWP [g]
           oe_output(oepix,oelin,3) = xr1(4) ! WIND
           oe_output(oepix,oelin,4) = rrr ! RR [mm/hr] 
           oe_output(oepix,oelin,5) = ch1 ! CHI SQUARED
           oe_output(oepix,oelin,6) = xr(7)  !SST [K]
           oe_output(oepix,oelin,8) = 0.0 !final_iwp
           oe_output(oepix,oelin,9) = 0.0 !final_rwp
           poster(oepix,oelin,1:nvar)  = sig1(:)  ! posterior errors - no matter retr type
           Amat_out(oepix,oelin,:)     = Amatrix(:)
           save_iter(oepix,oelin)      = n1
           Tb_diff(oepix,oelin,:) = oe_tbs(:)-test_out(:) 
           if (10**xr1(5) < lwmaxy) then !  save NR vars, cycle
              save_sw(oepix,oelin)        = swic
              last_oeout(:) = xr(:) !save for fg in next pixel (NR only)
              loeop = oepix !save last converged pix#
              loeos = oelin !save last converged scan#
              cycle 
           endif
           if ( 10**xr1(5) >= lwmaxy .and. ch1 < 1.0) then
              ! w/ no scattering, light drizzle... dont iterate further!
              rwpdif = (10**xr1(5) - lwmaxy)*1000 ! LW residual [g]
              frp = rwpdif-sqrt(rwpdif)
              if(frp<0) frp=0.0 ! possible if 0<rwpdif<1
              rrr = rcoef(drdx,1) + rcoef(drdx,2) * (frp*.001)
              save_sw(oepix,oelin) = 3 ! code for regression based RR, nonscattering
              if(rrr<0.0) rrr=0.0 !just in case
              oe_output(oepix,oelin,4) = rrr ! RR [mm/hr] 
              oe_output(oepix,oelin,2) = 10**xr1(5)*1000 - frp !CLWP
              oe_output(oepix,oelin,9) = frp !RWP
        !print*,rrr,save_sw(oepix,oelin),ch1,frp
              cycle
           endif
           if ( 10**xr1(5) >= lwmaxy .and. ch1 >= 1.0) then
              ! RWP, CWP priors for rain retrieval [in kg and log(kg)]
              rwp_pass = (10**xr1(5)-lwmaxy)-.001*sqrt(1000*(10**xr1(5)-lwmaxy)) !essentially RWP apriori
              ! convert rwp_pass to coeff on mean RWC profile:
              xa(2) = rwp_pass/sum(rwpm(drdx,:)) !voila: rwc1 prior
              lwp_pass = alog10(lwmaxy) ! CLWP left (in log)
              !lwp_pass = alog10(10**xr1(5)-rwp_pass) ! CLWP left (in log)
           ! NO cycle... this answer will be used if rain retrievals
           ! both fail!
           endif
        endif 

        ! however, no need to run rain retrievals if too much high freq
        ! scattering (defined as mean of 166V,H,190)
        if ((sum(oe_tbs(10:13),mask=mk)-sum(test_out(10:13),mask=mk))/3.0 < -8) then
           oe_output(oepix,oelin,:) = miss_flt
           poster(oepix,oelin,:)  = miss_flt  ! posterior errors
           save_iter(oepix,oelin) = n1
           save_sw(oepix,oelin)   = 4 ! code for too much scattering
           Amat_out(oepix,oelin,:)= Amatrix(:)
           cycle
        endif
        ! prep to run rain retrievals
        if (.not.flag .or. (flag.and.chsq>=chout).or. final_tpw>=75) then 
        !print*,'woulda used gprof' -- CAN TURN ON LATER AS SENSITIVITY
        !STUDY
           rwp_pass = .1 !? something nominal
           !rwp_pass = rwp(oepix,oelin) ! gprof rwp as prior?
           xa(2) = rwp_pass/sum(rwpm(drdx,:)) !voila: rwc1 prior
           lwp_pass = alog10(lwmaxy) ! as good a spot as any
           ! set output vars in case rain retrs fail too
           oe_output(oepix,oelin,:) = miss_flt
           poster(oepix,oelin,:)  = miss_flt  ! posterior errors
           save_iter(oepix,oelin) = -9
           save_sw(oepix,oelin)   = 6 ! code for getting to rain retr, not successful!
           Amat_out(oepix,oelin,:)= 0
        endif

        ! by here, should've 'cycled' if NR or drizzle or fwd model way
        ! off. time to set params for rain retrievals and run both.

        xa(1) = 0.01 !xr(1) !!take previous EOF1?
        !xa(2) = .03 !RWC1 EOF starting point --trying new thing (above)
        !xa(3) = 0.01 !-.01 !RWC2 EOF starting point --not used!
        xa(3) = rwp_pass   ! now used to pass through RWPdif (kinda cheating)
        xa(4) = pwindmag !xr(4) ! take previous wind?
        xa(5) = lwp_pass !xr1(5) 
        xa(6) = 0.09 !PIWC1 EOF starting point
        xa(7) = xr1(7)  !SST (still not retrieved)
        !  setting apriori covariance matrix:
        sa(:,:) = 0.0 ! re-initialize just to make sure
        sa(1,1) = (sqrt(mr_sigmas(ssdex,1))*0.4)**2 ! fraction of sigma value, WV EOF1
        sa(2,2) = 0.5**2 !madeup
        !sa(2,2) = 1.5*drizsa(1,1,drdx)  ! will interpolate later
        sa(3,3) = noretrieval !drizsa(2,2,drdx) 
        sa(4,4) = noretrieval !  NO wind retrieval in rain
        sa(5,5) = 0.01**2 ! CLWP tightly constrained
        sa(6,6) = 0.5**2 !madeup
        !sa(6,6) = drizsa(3,3,drdx) 
        sa(7,7) = noretrieval
        !sa(2,6) = drizsa(1,3,drdx) 
        !sa(6,2) = drizsa(3,1,drdx)
        !print*,'drizsa: ',drizsa

        ! set new min/max values for retr vars (index is by nretvar!!)
        xmin(2) = 0.0 !drzmi1(drdx)  ! rwc1
        xmax(2) = 4.0 !drzma1(drdx)
        xmin(nretvar-1) = x_min(5) ! lwp
        xmax(nretvar-1) = alog10(lwmaxy*1.5) !x_max(5)
        xmin(nretvar) = 0.0 !drzmi3(drdx) ! iwc1
        xmax(nretvar) = 4.0 !drzma3(drdx)
        
        ! wv profile from nr retrieval? ----NOPE
        mrmp(:) = apecmr(:) 

       !print*,'xa2: ',xa(2)
        swic = 1 !stratiform first
        call opt_est(swic,nz,noretrieval, &
                      oe_tbs,xa,sa,xr,xmax,xmin, &
                      flag,Amatrix,chsq,sigma,nit)
        ch2    = chsq
        A2(:)  = Amatrix(:)
        !print*,'rr2: ',saver,swic,sum(A2(:)),nit
        !print*,'++ ',xr(2),xr(6),final_iwp
        xr2(:) = xr(:) ! save output to *2
        fl2    = flag
        if (chsq>chout .or. final_tpw>75) fl2=0
        n2     = nit
        sig2(:)= sigma(:)
        saver2 = saver ! rr out of edd is global var...
        saverwp2 = final_rwp
        saveiwp2 = final_iwp
        savetp2 = final_tpw

        swic = 2 ! set to convective DSD -- try again
        call opt_est(swic,nz,noretrieval, &
                      oe_tbs,xa,sa,xr,xmax,xmin, &
                      flag,Amatrix,chsq,sigma,nit)
        ch3    = chsq
        A3(:)  = Amatrix(:)
        !print*,'rr3: ',saver,swic,sum(A3(:)),nit
        !print*,'++ ',xr(2),xr(6),final_iwp
        xr3(:) = xr(:) ! save output to *3
        fl3    = flag
        if (chsq>chout .or. final_tpw>75) fl3=0
        n3     = nit
        sig3(:)= sigma(:)
        saver3 = saver ! rr out of edd is global var...
        saverwp3 = final_rwp
        saveiwp3 = final_iwp
        savetp3 = final_tpw

        if (fl2+fl3==0) then   ! both failed... revert to NR solution (or lack thereof):
           if(fl1 .and. xr1(5)>=lwmaxy) then !(otherwise output given above stands)
              rwpdif = (10**xr1(5) - lwmaxy)*1000 ! LW residual [g]
              rrr = rcoef(drdx,1) + rcoef(drdx,2) * ((rwpdif-sqrt(rwpdif))*.001)
              !save_sw(oepix,oelin) = 5 ! code for reverting to NR rain
              if(rrr<0.0) rrr=0.0 !just in case
              oe_output(oepix,oelin,4) = rrr ! RR [mm/hr] 
              oe_output(oepix,oelin,2) = 10**xr1(5)*1000 -(rwpdif-sqrt(rwpdif))
              oe_output(oepix,oelin,9) = rwpdif-sqrt(rwpdif)
        !print*,'regn(bof)',oe_output(oepix,oelin,4),save_sw(oepix,oelin)
        !print*,'---',rwpdif-sqrt(rwpdif),10**xr2(5),ch2,ch3
           endif
           cycle
        endif
        ! strat or conv converged (with nonzero RR... if not revert to
        ! NR drizzle calc above)
        if ( (fl2 .or. fl3) .and. (saver2>0 .or. saver3>0) ) then
           if ((fl2.and.saver2>0.and..not.fl3) .or. &
               (fl3.and.fl2 .and. saver2>0.and.saver3==0) .or. &
               (fl3.and.fl2 .and.(saver2>0.and.saver3>0).and.ch2<=ch3) ) then
              oe_output(oepix,oelin,1) = savetp2 ! TPW
              oe_output(oepix,oelin,3) = xr(4) ! WIND
              oe_output(oepix,oelin,6) = xr(7)  !SST [K]
              oe_output(oepix,oelin,2) = 10**xr2(5) *1000 !CLWP [g]
              oe_output(oepix,oelin,4) = saver2/3.0 ! RR [mm/hr] 
              oe_output(oepix,oelin,5) = ch2 ! CHI SQUARED
              oe_output(oepix,oelin,8) = saveiwp2
              oe_output(oepix,oelin,9) = saverwp2
              poster(oepix,oelin,1:nvar)  = sig2(:)  ! posterior errors - no matter retr type
              Amat_out(oepix,oelin,:)     = A2(:)
              save_iter(oepix,oelin)      = n2
              Tb_diff(oepix,oelin,:)      = oe_tbs(:)-test_out2(:) 
              save_sw(oepix,oelin) = 1 ! code for reverting to NR rain
        !print*,'regn--',oe_output(oepix,oelin,4),save_sw(oepix,oelin)
        !print*,'---',saverwp2,1000*10**xr2(5),ch2,ch3
              cycle
           endif
           if ((fl3.and.saver3>0.and..not.fl2) .or. &
               (fl3.and.fl2 .and. saver3>0.and.saver2==0) .or. &
               (fl3.and.fl2 .and.(saver3>0.and.saver2>0).and.ch3<=ch2) ) then
              oe_output(oepix,oelin,1) = savetp3 ! TPW
              oe_output(oepix,oelin,3) = xr(4) ! WIND
              oe_output(oepix,oelin,6) = xr(7)  !SST [K]
              oe_output(oepix,oelin,2) = 10**xr3(5) *1000 !CLWP [g]
              oe_output(oepix,oelin,4) = saver3/3.0 ! RR [mm/hr] 
              oe_output(oepix,oelin,5) = ch3 ! CHI SQUARED
              oe_output(oepix,oelin,8) = saveiwp3
              oe_output(oepix,oelin,9) = saverwp3
              poster(oepix,oelin,1:nvar)  = sig3(:)  ! posterior errors - no matter retr type
              Amat_out(oepix,oelin,:)     = A3(:)
              save_iter(oepix,oelin)      = n3
              Tb_diff(oepix,oelin,:)      = oe_tbs(:)-test_out3(:) 
              save_sw(oepix,oelin) = 2 ! code for reverting to NR rain
        !print*,'regn--',oe_output(oepix,oelin,4),save_sw(oepix,oelin)
        !print*,'---',saverwp3,1000*10**xr3(5),ch2,ch3
              cycle
           endif
        endif

              enddo !pix loop
            enddo   !scan loop
            avg_iter = 9.9 !real(tot_iter) / real(run_count) ! iterations per retrieval -- wrong :)
            ! avg iterations kind of loses meaning  with 3-step
            ! retrievals going on, can revisit later

          call cleanup_crtm

        return
      end subroutine run_csu1dvar
!--------------------------------------------------- 


      subroutine noscatter(switch,X,TBOUT)

!     ROUTINES CALLED:
!     run_crtm:  calls crtm RT model

      integer :: switch
      INTEGER :: i,c !,j,k,m,n,z
      REAL :: TEMPAV(nz)
      REAL :: BTEMP,WIND,LiWaPa
      REAL :: IWP, ciwc(nz), clwc(nz),rwc(nz),piwc(nz),rawapa,PIWP
      REAL :: X(nvar),TBOUT(nch)
      real :: teebs(nch),tpw_lyr(nz)
      real :: teebc(nch),teebed(nch),teebr(nch),teed(nch), emu,edn
      real :: mixrat(nz),  clwced(nz)
      integer :: cldbindex, cldtindex
      real :: eco(3),reo(3) ! EOF coefficients
      real :: zeros(nz),rwpd(nz)

      zeros(:) = 0.0
      eco(:) = 0.0 ! initialize, then define eco(1:npc) below
      reo(:) = 0.0
      rwc(:) = 0.0
      piwc(:) = 0.0 !snow water, give or take
      rawapa = 0.0
      !final_rwp = 0.0
      WIND      = X(4)
      BTEMP     = X(7) ! ///JUST SST////
      IWP       = 0.0 !init
      LiWaPa    = 10**(X(5)) !LWP WAS IN LOG FORM; INVERT IT
      if(switch==0) then ! reminder of variable order: EOFC 1-3, WIND, LWP, IWP, SST.
        eco(1:npc)= X(1:npc)
      endif
      if(switch>0) then ! reminder of variable order: EOFC 1, DEOFC1-2, WIND, LWP, DEOFC3, SST.
        eco(1)    = X(1)
        eco(2:3)  = 0.0 ! 
        reo(1)    = X(2)
        reo(2)    = 0.0 !X(3)
        reo(3)    = X(6)
      endif

      cldbindex = nz-3 ! 16-3 = 13 -- 925-900mb
      if(btemp<280) cldbindex = nz-4 ! 16-4 = 12 -- 900-850mb
      cldtindex = nz-5 ! P-levels! look at def file for corresponding heights!
        ! change for ver3.0: cld layer 925-800mb cld layer (900-800 if sst<280K)
      clwc(:) = 0.0 ! initialize
      tempav(:) = apect(:) ! calculated above!
      ciwc = 0.0

      do i = 1, nz
        mixrat(i) = mrmp(i) + &
                  eco(1) * peofs(i,1) + & 
                  eco(2) * peofs(i,2) + & 
                  eco(3) * peofs(i,3) 
        if(mixrat(i).lt.0) mixrat(i)=0.0 ! cant have <zero mass!
      enddo

      if(switch>0) then
        do i = 3, nz !shifting CS data (missing last 2 p levs) down 2 levs
          !rwc(i) = rwpd(i) + reo(1)*rwpe(drdx,i-2,1)+reo(2)*rwpe(drdx,i-2,2)
          !rwc(i) = rwpd(i) + reo(1)*rwpm(drdx,i-2) !+reo(2)*rwpe(drdx,i-2,2)
          rwc(i) = reo(1)*rwpm(drdx,i-2) !+reo(2)*rwpe(drdx,i-2,2)
          if(rwc(i)<0.000001) rwc(i)=0.000001 ! round down to zero / neg -> 0
          !piwc(i) = piwpm(drdx,i-2) + reo(3)*piwpe(drdx,i-2,1)
          piwc(i) = piwpm(drdx,i-2) + reo(3)*piwpm(drdx,i-2)
          if(piwc(i)<.0001) piwc(i)=0.0 ! round down to zero / neg -> 0
          if(i<cldtindex) rwc(i)=0.0 !no RWC above designated cld top!
        enddo
        rawapa    = sum(rwc(:))
        PIWP      = sum(piwc(:))
      endif

      clwc(cldtindex:cldbindex) = liwapa/(cldbindex-cldtindex+1) !kg/m^2/layer

      final_lwp = liwapa*1000.0
      final_rwp = rawapa*1000.0
      final_iwp = piwp*1000.0
      final_tpw=sum(1.0/(1000.0*9.81) *mixrat(:)*(dp(:))*100.0)
      ! tpw_layer = 1/(rho_water*g) * MR * dP
      ! UNITS:  mm = (m^3/kg)(s^2/m)(g/kg)(kg/m*s^2)  -- need P in Pa!

!----- pass lyr temperature, wv mixing ratio,
!       LWP/IWP per layer, SST, wind speed

      if(switch==0) then
        call run_crtm(switch,tempav,mixrat,clwc,rwc,piwc,btemp,wind,teebs)
        tbout(:) = teebs(:)
      endif

      if(switch>0) then
        call run_crtm(switch,tempav,mixrat,clwc,zeros,zeros,btemp,wind,teebc)
        !print*,'RUNNING EDD CONTROL'
        emu = 1.0
        edn = 1.0 !shouldnt matter what mu/Dn are for control (0 rwc)
        call run_eddington(tempav,mixrat,clwc,zeros,zeros,btemp,emu,edn,wind,teebr)

        if(switch==1) then
          emu = 9 
          edn = 0.75 
          if(abs(dlat) < latdsd) then
            emu = 7
            edn = 0.83
          endif
          !print*,'RUNNING EDD RAIN STRAT'
          call run_eddington(tempav,mixrat,clwc,piwc,rwc, &
                    btemp,emu,edn,wind,teebed)
        endif
        if(switch==2) then
          emu = -1 
          edn = 1.8
          if(abs(dlat) < latdsd) then
            emu = 0.5
            edn = 1.60
          endif
          !print*,'RUNNING EDD RAIN CONV'
          call run_eddington(tempav,mixrat,clwc,piwc,rwc, &
                    btemp,emu,edn,wind,teebed)
        endif
        do c = 1,nch
          teed(c) = teebed(c) - teebr(c)
          tbout(c) = teebc(c) + teed(c)  !sim tbs = CRTM + (EddRain-EddCtrl)
        enddo
        !print*,'Tb c: ',teebc
        !print*,'Tb ed: ',teebr
        !print*,'Tb red: ',teebed
        !if(final_rwp>0) then
        ! print*,'reo/rr ',reo(1),saver,switch
        ! print*,'Tb ed diff: ',teed
        !endif
        !print*,'Tb rain: ',tbout
        !stop
        test_out(:) = tbout(:) !ie overwrites CRTM output if switch on
        if(switch==1) test_out2(:) = tbout(:)
        if(switch==2) test_out3(:) = tbout(:)
      endif

        !print*,'tbout: ',tbout(:)
        
      RETURN
      END subroutine noscatter     
!-----------------------------------------------------------

      SUBROUTINE opt_est(swtch,nz,noretrieval,&
                          y,xa2,sa2,x2,xmax,xmin,&
                          flag,Amatrix2,csq,sigma2,niter)

      implicit none

      integer :: swtch  ! IF swtch = ...
      integer :: nz   !Number of layers (nz) in atmosphere.  
      logical :: flag !flag indicates whether convergence was reached 

      ! Variables used in optimal estimation framework .
      integer :: n,niter ! defined above too!!
      integer    a,b,c,d,i,j,count,index,checkpol,end_flag
      real noretrieval
      real x(nretvar),x2(nvar),xnext(nretvar),xprime2(nvar)
      real xa(nretvar),xa2(nvar),xmax(nretvar),xmin(nretvar)
      real sa(nretvar,nretvar),sa2(nvar,nvar),sa_i(nretvar,nretvar)
      real K(nch,nretvar),K_t(nretvar,nch)
      real sx(nretvar,nretvar) ! posteriori error matrix
      real sx_i(nretvar,nretvar)
      real sigma2(nvar)
      real AP(nretvar,nretvar),Amatrix2(nvar)
      real F(nch),Fout(nch),Foutprime(nch) ! changed!
      real Fprime(nch),Fdbprime(nch),y(nch),dF(nch),dx(nretvar)
      real sum1(nretvar,1)         !xa-x
      real sum2(nch,1)           !y-F
      real sum3(nretvar,1)         !prod4+prod3
      real prod1(nretvar,nch)  !K_t*sy_i
      real prod2(nretvar,nretvar)!K_t*sy_i*K
      real prod3(nretvar)        !sa_i*(sum1)
      real prod4(nretvar)        !prod1*sum2
      real prod5(nretvar)        !sx*sum3
      real xdiff(nretvar,1)        !xnext-x
      real xdiff_t(1,nretvar)    !transposed
      real prod6(1,nretvar)      !xdiff_t*sx_i
      real prod7(1,1)            !prod6*xdiff 
      real sum4(nch,1)           !F-y
      real sum4_t(1,nch)       !(F-y) transpose
      real sum5(nretvar,1)         !x-xa
      real sum5_t(1,nretvar)     !(x-xa) transpose
      real prod8(1,nch)        !sum4_t*sy_i
      real prod9(1,1)              !prod8*sum2
      real prod10(1,nretvar)     !sum5_t*sa_i
      real prod11(1,1)           !prod10*sum5 
      real chisqtot(1,1)         !chi squared (prod9 + prod11) 
      real csq                 !cost function (apriori + measurement fit)
      real sum6(nretvar,nretvar) !(IM-AP)
      real prod12(nretvar,nch) !Sx*K_t
      real prod13(nretvar,nch) !Sx*K_t*Sy_i
      real prod14(nretvar)       !prod13*sum2 - contribution from obs
      real prod15(nretvar)       !sum6*sum1 - contribution from apriori
      real prod16(nretvar)       !prod14+prod16 = xnext-x

      real Fold(nch)           !old fwd modeled Tbs
      real diffF(nch,1),diffF_t(1,nch)  !F-F_old, transpose
      real Sdy(nch,nch),Sdy_i(nch,nch) ! (Rodgers Eq5.27) 
      real Ksa(nch,nretvar), Ksakt(nch,nch) ! K*Sa, K*Sa*K_t
      real Sdysum(nch,nch), Sdysum_i(nch,nch) ! (K*Sa*K_t+Sy), *_i
      real almost(nch,nch)! Sy*(Ksakt+Sy)_i (most of 5.27)
      real almost3(1,nch) ! diffF_t*Sdy_t
      real*8 disq(1,1)   !di squared (Eq5.33), used for conv test

      integer :: good(nretvar)
      real :: psdiff ! pix/scan distance from last convergence reached 

      real :: savejake(nch,nretvar) ! experimental!
      !integer :: moditer
      real :: osy(nch,nch),osy_i(nch,nch) ! opt_est-specific Sy,Sy_i
      real :: rosy(nch,nch) ! opt_est-specific Sy_rain
      real :: rref1,rref2,rcs(3),intf
      real*8 :: sg1(nch),sg2(nch),r1(nch,nch),r2(nch,nch)
      real :: rwadj(nch), rwad(nch,3,2)
        rwad(:,:,:) = 0.0 ! rwp analysis adjustments, rlo/mi/hi, stra/conv--10-166
        ! values from analysis (../rwp/), in K^2
        rwad(1:11,1,1) = (/0.0,0.0,0.0,0.0,0.0,0.0,0.1,0.0,0.1,0.0,0.0/)
        rwad(1:11,2,1) = (/0.0,0.0,0.0,0.1,0.0,0.1,0.6,0.2,0.1,0.2,0.1/)
        rwad(1:11,3,1) = (/0.0,0.0,0.1,0.3,0.1,0.5,2.2,0.6,0.1,0.5,0.4/)
        rwad(1:11,1,2) = (/0.0,0.0,0.0,0.0,0.0,0.0,0.2,0.0,0.0,0.0,0.0/)
        rwad(1:11,2,2) = (/0.0,0.0,0.1,0.3,0.1,0.1,1.2,0.1,0.1,0.0,0.0/)
        rwad(1:11,3,2) = (/0.1,0.2,0.4,1.3,0.3,0.3,3.5,0.4,0.1,0.1,0.1/)
        ! inits
        osy(:,:) = 0.0
        rosy(:,:) = 0.0

!-----Declare Apriori errors 
      sa(:,:) = 0.0
      b = 1
      do a=1,nvar ! less complex if not using covariances...
        if(sa2(a,a).ne.noretrieval) then 
          good(b) = a
          b=b+1
        endif
      enddo
      do a=1,nretvar
       do c=1,nretvar
         sa(a,c)=sa2(good(a),good(c))
       enddo
      enddo
        
!-----GENERATE FIRST GUESS FORWARD MODEL
! - - instead of using a priori as first guess, use last (good) retrieval
      psdiff = sqrt(real(loeop-oepix)**2 + real(loeos-oelin)**2)
      if(psdiff<=4.0 .AND. swtch==0) then !non-raining ONLY
        b=1
        do a = 1, nvar
          x2(a) = xa2(a) ! should be a priori
          if (sa2(a,a).NE.noretrieval) then
            x(b)  = last_oeout(a) ! use last 'good' OE pixel output
            xa(b) = xa2(a)
            x2(a) = last_oeout(a)
            b=b+1
          endif
        enddo
      else
        b=1
        do a = 1, nvar
          x2(a) = xa2(a)
          if (sa2(a,a).NE.noretrieval) then
            x(b)  = xa2(a)  ! simply use a priori values
            xa(b) = xa2(a)
            b=b+1
          endif
        enddo
      endif 

!-----NEWTONIAN ITERATION-----
      end_flag=0
      csq = 0.0 ! initialize
      niter = 0
      do n = 1, nn
        !write(*,*)'-- X: ',x2(1:6),n,swtch!,chisq !2/3/6 drz PCs
        !if(swtch>0)write(*,*)'X: ',x2(1:2),x2(5:6),n!,swtch!,chisq !2/3/6 drz PCs
        !write(*,*) x2(:),n
        call noscatter(swtch,x2,Fout)
        !if(swtch>=0)write(*,*)' Fout:',Fout !,n
        F(:) = Fout(:)

!-----CALCULATE JACOBIAN K-----
       !moditer = mod(n,2)
       !print*,'n mod 2: ',moditer
       ! uncomment to speed up code slightly (some downsides!)
       !if(moditer.eq.1) then
        do a = 1, nvar
          xprime2(a) = x2(a)              !Initialize xprime
        enddo
        a = 1
        do d = 1, nvar
         if (sa2(d,d).NE.noretrieval) then
          xprime2(d) = x2(d) + 0.02*x2(d)!//2% perturb to get a slope// 
          if(x2(d).eq.0) xprime2(d) = x2(d) + 0.02*xa2(d)
          dx(a) = xprime2(d) - x2(d)
          ! Forward model call -- EVERY TIME! (major resource hog!)
          call noscatter(swtch,xprime2,Foutprime)
          Fprime(:) = Foutprime(:) 
          dF(:) = Fprime(:) - F(:)
          K(:,a) = dF(:)/dx(a) ! Jacobian calculation!
          a = a + 1
          xprime2(d) = x2(d) ! reset for next iteration
         endif
        enddo
        savejake(:,:) = K(:,:)
        !if(swtch>0) then
          !print*,'=== jake dreo3 mm ',oepix,oelin
          !print*,savejake(:,2)
          !print*,savejake(:,3)
          !print*,savejake(:,5)
        !  print*,'---',minval(savejake(:,2)),maxval(savejake(:,2))
        !  print*,'---',minval(savejake(:,3)),maxval(savejake(:,3))
        !  print*,'---',minval(savejake(:,5)),maxval(savejake(:,5))
        !endif
       !endif
       !K(:,:) = savejake(:,:)
!          print*,'xpr2:',xprime2(1),xprime2(3),
!     >       xprime2(2),xprime2(4)

        osy(:,:)   = sy(:,:) !make sure NR Sy is used unless raining retrieval
        osy_i(:,:) = sy_i(:,:)
! new for rain retrieval: find Sy_rain and inv Sy_rain, add to Sy
! Sy(RWP) = sgn[*,*,0] * (RWP*difl[*,*])^2  [RWP<rlo]
! Sy(RWP) = syl[*,*] + sgn[*,*,1] * ((RWP-rlo)*difm[*,*])^2 [rlo<RWP<rme]
! Sy(RWP) = sym[*,*] + sgn[*,*,2] * ((RWP-rme)*difh[*,*])^2  [rme<RWP]
! where RWP is in [ g/m^2 ] !
! int factor = (RWP-rref1)/(rref2-rref1)  [rref1<RWP<rref2]
!  then f(RWP) = (1-intf)*f(rref1)+intf*f(rref2)  -- for sigmas and r value
        if(swtch>0) then
         if(abs(dlat) > latdsd) then
          if(swtch==1) rcs(:) = (/rs_rlo,rs_rme,rs_rhi/)
          if(swtch==2) rcs(:) = (/rc_rlo,rc_rme,rc_rhi/)
          if(final_rwp<=rcs(1)) then
            rref1=0.0
            rref2=rcs(1)
            sg1(:) = 0.0
            r1(:,:) = 0.0
            if(swtch==1) r2(:,:) = rs_r(:,:,1)
            if(swtch==2) r2(:,:) = rc_r(:,:,1)
            if(swtch==1) sg2(:) = rs_sgm(:,1)
            if(swtch==2) sg2(:) = rc_sgm(:,1)
          endif
          if(final_rwp<=rcs(2).and.final_rwp>rcs(1)) then
            rref1=rcs(1)
            rref2=rcs(2)
            if(swtch==1) r1(:,:) = rs_r(:,:,1)
            if(swtch==2) r1(:,:) = rc_r(:,:,1)
            if(swtch==1) r2(:,:) = rs_r(:,:,2)
            if(swtch==2) r2(:,:) = rc_r(:,:,2)
            if(swtch==1) sg1(:) = rs_sgm(:,1)
            if(swtch==2) sg1(:) = rc_sgm(:,1)
            if(swtch==1) sg2(:) = rs_sgm(:,2)
            if(swtch==2) sg2(:) = rc_sgm(:,2)
          endif
          if(final_rwp<=rcs(3).and.final_rwp>rcs(2)) then
            rref1=rcs(2)
            rref2=rcs(3)
            if(swtch==1) r1(:,:) = rs_r(:,:,2)
            if(swtch==2) r1(:,:) = rc_r(:,:,2)
            if(swtch==1) r2(:,:) = rs_r(:,:,3)
            if(swtch==2) r2(:,:) = rc_r(:,:,3)
            if(swtch==1) sg1(:) = rs_sgm(:,2)
            if(swtch==2) sg1(:) = rc_sgm(:,2)
            if(swtch==1) sg2(:) = rs_sgm(:,3)
            if(swtch==2) sg2(:) = rc_sgm(:,3)
          endif
         else ! low lat sy_dsd construction:
          if(swtch==1) rcs(:) = (/rsl_rlo,rsl_rme,rsl_rhi/)
          if(swtch==2) rcs(:) = (/rcl_rlo,rcl_rme,rcl_rhi/)
          if(final_rwp<=rcs(1)) then
            rref1=0.0
            rref2=rcs(1)
            sg1(:) = 0.0
            r1(:,:) = 0.0
            if(swtch==1) r2(:,:) = rsl_r(:,:,1)
            if(swtch==2) r2(:,:) = rcl_r(:,:,1)
            if(swtch==1) sg2(:) = rsl_sgm(:,1)
            if(swtch==2) sg2(:) = rcl_sgm(:,1)
          endif
          if(final_rwp<=rcs(2).and.final_rwp>rcs(1)) then
            rref1=rcs(1)
            rref2=rcs(2)
            if(swtch==1) r1(:,:) = rsl_r(:,:,1)
            if(swtch==2) r1(:,:) = rcl_r(:,:,1)
            if(swtch==1) r2(:,:) = rsl_r(:,:,2)
            if(swtch==2) r2(:,:) = rcl_r(:,:,2)
            if(swtch==1) sg1(:) = rsl_sgm(:,1)
            if(swtch==2) sg1(:) = rcl_sgm(:,1)
            if(swtch==1) sg2(:) = rsl_sgm(:,2)
            if(swtch==2) sg2(:) = rcl_sgm(:,2)
          endif
          if(final_rwp<=rcs(3).and.final_rwp>rcs(2)) then
            rref1=rcs(2)
            rref2=rcs(3)
            if(swtch==1) r1(:,:) = rsl_r(:,:,2)
            if(swtch==2) r1(:,:) = rcl_r(:,:,2)
            if(swtch==1) r2(:,:) = rsl_r(:,:,3)
            if(swtch==2) r2(:,:) = rcl_r(:,:,3)
            if(swtch==1) sg1(:) = rsl_sgm(:,2)
            if(swtch==2) sg1(:) = rcl_sgm(:,2)
            if(swtch==1) sg2(:) = rsl_sgm(:,3)
            if(swtch==2) sg2(:) = rcl_sgm(:,3)
          endif
         endif

          if(final_rwp>rcs(3)) then ! changed, goe10 to use lo below lo
           if(abs(dlat) > latdsd) then
            if(swtch==1) rosy(:,:) = real(rs_syh(:,:)) ! dont interpolate above!
            if(swtch==2) rosy(:,:) = real(rc_syh(:,:)) ! dont interpolate above!
           else
            if(swtch==1) rosy(:,:) = real(rsl_syh(:,:)) ! dont interpolate above!
            if(swtch==2) rosy(:,:) = real(rcl_syh(:,:)) ! dont interpolate above!
           endif
          else

            intf = (final_rwp-rref1)/(rref2-rref1)
            do j = 1, nch
             do c = 1, nch

              rosy(j,c) = real((intf*r2(j,c)+(1-intf)*r1(j,c)) * &
                          (intf*sg2(j)+(1-intf)*sg1(j))   * &
                          (intf*sg2(c)+(1-intf)*sg1(c)) )

            if(real((intf*r2(j,c)+(1-intf)*r1(j,c)))>1.0 .or. &
             (real((intf*r2(j,c)+(1-intf)*r1(j,c)))<-1.0)) then
             print*,real((intf*r2(j,c)+(1-intf)*r1(j,c))),j,c
             print*,intf,r1(j,c),r2(j,2)
             print*,final_rwp,rref1,rref2
             stop
            endif

             enddo
            enddo
          endif

          if(final_rwp<=rcs(1)) then
            rwadj(:) = 0.0
          endif
          if(final_rwp>rcs(1).and.final_rwp<=rcs(2)) then
            if(swtch==1) rwadj(:) = rwad(:,1,1)
            if(swtch==2) rwadj(:) = rwad(:,1,2)
          endif
          if(final_rwp>rcs(2).and.final_rwp<=rcs(3)) then
            if(swtch==1) rwadj(:) = rwad(:,2,1)
            if(swtch==2) rwadj(:) = rwad(:,2,2)
          endif
          if(final_rwp>rcs(3)) then
            if(swtch==1) rwadj(:) = rwad(:,3,1)
            if(swtch==2) rwadj(:) = rwad(:,3,1)
          endif

            do j = 1,nch
             do c = 1,nch
              osy(j,c) = sy(j,c)+rosy(j,c)
             enddo
             osy(j,j) = osy(j,j) + rwadj(j) ! new for later goe10
            enddo
        !print*,'r ',final_rwp,rref1,rref2,intf
        !print*,r1(:,:)
        !print*,r2(:,:)
        !stop

        !print*,rs_rlo,rs_rme,rs_rhi
        !print*,rc_rlo,rc_rme,rc_rhi
        !print*,'Sy:   ',sy(:,:)
        !print*,'osy:  ',osy(:,:)
        !print*,'---- sqrt sy, rosy, osy'
          call inverse(osy,osy_i,nch)
        !do c=1,13 
        ! print*,'#-#',sy(c,c),rosy(c,c),osy(c,c)!,rs_sym(c,c)
        ! print*,sqrt(sy(c,c)),sqrt(rosy(c,c)),sqrt(osy(c,c))!,sqrt(rs_sym(c,c))
        !enddo
        !stop
        endif


!-----EVALUATE NEWTONIAN MATRICES-----   
        !write(*,*)'1 args',K
        K_t(:,:) = transpose(K(:,:))
        !write(*,*)'1 args',K_t
        !write(*,*)'2 args',sa
        call inverse(sa,sa_i,nretvar)
        !write(*,*)'2 args',sa
        !write(*,*)'2 args',sa_i
        !write(*,*)'3 args',K_t,osy_i,
        prod1 = matmul(K_t,osy_i)
        !write(*,*)'5 args',prod1
        prod2 = matmul(prod1,K)
        !write(*,*)'6 args',prod2
        sx_i = sa_i + prod2
        !write(*,*)'7 args',sx_i
        call inverse(sx_i,sx,nretvar)
        !write(*,*)'8 args',sx

!-----PERFORM NEWTONIAN STEP-----
        sum1(:,1) = xa(:) - x(:)
        sum2(:,1) = y(:) - F(:)
        prod3 = matmul(sa_i,sum1(:,1))
        !print*,'sa ',sa(:,:)
        !print*,'prod3 ',prod3(:)
        prod4 = matmul(prod1,sum2(:,1))
        sum3(:,1) = prod4(:) + prod3(:)
        !print*,'sum3 ',sum3(:,1)
        prod5 = matmul(sx(:,:),sum3(:,1))
        xnext = x + prod5
        !print*,'prod5 ',prod5(:)
        !print*,'xnex1 ',xnext(:)

!-----ERROR DIAGNOSTICS-----
        AP = matmul(sx,prod2) ! Sx * (K_t * Sy_i * K) -- yields DOF
        sum4(:,1) = F(:) - y(:)
        sum4_t = transpose(sum4)
        sum5(:,1) = x(:) - xa(:)
        sum5_t = transpose(sum5)
        prod8 = matmul(sum4_t(:,:),osy_i(:,:))
        !print*,'(F-y)_t',sum4_t(1,:)
        !print*,'sy',sy(:,:)
        !print*,'rosy',rosy(:,:)
        !print*,'osy',osy(:,:)
        !print*,'osy_i',osy_i(:,:)
        prod9 = matmul(prod8(:,:),sum4(:,:)) !(F-y)_t * osy_i * (F-y)
        prod10 = matmul(sum5_t(:,:),sa_i(:,:))
        prod11 = matmul(prod10(:,:),sum5(:,:))
        chisqtot(1,1) = prod9(1,1) + prod11(1,1)
        !Check limits on xnext
        do a=1,nretvar
          if (xnext(a) .gt. xmax(a)) then
            xnext(a) = xmax(a)
          else if (xnext(a) .lt. xmin(a)) then
            xnext(a) = xmin(a)
          endif
        end do
        !print*,'xnex2 ',xnext(:)
        
        ! Evaluate closeness of Xnext and X
        xdiff(:,1) = xnext(:) - x(:)
        !print*,'xdiff ',xdiff(:,1)
        !print*,'x ',x(:)
        xdiff_t(:,:) = transpose(xdiff(:,:))
        prod6(1,:) = matmul(xdiff_t(1,:),sx_i(:,:))
        prod7(:,:) = matmul(prod6(:,:),xdiff(:,:))

        b = 1
        do a = 1, nvar
         if (sa2(a,a).NE.noretrieval) then
          x(b)  = xnext(b)
          x2(a) = xnext(b)
          b = b + 1
         endif
        end do

! ***As final retrieved state is found, TBs need to be
! ***recomputed.
        !print*,'xnex ',xnext(:)
        call noscatter(swtch,x2,Fout)
        F(:) = Fout(:)
        !print*,'Fout: ',Fout
        !chisq = prod9(1,1) !chisqtot(1,1) - prod11(1,1) 
        !chisq = chisqtot(1,1) - prod11(1,1) 

        if (n .gt. 1) then
          diffF(:,1) = F(:) - Fold(:)
          !print*,'Fout: ',Fout
          !print*,'Fold: ',Fold
          !print*,'diffF: ',diffF
          diffF_t = transpose(diffF)

          Ksa = matmul(K,sa)
          Ksakt = matmul(Ksa,K_t)
          !Sdysum = Ksakt + Sy !!WRONG!
          Sdysum = Ksakt + osy  !K*Sa*KT + (Sy+Sy_rain)
          call inverse(Sdysum,Sdysum_i,nch)
          almost = matmul(osy,Sdysum_i) 
          Sdy = matmul(almost,osy) ! Se*(K*Sa*KT+Se)*Se
          call inverse(Sdy,Sdy_i,nch)
          almost3(1,:) = matmul(diffF_t(1,:),Sdy_i(:,:)) ! try this????
          disq = matmul(dble(almost3),dble(diffF))
          !if(disq(1,1)<0) then
          !      print*,'++++++++++++++',oelin,oepix
          ! print*,disq(1,1),swtch
          !  print*,K
          !  stop
          !endif
        csq = (chisqtot(1,1) - prod11(1,1))/float(nch)
        !print*,n,disq(1,1),prod9(1,1),prod11(1,1)
        !print*,chisq!,prod9(1,1),prod11(1,1)

        ! first is old convergence criteria, 5.33 is new.
        !print*,'disq',disq(1,1),n
          !if (prod7(1,1)<nretvar/conv_factor .and. swtch>0) then ! Eq 5.29 in Rodgers 2000 
          !  end_flag=1
          !  if(swtch>0) print*,swtch,'$$$5.29conv$$$$',oepix
          !endif
          if (disq(1,1) < dble(nch/conv_factor).and.disq(1,1)>0) then ! Eq 5.33 in Rodgers 2000 
             !if(swtch>0) print*,swtch,'normal convergence!?',oepix
             end_flag=1
          endif
          !if(n>9.and.chisq>0.0.and.chisq<2.0) then 
          !  end_flag=1
          !  if(swtch>0)print*,'other conv crit,ch,disq',chisq,disq(1,1),n
          !endif
        ! secondary convergence crit. for 10+ iter, low chisq
          if(csq<1.0 .and. swtch>0) then 
             end_flag=1
             !print*,swtch,'#**3rd crit met**#',oepix
          endif
          if(disq(1,1)>=nch/conv_factor.and.n>=(nn-2).and.csq<chout) then 
             end_flag=1
             !if(swtch>0)print*,swtch,'#second crit##',oepix,disq(1,1)
          endif
          !if(disq(1,1) < 0 .AND. swtch==0) then ! should be impossible...
          !  print*,'convergence error -- negative value not allowed!'
          !  print*,'state ',x2(1:5)
          !  print*,'diffF ',diffF(:,1)
          !  print*,'xdif ',xdiff(:,1)
          !  print*,oepix,oelin,disq,swtch
          !  end_flag=1 ! essentially pass it off to raining retrieval
          !endif
          if(csq<0) then
            print*,'negative chisq!',oepix,oelin,swtch
            print*,n,disq(1,1),prod9(1,1),prod11(1,1)
            !print*,'F-y:',sum4(:,:)
            !print*,'rwp: ',final_rwp
        !open(unit=3,file='syo.bin',access='stream',status='unknown')
        !write(unit=3)rosy
        !write(unit=3)sy_i
        !write(unit=3)osy
        !write(unit=3)osy_i
        !close(unit=3)
        !stop
            flag=0
          endif
          if(disq(1,1) < 0) then ! should be impossible...
            print*,'(&)negative disq! ',final_rwp,swtch
            !print*,'(&)negative disq! ',oepix,oelin,swtch
            !print*,n,disq(1,1),prod9(1,1),prod11(1,1)
            !end_flag=0 
            !csq = 99
            !return
          endif
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (end_flag==1) then
        !chisq = chisqtot(1,1) - prod11(1,1) 
          b = 1
          do a = 1,nvar 
           if (sa2(a,a).NE.noretrieval) then	  
            Amatrix2(a)=AP(b,b) 
            sigma2(a)  =sqrt(sx(b,b)) ! -- posterior error standard deviations
            b = b + 1
           else
            Amatrix2(a) = 0.0
            sigma2(a)   = miss_flt 
           endif
          enddo
          flag=1
          niter = n
          if(csq == 99) flag=0 !sign of negative disq, bad iteration!
          !if(swtch>0) print*,'+++ sumAmat: ',sum(Amatrix2(1:6)),swtch
          !if(swtch==0) print*,'--- sumAmat: ',sum(Amatrix2(:))
          !if(swtch>0) print*,'+++ Amat: ',Amatrix2(:)!1:2),Amatrix2(5:6)
          !print*,'f=1, chi/sw/px: ',chisq,swtch,oepix
          return
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !chisq = (chisqtot(1,1) - prod11(1,1))/float(nch)

        Fold(:) = Fout(:)

      end do  !!!/// end of do loop starting way up above

!      write(*,*) 'Solution not found in ',nn,' steps'
!      No solution found upon reaching max iterations?  Set the following to
!      -99.99.
      do a = 1,nvar
       Amatrix2(a) = -99.99
       sigma2(a)   = -99.99
      enddo
      csq      = -99.99
      flag       = 0
      niter      = n
      return

      end subroutine opt_est

!-----------------------------------------------------------
      end module csu1dvar
