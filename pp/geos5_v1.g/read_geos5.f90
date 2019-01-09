      module read_geos5

!------ contains routine to read in netCDF GEOS5 FP IT file,
!------  and assign model grid values to satellite pixels.
!------  D Duncan, CSU, 11/2/16, last updated 11/3/16

      use netcdf
      use pp_definition

      implicit none

      contains

      subroutine geos5 

      logical :: fexists

      !integer,parameter :: nz=16 ! 16 levels, need 17 values
      integer :: ncid,ncid2,ncid3,ncid4
      character(len=20)  :: ppversion

      integer, parameter :: glo=576, gla=361, glev=42 ! fixed geos5 field sizes
      real*8 :: glocz(glo), glacz(gla) 
      integer :: glocz_varid, glacz_varid
      real :: gsp(glo,gla), gt(glo,gla,glev), gq(glo,gla,glev), grh(glo,gla,glev) ! from 3d
      real :: gh(glo,gla,glev) ! new. from 3d
      integer :: sp_varid, t_varid, q_varid, rh_varid,  h_varid
      real :: gq2m(glo,gla), gt2m(glo,gla), gts(glo,gla), gu(glo,gla),gslp(glo,gla) ! from 2d
      real :: gv(glo,gla),gtpw(glo,gla),  snod(glo,gla),snof(glo,gla),gsi(glo,gla) ! from 2d
      integer :: q2m_varid, t2m_varid, ts_varid,  si_varid
      integer :: u_varid, v_varid, tpw_varid, slp_varid, snod_varid,snof_varid 
        ! secondary variables
      real :: gwmag(glo,gla),gwdir(glo,gla), gtwb(glo,gla)
      real :: gwvmr(glo,gla,nz),gairt(glo,gla,nz), ghi(glo,gla,nz)
      real :: smixr, avq, glod,glad, rh, pfrsno,psnod,psci
      real*8 :: dpi = 3.14159265359
      integer :: p,s, gglo, ggla, pmo,pma, z, ila, ilo
      integer :: sub(nz+1) = (/25,23,21,19,17,16,15,14,13,11,9,7,5,4,3,2,1/)
        ! vars for getting date/time right
      character(len=4) :: yr, tim3,tim2
      character(len=2) :: mos(12) = (/'01','02','03','04','05','06',&
                                      '07','08','09','10','11','12'/)
      character(len=2) :: das(31) = (/'01','02','03','04','05','06',&
      '07','08','09','10','11','12','13','14','15','16','17','18',&
      '19','20','21','22','23','24','25','26','27','28','29','30','31'/)
      character(len=2) :: hos3(8) = (/'00','03','06','09','12','15','18','21'/)
      character(len=2) :: hos2(24)= (/'00','01','02','03','04','05','06',&
        '07','08','09','10','11','12','13','14','15','16','17','18',&
        '19','20','21','22','23'/)
                                      
      integer :: midscan, hodex3,hodex2
      character(len=2) :: mo, da, ho3,ho2

        ! experimental, fitting pixel level to lsmask,snowgrid resolutions...
      real,parameter    :: srinc    = .03125, sppd = 16. !srinc= 1/2 spacing
      integer,parameter :: nlon     = 5760,  nlat = 2880
      real :: newlo, snowres= 0.64
      integer :: snowj,snowi,ilat,ilon
      
! read in input/output filenames, then get npix/nscans from input file
! determine closest time step to take from analysis (different for 3D, 2D fields!)
      midscan = nint(nscans/2.0)
      write(yr,'(I4)') stdtime(midscan,1) ! year, 4char
      mo = mos(stdtime(midscan,2)) !use as index
      da = das(stdtime(midscan,3)) !use as index
      hodex3 = nint((stdtime(midscan,4) + stdtime(midscan,5)/60.0)/3.0)+1 !3D hour, instantaneous!
      hodex2 = nint( stdtime(midscan,4) - 0.5 +stdtime(midscan,5)/60.0)+1 !2D hour, time-averaged
      if(hodex3.le.0) hodex3=1
      if(hodex2.le.0) hodex2=1
      ho3 = hos3(hodex3)
      ho2 = hos2(hodex2)

        print*,'midscan time: ',stdtime(midscan,2:5)
        print*,'3D date/time: ',yr//mo//da//'_'//ho3//'00'
        print*,'2D date/time: ',yr//mo//da//'_'//ho2//'30'
      gdir='../input/'
      !gdir='/cdata2/archive/GEOS5/'//trim(yr)//'/'//trim(yr(3:4))//mo//'/'
        print*,'geos5 dir: ',trim(gdir)
      gf1='GEOS.fpit.asm.inst3_3d_asm_Np.GEOS5124.'//yr//mo//da//'_'//ho3//'00.V01.nc4' !hardcoded, obviously
      gf2='GEOS.fpit.asm.tavg1_2d_slv_Nx.GEOS5124.'//yr//mo//da//'_'//ho2//'30.V01.nc4'
      gf3='GEOS.fpit.asm.tavg1_2d_lnd_Nx.GEOS5124.'//yr//mo//da//'_'//ho2//'30.V01.nc4'
      gf4='GEOS.fpit.asm.tavg1_2d_flx_Nx.GEOS5124.'//yr//mo//da//'_'//ho2//'30.V01.nc4'
        !print*,gf1,gf2,gf3,gf4

      inquire(file=trim(gdir)//trim(gf1),exist = fexists)
      if (fexists .eq. .false.) then
        write(6,*) 'Input file ',trim(gdir)//trim(gf1),' not found!'
      endif

!--- start opening and reading netcdf file
      call check(nf90_open(trim(gdir)//trim(gf1), nf90_nowrite, ncid))

!---  get variable IDs from input file
      call check(nf90_inq_varid(ncid,"QV",q_varid)) !specific humidity
      call check(nf90_inq_varid(ncid,"PS",sp_varid)) !surface pressure
      call check(nf90_inq_varid(ncid,"T",t_varid)) ! air temp
      call check(nf90_inq_varid(ncid,"H",h_varid)) ! height
      !call check(nf90_inq_varid(ncid,"RH",rh_varid)) ! relative humidity

!---  read variables
      call check(nf90_get_var(ncid,q_varid,gq))
      call check(nf90_get_var(ncid,sp_varid,gsp))
      call check(nf90_get_var(ncid,t_varid,gt))
      call check(nf90_get_var(ncid,h_varid,gh)) ! new
      !call check(nf90_get_var(ncid,rh_varid,grh))

!---  close input file
      call check(nf90_close(ncid))

!--- start opening and reading netcdf file
      call check(nf90_open(trim(gdir)//trim(gf2), nf90_NoWrite, ncid2))

!---  get variable IDs from input file
      call check(nf90_inq_varid(ncid2,"QV2M",q2m_varid)) !2m spec hum
      call check(nf90_inq_varid(ncid2,"T2M",t2m_varid)) !2m temp
      call check(nf90_inq_varid(ncid2,"TS",ts_varid)) !skin temp
      call check(nf90_inq_varid(ncid2,"U10M",u_varid)) 
      call check(nf90_inq_varid(ncid2,"V10M",v_varid)) 
      call check(nf90_inq_varid(ncid2,"TQV",tpw_varid)) !total_precipitable_water_vapor
      call check(nf90_inq_varid(ncid2,"SLP",slp_varid)) !SLP (Pa)

!---  read variables
      call check(nf90_get_var(ncid2,q2m_varid,gq2m))
      call check(nf90_get_var(ncid2,t2m_varid,gt2m))
      call check(nf90_get_var(ncid2,ts_varid,gts))
      call check(nf90_get_var(ncid2,u_varid,gu))
      call check(nf90_get_var(ncid2,v_varid,gv))
      call check(nf90_get_var(ncid2,tpw_varid,gtpw))
      call check(nf90_get_var(ncid2,slp_varid,gslp))

!---  close input file2
      call check(nf90_close(ncid2))

!--- start opening and reading netcdf file
      call check(nf90_open(trim(gdir)//trim(gf3), nf90_NoWrite, ncid3))
        
!---  get variable IDs from input file
      call check(nf90_inq_varid(ncid3,"SNODP",snod_varid)) !snow depth (m)
      call check(nf90_inq_varid(ncid3,"FRSNO",snof_varid)) !snow fraction (0.0-1.0)
      call check(nf90_inq_varid(ncid3,"lon",glocz_varid))
      call check(nf90_inq_varid(ncid3,"lat",glacz_varid))

!---  read variables
      call check(nf90_get_var(ncid3,snod_varid,snod))
      call check(nf90_get_var(ncid3,snof_varid,snof))
      call check(nf90_get_var(ncid3,glocz_varid,glocz))
      call check(nf90_get_var(ncid3,glacz_varid,glacz))

!---  close input file3
      call check(nf90_close(ncid3))

!--- open 4th file, which contains Sea Ice Fraction
      call check(nf90_open(trim(gdir)//trim(gf4), nf90_nowrite, ncid4))

!---  get variable IDs from input file
      call check(nf90_inq_varid(ncid4,"FRSEAICE",si_varid)) !sea ice fraction (0-1.0)

!---  read variables
      call check(nf90_get_var(ncid4,si_varid,gsi))

!---  close input file
      call check(nf90_close(ncid4))

!--- start opening and reading netcdf file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! now calculate quantities desired
      do ilo = 1, glo
        do ila = 1, gla

          gwmag(ilo,ila) = sqrt(gu(ilo,ila)**2 + gv(ilo,ila)**2) ! 10m wind
          gwdir(ilo,ila) = 0.0 ! initialize in case u10m=0.0
          if(gu(ilo,ila).ne.0.0) then
            gwdir(ilo,ila) = real(atan(gv(ilo,ila)/gu(ilo,ila))*180.0/dpi)
            if(gu(ilo,ila)>0) gwdir(ilo,ila)= 270-gwdir(ilo,ila)
            if(gu(ilo,ila)<0) gwdir(ilo,ila)= 90 -gwdir(ilo,ila)
          endif
        !print*,'mag,dir: ',gwmag(ilo,ila),gwdir(ilo,ila)

          do z = 1, nz !average temp/mixratio in layer (what crtm needs)
            gairt(ilo,ila,z) = sum(gt(ilo,ila,sub(z+1):sub(z)))/float(sub(z)-sub(z+1)+1)
            avq = sum(gq(ilo,ila,sub(z+1):sub(z)))/float(sub(z)-sub(z+1)+1)
            gwvmr(ilo,ila,z) = avq/(1.0-avq*.001)*1000.0 !g/kg
            if(gt(ilo,ila,sub(z+1))>400) gairt(ilo,ila,z)=gairt(ilo,ila,z-1)
            if(gq(ilo,ila,sub(z+1))>1) gwvmr(ilo,ila,z)=gwvmr(ilo,ila,z-1)
            ghi(ilo,ila,z) = gh(ilo,ila,sub(z)) ! height OF TOP OF LAYER
            if(ghi(ilo,ila,z)<0) ghi(ilo,ila,z)=0
            if(gh(ilo,ila,sub(z))>25000) ghi(ilo,ila,z)=gh(ilo,ila,sub(z-1))
          enddo
          
          smixr = gq2m(ilo,ila)/(1.0-gq2m(ilo,ila)*.001)*1000.0 ! sfc mixratio [g/kg]
          rh = (gsp(ilo,ila)*.01)/(6.112*exp(17.27*(gt2m(ilo,ila)-273.15)/(gt2m(ilo,ila)-273.15+237.3))) * &
                (smixr*.001*28.966/18.016)/(1-(smixr*.001*28.966/18.016)) *100.0 
          if(rh>100.0) rh=100.0

          gtwb(ilo,ila) = (gt2m(ilo,ila)-273.15) * atan(0.151977*(rh+8.314)**0.5) + &
            atan((gt2m(ilo,ila)-273.15)+rh) - atan(rh-1.676)+ &
            (0.0039184*(rh**1.5) *atan(0.0231*rh)) - 4.686 ! from stull 2011
          gtwb(ilo,ila) = gtwb(ilo,ila) + 273.15 ! into Kelvin from C
        !print*,'Twb/t2m/rh: ',gtwb(ilo,ila)-273.15,gt2m(ilo,ila)-273.15,rh
        enddo
      enddo
        
      ! allocate pix/scan output variables
      allocate(tprof(npixl,nscans,nz),wvmr(npixl,nscans,nz))
      allocate(height(npixl,nscans,nz))
      allocate(wind(npixl,nscans),winddir(npixl,nscans))
      allocate(slp(npixl,nscans),tcwv(npixl,nscans))
      allocate(twb(npixl,nscans),t2m(npixl,nscans))
      allocate(skint(npixl,nscans),sic(npixl,nscans))
      allocate(sfccode(npixl,nscans),snowci(npixl,nscans))
      sfccode(:,:) = 0 ! for now...
      snowci(:,:) = 0 ! for now...

          ! lon goes -180 to 180 by 0.625deg
          ! lat goes -90 to 90 by 0.50deg
      ! now with all necessary variables, assign to satellite pixels
      do s = 1, nscans
        do p = 1, npix

          !gglo = floor((lon(p,s)+180)/.625)+1  !GEOS5 grid index, lon
          !ggla = floor(abs(lat(p,s)+90)/.5)+1 !GEOS5 grid index, lat
          gglo = nint((lon(p,s)+180)/.625)+1  !GEOS5 grid index, lon
          ggla = nint(abs(lat(p,s)+90)/.5)+1 !GEOS5 grid index, lat
          glod = (  lon(p,s) - glocz(gglo) )/.625
          glad = (  lat(p,s) - (glacz(ggla))  )/.5
          !glod = ( (lon(p,s)+180) - (gglo*.625 - .625/2.0) )/.625
          !glad = (  lat(p,s)      - (ggla*.5-90- .5/2.0)  )/.5
          pmo = nint(abs(glod)/glod)
          pma = nint(abs(glad)/glad)  ! just sign
          !pma = -1*nint(abs(glad)/glad)  ! just sign
          if((gglo+pmo)<=0 .or. (gglo+pmo)>576) pmo=0
          if((gla+pma) <=0 .or. (ggla+pma)>361) pma=0
        ! test to see if lat/lon matches up with lat/lon from GEOS5...
        !print*,'G lo,la: ',glocz(gglo),glacz(ggla)
        !print*,'glod,glad ',glod*.625,glad*.5,pmo,pma
        !print*,'p/s lo,la: ',lon(p,s),lat(p,s),' ',p,s
        !print*,'lo,la dif: ',(glocz(gglo)+abs(glod)*(glocz(gglo+pmo)-glocz(gglo)))-lon(p,s),(glacz(ggla)+abs(glad)*(glacz(ggla+pma)-glacz(ggla)))-lat(p,s),p,s
        !stop

          winddir(p,s)= nint(gwdir(gglo,ggla)) ! no interpolation, can't average degrees!
          ! interpolation between grid points, in i and j directions
          wind(p,s)= gwmag(gglo,ggla) + & 
           abs(glod)*(gwmag(gglo+pmo,ggla)-gwmag(gglo,ggla))+abs(glad)*(gwmag(gglo,ggla+pma)-gwmag(gglo,ggla))
        !print*,'closest ',gwmag(gglo,ggla)
        !print*,'p la/lo, g la/lo: ',lon(p,s),lat(p,s),glocz(gglo),glacz(ggla)
        !print*,gwmag(gglo-1,ggla+1),gwmag(gglo,ggla+1),gwmag(gglo+1,ggla+1)
        !print*,gwmag(gglo-1,ggla),gwmag(gglo,ggla),gwmag(gglo+1,ggla)
        !print*,gwmag(gglo-1,ggla-1),gwmag(gglo,ggla-1),gwmag(gglo+1,ggla-1)
        !print*,gglo,gglo+pmo,ggla,ggla+pma
        !print*,glod,glad
        !print*,glocz(gglo-1),glocz(gglo),glocz(gglo+1)
        !print*,glacz(ggla+1)
        !print*,glacz(ggla)
        !print*,glacz(ggla-1)
        !stop
          slp(p,s) = 0.01 * (gslp(gglo,ggla) + &  !Pa to hPa
           abs(glod)*(gslp(gglo+pmo,ggla)-gslp(gglo,ggla))+abs(glad)*(gslp(gglo,ggla+pma)-gslp(gglo,ggla)))
          twb(p,s) = gtwb(gglo,ggla) + & 
           abs(glod)*(gtwb(gglo+pmo,ggla)-gtwb(gglo,ggla))+abs(glad)*(gtwb(gglo,ggla+pma)-gtwb(gglo,ggla))
          tcwv(p,s) = gtpw(gglo,ggla) + & 
           abs(glod)*(gtpw(gglo+pmo,ggla)-gtpw(gglo,ggla))+abs(glad)*(gtpw(gglo,ggla+pma)-gtpw(gglo,ggla))
          t2m(p,s)  = gt2m(gglo,ggla) + & 
           abs(glod)*(gt2m(gglo+pmo,ggla)-gt2m(gglo,ggla))+abs(glad)*(gt2m(gglo,ggla+pma)-gt2m(gglo,ggla))
          skint(p,s)= gts(gglo,ggla) + & 
           abs(glod)*(gts(gglo+pmo,ggla)-gts(gglo,ggla))+abs(glad)*(gts(gglo,ggla+pma)-gts(gglo,ggla))
          ! frseaice no missing value, just 0.0 over land
          sic(p,s)= nint(100.0 * (gsi(gglo,ggla) + & 
           abs(glod)*(gsi(gglo+pmo,ggla)-gsi(gglo,ggla))+abs(glad)*(gsi(gglo,ggla+pma)-gsi(gglo,ggla))))
          wvmr(p,s,:) = gwvmr(gglo,ggla,:) + & 
           abs(glod)*(gwvmr(gglo+pmo,ggla,:)-gwvmr(gglo,ggla,:))+abs(glad)*(gwvmr(gglo,ggla+pma,:)-gwvmr(gglo,ggla,:))
          tprof(p,s,:)= gairt(gglo,ggla,:) + & 
           abs(glod)*(gairt(gglo+pmo,ggla,:)-gairt(gglo,ggla,:))+abs(glad)*(gairt(gglo,ggla+pma,:)-gairt(gglo,ggla,:))
          height(p,s,:)= nint(ghi(gglo,ggla,:) + & 
           abs(glod)*(ghi(gglo+pmo,ggla,:)-ghi(gglo,ggla,:))+abs(glad)*(ghi(gglo,ggla+pma,:)-ghi(gglo,ggla,:)))

          ! designate snow index here...
          pfrsno = 0.0
          psnod  = 0.0
          pfrsno = snof(gglo,ggla) + &
           abs(glod)*(snof(gglo+pmo,ggla)-snof(gglo,ggla))+abs(glad)*(snof(gglo,ggla+pma)-snof(gglo,ggla))
          psnod = snod(gglo,ggla) !+ &
           !abs(glod)*(snof(gglo+pmo,ggla)-snof(gglo,ggla))+abs(glad)*(snof(gglo,ggla+pma)-snof(gglo,ggla))
          if(snof(gglo,ggla+pma)>1) then ! remove 1e15 mis_val from interpolation
            pfrsno = snof(gglo,ggla) + abs(glod)*(snof(gglo+pmo,ggla)-snof(gglo,ggla))
            if(pfrsno>0) psnod = snod(gglo,ggla) !+ abs(glod)*(snod(gglo+pmo,ggla)-snod(gglo,ggla))
          endif
          if(snof(gglo+pmo,ggla)>1) then ! remove 1e15 mis_val from interpolation
            pfrsno = snof(gglo,ggla) + abs(glad)*(snof(gglo,ggla+pma)-snof(gglo,ggla))
            if(pfrsno>0) psnod  = snod(gglo,ggla) !+ abs(glad)*(snod(gglo,ggla+pma)-snod(gglo,ggla))
          endif
          if(snof(gglo+pmo,ggla)>1 .and. snof(gglo,ggla+pma)>1) then
            pfrsno=snof(gglo,ggla)
            if(pfrsno>0) psnod =snod(gglo,ggla)
          endif
          if(snof(gglo,ggla)>1) pfrsno = 0.0
          
          psci = pfrsno * psnod
          if(psci<0.02) snowci(p,s)=0  !definition of snow emiss classes -- not set in stone!!
          if(psci>0.02 .and. psci<0.2) snowci(p,s)=1
          if(psci>=0.2 .and. psci<0.5) snowci(p,s)=2
          if(psci>=0.5) snowci(p,s)=3

          ! quality controls...
          if(wind(p,s)<0) wind(p,s) = 0.0
          if(twb(p,s)>t2m(p,s)) twb(p,s) = t2m(p,s)
          if(sic(p,s)<0) sic(p,s) = 0
          if(sic(p,s)>100) sic(p,s) = 100
          if(minval(wvmr(p,s,:))<0) then
           do z = 1, nz
            if(wvmr(p,s,z)<0.0) wvmr(p,s,z) = 0.0 ! no negative WV!
           enddo 
          endif
          if(maxval(wvmr(p,s,:))>32.7) then
           do z = 1, nz
            if(wvmr(p,s,z)>32.7) wvmr(p,s,z) = 32.7 ! not allowed by compression to int
           enddo 
          endif
          if(maxval(tprof(p,s,:))>327.0) then
           do z = 1, nz
            if(tprof(p,s,z)>327) tprof(p,s,z) = 327.0 !not allowed by compression to int
           enddo 
          endif
        enddo
      enddo
      write(6,*) 'almost done read_geos5...'

      sfcgrid = 0 ! initialize
      snowgrid = miss_byt ! initialize
      do p = 1,npix
        do s = 1, nscans

        ! first do snowgrid fakeout
        snowi = int(1+ float(ilon-1)/snowres)
        snowj = int(1+ float(ilat-1)/snowres)
        if(snowi>9000) snowi=9000
        if(snowj>4500) snowj=4500
        if(snowi<=1) snowi=1
        if(snowj<=1) snowj=1
        if(snowgrid(snowi,snowj).ne. miss_byt) snowgrid(snowi,snowj)=snowci(p,s)

        ! now do sfcgrid definition from lo_flag
            newlo = lon(p,s)
            if(newlo > 180) newlo = newlo - 360.

            if(lat(p,s) > -91.0 .and. newlo >= -180.0) then !check for good lat/lon
              ilat= nint(sppd* ( 90.0 + srinc - lat(p,s))) !lat-lons in
              ilon= nint(sppd* ( 180.0 + srinc + newlo)) !-180:180 space          
              if(ilon .eq. nlon+1) ilon = 1
              if(ilat .eq. nlat+1) ilat = nlat
              
              if(lo_flag(p,s,rz)>=0 .and. lo_flag(p,s,rz)<=2) sfcgrid(ilon,ilat) = 10 !0-2% land = ocean
              if(lo_flag(p,s,rz)>2 .and. lo_flag(p,s,rz)<=95) sfcgrid(ilon,ilat) = 30 !2-95% land = coast
              if(lo_flag(p,s,rz)>95) sfcgrid(ilon,ilat) = 20 !96-100% land = land
            endif

        enddo
      enddo

      write(6,*) 'done read_geos5'
      !stop 'stop at end of read_geos5'
      end subroutine geos5

!______________________________________________________________________

      subroutine check(status)
      integer, intent ( in) :: status
      if (status /= nf90_noerr) then 
         print *, trim(nf90_strerror(status))
         stop "Stopped"
      end if
      end subroutine check  

      subroutine chkStatus(status, prtMsg)
!---  this routine simply checks the status code and writes an Error
!---  message if status < 0,  These are fatal errors so process stops.      
      integer           :: status
      character(len=*)  :: prtMsg
!      write(6,*) status, prtMsg
      if(status < 0) then
           write(*,*) 'FATAL error from: ', trim(prtMsg)
           stop
      endif
      return
      end subroutine chkStatus

      end module read_geos5
