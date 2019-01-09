      module csu1dvar_subr
        ! removed output routines, include back later

       use define_csu1dvar
       use define_intam2
       use GPM_arraydef_D
       !use csu1dvar

       implicit none

       contains


      subroutine read_reynolds_sst
       implicit none

!-----------------------------------------------------------------
! This routine will read in a Reynolds 0.25X0.25 daily SST and Sea 
! ice file selected from the time as specified in stdtime.  If the 
! individual SST file is not available, the most recent file is used.
!-----------------------------------------------------------------
    
       integer(2),parameter :: xsize=1440, ysize=720
       logical :: fexists
       character(len=256):: tfile
       character(len=100) :: dir_sst 
       character(len=8)  :: cdate
       character(len=4)  :: cyear
       integer :: i, j, k, m, read_date
       integer :: iyear,imon,iday, iyr,imo,ida       
       
       integer :: rlun, ios
       real    :: sbuff(1440)
       integer(2) :: ibuff(1440), err, anom            
       integer :: numday(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)

       integer(2),allocatable :: fice(:,:)
       integer(2),allocatable :: fsst(:,:)

       dir_sst = 'binary/' 
       !dir_sst = '/xdata/drandel/gpm/ppingest/sstdata_0.25/' !hard-coded!
       !allocate (icegrid(xsize,ysize) )   
       allocate (sstgrid(xsize,ysize) )  
       allocate (fsst(xsize,ysize) )         !temp sst from file
       allocate (fice(xsize,ysize) )         !temp sea-ice from file
                
!--- Find the Reynolds SST filename for this date 
       read_date = first_good_date       
       !write(*,*) '   read_date = ', read_date
        
       do i = 1,30                !search back 30 days for file
       
         write(cdate,'(i8)') read_date
         read(cdate(1:4),'(i4)') iyear
         read(cdate(5:6),'(i2)') imon
         read(cdate(7:8),'(i2)') iday
             !print*,cdate,iyear,imon,iday 
         if(read_date.lt.20020601 .or. read_date.gt.20111004) then
!AVHRR+AMSR starts this date               
               tfile = trim(dir_sst) // cdate(1:4) // 
     >                 '/avhrr-only-v2.' // cdate
               inquire(file = tfile,exist = fexists)     
         else               
               tfile = trim(dir_sst) // cdate(1:4) // 
     >                 '/amsr-avhrr-v2.' // cdate
               inquire(file = tfile,exist = fexists)
         endif 
         if(i .eq. 5) then
            write(*,*)' WARNING: sst more than 5 days old'
         endif  
         if(fexists) exit              !file for requested date found

         if(iday-1 .gt. 0) then        !if no file, look back in time
            read_date = read_date - 1  !still within month
         else            
            imon = imon - 1               !look in previous month
            if(imon .lt. 1) then          !look in previous year
                imon = 12
                iyear = iyear - 1
            endif
            if(mod(iyear,4) .eq. 0) numday(2) = 29
            iday = numday(imon)
            read_date = iyear *10000 + imon*100 + iday           
         endif   
       enddo
       if(.not. fexists) then       
            write(*,*)'no SST file : ',trim(tfile)
       endif 

       call gprof_lun(rlun)
!      write(*,*)'   reading sst/ice, filename    : ',
!    >                  trim(tfile)
       open(rlun,file=trim(tfile), form='unformatted',
     >      status='old', access='sequential',convert='big_endian',
     >      iostat=ios)
       if(ios.ne.0) write(*,*)'error opening SST file'

!--- read in SST and ICE grids
       
       read(rlun,iostat=ios)iyr,imo,ida,
     >           ((fsst(i,j),i=1,xsize),j=ysize,1,-1) !flip SST for NP up
       !write(*,*)'writing out sst grid!: ',fsst
       read(rlun,iostat=ios)iyr,imo,ida,
     >           ((anom,i=1,xsize),j=1,ysize)  !don't use anom field
       read(rlun,iostat=ios)iyr,imo,ida,
     >           ((err,i=1,xsize),j=1,ysize)  !don't use err field
       if(ios.ne.0) write(*,*)'issue 3' 
       read(rlun,iostat=ios)iyr,imo,ida,
     >           ((fice(i,j),i=1,xsize),j=ysize,1,-1) !flip ICE for NP up
       if(ios.ne.0) write(*,*)'issue 4'
       
!--- close file
       close(rlun)
      
!--- shift and scale sst field to center on 0 degrees

       do j = 1,ysize
         do i = 1,720
           if(fsst(i+720,j) .ne. -999) then
               sstgrid(i,j) = fsst(i+720,j) / 100.
           else
               sstgrid(i,j) = miss_flt
           endif
           
           if(fsst(i,j) .ne. -999) then    
               sstgrid(i+720,j)= fsst(i,j) / 100.
           else
               sstgrid(i+720,j) = miss_flt
           endif
         enddo
       enddo
       
!--- shift ice field to center on 0 degrees
     
!       do j = 1,ysize
!         do i = 1,720
!           ibuff(i)     = fice(i,j)
!           fice(i,j)    = fice(i+720,j)
!           fice(i+720,j)= ibuff(i)
!         enddo
!       enddo
!                         
!!--- assign ice to icegrid logical grid       
!     
!       do i = 1,xsize
!         do j = 1, ysize
!            if(fice(i,j) .gt. 0 .and. fice(i,j).le. 100) then
!                icegrid(i,j) = .true.
!            else
!                icegrid(i,j) = .false.
!            endif
!         enddo
!       enddo

       deallocate (fsst, fice)

      return
      end subroutine read_reynolds_sst

!------------------------------------------------------------------------------
      subroutine assign_sst

      implicit none

      real               :: alat, alon, rincs
      real               :: scl, sst1, sst2
      real               :: slat, slat1, slat2
      real               :: slon, slon1, slon2
      integer            :: li, pi, isst1=0, isst2=0
      integer(2)         :: ilat, ilon, islat, islon, smaskres
      integer(2)         :: islat1, islat2, islon1, islon2
      integer(2)         :: i,j,ixsize,iysize


      smaskres = 4           ! ppd  (4 = 0.25 X 0.25) 
      rincs    = (1. / float(smaskres))  * 0.5 !1/2 of grid increment

      do lin = 1, nscans
        do pix = 1, npix
                !write(*,*)'sfc_type(pix,lin)=',sfc_type(pix,lin)
          !if(sfc_type(pix,lin) .ne. 1) then
                 sst(pix,lin) = miss_flt   ! 20=land, 30=coast, ice handled below
          !       write(*,*)'sst(pix,lin)=',sst(pix,lin)
          !       cycle
          !endif
          alat = lat(pix,lin)            !input lat
          alon = lon(pix,lin)            !input lon centered on 0 deg
          alon=alon + 180.0              !shift for grid counter 

          islat = nint(smaskres*(90.+rincs-alat)) !assign sst array elements
          slat = smaskres*(90.+rincs-alat)
          !if(islat .eq. 181*smaskres) islat = 180*smaskres !chk
          !forsouthpole
          islon = nint(smaskres*(rincs+alon))
          slon = smaskres*(rincs+alon)

          if(islon .eq. (smaskres*360+1)) islon = 1
          if ((islat .lt. 1) .or. (islat .gt. smaskres*180) .or. !boundary check 
     >        (islon .lt. 1) .or. (islon .gt. smaskres*360))  then
              sst(pix,lin) = miss_flt
              cycle
          else
              if(sstgrid(islon,islat) .lt. -1.8) then ! why this cutoff?
                  sst(pix,lin) = miss_flt
                  !sfc_type(pix,lin) = 0 !31  !assign
                  cycle
              else
                  slat1= slat - islat ! lat index fraction -.5 to .5
                  slon1= slon - islon ! lon index fraction -.5 to .5
                  if(slat1.eq.0.0.and.slon1.eq.0.0) then
                    sst(pix,lin) = sstgrid(islon,islat) + 273.15
                    cycle
                  endif
                  if(minval(sstgrid(islon-1:islon+1,islat-1:islat+1))
     >               .lt.-1.8 .or. islat.eq.1.or.islat.eq.720 .or.
     >               islon.eq.1.or.islon.eq.1440) then
                    sst(pix,lin) = sstgrid(islon,islat) + 273.15
                    cycle
                  endif

                  if(slat1.gt.0)
     >             slat2=sstgrid(islon,islat+1)-sstgrid(islon,islat)
                  if(slat1.lt.0)
     >             slat2=sstgrid(islon,islat-1)-sstgrid(islon,islat)
                  if(slon1.gt.0)
     >             slon2=sstgrid(islon+1,islat)-sstgrid(islon,islat)
                  if(slon1.lt.0)
     >             slon2=sstgrid(islon-1,islat)-sstgrid(islon,islat)

                  sst(pix,lin) = 273.15 + sstgrid(islon,islat) +
     >                           slat2*abs(slat1) + slon2*abs(slon1)
              endif
          endif

        enddo   !pi
      enddo     !li

      deallocate (sstgrid)       !free up the sstgrid

      return
      end subroutine assign_sst

      subroutine prep_oe
        implicit none

        real :: freq,temp,pres,rho
        integer :: bin,era_lun1,era_lun2,era_lun3,rsslun,ios,syl,i,j,k
        character(len=2) :: charnch 
        character(len=3) :: ver='vW' ! LUT version
        character(len=4) :: yr
        character(len=2) :: mos(12) = (/'1','2','3','4','5','6',
     >                                '7','8','9','10','11','12'/)
        character(len=2) :: das(31) = (/'1','2','3','4','5','6',
     >     '7','8','9','10','11','12','13','14','15','16','17',
     >   '18','19','20','21','22','23','24','25','26','27','28',
     >   '29','30','31'/)
        character(len=2) :: mo, da
        character(len=3) :: emis
        character(len=1) :: spc
        character(len=100):: sydir
        !sydir = '/tdata1/dduncan/oe/sy/'
        sydir = 'binary/'

      !allocate(oe_output(npix,nscans,9),screen(npix,nscans)) !last dimension hard-coded!
      !allocate(Tb_diff(npix,nscans,maxchans),save_slp(npix,nscans)) 
      !allocate(save_wdir(npix,nscans),save_sal(npix,nscans)) 
      !allocate(save_iter(npix,nscans),poster(npix,nscans,8)) 
      !allocate(toffsets(nch),s_toffsets(nbins,nch),s_sy(nch,nch,nbins))
      !allocate(sy(nch,nch),sy_i(nch,nch)) ! local currently, with nch in name.

      if(nch.ge.10) write(charnch,'(I2)'), nch !
      if(nch.le.9)  write(charnch,'(I1)'), nch !
      !charnch = '9' ! OVERRIDE (for now)
      write(spc,'(I1)') npc
        spc = '3' ! for testing -- CHANGE
      call gprof_lun(syl)
!      changed Sy reading 7/21/15 DD, again 8/27, again 10/04 (just
!      directory call instead of within current dir)

      emis = 'F' ! default
      !if(run_rss) emis = 'R'
      if(useanalysis) emis = trim(emis) // 'a'
      if(useanalysis.ne..true.) emis = trim(emis) // 'c' !climatology
      if(nretvar.gt.5) emis=trim(emis)//'s'
        !emis='Fc' !OVERRIDE FOR TESTING!
        
      sy_file = trim(sydir)//satcode// '.SY.'//trim(emis)//
     >  '.'//trim(spc)//'.'//trim(charnch)//'.v2.1.bin'
      open(unit=syl,file=trim(sy_file),
     >  access='stream',status='old',form='unformatted',iostat=ios)
       if(ios .ne. 0) write(*,*)' cant open sy_file',sy_file
        read(syl) sy
        read(syl) toffsets
        read(syl) s_sy
        read(syl) s_toffsets
      close(syl)

      mr_sigmas(:,1) = (/1.33,1.27,1.37,1.45,1.47,1.44,1.40,1.34,1.30,
     > 1.25,1.20,1.15,1.15,1.14,1.11,1.11,1.12,1.12,1.10,1.07,1.08,
     > 1.05,1.02,1.00,0.98,1.03,1.05,1.09,1.20,1.20,1.20,1.20,1.20/)
      mr_sigmas(:,2) = (/0.76,0.87,0.95,0.94,0.93,0.87,0.88,0.91,0.91,
     > 0.89,0.86,0.85,0.88,0.92,0.93,0.95,0.95,0.96,0.94,0.91,0.86,
     > 0.79,0.74,0.68,0.63,0.60,0.60,0.66,0.66,0.66,0.66,0.66,0.66/)
      mr_sigmas(:,3) = 0.40

      call gprof_lun(syl)
      open(unit=syl,file= 'binary/SA-EOF_LW.3.v2.0c.bin',
     >  access='stream',status='old',form='unformatted',iostat=ios)
       if(ios .ne. 0) write(*,*)' cant open eof-lw cov file'
        read(syl) eof_lw_cov ! 3 x nbins
      close(syl)
     
!--- allocate and read in EC LUT data
      allocate(era_wm(nlon,nlat), era_ws(nlon,nlat))
      allocate(mprof(nbins,nz),eofs(nbins,nz,6))!npc))
      allocate(peofs(nz,6))

      write(yr,'(I4)') stdtime(5,1)
      mo = mos(stdtime(5,2)) !use as index
      da = das(stdtime(5,3)) !use as index

      erafile = 'binary/eof_mr.03.v3.bin' ! v3 has consistency across bins
      call gprof_lun(era_lun1)
      open(unit=era_lun1,file=erafile,access='stream',
     >     status='old',form='unformatted',iostat = ios)
       if(ios .ne. 0) write(*,*)' cant open EOF LUT file'
        read(era_lun1) mprof
        read(era_lun1) eofs
      close(era_lun1)

      !erafile ='/tdata1/dduncan/oe/LUT/ws_n128x4_vC.'//
      erafile ='binary/ws_n128x4_vC.'//
     >    trim(mo)//'.bin'
       call gprof_lun(era_lun1)
      open(unit=era_lun1,file=erafile,access='stream',
     >     status='old',form='unformatted',iostat = ios)
       if(ios .ne. 0) write(*,*)' cant open EC winds file'
        read(era_lun1) era_wm
        read(era_lun1) era_ws
      close(era_lun1)

      call inverse(sy,sy_i,nch) ! new to prep_oe!!!!


      end subroutine prep_oe

      subroutine qual_czech
       !use define_oeorbit
       implicit none
!----- new 11/17/15, D Duncan
!      this will screen out anomalously high TPW
!       pixels in which it is likely raining lightly or LWP
!       is being traded off for more WV.

       integer :: lin,pix,ct,s1,p1
       integer,parameter :: sd=3,pd=5
       real    :: sigmatp,meantp
       real    :: meth=5.0,sith=6.5,tottp,eachtp((2*pd+1)*(2*sd+1))
       real    :: temptp(2*pd+1,2*sd+1), tempch(2*pd+1,2*sd+1)
       
       screen(:,:) = oe_output(:,:,1)
       do lin = 1+sd,nscans-sd
         do pix = 1+pd,npix-pd
           if(oe_output(pix,lin,1).lt.0) cycle
           temptp(:,:) = oe_output((pix-pd):(pix+pd),
     >       (lin-sd):(lin+sd),1)
           tempch(:,:) = oe_output((pix-pd):(pix+pd),
     >       (lin-sd):(lin+sd),5)
           ct=0
           tottp=0.0
           do s1=1,2*sd+1
             do p1=1,2*pd+1
               if(temptp(p1,s1).gt.0.and.tempch(p1,s1).lt.chis1
     >            .and.temptp(p1,s1).ne.oe_output(pix,lin,1)) then
                 ct=ct+1
                 tottp=tottp+temptp(p1,s1)
                 eachtp(ct) = temptp(p1,s1)
               endif
             enddo
           enddo
           if(ct.lt.15) cycle ! need enough to get real background sense
           meantp = tottp/real(ct)
           sigmatp = SQRT(SUM((eachtp(1:ct)-meantp)**2)/real(ct-1))
           if(oe_output(pix,lin,1).gt.(meantp+meth) .or.
     >        sigmatp.gt.sith) then
             screen(pix,lin) = -996
           endif 
         enddo
       enddo
       oe_output(:,:,1) = screen(:,:)

      end subroutine qual_czech

!----------------------------------------------------------------------
      subroutine gprof_lun(ilun)
!    
!** This routine gets an open logical unit number starting with 100
!
       integer  :: ilun
       logical  :: llun
!       
       do ilun = 100,201
         inquire(unit=ilun, opened=llun)
         if(.not. llun) exit
       enddo
       return
      end subroutine gprof_lun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
        implicit none 
        integer n
        real a(n,n), c(n,n)
        real*8 aa(n,n),cc(n,n)
        real*8 L(n,n), U(n,n), b(n), d(n), x(n)
        real*8 coeff
        integer i, j, k

        ! step 0: initialization for matrices L and U and b
        ! Fortran 90/95 aloows such operations on matrices
        L=0.0
        U=0.0
        b=0.0

        aa(:,:) = dble(a(:,:)) ! do double precision!
        ! step 1: forward elimination
        do k=1, n-1
           do i=k+1,n
              coeff=aa(i,k)/aa(k,k)
              L(i,k) = coeff
              do j=k+1,n
                 aa(i,j) = aa(i,j)-coeff*aa(k,j)
              end do
           end do
        end do

        ! Step 2: prepare L and U matrices 
        ! L matrix is a matrix of the elimination coefficient
        ! + the diagonal elements are 1.0
        do i=1,n
          L(i,i) = 1.0
        end do
        ! U matrix is the upper triangular part of A
        do j=1,n
          do i=1,j
            U(i,j) = aa(i,j)
          end do
        end do

        ! Step 3: compute columns of the inverse matrix C
        do k=1,n
          b(k)=1.0
          d(1) = b(1)
        ! Step 3a: Solve Ld=b using the forward substitution
          do i=2,n
            d(i)=b(i)
            do j=1,i-1
              d(i) = d(i) - L(i,j)*d(j)
            end do
          end do
        ! Step 3b: Solve Ux=d using the back substitution
          x(n)=d(n)/U(n,n)
          do i = n-1,1,-1
            x(i) = d(i)
            do j=n,i+1,-1
              x(i)=x(i)-U(i,j)*x(j)
            end do
            x(i) = x(i)/u(i,i)
          end do
        ! Step 3c: fill the solutions x(n) into column k of C
          do i=1,n
            cc(i,k) = x(i)
          end do
          b(k)=0.0
        end do
        c(:,:) = real(cc(:,:))
      end subroutine inverse

      subroutine import_sal
      ! opens/reads binary file of Aquarius monthly SSS data from 2014, 
      !  (smoothed and formatted offline) to 720x360 grid, starting at
      !  0E,90N with mis_val=-999.9, type real
      ! D Duncan, CSU, 2/18/16

        integer :: lun,ios
        character(len=80) sfile
        character(len=180) dir
        character(len=2) :: money
        character(len=2) :: mons(12) = (/'01','02','03','04','05','06',
     >                                   '07','08','09','10','11','12'/)

        money = mons(stdtime(500,2)) !use as index
        dir = 'binary/' 
        !dir = '/tdata1/dduncan/Aquarius_SSS/'
        sfile = trim(dir) // 'sss.'//money//'.v1.bin'

        call gprof_lun(lun)
        open(unit=lun,file=trim(sfile),access='stream',
     >     status='old',form='unformatted',iostat = ios)
         if(ios .ne. 0) write(*,*)' cant open sss file'
          read(lun) sss
        close(lun)

      end subroutine import_sal

      ! new, lwmax read routine
      subroutine read_lwmax
        integer :: ios,lun
        call gprof_lun(lun)
        open(unit=lun,file='binary/lwmaxi.bin',
     >    access='stream',status='old',form='unformatted',iostat=ios)
         if(ios .ne. 0) write(*,*)' cant open lwmax file'
          read(lun) lwmax
        close(lun)
      end subroutine read_lwmax
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! new,  read routine to read in SA and EOFs for drizz
      subroutine read_drizz

        integer :: ios,luna,lune,lunr,syl1,syl2,syl3,syl4
        character(len=100) :: rosy1,rosy2,rosy3,rosy4
        call gprof_lun(lune)

        open(unit=lune,file='binary/eofs.20145.RPr.v1.bin',
     >    access='stream',status='old',form='unformatted',iostat=ios)
         if(ios .ne. 0) write(*,*)' cant open drizz eof file'
          read(lune) rwpm   !11x14
          read(lune) piwpm  !11x14
          read(lune) rwpe   !11x14x2
          read(lune) piwpe  !11x14x2
        close(lune)

        call gprof_lun(luna)
        open(unit=luna,file='binary/sa.20145.RPr.bin',
     >    access='stream',status='old',form='unformatted',iostat=ios)
         if(ios .ne. 0) write(*,*)' cant open drizz sa file'
          read(luna) drizsa   !3x3x11
        close(luna)

        call gprof_lun(lunr)
        open(unit=lunr,file='binary/rrcoeff.linab.bin',
     >    access='stream',status='old',form='unformatted',iostat=ios)
         if(ios .ne. 0) write(*,*)' cant open rr coeff file'
          read(lunr) rcoef   !11x2 (a+bx coefficients)
        close(lunr)

! new, read in RWP Sy values
!   i.e. Sy(RWP) = sign*(sqrt(abs(Sy_m))+(RWP-RWP_ref)*dif)^2  where RWP is in kg
      call gprof_lun(syl1)
      rosy1 = 'binary/sy.281.8.150.L4.conv.bin'
      open(unit=syl1,file=trim(rosy1),
     >  access='stream',status='old',form='unformatted',iostat=ios)
       if(ios .ne. 0) write(*,*)' cant open sy_rc_file',trim(rosy1)
        read(syl1) rc_rlo ! anchor points for interp
        read(syl1) rc_rme
        read(syl1) rc_rhi
        !print*,'rch:',rc_rlo,rc_rme,rc_rhi
        read(syl1) rc_syl ! Sy at rlo
        read(syl1) rc_sym ! Sy at rme
        read(syl1) rc_syh ! Sy at rhi
        read(syl1) rc_sgm ! sigma values, dblarr(nch,3)
        read(syl1) rc_r !r values, dblarr(nch,nch,3)
      close(syl1)
      call gprof_lun(syl2)
      rosy2 = 'binary/sy.281.8.150.L4.stra.bin'
      open(unit=syl2,file=trim(rosy2),
     >  access='stream',status='old',form='unformatted',iostat=ios)
       if(ios .ne. 0) write(*,*)' cant open sy_rs_file',trim(rosy2)
        read(syl2) rs_rlo ! anchor points for interp
        read(syl2) rs_rme
        read(syl2) rs_rhi
        !print*,'rsh:',rs_rlo,rs_rme,rs_rhi
        read(syl2) rs_syl ! Sy at rlo
        read(syl2) rs_sym ! Sy at rme
        read(syl2) rs_syh ! Sy at rhi
        read(syl2) rs_sgm ! sigma values, dblarr(nch,3)
        read(syl2) rs_r   ! r values, dblarr(nch,nch,3)
      close(syl2)
      call gprof_lun(syl2)
      rosy3 = 'binary/sy.281.8.150.L4l.conv.bin'
      open(unit=syl3,file=trim(rosy3),
     >  access='stream',status='old',form='unformatted',iostat=ios)
       if(ios .ne. 0) write(*,*)' cant open sy_rcl_file',trim(rosy1)
        read(syl3) rcl_rlo ! anchor points for interp
        read(syl3) rcl_rme
        read(syl3) rcl_rhi
        !print*,'rcl:',rcl_rlo,rcl_rme,rcl_rhi
        read(syl3) rcl_syl ! Sy at rlo
        read(syl3) rcl_sym ! Sy at rme
        read(syl3) rcl_syh ! Sy at rhi
        read(syl3) rcl_sgm ! sigma values, dblarr(nch,3)
        read(syl3) rcl_r !r values, dblarr(nch,nch,3)
      close(syl3)
      call gprof_lun(syl4)
      rosy4 = 'binary/sy.281.8.150.L4l.stra.bin'
      open(unit=syl4,file=trim(rosy4),
     >  access='stream',status='old',form='unformatted',iostat=ios)
       if(ios .ne. 0) write(*,*)' cant open sy_rsl_file',trim(rosy4)
        read(syl4) rsl_rlo ! anchor points for interp
        read(syl4) rsl_rme
        read(syl4) rsl_rhi
        !print*,'rsl:',rsl_rlo,rsl_rme,rsl_rhi
        read(syl4) rsl_syl ! Sy at rlo
        read(syl4) rsl_sym ! Sy at rme
        read(syl4) rsl_syh ! Sy at rhi
        read(syl4) rsl_sgm ! sigma values, dblarr(nch,3)
        read(syl4) rsl_r   ! r values, dblarr(nch,nch,3)
      close(syl4)
      ! scale to right Sy matrix and invert within csu1dvar.f

      end subroutine read_drizz
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine calc_az(slat,ipi,ascdes,sc_orient,rot_ang)

!       calculating azimuthal angle of pixel, 
!        relative to due North, for GMI/TMI's footprint

!      slat = spacecraft lat
!      ipi  = position in scan (pixel number, 1-221)   
!      ascdes= ascending/descending pass (1=asc, -1=desc)    
!      sc_orient = spacecraft orientation (180=backward, 0=forward)

      real :: slat!,sc_orient
      integer :: nfovhigh, iatt, ascdes, ipi!, badorient=0
      real*8 :: satincl != 65.00000000 ! satellite inclination, relative to E
      real*8 :: scanangle != 140.00
      real*8 :: dpi, d2r, theta, gtrack, sc_orient
      real   :: rot_ang

      if(satcode.eq.'GMI') then
        nfovhigh  = 221 ! # FOVs in high freq chans, thus # pix in L1C/L1CR
        satincl   = 65.000000
        scanangle = 140.00
      endif
      if(satcode.eq.'TMI') then
        nfovhigh  = 104 ! # FOVs in high freq chans, thus # pix in L1C/L1CR
        satincl   = 35.000000
        scanangle = 130.00
      endif

      iatt = 0
      dpi = 4.0d+0*atan(1.0) ! define Pi
      d2r = dpi / 180.0d+0   ! degrees to radians

        ! take into account direction of spacecraft -- rear/forward
      if(sc_orient.eq.180.0) iatt=-1
      if(sc_orient.eq.0.0)   iatt=1
      !if(iatt.eq.0) then
      !  badorient = badorient+1 ! counter for pixels with orientation
        ! not forward or backward, if wanting to output later to file?
        ! algorithm will skip this
        ! pixel due to missing Tbs anyway, so let it go...
      !endif

        ! make sure sat_lat works in below code
      if(slat.gt.satincl) slat = satincl
      if(slat.lt.(-1.0*satincl)) slat = -1.0*satincl

        ! contribution to azimuthal angle by pixel position
      theta = ((-0.5d+0*SCANANGLE) 
     >   +(SCANANGLE*(ipi-1))/(NFOVHIGH-1.0d+0))
        ! which is different if spacecraft moving backwards!
      if(iatt .eq. -1) then
        theta=((-0.5d+0*scanangle -180.0)+
     >   (scanangle*(ipi-1))/(nfovhigh-1.0d+0))
      endif

        ! contribution of lat on azimuthal angle
      gtrack=ascdes*(90.0d+0-
     >   asin(cos(satincl*d2r)/cos(abs(slat)*d2r))/d2r)

        ! now find angle between footprint long axis (from sat dir)
        !  and North, clockwise. gotta rotate 90 degrees, too.
      rot_ang = (90.0d+0 - gtrack - theta)
      if(rot_ang .ge. 360.0) rot_ang = rot_ang-360.0

      end subroutine calc_az


      end module csu1dvar_subr
