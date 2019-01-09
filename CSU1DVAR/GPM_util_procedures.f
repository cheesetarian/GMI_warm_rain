      module GPM_util_procedures
!--- GPM utility module, modified not at all apart from use declaration below,
!dduncan csu 1/24/17
       use GPM_arraydef_D
       use define_intam2
      
       contains

!-----------------------------------------------------------------------

      subroutine GPM_calc_EIA_offset
       implicit none

!---  internal variables
       integer :: lin, pix, ch
      
       allocate (tbdelta_EIA(npixs,nscans,15))
       tbdelta_eia = 0.0

c     write(log_lun,*)'chanavail = ', chan_avail
            
       do lin = 1, nscans    
         do pix = 1, npixs 
           do ch = 1,15           
            if(chan_avail(ch) .and. eia(pix,lin,ch).ne.miss_flt) then
                  tbdelta_eia(pix,lin,ch) = db_nominal_eia(ch) - 
     >                                         eia(pix,lin,ch)
             endif
           enddo          
         enddo    !pix
       enddo    !lin
      
       return
       end subroutine GPM_calc_EIA_offset

!---------------------------------------------------------------------------

       subroutine GPM_output
       implicit none
       
       logical            :: backsearch=.true.
       
       integer            :: ou_lun, ios, lin, pix, ilyr, reccnt=0
       integer            :: begi, i
       integer            :: temptime(8)
       character(len=5)   :: corbitn
       character(len=128) :: blank = ' '
       
       character(len=20):: SpeciesDesc(nspecies)=
     >    ( /'Rain Water Content','Cloud Water Content',
     >       'Ice Water Content', 'Snow Water Content',
     >       'Grauple/Hail Content'/)  
       real  :: HghtTopLayers(nlyrs)= 
     >           (/.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.,
     >            5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10.,
     >            11., 12., 13., 14., 15., 16., 17., 18./)
       real :: TempDesc(ntemps)=(/270,273,276,279,282,285,288,291,294,
     >                              297,300,303/)
       
!--- create output structures
       
       type :: Date6         
          integer(kind=knd2):: year         !for file creation dates
          integer(kind=knd2):: month
          integer(kind=knd2):: day
          integer(kind=knd2):: hour
          integer(kind=knd2):: minute
          integer(kind=knd2):: second
       end type Date6
       
       type :: Date7         
          integer(kind=knd2):: year        !scan date includes millisecs
          integer(kind=knd2):: month
          integer(kind=knd2):: day
          integer(kind=knd2):: hour
          integer(kind=knd2):: minute
          integer(kind=knd2):: second
          integer(kind=knd2):: millisec
       end type Date7
       
       type :: OrbitHdr                     !400 bytes  per file
          character(len=12) :: Satellite
          character(len=12) :: Sensor
          character(len=12) :: PreProcessorVersion
          character(len=12) :: AlgorithmVersion          
          character(len=128):: ProfileDatabaseFile
          character(len=128):: RadiometerFile
          type(Date6)       :: FileCreationDate
          type(Date6)       :: GranuleStartDate
          type(Date6)       :: GranuleEndDate
          integer           :: GranuleNumber
          integer(kind=knd2):: NumScansGranule
          integer(kind=knd2):: NumPixelsScan
          integer(kind=knd1):: ProfStructFlag
          character(len=51) :: Spare51
       end type OrbitHdr
       
       type :: ProfileInfo                  !537,864 bytes  per file
          integer(kind=knd1) :: Nspecies        
          integer(kind=knd1) :: Ntemps 
          integer(kind=knd1) :: Nlyrs 
          integer(kind=knd1) :: Nprfs 
          character(len=20)  :: SpeciesDesc(nspecies)    
          real               :: HghtTopLayers(nlyrs)     
          real               :: TempDesc(ntemps)      
          real  :: ProfileClusters(nspecies,ntemps,nlyrs,nprfs)
           
       end type ProfileInfo
      
       type :: ScanHdr                     ! 28 bytes in ScnHdr per scan
          real               :: Sclat
          real               :: Sclon
          real               :: Scalt
          type(Date7)        :: ScanDate
          integer(kind=knd2) :: Spare
       end type ScanHdr

       type :: DataRec                    ! 84 bytes in the Datarec/ pixel
          integer(kind=knd1) :: PixelStatus
          integer(kind=knd1) :: QualityFlag
          integer(kind=knd1) :: L1CQualflag
          integer(kind=knd1) :: SurfaceTypeIndex
          integer(kind=knd1) :: TotalColWaterVaporIndex  !model tcwvi
          integer(kind=knd1) :: ProbabilityofPrecip
          integer(kind=knd2) :: Temp2MeterIndex          !model T2mi
          
	  integer(kind=knd2) :: CAPE                     !model derived CAPE
          integer(kind=knd1) :: SunglintAngle          
	  integer(kind=knd1) :: Spare1
          
          real               :: Latitude
          real               :: Longitude
          real               :: SurfacePrcp
          real               :: FrozenPrcp
          real               :: ConvectivePrcp

          real               :: RainWaterPath  
          real               :: CloudWaterPath
          real               :: IceWaterPath
                    
          real               :: MostLikelyPrcp
          real               :: Prcp1stTertial
          real               :: Prcp2ndTertial
	  	  
	  integer(kind=knd2) :: ProfileTemp2mIndex         !2 bytes
          integer(kind=knd2) :: ProfileNumber(nspecies)    !2 bytes*5 = 10
          real               :: ProfileScale(nspecies)     !4 bytes*5 = 20

       end type DataRec
         
!--- assign structure names
       
       type(OrbitHdr)    :: orbhdr
       type(ScanHdr)     :: scnhdr
       type(ProfileInfo) :: prfinfo
       type(DataRec)     :: rec
       
!---open output file

       call GPM_lun(ou_lun)
       open(unit=ou_lun,file=trim(output_file),
     >      access='stream',status='unknown',iostat = ios)  
       if(ios .ne. 0) then
           call GPM_reprt_err(30,ou_lun,trim(output_file))    !opening output file
       endif
D      write(log_lun,*)'   output file                  : ',
D    >                 trim(output_file)
       
!-- Fill orbhdr with program variables
       
       orbhdr%Satellite           = sat_name
       orbhdr%Sensor              = sensor_name
       orbhdr%PreProcessorVersion = trim(pp_version)
       orbhdr%AlgorithmVersion    = trim(alg_version)
       orbhdr%ProfileDatabaseFile = trim(pdbfile)
       orbhdr%RadiometerFile      = trim(radfile)
              
       begi = scan(radfile,'/',backsearch) + 1
       orbhdr%RadiometerFile = trim(radfile(begi:256))

       call date_and_time(values=temptime)          
       orbhdr%FileCreationDate%year   = temptime(1)
       orbhdr%FileCreationDate%month  = temptime(2)
       orbhdr%FileCreationDate%day    = temptime(3)
       orbhdr%FileCreationDate%hour   = temptime(5)
       orbhdr%FileCreationDate%minute = temptime(6)
       orbhdr%FileCreationDate%second = temptime(7)

       do i = 1,nscans                                               !chk for 1st good scans             
          if(stdtime(i,1).ge.1987 .and. stdtime(i,1).le. 2050 .and.
     >       stdtime(i,2).ge.1    .and. stdtime(i,2).le. 12   .and.
     >       stdtime(i,3).ge.1    .and. stdtime(i,3).le. 31)  then
                exit
          endif
       enddo
       if(i .gt. nscans) call GPM_reprt_err(31,i,blank)         !finding a good time

       orbhdr%GranuleStartDate%year   = stdtime(i,1)
       orbhdr%GranuleStartDate%month  = stdtime(i,2)
       orbhdr%GranuleStartDate%day    = stdtime(i,3)
       orbhdr%GranuleStartDate%hour   = stdtime(i,4)
       orbhdr%GranuleStartDate%minute = stdtime(i,5)
       orbhdr%GranuleStartDate%second = stdtime(i,6)
       
       orbhdr%GranuleEndDate%year     = stdtime(nscans,1)
       orbhdr%GranuleEndDate%month    = stdtime(nscans,2)
       orbhdr%GranuleEndDate%day      = stdtime(nscans,3)
       orbhdr%GranuleEndDate%hour     = stdtime(nscans,4)
       orbhdr%GranuleEndDate%minute   = stdtime(nscans,5)
       orbhdr%GranuleEndDate%second   = stdtime(nscans,6)
       
       orbhdr%GranuleNumber           = granule
       orbhdr%NumScansGranule         = nscans
       orbhdr%NumPixelsScan           = npixs
       orbhdr%ProfStructFlag          = profstructflag
       orbhdr%spare51(1:40)            = calfile(1:40)

!--- write out orbit header info and cluster info

       write(ou_lun, iostat=ios)  orbhdr
       if(ios .ne. 0) call GPM_reprt_err(32,ou_lun,trim(output_file))     !writing orb hdr
       
!--- write profile information
      
       prfinfo%nspecies        = nspecies 
       prfinfo%ntemps          = ntemps 
       prfinfo%nlyrs           = nlyrs 
       prfinfo%nprfs           = nprfs
       prfinfo%SpeciesDesc     = speciesdesc
       prfinfo%HghtTopLayers   = hghtToplayers
       prfinfo%TempDesc        = tempdesc
       
       prfinfo%ProfileClusters = miss_flt
              
       prfinfo%ProfileClusters(1:5,1:12,1:28, 1:40) = profileclusters    !raining in 1-40
       prfinfo%ProfileClusters(2:3,1:12,1:28,41:80) = profileclustersNR  !non raining in 41-80
       
       write(ou_lun,iostat=ios) prfinfo
       if(ios .ne. 0) call GPM_reprt_err(33,ou_lun,trim(output_file))  !err writing prfinfo
       
!--- loop over all scans

       do lin = 1, nscans                            !outputs all scans in file       
         scnhdr%sclat = sclat(lin)
         scnhdr%sclon = sclon(lin)
         scnhdr%scalt = scalt(lin)
         scnhdr%scandate%year     = stdtime(lin,1)
         scnhdr%scandate%month    = stdtime(lin,2) 
         scnhdr%scandate%day      = stdtime(lin,3)
         scnhdr%scandate%hour     = stdtime(lin,4)
         scnhdr%scandate%minute   = stdtime(lin,5)
         scnhdr%scandate%second   = stdtime(lin,6)
         scnhdr%scandate%millisec = 0              !no millisec information

         write(ou_lun, iostat=ios)  scnhdr            !write out scan header
         if(ios .ne. 0)  call GPM_reprt_err(34,ou_lun,trim(output_file)) !err writing scnhdr

!---   loop over all pixels
      
         rec%SurfacePrcp = 0        
         do pix = 1, npixs         
           reccnt = reccnt + 1
            
!---      Set status and index variables
                      
           rec%PixelStatus = pixel_status(pix,lin)

!---      Set L1Cqual flag

           rec%L1CQualflag =  L1CQualflag (pix,lin)

!---      Set quality Flag, 0 = normal, 1 = use with caution, 2 = use with extreme caution

           if(quality_flag(pix,lin) .eq. 0) then                                 !qual flag not set yet, also set in GPM_rain
              if(rec%L1Cqualflag .gt. 0)	   quality_flag(pix,lin) = 1	 !sglint,RFI,geolocate,warmload
              if(abs(rec%SunglintAngle).le.20) quality_flag(pix,lin) = 1	 !sunglint contamination
              if(sfccode(pix,lin) .ge.  8 .and.				   !all snow surfaces 
     >           sfccode(pix,lin) .le. 11)     quality_flag(pix,lin) = 1	 ! set to 1
              if(rec%L1Cqualflag .le. -125)    quality_flag(pix,lin) = 2	 !channels are missing unexpectedly
              if(sfcprcp(pix,lin).eq.miss_flt) quality_flag(pix,lin) = 2	 !surface rain is missing
	   endif
           rec%QualityFlag = quality_flag(pix,lin)

!---      more indexes

           rec%SurfaceTypeIndex = sfccode(pix,lin)

	   if(pptcwv(pix,lin) .eq. miss_flt) then
	       rec%TotalColWaterVaporIndex = miss_byt
	   else
               rec%TotalColWaterVaporIndex = nint(pptcwv(pix,lin))
	   endif

!---      Insert the precip probability

           rec%ProbabilityofPrecip  = nint(pop(pix,lin)*100)
           if(rec%ProbabilityofPrecip .gt. 100) then
              rec%ProbabilityofPrecip = 100
           endif
           
!---      other parameters
	   
	   if(ppT2m(pix,lin) .eq. miss_flt) then
	       rec%Temp2meterIndex    =  miss_int2
	   else
	       rec%Temp2meterIndex    = nint(ppT2m(pix,lin))
           endif
	   rec%CAPE                   = CAPE(pix,lin)                      
           rec%SunglintAngle          = sglinta(pix,lin)

!---      Set Latitude and Longitude
                                                                          
           rec%Latitude    = lat(pix,lin)
           rec%Longitude   = lon(pix,lin)

!---      Set Precipitation variables

           rec%SurfacePrcp    = sfcprcp(pix,lin)           
           rec%FrozenPrcp     = frzprcp(pix,lin)
           rec%ConvectivePrcp = cnvprcp(pix,lin)

!---      Set integrated Path variables

           rec%RainWaterPath          = rwp(pix,lin)
           rec%CloudWaterPath         = cwp(pix,lin)
           rec%IceWaterPath           = iwp(pix,lin)

!---      Set Precipitation Diagnostics
           
           rec%MostLikelyPrcp         = mlprcp(pix,lin)
           rec%Prcp1stTertial         = prcp1stT(pix,lin)
           rec%Prcp2ndTertial         = prcp2ndT(pix,lin)
           
!---      Set profile declaration variables
                      
           rec%ProfileTemp2mIndex     = prfT2mindex(pix,lin)
           rec%ProfileNumber          = prfnum(pix,lin,:)
           rec%ProfileScale           = prfscale(pix,lin,:)

!---      Write out record
           
           write(ou_lun, iostat=ios)  rec                                  !write pixel structure
           if(ios .ne. 0)call GPM_reprt_err(35,ou_lun,trim(output_file))  !err writing pix%rec
           	   
	          
         enddo  !pix
       enddo   !lin
                 
          
       close(ou_lun) 
D      write(log_lun,*)'   number of records written    : ',reccnt 
       return
       end subroutine GPM_output
 
 !--------------------------------------------------------------------------
                
      subroutine GPM_reprt_err(err,istat,cstat)
       implicit none
      
      
!-- this subroutine reports errors and writes calls output_error which
!-- will write out a file with only header information including error comment.
                    
      integer  :: err,istat
      character(len=128):: cstat
      character(len=31) :: fe=' FATAL ERROR code:'
      character(len=51) :: errmess
      
 8    format(a19,i3,a)         !output format for termination error lines         
 9    format(a19,i3,a,a)         !output format for termination error lines       
 10   format(a19,i3,a,i)         !output format for termination error lines         
      write(log_lun,*)     

!--- select error message depending on err code

      select case(err)
      case(1)
        errmess = ' Improper command line arguments'
        write(*,8)fe,err,errmess
      case(2)
        errmess = ' Opening Log file= '
        write(*,9) fe,err,errmess,trim(cstat)
      case(3)
        errmess = ' Retrieving Logical Unit '
        write(*,8) fe,err,errmess
      case(4)
        errmess = ' Fatal error from Preprocessor: '
        write(log_lun,8) fe,err,errmess
      case(10)
        errmess = ' Opening PP file '
        write(log_lun,9)fe,err,errmess,trim(cstat)  
      case(11)
        errmess = ' Reading PP file Orbit header'
        write(log_lun,8)fe,err,errmess     
      case(13)
        errmess = ' Reading PP scan header, scan='
        write(log_lun,10)fe,err,errmess,istat  
      case(14)
        errmess = ' Reading PP pixel record, pixel='
        write(log_lun,10)fe,err,errmess,istat  
      case(15)
        errmess = ' No pixels in requested area'
        write(log_lun,8)fe,err,errmess

      case(16)
        errmess = ' No Tbs passing quality control'
        write(log_lun,8)fe,err,errmess       
      case(17)
        errmess = ' No bins selected in Histogram'
        write(log_lun,8)fe,err,errmess        
      
      case(20)
        errmess = ' Opening Non-raining Profile Clusters'
        write(log_lun,9)fe,err,errmess,trim(cstat)        
      case(21)
        errmess = ' Reading Non-raining Profile Clusters'
        write(log_lun,9)fe,err,errmess,trim(cstat)  
      case(22)
        errmess = ' Opening Species Profile Clusters'
        write(log_lun,9)fe,err,errmess,trim(cstat)        
      case(23)
        errmess = ' Reading Species Profile Clusters'
        write(log_lun,9)fe,err,errmess,trim(cstat)

      case(24)
        errmess = ' Opening Wet Bulb rain/snow file'
        write(log_lun,9)fe,err,errmess,trim(cstat)
      case(25)
        errmess = ' Reading Wet Bulb rain/snow file'
        write(log_lun,9)fe,err,errmess,trim(cstat)
      
      case(27)
              errmess = ' Opening Channel sensititivy file: '
        write(log_lun,9)fe,err,errmess,trim(cstat)
      case(28)  
              errmess = ' Reading Channel sensititivy file: '
        write(log_lun,9)fe,err,errmess,trim(cstat)
      
      case(30)
        errmess = ' Opening output file='
        write(log_lun,9)fe,err,errmess,trim(cstat)          
      case(31)
        errmess = ' Finding 1st good time'
        write(log_lun,8)fe,err,errmess  
      case(32)
        errmess = ' Writing Orb header to='
        write(log_lun,9)fe,err,errmess,trim(cstat)  
      case(33)
        errmess = ' Writing prof structure info to='
        write(log_lun,9)fe,err,errmess,trim(cstat)        
      case(34)
        errmess = ' Writing Scan header to='
        write(log_lun,9)fe,err,errmess,trim(cstat)  
      case(35) 
        errmess = ' Writing Pixel to='
        write(log_lun,9)fe,err,errmess,trim(cstat)                
      
      case(40)
        errmess = ' Opening Profile Database'
        write(log_lun,8)fe,err,errmess        
      case(41)
        errmess = ' Reading Total Database Clusters'
        write(log_lun,10)fe,err,errmess
      case(42)
        errmess = ' Reading Satcode in Read Database'
        write(log_lun,9)fe,err,errmess
      case(43)
        errmess = ' Database satcode mismatch'
        write(log_lun,9)fe,err,errmess 
      case(44)
        errmess = ' Reading Database nominalEIAs'
        write(log_lun,9)fe,err,errmess      
      case(45)
        errmess = ' Reading Cluster header'
        write(log_lun,8)fe,err,errmess
      case(46)
        errmess = ' Max cluster exceeded in database '
        write(log_lun,10)fe,err,errmess
      case(47)
        errmess = ' Reading clustered database profiles '
        write(log_lun,8)fe,err,errmess
      
      case(50)
        errmess = ' Reading probability cutoff file '
        write(log_lun,9)fe,err,errmess,trim(cstat)
      case(51)
        errmess = ' Reading probability cutoff file header '
        write(log_lun,9)fe,err,errmess,trim(cstat)
      case(52)
        errmess = ' Reading probability cutoff file arrays'
        write(log_lun,9)fe,err,errmess,trim(cstat)
	
      end select
    
!--- call output routine to write out output file when error occurs

      call GPM_output_error(errmess)     
      
      return
      end subroutine GPM_reprt_err

!---------------------------------------------------------------------
 
      subroutine GPM_output_error(errmess)
       implicit none
       
!--- this routine will write out the output file with only the orbital header info.
!--- it's called whenever a fatal error is reported.  The orbital header is written
!--- out including the reason for the error in the 'spare' comments section.  NSCANS will
!--- be set to zero.


       character(len=51)  :: errmess
       logical            :: backsearch=.true.
       
       integer            :: ou_lun, ios, lin, pix, ilyr, reccnt=0
       integer            :: begi, i
       integer            :: temptime(8)

       character(len=5)   :: corbitn
       real :: tbdif

D      integer  :: cpu1, cpu2
D      real     :: cputime 
       
!--- create ouput structures

       type :: Date6         
          integer(kind=knd2):: year
          integer(kind=knd2):: month
          integer(kind=knd2):: day
          integer(kind=knd2):: hour
          integer(kind=knd2):: minute
          integer(kind=knd2):: second
       end type Date6
       
       type :: OrbitHdr                           !400 bytes  per file
          character(len=12) :: Satellite
          character(len=12) :: Sensor
          character(len=12) :: PreProcessorVersion
          character(len=12) :: AlgorithmVersion          
          character(len=128):: ProfileDatabaseFile          !***  this is new ***
          character(len=128):: RadiometerFile
          type(Date6)       :: FileCreationDate
          type(Date6)       :: GranuleStartDate
          type(Date6)       :: GranuleEndDate
          integer           :: GranuleNumber
          integer(kind=knd2):: NumScansGranule
          integer(kind=knd2):: NumPixelsScan
          integer(kind=knd1):: ProfStructFlag
          character(len=51) :: Spare51
       end type OrbitHdr
       
!--- assign structure name
       
       type(OrbitHdr) :: orbhdr
       
!---open output file

       call GPM_lun(ou_lun)
       open(unit=ou_lun,file=trim(output_file), access='stream',
     >                          status='unknown',iostat = ios)
       if(ios .ne. 0) then
         write(log_lun,*)'  Error opening ERROR output file : ',
     >                      trim(output_file)
         stop
       else
D         write(log_lun,*)' writing output file with errors : ',
D    >                      trim(output_file)
       endif
       
!-- Fill orbhdr with program variables
       
       orbhdr%Satellite           = sat_name
       orbhdr%Sensor              = sensor_name
       orbhdr%PreProcessorVersion = trim(pp_version)
       orbhdr%AlgorithmVersion    = trim(alg_version)
       orbhdr%ProfileDatabaseFile = trim(pdbfile)
       orbhdr%RadiometerFile      = trim(radfile)
              
       begi = scan(radfile,'/',backsearch) + 1
       orbhdr%RadiometerFile = trim(radfile(begi:256))

       call date_and_time(values=temptime)          
       orbhdr%FileCreationDate%year   = temptime(1)
       orbhdr%FileCreationDate%month  = temptime(2)
       orbhdr%FileCreationDate%day    = temptime(3)
       orbhdr%FileCreationDate%hour   = temptime(5)
       orbhdr%FileCreationDate%minute = temptime(6)
       orbhdr%FileCreationDate%second = temptime(7)

       orbhdr%GranuleStartDate%year   = miss_int       
       orbhdr%GranuleStartDate%month  = miss_int       
       orbhdr%GranuleStartDate%day    = miss_int
       orbhdr%GranuleStartDate%hour   = miss_int
       orbhdr%GranuleStartDate%minute = miss_int
       orbhdr%GranuleStartDate%second = miss_int
       orbhdr%GranuleEndDate%year     = miss_int
       orbhdr%GranuleEndDate%month    = miss_int
       orbhdr%GranuleEndDate%day      = miss_int
       orbhdr%GranuleEndDate%hour     = miss_int
       orbhdr%GranuleEndDate%minute   = miss_int
       orbhdr%GranuleEndDate%second   = miss_int
       
       orbhdr%GranuleNumber           = granule
       orbhdr%NumScansGranule         = 0                       !set scans and pixs to zero
       orbhdr%NumPixelsScan           = 0
       orbhdr%ProfStructFlag          = profstructflag
       orbhdr%Spare51                 = trim(errmess)

!--- write out orbit header info and quit


       write(log_lun,*)' writing orbhdr'

       write(ou_lun, iostat=ios)  orbhdr
       if(ios .ne. 0) write(log_lun,*)'Error writing error orb header'
       
!--- calc CPU since start and write out

D     call cpu_time(cputime)    !cpu from start of execution in seconds
D     cputime = cputime / 60.   !change to minutes
D     cpu1 = int(cputime)
D     cpu2 = nint((cputime - cpu1)*60)
D     write(log_lun,'(a,i4,a,i2.2)')'  CPU time (min:sec) = ', 
D    >         cpu1,':',cpu2

!--- close log file and halt execution

      close(log_lun)
      stop 1
            
      end subroutine GPM_output_error
            
!---------------------------------------------------------------------
        
      subroutine GPM_lun(ilun)
!    
!** This routine gets an open logical unit number starting with 100
!
       integer  :: ilun
       logical  :: llun
       character(len=128) :: blank=' '
!        
       do ilun = 100,201
         inquire(unit=ilun, opened=llun)         
         if(.not. llun) exit
       enddo
       if(ilun .eq. 201) call GPM_reprt_err(3,ilun,blank)     !problem getting lun
       return
      end subroutine GPM_lun            
 
 !-------------------------------------------------------------------------- 
            
      end module GPM_util_procedures
