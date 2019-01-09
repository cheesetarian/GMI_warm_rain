      module GPROF_2017_V1

       USE define_intam2
       USE GPM_arraydef_D
       USE GPM_read_procedures
       USE GPM_rain_procedures 
       USE GPM_util_procedures      	                                   ! 1/2016       
       implicit none

       contains

      subroutine run_gprof(ofile,lfile,adir)
       USE define_intam2
        ! changing GPROF from program to module, commenting out call to
        !  read preprocessor (taken care of elsewhere), and using slightly
        !  modified arraydef module, not running output routine
        ! D Duncan, CSU, 1/30/16
        ! Updated to GPROF_2014_V2 mid Feb 2016
        ! 4/2/16 -- calling GPM_output from main program, not here now
        ! Updated to GPROF_2017_V1 Jan 2017, dduncan, csu

!-------------------------------------------------------------------
!
!     This is a retrieval program for GPM which reads the GPM
!     database profile files, and runs the Baysian retrieval. This 
!     version is a combined S0 (Petty channel reduction) and S1 with
!     an argument switch  in the command line to decide which to run.
!
!     written in November, 2012
!         by Dave Randel / CSU
!
!     UPDATES:
!
!     Version B4 (11/30/2013) : Added Path variables read from the database
!                              in order for paths to be calculated even
!                              when layers are turned off.  Added precip
!                              diagnostics, improved the S0 retrieval, and
!                              reads from clustered database profiles.
!
!     Version 1 (2/15/2014)  : final release pre-GPM core launch, S0 disabled
!
!     Version 1-2 (4/2/2014) : includes the Delta Tb correction (new DBase format)
!                               includes S0 retrieval for SSMIS

!     Version 1-3 (5/15/2014): reads in the sensor and other scene dependent errors from
!                               a file  Files are in the ancillary directory (chansens)
!                               Also new from the preprocessor is the calibration file
!                               (if applied).   The calibration filename is also inserted
!                               into the first 40 orbhdr 'spare' bytes.

!     Version 1-4 (6/27/2014): modification of the rain retrieval to expand the number 
!                               of bins to look for matching profiles if the number of 
!                               significant profiles in the bayesian average is below a 
!                               specified threshold.  This increases CPU time.
!
!     Version 1-5 (12/15/2014): Rain_procedures use T2m instead of skint.  Removed temp
!				dependency on ocean expansion since V1-5 is used with new
!				Combined/GMI empirical database.   V1-5 database (GANAL 
!				ancillary, and autosnow data and preprocessors (V1412 and 
!				later) used T2m
!
!     Version 1-5 (8/1/2015) continued: added additional parameter (L1Cqualflag) from
!                               preprocessor to better define GPROF output qual (0-1-2)
!                               and also to allow for missing channels from sensors  
!
!     Version 2beta (11/2015)   Uses the Combined/GMI databases, removed all S0 support
!                               including command argument
!
!
!     GPROF2014_Version 2 (1/22/2016)   Version for PPS processing - no layers yet     
!
!     GPROF2014_Version 2 (3/10/2016)   Version for PPS processing with layers     
!                                       alg_version = '2014_V2_1603v
!
!     GPROF 2017 V1  (10/2016)  Version for GPM V5.  Includes rainrate
!     probability cutoffs for each bin.  New output format with less
!                               params.

!-------------------------------------------------------------------      
	           
       character(len=256) :: infile,ofile,lfile,adir
!--- Command line arguments - <infile><outfile><logfile><ancil directory><layerflag>

       integer            :: n_arguments, iargc, isfccode
       integer            :: isfc
       character(len=128) :: command,blank=' ', cstat= ' '
       character(len=1)   :: cpsflag
               
!--- Housekeeping variables
       integer ::  i, j, k, ios, istat, nbin

       integer  :: cpu1, cpu2
       real     :: cputime 
  		   
        
       !input_file = trim(infile)
       output_file= trim(ofile)
       log_file   = trim(lfile)
       dir_anc    = trim(adir)
!--- Get command line inputs
!       n_arguments = iargc() 
!       
!       if(n_arguments.ne.5) call GPM_reprt_err(1,n_arguments,blank)     !wrong arguments
!       call getarg(1, input_file)      !input filename
!       call getarg(2, output_file)     !output filename
!       call getarg(3, log_file)        !log filename
!       call getarg(4, dir_anc)         !ancillary directory
!       call getarg(5, cpsflag)         !layer switch   =  0 or 1  (1 = layers on)
              
!--- Identify hostname which is running job, write to log file
       command = 'hostname > ' // trim(log_file)
       call system(trim(command))	
       
!--- Open log file - this will stay open throughout
       call GPM_lun(log_lun)
       open(unit=log_lun,file=log_file,form='formatted',access='append',
     >      iostat = ios)
       if(ios .ne. 0) call GPM_reprt_err(2,log_lun,trim(log_file))      !opening log file
       
!--- Set algorithm version and get integer profile structure flag, retrieval type
!--- As of 4/2/2014  S0 is only functional for sensor=SSMIS.  If S0 requested and sensor
!--- is NOT SSMIS then an error is generated and procedure will terminate.

       cpsflag = '0' ! force layers off
       read(cpsflag,'(i1)')  profstructflag

       alg_version = '2017_V1_1611'      !Alg version
       
!--- startup printout information      
       write(log_lun,*)' input  file name    : ', trim(input_file)
       write(log_lun,*)' output file name    : ', trim(output_file)
       write(log_lun,*)' log file name       : ', trim(log_file)
       write(log_lun,*)' ancillary directory : ', trim(dir_anc)
       write(log_lun,*)' profile struc flag  : ', profstructflag 
       write(log_lun,*)' algorithm version   : ', alg_version
         
!--- define lat-lon subset area for retrieval

       latmin = -90.0; latmax = 90.0; lonmin = -180; lonmax = 180
!       latmin =-45; latmax = 0.0; lonmin = 160; lonmax = 180
       
D      if(latmin .ne. -90 .and. latmax .ne. 90) then
D          write(log_lun,'(a,4F8.2)')'  area subset applied ',
D    >                       latmin,latmax,lonmin,lonmax
D      endif

!--- set logical if Probability cutoffs are applied to the rain variables

c       probcut = .false.
       probcut = .true.
D      write(log_lun,*)'  Probability flag = ', probcut

!--- Read in preprocessor data file

!D      write(log_lun,*)' reading input preprocessor file'
!       call GPM_read_preprocessor
       
!--- SET EIA offset option

       if(trim(sensor_name).eq.'GMI')   EIAoffset = .false.     !set EIA offset option
       if(trim(sensor_name).eq.'TMI')   EIAoffset = .false.
       if(trim(sensor_name).eq.'SSMIS') EIAoffset = .true.
       if(trim(sensor_name).eq.'AMSR2') EIAoffset = .true.
       if(trim(sensor_name).eq.'SSMI')  EIAoffset = .true.
       write(log_lun,*)' EIA offset flag     : ', EIAoffset

!--- create sfccode histgram

      write(log_lun,*)' assigning surface type histogram'
       call GPM_hist_preprocessor       
       
!--- Initialize arrays

      write(log_lun,*)' starting GPM_init_arrays'
       call GPM_init_arrays
		 
!--- read in the profile clusters, will return array filled with missing
!--- if profstructflag is off

      write(log_lun,*)' starting prof clusters read routine'
       call GPM_read_profile_clusters

!--- read in the channel sensitivity table

      write(log_lun,*)' starting read channel sensitivity'
       call GPM_read_chansens
 
!--- read in the channel sensitivity table
      write(log_lun,*)' starting read rainsnowtwb'
       call GPM_read_rainsnowtwb                        

!--- read in the probability cutoff / rain fraction file

      write(log_lun,*)' starting GPM_read_probability_cutoff routine'
       call GPM_read_probability_cutoff

!--- loop over all surface classes

       do isfc = 1,2 !minsfccode, maxsfccode ! just running ocean for now!
       !do isfc = minsfccode, maxsfccode

        if(isfc .ne. 1) cycle                !only do type 1

         isfccode = isfc
         if(histpp(isfccode) .eq. 0) cycle          !don't do this surface type if no pixels

D        write(log_lun,*)
        write(log_lun,*)' reading profile database, surface code = ',
     >                       isfccode
	 call GPM_read_dbase(isfccode)

!--- Define EIA database offset from nominal for all pixels
         
         if(EIAoffset) then 
D            write(log_lun,*)'  starting EIA offset calculation'
             call GPM_calc_EIA_offset
         endif
	    
!---     GPM 2014 Rainfall retrieval

D        write(log_lun,*)'  starting rain retrieval'
         call GPM_rain(isfc)

       call cpu_time(cputime)    !cpu from start of execution in seconds
       cputime = cputime / 60.   !change to minutes
       cpu1 = int(cputime)
       cpu2 = nint((cputime - cpu1)*60)

       write(log_lun,'(a,i4,a,i2.2)')' CPU time (min:sec) = ',
     >        cpu1,':',cpu2

       enddo  !isfccode
       
D      write(log_lun,*)' finished rain retrieval'
D      write(log_lun,*)' number of bins requested, but not found: ',
D    >                   binnfcnt

!--- write out rainfall retrieval file

!D      write(log_lun,*)' writing scans'
       !call GPM_output
      			 
!--- Exit with dump of CPU time
      
       call cpu_time(cputime)    !cpu from start of execution in seconds
       cputime = cputime / 60.   !change to minutes
       cpu1 = int(cputime)
       cpu2 = nint((cputime - cpu1)*60) 

       write(log_lun,'(a,i4,a,i2.2)')' CPU time (min:sec) = ', 
     >        cpu1,':',cpu2
       write(log_lun,*)'Normal completion, no fatal errors'

!      stop 0
      end subroutine run_gprof
      END MODULE GPROF_2017_V1
