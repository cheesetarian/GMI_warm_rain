      module miemod

      use define_csu1dvar

      contains

      subroutine mie_rain(freqy,temp,lwc,Dn,mu,ksca,asca,gsca,RR)
      !subroutine mie_rain(freqy,temp,lwc,Nw,Dn,mu,ksca,asca,gsca)
      
!*****Changed 2015 by Sarah Ringerud to use normalized gamma DSD with 
!     input Nw and Dm (mass weighted mean diameter) following 
!     Bringi et al. 2003

!     Modified D Duncan Feb 2017 to accept LWC and calc Nw (though that can be
!     changed back), calculate everything as real*8, and allow for any mu value

!     Compute the extinction, absorption, asymmetry parameter and
!     backscatter for a given water content of rain in [g/m^3], and
!     a particle size distribution n(D) with shape param mu, median volume
!     diameter Dn.

!     Input:
!     freqy		frequency of radiation [GHz]
!     temp		temperature of particles [K]
!     lwc		water content of rain distribution [g/m^3]
!     Dn                median volume diameter of exponential DSD [mm]
!     mu                shape parameter of exp dsd
        ! [[ Dn/Dm = (3.67+mu)/(4+mu) ]]

!     Output:
!     ksca		extinction coefficient [1/km]
!     asca		single-scatter albedo []
!     gsca		asymmetry factor []
!     pbck		backscatter phase function/(4*pi) []

IMPLICIT NONE

INTEGER :: i      
REAL*8 :: freqy, temp, lwc,f1, f2,f3,f4, ef
REAL ::  ksca, asca, gsca, pbck
REAL*8 :: wavel,rad
REAL*8 :: Dm, mu, lam, num, x
REAL*8 :: mom6, gamm, dr, dD, D,Dn, Nw,Nw1,Nw2
REAL*8 ::  qext, qsca, asym, qbsca, bext, bsca, bsym, bq11      
REAL :: epsreal, epsimag
COMPLEX ::   ewat, cref
real*8 :: refit,refib,vel,rain 

REAL*8 :: dens =1.0 !density of liquid water in g/cm^3
REAL :: ss=0.0 !salinity of rain drops
!REAL*8 :: ss=0.0 !salinity of rain drops
! issue with doing 'watoptic' in double prec., revert to ss, epsreal/epsimag
! real*4
real :: RR
     
!     Assign some useful constants
      wavel=300./freqy
      
!     Begin by checking if hydrometeors of this species are present.
!     If not, set scattering parameters to zero and return.

      if(lwc .lt. 0.001) then
        ksca=0.
        asca=0.
        gsca=0.
        pbck=0.
        !print*,'near-zero rwc!',lwc
        return
      endif


!     If hydrometeors are present, initialize the scattering parameters

      bext=0.0
      bsca=0.0
      bsym=0.0
      bq11=0.0
      !mom6=0.

      refit=0.0
      refib=0.0
      rain =0.0

!     Loop over particle sizes:

!     increments of particle radius are 0.05 mm; the particle
!     size distribution is expressed as a particle number density,
!     num, per radius increment

!     Convert Nw (given log10(Nw))
      !Nw1 = 10**(Nw) ! convert from log
      gamm = exp(log_gamma(mu+4)) 
        !use log_gamma (natural log of gamma function!) 
        !because normal gamma function can't handle non-integers
      ef = (6.0/(4.0**4.0)) * ((4+mu)**(mu+4))/gamm
      Dm = Dn * (4+mu)/(3.67+mu)

        ! [g/m^3 / (g/cm^3*mm^4)]  -> [1.e3 mm^-1 m^-3]
      Nw2 =1.0e3*(4.0**4.0)*lwc / (epi*dens*(Dm**4)) ! eq 4
      !newNw = dlog10(Nw2) ! pass back through via global var
      !print*,'Nw,Dn,mu,rwc',Nw2,Dn,mu,lwc
      !print*,'Nw / Nw ',Nw,alog10(Nw2)
      
      dr = 0.03
      !dr = 0.05
      dD = dr*2.
      !do i = 0, 200
      do i = 0, 100 !200
        rad = 0.5*dr + dr*float(i)  ! in mm
        !rad = 0.5*dr + 0.05*float(i)  ! in mm
        D = 2.0*rad

        num = Nw2*ef*((D/Dm)**mu) *exp(-(4.0+mu)*D/Dm)  !Bringi et al 2003 eq 8a
        !num = Nw1*ef*((D/Dm)**mu) *exp(-(4.0+mu)*D/Dm)  !Bringi et al 2003 eq 8a

        if (num .lt. 1d-60) then !numerical instability?
         ksca=0.
         asca=0.
         gsca=0.
         pbck=0.
         RR = 0.0
         return
        endif

        x = 2.*epi*rad/wavel !size param

!       complex refractive index of liquid water
        call watoptic(real(freqy),real(temp),ss,epsreal,epsimag)
        ewat = cmplx(epsreal,epsimag)
        cref = csqrt(ewat)
        !print*,freqy,temp,epsreal,epsimag
        !print*,ewat,cref

!       call Mie program
        call mie_LJBtable_all(x,cref,qsca,qext,asym,qbsca)
        !print*,'x,cref',x,num!cref,i
        !print*,'sca,ext,asym,bsca: ',qsca,qext,asym,qbsca

!       integrate over particle size distribution;

        ! beta_ext = sum(num(D)*qext*pi*r^2*dD)
        bext=bext+num*qext*epi*rad*rad*dD*1e-3
        bsca=bsca+num*qsca*epi*rad*rad*dD*1e-3
        bsym=bsym+num*qsca*asym*epi*rad*rad*dD*1e-3
        bq11=bq11+num*qbsca*epi*rad*rad*dD*1e-3
        !mom6=mom6*D**6*dD*1e-3
        !print*,'bext,bsca,num: ',bext,bsca,num 
        !print*,'###qext,qsca: ',qext,qsca
        !print*,'epsr,epsi: ',epsreal,epsimag,ewat,cref 

        !refit = refit + epi*D*D*D*num*dD  ! since num is N(D), use D,dD
        !refib = refib + epi*D*D*num*dD
        !refit = refit + epi*rad*rad*rad*num*dD*1e-3 !eff rad numerator
        !refib = refib + epi*rad*rad*num*dD*1e-3 !eff rad denominator

        ! calculate velocity to get rain rate
        vel = 3.78*(D**0.67) ! D in mm, vel in m/s  [e.g. Bringi et al. 2004]
        rain = rain + (epi/6.0)*D*D*D*num*vel*dD*1e-6 
        ! m/s = mm^3 * mm^-1 m^-3 m/s mm * (1e-3)^2 (density is included with num)
        !print*,'-=-',rain*3600.0,vel,num!,bext

      end do

      !reffective = (refit/refib)
      !print*,'reff [mm]',reffective
!
!     check for distribution with very small extinction;
!     set parameters to zero to avoid numerical problems
      if( bext .gt. 1.e-15) then
        ksca=real(bext)
        asca=real(bsca/bext)
        gsca=real(bsym/bsca)
        pbck=real(bq11/bsca)
        RR = rain*3600.0  ! mm/s -> mm/hr
      else
        ksca=0.
        asca=0.
        gsca=0.
        pbck=0.
        RR = 0.0
      end if
        !print*,'mrain ',rain,RR,num,bext !bext may be nan

      return
      end subroutine mie_rain

      !subroutine mie_LJBtable_ALL(dielec_rlist_ALL,dielec_ilist_ALL,&
      subroutine mie_LJBtable_ALL(xin,crefin,qscaout,qextout,asymout,qbscaout)      
                      !crefindex_ALL,xlist_ALL,mietable_ALL_db, &
      implicit none
      !INCLUDE 'variables.inc'
      
      real*8     :: xin
      complex    :: crefin   
      real*8     :: qscaout,qextout,asymout,qbscaout      
      real*8     :: missing_mie_value=-99.99

!--------variables will be needed to apply 1D linear interpolation
!        dimensions are assigned as follows
      real*8  :: ra, rb, rc, rd, re      
      real*8 :: crefrcnt, creficnt, xcnt
      logical :: crefrfound, crefifound, xfound      
      real*8  :: r1,r2,i1,i2,x1,x2
            
      type  :: interpolation_structure
      	  real*8  :: qsca
	  real*8 :: qext
	  real*8  :: asym
	  real*8  :: qbsca
      end type interpolation_structure      
      type(interpolation_structure) :: bestcref_x1,bestcref_x2
      

!--- location closest cref real locations 
      crefrcnt=0
      crefrfound=.false.
      
      ra=1
      rb=ndielec_rlist_ALL
      rc=int((ra+rb)/2)

        !print*,'diel,rc ',dielec_rlist_ALL(rc),rc
      do while ((crefrcnt.lt.10).and.(.not.crefrfound))
        if (real(crefin).ge.dielec_rlist_ALL(rc)) then
	  !on the right-hand side
	  ra=rc
	  rb=rb
	  if (rc.eq.int((ra+rb)/2)) crefrfound=.true.
	  rc=int((ra+rb)/2)
	else
	  !on the left-hand side
	  ra=ra
	  rb=rc-1
	  if (rc.eq.int((ra+rb)/2)) crefrfound=.true.
	  rc=int((ra+rb)/2)
	endif
	crefrcnt=crefrcnt+1.0
      enddo!

      ra=rc-3
      rb=rc+3
      if (ra.lt.1) ra=1
      if (rb.gt.ndielec_rlist_ALL) rb=ndielec_rlist_ALL
      rd=10000.0
      re=0.0
      do rc=ra, rb, 1
        if (abs(real(crefin)-dielec_rlist_ALL(rc)).lt.rd) then
	  rd=abs(real(crefin)-dielec_rlist_ALL(rc))
	  re=rc
	endif
      enddo
      if (real(crefin).ge.dielec_rlist_ALL(re)) then
        ra=re
	rb=re+1
	if (rb.gt.ndielec_rlist_ALL) then
	  ra=ndielec_rlist_ALL-1
	  rb=ndielec_rlist_ALL
	  crefin=cmplx(dielec_rlist_ALL(rb)-1.0e-6,aimag(crefin))
	endif
      else
        ra=re-1
	rb=re
	if (ra.lt.1) then
	  ra=1
	  rb=2
	  crefin=cmplx(dielec_rlist_ALL(ra)+1.0e-6,aimag(crefin))
	endif
      endif
      r1=ra
      r2=rb
      

 
 
!-------location closest cref imaginary locations      
      creficnt=0
      crefifound=.false.
      
      ra=1
      rb=ndielec_ilist_ALL
      rc=int((ra+rb)/2)
      
      do while ((creficnt.lt.15).and.(.not.crefifound))
        if (aimag(crefin).ge.dielec_ilist_ALL(rc)) then
	  !on the right-hand side
	  ra=rc
	  rb=rb
	  if (rc.eq.int((ra+rb)/2)) crefifound=.true.
	  rc=int((ra+rb)/2)
	else
	  !on the left-hand side
	  ra=ra
	  rb=rc-1
	  if (rc.eq.int((ra+rb)/2)) crefifound=.true.
	  rc=int((ra+rb)/2)
	endif
	creficnt=creficnt+1.0
      enddo!

      ra=rc-3
      rb=rc+3
      if (ra.lt.1) ra=1
      if (rb.gt.ndielec_ilist_ALL) rb=ndielec_ilist_ALL
      rd=10000.0
      re=0.0
      do rc=ra, rb, 1
        if (abs(aimag(crefin)-dielec_ilist_ALL(rc)).lt.rd) then
	  rd=abs(aimag(crefin)-dielec_ilist_ALL(rc))
	  re=rc
	endif
      enddo
      if (aimag(crefin).ge.dielec_ilist_ALL(re)) then
        ra=re
	rb=re+1
	if (rb.gt.ndielec_ilist_ALL) then
	  ra=ndielec_ilist_ALL-1
	  rb=ndielec_ilist_ALL
	  crefin=cmplx(real(crefin),dielec_ilist_ALL(rb)-1.0e-6)
	endif
      else
        ra=re-1
	rb=re
	if (ra.lt.1) then
	  ra=1
	  rb=2
	  crefin=cmplx(real(crefin),dielec_ilist_ALL(ra)+1.0e-6)
	endif
      endif
      i1=ra
      i2=rb
      
     
     
!-------location closest xsize locations      
      xcnt=0
      xfound=.false.
      
      ra=1
      rb=nxlist_ALL
      rc=int((ra+rb)/2)
      
      do while ((xcnt.lt.10).and.(.not.xfound))
        if (xin.ge.xlist_ALL(rc)) then
	  !on the right-hand side
	  ra=rc
	  rb=rb
	  if (rc.eq.int((ra+rb)/2)) xfound=.true.
	  rc=int((ra+rb)/2)
	else
	  !on the left-hand side
	  ra=ra
	  rb=rc-1
	  if (rc.eq.int((ra+rb)/2)) xfound=.true.
	  rc=int((ra+rb)/2)
	endif
	xcnt=xcnt+1.0
      enddo!

      ra=rc-3
      rb=rc+3
      if (ra.lt.1) ra=1
      if (rb.gt.nxlist_ALL) rb=nxlist_ALL
      rd=10000.0
      re=0.0
      do rc=ra, rb, 1
        if (abs(xin-xlist_ALL(rc)).lt.rd) then
	  rd=abs(xin-xlist_ALL(rc))
	  re=rc
	endif
      enddo
      if (xin.ge.xlist_ALL(re)) then
        ra=re
	rb=re+1
	if (rb.gt.nxlist_ALL) then
	  ra=nxlist_ALL-1
	  rb=nxlist_ALL
	  xin=xlist_ALL(rb)-1.0e-6
	endif
      else
        ra=re-1
	rb=re
	if (ra.lt.1) then
	  ra=1
	  rb=2
	  xin=xlist_ALL(ra)+1.0e-6
	endif
      endif
      x1=ra
      x2=rb
      !print*,'x1,x2 ',x1,x2
      

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--------APPLY 2D (for crefr+crefi) + 1D(for x) Linear Interpolation!
!------For all variables at x=x1
      bestcref_x1%qsca=1.0/( (crefindex_ALL(r2,i1)%real_cref &
     		              -crefindex_ALL(r1,i1)%real_cref) &
     			    *(crefindex_ALL(r1,i2)%imag_cref &
     			      -crefindex_ALL(r1,i1)%imag_cref)) &
      * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x1)%qsca &
     	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
         + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x1)%qsca &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
              *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
     	  + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x1)%qsca &
          *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)) &
         + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x1)%qsca &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))


      bestcref_x1%qext=1.0/( (crefindex_ALL(r2,i1)%real_cref &
     		              -crefindex_ALL(r1,i1)%real_cref) &
     			    *(crefindex_ALL(r1,i2)%imag_cref &
     			      -crefindex_ALL(r1,i1)%imag_cref)) &
      * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x1)%qext &
     	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
         + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x1)%qext &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
              *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
     	  + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x1)%qext &
     	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)) &
         + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x1)%qext &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))
     
      bestcref_x1%asym=1.0/( (crefindex_ALL(r2,i1)%real_cref &
     		              -crefindex_ALL(r1,i1)%real_cref) &
     			    *(crefindex_ALL(r1,i2)%imag_cref &
     			      -crefindex_ALL(r1,i1)%imag_cref)) &
      * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x1)%asym &
     	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
         + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x1)%asym &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
              *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
     	  + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x1)%asym &
     	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)) &
         + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x1)%asym &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))
     
      bestcref_x1%qbsca=1.0/( (crefindex_ALL(r2,i1)%real_cref &
     		              -crefindex_ALL(r1,i1)%real_cref) &
     			    *(crefindex_ALL(r1,i2)%imag_cref &
     			      -crefindex_ALL(r1,i1)%imag_cref)) &
      * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x1)%qbsca &
     	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
        + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x1)%qbsca &
    	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
              *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
     	  + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x1)%qbsca &
     	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)) &
         + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x1)%qbsca &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))          
     

!------For all variables at x=x2
      bestcref_x2%qsca=1.0/( (crefindex_ALL(r2,i1)%real_cref &
     		              -crefindex_ALL(r1,i1)%real_cref) &
     			    *(crefindex_ALL(r1,i2)%imag_cref &
     			      -crefindex_ALL(r1,i1)%imag_cref)) &
      * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x2)%qsca &
    	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
         + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x2)%qsca &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
              *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
     	  + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x2)%qsca &
     	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)) &
         + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x2)%qsca &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))


      bestcref_x2%qext=1.0/( (crefindex_ALL(r2,i1)%real_cref &
     		              -crefindex_ALL(r1,i1)%real_cref) &
     			    *(crefindex_ALL(r1,i2)%imag_cref &
     			      -crefindex_ALL(r1,i1)%imag_cref)) &
      * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x2)%qext &
     	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
         + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x2)%qext &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
              *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
     	  + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x2)%qext &
     	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)) &
         + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x2)%qext &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))
     
      bestcref_x2%asym=1.0/( (crefindex_ALL(r2,i1)%real_cref &
     		              -crefindex_ALL(r1,i1)%real_cref) &
     			    *(crefindex_ALL(r1,i2)%imag_cref &
     			      -crefindex_ALL(r1,i1)%imag_cref)) &
      * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x2)%asym &
     	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
         + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x2)%asym &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
              *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
     	  + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x2)%asym &
     	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) & 
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)) &
         + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x2)%asym &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))
     
      bestcref_x2%qbsca=1.0/( (crefindex_ALL(r2,i1)%real_cref &
     		              -crefindex_ALL(r1,i1)%real_cref) &
     			    *(crefindex_ALL(r1,i2)%imag_cref &
     			      -crefindex_ALL(r1,i1)%imag_cref)) &
      * (  mietable_ALL_db(crefindex_ALL(r1,i1)%loc_cref,x2)%qbsca &
     	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
         + mietable_ALL_db(crefindex_ALL(r2,i1)%loc_cref,x2)%qbsca &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
              *(crefindex_ALL(r1,i2)%imag_cref-aimag(crefin))) &
     	  + mietable_ALL_db(crefindex_ALL(r1,i2)%loc_cref,x2)%qbsca &
     	     *( (crefindex_ALL(r2,i1)%real_cref-real(crefin)) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)) &
         + mietable_ALL_db(crefindex_ALL(r2,i2)%loc_cref,x2)%qbsca &
     	     *( (real(crefin)-crefindex_ALL(r1,i1)%real_cref) &
     	       *(aimag(crefin)-crefindex_ALL(r1,i1)%imag_cref)))          
     
!------Interpolate between x1 and x2

      qscaout=bestcref_x1%qsca+(bestcref_x2%qsca-bestcref_x1%qsca) &
     	     *(xin-xlist_ALL(x1))/(xlist_ALL(x2)-xlist_ALL(x1))

      qextout=bestcref_x1%qext+(bestcref_x2%qext-bestcref_x1%qext) &
     	     *(xin-xlist_ALL(x1))/(xlist_ALL(x2)-xlist_ALL(x1))

      asymout=bestcref_x1%asym+(bestcref_x2%asym-bestcref_x1%asym) &
     	     *(xin-xlist_ALL(x1))/(xlist_ALL(x2)-xlist_ALL(x1))
     
      qbscaout=bestcref_x1%qbsca+(bestcref_x2%qbsca-bestcref_x1%qbsca) &
     	     *(xin-xlist_ALL(x1))/(xlist_ALL(x2)-xlist_ALL(x1))          


        !print*,'things',xin,xlist_all(x1),xlist_all(x2)
        !print*,bestcref_x2%qext,bestcref_x1%qext
        !print*,bestcref_x2%qsca,bestcref_x1%qsca
        !print*,bestcref_x2%asym,bestcref_x1%asym
        !print*,bestcref_x2%qbsca,bestcref_x1%qbsca
        !print*,'qsc,qex,qbsca',qscaout,qextout,qbscaout
      if ((qscaout.lt.0).or.(qextout.lt.0).or.(qbscaout.lt.0)) then
	
	if (min(qscaout,qextout,qbscaout).lt.(missing_mie_value/2.0)) then
	  !the value need to be interpolated is closer to the NAN, instead of any valid numbers
	  qscaout=0.0/0.0
	  qextout=0.0/0.0
	  asymout=0.0/0.0
	  qbscaout=0.0/0.0		
	else
	  !the value need to be interpolated is closer to a nearby valid point, instead of NAN	  
	  if (  (bestcref_x1%qsca.gt.0).and. &
     	        (bestcref_x1%qext.gt.0).and. &
     	        (bestcref_x1%qbsca.gt.0)    )then
	    !if interpolation at x1 is good
	    qscaout=bestcref_x1%qsca
	    qextout=bestcref_x1%qext
	    asymout=bestcref_x1%asym
	    qbscaout=bestcref_x1%qbsca
	  else
	    if (  (bestcref_x2%qsca.gt.0).and. &
     	          (bestcref_x2%qext.gt.0).and. &
     		  (bestcref_x2%qbsca.gt.0)    )then
	      !if interpolation at x2 is good
	      qscaout=bestcref_x2%qsca
	      qextout=bestcref_x2%qext
	      asymout=bestcref_x2%asym
	      qbscaout=bestcref_x2%qbsca	        
	    else
	      !if interpolation is bad at both x1 and x2
	      qscaout=0.0/0.0
	      qextout=0.0/0.0
	      asymout=0.0/0.0
	      qbscaout=0.0/0.0	      
	    endif!for if interpolation at x2
	  endif!for if interpolation at x1
	endif!for if relative location of interpolated result	

      endif!for if results are unreasonable
      
      return
      end subroutine mie_LJBtable_ALL   

! ---------------------------------------------------------------------------------------

      subroutine watoptic(freqy, temp, salinity, epsreal, epsimag)
      implicit   none

!     Input & output variables; eps the dielectric constant, epsilon  
      real     freqy, temp, salinity
      real     epsreal, epsimag

!C     internal variables      
      real     freqhz, ctemp!, pi
      real     omega, epsstat, trelax, fac1

      real ::  epshigh= 4.90
      real ::  pi=3.141592654

      freqhz = freqy*1.0e+09

      ctemp  = temp - 273.16
      omega  = 2.*pi*freqhz

      epsstat = (87.134 - 1.949E-01 * ctemp &
              - 1.276E-02 * ctemp * ctemp &
              + 2.491E-04 * ctemp * ctemp * ctemp) &
              * (1.0 + 1.613E-05 * salinity * ctemp &
              - 3.656E-03 * salinity &
              + 3.210E-05 * salinity * salinity &
              - 4.232E-07 * salinity * salinity * salinity)


      trelax = (1.768E-11 - 6.086E-13 * ctemp &
             + 1.104E-14 * ctemp * ctemp &
             - 8.111E-17 * ctemp * ctemp * ctemp) &
             * (1.0 + 2.282E-05 * salinity * ctemp &
             - 7.638E-04  * salinity &
             - 7.760E-06  * salinity * salinity &
             + 1.105E-08  * salinity * salinity * salinity)

      fac1    = 1.0 + omega*omega*trelax*trelax
      epsreal = epshigh + (epsstat - epshigh) / fac1
      epsimag = ((epsstat - epshigh) * omega * trelax) / fac1

      return
      end subroutine watoptic
      

      end module miemod
