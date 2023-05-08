! v1.0 29/11/21
! v2.0 14/08/22 - Updated to read off angle phi from SIAF through fits_siaf_read_v2v3.f90 writing msa_ref_v2v3.ascii file

     include '../empt_modules.f90'
!-------------------------------------
      program calculate_derived_parameters
      
! calculates:
! - mean on-sky projection  slit pitch and open area over viable slitlet map

      use path_to_ref_files
      use config_parameters, only : n_filt
      implicit none
! --- MSA Parameters------------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
!--------------------------------
      integer i,j,k,ii,ncount
      real slw,slh
      real avslw,avslh      
      real xmc(4),ymc(4)
      real v2c(4),v3c(4)
      real xm(2),ym(2)
      real v2(2),v3(2)
      real v1x,v1y,v2x,v2y,v1norm,v2norm
      real av_xpitch,av_ypitch,av_xopen,av_yopen
      real phi,phi_r
      real xm_ref,ym_ref
      real v2_ref,v3_ref,AngV3
      real avxgap,avxpitch
      character*70 model_directory
      character*70 filename

      ref_path='./' !overriden nominal path because run out of /reference_files
 
       n_filt=4 !CLEAR
    
       call read_V2V3_ref(V2_ref,V3_ref,AngV3)
             
       phi=180.-AngV3
              
       phi_r=phi/57.2957795131

       print *
       print *,'AngV3, phi:',AngV3,phi

       print *,'V2_ref,V3_ref:',V2_ref,V3_ref
       
       call sky_to_msa_full(V2_ref,V3_ref,xm_ref,ym_ref)

       print *,'xm_ref,ym_ref:',xm_ref,ym_ref
 
       
! inverse sanity check

!        call msa_to_sky_full(xm_ref,ym_ref,V2_ref,V3_ref)
!  
!        print *,'V2_ref,V3_ref:',V2_ref,V3_ref
! 
!        print *


! get instrument model name for reference

      filename=trim(adjustl(ref_path))//'model.conf'
      
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) model_directory
      close(20)


      write(*,*)
      write(*,*) 'Instrument Model: ',model_directory

   
! initialize msa routines

      call load_msa_dimensions_bd

! calculate mean shutter pitches on sky

      avslw=0.
      avslh=0.
      
      ncount=0

      do k=1,4
      do i=1,365
      do j=1,171
            
      ncount=ncount+1

! get corner coordinates of shutter on MSA
      
      call shutter_to_msa_corners(k,i,j,xmc,ymc)

! project to V2,V3
      
      do ii=1,4
      call msa_to_sky_bd(xmc(ii),ymc(ii),v2c(ii),v3c(ii))
      enddo

      slw=(sqrt((v2c(3)-v2c(4))*(v2c(3)-v2c(4))+(v3c(3)-v3c(4))*(v3c(3)-v3c(4)))+sqrt((v2c(1)-v2c(2))*(v2c(1)-v2c(2))+(v3c(1)-v3c(2))*(v3c(1)-v3c(2))))/2.*1000.
      slh=(sqrt((v2c(3)-v2c(1))*(v2c(3)-v2c(1))+(v3c(3)-v3c(1))*(v3c(3)-v3c(1)))+sqrt((v2c(4)-v2c(2))*(v2c(4)-v2c(2))+(v3c(4)-v3c(2))*(v3c(4)-v3c(2))))/2.*1000.
      
      avslw=avslw+slw
      avslh=avslh+slh     
      
      enddo
      enddo
      enddo

      av_xpitch=avslw/float(ncount)
      av_ypitch=avslh/float(ncount)
 

! calculate mean shutter open area dimensions on sky

      avslw=0.
      avslh=0.
      
      ncount=0


      do k=1,4
      do i=1,365
      do j=1,171

      ncount=ncount+1
 
! get corner coordinates of open area on MSA
      
      call shutter_to_msa_corners_open(k,i,j,xmc,ymc)
 
! project to V2,V3
      
      do ii=1,4
      call msa_to_sky_bd(xmc(ii),ymc(ii),v2c(ii),v3c(ii))
      enddo

      slw=(sqrt((v2c(3)-v2c(4))*(v2c(3)-v2c(4))+(v3c(3)-v3c(4))*(v3c(3)-v3c(4)))+sqrt((v2c(1)-v2c(2))*(v2c(1)-v2c(2))+(v3c(1)-v3c(2))*(v3c(1)-v3c(2))))/2.*1000.
      slh=(sqrt((v2c(3)-v2c(1))*(v2c(3)-v2c(1))+(v3c(3)-v3c(1))*(v3c(3)-v3c(1)))+sqrt((v2c(4)-v2c(2))*(v2c(4)-v2c(2))+(v3c(4)-v3c(2))*(v3c(4)-v3c(2))))/2.*1000.
            
      avslw=avslw+slw
      avslh=avslh+slh
            
      enddo
      enddo
      enddo

      av_xopen=avslw/float(ncount)
      av_yopen=avslh/float(ncount)

      write(*,*)
      write(*,*) 'Average shutter x pitch [mas]:',av_xpitch
      write(*,*) 'Average shutter y pitch [mas]:',av_ypitch
      write(*,*) 'Average shutter x open [mas]: ',av_xopen
      write(*,*) 'Average shutter y open [mas]: ',av_yopen
      write(*,*)
      
! calculate NIRSpec rotation w.r.t. V2,V3 phi


! vector along MSA y axis

!       xm(1)=0.
!       ym(1)=0.
!       xm(2)=xm(1)
!       ym(2)=0.5
!       
!       v1x=xm(2)-xm(1)
!       v1y=ym(2)-ym(1)
!       
!       v1norm=sqrt(v1x*v1x+v1y*v1y)
!       
!       v1x=v1x/v1norm
!       v1y=v1y/v1norm
!       
! ! project to V2,v3
! 
!       do ii=1,2
!       call msa_to_sky_bd(xm(ii),ym(ii),v2(ii),v3(ii))
!       enddo
! 
!       v2x=v2(2)-v2(1)
!       v2y=v3(2)-v3(1)
! 
!       v2norm=sqrt(v2x*v2x+v2y*v2y)
!       
!       v2x=v2x/v2norm
!       v2y=v2y/v2norm
! 
!       phi_r=acos(v1x*v2x+v1y*v2y)
!       
!       phi=phi_r*57.29578
!       
!       write(*,*) 'phi: ',phi
!       write(*,*)

 
 
! calculate avxgap between opposing MSA quadrants

        avxpitch=(xpitch(1)+xpitch(2)+xpitch(3)+xpitch(4))/4.

        avxgap=0.
        ncount=0.
        
       do j=1,171
        call shutter_to_msa_bd(3,1,j,xm(1),ym(1))
        call shutter_to_msa_bd(1,365,j,xm(2),ym(2))
        avxgap=avxgap+(xm(1)-xm(2))/avxpitch
        ncount=ncount+1
        call shutter_to_msa_bd(4,1,j,xm(1),ym(1))
        call shutter_to_msa_bd(2,365,j,xm(2),ym(2))
        avxgap=avxgap+(xm(1)-xm(2))/avxpitch
        ncount=ncount+1
       enddo

       avxgap=avxgap/float(ncount)-1.
       
       print *,'avxgap: ',avxgap
       print *
       
       
       
     
! write out reference file
 
       open(10,file='derived_parameters.ascii',status='replace',access='sequential', form='formatted')

       write(10,'(a)') '# Derived Instrument Parameters'    
       write(10,'(a)') '# Instrument Model Reference: '//trim(adjustl(model_directory))   
       write(10,*) 
       write(10,'(a)') '*av_xpitch  !mas'    
       write(10,'(f7.2)') av_xpitch    
       write(10,*) 
       write(10,'(a)') '*av_ypitch  !mas'    
       write(10,'(f7.2)') av_ypitch    
       write(10,*) 
       write(10,'(a)') '*av_xopen  !mas'    
       write(10,'(f7.2)') av_xopen    
       write(10,*) 
       write(10,'(a)') '*av_yopen  !mas'    
       write(10,'(f7.2)') av_yopen    
       write(10,*) 
       write(10,'(a)') '*phi  !deg'    
       write(10,'(f9.5)') phi    
       write(10,*) 
       write(10,'(a)') '*phi_r  !rad'    
       write(10,'(f10.6)') phi_r    
       write(10,*) 
       write(10,'(a)') '*V2_ref  !arcsec - should match active SIAF entry - NRS_FULL_MSA: V2Ref'
       write(10,'(f10.5)') V2_ref    
       write(10,*) 
       write(10,'(a)') '*V3_ref  !arcsec - should match active SIAF entry - NRS_FULL_MSA: V3Ref'
       write(10,'(f10.5)') V3_ref    
       write(10,*) 
       write(10,'(a)') '*xm_ref  !mm - calculated from V2_ref,V3_ref'    
       write(10,'(f12.7)') xm_ref    
       write(10,*) 
       write(10,'(a)') '*ym_ref  !mm - calculated from V2_ref,V3_ref'     
       write(10,'(f12.7)') ym_ref    
       write(10,*) 
       write(10,'(a)') '*avxgap  !mm' 
       write(10,'(f8.3)') avxgap    
       write(10,*) 
         
       close(10)
     
      stop
      end

!---------------------------------------------------------------------------
 
include './front_transforms_updt.f90'

include './back_transforms_updt.f90'

include './shutter_routines_new.f90'

include './misc_routines_new.f90'

!-------------------------------------------------------------------------

      subroutine shutter_to_msa_corners_open(mq,i,j,xmc,ymc)

!  Input:  mq = MSA Quadrant [1,2,3,4]
!           i = shutter x coordinate [1,365]
!           j = shutter y coordinate [1,171]
!  Output: xm,ym = position of center of shutter on MSA (mm)
!
! Asumes prior call of load_msa_dimensions
!
      use path_to_ref_files
      implicit none
      real xmc(4),ymc(4)
      integer mq,i,j
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
!-------------------------
       real xmm,ymm,xmr,ymr,corot,sirot

      if ((mq.lt.1).or.(mq.gt.4)) stop 'mq out of bounds: 1'
      if ((i.lt.1).or.(i.gt.mx)) stop 'i out of bounds: 1'
      if ((j.lt.1).or.(j.gt.my)) then
      write(*,*) mq,i,j
      stop 'j out of bounds: 1'
      endif

      corot=cos(-rot(mq))
      sirot=sin(-rot(mq))

! lower left corner of shutter  - 1

      xmr=(float(i)-1.0)*xopen(mq)
      ymr=(float(j)-1.0)*yopen(mq)
      xmm=xmr*corot+ymr*sirot
      ymm=-xmr*sirot+ymr*corot
      xmc(1)=xmm+x0(mq)-xopen(mq)/2.
      ymc(1)=ymm+y0(mq)-yopen(mq)/2.

! lower right corner of shutter - 2

      xmr=(float(i))*xopen(mq)
      ymr=(float(j)-1.0)*yopen(mq)
      xmm=xmr*corot+ymr*sirot
      ymm=-xmr*sirot+ymr*corot
      xmc(2)=xmm+x0(mq)-xopen(mq)/2.
      ymc(2)=ymm+y0(mq)-yopen(mq)/2.


! upper left corner of shutter  - 3

      xmr=(float(i)-1.0)*xopen(mq)
      ymr=(float(j))*yopen(mq)
      xmm=xmr*corot+ymr*sirot
      ymm=-xmr*sirot+ymr*corot
      xmc(3)=xmm+x0(mq)-xopen(mq)/2.
      ymc(3)=ymm+y0(mq)-yopen(mq)/2.


! upper right corner of shutter - 4 

      xmr=(float(i))*xopen(mq)
      ymr=(float(j))*yopen(mq)
      xmm=xmr*corot+ymr*sirot
      ymm=-xmr*sirot+ymr*corot
      xmc(4)=xmm+x0(mq)-xopen(mq)/2.
      ymc(4)=ymm+y0(mq)-yopen(mq)/2.

      return
      end

!----------------------------------------------------------------------

      subroutine shutter_to_msa_corners(mq,i,j,xmc,ymc)

!  Input:  mq = MSA Quadrant [1,2,3,4]
!           i = shutter x coordinate [1,365]
!           j = shutter y coordinate [1,171]
!  Output: xm,ym = position of center of shutter on MSA (mm)
!
! Asumes prior call of load_msa_dimensions
!
      use path_to_ref_files
      implicit none
      real xmc(4),ymc(4)
      integer mq,i,j
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
!-------------------------
       real xmm,ymm,xmr,ymr,corot,sirot

      if ((mq.lt.1).or.(mq.gt.4)) stop 'mq out of bounds: 1'
      if ((i.lt.1).or.(i.gt.mx)) stop 'i out of bounds: 1'
      if ((j.lt.1).or.(j.gt.my)) then
      write(*,*) mq,i,j
      stop 'j out of bounds: 1'
      endif

      corot=cos(-rot(mq))
      sirot=sin(-rot(mq))

! lower left corner of shutter  - 1

      xmr=(float(i)-1.0)*xpitch(mq)
      ymr=(float(j)-1.0)*ypitch(mq)
      xmm=xmr*corot+ymr*sirot
      ymm=-xmr*sirot+ymr*corot
      xmc(1)=xmm+x0(mq)-xpitch(mq)/2.
      ymc(1)=ymm+y0(mq)-ypitch(mq)/2.

! lower right corner of shutter - 2

      xmr=(float(i))*xpitch(mq)
      ymr=(float(j)-1.0)*ypitch(mq)
      xmm=xmr*corot+ymr*sirot
      ymm=-xmr*sirot+ymr*corot
      xmc(2)=xmm+x0(mq)-xpitch(mq)/2.
      ymc(2)=ymm+y0(mq)-ypitch(mq)/2.


! upper left corner of shutter  - 3

      xmr=(float(i)-1.0)*xpitch(mq)
      ymr=(float(j))*ypitch(mq)
      xmm=xmr*corot+ymr*sirot
      ymm=-xmr*sirot+ymr*corot
      xmc(3)=xmm+x0(mq)-xpitch(mq)/2.
      ymc(3)=ymm+y0(mq)-ypitch(mq)/2.


! upper right corner of shutter - 4 

      xmr=(float(i))*xpitch(mq)
      ymr=(float(j))*ypitch(mq)
      xmm=xmr*corot+ymr*sirot
      ymm=-xmr*sirot+ymr*corot
      xmc(4)=xmm+x0(mq)-xpitch(mq)/2.
      ymc(4)=ymm+y0(mq)-ypitch(mq)/2.

      return
      end

!---------------------------------------------------------------------------

      subroutine sky_to_msa_full(theta_x,theta_y,xm,ym)
!
!  Input: theta_x,theta_y = angular incidence angles on OTE (arcsec)
!  Output:  xm,ym = position in MSA plane (mm)
!
      implicit none
      real theta_x,theta_y,xm,ym
      real xs,ys

      call ote_to_stop_full(theta_x,theta_y,xs,ys)
      call stop_to_msa_bd(xs,ys,xm,ym)

      return
      end

!---------------------------------------------------------------------------

      subroutine msa_to_sky_full(xm,ym,theta_x,theta_y)
!
!  Input:  xm,ym = position in MSA plane (mm)
!  Output: theta_x,theta_y = angular incidence angles on OTE (arcsec)
!
      implicit none
      real xm,ym,theta_x,theta_y
      real xs,ys

      call msa_to_stop_bd(xm,ym,xs,ys)
      call stop_to_ote_full(xs,ys,theta_x,theta_y)

      return
      end

!---------------------------------------------------------------------------

      subroutine ote_to_stop_full(theta_x,theta_y,xs,ys)
!
!  Input: theta_x,theta_y = angular incidence angles on OTE (arcsec)
!  Output: xs,ys = position in field stop (mm)
!
      use path_to_ref_files
      implicit none
      real xs,ys,theta_x,theta_y
      real*8 xin,yin,xp,yp,xout,yout
      real*8 gammax,gammay,theta,x0in,y0in,x0out,y0out
      real*8 A(6,6),B(6,6)
      integer i,j
      character*70 filename
      character*70 model_directory
      character*30 parameter_tag
      logical first
      save first,A,B,gammax,gammay,theta,x0in,y0in,x0out,y0out
      data first/.true./

      if (.not.first) goto 1000

      first=.false.

!       gammax = -2.296 ! f_ote=131 m when input [deg] output [m]
!       gammay = 2.2953
! 
!       theta = 0.2075/57.29578 !Convert to rad
! 
!       x0in = 0.10539  !V2V3 Position of NIRSpec FOV center?
!       y0in = -0.11913000025
! 
!       x0out = -5.18289805611e-07
!       y0out = -1.92704532397e-09
! 
!       A(1,1) = -6.764695218042509e-12
!       A(1,2) = -0.001181845439725266
!       A(1,3) = -0.004421942263252729
!       A(1,4) = -0.00015786783293089
!       A(1,5) = -8.101344739691103e-05
!       A(1,6) = -0.00223528778759241
!       A(2,1) = 0.9999913490363831
!       A(2,2) = 0.01985321796323403
!       A(2,3) = 0.01836568669267155
!       A(2,4) = 0.0008476551557795109
!       A(2,5) = -0.004583374077531843
!       A(3,1) = 0.01299262834412884
!       A(3,2) = 0.0003111090225528627
!       A(3,3) = 0.0006487537603577787
!       A(3,4) = 0.01831138758726425
!       A(4,1) = 0.01829981571350813
!       A(4,2) = 0.0008920392687617073
!       A(4,3) = 0.01402112565533553
!       A(5,1) = 0.0007670516309624491
!       A(5,2) = 0.009553952594322013
!       A(6,1) = 0.0004297874598029328
! 
! 
!       B(1,1) = -6.115026307974087e-10
!       B(1,2) = 0.9999869591970522
!       B(1,3) = 0.01481498123407221
!       B(1,4) = 0.0183504542886734
!       B(1,5) = 0.0009281781374761249
!       B(1,6) = 0.0002958977074083435
!       B(2,1) = -0.001179801005691817
!       B(2,2) = 0.01742851464230401
!       B(2,3) = 0.0003646472234911172
!       B(2,4) = 0.0008189742785870893
!       B(2,5) = -0.001483363225125878
!       B(3,1) = -0.005044248520802729
!       B(3,2) = 0.01828606785677567
!       B(3,3) = 0.0009112362883588386
!       B(3,4) = -0.009363002229886064
!       B(4,1) = -0.0001649665026415299
!       B(4,2) = 0.00081286189015533
!       B(4,3) = 0.002997688363766571
!       B(5,1) = -6.145531569187097e-05
!       B(5,2) = -0.007795457001419592
!       B(6,1) = -2.764978182656641e-05



      filename=trim(adjustl(ref_path))//'model.conf'
      
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) model_directory
      close(20)
      
!       print *,model_directory
      
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/CoordTransform/OTE.pcf'

!       print *,filename

      open(20,file=filename,status='old',access='sequential',form='formatted')

10    read(20,*,end=100) parameter_tag
      if (parameter_tag(1:1).ne.'*') goto 10

      if (parameter_tag.eq.'*Factor') then
       read(20,*) gammax,gammay
!        print *,parameter_tag
!        print *,gammax,gammay
       goto 10
      endif

      if (parameter_tag.eq.'*Rotation') then
       read(20,*) theta
!        print *,parameter_tag
!        print *,theta
       goto 10
      endif

      if (parameter_tag.eq.'*InputRotationCentre') then
       read(20,*) x0in,y0in
!        print *,parameter_tag
!        print *,x0in,y0in
       goto 10
      endif

      if (parameter_tag.eq.'*OutputRotationCentre') then
       read(20,*) x0out,y0out
!        print *,parameter_tag
!        print *,x0out,y0out
       goto 10
      endif

      if (parameter_tag.eq.'*xForwardCoefficients') then
       read(20,*) ((A(i,j),j=1,7-i),i=1,6)
!        print *,parameter_tag
!        do i=1,6
!        do j=1,7-i
!        print *,A(i,j)
!        enddo !j
!        enddo !i
       goto 10
      endif


      if (parameter_tag.eq.'*yForwardCoefficients') then
       read(20,*) ((B(i,j),j=1,7-i),i=1,6)
!        print *,parameter_tag
!        do i=1,6
!        do j=1,7-i
!        print *,B(i,j)
!        enddo !j
!        enddo !i
       goto 10
      endif
 
      goto 10
 
100      close(20)

       theta = theta/57.29578 !Convert to rad


1000  continue

! Convert input angles [''] to [deg]

!      xin=dble(theta_x)/3600.+x0in !Origin V2,V3 added to input!
!       yin=dble(theta_y)/3600.+y0in

      xin=dble(theta_x)/3600.
      yin=dble(theta_y)/3600.


! scale and rotate

      xp=gammax*( (xin-x0in)*cos(theta)+(yin-y0in)*sin(theta))+x0out
      yp=gammay*(-(xin-x0in)*sin(theta)+(yin-y0in)*cos(theta))+y0out

! distort

      xout=0.
      yout=0.
      do i=1,6
      do j=1,7-i
      xout=xout+A(i,j)*(xp**(i-1))*(yp**(j-1))
      yout=yout+B(i,j)*(xp**(i-1))*(yp**(j-1))
      enddo
      enddo

! convert stop position from [m] to [mm]

      xs=real(xout*1000.) !BD transforms give output in meters
      ys=real(yout*1000.)

      return
      end

!-------------------------------------------------------------------------------

      subroutine stop_to_ote_full(xs,ys,theta_x,theta_y)
!
!  Input:  xs,ys = position in field stop (mm)
!  Output: theta_x,theta_y = angular incidence angles on OTE (arcsec)
!
      use path_to_ref_files
      implicit none
      real xs,ys,theta_x,theta_y
      real*8 xin,yin,xp,yp,xout,yout
      real*8 gammax,gammay,theta,x0in,y0in,x0out,y0out
      real*8 C(6,6),D(6,6)
      integer i,j
      character*70 filename
      character*70 model_directory
      character*30 parameter_tag
      logical first
      save first,C,D,gammax,gammay,theta,x0in,y0in,x0out,y0out
      data first/.true./

      if (.not.first) goto 1000

      first=.false.

!       gammax = -2.296 ! f_ote=131 m when input [deg] output [m]
!       gammay = 2.2953  !??? Should there be a minus here ????
! 
!       theta = 0.2075/57.29578 !Convert to rad
! 
!       x0in = 0.10539
!       y0in = -0.11913000025
! 
!       x0out = -5.18289805611e-07
!       y0out = -1.92704532397e-09
! 
!       C(1,1) = 5.887781571463152e-12
!       C(1,2) = 0.001181872741672517
!       C(1,3) = 0.004381098668339454
!       C(1,4) = -0.000103026266160280
!       C(1,5) = -0.000162113939511591
!       C(1,6) = 0.002239212905820631
!       C(2,1) = 1.000010045405057
!       C(2,2) = -0.0198948306438118
!       C(2,3) = -0.01794430688622102
!       C(2,4) = 0.0007702083937083082
!       C(2,5) = 0.005557741128844107
!       C(3,1) = -0.01301049927765961
!       C(3,2) = 0.0007263338972753573
!       C(3,3) = 0.0002815687482939574
!       C(3,4) = -0.01831970054423948
!       C(4,1) = -0.0180620155086792
!       C(4,2) = 0.0007213604156713583
!       C(4,3) = -0.0121120145515583
!       C(5,1) = 0.0004148141635435355
!       C(5,2) = -0.009560875485648879
!       C(6,1) = 0.0005642557760050515
! 
! 
!       D(1,1) = 6.116356498319719e-10
!       D(1,2) = 1.00001443535109
!       D(1,3) = -0.01483102426265499
!       D(1,4) = -0.01798861538665233
!       D(1,5) = 0.0004181691144437354
!       D(1,6) = 0.0007218197164462481
!       D(2,1) = 0.001179828191001408
!       D(2,2) = -0.01747569411172255
!       D(2,3) = 0.0006738474640055891
!       D(2,4) = 0.0005905203096196509
!       D(2,5) = 0.001516029056689128
!       D(3,1) = 0.005008481986839269
!       D(3,2) = -0.01810347472067562
!       D(3,3) = 0.0001467468673671769
!       D(3,4) = 0.01135524127936449
!       D(4,1) = -9.635229102552406e-05
!       D(4,2) = 0.000606605966116013
!       D(4,3) = -0.002935907543077931
!       D(5,1) = -0.0002147771257013271
!       D(5,2) = 0.008775038955027625
!       D(6,1) = 3.639541372213451e-05


      filename=trim(adjustl(ref_path))//'model.conf'
      
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) model_directory
      close(20)
      
!       print *,model_directory
      
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/CoordTransform/OTE.pcf'

!       print *,filename

      open(20,file=filename,status='old',access='sequential',form='formatted')

10    read(20,*,end=100) parameter_tag
      if (parameter_tag(1:1).ne.'*') goto 10

      if (parameter_tag.eq.'*Factor') then
       read(20,*) gammax,gammay
!        print *,parameter_tag
!        print *,gammax,gammay
       goto 10
      endif

      if (parameter_tag.eq.'*Rotation') then
       read(20,*) theta
!        print *,parameter_tag
!        print *,theta
       goto 10
      endif

      if (parameter_tag.eq.'*InputRotationCentre') then
       read(20,*) x0in,y0in
!        print *,parameter_tag
!        print *,x0in,y0in
       goto 10
      endif

      if (parameter_tag.eq.'*OutputRotationCentre') then
       read(20,*) x0out,y0out
!        print *,parameter_tag
!        print *,x0out,y0out
       goto 10
      endif

      if (parameter_tag.eq.'*xBackwardCoefficients') then
       read(20,*) ((c(i,j),j=1,7-i),i=1,6)
!        print *,parameter_tag
!        do i=1,6
!        do j=1,7-i
!        print *,C(i,j)
!        enddo !j
!        enddo !i
       goto 10
      endif


      if (parameter_tag.eq.'*yBackwardCoefficients') then
       read(20,*) ((d(i,j),j=1,7-i),i=1,6)
!        print *,parameter_tag
!        do i=1,6
!        do j=1,7-i
!        print *,D(i,j)
!        enddo !j
!        enddo !i
       goto 10
      endif
 
      goto 10
 
100      close(20)

       theta = theta/57.29578 !Convert to rad

1000  continue

! Convert stop position from [mm] to [m]

      xout=dble(xs)/1000.
      yout=dble(ys)/1000.

! de-distort

      xp=0.
      yp=0.
      do i=1,6
      do j=1,7-i
      xp=xp+C(i,j)*(xout**(i-1))*(yout**(j-1))
      yp=yp+D(i,j)*(xout**(i-1))*(yout**(j-1))
      enddo
      enddo

! scale and rotate

      xin=(xp-x0out)/gammax*cos(theta)-(yp-y0out)*sin(theta)/gammay+x0in
      yin=(xp-x0out)/gammax*sin(theta)+(yp-y0out)*cos(theta)/gammay+y0in

! convert OTE angles from [deg] to ["]

!       theta_x=real((xin-x0in)*3600.) !Origin V2,V3 subtracted output!
!       theta_y=real((yin-y0in)*3600.)

      theta_x=real((xin)*3600.) 
      theta_y=real((yin)*3600.)


      return
      end


!-------------------------------------------------------------------------------


     subroutine read_V2V3_ref(V2_ref,V3_ref,AngV3)
     
     use path_to_ref_files
     implicit none
     real V2_ref,V3_ref,AngV3
     character*30 parameter_tag
     character*70 filename
 
      filename=trim(adjustl(ref_path))//'msa_ref_v2v3.ascii'

      open(19,file=filename,status='old',access='sequential',form='formatted')

10    read(19,*,end=100) parameter_tag
      if (parameter_tag(1:1).ne.'*') goto 10

      if (parameter_tag.eq.'*V2Ref') then
       read(19,*) V2_ref
       goto 10
      endif

      if (parameter_tag.eq.'*V3Ref') then
       read(19,*) V3_ref
       goto 10
      endif


      if (parameter_tag.eq.'*AngV3') then
       read(19,*) AngV3
       goto 10
      endif

 
      goto 10
 
100      close(19)

      end subroutine read_V2V3_ref
      
!---------------------------------------------------------------------------

