!-----------------------------------------------------------------
!
! Collection of subroutines defining spectrograph "back end" of NIRSPec Instrument Model
!
! Version 1.0 - 6 September 2016
! Version 1.1 - 16 July 2017 - Minor corrections to fitted MSA shutter pitches introduced
! Version 1.2 - 15 July 2018 - Cleaned up all compiler warnings
! Version 2.0 - 30 July 2021 - Modified to read instrument parameters in from Python model
! Version 2.5 - 11 August 2021 - Modified to use module for common ref_path and remove incl_path
! Version 2.6 - 11 May 2022 - Modified to increase length of filename to 120 to accomodate long model directory names
! Version 2.6a - 1 June 2022 - Modified to change PRISM lambda1,lambda2 to 0.51 and 5.51 µm
! Version 2.9 - 24 June 2022 - CaF2 index of refraction model changed to read Sellmeier coefficients off model
!-----------------------------------------------------------------
      subroutine msa_to_tam_bd(xm,ym,theta_x,theta_y)
!
!  Input: xm,ym = position on MSA (mm)
!  Output: theta_x,theta_y = angular incidence angles on TAM (rad)
!
      use path_to_ref_files
      implicit none
      real xm,ym,theta_x,theta_y
      real xin,yin,xp,yp,xout,yout,norm,axt,ayt
      real gammax,gammay,theta,x0in,y0in,x0out,y0out
      real A(6,6),B(6,6)
      integer i,j
      character*120 filename
      character*70 model_directory
      character*30 parameter_tag
      logical first
      save first,A,B,gammax,gammay,theta,x0in,y0in,x0out,y0out
      data first/.true./

      if (.not.first) goto 1000

      first=.false.



      filename=trim(adjustl(ref_path))//'model.conf'
      
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) model_directory
      close(20)
      
!       print *,model_directory
      
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/CoordTransform/Collimator.pcf'

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
       do i=1,6
       do j=1,7-i
       read(20,*) A(i,j)
       enddo !j
       enddo !i
!        print *,parameter_tag
!        do i=1,6
!        do j=1,7-i
!        print *,A(i,j)
!        enddo !j
!        enddo !i
       goto 10
      endif


      if (parameter_tag.eq.'*yForwardCoefficients') then
       do i=1,6
       do j=1,7-i
       read(20,*) B(i,j)
       enddo !j
       enddo !i
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

      xin=xm/1000. !BD transforms have input in meters
      yin=ym/1000.

      xp=gammax*( (xin-x0in)*cos(theta)+(yin-y0in)*sin(theta))+x0out
      yp=gammay*(-(xin-x0in)*sin(theta)+(yin-y0in)*cos(theta))+y0out

      xout=0.
      yout=0.
      do i=1,6
      do j=1,7-i
      xout=xout+A(i,j)*(xp**(i-1))*(yp**(j-1))
      yout=yout+B(i,j)*(xp**(i-1))*(yp**(j-1))
      enddo
      enddo

! Convert "cosine" to direction cosines

      norm=sqrt(xout*xout+yout*yout+1.)

      axt=xout/norm
      ayt=yout/norm

! Convert to angles [rad]

      theta_y=asin(ayt)
      theta_x=asin(axt/cos(theta_y))

! Best fit paraxial model for testing purposes [mm] in [rad] out
! (rotation = 7.61670945E-03 rad)

!      theta_x=  xm*1.57225842E-03 + ym*1.19756678E-05 -2.04010587E-03
!      theta_y= -xm*1.15144103E-05 + ym*1.51170103E-03 -0.286156178

      return
      end


!-----------------------------------------------------------------
      subroutine tam_to_msa_bd(theta_x,theta_y,xm,ym)
!
!  Input: theta_x,theta_y = angular incidence angles on TAM (rad)
!  Ouput: xm,ym = position on MSA (mm)
!
      use path_to_ref_files
      implicit none
      real xm,ym,theta_x,theta_y
      real xin,yin,xp,yp,xout,yout,norm,axt,ayt
      real gammax,gammay,theta,x0in,y0in,x0out,y0out
      real C(6,6),D(6,6)
      integer i,j
      character*120 filename
      character*70 model_directory
      character*30 parameter_tag
      logical first
      save first,C,D,gammax,gammay,theta,x0in,y0in,x0out,y0out
      data first/.true./

      if (.not.first) goto 1000

      first=.false.



      filename=trim(adjustl(ref_path))//'model.conf'
      
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) model_directory
      close(20)
      
!       print *,model_directory
      
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/CoordTransform/Collimator.pcf'

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
       do i=1,6
       do j=1,7-i
       read(20,*) C(i,j)
       enddo !j
       enddo !i
!        print *,parameter_tag
!        do i=1,6
!        do j=1,7-i
!        print *,C(i,j)
!        enddo !j
!        enddo !i
       goto 10
      endif


      if (parameter_tag.eq.'*yBackwardCoefficients') then
       do i=1,6
       do j=1,7-i
       read(20,*) D(i,j)
       enddo !j
       enddo !i
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

! Convert input angle [rad] to direction cosines

      axt=sin(theta_x)*cos(theta_y)
      ayt=sin(theta_y)

! Convert direction cosines to "cosine"

      norm=sqrt(1-axt*axt-ayt*ayt)

      xout=axt/norm
      yout=ayt/norm

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

      xm=xin*1000. !BD transforms give output in meters
      ym=yin*1000.

      return
      end

!-----------------------------------------------------------------

      subroutine tam_to_fpa_bd(theta_x,theta_y,xf,yf)
!
!  Input: theta_x,theta_y = angular incidence angles on TAM (rad)
!  Output: xf,yf = position on FPA (mm)
!
      use path_to_ref_files
      implicit none
      real xf,yf,theta_x,theta_y
      real xin,yin,xp,yp,xout,yout,norm,axt,ayt
      real gammax,gammay,theta,x0in,y0in,x0out,y0out
      real A(6,6),B(6,6)
      integer i,j
      character*120 filename
      character*70 model_directory
      character*30 parameter_tag
      logical first
      save first,A,B,gammax,gammay,theta,x0in,y0in,x0out,y0out
      data first/.true./

      if (.not.first) goto 1000

      first=.false.



      filename=trim(adjustl(ref_path))//'model.conf'
      
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) model_directory
      close(20)
      
!       print *,model_directory
      
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/CoordTransform/Camera.pcf'

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
       do i=1,6
       do j=1,7-i
       read(20,*) A(i,j)
       enddo !j
       enddo !i
!        print *,parameter_tag
!        do i=1,6
!        do j=1,7-i
!        print *,A(i,j)
!        enddo !j
!        enddo !i
       goto 10
      endif


      if (parameter_tag.eq.'*yForwardCoefficients') then
       do i=1,6
       do j=1,7-i
       read(20,*) B(i,j)
       enddo !j
       enddo !i
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

! Convert input angle to direction cosines

      axt=sin(theta_x)*cos(theta_y)
      ayt=sin(theta_y)

! Convert direction cosines to "cosine"

      norm=sqrt(1-axt*axt-ayt*ayt)

      xin=axt/norm
      yin=ayt/norm

! rotate and scale

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

      xf=xout*1000. !BD transforms have output in meters
      yf=yout*1000.

! Paraxial Model for testing purposes [rad] in [mm] out
!    (-7.64358789E-03 rad rotation and offset)

!      xf = 283.724457*theta_x -2.16871524*theta_y + 5.48343211E-02
!      yf = 2.24546051*theta_x +293.764771*theta_y   -83.8842926

      return
      end

!-----------------------------------------------------------------

  
      subroutine fpa_to_tam_bd(xf,yf,theta_x,theta_y)
!
!  Input: xf,yf = position on FPA (mm)      
!  Output: theta_x,theta_y = angular incidence angles on TAM (rad)C
!
      use path_to_ref_files
      implicit none
      real xf,yf,theta_x,theta_y
      real xin,yin,xp,yp,xout,yout,norm,axt,ayt
      real gammax,gammay,theta,x0in,y0in,x0out,y0out
      real C(6,6),D(6,6)
      integer i,j
      character*120 filename
      character*70 model_directory
      character*30 parameter_tag
      logical first
      save first,C,D,gammax,gammay,theta,x0in,y0in,x0out,y0out
      data first/.true./

      if (.not.first) goto 1000

      first=.false.



      filename=trim(adjustl(ref_path))//'model.conf'
      
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) model_directory
      close(20)
      
!       print *,model_directory
      
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/CoordTransform/Camera.pcf'

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
       do i=1,6
       do j=1,7-i
       read(20,*) C(i,j)
       enddo !j
       enddo !i
!        print *,parameter_tag
!        do i=1,6
!        do j=1,7-i
!        print *,C(i,j)
!        enddo !j
!        enddo !i
       goto 10
      endif


      if (parameter_tag.eq.'*yBackwardCoefficients') then
       do i=1,6
       do j=1,7-i
       read(20,*) D(i,j)
       enddo !j
       enddo !i
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

! de-distort

      xout=xf/1000. !BD transforms in meters
      yout=yf/1000.
 
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
      
      
! Convert "cosine" to direction cosines

      norm=sqrt(xin*xin+yin*yin+1)
      
      axt=xin/norm
      ayt=yin/norm
            
! Convert to angles [rad]   
      
      theta_y=asin(ayt)
      theta_x=asin(axt/cos(theta_y))
                
      return
      end


!-----------------------------------------------------------------

      subroutine load_msa_dimensions_bd

      use path_to_ref_files
      implicit none
      character*120 filename
      character*70 model_directory
      character*30 parameter_tag
      integer i
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
!-------------------------



      filename=trim(adjustl(ref_path))//'MSA.msa.ascii'
      
!       print *,filename

      open(20,file=filename,status='old',access='sequential',form='formatted')

10    read(20,*,end=100) parameter_tag
      if (parameter_tag(1:1).ne.'*') goto 10

! Q1
      if (parameter_tag.eq.'*Q1xref') then
       read(20,*) x0(1)
!        print *,parameter_tag
!        print *,x0(1)
       goto 10
      endif

      if (parameter_tag.eq.'*Q1yref') then
       read(20,*) y0(1)
!        print *,parameter_tag
!        print *,y0(1)
       goto 10
      endif

      if (parameter_tag.eq.'*Q1rot') then
       read(20,*) rot(1)
!        print *,parameter_tag
!        print *,rot(1)
       goto 10
      endif

      if (parameter_tag.eq.'*Q1xpitch') then
       read(20,*) xpitch(1)
!        print *,parameter_tag
!        print *,xpitch(1)
       goto 10
      endif

      if (parameter_tag.eq.'*Q1ypitch') then
       read(20,*) ypitch(1)
!        print *,parameter_tag
!        print *,ypitch(1)
       goto 10
      endif

      if (parameter_tag.eq.'*Q1xopen') then
       read(20,*) xopen(1)
!        print *,parameter_tag
!        print *,xopen(1)
       goto 10
      endif

      if (parameter_tag.eq.'*Q1yopen') then
       read(20,*) yopen(1)
!        print *,parameter_tag
!        print *,yopen(1)
       goto 10
      endif

! Q2

     if (parameter_tag.eq.'*Q2xref') then
       read(20,*) x0(2)
!        print *,parameter_tag
!        print *,x0(2)
       goto 10
      endif

      if (parameter_tag.eq.'*Q2yref') then
       read(20,*) y0(2)
!        print *,parameter_tag
!        print *,y0(2)
       goto 10
      endif

      if (parameter_tag.eq.'*Q2rot') then
       read(20,*) rot(2)
!        print *,parameter_tag
!        print *,rot(2)
       goto 10
      endif

      if (parameter_tag.eq.'*Q2xpitch') then
       read(20,*) xpitch(2)
!        print *,parameter_tag
!        print *,xpitch(2)
       goto 10
      endif

      if (parameter_tag.eq.'*Q2ypitch') then
       read(20,*) ypitch(2)
!        print *,parameter_tag
!        print *,ypitch(2)
       goto 10
      endif

      if (parameter_tag.eq.'*Q2xopen') then
       read(20,*) xopen(2)
!        print *,parameter_tag
!        print *,xopen(2)
       goto 10
      endif

      if (parameter_tag.eq.'*Q2yopen') then
       read(20,*) yopen(2)
!        print *,parameter_tag
!        print *,yopen(2)
       goto 10
      endif

! Q3
 
      if (parameter_tag.eq.'*Q3xref') then
       read(20,*) x0(3)
!        print *,parameter_tag
!        print *,x0(3)
       goto 10
      endif

      if (parameter_tag.eq.'*Q3yref') then
       read(20,*) y0(3)
!        print *,parameter_tag
!        print *,y0(3)
       goto 10
      endif

      if (parameter_tag.eq.'*Q3rot') then
       read(20,*) rot(3)
!        print *,parameter_tag
!        print *,rot(3)
       goto 10
      endif

      if (parameter_tag.eq.'*Q3xpitch') then
       read(20,*) xpitch(3)
!        print *,parameter_tag
!        print *,xpitch(3)
       goto 10
      endif

      if (parameter_tag.eq.'*Q3ypitch') then
       read(20,*) ypitch(3)
!        print *,parameter_tag
!        print *,ypitch(3)
       goto 10
      endif

      if (parameter_tag.eq.'*Q3xopen') then
       read(20,*) xopen(3)
!        print *,parameter_tag
!        print *,xopen(3)
       goto 10
      endif

      if (parameter_tag.eq.'*Q3yopen') then
       read(20,*) yopen(3)
!        print *,parameter_tag
!        print *,yopen(3)
       goto 10
      endif

! Q4

      if (parameter_tag.eq.'*Q4xref') then
       read(20,*) x0(4)
!        print *,parameter_tag
!        print *,x0(4)
       goto 10
      endif

      if (parameter_tag.eq.'*Q4yref') then
       read(20,*) y0(4)
!        print *,parameter_tag
!        print *,y0(4)
       goto 10
      endif

      if (parameter_tag.eq.'*Q4rot') then
       read(20,*) rot(4)
!        print *,parameter_tag
!        print *,rot(4)
       goto 10
      endif

      if (parameter_tag.eq.'*Q4xpitch') then
       read(20,*) xpitch(4)
!        print *,parameter_tag
!        print *,xpitch(4)
       goto 10
      endif

      if (parameter_tag.eq.'*Q4ypitch') then
       read(20,*) ypitch(4)
!        print *,parameter_tag
!        print *,ypitch(4)
       goto 10
      endif

      if (parameter_tag.eq.'*Q4xopen') then
       read(20,*) xopen(4)
!        print *,parameter_tag
!        print *,xopen(4)
       goto 10
      endif

      if (parameter_tag.eq.'*Q4yopen') then
       read(20,*) yopen(4)
!        print *,parameter_tag
!        print *,yopen(4)
       goto 10
      endif

      goto 10
 
100      close(20)

      do i=1,4
      x0(i)=x0(i)*1000. !convert to mm
      y0(i)=y0(i)*1000.
      enddo

      mx=365
      my=171

      xbar=3.16 !mm
      ybar=3.14

      return
      end

!------------------------------------------------------------------------

      subroutine msa_to_shutter_bd(xm,ym,mq,i,j,rx,ry)

!  Input: xm,ym = position on MSA (mm)
!  Output: mq = MSA Quadrant [1,2,3,4]
!           i = shutter x coordinate [1,365]
!           j = shutter y coordinate [1,171]
!          rx = residual in x [fraction of MSA pitch in x]
!          ry = residual in y [fraction of MSA pitch in y]
!
! Asumes prior call of load_msa_dimensions
!
      implicit none
      real xm,ym,rx,ry
      integer i,j,mq
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
!-------------------------
       real xmm,ymm,xmr,ymr,corot,sirot


! determine MSA quadrant from x and y signs

      mq=0
      i=0
      j=0
      rx=10.
      ry=10.

      if ((xm.lt.0).and.(ym.lt.0)) then
         mq=1
      endif

      if ((xm.lt.0).and.(ym.gt.0)) then
         mq=2
      endif

      if ((xm.gt.0).and.(ym.lt.0)) then
         mq=3
      endif

      if ((xm.gt.0).and.(ym.gt.0)) then
         mq=4
      endif

      if (mq.eq.0) return

      xmm=xm-(x0(mq)-xpitch(mq)/2.)
      ymm=ym-(y0(mq)-ypitch(mq)/2.)
      corot=cos(rot(mq))
      sirot=sin(rot(mq))

      xmr=xmm*corot+ymm*sirot
      ymr=-xmm*sirot+ymm*corot

      i=int(xmr/xpitch(mq))+1
      j=int(ymr/ypitch(mq))+1

      if ((i.lt.1).or.(i.gt.mx).or.(j.lt.1).or.(j.gt.my)) then
      mq=0
      i=0
      j=0
      return
      endif

      rx=xmr/xpitch(mq)-float(i-1)-0.5
      ry=ymr/ypitch(mq)-float(j-1)-0.5

        return
        end

!-----------------------------------------------------------------

      subroutine shutter_to_msa_bd(mq,i,j,xm,ym)

!  Input:  mq = MSA Quadrant [1,2,3,4]
!           i = shutter x coordinate [1,365]
!           j = shutter y coordinate [1,171]
!  Output: xm,ym = position of center of shutter on MSA (mm)
!
! Asumes prior call of load_msa_dimensions
!
      implicit none
      real xm,ym
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

      xmr=(float(i)-0.5)*xpitch(mq)
      ymr=(float(j)-0.5)*ypitch(mq)
      corot=cos(-rot(mq))
      sirot=sin(-rot(mq))

      xmm=xmr*corot+ymr*sirot
      ymm=-xmr*sirot+ymr*corot

      xm=xmm+x0(mq)-xpitch(mq)/2.
      ym=ymm+y0(mq)-ypitch(mq)/2.

      return
      end

!----------------------------------------------------------------------

      subroutine shutter_to_msa_span_bd(mq,i,j,xmc,ymc,xmt,ymt,xmb,ymb)
!  Input:  mq = MSA Quadrant [1,2,3,4]
!           i = shutter x coordinate [1,365]
!           j = shutter y coordinate [1,171]
!  Output: xmc,ymc = position of center of shutter on MSA (mm)
!  Output: xmt,ymt = position of top of shutter on MSA (mm)
!  Output: xmb,ymb = position of bottom of shutter on MSA (mm)
!
! Asumes prior call of load_msa_dimensions
!
      implicit none
      real xmc,ymc,xmt,ymt,xmb,ymb
      integer mq,i,j
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
!-------------------------
       real xmm,ymm,xmr,ymr,corot,sirot

      if ((mq.lt.1).or.(mq.gt.4)) then
      write(*,*) mq,i,j
      stop 'mq out of bounds: 2'
      endif
      if ((i.lt.1).or.(i.gt.mx)) stop 'i out of bounds: 2'
      if ((j.lt.1).or.(j.gt.my)) then
      write(*,*) mq,i,j
      stop 'j out of bounds: 2'
      endif

      corot=cos(-rot(mq))
      sirot=sin(-rot(mq))

! center of shutter

      xmr=(float(i)-0.5)*xpitch(mq)
      ymr=(float(j)-0.5)*ypitch(mq)

      xmm=xmr*corot+ymr*sirot
      ymm=-xmr*sirot+ymr*corot

      xmc=xmm+x0(mq)-xpitch(mq)/2.
      ymc=ymm+y0(mq)-ypitch(mq)/2.

! top edge of shutter

      xmr=(float(i)-0.5)*xpitch(mq)
      ymr=(float(j))*ypitch(mq)

      xmm=xmr*corot+ymr*sirot
      ymm=-xmr*sirot+ymr*corot

      xmt=xmm+x0(mq)-xpitch(mq)/2.
      ymt=ymm+y0(mq)-ypitch(mq)/2.

! bottom edge of shutter

      xmr=(float(i)-0.5)*xpitch(mq)
      ymr=(float(j)-1.0)*ypitch(mq)

      xmm=xmr*corot+ymr*sirot
      ymm=-xmr*sirot+ymr*corot

      xmb=xmm+x0(mq)-xpitch(mq)/2.
      ymb=ymm+y0(mq)-ypitch(mq)/2.


      return
      end

!----------------------------------------------------------------------

      subroutine rottran(nq,m,xx,yy)
!  Input:   nq = MSA Quadrant [1,2,3,4]
!           xx = x coordinates within MSA Quadrant nq [1,m]
!           yy = y coordinates within Quadrant nq [1,m]
!  Output: xx,yy = rotated and shifted array in general MSA plane coordinates (mm)
!
! Asumes prior call of load_msa_dimensions
!
      implicit none
      integer nq,m,i
      real xx(m),yy(m)
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
!-------------------------
      real corot,sirot,xmr,ymr

      corot=cos(-rot(nq))
      sirot=sin(-rot(nq))

      do i=1,m

       xmr=xx(i)*corot+yy(i)*sirot
       ymr=-xx(i)*sirot+yy(i)*corot

       xx(i)=xmr+x0(nq)-xpitch(nq)/2.
       yy(i)=ymr+y0(nq)-ypitch(nq)/2.

      enddo

      return
      end


!-----------------------------------------------------------------

      subroutine load_fpa_dimensions_bd

      use path_to_ref_files
      implicit none
      real xff,yff,xfr,yfr,corot,sirot
      character*120 filename
      character*70 model_directory
      character*30 parameter_tag
      integer i
!--- FPA Parameters-----------
      integer mx,my
      real pitch_x(2),pitch_y(2),x0(2),y0(2),rot(2)
      common/fpa_dims/mx,my,pitch_x,pitch_y,x0,y0,rot
!----------------------------


      filename=trim(adjustl(ref_path))//'model.conf'
      
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) model_directory
      close(20)
      
!       print *,model_directory
      
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/FPA.fpa'
      
 
!       print *,filename

      open(20,file=filename,status='old',access='sequential',form='formatted')

10    read(20,*,end=100) parameter_tag
      if (parameter_tag(1:1).ne.'*') goto 10

! SCA 1

      if (parameter_tag.eq.'*SCA491_PitchX') then
       read(20,*) pitch_x(1)
!        print *,parameter_tag
!        print *,pitch_x(1)
       goto 10
      endif

      if (parameter_tag.eq.'*SCA491_PitchY') then
       read(20,*) pitch_y(1)
!        print *,parameter_tag
!        print *,pitch_y(1)
       goto 10
      endif

      if (parameter_tag.eq.'*SCA491_PosX') then
       read(20,*) x0(1)
!        print *,parameter_tag
!        print *,x0(1)
       goto 10
      endif

      if (parameter_tag.eq.'*SCA491_PosY') then
       read(20,*) y0(1)
!        print *,parameter_tag
!        print *,y0(1)
       goto 10
      endif

      if (parameter_tag.eq.'*SCA491_RotAngle') then
       read(20,*) rot(1)
!        print *,parameter_tag
!        print *,rot(1)
       goto 10
      endif


! SCA 2

      if (parameter_tag.eq.'*SCA492_PitchX') then
       read(20,*) pitch_x(2)
!        print *,parameter_tag
!        print *,pitch_x(2)
       goto 10
      endif

      if (parameter_tag.eq.'*SCA492_PitchY') then
       read(20,*) pitch_y(2)
!        print *,parameter_tag
!        print *,pitch_y(2)
       goto 10
      endif

      if (parameter_tag.eq.'*SCA492_PosX') then
       read(20,*) x0(2)
!        print *,parameter_tag
!        print *,x0(2)
       goto 10
      endif

      if (parameter_tag.eq.'*SCA492_PosY') then
       read(20,*) y0(2)
!        print *,parameter_tag
!        print *,y0(2)
       goto 10
      endif

      if (parameter_tag.eq.'*SCA492_RotAngle') then
       read(20,*) rot(2)
!        print *,parameter_tag
!        print *,rot(2)
       goto 10
      endif


      goto 10
 
100      close(20)

! convert from m to mm

     do i=1,2
     x0(i)=x0(i)*1000.
     y0(i)=y0(i)*1000.
     pitch_x(i)=pitch_x(i)*1000.
     pitch_y(i)=pitch_y(i)*1000.
     enddo

        mx=2048
        my=2048


! shift SCA 2 pixel origin to same convention as SCA 1

      xff=-float(mx-1)*pitch_x(2)
      yff=-float(my-1)*pitch_y(2)

      corot=cos(rot(2))
      sirot=sin(rot(2))

      xfr= xff*corot-yff*sirot+x0(2)
      yfr= xff*sirot+yff*corot+y0(2)

      x0(2)=xfr
      y0(2)=yfr

      return
      end


!------------------------------------------------------------------------

      subroutine fpa_to_pixel_bd(xf,yf,nsca,i,j,rx,ry)
!  Input: xf,yf = position on FPA (mm)
!  Output: nsca = SCA number [1,2]
!           i = pixel x coordinate [1,2048]
!           j = pixel y coordinate [1,2048]
!          rx = residual in x [fraction of FPA pixel in x]
!          ry = residual in y [fraction of FPA pixel in y]
!
! Assumes prior call of load_fpa_dimensions_bd
!
      implicit none
      real xf,yf,rx,ry
      integer nsca,i,j
      real xff,yff,xfr,yfr,corot,sirot
!--- FPA Parameters-----------
      integer mx,my
      real pitch_x(2),pitch_y(2),x0(2),y0(2),rot(2)
      common/fpa_dims/mx,my,pitch_x,pitch_y,x0,y0,rot
!----------------------------

! determine FPA quadrant from x sign

      nsca=0
      i=0
      j=0
      rx=10.
      ry=10.

      if (xf.lt.0) then
         nsca=1
      else
         nsca=2
      endif

      xff=xf-(x0(nsca)-pitch_x(nsca)/2.)
      yff=yf-(y0(nsca)-pitch_y(nsca)/2.)

      corot=cos(rot(nsca))
      sirot=sin(rot(nsca))

      xfr= xff*corot+yff*sirot
      yfr=-xff*sirot+yff*corot

      i=int(xfr/pitch_x(nsca))+1
      j=int(yfr/pitch_y(nsca))+1


      if ((i.lt.1).or.(i.gt.mx).or.(j.lt.1).or.(j.gt.my)) then
      nsca=0
      i=0
      j=0
      return
      endif

      rx=xfr/pitch_x(nsca)-float(i-1)-0.5
      ry=yfr/pitch_y(nsca)-float(j-1)-0.5

        return
        end

!-----------------------------------------------------------------

      subroutine pixel_to_fpa_bd(nsca,i,j,xf,yf)
!  Input:  nsca = SCA number [1,2]
!           i = pixel x coordinate [1,2048]
!           j = pixel y coordinate [1,2048]
!  Output: xf,yf = position of center of pixel on FPA (mm)
!
! Asumes prior call of load_fpa_dimensions_bd
!
      implicit none
      real xf,yf
      integer nsca,i,j
      real xff,yff,xfr,yfr,corot,sirot
!--- FPA Parameters-----------
      integer mx,my
      real pitch_x(2),pitch_y(2),x0(2),y0(2),rot(2)
      common/fpa_dims/mx,my,pitch_x,pitch_y,x0,y0,rot
!----------------------------

      if ((nsca.lt.1).or.(nsca.gt.2)) stop 'nsca out of bounds: 3'
      if ((i.lt.1).or.(i.gt.mx)) stop 'i out of bounds: 3'
      if ((j.lt.1).or.(j.gt.my)) stop 'j out of bounds: 3'


      xfr=(float(i)-0.5)*pitch_x(nsca)
      yfr=(float(j)-0.5)*pitch_y(nsca)

      corot=cos(-rot(nsca))
      sirot=sin(-rot(nsca))

      xff= xfr*corot+yfr*sirot
      yff=-xfr*sirot+yfr*corot

      xf=xff+x0(nsca)-pitch_x(nsca)/2.
      yf=yff+y0(nsca)-pitch_y(nsca)/2.

      return
      end

!----------------------------------------------------------------------

      subroutine msa_to_fpa_disp(n_disp,lambda,xm,ym,xf,yf)
! Input:
!  n_disp = 0 Mirror
!           1 G140M Band I   R=1000
!           2 G235M Band II  R=1000
!           3 G395M Band III R=1000
!           4 G140H Band I   R=2700
!           5 G235H Band II  R=2700
!           6 G395H Band III R=2700
!           7 Prism R=100
!        lambda wavelength [micron]
!        xm,ym position of shutter in MSA plane [mm]
! Output: xf,yf position of dispersed and focused light in FPA plane [mm]
!
! Requires single call with dummy parameters to initialize correctly
!
! 04/12/18: revised filter cutoffs to 0.70, 0.97. 1,66, 2.87 µm

      use path_to_ref_files
      implicit none
      character*120 filename
      character*70 model_directory
      character*30 parameter_tag
!---- Disperser parameters------
      integer n_disp,n_disp_prev
      real lambda1,lambda2
      character*6 gname
      common/disp/lambda1,lambda2,gname
      real lambda,xm,ym,xf,yf,theta_in,gamma_in,theta_out,gamma_out
      real ax0,ay0,az0,ax1,ay1,az1,ax2,ay2,az2,ax3,ay3,az3
      real ax_in,ay_in,az_in,ax_out,ay_out,az_out,axt,ayt,azt
      real csx,snx,csy,sny,csz,snz,csg,sng
      real tilt_x,tilt_y,tilt_z,tilt_g,d,m
      real alpha,n_T,nn,sna,csa,xm_prev,ym_prev
      save n_disp_prev,d,tilt_x,tilt_y,tilt_z,tilt_g,m,alpha,xm_prev,ym_prev,csx,snx,csy,sny,csz,snz,csg,sng,ax_in,ay_in,az_in
      data n_disp_prev/-10/
      data xm_prev,ym_prev/10000000.,10000000./


      if (n_disp.eq.n_disp_prev) goto 1000

      n_disp_prev=n_disp


      filename=trim(adjustl(ref_path))//'model.conf'
      
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) model_directory
      close(20)
      
!       print *,model_directory
                  
      if (n_disp.eq.0) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_MIRROR.dis'
      lambda1=0.5 !µm
      lambda2=5.3 !µm
      m=0.
      endif
      
      if (n_disp.eq.1) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_G140M.dis'
      lambda1=0.97 !µm
      lambda2=1.8 !µm
      m=-1.
      endif

      if (n_disp.eq.2) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_G235M.dis'
      lambda1=1.66 !µm
      lambda2=3.1 !µm
      m=-1.
      endif
     
      if (n_disp.eq.3) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_G395M.dis'
      lambda1=2.87 !µm
      lambda2=5.3 !µm
      m=-1.
      endif
      

      if (n_disp.eq.4) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_G140H.dis'
      lambda1=0.97 !µm
      lambda2=1.8 !µm
      m=-1.
      endif

      if (n_disp.eq.5) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_G235H.dis'
      lambda1=1.66 !µm
      lambda2=3.1 !µm
      m=-1.
      endif
     
      if (n_disp.eq.6) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_G395H.dis'
      lambda1=2.87 !µm
      lambda2=5.3 !µm
      m=-1.
      endif
   
      if (n_disp.eq.7) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_PRISM.pri'
      lambda1=0.51 !µm
      lambda2=5.51 !µm
      gname='PRISM'
      m=0.
      endif
 
 
!       print *,filename
      
!       print *,'lambda1:',lambda1
!       print *,'lambda2:',lambda2
!       print *,'order:',m

      open(20,file=filename,status='old',access='sequential',form='formatted')

10    read(20,*,end=100) parameter_tag
      if (parameter_tag(1:1).ne.'*') goto 10

      if (parameter_tag.eq.'*GRATINGNAME') then
       read(20,*) gname
!        print *,parameter_tag
!        print *,gname
       goto 10
      endif

      if (parameter_tag.eq.'*TILTY') then
       read(20,*) tilt_g
!        print *,parameter_tag
!        print *,tilt_g
       goto 10
      endif

      if (parameter_tag.eq.'*THETAX') then
       read(20,*) tilt_x
!        print *,parameter_tag
!        print *,tilt_x
       goto 10
      endif

      if (parameter_tag.eq.'*THETAY') then
       read(20,*) tilt_y
!        print *,parameter_tag
!        print *,tilt_y
       goto 10
      endif

      if (parameter_tag.eq.'*THETAZ') then
       read(20,*) tilt_z
!        print *,parameter_tag
!        print *,tilt_z
       goto 10
      endif

      if (parameter_tag.eq.'*ANGLE') then
       read(20,*) alpha
!        print *,parameter_tag
!        print *,alpha
       goto 10
      endif

      if (parameter_tag.eq.'*GROOVEDENSITY') then
       read(20,*) d
!        print *,parameter_tag
!        print *,d
       goto 10
      endif

      goto 10
 
100      close(20)


! convert groove density to 1/micron and disperser tilt angles to rads

      d=d/1.0E6 ![1/micron]
      tilt_g= -tilt_g/57.29578 ![rad]
      tilt_x= -tilt_x/57.29578/3600. ![rad]
      tilt_y= -tilt_y/57.29578/3600. ![rad]
      tilt_z= -tilt_z/57.29578/3600. ![rad]
      alpha=alpha/57.29578 ![rad]

! calculate relevant sines and cosines

      csx=cos(tilt_x)
      snx=sin(tilt_x)

      csy=cos(tilt_y)
      sny=sin(tilt_y)

      csz=cos(tilt_z)
      snz=sin(tilt_z)

      csg=cos(tilt_g)
      sng=sin(tilt_g)

1000  continue

! same xm,ym input location - only new lambda? Disperser tilt settings unchanged

      if ((xm.eq.xm_prev).and.(ym.eq.ym_prev)) goto 2000

      xm_prev=xm
      ym_prev=ym

      call msa_to_tam_bd(xm,ym,theta_in,gamma_in)

      axt=sin(theta_in)*cos(gamma_in)
      ayt=sin(gamma_in)
      azt=sqrt(1.0-axt*axt-ayt*ayt)

!      write(*,*) 'in',axt,ayt,azt


      ax1= axt
      ay1= ayt*csx - azt*snx
      az1= ayt*snx + azt*csx

!      write(*,*) '+x',ax1,ay1,az1


      ax2= ax1*csy + az1*sny
      ay2= ay1
      az2=-ax1*sny + az1*csy

!      write(*,*) '+y',ax2,ay2,az2


      ax3= ax2*csz - ay2*snz
      ay3= ax2*snz + ay2*csz
      az3= az2

!      write(*,*) '+z',ax3,ay3,az3


      ax_in= ax3*csg + az3*sng
      ay_in= ay3
      az_in=-ax3*sng + az3*csg


2000  continue

!      write(*,*) '+g',ax_in,ay_in,az_in

! initialize ax_out,ay_out & az_out for good measure to avoid warnings

      ax_out=0.
      ay_out=0.
      az_out=0.

      if ((n_disp.ge.1).and.(n_disp.le.6)) then

        ax_out=-ax_in-m*lambda*d
        ay_out=-ay_in
        az_out=sqrt(1.-ax_out*ax_out-ay_out*ay_out)

      endif


      if (n_disp.eq.7) then      !prism

        nn=n_T(lambda)           !index of refraction for wavelength

! Raytrace through double pass prism

!      write(*,*) 'in',ax_in,ay_in,az_in

        ax1=ax_in/nn
        ay1=ay_in/nn
        az1=sqrt(1-ax1*ax1-ay1*ay1)

!      write(*,*) 's1',ax1,ay1,az1

        csa=cos(alpha)
        sna=sin(alpha)

        ax2= ax1*csa - az1*sna
        ay2= ay1
        az2=ax1*sna + az1*csa

!      write(*,*) '+a',ax2,ay2,az2

        ax3=-ax2
        ay3=-ay2
        az3=az2

!      write(*,*) 'r1',ax3,ay3,az3

        ax1= ax3*csa + az3*sna
        ay1= ay3
        az1=-ax3*sna + az3*csa

!      write(*,*) '-a',ax1,ay1,az1

        ax_out=ax1*nn
        ay_out=ay1*nn
        az_out=sqrt(1.-ax_out*ax_out-ay_out*ay_out)

!      write(*,*) 's2',ax_out,ay_out,az_out

      endif


      if(n_disp.eq.0) then ! target acquisition mirror

       ax_out=ax_in
       ay_out=ay_in
       az_out=az_in

      endif

! rotate output beam back to angular system

!      write(*,*) 'di',ax_out,ay_out,az_out

      ax3= ax_out*csg - az_out*sng
      ay3= ay_out
      az3= ax_out*sng + az_out*csg

!      write(*,*) '-g',ax3,ay3,az3

      ax2= ax3*csz + ay3*snz
      ay2=-ax3*snz + ay3*csz
      az2= az3

!      write(*,*) '-z',ax2,ay2,az2

      ax1= ax2*csy - az2*sny
      ay1= ay2
      az1= ax2*sny + az2*csy

!      write(*,*) '-y',ax1,ay1,az1

      ax0= ax1
      ay0= ay1*csx + az1*snx
      az0=-ay1*snx + az1*csx

!      write(*,*) '-x',ax0,ay0,az0

! convert direction cosines to real angles

      gamma_out=asin(ay0)
      theta_out=asin(ax0/cos(gamma_out))

!      write(*,*) 'out',theta_out,gamma_out

! convert output beam direction to focussed spot in FPA plane

      call tam_to_fpa_bd(theta_out,gamma_out,xf,yf)

!      write(*,*) 'fpa',xf,yf

      return
      end

!--------------------------------------------------------------------

      real function n_T(lambda)
      use path_to_ref_files
      implicit none
      real lambda
      character*120 filename
      character*70 model_directory
      character*30 parameter_tag
      logical first
      real S1,S2,S3,l12,l22,l32,A1,A2,A3,l2
      save first,S1,S2,S3,l12,l22,l32
      data first/.true./

      if (.not.first) goto 1000

      first=.false.

      filename=trim(adjustl(ref_path))//'model.conf'
      
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) model_directory
      close(20)

      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_PRISM.pri'

      open(20,file=filename,status='old',access='sequential',form='formatted')

10    read(20,*,end=100) parameter_tag
      if (parameter_tag(1:1).ne.'*') goto 10

      if (parameter_tag.eq.'*COEFFORMULA') then
       read(20,*) S1
       read(20,*) l12
       read(20,*) S2
       read(20,*) l22
       read(20,*) S3
       read(20,*) l32
       goto 100
      endif

      goto 10

100   close (20)

1000 continue

      l2=lambda*lambda
      A1=S1*l2/(l2-l12)
      A2=S2*l2/(l2-l22)
      A3=S3*l2/(l2-l32)

      n_T=sqrt(1. + A1 + A2 + A3)

      return
      end

!----------------------------------------------------------------------
! 
!       real function n_T(lambda)
!       implicit none
!       real lambda,T,S10,S11,S12,S13,S14,S20,S21,S22,S23,S24,S30,S31,S32,S33,S34,S1,S2,S3,l1,l2,l3,A1,A2,A3
!       real l10,l11,l12,l13,l14,l20,l21,l22,l23,l24,l30,l31,l32,l33,l34
! 
! ! CaF2 Index of refraction Sellmeier fit - Leviton, Frey & Madison Proc SPIE 2008
! 
!       T=35. ! assumed operational temperature
! 
!       S10= 1.04834
!       S11=-2.21666E-4
!       S12=-6.73446E-6
!       S13= 1.50138E-8
!       S14=-2.77255E-11
! 
!       S20=-3.32723E-3
!       S21= 2.34683E-4
!       S22= 6.55744E-6
!       S23=-1.47028E-8
!       S24= 2.75023E-11
! 
!       S30= 3.72693
!       S31= 1.49844E-2
!       S32=-1.47511E-4
!       S33= 5.54293E-7
!       S34=-7.17298E-10
! 
! 
!       l10= 7.94375E-2
!       l11=-2.20758E-4
!       l12= 2.07862E-6
!       l13=-9.60254E-9
!       l14= 1.31401E-11
! 
!       l20= 0.258039
!       l21=-2.12833E-3
!       l22= 1.20393E-5
!       l23=-3.06973E-8
!       l24= 2.79793E-11
! 
!       l30= 34.0169
!       l31= 6.26867E-2
!       l32=-6.14541E-4
!       l33= 2.31517E-6
!       l34=-2.99638E-9
! 
!       l1=l10 + l11*T + l12*T**2 + l13*T**3 + l14*T**4
!       l2=l20 + l21*T + l22*T**2 + l23*T**3 + l24*T**4
!       l3=l30 + l31*T + l32*T**2 + l33*T**3 + l34*T**4
! 
! 
!       S1=S10 + S11*T + S12*T**2 + S13*T**3 + S14*T**4
!       S2=S20 + S21*T + S22*T**2 + S23*T**3 + S24*T**4
!       S3=S30 + S31*T + S32*T**2 + S33*T**3 + S34*T**4
! 
! 
!       A1=S1*lambda/(lambda+l1)*lambda/(lambda-l1)
!       A2=S2*lambda/(lambda+l2)*lambda/(lambda-l2)
!       A3=S3*lambda/(lambda+l3)*lambda/(lambda-l3)
! 
! 
!       n_T=sqrt(1. + A1 + A2 + A3)
!       
!       return
!       end function n_T
! 
! 
! ---------------------------------------------------------
! 
!       real function nn_T(lambda)
!       implicit none
!       real lambda,T,S10,S11,S12,S13,S14,S20,S21,S22,S23,S24,S30,S31,S32,S33,S34,S1,S2,S3,l1,l2,l3,A1,A2,A3
! 
! ! CaF2 Index of refraction Sellmeier fit - W. J. Topf, Optical Engineering, 34(5) 1369-1373 (1995)
! 
!       T=35. ! assumed operational temperature
! 
!       S10= 0.5815944
!       S11= 7.75880E-5
!       S12=-8.02380E-7
!       S13= 1.68964E-9
!       S14=-1.35186E-12
! 
!       S20= 0.4622821
!       S21=-6.01264E-5
!       S22= 5.95505E-7
!       S23=-1.30859E-9
!       S24= 1.06935E-12
! 
!       S30= 3.8729462
!       S31= 5.33254E-4
!       S32=-4.77551E-7
!       S33= 1.23826E-8
!       S34=-1.10489E-11
! 
! 
!       l1=0.050263605
!       l2=0.1003909
!       l3=34.649040
! 
!       S1=S10 + S11*T + S12*T**2 + S13*T**3 + S14*T**4
!       S2=S20 + S21*T + S22*T**2 + S23*T**3 + S24*T**4
!       S3=S30 + S31*T + S32*T**2 + S33*T**3 + S34*T**4
! 
! 
!       A1=S1*lambda/(lambda+l1)*lambda/(lambda-l1)
!       A2=S2*lambda/(lambda+l2)*lambda/(lambda-l2)
!       A3=S3*lambda/(lambda+l3)*lambda/(lambda-l3)
! 
! 
!       nn_T=sqrt(1. + A1 + A2 + A3)
! 
!       return
!       end
! 
! ---------------------------------------------------------

      subroutine load_msa_dimensions_bd_ideal
      implicit none
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims_ideal/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
!-------------------------
      real av,xx,yy
      integer i

!  regularize and then rotate approach

! MSA Metrology

      mx=365
      my=171


      x0(1)=-0.04268648784229399*1000.
      y0(1)=-0.04185273674364475*1000.
      rot(1)=-0.00036267844990266 !rad
      xpitch(1)=104.97/1000. !mm
      ypitch(1)=203.95/1000.
      xopen(1)=0.076 !mm
      yopen(1)=0.175

      x0(2)=-0.04280705897012871*1000.
      y0(2)=0.007269539733788312*1000.
      rot(2)=-0.00025189433508284
      xpitch(2)=104.98/1000. !mm
      ypitch(2)=203.79/1000.
      xopen(2)=0.076 !mm
      yopen(2)=0.175

      x0(3)=0.00459114910151604*1000.
      y0(3)=-0.04185427700300028*1000.
      rot(3)=-0.00017227917516437
      xpitch(3)=104.98/1000. !mm
      ypitch(3)=203.96/1000.
      xopen(3)=0.076 !mm
      yopen(3)=0.175

      x0(4)=0.004601540416117963*1000.
      y0(4)=0.007260364673386541*1000.
      rot(4)=-0.00015238322996341
      xpitch(4)=104.98/1000. !mm
      ypitch(4)=203.78/1000.
      xopen(4)=0.076 !mm
      yopen(4)=0.175

! Regularized array for testing purposes

      av=(x0(1)+x0(2))/2.+(73.0697937-64.2787476)/2000.
      x0(1)=av
      x0(2)=av

      av=(x0(3)+x0(4))/2.+(7.91898108-0.527828157)/2000.
      x0(3)=av
      x0(4)=av

      av=(y0(1)+y0(3))/2.+(-8.44610691 + 11.4249506)/2000.
      y0(1)=av
      y0(3)=av

      av=(y0(2)+y0(4))/2.+(-14.4054279-6.45742750)/2000.
      y0(2)=av
      y0(4)=av



      av=(xpitch(1)+xpitch(2)+xpitch(3)+xpitch(4))/4.
      xpitch(1)=av
      xpitch(2)=av
      xpitch(3)=av
      xpitch(4)=av

      av=(ypitch(1)+ypitch(2)+ypitch(3)+ypitch(4))/4.
      ypitch(1)=av
      ypitch(2)=av
      ypitch(3)=av
      ypitch(4)=av

      av=(rot(1)+rot(2)+rot(3)+rot(4))/4.
      av=-70./3600./57.29578
!      write(*,*) 'MSA Zero Clocking Angle [rad]: ',av
      rot(1)=av
      rot(2)=av
      rot(3)=av
      rot(4)=av

      do i=1,4
      xx=   x0(i)*cos(av)-y0(i)*sin(av)
      yy=   x0(i)*sin(av)+y0(i)*cos(av)
      x0(i)=xx
      y0(i)=yy
      enddo


      xbar=3.16 !mm
      ybar=3.14

      return
      end

!------------------------------------------------------------------------

      subroutine shutter_to_msa_bd_ideal(mq,i,j,xm,ym)
!  Input:  mq = MSA Quadrant [1,2,3,4]
!           i = shutter x coordinate [1,365]
!           j = shutter y coordinate [1,171]
!  Output: xm,ym = position of center of shutter on MSA (mm)
!
! Asumes prior call of load_msa_dimensions_ideal
!
      implicit none
      real xm,ym
      integer mq,i,j
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims_ideal/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
!-------------------------
       real xmm,ymm,xmr,ymr,corot,sirot

      if ((mq.lt.1).or.(mq.gt.4)) stop 'mq out of bounds: 4'
      if ((i.lt.1).or.(i.gt.mx)) stop 'i out of bounds: 4'
      if ((j.lt.1).or.(j.gt.my)) stop 'j out of bounds: 4'


      xmr=(float(i)-0.5)*xpitch(mq)
      ymr=(float(j)-0.5)*ypitch(mq)
      corot=cos(-rot(mq))
      sirot=sin(-rot(mq))

      xmm=xmr*corot+ymr*sirot
      ymm=-xmr*sirot+ymr*corot

      xm=xmm+x0(mq)-xpitch(mq)/2.
      ym=ymm+y0(mq)-ypitch(mq)/2.

      return
      end

!----------------------------------------------------------

      real function msa_gap(mq)
! routine used for calculating intra-quadrant horizontal distances between shutters
! including vertical gap between quadrants
      implicit none
      integer mq
      msa_gap=0.
      if ((mq.ge.3).or.(mq.ge.4)) msa_gap=1.
      return
      end


!-------------------------------------------------------------------


        SUBROUTINE LININT(X,Y,NPTS,XP,YP)
! Linear interpolation routine
        IMPLICIT NONE
        REAL X,Y,XP,YP,X1,X2,Y1,Y2,A0,A1
        INTEGER NPTS,J,J1,J2
        DIMENSION X(NPTS),Y(NPTS)
      IF(XP.LE.X(2))THEN
         J1=1
         J2=2
      ELSE IF(XP.GE.X(NPTS-1))THEN
         J1=NPTS-1
         J2=NPTS
      ELSE
         J=1
         DO WHILE((XP.GE.X(J)).AND.(J.LT.NPTS))
            J=J+1
         END DO
         IF(J.EQ.1)THEN
            J1=J
            J2=J+1
         ELSE IF(J.EQ.NPTS)THEN
            J1=J-1
            J2=J
         ELSE
            J1=J-1
            J2=J
         ENDIF
      ENDIF

      X1=X(J1)
      Y1=(Y(J1))
      X2=X(J2)
      Y2=(Y(J2))

      A1= (Y1-Y2)/(X1-X2)
      A0=  Y1-X1*A1
      YP=(A0+A1*XP)
        RETURN
        END
!--------------------------------------------------------

      subroutine msa_to_fpa_disp2(n,n_disp,lambda,xm,ym,xf,yf)
! Version with spectral order parseable
! Input:
!  n=0,-1,-2 spectral order
!  n_disp = 0 Mirror
!           1 G140M Band I   R=1000
!           2 G235M Band II  R=1000
!           3 G395M Band III R=1000
!           4 G140H Band I   R=2700
!           5 G235H Band II  R=2700
!           6 G395H Band III R=2700
!           7 Prism R=100
!        lambda wavelength [micron]
!        xm,ym position of shutter in MSA plane [mm]
! Output: xf,yf position of dispersed and focused light in FPA plane [mm]
!
! Requires single call with dummy parameters to initialize correctly
!
      use path_to_ref_files
      implicit none
      character*120 filename
      character*70 model_directory
      character*30 parameter_tag
!---- Disperser parameters------
      integer n_disp,n_disp_prev,n
      real lambda1,lambda2
      character*6 gname
      common/disp/lambda1,lambda2,gname
      real lambda,xm,ym,xf,yf,theta_in,gamma_in,theta_out,gamma_out
      real ax0,ay0,az0,ax1,ay1,az1,ax2,ay2,az2,ax3,ay3,az3
      real ax_in,ay_in,az_in,ax_out,ay_out,az_out,axt,ayt,azt
      real csx,snx,csy,sny,csz,snz,csg,sng
      real tilt_x,tilt_y,tilt_z,tilt_g,d,m
      real alpha,n_T,nn,sna,csa,xm_prev,ym_prev
      save m,n_disp_prev,d,tilt_x,tilt_y,tilt_z,tilt_g,alpha,xm_prev,ym_prev,csx,snx,csy,sny,csz,snz,csg,sng,ax_in,ay_in,az_in
      data n_disp_prev/-10/
      data xm_prev,ym_prev/10000000.,10000000./

      if (n_disp.eq.n_disp_prev) goto 1000

      n_disp_prev=n_disp



      filename=trim(adjustl(ref_path))//'model.conf'
      
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) model_directory
      close(20)
      
!       print *,model_directory
                  
      if (n_disp.eq.0) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_MIRROR.dis'
      lambda1=0.5 !µm
      lambda2=5.3 !µm
      m=0.
      endif
      
      if (n_disp.eq.1) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_G140M.dis'
      lambda1=0.97 !µm
      lambda2=1.8 !µm
      m=-1.
      endif

      if (n_disp.eq.2) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_G235M.dis'
      lambda1=1.66 !µm
      lambda2=3.1 !µm
      m=-1.
      endif
     
      if (n_disp.eq.3) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_G395M.dis'
      lambda1=2.87 !µm
      lambda2=5.3 !µm
      m=-1.
      endif
      

      if (n_disp.eq.4) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_G140H.dis'
      lambda1=0.97 !µm
      lambda2=1.8 !µm
      m=-1.
      endif

      if (n_disp.eq.5) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_G235H.dis'
      lambda1=1.66 !µm
      lambda2=3.1 !µm
      m=-1.
      endif
     
      if (n_disp.eq.6) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_G395H.dis'
      lambda1=2.87 !µm
      lambda2=5.3 !µm
      m=-1.
      endif
   
      if (n_disp.eq.7) then
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/Description/disperser_PRISM.pri'
      lambda1=0.51 !µm
      lambda2=5.51 !µm
      gname='PRISM'
      m=0.
      endif
 
 
      print *,filename
      
      print *,'lambda1:',lambda1
      print *,'lambda2:',lambda2
      print *,'order:',m

      open(20,file=filename,status='old',access='sequential',form='formatted')

10    read(20,*,end=100) parameter_tag
      if (parameter_tag(1:1).ne.'*') goto 10

      if (parameter_tag.eq.'*GRATINGNAME') then
       read(20,*) gname
       print *,parameter_tag
       print *,gname
       goto 10
      endif

      if (parameter_tag.eq.'*TILTY') then
       read(20,*) tilt_g
       print *,parameter_tag
       print *,tilt_g
       goto 10
      endif

      if (parameter_tag.eq.'*THETAX') then
       read(20,*) tilt_x
       print *,parameter_tag
       print *,tilt_x
       goto 10
      endif

      if (parameter_tag.eq.'*THETAY') then
       read(20,*) tilt_y
       print *,parameter_tag
       print *,tilt_y
       goto 10
      endif

      if (parameter_tag.eq.'*THETAZ') then
       read(20,*) tilt_z
       print *,parameter_tag
       print *,tilt_z
       goto 10
      endif

      if (parameter_tag.eq.'*ANGLE') then
       read(20,*) alpha
       print *,parameter_tag
       print *,alpha
       goto 10
      endif

      if (parameter_tag.eq.'*GROOVEDENSITY') then
       read(20,*) d
       print *,parameter_tag
       print *,d
       goto 10
      endif

      goto 10
 
100      close(20)



! convert groove density to 1/micron and disperser tilt angles to rads

      d=d/1.0E6 ![1/micron]
      tilt_g= -tilt_g/57.29578 ![rad]
      tilt_x= -tilt_x/57.29578/3600. ![rad]
      tilt_y= -tilt_y/57.29578/3600. ![rad]
      tilt_z= -tilt_z/57.29578/3600. ![rad]
      alpha=alpha/57.29578 ![rad]

      m=float(n)

! calculate relevant sines and cosines

      csx=cos(tilt_x)
      snx=sin(tilt_x)

      csy=cos(tilt_y)
      sny=sin(tilt_y)

      csz=cos(tilt_z)
      snz=sin(tilt_z)

      csg=cos(tilt_g)
      sng=sin(tilt_g)

1000  continue

! same xm,ym input location - only new lambda? Disperser tilt settings unchanged

      if ((xm.eq.xm_prev).and.(ym.eq.ym_prev)) goto 2000

      xm_prev=xm
      ym_prev=ym

      call msa_to_tam_bd(xm,ym,theta_in,gamma_in)

      axt=sin(theta_in)*cos(gamma_in)
      ayt=sin(gamma_in)
      azt=sqrt(1.0-axt*axt-ayt*ayt)

!      write(*,*) 'in',axt,ayt,azt


      ax1= axt
      ay1= ayt*csx - azt*snx
      az1= ayt*snx + azt*csx

!      write(*,*) '+x',ax1,ay1,az1


      ax2= ax1*csy + az1*sny
      ay2= ay1
      az2=-ax1*sny + az1*csy

!      write(*,*) '+y',ax2,ay2,az2


      ax3= ax2*csz - ay2*snz
      ay3= ax2*snz + ay2*csz
      az3= az2

!      write(*,*) '+z',ax3,ay3,az3


      ax_in= ax3*csg + az3*sng
      ay_in= ay3
      az_in=-ax3*sng + az3*csg


2000  continue

!      write(*,*) '+g',ax_in,ay_in,az_in

! initialize ax_out,ay_out & az_out for good measure to avoid warnings

      ax_out=0.
      ay_out=0.
      az_out=0.


      if ((n_disp.ge.1).and.(n_disp.le.6)) then

        ax_out=-ax_in-m*lambda*d
        ay_out=-ay_in
        az_out=sqrt(1.-ax_out*ax_out-ay_out*ay_out)

      endif


      if (n_disp.eq.7) then      !prism

        nn=n_T(lambda)           !index of refraction for wavelength

! Raytrace through double pass prism

!      write(*,*) 'in',ax_in,ay_in,az_in

        ax1=ax_in/nn
        ay1=ay_in/nn
        az1=sqrt(1-ax1*ax1-ay1*ay1)

!      write(*,*) 's1',ax1,ay1,az1

        csa=cos(alpha)
        sna=sin(alpha)

        ax2= ax1*csa - az1*sna
        ay2= ay1
        az2=ax1*sna + az1*csa

!      write(*,*) '+a',ax2,ay2,az2

        ax3=-ax2
        ay3=-ay2
        az3=az2

!      write(*,*) 'r1',ax3,ay3,az3

        ax1= ax3*csa + az3*sna
        ay1= ay3
        az1=-ax3*sna + az3*csa

!      write(*,*) '-a',ax1,ay1,az1

        ax_out=ax1*nn
        ay_out=ay1*nn
        az_out=sqrt(1.-ax_out*ax_out-ay_out*ay_out)

!      write(*,*) 's2',ax_out,ay_out,az_out

      endif


      if(n_disp.eq.0) then ! target acquisition mirror

       ax_out=ax_in
       ay_out=ay_in
       az_out=az_in

      endif

! rotate output beam back to angular system

!      write(*,*) 'di',ax_out,ay_out,az_out

      ax3= ax_out*csg - az_out*sng
      ay3= ay_out
      az3= ax_out*sng + az_out*csg

!      write(*,*) '-g',ax3,ay3,az3

      ax2= ax3*csz + ay3*snz
      ay2=-ax3*snz + ay3*csz
      az2= az3

!      write(*,*) '-z',ax2,ay2,az2

      ax1= ax2*csy - az2*sny
      ay1= ay2
      az1= ax2*sny + az2*csy

!      write(*,*) '-y',ax1,ay1,az1

      ax0= ax1
      ay0= ay1*csx + az1*snx
      az0=-ay1*snx + az1*csx

!      write(*,*) '-x',ax0,ay0,az0

! convert direction cosines to real angles

      gamma_out=asin(ay0)
      theta_out=asin(ax0/cos(gamma_out))

!      write(*,*) 'out',theta_out,gamma_out

! convert output beam direction to focussed spot in FPA plane

      call tam_to_fpa_bd(theta_out,gamma_out,xf,yf)

!      write(*,*) 'fpa',xf,yf

      return
      end

!--------------------------------------------------------------------
