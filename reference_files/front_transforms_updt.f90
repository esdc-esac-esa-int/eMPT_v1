!-----------------------------------------------------------------
!
! Collection of subroutines defining FORE + OTE "front end" of NIRSPec Instrument Model
!
! Version 1.1 - 4. July 2017 - double precision introduced internally in transforms
! Version 1.2 - 16. July 2017 - reference wavelength of 2.5 micron introduced to
!                               msa_to_stop_clear and stop_to_msa_clear subroutines
! Version 1.3 - 15 July 2018 - cleaned up all compiler warnings
! Version 2.0 - 30 July 2021 - updated to read instrument parameters from python model
! Version 2.5 - 11 August 2021 - Modified to use module for common ref_path and remove incl_path
! Version 2.6 - 11 May 2022 - Modified to increase length of filename to 120 to accomodate long model directory names
! Version 2.7 - 24 June 2022 - Modified to handle filter-dependent FORE transform
! Version 2.8 - 17 September 2022 - msa_to_sky_rot modified to take actual phi_r value into account
!---------------------------------------------------------------------------

      subroutine sky_to_msa_bd(theta_x,theta_y,xm,ym)
!
!  Input: theta_x,theta_y = angular incidence angles on OTE (arcsec)
!  Output:  xm,ym = position in MSA plane (mm)
!
      implicit none
      real theta_x,theta_y,xm,ym
      real xs,ys

      call ote_to_stop_bd(theta_x,theta_y,xs,ys)
      call stop_to_msa_bd(xs,ys,xm,ym)

      return
      end

!---------------------------------------------------------------------------

      subroutine msa_to_sky_bd(xm,ym,theta_x,theta_y)
!
!  Input:  xm,ym = position in MSA plane (mm)
!  Output: theta_x,theta_y = angular incidence angles on OTE (arcsec)
!
      implicit none
      real xm,ym,theta_x,theta_y
      real xs,ys

      call msa_to_stop_bd(xm,ym,xs,ys)
      call stop_to_ote_bd(xs,ys,theta_x,theta_y)

      return
      end

!---------------------------------------------------------------------------

      subroutine msa_to_sky_rot_bd(xm,ym,theta_x,theta_y)
!
! outputs rotated V2V3 sky image with projected MSA horizontal/vertical
!
!  Input:  xm,ym = position in MSA plane (mm)
!  Output: theta_x,theta_y = angular incidence angles on OTE in rotated plane (arcsec)
!
      use derived_parameters, only : phi_r
      implicit none
      real xm,ym,theta_x,theta_y
      real xs,ys,theta_x0,theta_y0,theta0

      call load_derived_parameters

      call msa_to_stop_bd(xm,ym,xs,ys)
      call stop_to_ote_bd(xs,ys,theta_x,theta_y)

! rotate back

      theta_x0=theta_x
      theta_y0=theta_y
!       theta0 = -41.508/57.29578 !Convert to rad
       theta0 = -phi_r

      theta_x=theta_x0*cos(theta0)-theta_y0*sin(theta0)
      theta_y=theta_x0*sin(theta0)+theta_y0*cos(theta0)

! Linear fit for testing purposes

!      theta_x=-xm*2.55058885+1.39565850E-02
!      theta_y=ym*2.59485984+0.374654382


      return
      end

!---------------------------------------------------------------------------

      subroutine sky_to_msa_rot_bd(theta_x,theta_y,xm,ym)
!
!  Input: theta_x,theta_y = angular incidence angles on OTE (arcsec)
!  Output:  xm,ym = position in MSA plane (mm)
!
      implicit none
      real theta_x,theta_y,xm,ym
      real xs,ys,theta_x0,theta_y0,theta0


! counterrotate

      theta0 = 41.508/57.29578 !Convert to rad

      theta_x0=theta_x*cos(theta0)-theta_y*sin(theta0)
      theta_y0=theta_x*sin(theta0)+theta_y*cos(theta0)

      call ote_to_stop_bd(theta_x0,theta_y0,xs,ys)
      call stop_to_msa_bd(xs,ys,xm,ym)

      return
      end

!---------------------------------------------------------------------------

      subroutine stop_to_ote_bd(xs,ys,theta_x,theta_y)
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
       goto 10
      endif


      if (parameter_tag.eq.'*yBackwardCoefficients') then
       read(20,*) ((d(i,j),j=1,7-i),i=1,6)
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

      theta_x=real((xin-x0in)*3600.) !Origin V2,V3 subtracted output!
      theta_y=real((yin-y0in)*3600.)

!       theta_x=real((xin)*3600.) 
!       theta_y=real((yin)*3600.)


      return
      end


!-------------------------------------------------------------------------------

      subroutine ote_to_stop_bd(theta_x,theta_y,xs,ys)
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
       goto 10
      endif

      if (parameter_tag.eq.'*InputRotationCentre') then
       read(20,*) x0in,y0in
       goto 10
      endif

      if (parameter_tag.eq.'*OutputRotationCentre') then
       read(20,*) x0out,y0out
       goto 10
      endif

      if (parameter_tag.eq.'*xForwardCoefficients') then
       read(20,*) ((A(i,j),j=1,7-i),i=1,6)
       goto 10
      endif


      if (parameter_tag.eq.'*yForwardCoefficients') then
       read(20,*) ((B(i,j),j=1,7-i),i=1,6)
       goto 10
      endif
 
      goto 10
 
100      close(20)

       theta = theta/57.29578 !Convert to rad


1000  continue

! Convert input angles [''] to [deg]

      xin=dble(theta_x)/3600.+x0in !Origin V2,V3 added to input!
      yin=dble(theta_y)/3600.+y0in

!       xin=dble(theta_x)/3600.
!       yin=dble(theta_y)/3600.


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


!-----------------------------------------------------------------

      subroutine stop_to_msa_bd(xs,ys,xm,ym)
!
!  Input:  xs,ys = position in field stop (mm)
!  Output: xm,ym = position in MSA plane (mm)
!
      use path_to_ref_files
      use config_parameters, only : n_filt
      implicit none
      real xs,ys,xm,ym
      real*8 xin,yin,xp,yp,xout,yout
      real*8 gammax,gammay,theta,x0in,y0in,x0out,y0out
      real*8 A1(6,6),A2(6,6),B1(6,6),B2(6,6)
      real*8 lambda
      character*120 filename
      character*70 model_directory
      character*30 parameter_tag
      character*20 foretrans
      integer i,j
      logical first
      save first,A1,A2,B1,B2,gammax,gammay,theta,x0in,y0in,x0out,y0out,lambda
      data first/.true./

      if (.not.first) goto 1000

      first=.false.

      if (n_filt.eq.0) foretrans='Fore_F070LP.pcf'
      if (n_filt.eq.1) foretrans='Fore_F100LP.pcf'
      if (n_filt.eq.2) foretrans='Fore_F170LP.pcf'
      if (n_filt.eq.3) foretrans='Fore_F290LP.pcf'
      if (n_filt.eq.4) foretrans='Fore_CLEAR.pcf'
      if (n_filt.eq.5) foretrans='Fore_F110W.pcf'
      if (n_filt.eq.6) foretrans='Fore_F140X.pcf'

       filename=trim(adjustl(ref_path))//'model.conf'
      
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) model_directory
      close(20)     
      
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/CoordTransform/'//trim(adjustl(foretrans))

      open(20,file=filename,status='old',access='sequential',form='formatted')

10    read(20,*,end=100) parameter_tag
      if (parameter_tag(1:1).ne.'*') goto 10

      if (parameter_tag.eq.'*Wavelength') then
       read(20,*) lambda
       goto 10
      endif

      if (parameter_tag.eq.'*Factor') then
       read(20,*) gammax,gammay
       goto 10
      endif

      if (parameter_tag.eq.'*Rotation') then
       read(20,*) theta
       goto 10
      endif

      if (parameter_tag.eq.'*InputRotationCentre') then
       read(20,*) x0in,y0in
       goto 10
      endif

      if (parameter_tag.eq.'*OutputRotationCentre') then
       read(20,*) x0out,y0out
       goto 10
      endif

      if (parameter_tag.eq.'*xForwardCoefficients') then
       do i=1,6
       do j=1,7-i
       read(20,*) A1(i,j)
       enddo !j
       enddo !i
       do i=1,6
       do j=1,7-i
       read(20,*) A2(i,j)
       enddo !j
       enddo !i
       goto 10
      endif


      if (parameter_tag.eq.'*yForwardCoefficients') then
       do i=1,6
       do j=1,7-i
       read(20,*) B1(i,j)
       enddo !j
       enddo !i
       do i=1,6
       do j=1,7-i
       read(20,*) B2(i,j)
       enddo !j
       enddo !i
       goto 10
      endif
 
      goto 10
 
100      close(20)

      lambda=2.5E-6

      theta = theta/57.29578 !Convert to rad

1000  continue

! convert field stop position from [mm] to [m]

      xin=dble(xs)/1000.
      yin=dble(ys)/1000.

! scale and rotate

      xp=gammax*( (xin-x0in)*cos(theta)+(yin-y0in)*sin(theta))+x0out
      yp=gammay*(-(xin-x0in)*sin(theta)+(yin-y0in)*cos(theta))+y0out

! distort

      xout=0.
      yout=0.
      do i=1,6
      do j=1,7-i
      xout=xout+(A1(i,j)+A2(i,j)*lambda)*(xp**(i-1))*(yp**(j-1))
      yout=yout+(B1(i,j)+B2(i,j)*lambda)*(xp**(i-1))*(yp**(j-1))
      enddo
      enddo

! convert MSA position from [m] to [mm]

      xm=real(xout*1000.) !BD transforms give output in meters
      ym=real(yout*1000.)

      return
      end

!-----------------------------------------------------------------

      subroutine msa_to_stop_bd(xm,ym,xs,ys)
!
!  Input:  xm,ym = position in MSA plane (mm)
!  Output: xs,ys = position in field stop (mm)
!
      use path_to_ref_files
      use config_parameters, only : n_filt
      implicit none
      real xs,ys,xm,ym
      real*8 xin,yin,xp,yp,xout,yout
      real*8 gammax,gammay,theta,x0in,y0in,x0out,y0out
      real*8 C1(6,6),C2(6,6),D1(6,6),D2(6,6)
      real*8 lambda
      integer i,j
      character*120 filename
      character*70 model_directory
      character*30 parameter_tag
      character*20 foretrans
      logical first
      save first,C1,C2,D1,D2,gammax,gammay,theta,x0in,y0in,x0out,y0out,lambda
      data first/.true./

      if (.not.first) goto 1000

      first=.false.

      if (n_filt.eq.0) foretrans='Fore_F070LP.pcf'
      if (n_filt.eq.1) foretrans='Fore_F100LP.pcf'
      if (n_filt.eq.2) foretrans='Fore_F170LP.pcf'
      if (n_filt.eq.3) foretrans='Fore_F290LP.pcf'
      if (n_filt.eq.4) foretrans='Fore_CLEAR.pcf'
      if (n_filt.eq.5) foretrans='Fore_F110W.pcf'
      if (n_filt.eq.6) foretrans='Fore_F140X.pcf'

      filename=trim(adjustl(ref_path))//'model.conf'
            
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) model_directory
      close(20)
            
      filename=trim(adjustl(ref_path))//trim(adjustl(model_directory))//'/CoordTransform/'//trim(adjustl(foretrans))

      open(20,file=filename,status='old',access='sequential',form='formatted')

10    read(20,*,end=100) parameter_tag
      if (parameter_tag(1:1).ne.'*') goto 10

      if (parameter_tag.eq.'*Wavelength') then
       read(20,*) lambda
       goto 10
      endif

      if (parameter_tag.eq.'*Factor') then
       read(20,*) gammax,gammay
       goto 10
      endif

      if (parameter_tag.eq.'*Rotation') then
       read(20,*) theta
       goto 10
      endif

      if (parameter_tag.eq.'*InputRotationCentre') then
       read(20,*) x0in,y0in
       goto 10
      endif

      if (parameter_tag.eq.'*OutputRotationCentre') then
       read(20,*) x0out,y0out
       goto 10
      endif

      if (parameter_tag.eq.'*xBackwardCoefficients') then
       do i=1,6
       do j=1,7-i
       read(20,*) C1(i,j)
       enddo !j
       enddo !i
       do i=1,6
       do j=1,7-i
       read(20,*) C2(i,j)
       enddo !j
       enddo !i
       goto 10
      endif


      if (parameter_tag.eq.'*yBackwardCoefficients') then
       do i=1,6
       do j=1,7-i
       read(20,*) D1(i,j)
       enddo !j
       enddo !i
       do i=1,6
       do j=1,7-i
       read(20,*) D2(i,j)
       enddo !j
       enddo !i
       goto 10
      endif
 
      goto 10
 
100      close(20)

      theta = theta/57.29578 !Convert to rad

      lambda=2.5E-6


1000  continue

! convert msa position from [mm] to [m]

      xout=dble(xm)/1000.
      yout=dble(ym)/1000.

! de-distort

      xp=0.
      yp=0.
      do i=1,6
      do j=1,7-i
      xp=xp+(C1(i,j)+C2(i,j)*lambda)*(xout**(i-1))*(yout**(j-1))
      yp=yp+(D1(i,j)+D2(i,j)*lambda)*(xout**(i-1))*(yout**(j-1))
      enddo
      enddo

! scale and rotate

      xin=(xp-x0out)/gammax*cos(theta)-(yp-y0out)*sin(theta)/gammay+x0in
      yin=(xp-x0out)/gammax*sin(theta)+(yp-y0out)*cos(theta)/gammay+y0in

! convert stop position from [m] to [mm]

      xs=real(xin*1000.)
      ys=real(yin*1000.)

      return
      end

!-----------------------------------------------------------------

