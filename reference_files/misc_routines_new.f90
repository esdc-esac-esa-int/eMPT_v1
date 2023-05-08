! 
! Collection of sundry mixed subroutines employed by the eMPT
! Version 1.0 19 March 2022
! 
! NB!!! - Special general release version purged of references to fitsio
!         29 September 2022
!-----------------------------------------------------------------------

      subroutine load_derived_parameters

      use path_to_ref_files
      use derived_parameters
      implicit none
      character*80 filename
      character*30 parameter_tag
      logical,save :: first
      data first/.true./

      if (.not.first) return
      
      first=.false.

      filename=trim(adjustl(ref_path))//'derived_parameters.ascii'

      open(20,file=filename,status='old',access='sequential',form='formatted')

10    read(20,*,end=100) parameter_tag
      if (parameter_tag(1:1).ne.'*') goto 10

      if (parameter_tag.eq.'*av_xpitch') then
       read(20,*) av_xpitch
!        print *,parameter_tag,av_xpitch
       goto 10
      endif

      if (parameter_tag.eq.'*av_ypitch') then
       read(20,*) av_ypitch
!        print *,parameter_tag,av_ypitch
       goto 10
      endif

      if (parameter_tag.eq.'*av_xopen') then
       read(20,*) av_xopen
!        print *,parameter_tag,av_xopen
       goto 10
      endif

      if (parameter_tag.eq.'*av_yopen') then
       read(20,*) av_yopen
!        print *,parameter_tag,av_yopen
       goto 10
      endif

      if (parameter_tag.eq.'*phi') then
       read(20,*) phi
!        print *,parameter_tag,phi
       goto 10
      endif

      if (parameter_tag.eq.'*phi_r') then
       read(20,*) phi_r
!        print *,parameter_tag,phi_r
       goto 10
      endif

      if (parameter_tag.eq.'*xm_ref') then
       read(20,*) xm_ref
!        print *,parameter_tag,xm_ref
       goto 10
      endif

      if (parameter_tag.eq.'*ym_ref') then
       read(20,*) ym_ref
!        print *,parameter_tag,ym_ref
       goto 10
      endif


      if (parameter_tag.eq.'*avxgap') then
       read(20,*) avxgap
!        print *,parameter_tag,avxgap
       goto 10
      endif

      goto 10

100   close(20)
      
      return

      end subroutine load_derived_parameters

!-----------------------------------------------------------------

      subroutine m_list_size_read()
      use config_parameters, only : n_dither
      use m_list
      implicit none
      integer nm1,nm2,nm3,ip,i
      integer nspt

! open m_list file to work out maximum number of m entries per pointing

      open(9,file=m_list_file,status='old',access='sequential',form='formatted')
 
      read(9,*) nspt
      
      n_max_m=-9999

      if (n_dither.eq.1) then
     
      do ip=1,nspt
       read(9,*) 
       read(9,*) 
       read(9,*) 
       read(9,*) nm1
       n_max_m=max(n_max_m,nm1)
         do i=1,nm1
           read(9,*)
         enddo
      enddo !ip

     close(9)
     
     return
     
     endif !n_dither=1

     
     if (n_dither.eq.2) then
     
      do ip=1,nspt
       read(9,*) 
       read(9,*) 
       read(9,*) 
       read(9,*) nm1
       n_max_m=max(n_max_m,nm1)
         do i=1,nm1
           read(9,*)
         enddo
       read(9,*) nm2
       n_max_m=max(n_max_m,nm2)
         do i=1,nm2
           read(9,*)
         enddo
      enddo !ip

     close(9)
     
     return
     
     endif !n_dither=2


     if (n_dither.eq.3) then
     
      do ip=1,nspt
       read(9,*) 
       read(9,*) 
       read(9,*) 
       read(9,*) nm1
       n_max_m=max(n_max_m,nm1)
         do i=1,nm1
           read(9,*)
         enddo
       read(9,*) nm2
       n_max_m=max(n_max_m,nm2)
         do i=1,nm2
           read(9,*)
         enddo
       read(9,*) nm3
       n_max_m=max(n_max_m,nm3)
         do i=1,nm3
           read(9,*)
         enddo
      enddo !ip

     close(9)
     
     return
     
     endif !n_dither=3


      end subroutine m_list_size_read

!-----------------------------------------------------------------

      subroutine catalog_read()

! Routine to read target catalog

      use target_cat
      use config_parameters, only : catfile
      implicit none
      integer i,idum,pri
      real*8 xdum,ydum

! Open target catalog

      open(9,file=catfile,status='old',access='sequential', form='formatted')

      call skip_hashes(9)

! Determine number of targets and number of priority classes in catalog
 
      nclass=0
      ntarg=0
10    read(9,*,end=20) idum,xdum,ydum,pri
      ntarg=ntarg+1
      nclass=max(nclass,pri)
      goto 10
20    close(9)


! allocate target arrays

      allocate (id_t(ntarg),ra_t(ntarg),dec_t(ntarg),pri_t(ntarg))

! Open target catalog second time and read in 

      open(9,file=catfile,status='old',access='sequential', form='formatted')

      call skip_hashes(9)      
      do i=1,ntarg
      read(9,*)  id_t(i),ra_t(i),dec_t(i),pri_t(i)
      enddo

      close(9)

      return

      end subroutine catalog_read

!-----------------------------------------------------------------
 
      subroutine pointings_read()
      
! Read IPA-Output pointings from configuration file
 
      use config_parameters, only : confname
      use pointings_and_groups
      character*4 start_marker,end_marker,strng
      integer idum
      
      start_marker='#IL#'      
      end_marker='#IE#'      


! open configuration file first time to work out how many pointings are listed

      open(9,file=confname,status='old',access='sequential',form='formatted')
  
      call find_marker(9,start_marker)
       
      n_p=0
      
10    read(9,'(a4)') strng
      if (strng.eq.end_marker) goto 20
      if (strng(1:1).ne.'#') then
      backspace(9)
      n_p=n_p+1
      read(9,*)
      endif
      goto 10
20    close(9)


! allocate pointing arrays

     allocate (no_p(n_p),ra_p(n_p),dec_p(n_p),pa_ap_p(n_p))
 
! open configuration file a second time to read in pointings
 
       open(9,file=confname,status='old',access='sequential',form='formatted')
  
      call find_marker(9,start_marker)
         
      n_p=0
      
100   read(9,'(a4)') strng
      if (strng.eq.end_marker) goto 200
      if (strng(1:1).ne.'#') then
      backspace(9)
      n_p=n_p+1
      read(9,*) no_p(n_p),idum,ra_p(n_p),dec_p(n_p),pa_ap_p(n_p)
      endif
      goto 100
200   close(9) 
     
      return

      end subroutine pointings_read

!-----------------------------------------------------------------

      subroutine k_list_read()
      use pointings_and_groups, only : n_p
      use k_list
      implicit none
      integer npointings,n_k,ip,k,idum


! open k_list_raw file first time to work out maximum number of k entried per pointing

      open(9,file=k_list_file,status='old',access='sequential',form='formatted')
      read(9,*) npointings
      
      if (npointings.ne.n_p) then
      write (*,*) ' Error: Discrepancy in number of pointings in configuration and k_list files' 
      stop
      endif
 
      n_max_k=-9999
      
      do ip=1,n_p
       read(9,*) idum,idum
       read(9,*) n_k
       n_max_k=max(n_max_k,n_k)
         do k=1,n_k
           read(9,*)
         enddo
      enddo

     close(9)

! allocate k_list arrays now that n_max_k is known

     allocate (nk(n_p),id(n_p,n_max_k),id_cat(n_p,n_max_k))
     allocate (pri_k(n_p,n_max_k),kt(n_p,n_max_k),it(n_p,n_max_k),jt(n_p,n_max_k),rx(n_p,n_max_k),ry(n_p,n_max_k))
 
! open k_list_raw file a second time to read in entries
 
      open(9,file=k_list_file,status='old',access='sequential',form='formatted')
      read(9,*) idum
      
      do ip=1,n_p
       read(9,*) idum,idum
       read(9,*) nk(ip)
         do k=1,nk(ip)
           read(9,*) id(ip,k),id_cat(ip,k),pri_k(ip,k),kt(ip,k),it(ip,k),jt(ip,k),rx(ip,k),ry(ip,k)
         enddo
      enddo

     close(9)
   
      return

      end subroutine k_list_read

!-----------------------------------------------------------------

      subroutine k_list_mod_read()
      use pointings_and_groups, only : n_p
      use k_list
      implicit none
!       integer npointings,n_k,ip,k,idum
      integer n_k,ip,k,idum


! open k_list_rmod file first time to work out maximum number of k entries per pointing

      open(9,file=k_list_mod_file,status='old',access='sequential',form='formatted')
 
      read(9,*) n_p
      
!       if (npointings.ne.n_p) then
!       write (*,*) ' Error: Discrepancy in number of pointings in configuration and k_list files' 
!       stop
!       endif
      
       
      n_max_k=-9999
      
      do ip=1,n_p
       read(9,*) idum,idum
       read(9,*) n_k
       n_max_k=max(n_max_k,n_k)
         do k=1,n_k
           read(9,*)
         enddo
      enddo

     close(9)

! allocate k_list arrays now that n_max_k is known

     allocate (nk(n_p),id(n_p,n_max_k),id_cat(n_p,n_max_k))
     allocate (pri_k(n_p,n_max_k),kt(n_p,n_max_k),it(n_p,n_max_k),jt(n_p,n_max_k),rx(n_p,n_max_k),ry(n_p,n_max_k))
 
! open k_list_mod file a second time to read in entries
 
      open(9,file=k_list_mod_file,status='old',access='sequential',form='formatted')
      read(9,*) idum
      
      do ip=1,n_p
       read(9,*) idum,idum
       read(9,*) nk(ip)
         do k=1,nk(ip)
           read(9,*) id(ip,k),id_cat(ip,k),pri_k(ip,k),kt(ip,k),it(ip,k),jt(ip,k),rx(ip,k),ry(ip,k)
         enddo
      enddo

     close(9)
   
      return

      end subroutine k_list_mod_read

!-----------------------------------------------------------------

      subroutine dva_calc(angle_to_target,dva_mag)
      implicit none
      real angle_to_target,beta,dva_mag
      parameter (beta=30.0/2.99792458E5) !km/s
      
      dva_mag=sqrt(1.-beta*beta)/(1.+beta*cos(angle_to_target/57.295780)) !dva on-sky magnification factor
      
      return
      
      end subroutine dva_calc

!-------------------------------------------------------------

      subroutine make_spoiler_masks3(ra_p0,dec_p0,pa_ap0,c0,c1,c2)
      use config_parameters, only : dva_mag
      use derived_parameters, only : xm_ref,ym_ref,ax_refv2,ay_refv3,phi
      use target_cat

      implicit none
      real*8 ra_p0,dec_p0,ra_p1,dec_p1,ra_p2,dec_p2
      real pa_ap0,pa_v30
      real ax_refra_p,ay_refdec_p
      real xt(ntarg),yt(ntarg),xtv2(ntarg),ytv3(ntarg)
      real xm,ym
      real xx,yy,theta_x,theta_y,av_ypitch_msa
      integer it
      real rx,ry
      integer ks,is,js
      integer c0(4,365,171),c1(4,365,171),c2(4,365,171)
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4), xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen, xbar,ybar


      pa_v30=pa_ap0-(180.-phi)
      if (pa_v30.lt.0.) pa_v30=pa_v30+360.
      if (pa_v30.gt.360.) pa_v30=pa_v30-360.

! work out pointing reference in RA,Dec for given roll

      call roll_v3(pa_v30,ax_refv2,ay_refv3,ax_refra_p,ay_refdec_p) !RA,Dec

      call load_msa_dimensions_bd
      call load_fpa_dimensions_bd

       av_ypitch_msa=(ypitch(1)+ypitch(2)+ypitch(3)+ypitch(4))/4.


! calculate two nodded pointing coordinates

      call msa_to_sky_bd(xm_ref,ym_ref+av_ypitch_msa,xx,yy)
      call roll_v3(pa_v30,xx,yy,theta_x,theta_y)

      call coordtang_single(ra_p0,dec_p0,theta_x-ax_refra_p,theta_y-ay_refdec_p,ra_p2,dec_p2)

      call msa_to_sky_bd(xm_ref,ym_ref-av_ypitch_msa,xx,yy)
      call roll_v3(pa_v30,xx,yy,theta_x,theta_y)

      call coordtang_single(ra_p0,dec_p0,theta_x-ax_refra_p,theta_y-ay_refdec_p,ra_p1,dec_p1)


! null occupancy arrays


      c0=0
      c1=0
      c2=0


! fill out occupancy array for central dither pointing first

      call tangcoord(ra_p0,dec_p0,ntarg,ra_t,dec_t,xt,yt)

      do it=1,ntarg

      xt(it)=xt(it)*dva_mag+ax_refra_p
      yt(it)=yt(it)*dva_mag+ay_refdec_P
      
! Rotate to V2,V3 system orientation

      call roll_v3(-pa_v30,xt(it),yt(it),xtv2(it),ytv3(it))

! May now be projected to MSA.

       call sky_to_msa_bd(xtv2(it),ytv3(it),xm,ym)

! Identify which quadrant and shutter target falls into (ks=0 ~ misses MSA entirely)

      call msa_to_shutter_bd(xm,ym,ks,is,js,rx,ry)

! flag target presence in occupancy array 0

      if ((ks.ge.1).and.(ks.le.4)) then
      c0(ks,is,js)=c0(ks,is,js)+1
      endif

      enddo !nt

! fill out occupancy array for dither pointing 1

      call tangcoord(ra_p1,dec_p1,ntarg,ra_t,dec_t,xt,yt)

      do it=1,ntarg

      xt(it)=xt(it)*dva_mag+ax_refra_p
      yt(it)=yt(it)*dva_mag+ay_refdec_p

! Rotate to V2,V3 system orientation

      call roll_v3(-pa_v30,xt(it),yt(it),xtv2(it),ytv3(it))

! May now be projected to MSA.

       call sky_to_msa_bd(xtv2(it),ytv3(it),xm,ym)

! Identify which quadrant and shutter target falls into (ks=0 ~ misses MSA entirely)

      call msa_to_shutter_bd(xm,ym,ks,is,js,rx,ry)

! flag target presence in occupancy array 1

      if ((ks.ge.1).and.(ks.le.4)) then
      c1(ks,is,js)=c1(ks,is,js)+1
      endif

      enddo !nt

! fill out occupancy array for dither pointing 2

      call tangcoord(ra_p2,dec_p2,ntarg,ra_t,dec_t,xt,yt)

      do it=1,ntarg

      xt(it)=xt(it)*dva_mag+ax_refra_p
      yt(it)=yt(it)*dva_mag+ay_refdec_P

! Rotate to V2,V3 system orientation

      call roll_v3(-pa_v30,xt(it),yt(it),xtv2(it),ytv3(it))

! May now be projected to MSA.

       call sky_to_msa_bd(xtv2(it),ytv3(it),xm,ym)

! Identify which quadrant and shutter target falls into (ks=0 ~ misses MSA entirely)

      call msa_to_shutter_bd(xm,ym,ks,is,js,rx,ry)

! flag target presence in occupancy array 2

      if ((ks.ge.1).and.(ks.le.4)) then
      c2(ks,is,js)=c2(ks,is,js)+1
      endif

      enddo !nt

      return

      end subroutine make_spoiler_masks3

!-----------------------------------------------------------------
! 
!       subroutine read_fits_image(filename)
!       use fits_image
!       implicit none
!       integer status,unit,readwrite,blocksize,naxes(2),nfound
!       integer group
!       real nullval
!       logical anynull
!       character filename*100,comment
!       real cdelt1,cdelt2
!       
! !  The STATUS parameter must always be initialized.
!       status=0
! 
! !  Get an unused Logical Unit Number to use to open the FITS file.
!       call ftgiou(unit,status)
! 
! !  Open the FITS file 
!       readwrite=0
!       call ftopen(unit,filename,readwrite,blocksize,status)
! 
! !  Determine the size of the image.
!       call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
!             
! 
! !  Check that it found both NAXIS1 and NAXIS2 keywords.
!       if (nfound .ne. 2)then
!           print *,'READIMAGE failed to read the NAXISn keywords.'
!           return
!        end if
! 
!       nx=naxes(1)
!       ny=naxes(2)
!       
! 
! ! allocate image array
! 
!       allocate (image(nx,ny))
! 
! ! work out what type of wcs the file has 
! 
! ! read common always present parameters
! 
!       call ftgkyd(unit,'CRPIX1',crpix1,comment,status)
!       call ftgkyd(unit,'CRPIX2',crpix2,comment,status)
!       call ftgkyd(unit,'CRVAL1',crval1,comment,status)
!       call ftgkyd(unit,'CRVAL2',crval2,comment,status)
! 
! 
! ! Normal STScI-style matrix?
! 
!       call ftgkyd(unit,'CD1_1',cd1_1,comment,status)
!       
!       if (status.gt.0) then
!       status=0
!       goto 100 !wcs does not conform to CD1_1 convention
!       endif
! 
!       call ftgkyd(unit,'CD2_2',cd2_2,comment,status)
!       if (status.gt.0) then
!       status=0 
!       goto 100 !wcs does not conform to CD2_2 convention
!       endif
! 
!       call ftgkyd(unit,'CD1_2',cd1_2,comment,status)
!       if (status.gt.0) then
!       cd1_2=0. ! keyword not in header, set it
!       status=0
!       endif
!       
!       call ftgkyd(unit,'CD2_1',cd2_1,comment,status)
!       if (status.gt.0) then
!       cd2_1=0. ! keyword not in header, set it
!       status=0
!       endif
! 
!       goto 1000
! 
! 100   continue
! 
! ! Brant-style scaleless matrix?
! 
!       call ftgkyd(unit,'PC1_1',cd1_1,comment,status)
! 
!       if (status.gt.0) then
!       status=0
!       goto 200
!       endif
! 
!       call ftgkyd(unit,'PC1_2',cd1_2,comment,status)
!       call ftgkyd(unit,'PC2_1',cd2_1,comment,status)
!       call ftgkyd(unit,'PC2_2',cd2_2,comment,status)
! 
!       call ftgkye(unit,'CDELT1',cdelt1,comment,status)
!       call ftgkye(unit,'CDELT2',cdelt2,comment,status)
! 
! 
!       cd1_1=cd1_1*cdelt1
!       cd1_2=cd1_2*cdelt1
!       cd2_1=cd2_1*cdelt2
!       cd2_2=cd2_2*cdelt2
! 
!       goto 1000
! 
! 200 continue
! 
! ! Sandro-style no matrix at all?      
! 
!       call ftgkye(unit,'CDELT1',cdelt1,comment,status)
!       call ftgkye(unit,'CDELT2',cdelt2,comment,status)
! 
!       cd1_1=1.*cdelt1
!       cd1_2=0.*cdelt1
!       cd2_1=0.*cdelt2
!       cd2_2=1.*cdelt2
! 
! 1000 continue
! 
! 
! 
! !  Initialize variables
!       group=1
!       nullval=-999
!       
!        call ftg2de(unit,group,nullval,nx,nx,ny,image,anynull,status)
!        
!        
! !  The FITS file must always be closed before exiting the program. 
! !  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
!       call ftclos(unit, status)
!       call ftfiou(unit, status)
! 
! !  Check for any error, and if so print out error messages.
! !  The PRINTERROR subroutine is listed near the end of this file.
!       if (status .gt. 0) call printerror(status)
!       
!       return
!       
!       end subroutine read_fits_image
! 
!-----------------------------------------------------------------
! 
!       subroutine read_segmap(filename)
!       use fits_image
!       implicit none
!       integer status,unit,readwrite,blocksize,naxes(2),nfound
!       integer group
!       real nullval
!       logical anynull
!       character filename*100,comment
!       real cdelt1,cdelt2
!       
! !  The STATUS parameter must always be initialized.
!       status=0
! 
! !  Get an unused Logical Unit Number to use to open the FITS file.
!       call ftgiou(unit,status)
! 
! !  Open the FITS file 
!       readwrite=0
!       call ftopen(unit,filename,readwrite,blocksize,status)
! 
! !  Determine the size of the image.
!       call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
!             
! 
! !  Check that it found both NAXIS1 and NAXIS2 keywords.
!       if (nfound .ne. 2)then
!           print *,'READIMAGE failed to read the NAXISn keywords.'
!           return
!        end if
! 
!       nx=naxes(1)
!       ny=naxes(2)
!       
! 
! ! allocate image array
! 
!       allocate (segmap(nx,ny))
! 
! ! work out what type of wcs the file has 
! 
! ! read common always present parameters
! 
!       call ftgkyd(unit,'CRPIX1',crpix1,comment,status)
!       call ftgkyd(unit,'CRPIX2',crpix2,comment,status)
!       call ftgkyd(unit,'CRVAL1',crval1,comment,status)
!       call ftgkyd(unit,'CRVAL2',crval2,comment,status)
! 
! 
! ! Normal STScI-style matrix?
! 
!       call ftgkyd(unit,'CD1_1',cd1_1,comment,status)
!       
!       if (status.gt.0) then
!       status=0
!       goto 100 !wcs does not conform to CD1_1 convention
!       endif
! 
!       call ftgkyd(unit,'CD2_2',cd2_2,comment,status)
!       if (status.gt.0) then
!       status=0 
!       goto 100 !wcs does not conform to CD2_2 convention
!       endif
! 
!       call ftgkyd(unit,'CD1_2',cd1_2,comment,status)
!       if (status.gt.0) then
!       cd1_2=0. ! keyword not in header, set it
!       status=0
!       endif
!       
!       call ftgkyd(unit,'CD2_1',cd2_1,comment,status)
!       if (status.gt.0) then
!       cd2_1=0. ! keyword not in header, set it
!       status=0
!       endif
! 
!       goto 1000
! 
! 100   continue
! 
! ! Brant-style scaleless matrix?
! 
!       call ftgkyd(unit,'PC1_1',cd1_1,comment,status)
! 
!       if (status.gt.0) then
!       status=0
!       goto 200
!       endif
! 
!       call ftgkyd(unit,'PC1_2',cd1_2,comment,status)
!       call ftgkyd(unit,'PC2_1',cd2_1,comment,status)
!       call ftgkyd(unit,'PC2_2',cd2_2,comment,status)
! 
!       call ftgkye(unit,'CDELT1',cdelt1,comment,status)
!       call ftgkye(unit,'CDELT2',cdelt2,comment,status)
! 
! 
!       cd1_1=cd1_1*cdelt1
!       cd1_2=cd1_2*cdelt1
!       cd2_1=cd2_1*cdelt2
!       cd2_2=cd2_2*cdelt2
! 
!       goto 1000
! 
! 200 continue
! 
! ! Sandro-style no matrix at all?      
! 
!       call ftgkye(unit,'CDELT1',cdelt1,comment,status)
!       call ftgkye(unit,'CDELT2',cdelt2,comment,status)
! 
!       cd1_1=1.*cdelt1
!       cd1_2=0.*cdelt1
!       cd2_1=0.*cdelt2
!       cd2_2=1.*cdelt2
! 
! 1000 continue
! 
! 
! 
! !  Initialize variables
!       group=1
!       nullval=-999
!       
!        call ftg2dj(unit,group,nullval,nx,nx,ny,segmap,anynull,status)
!        
!        
! !  The FITS file must always be closed before exiting the program. 
! !  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
!       call ftclos(unit, status)
!       call ftfiou(unit, status)
! 
! !  Check for any error, and if so print out error messages.
! !  The PRINTERROR subroutine is listed near the end of this file.
!       if (status .gt. 0) call printerror(status)
!       
!       return
!       
!       end subroutine read_segmap
! 
!-----------------------------------------------------------------
! 
!       subroutine printerror(status)
! 
! !  This subroutine prints out the descriptive text corresponding to the
! !  error status value and prints out the contents of the internal
! !  error message stack generated by FITSIO whenever an error occurs.
! 
!       integer status
!       character errtext*30,errmessage*80
! 
! !  Check if status is OK (no error); if so, simply return
!       if (status .le. 0)return
! 
! !  The FTGERR subroutine returns a descriptive 30-character text string that
! !  corresponds to the integer error status number.  A complete list of all
! !  the error numbers can be found in the back of the FITSIO User's Guide.
!       call ftgerr(status,errtext)
!       print *,'FITSIO Error Status =',status,': ',errtext
! 
! !  FITSIO usually generates an internal stack of error messages whenever
! !  an error occurs.  These messages provide much more information on the
! !  cause of the problem than can be provided by the single integer error
! !  status value.  The FTGMSG subroutine retrieves the oldest message from
! !  the stack and shifts any remaining messages on the stack down one
! !  position.  FTGMSG is called repeatedly until a blank message is
! !  returned, which indicates that the stack is empty.  Each error message
! !  may be up to 80 characters in length.  Another subroutine, called
! !  FTCMSG, is available to simply clear the whole error message stack in
! !  cases where one is not interested in the contents.
!       call ftgmsg(errmessage)
!       do while (errmessage .ne. ' ')
!           print *,errmessage
!           call ftgmsg(errmessage)
!       end do
! 
!       return
! 
!       end subroutine printerror
! 
!------------------------------------------------------------------------

      subroutine tangcoord(cra,cdec,n,ra,dec,xt,yt)
!  Calculate xt,yt tangential coordinates [arcsec] from ra,dec [deg] and central cra,cdec [deg]
      integer n,i
      real*8 cra,cdec,ra(n),dec(n)
      real xt(n),yt(n)
      real*8 cra_r,cdec_r,r,xx,yy,dec_r,ra_r

! Scale from deg to rad

      cra_r=cra/57.2957795 !rad
      cdec_r=cdec/57.2957795


      do i=1,n

      ra_r=ra(i)/57.2957795
      dec_r=dec(i)/57.2957795

      r=sin(dec_r)*sin(cdec_r)+cos(dec_r)*cos(cdec_r)*cos(ra_r-cra_r)
      xx=(cos(dec_r)*sin(ra_r-cra_r))/r
      yy=(sin(dec_r)*cos(cdec_r)-cos(dec_r)*sin(cdec_r)*cos(ra_r-cra_r))/r

!  Scale from rad to arcsec

      xt(i)=real(206264.8062*xx)
      yt(i)=real(206264.8062*yy)

      enddo

      return
      end subroutine tangcoord

! --------------------------------------------------------------------

      subroutine tangcoord_single(cra,cdec,ra,dec,xt,yt)
!  Calculate xt,yt tangential coordinates [arcsec] from ra,dec [deg] and central cra,cdec [deg]
      real*8 cra,cdec,ra,dec
      real xt,yt
      real*8 cra_r,cdec_r,r,xx,yy,dec_r,ra_r

! Scale from deg to rad

      cra_r=cra/57.2957795 !rad
      cdec_r=cdec/57.2957795

      ra_r=ra/57.2957795
      dec_r=dec/57.2957795

      r=sin(dec_r)*sin(cdec_r)+cos(dec_r)*cos(cdec_r)*cos(ra_r-cra_r)
      xx=(cos(dec_r)*sin(ra_r-cra_r))/r
      yy=(sin(dec_r)*cos(cdec_r)-cos(dec_r)*sin(cdec_r)*cos(ra_r-cra_r))/r

!  Scale from rad to arcsec

      xt=real(206264.8062*xx)
      yt=real(206264.8062*yy)


      return
      end subroutine tangcoord_single

!--------------------------------------------------------------------------------

      subroutine coordtang_single(cra,cdec,xt,yt,ra,dec)
!  Calculate ra,dec [deg] from xt,yt tangential coordinates [arcsec] and central cra,cdec [deg]
      real*8 cra,cdec,ra,dec
      real xt,yt
      real*8 cra_r,cdec_r,xx,yy,dec_r,ra_r

! Scale from deg to rad

      cra_r=cra/57.2957795 !rad
      cdec_r=cdec/57.2957795


! Scale from arcsec to rad

      xx=dble(xt/206264.8062)
      yy=dble(yt/206264.8062)

      ra_r=cra_r+atan(xx/(cos(cdec_r)-yy*sin(cdec_r)))

      dec_r=asin((sin(cdec_r)+yy*cos(cdec_r))/sqrt(1+xx*xx+yy*yy))

! scale from rads to degrees

      ra=ra_r*57.2957795
      dec=dec_r*57.2957795



      return
      end subroutine coordtang_single

!--------------------------------------------------------------------------------

     subroutine roll_v3(pa_v3,x0,y0,xr,yr)

! Rotate MSA quadrants counterclockwise on sky by V3 Postion Angle (deg)

         implicit none
         real pa_v3,x0,y0,xr,yr
         real theta0

         theta0=pa_v3/57.2957795 !rad
         xr=x0*cos(theta0)+y0*sin(theta0)
         yr=-x0*sin(theta0)+y0*cos(theta0)

         return
         end subroutine roll_v3

!------------------------------------------------------------------------------------------

      subroutine xy_to_radec(x,y,ra,dec)
      use fits_image
      implicit none
      real*8 ra,dec
      real x,y
      real*8 cra_r,cdec_r,xx,yy,xr,yr,ra_r,dec_r

! Scale from deg to rad

      cra_r=crval1/57.2957795 !rad
      cdec_r=crval2/57.2957795


! Scale from deg to rad

      xx=(dble(x)-crpix1)/57.2957795
      yy=(dble(y)-crpix2)/57.2957795
      
      xr=(cd1_1*xx+cd1_2*yy)
      yr=(cd1_2*xx+cd2_2*yy)

      
      ra_r=cra_r+atan(xr/(cos(cdec_r)-yr*sin(cdec_r)))
      dec_r=asin((sin(cdec_r)+yr*cos(cdec_r))/sqrt(1+xr*xr+yr*yr))

! scale from rads to degrees

      ra=ra_r*57.2957795
      dec=dec_r*57.2957795

      return
      end subroutine xy_to_radec
      
!-----------------------------------------------------------------

      subroutine radec_to_xy(ra,dec,x,y)
      use fits_image
      implicit none
      real*8 ra,dec
      real x,y
      real*8 cra_r,cdec_r,r,xx,yy,dec_r,ra_r,xr,yr

! Scale from deg to rad

      cra_r=crval1/57.2957795 !rad
      cdec_r=crval2/57.2957795

      ra_r=ra/57.2957795
      dec_r=dec/57.2957795

      r=sin(dec_r)*sin(cdec_r)+cos(dec_r)*cos(cdec_r)*cos(ra_r-cra_r)
      xx=(cos(dec_r)*sin(ra_r-cra_r))/r
      yy=(sin(dec_r)*cos(cdec_r)-cos(dec_r)*sin(cdec_r)*cos(ra_r-cra_r))/r

!  Scale from rad to deg

      xx=xx*57.2957795
      yy=yy*57.2957795
      
      xr=( cd2_2*xx - cd1_2*yy)/(cd1_1*cd2_2-cd1_2*cd2_1) + crpix1
      yr=(-cd2_1*xx + cd1_1*yy)/(cd1_1*cd2_2-cd1_2*cd2_1) + crpix2
      
      x=real(xr)
      y=real(yr)
      
      return
      end subroutine radec_to_xy

!-----------------------------------------------------------------

      subroutine ij_to_radec(i,j,ra,dec)
      implicit none
      integer i,j
      real*8 ra,dec
      real x,y

      x=float(i)
      y=float(j)      

      call xy_to_radec(x,y,ra,dec)

      return
      end subroutine ij_to_radec
      
!-----------------------------------------------------------------

      subroutine radec_to_ij(ra,dec,i,j)
      implicit none
      integer i,j
      real*8 ra,dec
      real x,y

      call radec_to_xy(ra,dec,x,y)

      i=int(x+0.5)
      j=int(y+0.5)
 
      return
      end subroutine radec_to_ij
      
!-----------------------------------------------------------------

      subroutine ij_to_radec_box(i,j,ra,dec)
      implicit none
      integer i,j
      real*8 ra(5),dec(5)
      real x,y

      x=float(i)-0.5
      y=float(j)-0.5      

      call xy_to_radec(x,y,ra(5),dec(5))

      x=float(i)+0.5
      y=float(j)-0.5      

      call xy_to_radec(x,y,ra(4),dec(4))

      x=float(i)+0.5
      y=float(j)+0.5      

      call xy_to_radec(x,y,ra(3),dec(3))

      x=float(i)-0.5
      y=float(j)+0.5      

      call xy_to_radec(x,y,ra(2),dec(2))
      
      ra(1)=ra(5)
      dec(1)=dec(5)

!         call swap_order_d(ra,dec,5)

      return
      end subroutine ij_to_radec_box

!-----------------------------------------------------------------
 
    subroutine overlap_area_frac(xr,yr,xc,yc,frac,no)
    implicit none
    double precision xr(5),yr(5),xc(5),yc(5),xo(10),yo(10),a,a0,frac
    integer no

     frac=0.
     
     no=-1
     
     call slh(xr,yr,5,xc,yc,5,xo,yo,no)
     
     if (no.le.0) return
     

     call poly_area(xo,yo,no,a)
     
 
     a0=1.
     if (a.gt.0) call poly_area(xr,yr,5,a0)
  
     frac=a/a0
           
     return
     
     end subroutine overlap_area_frac

! -------------------------------------------------------- !

  subroutine slh(xr,yr,nr,xc,yc,nc,xo,yo,no)
  
! xr,yr vertices of reference polygon
! xc,yc vertices of clipping polygon
! xo,yo vertices of projection of xr,yr onto xc,yc
  
! NB polygon vertices xr,yr and xc,yc must be listed in counter-clockwise order !
  
  implicit none
  integer nr
  integer nc
  double precision xr(nr),yr(nr),xc(nc),yc(nc),xw(nr+nc),yw(nr+nc),xo(nr+nc),yo(nr+nc)  
  double precision x1,y1,x2,y2
  integer nw,no
  integer i

  nw=nc
  xw(1:nw)=xc(1:nw)
  yw(1:nw)=yc(1:nw)
  
  do i=1,nr-1
   
   x1=xr(i)
   y1=yr(i)
   x2=xr(i+1)
   y2=yr(i+1)

   call edge_clip(nc,nr,xw,yw,nw,x1,y1,x2,y2,xo,yo,no)
  
   nw=no
   xw(1:nw)=xo(1:nw)
   yw(1:nw)=yo(1:nw)
   
  enddo
  
  return
  
  end subroutine slh
  
! -------------------------------------------------------- !
 
    subroutine edge_clip(nc,nr,xi,yi,ni,x1,y1,x2,y2,xo,yo,no)
    implicit none
    integer :: nr
    integer :: nc
    double precision xi(ni),yi(ni),xo(nc+nr),yo(nc+nr)  
    double precision x1,y1,x2,y2,x0,y0,xx1,yy1,xx2,yy2
    integer ni,no,i,c
    logical inside
    
    c = 0 ! counter for the output polygon
 
    do i=1,ni-1 ! for each edge i of poly

      xx1=xi(i)
      yy1=yi(i)
      xx2=xi(i+1)
      yy2=yi(i+1)
      
      if (inside(xx1,yy1,x1,y1,x2,y2)) then
        if (inside(xx2,yy2,x1,y1,x2,y2)) then
         c=c+1
         xo(c)=xx2
         yo(c)=yy2
         
        else
        
         call intersection(xx1,yy1,xx2,yy2,x1,y1,x2,y2,x0,y0)
         c=c+1
         xo(c)=x0
         yo(c)=y0
         end if
         
        else
 
       if (inside(xx2,yy2,x1,y1,x2,y2)) then
          call intersection(xx1,yy1,xx2,yy2,x1,y1,x2,y2,x0,y0)
          c=c+1
          xo(c)=x0
          yo(c)=y0
          c=c+1
          xo(c)=xx2
          yo(c)=yy2         
        end if
      end if
    end do
     
    if (c .gt. 0) then
       if ((xo(1).ne.xo(c)).or.(yo(1).ne.yo(c))) then
        c=c+1
        xo(c)=xo(1)
        yo(c)=yo(1)
      end if
    end if

    no = c
    
    return
    
  end subroutine edge_clip
  
! -------------------------------------------------------- !

  subroutine intersection(xx1,yy1,xx2,yy2,x1,y1,x2,y2,x0,y0)
  implicit none
  double precision xx1,yy1,xx2,yy2,x1,y1,x2,y2,x0,y0
  double precision vx1,vy1,vx2,vy2,xyx,xyy
  double precision a
  double precision cross_product
      
    vx1=xx2-xx1
    vy1=yy2-yy1
    
    vx2=x2-x1
    vy2=y2-y1
    
    if (cross_product(vx1,vy1,vx2,vy2).eq.0.d0) then

      xyx=x1-xx1
      xyy=y1-yy1
 
       if (cross_product(xyx,xyy,vx1,vy1).eq.0.d0) then
 
        x0=xx2
        y0=yy2
        
       endif
       
     else
     
     xyx=x1-xx1
     xyy=y1-yy1
     
     a = cross_product(xyx,xyy,vx2,vy2)/cross_product(vx1,vy1,vx2,vy2)
      
     if ((a.gt.1.d0).or.(a.lt.0)) then
        ! no intersection
      else
        x0=xx1+a*vx1
        y0=yy1+a*vy1        
      end if
    end if

    return
    
    end subroutine intersection
    
! -------------------------------------------------------- !

  function inside(px,py,x1,y1,x2,y2)
  implicit none
  double precision px,py,x1,y1,x2,y2,vx1,vy1,vx2,vy2
  logical inside
  double precision cross_product
   vx1=x2-x1
   vy1=y2-y1
   vx2=px-x1
   vy2=py-y1
    if (cross_product(vx1,vy1,vx2,vy2).ge.0.d0) then
     inside=.true.
    else
     inside=.false.
    endif
   end function inside

! -------------------------------------------------------- !

  function cross_product(x1,y1,x2,y2)
    implicit none
    double precision x1,y1,x2,y2
    double precision :: cross_product

    cross_product=x1*y2-y1*x2
 
  end function cross_product

! -------------------------------------------------------- !

   subroutine poly_area(x,y,nb,area) 
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-04  Time: 12:24:06

IMPLICIT NONE
INTEGER :: nb
DOUBLE PRECISION :: x(nb)
DOUBLE PRECISION :: y(nb)
DOUBLE PRECISION :: area

!*****************************************************************

!   GIVEN A SEQUENCE OF NB POINTS (X(I),Y(I)),  polyarea COMPUTES THE AREA
! BOUNDED BY THE CLOSED POLYGONAL CURVE WHICH PASSES THROUGH THE POINTS IN
! THE ORDER THAT THEY ARE INDEXED.  THE FINAL POINT OF THE CURVE IS ASSUMED
! TO BE THE FIRST POINT GIVEN.  THEREFORE, IT NEED NOT BE LISTED AT THE END
! OF X AND Y.  THE CURVE IS NOT REQUIRED TO BE SIMPLE.  e.g. It may cross over
! itself.

!*****************************************************************

INTEGER  :: i, n, nm1
DOUBLE PRECISION     :: a

n = nb
IF (x(1) == x(n) .AND. y(1) == y(n)) n = n - 1

SELECT CASE (n)
  CASE (:2)
     area = 0.0

  CASE (3)
     area = 0.5*((x(2) - x(1))*(y(3) - y(1)) - (x(3) - x(1))*(y(2) - y(1)))

  CASE DEFAULT
    nm1 = n - 1
    a = x(1)*(y(2) - y(n)) + x(n)*(y(1) - y(nm1))
    DO  i = 2, nm1
      a = a + x(i)*(y(i+1) - y(i-1))
    END DO
    area = 0.5*a
END SELECT
!    area=abs(area)
RETURN
 end subroutine poly_area

!-----------------------------------------------------------
 
      SUBROUTINE PNPOLY(PX,PY,XX,YY,N,INOUT)                            

!                                                                       
!        SUBROUTINE PNPOLY                                              
!                                                                       
!        PURPOSE                                                        
!           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
!                                                                       
!        USAGE                                                          
!           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     
!                                                                       
!        DESCRIPTION OF THE PARAMETERS                                  
!           PX      - X-COORDINATE OF POINT IN QUESTION.                
!           PY      - Y-COORDINATE OF POINT IN QUESTION.                
!           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
!                     VERTICES OF POLYGON.                              
!           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
!                     VERTICES OF POLYGON.                              
!           N       - NUMBER OF VERTICES IN THE POLYGON.                
!           INOUT   - THE SIGNAL RETURNED:                              
!                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
!                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
!                      1 IF THE POINT IS INSIDE OF THE POLYGON.         
!                                                                       
!        REMARKS                                                        
!           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
!           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
!           OPTIONALLY BE INCREASED BY 1.                               
!           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
!           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
!           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
!           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
!           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
!           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
!           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
!                                                                       
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
!           NONE                                                        
!                                                                       
!        METHOD                                                         
!           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
!           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
!           POINT IS INSIDE OF THE POLYGON.                             
!                                                                       
!     ..................................................................
                                                                       
      IMPLICIT NONE
      REAL*8 X(10),Y(10),XX(N),YY(N),PX,PY
      INTEGER I,J,N,MAXDIM,INOUT                                  
      LOGICAL MX,MY,NX,NY                                               
!      INTEGER O                                                         
!      OUTPUT UNIT FOR PRINTED MESSAGES                                 
!      DATA O/6/                                                         
      MAXDIM=10                                                        
      IF(N.LE.MAXDIM)GO TO 6                                            
      WRITE(*,7)  N                                                      
7     FORMAT('0WARNING:',I5,' TOO GREAT FOR THIS VERSION OF PNPOLY. RESULTS INVALID')                                                 
      RETURN                                                            

6     DO I=1,N                                                        
      X(I)=XX(I)-PX                                                     
      Y(I)=YY(I)-PY 
      ENDDO                                                    


!       INOUT=-1                                                          
!       DO 2 I=1,N                                                        
!       J=1+MOD(I,N)                                                      
!       MX=X(I).GE.0.0                                                    
!       NX=X(J).GE.0.0                                                    
!       MY=Y(I).GE.0.0                                                    
!       NY=Y(J).GE.0.0                                                    
!       IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2       
!       IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3  
!       INOUT=-INOUT                                                      
!       GO TO 2                                                           
! 3     IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5                       
! 4     INOUT=0                                                           
!       RETURN                                                            
! 5     INOUT=-INOUT                                                      
! 2     CONTINUE                                                          
!       RETURN                                                            


      INOUT=-1                                                          
      DO I=1,N                                                        
      J=1+MOD(I,N)                                                      
      MX=X(I).GE.0.0                                                    
      NX=X(J).GE.0.0                                                    
      MY=Y(I).GE.0.0                                                    
      NY=Y(J).GE.0.0                                                    
      IF(.NOT.((MY.OR.NY).AND.(MX.OR.NX)).OR.(MX.AND.NX)) GO TO 2       
      IF(.NOT.(MY.AND.NY.AND.(MX.OR.NX).AND..NOT.(MX.AND.NX))) GO TO 3  
      INOUT=-INOUT                                                      
      GO TO 2                                                           
3     IF( (Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I)).LT.0. ) GOTO 2               
      IF( (Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I)).EQ.0. ) GOTO 4                    
      IF( (Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I)).GT.0. ) GOTO 5                      
4     INOUT=0                                                           
      RETURN                                                            
5     INOUT=-INOUT
2     CONTINUE
      ENDDO                                                                                                             
      RETURN                                                            

      END SUBROUTINE PNPOLY     

!-----------------------------------------------------------------

      subroutine swap_order(x,y,n)
      implicit none
      integer n
      integer i
      real x(n),y(n),x1(n),y1(n)
      
      do i=1,n
      x1(i)=x(n-i+1)
      y1(i)=y(n-i+1)
      enddo
      
      do i=1,n
      x(i)=x1(i)
      y(i)=y1(i)
      enddo
      
      return
      end subroutine swap_order
      
!-----------------------------------------------------------------
      
      subroutine swap_order_d(x,y,n)
      implicit none
      integer n
      integer i
      real*8 x(n),y(n),x1(n),y1(n)
      
      do i=1,n
      x1(i)=x(n-i+1)
      y1(i)=y(n-i+1)
      enddo
      
      do i=1,n
      x(i)=x1(i)
      y(i)=y1(i)
      enddo
      
      return
      end subroutine swap_order_d
      
!-----------------------------------------------------------------

       subroutine check_sense(n,x,y)
       implicit none
       integer n
       real*8 x(5),y(5)
       real*8 xc,yc
       integer i
       real*8 vx(5),vy(5),cross_product,a
       integer sgn
      
       xc=0.
       yc=0.
       do i=1,n
       xc=xc+x(i)
       yc=yc+y(i)
       enddo
       xc=xc/float(n)
       yc=yc/float(n)
       
       do i=1,5
       vx(i)=x(i)-xc
       vy(i)=y(i)-yc
       enddo
 
       do i=1,4
       a=cross_product(vx(i),vy(i),vx(i+1),vy(i+1))
       if (a.gt.0.) sgn=1
       if (a.eq.0.) sgn=0
       if (a.lt.0.) sgn=-1
       print *,sgn
       enddo
       
       return
       end subroutine check_sense
       
!-----------------------------------------------------------------

       subroutine check_sense_r(n,x,y)
       implicit none
       integer n
       real x(5),y(5)
       real xc,yc
       integer i
       real*8 vx(5),vy(5),cross_product,a
       integer sgn
      
       xc=0.
       yc=0.
       do i=1,n
       xc=xc+x(i)
       yc=yc+y(i)
       enddo
       xc=xc/float(n)
       yc=yc/float(n)
       
       do i=1,5
       vx(i)=dble(x(i)-xc)
       vy(i)=dble(y(i)-yc)
       enddo
 
       do i=1,4
       a=cross_product(vx(i),vy(i),vx(i+1),vy(i+1))
       if (a.gt.0.) sgn=1
       if (a.eq.0.) sgn=0
       if (a.lt.0.) sgn=-1
       print *,sgn
       enddo
       
       return
       end subroutine check_sense_r

!-----------------------------------------------------------------

      subroutine sort(n,arr)
      INTEGER n,M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      REAL a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)stop 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      end subroutine sort

!-------------------------------------------------------------------------------
      
      SUBROUTINE indexn(n,arr,indx)
      INTEGER n,indx(n)
      INTEGER arr(n)
      integer, parameter :: M=7,NSTACK=50
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      INTEGER a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)stop 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END SUBROUTINE indexn

!-------------------------------------------------------------------------------

      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n)
      REAL arr(n)
      integer, parameter :: M=7,NSTACK=50
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)stop 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END SUBROUTINE indexx

!-------------------------------------------------------------------------------

      subroutine skip_hashes(lu)
      implicit none
      integer lu
      character*1 hash
10    read(lu,11) hash
11    format((a1))
      if (hash.ne.'#') then
      backspace(lu)
      goto 12
      endif
      goto 10
12    return
      end subroutine skip_hashes

!---------------------------------------------------------------------------

      subroutine find_marker(lu,marker)
      implicit none
      integer lu
      character*4 marker,strng
10      read(lu,*) strng
      if (strng.eq.marker) goto 20
      goto 10
20    return
      end subroutine find_marker

!-------------------------------------------------------------

      subroutine set_rainbow_color_fine(nclass,iclass)
      implicit none
      integer n1,n2,nclass,iclass
      integer m
      real x1,x2,y1,y2


      if ((nclass.eq.1).and.(iclass.ne.0)) then
      call pgsci(2)
      return
      endif

      
      if ((iclass.eq.0).or.(iclass.gt.nclass)) then
      call pgsci(15)
      return
      endif



      n1=35
      n2=99


      x1=float(n1)
      x2=float(n2)
      y2=float(nclass)
      y1=float(nclass-abs(iclass)+1)

      m=int(n1+(y1-1.)/(y2-1.)*(x2-x1)+0.5)


      call pgsci(m)

      return
      end subroutine set_rainbow_color_fine

!--------------------------------------------
     
      subroutine set_rainbow_lut
      implicit none
      integer i
      real l1,l2,l,r,g,b
      

      l1=450.
      l2=700.
      
      do i=20,99
      l=l1+(l2-l1)/float(99-20)*float(i-20)  
!       call rainbow_rgb(l,r,g,b)
      call rainbow_rgb2(l,r,g,b)
      call pgscr(i,r,g,b)
      enddo
      
      return 
      end subroutine set_rainbow_lut

!--------------------------------------------
      
      subroutine rainbow_rgb(l,r,g,b)
! rainbow 400 nm to 700 nm
      implicit none
      real l,r,g,b,t

          r=0.0
          g=0.0
          b=0.0
   
          if ((l.ge.400.0).and.(l.lt.410.0)) then
          t=(l-400.0)/(410.0-400.0)
          r=(0.33*t)-(0.20*t*t)  
          else if ((l.ge.410.0).and.(l.lt.475.0)) then
          t=(l-410.0)/(475.0-410.0)
          r=0.14-(0.13*t*t)  
          else if ((l.ge.545.0).and.(l.lt.595.0)) then
          t=(l-545.0)/(595.0-545.0)
          r=(1.98*t)-(t*t)
          else if ((l.ge.595.0).and.(l.lt.650.0)) then 
          t=(l-595.0)/(650.0-595.0)
          r=0.98+(0.06*t)-(0.40*t*t)
          else if ((l.ge.650.0).and.(l.lt.700.0)) then
          t=(l-650.0)/(700.0-650.0)
          r=0.65-(0.84*t)+(0.20*t*t)
          endif
      
          if ((l.ge.415.0).and.(l.lt.475.0)) then
          t=(l-415.0)/(475.0-415.0)
          g=(0.80*t*t)
          else if ((l.ge.475.0).and.(l.lt.590.0)) then 
          t=(l-475.0)/(590.0-475.0)
          g=0.8+(0.76*t)-(0.80*t*t)
          else if ((l.ge.585.0).and.(l.lt.639.0)) then
          t=(l-585.0)/(639.0-585.0)
          g=0.84-(0.84*t)   
          endif        
      
          if ((l.ge.400.0).and.(l.lt.475.0)) then
          t=(l-400.0)/(475.0-400.0)
          b=(2.20*t)-(1.50*t*t)
          else if ((l.ge.475.0).and.(l.lt.560.0)) then
          t=(l-475.0)/(560.0-475.0)
          b=0.7-(t)+(0.30*t*t)
          endif
   
      return
      end subroutine rainbow_rgb
     
!--------------------------------------------------------------------
      
      subroutine rainbow_rgb2(l,r,g,b)
! rainbow 400 nm to 700 nm
      implicit none
      real l,r,g,b

          r=0.0
          g=0.0
          b=0.0
   
          if ((l.le.700.0).and.(l.gt.650.0)) then
           r=1.
           g=(700.-l)/50.
           b=0.
          else if ((l.le.650.0).and.(l.gt.600.0)) then
           r=(l-600.)/50.
           g=1.
           b=0.
          else if ((l.le.600.0).and.(l.gt.550.0)) then
           r=0.
           g=1.
           b=(600.-l)/50.
          else if ((l.le.550.0).and.(l.gt.500.0)) then 
           r=0.
           g=(l-500.)/50.
           b=1.
          else if ((l.le.500.0).and.(l.gt.450.0)) then
           r=(500.-l)/50.
           g=0.
           b=1.  
          else if ((l.le.450.0).and.(l.gt.400.0)) then
           r=1.
           g=0.
           b=(l-400.)/50.  
          endif
   
      return
      end subroutine rainbow_rgb2

!--------------------------------------------------------------------

      subroutine readem(filename,x,y,n)
      implicit none
      integer n
      real x(n),y(n)
      character*11 filename
      open(09,file=filename,status='old',access='sequential', form='formatted')

      n=0
10    n=n+1
      read(9,*,end=20) x(n),y(n)
      goto 10
20    close(9)
      n=n-1

      return
      end subroutine readem

!------------------------------------------------------------------------
