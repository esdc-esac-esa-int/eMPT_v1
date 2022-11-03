! v1.7 rewritten to handle filter-dependent FORE transformation
      
     include 'empt_modules.f90'

!-------------------------------------------------------

      program k_make     
      use path_to_ref_files
      use config_parameters
      use target_cat
      use derived_parameters
      use pointings_and_groups, only : n_p,no_p,ra_p,dec_p,pa_ap_p

      implicit none
      real t_start,t_stop
      integer i_p,i_t
      integer, allocatable :: kt(:),it(:),jt(:),it_ns(:)
      real, allocatable :: xt(:),yt(:),rx(:),ry(:)
      real xm,ym,xtv2,ytv3
      real ccpa_v3
      real ax_refra_p,ay_refdec_p
      integer ns,i_s
      
      character*2 rn
      character*70 k_list_file
! - slipmap
      integer slitmap(4,365,171)
      common/slit_map/slitmap
! - spectral parameters
!      integer n_disp
      real lambda1,lambda2
      character*6 gname
      common/disp/lambda1,lambda2,gname
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar


      call cpu_time(t_start)

      write(*,*)
      write(*,*) '------------------------------------------------------------------------------------------'
      write(*,*) '--- k_make -- eMPT Slit Assignment k-list Generation Module Version 1.7 24 June 2022 -----'
      write(*,*) '------------------------------------------------------------------------------------------'
      write(*,*)

      call load_derived_parameters

      call file_entry_k_make()

! create output directory 

      write(rn,'(I2)') ntile
      rn=trim(adjustl(rn))
      if (ntile.lt.10) rn='0'//rn(1:1)

      call system('mkdir trial_'//rn(1:2))
      
      call system('mkdir trial_'//rn(1:2)//'/k_make_output')

      call catalog_read()

      call pointings_read()


       write(*,*)
       write(*,'(a,a2)')    ' Trial Identifier No.:    ',rn
       write(*,'(a,a)')     ' Configuration File:      ',confname
       write(*,'(a,a)')     ' Disperser:               ',disperser
       write(*,'(a,f4.2)')  ' S_0:                     ',sthresh
       write(*,'(a,2f7.3)') ' Acceptance Zone:         ',raccx,raccy
       write(*,'(a,a)') ' Target Catalog:             ',catfile
       write(*,'(a,a)') ' Reference Image:            ',imname
       write(*,'(a,a)') ' Segmentation Map:           ',segname
       write(*,'(a,f10.2,f7.3)') ' Angle_to_target,dva_mag:',angle_to_target,dva_mag
       write(*,'(a,i0)') ' n_dither:                    ',n_dither
       write(*,'(a,i0)') ' Total number of targets in input catalog:        ',ntarg
       write(*,'(a,i0)') ' Number of Priority Classes in input catalog:     ',nclass
       write(*,'(a,i0)') ' Total number of pointings in configuration file: ',n_p
       write(*,*)

! Initialize mode-specific subroutines (needed before shutter_values and slitmap construction)

      call load_msa_dimensions_bd
      call load_fpa_dimensions_bd
      call msa_to_fpa_disp(n_disp,1.,xm,ym,xm,ym)

! set up viable shutter map

      call read_shutter_values
      call read_base_msamap
      call make_slitlet_map
      call initialize_slitmap(n_disp,sthresh,1,1)

! create and open k_list output file 

      k_list_file='trial_'//rn(1:2)//'/k_make_output/k_list_raw.txt'

      open(9,file=k_list_file,status='replace',access='sequential', form='formatted')
      write(9,*) n_p

! allocate target k_list arrays

      allocate (xt(ntarg),yt(ntarg),kt(ntarg),it(ntarg),jt(ntarg),rx(ntarg),ry(ntarg),it_ns(ntarg))


      do i_p=1,n_p ! all pointings

      write(*,'(a,i0,a,i0)') ' Processing k-space list for pointing ',i_p,' of ',n_p

      ccpa_v3=pa_ap_p(i_p)-(180.-phi)
      if (ccpa_v3.lt.0.) ccpa_v3=ccpa_v3+360.
      if (ccpa_v3.gt.360.) ccpa_v3=ccpa_v3-360.

! calculate ra,dec reference point for pointing roll angle

      call roll_v3(ccpa_v3,ax_refv2,ay_refv3,ax_refra_p,ay_refdec_p) !RA,Dec


      kt=0
      it=0
      jt=0

      call tangcoord(ra_p(i_p),dec_p(i_p),ntarg,ra_t,dec_t,xt,yt)

      ns=0

      do i_t=1,ntarg ! all targets
      

! does target have a positive priority class of interest?

      if (pri_t(i_t).gt.0) then

      xt(i_t)=xt(i_t)*dva_mag+ax_refra_p !apply dva magnification and add origin at ref position
      yt(i_t)=yt(i_t)*dva_mag+ay_refdec_p 

! Rotate to V2,V3 system orientation

      call roll_v3(-ccpa_v3,xt(i_t),yt(i_t),xtv2,ytv3)

! May now be projected to MSA.

       call sky_to_msa_bd(xtv2,ytv3,xm,ym)
       
 
! Identify which quadrant and shutter targets fall in (kt(i)=0 ~ misses MSA entirely)

      call msa_to_shutter_bd(xm,ym,kt(i_t),it(i_t),jt(i_t),rx(i_t),ry(i_t))


! target in MSA quadrant?

      if ((kt(i_t).ge.1).and.(kt(i_t).le.4)) then

! In viable slitlet?

       if (slitmap(kt(i_t),it(i_t),jt(i_t)).eq.4) then

! In Acceptance zone?

       if ((abs(rx(i_t)).le.raccx).and.(abs(ry(i_t)).le.raccy)) then


! Target is viable, save target catalog index

       ns=ns+1


       it_ns(ns)=i_t

      endif !In Acceptance Zone
      endif !In Viable Slitlet
      endif !In Quadrant
      endif !Valid priority class 

      enddo ! i_t target loop


      write(9,*) i_p,no_p(i_p)

      write(9,*) ns

      do i_s=1,ns
      write(9,211) it_ns(i_s),id_t(it_ns(i_s)),pri_t(it_ns(i_s)),kt(it_ns(i_s)),it(it_ns(i_s)),jt(it_ns(i_s)),rx(it_ns(i_s)),ry(it_ns(i_s))
211   format(I10,2x,I10,x,I3,3x,I1,x,I3,x,I3,x,f6.3,x,f6.3)
      enddo

      enddo !i_p pointing loop

      close(9)

     write(*,*)


! Append parameters to configuration file


      open(11,file=confname,status='old', access='sequential', form='formatted')

      call find_marker(11,'#KM#')

      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- SECTION AUTOMATICALLY APPENDED BY THE K_MAKE MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a)') '# The generated k-list file is located in:'
      write(11,'(a)') '#'
      write(11,'(a)') '#    ./trial_'//rn(1:2)//'/k_make_output/k_list_raw.txt'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- END OF SECTION APPENDED BY THE K_MAKE MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a)') '#KC# -- Marker do not delete or move this line'

      close(11)

      write(*,*)
      call cpu_time(t_stop)
      write(*,'(a)')      '-------------------------------------------------'
      write(*,'(a,f7.1)') ' CPU time used by k_make module [s]:', t_stop-t_start
      write(*,'(a)')      '-------------------------------------------------'
      write(*,*)


      end program k_make
      
!------------------------------------------------------------
 
 
include './reference_files/front_transforms_updt.f90'

include './reference_files/back_transforms_updt.f90'

include './reference_files/shutter_routines_new.f90'

include './reference_files/misc_routines_new.f90'

!------------------------------------------------------------

      subroutine file_entry_k_make()
      use config_parameters
      use derived_parameters
      implicit none
      integer num_arg,iargc

      num_arg=iargc()       
      
      if (num_arg.eq.1) then
      call getarg(num_arg,confname)
      confname=adjustl(confname)
      goto 101
      end if

100   write(*,*)
      write(*,'(a)',advance="no") ' Enter configuration file name: '
      read(*,'(a)') confname
      confname=adjustl(confname)

101   open(9,file=confname,status='old',access='sequential',form='formatted',err=102)
      goto 103
102      write(*,*)
         write(*,*) 'Error: Configuration file not found'
      close(9)
      goto 100
103   continue

      call skip_hashes(9)
      read(9,*) ntile
      if ((ntile.lt.0.).or.(ntile.gt.99.)) then
      write(*,*) 'Error: Run ID out of range in configuration file (00-99)'
      stop
      endif

      call skip_hashes(9)

      read(9,'(a)') disperser
      disperser=trim(adjustl(disperser))
      n_disp=0
      if (disperser.eq.'PRISM/CLEAR') then
      n_disp=7
      n_filt=4
      endif
      if (disperser.eq.'G140M/F070LP') then
      n_disp=1
      n_filt=0
      endif
      if (disperser.eq.'G140M/F100LP') then
      n_disp=1
      n_filt=1
      endif
      if (disperser.eq.'G235M/F170LP') then
      n_disp=2
      n_filt=2
      endif
      if (disperser.eq.'G395M/F290LP') then
      n_disp=3
      n_filt=3
      endif
      if (disperser.eq.'G140H/F070LP') then
      n_disp=4
      n_filt=0
      endif
      if (disperser.eq.'G140H/F100LP') then
      n_disp=4
      n_filt=1
      endif
      if (disperser.eq.'G235H/F170LP') then
      n_disp=5
      n_filt=2
      endif
      if (disperser.eq.'G395H/F290LP') then
      n_disp=6
      n_filt=3
      endif
      if (n_disp.eq.0) then
      write(*,*)
      write(*,*) 'Error: Incorrect specification of disperser/filter in configuration file'
      write(*,*)
      stop
      endif
 
 
      call skip_hashes(9)
      read(9,*) raccx,raccy
      if ((raccx.gt.0.5).or.(raccy.gt.0.5).or.(raccx.lt.0.).or.(raccy.lt.0)) then
      write(*,*)
      write(*,*) 'Error: Acceptance Zone out of range in configuration file'
      write(*,*)
      stop
      endif

      call skip_hashes(9)
      read(9,*) sthresh
      if (sthresh.lt.0.) then
      write(*,*)
      write(*,*) 'Error: Vertical spectral separation threshold out of range in configuration file'
      write(*,*)
      stop
      endif


      call skip_hashes(9)
      read(9,'(a70)') catfile
      catfile=adjustl(catfile)
      open(10,file=trim(catfile),status='old',access='sequential',form='formatted',err=200)
      close(10)
      goto 201
200   write(*,*) 'Error: Target Input Catalog file not found'
      write(*,*) catfile
      close(10)
      stop
201   continue

      call skip_hashes(9)
      read(9,*) imname
      
      call skip_hashes(9)
      read(9,*) segname
      
      call skip_hashes(9)
      read(9,*) !cra,cdec,cpa_ap

! work out pointing reference point coordinates in V2,V3 system

      call msa_to_sky_bd(xm_ref,ym_ref,ax_refv2,ay_refv3) !v2,V3

      call skip_hashes(9)
      read(9,*) !szone_x,szone_y

      call skip_hashes(9)
      read(9,*) angle_to_target
      call dva_calc(angle_to_target,dva_mag)

      call skip_hashes(9)
      read(9,*) !max_c_score1
 
      call skip_hashes(9)
      read(9,*) !max_c_score2

      call skip_hashes(9)
      read(9,*) n_dither
      if ((n_dither.lt.1).or.(n_dither.gt.3)) then
      write(*,*)
      write(*,*)'Error: Number of dithered pointings specified in configuration file out of range'
      write(*,*)'       Allowed values: 1, 2 or 3'
      write(*,*)
      stop
      endif
 
      close(9)

      return

      end subroutine file_entry_k_make

!------------------------------------------------------------
