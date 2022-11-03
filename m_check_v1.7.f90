!v1.7 incresed ID field to 13 digits
           
     include 'empt_modules.f90'

!-------------------------------------------------------

      program m_check
      use path_to_ref_files
      use config_parameters
      use derived_parameters
      use target_cat, only : ntarg,nclass
      use pointings_and_groups, only : n_p,ra_p,dec_p,pa_ap_p
      use m_list
      implicit none
      real t_start,t_stop
      character*2 rn
      character*5 ipap,jpap
      character*70 outdir_1
      character*70 outdir_2
      character*70 outdir_3
      character*120 filename
      character*70 spare_slit_file
      character*80 s_infile

      integer, parameter :: max_dither=3
      integer nn_p(max_dither),nm(max_dither),ip_m(max_dither)
      integer nsm1,nsm2,nsm3
      integer i
 
      integer im,ip,idum
      integer, allocatable :: id(:,:),id_cat(:,:),pri_m(:,:),kt(:,:),it(:,:),jt(:,:)
      real, allocatable :: rx(:,:),ry(:,:)
      
      integer, parameter :: max_spare_m=250 !Maximum number of non-overlapping unused slitlets on MSA
      integer ksm1(max_spare_m),ism1(max_spare_m),jsm1(max_spare_m)
      integer ksm2(max_spare_m),ism2(max_spare_m),jsm2(max_spare_m)
      integer ksm3(max_spare_m),ism3(max_spare_m),jsm3(max_spare_m)


      integer ntt
      integer, allocatable :: ktt(:),itt(:),jtt(:),idtt(:)

      call cpu_time(t_start)


      write(*,*)
      write(*,*) '-------------------------------------------------------------------------------'
      write(*,*) '--- m_check - final m-list Spoiler Check Module Version 1.7  24 October 2022 --'
      write(*,*) '-------------------------------------------------------------------------------'
      write(*,*)

      call load_derived_parameters

      call file_entry_m_check()

      write(rn,'(i2)') ntile
      rn=trim(adjustl(rn))
      if (ntile.lt.10) rn='0'//rn(1:1)

      call system('mkdir trial_'//rn(1:2)//'/m_check_output')

      call pointings_read()

      call catalog_read()


      if (n_dither.eq.1) m_list_file='trial_'//rn(1:2)//'/m_make_output/single_m_list.txt'
      if (n_dither.eq.2) m_list_file='trial_'//rn(1:2)//'/m_make_output/pair_m_list.txt'
      if (n_dither.eq.3) m_list_file='trial_'//rn(1:2)//'/m_make_output/triple_m_list.txt'
   
      call m_list_size_read()

       write(*,*)
       write(*,'(a,a2)') ' Trial Identifier No.:        ',rn
       write(*,'(a,a)')  ' Configuration File:          ',confname
       write(*,'(a,i0)') ' n_dither:                    ',n_dither
       write(*,'(a,i0)') ' Number of targets in input catalog:               ',ntarg
       write(*,'(a,i0)') ' Number of Priority Classes in input catalog:      ',nclass
       write(*,'(a,i0)') ' Total number of pointings in configuration file:  ',n_p
       write(*,'(a,a)')  ' m_list file:                                      ',m_list_file
       write(*,'(a,i0)') ' Maximum number of targets in m_list per pointing: ',n_max_m
       write(*,'(a,i0)') ' Target optimal single/pair/triple pointing:       ',mtarget
       write(*,*)


! open and read in m_list file

      allocate (id(max_dither,n_max_m),id_cat(max_dither,n_max_m),pri_m(max_dither,n_max_m))
      allocate (kt(max_dither,n_max_m),it(max_dither,n_max_m),jt(max_dither,n_max_m))
      allocate (rx(max_dither,n_max_m),ry(max_dither,n_max_m))
 
 !---------------------------------------------------------------

      if (n_dither.eq.1) then


      open(11,file=m_list_file, status='old',access='sequential', form='formatted')

      read(11,*) !nsingles

      do i=1,mtarget

      read(11,*)
      read(11,*) ip_m(1)
      read(11,*) nn_p(1)

      ip=1

      read(11,*) nm(ip)
      do im=1,nm(ip)
      read(11,*) id(ip,im),id_cat(ip,im),pri_m(ip,im),kt(ip,im),it(ip,im),jt(ip,im),rx(ip,im),ry(ip,im)
      enddo !im loop

      enddo !i

      close(11)

      write(*,*)
      write(*,*) 'm-list file located and read in'

      write(ipap,"(i0)") nn_p(1)
      spare_slit_file='./trial_'//rn(1:2)//'/m_pick_output/pointing_'//trim(adjustl(ipap))//'/sky_shutters.txt'

      nsm1=0

      if (.not.skypad) goto 7

      open(11,file=spare_slit_file, status='old',access='sequential', form='formatted', err=7)

      
      read(11,*) !skip header
5     nsm1=nsm1+1
      read(11,*,end=6) idum,ksm1(nsm1),ism1(nsm1),jsm1(nsm1)
      goto 5
6     nsm1=nsm1-1
      close(11)

      write(*,*)
      write(*,*) 'Sky shutter list located and read in'


7     continue
      
      
      
      endif !ndither=1

!---------------------------------------------------------

      if (n_dither.eq.2) then

      open(11,file=m_list_file, status='old',access='sequential', form='formatted')

      read(11,*) !npairs

      do i=1,mtarget

      read(11,*)
      read(11,*) ip_m(1),ip_m(2)
      read(11,*) nn_p(1),nn_p(2)

      ip=1

      read(11,*) nm(ip)
      do im=1,nm(ip)
      read(11,*) id(ip,im),id_cat(ip,im),pri_m(ip,im),kt(ip,im),it(ip,im),jt(ip,im),rx(ip,im),ry(ip,im)
      enddo !im loop

      ip=2

      read(11,*) nm(ip)
      do im=1,nm(ip)
      read(11,*) id(ip,im),id_cat(ip,im),pri_m(ip,im),kt(ip,im),it(ip,im),jt(ip,im),rx(ip,im),ry(ip,im)
      enddo !im loop

      enddo !i

      close(11)

      write(*,*)
      write(*,*) 'm-list file located and read in'

      write(ipap,"(i0)") nn_p(1)
      write(jpap,"(i0)") mtarget
      spare_slit_file='./trial_'//rn(1:2)//'/m_pick_output/pair_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/sky_shutters.txt'

      nsm1=0

      if (.not.skypad) goto 17


      open(11,file=spare_slit_file, status='old',access='sequential', form='formatted', err=17)
      
      read(11,*) !skip header
15     nsm1=nsm1+1
      read(11,*,end=16) idum,ksm1(nsm1),ism1(nsm1),jsm1(nsm1)
      goto 15
16     nsm1=nsm1-1
      close(11)

17     continue

      write(ipap,"(i0)") nn_p(2)
      write(jpap,"(i0)") mtarget
      spare_slit_file='./trial_'//rn(1:2)//'/m_pick_output/pair_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/sky_shutters.txt'

      nsm2=0

      if (.not.skypad) goto 27


      open(11,file=spare_slit_file, status='old',access='sequential', form='formatted', err=27)
      
      read(11,*) !skip header
25     nsm2=nsm2+1
      read(11,*,end=26) idum,ksm2(nsm2),ism2(nsm2),jsm2(nsm2)
      goto 25
26     nsm2=nsm2-1
      close(11)


      write(*,*)
      write(*,*) 'Sky shutter lists located and read in'


27     continue



      endif

!---------------------------------------------------------

      if (n_dither.eq.3) then

 
      open(11,file=m_list_file, status='old',access='sequential', form='formatted')

      read(11,*) !ntriples

      do i=1,mtarget

      read(11,*)
      read(11,*) ip_m(1),ip_m(2),ip_m(3)
      read(11,*) nn_p(1),nn_p(2),nn_p(3)

      ip=1

      read(11,*) nm(ip)
      do im=1,nm(ip)
      read(11,*) id(ip,im),id_cat(ip,im),pri_m(ip,im),kt(ip,im),it(ip,im),jt(ip,im),rx(ip,im),ry(ip,im)
      enddo !im loop

      ip=2

      read(11,*) nm(ip)
      do im=1,nm(ip)
      read(11,*) id(ip,im),id_cat(ip,im),pri_m(ip,im),kt(ip,im),it(ip,im),jt(ip,im),rx(ip,im),ry(ip,im)
      enddo !im loop

       ip=3

      read(11,*) nm(ip)
      do im=1,nm(ip)
      read(11,*) id(ip,im),id_cat(ip,im),pri_m(ip,im),kt(ip,im),it(ip,im),jt(ip,im),rx(ip,im),ry(ip,im)
      enddo !im loop

      enddo !i

      close(11)

      write(*,*)
      write(*,*) 'm-list file located and read in'

      write(ipap,"(i0)") nn_p(1)
      write(jpap,"(i0)") mtarget
      spare_slit_file='./trial_'//rn(1:2)//'/m_pick_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/sky_shutters.txt'

      nsm1=0

      if (.not.skypad) goto 37


      open(11,file=spare_slit_file, status='old',access='sequential', form='formatted', err=37)
      
      read(11,*) !skip header
35     nsm1=nsm1+1
      read(11,*,end=36) idum,ksm1(nsm1),ism1(nsm1),jsm1(nsm1)
      goto 35
36     nsm1=nsm1-1
      close(11)

37     continue

      write(ipap,"(i0)") nn_p(2)
      write(jpap,"(i0)") mtarget
      spare_slit_file='./trial_'//rn(1:2)//'/m_pick_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/sky_shutters.txt'

      nsm2=0

      if (.not.skypad) goto 47

      open(11,file=spare_slit_file, status='old',access='sequential', form='formatted', err=47)
      
      read(11,*) !skip header
45     nsm2=nsm2+1
      read(11,*,end=46) idum,ksm2(nsm2),ism2(nsm2),jsm2(nsm2)
      goto 45
46     nsm2=nsm2-1
      close(11)


47     continue

      write(ipap,"(i0)") nn_p(3)
      write(jpap,"(i0)") mtarget
      spare_slit_file='./trial_'//rn(1:2)//'/m_pick_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/sky_shutters.txt'

      nsm3=0

      if (.not.skypad) goto 57

      open(11,file=spare_slit_file, status='old',access='sequential', form='formatted',err=57)
      
      read(11,*) !skip header
55     nsm3=nsm3+1
      read(11,*,end=56) idum,ksm3(nsm3),ism3(nsm3),jsm3(nsm3)
      goto 55
56     nsm3=nsm3-1
      close(11)


      write(*,*)
      write(*,*) 'Sky shutter lists located and read in'


57     continue

      endif


!---------------------------------------------------------

! create output directories

      if (n_dither.eq.1) then
      write(ipap,"(i0)") nn_p(1)
      outdir_1='trial_'//rn(1:2)//'/m_check_output/pointing_'//trim(adjustl(ipap))
      call system( 'mkdir '//outdir_1)
      endif

      if (n_dither.eq.2) then
      write(jpap,"(i0)") mtarget
      call system( 'mkdir trial_'//rn(1:2)//'/m_check_output/pair_'//trim(adjustl(jpap)))
      write(ipap,"(i0)") nn_p(1)
      outdir_1='trial_'//rn(1:2)//'/m_check_output/pair_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))
      call system( 'mkdir '//outdir_1)
      write(ipap,"(i0)") nn_p(2)
      outdir_2='trial_'//rn(1:2)//'/m_check_output/pair_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))
      call system( 'mkdir '//outdir_2)
      endif

      if (n_dither.eq.3) then
      write(jpap,"(i0)") mtarget
      call system( 'mkdir trial_'//rn(1:2)//'/m_check_output/triple_'//trim(adjustl(jpap)))
      write(ipap,"(i0)") nn_p(1)
      outdir_1='trial_'//rn(1:2)//'/m_check_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))
      call system( 'mkdir '//outdir_1)
      write(ipap,"(i0)") nn_p(2)
      outdir_2='trial_'//rn(1:2)//'/m_check_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))
      call system( 'mkdir '//outdir_2)
      write(ipap,"(i0)") nn_p(3)
      outdir_3='trial_'//rn(1:2)//'/m_check_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))
      call system( 'mkdir '//outdir_3)
      endif


!      make panel plots
       

       allocate (ktt(n_max_m),itt(n_max_m),jtt(n_max_m),idtt(n_max_m))


      do ip=n_dither,1,-1


      ntt=nm(ip)
      
      do i=1,nm(ip)
      ktt(i)=kt(ip,i)
      itt(i)=it(ip,i)
      jtt(i)=jt(ip,i)
      idtt(i)=id_cat(ip,i)
      enddo
     
       if (ip.eq.1) filename=trim(adjustl(outdir_1))//'/slitlet_panel_plot.ps/cps'
       if (ip.eq.2) filename=trim(adjustl(outdir_2))//'/slitlet_panel_plot.ps/cps'
       if (ip.eq.3) filename=trim(adjustl(outdir_3))//'/slitlet_panel_plot.ps/cps'


       if (ip.eq.1) call draw_panels_on_sky(ra_p(ip_m(ip)),dec_p(ip_m(ip)),pa_ap_p(ip_m(ip)),ntt,idtt,ktt,itt,jtt,nsm1,ksm1,ism1,jsm1,filename)       
       if (ip.eq.2) call draw_panels_on_sky(ra_p(ip_m(ip)),dec_p(ip_m(ip)),pa_ap_p(ip_m(ip)),ntt,idtt,ktt,itt,jtt,nsm2,ksm2,ism2,jsm2,filename)
       if (ip.eq.3) call draw_panels_on_sky(ra_p(ip_m(ip)),dec_p(ip_m(ip)),pa_ap_p(ip_m(ip)),ntt,idtt,ktt,itt,jtt,nsm3,ksm3,ism3,jsm3,filename)

       enddo

! append to configuration file 

      open(11,file=confname,status='old',access='sequential', form='formatted')

      call find_marker(11,'#MC#')

      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- SECTION AUTOMATICALLY APPENDED BY THE M_CHECK MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# Refer to the close-up slitlet proximity maps in:'
      write(11,'(a)') '#'
        do ip=1,n_dither
         if (ip.eq.1) filename='./'//trim(adjustl(outdir_1))//'/slitlet_panel_plot.ps'
         if (ip.eq.2) filename='./'//trim(adjustl(outdir_2))//'/slitlet_panel_plot.ps'
         if (ip.eq.3) filename='./'//trim(adjustl(outdir_3))//'/slitlet_panel_plot.ps'
       write(11,'(a)') '#   '//filename 
        enddo
       write(11,'(a)') '#'
      write(11,'(a)') '# to gauge the contamination levels of each target at each final pointing.'
      write(11,'(a)') '#'
      write(11,'(a)') '# Targets deemed to be excessively contaminated may be excluded from observation by' 
      write(11,'(a)') '# changing the sign of their Priority Class assignments to be negative in the'
      write(11,'(a)') '# input catalog, and running the ipa and subsequent modules anew.'
      write(11,'(a)') '#'
      write(11,'(a)') '#'  
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- END OF SECTION APPENDED BY THE M_CHECK MODULE ----------'
      write(11,'(a)') '#'


 
       write(*,*)
       call cpu_time(t_stop)
       write(*,'(a)')      '-------------------------------------------------'
       write(*,'(a,f7.1)') ' CPU time used by m_check module [s]:',t_stop-t_start
       write(*,'(a)')      '-------------------------------------------------'
       write(*,*)

      end program m_check

!------------------------------------------------------------

include './reference_files/front_transforms_updt.f90'

include './reference_files/back_transforms_updt.f90'

include './reference_files/shutter_routines_new.f90'

include './reference_files/misc_routines_new.f90'

!------------------------------------------------------------

      subroutine file_entry_m_check()
      use config_parameters
      use derived_parameters
      use target_cat, only : nclass
      implicit none
      integer num_arg,iargc
      character*70 skystring

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
      read(9,*) !disperser

      call skip_hashes(9)
      read(9,*) !raccx,raccy

      call skip_hashes(9)
      read(9,*) !sthresh


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
      read(9,*) !ref image name
      
      call skip_hashes(9)
      read(9,*) !seg map name
      
      call skip_hashes(9)
      read(9,*) !cra,cdec,cpa_ap

 
 ! work out pointing reference point coordinates in V2,V3 on-sky systems

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
  
      call find_marker(9,'#KC#')

      call skip_hashes(9)
      read(9,*) nclass

 
      call skip_hashes(9)
      if (n_dither.eq.1) then
      read(9,*) !(order_matrix(i,1), i=1,nclass)
      endif !n_dither=1


      if (n_dither.eq.2) then
      read(9,*) !(order_matrix(i,1), i=1,nclass)
      call skip_hashes(9)
      read(9,*) !(order_matrix(i,2), i=1,nclass)
      endif !n_dither=2


      if (n_dither.eq.3) then
      read(9,*) !(order_matrix(i,1), i=1,nclass)
      call skip_hashes(9)
      read(9,*) !(order_matrix(i,2), i=1,nclass)
      call skip_hashes(9)
      read(9,*) !(order_matrix(i,3), i=1,nclass)
     endif !n_dither=3


      if (n_dither.gt.1) then
      call skip_hashes(9)
      read(9,*) !min_dx,max_dx,min_dy,max_dy
      endif


      call skip_hashes(9)
      read(9,*) !(weights(i),i=1,nclass)

      call skip_hashes(9)
      read(9,*) mtarget !optimal single, pair or triple pointing number

      call skip_hashes(9)
      skypad=.false.
      read(9,*) skystring
      skystring=trim(adjustl(skystring))
      if ((skystring.eq.'Y').or.(skystring.eq.'y')) skypad=.true.


      close(9)

      return

      end subroutine file_entry_m_check

!-----------------------------------------------------------------

      subroutine draw_panels_on_sky(ccra,ccdec,ccpa_ap,ntt,idtt,ktt,itt,jtt,ns,ks,is,js,sky_plot)
      use config_parameters, only : dva_mag
      use derived_parameters, only : phi,ax_refv2,ay_refv3
      use target_cat, only : ntarg,nclass,id_t,ra_t,dec_t,pri_t
      implicit none
      real ccpa_ap,ccpa_v3
      real*8 ccra,ccdec
      real ax_refra_p,ay_refdec_p
      integer ns
      integer, parameter :: max_spare_m=250 !Maximum number of non-overlapping unused slitlets on MSA
      integer ks(max_spare_m),is(max_spare_m),js(max_spare_m)
      integer ntt
      integer ktt(ntt),itt(ntt),jtt(ntt),idtt(ntt)
      real ax,ay,axr,ayr
      integer i,ic
      real xx,yy,xp1,xp2,yp1,yp2
      real x1(5),y1(5)
      real av_xpitch_m,av_ypitch_m
      integer kk,ii,jj,it,iq,jq
      real fx1,fx2,fy1,fy2
      character*13 tnum
      integer ntnum
      character*120 sky_plot
      character*2 rn
      common/run/rn
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar

       
      call load_msa_dimensions_bd

      av_xpitch_m=(xpitch(1)+xpitch(2)+xpitch(3)+xpitch(4))/4.
      av_ypitch_m=(ypitch(1)+ypitch(2)+ypitch(3)+ypitch(4))/4.     
      
      ccpa_v3=ccpa_ap-(180.-phi)
      if (ccpa_v3.lt.0.) ccpa_v3=ccpa_v3+360.
      if (ccpa_v3.gt.360.) ccpa_v3=ccpa_v3-360.

! work out pointing reference position for present roll angle

      call roll_v3(ccpa_v3,ax_refv2,ay_refv3,ax_refra_p,ay_refdec_p) !RA,Dec


! Plot target catalog and MSA FOV and on sky

        CALL PGBEGIN(0,sky_plot,1,1)
        CALL PGSLS(1)
        CALL PGSLW(1)
        CALL PGSCF(1)
        CALL PGSCH(1.)

      call set_rainbow_lut

      it=0

      do jq=10,1,-1

      do iq=1,25
      
      it=it+1
    
      if (it.gt.ntt+ns) goto 1000
      
      fx1=float(iq-1)*1/25.
      fx2=float(iq)*1/25.
      fy1=float(jq-1)*1/10.
      fy2=fy1+1/11.

      if (it.le.ntt) then

      kk=ktt(it)
      ii=itt(it)
      jj=jtt(it)
 
      x1(1)=(float(ii)-0.5)*xpitch(kk)
      y1(1)=(float(jj)-0.5)*ypitch(kk)
      call rottran(kk,1,x1,y1) !handles intra-quadrant rotations

       call msa_to_sky_rot_bd(x1(1),y1(1),xx,yy) 

       x1(1)=xx
       y1(1)=yy

       endif
       
      if (it.gt.ntt) then

      kk=ks(it-ntt)
      ii=is(it-ntt)
      jj=js(it-ntt)
 
      x1(1)=(float(ii)-0.5)*xpitch(kk)
      y1(1)=(float(jj)-0.5)*ypitch(kk)
      call rottran(kk,1,x1,y1) !handles intra-quadrant rotations

       call msa_to_sky_rot_bd(x1(1),y1(1),xx,yy) 

       x1(1)=xx
       y1(1)=yy

       endif
             

        xp1=x1(1)-3.0*av_xpitch_m*2.55 !nominal plate scale in arcsec/mm
        xp2=x1(1)+3.0*av_xpitch_m*2.55
        yp1=y1(1)-3.5*av_ypitch_m*2.59
        yp2=y1(1)+3.5*av_ypitch_m*2.59


        CALL PGVPORT(fx1,fx2,fy1,fy2)
       CALL PGWNAD(xp2,xp1,yp1,yp2) !flipped x-direction
!         CALL PGWNAD(xp1,xp2,yp2,yp1) !rotated to STScI APT convention 
        call pgsci(1)
        CALL PGSLW(1)
        call PGBOX('BC',0.,0,'BC',0.,0)

! draw target slitlets    

      if (it.le.ntt) then

      call pgsci(15)
      call pgslw(1)

      call draw_shutter(ktt(it),itt(it),jtt(it))
      call draw_shutter(ktt(it),itt(it),jtt(it)+1)
      call draw_shutter(ktt(it),itt(it),jtt(it)-1)
      
      call pgnumb(idtt(it),0,1,tnum,ntnum)
      xx=(xp1+xp2)/2.
      yy=yp1+(yp2-yp1)*0.05
!        yy=yp2+(yp1-yp2)*0.05
      call pgsch(0.3)
      call pgsci(1)
      call pgptext(xx,yy,0.,0.5,tnum(1:ntnum))
      
      endif

! draw sky slitlets    

      if (it.gt.ntt) then

      call pgsci(15)
      call pgslw(1)
     
      call draw_shutter(ks(it-ntt),is(it-ntt),js(it-ntt))
      call draw_shutter(ks(it-ntt),is(it-ntt),js(it-ntt)+1)
      call draw_shutter(ks(it-ntt),is(it-ntt),js(it-ntt)-1)

      call pgnumb(it-ntt,0,1,tnum,ntnum)
      tnum='B'//tnum(1:ntnum)
      ntnum=ntnum+1
      xx=(xp1+xp2)/2.
      yy=yp1+(yp2-yp1)*0.05
!       yy=yp2+(yp1-yp2)*0.05
      call pgsch(0.3)
      call pgsci(1)
      call pgptext(xx,yy,0.,0.5,tnum(1:ntnum))

      endif

! overlay targets

! priority class 0 and below targets first

      call pgslw(1)

      ic=0
      do i=1,ntarg
      if (pri_t(i).le.ic) then
      call tangcoord_single(ccra,ccdec,ra_t(i),dec_t(i),ax,ay)
      ax=ax*dva_mag+ax_refra_p
      ay=ay*dva_mag+ay_refdec_p
      call roll_v3(-ccpa_v3+phi,ax,ay,axr,ayr)
      ax=axr
      ay=ayr
      if ((ax.ge.xp1).and.(ax.le.xp2).and.(ay.ge.yp1).and.(ay.le.yp2)) then
!      call set_rainbow_color_fine(nclass,pri_t(i))
      call pgsci(15)
      call pgsch(0.5)
      call pgpoint(1,ax,ay,17)
      call pgsci(1)
      call pgnumb(id_t(i),0,1,tnum,ntnum)            ! label sub_cat id number
      call pgsch(.15)
      call pgptext(ax,ay+0.1,0.,0.5,tnum(1:ntnum))
      endif
      endif
      enddo

! all other priority classes

      call pgslw(1)

      do ic=nclass,1,-1
      do i=1,ntarg
      if (pri_t(i).eq.ic) then
      call tangcoord_single(ccra,ccdec,ra_t(i),dec_t(i),ax,ay)
      ax=ax*dva_mag+ax_refra_p
      ay=ay*dva_mag+ay_refdec_p
      call roll_v3(-ccpa_v3+phi,ax,ay,axr,ayr)
      ax=axr
      ay=ayr
      if ((ax.ge.xp1).and.(ax.le.xp2).and.(ay.ge.yp1).and.(ay.le.yp2)) then
      call set_rainbow_color_fine(nclass,pri_t(i))
      call pgsch(0.5)
      call pgpoint(1,ax,ay,17)
      call pgsci(1)
      call pgnumb(id_t(i),0,1,tnum,ntnum)            ! label sub_cat id number
      call pgsch(.15)
      call pgptext(ax,ay+0.1,0.,0.5,tnum(1:ntnum))
      endif
      endif
      enddo
      enddo


      enddo !iq
      enddo !jq
      
1000 continue

      call pgend

      return

      end subroutine draw_panels_on_sky

!---------------------------------------------------------------------------    

      subroutine draw_shutter(kk,ii,jj)
      implicit none
      integer kk,ii,jj,i
      real x1(5),y1(5),x1r(5),y1r(5),xx,yy
      real av_xpitch_m,av_ypitch_m
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar

      av_xpitch_m=(xpitch(1)+xpitch(2)+xpitch(3)+xpitch(4))/4.
      av_ypitch_m=(ypitch(1)+ypitch(2)+ypitch(3)+ypitch(4))/4.


      x1(1)=float(ii-1)*xpitch(kk)
      y1(1)=float(jj-1)*ypitch(kk)


      x1(2)=float(ii-1)*xpitch(kk)
      y1(2)=float(jj)*ypitch(kk)

      x1(3)=float(ii)*xpitch(kk)
      y1(3)=float(jj)*ypitch(kk)

      x1(4)=float(ii)*xpitch(kk)
      y1(4)=float(jj-1)*ypitch(kk)


      x1(1)=x1(1)+(xpitch(kk)-xopen(kk))/2.
      y1(1)=y1(1)+(ypitch(kk)-yopen(kk))/2.

      x1(2)=x1(2)+(xpitch(kk)-xopen(kk))/2.
      y1(2)=y1(2)-(ypitch(kk)-yopen(kk))/2.

      x1(3)=x1(3)-(xpitch(kk)-xopen(kk))/2.
      y1(3)=y1(3)-(ypitch(kk)-yopen(kk))/2.

      x1(4)=x1(4)-(xpitch(kk)-xopen(kk))/2.
      y1(4)=y1(4)+(ypitch(kk)-yopen(kk))/2.


      call rottran(kk,4,x1,y1) !handles intra-quadrant rotations
      
      do i=1,4
!        call msa_to_sky_bd(x1(i),y1(i),xx,yy) !v2,v3
       call msa_to_sky_rot_bd(x1(i),y1(i),xx,yy) !stays in msa orientation
       x1(i)=xx
       y1(i)=yy
      enddo

       x1(5)=x1(1)
       y1(5)=y1(1)

        do i=1,5
!         call roll_v3(pa_v3,x1(i),y1(i),x1r(i),y1r(i)) !ra,dec
        x1r(i)=x1(i)
        y1r(i)=y1(i)
        enddo

!        call pgpoly(5,x1r,y1r)
      call pgslw(1)
      call pgline(5,x1r,y1r) 

             
      x1(1)=float(ii-1)*xpitch(kk)
      y1(1)=float(jj-1)*ypitch(kk)+av_ypitch_m


      x1(2)=float(ii-1)*xpitch(kk)
      y1(2)=float(jj)*ypitch(kk)+av_ypitch_m

      x1(3)=float(ii)*xpitch(kk)
      y1(3)=float(jj)*ypitch(kk)+av_ypitch_m

      x1(4)=float(ii)*xpitch(kk)
      y1(4)=float(jj-1)*ypitch(kk)+av_ypitch_m

      x1(1)=x1(1)+(xpitch(kk)-xopen(kk))/2.
      y1(1)=y1(1)+(ypitch(kk)-yopen(kk))/2.

      x1(2)=x1(2)+(xpitch(kk)-xopen(kk))/2.
      y1(2)=y1(2)-(ypitch(kk)-yopen(kk))/2.

      x1(3)=x1(3)-(xpitch(kk)-xopen(kk))/2.
      y1(3)=y1(3)-(ypitch(kk)-yopen(kk))/2.

      x1(4)=x1(4)-(xpitch(kk)-xopen(kk))/2.
      y1(4)=y1(4)+(ypitch(kk)-yopen(kk))/2.


      call rottran(kk,4,x1,y1) !handles intra-quadrant rotations
      
      do i=1,4
!        call msa_to_sky_bd(x1(i),y1(i),xx,yy) !v2,v3
       call msa_to_sky_rot_bd(x1(i),y1(i),xx,yy) !v2,v3
       x1(i)=xx
       y1(i)=yy
      enddo

       x1(5)=x1(1)
       y1(5)=y1(1)

        do i=1,5
!         call roll_v3(pa_v3,x1(i),y1(i),x1r(i),y1r(i)) !ra,dec
        x1r(i)=x1(i)
        y1r(i)=y1(i)
        enddo


!        call pgpoly(5,x1r,y1r)
      call pgslw(1)
      call pgline(5,x1r,y1r) 

             
      x1(1)=float(ii-1)*xpitch(kk)
      y1(1)=float(jj-1)*ypitch(kk)-av_ypitch_m


      x1(2)=float(ii-1)*xpitch(kk)
      y1(2)=float(jj)*ypitch(kk)-av_ypitch_m

      x1(3)=float(ii)*xpitch(kk)
      y1(3)=float(jj)*ypitch(kk)-av_ypitch_m

      x1(4)=float(ii)*xpitch(kk)
      y1(4)=float(jj-1)*ypitch(kk)-av_ypitch_m

      x1(1)=x1(1)+(xpitch(kk)-xopen(kk))/2.
      y1(1)=y1(1)+(ypitch(kk)-yopen(kk))/2.

      x1(2)=x1(2)+(xpitch(kk)-xopen(kk))/2.
      y1(2)=y1(2)-(ypitch(kk)-yopen(kk))/2.

      x1(3)=x1(3)-(xpitch(kk)-xopen(kk))/2.
      y1(3)=y1(3)-(ypitch(kk)-yopen(kk))/2.

      x1(4)=x1(4)-(xpitch(kk)-xopen(kk))/2.
      y1(4)=y1(4)+(ypitch(kk)-yopen(kk))/2.

      call rottran(kk,4,x1,y1) !handles intra-quadrant rotations
      
      do i=1,4
!        call msa_to_sky_bd(x1(i),y1(i),xx,yy) !v2,v3
       call msa_to_sky_rot_bd(x1(i),y1(i),xx,yy) !v2,v3
       x1(i)=xx
       y1(i)=yy
      enddo

       x1(5)=x1(1)
       y1(5)=y1(1)

        do i=1,5
!         call roll_v3(pa_v3,x1(i),y1(i),x1r(i),y1r(i)) !ra,dec
        x1r(i)=x1(i)
        y1r(i)=y1(i)
        enddo
             

!        call pgpoly(5,x1r,y1r)
      call pgslw(1)
      call pgline(5,x1r,y1r) 
      
      
      return
      
      end subroutine draw_shutter
      
!------------------------------------------------------------------------------------------
