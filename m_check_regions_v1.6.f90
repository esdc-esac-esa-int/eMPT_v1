! v1.7 id string length increased to 13

      include 'empt_modules.f90'

!-------------------------------------------------------------------------

      program m_check_regions
      use path_to_ref_files
      use config_parameters, only : ntile,confname,n_dither,mtarget,skypad
      use pointings_and_groups, only : n_p,ra_p,dec_p,pa_ap_p
      use m_list
      use target_cat, only : ntarg,nclass
      use derived_parameters
      implicit none
      integer ip,idum,i
      character*5 ipap,jpap
      character*120 filename,filename2
      character*70 outdir_1
      character*70 outdir_2
      character*70 outdir_3
      character*80 s_infile
      character*120 sc_infile
      logical shutter_scores

      character*2 rn
      real t_start,t_stop


      integer, parameter :: max_dither=3
      
      integer, allocatable :: id(:,:),id_cat(:,:),pri_m(:,:)
      integer, allocatable :: kt(:,:),it(:,:),jt(:,:)
      
      integer, allocatable :: sh_sc1(:,:),sh_sc2(:,:),sh_sc3(:,:)

      real, allocatable :: rx(:,:),ry(:,:)

      integer :: im
      integer, allocatable :: nn_p(:),nm(:),ip_m(:)

      integer :: ntt
      integer, allocatable :: ktt(:),itt(:),jtt(:),idtt(:)

      integer :: nsm1,nsm2,nsm3      
      integer, parameter :: max_spare_m=250 !Maximum number of non-overlapping unused slitlets on MSA
      integer :: ksm1(max_spare_m),ism1(max_spare_m),jsm1(max_spare_m)
      integer :: ksm2(max_spare_m),ism2(max_spare_m),jsm2(max_spare_m)
      integer :: ksm3(max_spare_m),ism3(max_spare_m),jsm3(max_spare_m)
      integer spsh_sc1(5,max_spare_m),spsh_sc2(5,max_spare_m),spsh_sc3(5,max_spare_m)
      

! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar

      call cpu_time(t_start)


      write(*,*)
      write(*,*) '--------------------------------------------------------------------------------------------'
      write(*,*) '--- m_check_regions - additional m-list Spoiler Check Module Version 1.7 4 August 2022 --=--'
      write(*,*) '--------------------------------------------------------------------------------------------'
      write(*,*)

     call load_derived_parameters

!    read in configuration file with pointing listing

      call file_entry_m_check_regions()


! create output directory

      write(rn,'(i2)') ntile
      rn=trim(adjustl(rn))
      if (ntile.lt.10) rn='0'//rn(1:1)

      call system('mkdir trial_'//rn(1:2)//'/m_check_output')

      call catalog_read()

      call pointings_read()


      if (n_dither.eq.1) m_list_file='trial_'//rn(1:2)//'/m_make_output/single_m_list.txt'
      if (n_dither.eq.2) m_list_file='trial_'//rn(1:2)//'/m_make_output/pair_m_list.txt'
      if (n_dither.eq.3) m_list_file='trial_'//rn(1:2)//'/m_make_output/triple_m_list.txt'
   
      call m_list_size_read()

      allocate (id(max_dither,n_max_m),id_cat(max_dither,n_max_m),pri_m(max_dither,n_max_m))
      allocate (kt(max_dither,n_max_m),it(max_dither,n_max_m),jt(max_dither,n_max_m))
      allocate (sh_sc1(5,n_max_m),sh_sc2(5,n_max_m),sh_sc3(5,n_max_m))
      allocate (rx(max_dither,n_max_m),ry(max_dither,n_max_m))
      allocate (nn_p(n_p),nm(max_dither),ip_m(max_dither))
      allocate (ktt(n_max_m),itt(n_max_m),jtt(n_max_m),idtt(n_max_m))


! initialize shutter scores so all plotted green in case shutter score file not present

      sh_sc1=0
      sh_sc2=0
      sh_sc3=0
      
      spsh_sc1=0
      spsh_sc2=0
      spsh_sc3=0

      sh_sc1(1,:)=1 !central shutters
      sh_sc2(1,:)=1
      sh_sc3(1,:)=1

      
      shutter_scores=.true.

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
      s_infile='./trial_'//rn(1:2)//'/m_pick_output/pointing_'//trim(adjustl(ipap))//'/sky_shutters.txt'

      nsm1=0

      if (.not.skypad) goto 7

      open(11,file=s_infile, status='old',access='sequential', form='formatted', err=7)

      
      read(11,*) !skip header
5     nsm1=nsm1+1
      read(11,*,end=6) idum,ksm1(nsm1),ism1(nsm1),jsm1(nsm1)
      goto 5
6     nsm1=nsm1-1

      write(*,*)
      write(*,*) 'Sky shutter list located and read in'

7     continue

      close(11)

      sc_infile='./trial_'//rn(1:2)//'/m_check_output/pointing_'//trim(adjustl(ipap))//'/slitlet_contamination_scores.txt'
      
      ip=1

      open(11,file=sc_infile, status='old',access='sequential', form='formatted',err=1010)
      call skip_hashes(11)
      do im=1,nm(ip)
      read(11,*) idum,idum,idum,idum,idum,idum,idum,sh_sc1(1,im),sh_sc1(3,im),sh_sc1(2,im),idum,idum,idum,idum,idum,sh_sc1(4,im),idum,idum,idum,idum,sh_sc1(5,im)
      enddo     
      if (nsm1.le.0) goto 8
      call skip_hashes(11)
      do im=1,nsm1
      read(11,*) idum,idum,idum,idum,spsh_sc1(1,im),spsh_sc1(3,im),spsh_sc1(2,im),idum,idum,idum,idum,idum,spsh_sc1(4,im),idum,idum,idum,idum,spsh_sc1(5,im)
      enddo
      
8      close(11)
 
      write(*,*)
      write(*,*) 'Shutter scores located and read in'
      write(*,*)

      goto 1011
1010  close(11)
      write(*,*)
      !write(*,*) ' Shutter score file not found - shutter regions will not be color-coded'  
      !write(*,*) '      Run m_check_shutter_scores module first if needed     ' 
      write(*,*)
      shutter_scores=.false.
1011  continue

        
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
      s_infile='./trial_'//rn(1:2)//'/m_pick_output/pair_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/sky_shutters.txt'

      nsm1=0

      if (.not.skypad) goto 17

      open(11,file=s_infile, status='old',access='sequential', form='formatted', err=17)
      
      read(11,*) !skip header
15     nsm1=nsm1+1
      read(11,*,end=16) idum,ksm1(nsm1),ism1(nsm1),jsm1(nsm1)
      goto 15
16     nsm1=nsm1-1

17     continue

      close(11)

      write(ipap,"(i0)") nn_p(2)
      write(jpap,"(i0)") mtarget
      s_infile='./trial_'//rn(1:2)//'/m_pick_output/pair_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/sky_shutters.txt'

      nsm2=0

      if (.not.skypad) goto 27

      open(11,file=s_infile, status='old',access='sequential', form='formatted', err=27)
      
      read(11,*) !skip header
25     nsm2=nsm2+1
      read(11,*,end=26) idum,ksm2(nsm2),ism2(nsm2),jsm2(nsm2)
      goto 25
26     nsm2=nsm2-1

      write(*,*)
      write(*,*) 'Sky shutter lists located and read in'

27     continue

      close(11)

      ip=1
      
      write(ipap,"(i0)") nn_p(ip)
      write(jpap,"(i0)") mtarget
      sc_infile='./trial_'//rn(1:2)//'/m_check_output/pair_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/slitlet_contamination_scores.txt'

      open(11,file=sc_infile, status='old',access='sequential', form='formatted',err=1020)
      call skip_hashes(11)
      do im=1,nm(ip)
      read(11,*) idum,idum,idum,idum,idum,idum,idum,sh_sc1(1,im),sh_sc1(3,im),sh_sc1(2,im),idum,idum,idum,idum,idum,sh_sc1(4,im),idum,idum,idum,idum,sh_sc1(5,im)
      enddo
      if (nsm1.le.0) goto 28
      call skip_hashes(11)
      do im=1,nsm1
      read(11,*) idum,idum,idum,idum,spsh_sc1(1,im),spsh_sc1(3,im),spsh_sc1(2,im),idum,idum,idum,idum,idum,spsh_sc1(4,im),idum,idum,idum,idum,spsh_sc1(5,im)
      enddo
28     close(11)

      goto 1021
1020  close(11)
      write(*,*)
      !write(*,*) ' Shutter score file not found - shutter regions will not be color-coded'  
      !write(*,*) '      Run m_check_shutter_scores module first if needed     ' 
      write(*,*)
      shutter_scores=.false.
1021  continue



      ip=2

      write(ipap,"(i0)") nn_p(ip)
      write(jpap,"(i0)") mtarget
      sc_infile='./trial_'//rn(1:2)//'/m_check_output/pair_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/slitlet_contamination_scores.txt'
      open(11,file=sc_infile, status='old',access='sequential', form='formatted',err=1022)
      call skip_hashes(11)
      do im=1,nm(ip)
      read(11,*) idum,idum,idum,idum,idum,idum,idum,sh_sc2(1,im),sh_sc2(3,im),sh_sc2(2,im),idum,idum,idum,idum,idum,sh_sc2(4,im),idum,idum,idum,idum,sh_sc2(5,im)
      enddo
      if (nsm2.le.0) goto 29
      call skip_hashes(11)
      do im=1,nsm2
      read(11,*) idum,idum,idum,idum,spsh_sc2(1,im),spsh_sc2(3,im),spsh_sc2(2,im),idum,idum,idum,idum,idum,spsh_sc2(4,im),idum,idum,idum,idum,spsh_sc2(5,im)
      enddo
29     close(11)


      write(*,*)
      write(*,*) 'Shutter scores located and read in'
      write(*,*)
 
      goto 1023
1022  close(11)
      write(*,*)
      !write(*,*) ' Shutter score file not found - shutter regions will not be color-coded'  
      !write(*,*) '      Run m_check_shutter_scores module first if needed     ' 
      write(*,*)
      shutter_scores=.false.
1023  continue


      endif !n_dither=2

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
      s_infile='./trial_'//rn(1:2)//'/m_pick_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/sky_shutters.txt'

      nsm1=0

      if (.not.skypad) goto 37

      open(11,file=s_infile, status='old',access='sequential', form='formatted', err=37)
      
      read(11,*) !skip header
35     nsm1=nsm1+1
      read(11,*,end=36) idum,ksm1(nsm1),ism1(nsm1),jsm1(nsm1)
      goto 35
36     nsm1=nsm1-1

37     continue

       close(11)

      write(ipap,"(i0)") nn_p(2)
      write(jpap,"(i0)") mtarget
      s_infile='./trial_'//rn(1:2)//'/m_pick_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/sky_shutters.txt'

      nsm2=0

      if (.not.skypad) goto 47

      open(11,file=s_infile, status='old',access='sequential', form='formatted', err=47)
      
      read(11,*) !skip header
45     nsm2=nsm2+1
      read(11,*,end=46) idum,ksm2(nsm2),ism2(nsm2),jsm2(nsm2)
      goto 45
46     nsm2=nsm2-1

47     continue

       close(11)

      write(ipap,"(i0)") nn_p(3)
      write(jpap,"(i0)") mtarget
      s_infile='./trial_'//rn(1:2)//'/m_pick_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/sky_shutters.txt'

      nsm3=0

      if (.not.skypad) goto 57

      open(11,file=s_infile, status='old',access='sequential', form='formatted',err=57)
      
      read(11,*) !skip header
55     nsm3=nsm3+1
      read(11,*,end=56) idum,ksm3(nsm3),ism3(nsm3),jsm3(nsm3)
      goto 55
56     nsm3=nsm3-1

      write(*,*)
      write(*,*) 'Sky shutter lists located and read in'

57     continue

      close(11)

      ip=1

      write(ipap,"(i0)") nn_p(ip)
      write(jpap,"(i0)") mtarget
      sc_infile='./trial_'//rn(1:2)//'/m_check_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/slitlet_contamination_scores.txt'

      open(11,file=sc_infile, status='old',access='sequential', form='formatted',err=1030)
      call skip_hashes(11)
      do im=1,nm(ip)
      read(11,*) idum,idum,idum,idum,idum,idum,idum,sh_sc1(1,im),sh_sc1(3,im),sh_sc1(2,im),idum,idum,idum,idum,idum,sh_sc1(4,im),idum,idum,idum,idum,sh_sc1(5,im)
      enddo
      if (nsm1.le.0) goto 58
      call skip_hashes(11)
      do im=1,nsm1
      read(11,*) idum,idum,idum,idum,spsh_sc1(1,im),spsh_sc1(3,im),spsh_sc1(2,im),idum,idum,idum,idum,idum,spsh_sc1(4,im),idum,idum,idum,idum,spsh_sc1(5,im)
!        print *,im,(spsh_sc1(j,im),j=1,5)      
      enddo

58      close(11)

      goto 1031
1030  close(11)
      write(*,*)
      !write(*,*) ' Shutter score file not found - shutter regions will not be color-coded'  
      !write(*,*) '      Run m_check_shutter_scores module first if needed     ' 
      write(*,*)
      shutter_scores=.false.
1031  continue


      ip=2

      write(ipap,"(i0)") nn_p(ip)
      write(jpap,"(i0)") mtarget
      sc_infile='./trial_'//rn(1:2)//'/m_check_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/slitlet_contamination_scores.txt'

      open(11,file=sc_infile, status='old',access='sequential', form='formatted',err=1032)
      call skip_hashes(11)
      do im=1,nm(ip)
      read(11,*) idum,idum,idum,idum,idum,idum,idum,sh_sc2(1,im),sh_sc2(3,im),sh_sc2(2,im),idum,idum,idum,idum,idum,sh_sc2(4,im),idum,idum,idum,idum,sh_sc2(5,im)
      enddo
      if (nsm2.le.0) goto 59
      call skip_hashes(11)
      do im=1,nsm2
      read(11,*) idum,idum,idum,idum,spsh_sc2(1,im),spsh_sc2(3,im),spsh_sc2(2,im),idum,idum,idum,idum,idum,spsh_sc2(4,im),idum,idum,idum,idum,spsh_sc2(5,im)
      enddo
59      close(11)

      goto 1033
1032  close(11)
      write(*,*)
      !write(*,*) ' Shutter score file not found - shutter regions will not be color-coded'  
      !write(*,*) '      Run m_check_shutter_scores module first if needed     ' 
      write(*,*)
      shutter_scores=.false.
1033  continue



      ip=3

      write(ipap,"(i0)") nn_p(ip)
      write(jpap,"(i0)") mtarget
      sc_infile='./trial_'//rn(1:2)//'/m_check_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))//'/slitlet_contamination_scores.txt'

      open(11,file=sc_infile, status='old',access='sequential', form='formatted',err=1034)
      call skip_hashes(11)
      do im=1,nm(ip)
      read(11,*) idum,idum,idum,idum,idum,idum,idum,sh_sc3(1,im),sh_sc3(3,im),sh_sc3(2,im),idum,idum,idum,idum,idum,sh_sc3(4,im),idum,idum,idum,idum,sh_sc3(5,im)
      enddo
      if (nsm3.le.0) goto 60
      call skip_hashes(11)
      do im=1,nsm3
      read(11,*) idum,idum,idum,idum,spsh_sc3(1,im),spsh_sc3(3,im),spsh_sc3(2,im),idum,idum,idum,idum,idum,spsh_sc3(4,im),idum,idum,idum,idum,spsh_sc3(5,im)
      enddo
60      close(11)

      write(*,*)
      write(*,*) 'Shutter scores located and read in'
      write(*,*)

      goto 1035
1034  close(11)
      write(*,*)
      !write(*,*) ' Shutter score file not found - shutter regions will not be color-coded'  
      !write(*,*) '      Run m_check_shutter_scores module first if needed     ' 
      write(*,*)
      shutter_scores=.false.
1035  continue
 

      endif !n_dither=3


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


! Write out DS9 target region files
       
       do ip=n_dither,1,-1

      ntt=nm(ip)
      
      do i=1,nm(ip)
      ktt(i)=kt(ip,i)
      itt(i)=it(ip,i)
      jtt(i)=jt(ip,i)
      idtt(i)=id_cat(ip,i)
      enddo
     
       if (ip.eq.1) filename=trim(adjustl(outdir_1))//'/ds9_slit_regions.reg'
       if (ip.eq.2) filename=trim(adjustl(outdir_2))//'/ds9_slit_regions.reg'
       if (ip.eq.3) filename=trim(adjustl(outdir_3))//'/ds9_slit_regions.reg'

       if (ip.eq.1) filename2=trim(adjustl(outdir_1))//'/ds9_targ_regions.reg'
       if (ip.eq.2) filename2=trim(adjustl(outdir_2))//'/ds9_targ_regions.reg'
       if (ip.eq.3) filename2=trim(adjustl(outdir_3))//'/ds9_targ_regions.reg'


       if (ip.eq.1) call write_regions(ra_p(ip_m(ip)),dec_p(ip_m(ip)),pa_ap_p(ip_m(ip)),ntt,idtt,ktt,itt,jtt,nsm1,ksm1,ism1,jsm1,filename,filename2,sh_sc1,spsh_sc1,shutter_scores)       
       if (ip.eq.2) call write_regions(ra_p(ip_m(ip)),dec_p(ip_m(ip)),pa_ap_p(ip_m(ip)),ntt,idtt,ktt,itt,jtt,nsm2,ksm2,ism2,jsm2,filename,filename2,sh_sc2,spsh_sc2,shutter_scores)
       if (ip.eq.3) call write_regions(ra_p(ip_m(ip)),dec_p(ip_m(ip)),pa_ap_p(ip_m(ip)),ntt,idtt,ktt,itt,jtt,nsm3,ksm3,ism3,jsm3,filename,filename2,sh_sc3,spsh_sc3,shutter_scores)

       enddo


       write(*,*)
       call cpu_time(t_stop)
       write(*,'(a)')      '---------------------------------------------------------'
       write(*,'(a,f7.1)') ' CPU time used by m_check_regions module [s]:',t_stop-t_start
       write(*,'(a)')      '---------------------------------------------------------'
       write(*,*)


      end program m_check_regions

!------------------------------------------------------------

 
include './reference_files/front_transforms_updt.f90'

include './reference_files/back_transforms_updt.f90'

include './reference_files/shutter_routines_new.f90'

include './reference_files/misc_routines_new.f90'

!------------------------------------------------------------

      subroutine file_entry_m_check_regions()
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
      read(9,*) !imname
      
      
      call skip_hashes(9)
      read(9,*) !segname

      
      call skip_hashes(9)
      read(9,*) !cra,cdec,cpa_ap


      call skip_hashes(9)
      read(9,*) !szone_x,szone_y


! work out pointing reference point coordinates in V2,V3 on-sky system

      call msa_to_sky_bd(xm_ref,ym_ref,ax_refv2,ay_refv3) !v2,V3

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

      end subroutine file_entry_m_check_regions

!------------------------------------------------------------------------

      subroutine write_regions(ccra,ccdec,ccpa_ap,ntt,idtt,ktt,itt,jtt,ns,ks,is,js,reg_file,reg_file2,sh_sc,spsh_sc,shutter_scores)
      use derived_parameters, only : phi
      use target_cat, only : ntarg,ra_t,dec_t,id_t
      use m_list, only : n_max_m
      implicit none
      real ccpa_ap,ccpa_v3
      real*8 ccra,ccdec
      integer kk,ii,jj,it
      character*120 reg_file
      character*120 reg_file2
      character*13 label
      logical shutter_scores

      integer ntt
      integer ktt(n_max_m),itt(n_max_m),jtt(n_max_m),idtt(n_max_m)
      integer sh_sc(5,n_max_m)
      integer spsh_sc(5,n_max_m)

      integer ns
      integer, parameter :: max_spare_m=250 !Maximum number of non-overlapping unused slitlets on MSA
      integer ks(max_spare_m),is(max_spare_m),js(max_spare_m)

! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar

 
       open(11,file=reg_file,status='replace',access='sequential', form='formatted')    

       write(11,'(a)') '# Region file format: DS9 version 4.1 written by NIRSpec IPA/eMPT'
       write(11,'(a)') '# Output by NIRSpec IPA/eMPT'
       write(11,'(a)') '#'
       write(11,'(a)') '#'
      if (shutter_scores) then
       write(11,'(a)') '# Shutter content (additive) scoring and color scheme:'
       write(11,'(a)') '#'
       write(11,'(a)') '# Shutters nominally containing intended targets:'
       write(11,'(a)') '#'
       write(11,'(a)') '#  0: No object present in shutter                - yellow'
       write(11,'(a)') '#  1: Intended target present in shutter          - green'
       write(11,'(a)') '#  4: One non-target object present in shutter    - red'
       write(11,'(a)') '#  8: Two non-target objects present in shutter   - red'
       write(11,'(a)') '# 16: Three non-target objects present in shutter - red'      
       write(11,'(a)') '#'
       write(11,'(a)') '# Shutters not containing intended targets:'
       write(11,'(a)') '#'
       write(11,'(a)') '#  0: No object present in shutter                - green'
       write(11,'(a)') '#  2: Spill-over from target present in shutter   - yellow'     
       write(11,'(a)') '#  4: One non-target object present in shutter    - red'
       write(11,'(a)') '#  8: Two non-target objects present in shutter   - red' 
       write(11,'(a)') '# 16: Three non-target objects present in shutter - red'       
       write(11,'(a)') '#'
       write(11,'(a)') '# Spare slitlets:'
       write(11,'(a)') '#'
       write(11,'(a)') '#  0: No object present in shutter                - green'   
       write(11,'(a)') '#  4: Object present in shutter                   - red'
       write(11,'(a)') '#'
       write(11,'(a)') '# 32: Empty shutter with anomalous faint or high flux - cyan'
       write(11,'(a)') '#'
      endif !shutter_scores
       write(11,'(a)') '#'
       write(11,'(a)') 'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=0 move=0 delete=1 rotate=0 include=1 source=1'
       write(11,'(a)') 'fk5'


      call load_msa_dimensions_bd
           
    
      ccpa_v3=ccpa_ap-(180.-phi)
      if (ccpa_v3.lt.0.) ccpa_v3=ccpa_v3+360.
      if (ccpa_v3.gt.360.) ccpa_v3=ccpa_v3-360.


       do it=1,ntt

! shutter kk,ii,jj
  
       kk=ktt(it)
       ii=itt(it)
       jj=jtt(it)
       
      write(label,"(i0)") idtt(it)
      label=trim(adjustl(label))
       
      call write_shutter_region_label(ccra,ccdec,ccpa_v3,kk,ii,jj,label,sh_sc(1,it))  !1   
      call write_shutter_region(ccra,ccdec,ccpa_v3,kk,ii,jj+1,sh_sc(2,it))     !2
      call write_shutter_region(ccra,ccdec,ccpa_v3,kk,ii,jj-1,sh_sc(3,it))     !3
      call write_shutter_region(ccra,ccdec,ccpa_v3,kk,ii,jj+2,sh_sc(4,it))     !4 
      call write_shutter_region(ccra,ccdec,ccpa_v3,kk,ii,jj-2,sh_sc(5,it))     !5
           
      enddo !it


! spare slitlets

       do it=1,ns

! shutter kk,ii,jj
  
       kk=ks(it)
       ii=is(it)
       jj=js(it)
       
      write(label,"(i0)") it
      label='B'//trim(adjustl(label))
      label=trim(adjustl(label))
      call write_shutter_region_label_spare(ccra,ccdec,ccpa_v3,kk,ii,jj,label,spsh_sc(1,it))      
      call write_shutter_region(ccra,ccdec,ccpa_v3,kk,ii,jj+1,spsh_sc(2,it))      
      call write_shutter_region(ccra,ccdec,ccpa_v3,kk,ii,jj-1,spsh_sc(3,it))      
      call write_shutter_region(ccra,ccdec,ccpa_v3,kk,ii,jj+2,spsh_sc(4,it))      
      call write_shutter_region(ccra,ccdec,ccpa_v3,kk,ii,jj-2,spsh_sc(5,it))      
      
      enddo !it

      
      close(11)

! Write catalog region file


       open(11,file=reg_file2,status='replace',access='sequential', form='formatted')    
 
       write(11,'(a)') '# Region file format: DS9 version 4.1 written by NIRSpec IPA/eMPT'
       write(11,'(a)') '# Output by NIRSpec IPA/eMPT'
       write(11,'(a)') 'global color=red dashlist=8 3 width=1 font="helvetica 8 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=0 move=0 delete=1 rotate=0 include=1 source=1'
       write(11,'(a)') 'fk5'


      do it=1,ntarg
       write(11,'(a,f11.7,a,f11.7,a,i0,a)') 'point(',ra_t(it),',',dec_t(it),') # point=cross text={',id_t(it),'}'
      enddo
  
      close(11)
    
      return
      
      end subroutine write_regions
      
!-----------------------------------------------------------------------------------------------    
       
      subroutine write_shutter_region(ccra,ccdec,ccpa_v3,kk,ii,jj,isc)
      use config_parameters, only : dva_mag
      use derived_parameters, only : ax_refv2,ay_refv3
      implicit none
      real*8 ccra,ccdec
      real ccpa_v3
      real x1(4),y1(4)
      real*8 ra_p(4),dec_p(4)
      real a_x,a_y,theta_x,theta_y
      integer kk,ii,jj,i,isc
      character*8 color
 ! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar

      if (isc.eq.0) color='green'
      if (isc.eq.2) color='yellow' 
      if (isc.ge.4) color='red'
      if (isc.eq.32) color='cyan'

      
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
        call msa_to_sky_bd(x1(i),y1(i),theta_x,theta_y) !v2,v3 arcsec 
        theta_x=(theta_x-ax_refv2)/dva_mag
        theta_y=(theta_y-ay_refv3)/dva_mag
        call roll_v3(ccpa_v3,theta_x,theta_y,a_x,a_y) 
        call coordtang_single(ccra,ccdec,a_x,a_y,ra_p(i),dec_p(i)) !ra,dec 
      enddo
      
 
      write(11,'(a,f12.7,a,f12.7,a,f12.7,a,f12.7,a,f12.7,a,f12.7,a,f12.7,a,f12.7,a)') 'polygon(',ra_p(1),',',dec_p(1),',',ra_p(2),',',dec_p(2),',',ra_p(3),',',dec_p(3),',',ra_p(4),',',dec_p(4),') # color={'//trim(color)//'}'

      
      return
      end subroutine write_shutter_region

!-----------------------------------------------------------------------------------------------    
      
      subroutine write_shutter_region_label(ccra,ccdec,ccpa_v3,kk,ii,jj,label,isc)
      use config_parameters, only : dva_mag
      use derived_parameters, only : ax_refv2,ay_refv3
      implicit none
      real*8 ccra,ccdec
      real ccpa_v3
      real x1(4),y1(4)
      real*8 ra_p(4),dec_p(4)
      real a_x,a_y,theta_x,theta_y
      integer kk,ii,jj,i,isc
      character*13 label
      character*8 color
 ! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
      
      if (isc.eq.1) color='green'
      if (isc.eq.0) color='yellow'
      if (isc.gt.1) color='red'
      
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
        call msa_to_sky_bd(x1(i),y1(i),theta_x,theta_y) !arcsec shutter center V
        theta_x=(theta_x-ax_refv2)/dva_mag
        theta_y=(theta_y-ay_refv3)/dva_mag
        call roll_v3(ccpa_v3,theta_x,theta_y,a_x,a_y) 
        call coordtang_single(ccra,ccdec,a_x,a_y,ra_p(i),dec_p(i)) !ra,dec 
      enddo
      

      write(11,'(a,f12.7,a,f12.7,a,f12.7,a,f12.7,a,f12.7,a,f12.7,a,f12.7,a,f12.7,a)') 'polygon(',ra_p(1),',',dec_p(1),',',ra_p(2),',',dec_p(2),',',ra_p(3),',',dec_p(3),',',ra_p(4),',',dec_p(4),') # text={'//trim(label)//'} color={'//trim(color)//'}'

      
      return
      end subroutine write_shutter_region_label

!-----------------------------------------------------------------------------------------------    
     
      subroutine write_shutter_region_label_spare(ccra,ccdec,ccpa_v3,kk,ii,jj,label,isc)
      use config_parameters, only : dva_mag
      use derived_parameters, only : ax_refv2,ay_refv3
      implicit none
      real ccpa_v3
      real x1(4),y1(4)
      real*8 ccra,ccdec
      real*8 ra_p(4),dec_p(4)
      real a_x,a_y,theta_x,theta_y
      integer kk,ii,jj,i,isc
      character*13 label
      character*8 color
 ! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
      
      if (isc.eq.0) color='green'
      if (isc.gt.0) color='red'
      
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
        call msa_to_sky_bd(x1(i),y1(i),theta_x,theta_y) !arcsec shutter center V
        theta_x=(theta_x-ax_refv2)/dva_mag
        theta_y=(theta_y-ay_refv3)/dva_mag
        call roll_v3(ccpa_v3,theta_x,theta_y,a_x,a_y) 
        call coordtang_single(ccra,ccdec,a_x,a_y,ra_p(i),dec_p(i)) !ra,dec 
      enddo
      

      write(11,'(a,f12.7,a,f12.7,a,f12.7,a,f12.7,a,f12.7,a,f12.7,a,f12.7,a,f12.7,a)') 'polygon(',ra_p(1),',',dec_p(1),',',ra_p(2),',',dec_p(2),',',ra_p(3),',',dec_p(3),',',ra_p(4),',',dec_p(4),') # text={'//trim(label)//'} color={'//trim(color)//'}'

      
      return
      end subroutine write_shutter_region_label_spare

!-----------------------------------------------------------------------------------------------    
