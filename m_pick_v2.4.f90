! V2.0  max_failo increased to 50  19/03/22  
! v2.1 pointing_summary.txt output file changed to reflect MPT roll-assigment approach 
! v2.1 rewritten to handle filter-dependent FORE transformation
! v2.2 increased ID field to 13 digits throughout   
! v2.3 text in summary file changed to emphasize ASSIGNED PA_AP 
! v2.4 swapped nodded pointing order in summary table

     include 'empt_modules.f90'

!-------------------------------------------------------

      program m_pick
      use path_to_ref_files
      use config_parameters
      use derived_parameters
      use target_cat
      use pointings_and_groups, only : n_p,no_p,ra_p,dec_p,pa_ap_p
      use m_list
      use histograms
      implicit none
      real t_start,t_stop
      character*2 rn
      character*5 ipap,jpap
      character*70 outdir_1
      character*70 outdir_2
      character*70 outdir_3
      character*70 filename
      character*70 outfile
      character*70 k_list_mod_file
      integer nclass_conf
      integer, allocatable :: nsingles,npairs,ntriples
      integer, allocatable :: nsinglecount,npaircount,ntriplecount
      integer, allocatable :: ip1,ip2,ip3
      integer, allocatable :: np1,np2,np3
      integer, allocatable :: im1,im2,im3
      integer, allocatable :: nm1,nm2,nm3
      integer, allocatable :: id01(:),id1(:),cl1(:)
      integer, allocatable ::  kt1(:),it1(:),jt1(:)
      real, allocatable ::  rx1(:),ry1(:)
      real*8, allocatable ::  ra_t1(:),dec_t1(:)
      integer, allocatable ::  id02(:),id2(:),cl2(:)
      integer, allocatable ::  kt2(:),it2(:),jt2(:)
      real, allocatable ::  rx2(:),ry2(:)
      real*8, allocatable ::  ra_t2(:),dec_t2(:)
      integer, allocatable ::  id03(:),id3(:),cl3(:)
      integer, allocatable ::  kt3(:),it3(:),jt3(:)
      real, allocatable ::  rx3(:),ry3(:)
      real*8, allocatable ::  ra_t3(:),dec_t3(:)

      integer  i,idum,ip
      integer ik,nk,nk1,nk2,nk3,prcl
      real*8 ra_p1,dec_p1
      real*8 ra_p2,dec_p2
      real*8 ra_p3,dec_p3
      real pa_ap_p1,pa_ap_p2,pa_ap_p3
      real*8 cran1,cdecn1,cran2,cdecn2


      integer, allocatable ::  ktn1(:),itn1(:),jtn1(:)
      real, allocatable ::  rxn1(:),ryn1(:)
      integer, allocatable ::  ktn2(:),itn2(:),jtn2(:)
      real, allocatable ::  rxn2(:),ryn2(:)
      integer nsm1,nsm2,nsm3
 
      integer, allocatable :: id0_tmp(:),id_tmp(:),cl_tmp(:),kt_tmp(:),it_tmp(:),jt_tmp(:)
     
      integer, parameter :: max_spare_m=250 !Maximum number of non-overlapping unused slitlets on MSA
      integer ksm1(max_spare_m),ism1(max_spare_m),jsm1(max_spare_m)
      integer ksm2(max_spare_m),ism2(max_spare_m),jsm2(max_spare_m)
      integer ksm3(max_spare_m),ism3(max_spare_m),jsm3(max_spare_m)
     
! - spectral parameters      
!      integer n_disp
      real lambda1,lambda2
      character*6 gname  
      common/disp/lambda1,lambda2,gname
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
! - prism seperations
      integer prism_sep_p(4,365,171),prism_sep_m(4,365,171)
      common/prism_length/prism_sep_p,prism_sep_m


      real*8 ra_p0,dec_p0
      real pa_ap_p0


      call cpu_time(t_start)

! Note priority class for each target assumed to be the same for all n_dither dithers
! assigned priority class in order of increasing dither

      write(*,*)
      write(*,*) '-----------------------------------------------------------------------------------------'
      write(*,*) '--- m_pick -- eMPT General Dithered m-list Parsing Module Version 2.4 29 October 2022 ---'
      write(*,*) '-----------------------------------------------------------------------------------------'
      write(*,*)

      call load_derived_parameters

      call file_entry_m_pick()
      
      write(rn,'(i2)') ntile
      rn=trim(adjustl(rn))
      if (ntile.lt.10) rn='0'//rn(1:1)

      call system('mkdir trial_'//rn(1:2)//'/m_pick_output')

      if (n_disp.eq.7) call read_prism_separations


      call pointings_read()
      
      
      if (n_dither.eq.1) m_list_file='trial_'//rn(1:2)//'/m_make_output/single_m_list.txt'
      if (n_dither.eq.2) m_list_file='trial_'//rn(1:2)//'/m_make_output/pair_m_list.txt'
      if (n_dither.eq.3) m_list_file='trial_'//rn(1:2)//'/m_make_output/triple_m_list.txt'
   
      call m_list_size_read()
      nclass_conf=nclass

      allocate (ip1,ip2,ip3,np1,np2,np3,im1,im2,im3,nm1,nm2,nm3)
      allocate (id01(n_max_m),id1(n_max_m),cl1(n_max_m))
      allocate (id02(n_max_m),id2(n_max_m),cl2(n_max_m))
      allocate (id03(n_max_m),id3(n_max_m),cl3(n_max_m))


      k_list_mod_file='trial_'//rn(1:2)//'/k_clean_output/k_list_mod.txt'

      call catalog_read()
      nclass=nclass_conf

      
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

! Make catalog Priority Class histogram

      allocate (hist_k1(nclass),hist_k2(nclass),hist_k3(nclass))
      allocate (hist_m1(nclass),hist_m2(nclass),hist_m3(nclass))
      allocate (hist_cat(nclass),hist_com(nclass))
      allocate (avexp(nclass))
      
      hist_cat=0
      
      do i=1,ntarg
      if (pri_t(i).gt.0) hist_cat(pri_t(i))=hist_cat(pri_t(i))+1
      enddo

! Null k_list and m_list histograms

      hist_m1=0
      hist_m2=0
      hist_m3=0
      hist_k1=0
      hist_k2=0
      hist_k3=0
 
      write(*,*)
      write(*,*) 'Constructing catalog(s) of observed targets'

! Read in observed target lists from target optimal m-list and recover target coordinates from input catalog

!***************************************************

      if (n_dither.eq.1) then

      allocate (nsingles,nsinglecount)
 
      open(10,file=m_list_file,status='old',access='sequential', form='formatted')


      read(10,*) nsingles

      allocate (kt1(n_max_m),it1(n_max_m),jt1(n_max_m))
      allocate (rx1(n_max_m),ry1(n_max_m))
      allocate (ra_t1(n_max_m),dec_t1(n_max_m))


      do i=1,mtarget

      read(10,*) nsinglecount
      read(10,*) ip1
      read(10,*) np1


      read(10,*) nm1
      do im1=1,nm1
      read(10,*) id01(im1),id1(im1),cl1(im1),kt1(im1),it1(im1),jt1(im1),rx1(im1),ry1(im1)
      ra_t1(im1)=ra_t(id01(im1))
      dec_t1(im1)=dec_t(id01(im1))
      if (id1(im1).ne.id_t(id01(im1))) then
      write(*,*) 'ID mismatch between m_list and target list:',im1,id1(im1),id_t(id01(im1))
      stop
      endif
      enddo!im1

      enddo!i

      close(10)


 
! fill in m_list Priority Class histogram

      do im1=1,nm1
      if (cl1(im1).gt.0) hist_m1(cl1(im1))=hist_m1(cl1(im1))+1
      enddo

 
      endif !n_dither=1

!***************************************************

      if (n_dither.eq.2) then

       open(10,file=m_list_file,status='old',access='sequential', form='formatted')

       allocate (npairs,npaircount)

      read(10,*) npairs

      allocate (kt1(n_max_m),it1(n_max_m),jt1(n_max_m))
      allocate (rx1(n_max_m),ry1(n_max_m))
      allocate (ra_t1(n_max_m),dec_t1(n_max_m))
      allocate (kt2(n_max_m),it2(n_max_m),jt2(n_max_m))
      allocate (rx2(n_max_m),ry2(n_max_m))
      allocate (ra_t2(n_max_m),dec_t2(n_max_m))

      do i=1,mtarget

      read(10,*) npaircount
      read(10,*) ip1,ip2
      read(10,*) np1,np2

 
      read(10,*) nm1
      do im1=1,nm1
      read(10,*) id01(im1),id1(im1),cl1(im1),kt1(im1),it1(im1),jt1(im1),rx1(im1),ry1(im1)
      ra_t1(im1)=ra_t(id01(im1))
      dec_t1(im1)=dec_t(id01(im1))
      if (id1(im1).ne.id_t(id01(im1))) then
      write(*,*) 'ID mismatch between m_list and target list:',im1,id1(im1),id_t(id01(im1))
      stop
      endif
      enddo!im1

      read(10,*) nm2
      do im2=1,nm2
      read(10,*) id02(im2),id2(im2),cl2(im2),kt2(im2),it2(im2),jt2(im2),rx2(im2),ry2(im2)
      ra_t2(im2)=ra_t(id02(im2))
      dec_t2(im2)=dec_t(id02(im2))
      if (id2(im2).ne.id_t(id02(im2))) then
      write(*,*) 'ID mismatch between m_list and target list:',im2,id2(im2),id_t(id02(im2))
      stop
      endif
      enddo!im2

      enddo !i

      close(10)

! fill in m_list Priority Class histograms

      do im1=1,nm1
      if (cl1(im1).gt.0) hist_m1(cl1(im1))=hist_m1(cl1(im1))+1
      enddo

      do im2=1,nm2
      if (cl2(im2).gt.0) hist_m2(cl2(im2))=hist_m2(cl2(im2))+1
      enddo

      endif !n_dither=2

!***************************************************

      if (n_dither.eq.3) then


       open(10,file=m_list_file,status='old',access='sequential', form='formatted')

      allocate(ntriples,ntriplecount)

      read(10,*) ntriples

      allocate (kt1(n_max_m),it1(n_max_m),jt1(n_max_m))
      allocate (rx1(n_max_m),ry1(n_max_m))
      allocate (ra_t1(n_max_m),dec_t1(n_max_m))
      allocate (kt2(n_max_m),it2(n_max_m),jt2(n_max_m))
      allocate (rx2(n_max_m),ry2(n_max_m))
      allocate (ra_t2(n_max_m),dec_t2(n_max_m))
      allocate (kt3(n_max_m),it3(n_max_m),jt3(n_max_m))
      allocate (rx3(n_max_m),ry3(n_max_m))
      allocate (ra_t3(n_max_m),dec_t3(n_max_m))

      do i=1,mtarget

      read(10,*) ntriplecount
      read(10,*) ip1,ip2,ip3
      read(10,*) np1,np2,np3

      read(10,*) nm1
      do im1=1,nm1
      read(10,*) id01(im1),id1(im1),cl1(im1),kt1(im1),it1(im1),jt1(im1),rx1(im1),ry1(im1)
      ra_t1(im1)=ra_t(id01(im1))
      dec_t1(im1)=dec_t(id01(im1))
      if (id1(im1).ne.id_t(id01(im1))) then
      write(*,*) 'ID mismatch between m_list and target list:',im1,id1(im1),id_t(id01(im1))
      stop
      endif
      enddo!im1

      read(10,*) nm2
      do im2=1,nm2
      read(10,*) id02(im2),id2(im2),cl2(im2),kt2(im2),it2(im2),jt2(im2),rx2(im2),ry2(im2)
      ra_t2(im2)=ra_t(id02(im2))
      dec_t2(im2)=dec_t(id02(im2))
      if (id2(im2).ne.id_t(id02(im2))) then
      write(*,*) 'ID mismatch between m_list and target list:',im2,id2(im2),id_t(id02(im2))
      stop
      endif
      enddo!im2

      read(10,*) nm3
      do im3=1,nm3
      read(10,*) id03(im3),id3(im3),cl3(im3),kt3(im3),it3(im3),jt3(im3),rx3(im3),ry3(im3)
      ra_t3(im3)=ra_t(id03(im3))
      dec_t3(im3)=dec_t(id03(im3))
      if (id3(im3).ne.id_t(id03(im3))) then
      write(*,*) 'ID mismatch between m_list and target list:',im3,id3(im3),id_t(id03(im3))
      stop
      endif
      enddo!im3

      enddo !i

      close(10)

! fill in m_list Priority Class histograms

      do im1=1,nm1
      if (cl1(im1).gt.0) hist_m1(cl1(im1))=hist_m1(cl1(im1))+1
      enddo

      do im2=1,nm2
      if (cl2(im2).gt.0) hist_m2(cl2(im2))=hist_m2(cl2(im2))+1
      enddo

      do im3=1,nm3
      if (cl3(im3).gt.0) hist_m3(cl3(im3))=hist_m3(cl3(im3))+1
      enddo


      endif !n_dither=3

!***************************************************

      write(*,*)
      write(*,*) 'Creating output directories'


! create output directories

      if (n_dither.eq.1) then
      write(ipap,"(i0)") np1
      outdir_1='trial_'//rn(1:2)//'/m_pick_output/pointing_'//trim(adjustl(ipap))
      call system( 'mkdir '//outdir_1)
      endif

      if (n_dither.eq.2) then
      write(jpap,"(i0)") npaircount
      call system( 'mkdir trial_'//rn(1:2)//'/m_pick_output/pair_'//trim(adjustl(jpap)))
      write(ipap,"(i0)") np1
      outdir_1='trial_'//rn(1:2)//'/m_pick_output/pair_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))
      call system( 'mkdir '//outdir_1)
      write(ipap,"(i0)") np2
      outdir_2='trial_'//rn(1:2)//'/m_pick_output/pair_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))
      call system( 'mkdir '//outdir_2)
      endif

      if (n_dither.eq.3) then
      write(jpap,"(i0)") ntriplecount
      call system( 'mkdir trial_'//rn(1:2)//'/m_pick_output/triple_'//trim(adjustl(jpap)))
      write(ipap,"(i0)") np1
      outdir_1='trial_'//rn(1:2)//'/m_pick_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))
      call system( 'mkdir '//outdir_1)
      write(ipap,"(i0)") np2
      outdir_2='trial_'//rn(1:2)//'/m_pick_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))
      call system( 'mkdir '//outdir_2)
      write(ipap,"(i0)") np3
      outdir_3='trial_'//rn(1:2)//'/m_pick_output/triple_'//trim(adjustl(jpap))//'/pointing_'//trim(adjustl(ipap))
      call system( 'mkdir '//outdir_3)
      endif

!***************************************************

! fill in k-list histogram(s)

! read in pointing k space list for pointing ip1 and fill in k_list histogram

       open(9,file=k_list_mod_file,status='old',access='sequential', form='formatted')

      read(9,*) nk

      if (nk.ne.n_p) then
      write(*,*) 'Inconsistency between number of pointings specified in configuration and k-space files'
      stop
      endif

! skip up to pointing ip1

      do ip=1,ip1-1
       read(9,*) !ip,np
       read(9,*) nk
        do ik=1,nk
          read(9,*) idum
        enddo !ik loop
      enddo! ip loop

! read pointing ip1

      read(9,*) !ip,np
      read(9,*) nk1

       do ik=1,nk1
        read(9,*) idum,idum,prcl
        if (prcl.gt.0) hist_k1(prcl)=hist_k1(prcl)+1
       enddo !ik loop

       close(9)


      if (n_dither.gt.1) then

! read in pointing k space list for pointing ip2 and fill in k_list histogram

        open(9,file=k_list_mod_file,status='old',access='sequential', form='formatted')

      read(9,*) idum

! skip up to pointing ip2

      do ip=1,ip2-1
       read(9,*) !ip,np
       read(9,*) nk
        do ik=1,nk
          read(9,*) idum
        enddo !ik loop
      enddo! ip loop

! read pointing ip2

      read(9,*) !ip,np
      read(9,*) nk2

       do ik=1,nk2
        read(9,*) idum,idum,prcl
        if (prcl.gt.0) hist_k2(prcl)=hist_k2(prcl)+1
       enddo !ik loop

       close(9)

       endif !n_dither>1


! read in pointing k space list for pointing ip3 and fill in k_list histogram

      if (n_dither.gt.2) then

        open(9,file=k_list_mod_file,status='old',access='sequential', form='formatted')

      read(9,*) idum

! skip up to pointing ip3



      do ip=1,ip3-1
       read(9,*) !ip,np
       read(9,*) nk
        do ik=1,nk
          read(9,*) idum
        enddo !ik loop
      enddo! ip loop

! read pointing ip3

      read(9,*) !ip,np
      read(9,*) nk3

       do ik=1,nk3
        read(9,*) idum,idum,prcl
        if (prcl.gt.0) hist_k3(prcl)=hist_k3(prcl)+1
       enddo !ik loop

       close(9)

       endif !n_dither>2

!***************************************************


!   Recover pointing 1 ra,dec,pa_ap

      do i=1,n_p
      if (no_p(i).eq.np1) then
      ra_p1=ra_p(i)
      dec_p1=dec_p(i)
      pa_ap_p1=pa_ap_p(i)
      goto 111
      endif
      enddo
111   continue


! write out pointing 1 observed target catalog to file

      filename=trim(outdir_1)//'/observed_targets.cat'
      open(10,file=filename,status='replace',access='sequential', form='formatted')
      write(10,101) '# No   No_sub      No_cat    Pr    RA[deg]     Dec[deg]'
      do im1=1,nm1
      write(10,100) im1,id01(im1),id1(im1),cl1(im1),ra_t1(im1),dec_t1(im1)
      enddo
      close(10)

100   format(i4,x,i8,x,i12,x,i3,f13.7,x,f12.7)
101   format(a)

! redo shutter assignments for pointing 1 nods


      allocate(ktn1(n_max_m),itn1(n_max_m),jtn1(n_max_m))
      allocate(rxn1(n_max_m),ryn1(n_max_m))
      allocate(ktn2(n_max_m),itn2(n_max_m),jtn2(n_max_m))
      allocate(rxn2(n_max_m),ryn2(n_max_m))

      call calc_nods(ra_p1,dec_p1,pa_ap_p1,nm1,ra_t1,dec_t1,cran1,cdecn1,ktn1,itn1,jtn1,rxn1,ryn1,cran2,cdecn2,ktn2,itn2,jtn2,rxn2,ryn2)

! sanity check

      do i=1,nm1
      if ((kt1(i).ne.ktn1(i)).or.(kt1(i).ne.ktn2(i))) then
      write(*,*) 'k discrepancy P1 - target has moved out of shutter between nods',i
!       stop
      endif
      if ((it1(i).ne.itn1(i)).or.(it1(i).ne.itn2(i))) then
      write(*,*) 'i discrepancy P1 - target has moved out of shutter between nods',i
!       stop
      endif
      if ((jt1(i).ne.jtn1(i)-1).or.(jt1(i).ne.jtn2(i)+1)) then
      write(*,*) 'j discrepancy P1 - target has moved out of shutter between nods',i
!       stop
      endif
      enddo

      call draw_collapsed_shutter(outdir_1,nm1,rx1,ry1,rxn1,ryn1,rxn2,ryn2,raccx,raccy,cl1,id1)

      write(*,*)
      write(*,*) 'Placing Sky Background Slitlets in un-used shutters (patience please)'
      write(*,*)

      ra_p0=ra_p1
      dec_p0=dec_p1
      pa_ap_p0=pa_ap_p1

      call make_slitmask(outdir_1,nm1,kt1,it1,jt1,nsm1,ksm1,ism1,jsm1,ra_p0,dec_p0,pa_ap_p0,skypad)

      write(*,'(a,i0,a,i0,a)') ' ',nsm1,' Sky Background slitlets added to the IPA Pointing ',np1,' MSA Mask'


! append sky slitlets to target slit list for pointing 1

      if (skypad) then
 
      allocate (id0_tmp(n_max_m),id_tmp(n_max_m),cl_tmp(n_max_m),kt_tmp(n_max_m),it_tmp(n_max_m),jt_tmp(n_max_m))
 
      do i=1,nm1
       kt_tmp(i)=kt1(i)
       it_tmp(i)=it1(i)
       jt_tmp(i)=jt1(i)
       cl_tmp(i)=cl1(i)
       id_tmp(i)=id1(i)
       id0_tmp(i)=id01(i)
      enddo

      deallocate (id01,id1,cl1,kt1,it1,jt1)
      allocate (id01(nm1+nsm1),id1(nm1+nsm1),cl1(nm1+nsm1),kt1(nm1+nsm1),it1(nm1+nsm1),jt1(nm1+nsm1))
 
      do i=1,nm1
       kt1(i)=kt_tmp(i)
       it1(i)=it_tmp(i)
       jt1(i)=jt_tmp(i)
       cl1(i)=cl_tmp(i)
       id1(i)=id_tmp(i)
       id01(i)=id0_tmp(i)
      enddo

      do i=1,nsm1
       kt1(nm1+i)=ksm1(i)
       it1(nm1+i)=ism1(i)
       jt1(nm1+i)=jsm1(i)
       cl1(nm1+i)=0
       id1(nm1+i)=0
       id01(nm1+i)=0
      enddo

      deallocate (id0_tmp,id_tmp,cl_tmp,kt_tmp,it_tmp,jt_tmp)

      endif !skypad

      call draw_msa_on_sky(ra_p1,dec_p1,pa_ap_p1,nm1,id1,ra_t1,dec_t1,cl1,outdir_1)

      call write_summary(outdir_1,gname,np1,ra_p1,dec_p1,pa_ap_p1,nm1,id01,id1,cl1,kt1,it1,jt1,rx1,ry1,nsm1,hist_m1,hist_k1,hist_cat,cran1,cdecn1,ktn1,itn1,jtn1,rxn1,ryn1,cran2,cdecn2,ktn2,itn2,jtn2,rxn2,ryn2)

!--------------

      if (n_dither.gt.1) then


!   Recover pointing 2 ra,dec,pa_ap

      do i=1,n_p
      if (no_p(i).eq.np2) then
      ra_p2=ra_p(i)
      dec_p2=dec_p(i)
      pa_ap_p2=pa_ap_p(i)
      goto 222
      endif
      enddo
222   continue


! write out pointing 2 observed target catalog to file

      filename=trim(outdir_2)//'/observed_targets.cat'
      open(10,file=filename,status='replace',access='sequential', form='formatted')
      write(10,101) '# No No_sub No_cat  Pr    RA[deg]     Dec[deg]'
      do im2=1,nm2
      write(10,100) im2,id02(im2),id2(im2),cl2(im2),ra_t2(im2),dec_t2(im2)
      enddo
      close(10)

! redo shutter assignments for pointing 2 nods

      call calc_nods(ra_p2,dec_p2,pa_ap_p2,nm2,ra_t2,dec_t2,cran1,cdecn1,ktn1,itn1,jtn1,rxn1,ryn1,cran2,cdecn2,ktn2,itn2,jtn2,rxn2,ryn2)

! sanity check

      do i=1,nm2
      if ((kt2(i).ne.ktn1(i)).or.(kt2(i).ne.ktn2(i))) then
      write(*,*) 'k discrepancy P2 - target has moved out of shutter between nods',i
!       stop
      endif
      if ((it2(i).ne.itn1(i)).or.(it2(i).ne.itn2(i))) then
      write(*,*) 'i discrepancy P2 - target has moved out of shutter between nods',i
!       stop
      endif
      if ((jt2(i).ne.jtn1(i)-1).or.(jt2(i).ne.jtn2(i)+1)) then
      write(*,*) 'j discrepancy P2 - target has moved out of shutter between nods',i
!       stop
      endif
      enddo

      call draw_collapsed_shutter(outdir_2,nm2,rx2,ry2,rxn1,ryn1,rxn2,ryn2,raccx,raccy,cl2,id2)


      ra_p0=ra_p2
      dec_p0=dec_p2
      pa_ap_p0=pa_ap_p2

      call make_slitmask(outdir_2,nm2,kt2,it2,jt2,nsm2,ksm2,ism2,jsm2,ra_p0,dec_p0,pa_ap_p0,skypad)

      write(*,'(a,i0,a,i0,a)') ' ',nsm2,' Sky Background slitlets added to the IPA Pointing ',np2,' MSA Mask'

! append sky slitlets to target slit list

      if (skypad) then
 
      allocate (id0_tmp(n_max_m),id_tmp(n_max_m),cl_tmp(n_max_m),kt_tmp(n_max_m),it_tmp(n_max_m),jt_tmp(n_max_m))
 
      do i=1,nm2
       kt_tmp(i)=kt2(i)
       it_tmp(i)=it2(i)
       jt_tmp(i)=jt2(i)
       cl_tmp(i)=cl2(i)
       id_tmp(i)=id2(i)
       id0_tmp(i)=id02(i)
      enddo

      deallocate (id02,id2,cl2,kt2,it2,jt2)
      allocate (id02(nm2+nsm2),id2(nm2+nsm2),cl2(nm2+nsm2),kt2(nm2+nsm2),it2(nm2+nsm2),jt2(nm2+nsm2))
 
      do i=1,nm2
       kt2(i)=kt_tmp(i)
       it2(i)=it_tmp(i)
       jt2(i)=jt_tmp(i)
       cl2(i)=cl_tmp(i)
       id2(i)=id_tmp(i)
      id02(i)=id0_tmp(i)
      enddo

      do i=1,nsm2
       kt2(nm2+i)=ksm2(i)
       it2(nm2+i)=ism2(i)
       jt2(nm2+i)=jsm2(i)
       cl2(nm2+i)=0
       id2(nm2+i)=0
       id02(nm2+i)=0
      enddo

      deallocate (id0_tmp,id_tmp,cl_tmp,kt_tmp,it_tmp,jt_tmp)

      endif !skypad


      call draw_msa_on_sky(ra_p2,dec_p2,pa_ap_p2,nm2,id2,ra_t2,dec_t2,cl2,outdir_2)


      call write_summary(outdir_2,gname,np2,ra_p2,dec_p2,pa_ap_p2,nm2,id02,id2,cl2,kt2,it2,jt2,rx2,ry2,nsm2,hist_m2,hist_k2,hist_cat,cran1,cdecn1,ktn1,itn1,jtn1,rxn1,ryn1,cran2,cdecn2,ktn2,itn2,jtn2,rxn2,ryn2)

      endif! n_dither>1


!--------------

      if (n_dither.gt.2) then

!         Recover pointing 3 ra,dec,pa_ap

      do i=1,n_p
      if (no_p(i).eq.np3) then
      ra_p3=ra_p(i)
      dec_p3=dec_p(i)
      pa_ap_p3=pa_ap_p(i)
      goto 333
      endif
      enddo
333   continue

!     write(*,*)
!     write(*,*) 'Dither pointing 3:',np3
!     write(*,*) 'Coordinates: ',ra_p3,dec_p3,pa_ap_p3
!     write(*,*) 'Number of targets: ',nm3

! write out pointing 3 observed target catalog to file

      filename=trim(outdir_3)//'/observed_targets.cat'
      open(10,file=filename,status='replace',access='sequential', form='formatted')
      write(10,101) '# No No_sub No_cat  Pr    RA[deg]     Dec[deg]'
      do im3=1,nm3
      write(10,100) im3,id03(im3),id3(im3),cl3(im3),ra_t3(im3),dec_t3(im3)
      enddo
      close(10)

! redo shutter assignments for pointing 3 nods

      call calc_nods(ra_p3,dec_p3,pa_ap_p3,nm3,ra_t3,dec_t3,cran1,cdecn1,ktn1,itn1,jtn1,rxn1,ryn1,cran2,cdecn2,ktn2,itn2,jtn2,rxn2,ryn2)

! sanity check

      do i=1,nm3
      if ((kt3(i).ne.ktn1(i)).or.(kt3(i).ne.ktn2(i))) then
      write(*,*) 'k discrepancy P3 - target has moved out of shutter between nods',i
!       stop
      endif
      if ((it3(i).ne.itn1(i)).or.(it3(i).ne.itn2(i))) then
      write(*,*) 'i discrepancy P3 - target has moved out of shutter between nods',i
!       stop
      endif
      if ((jt3(i).ne.jtn1(i)-1).or.(jt3(i).ne.jtn2(i)+1)) then
      write(*,*) 'j discrepancy P3 - target has moved out of shutter between nods',i
!       stop
      endif
      enddo

      call draw_collapsed_shutter(outdir_3,nm3,rx3,ry3,rxn1,ryn1,rxn2,ryn2,raccx,raccy,cl3,id3)

      ra_p0=ra_p3
      dec_p0=dec_p3
      pa_ap_p0=pa_ap_p3

      call make_slitmask(outdir_3,nm3,kt3,it3,jt3,nsm3,ksm3,ism3,jsm3,ra_p0,dec_p0,pa_ap_p0,skypad)

      write(*,'(a,i0,a,i0,a)') ' ',nsm3,' Sky Background slitlets added to the IPA Pointing ',np3,' MSA Mask'


! append sky slitlets to target slit list

      if (skypad) then
 
      allocate (id0_tmp(n_max_m),id_tmp(n_max_m),cl_tmp(n_max_m),kt_tmp(n_max_m),it_tmp(n_max_m),jt_tmp(n_max_m))
 
      do i=1,nm3
       kt_tmp(i)=kt3(i)
       it_tmp(i)=it3(i)
       jt_tmp(i)=jt3(i)
       cl_tmp(i)=cl3(i)
       id_tmp(i)=id3(i)
       id0_tmp(i)=id03(i)
      enddo

      deallocate (id03,id3,cl3,kt3,it3,jt3)
      allocate (id03(nm3+nsm3),id3(nm3+nsm3),cl3(nm3+nsm3),kt3(nm3+nsm3),it3(nm3+nsm3),jt3(nm3+nsm3))
 
      do i=1,nm3
       kt3(i)=kt_tmp(i)
       it3(i)=it_tmp(i)
       jt3(i)=jt_tmp(i)
       cl3(i)=cl_tmp(i)
       id3(i)=id_tmp(i)
      id03(i)=id0_tmp(i)
      enddo

      do i=1,nsm3
       kt3(nm3+i)=ksm3(i)
       it3(nm3+i)=ism3(i)
       jt3(nm3+i)=jsm3(i)
       cl3(nm3+i)=0
       id3(nm3+i)=0
       id03(nm3+i)=0
      enddo

      deallocate (id0_tmp,id_tmp,cl_tmp,kt_tmp,it_tmp,jt_tmp)

      endif !skypad

      call draw_msa_on_sky(ra_p3,dec_p3,pa_ap_p3,nm3,id3,ra_t3,dec_t3,cl3,outdir_3)

      call write_summary(outdir_3,gname,np3,ra_p3,dec_p3,pa_ap_p3,nm3,id03,id3,cl3,kt3,it3,jt3,rx3,ry3,nsm3,hist_m3,hist_k3,hist_cat,cran1,cdecn1,ktn1,itn1,jtn1,rxn1,ryn1,cran2,cdecn2,ktn2,itn2,jtn2,rxn2,ryn2)

 
      endif !n_dither>2



      write(*,*)
      write(*,*)'Analyzing combined observed target sample'


      if (n_dither.eq.1) then
      write(ipap,"(i0)") np1
      outfile='trial_'//rn(1:2)//'/m_pick_output/pointing_'//trim(adjustl(ipap))//'/target_list.txt'
      endif
      
      

      if (n_dither.eq.2) then
      write(jpap,"(i0)") npaircount
      outfile='trial_'//rn(1:2)//'/m_pick_output/pair_'//trim(adjustl(jpap))//'/combined_target_list.txt'
      endif


      if (n_dither.eq.3) then
      write(jpap,"(i0)") ntriplecount
      outfile='trial_'//rn(1:2)//'/m_pick_output/triple_'//trim(adjustl(jpap))//'/combined_target_list.txt'
      endif


      call common_target_list(nm1,id01,id1,cl1,nm2,id02,id2,cl2,nm3,id03,id3,cl3,np1,np2,np3,outfile)


! append to configuration file

      open(11,file=confname,status='old',access='sequential', form='formatted')

      call find_marker(11,'#MP#')

      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- SECTION AUTOMATICALLY APPENDED BY THE M_PICK MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a)') '#'



      if (n_dither.eq.1) then
      write(11,'(a)') '# Selected single optimal pointing: '
      write(11,'(a)') '#'
      write(11,1120) '# IPA Pointing No:      ',np1
      write(11,1130) '# Coordinates and actual MSA Roll: ',ra_p1,dec_p1,pa_ap_p1
      write(11,1140) '# Number of targets:      ',nm1
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# Target breakdown by Priority Class:'
      write(11,'(a)') '#'
      write(11,'(a)') '#    Pri   Number'
      do i=1,nclass
      write(11,1160) '#',i,hist_com(i)
1160  format(a,3x,i3,2x,i5)
      enddo !i
      write(11,'(a)') '#'
      write(11,'(a,i0,a,i0,a)') '# ',nsm1,' Sky Background slitlets added to the IPA Pointing ',np1,' MSA Mask'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# The full output from the m_pick module can be found in the directory:'
      write(11,'(a)') '#'
      write(11,'(a)') '#    ./'//trim(adjustl(outdir_1))//'/'
      write(11,'(a)') '#'
      write(11,'(a)') '#  and in the file:'
      write(11,'(a)') '#'
      write(11,'(a)') '#    ./trial_'//rn(1:2)//'/m_pick_output/pointing_'//trim(adjustl(ipap))//'/target_list.txt'
      endif !n_dither=1


1120  format(a,I4)
1130  format(a,f13.7,3x,f11.7,3x,f10.6)
1140  format(a,I4)
1100  format(a)

      if (n_dither.eq.2) then
      write(11,'(a,i0)') '# Selected optimal dithered pointing pair: ',mtarget
      write(11,'(a)') '#'
      write(11,1120) '# IPA Pointing No:      ',np1
      write(11,1130) '# Coordinates and Roll: ',ra_p1,dec_p1,pa_ap_p1
      write(11,1140) '# Number of targets:      ',nm1
      write(11,'(a)') '#'
      write(11,1120) '# IPA Pointing No:      ',np2
      write(11,1130) '# Coordinates and Roll: ',ra_p2,dec_p2,pa_ap_p2
      write(11,1140) '# Number of targets:      ',nm2
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,1100) '# Combined Sample:'
      write(11,'(a)') '#'
      write(11,1180) '# Total number of spectra: ',n_spec
      write(11,1180) '# Number of unique targets:',n_com
      write(11,'(a)') '#'
      write(11,1100) '# Unique target breakdown by Priority Class:'
      write(11,'(a)') '#'
      write(11,1100) '#    Pri  Number   AvExp'
       do i=1,nclass
        write(11,1170) '#',i,hist_com(i),avexp(i)
       enddo !i
      write(11,'(a)') '#'

      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a,i0,a,i0,a)') '# ',nsm1,' Sky Background slitlets added to the IPA Pointing ',np1,' MSA Mask'
      write(11,'(a,i0,a,i0,a)') '# ',nsm2,' Sky Background slitlets added to the IPA Pointing ',np2,' MSA Mask'
      write(11,'(a)') '#'
      write(11,'(a)') '#'


      write(11,1100) '# The full output from the m_pick module can be found in the directories:'
      write(11,'(a)') '#'
      write(11,'(a)') '#    ./'//trim(adjustl(outdir_1))//'/'
      write(11,'(a)') '#    ./'//trim(adjustl(outdir_2))//'/'
      write(11,'(a)') '#'
      write(11,1100) '#  and in the file:'
      write(11,'(a)') '#'
      write(11,1100) '#    ./trial_'//rn(1:2)//'/m_pick_output/pair_'//trim(adjustl(jpap))//'/combined_target_list.txt'
      endif !n_dither=2



      if (n_dither.eq.3) then
      write(11,'(a,i0)') '# Selected optimal dithered pointing triple: ',mtarget
      write(11,'(a)') '#'
      write(11,1100) '# Dither Pointing 1:'
      write(11,1120) '# IPA Pointing No:      ',np1
      write(11,1130) '# Coordinates and Roll: ',ra_p1,dec_p1,pa_ap_p1
      write(11,1140) '# Number of targets:      ',nm1
      write(11,'(a)') '#'
      write(11,1100) '# Dither Pointing 2:'
      write(11,1120) '# IPA Pointing No:      ',np2
      write(11,1130) '# Coordinates and Roll: ',ra_p2,dec_p2,pa_ap_p2
      write(11,1140) '# Number of targets:      ',nm2
      write(11,'(a)') '#'
      write(11,1100) '# Dither Pointing 3:'
      write(11,1120) '# IPA Pointing No:      ',np3
      write(11,1130) '# Coordinates and Roll: ',ra_p3,dec_p3,pa_ap_p3
      write(11,1140) '# Number of targets:      ',nm3
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,1100) '# Combined Sample:'
      write(11,'(a)') '#'
      write(11,1180) '# Total number of spectra: ',n_spec
      write(11,1180) '# Number of unique targets:',n_com
      write(11,'(a)') '#'
      write(11,1100) '# Unique target breakdown by Priority Class:'
      write(11,'(a)') '#'
      write(11,1100) '#   Class  Number  AvExp'
       do i=1,nclass
        write(11,1170) '#',i,hist_com(i),avexp(i)
       enddo !i
      write(11,'(a)') '#'
1170  format(a,3x,i3,2x,i5,5x,f4.2)
1180  format(a,i5)


      write(11,'(a)') '#'
      write(11,'(a,i0,a,i0,a)') '# ',nsm1,' Sky Background slitlets added to the IPA Pointing ',np1,' MSA Mask'
      write(11,'(a,i0,a,i0,a)') '# ',nsm2,' Sky Background slitlets added to the IPA Pointing ',np2,' MSA Mask'
      write(11,'(a,i0,a,i0,a)') '# ',nsm3,' Sky Background slitlets added to the IPA Pointing ',np3,' MSA Mask'
      write(11,'(a)') '#'
      write(11,'(a)') '#'


      write(11,1100) '# The full output from the m_pick module can be found in the directories:'
      write(11,'(a)') '#'
      write(11,'(a)') '#    ./'//trim(adjustl(outdir_1))//'/'
      write(11,'(a)') '#    ./'//trim(adjustl(outdir_2))//'/'
      write(11,'(a)') '#    ./'//trim(adjustl(outdir_3))//'/'
      write(11,'(a)') '#'
      write(11,1100) '#  and in the file:'
      write(11,'(a)') '#'
      write(11,1100) '#    ./trial_'//rn(1:2)//'/m_pick_output/triple_'//trim(adjustl(jpap))//'/combined_target_list.txt'
      endif !n_dither=3


      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- END OF SECTION APPENDED BY THE M_PICK MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a)') '#  -- Changing the preferred choice among the optimal pointings requires re-running the m_pick and subsequent modules ---'
      write(11,'(a)') '#'
      write(11,'(a)') '#MC# -- Marker do not delete or move this line'

       close(11)

! plot spectra for all dispersers

      write(*,*)'Generating spectral trace plots'
      write(*,*)


      call plot_spectra(outdir_1,7,nm1,id1,cl1,kt1,it1,jt1,nsm1)
      call plot_spectra(outdir_1,1,nm1,id1,cl1,kt1,it1,jt1,nsm1)
      call plot_spectra(outdir_1,2,nm1,id1,cl1,kt1,it1,jt1,nsm1)
      call plot_spectra(outdir_1,3,nm1,id1,cl1,kt1,it1,jt1,nsm1)
      call plot_spectra(outdir_1,4,nm1,id1,cl1,kt1,it1,jt1,nsm1)
      call plot_spectra(outdir_1,5,nm1,id1,cl1,kt1,it1,jt1,nsm1)
      call plot_spectra(outdir_1,6,nm1,id1,cl1,kt1,it1,jt1,nsm1)

      if (n_dither.gt.1) then
      call plot_spectra(outdir_2,7,nm2,id2,cl2,kt2,it2,jt2,nsm2)
      call plot_spectra(outdir_2,1,nm2,id2,cl2,kt2,it2,jt2,nsm2)
      call plot_spectra(outdir_2,2,nm2,id2,cl2,kt2,it2,jt2,nsm2)
      call plot_spectra(outdir_2,3,nm2,id2,cl2,kt2,it2,jt2,nsm2)
      call plot_spectra(outdir_2,4,nm2,id2,cl2,kt2,it2,jt2,nsm2)
      call plot_spectra(outdir_2,5,nm2,id2,cl2,kt2,it2,jt2,nsm2)
      call plot_spectra(outdir_2,6,nm2,id2,cl2,kt2,it2,jt2,nsm2)
      endif

      if (n_dither.gt.2) then
      call plot_spectra(outdir_3,7,nm3,id3,cl3,kt3,it3,jt3,nsm3)
      call plot_spectra(outdir_3,1,nm3,id3,cl3,kt3,it3,jt3,nsm3)
      call plot_spectra(outdir_3,2,nm3,id3,cl3,kt3,it3,jt3,nsm3)
      call plot_spectra(outdir_3,3,nm3,id3,cl3,kt3,it3,jt3,nsm3)
      call plot_spectra(outdir_3,4,nm3,id3,cl3,kt3,it3,jt3,nsm3)
      call plot_spectra(outdir_3,5,nm3,id3,cl3,kt3,it3,jt3,nsm3)
      call plot_spectra(outdir_3,6,nm3,id3,cl3,kt3,it3,jt3,nsm3)
      endif


      write(*,*)
      call cpu_time(t_stop)
      write(*,'(a)')      '-------------------------------------------------'
      write(*,'(a,f7.1)') ' CPU time used by m_pick module [s]:',t_stop-t_start
      write(*,'(a)')      '-------------------------------------------------'
      write(*,*)

         
      end program m_pick
      
!-------------------------------------------------------
 
 
include './reference_files/front_transforms_updt.f90'

include './reference_files/back_transforms_updt.f90'

include './reference_files/shutter_routines_new.f90'

include './reference_files/misc_routines_new.f90'


!------------------------------------------------------------

      subroutine file_entry_m_pick()
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

      call skip_hashes(9)
      read(9,*) sthresh


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
       read(9,*) cra,cdec,cpa_ap

! work out pointing reference point coordinates in on-sky systems

      call msa_to_sky_bd(xm_ref,ym_ref,ax_refv2,ay_refv3) !v2,V3


      call skip_hashes(9)
      read(9,*) !szone_x,szone_y

      call skip_hashes(9)
      read(9,*) angle_to_target
      call dva_calc(angle_to_target,dva_mag)

      call skip_hashes(9)
      read(9,*) !max_c_scores1
 
      call skip_hashes(9)
      read(9,*) !max_c_scores2

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

      end subroutine file_entry_m_pick

!------------------------------------------------------------

      subroutine calc_nods(ccra,ccdec,ccpa_ap,nt,ra_t,dec_t,cran1,cdecn1,ktn1,itn1,jtn1,rxn1,ryn1,cran2,cdecn2,ktn2,itn2,jtn2,rxn2,ryn2)
      use derived_parameters, only : ax_refv2,ay_refv3,xm_ref,ym_ref,phi     
      use config_parameters, only : dva_mag
      implicit none
      integer nt
      real*8 ccra,ccdec
      real ccpa_ap,ccpa_v3
      real*8 ra_t(nt),dec_t(nt)
      real xx,yy,xm,ym
      real xt(nt),yt(nt)
      real xtv2(nt),ytv3(nt)
      real theta_x,theta_y
      real*8 cran1,cdecn1,cran2,cdecn2
      real av_ypitch_msa
      real ax_refra_p,ay_refdec_p
      integer i
      integer ktn1(nt),itn1(nt),jtn1(nt)
      real rxn1(nt),ryn1(nt)
      integer ktn2(nt),itn2(nt),jtn2(nt)
      real rxn2(nt),ryn2(nt)
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar

       ccpa_v3=ccpa_ap-(180.-phi)
      if (ccpa_v3.lt.0.) ccpa_v3=ccpa_v3+360.
      if (ccpa_v3.gt.360.) ccpa_v3=ccpa_v3-360.

! work out pointing reference point in RA,Dec for roll in question

      call roll_v3(ccpa_v3,ax_refv2,ay_refv3,ax_refra_p,ay_refdec_p) !RA,Dec

! Initialize mode-specific subroutines

      call load_msa_dimensions_bd
      call load_fpa_dimensions_bd

       av_ypitch_msa=(ypitch(1)+ypitch(2)+ypitch(3)+ypitch(4))/4.


! calculate nodded 3-slitlet pointings
!

      call msa_to_sky_bd(xm_ref,ym_ref+av_ypitch_msa,xx,yy)
      call roll_v3(ccpa_v3,xx,yy,theta_x,theta_y)

      call coordtang_single(ccra,ccdec,theta_x-ax_refra_p,theta_y-ay_refdec_p,cran2,cdecn2)

      call msa_to_sky_bd(xm_ref,ym_ref-av_ypitch_msa,xx,yy)
      call roll_v3(ccpa_v3,xx,yy,theta_x,theta_y)

      call coordtang_single(ccra,ccdec,theta_x-ax_refra_p,theta_y-ay_refdec_p,cran1,cdecn1)

!     write(*,*) 'RA, Dec of Nodded Pointings:'
!     write(*,99)  ' Nod 1:',cran1,cdecn1 !-1 shutter down
!     write(*,99)  ' Nod 2:',cran2,cdecn2 !+1 shutter up
!99    format(a7,2x,f13.7,3x,f11.7)
!     write(*,*)
!     write(*,*)


! Calculate shutter MSA locations for top nod

      call tangcoord(cran1,cdecn1,nt,ra_t,dec_t,xt,yt)

      do i=1,nt

      xt(i)=(xt(i)*dva_mag+ax_refra_p)
      yt(i)=(yt(i)*dva_mag+ay_refdec_p)


! Rotate to V2,V3 system orientation

      call roll_v3(-ccpa_v3,xt(i),yt(i),xtv2(i),ytv3(i))

! May now be projected to MSA.

       call sky_to_msa_bd(xtv2(i),ytv3(i),xm,ym)

! Identify which quadrant and shutter targets fall in (kt(i)=0 ~ misses MSA entirely)

      call msa_to_shutter_bd(xm,ym,ktn1(i),itn1(i),jtn1(i),rxn1(i),ryn1(i))

      enddo


! Calculate shutter MSA locations for bottom nod

      call tangcoord(cran2,cdecn2,nt,ra_t,dec_t,xt,yt)

      do i=1,nt

      xt(i)=(xt(i)*dva_mag+ax_refra_p)
      yt(i)=(yt(i)*dva_mag+ay_refdec_p)


! Rotate to V2,V3 system orientation

      call roll_v3(-ccpa_v3,xt(i),yt(i),xtv2(i),ytv3(i))

! May now be projected to MSA.

       call sky_to_msa_bd(xtv2(i),ytv3(i),xm,ym)

! Identify which quadrant and shutter targets fall in (kt(i)=0 ~ misses MSA entirely)

      call msa_to_shutter_bd(xm,ym,ktn2(i),itn2(i),jtn2(i),rxn2(i),ryn2(i))

      enddo


      return

      end subroutine calc_nods

!----------------------------------------------------------------------------------------

      subroutine draw_collapsed_shutter(outdir,n,rx0,ry0,rx1,ry1,rx2,ry2,raccx,raccy,cl,id)

      use derived_parameters, only : av_xpitch,av_ypitch,av_xopen,av_yopen
      implicit none
      integer n,i,j
      integer cl(n),id(n)
      real rx0(n),ry0(n),rx1(n),ry1(n),rx2(n),ry2(n)
      real raccx,raccy
      real av_ypitch_as,av_xpitch_as,av_xopen_as,av_yopen_as
      real x1(5),y1(5),xp1,xp2,yp1,yp2,xx,yy
      character*13 tnum
      integer ntnum
      character*70 filename
      character*70 outdir

      filename=trim(outdir)//'/collapsed_shutter.ps/cps'

        av_xpitch_as=av_xpitch/1000.
        av_ypitch_as=av_ypitch/1000.
        av_xopen_as=av_xopen/1000.
        av_yopen_as=av_yopen/1000.


        CALL PGBEGIN(0,filename,1,1)
!       CALL PGBEGIN(0,'?',1,1)
        CALL PGSLS(1)
        CALL PGSLW(1)
        CALL PGSCF(1)
        CALL PGSCH(1.)

        xp1=-av_xpitch_as*0.7
        xp2=+av_xpitch_as*0.7
        yp1=av_ypitch_as*0.7
        yp2=-av_ypitch_as*0.7

        CALL PGVPORT(.0,1.,.0,1.)
        CALL PGWNAD(xp2,xp1,yp1,yp2)
        CALL PGSLW(1)

! draw shutter outline

      do i=-1,1
      do j=-1,1


      x1(1)=-av_xpitch_as/2.+float(i)*av_xpitch_as
      y1(1)=-av_ypitch_as/2.+float(j)*av_ypitch_as
      x1(2)=av_xpitch_as/2.+float(i)*av_xpitch_as
      y1(2)=-av_ypitch_as/2.+float(j)*av_ypitch_as
      x1(3)=av_xpitch_as/2.+float(i)*av_xpitch_as
      y1(3)=av_ypitch_as/2.+float(j)*av_ypitch_as
      x1(4)=-av_xpitch_as/2.+float(i)*av_xpitch_as
      y1(4)=av_ypitch_as/2.+float(j)*av_ypitch_as
      x1(5)=x1(1)
      y1(5)=y1(1)

      call pgsci(15)
      call pgpoly(5,x1,y1)

      call pgsci(1)
      call pgline(5,x1,y1)

      x1(1)=-av_xopen_as/2.+float(i)*av_xpitch_as
      y1(1)=-av_yopen_as/2.+float(j)*av_ypitch_as
      x1(2)=av_xopen_as/2.+float(i)*av_xpitch_as
      y1(2)=-av_yopen_as/2.+float(j)*av_ypitch_as
      x1(3)=av_xopen_as/2.+float(i)*av_xpitch_as
      y1(3)=av_yopen_as/2.+float(j)*av_ypitch_as
      x1(4)=-av_xopen_as/2.+float(i)*av_xpitch_as
      y1(4)=av_yopen_as/2.+float(j)*av_ypitch_as
      x1(5)=x1(1)
      y1(5)=y1(1)

      call pgsci(7)
      call pgpoly(5,x1,y1)

      call pgsci(15)
      call pgline(5,x1,y1)
      call pgsci(1)

      enddo
      enddo


! draw acceptance zone


      x1(1)=-raccx*av_xpitch_as
      y1(1)=+raccy*av_ypitch_as
      x1(2)=+raccx*av_xpitch_as
      y1(2)=+raccy*av_ypitch_as
      x1(3)=+raccx*av_xpitch_as
      y1(3)=-raccy*av_ypitch_as
      x1(4)=-raccx*av_xpitch_as
      y1(4)=-raccy*av_ypitch_as
      x1(5)=x1(1)
      y1(5)=y1(1)
      call pgsci(2)
      call pgline(5,x1,y1)


      call pgscr(20,0./255.,139./255.,0./255.) !Green 4

      call pgscr(21,192./255.,255./255.,62./255.) !Olivedrab 1

      call pgscr(22,154./255.,255./255.,154./255.) !Palegreen 1


      do i=n,1,-1
      xx=rx1(i)*av_xpitch_as
      yy=ry1(i)*av_ypitch_as
      call pgsci(21)
      call pgpoint(1,xx,yy,-16)
      enddo

      do i=n,1,-1
      xx=rx2(i)*av_xpitch_as
      yy=ry2(i)*av_ypitch_as
      call pgsci(22)
      call pgpoint(1,xx,yy,-16)
      enddo

      do i=n,1,-1
      xx=rx0(i)*av_xpitch_as
      yy=ry0(i)*av_ypitch_as
      call pgsci(20)
      if (cl(i).eq.1) call pgsci(2)
      call pgpoint(1,xx,yy,-16)
      enddo

      call pgsci(1)
      call pgsch(0.3)

      do i=n,1,-1
      xx=rx0(i)*av_xpitch_as
      yy=ry0(i)*av_ypitch_as
      call pgnumb(id(i),0,1,tnum,ntnum)
      call pgptext(xx,yy+0.0015,0.,0.5,tnum(1:ntnum))
      enddo


      call pgend

      return
      end subroutine draw_collapsed_shutter

!--------------------------------------------------------------------

      subroutine make_slitmask(outdir,nt,kt,it,jt,nsm,ksm,ism,jsm,ra_p0,dec_p0,pa_ap_p0,skypad)
      use path_to_ref_files
      use config_parameters, only : n_disp,sthresh
      use derived_parameters, only : avxgap
      implicit none
      integer nt
      integer kt(nt),it(nt),jt(nt)
      integer nsm
      integer, parameter :: max_spare_m=250 !Maximum number of non-overlapping unu
      integer ksm(max_spare_m),ism(max_spare_m),jsm(max_spare_m)
      real*8 ra_p0,dec_p0
      real pa_ap_p0
      logical skypad
      integer im,ii,jj,j1,j2,dj,k,i,j,k2,ntout
      real xx,yy
      real i1,i2
      character*70 filename
      character*70 outdir
! - slitmap
      integer slitmap(4,365,171)
      common/slit_map/slitmap
! - shutter values
      real shval(4,365,171),msa_gap
      common/shutter_values/shval
 ! - spectral parameters
!      integer n_disp
      real lambda1,lambda2
      character*6 gname
      common/disp/lambda1,lambda2,gname
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
! - prism seperations
      integer prism_sep_p(4,365,171),prism_sep_m(4,365,171)
      common/prism_length/prism_sep_p,prism_sep_m
      real sep_m,sep_p

! Compiler warning avoidance

      k2=0

      call load_msa_dimensions_bd
      call load_fpa_dimensions_bd
      call msa_to_fpa_disp(n_disp,1.,xx,yy,xx,yy) !initialize disperser

      call read_shutter_values
      if (n_disp.eq.7) call read_prism_separations
      call read_base_msamap
      call make_slitlet_map

      filename=trim(outdir)//'/slitlet_usage_stats.txt'
      open(10,file=filename,status='replace',access='sequential', form='formatted')


      write(10,*)
      write(10,*) 'Slitmap statistics at start:'
      write(10,*)

      call slitmap_stats(10)

      call initialize_slitmap(n_disp,sthresh,1,1)

      filename=trim(outdir)//'/slitmap_start.ps/cps'
      call plot_slitmap(filename)

      write(10,*)
      write(10,*) 'Slitmap after FO overlap and PRISM truncation purge:'
      write(10,*)

      call slitmap_stats(10)

! March through accepted targets

      ntout=0
      dj=int(2.*sthresh)


      do im=1,nt

      k=kt(im)
      i=it(im)
      j=jt(im)

! Mark commanded open target slitlets in slitmap

      if (slitmap(k,i,j).gt.1) slitmap(k,i,j)=-1

      if (j+1.le.171) then
      if (slitmap(k,i,j+1).gt.1) slitmap(k,i,j+1)=-1
      endif

      if (j-1.ge.1) then
      if (slitmap(k,i,j-1).gt.1) slitmap(k,i,j-1)=-1
      endif

! Map out shutters excluded due to target overlaps in slitmap

      j1=max(j-dj,1)
      j2=min(j+dj,171)
      if (k.eq.1) k2=3
      if (k.eq.3) k2=1
      if (k.eq.2) k2=4
      if (k.eq.4) k2=2

      do ii=1,365
      do jj=j1,j2

      if (abs(shval(k,i,j)-shval(k,ii,jj)).le.sthresh) then

       if (n_disp.eq.7) then
       i1=float(i)+msa_gap(k)*(avxgap+365.)
       i2=float(ii)+msa_gap(k)*(avxgap+365.)
         sep_p=float(prism_sep_p(k,i,j))
         sep_m=float(prism_sep_m(k,i,j))
       if (((i2-i1).le.sep_m).or.((i2-i1).ge.sep_p)) goto 323
       endif

        if ((slitmap(k,ii,jj).ge.2).and.(slitmap(k,ii,jj).le.4)) then
        slitmap(k,ii,jj)=10
        ntout=ntout+1
        endif

323   continue
      endif


      if (abs(shval(k,i,j)-shval(k2,ii,jj)).le.sthresh) then

       if (n_disp.eq.7) then
       i1=float(i)+msa_gap(k)*(avxgap+365.)
       i2=float(ii)+msa_gap(k2)*(avxgap+365.)
         sep_p=float(prism_sep_p(k,i,j))
         sep_m=float(prism_sep_m(k,i,j))
       if (((i2-i1).le.sep_m).or.((i2-i1).ge.sep_p)) goto 324
       endif


        if ((slitmap(k2,ii,jj).ge.2).and.(slitmap(k2,ii,jj).le.4)) then
        slitmap(k2,ii,jj)=10
        ntout=ntout+1
        endif

324   continue
      endif

      enddo !jj
      enddo !ii

      enddo !im

!     write(*,*) 'Vacant slitlets claimed:',ntout
!     write(*,*)

      write(10,*)
      write(10,*) 'Slitmap after target placement:'
      write(10,*)

      call slitmap_stats(10)

      nsm=0

      if (skypad) then

      call use_spare_shutters(ra_p0,dec_p0,pa_ap_p0,nsm,ksm,ism,jsm)

      write(10,*)
      write(10,*) 'Slitmap after placement of spare sky slitlets:'
      write(10,*)

      call slitmap_stats(10)

      endif !skypad

      write(10,*)

      close(10)

      filename=trim(outdir)//'/slitmap_end.ps/cps'
      call plot_slitmap(filename)

      filename=trim(outdir)//'/shutter_mask.csv'
      call make_csv_file(filename)


      return
      end subroutine make_slitmask

!------------------------------------------------------------------------------------

      subroutine draw_msa_on_sky(ccra,ccdec,ccpa_ap,nt,id,rat,dect,pri,outdir)
      use config_parameters, only : dva_mag,cra,cdec
      use target_cat, only : nclass
      use derived_parameters, only : ax_refv2,ay_refv3,phi     
      implicit none
      real*8 rat(nt),dect(nt)
      integer pri(nt),id(nt)
      integer nt
      real ccpa_ap,ccpa_v3
      real*8 ccra,ccdec
      real ax_refra_p,ay_refdec_p
      real ax,ay,ax_0,ay_0
      integer i,k
      real xx,yy,xp1,xp2,yp1,yp2
      real x1(5),y1(5),x2(5),y2(5),x3(5),y3(5),x4(5),y4(5)
      real x1r(5),y1r(5),x2r(5),y2r(5),x3r(5),y3r(5),x4r(5),y4r(5)
      real r0,twopi,ang,xxs(50),yys(50)
      character*13 tnum
      integer ntnum
      character*70 outdir
      character*70 filename
      character*2 rn
      common/run/rn
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
      
      filename=trim(outdir)//'/fov_on_sky.ps/cps'

       ccpa_v3=ccpa_ap-(180.-phi)
      if (ccpa_v3.lt.0.) ccpa_v3=ccpa_v3+360.
      if (ccpa_v3.gt.360.) ccpa_v3=ccpa_v3-360.


! work out pointing reference point for present roll angle
 
      call roll_v3(ccpa_v3,ax_refv2,ay_refv3,ax_refra_p,ay_refdec_p) !RA,Dec

!       print *
!       print *,'ax_refv2,ay_refv3:',ax_refv2,ay_refv3
!  
!       print *,'ax_refra_p,ay_refdec_p:',ax_refra_p,ay_refdec_p
!       print *  



! Determine MSA Quadrant Outlines

! Quadrant 1

      x1(1)=0.
      y1(1)=0.
      x1(2)=x1(1)
      y1(2)=y1(1)+float(my)*ypitch(1)
      x1(3)=x1(2)+float(mx)*xpitch(1)
      y1(3)=y1(2)
      x1(4)=x1(3)
      y1(4)=y1(1)
      call rottran(1,4,x1,y1)

      do i=1,4
       call msa_to_sky_bd(x1(i),y1(i),xx,yy)
       x1(i)=xx
       y1(i)=yy
      enddo

      x1(5)=x1(1)
      y1(5)=y1(1)


! Quadrant 2

      x2(1)=0.
      y2(1)=0.
      x2(2)=x2(1)
      y2(2)=y2(1)+float(my)*ypitch(2)
      x2(3)=x2(2)+float(mx)*xpitch(2)
      y2(3)=y2(2)
      x2(4)=x2(3)
      y2(4)=y2(1)
      call rottran(2,4,x2,y2)

      do i=1,4
       call msa_to_sky_bd(x2(i),y2(i),xx,yy)
       x2(i)=xx
       y2(i)=yy
      enddo

      x2(5)=x2(1)
      y2(5)=y2(1)

! Quadrant 3

      x3(1)=0.
      y3(1)=0.
      x3(2)=x3(1)
      y3(2)=y3(1)+float(my)*ypitch(3)
      x3(3)=x3(2)+float(mx)*xpitch(3)
      y3(3)=y3(2)
      x3(4)=x3(3)
      y3(4)=y3(1)
      call rottran(3,4,x3,y3)

      do i=1,4
       call msa_to_sky_bd(x3(i),y3(i),xx,yy)
       x3(i)=xx
       y3(i)=yy
      enddo

      x3(5)=x3(1)
      y3(5)=y3(1)

! Quadrant 4

      x4(1)=0.
      y4(1)=0.
      x4(2)=x4(1)
      y4(2)=y4(1)+float(my)*ypitch(4)
      x4(3)=x4(2)+float(mx)*xpitch(4)
      y4(3)=y4(2)
      x4(4)=x4(3)
      y4(4)=y4(1)
      call rottran(4,4,x4,y4) !minor quadrant-dependent rotation in xm,y

      do i=1,4
       call msa_to_sky_bd(x4(i),y4(i),xx,yy) !V2,V3
       x4(i)=xx
       y4(i)=yy
      enddo

      x4(5)=x4(1)
      y4(5)=y4(1)


! Rotate MSA quadrant corners to V3 roll angle on sky

        do i=1,5
        call roll_v3(ccpa_v3,x1(i),y1(i),x1r(i),y1r(i)) !ra,dec
        call roll_v3(ccpa_v3,x2(i),y2(i),x2r(i),y2r(i))
        call roll_v3(ccpa_v3,x3(i),y3(i),x3r(i),y3r(i))
        call roll_v3(ccpa_v3,x4(i),y4(i),x4r(i),y4r(i))
        enddo


! Plot MSA FOV and target catalog on sky

        CALL PGBEGIN(0,filename,1,1)
        CALL PGSLS(1)
        CALL PGSLW(1)
        CALL PGSCF(2)
        CALL PGSCH(1.)


        xp1=-55.*2.516 !nominal plate scale in arcsec/mm
        xp2= 55.*2.516
        yp1=-50.*2.516
        yp2= 50.*2.516


        xp1=xp1*1.5+ax_refra_p !
        xp2=xp2*1.5+ax_refra_p
        yp1=yp1*1.5+ay_refdec_p
        yp2=yp2*1.5+ay_refdec_p
        
!         print *,'xp1,xp2,yp1,yp2',xp1,xp2,yp1,yp2
        

        CALL PGVPORT(.0,1.,.0,1.)
        CALL PGWNAD(xp2,xp1,yp1,yp2) !flipped x-direction
        CALL PGSLW(1)



! Draw MSA quadrant outlines

      call pgsch(1.0)
      call pgscf(2)
      call pgslw(1)

      call pgsci(14)
      call pgline(5,x1r,y1r)
      xx=(x1r(1)+x1r(3))/2.
      yy=(y1r(1)+y1r(3))/2.
      call pgsci(15)
      call pgptext(xx,yy,0.,0.5,'Q1')

      call pgsci(14)
      call pgline(5,x2r,y2r)
      xx=(x2r(1)+x2r(3))/2.
      yy=(y2r(1)+y2r(3))/2.
      call pgsci(15)
      call pgptext(xx,yy,0.,0.5,'Q2')

      call pgsci(14)
      call pgline(5,x3r,y3r)
      xx=(x3r(1)+x3r(3))/2.
      yy=(y3r(1)+y3r(3))/2.
      call pgsci(15)
      call pgptext(xx,yy,0.,0.5,'Q3')

      call pgsci(14)
      call pgline(5,x4r,y4r)
      xx=(x4r(1)+x4r(3))/2.
      yy=(y4r(1)+y4r(3))/2.
      call pgsci(15)
      call pgptext(xx,yy,0.,0.5,'Q4')


! draw fixed slits

      call draw_slits_sky_rot(ccpa_v3)

      call pgsci(1)

      call draw_v2v3(ccpa_v3,ax_refra_p,ay_refdec_p)

      call draw_compass(0.,ax_refra_p,ay_refdec_p)

! calculate offset of catalog center

      call tangcoord_single(cra,cdec,ccra,ccdec,ax_0,ay_0)

!     write(*,*) 'ra,dec offset from catalog center ["]:',ax_0,ay_0,sqrt(ax_0*ax_0+ay_0*ay_0)


! draw offset over-booking zone

         r0=180.
         twopi=6.283185
         do i=1,49
         ang=twopi*float(i)/49.
         xxs(i)=-ax_0+r0*cos(ang)+ax_refra_p
         yys(i)=-ay_0+r0*sin(ang)+ay_refdec_p
         enddo
         xxs(50)=xxs(1)
         yys(50)=yys(1)

         call pgslw(1)
         call pgsci(14)
         call pgline(50,xxs,yys)
         call pgslw(1)



      call set_rainbow_lut


      call pgsch(0.4)
      call pgslw(2)

      do i=1,nt
      call tangcoord_single(ccra,ccdec,rat(i),dect(i),ax,ay)
      ax=(ax*dva_mag+ax_refra_p)
      ay=(ay*dva_mag+ay_refdec_p)
      k=pri(i)
      call set_rainbow_color_fine(nclass,k)
      call pgpoint(1,ax,ay,-12)
      call pgsci(1)
      call pgnumb(id(i),0,1,tnum,ntnum) ! label target ID
 !     call pgnumb(i,0,1,tnum,ntnum)            ! label output target num
      call pgptext(ax,ay+2.5,0.,0.5,tnum(1:ntnum))
      enddo


      call pgend

      return

      end subroutine draw_msa_on_sky

!--------------------------------------------

      subroutine draw_compass(roll,ax_refra_p,ay_refdec_p) !x flipped version!
      implicit none
      real roll !deg
      real x0,y0,xx(2),yy(2),len,roll_r
      real ax_refra_p,ay_refdec_p
      len=5.
      x0=60.*2.516+ax_refra_p
      y0=50.*2.516+ay_refdec_p
      roll_r=roll/57.2957795
      xx(1)=x0
      yy(1)=y0
      xx(2)=xx(1)+len*sin(roll_r)
      yy(2)=yy(1)+len*cos(roll_r)
      call pgline(2,xx,yy)
      call pgsch(0.5)
      call pgscf(1)
      xx(2)=xx(1)+(len+4.)*sin(roll_r)
      yy(2)=yy(1)+(len+4.)*cos(roll_r)
      call pgptxt(xx(2),yy(2),roll,0.5,'N')
      xx(2)=xx(1)+len*cos(roll_r)
      yy(2)=yy(1)-len*sin(roll_r)
      call pgline(2,xx,yy)
      xx(2)=xx(1)+(len+4.)*cos(roll_r)
      yy(2)=yy(1)-(len+4.)*sin(roll_r)
      call pgptxt(xx(2),yy(2),roll,0.5,'E')

      return
      end subroutine draw_compass

!---------------------------------------------------------------------------

      subroutine draw_slits_sky_rot(pa_v3)
      real pa_v3
      real xx(60),yy(60),xd,yd,xr(60),yr(60)
      real x00,y00,r0,twopi,ang,w,h
      integer i
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar


! IFU Aperture

         X00=0.04125475*1000.
         Y00=(6.99E-06)*1000.
         r0=1.023
         twopi=6.283185
         do i=1,59
         ang=twopi*float(i)/59.
         xx(i)=X00+r0*cos(ang)
         yy(i)=Y00+r0*sin(ang)
      call msa_to_sky_bd(xx(i),yy(i),xd,yd)
      call roll_v3(pa_v3,xd,yd,xr(i),yr(i))
         enddo
         xr(60)=xr(1)
         yr(60)=yr(1)

         call pgslw(1)
         call pgline(60,xr,yr)
         call pgslw(1)


        xx(2)=X00+0.66
        yy(2)=Y00+0.66
        xx(3)=xx(2)
        yy(3)=yy(2)-0.66*2.
        xx(4)=xx(3)-0.66*2.
        yy(4)=yy(3)
        xx(5)=xx(4)
        yy(5)=yy(4)+0.66*2.
        xx(1)=xx(5)
        yy(1)=yy(5)

        do i=1,5
        call msa_to_sky_bd(xx(i),yy(i),xd,yd)
        call roll_v3(pa_v3,xd,yd,xr(i),yr(i))
        enddo

        call pgline(5,xr,yr)


        yy(1)=yy(2)
        yy(2)=yy(3)


        xx(1)=xx(4)


        do i=1,29
        xx(1)=xx(1)+2.*0.66/30.
        xx(2)=xx(1)
         do j=1,2
         call msa_to_sky_bd(xx(j),yy(j),xd,yd)
         call roll_v3(pa_v3,xd,yd,xr(j),yr(j))
         enddo
        call pgline(2,xr,yr)
        enddo




! SLIT_A_200_1


         X00=0.02697243*1000. !mm
         Y00=-0.002716702*1000.
         w=(8.135E-05)*1000.
         h=0.00127117*1000.

         xx(1)=X00-w/2.
         yy(1)=Y00-h/2.
         xx(2)=xx(1)+w
         yy(2)=yy(1)
         xx(3)=xx(2)
         yy(3)=yy(2)+h
         xx(4)=xx(3)-w
         yy(4)=yy(3)
         xx(5)=xx(1)
         yy(5)=yy(1)

        do i=1,5
        call msa_to_sky_bd(xx(i),yy(i),xd,yd)
        call roll_v3(pa_v3,xd,yd,xr(i),yr(i))
        enddo

         call pgline(5,xr,yr)


! SLIT_A_400

         X00=0.02966518*1000.
         Y00=0.0002851795*1000.
         w=0.00015924*1000.
         h=0.00144982*1000.

         xx(1)=X00-w/2.
         yy(1)=Y00-h/2.
         xx(2)=xx(1)+w
         yy(2)=yy(1)
         xx(3)=xx(2)
         yy(3)=yy(2)+h
         xx(4)=xx(3)-w
         yy(4)=yy(3)
         xx(5)=xx(1)
         yy(5)=yy(1)

        do i=1,5
        call msa_to_sky_bd(xx(i),yy(i),xd,yd)
        call roll_v3(pa_v3,xd,yd,xr(i),yr(i))
        enddo


         call pgline(5,xr,yr)



! SLIT_A_200_2

         X00=0.03466442*1000.
         Y00=-0.001256605*1000.
         w=(8.07E-05)*1000.
         h=0.00126707*1000.

         xx(1)=X00-w/2.
         yy(1)=Y00-h/2.
         xx(2)=xx(1)+w
         yy(2)=yy(1)
         xx(3)=xx(2)
         yy(3)=yy(2)+h
         xx(4)=xx(3)-w
         yy(4)=yy(3)
         xx(5)=xx(1)
         yy(5)=yy(1)

        do i=1,5
        call msa_to_sky_bd(xx(i),yy(i),xd,yd)
        call roll_v3(pa_v3,xd,yd,xr(i),yr(i))
        enddo


         call pgline(5,xr,yr)


! 1600x1600 square hole

         X00=0.02866519*1000.
         Y00=0.001557565*1000.
         w=0.00063134*1000.
         h=0.00062011*1000.

         xx(1)=X00-w/2.
         yy(1)=Y00-h/2.
         xx(2)=xx(1)+w
         yy(2)=yy(1)
         xx(3)=xx(2)
         yy(3)=yy(2)+h
         xx(4)=xx(3)-w
         yy(4)=yy(3)
         xx(5)=xx(1)
         yy(5)=yy(1)

        do i=1,5
        call msa_to_sky_bd(xx(i),yy(i),xd,yd)
        call roll_v3(pa_v3,xd,yd,xr(i),yr(i))
        enddo

         call pgline(5,xr,yr)


! SLIT_B_200

         X00=-0.03467526*1000.
         Y00=0.00275558*1000.
         w=(7.926E-05)*1000.
         h=0.00127065*1000.

         xx(1)=X00-w/2.
         yy(1)=Y00-h/2.
         xx(2)=xx(1)+w
         yy(2)=yy(1)
         xx(3)=xx(2)
         yy(3)=yy(2)+h
         xx(4)=xx(3)-w
         yy(4)=yy(3)
         xx(5)=xx(1)
         yy(5)=yy(1)

        do i=1,5
        call msa_to_sky_bd(xx(i),yy(i),xd,yd)
        call roll_v3(pa_v3,xd,yd,xr(i),yr(i))
        enddo


         call pgline(5,xr,yr)

      return
      end subroutine draw_slits_sky_rot

!---------------------------------------------------------------------------

      subroutine draw_v2v3(roll,ax_refra_p,ay_refdec_p) !x flipped version!
      implicit none
      real roll !deg
      real ax_refra_p,ay_refdec_p
      real x0,y0,xx(2),yy(2),len,roll_r
      len=5.
      x0=60.*2.516
      y0=-60.*2.516
      roll_r=roll/57.2957795
      xx(1)=x0+ax_refra_p
      yy(1)=y0+ay_refdec_p
      xx(2)=xx(1)+len*sin(roll_r)
      yy(2)=yy(1)+len*cos(roll_r)
      call pgline(2,xx,yy)
      call pgsch(0.5)
      call pgscf(1)
      xx(2)=xx(1)+(len+5.)*sin(roll_r)
      yy(2)=yy(1)+(len+5.)*cos(roll_r)
      call pgptxt(xx(2),yy(2),roll,0.5,'V3')
      xx(2)=xx(1)+len*cos(roll_r)
      yy(2)=yy(1)-len*sin(roll_r)
      call pgline(2,xx,yy)
      xx(2)=xx(1)+(len+5.)*cos(roll_r)
      yy(2)=yy(1)-(len+5.)*sin(roll_r)
      call pgptxt(xx(2),yy(2),roll,0.5,'V2')

      return
      end subroutine draw_v2v3

!------------------------------------------------------------------------------------------
 
      subroutine write_summary(outdir,gname,np,ccra,ccdec,ccpa_ap,nm,id0,id,pri,k0,i0,j0,rx0,ry0,nsm,hist_m,hist_k,hist_cat,cran1,cdecn1,k1,i1,j1,rx1,ry1,cran2,cdecn2,k2,i2,j2,rx2,ry2)     
      use config_parameters, only : n_dither,raccx,raccy,sthresh,ntile,mtarget,disperser,cpa_ap
      use target_cat, only : ntarg,nclass
      use m_list, only : n_max_m 
      use derived_parameters, only : phi
      implicit none
      real*8 ccra,ccdec,cran1,cdecn1,cran2,cdecn2
      real ccpa_ap,ccpa_v3
      character*70 outdir,filename
      integer np,nm,nsm
      integer id0(nm+nsm),id(nm+nsm),pri(nm+nsm)
      integer k0(nm+nsm),i0(nm+nsm),j0(nm+nsm)
      integer k1(n_max_m),i1(n_max_m),j1(n_max_m)
      integer k2(n_max_m),i2(n_max_m),j2(n_max_m)
      real rx0(n_max_m),ry0(n_max_m)
      real rx1(n_max_m),ry1(n_max_m)
      real rx2(n_max_m),ry2(n_max_m)
      integer i,ii,k
      integer hist_k(nclass)
      integer hist_m(nclass)
      integer hist_cat(nclass)
      integer nk
      character*6 gname
      character*2 rn

      write(rn,'(i2)') ntile
      rn=trim(adjustl(rn))
      if (ntile.lt.10) rn='0'//rn(1:1)

      ccpa_v3=ccpa_ap-(180.-phi)
      if (ccpa_v3.lt.0.) ccpa_v3=ccpa_v3+360.
      if (ccpa_v3.gt.360.) ccpa_v3=ccpa_v3-360.
      
!       print *,ccpa_v3


      filename=trim(outdir)//'/pointing_summary.txt'

      open(10,file=filename,status='replace',access='sequential', form='formatted')



      write(10,*)
      write(10,*)
      write(10,*) '------- Summary -----------'
      write(10,*)


      write(10,'(a,a)') ' Trial: ',rn
      write(10,*)

      if (n_dither.eq.1) write(10,'(a)') ' Single Pointing'
      if (n_dither.eq.2) write(10,'(a,i4)') ' Part of paired dithered pointings:',mtarget
      if (n_dither.eq.3) write(10,'(a,i4)') ' Part of triple dithered pointings:',mtarget

      write(10,*)
      write(10,*) 'IPA Pointing No:',np
      write(10,*)
      write(10,*) 'Pointing information:'
      write(10,*)
      write(10,*) 'RA, Dec of Central Pointing:'
      write(10,200) ' Nod 0:',ccra,ccdec
200   format(a7,2x,f13.7,3x,f11.7)

      write(10,*)
      write(10,*) 'Official Assigned APT/MPT roll angle:'
      write(10,210) ' PA_AP: ',cpa_ap
      write(10,*)
      write(10,*) 'Actual MSA roll angle:'
      write(10,210) ' PA_AP: ',ccpa_ap
      write(10,210) ' PA_V3: ',ccpa_v3
      write(10,*)
210   format(A8,2x,f10.6)

      write(10,*) 'RA, Dec of Nodded Pointings:'
      write(10,200)  ' Nod 1:',cran1,cdecn1 !-1 shutter down
      write(10,200)  ' Nod 2:',cran2,cdecn2 !+1 shutter up
      write(10,*)
      write(10,204) ' MSA mask intended for: ',disperser
      write(10,*)
      write(10,201) ' Acceptance Zone Thresholds:',raccx,raccy
      write(10,202) ' Acceptance Zone Open Area Filling Factor:',4.*raccx*raccy
201   format(a,x,2f6.3)
202   format(a,x,f6.3)
      write(10,*)
      write(10,203) 'Overlap Acceptance Threshold:',sthresh
203   format(a30,f6.3)
204   format(a,a)
      write(10,*)


!      calculate nk

      nk=0
      do i=1,nclass
      nk=nk+hist_k(i)
      enddo


      write(10,*) 'Total number of targets in input catalog:    ',ntarg
      write(10,*) 'Number of targets in Viable Slitlets:        ',nk
      write(10,*) 'Number of accepted targets:                  ',nm
      write(10,*) 'Number of added spare slitlets:              ',nsm
      write(10,*)



      write(10,*)
      write(10,*) 'Target breakdown by Priority Class:'
      write(10,*)
      write(10,219) '   PrCl   In Catalog   In Slitlets   Accepted   % Accepted'
219   format(a)

      do i=1,nclass
      write(10,220) i,hist_cat(i),hist_k(i),hist_m(i),float(hist_m(i))/float(hist_k(i))*100.
      enddo
220   format(3x,i2,6x,i6,7x,i5,8x,i4,9x,f6.2)
      write(10,*)



! write out selected targets and intra-shutter locations

      write(10,*)
      write(10,*) 'Accepted targets, assigned slitlets and intra-shutter locations:'
      write(10,*)
      write(10,*) ' No    ID_sub        ID_cat  pri  k0 i0  j0   rx0    ry0    k2 i2  j2   rx2    ry2    k1 i1  j1   rx1    ry1'

      do i=1,nm
       write(10,2000) i,id0(i),id(i),pri(i),k0(i),i0(i),j0(i),rx0(i),ry0(i),k2(i),i2(i),j2(i),rx2(i),ry2(i),k1(i),i1(i),j1(i),rx1(i),ry1(i)
2000  format(I4,2x,I8,2x,I12,x,I3,3x,I1,x,I3,x,I3,x,f6.3,x,f6.3,3x,I1,x,I3,x,I3,x,f6.3,x,f6.3,3x,I1,x,I3,x,I3,x,f6.3,x,f6.3)
      enddo
      write(10,*)
      write(10,*)
      write(10,2010) ' Sky Background Slitlets added to the MSA Mask:'
      write(10,*)
2010  format(a)

      if (nsm.gt.0) then
      ii=0
      write(10,*) '      k   i   j'
      do k=1,4
      do i=1,nsm
      if (k0(nm+i).eq.k) then
      ii=ii+1
      write(10,2020) ii,k0(nm+i),i0(nm+i),j0(nm+i)
      endif
      enddo !i
      enddo !k
      endif !nsm>0

2020  format(i4,3x,i1,x,i3,x,i3)

      close(10)

! write out special file with empty sky slitlets

      if (nsm.gt.0) then


      filename=trim(outdir)//'/sky_shutters.txt'

      open(10,file=filename,status='replace',access='sequential', form='formatted')


      ii=0
      write(10,*) '#     k   i   j'
      do k=1,4
      do i=1,nsm
      if (k0(nm+i).eq.k) then
      ii=ii+1
      write(10,2020) ii,k0(nm+i),i0(nm+i),j0(nm+i)
      endif
      enddo !i
      enddo !k

      close(10)
      endif !nsm>0




      return

      end subroutine write_summary

!------------------------------------------------------------------------------------------
 
      subroutine use_spare_shutters(ra_p0,dec_p0,pa_ap_p0,nsm,ksm,ism,jsm)
      implicit none
      real*8 ra_p0,dec_p0
      real pa_ap_p0
      integer c0(4,365,171),c1(4,365,171),c2(4,365,171)
!       integer, parameter :: max_dup=5
!       integer ids0(4,365,171,max_dup),ids1(4,365,171,max_dup)
!       integer ids2(4,365,171,max_dup)
      integer ns


      integer, allocatable :: ks(:),is(:),js(:)

      integer, parameter :: max_targ_k=10000 
      integer kss(max_targ_k),iss(max_targ_k),jss(max_targ_k),nss
 
      integer nsm
      integer, parameter :: max_spare_m=250 !Maximum final number of non-overlapping unused slitlets
      integer ksm(max_spare_m),ism(max_spare_m),jsm(max_spare_m)
      integer i_m(max_spare_m)

      real u
      integer istore
      integer i,j,k
! - slitmap
      integer slitmap(4,365,171)
      common/slit_map/slitmap
! - spectral parameters
!      integer n_disp
      real lambda1,lambda2
      character*6 gname
      common/disp/lambda1,lambda2,gname
          
      call make_spoiler_masks3(ra_p0,dec_p0,pa_ap_p0,c0,c1,c2)
 


! run though spare slitlets first to determine value of ns

      ns=0

      do k=1,4
       do i=1,365
        do j=1,171


       if (slitmap(k,i,j).eq.4) then

       
       if ((c0(k,i,j).gt.0).or.(c0(k,i,j+1).gt.0).or.(c0(k,i,j-1).gt.0).or.(c1(k,i,j).gt.0).or.(c1(k,i,j+1).gt.0).or.(c1(k,i,j-1).gt.0).or.(c2(k,i,j).gt.0).or.(c2(k,i,j+1).gt.0).or.(c2(k,i,j-1).gt.0)) then
       goto 66 !slitlet is contaminated skip entirely
       endif !slitlet contaminated
       
        ns=ns+1
               
       endif! slitmap(k,i,j)=4
       
66    continue

        enddo! j
       enddo! i
      enddo! k


     allocate (ks(ns),is(ns),js(ns))


! run through spare slitlets second time and capture them

      ns=0
      do k=1,4
       do i=1,365
        do j=1,171


       if (slitmap(k,i,j).eq.4) then


       
       if ((c0(k,i,j).gt.0).or.(c0(k,i,j+1).gt.0).or.(c0(k,i,j-1).gt.0).or.(c1(k,i,j).gt.0).or.(c1(k,i,j+1).gt.0).or.(c1(k,i,j-1).gt.0).or.(c2(k,i,j).gt.0).or.(c2(k,i,j+1).gt.0).or.(c2(k,i,j-1).gt.0)) then
       goto 666 !slitlet is contaminated skip entirely
       endif !slitlet contaminated
       
        ns=ns+1
                 
        ks(ns)=k
        is(ns)=i
        js(ns)=j
       
       endif! slitmap(k,i,j)=4
       
666    continue

        enddo! j
       enddo! i
      enddo! k

! randomize full spare shutter list to avoid quadrant bias in Arribas algorithm

      do i=ns,1,-1
      call random_number(u)
      j = 1 + floor(i*u)

      istore=ks(i)
      ks(i)=ks(j)
      ks(j)=istore

      istore=is(i)
      is(i)=is(j)
      is(j)=istore


      istore=js(i)
      js(i)=js(j)
      js(j)=istore

      enddo

! limit to first max_targ_k randomized empty slitlets if ns > max_targ_k to shorten selection

     nss=min(ns,max_targ_k)

     do i=1,nss
     kss(i)=ks(i)
     iss(i)=is(i)
     jss(i)=js(i)
     enddo


      call arribas_m_pick_test(nss,kss,iss,jss,nsm,i_m)

      do i=1,nsm
      ksm(i)=kss(i_m(i))
      ism(i)=iss(i_m(i))
      jsm(i)=jss(i_m(i))
      enddo

      return
      end subroutine use_spare_shutters

!---------------------------------------------------------------

      subroutine arribas_m_pick_test(nk,kt,it,jt,nm,i_m)

! determine nm non-overlapping targets from nk input
! This m_make version updates the viable slitlet map
! NB: shuttervalues already loaded in call to initialize_slitmap

      use config_parameters, only : n_disp,sthresh
      use derived_parameters, only : avxgap
      implicit none
      integer nk
      integer kt(nk),it(nk),jt(nk)
      integer, parameter :: max_spare_m=250 !Maximum final number of non-overlapping unused slitlets
      integer nm,i_m(max_spare_m)
      integer ovlap(nk,nk),tovlap(nk)
      real shv_i,shv_j
      integer n_purge,max_tovlap,idel,idum,ntout
      integer k,i,j,k2,j1,j2,m,ii,jj,dj
      real i1,i2
      real msa_gap
! - slitmap
      integer slitmap(4,365,171)
      common/slit_map/slitmap
! - shutter values
      real shval(4,365,171)
      common/shutter_values/shval
! - spectral parameters
!      integer n_disp
      real lambda1,lambda2
      character*6 gname
      common/disp/lambda1,lambda2,gname
! - prism seperations
      integer prism_sep_p(4,365,171),prism_sep_m(4,365,171)
      common/prism_length/prism_sep_p,prism_sep_m
      real sep_m,sep_p
 
      idum=-1

! Compiler warning avoidance

      idel=0
      sep_p=0.
      sep_m=0.
      ii=0
      i1=0

      k2=0 !initialzation to avoid compiler warning

      if (nk.eq.0) then
      nm=0
      return
      endif

! Fill in overlap matrix 

      do i=1,nk
      shv_i=shval(kt(i),it(i),jt(i))
       if (n_disp.eq.7) then
         i1=float(it(i))+msa_gap(kt(i))*(avxgap+365.)
         sep_p=float(prism_sep_p(kt(i),it(i),jt(i)))
         sep_m=float(prism_sep_m(kt(i),it(i),jt(i)))
        endif
      do j=1,nk
      ovlap(i,j)=0
      k=kt(i)
      if (k.eq.1) k2=3
      if (k.eq.3) k2=1
      if (k.eq.2) k2=4
      if (k.eq.4) k2=2
      if (.not.((kt(j).eq.k).or.(kt(j).eq.k2))) goto 6666
      shv_j=shval(kt(j),it(j),jt(j))
      if (abs(shv_i-shv_j).lt.sthresh) ovlap(i,j)=1
       if (n_disp.eq.7) then
         i2=float(it(j))+msa_gap(kt(j))*(avxgap+365.)
         if (((i2-i1).le.sep_m).or.((i2-i1).ge.sep_p)) ovlap(i,j)=0
        endif
6666  continue !target pair i,j not in same or adjacent quadrants - no overlap
      enddo
      enddo


! Eliminate targets in overlap matrix situated in previously accounted for slitlets
! 
!       nk_p=0
!       do i=1,nk
!       if (slitmap(kt(i),it(i),jt(i)).ne.4) then
!        nk_p=nk_p+1
!        do ii=1,nk
!        ovlap(ii,i)=0
!        ovlap(i,ii)=0
!        enddo
!       endif
!       enddo
! 
!       if (nk_p.gt.0) print *,'Nonvalid slitlets:', nk_p
! 
!       if (nk_p-nk.eq.0) then
!       nm=0
!     write(*,*)
!     write(*,*) 'No targets in non-occupied slitlets'
!     write(*,*)
!     return !No further processing needed
!       return
!       endif

      n_purge=0

! calculate marginals of overlap matrix and find largest value

      max_tovlap=-1

      do i=1,nk
      tovlap(i)=0
       do ii=1,nk
        tovlap(i)=tovlap(i)+ovlap(i,ii)
       enddo
      max_tovlap=max(max_tovlap,tovlap(i))
      if (tovlap(i).eq.max_tovlap) idel=i
      enddo

      if (max_tovlap.eq.1) goto 1200 !no overlaps to start with

! Main overlap-eliminating loop begin

1100  continue 

! Last maximum overlapping target idel to be nulled in overlap matrix

! Prepare to find next value of idel

      max_tovlap=-1

! Decrease overlap count of targets affected by target idel by one 

      do i=1,nk
      if ((ovlap(i,idel).gt.0).and.(i.ne.idel)) tovlap(i)=tovlap(i)-1
      max_tovlap=max(max_tovlap,tovlap(i))     
      if (tovlap(i).eq.max_tovlap) ii=i
      enddo

! Null idel rows and columns in overlap matrix

      ovlap(:,idel)=0 !old idel
      ovlap(idel,:)=0
      ovlap(idel,idel)=0
      tovlap(idel)=0 
     
      n_purge=n_purge+1
      
      idel=ii !new idel

      if (max_tovlap.eq.1) goto 1200 !no more overlaps

      goto 1100

1200  continue


! re-order target list and do sanity check

       nm=0
       do i=1,nk
       if (ovlap(i,i).eq.1) then
       nm=nm+1
       i_m(nm)=i
       endif
       enddo


!     write(*,*) 'Non-overlapping prism targets in non-occupied slitlets:'
!     do i=1,nm
!      write(*,*) i,tovlap(i_m(i)),i_m(i),kt(i_m(i)),it(i_m(i)),jt(i_m(i))
!     enddo
!     write(*,*)


! Fill in reshuffled overlap matrix

!        ovlap=0
! 
! 
!       do i=1,nm
!       do j=1,nm
!       ovlap(i,j)=0
!       k=kt(i_m(i))
!       if (k.eq.1) k2=3
!       if (k.eq.3) k2=1
!       if (k.eq.2) k2=4
!       if (k.eq.4) k2=2
!       if (.not.((kt(i_m(j)).eq.k).or.(kt(i_m(j)).eq.k2))) goto 7777
!       shv_i=shval(kt(i_m(i)),it(i_m(i)),jt(i_m(i)))
!       shv_j=shval(kt(i_m(j)),it(i_m(j)),jt(i_m(j)))
!       if (abs(shv_i-shv_j).lt.sthresh) ovlap(i,j)=1
!       if (n_disp.eq.7) then
!          sep_p=float(prism_sep_p(kt(i_m(i)),it(i_m(i)),jt(i_m(i))))
!          sep_m=float(prism_sep_m(kt(i_m(i)),it(i_m(i)),jt(i_m(i))))
!          i1=float(it(i_m(i)))+msa_gap(kt(i_m(i)))*(avxgap+365.)
!          i2=float(it(i_m(j)))+msa_gap(kt(i_m(j)))*(avxgap+365.)
!          if (((i2-i1).le.sep_m).or.((i2-i1).ge.sep_p)) ovlap(i,j)=0
!       endif
! 7777  continue
!       enddo
!       enddo
! 
! 
!        do i=1,nm
!        tovlap(i)=0
!        do j=1,nk
!        tovlap(i)=tovlap(i)+ovlap(i,j)
!        enddo
!        enddo
! 
! 
!     write(*,*) 'Non-overlapping prism targets in non-occupied slitlets:'
!     do i=1,nm
!      write(*,*) i,tovlap(i),i_m(i),kt(i_m(i)),it(i_m(i)),jt(i_m(i))
!     enddo
!     write(*,*)


! Mark commanded open target slitlets in slitmap

      do m=1,nm

      k=kt(i_m(m))
      i=it(i_m(m))
      j=jt(i_m(m))

      if (slitmap(k,i,j).gt.1) slitmap(k,i,j)=-1

      if (j+1.le.171) then
      if (slitmap(k,i,j+1).gt.1) slitmap(k,i,j+1)=-1
      endif

      if (j-1.ge.1) then
      if (slitmap(k,i,j-1).gt.1) slitmap(k,i,j-1)=-1
      endif

! Map out shutters excluded due to target overlaps in slitmap

      ntout=0
      dj=int(2.*sthresh)
      j1=max(j-dj,1)
      j2=min(j+dj,171)
      if (k.eq.1) k2=3
      if (k.eq.3) k2=1
      if (k.eq.2) k2=4
      if (k.eq.4) k2=2

      do ii=1,365
      do jj=j1,j2

      if (abs(shval(k,i,j)-shval(k,ii,jj)).le.sthresh) then !shutter cla

       if (n_disp.eq.7) then
         sep_p=float(prism_sep_p(k,i,j))
         sep_m=float(prism_sep_m(k,i,j))
         i1=float(i)+msa_gap(k)*(avxgap+365.)
         i2=float(ii)+msa_gap(k)*(avxgap+365.)
         if (((i2-i1).le.sep_m).or.((i2-i1).ge.sep_p)) goto 323 ! except perhaps prism
       endif

        if ((slitmap(k,ii,jj).ge.2).and.(slitmap(k,ii,jj).le.4)) then
        slitmap(k,ii,jj)=10
        ntout=ntout+1
        endif

323   continue
      endif


      if (abs(shval(k,i,j)-shval(k2,ii,jj)).le.sthresh) then !repeat for adjacent quadrant

       if (n_disp.eq.7) then
         sep_p=float(prism_sep_p(k,i,j))
         sep_m=float(prism_sep_m(k,i,j))
         i1=float(i)+msa_gap(k)*(avxgap+365.)
         i2=float(ii)+msa_gap(k2)*(avxgap+365.)
         if (((i2-i1).le.sep_m).or.((i2-i1).ge.sep_p)) goto 324
       endif


        if ((slitmap(k2,ii,jj).ge.2).and.(slitmap(k2,ii,jj).le.4)) then
        slitmap(k2,ii,jj)=10
        ntout=ntout+1
        endif

324   continue
      endif

      enddo
      enddo

      enddo

!     write(*,*) 'Vacant slitlets claimed:',ntout
!     write(*,*)


      return
      end subroutine arribas_m_pick_test

!-----------------------------------------------------------------

      subroutine common_target_list(nm1,id01,id1,cl1,nm2,id02,id2,cl2,nm3,id03,id3,cl3,np1,np2,np3,outfile)
      use config_parameters, only : n_dither
      use target_cat, only : nclass
      use m_list, only : n_max_m
      use histograms, only : hist_com,avexp,n_com,n_spec
      implicit none
      integer ncl(nclass)
      integer id01(n_max_m),id1(n_max_m),cl1(n_max_m)
      integer id02(n_max_m),id2(n_max_m),cl2(n_max_m)
      integer id03(n_max_m),id3(n_max_m),cl3(n_max_m)
      integer nm1,nm2,nm3,nt,np1,np2,np3
      integer idt(3*n_max_m),clt(3*n_max_m)
      integer id0t(3*n_max_m)
      integer id(3,n_max_m)
      integer id0(3,n_max_m)
      logical d(3,3*n_max_m)
      integer im,i1,i2,i3,it,ip,ic
      character*70 outfile
      integer nm(3)
      integer nexp,maxexp,minexp
      integer n_tot,n_single,n_double,n_triple,n_cl
      real frac_triple,frac_double,frac_single
      integer ntnum


! Compiler warning avoidance

      ntnum=0

      if (n_dither.eq.1) then
      nm2=0
      nm3=0
      endif

      if (n_dither.eq.2) then
      nm3=0
      endif


! construct combined target list

      nt=0

      do im=1,nclass

! grab class m targets in list 1 and add to total list

      do i1=1,nm1
        if (cl1(i1).eq.im) then
        nt=nt+1
        idt(nt)=id1(i1)
        id0t(nt)=id01(i1)
        clt(nt)=cl1(i1)
        endif
      enddo !i1


! grab further unique targets in list 2 and add to total list

       do i2=1,nm2
        if (cl2(i2).eq.im) then
          do it=1,nt
          if (id2(i2).eq.idt(it)) goto 100 !target already in list, skip
          enddo !it
          nt=nt+1
          idt(nt)=id2(i2)
          id0t(nt)=id02(i2)
          clt(nt)=cl2(i2)
100       continue
        endif
      enddo !i2

! grab further unique targets in list 3 and add to total list

       do i3=1,nm3
        if (cl3(i3).eq.im) then
          do it=1,nt
          if (id3(i3).eq.idt(it)) goto 200 !target already in list, skip
          enddo !it
          nt=nt+1
          idt(nt)=id3(i3)
          id0t(nt)=id03(i3)
          clt(nt)=cl3(i3)
200       continue
        endif
      enddo !i3

      enddo !m


      nm(1)=nm1
      do im=1,nm1
      id(1,im)=id1(im)
      id0(1,im)=id01(im)
      enddo

      nm(2)=nm2
      do im=1,nm2
      id(2,im)=id2(im)
      id0(2,im)=id02(im)
      enddo

      nm(3)=nm3
      do im=1,nm3
      id(3,im)=id3(im)
      id0(3,im)=id03(im)
      enddo

      do ip=1,3
      do it=1,3*n_max_m
      d(ip,it)=.false.
      enddo
      enddo

! map out which targets are in which lists

      do it=1,nt

      do ip=1,n_dither

      d(ip,it)=.false.

      do im=1,nm(ip)
      if (id(ip,im).eq.idt(it)) d(ip,it)=.true.
      enddo

      enddo

      enddo


      open(10,file=outfile,status='replace',access='sequential', form='formatted')
      write(10,*)

      if (n_dither.gt.1) then
      write(10,*) 'Combined target list:'
      if (n_dither.eq.2) write(10,*) '  No   No_sub        No_cat  Pri P1 P2'
      if (n_dither.eq.3) write(10,*) '  No   No_sub        No_cat  Pri P1 P2 P3'
      do it=1,nt
      write(10,1000) it,id0t(it),idt(it),clt(it),(d(ip,it),ip=1,n_dither)
      enddo
      endif !n_dither>1

      if (n_dither.eq.1) then
      write(10,*) 'Target list:'
      write(10,*) '  No   No_sub        No_cat  Pri'
      do it=1,nt
      write(10,1001) it,id0t(it),idt(it),clt(it)
      enddo
      endif !n_dither=1


1000  format(i5,x,i8,x,i13,x,i3,2x,3L2)
1001  format(i5,x,i8,x,i13,x,i3)


! write out targets in each catalog

      write(10,*)

      if (n_dither.eq.1) then
      write(10,*) 'Number of targets:',nt
      write(10,*)
      endif !n_dither=1



      if (n_dither.gt.1) then

      write(*,*)
      write(*,*) 'Number of targets per dithered pointing:'
      im=0
      do ip=1,n_dither
      if (ip.eq.1) write(*,'(a,i0,a,i0)') ' IPA Pointing ',np1,': ',nm(ip)
      if (ip.eq.2) write(*,'(a,i0,a,i0)') ' IPA Pointing ',np2,': ',nm(ip)
      if (ip.eq.3) write(*,'(a,i0,a,i0)') ' IPA Pointing ',np3,': ',nm(ip)
      im=im+nm(ip)
      enddo
      write(*,*)
      write(*,*) 'Total number of spectra: ',im
      write(*,*) 'Number of unique targets:',nt
      write(*,*)

      write(10,*)
      write(10,*)
      write(10,*) 'Number of targets per dithered pointing:'
      im=0
      do ip=1,n_dither
      if (ip.eq.1) write(10,*) 'Dither Pointing 1:',nm(ip)
      if (ip.eq.2) write(10,*) 'Dither Pointing 2:',nm(ip)
      if (ip.eq.3) write(10,*) 'Dither Pointing 3:',nm(ip)
      im=im+nm(ip)
      enddo
      write(10,*)
      write(10,*) 'Total number of spectra: ',im
      write(10,*) 'Number of unique targets:',nt
      write(10,*)


      endif !n_dither>1


! determine maximum multiple exposure rate

      maxexp=-666
      minexp=666

       do it=1,nt
      nexp=0
       do ip=1,n_dither
       if (d(ip,it)) nexp=nexp+1
       enddo
      maxexp=max(maxexp,nexp)
      minexp=min(minexp,nexp)
      enddo


      if (n_dither.gt.1) then

      write(10,*) 'Minimum exposure:',minexp
      write(10,*) 'Maximum exposure:',maxexp

! Breakdown by priority class

      write(10,*)
      write(10,*) 'Breakdown by priority Class:'

      endif !n_dither>1

      do ic=1,nclass

      n_cl=ic

      n_tot=0
      n_single=0
      n_double=0
      n_triple=0


      do it=1,nt

      if (clt(it).eq.n_cl) then

      n_tot=n_tot+1

       nexp=0
       do ip=1,n_dither
       if (d(ip,it)) nexp=nexp+1
       enddo

      if (nexp.eq.0) stop 'Null Exposure!'
      if (nexp.eq.1) n_single=n_single+1
      if (nexp.eq.2) n_double=n_double+1
      if (nexp.eq.3) n_triple=n_triple+1

      endif

      enddo !it

      frac_single=0.
      frac_double=0.
      frac_triple=0.

      if (n_tot.gt.0) then
      frac_single=float(n_single)/float(n_tot)*100.
      frac_double=float(n_double)/float(n_tot)*100.
      frac_triple=float(n_triple)/float(n_tot)*100.
      endif

       ncl(ic)=n_tot
       avexp(ic)=(frac_single+2.*frac_double+3.*frac_triple)/100.


!     if (ic.gt.1) then

      if (n_dither.gt.1) then

      write(10,*)
      write(10,*)
      write(10,98) ' Target Class:',n_cl
      write(10,*)
      if (n_dither.gt.1) then
      write(10,99) ' Unique targets in Class:',n_tot
      write(10,*)
      if (n_dither.gt.2) write(10,99) ' Triple Exposed Targets:  ',n_triple,frac_triple,' %'
      if (n_dither.gt.1) write(10,99) ' Double Exposed Targets:  ',n_double,frac_double,' %'
      write(10,99) ' Single Exposed Targets:  ',n_single,frac_single,' %'
      endif
      write(10,*) ' Checksum:',n_tot-n_single-n_double-n_triple

      endif !n_dither>1

!     endif !ic

      enddo !ic

98    format(a,I5)
99    format(a,I5,f6.1,a2)


      if (n_dither.gt.1) then

! same for combined Sample

      n_tot=0
      n_single=0
      n_double=0
      n_triple=0

      do it=1,nt

      n_tot=n_tot+1

       nexp=0
       do ip=1,n_dither
       if (d(ip,it)) nexp=nexp+1
       enddo

      if (nexp.eq.0) stop 'Null Exposure!'
      if (nexp.eq.1) n_single=n_single+1
      if (nexp.eq.2) n_double=n_double+1
      if (nexp.eq.3) n_triple=n_triple+1

      enddo !it


      if (n_tot.gt.0) then
      frac_single=float(n_single)/float(n_tot)*100.
      frac_double=float(n_double)/float(n_tot)*100.
      frac_triple=float(n_triple)/float(n_tot)*100.
      endif


      write(10,*)
      write(10,*)
      write(10,*)
      write(10,96) ' Combined Sample:'
      write(10,*)
      write(10,97) ' Unique targets all Classes:',n_tot
      write(10,*)
      write(10,99)
      if (n_dither.gt.2) write(10,99) ' Triple Exposed Targets:  ',n_triple,frac_triple,' %'
      if (n_dither.gt.1) write(10,99) ' Double Exposed Targets:  ',n_double,frac_double,' %'
      write(10,99) ' Single Exposed Targets:  ',n_single,frac_single,' %'
 96    format(a,I5)
 97    format(a,I5)
      write(10,*) ' Checksum:',n_tot-n_single-n_double-n_triple
      write(10,*)


      write(10,*)
      write(10,*) 'Summary:'
      write(10,*)
      write(10,*) 'Unique target breakdown by Priority Class:'
      write(10,*)
      write(10,93) '  Class  Number  AvExp'

      ntnum=0

      do ic=1,nclass
      n_cl=ic
      write(10,95) n_cl,ncl(ic),avexp(ic)
      ntnum=ntnum+ncl(ic)
      hist_com(ic)=ncl(ic) !load commons for use outside subroutine
      enddo !ic


      write(10,94) '  Tot',nt

      write(10,*) 'Checksum:',nt-ntnum
      write(10,*)
95    format(i5,2x,i6,2x,f6.2)
94    format(a5,2x,i6)
93    format(a22)



      write(*,*) 'Unique target breakdown by Priority Class:'
      write(*,*)
      write(*,93) '  Class  Number  AvExp'

      ntnum=0

      do ic=1,nclass
      n_cl=ic
      write(*,95) n_cl,ncl(ic),avexp(ic)
      ntnum=ntnum+ncl(ic)
      enddo !ic

      write(*,94) 'Tot',nt

      write(*,*)

      endif !n_dither>1


      if (n_dither.eq.1) then

      write(*,*)
      write(*,*) 'Target breakdown by Priority Class:'
      write(*,*)
      write(*,*) '  Pri    Number'
      ntnum=0

      do ic=1,nclass
      n_cl=ic
      write(*,87) n_cl,ncl(ic)
      ntnum=ntnum+ncl(ic)
      enddo !ic
87    format(i5,2x,i6)
      write(*,94) 'Tot',nt
      write(*,*)


      write(10,*)
      write(10,*) 'Target breakdown by Priority Class:'
      write(10,*)
      write(10,*) '  Pri    Number'
      ntnum=0

      do ic=1,nclass
      n_cl=ic
      write(10,87) n_cl,ncl(ic)
      ntnum=ntnum+ncl(ic)
      enddo !ic
       write(10,*)
!     write(10,94) 'Tot',nt
      write(10,*)
      write(10,*)


      endif !n_dither=1

      close(10)

! load common statistics for use outside subroutine

      n_spec=nm1+nm2+nm3
      n_com=ntnum

      do ic=1,nclass
      hist_com(ic)=ncl(ic)
      enddo


      return

      end subroutine common_target_list

!------------------------------------------------------------

      subroutine plot_spectra(outdir,n_dispx,nm,id,cl,kt,it,jt,nsm)
! Label moved to left side
      use path_to_ref_files
      use target_cat, only : nclass
!       use m_list, only : n_max_m
      implicit none
      integer n_dispx
      integer nm,nsm
      integer id(nm+nsm),cl(nm+nsm)
      integer kt(nm+nsm),it(nm+nsm),jt(nm+nsm)
      integer, parameter :: max_failo=50 !maximum number of failed open shutters
      integer nfailo,ko(max_failo),io(max_failo),jo(max_failo)
      common/failed_open/nfailo,ko,io,jo
      real x1(5),x2(5),y1(5),y2(5),xmb,ymb,xmt,ymt,xm,ym,xp1,yp1,xp2,yp2
      real xf(1000),yf(1000)
!       logical skypad
      integer i,j,k
      character*70 outdir
      character*70 filename
      character*50 header
      character*8 tnum
      integer ntnum
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
! - spectral parameters
      integer nwave
!      integer n_disp
      real lambda1,lambda2,wave
      character*6 gname
      common/disp/lambda1,lambda2,gname


      call msa_to_fpa_disp(n_dispx,wave,xmb,ymb,xmt,ymt)

      filename=trim(outdir)//'/'//gname(1:5)//'_spectra.ps/cps'


! Plot spectra

      nwave=100

        xp1=-45.
        xp2=45.
        yp1=-20.
        yp2=20.

        CALL PGBEGIN(0,filename,1,1)
        CALL PGSLS(1)
        CALL PGSLW(1)
        CALL PGSCF(2)
        CALL PGSCH(1.)

        CALL PGVPORT(.0,1.,.0,1.)
        CALL PGWNAD(xp1,xp2,yp1,yp2)
        CALL PGSLW(1)

! draw FPA outline

      call load_fpa_dimensions_bd

      call pixel_to_fpa_bd(1,1,1,x1(1),y1(1))
      call pixel_to_fpa_bd(1,1,2048,x1(2),y1(2))
      call pixel_to_fpa_bd(1,2048,2048,x1(3),y1(3))
      call pixel_to_fpa_bd(1,2048,1,x1(4),y1(4))
      x1(5)=x1(1)
      y1(5)=y1(1)
      call pixel_to_fpa_bd(2,1,1,x2(1),y2(1))
      call pixel_to_fpa_bd(2,1,2048,x2(2),y2(2))
      call pixel_to_fpa_bd(2,2048,2048,x2(3),y2(3))
      call pixel_to_fpa_bd(2,2048,1,x2(4),y2(4))
      x2(5)=x2(1)
      y2(5)=y2(1)

      call pgsci(15)
      call pgline(5,x1,y1)
      call pgline(5,x2,y2)

      call load_msa_dimensions_bd


! intialize disperser

      call msa_to_fpa_disp(n_dispx,wave,xmb,ymb,xmt,ymt)
      call pgslw(1)
!      if (n_disp.eq.1) then
!      lambda1=0.7
!      lambda2=1.4
!      endif

      call read_shutter_values
      call read_base_msamap
      call make_slitlet_map


! set colors

      call set_rainbow_lut


! Plot failed open spectra


      do i=1,nfailo

! central shutter of slitlet

      call shutter_to_msa_span_bd(ko(i),io(i),jo(i),xm,ym,xmt,ymt,xmb,ymb)

      do j=1,nwave
      wave=lambda1+(lambda2-lambda1)/float(nwave-1)*float(j-1)
      call msa_to_fpa_disp(n_dispx,wave,xmt,ymt,xf(j),yf(j))
      enddo

      call pgsci(6)
      call pgline(nwave,xf,yf)

      do j=1,nwave
      wave=lambda1+(lambda2-lambda1)/float(nwave-1)*float(j-1)
      call msa_to_fpa_disp(n_dispx,wave,xmb,ymb,xf(j),yf(j))
      enddo

      call pgline(nwave,xf,yf)

      enddo



      do i=nm+nsm,1,-1 !all spectra bottom of list first

! determine priority class color

      k=cl(i)
      call set_rainbow_color_fine(nclass,k)

      if (cl(i).eq.0) call pgsci(15) !empty shutters

!  top shutter of slitlet

      call shutter_to_msa_span_bd(kt(i),it(i),jt(i)+1,xm,ym,xmt,ymt,xmb,ymb)

      do j=1,nwave
      wave=lambda1+(lambda2-lambda1)/float(nwave-1)*float(j-1)
      call msa_to_fpa_disp(n_dispx,wave,xmt,ymt,xf(j),yf(j))
      enddo

      call pgline(nwave,xf,yf)

! bottom shutter of slitlet

      call shutter_to_msa_span_bd(kt(i),it(i),jt(i)-1,xm,ym,xmt,ymt,xmb,ymb)

      do j=1,nwave
      wave=lambda1+(lambda2-lambda1)/float(nwave-1)*float(j-1)
      call msa_to_fpa_disp(n_dispx,wave,xmb,ymb,xf(j),yf(j))
      enddo

      call pgline(nwave,xf,yf)

! central shutter of slitlet

      call shutter_to_msa_span_bd(kt(i),it(i),jt(i),xm,ym,xmt,ymt,xmb,ymb)

      do j=1,nwave
      wave=lambda1+(lambda2-lambda1)/float(nwave-1)*float(j-1)
      call msa_to_fpa_disp(n_dispx,wave,xm,ym,xf(j),yf(j))
      enddo

      call pgline(nwave,xf,yf)

      if (cl(i).ge.1) then
      call pgnumb(id(i),0,1,tnum,ntnum) !label catalog ID
!       call pgnumb(i,0,1,tnum,ntnum) !label output catalog number
      call pgsci(1)
      call pgsch(0.4)
      call pgscf(1)
      call pgptext(xf(1)-0.12,yf(1)-0.2,0.,1.0,tnum(1:ntnum))
      call pgscf(2)
      call pgsch(1.0)
      endif

      enddo


      call pgnumb(nm,0,1,tnum,ntnum)
       header=gname//' '//tnum(1:ntnum)//' Target Spectra'
!     if (np.eq.1) header=
!    #gname//' '//tnum(1:ntnum)//' Overlapping Spectra'
!     if (np.eq.2) header=
!    #gname//' '//tnum(1:ntnum)//' Non-Overlapping Spectra'
      call pgsci(1)
      call pgptext(0.,yp2+1.3,0.,0.5,header)

      call pgnumb(nsm,0,1,tnum,ntnum)
      header=tnum(1:ntnum)//' Sky Background Spectra'
      call pgptext(0.,yp2-0.8,0.,0.5,header)


      call pgend

      return
      end subroutine plot_spectra

!-------------------------------------------------------------------
