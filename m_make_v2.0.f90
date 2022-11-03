 
! V1.8 target ID output increased from 5 to 7 digits
! v1.9 rewritten to handle filter-dependent FORE transformation      
! v2.0 increased ID field to I13 in m_list file format statement 200

     include 'empt_modules.f90'

!-------------------------------------------------------

      program m_make
      use path_to_ref_files
      use config_parameters
      use derived_parameters
      use target_cat, only : nclass
      use pointings_and_groups, only : n_p,no_p,ra_p,dec_p,pa_ap_p
      use k_list
      use order_and_weights
      implicit none
      character*70 outfile
      character*2 rn
      real*8 av_ra_p,av_dec_p
      real pa_v3_p0
      integer npointings
      real, allocatable :: tx_p(:),ty_p(:)
      real ax,ay,aax,aay
      logical l1,l2,l3
      real dx12,dy12,dx13,dy13,dx23,dy23
      integer nm
      integer, allocatable :: i_m(:)
      real xx,yy
      integer ip,ik,ip1,ip2,ip3,i_c
      integer npairs,npaircount,ntriplecount,ntriples,n_possible
      integer nsingles,nsinglecount
      integer, allocatable :: id1(:),id2(:),id3(:)
      integer, allocatable :: idc1(:),idp1(:),idu1(:)
      integer, allocatable :: idc2(:),idp2(:),idu2(:)
      integer, allocatable :: idc3(:),idp3(:),idu3(:)
      integer ik1,ik2,ik3
      integer nk1,nk2,nk3
      integer nu1,nu2,nu3,ncom1,ncom2,ncom3,np1,np2,np3
      integer nkc
      integer, allocatable :: ktc(:),itc(:),jtc(:)
      integer, allocatable :: i_nkc(:)
      integer, allocatable :: i_mc(:)
      integer im,ic,nmc,ikk,no,nn,jf
      integer, allocatable :: ipp1(:),ipp2(:),ipp3(:)

      real t_start,t_stop

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
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
! - prism seperations
      integer prism_sep_p(4,365,171),prism_sep_m(4,365,171)
      common/prism_length/prism_sep_p,prism_sep_m

   
      call cpu_time(t_start)
      
      nn=0 !initialization to avoid compiler warning
      ikk=0 !initialization to avoid compiler warning

 
      write(*,*)
      write(*,*) '-------------------------------------------------------------------------------------------'
      write(*,*) '--- m_make -- eMPT Arribas Algorithm m-list Generation Module Version 2.0 4 August 2022 ---'
      write(*,*) '-------------------------------------------------------------------------------------------'
      write(*,*)

      call load_derived_parameters

      call file_entry_m_make()
           
      write(rn,'(i2)') ntile
      rn=trim(adjustl(rn))
      if (ntile.lt.10) rn='0'//rn(1:1)
      
      k_list_mod_file='trial_'//rn(1:2)//'/k_clean_output/k_list_mod.txt'

      call pointings_read()
      
      npointings=n_p

      if (npointings.ne.n_p) then
      write (*,*)
      write (*,*) ' ******************************************************************************' 
      write (*,*) ' Error: Discrepancy in number of pointings in configuration and k_list_mod file' 
      write (*,*) ' ******************************************************************************' 
      write (*,*)
      stop
      endif

      
      call k_list_mod_read()

       write(*,*)
       write(*,'(a,a2)') ' Trial Identifier No.:        ',rn
       write(*,'(a,a)')  ' Configuration File:          ',confname
       write(*,'(a,i0)') ' n_dither:                    ',n_dither
       write(*,'(a,i0)') ' Number of Priority Classes in input catalog:      ',nclass
       write(*,'(a,i0)') ' Total number of pointings in configuration file:  ',n_p
       write(*,'(a,a)')  ' Modified input k_list file:                       ',k_list_mod_file
       write(*,'(a,i0)') ' Maximum number of targets in k_list per pointing: ',n_max_k
       write(*,*)

 
      call system('mkdir trial_'//rn(1:2)//'/m_make_output')

      if (n_dither.eq.1) outfile='trial_'//rn(1:2)//'/m_make_output/single_m_list.txt'
      if (n_dither.eq.2) outfile='trial_'//rn(1:2)//'/m_make_output/pair_m_list.txt'
      if (n_dither.eq.3) outfile='trial_'//rn(1:2)//'/m_make_output/triple_m_list.txt'

      allocate (i_m(n_max_k))
      allocate (ktc(n_max_k),itc(n_max_k),jtc(n_max_k)) 
      allocate (i_nkc(n_max_k))
      allocate (i_mc(n_max_k))


! create tangential coordinates of pointing positions if required for legal pair or triple assessment

      if (n_dither.gt.1) then

! calculate mean ra and dec

      av_ra_p=0.
      av_dec_p=0.

      do ip=1,npointings
      av_ra_p=av_ra_p+ra_p(ip)
      av_dec_p=av_dec_p+dec_p(ip)
      enddo

      av_ra_p=av_ra_p/float(npointings)
      av_dec_p=av_dec_p/float(npointings)

      allocate (tx_p(n_p),ty_p(n_p))      

! calculate relative pointing positions in MSA plane orientation

      do ip=1,npointings

      pa_v3_p0=pa_ap_p(ip)-(180.-phi)
      if (pa_v3_p0.lt.0.) pa_v3_p0=pa_v3_p0+360.
      if (pa_v3_p0.gt.360.) pa_v3_p0=pa_v3_p0-360.

      call tangcoord_single(av_ra_p,av_dec_p,ra_p(ip),dec_p(ip),aax,aay)

! Rotate to V2,V3

      call roll_v3(-pa_v3_p0,aax,aay,ax,ay) !V2,V3

! Rotate to MSA coordinate system orientation

      tx_p(ip)= ax*cos(phi_r)+ay*sin(phi_r) !xm,ym as
      ty_p(ip)=-ax*sin(phi_r)+ay*cos(phi_r)

      tx_p(ip)=-tx_p(ip)*1000./av_xpitch !mas
      ty_p(ip)=-ty_p(ip)*1000./av_ypitch

      enddo

      deallocate (ra_p,dec_p,pa_ap_p)

      endif !n_dither>1


! Initialize mode-specific subroutines (needed before shutter_values and slitmap construction)

      call load_msa_dimensions_bd
      call load_fpa_dimensions_bd
      call msa_to_fpa_disp(n_disp,1.,xx,yy,xx,yy)

! read in shutter values and set up viable shutter map

      call read_shutter_values
      if (n_disp.eq.7) call read_prism_separations
      call read_base_msamap



!*************************************

       if (n_dither.eq.1) then

       nsingles=npointings

       allocate (ipp1(nsingles)) 
       allocate (id1(n_max_k))
       allocate (idu1(n_max_k))

       write(*,*)
       write(*,110) ' Number of single pointings: ',nsingles
       write(*,110) ' Number of Priority Classes: ',nclass
       write(*,*)


!      open output m-space target list file

       open(10,file=outfile,status='replace',access='sequential', form='formatted')

       write(10,*) nsingles


! Begin pointing loop

      do nsinglecount=1,nsingles

      ip1=nsinglecount


      write(10,*) nsinglecount
      write(10,*) ip1  !pointing list sequence numbers
      write(10,*) no_p(ip1) !IPA Pointing list sequence number


      write(*,88) ' Processing single pointing ',nsinglecount,' of ',nsingles
88    format(a,i0,a,i0)

      nk1=nk(ip1)

      do ik1=1,nk1
      id1(ik1)=id_cat(ip1,ik1)
      enddo


! move k-list 1 into unique part (nu1)

      nu1=0

      do ik1=1,nk1
      nu1=nu1+1
      idu1(nu1)=ik1
      enddo !ik1


! Apply Arribas algorithm to k-list 1 in order of priority classes

! reset slitlet map for k-list 1

      call make_slitlet_map

      call initialize_slitmap(n_disp,sthresh,1,1)

      nm=0

      do no=1,nclass

      ic=place_order(no,1)
      nn=nu1

      nkc=0
      do ik=1,nn
      ikk=idu1(ik)
      if (pri_k(ip1,ikk).eq.ic) then
      nkc=nkc+1
      ktc(nkc)=kt(ip1,ikk)
      itc(nkc)=it(ip1,ikk)
      jtc(nkc)=jt(ip1,ikk)
      i_nkc(nkc)=ikk !full k-space input catalog reference
      endif
      enddo

      call arribas_m_make_test(nkc,ktc,itc,jtc,nmc,i_mc)

! Gather placed objects for pointing in sequence

      do im=1,nmc
      nm=nm+1
      i_m(nm)=i_nkc(i_mc(im))
      enddo

      enddo !no placement sequence order loop


      write(10,*) nm
!     nm1=nm

      do im=1,nm
       write(10,200) id(ip1,i_m(im)),id_cat(ip1,i_m(im)),pri_k(ip1,i_m(im)),kt(ip1,i_m(im)),it(ip1,i_m(im)),jt(ip1,i_m(im)),rx(ip1,i_m(im)),ry(ip1,i_m(im))
      enddo

      enddo ! singlecount

      close(10)

      write(*,*)

      endif !n_dither=1

!*************************************


       if (n_dither.eq.2) then


       n_possible=0
       do ip1=1,npointings
       do ip2=1,ip1-1
       n_possible=n_possible+1
       enddo
       enddo


! first pass through to establish number of legal pairs

       npairs=0

       do ip1=1,npointings !march through ip1

       do ip2=1,ip1-1 !march through ip2 for each ip1

! check legality of ip2 versus ip1

       dy12=abs(ty_p(ip1)-ty_p(ip2))

       if (dy12.gt.max_dy) goto 67 !ip2 too far away from ip1 in y, next

       dx12=abs(tx_p(ip1)-tx_p(ip2))

       if (dx12.gt.max_dx) goto 67 !ip2 too far away from ip1 in x, next

       if ((dy12.lt.min_dy).and.(dx12.le.min_dx)) goto 67 !ip2 in same r


! ip1 and ip2 are a mutually legal pairs at this point

       npairs=npairs+1


67     continue !next ip2
       enddo !ip2

       enddo !ip1


       write(*,*)
       write(*,111) ' Number of pointings:        ',npointings
       write(*,111) ' Number of possible pairs:   ',n_possible
       write(*,111) ' Number of legal pairs:      ',npairs
       write(*,111) ' Number of Priority Classes: ',nclass
111    format(a,i0)
       write(*,*)

      if (npairs.lt.1) then
      write(*,*)
      write(*,*) '*****************************************************'
      write(*,*)
      write(*,*) '          No legal dither pairs found!'
      write(*,*) '             Program m_make stopped'
      write(*,*)
      write(*,*) '  Either relax the legal pair proximity rules in the'
      write(*,*) '  configuration file and run m_make again, and/or'
      write(*,*) '  add additional less optimal pointings from the IPA'
      write(*,*) '  output to the configuration file for processing'
      write(*,*) '  and run k_make followed by m_clean and m_make again'
      write(*,*)
      write(*,*) '*****************************************************'
      write(*,*)
      stop
      endif !no legal pairs



! 2nd pass through to fill out ipp1 and ipp2 arrays

       allocate (ipp1(npairs),ipp2(npairs)) 

       npairs=0

       do ip1=1,npointings !march through ip1

       do ip2=1,ip1-1 !march through ip2 for each ip1

! check legality of ip2 versus ip1

       dy12=abs(ty_p(ip1)-ty_p(ip2))

       if (dy12.gt.max_dy) goto 670 !ip2 too far away from ip1 in y, next

       dx12=abs(tx_p(ip1)-tx_p(ip2))

       if (dx12.gt.max_dx) goto 670 !ip2 too far away from ip1 in x, next

       if ((dy12.lt.min_dy).and.(dx12.le.min_dx)) goto 670 !ip2 in same r


! ip1 and ip2 are a mutually legal pairs at this point

       npairs=npairs+1

       ipp1(npairs)=ip1 !target list sequence numbers
       ipp2(npairs)=ip2


670    continue !next ip2
       enddo !ip2

       enddo !ip1


!      open output m-space target list file

       open(10,file=outfile,status='replace',access='sequential', form='formatted')

       write(10,*) npairs


! Begin pointing pair loop

      allocate (id1(n_max_k),id2(n_max_k))
      allocate (idc1(n_max_k),idu1(n_max_k))
      allocate (idc2(n_max_k),idu2(n_max_k))

      do npaircount=1,npairs

      ip1=ipp1(npaircount)
      ip2=ipp2(npaircount)


      write(10,*) npaircount
      write(10,*) ip1,ip2   !pointing list sequence numbers
      write(10,*) no_p(ip1),no_p(ip2) !IPA Pointing list sequence numbers



      write(*,99) ' Processing pointing pair ',npaircount,' of ',npairs
99    format(a,i0,a,i0)

      nk1=nk(ip1)
      nk2=nk(ip2)

      do ik1=1,nk1
      id1(ik1)=id_cat(ip1,ik1)
      enddo

      do ik2=1,nk2
      id2(ik2)=id_cat(ip2,ik2)
      enddo

! sort k-list 1 into common (ncom1) and unique parts (nu1)

      ncom1=0
      nu1=0

      do ik1=1,nk1


      l2=.false.
      do ik2=1,nk2
        if (id1(ik1).eq.id2(ik2)) then
         l2=.true.
         goto 555
        endif
       enddo !ik2

555    continue


        if (l2) then
         ncom1=ncom1+1
         idc1(ncom1)=ik1
         goto 558
        endif


        if (.not.l2) then
         nu1=nu1+1
         idu1(nu1)=ik1
         goto 558
        endif

558   continue

      enddo !ik1


! sort k-list 2 into common (ncom2), partial (np2) and unique parts (nu2)

      ncom2=0
      nu2=0

       do ik2=1,nk2

      l1=.false.
      do ik1=1,nk1
        if (id2(ik2).eq.id1(ik1)) then
         l1=.true.
         goto 559
        endif
       enddo !ik1

559    continue


        if (l1) then
         ncom2=ncom2+1
         idc2(ncom2)=ik2
         goto 560
        endif


        if (.not.l1) then
         nu2=nu2+1
         idu2(nu2)=ik2
         goto 560
        endif



560   continue
      enddo !ik2



! Apply Arribas algorithm to partitoned k-list 1 in order specified by place_order sequence

! reset slitlet map for k-list 1

      call make_slitlet_map

      call initialize_slitmap(n_disp,sthresh,1,1)

      nm=0

      do no=1,2*nclass

      ic=place_order(no,1) !priority class set

      jf=place_order(no,2) !k_list portion set 1=com,2=unique
      if (jf.eq.1) nn=ncom1
      if (jf.eq.2) nn=nu1

      nkc=0
      do ik=1,nn
      if (jf.eq.1) ikk=idc1(ik)
      if (jf.eq.2) ikk=idu1(ik)
      if (pri_k(ip1,ikk).eq.ic) then
      nkc=nkc+1
      ktc(nkc)=kt(ip1,ikk)
      itc(nkc)=it(ip1,ikk)
      jtc(nkc)=jt(ip1,ikk)
      i_nkc(nkc)=ikk !full k-space input catalog reference
      endif
      enddo

      call arribas_m_make_test(nkc,ktc,itc,jtc,nmc,i_mc)

! Gather placed objects for pointing in sequence

      do im=1,nmc
      nm=nm+1
      i_m(nm)=i_nkc(i_mc(im))
      enddo

      enddo !no placement sequence order loop


      write(10,*) nm

      do im=1,nm
       write(10,200) id(ip1,i_m(im)),id_cat(ip1,i_m(im)),pri_k(ip1,i_m(im)),kt(ip1,i_m(im)),it(ip1,i_m(im)),jt(ip1,i_m(im)),rx(ip1,i_m(im)),ry(ip1,i_m(im))
      enddo


! Apply Arribas algorithm to partitoned k-list 2 in order specified by place_order sequence

! reset slitlet map for k-list 2

      call make_slitlet_map

      call initialize_slitmap(n_disp,sthresh,1,1)

      nm=0

      do no=1,2*nclass

      ic=place_order(no,1) !priority class set

      jf=place_order(no,2) !k_list portion set 1=com,2=part,3=unique
      if (jf.eq.1) nn=ncom2
      if (jf.eq.2) nn=nu2

      nkc=0
      do ik=1,nn
      if (jf.eq.1) ikk=idc2(ik)
      if (jf.eq.2) ikk=idu2(ik)
      if (pri_k(ip2,ikk).eq.ic) then
      nkc=nkc+1
      ktc(nkc)=kt(ip2,ikk)
      itc(nkc)=it(ip2,ikk)
      jtc(nkc)=jt(ip2,ikk)
      i_nkc(nkc)=ikk !full k-space input catalog reference
      endif
      enddo

      call arribas_m_make_test(nkc,ktc,itc,jtc,nmc,i_mc)

! Gather placed objects for pointing in sequence

      do im=1,nmc
      nm=nm+1
      i_m(nm)=i_nkc(i_mc(im))
      enddo

      enddo !no placement sequence order loop


      write(10,*) nm
!     nm1=nm

      do im=1,nm
       write(10,200) id(ip2,i_m(im)),id_cat(ip2,i_m(im)),pri_k(ip2,i_m(im)),kt(ip2,i_m(im)),it(ip2,i_m(im)),jt(ip2,i_m(im)),rx(ip2,i_m(im)),ry(ip2,i_m(im))
      enddo


      enddo ! npaircount

      close(10)

      write(*,*)

      endif !n_dither=2



!*************************************


       if (n_dither.eq.3) then

       n_possible=0
       do ip1=1,npointings !march through ip1
       do ip2=1,ip1-1 !march through ip2 for each ip1
       do ip3=1,ip2-1 !march through ip3 for each legal pair ip1 and ip2
       n_possible=n_possible+1
       enddo
       enddo
       enddo

! first run through to establish number of legal triples 

       ntriples=0

       do ip1=1,npointings !march through ip1

       do ip2=1,ip1-1 !march through ip2 for each ip1

! check legality of ip2 versus ip1

       dy12=abs(ty_p(ip1)-ty_p(ip2))

       if (dy12.gt.max_dy) goto 77 !ip2 too far away from ip1 in y, next

       dx12=abs(tx_p(ip1)-tx_p(ip2))

       if (dx12.gt.max_dx) goto 77 !ip2 too far away from ip1 in x, next

       if ((dy12.lt.min_dy).and.(dx12.le.min_dx)) goto 77 !ip2 in same r


! ip1 and ip2 mutually legal at this point,

       do ip3=1,ip2-1 !march through ip3 for each legal pair ip1 and ip2

! check legality of ip3 against ip1

       dy13=abs(ty_p(ip1)-ty_p(ip3))

       if (dy13.gt.max_dy) goto 66 !ip3 too far away from ip1 in y, next

       dx13=abs(tx_p(ip1)-tx_p(ip3))

       if (dx13.gt.max_dx) goto 66 !ip3 too far away from ip1 in x, next

       if ((dy13.lt.min_dy).and.(dx13.le.min_dx)) goto 66 !ip3 in same r


! check legality of ip3 against ip2

       dy23=abs(ty_p(ip2)-ty_p(ip3))

       if (dy23.gt.max_dy) goto 66 !ip3 too far away from ip2 in y, next

       dx23=abs(tx_p(ip2)-tx_p(ip3))

       if (dx23.gt.max_dx) goto 66 !ip3 too far away from ip2 in x, next

       if ((dy23.lt.min_dy).and.(dx23.le.min_dx)) goto 66 !ip3 in same r


! ip1, ip2 and ip3 are a mutually legal triple at this point

       ntriples=ntriples+1


66     continue !next ip3
       enddo !ip3

77     continue !next ip2
       enddo !ip2

       enddo !ip1

       write(*,*)
       write(*,110) ' Number of pointings:        ',npointings
       write(*,110) ' Number of possible triples: ',n_possible
       write(*,110) ' Number of legal triples:    ',ntriples
       write(*,110) ' Number of Priority Classes: ',nclass
110    format(a,i0)
       write(*,*)



      if (ntriples.lt.1) then
      write(*,*)
      write(*,*) '*****************************************************'
      write(*,*)
      write(*,*) '          No legal dither triples found!'
      write(*,*) '              Program m_make stopped'
      write(*,*)
      write(*,*) ' Either relax the legal triple proximity rules in'
      write(*,*) ' the configuration file and run m_make again, and/or'
      write(*,*) ' add additional less optimal pointings from the IPA'
      write(*,*) ' output to the configuration file for processing'
      write(*,*) ' and run k_make followed by m_make again'
      write(*,*)
      write(*,*) '*****************************************************'
      write(*,*)
      stop
      endif !no legal triple



! second run through to fill out ipp1, ipp2, 1pp3 lists

       allocate (ipp1(ntriples),ipp2(ntriples),ipp3(ntriples)) 

       ntriples=0

       do ip1=1,npointings !march through ip1

       do ip2=1,ip1-1 !march through ip2 for each ip1

! check legality of ip2 versus ip1

       dy12=abs(ty_p(ip1)-ty_p(ip2))

       if (dy12.gt.max_dy) goto 770 !ip2 too far away from ip1 in y, next

       dx12=abs(tx_p(ip1)-tx_p(ip2))

       if (dx12.gt.max_dx) goto 770 !ip2 too far away from ip1 in x, next

       if ((dy12.lt.min_dy).and.(dx12.le.min_dx)) goto 770 !ip2 in same r


! ip1 and ip2 mutually legal at this point,

       do ip3=1,ip2-1 !march through ip3 for each legal pair ip1 and ip2

! check legality of ip3 against ip1

       dy13=abs(ty_p(ip1)-ty_p(ip3))

       if (dy13.gt.max_dy) goto 660 !ip3 too far away from ip1 in y, next

       dx13=abs(tx_p(ip1)-tx_p(ip3))

       if (dx13.gt.max_dx) goto 660 !ip3 too far away from ip1 in x, next

       if ((dy13.lt.min_dy).and.(dx13.le.min_dx)) goto 660 !ip3 in same r

! check legality of ip3 against ip2

       dy23=abs(ty_p(ip2)-ty_p(ip3))

       if (dy23.gt.max_dy) goto 660 !ip3 too far away from ip2 in y, next

       dx23=abs(tx_p(ip2)-tx_p(ip3))

       if (dx23.gt.max_dx) goto 660 !ip3 too far away from ip2 in x, next

       if ((dy23.lt.min_dy).and.(dx23.le.min_dx)) goto 660 !ip3 in same r


! ip1, ip2 and ip3 are a mutually legal triple at this point

       ntriples=ntriples+1

       ipp1(ntriples)=ip1 !target list sequence numbers
       ipp2(ntriples)=ip2
       ipp3(ntriples)=ip3


660    continue !next ip3
       enddo !ip3

770    continue !next ip2
       enddo !ip2

       enddo !ip1


!      open output m-space target list file

       open(10,file=outfile,status='replace',access='sequential', form='formatted')

       write(10,*) ntriples


! Begin pointing triple loop

      allocate (id1(n_max_k),id2(n_max_k),id3(n_max_k))
      allocate (idc1(n_max_k),idp1(n_max_k),idu1(n_max_k))
      allocate (idc2(n_max_k),idp2(n_max_k),idu2(n_max_k))
      allocate (idc3(n_max_k),idp3(n_max_k),idu3(n_max_k))

      do ntriplecount=1,ntriples

      ip1=ipp1(ntriplecount)
      ip2=ipp2(ntriplecount)
      ip3=ipp3(ntriplecount)


      write(10,*) ntriplecount
      write(10,*) ip1,ip2,ip3   !pointing list sequence numbers
      write(10,*) no_p(ip1),no_p(ip2),no_p(ip3) !IPA Pointing list sequence



      write(*,100) ' Processing triple pointing ',ntriplecount,' of ',ntriples
!     write(*,*) ip1,ip2,ip3
100   format(a,i0,a,i0)

      nk1=nk(ip1)
      nk2=nk(ip2)
      nk3=nk(ip3)

      do ik1=1,nk1
      id1(ik1)=id_cat(ip1,ik1)
      enddo

      do ik2=1,nk2
      id2(ik2)=id_cat(ip2,ik2)
      enddo

      do ik3=1,nk3
      id3(ik3)=id_cat(ip3,ik3)
      enddo


! sort k-list 1 into common (ncom1), partial (np1) and unique parts (nu1)

      ncom1=0
      np1=0
      nu1=0

      do ik1=1,nk1


      l2=.false.
      do ik2=1,nk2
        if (id1(ik1).eq.id2(ik2)) then
         l2=.true.
         goto 666
        endif
       enddo !ik2

666    continue


      l3=.false.
      do ik3=1,nk3
       if (id1(ik1).eq.id3(ik3)) then
        l3=.true.
        goto 667
        endif
       enddo !ik3

 667   continue


        if (l2.and.l3) then
         ncom1=ncom1+1
         idc1(ncom1)=ik1
         goto 668
        endif

        if (l2.or.l3) then
         np1=np1+1
         idp1(np1)=ik1 !sequential number in whole k-list
         goto 668
       endif

        if ((.not.l2).and.(.not.l3)) then
         nu1=nu1+1
         idu1(nu1)=ik1
         goto 668
        endif

668   continue

      enddo !ik1


! sort k-list 2 into common (ncom2), partial (np2) and unique parts (nu2)

      ncom2=0
      np2=0
      nu2=0

       do ik2=1,nk2

      l1=.false.
      do ik1=1,nk1
        if (id2(ik2).eq.id1(ik1)) then
         l1=.true.
         goto 766
        endif
       enddo !ik1

766    continue

      l3=.false.
      do ik3=1,nk3
       if (id2(ik2).eq.id3(ik3)) then
        l3=.true.
        goto 767
        endif
       enddo !ik3

767   continue

        if (l1.and.l3) then
         ncom2=ncom2+1
         idc2(ncom2)=ik2
         goto 768
        endif

        if (l1.or.l3) then
         np2=np2+1
         idp2(np2)=ik2 !sequential number in whole k-list
         goto 768
       endif

        if ((.not.l1).and.(.not.l3)) then
         nu2=nu2+1
         idu2(nu2)=ik2
         goto 768
        endif



768   continue
      enddo !ik2


! sort k-list 3 into common (ncom3), partial (np3) and unique parts (nu3)

      ncom3=0
      np3=0
      nu3=0

       do ik3=1,nk3

      l1=.false.
      do ik1=1,nk1
        if (id3(ik3).eq.id1(ik1)) then
         l1=.true.
         goto 866
        endif
       enddo !ik1

866    continue

      l2=.false.
      do ik2=1,nk2
       if (id3(ik3).eq.id2(ik2)) then
        l2=.true.
        goto 867
        endif
       enddo !ik2

867   continue

        if (l1.and.l2) then
         ncom3=ncom3+1
         idc3(ncom3)=ik3
         goto 868
        endif

        if (l1.or.l2) then
         np3=np3+1
         idp3(np3)=ik3 !sequential number in whole k-list
         goto 868
       endif

        if ((.not.l1).and.(.not.l2)) then
         nu3=nu3+1
         idu3(nu3)=ik3
         goto 868
        endif

868   continue
      enddo !ik3


!     write(*,*) nk1,ncom1,np1,nu1,ncom1+np1+nu1-nk1
!     write(*,*) nk2,ncom2,np2,nu2,ncom2+np2+nu2-nk2
!     write(*,*) nk3,ncom3,np3,nu3,ncom3+np3+nu3-nk3

!     do ik=1,ncom1
!     write(*,*) idc1(ik),idc2(ik)
!     enddo


! Apply Arribas algorithm to partitoned k-list 1 in order specified by place_order sequence

! reset slitlet map for k-list 1

      call make_slitlet_map

      call initialize_slitmap(n_disp,sthresh,1,1)

      nm=0

      do no=1,3*nclass

      ic=place_order(no,1) !priority class set

      jf=place_order(no,2) !k_list portion set 1=com,2=part,3=unique
      if (jf.eq.1) nn=ncom1
      if (jf.eq.2) nn=np1
      if (jf.eq.3) nn=nu1

      nkc=0
      do ik=1,nn
      if (jf.eq.1) ikk=idc1(ik)
      if (jf.eq.2) ikk=idp1(ik)
      if (jf.eq.3) ikk=idu1(ik)
      if (pri_k(ip1,ikk).eq.ic) then
      nkc=nkc+1
      ktc(nkc)=kt(ip1,ikk)
      itc(nkc)=it(ip1,ikk)
      jtc(nkc)=jt(ip1,ikk)
      i_nkc(nkc)=ikk !full k-space input catalog reference
      endif
      enddo

      call arribas_m_make_test(nkc,ktc,itc,jtc,nmc,i_mc)

! Gather placed objects for pointing in sequence

      do im=1,nmc
      nm=nm+1
      i_m(nm)=i_nkc(i_mc(im))
      enddo

      enddo !no placement sequence order loop


      write(10,*) nm
!     nm1=nm

      do im=1,nm
       write(10,200) id(ip1,i_m(im)),id_cat(ip1,i_m(im)),pri_k(ip1,i_m(im)),kt(ip1,i_m(im)),it(ip1,i_m(im)),jt(ip1,i_m(im)),rx(ip1,i_m(im)),ry(ip1,i_m(im))
200   format(I7,2x,I13,x,I3,3x,I1,x,I3,x,I3,x,f6.3,x,f6.3)
      enddo


! Apply Arribas algorithm to partitoned k-list 2 in order specified by place_order sequence

! reset slitlet map for k-list 2

      call make_slitlet_map

      call initialize_slitmap(n_disp,sthresh,1,1)

      nm=0

      do no=1,3*nclass

      ic=place_order(no,1) !priority class set

      jf=place_order(no,2) !k_list portion set 1=com,2=part,3=unique
      if (jf.eq.1) nn=ncom2
      if (jf.eq.2) nn=np2
      if (jf.eq.3) nn=nu2

      nkc=0
      do ik=1,nn
      if (jf.eq.1) ikk=idc2(ik)
      if (jf.eq.2) ikk=idp2(ik)
      if (jf.eq.3) ikk=idu2(ik)
      if (pri_k(ip2,ikk).eq.ic) then
      nkc=nkc+1
      ktc(nkc)=kt(ip2,ikk)
      itc(nkc)=it(ip2,ikk)
      jtc(nkc)=jt(ip2,ikk)
      i_nkc(nkc)=ikk !full k-space input catalog reference
      endif
      enddo

      call arribas_m_make_test(nkc,ktc,itc,jtc,nmc,i_mc)

! Gather placed objects for pointing in sequence

      do im=1,nmc
      nm=nm+1
      i_m(nm)=i_nkc(i_mc(im))
      enddo

      enddo !no placement sequence order loop


      write(10,*) nm
!     nm1=nm

      do im=1,nm
       write(10,200) id(ip2,i_m(im)),id_cat(ip2,i_m(im)),pri_k(ip2,i_m(im)),kt(ip2,i_m(im)),it(ip2,i_m(im)),jt(ip2,i_m(im)),rx(ip2,i_m(im)),ry(ip2,i_m(im))
      enddo


! Apply Arribas algorithm to partitoned k-list 3 in order specified by place_order sequence

! reset slitlet map for k-list 3

      call make_slitlet_map

      call initialize_slitmap(n_disp,sthresh,1,1)

      nm=0

      do no=1,3*nclass

      ic=place_order(no,1) !priority class set

      jf=place_order(no,2) !k_list portion set 1=com,2=part,3=unique
      if (jf.eq.1) nn=ncom3
      if (jf.eq.2) nn=np3
      if (jf.eq.3) nn=nu3

      nkc=0
      do ik=1,nn
      if (jf.eq.1) ikk=idc3(ik)
      if (jf.eq.2) ikk=idp3(ik)
      if (jf.eq.3) ikk=idu3(ik)
      if (pri_k(ip3,ikk).eq.ic) then
      nkc=nkc+1
      ktc(nkc)=kt(ip3,ikk)
      itc(nkc)=it(ip3,ikk)
      jtc(nkc)=jt(ip3,ikk)
      i_nkc(nkc)=ikk !full k-space input catalog reference
      endif
      enddo

      call arribas_m_make_test(nkc,ktc,itc,jtc,nmc,i_mc)

! Gather placed objects for pointing in sequence

      do im=1,nmc
      nm=nm+1
      i_m(nm)=i_nkc(i_mc(im))
      enddo

      enddo !no placement sequence order loop


      write(10,*) nm
!     nm1=nm

      do im=1,nm
       write(10,200) id(ip3,i_m(im)),id_cat(ip3,i_m(im)),pri_k(ip3,i_m(im)),kt(ip3,i_m(im)),it(ip3,i_m(im)),jt(ip3,i_m(im)),rx(ip3,i_m(im)),ry(ip3,i_m(im))
      enddo


      enddo ! ntriplecount

      close(10)

      write(*,*)

      endif !n_dither=3

      open(11,file=confname,status='old',access='sequential',form='formatted')

      call find_marker(11,'#MM#')

      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- SECTION AUTOMATICALLY APPENDED BY THE M_MAKE MODULE ----------'
      write(11,'(a)') '#'

      write(11,'(a)') '# The generated m-list file is located in:'
      write(11,'(a)') '#'
      write(11,'(a)') '#   ./'//outfile
      write(11,'(a)') '#'

       write(11,'(a,i0)') '# Number of pointings:        ',npointings

       if (n_dither.eq.2) then
       write(11,'(a,i0)') '# Number of possible pairs:   ',n_possible
       write(11,'(a,i0)') '# Number of legal pairs:      ',npairs
       endif

       if (n_dither.eq.3) then
       write(11,'(a,i0)') '# Number of possible triples: ',n_possible
       write(11,'(a,i0)') '# Number of legal triples:    ',ntriples
       endif

      write(11,'(a)') '#'



      write(11,'(a)')  '# The following default parameters should be reviewed and modified as needed before running the m_sort module:'



      write(11,'(a)') '#'


      write(11,'(a)') '# Figure of Merit Weights for each Priority Class in increasing order'
       write(11,*) (0.5**(i_c-1),i_c=1,nclass)
!      write(11,*) (0.002**(i_c-1),i_c=1,nclass)
!       write(11,*) 1.,(0.,i_c=1,nclass-1)
!      write(11,*) 1.0,0.670475,0.294855,0.149834,0.082171,0.048325,0.013295,0.000961
      
      
    
      


      write(11,'(a)') '#'


      write(11,'(a)') '# ----- END OF SECTION APPENDED BY THE M_MAKE MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a)') '#  -- Modifying any parameter in the above section requires re-running the m_sort and following modules ---'
      write(11,'(a)') '#'
      write(11,'(a)') '#MS# -- Marker do not delete or move this line'
 
      close(11)



       write(*,*)
       call cpu_time(t_stop)
       write(*,'(a)')      '------------------------------------------------'
       write(*,'(a,f7.1)') ' CPU time used by m_make module [s]:',t_stop-t_start
       write(*,'(a)')      '------------------------------------------------'

       write(*,*)
         
      end program m_make
      
!-------------------------------------------------------

 
include './reference_files/front_transforms_updt.f90'

include './reference_files/back_transforms_updt.f90'

include './reference_files/shutter_routines_new.f90'

include './reference_files/misc_routines_new.f90'

!------------------------------------------------------------

      subroutine file_entry_m_make()
      use derived_parameters, only : phi
      use config_parameters
      use order_and_weights
      use target_cat, only : nclass
      implicit none
      integer num_arg,iargc,i,j,io
      integer checksum1,checksum2

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
      read(9,*) !raccx,raccy

      call skip_hashes(9)
      read(9,*) sthresh


      call skip_hashes(9)
      read(9,*) !catfile

      call skip_hashes(9)
      read(9,*) !ref image name
      
      call skip_hashes(9)
      read(9,*) !seg map name
      
      call skip_hashes(9)
      read(9,*) !cra,cdec,cpa_ap

!       if ((cpa_ap.lt.0.).or.(cpa_ap.gt.360.)) then
!       write(*,*)
!       write(*,*) 'Error: Nominal Roll Angle out of range in configuration file'
!       write(*,*)
!       stop
!       endif
! 
!       cpa_v3=cpa_ap-(180.-phi)
!       if (cpa_v3.lt.0.) cpa_v3=cpa_v3+360.
!       if (cpa_v3.gt.360.) cpa_v3=cpa_v3-360.


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

      allocate (order_matrix(nclass,n_dither),place_order(nclass*n_dither,2))
 
      call skip_hashes(9)

      if (n_dither.eq.1) then
      read(9,*) (order_matrix(i,1), i=1,nclass)

      do io=1,nclass
      do i=1,nclass
      if (order_matrix(i,1).eq.io) then
      place_order(io,1)=i
      goto 997
      endif
      enddo !i
997   continue
      enddo !io

      endif !n_dither=1


      if (n_dither.eq.2) then
      read(9,*) (order_matrix(i,1), i=1,nclass)
      call skip_hashes(9)
      read(9,*) (order_matrix(i,2), i=1,nclass)

      do io=1,2*nclass
      do i=1,nclass
      do j=1,2
      if (order_matrix(i,j).eq.io) then
      place_order(io,1)=i
      place_order(io,2)=j
      goto 999
      endif
      enddo !j
      enddo !i
999   continue
      enddo !io

      checksum1=0
      do i=1,2*nclass
      checksum1=checksum1+i
      enddo

      checksum2=0
      do i=1,nclass
      do j=1,2
      checksum2=checksum2+order_matrix(i,j)
      enddo
      enddo

      if (checksum2.ne.checksum1) then
      write(*,*)
      write(*,*) 'Error in specification of placement order matrix in configuration file!'
      write(*,*)
      stop
      endif

      endif !n_dither=2


      if (n_dither.eq.3) then
      read(9,*) (order_matrix(i,1), i=1,nclass)
      call skip_hashes(9)
      read(9,*) (order_matrix(i,2), i=1,nclass)
      call skip_hashes(9)
      read(9,*) (order_matrix(i,3), i=1,nclass)

      do io=1,3*nclass
      do i=1,nclass
      do j=1,3
      if (order_matrix(i,j).eq.io) then
      place_order(io,1)=i
      place_order(io,2)=j
      goto 777
      endif
      enddo !j
      enddo !i
777   continue
      enddo !io


      checksum1=0
      do i=1,3*nclass
      checksum1=checksum1+i
      enddo

      checksum2=0
      do i=1,nclass
      do j=1,3
      checksum2=checksum2+order_matrix(i,j)
      enddo
      enddo

      if (checksum2.ne.checksum1) then
      write(*,*)
      write(*,*) 'Error: Incorrect specification of placement order matrix in configuration file!'
      write(*,*)
      stop
      endif

      endif !n_dither=3



      if (n_dither.gt.1) then
      call skip_hashes(9)
      read(9,*) min_dx,max_dx,min_dy,max_dy
      endif


 
      close(9)

      return

      end subroutine file_entry_m_make

!------------------------------------------------------------------------
 
      subroutine arribas_m_make_test(nk,kt,it,jt,nm,i_m)

! determine nm non-overlapping targets from nk input
! This m_make version updates the viable slitlet map
! NB: shuttervalues already loaded in call to initialize_slitmap

      use config_parameters, only : n_disp,sthresh
      use derived_parameters, only : avxgap
      use k_list, only : n_max_k
      implicit none
      integer nk
      integer kt(nk),it(nk),jt(nk)
      integer nm,i_m(n_max_k)
      integer ovlap(nk,nk),tovlap(nk),max_targ(nk)
      real shv_i,shv_j
      integer n_purge,max_tovlap,n_max,idel,idum,ntout
      integer k,i,j,k2,j1,j2,m,ii,jj,nk_p,dj
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

      nk_p=0
      do i=1,nk
      if (slitmap(kt(i),it(i),jt(i)).ne.4) then
       nk_p=nk_p+1
        do ii=1,nk
         ovlap(ii,i)=0
         ovlap(i,ii)=0
        enddo
      endif
      enddo

      if (nk_p-nk.eq.0) then
      nm=0
!     write(*,*)
!     write(*,*) 'No targets in non-occupied slitlets'
!     write(*,*)
!     return !No further processing needed
      return
      endif


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

! Fill in reshuffled overlap matrix

       ovlap=0


      do i=1,nm
      do j=1,nm
      ovlap(i,j)=0
      k=kt(i_m(i))
      if (k.eq.1) k2=3
      if (k.eq.3) k2=1
      if (k.eq.2) k2=4
      if (k.eq.4) k2=2
      if (.not.((kt(i_m(j)).eq.k).or.(kt(i_m(j)).eq.k2))) goto 7777
      shv_i=shval(kt(i_m(i)),it(i_m(i)),jt(i_m(i)))
      shv_j=shval(kt(i_m(j)),it(i_m(j)),jt(i_m(j)))
      if (abs(shv_i-shv_j).lt.sthresh) ovlap(i,j)=1
      if (n_disp.eq.7) then
         sep_p=float(prism_sep_p(kt(i_m(i)),it(i_m(i)),jt(i_m(i))))
         sep_m=float(prism_sep_m(kt(i_m(i)),it(i_m(i)),jt(i_m(i))))
       i1=float(it(i_m(i)))+msa_gap(kt(i_m(i)))*(avxgap+365.)
       i2=float(it(i_m(j)))+msa_gap(kt(i_m(j)))*(avxgap+365.)
       if (((i2-i1).le.sep_m).or.((i2-i1).ge.sep_p)) ovlap(i,j)=0
      endif
7777  continue
      enddo
      enddo


       do i=1,nm
       tovlap(i)=0
       do j=1,nk
       tovlap(i)=tovlap(i)+ovlap(i,j)
       enddo
       enddo


!     write(*,*) 'Non-overlapping prism targets in non-occupied slitlets:'
!     do i=1,nm
!      write(*,*) i,tovlap(i_m(i)),i_m(i),kt(i_m(i)),it(i_m(i)),jt(i_m(i))
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


      if (abs(shval(k,i,j)-shval(k2,ii,jj)).le.sthresh) then !repeat for

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
      end subroutine arribas_m_make_test

!------------------------------------------------------------------------
 
