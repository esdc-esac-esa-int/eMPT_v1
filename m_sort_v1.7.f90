! V1.7 fixed sky shutter flag write to config file inconsistency  


      include 'empt_modules.f90'

!-------------------------------------------------------

      program m_sort
      use path_to_ref_files
      use config_parameters
      use target_cat, only : nclass
      use pointings_and_groups, only : n_p
      use m_list
      use order_and_weights, only : weights
      implicit none
      real t_start,t_stop
      character*2 rn
      character*70 filename
      integer, allocatable :: nsingles,npairs,ntriples
      integer, allocatable :: nsinglecount,npaircount,ntriplecount
      integer, allocatable :: ip1,ip2,ip3
      integer, allocatable :: np1,np2,np3
      integer, allocatable :: im1,im2,im3
      integer, allocatable :: nm1,nm2,nm3
      integer, allocatable :: id1(:),id2(:),id3(:)
      integer, allocatable :: cl1(:),cl2(:),cl3(:)
      integer, allocatable :: ipp1(:),ipp2(:),ipp3(:)
      integer, allocatable :: npp1(:),npp2(:),npp3(:)
      integer, allocatable :: indx(:)
      integer, allocatable :: n_pair_best,n_triple_best
      real, allocatable :: m_fom1(:),m_fom2(:),m_fom3(:)
      real, allocatable :: m_fom(:)
      real, allocatable :: fom1,fom2,fom3
      real m_fom_best
      integer nstop
      integer  i,idum

      call cpu_time(t_start)


      write(*,*)
      write(*,*) '-----------------------------------------------------------------------------------------'
      write(*,*) '--- m_sort -- eMPT General Dithered m-list Sorting Module Version 1.7 26 September 2022 -'
      write(*,*) '-----------------------------------------------------------------------------------------'
      write(*,*)

      call file_entry_m_sort()
      
      write(rn,'(i2)') ntile
      rn=trim(adjustl(rn))
      if (ntile.lt.10) rn='0'//rn(1:1)

      call system('mkdir trial_'//rn(1:2)//'/m_sort_output')
      
      if (n_dither.eq.1) m_list_file='trial_'//rn(1:2)//'/m_make_output/single_m_list.txt'
      if (n_dither.eq.2) m_list_file='trial_'//rn(1:2)//'/m_make_output/pair_m_list.txt'
      if (n_dither.eq.3) m_list_file='trial_'//rn(1:2)//'/m_make_output/triple_m_list.txt'
   
      call m_list_size_read()
      
       write(*,*)
       write(*,'(a,a2)') ' Trial Identifier No.:        ',rn
       write(*,'(a,a)')  ' Configuration File:          ',confname
       write(*,'(a,i0)') ' n_dither:                    ',n_dither
       write(*,'(a,i0)') ' Number of Priority Classes in input catalog:      ',nclass
       write(*,'(a,i0)') ' Total number of pointings in configuration file:  ',n_p
       write(*,'(a,a)')  ' m_list file:                                      ',m_list_file
       write(*,'(a,i0)') ' Maximum number of targets in m_list per pointing: ',n_max_m
       write(*,*)


!***************************************************


       if (n_dither.eq.1) then

!      open output m-space target list file

       open(10,file=m_list_file,status='old',access='sequential', form='formatted')

      allocate (nsingles,nsinglecount)

      read(10,*) nsingles

       write(*,*)
       write(*,'(a,i0)') ' Number of single pointings:       ',nsingles

      allocate (ip1,np1,im1,nm1,fom1)
      allocate (ipp1(nsingles),npp1(nsingles))
      allocate (id1(n_max_m),cl1(n_max_m))
      allocate (m_fom1(nsingles))
      allocate (m_fom(nsingles))
      allocate (indx(nsingles))

      write(*,*)
      write(*,*)

! Read and process each m-space single pointings in turn, saving pointing combinations for k-space calculation

      do i=1,nsingles

      read(10,*) nsinglecount
      read(10,*) ip1

      ipp1(nsinglecount)=ip1

      read(10,*) np1

      npp1(nsinglecount)=np1

      write(*,114) ' Processing m-space single pointing ',nsinglecount,' of ',nsingles
114   format(a,i0,a,i0)

      read(10,*) nm1
      do im1=1,nm1
      read(10,*) idum,id1(im1),cl1(im1)
      enddo!im1

       call single_sort(n_max_m,nm1,id1,cl1,fom1)

       m_fom1(nsinglecount)=fom1

      enddo!i





      do i=1,nsingles
      m_fom(i)=(m_fom1(i))
      enddo!i


      call indexx(nsingles,m_fom,indx)

      write(*,*)
      write(*,*) 'FOM2: Total weighted number of targets observed'
      write(*,*) 'Highest FOM2 single pointings:'
      write(*,*) '          Pointing   FOM2'
      nstop=max(1,nsingles-9)
      do i=nsingles,nstop,-1
      write(*,2001) indx(i),npp1(indx(i)),m_fom(indx(i))
      enddo
2001  format(3x,i4,3x,i4,3x,f9.5)

      write(*,*)

      write(*,*) 'Lowest FOM2 single pointings:'
      write(*,*) '          Pointing   FOM2'
      nstop=min(nsingles,10)
      do i=1,nstop
      write(*,2001) indx(i),npp1(indx(i)),m_fom(indx(i))
      enddo

! write ranked pointings out to file

      filename='trial_'//rn(1:2)//'/m_sort_output/single_list_fom2.txt'

       open(10,file=filename,status='replace',access='sequential', form='formatted')
      write(10,'(a)') '#         Pointing  FOM2'
      do i=nsingles,1,-1
      write(10,2001) indx(i),npp1(indx(i)),m_fom(indx(i))
      enddo
      close(10)

      write(*,*)

      filename= 'trial_'//rn(1:2)//'/m_sort_output/single_plot_fom2.ps/cps'



! append configuration file

      open(11,file=confname,status='old',access='sequential', form='formatted')

      call find_marker(11,'#MS#')

      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- SECTION AUTOMATICALLY APPENDED BY THE M_SORT MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a)') '# The figure-of-merit ranked single pointing list is located in:'
      write(11,'(a)') '#'
      write(11,'(a)') '#   ./trial_'//rn(1:2)//'/m_sort_output/single_list_fom2.txt'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# The following default parameters should be reviewed and modified as needed before running the m_pick module:'
      write(11,'(a)') '#'
      write(11,'(a)') '# Top ranked (FOM2) pointings:'
      write(11,'(a)') '#         Pointing   FOM2'
      i=nsingles
      write(11,2001) indx(i),npp1(indx(i)),m_fom(indx(i))
      nstop=max(1,nsingles-9)
      do i=nsingles-1,nstop,-1
      write(11,'(a1,3x,i4,3x,i4,3x,f9.5)') '#',indx(i),npp1(indx(i)),m_fom(indx(i))
      enddo
      write(11,'(a)') '#'
      write(11,'(a)') '# Automatically fill out final MSA mask with Sky Background Slitlets?'
      write(11,'(a)') '# [Y or N [Default] ]'
      write(11,'(a)') '  N'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- END OF SECTION APPENDED BY THE M_SORT MODULE ----------'

      write(11,'(a)') '#'
      write(11,'(a)') '#  -- Modifying any parameter in the above section requires re-running the m_pick module ---'
      write(11,'(a)') '#'
      write(11,'(a)') '#MP# -- Marker do not delete or move this line'

       close(11)


      endif !n_dither=1


!***************************************************


       if (n_dither.eq.2) then

!      open output m-space target list file

       allocate (npairs,npaircount)

       open(10,file=m_list_file,status='old',access='sequential', form='formatted')

      read(10,*) npairs


       write(*,*)
       write(*,'(a,i0)') ' Number of legal pairs:     ',npairs

      allocate (ip1,ip2,np1,np2,im1,im2,nm1,nm2,fom1,fom2)
      allocate (ipp1(npairs),npp1(npairs),ipp2(npairs),npp2(npairs))
      allocate (id1(n_max_m),cl1(n_max_m),id2(n_max_m),cl2(n_max_m))
      allocate (m_fom1(npairs),m_fom2(npairs))
      allocate (m_fom(npairs))
      allocate (n_pair_best)
      allocate (indx(npairs))

      write(*,*)
      write(*,*)

! Read and process each m-space triple pointings in turn, saving pointing combinations for k-space calculation

      do i=1,npairs

      read(10,*) npaircount
      read(10,*) ip1,ip2

      ipp1(npaircount)=ip1
      ipp2(npaircount)=ip2

      read(10,*) np1,np2

      npp1(npaircount)=np1
      npp2(npaircount)=np2

      write(*,112) ' Processing m-space pointing pair ',npaircount,' of ',npairs
112   format(a,i0,a,i0)

      read(10,*) nm1
      do im1=1,nm1
      read(10,*) idum,id1(im1),cl1(im1)
      enddo!im1

      read(10,*) nm2
      do im2=1,nm2
      read(10,*) idum,id2(im2),cl2(im2)
      enddo!im2


       call pair_sort(n_max_m,nm1,id1,cl1,nm2,id2,cl2,fom2,fom1)

       m_fom2(npaircount)=fom2
       m_fom1(npaircount)=fom1


      enddo!i



      m_fom_best=-10.
      do i=1,npairs
      m_fom(i)=(2.*m_fom2(i)+1.*m_fom1(i))/(m_fom2(i)+m_fom1(i))
     
      if (m_fom(i).gt.m_fom_best) then
      m_fom_best=m_fom(i)
      n_pair_best=i
      endif

      enddo!i


      call indexx(npairs,m_fom,indx)

      write(*,*)
      write(*,*) 'FOM1: Average exposure per weighted target'
      write(*,*) 'Highest FOM1 paired pointings:'
      write(*,*) ' Pair    P2  P1     FOM1'
      nstop=max(1,npairs-9)
      do i=npairs,nstop,-1
      write(*,1007) indx(i),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      enddo
1007  format(i5,3x,2i4,3x,f9.5)
      write(*,*)

      write(*,*) 'Lowest FOM1 paired pointings:'
      write(*,*) ' Pair    P2  P1     FOM1'
      nstop=min(npairs,10)
      do i=1,nstop
      write(*,1007) indx(i),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      enddo

      write(*,*)

! write ranked pairs out to file for fom1

      filename='trial_'//rn(1:2)//'/m_sort_output/pair_list_fom1.txt'

       open(10,file=filename,status='replace',access='sequential', form='formatted')
      write(10,'(a)') '# Pair    P2  P1     FOM1'
      do i=npairs,1,-1
      write(10,1007) indx(i),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      enddo
      close(10)

      write(*,*)

! append configuration file

      open(11,file=confname,status='old',access='sequential', form='formatted')

      call find_marker(11,'#MS#')

      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- SECTION AUTOMATICALLY APPENDED BY THE M_SORT MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a)') '# The figure-of-merit ranked pointing pair lists are located in:'
      write(11,'(a)') '#'
      write(11,'(a)') '#   ./trial_'//rn(1:2)//'/m_sort_output/pair_list_fom1.txt'
      write(11,'(a)') '#   ./trial_'//rn(1:2)//'/m_sort_output/pair_list_fom2.txt'
      write(11,'(a)') '#'
      write(11,'(a)') '# Top ranked (FOM1) pointing pairs:'
      write(11,'(a)') '# Pair    P2  P1     FOM1'
      i=npairs
      write(11,1007) indx(i),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      nstop=max(1,npairs-9)
      do i=npairs-1,nstop,-1
      write(11,'(a1,i5,3x,2i4,3x,f9.5)') '#',indx(i),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      enddo
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# Automatically fill out final MSA mask with Sky Background Slitlets?'
      write(11,'(a)') '# [Y or N [Default] ]'
      write(11,'(a)') '  N'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- END OF SECTION APPENDED BY THE M_SORT MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a)') '#  -- Modifying any parameter in the above section requires re-running the m_pick module ---'
      write(11,'(a)') '#'
      write(11,'(a)') '#MP# -- Marker do not delete or move this line'

       close(11)



! repeat fom calculation for fom2


      do i=1,npairs
      m_fom(i)=(m_fom2(i)+m_fom1(i))
      enddo!i


      call indexx(npairs,m_fom,indx)

      write(*,*)
      write(*,*) 'FOM2: Total weighted number of targets observed'
      write(*,*) 'Highest FOM2 paired pointings:'
      write(*,*) ' Pair    P2  P1     FOM2'
      nstop=max(1,npairs-9)
      do i=npairs,nstop,-1
      write(*,2007) indx(i),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      enddo
2007  format(i5,3x,2i4,3x,f9.5)

      write(*,*)


      write(*,*) 'Lowest FOM2 paired pointings:'
      write(*,*) ' Pair    P2  P1     FOM2'
      nstop=min(npairs,10)
      do i=1,nstop
      write(*,2007) indx(i),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      enddo

      write(*,*)

      filename='trial_'//rn(1:2)//'/m_sort_output/pair_list_fom2.txt'

       open(10,file=filename,status='replace',access='sequential', form='formatted')
      write(10,'(a)') '# Pair    P2  P1     FOM2'
      do i=npairs,1,-1
      write(10,2007) indx(i),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      enddo
      close(10)

      write(*,*)


      endif !n_dither=2


!***************************************************


       if (n_dither.eq.3) then


!      open output m-space target list file

       open(10,file=m_list_file,status='old',access='sequential', form='formatted')

      allocate(ntriples,ntriplecount)

      read(10,*) ntriples


       write(*,*)
       write(*,'(a,i4)') ' Number of legal triples:   ',ntriples

      allocate (ip1,ip2,ip3,np1,np2,np3,im1,im2,im3,nm1,nm2,nm3,fom1,fom2,fom3)
      allocate (ipp1(ntriples),npp1(ntriples),ipp2(ntriples),npp2(ntriples),ipp3(ntriples),npp3(ntriples))
      allocate (id1(n_max_m),cl1(n_max_m),id2(n_max_m),cl2(n_max_m),id3(n_max_m),cl3(n_max_m))
      allocate (m_fom1(ntriples),m_fom2(ntriples),m_fom3(ntriples))
      allocate (m_fom(ntriples))
      allocate (n_triple_best)
      allocate (indx(ntriples))

      write(*,*)
      write(*,*)

! Read and process each m-space triple pointings in turn, saving pointing combinations for k-space calculation

      do i=1,ntriples

      read(10,*) ntriplecount
      read(10,*) ip1,ip2,ip3

      ipp1(ntriplecount)=ip1
      ipp2(ntriplecount)=ip2
      ipp3(ntriplecount)=ip3

      read(10,*) np1,np2,np3

      npp1(ntriplecount)=np1
      npp2(ntriplecount)=np2
      npp3(ntriplecount)=np3

      write(*,110) ' Processing m-space triple pointing ',ntriplecount,' of ',ntriples
110   format(a,i0,a,i0)


      read(10,*) nm1
      do im1=1,nm1
      read(10,*) idum,id1(im1),cl1(im1)
      enddo!im1

      read(10,*) nm2
      do im2=1,nm2
      read(10,*) idum,id2(im2),cl2(im2)
      enddo!im2

      read(10,*) nm3
      do im3=1,nm3
      read(10,*) idum,id3(im3),cl3(im3)
      enddo!im3

       call triple_sort(n_max_m,nm1,id1,cl1,nm2,id2,cl2,nm3,id3,cl3,fom3,fom2,fom1)

       m_fom3(ntriplecount)=fom3
       m_fom2(ntriplecount)=fom2
       m_fom1(ntriplecount)=fom1

      enddo!i




      m_fom_best=-10.
      do i=1,ntriples
      m_fom(i)=(3.*m_fom3(i)+2.*m_fom2(i)+1.*m_fom1(i))/(m_fom3(i)+m_fom2(i)+m_fom1(i))

      if (m_fom(i).gt.m_fom_best) then
      m_fom_best=m_fom(i)
      n_triple_best=i
      endif

      enddo!i


      call indexx(ntriples,m_fom,indx)
      
      
      write(*,*)
      write(*,*) 'FOM1: Average exposure per weighted target'
      write(*,*) 'Highest FOM1 triple pointings:'
      write(*,*) 'Triple   P3  P2  P1     FOM1'
      nstop=max(1,ntriples-9)
      do i=ntriples,nstop,-1
!     write(*,1000) indx(i),(n_triple(indx(i),j),j=3,1,-1),
!    #m_fom(indx(i))
      write(*,1001) indx(i),npp3(indx(i)),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      enddo
1001  format(i5,3x,3i4,3x,f9.5)
      write(*,*)

      write(*,*) 'Lowest FOM1 triple pointings:'
      write(*,*) 'Triple   P3  P2  P1     FOM1'
      nstop=min(ntriples,10)
      do i=1,nstop
      write(*,1001) indx(i),npp3(indx(i)),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      enddo

      write(*,*)

      filename='trial_'//rn(1:2)//'/m_sort_output/triple_list_fom1.txt'

       open(10,file=filename,status='replace',access='sequential', form='formatted')
      write(10,2004) '# Triple  P3  P2  P1     FOM1'
      do i=ntriples,1,-1
      write(10,1001) indx(i),npp3(indx(i)),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      enddo
2004  format(a29)
      close(10)
      write(*,*)

      filename='trial_'//rn(1:2)//'/m_sort_output/triple_plot_fom1.ps/cps'



! append configuration file

      open(11,file=confname,status='old',access='sequential', form='formatted')

      call find_marker(11,'#MS#')

      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- SECTION AUTOMATICALLY APPENDED BY THE M-SORT MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a)') '# The figure-of-merit ranked pointing triple lists are located in:'
      write(11,'(a)') '#'
      write(11,'(a)') '#   ./trial_'//rn(1:2)//'/m_sort_output/triple_list_fom1.txt'
      write(11,'(a)') '#   ./trial_'//rn(1:2)//'/m_sort_output/triple_list_fom2.txt'
      write(11,'(a)') '#'
      write(11,'(a)') '# Top ranked (FOM1) triple pointings:'
      write(11,'(a)') '# Triple  P3  P2  P1     FOM1'
      i=ntriples
      write(11,1001) indx(i),npp3(indx(i)),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      nstop=max(1,ntriples-9)
      do i=ntriples-1,nstop,-1
      write(11,'(a1,i5,3x,3i4,3x,f9.5)') '#',indx(i),npp3(indx(i)),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      enddo
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# Automatically fill out final MSA mask with Sky Background Slitlets?'
      write(11,'(a)') '# [Y or N [Default] ]'
      write(11,'(a)') '  N'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- END OF SECTION APPENDED BY THE M_SORT MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a)') '#  -- Modifying any parameter in the above section requires re-running the m_pick module ---'
      write(11,'(a)') '#'
      write(11,'(a)') '#MP# -- Marker do not delete or move this line'


       close(11)

! repeat exercise for fom2


      do i=1,ntriples
      m_fom(i)=(m_fom3(i)+m_fom2(i)+m_fom1(i))
      enddo!i


      call indexx(ntriples,m_fom,indx)

      write(*,*)
      write(*,*) 'FOM2: Total weighted number of targets observed'
      write(*,*) 'Highest FOM2 triple pointings:'
      write(*,*) 'Triple   P3  P2  P1     FOM2'
      nstop=max(1,ntriples-9)
      do i=ntriples,nstop,-1
      write(*,2000) indx(i),npp3(indx(i)),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      enddo
2000  format(i5,3x,3i4,3x,f9.5)

      write(*,*)


      write(*,*) 'Lowest FOM2 triple pointings:'
      write(*,*) 'Triple   P3  P2  P1     FOM2'
      nstop=min(ntriples,10)
      do i=1,nstop
      write(*,2000) indx(i),npp3(indx(i)),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      enddo

      write(*,*)
      filename='trial_'//rn(1:2)//'/m_sort_output/triple_list_fom2.txt'

       open(10,file=filename,status='replace',access='sequential', form='formatted')
      write(10,2004) '# Triple  P3  P2  P1     FOM2'
      do i=ntriples,1,-1
      write(10,2000) indx(i),npp3(indx(i)),npp2(indx(i)),npp1(indx(i)),m_fom(indx(i))
      enddo
      close(10)
      write(*,*)


      endif !n_dither=3



      write(*,*)
      call cpu_time(t_stop)
      write(*,'(a)')      '----------------------------------------------'
      write(*,'(a,f7.1)') ' CPU time used by m_sort module [s]:',t_stop-t_start
      write(*,'(a)')      '----------------------------------------------'
      write(*,*)

         
      end program m_sort
      
!-------------------------------------------------------
 
 
include './reference_files/front_transforms_updt.f90'

include './reference_files/back_transforms_updt.f90'

include './reference_files/shutter_routines_new.f90'

include './reference_files/misc_routines_new.f90'

!------------------------------------------------------------

      subroutine file_entry_m_sort()
      use config_parameters
      use order_and_weights, only : weights
      use target_cat, only : nclass
      implicit none
      integer num_arg,iargc,i

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
      read(9,*) !catfile

      call skip_hashes(9)
      read(9,*) !ref image name
      
      call skip_hashes(9)
      read(9,*) !seg map name
      
      call skip_hashes(9)
      read(9,*) !cra,cdec,cpa_ap

      call skip_hashes(9)
      read(9,*) !szone_x,szone_y

      call skip_hashes(9)
      read(9,*) !angle_to_target

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

      allocate (weights(nclass))

      call skip_hashes(9)
      read(9,*) (weights(i),i=1,nclass)

      close(9)

      return

      end subroutine file_entry_m_sort

!------------------------------------------------------------

      subroutine single_sort(n_max,n1,id1,cl1,fom1)
      use order_and_weights, only : weights
      use target_cat, only : nclass
      implicit none
      integer n_max
      integer id1(n_max)
      integer cl1(n_max)
      integer n1
      integer i1
      real fom1
      integer nt,it,idt(n_max),clt(n_max)
      integer m
      logical l1(n_max)
      integer hist1(nclass)
      integer ncheck

! construct combined target list

      nt=0

      do m=1,nclass

! grab class m targets in list 1 and add to total list

      do i1=1,n1
        if (cl1(i1).eq.m) then
        nt=nt+1
        idt(nt)=id1(i1)
        clt(nt)=cl1(i1)
        endif
      enddo !i1

      enddo !m


!     write(*,*) 'combined list',nt
!     do it=1,nt
!     write(*,*) it,idt(it),clt(it)
!     enddo

! map out which targets are in which lists

      do it=1,nt

      l1(it)=.false.


        do i1=1,n1
         if (idt(it).eq.id1(i1)) l1(it)=.true.
        enddo! n1


      enddo! it

! single cover histograms

      do m=1,nclass

       hist1(m)=0

       do it=1,nt

       if (clt(it).eq.m) then
       hist1(m)=hist1(m)+1
       endif !m

       enddo !it

      enddo !m

! calculate checksum

      ncheck=0

      do m=1,nclass
      ncheck=ncheck+hist1(m)
      enddo!m

      ncheck=ncheck-nt

      if (ncheck.ne.0) stop 'Error in  checksum'

! calculate figure of merits

      fom1=0

      do m=1,nclass
      fom1=fom1+weights(m)*float(hist1(m))
      enddo! m

      return

      end subroutine single_sort

!-------------------------------------------------------------

      subroutine pair_sort(n_max,n1,id1,cl1,n2,id2,cl2,fom2,fom1)
      use order_and_weights, only : weights
      use target_cat, only : nclass
      implicit none
      integer n_max
      integer id1(n_max),id2(n_max)
      integer cl1(n_max),cl2(n_max)
      integer n1,n2
      integer i1,i2
      real fom2,fom1
      integer nt,it,idt(2*n_max),clt(2*n_max)
      integer m,ne
      logical l1(2*n_max),l2(2*n_max)
      integer hist1(nclass),hist2(nclass)
      integer ncheck

! construct combined target list

      nt=0

      do m=1,nclass

! grab class m targets in list 1 and add to total list

      do i1=1,n1
        if (cl1(i1).eq.m) then
        nt=nt+1
        idt(nt)=id1(i1)
        clt(nt)=cl1(i1)
        endif
      enddo !i1

! grab further unique targets in list 2 and add to total list

       do i2=1,n2
        if (cl2(i2).eq.m) then
          do it=1,nt
          if (id2(i2).eq.idt(it)) goto 100 !target already in list, skip
          enddo !it
          nt=nt+1
          idt(nt)=id2(i2)
          clt(nt)=cl2(i2)
100       continue
        endif
      enddo !i2

      enddo !m



! map out which targets are in which lists

      do it=1,nt

      l1(it)=.false.
      l2(it)=.false.

        do i1=1,n1
         if (idt(it).eq.id1(i1)) l1(it)=.true.
        enddo! n1

        do i2=1,n2
         if (idt(it).eq.id2(i2)) l2(it)=.true.
        enddo! n2


      enddo! it

! make double and single cover histograms

      do m=1,nclass

       hist2(m)=0
       hist1(m)=0

       do it=1,nt

       if (clt(it).eq.m) then
       ne=0
       if (l1(it)) ne=ne+1
       if (l2(it)) ne=ne+1

       if (ne.eq.0) stop 'No Match!'
       if (ne.eq.2) hist2(m)=hist2(m)+1
       if (ne.eq.1) hist1(m)=hist1(m)+1

       endif !m

       enddo !it

      enddo !m

! calculate checksum

      ncheck=0

      do m=1,nclass
      ncheck=ncheck+hist1(m)
      ncheck=ncheck+hist2(m)
      enddo!m

      ncheck=ncheck-nt

      if (ncheck.ne.0) stop 'Error in sorting checksum'

! calculate figure of merits

      fom2=0
      fom1=0

      do m=1,nclass
      fom2=fom2+weights(m)*float(hist2(m))
      fom1=fom1+weights(m)*float(hist1(m))
      enddo! m

      return

      end subroutine pair_sort

!-------------------------------------------------------------

      subroutine triple_sort(n_max,n1,id1,cl1,n2,id2,cl2,n3,id3,cl3,fom3,fom2,fom1)
      use order_and_weights, only : weights
      use target_cat, only : nclass
      implicit none
      integer n_max
      integer id1(n_max),id2(n_max),id3(n_max)
      integer cl1(n_max),cl2(n_max),cl3(n_max)
      integer n1,n2,n3
      integer i1,i2,i3
      real fom3,fom2,fom1
      integer nt,it,idt(3*n_max),clt(3*n_max)
      integer m,ne
      logical l1(3*n_max),l2(3*n_max),l3(3*n_max)
      integer hist1(nclass),hist2(nclass),hist3(nclass)
      integer ncheck

! construct combined target list

      nt=0

      do m=1,nclass

! grab class m targets in list 1 and add to total list

      do i1=1,n1
        if (cl1(i1).eq.m) then
        nt=nt+1
        idt(nt)=id1(i1)
        clt(nt)=cl1(i1)
        endif
      enddo !i1

! grab further unique targets in list 2 and add to total list

       do i2=1,n2
        if (cl2(i2).eq.m) then
          do it=1,nt
          if (id2(i2).eq.idt(it)) goto 100 !target already in list, skip
          enddo !it
          nt=nt+1
          idt(nt)=id2(i2)
          clt(nt)=cl2(i2)
100       continue
        endif
      enddo !i2

! grab further unique targets in list 3 and add to total list

       do i3=1,n3
        if (cl3(i3).eq.m) then
          do it=1,nt
          if (id3(i3).eq.idt(it)) goto 200 !target already in list, skip
          enddo !it
          nt=nt+1
          idt(nt)=id3(i3)
          clt(nt)=cl3(i3)
200       continue
        endif
      enddo !i3

      enddo !m


! map out which targets are in which lists

      do it=1,nt

      l1(it)=.false.
      l2(it)=.false.
      l3(it)=.false.

        do i1=1,n1
         if (idt(it).eq.id1(i1)) l1(it)=.true.
        enddo! n1

        do i2=1,n2
         if (idt(it).eq.id2(i2)) l2(it)=.true.
        enddo! n2

        do i3=1,n3
         if (idt(it).eq.id3(i3)) l3(it)=.true.
        enddo! n3

      enddo! it

! make triple, double and single cover histograms

      do m=1,nclass

       hist3(m)=0
       hist2(m)=0
       hist1(m)=0

       do it=1,nt

       if (clt(it).eq.m) then
       ne=0
       if (l1(it)) ne=ne+1
       if (l2(it)) ne=ne+1
       if (l3(it)) ne=ne+1

       if (ne.eq.0) stop 'No Match!'
       if (ne.eq.3) hist3(m)=hist3(m)+1
       if (ne.eq.2) hist2(m)=hist2(m)+1
       if (ne.eq.1) hist1(m)=hist1(m)+1

       endif !m

       enddo !it

      enddo !m

! calculate checksum

      ncheck=0

      do m=1,nclass
      ncheck=ncheck+hist1(m)
      ncheck=ncheck+hist2(m)
      ncheck=ncheck+hist3(m)
      enddo!m

      ncheck=ncheck-nt

      if (ncheck.ne.0) stop 'Error in sorting checksum'

! calculate figure of merits

      fom3=0
      fom2=0
      fom1=0

      do m=1,nclass
      fom3=fom3+weights(m)*float(hist3(m))
      fom2=fom2+weights(m)*float(hist2(m))
      fom1=fom1+weights(m)*float(hist1(m))
      enddo! m

      return

      end subroutine triple_sort

!-------------------------------------------------------------
