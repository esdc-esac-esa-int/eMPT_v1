          
     include 'empt_modules.f90'

!-------------------------------------------------------

      program k_clean     
      use path_to_ref_files
      use config_parameters
      use target_cat
      use derived_parameters
      use pointings_and_groups, only : n_p,no_p,ra_p,dec_p,pa_ap_p
      use k_list

      implicit none
      real t_start,t_stop
      
      character*2 rn

      integer c0(4,365,171),c1(4,365,171),c2(4,365,171)
      integer i_c,ik
      integer pcount,npointings,nspoil
      integer maxscore1,maxscore23,maxscore45
      integer ip,kk,ii,jj
      integer cmax
      integer score_1,score_2,score_3,score_4,score_5


      call cpu_time(t_start)


      write(*,*)
      write(*,*) '----------------------------------------------------------------------------------------'
      write(*,*) '--- k_clean -- eMPT Contaminated Target Elimination Module Version 1.7 4 December 2021 -'
      write(*,*) '----------------------------------------------------------------------------------------'
      write(*,*)

      call load_derived_parameters
      call file_entry_k_clean()

! create output directory 

      write(rn,'(i2)') ntile
      rn=trim(adjustl(rn))
      if (ntile.lt.10) rn='0'//rn(1:1)
      
      k_list_file='trial_'//rn(1:2)//'/k_make_output/k_list_raw.txt'

      call system('mkdir trial_'//rn(1:2)//'/k_clean_output')
      
      k_list_mod_file='trial_'//rn(1:2)//'/k_clean_output/k_list_mod.txt'
      
      call catalog_read()      

       write(*,*)
       write(*,*) ' Reading target input catalog'

       call pointings_read()

       write(*,*)
       write(*,*) ' Reading pointings to be processed'

       npointings=n_p

       call k_list_read()

      if (npointings.ne.n_p) then
      write (*,*)
      write (*,*) ' ******************************************************************************' 
      write (*,*) ' Error: Discrepancy in number of pointings in configuration and k_list_raw file' 
      write (*,*) ' ******************************************************************************' 
      write (*,*)
      stop
      endif


       write(*,*)
       write(*,'(a,a2)') ' Trial Identifier No.:        ',rn
       write(*,'(a,a)')  ' Configuration File:          ',confname
       write(*,'(a,a)')  ' Target Catalog:              ',catfile
!        write(*,'(a,a)')  ' Reference Image:             ',imname
!        write(*,'(a,a)')  ' Segmentation Map:            ',segname
       write(*,'(a,i0)') ' n_dither:                    ',n_dither
       write(*,'(a,i0)') ' Total number of targets in input catalog:         ',ntarg
       write(*,'(a,i0)') ' Number of Priority Classes in input catalog:      ',nclass
       write(*,'(a,i0)') ' Total number of pointings in configuration file:  ',n_p
       write(*,'(a,a)')  ' Input k_list file:                                ',k_list_file
       write(*,'(a,i0)') ' Maximum number of targets in k_list per pointing: ',n_max_k
       write(*,'(a,i0,2x,i0,2x,i0)') ' Maximum allowable contamination scores for Priority Class 1 targets:         ',max_c_score1_1,max_c_score23_1,max_c_score45_1
       write(*,'(a,i0,2x,i0,2x,i0)') ' Maximum allowable contamination scores for all other Priority Class targets: ',max_c_score1_2,max_c_score23_2,max_c_score45_2
       write(*,*)




! loop through all pointings


 
      pcount=0

      write(*,*)

      do ip=1,npointings

      pcount=pcount+1

      write(*,100) ' Processing pointing ',pcount,' out of ',npointings
100   format(a,i0,a,i0)

      
       call make_spoiler_masks3(ra_p(ip),dec_p(ip),pa_ap_p(ip),c0,c1,c2)


! run through targets in k_list for pointing  ip

      nspoil=0

 
      do ik=1,nk(ip)

      kk=kt(ip,ik)
      ii=it(ip,ik)
      jj=jt(ip,ik)

      score_1=0
      score_2=0
      score_3=0
      score_4=0
      score_5=0
      

      if (pri_k(ip,ik).lt.0) goto 777 !target excluded from consideration as target in k-list


! Shutter 1

     if (c0(kk,ii,jj).ge.1) score_1=1

     if ((c0(kk,ii,jj).gt.1).or.(c1(kk,ii,jj+1).gt.1).or.(c1(kk,ii,jj-1).gt.1)) then

     cmax=max(c0(kk,ii,jj),c2(kk,ii,jj-1),c1(kk,ii,jj+1))
     
     score_1=score_1+4*(cmax-1)
     
     endif
     
! Shutter 2

     if ((c0(kk,ii,jj+1).gt.0).or.(c2(kk,ii,jj).gt.0)) then

     cmax=max(c0(kk,ii,jj+1),c2(kk,ii,jj))
     
     score_2=4*cmax
     
     endif
     
! Shutter 3

     if ((c0(kk,ii,jj-1).gt.0).or.(c1(kk,ii,jj).gt.0)) then

     cmax=max(c0(kk,ii,jj-1),c1(kk,ii,jj))
     
     score_3=4*cmax
     
     endif

! Shutter 4

     if (c2(kk,ii,jj+1).gt.0) then

     cmax=c2(kk,ii,jj+1)
     
     score_4=4*cmax
     
     endif

! Shutter 5

     if (c1(kk,ii,jj-1).gt.0) then

     cmax=c1(kk,ii,jj-1)
     
     score_5=4*cmax
     
     endif


      maxscore1=max_c_score1_2
      maxscore23=max_c_score23_2
      maxscore45=max_c_score45_2


     if (pri_k(ip,ik).eq.1) then
      maxscore1=max_c_score1_1
      maxscore23=max_c_score23_1
      maxscore45=max_c_score45_1
     endif

     if ((score_1.gt.maxscore1).or.(score_2.gt.maxscore23).or.(score_3.gt.maxscore23).or.(score_4.gt.maxscore45).or.(score_5.gt.maxscore45)) then
      if (pri_k(ip,ik).gt.0)  then
       pri_k(ip,ik)=-pri_k(ip,ik) !make priority class negative for contaminated target ik of pointing ip
       nspoil=nspoil+1
      endif
     endif

777   continue

      enddo! nk

       write(*,220) nspoil,' out of ',nk(ip),' targets flagged as contaminated (',float(nspoil)/float(nk(ip))*100.,'%) and removed from mask'
220   format(x,i0,a,i0,a,f5.2,a,i0,a)

      enddo !ip

! write out modified k-list file

       open(9,file=k_list_mod_file,  status='replace',access='sequential', form='formatted')

      write(9,*) npointings

      do ip=1,npointings

      write(9,*) ip,no_p(ip)

      write(9,*) nk(ip)

      do ik=1,nk(ip)
      write(9,211) id(ip,ik),id_cat(ip,ik),pri_k(ip,ik),kt(ip,ik),it(ip,ik),jt(ip,ik),rx(ip,ik),ry(ip,ik)
      enddo !ik loop

      enddo !ip pointing loop

      close(9)

211   format(I10,2x,I10,x,I3,3x,I1,x,I3,x,I3,x,f6.3,x,f6.3)
      

! Append parameters to configuration file


      open(11,file=confname,status='old', access='sequential', form='formatted')

      call find_marker(11,'#KC#')

      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- SECTION AUTOMATICALLY APPENDED BY THE K_CLEAN MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a)') '# The modified k-list file with contaminated targets marked with negative Priority Class designations is located in:'
      write(11,'(a)') '#'
      write(11,'(a)') '#    ./trial_'//rn(1:2)//'/k_clean_output/k_list_mod.txt'
      write(11,'(a)') '#'     
      write(11,'(a)') '#'
      write(11,'(a)') '#'

      
      write(11,'(a)') '# The following suggested default parameters should be reviewed and modified as needed before running the m_make module:'
      write(11,'(a)') '#'

      write(11,'(a)') '# Number of Priority Classes present in Input Catalog:'
      write(11,*) nclass
      write(11,'(a)') '#'
      write(11,'(a)') '#'

      if (n_dither.eq.1) then

      write(11,'(a)') '# Order in which segmented k-lists are to be fed to the Arribas algorithm'


      write(11,'(a)') '# Default sequence in strict order of Priority Class:'
      write(11,*) (i_c,i_c=1,nclass)
      write(11,'(a)') '#'

      endif !n_dither=1


      if (n_dither.eq.2) then

      write(11,'(a)') '# Order in which segmented k-lists are to be fed to the Arribas algorithm'
      write(11,'(a)') '# and their spectra attempted placed on the FPA without overlap'
      write(11,'(a)') '#'
      write(11,'(a)') '# Example: Pr. Class 1-4 Commons first, followed by Pr. Class 1-4 '
      write(11,'(a)') '# Uniques, and finally Pr.Class 5'
      write(11,'(a)') '# fillers minimizing commonality'
      write(11,'(a)') '#'
      write(11,'(a)') '# Priority Class:               1      2     3      4      5 '
      write(11,'(a)') '#                           -------------------------------------'
      write(11,'(a)') '# Fully Common Targets:     |   1  |   2  |  3  |   4   |  9  | '
      write(11,'(a)') '# Unique Targets:           |   5  |   6  |  7  |   8   | 10  | '
      write(11,'(a)') '#                           -------------------------------------'
      write(11,'(a)') '#'



      write(11,'(a)') '# Sequence optimizing target commonality between dithers with highest priority class filler targets placed last'
      write(11,*) (i_c,i_c=1,nclass-1),2*nclass-1
      write(11,*) (i_c+nclass-1,i_c=1,nclass-1),2*nclass


      write(11,'(a)') '#'

      endif !n_dither=2


      if (n_dither.eq.3) then

      write(11,'(a)') '# Order in which segmented k-lists are to be fed to the Arribas algorithm'
      write(11,'(a)') '# and their spectra attempted placed on the FPA without overlap'
      write(11,'(a)') '#'
      write(11,'(a)') '# Example: Pr. Class 1-4 Commons first, followed by Pr. Class 1-4 '
      write(11,'(a)') '# Partials, followed by Pr. Class 1-4 Uniques, and finally Pr.Class 5'
      write(11,'(a)') '# fillers minimizing commonality'
      write(11,'(a)') '#'
      write(11,'(a)') '# Priority Class:               1      2      3      4      5 '
      write(11,'(a)') '#                           -------------------------------------'
      write(11,'(a)') '# Fully Common Targets:     |   1  |   2  |  3  |   4   |  15  | '
      write(11,'(a)') '# Partially Common Targets: |   5  |   6  |  7  |   8   |  14  | '
      write(11,'(a)') '# Unique Targets:           |   9  |  10  | 11  |  12   |  13  | '
      write(11,'(a)') '#                           -------------------------------------'
      write(11,'(a)') '#'



      write(11,'(a)') '# Sequence optimizing target commonality between dithers with highest priority class filler targets placed last'
      write(11,*) (i_c,i_c=1,nclass-1),3*nclass-2
      write(11,*) (i_c+nclass-1,i_c=1,nclass-1),3*nclass-1
      write(11,*) (i_c+2*(nclass-1),i_c=1,nclass-1),3*nclass

      endif !n_dither=3

     
      if (n_dither.gt.1) then

      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# Limits defining minimum and maximum separations between legal multiple pointing dithers'
      write(11,'(a)') '# min_dx, max_dx, min_dy max_dy [shutter X and Y facets]'
      write(11,'(a)') '#'
      write(11,'(a)') '# Rules: 1) Legal dithers must be separated by no more than max_dx in X'
      write(11,'(a)') '#        2) Legal dithers must be separated by no more than max_dy in Y'
      write(11,'(a)') '#        3) If two dithers are closer together than min_dy in Y, then'
      write(11,'(a)') '#           they must be separated by more than min_dx in X'
      write(11,'(a)') '#'
      write(11,'(a)') '# To disable set to: 0. 1000. 0. 1000.'
      write(11,'(a)') '# [Default: 3.0 15.0 1.0 8.0]'
      if (n_p.gt.n_dither) write(11,*) 3.0,15.0,1.0,8.0
      if (n_p.le.n_dither) then
      write(11,*) 0.,1000.,0.,1000.0
      write(11,'(a)') '# Disabled due to no need for dither censoring'
      endif !n_p<=n_dither

      endif !n_dither>1

      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- END OF SECTION APPENDED BY THE K_CLEAN MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a)') '#  -- Modifying any parameter in the above section requires re-running the m_make and following modules ---'
      write(11,'(a)') '#'
      write(11,'(a)') '#MM# -- Marker do not delete or move this line'

      close(11)

      write(*,*)
      call cpu_time(t_stop)
      write(*,'(a)')      '------------------------------------------------'
      write(*,'(a,f7.1)') ' CPU time used by k_clean module [s]:', t_stop-t_start
      write(*,'(a)')      '------------------------------------------------'
      write(*,*)


      end program k_clean
      
!------------------------------------------------------------
 
 
include './reference_files/front_transforms_updt.f90'

include './reference_files/back_transforms_updt.f90'

include './reference_files/shutter_routines_new.f90'

include './reference_files/misc_routines_new.f90'

!------------------------------------------------------------

      subroutine file_entry_k_clean()
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
! 

! work out pointing reference point coordinates in V2,V3 on-sky systems

      call msa_to_sky_bd(xm_ref,ym_ref,ax_refv2,ay_refv3) !v2,V3


      call skip_hashes(9)
      read(9,*) !szone_x,szone_y

      call skip_hashes(9)
      read(9,*) angle_to_target
      call dva_calc(angle_to_target,dva_mag)

      call skip_hashes(9)
      read(9,*) max_c_score1_1,max_c_score23_1,max_c_score45_1
 
      call skip_hashes(9)
      read(9,*) max_c_score1_2,max_c_score23_2,max_c_score45_2

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

      end subroutine file_entry_k_clean

!------------------------------------------------------------
