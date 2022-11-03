! v2.0 rewriten to include MPT method of roll designation
! v2.1 rewritten to handle filter-dependent FORE transformation
! v2.2 increased reachable target id to I10 in screen output
! v2.3 increased digits of dva_mag output, fixed imname and segname output
! v2.4 fixed allocation/deallocation problem with pa arrays when looping 
      
     include 'empt_modules.f90'

!-------------------------------------------------------
      program ipa      
      use path_to_ref_files
      use config_parameters
      use target_cat
      use derived_parameters
      use pointings_and_groups

      implicit none
      
! --- optimal pointing variables
      integer, parameter :: msx=7 ! vector space sampling as odd fraction of acceptance zone in x
      integer, parameter :: msy=13 ! vector space sampling as odd fraction of acceptance zone in x
      real :: dx,dy
      integer :: nx,ny
      integer, allocatable :: smap(:,:) !digital shift vector array
      real, allocatable :: tx(:),ty(:) ! tangential coordinates of targets
      real, allocatable :: sx(:),sy(:) ! v2,v3 coordinates of viable slitlets (mas)
      real*8, allocatable :: ra_f(:),dec_f(:)
      real, allocatable :: pa_v3_f(:)
      integer, allocatable :: n_targs_f(:),ids_f(:,:),id_list(:)
      real ax_refra_p,ay_refdec_p
      real apa_v3,roll_corr

      
! --- input and output files
      character*70 list_peakfile
      character*70 group_pointing_file
      character*70 sky_plot
      character*70 swarm_plot
      character*70 pointing_map_plot
      character*70 digital_map_plot
      character*2 rn
! - spectral parameters
!      integer :: n_disp
      real :: lambda1,lambda2
      character*6 gname
      common/disp/lambda1,lambda2,gname
! --- MSA Parameters-------
      integer :: mx,my
      real :: x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
! --- slitmap --------------
      integer :: slitmap(4,365,171)
      common/slit_map/slitmap
! - prism seperations
      integer prism_sep_p(4,365,171),prism_sep_m(4,365,171)
      common/prism_length/prism_sep_p,prism_sep_m

      real :: t_start,t_stop
      integer i,j,k,ii,jj,ih,jh,iho,jho,ig,jg
      integer nt,ns
      integer n_vectors
      real ax,ay,xm,ym,aax,aay
      real dist_x,dist_y
      integer max_peak,min_peak,n_peaks,peak_val
      integer ip
      real dra,ddec
      real dlim_x,dlim_y,dax,day,ss,cc,pa_v3_tmp
      real*8 ra_tmp,dec_tmp
      integer it_reach,n_pass
      logical group_output
      integer initial_dig_depth,min_keep_depth
      real*8 ra_old,dec_old,ra_diff,dec_diff
      integer loopmax
     
      call cpu_time(t_start)


      write(*,*)
      write(*,*) '-------------------------------------------------------------------------------'
      write(*,*) '--- ipa - Initial Pointing Algorithm Module - Version 2.4 -- 1 September 2022 -'
      write(*,*) '-------------------------------------------------------------------------------'
      write(*,*)

      call load_derived_parameters
      dx=av_xopen/float(msx) ! vector space sampling in mas in x
      dy=av_yopen/float(msy) ! vector space sampling in mas in y  


      group_output=.false.

      initial_dig_depth=3
      min_keep_depth=2
      min_numb_groups=1
      min_pointings_in_group=1

! read trial configuration file

      call file_entry_ipa()
      
      n_filt=2


! read ipa configuration file parameters

     open(9,file=trim(ref_path)//'ipa.conf',status='old',access='sequential',form='formatted',err=10)
     call skip_hashes(9)
     read(9,*) group_output
     call skip_hashes(9)
     read(9,*) initial_dig_depth
     call skip_hashes(9)
     read(9,*) min_numb_groups
     read(9,*) min_pointings_in_group
     if (min_pointings_in_group.le.0) min_pointings_in_group=n_dither
     call skip_hashes(9)
     read(9,*) min_keep_depth
    close(9)
   
    goto 20
10  write(*,*) 'Error: IPA Configuration file not found'
    stop
20 continue     

! work out actual roll angle of central pointing

      call get_cat_ref_position()

      call get_mpt_roll_offset(cra,cdec,roll_corr)
      
      apa_v3=cpa_v3+roll_corr

! work out pointing reference point for present roll angle

      call roll_v3(apa_v3,ax_refv2,ay_refv3,ax_refra_p,ay_refdec_p) !RA,Dec

       write(*,*)
       write(*,'(a,I0)') ' Trial Identifier No.:        ',ntile
       write(*,'(a,a)') ' Configuration File:          ',confname
       write(*,'(a,a)') ' Disperser:                   ',disperser
       write(*,'(a,f4.2)') ' S_0:                         ',sthresh
       write(*,'(a,2f7.3)') ' Acceptance Zone:           ',raccx,raccy
       write(*,'(a,a)')' Target Catalog:              ',catfile
       write(*,'(a,2f13.7)') ' Catalog Median Reference Position:  ',ra_c,dec_c
       write(*,'(a,a)') ' Reference Image:            ',imname
       write(*,'(a,a)') ' Segmentation Map:           ',segname
       write(*,'(a,2f13.7,f10.4)') ' Nominal Pointing & Roll:  ',cra,cdec,cpa_ap
       write(*,'(a,f10.4)') ' Central Roll Correction:  ',roll_corr


       write(*,'(a,2f10.2)') ' Search Zone:            ',szone_x,szone_y
       write(*,'(a,f10.2,f9.5)') ' Angle_to_target,dva_mag:',angle_to_target,dva_mag
       write(*,'(a,i0)') ' n_dither:                    ',n_dither
       write(*,*)
       if (group_output) then
        write(*,'(a)') ' IPA Mode: Pointings Grouped'
        write(*,'(a,i0)') ' Initial dig depth::       ',initial_dig_depth
        write(*,'(a,i0)') ' Min pointings in group:   ',min_pointings_in_group
        write(*,'(a,i0)') ' Min number of groups:     ',min_numb_groups
       else
        write(*,'(a)') ' IPA Mode: Pointings Not Grouped'
        write(*,'(a,i0)') ' Initial dig depth::       ',initial_dig_depth
        write(*,'(a,i0)') ' Min below max to save:    ',min_keep_depth
       endif
        write(*,*)

 
! create output directory 

      write(rn,'(I2)') ntile
      rn=trim(adjustl(rn))
      if (ntile.lt.10) rn='0'//rn(1:1)
      
      call system('mkdir trial_'//rn(1:2))
      call system('mkdir trial_'//rn(1:2)//'/ipa_output')

! name output files

      sky_plot='trial_'//rn(1:2)//'/ipa_output/input_catalog_on_sky.ps/cps'
      swarm_plot='trial_'//rn(1:2)//'/ipa_output/on_sky_map_start.ps/cps'
      digital_map_plot='trial_'//rn(1:2)//'/ipa_output/digital_peak_map.ps/cps'
      group_pointing_file='trial_'//rn(1:2)//'/ipa_output/grouped_optimal_pointings.txt'
      pointing_map_plot='trial_'//rn(1:2)//'/ipa_output/map_optimal_pointings.ps/cps'

! read in priority 1 targets from catalog
     
       call cat_pri_one_read()


      if (ntarg.le.0) then
      write(*,*)
      write(*,*) 'Error: No Targets in Input Catalog'
      write(*,*)
      stop
      endif


      if (ntarg_1.le.0) then
      write(*,*)
      write(*,*) 'Error: No Priority Class 1 Targets in Input Catalog'
      write(*,*)
      stop
      endif

      if (nclass.le.0) then
      write(*,*)
      write(*,*) 'Error: No Priority Class Designations in Input Catalog'
      write(*,*)
      stop
      endif

       print *
       print *,'Targets in Input Catalog:  ',ntarg
       print *,'Priority 1 Targets:        ',ntarg_1
       print *,'Number of Priority Classes:',nclass
       


! initialize mode-specific instrument model subroutines

      call load_msa_dimensions_bd
      call load_fpa_dimensions_bd
      call msa_to_fpa_disp(n_disp,2.,xm,ym,ax,ay)

! allocate priority 1 target coordinate list and reachable table

      nt=ntarg_1
      allocate (tx(nt),ty(nt),targ_poss(nt))


! construct priority 1 coordinate swarm on sky_plot at nominal pointing

      do i=1,nt

      call tangcoord_single(cra,cdec,ra_t(i),dec_t(i),ax,ay) !Ra,Dec

      aax=ax*dva_mag+ax_refra_p !apply dva magnification and shift to origin at ref position
      aay=ay*dva_mag+ay_refdec_p


! Rotate to V2,V3

      call roll_v3(-apa_v3,aax,aay,ax,ay) !V2,V3

! Rotate to MSA coordinate system orientation

      tx(i)= ax*cos(phi_r)+ay*sin(phi_r) !xm,ym as
      ty(i)=-ax*sin(phi_r)+ay*cos(phi_r)

      tx(i)=tx(i)*1000. !mas
      ty(i)=ty(i)*1000.

      enddo !nt


! construct viable slitlet swarm on sky

      call read_shutter_values
      if (n_disp.eq.7) call read_prism_separations
      call read_base_msamap
      call make_slitlet_map
      call initialize_slitmap(n_disp,sthresh,1,1)

! determine how many viable slitlets are on msa

      ns=0

      do k=1,4
      do i=1,365
      do j=1,171
       if (slitmap(k,i,j).eq.4) ns=ns+1   
      enddo !j
      enddo !i
      enddo !k

      print *
      print *,'Number Viable Slitlets on MSA: ',ns

! allocate viable slitlet coordinate list
      
      allocate (sx(ns),sy(ns))
 
 ! load viable slitlet coordinate list
  
      ii=0
  
      do k=1,4
      do i=1,365
      do j=1,171

      if (slitmap(k,i,j).ne.4) cycle
      
      ii=ii+1
      
      call shutter_to_msa_bd(k,i,j,xm,ym)

      call msa_to_sky_bd(xm,ym,ax,ay) !V2,V3

! rotate to orientation of MSA system

      aax= ax*cos(phi_r)+ay*sin(phi_r)
      aay=-ax*sin(phi_r)+ay*cos(phi_r)

      sx(ii)=aax*1000. !mas
      sy(ii)=aay*1000.
        
      enddo !j
      enddo !i
      enddo !k

 ! Plot initial slitlet and target swarms on sky

      call plot_intial_swarms(nt,tx,ty,id_t,ns,sx,sy,szone_x,szone_y,swarm_plot)


! ------- construct digital shift vector array --------------------------------------


! calculate required size of digital shift array, allocate and null it

       nx=int(szone_x*1000./dx) ! shift vector array -nx:nx in x
       ny=int(szone_y*1000./dy) ! shift vector array -ny:ny in y

       allocate (smap(-nx:nx,-ny:ny))
       smap=0


! Fill in shift vector plane with all permissable target-slitlet combinations

        n_vectors=0

        do i=1,nt !run through all targets

        targ_poss(i)=.false. !non-shiftable target assumed initially

        do j=1,ns !run through all viable slitlets

! shift vector needed to move target i to center of viable slitlet j

        dist_x=sx(j)-tx(i)
        dist_y=sy(j)-ty(i)
        
! Is shift vector within allowable shift vector range?

        if ((abs(dist_x).lt.szone_x*1000.).and.(abs(dist_y).lt.szone_y*1000.)) then

        targ_poss(i)=.true. !target is shiftable

! Count peaks in shift vector plane

        n_vectors=n_vectors+1 !number of peaks in shift vector plane

! round shift vector to nearest digitized (signed) pixel

        ih=nint(dist_x/dx)
        jh=nint(dist_y/dy)

! Increment shift vector plane by one unit over acceptance zone footprint in x and y

        do ii=-(msx-1)/2,(msx-1)/2
        do jj=-(msy-1)/2,(msy-1)/2
        iho=ih+ii
        jho=jh+jj
         if ((abs(iho).le.nx).and.(abs(jho).le.ny)) then
          smap(iho,jho)=smap(iho,jho)+1
         endif
        enddo
        enddo

        endif !in allowable shift zone end

        enddo    !j loop end
        enddo    !i loop end


! -------  Digital shift vector array calculation complete ---------------------------     


! determine number of (not necessarily simultaneously) reachable targets

        n_reach=0
        do i=1,nt
        if (targ_poss(i)) n_reach=n_reach+1
        enddo

      if (n_reach.le.0) then
      write(*,*)
      write(*,*) 'Error: No Priority Class 1 Targets Reachable in Search Box at Nominal Pointing'
      write(*,'(a,a)') ' Consult figure: ',swarm_plot
      write(*,*)
      stop
      endif


      write(*,*)
      write(*,'(a,i0)') ' Total Number of reachable Priority Class 1 targets within Shift Range: ',n_reach
      
      write(*,*)
      write(*,*) 'Reachable targets:'
      write(*,*) '  ID_cat'
        
        do i=1,nt
         if (targ_poss(i)) then 
          write(*,'(i10)') id_t(i)
         endif
        enddo

       write(*,*)
       write(*,'(a,i0)') ' Total number of shift vectors in shift array: ',n_vectors
       write(*,'(a,2f7.1)') ' Shift Array digitization [mas]:',dx,dy


! Single Reachable Priority 1 target shortcut start ----------------

    if (n_reach.eq.1) then
    
       ra_diff=0.0001/3600./cos(cdec/57.2957795) ! mas
       dec_diff=0.0001/3600. !1 mas
       loopmax=10

       min_numb_groups=1
 
        n_pass=1
        max_peak=1
        call plot_digital_shift_map(nx,ny,smap,dx,dy,max_peak,digital_map_plot)
        deallocate (smap)

! figure out which target is the reachable one

        do i=1,nt
        if (targ_poss(i)) it_reach=i
        enddo

! figure out how many shutters within the search zone the lone target can be placed in

        n_peaks=0
        
        do j=1,ns

          dist_x=sx(j)-tx(it_reach)
          dist_y=sy(j)-ty(it_reach)

           if ((abs(dist_x).lt.szone_x*1000.).and.(abs(dist_y).lt.szone_y*1000.)) n_peaks=n_peaks+1
           
        enddo 
        
! allocate pointing save lists

       allocate (ra_p(n_peaks),dec_p(n_peaks),pa_v3_p(n_peaks),pa_ap_p(n_peaks),n_targs_p(n_peaks),ids_p(n_peaks,1))
 
       pa_v3_p(:)=apa_v3
       pa_ap_p(:)=0.
       
       ids_p=0 

! go through through legal shifts a second time to calculate exact pointing values            

        n_peaks=0
        
        do j=1,ns

          dist_x=sx(j)-tx(it_reach)
          dist_y=sy(j)-ty(it_reach)

           if ((abs(dist_x).lt.szone_x*1000.).and.(abs(dist_y).lt.szone_y*1000.)) then
        
            n_peaks=n_peaks+1
 
            ax=dist_x/1000. !as
            ay=dist_y/1000.
            
            aax= ax*cos(phi_r)-ay*sin(phi_r) !rotate offset to v2,v3 orientation
            aay= ax*sin(phi_r)+ay*cos(phi_r)

            call roll_v3(pa_v3_p(n_peaks),aax,aay,ax,ay) !ra,dec arcsec orientation

           ax=-ax
           ay=-ay

! calculate matching pointing coordinate 

            call coordtang_single(cra,cdec,ax,ay,ra_old,dec_old)

 ! attach individual roll to pointing and do mac loopmax iterations on pointing
 
        do ii=1,loopmax

        call get_mpt_roll_offset(ra_old,dec_old,roll_corr)      
         pa_v3_p(n_peaks)=cpa_v3+roll_corr
 
        call roll_v3(pa_v3_p(n_peaks),aax,aay,ax,ay) !ra,dec arcsec orientation
        call coordtang_single(cra,cdec,-ax,-ay,ra_p(n_peaks),dec_p(n_peaks))
        
        if ( (abs(ra_old-ra_p(n_peaks)).lt.ra_diff).and.(abs(dec_old-dec_p(n_peaks)).lt.dec_diff)) exit
        
        ra_old=ra_p(n_peaks)
        dec_old=dec_p(n_peaks)
        
        enddo !ii
               
            n_targs_p(n_peaks)=1
            ids_p(n_peaks,1)=id_t(it_reach)

! fill in pa_ap_p array

           pa_ap_p(n_peaks)=pa_v3_p(n_peaks)+(180.-phi)
           if (pa_ap_p(n_peaks).gt.360.) pa_ap_p(n_peaks)=pa_ap_p(n_peaks)-360.
           if (pa_ap_p(n_peaks).lt.0.) pa_ap_p(n_peaks)=pa_ap_p(n_peaks)+360.
 
           endif
        
          enddo !ns

           deallocate (tx,ty,sx,sy)

           max_targs=1
           min_targs=0        
  
           write(*,*)
           write(*,'(a,i0,a,i0,a)') ' Pointings providing simultaneous coverage from ',max_targs,' and down to ',max(min_targs,1),' targets will be sorted and saved'
 
       
       n_p=n_peaks
       
       np_max_targs=n_peaks
       np_min_targs=0


       write(*,*)
       write(*,'(a,i0,a,i0)') ' Final number of identified ',max_targs,' target coverage pointings: ',np_max_targs
         
        goto 10001
        
        endif

! Single Reachable Priority 1 target shortcut end -------------------------

        deallocate (tx,ty,sx,sy)
 
! determine maximum target overlap hits in digital shift vector array ------

        max_peak=-666
        do i=-nx,nx
        do j=-ny,ny
        max_peak=max(max_peak,smap(i,j))
        enddo
        enddo

        call plot_digital_shift_map(nx,ny,smap,dx,dy,max_peak,digital_map_plot)

        min_peak=max(max_peak-initial_dig_depth,1) ! initial starting assumption that single pass through is sufficient to locate maximum group with more than n_dither pointings
 
        ng_max_ok=min_numb_groups+1 ! initial starting assumption that single pass through is sufficient to locate maximum group with more than n_dither pointings

        n_pass=0
        

!-----------------------------------------------------------------------------------

1111   n_pass=n_pass+1 ! Entry point for additional deeper passes through digital shift array

!-----------------------------------------------------------------------------------

       if (ng_max_ok.le.min_numb_groups) then
       write(*,*)
       write(*,'(a,i0,a)') ' Pass ',n_pass,' through convolved shift array initiated'
       endif

       print *
       write(*,'(a,I0)') ' Maximum peak in convolved shift array:  ',max_peak

       if (ng_max_ok.le.min_numb_groups) min_peak=max(min_peak-1,1) !dig one unit deeper on each additional pass

       write(*,'(a,I0)') ' Minimum amplitude peaks to be explored: ',min_peak

 
 ! determine number of digital peaks worth fine-tuning
 
        n_peaks=0
        do i=-nx,nx
        do j=-ny,ny
        if (smap(i,j).lt.min_peak) cycle
        n_peaks=n_peaks+1
        enddo
        enddo
       
       write(*,'(a,I0)') ' Total number of digital peaks to be fine-tuned: ',n_peaks
       
! allocate pointing save lists

       allocate (ra_p(n_peaks),dec_p(n_peaks),pa_v3_p(n_peaks),n_targs_p(n_peaks),ids_p(n_peaks,2*max_peak),id_list(2*max_peak))
       ids_p(:,:)=0 
  
 
 ! fill out starting point digital peak-derived pointing save lists

        ii=0
        
        do jj=max_peak,min_peak,-1
        
        do i=-nx,nx
        do j=-ny,ny
        
        if (smap(i,j).ne.jj) cycle
        
        ii=ii+1

 ! get digitized shift vector of pixel i,j and calculate approximate ra,dec of peak

        dist_x=float(i)*dx/1000. !xm,ym
        dist_y=float(j)*dy/1000.

        ax= dist_x*cos(-phi_r)+dist_y*sin(-phi_r) !V2,V3 as orientation
        ay=-dist_x*sin(-phi_r)+dist_y*cos(-phi_r)

        call roll_v3(apa_v3,ax,ay,aax,aay) !ra,dec arcsec orientation

        aax=-aax
        aay=-aay

        call coordtang_single(cra,cdec,aax,aay,ra_p(ii),dec_p(ii))
       
        n_targs_p(ii)=smap(i,j)
        
 ! attach individual roll to pointing
 
        call get_mpt_roll_offset(ra_p(ii),dec_p(ii),roll_corr)      
        pa_v3_p(ii)=cpa_v3+roll_corr
        
               
        enddo !j
        enddo !i
        
        enddo !jj

       
 
! -------  Convert in place to fine-tune extracted list of Digital Peaks  ---------------  
 
           print *
 
           max_targs=-666

           ii=0
 
           do ip=1,n_peaks
           
            call fine_tune_pointing_roll(ra_p(ip),dec_p(ip),pa_v3_p(ip),n_targs_p(ip),id_list,max_peak)
     
           do j=1,n_targs_p(ip)
           ids_p(ip,j)=id_list(j)
           enddo

           ii=ii+1
           
           max_targs=max(max_targs,n_targs_p(ip))
           
           enddo !ip n_peaks 


           write(*,'(a,i0,a,i0,a)') ' ',ii,' peaks out of ',n_peaks,' processed'
           write(*,*)
           write(*,'(a,i0)') ' Maximum simultaneous Priority 1 Target Coverage: ',max_targs
 
 
 ! -------  Fine-tuning of list of Digital Peaks completed  ---- 
  
           if (group_output) then

!            if (ng_max_ok.gt.min_numb_groups) min_targs=max_targs !initial starting point     
!            if (ng_max_ok.le.min_numb_groups) min_targs=max(min_targs-1,1) !go one step deeper on every subsequent pass 
             if (ng_max_ok.gt.min_numb_groups) then
              min_targs=max_targs !initial starting point
             else    
              min_targs=max(min_targs-1,1) !go one step deeper on every subsequent pass 
             endif
          else
            
            min_targs=max(max_targs-min_keep_depth,1)
           
           endif !group_output
           
           write(*,*)
           write(*,'(a,i0,a,i0,a)') ' Pointings providing simultaneous coverage from ',max_targs,' and down to ',min_targs,' targets will be saved'
           
           print *
           print *, 'Eliminating low-coverage and duplicate pointings from fine-tuned list (patience please)'


       dlim_x=1.0*av_xopen/1000. !as
       dlim_y=1.0*av_yopen/1000.

!        ss=sin(phi_r-cpa_v3/57.2957795)
!        cc=cos(phi_r-cpa_v3/57.2957795)

       
       do i=1,n_peaks 

! only care about peaks above min_targs, set n_targs_p(i)=0 if below min_targs
       
       if (n_targs_p(i).lt.min_targs) then
        n_targs_p(i)=0
         cycle
       endif
       
        do j=1,i-1 ! all previously examined pointings before i

        if (n_targs_p(j).lt.min_targs) cycle
 
        call tangcoord_single(ra_p(j),dec_p(j),ra_p(i),dec_p(i),dra,ddec)  !ra,dec

! Rotate to MSA coordinate system orientation

       ss=sin(phi_r-pa_v3_p(i)/57.2957795)
       cc=cos(phi_r-pa_v3_p(i)/57.2957795)

          dax= dra*cc+ddec*ss
          day=-dra*ss+ddec*cc

         if ((abs(dax).gt.dlim_x).or.(abs(day).gt.dlim_y)) cycle !pointing j not spatially coincident with pointing i -> next j

! Pointings i and j coincide spatially, does j cover a greater or equal number of targets than trial pointing i? -> if yes, pointing i uninteresting, next j

         if (n_targs_p(j).ge.n_targs_p(i)) then
          n_targs_p(i)=0  
         cycle 
         endif   
         
! Pointings i and j coincide spatially, and trial pointing i has higher coverage than 
! previous pointing j at same location -> swap pointing i with j. New trial pointing i now uninteresting


          ra_tmp=ra_p(j)
          dec_tmp=dec_p(j)
          pa_v3_tmp=pa_v3_p(j)

         ra_p(j)=ra_p(i)
         dec_p(j)=dec_p(i)
         pa_v3_p(j)=pa_v3_p(i)
         n_targs_p(j)=n_targs_p(i)
         do ii=1,n_targs_p(i)
           ids_p(j,ii)=ids_p(i,ii)
         enddo     

          ra_p(i)=ra_tmp
          dec_p(i)=dec_tmp
          pa_v3_p(i)=pa_v3_tmp
          n_targs_p(i)=0            
        
       enddo !j
       
       enddo !i



! count number of unique non-zero fine-tined peaks of interest 

       np_min_targs=0 ! number of pointings having at least min_targs coverage

       do ip=1,n_peaks
       if (n_targs_p(ip).ge.min_targs) np_min_targs=np_min_targs+1
        enddo

       n_p=np_min_targs ! new final number of saved fine-tuned pointings

       write(*,*)
       write(*,'(a,i0,a,i0)') ' Final number of identified ',min_targs,' or greater target coverage pointings: ',np_min_targs
 

! purge final pointing list of uninteresting nulled pointings and compact list

! allocate temporary storage

       allocate (ra_f(n_p),dec_f(n_p),pa_v3_f(n_p),n_targs_f(n_p),ids_f(n_p,max_targs))

       ii=0

       do peak_val=max_targs,min_targs,-1

       do ip=1,n_peaks
    
       if (n_targs_p(ip).eq.peak_val) then
        ii=ii+1
        ra_f(ii)=ra_p(ip)
        dec_f(ii)=dec_p(ip)
        pa_v3_f(ii)=pa_v3_p(ip)
        n_targs_f(ii)=n_targs_p(ip)
        do j=1,n_targs_p(ip)
         ids_f(ii,j)=ids_p(ip,j)
        enddo
       endif    

       enddo !ip
 
       enddo !peak_val
                     
       deallocate (ra_p,dec_p,pa_v3_p,n_targs_p,ids_p)
             
       allocate (ra_p(n_p),dec_p(n_p),pa_v3_p(n_p),pa_ap_p(n_p),n_targs_p(n_p),ids_p(n_p,max_targs))
     
       ra_p=ra_f
       dec_p=dec_f
       pa_v3_p=pa_v3_f
       n_targs_p=n_targs_f
       ids_p=ids_f
       
        deallocate (ra_f,dec_f,pa_v3_f,n_targs_f,ids_f)

! fill in pa_ap_p array

     do ip=1,n_p
      pa_ap_p(ip)=pa_v3_p(ip)+(180.-phi)
      if (pa_ap_p(ip).gt.360.) pa_ap_p(ip)=pa_ap_p(ip)-360.
      if (pa_ap_p(ip).lt.0.) pa_ap_p(ip)=pa_ap_p(ip)+360.
     enddo

        
!         write(*,*)
!         do ip=1,n_p
!         write(*,*) ra_p(ip),dec_p(ip),n_targs_p(ip),(ids_p(ip,j),j=1,n_targs_p(ip))
!         enddo
      

10001    continue ! re-entry jump point for special n_reach=1 case


! group final pointings by targets covered

       if (group_output) then

       write(*,*)
       write(*,*) 'Grouping optimal pointings per targets covered'
       write(*,*)

       call group_pointings()


      do i=1,n_groups
      write(*,'(a,i0,a,i0,a,i0,a,i0,a)') ' Pass ',n_pass,' - Group ',i,' covers ',gr_tcount(gr_indx(i)),' targets simultaneously in ',gr_pcount(gr_indx(i)),' pointings'

      do j=1,n_groups
      if (gr_subs(gr_indx(i),gr_indx(j)).ne.0) write(*,'(a,i0,a,i0,a,i0,a)') '          Group ',i,' is a subset of group ',j,', adding a further ',gr_pcount(gr_indx(j)),' pointings'
      enddo
      
      if (gr_ptcount(gr_indx(i)).gt.gr_pcount(gr_indx(i))) write(*,'(a,i0,a,i0,a)') '          The targets of group ',i,' are covered in a total of ',gr_ptcount(gr_indx(i)),' pointings'
      
      write(*,'(a)') '         ID_cat'        
        do j=1,n_targs_p(gr_targs(gr_indx(i)))
        write(*,*) ids_p(gr_targs(gr_indx(i)),j)
        enddo !j

      write(*,*)


      enddo !i
 
! determine from scratch whether there is a sufficient number of groups containing a sufficient number of pointings  
! exceeding the specified value of n_dither

      max_ok_gr_targs=-9999
      min_ok_gr_targs=9999
      do i=1,n_groups
       if (gr_ptcount(gr_indx(i)).ge.min_pointings_in_group) then
        max_ok_gr_targs=max(max_ok_gr_targs,gr_tcount(gr_indx(i)))
        min_ok_gr_targs=min(min_ok_gr_targs,gr_tcount(gr_indx(i)))
       endif
      enddo

! check whether max_ok_gr_targs and min_ok_gr_targs differ by more than one and if so whether min_ok_gr_targs can be increased

    if (max_ok_gr_targs-min_ok_gr_targs.ge.1) then
    
       do j=max_ok_gr_targs-1,min_ok_gr_targs,-1
       ng_min_ok=0
        do i=1,n_groups
         if (gr_tcount(gr_indx(i)).eq.j) ng_min_ok=ng_min_ok+1
        enddo !i
       if (ng_min_ok.gt.1) then
        min_ok_gr_targs=j
        exit
       endif
       enddo !j
    
    endif      

     ng_min_ok=0
     ng_max_ok=0
     ng_all_ok=0


      do i=1,n_groups
       if (gr_ptcount(gr_indx(i)).ge.min_pointings_in_group) then
        if (gr_tcount(gr_indx(i)).eq.max_ok_gr_targs) ng_max_ok=ng_max_ok+1
        if (gr_tcount(gr_indx(i)).eq.min_ok_gr_targs) ng_min_ok=ng_min_ok+1
        if (gr_tcount(gr_indx(i)).ge.min_ok_gr_targs) ng_all_ok=ng_all_ok+1
       endif
      enddo

      if (min_ok_gr_targs.ge.max_ok_gr_targs) ng_min_ok=0 !To handle only max_ok_gr_targs situation


     if (ng_all_ok.ge.min_numb_groups) then 
      write(*,*)  
      write(*,'(a,i0,a,i0,a,i0,a)') ' ',ng_max_ok,' groups covering ',max_ok_gr_targs,' targets at ',min_pointings_in_group,' or more pointings were found' 
      if (ng_min_ok.gt.0) write(*,'(a,i0,a,i0,a,i0,a)') ' ',ng_all_ok,' groups covering ',min_ok_gr_targs,' or more targets at ',min_pointings_in_group,' or more pointings were found' 
     endif

     
     if (ng_all_ok.lt.min_numb_groups) then !loop back for deeper probe of digital peak map
     
     deallocate (ra_p,dec_p,pa_v3_p,pa_ap_p,n_targs_p,ids_p,id_list)
     deallocate (gr_flag,gr_targs,gr_pcount,gr_tcount,gr_indx,gr_inv,gr_subs,gr_ptcount)
    
         
      write(*,*)  
      write(*,'(a,i0,a)') ' No group with at least ',min_pointings_in_group,' pointings was found'        

      if (min_targs.le.1) then
        write(*,*)
        write(*,'(a)') ' Error: All digital peaks in the shift vector plane have been fully explored '
        write(*,'(a)') ' Algorithm cannot proceed further'
        write(*,*)
        stop
      endif

 

     goto 1111
      endif


     
! we have identified at least min_numb_groups containing at least min_targets_in_group pointings covering the same group of targets and are good to go


! print out final pointing groups to file

      call write_groups_to_file(group_pointing_file)


      endif ! plot_output=.true.


! fill in pa_ap_p array

     do ip=1,n_p
      pa_ap_p(ip)=pa_v3_p(ip)+(180.-phi)
      if (pa_ap_p(ip).gt.360.) pa_ap_p(ip)=pa_ap_p(ip)-360.
      if (pa_ap_p(ip).lt.0.) pa_ap_p(ip)=pa_ap_p(ip)+360.
     enddo


! plot final optimal pointings in shift plane


      call plot_pointing_map(pointing_map_plot,nx,ny,dx,dy)

! plot full catalog on sky

      call draw_cat_on_sky2(sky_plot)

! open  configuration file 

      open(11,file=confname,status='old',access='sequential', form='formatted')

! Append peak IPA output summary and pointing coordinate list to configuration file

      call find_marker(11,'#IP#')

      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- SECTION AUTOMATICALLY APPENDED BY THE IPA MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a,i0)') '# Total number of observable Priority Class 1 targets in range: ',n_reach
    
      write(11,'(a)') '#'
      write(11,'(a,i0)') '# Maximum number of Priority Class 1 targets observable simultaneously in a single pointing: ',max_targs
      write(11,'(a)') '#'
      if (n_dither.gt.1) then
       write(11,'(a,i0,a,i0)') '# Maximum number of Priority Class 1 targets observable simultaneously in ',n_dither,' or more dithered pointings: ',max_ok_gr_targs
       write(11,'(a)') '#'
      endif

      write(11,'(a)') '#'

      write(11,'(a)') '# A detailed sorted list and map of all optimal pointings identified by the IPA module can be found in the directory:'
      write(11,'(a)') '#'
      write(11,'(a)') '#   ./trial_'//rn(1:2)//'/ipa_output/'
      write(11,'(a)') '#'
     
     if (group_output) then
     
     if (ng_all_ok.ge.min_numb_groups) then 
      write(11,'(a)') '#'
      write(11,'(a,i0,a,i0,a,i0,a)') '#   ',ng_max_ok,' groups covering ',max_ok_gr_targs,' targets at ',min_pointings_in_group,' or more pointings were found' 
      if (ng_min_ok.gt.0) write(11,'(a,i0,a,i0,a,i0,a)') '#   ',ng_all_ok,' groups covering ',min_ok_gr_targs,' or more targets at ',min_pointings_in_group,' or more pointings were found' 
      write(11,'(a)') '#'
     endif
          
     
      write(11,'(a)') '# List of prime candidate optimal pointings automatically generated by the IPA module:'
      write(11,'(a)') '#     (selected as all pointings achieving the highest possible simultaneos priority '
      write(11,'(a)') '#     class 1 target coverage and belonging and groups of pointings covering the same '
      write(11,'(a)') '#     targets at n_dither or more pointings)'
      write(11,'(a)') '#'
      write(11,'(a)') '#  Consult ./trial_'//rn(1:2)//'/ipa_output/grouped_optimal_poinings.txt for lists of'
      write(11,'(a)') '#      which priority 1 targets are covered by the individual pointings of each group'
      write(11,'(a)') '#'
      
      endif !group_output=.true

      write(11,'(a)') '#IL# -- Marker. Do not delete or move this line'
      
      write(11,'(a)') '#'

 
      write(11,'(a)') '# Pointing'
      write(11,'(a)') '#  No. Coverage   RA [deg]       Dec [deg]     pa_ap [deg]   Group No.'

       if (group_output) then

       allocate (p_free(n_p))
       
       p_free=.true.

       do ig=1,n_groups ! add superset groupings first      


       if (gr_ptcount(gr_indx(ig)).lt.min_pointings_in_group) cycle !too few pointings in group -> skip group
       
 
        do ip=1,n_p !primary group pointings
         
          if (gr_flag(ip).eq.gr_indx(ig)) then
          write(11,'(a,i5,i6,2x,3f15.9,7x,i0)') ' ',ip,n_targs_p(ip),ra_p(ip),dec_p(ip),pa_ap_p(ip),ig
           p_free(ip)=.false.
         endif
       
        enddo !ip


       if (gr_ptcount(gr_indx(ig)).gt.gr_pcount(gr_indx(ig))) then !pointings of super-set groups need to be added to group
       
        do jg=1,n_groups
       
         if (gr_subs(gr_indx(ig),gr_indx(jg)).gt.0) then !super-set gr_indx(jg) needs to be included
         
           do ip=1,n_p !finde pointings belonging to super-set gr_indx(jg) and print out 

             if (gr_flag(ip).eq.gr_indx(jg)) then
               
              if (p_free(ip)) then
              write(11,'(a,i5,i6,2x,3f15.9,7x,i0)') ' ',ip,n_targs_p(ip),ra_p(ip),dec_p(ip),pa_ap_p(ip),jg
              p_free(ip)=.false.
              else
              write(11,'(a,i5,i6,2x,3f15.9,7x,i0)') '#',ip,n_targs_p(ip),ra_p(ip),dec_p(ip),pa_ap_p(ip),jg
              endif !p_free
        
            endif !super-set jq

          enddo !ip
         
         endif !super-set group pointings need to be output
        
        enddo !jg
        
         endif 
                     
         write(11,'(a)') '#'
      
       enddo !ig

       else
       
        ig=1
        do ip=1,n_p !primary group pointings       
          write(11,'(a,i5,i6,2x,3f15.9,7x,i0)') ' ',ip,n_targs_p(ip),ra_p(ip),dec_p(ip),pa_ap_p(ip),ig      
        enddo !ip
        
        endif !group_output=.true.     

      write(11,'(a)') '#'
 
      write(11,'(a)') '#IE# -- Marker. Do not delete or move this line'
      
      write(11,'(a)') '#'

 
      write(11,'(a)') '#'
      write(11,'(a)') '# CAUTION: The above automatically selected candidate pointing list is not necessarily optimal depending on the task at hand'
      write(11,'(a)') '#'
      write(11,'(a)') '#             Consult the complete IPA output in ./trial_'//rn(1:2)//'/ipa_output/'
      write(11,'(a)') '#'
      write(11,'(a)') '# ----- END OF SECTION APPENDED BY THE IPA MODULE ----------'
      write(11,'(a)') '#'
      write(11,'(a)') '#  -- Modifying the above list of optimal pointings requires re-running the k_make and subsequent modules ---'
2222  write(11,'(a)') '#'
      write(11,'(a)') '#KM# -- Marker do not delete or move this line'
 
      close(11)

! done writing out to configuration, file close up shop.
     
       print *
       call cpu_time(t_stop)
       write(*,'(a)')      '----------------------------------------------'
       write(*,'(a,f7.1)') ' CPU time used by ipa module [s]:',t_stop-t_start
       write(*,'(a)')      '----------------------------------------------'
       write(*,*)
       
       end program ipa

!------------------------------------------------------------
 
include './reference_files/front_transforms_updt.f90'

include './reference_files/back_transforms_updt.f90'

include './reference_files/shutter_routines_new.f90'

include './reference_files/misc_routines_new.f90'

!------------------------------------------------------------

      subroutine file_entry_ipa()
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
      if ((raccx.gt.0.5).or.(raccy.gt.0.5)) then
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

      imname=''
      segname=''

      call skip_hashes(9)
      read(9,*) !ref image name
      
      call skip_hashes(9)
      read(9,*) !seg map name
      

      call skip_hashes(9)
      read(9,*) cra,cdec,cpa_ap

      if ((cpa_ap.lt.0.).or.(cpa_ap.gt.360.)) then
      write(*,*)
      write(*,*) 'Error: Nominal Roll Angle out of range in configuration file'
      write(*,*)
      stop
      endif

      cpa_v3=cpa_ap-(180.-phi)
      if (cpa_v3.lt.0.) cpa_v3=cpa_v3+360.
      if (cpa_v3.gt.360.) cpa_v3=cpa_v3-360.

! work out pointing reference point coordinates in V2,V3 on-sky systems

      call msa_to_sky_bd(xm_ref,ym_ref,ax_refv2,ay_refv3) !v2,V3


      if ((cra.lt.0.).or.(cra.gt.360.)) then
      write(*,*)
      write(*,*) 'Error: Nominal Right Ascension out of range in configuration file'
      write(*,*)
      stop
      endif

      if ((cdec.lt.-90.).or.(cdec.gt.90.)) then
      write(*,*)
      write(*,*) 'Error: Nominal Declination out of range in configuration file'
      write(*,*)
      stop
      endif

      call skip_hashes(9)

      read(9,*) szone_x,szone_y
      if ((szone_x.le.0).or.(szone_y.le.0.)) then
      write(*,*)
      write(*,*) 'Error: Search Zone specification negative or zero in configuration file'
      write(*,*)
      stop
      endif


      call skip_hashes(9)
      read(9,*) angle_to_target
      call dva_calc(angle_to_target,dva_mag)


      call skip_hashes(9)
      read(9,*) !max_c_score1 !

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

      end subroutine file_entry_ipa

!-------------------------------------------------------------

      subroutine cat_pri_one_read()

! Routine to read target catalog

      use target_cat
      use config_parameters, only : catfile
      implicit none
      integer i,idum,pri
      real*8 xdum,ydum

! Open target catalog

      open(9,file=catfile,status='old',access='sequential', form='formatted')

      call skip_hashes(9)

! Determine number of priority 1 targets and number of priority classes in catalog
 
      nclass=0
      ntarg=0
      ntarg_1=0
10    read(9,*,end=20) idum,xdum,ydum,pri
      ntarg=ntarg+1
      if (pri.eq.1) ntarg_1=ntarg_1+1
      nclass=max(nclass,pri)
!       print *,idum,nclass,pri
      goto 10
20    close(9)

      allocate (id_t(ntarg_1),ra_t(ntarg_1),dec_t(ntarg_1),pri_t(ntarg_1))

! Open target catalog again

      open(9,file=catfile,status='old',access='sequential', form='formatted')

      call skip_hashes(9)
      
      i=0
30    read(9,*,end=40) idum,xdum,ydum,pri

      if (pri.eq.1) then
       i=i+1
       id_t(i)=idum
       ra_t(i)=xdum
       dec_t(i)=ydum
       pri_t(i)=pri
      endif

      goto 30
40    close(9)

      return

      end subroutine cat_pri_one_read

!-----------------------------------------------------------------
      
      subroutine draw_cat_on_sky2(sky_plot)
      use derived_parameters, only : ax_refv2,ay_refv3
      use config_parameters
      use target_cat
      implicit none
!       integer ntarg,nclass
!       real*8 rat(ntarg),dect(ntarg)
!       integer pri(ntarg),id(ntarg)
!       real cpa_ap,cpa_v3
!       real*8 cra,cdec
      real*8 ra_min,ra_max,dec_min,dec_max
      real ax_refra_p,ay_refdec_p
      real ax,ay
      real thetax_ref,thetay_ref,ax_0,ay_0
      integer i,ic
      real xx,yy,xp1,xp2,yp1,yp2
      real x1(5),y1(5),x2(5),y2(5),x3(5),y3(5),x4(5),y4(5)
      real x1r(5),y1r(5),x2r(5),y2r(5),x3r(5),y3r(5),x4r(5),y4r(5)
      real r0,twopi,ang,xxs(50),yys(50)
      character*70 sky_plot
!       character*70 catname
      character*2 rn
      common/run/rn
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar

! work out pointing reference point for present roll angle

      call roll_v3(cpa_v3,ax_refv2,ay_refv3,ax_refra_p,ay_refdec_p) !RA,Dec


! Input full target catalog

       deallocate (id_t,ra_t,dec_t,pri_t) ! from previous pc=1 only read

       call catalog_read()

       ra_max=-999.
       ra_min=999.
       dec_max=-999.
       dec_min=999.


!        do i=1,ntarg
!        call tangcoord_single(cra,cdec,ra_t(i),dec_t(i),ax,ay)
!        ra_max=max(ra_max,ra_t(i))
!        ra_min=min(ra_min,ra_t(i))
!        dec_max=max(dec_max,dec_t(i))
!        dec_min=min(dec_min,dec_t(i))
!        enddo

       do i=1,ntarg
       call tangcoord_single(cra,cdec,ra_t(i),dec_t(i),ax,ay)
       ra_max=max(ra_max,ax)
       ra_min=min(ra_min,ax)
       dec_max=max(dec_max,ay)
       dec_min=min(dec_min,ay)
       enddo



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
        call roll_v3(cpa_v3,x1(i),y1(i),x1r(i),y1r(i)) !ra,dec
        call roll_v3(cpa_v3,x2(i),y2(i),x2r(i),y2r(i))
        call roll_v3(cpa_v3,x3(i),y3(i),x3r(i),y3r(i))
        call roll_v3(cpa_v3,x4(i),y4(i),x4r(i),y4r(i))
        enddo


! Plot target catalog and MSA FOV and on sky

        CALL PGBEGIN(0,sky_plot,1,1)
!        CALL PGBEGIN(0,'?',1,1)
        CALL PGSLS(1)
        CALL PGSLW(1)
        CALL PGSCF(2)
        CALL PGSCH(1.)



       xp1=ra_min
       xp2=ra_max
       yp1=dec_min
       yp2=dec_max
       
        xp1=xp1*1.25
        xp2=xp2*1.25
        yp1=yp1*1.25
        yp2=yp2*1.25


      xp1=xp1+ax_refra_p
      xp2=xp2+ax_refra_p
      yp1=yp1+ay_refdec_p
      yp2=yp2+ay_refdec_p


!        print *,'ntarg',ntarg
!        print *,'ax_refra_p,ay_refdec_p',ax_refra_p,ay_refdec_p
!        print *,'xp1,xp2,yp1,yp2',xp1,xp2,yp1,yp2



        CALL PGVPORT(.0,1.,.0,1.)
        CALL PGWNAD(xp2,xp1,yp1,yp2) !flipped x-direction
        CALL PGSLW(1)

      call set_rainbow_lut

      call pgsch(0.2)
      call pgslw(2)

! priority class 0 targets first

      ic=0
      call set_rainbow_color_fine(nclass,ic)
      do i=1,ntarg
      if (pri_t(i).eq.ic) then
      call tangcoord_single(cra,cdec,ra_t(i),dec_t(i),ax,ay)
      ax=(ax*dva_mag+ax_refra_p)
      ay=(ay*dva_mag+ay_refdec_p)
      call pgpoint(1,ax,ay,-12)
      endif
      enddo

! other positive priority class targets

      call pgsch(0.2)
      call pgslw(2)

      do ic=nclass,1,-1
      call set_rainbow_color_fine(nclass,ic)
      do i=1,ntarg
      if (pri_t(i).eq.ic) then
      call tangcoord_single(cra,cdec,ra_t(i),dec_t(i),ax,ay)
      ax=(ax*dva_mag+ax_refra_p)
      ay=(ay*dva_mag+ay_refdec_p)
      call pgpoint(1,ax,ay,-12)
      endif
      enddo
      enddo


! zero and negative priority class targets

      call pgsch(0.2)
      call pgslw(2)

      do ic=nclass,1,-1
      call pgsci(15)
      do i=1,ntarg
      if (pri_t(i).le.0) then
      call tangcoord_single(cra,cdec,ra_t(i),dec_t(i),ax,ay)
      ax=(ax*dva_mag+ax_refra_p)
      ay=(ay*dva_mag+ay_refdec_p)
      call pgpoint(1,ax,ay,-12)
!     call pgsci(1)
!     call pgnumb(id(i_nss(i)),0,1,tnum,ntnum) ! label target ID
!     call pgnumb(i,0,1,tnum,ntnum)            ! label output target number
!     call pgptext(ax,ay+2.5,0.,0.5,tnum(1:ntnum))
      endif
      enddo
      enddo



! Draw MSA quadrant outlines

      call pgsch(1.0)
      call pgscf(2)
      call pgslw(2)

      call pgsci(1)
      call pgline(5,x1r,y1r)
      xx=(x1r(1)+x1r(3))/2.
      yy=(y1r(1)+y1r(3))/2.
      call pgslw(2)
      call pgptext(xx,yy,0.,0.5,'Q1')
      call pgslw(2)

      call pgline(5,x2r,y2r)
      xx=(x2r(1)+x2r(3))/2.
      yy=(y2r(1)+y2r(3))/2.
      call pgslw(2)
      call pgptext(xx,yy,0.,0.5,'Q2')
      call pgslw(2)

      call pgline(5,x3r,y3r)
      xx=(x3r(1)+x3r(3))/2.
      yy=(y3r(1)+y3r(3))/2.
      call pgslw(2)
      call pgptext(xx,yy,0.,0.5,'Q3')
      call pgslw(2)

      call pgline(5,x4r,y4r)
      xx=(x4r(1)+x4r(3))/2.
      yy=(y4r(1)+y4r(3))/2.
      call pgslw(2)
      call pgptext(xx,yy,0.,0.5,'Q4')
      call pgslw(2)

      call pgsci(1)

! draw fixed slits

      call draw_slits_sky_rot(cpa_v3)

      call draw_v2v3(cpa_v3,xp2-0.04*(xp2-xp1),yp1+0.2*(yp2-yp1))

      call draw_compass(0.,xp2-0.04*(xp2-xp1),yp2-0.2*(yp2-yp1))

! draw over-booking zone

         r0=180.
         twopi=6.283185
         do i=1,49
         ang=twopi*float(i)/49.
         xxs(i)=ax_refra_p+r0*cos(ang)
         yys(i)=ay_refdec_p+r0*sin(ang)
         enddo
         xxs(50)=xxs(1)
         yys(50)=yys(1)

         call pgslw(1)
         call pgsci(1)
         call pgline(50,xxs,yys)
         call pgslw(1)

      call pgend

      return

      end subroutine draw_cat_on_sky2

!--------------------------------------------!---------------------------

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

!-----------------------------------------------------------------------

      subroutine draw_compass(roll,x0,y0) !x flipped version!
      implicit none
      real roll !deg
      real x0,y0,xx(2),yy(2),len,roll_r
      len=x0/12.
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

      subroutine draw_v2v3(roll,x0,y0) !x flipped version!
      implicit none
      real roll !deg
      real x0,y0,xx(2),yy(2),len,roll_r
      len=x0/12.
      roll_r=roll/57.2957795
      xx(1)=x0
      yy(1)=y0
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

!-------------------------------------------------------------------------

      subroutine plot_intial_swarms(nt,tx,ty,id_t,ns,sx,sy,szone_x,szone_y,swarm_plot)
        use derived_parameters, only : ax_refv2,ay_refv3,phi_r
        implicit none
        integer nt,ns
        integer id_t(nt)
        real tx(nt),ty(nt),sx(ns),sy(ns),szone_x,szone_y
        character*70 swarm_plot
        integer i
        character*10 targ_id
        integer ntnum
        real xx(5),yy(5)
        real xm_refv2v3,ym_refv2v3

! work out angular v2v3 system origin in MSA orientation
        
      xm_refv2v3= ax_refv2*cos(phi_r)+ay_refv3*sin(phi_r) 
      ym_refv2v3=-ax_refv2*sin(phi_r)+ay_refv3*cos(phi_r)
      
!       print *
!       print *,'ax_refv2,ay_refv3:',ax_refv2,ay_refv3
!  
!       print *,'xm_refv2v3,ym_refv3:',xm_refv2v3,ym_refv2v3
!       print *  
      
      xm_refv2v3=xm_refv2v3*1000. !mas
      ym_refv2v3=ym_refv2v3*1000.
                          
        CALL PGBEGIN(0,swarm_plot,1,1)
        CALL PGSLS(1)
        CALL PGSLW(1)
        CALL PGSCF(2)
        CALL PGSCH(1.)


        CALL PGVPORT(.0,1.,.0,1.)
        call PGADVANCE
        CALL PGWNAD(150000.+xm_refv2v3,-150000.+xm_refv2v3,-150000.+ym_refv2v3,150000.+ym_refv2v3)
        CALL PGSLW(2)

        call pgscr(17,0./256.,206./256.,38./256)
        call pgsci(17)
        call pgsch(0.01)

        call pgpoint(ns,sx,sy,2)

         call pgsch(.5)

         do i=1,nt
         call pgsci(2)
         call pgpoint(1,tx(i),ty(i),18)
         xx(1)=tx(i)-szone_x*1000.
         yy(1)=ty(i)+szone_y*1000.
         xx(2)=tx(i)+szone_x*1000.
         yy(2)=yy(1)
         xx(3)=xx(2)
         yy(3)=ty(i)-szone_y*1000.
         xx(4)=xx(1)
         yy(4)=yy(3)
         xx(5)=xx(1)
         yy(5)=yy(1)
         call pgsci(4)
         call pgline(5,xx,yy)
         call pgsci(1)
         call pgnumb(id_t(i),0,1,targ_id,ntnum)
         call pgptext(tx(i),ty(i)+4000.,0.0,0.5,targ_id)
         enddo

        call pgend
        end subroutine plot_intial_swarms

!-------------------------------------------------------------------------

      subroutine plot_digital_shift_map(nx,ny,smap,dx,dy,max_peak,map_plot)
      implicit none
      integer nx,ny
      integer smap(-nx:nx,-ny:ny)
      real dx,dy
      integer max_peak
      character*70 map_plot
      real :: pmap(2*nx+1,2*ny+1)  
      real tr(6)
      integer i,j,ii,jj
    

        tr(1)=0.
        tr(2)=dx
        tr(3)=0.
        tr(4)=0.
        tr(5)=0.
        tr(6)=dy

        ii=0
        jj=0
        do i=-nx,nx
        ii=i+nx+1
        do j=-ny,ny
        jj=j+ny+1
        pmap(ii,jj)=float(smap(i,j))
        enddo
        enddo
       

        CALL PGBEGIN(0,map_plot,1,1)
        CALL PGSLS(1)
        CALL PGSLW(1)
        CALL PGSCF(2)
        CALL PGSCH(1.)

        CALL PGVPORT(.0,1.,0.05,1.)
        call PGADVANCE
        CALL PGWNAD(float(2*nx+1)*dx,0.,0.,float(2*ny+1)*dy)

      CALL PALETT(2,1.0,0.5)

        call pgimag(pmap,2*nx+1,2*ny+1,1,2*nx+1,1,2*ny+1,0.,float(max_peak),tr)

      CALL PGSCH(1.0)
      call pgscf(2)
      CALL PGWEDG('BI',0.5,2.0,0.,float(max_peak),' ')

        call pgend
        
        return
         end subroutine plot_digital_shift_map

!-------------------------------------------------------------------------
   
      SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
!-----------------------------------------------------------------------
! Set a "palette" of colors in the range of color indices used by
! PGIMAG.
!-----------------------------------------------------------------------
      INTEGER TYPE
      REAL CONTRA, BRIGHT
!
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)
!
      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/
!
      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/
!
      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/
!
      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/
!
      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
!
      IF (TYPE.EQ.1) THEN
!        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
!        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
!        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
!        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
!        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
        END SUBROUTINE PALETT
        
!-----------------------------------------------------------------------

     subroutine fine_tune_pointing_roll(ra,dec,pa_v3,n_list,id_list,max_peak)
     use target_cat, only : ntarg_1
     implicit none
     integer max_peak
     real*8 ra,dec
     real pa_v3
     integer n_list
     integer id_list(2*max_peak)

     integer nt    
!      integer idt(ntarg_1),kt(ntarg_1),it(ntarg_1),jt(ntarg_1)      
!      real rxt(ntarg_1),ryt(ntarg_1)
!      integer nr,idr(ntarg_1)

     integer idt(2*max_peak),kt(2*max_peak),it(2*max_peak),jt(2*max_peak)      
     real rxt(2*max_peak),ryt(2*max_peak)
     integer nr,idr(2*max_peak)
     
     real*8 ra_trim,dec_trim
     real ra_shift,dec_shift,pa_v3_trim
     
     real ra_tol,dec_tol
     
     integer i, n_max

     ra_tol= 5./1000. !as
     dec_tol=5./1000. !as
     n_max=5


! First iteration starting from digital peak RA,Dec

     call project_targets_roll(ra,dec,pa_v3,nt,idt,kt,it,jt,rxt,ryt) 
          
     call purge_overlaps(nt,idt,kt,it,jt,rxt,ryt,nr,idr)
     
     call fine_tune_bx_roll(ra,dec,pa_v3,nt,rxt,ryt,ra_trim,dec_trim,pa_v3_trim,ra_shift,dec_shift)  

     ra=ra_trim
     dec=dec_trim

! Continue iteration until shift correction below ra_tol and dec_tol

     do i=1,n_max

     call project_targets_roll(ra,dec,pa_v3,nt,idt,kt,it,jt,rxt,ryt) 
   
     call purge_overlaps(nt,idt,kt,it,jt,rxt,ryt,nr,idr)
          
     call fine_tune_bx_roll(ra,dec,pa_v3,nt,rxt,ryt,ra_trim,dec_trim,pa_v3_trim,ra_shift,dec_shift)  
 
     ra=ra_trim
     dec=dec_trim
     pa_v3=pa_v3_trim

     if ((ra_shift.lt.ra_tol).and.(dec_shift.lt.dec_tol)) exit
 
     enddo
 
!      if (i.eq.n_max) print *,'Max target group centering iterations exceeded'
         
     n_list=nt

! Load target IDs of covered targets
     
     do i=1,nt
      id_list(i)=idt(i)
     enddo
         
     return
     
     end subroutine fine_tune_pointing_roll

!-----------------------------------------------------------------------

      subroutine project_targets_roll(pra,pdec,ppa_v3,n_out,id_out,k_out,i_out,j_out,rx_out,ry_out)
      use derived_parameters
      use config_parameters
      use target_cat
      implicit none
      integer n_out,n_in
      real*8 pra,pdec
      real ppa_v3
      real ax_refra_p,ay_refdec_p
      integer id_out(ntarg_1),k_out(ntarg_1),i_out(ntarg_1),j_out(ntarg_1)      
      real rx_out(ntarg_1),ry_out(ntarg_1)
      integer slitmap(4,365,171)
      common/slit_map/slitmap
      real ax,ay,aax,aay,xm,ym,rx,ry      
      integer i,kk,ii,jj
      
      n_in=ntarg_1
      n_out=0

! work out pointing reference point for present roll angle

      call roll_v3(ppa_v3,ax_refv2,ay_refv3,ax_refra_p,ay_refdec_p) !RA,Dec
      
      do i=1,n_in
       
       call tangcoord_single(pra,pdec,ra_t(i),dec_t(i),ax,ay)  !ra,dec
       ax=ax*dva_mag+ax_refra_p
       ay=ay*dva_mag+ay_refdec_p
       
       call roll_v3(-ppa_v3,ax,ay,aax,aay) !v2,v3
       call sky_to_msa_bd(aax,aay,xm,ym)
       call msa_to_shutter_bd(xm,ym,kk,ii,jj,rx,ry)
        if (kk.ne.0) then !in shutter
        if (slitmap(kk,ii,jj).eq.4) then !in viable slitlet
        if ((abs(rx).le.raccx).and.(abs(ry).le.raccy)) then !in acceptance

         n_out=n_out+1   !objects in viable slitlets
         id_out(n_out)=id_t(i)
         k_out(n_out)=kk
         i_out(n_out)=ii
         j_out(n_out)=jj
         rx_out(n_out)=rx
         ry_out(n_out)=ry

        endif
        endif
        endif

       enddo

       return
       end subroutine project_targets_roll
       
!-------------------------------------------------------------------------

      subroutine arribas_ipa_test(nk,kt,it,jt,nm,i_m)

! determine nm non-overlapping targets from nk input
! This IPA version does not update viable slitlet map
! NB: shuttervalues already loaded in call to initialize_slitmap
      use derived_parameters,only : avxgap
      use config_parameters, only : n_disp,sthresh
      use target_cat, only : ntarg_1
      implicit none
      integer nk
      integer kt(ntarg_1),it(ntarg_1),jt(ntarg_1)
      integer nm,i_m(ntarg_1)
      integer ovlap(nk,nk),tovlap(nk),max_targ(nk)
      real shv_i,shv_j
      integer n_purge,max_tovlap,n_max,idel,idum
      integer k,i,j,k2,ii
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


      return
      end subroutine arribas_ipa_test

!------------------------------------------------------------------------
  
      subroutine fine_tune_bx_roll(cra_in,cdec_in,cpa_v3_in,nt,rxt,ryt,cra_out,cdec_out,cpa_v3_out,ra_shift,dec_shift)     
      use config_parameters, only : raccx,raccy,cpa_v3
      use derived_parameters, only : av_xpitch,av_ypitch,phi_r
      implicit none
      real*8 cra_in,cdec_in,cra_out,cdec_out
      real cpa_v3_in,cpa_v3_out,roll_corr
      integer nt,it
      real rxt(nt),ryt(nt)
      real max_rx,min_rx,max_ry,min_ry
      real rx_shift,ry_shift,v2_shift,v3_shift,ra_shift,dec_shift

      max_rx=-666.
      min_rx=666.
      max_ry=-666.
      min_ry=666.
            
      do it=1,nt      
            max_rx=max(max_rx,rxt(it))
            min_rx=min(min_rx,rxt(it))
            max_ry=max(max_ry,ryt(it))
            min_ry=min(min_ry,ryt(it))      
      enddo !it

      rx_shift=-(max_rx+min_rx)/2.*av_xpitch/1000. !as on MSA
      ry_shift=(max_ry+min_ry)/2.*av_ypitch/1000. !as on MSA

! rotate shift vector from MSA orientation to V2,V3 orientation

      v2_shift= rx_shift*cos(-phi_r)+ry_shift*sin(-phi_r) 
      v3_shift=-rx_shift*sin(-phi_r)+ry_shift*cos(-phi_r)

! rotate shift vector from  V2,V3 orientation to RA,Dec orientation
      
      call roll_v3(cpa_v3_in,v2_shift,v3_shift,ra_shift,dec_shift) !V2,V3

! calculate updated RA,Dec
      
      call coordtang_single(cra_in,cdec_in,ra_shift,dec_shift,cra_out,cdec_out)
      
      call get_mpt_roll_offset(cra_out,cdec_out,roll_corr)
      
      cpa_v3_out=cpa_v3+roll_corr
           
      return
      
      end subroutine fine_tune_bx_roll

!------------------------------------------------------------------------

      subroutine purge_overlaps(n,id,k,i,j,rx,ry,n_o,id_o)
      use config_parameters
      use target_cat, only : ntarg_1
      implicit none
      integer n
      integer id(ntarg_1),k(ntarg_1),i(ntarg_1),j(ntarg_1)      
      real rx(ntarg_1),ry(ntarg_1)
      integer id_tmp(ntarg_1),k_tmp(ntarg_1),i_tmp(ntarg_1),j_tmp(ntarg_1)      
      real rx_tmp(ntarg_1),ry_tmp(ntarg_1)
      integer nm,i_m(ntarg_1),n_o,id_o(ntarg_1)
      integer m,i_i,i_o

!      call arribas_ipa(n,k,i,j,nm,i_m)
      call arribas_ipa_test(n,k,i,j,nm,i_m)
      
      n_o=0
      do i_i=1,n
       do i_o=1,nm
       if (id(i_i).eq.id(i_m(i_o))) goto 1000
       enddo
       n_o=n_o+1
       id_o(n_o)=id(i_i)
1000   continue
      enddo
      
  
      do m=1,nm
      id_tmp(m)=id(i_m(m))
      k_tmp(m)=k(i_m(m))
      i_tmp(m)=i(i_m(m))
      j_tmp(m)=j(i_m(m))
      rx_tmp(m)=rx(i_m(m))
      ry_tmp(m)=ry(i_m(m))
      enddo

      do m=1,nm
      id(m)=id_tmp(m)
      k(m)=k_tmp(m)
      i(m)=i_tmp(m)
      j(m)=j_tmp(m)
      rx(m)=rx_tmp(m)
      ry(m)=ry_tmp(m)
      enddo

      n=nm
          
       return
       end subroutine purge_overlaps

!-------------------------------------------------------------------------

      subroutine group_pointings()

! NB: The original IPA pointing numbers after renumbering in main ipa program remain static. 
! the initial arrived at group designations are also maintained but ranked by number of targets in group
! (gr_targs) followed by number of pointings in group (gr_pcount) by means of the indices gr_indx 
! and its inverse gr_inv 

      use config_parameters, only : cpa_ap,n_dither
      use target_cat, only : n_reach      
      use pointings_and_groups
      implicit none
      integer :: ip1,ip2,it1,it2,id1,id2
      integer :: peak_val,n_groups_carry
      integer ii,i,j
      integer, allocatable ::gr_tmp(:)
      
      call max_group_calc(n_reach,max_targs,min_targs,max_groups)  !scope out largest mathematically possible number of groups    

!       print *,' Max groups:',max_groups
      
      allocate (gr_flag(n_p))
      allocate (gr_targs(max_groups),gr_ptcount(max_groups),gr_pcount(max_groups),gr_pcountr(max_groups),gr_tcount(max_groups))
      allocate (gr_tmp(max_groups),gr_indx(max_groups),gr_inv(max_groups))
 

      n_groups=0 ! total number of groups found
      n_groups_carry=0
      gr_flag=0 ! number of group pointing i belongs to
      gr_targs=0 ! pointer to number of IPA pointing identifying group member


      do peak_val=max_targs,min_targs,-1 ! count down simultaneous target coverage

      do ip1=1,n_p ! initially assume each pointing target set ip1 to constitute a group

      if (n_targs_p(ip1).ne.peak_val) goto 888 !ip1 not of peak under consideration
      if (gr_flag(ip1).ne.0) goto 888  !ip1 already assigned a group, next ip1

      n_groups=n_groups+1 ! pointing ip1 constitutes new group 

      gr_targs(n_groups)=ip1 !target list of new group n_groups defined by that of pointing ip

      gr_flag(ip1)=n_groups ! pointing ip belongs to new group being considered

! check if any other pointings have same target list as ip1

      do ip2=1,n_p

      if (ip2.eq.ip1) goto 777 !no need to compare ip1 with self, next i
      if (n_targs_p(ip2).ne.peak_val) goto 777 !ip2 not top peak under consideration
      if (gr_flag(ip2).ne.0) goto 777 !pointing ip2 already classified

! pointing ip2 is a possible match to ip1 - compare target lists in detail

      do it1=1,n_targs_p(ip1) !run through all targets of pointing ip1

      id1=ids_p(ip1,it1)

      do it2=1,n_targs_p(ip2) !run through all targets of pointing ip2

      id2=ids_p(ip2,it2)

      if (id2.eq.id1) goto 666 !id1 is in ip2 target list, next it1

      enddo !it2

      goto 777 !!id1 is not in ip2 target list, next ip2

666   continue

      enddo !it1

! all nt1 and nt2 targets match, ip2 is in ip1 group, next ip2

      gr_flag(ip2)=n_groups ! mark ip2 as belonging to group defined by ip1

777   continue

      enddo !ip2

! id1 not in ip2 target list -> ip2 not in ip1 group

888   continue

      enddo !ip1

      if (peak_val.gt.0) write(*,850) ' Number of different ',peak_val,' simultaneous target groupings: ',n_groups-n_groups_carry
850   format(a,i0,a,i0)

      n_groups_carry=n_groups

      enddo !peak_val
            
      write(*,*)
   
      do i=1,n_groups
       gr_tcount(i)=n_targs_p(gr_targs(i)) !number of targets covered in group i
       gr_pcount(i)=0                      !number of pointings in group i
       do j=1,n_p
        if (gr_flag(j).eq.i) then
        gr_pcount(i)=gr_pcount(i)+1
        endif
       enddo
!      write(*,'(a,i0,a,i0,a,i0,a)') ' Group: ',i,' spans ',gr_tcount(i),' simultaneous targets at ',gr_pcount(i),' pointings'
      enddo


        gr_pcountr=float(gr_pcount) !index routine wants real array


! renumber group numbers so largest coverage ones first

       call indexx(n_groups,gr_pcountr,gr_tmp) !gr_tmp group list ranked per decreasing pointings gr_pcount
       
      deallocate (gr_pcountr)


! rank groups per number of targets covered and then number of pointings in group

       ii=0
       
       do peak_val=max_targs,min_targs,-1
       
       do j=n_groups,1,-1
       if (gr_tcount(gr_tmp(j)).eq.peak_val) then
       ii=ii+1
       gr_indx(ii)=gr_tmp(j)
       endif
       enddo !j
       
       enddo !peak_val
             
       deallocate (gr_tmp)

      do i=1,n_groups
      gr_inv(gr_indx(i))=i
      enddo

      call map_subgroups()

! calculate total number of pointings in group including pointings of groups of which a group is a subset

      do i=1,n_groups      
     
      gr_ptcount(gr_indx(i))=gr_pcount(gr_indx(i))
     
        do j=1,n_groups 
         if (gr_subs(gr_indx(i),gr_indx(j)).ne.0) gr_ptcount(gr_indx(i))=gr_ptcount(gr_indx(i))+gr_pcount(gr_indx(j))
        enddo !j
        
      enddo !i
      
      return
      
      end subroutine group_pointings

!-------------------------------------------------------------------------------\

      subroutine map_subgroups()
      use pointings_and_groups, only : n_groups,gr_indx,gr_tcount,gr_targs,ids_p,gr_subs,gr_inv
      implicit none
      integer i,j,ii,jj
      
      allocate (gr_subs(n_groups,n_groups))
      
      gr_subs=0

        do i=1,n_groups
              
        do j=1,n_groups

        if (i.eq.j) then !group is always a subset of it self
         gr_subs(gr_indx(i),gr_indx(j))=0
         goto 100
        endif		



! is group i as subset of group j?

		if (gr_tcount(gr_indx(j)).lt.gr_tcount(gr_indx(i))) then 
! group J contains fewer targets than group i -> i cannot be a subset of j
		 gr_subs(gr_indx(i),gr_indx(j))=0
		 goto 100 !next j
		endif

! check whether all targets in group i are present in group j

		do ii=1,gr_tcount(gr_indx(i)) !target ii in group i

		  do jj=1,gr_tcount(gr_indx(j)) !target jj in group j

		   if (ids_p(gr_targs(gr_indx(i)),ii).eq.ids_p(gr_targs(gr_indx(j)),jj)) goto 200 !target ii is in group j list next ii

		  enddo !jj

! target ii is not in group j list -> group i is not a subset of group j 

 		   gr_subs(gr_indx(i),gr_indx(j))=0
		   goto 100 !next j

200		enddo !ii

 		   gr_subs(gr_indx(i),gr_indx(j))=1		
		
100		enddo !j
		
		enddo !i
		
       return
       end subroutine map_subgroups

!-------------------------------------------------------------------------------

      subroutine write_groups_to_file(group_pointing_file)
      use pointings_and_groups
      use config_parameters, only : cpa_ap
      implicit none
      character*70 group_pointing_file
      real*8 ra_min,ra_max,dec_min,dec_max,av_ra,av_dec,d_ra,d_dec
      integer i,j,ntot,ip,ig,jg
    
       open(9,file=group_pointing_file,status='replace',access='sequential', form='formatted')

      do ig=1,n_groups
      
      write(9,*)      
      write(9,*)
      write(9,'(a,i0)') 'Group Number: ',ig
      write(9,'(a,i0)') 'Number of common targets in Group: ',gr_tcount(gr_indx(ig))
      write(9,'(a,i0)') 'Number of pointings in Group: ',gr_pcount(gr_indx(ig))
       do jg=1,n_groups
      if (gr_subs(gr_indx(ig),gr_indx(jg)).ne.0) write(9,'(a,i0,a,i0,a,i0,a)') 'Group ',ig,' is a subset of group ',jg,', spanning ',gr_pcount(gr_indx(jg)),' additional pointings'
      enddo !j 
      if (gr_ptcount(gr_indx(ig)).gt.gr_pcount(gr_indx(ig))) write(9,'(a,i0)') 'Total number of pointings in group: ',gr_ptcount(gr_indx(ig))
      
      
      write(9,*) 'Targets in group: '
      write(9,*) '        ID_cat'
       do j=1,n_targs_p(gr_targs(gr_indx(ig)))
        write(9,*) ids_p(gr_targs(gr_indx(ig)),j)
       enddo !j



      write(9,*) 'Pointing'
      write(9,*) '  No. Coverage   RA [deg]       Dec [deg]     pa_ap [deg]'
 

      ra_max=-9999.
      ra_min=+9999.
      dec_max=-9999.
      dec_min=+9999.

       do ip=1,n_p
        if (gr_flag(ip).eq.gr_indx(ig)) then
         write(9,750) ip,n_targs_p(ip),ra_p(ip),dec_p(ip),pa_ap_p(ip)
         ra_min=min(ra_min,ra_p(ip))
         ra_max=max(ra_max,ra_p(ip))
         dec_min=min(dec_min,dec_p(ip))
         dec_max=max(dec_max,dec_p(ip))
        endif
       enddo !ip

       if (gr_ptcount(gr_indx(ig)).gt.gr_pcount(gr_indx(ig))) then !pointings of super-set groups need to be included
       
        do jg=1,n_groups       
         if (gr_subs(gr_indx(ig),gr_indx(jg)).gt.0) then !super-set gr_indx(jg) needs to be included        
           do ip=1,n_p !finde pointings belonging to super-set gr_indx(jg) and print out 
             if (gr_flag(ip).eq.gr_indx(jg)) then 
                 write(9,750) ip,n_targs_p(ip),ra_p(ip),dec_p(ip),pa_ap_p(ip)
                 ra_min=min(ra_min,ra_p(ip))
                 ra_max=max(ra_max,ra_p(ip))
                 dec_min=min(dec_min,dec_p(ip))
                 dec_max=max(dec_max,dec_p(ip))
             endif !super-set jq
            enddo !ip       
         endif !super-set group pointings need to be output       
        enddo !jg
        
        endif !pointings of super-set groups need to be included





      av_ra=(ra_max+ra_min)/2.
      av_dec=(dec_max+dec_min)/2.
      d_ra=(ra_max-ra_min)*cos(av_dec/57.295779)*3600./2.
      d_dec=(dec_max-dec_min)*3600./2.

      write(9,'(a)') 'Enclosing envelope of group:'
      write(9,760) 'RA:  ',av_ra,' [deg] +/- ',d_ra,' [as]'
      write(9,760) 'Dec: ',av_dec,' [deg] +/- ',d_dec,' [as]'

      enddo !ig

750   format(I5,I6,2x,3f15.9)
760   format(a,f15.9,a,f6.1,a)


 
! sanity check

      ntot=0
      do i=1,n_groups
      ntot=ntot+gr_pcount(i)
    enddo


! Sanity check

      write(9,*)
      write(9,*) 'Checksum:',ntot-n_p


      close(9)

      return
      
      end subroutine write_groups_to_file
 
!-------------------------------------------------------------------------------

      subroutine plot_pointing_map(pointing_map_plot,nx,ny,dx,dy)
      use config_parameters, only : cra,cdec,cpa_v3,szone_x,szone_y,n_dither
      use derived_parameters, only : av_xpitch,av_ypitch,phi_r
      use pointings_and_groups
      implicit none
      character*70 pointing_map_plot
      integer nx,ny
      real dx,dy
      integer nnx,nny,m,ii,ip,ig,jg
      real xx(2),yy(2)
      real ax,ay,aax,aay
      real tx_peak,ty_peak
      real xm,ym
      character*6 tnum
      integer ntnum
       integer, allocatable :: plot_order(:),gr_order
!        integer, allocatable :: p_free(:)
       logical, allocatable :: g_free(:)

        CALL PGBEGIN(0,pointing_map_plot,1,1)
        CALL PGSLS(1)
        CALL PGSLW(1)
        CALL PGSCF(1)
        CALL PGSCH(1.)

        CALL PGVPORT(.0,1.,0.05,0.95)
        call PGADVANCE
        CALL PGWNAD(float(nx)*dx,-float(nx)*dx,-float(ny)*dy,float(ny)*dy)
        CALL PGBOX('BC',0.0,0,'BC',0.0,0)

! draw shutter facet grid for scale

      call pgslw(1)
      call pgsci(15)

      nnx=int(szone_x*1000./av_xpitch)+10
      nny=int(szone_y*1000./av_ypitch)+10

      xx(1)=-float(nnx)*av_xpitch
      yy(1)=-float(nny)*av_ypitch
      yy(2)=+float(nny)*av_ypitch

      do ii=-nnx,nnx
      xx(1)=xx(1)+av_xpitch
      xx(2)=xx(1)
      call pgline(2,xx,yy)
      enddo

      yy(1)=-float(nny)*av_ypitch
      xx(1)=-float(nnx)*av_xpitch
      xx(2)=+float(nnx)*av_xpitch

      do ii=-nny,nny
      yy(1)=yy(1)+av_ypitch
      yy(2)=yy(1)
      call pgline(2,xx,yy)
      enddo


       allocate (p_free(n_p))      
       allocate (plot_order(n_groups))
       allocate (g_free(n_groups))



! work out plot sequence of max_ok_gr_targs annd mon_ok_gr_targ groups 

       
       g_free=.true.
       
        gr_order=0
       
        do ig=1,n_groups ! add superset groupings first      

       if (gr_ptcount(gr_indx(ig)).lt.min_pointings_in_group) cycle !too few pointings in group -> skip group

       if (gr_tcount(gr_indx(ig)).ne.max_ok_gr_targs) cycle !not max group of interest  -> skip group
      
       if (g_free(gr_indx(ig))) then
        
          gr_order=gr_order+1
        
          plot_order(gr_indx(ig))=gr_order
        
          g_free(gr_indx(ig))=.false.
        
        endif
        
        enddo !ig


        
        gr_order=0
       
        do ig=1,n_groups ! add superset groupings first      

       if (gr_ptcount(gr_indx(ig)).lt.min_pointings_in_group) cycle !too few pointings in group -> skip group

       if (gr_tcount(gr_indx(ig)).ne.min_ok_gr_targs) cycle !not min group of interest  -> skip group
       
       if (g_free(gr_indx(ig))) then
        
          gr_order=gr_order+1
        
          plot_order(gr_indx(ig))=gr_order
        
          g_free(gr_indx(ig))=.false.
        
        endif
        
        enddo !ig
       

       
!        do ig=1,n_groups
! 
!        if (gr_ptcount(gr_indx(ig)).lt.min_pointings_in_group) cycle !too few pointings in group -> skip group
! 
!        if (gr_tcount(gr_indx(ig)).ne.max_ok_gr_targs) cycle !not min group of interest  -> skip group
!         
!         print *,gr_indx(ig),plot_order(gr_indx(ig))
!        
!        enddo
       

       
!         do ig=1,n_groups
! 
!        if (gr_ptcount(gr_indx(ig)).lt.min_pointings_in_group) cycle !too few pointings in group -> skip group
! 
!        if (gr_tcount(gr_indx(ig)).ne.min_ok_gr_targs) cycle !not min group of interest  -> skip group
!         
!         print *,gr_indx(ig),plot_order(gr_indx(ig))
!        
!        enddo      
       
! plot max_ok_gr_targs groups first
           
       
       p_free=.true.

       do ig=1,n_groups    

       if (gr_ptcount(gr_indx(ig)).lt.min_pointings_in_group) cycle !too few pointings in group -> skip group

       if (gr_tcount(gr_indx(ig)).ne.max_ok_gr_targs) cycle !not group of interest  -> skip group
 
 
         do ip=1,n_p !primary group pointings
         
         if (.not.p_free(ip)) cycle
         
          if (gr_flag(ip).eq.gr_indx(ig)) then
      
            call tangcoord_single(cra,cdec,ra_p(ip),dec_p(ip),ax,ay) 

         aax=ax
         aay=ay

         call roll_v3(-cpa_v3,aax,aay,ax,ay) !V2,V3

         tx_peak= ax*cos(phi_r)+ay*sin(phi_r) !xm,ym as
         ty_peak=-ax*sin(phi_r)+ay*cos(phi_r)

         tx_peak=-tx_peak*1000. !mas
         ty_peak=-ty_peak*1000.


       call set_color_shade_red(ng_max_ok,plot_order(gr_indx(ig)))

         call pgslw(1)
         call pgsch(0.5)
         call pgpoint(1,tx_peak,ty_peak,-16) !circle
         call pgsch(0.15)
         call pgsci(1) !Black in .ps
         call pgnumb(ip,0,1,tnum,ntnum)
         call pgptext(tx_peak,ty_peak-szone_y/330.,0.,0.5,tnum(1:ntnum))


       p_free(ip)=.false.

       endif !superset ig
       
        enddo !ip


       if (gr_ptcount(gr_indx(ig)).gt.gr_pcount(gr_indx(ig))) then !pointings of super-set groups need to be included
       
        do jg=1,n_groups
       
         if (gr_subs(gr_indx(ig),gr_indx(jg)).gt.0) then !super-set gr_indx(jg) needs to be included
         
           do ip=1,n_p !finde pointings belonging to super-set gr_indx(jg) and print out 

         if (.not.p_free(ip)) cycle

             if (gr_flag(ip).eq.gr_indx(jg)) then
               
            call tangcoord_single(cra,cdec,ra_p(ip),dec_p(ip),ax,ay) 

         aax=ax
         aay=ay

         call roll_v3(-cpa_v3,aax,aay,ax,ay) !V2,V3

         tx_peak= ax*cos(phi_r)+ay*sin(phi_r) !xm,ym as
         ty_peak=-ax*sin(phi_r)+ay*cos(phi_r)

         tx_peak=-tx_peak*1000. !mas
         ty_peak=-ty_peak*1000.


       call set_color_shade_red(ng_max_ok,plot_order(gr_indx(ig)))

         call pgslw(1)
         call pgsch(0.5)
         call pgpoint(1,tx_peak,ty_peak,-16) !circle
         call pgsch(0.15)
         call pgsci(1) !Black in .ps
         call pgnumb(ip,0,1,tnum,ntnum)
         call pgptext(tx_peak,ty_peak-szone_y/330.,0.,0.5,tnum(1:ntnum))


       p_free(ip)=.false.

        
            endif !super-set jq

          enddo !ip
         
         endif !super-set group pointings need to be output
        
        enddo !jg
        
         endif 
         
      
       enddo !ig


! repeat for min_ok_gr_targs groups 

      do ig=1,n_groups    

       if (gr_ptcount(gr_indx(ig)).lt.min_pointings_in_group) cycle !too few pointings in group -> skip group

       if (gr_tcount(gr_indx(ig)).ne.min_ok_gr_targs) cycle !not group of interest  -> skip group
 
 
         do ip=1,n_p !primary group pointings
 
          if (.not.p_free(ip)) cycle        
 
          if (gr_flag(ip).eq.gr_indx(ig)) then
      
            call tangcoord_single(cra,cdec,ra_p(ip),dec_p(ip),ax,ay) 

         aax=ax
         aay=ay

         call roll_v3(-cpa_v3,aax,aay,ax,ay) !V2,V3

         tx_peak= ax*cos(phi_r)+ay*sin(phi_r) !xm,ym as
         ty_peak=-ax*sin(phi_r)+ay*cos(phi_r)

         tx_peak=-tx_peak*1000. !mas
         ty_peak=-ty_peak*1000.


       call set_color_shade_green(ng_min_ok,plot_order(gr_indx(ig)))

         call pgslw(1)
         call pgsch(0.5)
         call pgpoint(1,tx_peak,ty_peak,-16) !circle
         call pgsch(0.15)
         call pgsci(1) !Black in .ps
         call pgnumb(ip,0,1,tnum,ntnum)
         call pgptext(tx_peak,ty_peak-szone_y/330.,0.,0.5,tnum(1:ntnum))

       p_free(ip)=.false.

       endif !superset ig
       
        enddo !ip


       if (gr_ptcount(gr_indx(ig)).gt.gr_pcount(gr_indx(ig))) then !pointings of super-set groups need to be included
       
        do jg=1,n_groups
       
         if (gr_subs(gr_indx(ig),gr_indx(jg)).gt.0) then !super-set gr_indx(jg) needs to be included
         
           do ip=1,n_p !finde pointings belonging to super-set gr_indx(jg) and print out 

          if (.not.p_free(ip)) cycle
 
             if (gr_flag(ip).eq.gr_indx(jg)) then
               
            call tangcoord_single(cra,cdec,ra_p(ip),dec_p(ip),ax,ay) 

         aax=ax
         aay=ay

         call roll_v3(-cpa_v3,aax,aay,ax,ay) !V2,V3

         tx_peak= ax*cos(phi_r)+ay*sin(phi_r) !xm,ym as
         ty_peak=-ax*sin(phi_r)+ay*cos(phi_r)

         tx_peak=-tx_peak*1000. !mas
         ty_peak=-ty_peak*1000.


       call set_color_shade_green(ng_min_ok,plot_order(gr_indx(ig)))

         call pgslw(1)
         call pgsch(0.5)
         call pgpoint(1,tx_peak,ty_peak,-16) !circle
         call pgsch(0.15)
         call pgsci(1) !Black in .ps
         call pgnumb(ip,0,1,tnum,ntnum)
         call pgptext(tx_peak,ty_peak-szone_y/330.,0.,0.5,tnum(1:ntnum))


       p_free(ip)=.false.

        
            endif !super-set jq

          enddo !ip
         
         endif !super-set group pointings need to be output
        
        enddo !jg
        
         endif 
         
      
       enddo !ig


      deallocate(p_free)


      call pgsch(1.0)
      call pgscf(2)
      call pgsci(1)
      xm=0.
      ym=szone_y*1000.*1.05
      call pgnumb(max_targs,0,1,tnum,ntnum)
      call pgptext(xm,ym,0.,0.5,'Maximum number of simultaneously observable targets: '//tnum(1:ntnum))

      call pgend
      
      return
      
      end subroutine plot_pointing_map
      
!-------------------------------------------------------

       subroutine set_color_shade_red(n,i)
       implicit none
       integer n,i
       real xr,xg,xb
       integer xr_s,xg_s,xb_s
       integer xr_e,xg_e,xb_e
       
       xr_s=220
       xg_s=20
       xb_s=60

       if (n.eq.1) then
       call pgscr(20,xr_s/255.,xg_s/255.,xb_s/255.)       
       call pgsci(20)
       return
       endif
             
       xr_e=255
       xg_e=250
       xb_e=0
             
       xr=(float(xr_s)+float(i-1)/float(n-1)*(float(xr_e-xr_s)))/255.
       xg=(float(xg_s)+float(i-1)/float(n-1)*(float(xg_e-xg_s)))/255.
       xb=(float(xb_s)+float(i-1)/float(n-1)*(float(xb_e-xb_s)))/255.

       call pgscr(20,xr,xg,xb)       
       call pgsci(20)

       return
       
       end subroutine set_color_shade_red
       
!-------------------------------------------------------

       subroutine set_color_shade_green(n,i)
       implicit none
       integer n,i
       real xr,xg,xb
       integer xr_s,xg_s,xb_s
       integer xr_e,xg_e,xb_e
       
       xr_s=0
       xg_s=139
       xb_s=34

       if (n.eq.1) then
       call pgscr(20,xr_s/255.,xg_s/255.,xb_s/255.)       
       call pgsci(20)
       return
       endif
       
       xr_e=173
       xg_e=255
       xb_e=47
         
       xr=(float(xr_s)+float(i-1)/float(n-1)*(float(xr_e-xr_s)))/255.
       xg=(float(xg_s)+float(i-1)/float(n-1)*(float(xg_e-xg_s)))/255.
       xb=(float(xb_s)+float(i-1)/float(n-1)*(float(xb_e-xb_s)))/255.


       call pgscr(20,xr,xg,xb)       
       call pgsci(20)

       return
       
       end subroutine set_color_shade_green

!-------------------------------------------------------
 
      subroutine max_group_calc(n_reach,max_targs,min_targs,max_groups)
! calculates maximum number of groups of max_targs down to min_targs members drawn from n_reach members  
      implicit none
      integer n_reach,max_targs,min_targs
      integer n_groups,max_groups
      integer k,n,m
      real r,r1,r2,r3
      real, parameter ::  e=2.718281828
      real, parameter :: pi=3.141592654

      if (n_reach.eq.1) then
       n_groups=1
       max_groups=1
       return
      endif

        max_groups=0

        do k=max_targs,min_targs,-1

        n=n_reach      
        m=n     
        r1=log(2.*pi*float(m))/2.+float(m)*log(float(m)/e)
        m=k     
        r2=log(2.*pi*float(m))/2.+float(m)*log(float(m)/e)      
        m=n-k      
        r3=log(2.*pi*float(m))/2.+float(m)*log(float(m)/e)
        if (r1-r2-r3.gt.7.0) then
         n_groups=1000
        else
         n_groups=int(exp(r1-r2-r3))
        endif       
              
        max_groups=max_groups+n_groups
        
        enddo !k
          
      return
      
      end subroutine max_group_calc


!-----------------------------------------------------
        subroutine get_mpt_roll_offset(ra_p,dec_p,roll_corr)
        use target_cat, only : ra_c,dec_c
        implicit none
        real*8 ra_p,dec_p,d_ra,dra,decc,decp,d_roll,x,y,d_roll2
        real*8 bc,bp
        real roll_corr

        dra=(ra_p-ra_c)/57.2957795131
        decc=dec_c/57.2957795131
        decp=dec_p/57.2957795131
         

        x=cos(decp)*sin(dra)/(cos(decc)*sin(decp)-sin(decc)*cos(decp)*cos(dra))

        y=-cos(decc)*sin(dra)/(cos(decp)*sin(decc)-sin(decp)*cos(decc)*cos(dra))


        bc=atan(x)*57.2957795131
        
        bp=atan(y)*57.2957795131

        d_roll=-atan((x-y)/(1.+x*y))*57.2957795131
        
        
        d_roll2=atan((sin(dra)*(sin(decc)+sin(decp)))/(cos(decc)*cos(decp)+cos(dra)*(1.+sin(decc)*sin(decp))))*57.2957795131
        
!         write(*,'(a,f13.7)') 'bearing c [deg]:',bc
!         write(*,'(a,f13.7)') 'bearing p [deg]:',bp 
!         write(*,'(a,3f13.7)')  'Roll diff [deg]:',bp-bc,d_roll,d_roll2
        
        roll_corr=real(d_roll2)
        
        return
        
        end subroutine get_mpt_roll_offset
        
!----------------------------------------------------------------------------

      subroutine get_cat_ref_position()
      use config_parameters, only : catfile
      use target_cat, only : ra_c,dec_c
      implicit none
      integer idum,ncat,i,imed
      real*8,allocatable ::ra(:),dec(:)
      integer, allocatable :: indx(:),idno(:)
      integer ih,ihm,id,idm
      real*8 mmra,mmdec
      real fhs,fds,sgn
        
      open(20,file=catfile,status='old',access='sequential', form='formatted')
      call skip_hashes(20)                 
      ncat=0
10    read(20,*,END=11) idum
      ncat=ncat+1
      goto 10
11    continue      
      close(20)
 
      allocate(ra(ncat),dec(ncat),indx(ncat),idno(ncat))

      open(20,file=catfile,status='old',access='sequential', form='formatted')     
      call skip_hashes(20)                 
      do i=1,ncat
      read(20,*) idno(i),ra(i),dec(i)
      enddo     
      close(20)

      imed=int(float(ncat)/2.+1.)

      call indexdx(ncat,ra,indx)

      ra_c=ra(indx(imed))

      call indexdx(ncat,dec,indx)

      dec_c=dec(indx(imed))

      deallocate(ra,dec,indx,idno)

    
      print *
      write(*,'(a,2f13.7)') 'Median RA,Dec:',ra_c,dec_c

      mmra=ra_c/15.
      ih=int(mmra)
      mmra=(mmra-float(ih))*60.
      ihm=int(mmra)
      fhs=(mmra-float(ihm))*60.

      sgn=sign(1.,real(dec_c)) 
      mmdec=abs(dec_c)     
      id=int(mmdec)
      mmdec=(mmdec-float(id))*60.
      idm=int(mmdec)
      fds=(mmdec-float(idm))*60.
      id=id*int(sgn)

      print *     
      write(*,'(a,i3,x,i2,x,f6.3,2x,i3,x,i2,x,f6.3)') 'Median RA,Dec: ',ih,ihm,fhs,id,idm,fds
      print *     

      return
      
      end subroutine get_cat_ref_position
      
      
!------------------------------------------------------------------------------

      SUBROUTINE indexdx(n,arr,indx)
      INTEGER n,indx(n)
      REAL*8 arr(n)
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
      END SUBROUTINE indexdx

!-------------------------------------------------------------------------------

