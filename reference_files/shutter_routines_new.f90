!
! Collection of MSA-specific subroutines employed by the eMPT
! Version 1.0 19 March 2022 - parameter max_failo increased to 50
!
!-------------------------------------------------------------

      subroutine read_base_msamap

! NB! New revised coding scheme for input esa_msa_map:
!
!     0: Failed closed
!     1: Failed open
!     2: Viable 1-slitlet
!     3: Viable 2-slitlet
!     4: Viable 3-slitlet

      use path_to_ref_files
      implicit none
      character*73 filename
      integer msamap(4,365,171)
      integer k,i,j
      integer nfailo
      integer, parameter :: max_failo=50 !maximum number of failed open shutters
      integer ko(max_failo),io(max_failo),jo(max_failo)
      integer nfailc,nok
      common/msa_map/msamap
      common/failed_open/nfailo,ko,io,jo


      filename=trim(adjustl(ref_path))//'esa_msa_map.dat'

      open(12,file=filename,status='old',access='sequential', form='formatted')
      do k=1,4
        read(12,*) ((msamap(k,i,j),i=1,365),j=1,171)
      enddo
      close(12)

! Gather failed and functional shutters

      nfailc=0
      nfailo=0
      nok=0

      do k=1,4
      do i=1,365
      do j=1,171


! count failed closed

      if (msamap(k,i,j).eq.0) nfailc=nfailc+1

! count and save failed opens

      if (msamap(k,i,j).eq.1) then
        nfailo=nfailo+1
        ko(nfailo)=k
        io(nfailo)=i
        jo(nfailo)=j
      endif

! count operational shutters

      if ((msamap(k,i,j).ge.2).and.(msamap(k,i,j).le.4)) nok=nok+1

      enddo
      enddo
      enddo

!      nfailo=0 ! uncomment to turn off failed open spectral overlap

!     write(*,*)
!     write(*,*) 'Bare MSA Statistics at start:'
!     write(*,*)
!     write(*,*) 'Number of failed closed shutters:',nfailc
!     write(*,*) 'Number of failed open shutters:  ',nfailo
!     write(*,*) 'Number of functional shutters:   ',nok
!     write(*,*) 'Shutter count checksum:',
!    #4*365*171-nfailc-nfailo-nok
!     write(*,*)
!     write(*,*) 'List of failed Open shutters:'
!      do i=1,nfailo
!       write(*,*) ko(i),io(i),jo(i)
!      enddo
!     write(*,*)

      return
      end subroutine read_base_msamap

!---------------------------------------------------------------------------

      subroutine make_slitlet_map

! NB! New revised coding scheme for input esa_msa_map:
!
!     0: Failed closed
!     1: Failed open
!     2: Functional Shutter

! NB! new code for localized slitmap
!
!     0: 0-slitlet (3 vertically adjacent failed open or closed shutters)
!     1: Failed open
!     2: Viable 1-slitlet
!     3: Viable 2-slitlet
!     4: Viable 3-slitlet
!     5: Excluded due to failed open overlap
!     7: Excluded due truncated prism spectra
!    10: Excluded due placed target overlap
!    -1: Commanded open shutters

      implicit none
      integer msamap(4,365,171),slitmap(4,365,171)
      integer k,i,j
      integer nfailo
      integer, parameter :: max_failo=50 !maximum number of failed open shutters
      integer ko(max_failo),io(max_failo),jo(max_failo)
      common/msa_map/msamap
      common/failed_open/nfailo,ko,io,jo
      common/slit_map/slitmap
      integer c0,c1,c2,cs

! create bare shutter map including 3, 2 and 1 slitlets

      do k=1,4
      do i=1,365
      do j=1,171

      if (msamap(k,i,j).eq.1) then !failed open - don't use at all
      slitmap(k,i,j)=1
      goto 100
      endif

      cs=msamap(k,i,j)
      if (cs.gt.1) c0=1 !OK shutter
      if (cs.le.1) c0=0 !Failed closed or open

      if (j+1.gt.171) then
      c1=0
      else
      cs=msamap(k,i,j+1)
      if (cs.gt.1) c1=1 !OK shutter
      if (cs.le.1) c1=0 !Failed closed or open
      endif

      if (j-1.lt.1) then
      c2=0
      else
      cs=msamap(k,i,j-1)
      if (cs.gt.1) c2=1 !OK shutter
      if (cs.le.1) c2=0 !Failed closed or open
      endif


      cs=c0+c1+c2

      if (cs.eq.0) slitmap(k,i,j)=0
      if (cs.ge.1) slitmap(k,i,j)=cs+1

100   continue

      enddo
      enddo
      enddo

!     n0slit=0
!     n1slit=0
!     n2slit=0
!     n3slit=0
!     nfopen=0
!
!     do k=1,4
!     do i=1,365
!     do j=1,171
!
!C count failed open shutters
!
!     if (slitmap(k,i,j).eq.1) nfopen=nfopen+1
!
!C count 0 slitlets
!
!     if (slitmap(k,i,j).eq.0) n0slit=n0slit+1
!
!C count 1 slitlets
!
!     if (slitmap(k,i,j).eq.2) n1slit=n1slit+1
!
!C count 2 slitlets
!
!     if (slitmap(k,i,j).eq.3) n2slit=n2slit+1
!
!C count 3 slitlets
!
!     if (slitmap(k,i,j).eq.4) n3slit=n3slit+1
!
!     enddo
!     enddo
!     enddo
!
!
!     write(*,*)
!     write(*,*) 'Viable slitlet Statistics at start:'
!     write(*,*)
!     write(*,*) 'Number of 3 slitlets:          ',n3slit
!     write(*,*) 'Number of 2 slitlets:          ',n2slit
!     write(*,*) 'Number of 1 slitlets:          ',n1slit
!     write(*,*) 'Number of 0 slitlets:          ',n0slit
!     write(*,*) 'Number of Failed Open Shutters:',nfopen
!     write(*,*) 'Shutter Count Check:',
!    #4*365*171-n0slit-n1slit-n2slit-n3slit-nfopen
!     write(*,*)

      return
      end subroutine make_slitlet_map

!---------------------------------------------------------------------------

      subroutine initialize_slitmap(n_disp,sthresh,trflag,opflag)
!
! set trflag=0 to allow detector-gap truncated prism spectra
! set opflag=0 to allow overlap with spectra from failed open shutters
!
!
! NB! New revised coding scheme for input msa_map:
!
!     0: Failed closed
!     1: Failed open
!     2: Functional Shutter

! NB! new code for localized slitmap
!
!     0: 0-slitlet (3 vertically adjacent failed open or closed shutters)
!     1: Failed open
!     2: Viable 1-slitlet
!     3: Viable 2-slitlet
!     4: Viable 3-slitlet
!     5: Excluded due to failed open overlap
!     7: Excluded due truncated prism spectra
!    10: Excluded due placed target overlap
!    -1: Commanded open shutters

      use derived_parameters, only : avxgap
      implicit none
      real sthresh
      integer trflag,opflag
      integer msamap(4,365,171),slitmap(4,365,171)
      real shval(4,365,171),msa_gap,pix
      integer k,i,j,m
      integer nfailo
      integer, parameter :: max_failo=50 !maximum number of failed open shutters
      integer ko(max_failo),io(max_failo),jo(max_failo)
      common/msa_map/msamap
      common/failed_open/nfailo,ko,io,jo
      common/slit_map/slitmap
      common/shutter_values/shval
      real xg1,yg1,xg2,yg2,xg3,yg3,xg4,yg4,xm,ym,xs1,ys1,xs2,ys2
      integer ntrunc,nkillo,j1,j2,k2,ii,jj,dj
      real i1,i2
! - spectral parameters
      integer n_disp
      real lambda1,lambda2
      character*6 gname
      common/disp/lambda1,lambda2,gname
! - prism seperations
      integer prism_sep_p(4,365,171),prism_sep_m(4,365,171)
      common/prism_length/prism_sep_p,prism_sep_m
      real sep_m,sep_p
      
      call load_derived_parameters
      pix=0.018 !mm

      k2=0 !initialzation to avoid compiler warning


! Exclude viable 1, 2 & 3 slitlets leading to prism spectra truncation


      ntrunc=0

!          goto 666

      if ((n_disp.eq.7).and.(trflag.ne.0)) then

!  calculate detector gap edges (including reference pixels)

      call pixel_to_fpa_bd(1,2044,5,xg1,yg1)
      call pixel_to_fpa_bd(2,5,5,xg2,yg2)
      call pixel_to_fpa_bd(1,2044,2044,xg3,yg3)
      call pixel_to_fpa_bd(2,5,2044,xg4,yg4)

      xg1=(xg1+xg3)/2.
      yg1=(yg1+yg2)/2.
      xg2=(xg2+xg4)/2.
      yg2=(yg3+yg4)/2.

      xg1=xg1+pix/2.
      yg1=yg1+pix/2.
      xg2=xg2-pix/2.
      yg2=yg2-pix/2.



      do k=1,4
       do j=1,171
        do i=1,365
        if (slitmap(k,i,j).ne.7) then
        if ((slitmap(k,i,j).ge.2).and.(slitmap(k,i,j).le.4)) then
        call shutter_to_msa_bd(k,i,j,xm,ym)
        call msa_to_fpa_disp(n_disp,lambda1,xm,ym,xs1,ys1)
        call msa_to_fpa_disp(n_disp,lambda2,xm,ym,xs2,ys2)
        if (.not.((xs2.lt.xg1).or.(xs1.gt.xg2))) then
         slitmap(k,i,j)=7
         ntrunc=ntrunc+1
        endif
        endif
        endif
         enddo
        enddo
       enddo

      endif


666 continue

! Exclude viable 1,2 & 3 slitlets due to failed opens

      nkillo=0

      if (opflag.ne.0) then

      do m=1,nfailo

      k=ko(m)
      i=io(m)
      j=jo(m)
      dj=int(2.*sthresh)
      j1=max(j-dj,1)
      j2=min(j+dj,171)
      if (k.eq.1) k2=3
      if (k.eq.3) k2=1
      if (k.eq.2) k2=4
      if (k.eq.4) k2=2

      do ii=1,365
      do jj=j1,j2

      if (abs(shval(k,i,j)-shval(k,ii,jj)).le.(sthresh-1.))then

       if (n_disp.eq.7) then
       sep_p=float(prism_sep_p(k,i,j))
       sep_m=float(prism_sep_m(k,i,j))
       i1=float(i)+msa_gap(k)*(avxgap+365.)
       i2=float(ii)+msa_gap(k)*(avxgap+365.)
       if (((i2-i1).le.sep_m).or.((i2-i1).ge.sep_p)) goto 123
       endif

       if ((slitmap(k,ii,jj).ge.2).and.(slitmap(k,ii,jj).le.4)) then
       nkillo=nkillo+1
       slitmap(k,ii,jj)=5
       endif

123   continue
      endif

      if (abs(shval(k,i,j)-shval(k2,ii,jj)).lt.(sthresh-1.)) then

       if (n_disp.eq.7) then
       sep_p=float(prism_sep_p(k,i,j))
       sep_m=float(prism_sep_m(k,i,j))
       i1=float(i)+msa_gap(k)*(avxgap+365.)
       i2=float(ii)+msa_gap(k2)*(avxgap+365.)
       if (((i2-i1).le.sep_m).or.((i2-i1).ge.sep_p)) goto 124
       endif

       if ((slitmap(k2,ii,jj).ge.2).and.(slitmap(k2,ii,jj).le.4)) then
       nkillo=nkillo+1
       slitmap(k2,ii,jj)=5
       endif

124   continue
      endif

      enddo
      enddo

      enddo

      endif

      return
      end subroutine initialize_slitmap

!---------------------------------------------------------------------------

      subroutine slitmap_stats(lu)

      implicit none
      integer lu
      integer slitmap(4,365,171)
      common/slit_map/slitmap
      integer ntrunc,nkillo,nfopen,n3slit,n2slit,n1slit,n0slit,ntover,ncopen
      integer i,k,j


      ntrunc=0
      nkillo=0
      nfopen=0
      n3slit=0
      n2slit=0
      n1slit=0
      n0slit=0
      ntover=0
      ncopen=0

      do k=1,4
      do i=1,365
      do j=1,171

! count viable 3-slitlets

      if (slitmap(k,i,j).eq.4) n3slit=n3slit+1

! count viable 2-slitlets

      if (slitmap(k,i,j).eq.3) n2slit=n2slit+1

! count viable 1-slitlets

      if (slitmap(k,i,j).eq.2) n1slit=n1slit+1

! count viable 0-slitlets

      if (slitmap(k,i,j).eq.0) n0slit=n0slit+1

! count failed open

      if (slitmap(k,i,j).eq.1) nfopen=nfopen+1


! recount FO overlap vetoed slitlets after PRISM truncation purge

      if (slitmap(k,i,j).eq.5) nkillo=nkillo+1

! recount PRISM truncation vetoed slitlets

      if (slitmap(k,i,j).eq.7) ntrunc=ntrunc+1

! count accepted targets overlap  vetoed slitlets

      if (slitmap(k,i,j).eq.10) ntover=ntover+1

! count commanded open shutters

      if (slitmap(k,i,j).eq.-1) ncopen=ncopen+1


      enddo
      enddo
      enddo



      write(lu,*) 'Slitlets censored due to FO ovlp:',nkillo
      write(lu,*) 'Slitlets censored by PRISM Trunc:',ntrunc
      write(lu,*) 'Slitlets censored by target ovlp:',ntover
      write(lu,*) 'Shutters commanded open:         ',ncopen
      write(lu,*) 'Remaining viable 3-slitlets:     ',n3slit
      write(lu,*) 'Remaining viable 2-slitlets:     ',n2slit
      write(lu,*) 'Remaining viable 1-slitlets:     ',n1slit
      write(lu,*) 'Number 0-slitlets:               ',n0slit
      write(lu,*) 'Number of failed open shutters:  ',nfopen
      write(lu,*) 'Shutter count checksum:',4*365*171-nkillo-ntrunc-nfopen-n3slit-n2slit-n1slit-n0slit-ntover-ncopen


      return
      end subroutine slitmap_stats

!-------------------------------------------------------------------------

      subroutine read_shutter_values
      use path_to_ref_files
! choice of disperser parsed through common gname string
      implicit none
      character*79 filename
      real shval(4,365,171)
      integer k,i,j
      common/shutter_values/shval
! - spectral parameters
!     integer n_disp
      real lambda1,lambda2
      character*6 gname
      common/disp/lambda1,lambda2,gname

      filename=trim(adjustl(ref_path))//'th_shval_pj_'//gname(1:5)//'.dat'
!       filename=trim(ref_path)//'bf_shval_pf_'//gname(1:5)//'.dat'

      open(09,file=filename,status='old',access='sequential', form='formatted')
      do k=1,4
        read(09,*) ((shval(k,i,j),i=1,365),j=1,171)
      enddo
      close(9)


      return
      end subroutine read_shutter_values

!-------------------------------------------------------------

      subroutine read_prism_separations
      use path_to_ref_files
      implicit none
      character*79 filename
      integer k,i,j
! - prism seperations
      integer prism_sep_p(4,365,171),prism_sep_m(4,365,171)
      common/prism_length/prism_sep_p,prism_sep_m

      filename=trim(adjustl(ref_path))//'prism_sep.dat'

      open(9,file=filename,status='old',access='sequential', form='formatted')
      do k=1,4
        read(9,*) ((prism_sep_p(k,i,j),i=1,365),j=1,171)
      enddo
      do k=1,4
        read(9,*) ((prism_sep_m(k,i,j),i=1,365),j=1,171)
      enddo
      close(9)

      return
      end subroutine read_prism_separations

!-------------------------------------------------------------

      subroutine make_csv_file(filename)
      implicit none
      integer msamap(4,365,171)
      common/msa_map/msamap
      integer slitmap(4,365,171)
      common/slit_map/slitmap
      character*683 line !2*2*171-1  two y-rows
      character*1 shut_stat
      integer ir,jr,ii,jj,kk
      character*70 filename
!     character*2 rn
!     common/run/rn


! NB! Revised updated version that avoids issue of failed closed shutters being designated as
! viable 2-slitlets in slitmask array

! NB! New revised coding scheme for msamap:
!
!     0: Failed closed
!     1: Failed open
!     2: functional


! NB! new code slitmap
!
!     0: Failed closed
!     1: Failed open
!     2: Viable 1-slitlet
!     3: Viable 2-slitlet
!     4: Viable 3-slitlet
!     5: Excluded due to failed open overlap
!     7: Excluded due truncated prism spectra
!    10: Excluded due placed target overlap
!    -1: Commanded open shutters

      open(9,file=filename,status='replace',access='sequential',form='formatted')

      line='# This CSV indicates which shutters should be open/closed on the MSA - created by ESA NIRSpec Team'

      write(9,'(a)') line


! write msa left half q1 and q2


      do ir=1,365

       jj=0
       do jr=1,341,2

        kk=1
        ii=ir
        jj=jj+1

        if (msamap(kk,ii,jj).eq.0) shut_stat='x'
        if (msamap(kk,ii,jj).eq.1) shut_stat='s'
        if (msamap(kk,ii,jj).ge.2) shut_stat='1'
        if (slitmap(kk,ii,jj).eq.-1) shut_stat='0'

        line(jr:jr)=shut_stat

       enddo

       jj=0
       do jr=343,683,2

        kk=2
        ii=ir
        jj=jj+1

        if (msamap(kk,ii,jj).eq.0) shut_stat='x'
        if (msamap(kk,ii,jj).eq.1) shut_stat='s'
        if (msamap(kk,ii,jj).ge.2) shut_stat='1'
        if (slitmap(kk,ii,jj).eq.-1) shut_stat='0'

        line(jr:jr)=shut_stat

       enddo

! write seperation commas

      shut_stat=','
      do jr=2,682,2
      line(jr:jr)=shut_stat
      enddo

      write(9,'(a)') line

      enddo


! write msa right half q3 and q4

      do ir=366,730

       jj=0
       do jr=1,341,2

        kk=3
        ii=ir-365
        jj=jj+1

        if (msamap(kk,ii,jj).eq.0) shut_stat='x'
        if (msamap(kk,ii,jj).eq.1) shut_stat='s'
        if (msamap(kk,ii,jj).ge.2) shut_stat='1'
        if (slitmap(kk,ii,jj).eq.-1) shut_stat='0'

        line(jr:jr)=shut_stat

       enddo

       jj=0
       do jr=343,683,2

        kk=4
        ii=ir-365
        jj=jj+1

        if (msamap(kk,ii,jj).eq.0) shut_stat='x'
        if (msamap(kk,ii,jj).eq.1) shut_stat='s'
        if (msamap(kk,ii,jj).ge.2) shut_stat='1'
        if (slitmap(kk,ii,jj).eq.-1) shut_stat='0'

        line(jr:jr)=shut_stat

       enddo

! write seperation commas

      shut_stat=','
      do jr=2,682,2
      line(jr:jr)=shut_stat
      enddo

      write(9,'(a)') line

      enddo


      close(9)

      return
      end subroutine make_csv_file

!------------------------------------------------------------------------------------------

      subroutine plot_slitmap(filename)
      implicit none
      integer msamap(4,365,171),slitmap(4,365,171)
      integer k,i,j
      integer nfailo
      integer, parameter :: max_failo=50 !maximum number of failed open shutters
      integer ko(max_failo),io(max_failo),jo(max_failo)
      integer nplug
      integer, parameter :: max_plug=20000 !maximum number of plugged shutters
      real xplug(max_plug),yplug(max_plug)
      real x1,x2,y1,y2,xd(5),yd(5),xm,ym
      common/msa_map/msamap
      common/failed_open/nfailo,ko,io,jo
      common/slit_map/slitmap
      character*70 filename
! - spectral parameters
!       real lambda1,lambda2
!       character*6 gname
!       common/disp/lambda1,lambda2,gname
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar

! NB! new code for localized msamap
!
!     0: Failed closed
!     1: Failed open
!     2: Viable 1-slitlet
!     3: Viable 2-slitlet
!     4: Viable 3-slitlet
!     5: Excluded due to failed open overlap
!     7: Excluded due truncated prism spectra
!    10: Excluded due placed target overlap
!    -1: Commanded open shutters


        CALL PGBEGIN(0,filename,1,1)
        CALL PGSLS(1)
        CALL PGSLW(1)
        CALL PGSCF(2)
        CALL PGSCH(1.)

        x1=-55.
        x2=55.
        y1=-50.
        y2=50.

        CALL PGVPORT(.0,1.,.0,1.)
        call PGADVANCE
!        CALL PGWNAD(x1,x2,y1,y2)
        CALL PGWNAD(x2,x1,y2,y1) !Rotated to match STScI convention
        CALL PGSLW(2)

! MSA Quadrant 4

        nplug=0
        do i=1,365
        do j=1,171
         if (msamap(4,i,j).eq.0) then
          nplug=nplug+1
          xplug(nplug)=float(i)
          yplug(nplug)=float(j)
         endif
        enddo
        enddo

       call draw_msa_quad(4,xplug,yplug,nplug)


! MSA Quadrant 2


        nplug=0
        do i=1,365
        do j=1,171
         if (msamap(2,i,j).eq.0) then
          nplug=nplug+1
          xplug(nplug)=float(i)
          yplug(nplug)=float(j)
         endif
        enddo
        enddo

        call draw_msa_quad(2,xplug,yplug,nplug)


! MSA Quadrant 3

        nplug=0
        do i=1,365
        do j=1,171
         if (msamap(3,i,j).eq.0) then
          nplug=nplug+1
          xplug(nplug)=float(i)
          yplug(nplug)=float(j)
         endif
        enddo
        enddo


        call draw_msa_quad(3,xplug,yplug,nplug)


! MSA Quadrant 1

        nplug=0
        do i=1,365
        do j=1,171
         if (msamap(1,i,j).eq.0) then
          nplug=nplug+1
          xplug(nplug)=float(i)
          yplug(nplug)=float(j)
         endif
        enddo
        enddo

        call draw_msa_quad(1,xplug,yplug,nplug)


! Draw fixed slits

      call draw_fixed_slits

! Map out various classes of shutter per code

! Draw failed open shutters (code 1, red)

      call pgsci(2)

      do k=1,4
      do i=1,365
      do j=1,171

      if (msamap(k,i,j).eq.1) then
        call shutter_to_msa_bd(k,i,j,xm,ym)
        xd(1)=xm-0.5*xpitch(k)
        yd(1)=ym-0.5*ypitch(k)
        xd(2)=xd(1)+xpitch(k)
        yd(2)=yd(1)
        xd(3)=xd(2)
        yd(3)=yd(1)+ypitch(k)
        xd(4)=xd(1)
        yd(4)=yd(3)
        xd(5)=xd(1)
        yd(5)=yd(1)
        call pgpoly(5,xd,yd)
       endif

      enddo
      enddo
      enddo


! Draw slitlets containing final targets (code -1, blue)

      call pgsci(4)

      do k=1,4
      do i=1,365
      do j=1,171

      if (slitmap(k,i,j).eq.-1) then
        call shutter_to_msa_bd(k,i,j,xm,ym)
        xd(1)=xm-0.5*xpitch(k)
        yd(1)=ym-0.5*ypitch(k)
        xd(2)=xd(1)+xpitch(k)
        yd(2)=yd(1)
        xd(3)=xd(2)
        yd(3)=yd(1)+ypitch(k)
        xd(4)=xd(1)
        yd(4)=yd(3)
        xd(5)=xd(1)
        yd(5)=yd(1)
        call pgpoly(5,xd,yd)
       endif

      enddo
      enddo
      enddo


! Mark killed slitlets due to failed open shutters (code 5, red)

      call pgsci(2)
      call pgsch(0.01)

      do k=1,4
      do i=1,365
      do j=1,171

      if (slitmap(k,i,j).eq.5) then
      call shutter_to_msa_bd(k,i,j,xm,ym)
      call pgpoint(1,xm,ym,-1)
      endif

      enddo
      enddo
      enddo


! Mark target overlap shutters (code 10, pale orange)

      call pgscr(16,255./256.,165./256.,0./256.)
      call pgsci(16)
      call pgsch(0.01)

      do k=1,4
      do i=1,365
      do j=1,171
       if (slitmap(k,i,j).eq.10) then
       call shutter_to_msa_bd(k,i,j,xm,ym)
       call pgpoint(1,xm,ym,-1)
       endif
      enddo
      enddo
      enddo

! Mark leftover viable 3-shutters (code 4, green)

      call pgscr(17,0./256.,206./256.,38./256)
      call pgsci(17)
      call pgsch(0.01)

      do k=1,4
      do i=1,365
      do j=1,171
       if (slitmap(k,i,j).eq.4) then
       call shutter_to_msa_bd(k,i,j,xm,ym)
       call pgpoint(1,xm,ym,-1)
       endif
      enddo
      enddo
      enddo

! Mark leftover viable 2-shutters (code 3, light green)

      call pgscr(18,151./256.,234./256.,42./256.)
      call pgsci(18)
      call pgsch(0.01)

      do k=1,4
      do i=1,365
      do j=1,171
       if (slitmap(k,i,j).eq.3) then
       call shutter_to_msa_bd(k,i,j,xm,ym)
       call pgpoint(1,xm,ym,-1)
       endif
      enddo
      enddo
      enddo

! Mark leftover viable 1-shutters (code 2, yellow green)

      call pgscr(19,218./256.,240./256.,30./256.)
      call pgsci(19)
      call pgsch(0.01)

      do k=1,4
      do i=1,365
      do j=1,171
       if (slitmap(k,i,j).eq.2) then
       call shutter_to_msa_bd(k,i,j,xm,ym)
       call pgpoint(1,xm,ym,-1)
       endif
      enddo
      enddo
      enddo

! Mark PRISM Truncated shutters (code 7, light grey)

      call pgscr(20,211./256.,211./256.,211./256.)
      call pgsci(20)
      call pgsch(0.01)

      do k=1,4
      do i=1,365
      do j=1,171
       if (slitmap(k,i,j).eq.7) then
       call shutter_to_msa_bd(k,i,j,xm,ym)
       call pgpoint(1,xm,ym,-1)
       endif
      enddo
      enddo
      enddo


      call pgend

      return
      end subroutine plot_slitmap

!--------------------------------------------------------------

      subroutine draw_fixed_slits
      real xx(60),yy(60)
      real x00,y00,r0,twopi,ang,w,h
      integer i
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
!-------------------------

! IFU Aperture


         X00=0.04125475*1000.
         Y00=(6.99E-06)*1000.
         r0=1.023
         twopi=6.283185
         do i=1,59
         ang=twopi*float(i)/59.
         xx(i)=X00+r0*cos(ang)
         yy(i)=Y00+r0*sin(ang)
         enddo
         xx(60)=xx(1)
         yy(60)=yy(1)

         call pgslw(1)
         call pgline(60,xx,yy)
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

         call pgline(5,xx,yy)

        yy(1)=yy(2)
        yy(2)=yy(3)


        xx(1)=xx(4)

        do i=1,29
        xx(1)=xx(1)+2.*0.66/30.
        xx(2)=xx(1)
        call pgline(2,xx,yy)
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

         call pgsci(1)
         call pgline(5,xx,yy)
         call pgsci(1)


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

         call pgsci(1)
         call pgline(5,xx,yy)
         call pgsci(1)



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

         call pgsci(1)
         call pgline(5,xx,yy)
         call pgsci(1)


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

         call pgsci(1)
         call pgline(5,xx,yy)
         call pgsci(1)


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

         call pgsci(1)
         call pgline(5,xx,yy)
         call pgsci(1)

      return
      end subroutine draw_fixed_slits

! --------------------------------------------------------------------

      subroutine draw_msa_quad(nq,xplug,yplug,nplug)
      implicit none
      integer nplug
      real xplug(nplug),yplug(nplug)
      real xx(101),yy(101)
      integer i,j,m,mm,nq
      real x1,y1,x2,y2
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
!-------------------------
      real xgrid,ygrid

      xgrid=(xpitch(nq)-xopen(nq))/2.
      ygrid=(ypitch(nq)-yopen(nq))/2.

      m=101 ! = mm*4+1
      mm=(m-1)/4

! draw shutter grid

        call pgsci(15)
        call pgslw(1)

        x1=0.
        x2=x1+float(mx)*xpitch(nq)

        do i=1,my+1
        y1=float(i-1)*ypitch(nq)-ygrid
        y2=y1+2.*ygrid

        xx(1)=x1
        yy(1)=y1
        do j=2,(m+1)/2
        xx(j)=x1+(x2-x1)/float((m-3)/2)*float(j-2)
        yy(j)=y2
        xx(m+2-j)=xx(j)
        yy(m+2-j)=y1
        enddo


        call rottran(nq,m,xx,yy)
        call pgpoly(m,xx,yy)

        enddo

        y1=0.
        y2=float(my)*ypitch(nq)

        do i=1,mx+1

        x1=float(i-1)*xpitch(nq)-xgrid
        x2=x1+2.*xgrid

        xx(1)=x1
        yy(1)=y1
        do j=2,(m+1)/2
        yy(j)=y1+(y2-y1)/float((m-3)/2)*float(j-2)
        xx(j)=x2
        yy(m+2-j)=yy(j)
        xx(m+2-j)=x1
        enddo

        call rottran(nq,m,xx,yy)
        call pgpoly(m,xx,yy)

        enddo


        call pgsci(1)

! draw edges

         x1=0.
         y1=0.
         x2=x1+float(mx)*xpitch(nq)
         y2=y1+float(my)*ypitch(nq)

         xx(1)=x1
         yy(1)=y1
         do j=2,mm+1
         xx(j)=x1
         yy(j)=y1+(y2-y1)/float(mm)*float(j-1)
         xx(j+mm)=x1+(x2-x1)/float(mm)*float(j-1)
         yy(j+mm)=y2
         xx(j+2*mm)=x2
         yy(j+2*mm)=y2-(y2-y1)/float(mm)*float(j-1)
         xx(j+3*mm)=x2-(x2-x1)/float(mm)*float(j-1)
         yy(j+3*mm)=y1
         enddo

        call rottran(nq,m,xx,yy)
        call pgline(m,xx,yy)


         x1=x1-xbar
         y1=y1-ybar
         x2=x2+xbar
         y2=y2+ybar

         xx(1)=x1
         yy(1)=y1
         do j=2,mm+1
         xx(j)=x1
         yy(j)=y1+(y2-y1)/float(mm)*float(j-1)
         xx(j+mm)=x1+(x2-x1)/float(mm)*float(j-1)
         yy(j+mm)=y2
         xx(j+2*mm)=x2
         yy(j+2*mm)=y2-(y2-y1)/float(mm)*float(j-1)
         xx(j+3*mm)=x2-(x2-x1)/float(mm)*float(j-1)
         yy(j+3*mm)=y1
         enddo

        call rottran(nq,m,xx,yy)
        call pgline(m,xx,yy)

        call pgsci(1)
        call pgslw(1)

! draw plugs

        do i=1,nplug
        xx(1)=(xplug(i)-1.)*xpitch(nq)
        yy(1)=(yplug(i)-1.)*ypitch(nq)
        xx(2)=xx(1)+xpitch(nq)
        yy(2)=yy(1)
        xx(3)=xx(2)
        yy(3)=yy(1)+ypitch(nq)
        xx(4)=xx(1)
        yy(4)=yy(3)
        xx(5)=xx(1)
        yy(5)=yy(1)
        call rottran(nq,m,xx,yy)
        call pgpoly(5,xx,yy)
        enddo



        return
        end subroutine draw_msa_quad

!-----------------------------------------------------

      subroutine draw_slits_on_sky
      real xx(60),yy(60)
      real x00,y00,r0,twopi,ang,w,h
      integer i
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
!-------------------------


! IFU Aperture


         X00=0.04125475*1000.
         Y00=(6.99E-06)*1000.
         r0=1.023
         twopi=6.283185
         do i=1,59
         ang=twopi*float(i)/59.
         xx(i)=X00+r0*cos(ang)
         yy(i)=Y00+r0*sin(ang)
         enddo
         xx(60)=xx(1)
         yy(60)=yy(1)

         call pgslw(1)
         call tran_pgline(60,xx,yy)
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

         call tran_pgline(5,xx,yy)

        yy(1)=yy(2)
        yy(2)=yy(3)


        xx(1)=xx(4)

        do i=1,29
        xx(1)=xx(1)+2.*0.66/30.
        xx(2)=xx(1)
        call tran_pgline(2,xx,yy)
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

         call pgsci(1)
         call tran_pgline(5,xx,yy)
         call pgsci(1)


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

         call pgsci(1)
         call tran_pgline(5,xx,yy)
         call pgsci(1)



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

         call pgsci(1)
         call tran_pgline(5,xx,yy)
         call pgsci(1)


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

         call pgsci(1)
         call tran_pgline(5,xx,yy)
         call pgsci(1)


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

         call pgsci(1)
         call tran_pgline(5,xx,yy)
         call pgsci(1)

      Return
      end subroutine draw_slits_on_sky

!-----------------------------------

      subroutine draw_msa_quad_on_sky(nq,xplug,yplug,nplug)
      implicit none
      integer nplug
      real xplug(nplug),yplug(nplug)
      real xx(101),yy(101)
      integer i,j,m,mm,ii,jj,mq,nq
      real x1,y1,x2,y2
      real map(365,171)
! --- MSA Parameters-------
      integer mx,my
      real x0(4),y0(4),rot(4),xpitch(4),ypitch(4),xopen(4),yopen(4),xbar,ybar
      common/msa_dims/mx,my,x0,y0,rot,xpitch,ypitch,xopen,yopen,xbar,ybar
!-------------------------
      real xgrid,ygrid
      integer nq_target(8000),i_target(8000),j_target(8000),shut_flag
      integer ntarg
      common/targets/ntarg,nq_target,i_target,j_target,shut_flag

      xgrid=(xpitch(nq)-xopen(nq))/2.
      ygrid=(ypitch(nq)-yopen(nq))/2.


      m=101 ! = mm*4+1
      mm=(m-1)/4

! draw shutter grid

        call pgsci(15)
        call pgslw(1)

        x1=0.
        x2=x1+float(mx)*xpitch(nq)

        do i=1,my+1
        y1=float(i-1)*ypitch(nq)-ygrid
        y2=y1+2.*ygrid

        xx(1)=x1
        yy(1)=y1
        do j=2,(m+1)/2
        xx(j)=x1+(x2-x1)/float((m-3)/2)*float(j-2)
        yy(j)=y2
        xx(m+2-j)=xx(j)
        yy(m+2-j)=y1
        enddo

!        do j=1,m
!        write(*,*) j,xx(j),yy(j)
!        enddo

        call rottran(nq,m,xx,yy)
        call tran_pgpoly(m,xx,yy)

        enddo

        y1=0.
        y2=float(my)*ypitch(nq)

        do i=1,mx+1

        x1=float(i-1)*xpitch(nq)-xgrid
        x2=x1+2.*xgrid

        xx(1)=x1
        yy(1)=y1
        do j=2,(m+1)/2
        yy(j)=y1+(y2-y1)/float((m-3)/2)*float(j-2)
        xx(j)=x2
        yy(m+2-j)=yy(j)
        xx(m+2-j)=x1
        enddo

        call rottran(nq,m,xx,yy)
        call tran_pgpoly(m,xx,yy)

        enddo



        call pgsci(1)

! draw edges


         x1=0.
         y1=0.
         x2=x1+float(mx)*xpitch(nq)
         y2=y1+float(my)*ypitch(nq)

         xx(1)=x1
         yy(1)=y1
         do j=2,mm+1
         xx(j)=x1
         yy(j)=y1+(y2-y1)/float(mm)*float(j-1)
         xx(j+mm)=x1+(x2-x1)/float(mm)*float(j-1)
         yy(j+mm)=y2
         xx(j+2*mm)=x2
         yy(j+2*mm)=y2-(y2-y1)/float(mm)*float(j-1)
         xx(j+3*mm)=x2-(x2-x1)/float(mm)*float(j-1)
         yy(j+3*mm)=y1
         enddo

        call rottran(nq,m,xx,yy)
        call tran_pgline(m,xx,yy)


         x1=x1-xbar
         y1=y1-ybar
         x2=x2+xbar
         y2=y2+ybar

         xx(1)=x1
         yy(1)=y1
         do j=2,mm+1
         xx(j)=x1
         yy(j)=y1+(y2-y1)/float(mm)*float(j-1)
         xx(j+mm)=x1+(x2-x1)/float(mm)*float(j-1)
         yy(j+mm)=y2
         xx(j+2*mm)=x2
         yy(j+2*mm)=y2-(y2-y1)/float(mm)*float(j-1)
         xx(j+3*mm)=x2-(x2-x1)/float(mm)*float(j-1)
         yy(j+3*mm)=y1
         enddo

        call rottran(nq,m,xx,yy)
        call tran_pgline(m,xx,yy)

        call pgsci(1)
        call pgslw(1)

! draw plugs

        if (shut_flag.eq.0) then

        do i=1,nplug
        xx(1)=(xplug(i)-1.)*xpitch(nq)
        yy(1)=(yplug(i)-1.)*ypitch(nq)
        xx(2)=xx(1)+xpitch(nq)
        yy(2)=yy(1)
        xx(3)=xx(2)
        yy(3)=yy(1)+ypitch(nq)
        xx(4)=xx(1)
        yy(4)=yy(3)
        xx(5)=xx(1)
        yy(5)=yy(1)
        call rottran(nq,m,xx,yy)
        call tran_pgpoly(5,xx,yy)
        enddo

        endif

        if (shut_flag.eq.1) then

        do i=1,mx
        do j=1,my
        map(i,j)=1.
        enddo
        enddo

        do i=1,ntarg
        mq=nq_target(i)
        if (mq.eq.nq) then
        ii=i_target(i)
        jj=j_target(i)
        map(ii,jj)=0.
        map(ii,jj-1)=0
        map(ii,jj+1)=0
        endif
        enddo

        do i=1,mx
        do j=1,my

        if (map(i,j).eq.1) then
        xx(1)=float(i-1)*xpitch(nq)
        yy(1)=float(j-1)*ypitch(nq)
        xx(2)=xx(1)+xpitch(nq)
        yy(2)=yy(1)
        xx(3)=xx(2)
        yy(3)=yy(1)+ypitch(nq)
        xx(4)=xx(1)
        yy(4)=yy(3)
        xx(5)=xx(1)
        yy(5)=yy(1)
        call rottran(nq,m,xx,yy)
        call tran_pgpoly(5,xx,yy)
        endif

        enddo
        enddo



        endif


        return
        end subroutine draw_msa_quad_on_sky

! --------------------------------------------------------------------

      subroutine tran_pgline(n,x,y)
      implicit none
      integer n,i
      real x(n),y(n),xt(n),yt(n)
      do i=1,n
      call msa_to_sky_rot_bd(x(i),y(i),xt(i),yt(i))
      enddo
      call pgline(n,xt,yt)
      return
      end subroutine tran_pgline

! --------------------------------------------------------------------

      subroutine tran_pgpoly(n,x,y)
      implicit none
      integer n,i
      real x(n),y(n),xt(n),yt(n)
      do i=1,n
      call msa_to_sky_rot_bd(x(i),y(i),xt(i),yt(i))
      enddo
      call pgpoly(n,xt,yt)
      return
      end subroutine tran_pgpoly

! --------------------------------------------------------------------

