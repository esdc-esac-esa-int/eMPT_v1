! v1.0 10/05/22

     include '../empt_modules.f90'

!-------------------------------------
      program msa_fits_read
      use path_to_ref_files
      implicit none
      character*120 filename_chk,filename_sht,filename_vig,filename
      integer check(4,365,171),short(4,365,171),fstop(4,365,171),esa_map(4,365,171)
      integer k,i,j


      ref_path='./' !overriden nominal path because run out of /reference_files

! get instrument model name for reference

      filename=trim(adjustl(ref_path))//'model.conf'
      
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) !model_directory
      call skip_hashes(20)
      read(20,*) filename_chk
      read(20,*) filename_vig
      read(20,*) filename_sht
      close(20)


! flight
      
      call readheader(filename_chk)
      call readtable(filename_chk,check)

      call readheader(filename_sht)
      call readtable(filename_sht,short)

      call readheader(filename_vig)
      call readtable(filename_vig,fstop)

        esa_map(:,:,:)=0


		do k=1,4
		do i=1,365
		do j=1,171

		if (check(k,i,j).eq.0)   esa_map(k,i,j)=2 !(functional)
		if (check(k,i,j).eq.2)   esa_map(k,i,j)=0 !(failed closed - permanent)
		if (check(k,i,j).eq.4)   esa_map(k,i,j)=0 !(failed closed - intermittent)
		if (check(k,i,j).eq.8)   esa_map(k,i,j)=0 !(failed closed - vignetted)

		if ((short(k,i,j).ge.16).and.(short(k,i,j).le.255))   esa_map(k,i,j)=0 !(failed closed - masked off)

		if (check(k,i,j).eq.256) esa_map(k,i,j)=1 !(failed open - permanent)
		if (check(k,i,j).eq.258) esa_map(k,i,j)=1 !(failed open - permanent)
		if (check(k,i,j).eq.512) esa_map(k,i,j)=1 !(failed open - intermittent)

		if (fstop(k,i,j).eq.8)   esa_map(k,i,j)=0 !(failed closed - vignetted)

		enddo !i
		enddo !j
		enddo !k
  
  
        open(12,file='esa_msa_map.dat',status='replace',access='sequential', form='formatted')
        do k=1,4
        write(12,*) ((esa_map(k,i,j),i=1,365),j=1,171)
        enddo       
        close(12)

      end program msa_fits_read

!---------------------------------------------------------------------------
 
include './front_transforms_updt.f90'

include './back_transforms_updt.f90'

include './shutter_routines_new.f90'

include './misc_routines_new.f90'

! *************************************************************************

      subroutine readtable(filename,statcode)

      integer status,unit,readwrite,blocksize,hdutype,ntable,naxes(2)
      integer ncol,col_stat,col_no
      integer statcode(4,365,171)
      character filename*120,nullstr*1,name*8,ttype(500)*20,tform(500)*1
      logical anynull
      integer i,j,k,w

!  The STATUS parameter must always be initialized.
      status=0

!  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

!  Open the FITS file
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)
      
      print *

         ntable=2

!  Move to the next extension
          call ftmahd(unit,ntable,hdutype,status)

          print *,' '
          if (hdutype .eq. 1) then
              print *,'Reading ASCII table in HDU ',ntable
          else if (hdutype .eq. 2) then
              print *,'Reading binary table in HDU ',ntable
          end if


!  Read the TTYPEn keywords, which give the names of the columns
          call ftgkns(unit,'TTYPE',1,500,ttype,ncol,status)


!  Read the TFORMn keywords, which give the types of the columns
          call ftgkns(unit,'TFORM',1,500,tform,ncol,status)
 
          write(*,*) 'Number of columns',ncol
          do i=1,ncol
          write(*,'(i4,3x,a20,3x,a1)') i,ttype(i),tform(i)
          enddo
! 
          print *
          call ftgcno(unit,.false.,'NO',col_no,status)
          print *,'NO number          ',col_no,tform(col_no)
          call ftgcno(unit,.false.,'STATUS',col_stat,status)
          print *,'STATUS column number   ',col_stat,tform(col_stat)
          print *
 
            ntable=1
            do k=1,4
            ntable=ntable+1
            call ftmahd(unit,ntable,hdutype,status)
 
             if (hdutype .eq. 1) then
              print *,'Reading ASCII table in HDU ',ntable
            else if (hdutype .eq. 2) then
              print *,'Reading binary table in HDU ',ntable
            end if
            print *,' '

            do i=1,365
            do j=1,171

            w=i+(j-1)*365

            call ftgcvj(unit,col_stat,w,1,1,0,statcode(k,i,j),anynull,status)

            enddo
            enddo
            enddo

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit,status)
      call ftfiou(unit,status)

!  Check for any error, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      if (status .gt. 0)call printerror(status)
      
      end subroutine readtable 

! *************************************************************************
      subroutine readheader(filename)

!  Print out all the header keywords in all extensions of a FITS file

      integer status,unit,readwrite,blocksize,nkeys,nspace,hdutype,i,j
      character filename*120,record*80

!  The STATUS parameter must always be initialized.
      status=0

!  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

!     name of FITS file 
!      filename='ATESTFILEZ.FITS'

!     open the FITS file, with read-only access.  The returned BLOCKSIZE
!     parameter is obsolete and should be ignored. 
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)

      j = 0
100   continue
      j = j + 1

      print *,'Header listing for HDU', j

!  The FTGHSP subroutine returns the number of existing keywords in the
!  current header data unit (CHDU), not counting the required END keyword,
      call ftghsp(unit,nkeys,nspace,status)

!  Read each 80-character keyword record, and print it out.
      do i = 1, nkeys
          call ftgrec(unit,i,record,status)
          print *,record
      end do

!  Print out an END record, and a blank line to mark the end of the header.
      if (status .eq. 0)then
          print *,'END'
          print *,' '
      end if

!  Try moving to the next extension in the FITS file, if it exists.
!  The FTMRHD subroutine attempts to move to the next HDU, as specified by
!  the second parameter.   This subroutine moves by a relative number of
!  HDUs from the current HDU.  The related FTMAHD routine may be used to
!  move to an absolute HDU number in the FITS file.  If the end-of-file is
!  encountered when trying to move to the specified extension, then a
!  status = 107 is returned.
      call ftmrhd(unit,1,hdutype,status)

      if (status .eq. 0)then
!         success, so jump back and print out keywords in this extension
          go to 100

      else if (status .eq. 107)then
!         hit end of file, so quit
          status=0
      end if

!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)

!  Check for any error, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      if (status .gt. 0)call printerror(status)
      end

! *************************************************************************
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
!       end
! 
! *************************************************************************
