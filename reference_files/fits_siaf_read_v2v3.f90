      program fits_table_read
      implicit none
      character*80 siafname,filename
      real V2Ref,V3Ref,AngV3


! get instrument model siaf name for reference

      filename='model.conf'
      
      open(20,file=filename,status='old',access='sequential',form='formatted')
      call skip_hashes(20)
      read(20,*) !model_directory
      call skip_hashes(20)
      read(20,*) 
      read(20,*) !MSA tables
      read(20,*) 
      call skip_hashes(20)
      read(20,*) siafname
      close(20)

      call readtable(siafname,V2Ref,V3Ref,AngV3)
      
      print *,'V2Ref=',V2Ref
      print *,'V3Ref=',V3Ref
      print *,'V3Ang=',AngV3
      print *

      open(20,file='msa_ref_v2v3.ascii',status='replace',access='sequential',form='formatted')
      write(20,'(a)') '# NRS_MSA_FULL V2,V3 reference point'
      write(20,*)
      write(20,'(a)') '*V2Ref'
      write(20,'(f12.5)') V2Ref
      write(20,*)
      write(20,'(a)') '*V3Ref'
      write(20,'(f12.5)') V3Ref     
      write(20,*)
      write(20,'(a)') '*AngV3'
      write(20,'(f12.5)') AngV3    
      close(20)

      print *,'Reference values written to msa_ref_v2v3.ascii file'

      end

! *************************************************************************

      subroutine readtable(filename,V2Ref,V3Ref,AngV3)
      implicit none
      integer status,unit,readwrite,blocksize,hdutype,ntable,naxes(2),ncol,nfound
      character filename*80,nullstr*1,name*8,ttype(500)*20,tform(500)*3
      logical anynull
      character*20 ap_name(60)
      real c1xm(60),c1ym(60),c1v2(60),c1v3(60)
      real c2xm(60),c2ym(60),c2v2(60),c2v3(60)
      real c3xm(60),c3ym(60),c3v2(60),c3v3(60)
      real c4xm(60),c4ym(60),c4v2(60),c4v3(60)
      real rxm(60),rym(60),rv2(60),rv3(60)
      real ang(60)
      integer i
      real V2Ref,V3Ref,AngV3

!  The STATUS parameter must always be initialized.
      status=0


!  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

!  Open the FITS file
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)
      
!       call readheader(filename)

         ntable=2

!  Move to the next extension
          call ftmahd(unit,ntable,hdutype,status)

          print *,' '
          if (hdutype .eq. 1) then
              print *,'Reading ASCII table in HDU ',ntable
          else if (hdutype .eq. 2) then
              print *,'Reading binary table in HDU ',ntable
          end if

!  Determine the size of the image.
      call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
      write(*,*) 'Number of rows   ',naxes(2)


!  Read the TTYPEn keywords, which give the names of the columns
          call ftgkns(unit,'TTYPE',1,500,ttype,ncol,status)


!  Read the TFORMn keywords, which give the types of the columns
          call ftgkns(unit,'TFORM',1,500,tform,ncol,status)
 
          write(*,*) 'Number of columns',ncol
!           do i=1,ncol
!           write(*,'(i4,3x,a20,3x,a3)') i,ttype(i),tform(i)
!           enddo

 
!  Read the data, one row at a time, and print them out

           do i=1,naxes(2)

  
               call ftgcvs(unit,1,i,1,1,' ',ap_name(i),anynull,status)

              call ftgcve(unit,2,i,1,1,0,c1xm(i),anynull,status)
              call ftgcve(unit,3,i,1,1,0,c1ym(i),anynull,status)
              call ftgcve(unit,4,i,1,1,0,c2xm(i),anynull,status)
              call ftgcve(unit,5,i,1,1,0,c2ym(i),anynull,status)
              call ftgcve(unit,6,i,1,1,0,c3xm(i),anynull,status)
              call ftgcve(unit,7,i,1,1,0,c3ym(i),anynull,status)
              call ftgcve(unit,8,i,1,1,0,c4xm(i),anynull,status)
              call ftgcve(unit,9,i,1,1,0,c4ym(i),anynull,status)
              call ftgcve(unit,10,i,1,1,0,rxm(i),anynull,status)
              call ftgcve(unit,11,i,1,1,0,rym(i),anynull,status)
              call ftgcve(unit,12,i,1,1,0,c1v2(i),anynull,status)
              call ftgcve(unit,13,i,1,1,0,c1v3(i),anynull,status)
              call ftgcve(unit,14,i,1,1,0,c2v2(i),anynull,status)
              call ftgcve(unit,15,i,1,1,0,c2v3(i),anynull,status)
              call ftgcve(unit,16,i,1,1,0,c3v2(i),anynull,status)
              call ftgcve(unit,17,i,1,1,0,c3v3(i),anynull,status)
              call ftgcve(unit,18,i,1,1,0,c4v2(i),anynull,status)
              call ftgcve(unit,19,i,1,1,0,c4v3(i),anynull,status)
              call ftgcve(unit,20,i,1,1,0,rv2(i),anynull,status)
              call ftgcve(unit,21,i,1,1,0,rv3(i),anynull,status)
              call ftgcve(unit,22,i,1,1,0,ang(i),anynull,status)
           
           enddo

 

      write(*,*) 
!       do i=1,naxes(2)
!          write(*,*) i,ap_name(i),c1xm(i),c1ym(i),c2xm(i),c2ym(i),c3xm(i),c3ym(i),c4xm(i),c4ym(i),rxm(i),rym(i), &
!          c1v2(i),c1v3(i),c2v2(i),c2v3(i),c3v2(i),c3v3(i),c4v2(i),c4v3(i),rv2(i),rv3(i),ang(i)
!          write(10,*) ap_name(i),c1xm(i),c1ym(i),c2xm(i),c2ym(i),c3xm(i),c3ym(i),c4xm(i),c4ym(i),rxm(i),rym(i), &
!          c1v2(i),c1v3(i),c2v2(i),c2v3(i),c3v2(i),c3v3(i),c4v2(i),c4v3(i),rv2(i),rv3(i),ang(i)
!  
!         write(*,*) i,ap_name(i),c1v2(i),c1v3(i),c2v2(i),c2v3(i),c3v2(i),c3v3(i),c4v2(i),c4v3(i),rv2(i),rv3(i),ang(i)
!  
!        enddo

!       write(*,*) 
!       do i=38,41
!       write(*,*) ap_name(i),c1v2(i),c1v3(i),c2v2(i),c2v3(i),c3v2(i),c3v3(i),c4v2(i),c4v3(i) 
!       enddo
!       write(*,*) rv2(37),rv3(37),rxm(37)*1000.,rym(37)*1000.
!       write(*,*) 
!       do i=43,46
!       write(*,*) ap_name(i),c1v2(i),c1v3(i),c2v2(i),c2v3(i),c3v2(i),c3v3(i),c4v2(i),c4v3(i) 
!       enddo
!       write(*,*) rv2(42),rv3(42),rxm(42)*1000.,rym(42)*1000.


        V2Ref=rv2(37)
        V3ref=rv3(37)
        AngV3=ang(37)

 
!  The FITS file must always be closed before exiting the program. 
!  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit,status)
      call ftfiou(unit,status)

!  Check for any error, and if so print out error messages.
!  The PRINTERROR subroutine is listed near the end of this file.
      if (status .gt. 0)call printerror(status)

      return
      end

! *************************************************************************

      subroutine printerror(status)

!  This subroutine prints out the descriptive text corresponding to the
!  error status value and prints out the contents of the internal
!  error message stack generated by FITSIO whenever an error occurs.

      integer status
      character errtext*30,errmessage*80

!  Check if status is OK (no error); if so, simply return
      if (status .le. 0)return

!  The FTGERR subroutine returns a descriptive 30-character text string that
!  corresponds to the integer error status number.  A complete list of all
!  the error numbers can be found in the back of the FITSIO User's Guide.
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

!  FITSIO usually generates an internal stack of error messages whenever
!  an error occurs.  These messages provide much more information on the
!  cause of the problem than can be provided by the single integer error
!  status value.  The FTGMSG subroutine retrieves the oldest message from
!  the stack and shifts any remaining messages on the stack down one
!  position.  FTGMSG is called repeatedly until a blank message is
!  returned, which indicates that the stack is empty.  Each error message
!  may be up to 80 characters in length.  Another subroutine, called
!  FTCMSG, is available to simply clear the whole error message stack in
!  cases where one is not interested in the contents.
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
      end do
      end

! *************************************************************************

      subroutine readheader(filename)

!  Print out all the header keywords in all extensions of a FITS file

      integer status,unit,readwrite,blocksize,nkeys,nspace,hdutype,i,j
      character filename*80,record*80

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

!-------------------------------------------------------------------------------

      subroutine skip_hashes(lu)
      implicit none
      integer lu
      character*1 hash
10    read(lu,11) hash
11    format((a1))
      if (hash.ne.'#') then
      backspace(lu)
      goto 12
      endif
      goto 10
12    return
      end subroutine skip_hashes

!---------------------------------------------------------------------------
