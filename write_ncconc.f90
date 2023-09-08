!***********************************************************************
!* Copyright 2012,2013                                                *
!* Jerome Brioude, Delia Arnold, Andreas Stohl, Wayne Angevine,       *
!* John Burkhart, Massimo Cassiani, Adam Dingwell, Richard C Easter, Sabine Eckhardt,*
!* Stephanie Evan, Jerome D Fast, Don Morton, Ignacio Pisso,          *
!* Petra Seibert, Gerard Wotawa, Caroline Forster, Harald Sodemann,   *
!* Adam Dingwell                                                      *
!*                                                                    *
!* This file is part of FLEXPART WRF                                  *
!                                                                     *
! FLEXPART is free software: you can redistribute it and/or modify    *
! it under the terms of the GNU General Public License as published by*
! the Free Software Foundation, either version 3 of the License, or   *
! (at your option) any later version.                                 *
!                                                                     *
! FLEXPART is distributed in the hope that it will be useful,         *
! but WITHOUT ANY WARRANTY; without even the implied warranty of      *
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
! GNU General Public License for more details.                        *
!                                                                     *
! You should have received a copy of the GNU General Public License   *
! along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.   *
!**********************************************************************

subroutine write_ncconc(itime,outnum,ks,kp,nage,tot_mu_scalar,nesting_level)
  
  !*****************************************************************************
  !                                                                            *
  !  This routine writes concentration, mixing ratio and deposition fields     *
  !  to a netcdf file defined by flex_ncheader.                                *
  !                                                                            *
  !  flex_ncheader is called from within write_ncconc when it's time for a new *
  !  output file.                                                              *
  !                                                                            *
  !  write_ncconc should be called by concoutput_irreg and concoutput_reg      *
  !  it is separate from the binary and ascii output routines to avoid mixing  *
  !  of sparse and full grid approaches.  Netcdf will output the full grid.    *
  !                                                                            *
  !      Author: A. Dingwell                                                   *
  !                                                                            *
  !      29 May 2013                                                           *
  !                                                                            *
  ! Modifications:                                                             *
  ! June 5 2013: J. Brioude: compression using deflate level, optimization of  *
  !  the writing procedure. bug fixes for backtrajectory mode                  *
  ! Feb 2014: A. Dingwell: bug fix
  !*****************************************************************************
  
  use point_mod
  use outg_mod
  use par_mod
  use com_mod

  implicit none

  include 'netcdf.inc'

  real    :: outnum         ! Number of samples for each concentration calculation
  integer :: itime          ! Current simulation time [s]
  integer :: ks,kp,nage     ! species, maxpointspec_act and nageclass indices resp.
  real    :: tot_mu_scalar  ! total mass per source and species (backward)
                            ! or unity (forward).  Should probably be sent as
                            ! tot_mu(ks,kp) from concoutput*.f90
  integer :: nesting_level  ! 0 for main (mother) grid, 1 for nest (child)

  real(kind=dp) :: jul          ! Julian date
  integer   :: jjjjmmdd,ihmmss  ! date & time as integer
  character :: adate*8,atime*6  ! date and time strings, used for filename

  integer :: ncid           ! Pointer to netcdf file, depends on nesting level
  integer :: grid_nx,grid_ny! outgrid dimensions, depend on the nesting level
  integer :: ncret          ! Netcdf:  return code
  integer :: ix,jy,kz       ! iterators
  character :: datestr*15   ! For the Times variable
  integer :: deflate_level=5 ! compression level
    
  if (option_verbose.ge.1) then
    write(*,*)'write_ncconc: writing netcdf output for: domain,kp,nage =',&
      nesting_level+1,kp,nage
  endif

  ! Determine which nest/outfile we are writing to
  !***********************************************
  if (nesting_level.eq.0) then
    ncid    = ncout
    grid_nx = numxgrid
    grid_ny = numygrid
  elseif (nesting_level.eq.1) then
    ncid    = ncoutn
    grid_nx = numxgridn
    grid_ny = numygridn
  else
    write(*,*) '***write_ncconc error: nesting level must be 0 or 1'
    ! Note for future development: If additional output nests are to be
    ! supported for netcdf output, modification must be made here as well as in
    ! the respective nesting_level if-block in write_ncheader
  endif
  ! Update/Initialize record index
  !*******************************
   if ((ks.eq.1).and.(kp.eq.1).and.(nage.eq.1)) then
!   print*,'ncirec',ncirec,ncnumrec
  if (nesting_level.eq.0) then  ! Only update for first domain
    if (itime.eq.loutstep) then  ! first output
      ncirec = 1  ! initialize record index
    elseif (ncirec.eq.ncnumrec) then  ! new file
!      print*,'file is closing'
      ncirec = 1  ! reset record index
      ncret=nf_close(ncid)      ! close the old file
      call check_ncerror(ncret)
!      print*,'file is closed'
    else
      ncirec=ncirec+1 ! move on to next record
    endif
  endif
!   print*,'ncirec',ncirec,ncnumrec
  endif

  ! Check if it's time to create a new netcdf file
  !***********************************************
  if (ncirec.eq.1) then         ! First output in current file?
!   write(*,*) 'itime=',itime
   if ((ks.eq.1).and.(kp.eq.1).and.(nage.eq.1)) then
!    if (itime.ne.loutstep) then ! Not the first output file?
!      ncret=nf_close(ncid)      ! close the old file
!      call check_ncerror(ncret)
!     print*,'file is closed'
!    endif
!   call write_ncheader(itime,nesting_level)  ! Create new file
    if (option_verbose.ge.1) &
      write(*,*)'write_ncconc: calling write_ncinfo'
    call write_ncinfo(itime,nesting_level)  ! Create new file
    ! Reassign file handle to the newly created file:
     endif
    if (nesting_level.eq.0) ncid=ncout
    if (nesting_level.eq.1) ncid=ncoutn
  endif

  if (option_verbose.ge.50) &
    write(*,*) 'ncid,nccovid=',ncid,nccovid


  ! Create output for the current record index
  !*******************************************
  jul=bdate+real(itime,kind=dp)/86400._dp
  call caldate(jul,jjjjmmdd,ihmmss)
  write(adate,'(i8.8)') jjjjmmdd
  write(atime,'(i6.6)') ihmmss
  
   if ((ks.eq.1).and.(kp.eq.1).and.(nage.eq.1)) then

  if (option_verbose.ge.10) write(*,*)'write_ncconc: record index',ncirec
  write(datestr,'(I8.8,A1,I6.6)') jjjjmmdd,'_',ihmmss
  ncret = nf_put_vara_text(ncid,ncrecvid,(/1,ncirec/),(/15,1/),datestr)
  call check_ncerror(ncret)
   endif

   if ((iout.eq.1).or.(iout.eq.3).or.(iout.eq.5)) then ! concentration
    if (option_verbose.ge.1)  &
      write(*,*)'write_ncconc: concentration output',kp,nage,ncirec,nccovid,ncid

    write(*,*) "[JK] Test2: ", area(ix,jy)*hmix(ix,iy,1,2), volume(ix,iy,kz)

    do kz=1,numzgrid
    do jy=0,grid_ny-1
    do ix=0,grid_nx-1
      grid2(ix,jy,kz,kp,nage)= grid(ix,jy,kz)*factor3d(ix,jy,kz)/tot_mu_scalar
    enddo ! ix=1,grid_nx-1
    enddo ! jy=1,grid_ny-1
    enddo ! kz=1,numzgrid

    if (kp.eq.maxpointspec_act .and. nage.eq.nageclass) then
     if (ldirect.eq.-1) then 
      ncret = nf_put_vara_real(ncid,nccovid, &
        (/1,1,1,1,1,ncirec/),(/grid_nx,grid_ny,numzgrid,kp,nage,1/), &
        grid2(0:grid_nx-1,0:grid_ny-1,1:numzgrid,1:kp,1:nage))
      call check_ncerror(ncret)
     else
      if (kp.gt.1) then
       ncret = nf_put_vara_real(ncid,nccovid, &
        (/1,1,1,1,ks,1,ncirec/),(/grid_nx,grid_ny,numzgrid,kp,1,nage,1/), &
        grid2(0:grid_nx-1,0:grid_ny-1,1:numzgrid,1:kp,1:nage)) 
      else
       ncret = nf_put_vara_real(ncid,nccovid, &
        (/1,1,1,ks,1,ncirec/),(/grid_nx,grid_ny,numzgrid,1,nage,1/), &
        grid2(0:grid_nx-1,0:grid_ny-1,1:numzgrid,1,1:nage)) 
      endif
      call check_ncerror(ncret)
     endif
    endif
!    do kz=1,numzgrid
!    do jy=0,grid_ny-1
!    do ix=0,grid_nx-1
!      ncret = nf_put_vara_real(ncid,nccovid, &
!        (/ix+1,jy+1,kz,kp,nage,ncirec/),(/1,1,1,1,1,1/), &
!        grid(ix,jy,kz)*factor3d(ix,jy,kz)/tot_mu_scalar)
!      call check_ncerror(ncret)
!    enddo ! ix=1,grid_nx-1
!    enddo ! jy=1,grid_ny-1
!    enddo ! kz=1,numzgrid
   endif ! concentraion

   if ((iout.eq.2).or.(iout.eq.3)) then  ! mixing ratio
    if (option_verbose.ge.1)write(*,*)'write_ncconc: mixing ratio output'

    write(*,*) "[JK] Test3: ", area(ix,jy)*hmix(ix,iy,1,2), volume(ix,iy,kz)
    do kz=1,numzgrid
    do jy=0,grid_ny-1
    do ix=0,grid_nx-1
!      grid3(ix,jy,kz,kp,nage)= 1.e12*grid(ix,jy,kz)/volume(ix,jy,kz)/outnum*  &
!        weightair/weightmolar(ks)/densityoutgrid(ix,jy,kz)
!   // JK: trying to get mean mixing ratio in boundary layer only
      write(*,*) "[JK] Test: ", area(ix,jy)*hmix(ix,iy,1,2), volume(ix,iy,kz)

      grid3(ix,jy,kz,kp,nage)= 1.e12*grid(ix,jy,kz)/(area(ix,jy)*hmix(ix,iy,1,2))/outnum*  &
        weightair/weightmolar(ks)/densityoutgrid(ix,jy,kz)
    enddo ! ix=1,grid_nx-1
    enddo ! jy=1,grid_ny-1
    enddo ! kz=1,numzgrid
    if (kp.eq.maxpointspec_act .and. nage.eq.nageclass) then
     if (ldirect.eq.-1) then
      ncret = nf_put_vara_real(ncid,ncravid, &
        (/1,1,1,1,1,ncirec/),(/grid_nx,grid_ny,numzgrid,kp,nage,1/), &
        grid3(0:grid_nx-1,0:grid_ny-1,1:numzgrid,1:kp,1:nage))
      call check_ncerror(ncret)
     else
      if (kp.gt.1) then
       ncret = nf_put_vara_real(ncid,ncravid, &
        (/1,1,1,1,ks,1,ncirec/),(/grid_nx,grid_ny,numzgrid,kp,1,nage,1/), &
        grid3(0:grid_nx-1,0:grid_ny-1,1:numzgrid,1:kp,1:nage))
      else
       ncret = nf_put_vara_real(ncid,ncravid, &
        (/1,1,1,ks,1,ncirec/),(/grid_nx,grid_ny,numzgrid,1,nage,1/), &
        grid3(0:grid_nx-1,0:grid_ny-1,1:numzgrid,1,1:nage))
      endif
      call check_ncerror(ncret)
     endif
    endif

!    do kz=1,numzgrid
!    do jy=0,grid_ny-1
!    do ix=0,grid_nx-1
!      ncret = nf_put_vara_real(ncid,ncravid, &
!        (/ix+1,jy+1,kz,kp,nage,ncirec/),(/1,1,1,1,1,1/), &
!        1.e12*grid(ix,jy,kz)/volume(ix,jy,kz)/outnum*  &
!        weightair/weightmolar(ks)/densityoutgrid(ix,jy,kz))
!      call check_ncerror(ncret)
!    enddo ! ix=1,grid_nx-1
!    enddo ! jy=1,numygrid-1
!    enddo ! kz=1,numzgrid
     endif ! mixing ratio

  if ((ldirect.eq.1).and.(WETDEP)) then ! WETDEP
    if (option_verbose.ge.1)write(*,*)'write_ncconc: wet deposition output'
    do jy=0,grid_ny-1
    do ix=0,grid_nx-1
    if (nesting_level.eq.0)  wetgrid2(ix,jy,kp,nage)=1.e12*wetgrid(ix,jy)/area(ix,jy)
    if (nesting_level.eq.1)  wetgrid2(ix,jy,kp,nage)=1.e12*wetgrid(ix,jy)/arean(ix,jy)
    enddo ! ix=1,grid_nx-1
    enddo ! jy=1,grid_ny-1
    if (kp.eq.maxpointspec_act .and. nage.eq.nageclass) then
  if (ldirect.eq.-1) then
      ncret = nf_put_vara_real(ncid,ncwdvid, &
        (/1,1,1,1,ncirec/),(/grid_nx,grid_ny,kp,nage,1/), &
        wetgrid2(0:grid_nx-1,0:grid_ny-1,1:kp,1:nage))
      call check_ncerror(ncret)
  else
    if (kp.gt.1) then
      ncret = nf_put_vara_real(ncid,ncwdvid, &
        (/1,1,1,ks,1,ncirec/),(/grid_nx,grid_ny,kp,1,nage,1/), &
        wetgrid2(0:grid_nx-1,0:grid_ny-1,1:kp,1:nage))
    else
      ncret = nf_put_vara_real(ncid,ncwdvid, &
        (/1,1,ks,1,ncirec/),(/grid_nx,grid_ny,1,nage,1/), &
        wetgrid2(0:grid_nx-1,0:grid_ny-1,1,1:nage))
    endif
      call check_ncerror(ncret)
    endif
  endif
!    do jy=0,grid_ny-1
!    do ix=0,grid_nx-1
!      ncret = nf_put_vara_real(ncid,ncwdvid, &
!        (/ix+1,jy+1,kp,nage,ncirec/),(/1,1,1,1,1/), &
!        1.e12*wetgrid(ix,jy)/area(ix,jy))
!      call check_ncerror(ncret)
!    enddo ! ix=1,grid_nx-1
!    enddo ! jy=1,numygrid-1
  endif ! WETDEP

  if ((ldirect.eq.1).and.(DRYDEP)) then ! DRYDEP
    if (option_verbose.ge.1)write(*,*)'write_ncconc: dry deposition output'
    do jy=0,grid_ny-1
    do ix=0,grid_nx-1
    if (nesting_level.eq.0)  drygrid2(ix,jy,kp,nage)=1.e12*drygrid(ix,jy)/area(ix,jy)
    if (nesting_level.eq.1)  drygrid2(ix,jy,kp,nage)=1.e12*drygrid(ix,jy)/arean(ix,jy)
    enddo ! ix=1,grid_nx-1
    enddo ! jy=1,grid_ny-1
    if (kp.eq.maxpointspec_act .and. nage.eq.nageclass) then
  if (ldirect.eq.-1) then
      ncret = nf_put_vara_real(ncid,ncddvid, &
        (/1,1,1,1,ncirec/),(/grid_nx,grid_ny,kp,nage,1/), &
        drygrid2(0:grid_nx-1,0:grid_ny-1,1:kp,1:nage))
      call check_ncerror(ncret)
  else
    if (kp.gt.1) then
      ncret = nf_put_vara_real(ncid,ncddvid, &
        (/1,1,1,ks,1,ncirec/),(/grid_nx,grid_ny,kp,1,nage,1/), &
        drygrid2(0:grid_nx-1,0:grid_ny-1,1:kp,1:nage))
    else
      ncret = nf_put_vara_real(ncid,ncddvid, &
        (/1,1,ks,1,ncirec/),(/grid_nx,grid_ny,1,nage,1/), &
        drygrid2(0:grid_nx-1,0:grid_ny-1,1,1:nage))
    endif
      call check_ncerror(ncret)
    endif
  endif

!    do jy=0,grid_ny-1
!    do ix=0,grid_nx-1
!      ncret = nf_put_vara_real(ncid,ncddvid, &
!        (/ix+1,jy+1,kp,nage,ncirec/),(/1,1,1,1,1/), &
!        1.e12*drygrid(ix,jy)/area(ix,jy))
!      call check_ncerror(ncret)
!    enddo ! ix=1,grid_nx-1
!    enddo ! jy=1,numygrid-1
  endif ! DRYDEP

   ncret=nf_sync(ncid)
   call check_ncerror(ncret)

end subroutine write_ncconc
