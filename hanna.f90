!**********************************************************************
! Copyright 1998,1999,2000,2001,2002,2005,2007,2008,2009,2010         *
! Andreas Stohl, Petra Seibert, A. Frank, Gerhard Wotawa,             *
! Caroline Forster, Sabine Eckhardt, John Burkhart, Harald Sodemann   *
!                                                                     *
! This file is part of FLEXPART.                                      *
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

subroutine hanna(z, &
   ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw, &
   sigw,dsigwdz,dsigw2dz)

  !                 i
  !*****************************************************************************
  !                                                                            *
  !   Computation of \sigma_i and \tau_L based on the scheme of Hanna (1982)   *
  !                                                                            *
  !   Author: A. Stohl                                                         *
  !                                                                            *
  !   4 December 1997                                                          *
  !                                                                            *
  !*****************************************************************************
  !                                                                            *
  ! Variables:                                                                 *
  ! dsigwdz [1/s]     vertical gradient of sigw                                *
  ! ol [m]            Obukhov length                                           *
  ! sigu, sigv, sigw  standard deviations of turbulent velocity fluctuations   *
  ! tlu [s]           Lagrangian time scale for the along wind component.      *
  ! tlv [s]           Lagrangian time scale for the cross wind component.      *
  ! tlw [s]           Lagrangian time scale for the vertical wind component.   *
  ! ust, ustar [m/s]  friction velocity                                        *
  ! wst, wstar [m/s]  convective velocity scale                                *
  !                                                                            *
  !*****************************************************************************

  use par_mod
  use com_mod
!  use hanna_mod

  implicit none

  real :: corr,z
  real :: ust,wst,ol,h,zeta,sigu,sigv,tlu,tlv,tlw
  real :: sigw,dsigwdz,dsigw2dz



  !**********************
  ! 1. Neutral conditions
  !**********************

  if (h/abs(ol).lt.1.) then
    ust=max(1.e-4,ust)
    corr=z/ust
    sigu=1.e-2+2.0*ust*exp(-3.e-4*corr)
    sigw=1.3*ust*exp(-2.e-4*corr)
    dsigwdz=-2.e-4*sigw
    sigw=sigw+1.e-2
    sigv=sigw
    tlu=0.5*z/sigw/(1.+1.5e-3*corr)
    tlv=tlu
    tlw=tlu
    !write (*,*) 'PBLH:', h


  !***********************
  ! 2. Unstable conditions
  !***********************

  else if (ol.lt.0.) then
    !write (*,*) 'PBLH:', h


  ! Determine sigmas
  !*****************

    sigu=1.e-2+ust*(12-0.5*h/ol)**0.33333
    sigv=sigu
    sigw=sqrt(1.2*wst**2*(1.-.9*zeta)*zeta**0.66666+ &
         (1.8-1.4*zeta)*ust**2)+1.e-2
    dsigwdz=0.5/sigw/h*(-1.4*ust**2+wst**2* &
         (0.8*max(zeta,1.e-3)**(-.33333)-1.8*zeta**0.66666))


  ! Determine average Lagrangian time scale
  !****************************************

    tlu=0.15*h/sigu
    tlv=tlu
    if (z.lt.abs(ol)) then
      tlw=0.1*z/(sigw*(0.55-0.38*abs(z/ol)))
    else if (zeta.lt.0.1) then
      tlw=0.59*z/sigw
    else
      tlw=0.15*h/sigw*(1.-exp(-5*zeta))
    endif


  !*********************
  ! 3. Stable conditions
  !*********************

  else
    sigu=1.e-2+2.*ust*(1.-zeta)
    sigv=1.e-2+1.3*ust*(1.-zeta)
    sigw=sigv
    dsigwdz=-1.3*ust/h
    tlu=0.15*h/sigu*(sqrt(zeta))
    tlv=0.467*tlu
    tlw=0.1*h/sigw*zeta**0.8
  endif


  tlu=max(10.,tlu)
  tlv=max(10.,tlv)
  tlw=max(30.,tlw)

  if (dsigwdz.eq.0.) dsigwdz=1.e-10

end subroutine hanna
