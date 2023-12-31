!***********************************************************************
!* Copyright 2012,2013                                                *
!* Jerome Brioude, Delia Arnold, Andreas Stohl, Wayne Angevine,       *
!* John Burkhart, Massimo Cassiani, Adam Dingwell, Richard C Easter, Sabine Eckhardt,*
!* Stephanie Evan, Jerome D Fast, Don Morton, Ignacio Pisso,          *
!* Petra Seibert, Gerard Wotawa, Caroline Forster, Harald Sodemann,   *
!* Adam Dingwell,                                                      *
!*                                                                     *
!* This file is part of FLEXPART WRF                                   *
!*                                                                     *
!* FLEXPART is free software: you can redistribute it and/or modify    *
!* it under the terms of the GNU General Public License as published by*
!* the Free Software Foundation, either version 3 of the License, or   *
!* (at your option) any later version.                                 *
!*                                                                     *
!* FLEXPART is distributed in the hope that it will be useful,         *
!* but WITHOUT ANY WARRANTY; without even the implied warranty of      *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
!* GNU General Public License for more details.                        *
!*                                                                     *
!* You should have received a copy of the GNU General Public License   *
!* along with FLEXPART.  If not, see <http://www.gnu.org/licenses/>.   *
!***********************************************************************
!*******************************************************************************
!        Include file for particle diffusion model FLEXPART                    *
!        This file contains a global common block used by FLEXPART             *
!                                                                              *
!        Author: A. Stohl                                                      *
!                                                                              *
!        June 1996                                                             *
!                                                                              *
!        Last update: 9 August 2000                                            *
!                                                                              *
!*******************************************************************************

module com_mod

  use par_mod, only: dp, numpath, maxnests, maxageclass, maxspec, ni, &
       numclass, nymax, nxmax, maxcolumn, maxwf, nzmax, nxmaxn, nymaxn, &
       maxreceptor, maxpart, maxrand, nwzmax, nuvzmax,maxomp,maxrandomp, &
       maxrand3,max_stream
  implicit none

  !****************************************************************
  ! Variables defining where FLEXPART input/output files are stored
  !****************************************************************

  character :: path(numpath+2*maxnests)*120
  integer :: length(numpath+2*maxnests)
  character :: inputname*200
  integer   :: len2
  integer   :: ncout,ncoutn

  ! path                    path names needed for trajectory model
  ! length                  length of path names needed for trajectory model
  ! ncout                   file handle of netcdf output file
  ! ncoutn                  file handle of netcdf output nest file


  !********************************************************
  ! Variables defining the general model run specifications
  !********************************************************

  integer :: ibdate,ibtime,iedate,ietime
  real(kind=dp) :: bdate,edate


  ! ibdate                  beginning date (YYYYMMDD)
  ! ibtime                  beginning time (HHMISS)
  ! iedate                  ending date (YYYYMMDD)
  ! ietime                  ending time (HHMISS)
  ! bdate                   beginning date of simulation (julian date)
  ! edate                   ending date of simulation (julian date)


  integer :: ldirect,ideltas

  ! ldirect                 1 for forward, -1 for backward simulation
  ! ideltas                 length of trajectory loop from beginning to
  !                    ending date (s)

  integer :: loutstep,loutaver,loutsample,method,lsynctime
  real :: outstep
  integer :: ncirec,ncnumrec

  ! loutstep [s]            gridded concentration output every loutstep seconds
  ! loutaver [s]            concentration output is an average over [s] seconds
  ! loutsample [s]          sampling interval of gridded concentration output
  ! lsynctime [s]           synchronisation time of all particles
  ! method                  indicator which dispersion method is to be used
  ! outstep = real(abs(loutstep))
  ! ncired                  record index for netcdf files
  ! ncnumrec                maximum record length of netcdf output files

  real :: ctl,fine,dt_conv
  integer :: ifine,iout,ipout,ipin,iflux,mdomainfill
  integer :: mquasilag,nested_output,ind_source,ind_receptor
  integer :: ind_rel,ind_samp,ioutputforeachrelease,linit_cond
  integer :: turb_option, method_w_terrain_correction
  integer :: sfc_option
  integer :: lu_option
  integer :: wind_option,time_option
  integer :: outgrid_option
  integer :: numpoint_option,option_verbose
  integer :: iouttype
  integer :: cblflag ! added by mc for CBL skewed pdf calculations

  logical :: turbswitch

  ! ctl           factor, by which time step must be smaller than Lagrangian time scale
  ! fine          real(ifine)
  ! ifine         reduction factor for time step used for vertical wind
  !     Langevin equation for the vertical wind component
  ! iout          output options: 1 conc. output (ng/m3), 2 mixing ratio (pptv), 3 both
  ! ipout         particle dump options: 0 no, 1 every output interval, 2 only at end
  ! ipin          read in particle positions from dumped file from a previous run
  ! iflux         flux calculation options: 1 calculation of fluxes, 2 no fluxes
  ! mdomainfill   0: normal run
  !               1: particles are initialized according to atmospheric mass distribution
  ! mquasilag     0: normal run
  !               1: Particle position output is produced in a condensed format and particles are numbered
  ! nested_output 0: no
  !               1: yes
  ! ind_source    switches between different units for concentrations at the source
  !               NOTE that in backward simulations the release of computational particles
  !               takes place at the "receptor" and the sampling of particles at the "source".
  !               1: mass units
  !               2: mass mixing ratio units
  ! ind_receptor  switches between different units for FLEXPART concentration at the receptor
  !               1: mass units
  !               2: mass mixing ratio units
  ! ind_rel, ind_samp  are used within the code to change between mass and mass-mix (see readcommand.f)
  ! ioutputforeachrelease   Should each release be a seperate output field?
  ! linit_cond    switch on the output of sensitivity to initial conditions for backward runs
  !               0: no
  !               1: mass unit
  !               2: mass mixing ratio unit
  ! turboption    
  ! method_w_terrain_correction
  ! sfc_option
  ! lu_option     0: Old land use (MODIS IGBP distributed with FLEXPART)
  !               1: Get land use from WRF output
  ! turbswitch    determines how the Markov chain is formulated


  integer :: mintime,itsplit

  ! mintime                 minimum time step to be used by FLEXPART
  ! itsplit                 time constant for splitting particles

  integer :: lsubgrid,lconvection,lagespectra

  ! lsubgrid     1 if subgrid topography parameterization switched on, 2 if not
  ! lconvection  1 if convection parameterization switched on, 0 if not
  ! lagespectra  1 if age spectra calculation switched on, 2 if not


  integer :: nageclass,lage(maxageclass)

  ! nageclass               number of ageclasses for the age spectra calculation
  ! lage [s]                ageclasses for the age spectra calculation


  logical :: gdomainfill

  ! gdomainfill             .T., if domain-filling is global, .F. if not



  !*********************************************************************
  ! Variables defining the release locations, released species and their
  ! properties, etc.
  !*********************************************************************

  !change Sabine Eckhardt, only save the first 1000 identifier for releasepoints
  character :: compoint(2001)*45
  integer :: numpoint
  !sec, now dynamically allocated:
  ! ireleasestart(maxpoint),ireleaseend(maxpoint)
  !      real xpoint1(maxpoint),ypoint1(maxpoint)
  !real xpoint2(maxpoint),ypoint2(maxpoint)
  !real zpoint1(maxpoint),zpoint2(maxpoint)
  !integer*2 kindz(maxpoint)
  integer :: specnum(maxspec)
  real :: decay(maxspec)
  real :: weta(maxspec),wetb(maxspec)
  real :: reldiff(maxspec),henry(maxspec),f0(maxspec)
  real :: density(maxspec),dquer(maxspec),dsigma(maxspec)
  real :: vsetaver(maxspec),cunningham(maxspec),weightmolar(maxspec)
  real :: vset(maxspec,ni),schmi(maxspec,ni),fract(maxspec,ni)
  real :: ri(5,numclass),rac(5,numclass),rcl(maxspec,5,numclass)
  real :: rgs(maxspec,5,numclass),rlu(maxspec,5,numclass)
  real :: rm(maxspec),dryvel(maxspec),kao(maxspec),ohreact(maxspec)
  ! se  it is possible to associate a species with a second one to make transfer from gas to aerosol
  integer :: spec_ass(maxspec)

  real :: area_hour(maxspec,24),point_hour(maxspec,24)
  real :: area_dow(maxspec,7),point_dow(maxspec,7)

  !integer npart(maxpoint)
  integer :: nspec,maxpointspec_act,link(maxspec)
  character(len=10) :: species(maxspec)


  ! compoint                comment, also "name" of each starting point
  ! numpoint                actual number of trajectory starting/ending points
  ! ireleasestart,ireleaseend [s] starting and ending time of each release
  ! xmass(point,spec)       total mass emitted
  ! xpoint1,ypoint1         lower left coordinates of release area
  ! xpoint2,ypoint2         upper right coordinates of release area
  ! zpoint1,zpoint2         min./max. z-coordinates of release points
  ! kindz                   1: zpoint is in m agl, 2: zpoint is in m asl
  ! npart                   number of particles per release point
  ! nspec                   number of different species allowed for one release
  ! maxpointspec_act        number of releaspoints for which a different output shall be created
  ! species                 name of species
  ! decay                   decay constant of radionuclide

  ! WET DEPOSITION
  ! weta, wetb              parameters for determining wet scavenging coefficients

  ! GAS DEPOSITION
  ! reldiff                 diffusivitiy of species relative to diff. of H2O
  ! henry [M/atm]           Henry constant
  ! f0                      reactivity relative to that of O3
  ! ri [s/m]                stomatal resistance
  ! rcl [s/m]               lower canopy resistance
  ! rgs [s/m]               ground resistance
  ! rlu [s/m]               leaf cuticular resistance
  ! rm [s/m]                mesophyll resistance
  ! dryvel [m/s]            constant dry deposition velocity

  ! PARTICLE DEPOSITION
  ! density [kg/m3]         density of particles
  ! dquer [m]               mean diameter of particles
  ! dsigma                  dsigma=10 or dsigma=0.1 means that 68% of the
  !                    mass are between 0.1*dquer and 10*dquer

  ! fract                   mass fraction of each diameter interval
  ! vset [m/s]              gravitational settling velocity in ni intervals
  ! cunningham              Cunningham slip correction (strictly valid only near surface)
  ! vsetaver [m/s]          average gravitational settling velocity
  ! schmi                   Schmidt number**2/3 of each diameter interval
  ! weightmolar [g/mol]     molecular weight

  ! TIME VARIATION OF EMISSION
  ! area_hour, point_hour   daily variation of emission strengths for area and point sources
  ! area_dow, point_dow     day-of-week variation of emission strengths for area and point sources



  !**********************************************************
  ! Variables used for domain-filling trajectory calculations
  !**********************************************************

  integer :: nx_we(2),ny_sn(2)
  integer :: numcolumn
  integer :: numcolumn_we(2,0:nymax-1),numcolumn_sn(2,0:nxmax-1)
  real :: zcolumn_we(2,0:nymax-1,maxcolumn)
  real :: zcolumn_sn(2,0:nxmax-1,maxcolumn)
  real :: xmassperparticle
  real :: acc_mass_we(2,0:nymax-1,maxcolumn)
  real :: acc_mass_sn(2,0:nxmax-1,maxcolumn)

  ! nx_we(2)                x indices of western and eastern boundary of domain-filling
  ! ny_sn(2)                y indices of southern and northern boundary of domain-filling
  ! numcolumn_we            number of particles to be released within one column
  !                    at the western and eastern boundary surfaces
  ! numcolumn_sn            same as numcolumn_we, but for southern and northern domain boundary
  ! numcolumn               maximum number of particles to be released within a single
  !                    column
  ! zcolumn_we              altitudes where particles are to be released
  !                    at the western and eastern boundary surfaces
  ! zcolumn_sn              same as zcolumn_we, but for southern and northern domain boundary
  ! xmassperparticle        air mass per particle in the domain-filling traj. option
  ! acc_mass_we             mass that has accumulated at the western and eastern boundary;
  !                    if it exceeds xmassperparticle, a particle is released and
  !                    acc_mass_we is reduced accordingly
  ! acc_mass_sn             same as acc_mass_we, but for southern and northern domain boundary



  !******************************************************************************
  ! Variables associated with the ECMWF meteorological input data ("wind fields")
  !******************************************************************************

  integer :: numbwf,wftime(maxwf),lwindinterv,wfdt(maxwf)
  character(len=255) :: wfname(maxwf),wfspec(maxwf)

  ! lwindinterv [s]         Interval between wind fields currently in memory
  ! numbwf                  actual number of wind fields
  ! wftime(maxwf) [s]       times relative to beginning time of wind fields
  ! wfname(maxwf)           file names of wind fields
  ! wfspec(maxwf)           specifications of wind field file, e.g. if on hard
  !                    disc or on tape

  integer :: memtime(2),memind(2)

  ! memtime [s]             validation times of wind fields in memory
  ! memind                  pointer to wind field, in order to avoid shuffling
  !                    of wind fields



  !****************************************************************************
  ! Variables defining actual size and geographical location of the wind fields
  !****************************************************************************

  integer :: nx,ny,nxmin1,nymin1,nxfield,nuvz,nwz,nz,nmixz,nlev_ec
  real :: dx,dy,xlon0,ylat0,dxconst,dyconst,height(nzmax)
  integer :: add_sfc_level,itimeold
  integer :: map_proj_id, map_proj_method
  real :: xmet0,ymet0,coefdx
  real :: heightmid(nzmax)
  real :: xlon2d(0:nxmax-1,0:nymax-1),ylat2d(0:nxmax-1,0:nymax-1)
  real :: map_stdlon, map_truelat1, map_truelat2
  integer :: mp_physics

  ! nx,ny,nz                actual dimensions of wind fields in x,y and z
  !                    direction, respectively
  ! nxmin1,nymin1           nx-1, ny-1, respectively
  ! nuvz,nwz                vertical dimension of original ECMWF data
  ! nxfield                 same as nx for limited area fields,
  !                    but for global fields nx=nxfield+1
  ! nmixz                   number of levels up to maximum PBL height (3500 m)

  ! nuvz is used for u,v components
  ! nwz is used for w components (staggered grid)
  ! nz is used for the levels in transformed coordinates (terrain-following Cartesian
  ! coordinates)

  ! nlev_ec  number of levels ECMWF model
  ! dx                      grid distance in x direction
  ! dy                      grid distance in y direction
  ! dxconst,dyconst         auxiliary variables for utransform,vtransform
  ! height                  heights of all levels
  ! xlon0                   geographical longitude and
  ! ylat0                   geographical latitude of lower left grid point



  !*************************************************
  ! Variables used for vertical model discretization
  !*************************************************

  real :: akm(nwzmax),bkm(nwzmax)
  real :: akz(nuvzmax),bkz(nuvzmax)
!  real :: aknew(nzmax),bknew(nzmax)
  real :: eta_w_wrf(nzmax), eta_u_wrf(nzmax), p_top_wrf

  ! akm,bkm: coeffizients which regulate vertical discretization of ecmwf model
  !     (at the border of model layers)
  ! akz,bkz: model discretization coeffizients at the centre of the layers
  ! aknew,bknew model discretization coeffizients at the interpolated levels



  ! Fixed fields, unchangeable with time
  !*************************************

  real :: oro(0:nxmax-1,0:nymax-1)
  real :: excessoro(0:nxmax-1,0:nymax-1)
  real :: lsm(0:nxmax-1,0:nymax-1)

  ! oro [m]              orography of the ECMWF model
  ! excessoro            excess orography mother domain
  ! lsm                  land sea mask of the ECMWF model
  ! 3d fields
  !**********

  real :: uu(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: vv(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: uupol(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: vvpol(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: ww(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: tt(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: qv(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: pv(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: rho(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: drhodz(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: tth(0:nxmax-1,0:nymax-1,nuvzmax,2)
  real :: qvh(0:nxmax-1,0:nymax-1,nuvzmax,2)
  real :: pph(0:nxmax-1,0:nymax-1,nuvzmax,2)
  real :: zzh(0:nxmax-1,0:nymax-1,nwzmax+1,2)
  real :: tkeh(0:nxmax-1,0:nymax-1,nuvzmax,2)
  real :: tke(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: ptth(0:nxmax-1,0:nymax-1,nuvzmax,2)
  real :: ptt(0:nxmax-1,0:nymax-1,nzmax,2)
  real :: u_wrf(0:nxmax-1,0:nymax-1,nuvzmax,2)
  real :: v_wrf(0:nxmax-1,0:nymax-1,nuvzmax,2)
  real :: w_wrf(0:nxmax-1,0:nymax-1,nwzmax,2)
! CDA added for new wet depo scheme
  integer :: icloudbot(0:nxmax-1,0:nymax-1,2)
  integer :: icloudthck(0:nxmax-1,0:nymax-1,2)


!  real :: pplev(0:nxmax-1,0:nymax-1,nuvzmax,2)
! CDA Commented old cloud variables
!  integer(kind=1) :: clouds(0:nxmax-1,0:nymax-1,nzmax,2)
!  integer :: cloudsh(0:nxmax-1,0:nymax-1,2)

  ! uu,vv,ww [m/2]       wind components in x,y and z direction
  ! uupol,vvpol [m/s]    wind components in polar stereographic projection
  ! tt [K]               temperature data
  ! qv                   specific humidity data
  ! pv (pvu)             potential vorticity
  ! rho [kg/m3]          air density
  ! drhodz [kg/m2]       vertical air density gradient
  ! tth,qvh              tth,qvh on original eta levels
  ! clouds:   no cloud, no precipitation   0
  !      cloud, no precipitation      1
  !      rainout  conv/lsp dominated  2/3
  !      washout  conv/lsp dominated  4/5
  ! pplev for the GFS version

  ! 2d fields
  !**********

  real :: xlanduse(0:nxmax-1,0:nymax-1,numclass)  ! no interpolation!
  real :: ps(0:nxmax-1,0:nymax-1,1,2)
  real :: sd(0:nxmax-1,0:nymax-1,1,2)
  real :: msl(0:nxmax-1,0:nymax-1,1,2)
  real :: tcc(0:nxmax-1,0:nymax-1,1,2)
  real :: u10(0:nxmax-1,0:nymax-1,1,2)
  real :: v10(0:nxmax-1,0:nymax-1,1,2)
  real :: tt2(0:nxmax-1,0:nymax-1,1,2)
  real :: m_x(0:nxmax-1,0:nymax-1,1)
  real :: m_y(0:nxmax-1,0:nymax-1,1)
  real :: td2(0:nxmax-1,0:nymax-1,1,2)
  real :: lsprec(0:nxmax-1,0:nymax-1,1,2)
  real :: convprec(0:nxmax-1,0:nymax-1,1,2)
  real :: sshf(0:nxmax-1,0:nymax-1,1,2)
  real :: ssr(0:nxmax-1,0:nymax-1,1,2)
  real :: surfstr(0:nxmax-1,0:nymax-1,1,2)
  real :: ustar(0:nxmax-1,0:nymax-1,1,2)
  real :: wstar(0:nxmax-1,0:nymax-1,1,2)
  real :: hmix(0:nxmax-1,0:nymax-1,1,2)
  real :: tropopause(0:nxmax-1,0:nymax-1,1,2)
  real :: oli(0:nxmax-1,0:nymax-1,1,2)
  real :: diffk(0:nxmax-1,0:nymax-1,1,2)

  ! xlanduse [0-1]       area fractions in percent
  ! ps                   surface pressure
  ! sd                   snow depth
  ! msl                  mean sea level pressure
  ! tcc                  total cloud cover
  ! u10                  10 meter u
  ! v10                  10 meter v
  ! tt2                  2 meter temperature
  ! td2                  2 meter dew point
  ! lsprec [mm/h]        large scale total precipitation
  ! convprec [mm/h]      convective precipitation
  ! sshf                 surface sensible heat flux
  ! ssr                  surface solar radiation
  ! surfstr              surface stress
  ! ustar [m/s]          friction velocity
  ! wstar [m/s]          convective velocity scale
  ! hmix  [m]            mixing height
  ! tropopause [m]       altitude of thermal tropopause
  ! oli [m]              inverse Obukhov length (1/L)
  ! diffk [m2/s]         diffusion coefficient at reference height


  real :: vdep(0:nxmax-1,0:nymax-1,maxspec,2)

  ! vdep [m/s]           deposition velocities


  !********************************************************************
  ! Variables associated with the ECMWF input data (nested wind fields)
  !********************************************************************

  ! NOTE: all nested variables have the same name as the variables used
  ! for the mother domain, except with a 'n' appended at the end
  !********************************************************************

  integer :: numbnests

  ! numbnests    number of nested grids

  character(len=255) :: wfnamen(maxnests,maxwf)
  character(len=18) :: wfspecn(maxnests,maxwf)

  ! wfnamen      nested wind field names
  ! wfspecn      specifications of wind field file, e.g. if on hard
  !         disc or on tape


  !*********************************************************************
  ! Variables characterizing size and location of the nested wind fields
  !*********************************************************************

  integer :: nxn(maxnests),nyn(maxnests)
  real :: dxn(maxnests),dyn(maxnests),xlon0n(maxnests),ylat0n(maxnests)
  real :: xmet0n(maxnests),ymet0n(maxnests)
      real :: xlon2dn(0:nxmaxn-1,0:nymaxn-1,maxnests)
      real :: ylat2dn(0:nxmaxn-1,0:nymaxn-1,maxnests)
      integer :: m_grid_id(0:maxnests), m_parent_grid_id(0:maxnests)
      integer :: m_parent_grid_ratio(0:maxnests)
      integer :: l_parent_nest_id(0:maxnests)
      integer :: i_parent_start(0:maxnests), j_parent_start(0:maxnests)
      integer :: mp_physicsn(maxnests)


  ! nxn,nyn      actual dimensions of nested wind fields in x and y direction
  ! dxn,dyn      grid distances in x,y direction for the nested grids
  ! xlon0n       geographical longitude of lower left grid point of nested wind fields
  ! ylat0n       geographical latitude of lower left grid point of nested wind fields


  ! Nested fields, unchangeable with time
  !**************************************

  real :: oron(0:nxmaxn-1,0:nymaxn-1,maxnests)
  real :: excessoron(0:nxmaxn-1,0:nymaxn-1,maxnests)
  real :: lsmn(0:nxmaxn-1,0:nymaxn-1,maxnests)

  ! 3d nested fields
  !*****************

  real :: uun(0:nxmaxn-1,0:nymaxn-1,nzmax,2,maxnests)
  real :: vvn(0:nxmaxn-1,0:nymaxn-1,nzmax,2,maxnests)
  real :: wwn(0:nxmaxn-1,0:nymaxn-1,nzmax,2,maxnests)
  real :: ttn(0:nxmaxn-1,0:nymaxn-1,nzmax,2,maxnests)
  real :: qvn(0:nxmaxn-1,0:nymaxn-1,nzmax,2,maxnests)
  real :: pvn(0:nxmaxn-1,0:nymaxn-1,nzmax,2,maxnests)
!  integer(kind=1) :: cloudsn(0:nxmaxn-1,0:nymaxn-1,0:nzmax,2,maxnests)
!  integer :: cloudsnh(0:nxmaxn-1,0:nymaxn-1,2,maxnests)

! CDA added for new wet depo scheme
  integer :: icloudbotn(0:nxmaxn-1,0:nymaxn-1,2,maxnests)
  integer :: icloudthckn(0:nxmaxn-1,0:nymaxn-1,2,maxnests)

  real :: rhon(0:nxmaxn-1,0:nymaxn-1,nzmax,2,maxnests)
  real :: drhodzn(0:nxmaxn-1,0:nymaxn-1,nzmax,2,maxnests)
  real :: tthn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,2,maxnests)
  real :: qvhn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,2,maxnests)
      real :: pphn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,2,maxnests)
      real :: zzhn(0:nxmaxn-1,0:nymaxn-1,nwzmax+1,2,maxnests)
      real :: tkehn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,2,maxnests)
      real :: tken(0:nxmaxn-1,0:nymaxn-1,nzmax,2,maxnests)
      real :: ptthn(0:nxmaxn-1,0:nymaxn-1,nuvzmax,2,maxnests)
      real :: pttn(0:nxmaxn-1,0:nymaxn-1,nzmax,2,maxnests)

      real :: un_wrf(0:nxmaxn-1,0:nymaxn-1,nuvzmax,2,maxnests)
      real :: vn_wrf(0:nxmaxn-1,0:nymaxn-1,nuvzmax,2,maxnests)
      real :: wn_wrf(0:nxmaxn-1,0:nymaxn-1,nwzmax,2,maxnests)

  ! 2d nested fields
  !*****************

  real :: xlandusen(0:nxmaxn-1,0:nymaxn-1,numclass,maxnests)  ! no interpolation!
  real :: psn(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: sdn(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: msln(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: tccn(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: u10n(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: v10n(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: tt2n(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real(kind=4) :: m_yn(0:nxmaxn-1,0:nymaxn-1,1,maxnests)
  real(kind=4) :: m_xn(0:nxmaxn-1,0:nymaxn-1,1,maxnests)
  real :: td2n(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: lsprecn(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: convprecn(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: sshfn(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: ssrn(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: surfstrn(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: ustarn(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: wstarn(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: hmixn(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: tropopausen(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: olin(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: diffkn(0:nxmaxn-1,0:nymaxn-1,1,2,maxnests)
  real :: vdepn(0:nxmaxn-1,0:nymaxn-1,maxspec,2,maxnests)


  !*************************************************
  ! Certain auxiliary variables needed for the nests
  !*************************************************

  real :: xresoln(0:maxnests),yresoln(0:maxnests)

  ! xresoln, yresoln   Factors by which the resolutions in the nests
  !               are enhanced compared to mother grid

  real :: xln(maxnests),yln(maxnests),xrn(maxnests),yrn(maxnests)

  ! xln,yln,xrn,yrn    Corner points of nested grids in grid coordinates
  !               of mother grid


  !******************************************************
  ! Variables defining the polar stereographic projection
  !******************************************************

  logical :: xglobal,sglobal,nglobal
  real :: switchnorthg,switchsouthg

  !xglobal             T for global fields, F for limited area fields
  !sglobal             T if domain extends towards south pole
  !nglobal             T if domain extends towards north pole
  !switchnorthg,switchsouthg   same as parameters switchnorth,
  !                    switchsouth, but in grid units

  real :: southpolemap(9),northpolemap(9)

  !southpolemap,northpolemap   define stereographic projections
  !                    at the two poles


  !******************
  ! Landuse inventory
  ! Sabine Eckhardt Dec 06: change to new landuse inventary - 11 classes, 1200 x 600 global
  !******************

  integer(kind=1) :: landinvent(1200,600,6)
  real :: z0(numclass)

  ! landinvent          landuse inventory (numclass=11 classes)
  !                     The last dimension holds lu_categories (pos 1-3) with
  !                     corresponding fractions (pos 4-6)
  ! z0                  roughness length for the landuse classes

  !**************************************************************************
  ! Variables characterizing the output grid and containing the model results
  !**************************************************************************

  integer :: numxgrid,numygrid,numzgrid
  real :: dxout,dyout,outlon0,outlat0,xoutshift,youtshift
      real :: dyoutl, dxoutl,out_xm0,out_ym0
  real :: dxoutln,dyoutln
  integer :: numxgridn,numygridn
  real :: dxoutn,dyoutn,outlon0n,outlat0n,xoutshiftn,youtshiftn
      real :: outgrid_swlon,outgrid_swlat,outgrid_nelon,outgrid_nelat
      real :: outgridn_swlon,outgridn_swlat,outgridn_nelon,outgridn_nelat
  real :: out_xm0n,out_ym0n
  !real outheight(maxzgrid),outheighthalf(maxzgrid)
  logical :: DEP,DRYDEP,DRYDEPSPEC(maxspec),WETDEP,OHREA,ASSSPEC
  integer nccovid,ncravid,ncwdvid,ncddvid,ncrecvid,ncrcovid,ncrmivid

  ! numxgrid,numygrid       number of grid points in x,y-direction
  ! numxgridn,numygridn     number of grid points in x,y-direction for nested output grid
  ! numzgrid                number of vertical levels of output grid
  ! dxout,dyout             grid distance of output grid
  ! dxoutn,dyoutn           grid distance of nested output grid
  ! outlon0,outlat0         lower left corner of output grid
  ! outlon0n,outlat0n       lower left corner of nested output grid
  ! xoutshift,youtshift     xlon0-outlon0, ylat0-outlat0
  ! xoutshiftn,youtshiftn   xlon0-outlon0n, ylat0-outlat0n
  ! out_xm0,out_ym0         lower left corner of output grid in grid metres
  ! outheight [m]           upper levels of the output grid
  ! outheighthalf [m]       half (middle) levels of the output grid cells
  ! DEP                     .true., if either dry or wet depos. is switched on
  ! DRYDEP                  .true., if dry deposition is switched on
  ! DRYDEPSPEC              .true., if dry deposition is switched on for that species
  ! WETDEP                  .true., if wet deposition is switched on
  ! OHREA                   .true., if OH reaction is switched on
  ! ASSSPEC                 .true., if there are two species asscoiated
  !                    (i.e. transfer of mass between these two occurs
  ! nccovid                 netcdf variable ID of concentration grid
  ! ncravid                 netcdf variable ID of mixing ratio grid
  ! ncddvid                 netcdf variable ID of drydep grid
  ! ncwdvid                 netcdf variable ID of wetdep grid
  ! ncrcovid                netcdf variable ID of receptor concentration
  ! ncrmivid                netcdf variable ID of receptor mixing ratio

  !  if output for each releasepoint shall be created maxpointspec=number of releasepoints
  !  else maxpointspec is 1 -> moved to unc_mod
  !  the OUTGRID is moved to the module outg_mod
  !******************************************************************************

  !real gridunc(0:maxxgrid-1,0:maxygrid-1,maxzgrid,maxspec,
  !    +             maxpointspec_act,nclassunc,maxageclass)
  !real griduncn(0:maxxgridn-1,0:maxygridn-1,maxzgrid,maxspec,
  !    +              maxpointspec_act,nclassunc,maxageclass)
  !real wetgridunc(0:maxxgrid-1,0:maxygrid-1,maxspec,
  !    +                maxpointspec_act,nclassunc,maxageclass)
  !real wetgriduncn(0:maxxgridn-1,0:maxygridn-1,maxspec,
  !    +ct                 maxpointspec,nclassunc,maxageclass)
  !real drygridunc(0:maxxgrid-1,0:maxygrid-1,maxspec,maxpointspec,
  !    +                nclassunc,maxageclass)
  !real drygriduncn(0:maxxgridn-1,0:maxygridn-1,maxspec,
  !    +                 maxpointspec,nclassunc,maxageclass)

  !real oroout(0:maxxgrid-1,0:maxygrid-1)
  !real orooutn(0:maxxgridn-1,0:maxygridn-1)
  !     real area(0:maxxgrid-1,0:maxygrid-1)
  !real arean(0:maxxgridn-1,0:maxygridn-1)
  !real volume(0:maxxgrid-1,0:maxygrid-1,maxzgrid)
  !real volumen(0:maxxgridn-1,0:maxygridn-1,maxzgrid)

  !real areaeast(0:maxxgrid-1,0:maxygrid-1,maxzgrid)
  !real areanorth(0:maxxgrid-1,0:maxygrid-1,maxzgrid)


  ! gridunc,griduncn        uncertainty of outputted concentrations
  ! wetgridunc,wetgriduncn  uncertainty of accumulated wet deposited mass on output grid
  ! drygridunc,drygriduncn  uncertainty of accumulated dry deposited mass on output grid
  ! oroout,orooutn [m]      height of model topography at output grid
  ! area,arean [m2]         area of each grid cell
  ! volume,volumen [m3]     volume of each grid cell
  ! ... field names with n at the end indicate a nested output grid


  !***********************************
  ! Variables defining receptor points
  !***********************************

  real :: xreceptor(maxreceptor),yreceptor(maxreceptor)
  real :: receptor_lon(maxreceptor),receptor_lat(maxreceptor)
  real :: receptorarea(maxreceptor)
  real :: creceptor(maxreceptor,maxspec)
  character(len=16) :: receptorname(maxreceptor)
  integer :: numreceptor

  ! xreceptor,yreceptor     receptor position
  ! creceptor               concentrations at receptor points
  ! receptorarea            area of 1*1 grid cell at receptor point



  !***************************************
  ! Variables characterizing each particle
  !***************************************

  integer :: numpart
  integer :: numparticlecount

  real,allocatable, dimension (:,:) :: xmass1
! real,allocatable, dimension (:,:) :: drydep1
  real,allocatable, dimension(:) :: xmasssave
  integer,allocatable, dimension (:) :: itra1
  integer,allocatable, dimension (:) :: npoint,nclass
  integer,allocatable, dimension (:) :: idt,itramem,itrasplit
  real(kind=dp),allocatable, dimension (:) :: xtra1,ytra1
  real,allocatable, dimension (:) :: ztra1
  real,allocatable, dimension (:) :: uap,ucp,uzp
  real,allocatable, dimension (:) :: us,vs,ws
  integer(kind=2),allocatable, dimension (:) :: cbt
  
  ! numpart                       actual number of particles in memory
  ! itra1 (particles) [s]         temporal positions of the particles
  ! npoint(particles)             indicates the release point of each particle
  ! nclass (particles)            one of nclassunc classes to which the particle is attributed
  ! itramem (particles) [s]       memorized release times of the particles
  ! itrasplit (particles) [s]     next time when particle is to be split into two
  ! idt(part) [s]                 time step to be used for next integration
  ! numparticlecount              counts the total number of particles that have been released
  ! xtra1,ytra1,ztra1             spatial positions of the particles
  ! xmass1(npartot,species) [kg]  particle masses



  !*******************************************************
  ! Info table on available chemical species/radionuclides
  !*******************************************************

  !character*10 specname(maxtable)
  real :: decaytime(maxspec),wetscava(maxspec),wetscavb(maxspec)
  real :: drydiff(maxspec),dryhenry(maxspec),dryactiv(maxspec)
  real :: partrho(maxspec),partmean(maxspec),partsig(maxspec)
  real :: dryvelo(maxspec),weightmol(maxspec)

  ! specname            Name of chemical species/radionuclide
  ! decaytime           Half time of radionuclides
  ! wetscava, wetscavb  Parameters for calculating scavenging coefficients
  ! drydiff             diffusivitiy of species relative to diff. of H2O
  ! dryhenry [M/atm]    Henry constant
  ! dryactiv            reactivity relative to that of O3
  ! partrho [kg/m3]     density of particles
  ! partmean [m]        mean diameter of particles
  ! partsig [m]         mean stand. deviation of particle diameter
  ! dryvelo [cm/s]      constant dry deposition velocity
  ! weightmol [g/mol]   molecular weight
  ! ohreact             OH reaction rate



  !********************
  ! Random number field
  !********************

! real :: rannumb(maxrandomp)
!#if newrandomgen.EQ.1               ! added by mc: if not commented use a
!preprocessor to select 
   real :: rannumb(maxrand3)           !modified by mc^M
   real :: uniform_rannumb(maxrand3)   !modified by mc^M

!#else                                     ! if not commented use a preprocessor
!    real :: rannumb(maxrandomp)           !modified by mc^M
!   real :: uniform_rannumb(maxrandomp)   !modified by mc^M
!  integer :: nan_count(MAX_STREAM)=0    !modified by mc
!#endif                                  ! if not commented use a preprocessor



end module com_mod
