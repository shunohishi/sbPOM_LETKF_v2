!--------------------------------------------------------
! Variable Name |
!--------------------------------------------------------

module mod_varname

contains

  subroutine varname2d(ivar,varname,long_name,units_name)

    implicit none

    !---IN
    integer,intent(in) :: ivar

    !---OUT
    character(*),intent(out) :: varname,long_name,units_name

    if(ivar == 1)then
       varname="lhf"
       long_name="Latent heat flux"
       units_name="W m^-2"       
    else if(ivar == 2)then
       varname="shf"
       long_name="Sensible heat flux"
       units_name="W m^-2"       
    else if(ivar == 3)then
       varname="lwr"
       long_name="Net longwave radiation"
       units_name="W m^-2"       
    else if(ivar == 4)then
       varname="swr"
       long_name="Net shortwave radiation"
       units_name="W m^-2"
    else if(ivar == 5)then
       varname="windu"
       long_name="Surface zonal wind"
       units_name="m s^-1"
    else if(ivar == 6)then
       varname="windv"
       long_name="Surface meridional wind"
       units_name="m s^-1"       
    else if(ivar == 7)then
       varname="winds"
       long_name="Surface wind speed"
       units_name="m s^-1"       
    else if(ivar == 8)then
       varname="tauu"
       long_name="Surface zonal wind stress"
       units_name="N m^-2"       
    else if(ivar == 9)then
       varname="tauv"
       long_name="Surface meridional wind stress"
       units_name="N m^-2"              
    else if(ivar == 10)then
       varname="taus"
       long_name="Surface wind stress magnitude"
       units_name="N m^-2"                     
    else if(ivar == 11)then
       varname="qa"
       long_name="Surface air specific humidity"
       units_name="g kg^-1"                            
    else if(ivar == 12)then
       varname="qs"
       long_name="Surface saturated specific humidity"
       units_name="g kg^-1"                                   
    else if(ivar == 13)then
       varname="ta"
       long_name="Surface air temperature"
       units_name="degree C"
    else if(ivar == 14)then
       varname="tsfc"
       long_name="Contributions from latent and sensible heat fluxes and net longwave radiation to SST"
       units_name="degree C day^-1"
    else if(ivar == 15)then
       varname="ssfc"
       long_name="Contributions from freshwater flux to SSS"
       units_name="day^-1"
    else
       write(*,*) "***Error: ivar in 2D ==> ",ivar
       stop
    end if

  end subroutine varname2d

  !-------------------------------------------------------

  subroutine varname3d(ivar,varname,long_name,units_name)

    implicit none

    !---IN
    integer,intent(in) :: ivar

    !---OUT
    character(*),intent(out) :: varname,long_name,units_name

    if(ivar == 1)then
       varname="t"
       long_name="Potential temperature"
       units_name="degree C"
    else if(ivar == 2)then
       varname="s"
       long_name="Salinity"
       units_name="1"
    else if(ivar == 3)then
       varname="u"
       long_name="Zonal velocity"
       units_name="m s^-1"
    else if(ivar == 4)then
       varname="v"
       long_name="Meridional velocity"
       units_name="m s^-1"
    else if(ivar == 5)then
       varname="w"
       long_name="Vertical velocity in sigma coordinates"
       units_name="m s^-1"
    else if(ivar == 6)then
       varname="wr"
       long_name="Vertical velocity in z coordinates"
       units_name="m s^-1"
    else
       write(*,*) "***Error: ivar in 3D ==> ",ivar
       stop
    end if

  end subroutine varname3d

  !----------------------------------------------------------
  
  subroutine varname_ens(ivar,varname,long_name,units_name)

    implicit none

    !---IN
    integer,intent(in) :: ivar

    !---OUT
    character(*),intent(out) :: varname,long_name,units_name

    if(ivar == 1)then
       varname="el"
       long_name="Sea surface height"
       units_name="meter"
    else if(ivar == 2)then
       varname="u"
       long_name="Surface zonal velocity"
       units_name="m s^-1"
    else if(ivar == 3)then
       varname="v"
       long_name="Surface meridional velocity"
       units_name="m s^-1"
    else if(ivar == 4)then
       varname="t"
       long_name="Sea surface temperature"
       units_name="degree C"
    else if(ivar == 5)then
       varname="s"
       long_name="Sea surface salinity"
       units_name="1"
    else
       write(*,*) "***Error: ivar in Ensemble ==> ",ivar
       stop
    end if

  end subroutine varname_ens

  !----------------------------------------------------------

  subroutine varname_mld(ivar,varname,long_name,units_name)

    implicit none

    !---IN
    integer,intent(in) :: ivar

    !---OUT
    character(*),intent(out) :: varname,long_name,units_name

    if(ivar == 1)then
       varname="mld"
       long_name="Mixed-layer depth"
       units_name="meter"
    else if(ivar == 2)then
       varname="mld_ent"
       long_name="Mixed-layer depth used in entrainment and detrainment estimation"
       units_name="meter"
    else if(ivar == 3)then
       varname="dhdt"
       long_name="Mixed-layer depth tendency"
       units_name="m day^-1"
    else if(ivar == 4)then
       varname="delta_t"
       long_name="Difference between mixed-layer mean temperature and temperature at the base of mixed layer"
       units_name="degree C"
    else if(ivar == 5)then
       varname="delta_s"
       long_name="Difference between mixed-layer mean salinity and salinity at the base of mixed layer"
       units_name="1"
    else
       write(*,*) "***Error: ivar in MLD ==> ",ivar
       stop
    end if

  end subroutine varname_mld
  
  !----------------------------------------------------------

  subroutine varname_mlt(ivar,varname,long_name,units_name)

    implicit none

    !---IN
    integer,intent(in) :: ivar

    !---OUT
    character(*),intent(out) :: varname,long_name,units_name

    if(ivar == 1)then
       varname="dtdt"
       long_name="Mixed-layer temperature tendency"
       units_name="degree C day^-1"
    else if(ivar == 2)then !***2D***
       varname="tsfc"
       long_name="Contribution of non-shortwave surface heat fluxes to mixed-layer temperature tendency"
       units_name="degree C day^-1"
    else if(ivar == 3)then
       varname="qz"
       long_name="Contribution of shortwave radiation to mixed-layer temperature tendency"
       units_name="degree C day^-1"
    else if(ivar == 4)then
       varname="txadv"
       long_name="Contribution of zonal advection to mixed-layer temperature tendency"
       units_name="degree C day^-1"
    else if(ivar == 5)then
       varname="tyadv"
       long_name="Contribution of meridional advection to mixed-layer temperature tendency"
       units_name="degree C day^-1"
    else if(ivar == 6)then
       varname="tzadv"
       long_name="Contribution of vertical advection to mixed-layer temperature tendency"
       units_name="degree C day^-1"
    else if(ivar == 7)then
       varname="txdif"
       long_name="Contribution of zonal diffusion to mixed-layer temperature tendency"
       units_name="degree C day^-1"
    else if(ivar == 8)then
       varname="tydif"
       long_name="Contribution of meridional diffusion to mixed-layer temperature tendency"
       units_name="degree C day^-1"
    else if(ivar == 9)then
       varname="tzdif"
       long_name="Contribution of vertical diffusion to mixed-layer temperature tendency"
       units_name="degree C day^-1"
    else if(ivar == 10)then
       varname="tent"
       long_name="Contribution of detrainment and entrainment to mixed-layer temperature tendency"
       units_name="degree C day^-1"
    else if(ivar == 11)then
       varname="tiau"
       long_name="Contribution of analysis increment to mixed-layer temperature tendency"
       units_name="degree C day^-1"    
    else if(ivar == 12)then
       varname="tres"
       long_name="Residual of mixed-layer temperature budget"
       units_name="degree C day^-1"    
    else
       write(*,*) "***Error: ivar in MLT ==> ",ivar
       stop
    end if

  end subroutine varname_mlt

  !---------------------------------------------------------------

  subroutine varname_mls(ivar,varname,long_name,units_name)

    implicit none

    !---IN
    integer,intent(in) :: ivar

    !---OUT
    character(*),intent(out) :: varname,long_name,units_name

    if(ivar == 1)then
       varname="dsdt"
       long_name="Mixed-layer salinity tendency"
       units_name="day^-1"
    else if(ivar == 2)then !***2D***
       varname="ssfc"
       long_name="Contribution of freshwater flux to mixed-layer salinity tendency"
       units_name="day^-1"
    else if(ivar == 3)then
       varname="sxadv"
       long_name="Contribution of zonal advection to mixed-layer salinity tendency"
       units_name="day^-1"
    else if(ivar == 4)then
       varname="syadv"
       long_name="Contribution of meridional advection to mixed-layer salinity tendency"
       units_name="day^-1"
    else if(ivar == 5)then
       varname="szadv"
       long_name="Contribution of vertical advection to mixed-layer salinity tendency"
       units_name="day^-1"
    else if(ivar == 6)then
       varname="sxdif"
       long_name="Contribution of zonal diffusion to mixed-layer salinity tendency"
       units_name="day^-1"
    else if(ivar == 7)then
       varname="sydif"
       long_name="Contribution of meridional diffusion to mixed-layer salinity tendency"
       units_name="day^-1"
    else if(ivar == 8)then
       varname="szdif"
       long_name="Contribution of vertical diffusion to mixed-layer salinity tendency"
       units_name="day^-1"
    else if(ivar == 9)then
       varname="sent"
       long_name="Contribution of detrainment and entrainment to mixed-layer salinity tendency"
       units_name="day^-1"
    else if(ivar == 10)then
       varname="siau"
       long_name="Contribution of analysis increment to mixed-layer salinity tendency"
       units_name="day^-1"    
    else if(ivar == 11)then
       varname="snudge"
       long_name="Contribution of nudging toward observational climatology to mixed-layer salinity tendency"
       units_name="day^-1"    
    else if(ivar == 12)then
       varname="sres"
       long_name="Residual of mixed-layer salinity budget"
       units_name="day^-1"    
    else
       write(*,*) "***Error: ivar in MLS ==> ",ivar
       stop
    end if

  end subroutine varname_mls

end module mod_varname
