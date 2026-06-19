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
       long_name="Temperature"
       units_name="degree C"
    else if(ivar == 2)then
       varname="s"
       long_name="Salinity"
       units_name="-"
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
       long_name="Vertical velocity (sigma coordinate)"
       units_name="m s^-1"
    else if(ivar == 6)then
       varname="wr"
       long_name="Vertical velocity (z coordinate)"
       units_name="m s^-1"
    else
       write(*,*) "***Error: ivar in 3D ==> ",ivar
       stop
    end if

  end subroutine varname3d

  !----------------------------------------------------------

  subroutine varname_mlt(ivar,varname,long_name,units_name)

    implicit none

    !---IN
    integer,intent(in) :: ivar

    !---OUT
    character(*),intent(out) :: varname,long_name,units_name

    if(ivar == 1)then
       varname="dtdt"
       long_name="MLT tendency"
       units_name="degree C day^-1"
    else if(ivar == 2)then !***2D***
       varname="tsfc"
       long_name="Contributions from latent and sensible heat fluxes and net longwave radiation to MLT"
       units_name="degree C day^-1"
    else if(ivar == 3)then
       varname="qz"
       long_name="Contributions from net shortwave radiation to MLT"
       units_name="degree C day^-1"
    else if(ivar == 4)then
       varname="txadv"
       long_name="Zonal MLT advection"
       units_name="degree C day^-1"
    else if(ivar == 5)then
       varname="tyadv"
       long_name="Meridional MLT advection"
       units_name="degree C day^-1"
    else if(ivar == 6)then
       varname="tzadv"
       long_name="Vertical MLT advection"
       units_name="degree C day^-1"
    else if(ivar == 7)then
       varname="txdif"
       long_name="Zonal MLT diffusion"
       units_name="degree C day^-1"
    else if(ivar == 8)then
       varname="tydif"
       long_name="Meridional MLT diffusion"
       units_name="degree C day^-1"
    else if(ivar == 9)then
       varname="tzdif"
       long_name="Vertical MLT diffusion"
       units_name="degree C day^-1"
    else if(ivar == 10)then
       varname="tent"
       long_name="MLT detrainment and entrainment"
       units_name="degree C day^-1"
    else if(ivar == 11)then
       varname="tiau"
       long_name="MLT analysis increment"
       units_name="degree C day^-1"    
    else if(ivar == 12)then
       varname="tres"
       long_name="MLT residual"
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
       long_name="MLS tendency"
       units_name="day^-1"
    else if(ivar == 2)then !***2D***
       varname="ssfc"
       long_name="Contributions from freshwater flux to MLS"
       units_name="day^-1"
    else if(ivar == 3)then
       varname="sxadv"
       long_name="Zonal MLS advection"
       units_name="day^-1"
    else if(ivar == 4)then
       varname="syadv"
       long_name="Meridional MLS advection"
       units_name="day^-1"
    else if(ivar == 5)then
       varname="szadv"
       long_name="Vertical MLS advection"
       units_name="day^-1"
    else if(ivar == 6)then
       varname="sxdif"
       long_name="Zonal MLS diffusion"
       units_name="day^-1"
    else if(ivar == 7)then
       varname="sydif"
       long_name="Meridional MLS diffusion"
       units_name="day^-1"
    else if(ivar == 8)then
       varname="szdif"
       long_name="Vertical MLS diffusion"
       units_name="day^-1"
    else if(ivar == 9)then
       varname="sent"
       long_name="MLS detrainment and entrainment"
       units_name="day^-1"
    else if(ivar == 10)then
       varname="siau"
       long_name="MLS analysis increment"
       units_name="day^-1"    
    else if(ivar == 11)then
       varname="snudge"
       long_name="MLS nudging"
       units_name="day^-1"    
    else if(ivar == 12)then
       varname="sres"
       long_name="MLS residual"
       units_name="day^-1"    
    else
       write(*,*) "***Error: ivar in MLS ==> ",ivar
       stop
    end if

  end subroutine varname_mls

  !----------------------------------------------------------

  subroutine varname_mld(ivar,varname,long_name,units_name)

    implicit none

    !---IN
    integer,intent(in) :: ivar

    !---OUT
    character(*),intent(out) :: varname,long_name,units_name

    if(ivar == 1)then
       varname="mld"
       long_name="MLD"
       units_name="meter"
    else if(ivar == 2)then
       varname="mld_ent"
       long_name="MLD of detrainment and estrainment"
       units_name="meter"
    else if(ivar == 3)then
       varname="dhdt"
       long_name="MLD tendency"
       units_name="m day^-1"
    else if(ivar == 4)then
       varname="delta_t"
       long_name="Temperature difference between ML average and at the base of ML"
       units_name="degree C"
    else if(ivar == 5)then
       varname="delta_s"
       long_name="Salinity difference between ML average and at the base of ML"
       units_name="-"
    else
       write(*,*) "***Error: ivar in MLS ==> ",ivar
       stop
    end if

  end subroutine varname_mld

  !-------------------------------------------------------

  subroutine varname_ens(ivar,varname,long_name,units_name)

    implicit none

    !---IN
    integer,intent(in) :: ivar

    !---OUT
    character(*),intent(out) :: varname,long_name,units_name

    if(ivar == 1)then
       varname="h"
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
       units_name="-"
    else
       write(*,*) "***Error: ivar in Ensemble ==> ",ivar
       stop
    end if

  end subroutine varname_ens

end module mod_varname
