!-------------------------------------------------------------------------
! Define variable |
!-------------------------------------------------------------------------
subroutine define_var_netcdf(ncid,ndim,dim,varid,name,long_name,units_name)

  use netcdf
  implicit none

  integer status

  integer,intent(in) :: ncid
  integer,intent(in) :: ndim,dim(ndim)
  integer,intent(inout) :: varid

  character(*),intent(in) :: name,long_name,units_name

  !*** Define variable: Real Only ***
  status=nf90_def_var(ncid,trim(name),nf90_float,dim,varid)
  call handle_error_netcdf('nf90_def_var:'//trim(name),status,nf90_noerr)

  !***Deflate level: 5***
  status=nf90_def_var_deflate(ncid,varid,shuffle=1,deflate=1,deflate_level=5)
  call handle_error_netcdf('nf90_def_var_deflate:'//trim(name),status,nf90_noerr)

  !*** Name description ***
  status=nf90_put_att(ncid,varid,"long_name",trim(long_name))
  call handle_error_netcdf('nf90_put_att:'//trim(long_name),status,nf90_noerr)

  !*** Units description ***
  status=nf90_put_att(ncid,varid,"units",trim(units_name))
  call handle_error_netcdf('nf90_put_att:'//trim(units_name),status,nf90_noerr)

end subroutine define_var_netcdf
!_______________________________________________________________________
subroutine write_netcdf_var_sngl(ncid,im,jm,km,varname,glb)

  use mpi
  use netcdf
  use common_pom_var, only: r_sngl
  implicit none

  integer status,varid
  
  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,km

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(in) :: glb(im,jm,km)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_put_var(ncid,varid,glb,(/1,1,1/),(/im,jm,km/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)
  
end subroutine write_netcdf_var_sngl
!_______________________________________________________________________
subroutine write_pnetcdf_var_sngl(ncid,im,jm,km,varname,glb)

  use mpi
  use netcdf
  use common_pom_var, only: r_sngl,i_global,j_global
  implicit none

  integer status,varid
  
  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,km

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(in) :: glb(im,jm,km)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_var_par_access(ncid,varid,nf90_collective)
  call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  status=nf90_put_var(ncid,varid,glb,(/i_global(1),j_global(1),1/),(/im,jm,km/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)
  
end subroutine write_pnetcdf_var_sngl
!_______________________________________________________________________
subroutine write_netcdf_var0d_time_sngl(ncid,it,varname,glb)

  use mpi
  use netcdf
  use common_pom_var, only: r_sngl
  implicit none

  integer status,varid
  
  integer,intent(in) :: ncid
  integer,intent(in) :: it

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(in) :: glb(1)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_put_var(ncid,varid,glb,(/it/),(/1/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)
  
end subroutine write_netcdf_var0d_time_sngl
!_______________________________________________________________________
subroutine write_pnetcdf_var0d_time_sngl(ncid,it,varname,glb)

  use mpi
  use netcdf
  use common_pom_var, only: r_sngl
  implicit none

  integer status,varid
  
  integer,intent(in) :: ncid
  integer,intent(in) :: it

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(in) :: glb(1)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_var_par_access(ncid,varid,nf90_collective)
  call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  status=nf90_put_var(ncid,varid,glb,(/it/),(/1/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)
  
end subroutine write_pnetcdf_var0d_time_sngl
!_______________________________________________________________________
subroutine write_netcdf_var1d_time_sngl(ncid,im,it,varname,glb)

  use mpi
  use netcdf
  use common_pom_var, only: r_sngl
  implicit none

  integer status,varid
  
  integer,intent(in) :: ncid
  integer,intent(in) :: im,it

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(in) :: glb(im)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_put_var(ncid,varid,glb,(/1,it/),(/im,1/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)
  
end subroutine write_netcdf_var1d_time_sngl
!_______________________________________________________________________
subroutine write_pnetcdf_var1d_time_sngl(ncid,im,it,varname,glb)

  use mpi
  use netcdf
  use common_pom_var, only: r_sngl
  implicit none

  integer status,varid
  
  integer,intent(in) :: ncid
  integer,intent(in) :: im,it

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(in) :: glb(im)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_var_par_access(ncid,varid,nf90_collective)
  call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  status=nf90_put_var(ncid,varid,glb,(/1,it/),(/im,1/))
  call handle_error_netcdf('nf90_put_var:'//trim(varname),status,nf90_noerr)
  
end subroutine write_pnetcdf_var1d_time_sngl
!_______________________________________________________________________
subroutine write_pnetcdf_var1d_sngl(ncid,im,varname,glb)

  use mpi
  use netcdf
  use common_pom_var, only: r_sngl
  implicit none

  integer status,varid
  
  integer,intent(in) :: ncid
  integer,intent(in) :: im

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(in) :: glb(im)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_var_par_access(ncid,varid,nf90_collective)
  call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  status=nf90_put_var(ncid,varid,glb)
  call handle_error_netcdf('nf90_put_var:'//trim(varname),status,nf90_noerr)
  
end subroutine write_pnetcdf_var1d_sngl
!_______________________________________________________________________
subroutine write_netcdf_var2d_time_sngl(ncid,im,jm,it,varname,glb)

  use mpi
  use netcdf
  use common_pom_var, only: r_sngl
  implicit none

  integer status,varid
  
  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,it

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(in) :: glb(im,jm)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_put_var(ncid,varid,glb,(/1,1,it/),(/im,jm,1/))
  call handle_error_netcdf('nf90_put_var:'//trim(varname),status,nf90_noerr)
  
end subroutine write_netcdf_var2d_time_sngl
!_______________________________________________________________________
subroutine write_pnetcdf_var2d_time_sngl(ncid,im,jm,it,varname,glb)

  use mpi
  use netcdf
  use common_pom_var, only: r_sngl,i_global,j_global
  implicit none

  integer status,varid
  
  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,it

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(in) :: glb(im,jm)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_var_par_access(ncid,varid,nf90_collective)
  call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  status=nf90_put_var(ncid,varid,glb, &
       & (/i_global(1),j_global(1),it/),(/im,jm,1/))
  call handle_error_netcdf('nf90_put_var:'//trim(varname),status,nf90_noerr)
  
end subroutine write_pnetcdf_var2d_time_sngl
!_______________________________________________________________________
subroutine write_netcdf_var3d_time_sngl(ncid,im,jm,km,it,varname,glb)

  use mpi
  use netcdf
  use common_pom_var, only: r_sngl
  implicit none

  integer status,varid
  
  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,km,it

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(in) :: glb(im,jm,km)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_put_var(ncid,varid,glb,(/1,1,1,it/),(/im,jm,km,1/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)
  
end subroutine write_netcdf_var3d_time_sngl
!_______________________________________________________________________
subroutine write_pnetcdf_var3d_time_sngl(ncid,im,jm,km,it,varname,glb)

  use mpi
  use netcdf
  use common_pom_var, only: r_sngl,i_global,j_global
  implicit none

  integer status,varid
  
  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,km,it

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(in) :: glb(im,jm,km)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_var_par_access(ncid,varid,nf90_collective)
  call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  status=nf90_put_var(ncid,varid,glb, &
       & (/i_global(1),j_global(1),1,it/),(/im,jm,km,1/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)
  
end subroutine write_pnetcdf_var3d_time_sngl
!_______________________________________________________________________
subroutine write_output_netcdf

  use mpi
  use netcdf
  use common_pom_var
  implicit none

  integer status  
  integer ncid,varid
  integer x_dimid,y_dimid,z_dimid,time_dimid

  integer ndim
  integer,allocatable :: dim(:)

  integer nprint

  integer ierr
  
  character(120) netcdf_out_file

  real(kind = r_sngl),allocatable :: glb2d(:,:),glb3d(:,:,:)

  if(.not. lpnetcdf) &
       & allocate(glb2d(im_global,jm_global),glb3d(im_global,jm_global,kb))
  
  !--- create netcdf file ---------------------------------------------------------------
  nprint=nint(iint/dble(iprint))
  write(netcdf_out_file,'(''out/'',a,''.nc'')') trim(netcdf_file)
  if(my_task == master_task) &
       & write(*,'(/''writing file '',a)') trim(netcdf_out_file)

  if(nprint == 1)then

     if(lpnetcdf)then
        status=nf90_create(trim(netcdf_out_file),ior(nf90_netcdf4,nf90_mpiio),ncid, &
             & comm=mpi_comm_world,info=mpi_info_null)
        call handle_error_netcdf('nf90_create: '//trim(netcdf_out_file),status,nf90_noerr)
     else if(my_task == master_task)then
        status=nf90_create(trim(netcdf_out_file),nf90_netcdf4,ncid)
        call handle_error_netcdf('nf90_create: '//trim(netcdf_out_file),status,nf90_noerr)
     end if

     if(lpnetcdf .or. my_task == master_task)then

        !     define global attributes
        status=nf90_put_att(ncid,nf90_global,"title",trim(title))
        call handle_error_netcdf('nf90_put_att: title',status,nf90_noerr)
        status=nf90_put_att(ncid,nf90_global,"description","output file")
        call handle_error_netcdf('nf90_put_att: description',status,nf90_noerr)

        !     define dimensions
        status=nf90_def_dim(ncid,"x",im_global,x_dimid)
        call handle_error_netcdf('nf90_def_dim: x',status,nf90_noerr)

        status=nf90_def_dim(ncid,"y",jm_global,y_dimid)
        call handle_error_netcdf('nf90_def_dim: y',status,nf90_noerr)

        status=nf90_def_dim(ncid,"z",kb,z_dimid)
        call handle_error_netcdf('nf90_def_dim: z',status,nf90_noerr)

        status=nf90_def_dim(ncid,"time",nf90_unlimited,time_dimid)
        call handle_error_netcdf('nf90_def_dim: time',status,nf90_noerr)

        !     define variables and their attributes
        !1D
        ndim=1
        allocate(dim(ndim))
        dim(1)=time_dimid
        call define_var_netcdf(ncid,ndim,dim,varid,"time","time","days since "//time_start)

        deallocate(dim)

        !2D
        ndim=2
        allocate(dim(ndim))
        dim(1)=x_dimid
        dim(2)=y_dimid

        call define_var_netcdf(ncid,ndim,dim,varid,"east_u","longitude of u points","degree E")
        call define_var_netcdf(ncid,ndim,dim,varid,"east_v","longitude of v points","degree E")
        call define_var_netcdf(ncid,ndim,dim,varid,"east_e","longitude of elevation points","degree E")

        call define_var_netcdf(ncid,ndim,dim,varid,"north_u","latitude of u points","degree N")
        call define_var_netcdf(ncid,ndim,dim,varid,"north_v","latitude of v points","degree N")
        call define_var_netcdf(ncid,ndim,dim,varid,"north_e","latitude of elevation points","degree N")

        call define_var_netcdf(ncid,ndim,dim,varid,"h","undisturbed water depth","meter")

        call define_var_netcdf(ncid,ndim,dim,varid,"dum","u surface mask [u]","-")
        call define_var_netcdf(ncid,ndim,dim,varid,"dvm","v surface mask [v]","-")
        call define_var_netcdf(ncid,ndim,dim,varid,"fsm","elevation surface mask [el]","-")

        deallocate(dim)

        !3D
        ndim=3
        allocate(dim(ndim))
        dim(1)=x_dimid
        dim(2)=y_dimid
        dim(3)=z_dimid

        call define_var_netcdf(ncid,ndim,dim,varid,"z_w","sigma coordinate of w points","meter")
        call define_var_netcdf(ncid,ndim,dim,varid,"z_e","sigmacoordinate of tsuv points","meter")

        deallocate(dim)
        
        !2D+Time
        ndim=3
        allocate(dim(ndim))
        dim(1)=x_dimid
        dim(2)=y_dimid
        dim(3)=time_dimid

        call define_var_netcdf(ncid,ndim,dim,varid,"el","surface elevation [el]","meter")
        call define_var_netcdf(ncid,ndim,dim,varid,"lhf","latent heat flux [el]","W/m^2")
        call define_var_netcdf(ncid,ndim,dim,varid,"shf","sensible heat flux [el]","W/m^2")
        call define_var_netcdf(ncid,ndim,dim,varid,"lwr","longwave radiation [el]","W/m^2")
        call define_var_netcdf(ncid,ndim,dim,varid,"swr","shortwave radiation [el]","W/m^2")

        call define_var_netcdf(ncid,ndim,dim,varid,"windu","zonal wind [el]","m/s")
        call define_var_netcdf(ncid,ndim,dim,varid,"windv","meridional wind [el]","m/s")
        call define_var_netcdf(ncid,ndim,dim,varid,"winds","wind speed [el]","m/s")

        call define_var_netcdf(ncid,ndim,dim,varid,"tauu","zonal wind stress [el]","N/m^2")
        call define_var_netcdf(ncid,ndim,dim,varid,"tauv","meriidonal wind stress [el]","N/m^2")
        call define_var_netcdf(ncid,ndim,dim,varid,"taus","magnitude of wind stress [el]","N/m^2")

        call define_var_netcdf(ncid,ndim,dim,varid,"qa","air specific humidity [el]","g/kg")
        call define_var_netcdf(ncid,ndim,dim,varid,"qs","surface saturated specific humidity [el]","g/kg")
        call define_var_netcdf(ncid,ndim,dim,varid,"ta","air temperature [el]","degree C")

        if(issf == 1)then
           call define_var_netcdf(ncid,ndim,dim,varid,"evap","evaporation [el]","mm/day")
           call define_var_netcdf(ncid,ndim,dim,varid,"prep","precipitation [el]","mm/day")
           call define_var_netcdf(ncid,ndim,dim,varid,"river","river discharge [el]","mm/day")

           call define_var_netcdf(ncid,ndim,dim,varid,"eflux","evaporation flux [el]","m/s")
           call define_var_netcdf(ncid,ndim,dim,varid,"pflux","precipitation flux [el]","m/s")
           call define_var_netcdf(ncid,ndim,dim,varid,"rflux","river discharge flux [el]","m/s")        
        end if

        if(budget == 1)then
           call define_var_netcdf(ncid,ndim,dim,varid,"tsfc","surface heat flux forcing [el]","degree C/day")
           call define_var_netcdf(ncid,ndim,dim,varid,"ssfc","surface freshwater flux forcing [el]","1/day")
        end if

        deallocate(dim)

        !3D+Time
        ndim=4
        allocate(dim(ndim))
        dim(1)=x_dimid
        dim(2)=y_dimid
        dim(3)=z_dimid
        dim(4)=time_dimid

        call define_var_netcdf(ncid,ndim,dim,varid,"u","zonal velocity [u,zz]","m/s")
        call define_var_netcdf(ncid,ndim,dim,varid,"v","meridional velocity [v,zz]","m/s")
        call define_var_netcdf(ncid,ndim,dim,varid,"w","vertical velocity (sigma-cordinate) [el,z]","m/s")
        call define_var_netcdf(ncid,ndim,dim,varid,"wr","vertical velocity (z-cordinate) [el,z]","m/s")

        call define_var_netcdf(ncid,ndim,dim,varid,"t","potential temperature [el,zz]","degree C")
        call define_var_netcdf(ncid,ndim,dim,varid,"s","salinity [el,zz]","-")

        if(budget == 1)then

           call define_var_netcdf(ncid,ndim,dim,varid,"dtdt","temperature tendency [el,zz]","degree C/day")

           call define_var_netcdf(ncid,ndim,dim,varid,"txadv","zonal temperature advection [el,zz]","degree C/day")
           call define_var_netcdf(ncid,ndim,dim,varid,"tyadv","meridional temperature advection [el,zz]","degree C/day")
           call define_var_netcdf(ncid,ndim,dim,varid,"tzadv","vertical temperature advection [el,zz]","degree C/day")

           call define_var_netcdf(ncid,ndim,dim,varid,"txdif","zonal temperature diffusion [el,zz]","degree C/day")
           call define_var_netcdf(ncid,ndim,dim,varid,"tydif","meridional temperature diffusion [el,zz]","degree C/day")
           call define_var_netcdf(ncid,ndim,dim,varid,"tzdif","vertical temperature diffusion [el,zz]","degree C/day")

           call define_var_netcdf(ncid,ndim,dim,varid,"qz","shortwave penetration [el,zz]","degree C/day")

           call define_var_netcdf(ncid,ndim,dim,varid,"tnudge","temperature nudging [el,zz]","degree C/day")
!           call define_var_netcdf(ncid,ndim,dim,varid,"tres","temperature residual [el,zz]","degree C/day")

           call define_var_netcdf(ncid,ndim,dim,varid,"dsdt","salinity tendency [el,zz]","1/day")

           call define_var_netcdf(ncid,ndim,dim,varid,"sxadv","zonal salinity advection [el,zz]","1/day")
           call define_var_netcdf(ncid,ndim,dim,varid,"syadv","meridional salinity advection [el,zz]","1/day")
           call define_var_netcdf(ncid,ndim,dim,varid,"szadv","vertical salinity advection [el,zz]","1/day")


           call define_var_netcdf(ncid,ndim,dim,varid,"sxdif","zonal salinity diffusion [el,zz]","1/day")
           call define_var_netcdf(ncid,ndim,dim,varid,"sydif","meridional salinity diffusion [el,zz]","1/day")
           call define_var_netcdf(ncid,ndim,dim,varid,"szdif","vertical salinity diffusion [el,zz]","1/day")

           call define_var_netcdf(ncid,ndim,dim,varid,"snudge","salinity nudging [el,zz]","1/day")
!           call define_var_netcdf(ncid,ndim,dim,varid,"sres","Salinity residual [el,zz]","degree C/day")

           call define_var_netcdf(ncid,ndim,dim,varid,"aam","Am (Ah=0.2*Am) [el,zz]","m^2/s")
           call define_var_netcdf(ncid,ndim,dim,varid,"km","vertical kinematic viscosity [el,z]","m^2/s")
           call define_var_netcdf(ncid,ndim,dim,varid,"kh","vertical diffusivity [el,z]","m^2/s")

        end if

        if(budget == 1 .and. assim == 2)then
           call define_var_netcdf(ncid,ndim,dim,varid,"tiau","temperature analysis increment [el]","degree C/day")
           call define_var_netcdf(ncid,ndim,dim,varid,"siau","salinity analysis increment [el]","1/day")
        end if

        if(budget == 1 .and. lroff)then
           call define_var_netcdf(ncid,ndim,dim,varid,"troff","temperature rounding value [el]","degree C/day")
           call define_var_netcdf(ncid,ndim,dim,varid,"sroff","salinity rounding value [el]","1/day")
        end if

        deallocate(dim)

        !End definittion
        status=nf90_enddef(ncid)
        call handle_error_netcdf('nf90_enddef',status,nf90_noerr)
        status=nf90_close(ncid)
        call handle_error_netcdf('nf90_close',status,nf90_noerr)
        
     end if !lpnetcdf or mytask==master_task

  end if !nprint = 1

  call MPI_Barrier(pom_comm,ierr)
  
  !--- write data ------------------------------------------

  !Open
  if(lpnetcdf)then
     status=nf90_open(trim(netcdf_out_file),nf90_write,ncid, &
          & comm=mpi_comm_world,info=mpi_info_null)
     call handle_error_netcdf('nf90_open',status,nf90_noerr)
  else if(my_task == master_task)then
     status=nf90_open(trim(netcdf_out_file),nf90_write,ncid)
     call handle_error_netcdf('nf90_open',status,nf90_noerr)
  end if

  !Time
  if(lpnetcdf)then
     call write_pnetcdf_var0d_time_sngl(ncid,nprint,"time",real(time-0.5d0))
  else if(my_task == master_task)then
     call write_netcdf_var1d_time_sngl(ncid,1,nprint,"time",real(time-0.5d0))
  end if
  
  if(iint/iprint == 1)then
     
     !2D
     if(lpnetcdf)then

        call write_pnetcdf_var_sngl(ncid,im,jm,1,"east_u",real(east_u))

        call write_pnetcdf_var_sngl(ncid,im,jm,1,"east_v",real(east_v))

        call write_pnetcdf_var_sngl(ncid,im,jm,1,"east_e",real(east_e))

        call write_pnetcdf_var_sngl(ncid,im,jm,1,"north_u",real(north_u))

        call write_pnetcdf_var_sngl(ncid,im,jm,1,"north_v",real(north_v))

        call write_pnetcdf_var_sngl(ncid,im,jm,1,"north_e",real(north_e))

        call write_pnetcdf_var_sngl(ncid,im,jm,1,"h",real(h))

        call write_pnetcdf_var_sngl(ncid,im,jm,1,"fsm",real(fsm))

        call write_pnetcdf_var_sngl(ncid,im,jm,1,"dum",real(dum))

        call write_pnetcdf_var_sngl(ncid,im,jm,1,"dvm",real(dvm))
        
     else

        call merge2d(glb2d,east_u)
        if(my_task == master_task) &
             & call write_netcdf_var_sngl(ncid,im_global,jm_global,1,"east_u",real(glb2d))

        call merge2d(glb2d,east_v)
        if(my_task == master_task) &
             & call write_netcdf_var_sngl(ncid,im_global,jm_global,1,"east_v",real(glb2d))

        call merge2d(glb2d,east_e)
        if(my_task == master_task) &
             & call write_netcdf_var_sngl(ncid,im_global,jm_global,1,"east_e",real(glb2d))

        call merge2d(glb2d,north_u)
        if(my_task == master_task) &
             & call write_netcdf_var_sngl(ncid,im_global,jm_global,1,"north_u",real(glb2d))

        call merge2d(glb2d,north_v)
        if(my_task == master_task) &
             & call write_netcdf_var_sngl(ncid,im_global,jm_global,1,"north_v",real(glb2d))

        call merge2d(glb2d,north_e)
        if(my_task == master_task) &
             & call write_netcdf_var_sngl(ncid,im_global,jm_global,1,"north_e",real(glb2d))

        call merge2d(glb2d,h)
        if(my_task == master_task) &
             & call write_netcdf_var_sngl(ncid,im_global,jm_global,1,"h",real(glb2d))

        call merge2d(glb2d,fsm)
        if(my_task == master_task) &
             & call write_netcdf_var_sngl(ncid,im_global,jm_global,1,"fsm",real(glb2d))

        call merge2d(glb2d,dum)
        if(my_task == master_task) &
             & call write_netcdf_var_sngl(ncid,im_global,jm_global,1,"dum",real(glb2d))

        call merge2d(glb2d,dvm)
        if(my_task == master_task) &
             & call write_netcdf_var_sngl(ncid,im_global,jm_global,1,"dvm",real(glb2d))

     end if !lpnetcdf

     !3D
     if(lpnetcdf)then
        call write_pnetcdf_var_sngl(ncid,im,jm,kb,"z_w",real(z))
        call write_pnetcdf_var_sngl(ncid,im,jm,kb,"z_e",real(zz))
     else
        call merge3d(glb3d,z)
        if(my_task == master_task) &
             & call write_netcdf_var_sngl(ncid,im_global,jm_global,kb,"z_w",real(glb3d))
        call merge3d(glb3d,zz)
        if(my_task == master_task) &
             & call write_netcdf_var_sngl(ncid,im_global,jm_global,kb,"z_e",real(glb3d))
     end if
     
  end if !iint/iprint
  
  !2D+Time
  if(lpnetcdf)then
     
     if(idave == 1)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"el",real(el_dave))
     else if(idave == 2)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"el",real(el))
     end if

     if(idave == 1)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"lhf",real(lhf_dave))
     else if(idave == 2)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"lhf",real(lhf))
     end if

     if(idave == 1)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"shf",real(shf_dave))
     else if(idave == 2)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"shf",real(shf))
     end if

     if(idave == 1)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"lwr",real(lwr_dave))
     else if(idave == 2)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"lwr",real(lwr))
     end if

     if(idave == 1)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"swr",real(swr_dave))
     else if(idave == 2)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"swr",real(swr))
     end if

     if(idave == 1)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"windu",real(windu_dave))
     else if(idave == 2)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"windu",real(windu))
     end if

     if(idave == 1)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"windv",real(windv_dave))
     else if(idave == 2)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"windv",real(windv))
     end if
     
     if(idave == 1)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"winds",real(winds_dave))
     else if(idave == 2)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"winds",real(winds))
     end if

     if(idave == 1)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"tauu",real(tauu_dave))
     else if(idave == 2)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"tauu",real(tauu))
     end if

     if(idave == 1)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"tauv",real(tauv_dave))
     else if(idave == 2)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"tauv",real(tauv))
     end if
     
     if(idave == 1)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"taus",real(taus_dave))
     else if(idave == 2)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"taus",real(taus))
     end if

     if(idave == 1)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"qa",real(qa_dave))
     else if(idave == 2)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"qa",real(qa))
     end if

     if(idave == 1)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"qs",real(qs_dave))
     else if(idave == 2)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"qs",real(qs))
     end if

     if(idave == 1)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"ta",real(ta_dave))
     else if(idave == 2)then
        call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"ta",real(ta))
     end if
     
     if(issf == 1)then

        if(idave == 1)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"evap",real(evap_dave))
        else if(idave == 2)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"evap",real(evap))
        end if

        if(idave == 1)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"prep",real(prep_dave))
        else if(idave == 2)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"prep",real(prep))
        end if


        if(idave == 1)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"river",real(river_dave))
        else if(idave == 2)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"river",real(river))
        end if

        if(idave == 1)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"eflux",real(eflux_dave))
        else if(idave == 2)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"eflux",real(eflux))
        end if
        
        if(idave == 1)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"pflux",real(pflux_dave))
        else if(idave == 2)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"pflux",real(pflux))
        end if

        if(idave == 1)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"rflux",real(rflux_dave))
        else if(idave == 2)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"rflux",real(rflux))
        end if
        
     end if !issf

     if(budget == 1)then

        if(idave == 1)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"tsfc",real(tsfc_dave))
        else if(idave == 2)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"tsfc",real(tsfc))
        end if

        if(idave == 1)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"ssfc",real(ssfc_dave))
        else if(idave == 2)then
           call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"ssfc",real(ssfc))
        end if
        
     end if !budget
     
  else
  
     if(idave == 1)then
        call merge2d(glb2d,el_dave)
     else if(idave == 2)then
        call merge2d(glb2d,el)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"el",real(glb2d))

     if(idave == 1)then
        call merge2d(glb2d,lhf_dave)
     else if(idave == 2)then
        call merge2d(glb2d,lhf)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"lhf",real(glb2d))

     if(idave == 1)then
        call merge2d(glb2d,shf_dave)
     else if(idave == 2)then
        call merge2d(glb2d,shf)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"shf",real(glb2d))

     if(idave == 1)then
        call merge2d(glb2d,lwr_dave)
     else if(idave == 2)then
        call merge2d(glb2d,lwr)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"lwr",real(glb2d))

     if(idave == 1)then
        call merge2d(glb2d,swr_dave)
     else if(idave == 2)then
        call merge2d(glb2d,swr)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"swr",real(glb2d))

     if(idave == 1)then
        call merge2d(glb2d,windu_dave)
     else if(idave == 2)then
        call merge2d(glb2d,windu)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"windu",real(glb2d))

     if(idave == 1)then
        call merge2d(glb2d,windv_dave)
     else if(idave == 2)then
        call merge2d(glb2d,windv)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"windv",real(glb2d))

     if(idave == 1)then
        call merge2d(glb2d,winds_dave)
     else if(idave == 2)then
        call merge2d(glb2d,winds)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"winds",real(glb2d))

     if(idave == 1)then
        call merge2d(glb2d,tauu_dave)
     else if(idave == 2)then
        call merge2d(glb2d,tauu)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"tauu",real(glb2d))

     if(idave == 1)then
        call merge2d(glb2d,tauv_dave)
     else if(idave == 2)then
        call merge2d(glb2d,tauv)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"tauv",real(glb2d))

     if(idave == 1)then
        call merge2d(glb2d,taus_dave)
     else if(idave == 2)then
        call merge2d(glb2d,taus)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"taus",real(glb2d))

     if(idave == 1)then
        call merge2d(glb2d,qa_dave)
     else if(idave == 2)then
        call merge2d(glb2d,qa)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"qa",real(glb2d))

     if(idave == 1)then
        call merge2d(glb2d,qs_dave)
     else if(idave == 2)then
        call merge2d(glb2d,qs)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"qs",real(glb2d))

     if(idave == 1)then
        call merge2d(glb2d,ta_dave)
     else if(idave == 2)then
        call merge2d(glb2d,ta)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"ta",real(glb2d))

     if(issf == 1)then

        if(idave == 1)then
           call merge2d(glb2d,evap_dave)
        else if(idave == 2)then
           call merge2d(glb2d,evap)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"evap",real(glb2d))

        if(idave == 1)then
           call merge2d(glb2d,prep_dave)
        else if(idave == 2)then
           call merge2d(glb2d,prep)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"prep",real(glb2d))

        if(idave == 1)then
           call merge2d(glb2d,river_dave)
        else if(idave == 2)then
           call merge2d(glb2d,river)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"river",real(glb2d))

        if(idave == 1)then
           call merge2d(glb2d,eflux_dave)
        else if(idave == 2)then
           call merge2d(glb2d,eflux)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"eflux",real(glb2d))

        if(idave == 1)then
           call merge2d(glb2d,pflux_dave)
        else if(idave == 2)then
           call merge2d(glb2d,pflux)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"pflux",real(glb2d))

        if(idave == 1)then
           call merge2d(glb2d,rflux_dave)
        else if(idave == 2)then
           call merge2d(glb2d,rflux)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"rflux",real(glb2d))

     end if !issf

     if(budget == 1)then

        if(idave == 1)then
           call merge2d(glb2d,tsfc_dave)
        else if(idave == 2)then
           call merge2d(glb2d,tsfc)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"tsfc",real(glb2d))

        if(idave == 1)then
           call merge2d(glb2d,ssfc_dave)
        else if(idave == 2)then
           call merge2d(glb2d,ssfc)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"ssfc",real(glb2d))

     end if !budget
     
  end if !lpnetcdf


  if(lpnetcdf)then

     !3D+Time
     if(idave == 1)then
        call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"u",real(u_dave))
     else if(idave == 2)then
        call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"u",real(u))
     end if
     
     if(idave == 1)then
        call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"v",real(v_dave))
     else if(idave == 2)then
        call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"v",real(v))
     end if

     if(idave == 1)then
        call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"w",real(w_dave))
     else if(idave == 2)then
        call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"w",real(w))
     end if

     if(idave == 1)then
        call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"wr",real(wr_dave))
     else if(idave == 2)then
        call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"wr",real(wr))
     end if
     
     if(idave == 1)then
        call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"t",real(t_dave))
     else if(idave == 2)then
        call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"t",real(t))
     end if

     if(idave == 1)then
        call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"s",real(s_dave))
     else if(idave == 2)then
        call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"s",real(s))
     end if
     
     if(budget == 1)then

        !dT/dt
        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"dtdt",real(dtdt_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"dtdt",real(dtdt))
        end if
        
        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"txadv",real(txadv_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"txadv",real(txadv))
        end if
        
        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"tyadv",real(tyadv_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"tyadv",real(tyadv))
        end if

        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"tzadv",real(tzadv_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"tzadv",real(tzadv))
        end if
        
        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"txdif",real(txdif_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"txdif",real(txdif))
        end if
        
        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"tydif",real(tydif_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"tydif",real(tydif))
        end if

        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"tzdif",real(tzdif_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"tzdif",real(tzdif))
        end if
        
        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"qz",real(qz_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"qz",real(qz))
        end if
        
        !dS/dt
        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"dsdt",real(dsdt_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"dsdt",real(dsdt))
        end if
        
        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"sxadv",real(sxadv_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"sxadv",real(sxadv))
        end if

        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"syadv",real(syadv_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"syadv",real(syadv))
        end if

        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"szadv",real(szadv_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"szadv",real(szadv))
        end if

        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"sxdif",real(sxdif_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"sxdif",real(sxdif))
        end if
        
        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"sydif",real(sydif_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"sydif",real(sydif))
        end if

        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"szdif",real(szdif_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"szdif",real(szdif))
        end if
        
        !Viscoity/Diffusivity
        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"aam",real(aam_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"aam",real(aam))
        end if
        
        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"km",real(km_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"km",real(km))
        end if

        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"kh",real(kh_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"kh",real(kh))
        end if
        
     end if !budget

     if(budget == 1 .and. assim == 2)then

        call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"tiau",real(t_iau))
        call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"siau",real(s_iau))

     end if !budget & assim

     if(budget == 1 .and. (ts_nudge /= 0.d0 .or. ti_nudge /= 0.d0 .or. ss_nudge /= 0.d0 .or. si_nudge /= 0.d0))then

        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"tnudge",real(tnudge_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"tnudge",real(tnudge))
        end if
        
        if(idave == 1)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"snudge",real(snudge_dave))
        else if(idave == 2)then
           call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"snudge",real(snudge))
        end if
        
     end if
     
  else
  
     !3D+Time
     if(idave == 1)then
        call merge3d(glb3d,u_dave)
     else if(idave == 2)then
        call merge3d(glb3d,u)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"u",real(glb3d))

     if(idave == 1)then
        call merge3d(glb3d,v_dave)
     else if(idave == 2)then
        call merge3d(glb3d,v)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"v",real(glb3d))

     if(idave == 1)then
        call merge3d(glb3d,w_dave)
     else if(idave == 2)then
        call merge3d(glb3d,w)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"w",real(glb3d))

     if(idave == 1)then
        call merge3d(glb3d,wr_dave)
     else if(idave == 2)then
        call merge3d(glb3d,wr)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"wr",real(glb3d))

     if(idave == 1)then
        call merge3d(glb3d,t_dave)
     else if(idave == 2)then
        call merge3d(glb3d,t)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"t",real(glb3d))

     if(idave == 1)then
        call merge3d(glb3d,s_dave)
     else if(idave == 2)then
        call merge3d(glb3d,s)
     end if
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"s",real(glb3d))

     if(budget == 1)then

        !dT/dt
        if(idave == 1)then
           call merge3d(glb3d,dtdt_dave)
        else if(idave == 2)then
           call merge3d(glb3d,dtdt)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"dtdt",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,txadv_dave)
        else if(idave == 2)then
           call merge3d(glb3d,txadv)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"txadv",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,tyadv_dave)
        else if(idave == 2)then
           call merge3d(glb3d,tyadv)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"tyadv",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,tzadv_dave)
        else if(idave == 2)then
           call merge3d(glb3d,tzadv)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"tzadv",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,txdif_dave)
        else if(idave == 2)then
           call merge3d(glb3d,txdif)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"txdif",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,tydif_dave)
        else if(idave == 2)then
           call merge3d(glb3d,tydif)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"tydif",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,tzdif_dave)
        else if(idave == 2)then
           call merge3d(glb3d,tzdif)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"tzdif",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,qz_dave)
        else if(idave == 2)then
           call merge3d(glb3d,qz)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"qz",real(glb3d))

        !dS/dt
        if(idave == 1)then
           call merge3d(glb3d,dsdt_dave)
        else if(idave == 2)then
           call merge3d(glb3d,dsdt)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"dsdt",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,sxadv_dave)
        else if(idave == 2)then
           call merge3d(glb3d,sxadv)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"sxadv",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,syadv_dave)
        else if(idave == 2)then
           call merge3d(glb3d,syadv)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"syadv",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,szadv_dave)
        else if(idave == 2)then
           call merge3d(glb3d,szadv)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"szadv",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,sxdif_dave)
        else if(idave == 2)then
           call merge3d(glb3d,sxdif)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"sxdif",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,sydif_dave)
        else if(idave == 2)then
           call merge3d(glb3d,sydif)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"sydif",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,szdif_dave)
        else if(idave == 2)then
           call merge3d(glb3d,szdif)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"szdif",real(glb3d))

        !Viscoity/Diffusivity
        if(idave == 1)then
           call merge3d(glb3d,aam_dave)
        else if(idave == 2)then
           call merge3d(glb3d,aam)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"aam",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,km_dave)
        else if(idave == 2)then
           call merge3d(glb3d,km)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"km",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,kh_dave)
        else if(idave == 2)then
           call merge3d(glb3d,kh)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"kh",real(glb3d))

     end if

     if(budget == 1 .and. assim == 2)then

        call merge3d(glb3d,t_iau)
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"tiau",real(glb3d))

        call merge3d(glb3d,s_iau)
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"siau",real(glb3d))

     end if

     if(budget == 1 .and. (ts_nudge /= 0.d0 .or. ti_nudge /= 0.d0 .or. ss_nudge /= 0.d0 .or. si_nudge /= 0.d0))then

        if(idave == 1)then
           call merge3d(glb3d,tnudge_dave)
        else if(idave == 2)then
           call merge3d(glb3d,tnudge)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"tnudge",real(glb3d))

        if(idave == 1)then
           call merge3d(glb3d,snudge_dave)
        else if(idave == 2)then
           call merge3d(glb3d,snudge)
        end if
        if(my_task == master_task) &
             & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"snudge",real(glb3d))

     end if

  end if !lpnetcdf
     
  !     close file
  if(lpnetcdf .or. my_task == master_task)then
     status=nf90_close(ncid)
     call handle_error_netcdf('nf90_close: output',status,nf90_noerr)
  end if

  if(.not. lpnetcdf) deallocate(glb2d,glb3d)
  
end subroutine write_output_netcdf

!_______________________________________________________________________
subroutine write_restart_netcdf

  use mpi
  use netcdf
  use common_pom_var  
  implicit none

  integer,parameter :: nprint=1
  
  integer status,ncid,varid
  integer time_dimid,x_dimid,y_dimid,z_dimid

  integer ndim
  integer,allocatable :: dim(:)

  integer ierr
  
  real(kind = r_sngl),allocatable :: glb2d(:,:),glb3d(:,:,:)

  character(120) netcdf_out_file  

  if(.not. lpnetcdf) &
       & allocate(glb2d(im_global,jm_global),glb3d(im_global,jm_global,kb))
  
  !--- create netcdf restart file ------------------------------------------------------
  if(mod(iint,irestart) == 0)then
     netcdf_out_file="out/"//trim(write_rst_file)//".nc"
  else
     netcdf_out_file="out/"//trim(write_rst_file)//".nc"     
  end if

  if(my_task == master_task) &
       & write(*,'(/''writing file '',a)') trim(netcdf_out_file)

  if(lpnetcdf .or. my_task == master_task)then      

     if(lpnetcdf)then
        status=nf90_create(trim(netcdf_out_file),ior(nf90_netcdf4,nf90_mpiio),ncid, &
             & comm=mpi_comm_world,info=mpi_info_null)
        call handle_error_netcdf('nf90_create: '//trim(netcdf_out_file),status,nf90_noerr)
     else
        status=nf90_create(trim(netcdf_out_file),nf90_netcdf4,ncid)
        call handle_error_netcdf('nf90_create',status,nf90_noerr)
     end if

     status=nf90_put_att(ncid,nf90_global,"title",trim(title))
     call handle_error_netcdf('nf90_put_att: title',status,nf90_noerr)
     status=nf90_put_att(ncid,nf90_global,"description","restart file")
     call handle_error_netcdf('nf90_put_att: description',status,nf90_noerr)

     !     define dimensions
     status=nf90_def_dim(ncid,"x",im_global,x_dimid)
     call handle_error_netcdf('nf90_def_dim: x',status,nf90_noerr)

     status=nf90_def_dim(ncid,"y",jm_global,y_dimid)
     call handle_error_netcdf('nf90_def_dim: y',status,nf90_noerr)

     status=nf90_def_dim(ncid,"z",kb,z_dimid)
     call handle_error_netcdf('nf90_def_dim: z',status,nf90_noerr)

     status=nf90_def_dim(ncid,"time",nf90_unlimited,time_dimid)
     call handle_error_netcdf('nf90_def_dim: time',status,nf90_noerr)

     !     define variables and their attributes
     !1D
     ndim=1
     call define_var_netcdf(ncid,ndim,time_dimid,varid,"time","time","days since "//time_start)

     !2D+Time
     ndim=3
     allocate(dim(ndim))
     dim(1)=x_dimid
     dim(2)=y_dimid
     dim(3)=time_dimid

     call define_var_netcdf(ncid,ndim,dim,varid,"wubot","zonal momentum flux at the bottom [u]","m^2/s^2")
     call define_var_netcdf(ncid,ndim,dim,varid,"wvbot","meridional momentum flux at the bottom [v]","m^2/s^2")

     call define_var_netcdf(ncid,ndim,dim,varid,"aam2d","vertical averaged aam [el]","m^2/s")

     call define_var_netcdf(ncid,ndim,dim,varid,"ua","vertical averaged u [u]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"uab","vertical averaged u at time -dt [u]","m/s")

     call define_var_netcdf(ncid,ndim,dim,varid,"va","vertical averaged v [v]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"vab","vertical averaged v at time -dt [v]","m/s")

     call define_var_netcdf(ncid,ndim,dim,varid,"el","surface elevation in external mode [el]","meter")
     call define_var_netcdf(ncid,ndim,dim,varid,"elb","surface elevation in external mode at time -dt [el]","meter")

     call define_var_netcdf(ncid,ndim,dim,varid,"et","surface elevation in internal mode [el]","meter")
     call define_var_netcdf(ncid,ndim,dim,varid,"etb","surface elevation in internal mode at time -dt [el]","meter")

     call define_var_netcdf(ncid,ndim,dim,varid,"egb","surface elevation for pressure gradient at time -dt [el]","meter")

     call define_var_netcdf(ncid,ndim,dim,varid,"utb","ua time averaged over dti [u]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"vtb","va time averaged over dti [v]","m/s")

     call define_var_netcdf(ncid,ndim,dim,varid,"adx2d","vertical integrated advx [u]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"ady2d","vertical integrated advy [v]","m/s")

     call define_var_netcdf(ncid,ndim,dim,varid,"advua","sum of 2nd, 3rd, and 4th terms in Eq. (18) [u]","-")
     call define_var_netcdf(ncid,ndim,dim,varid,"advva","sum of 2nd, 3rd, and 4th terms in Eq. (19) [v]","-")

     deallocate(dim)

     !3D+Time
     ndim=4
     allocate(dim(ndim))
     dim(1)=x_dimid
     dim(2)=y_dimid
     dim(3)=z_dimid
     dim(4)=time_dimid

     call define_var_netcdf(ncid,ndim,dim,varid,"u","zonal velocity [u]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"ub","zonal velocity at time -dt [u]","m/s")

     call define_var_netcdf(ncid,ndim,dim,varid,"v","meridional velocity [v]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"vb","meridional velocity at time -dt [v]","m/s")

     call define_var_netcdf(ncid,ndim,dim,varid,"w","sigma-velocity [el,z]","m/s")

     call define_var_netcdf(ncid,ndim,dim,varid,"t","potential temperature [el,zz]","degree C")
     call define_var_netcdf(ncid,ndim,dim,varid,"tb","potential temperature at time -dt [el,zz]","degree C")

     call define_var_netcdf(ncid,ndim,dim,varid,"s","salinity [el,zz]","-")
     call define_var_netcdf(ncid,ndim,dim,varid,"sb","salinity at time -dt [el,zz]","-")

     call define_var_netcdf(ncid,ndim,dim,varid,"rho","(density-1000)/rhoref [el,zz]","-")

     call define_var_netcdf(ncid,ndim,dim,varid,"km","vertical kinematic viscosity [el,z]","m^2/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"kh","vertical diffusivity [el,z]","m^2/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"kq","kq [el,z]","m^2/s")

     call define_var_netcdf(ncid,ndim,dim,varid,"l","turbulence length scale [el,z]","-")

     call define_var_netcdf(ncid,ndim,dim,varid,"q2","twice the turbulent kinetic energy [el,z]","m^2/s^2")
     call define_var_netcdf(ncid,ndim,dim,varid,"q2b","twice the turbulent kinetic energy at time -dt [el,z]","m^2/s^2")
     call define_var_netcdf(ncid,ndim,dim,varid,"q2l","q2 x l [el,z]","m^3/s^2")
     call define_var_netcdf(ncid,ndim,dim,varid,"q2lb","q2 x l at time -dt [el,z]","m^3/s^2")

     call define_var_netcdf(ncid,ndim,dim,varid,"aam","horizontal kinematic viscosity [el,zz]","m^2/s")

     deallocate(dim)

     !     end definitions
     status=nf90_enddef(ncid)
     call handle_error_netcdf('nf90_enddef',status,nf90_noerr)
     status=nf90_close(ncid)
     call handle_error_netcdf('nf90_close',status,nf90_noerr)
     
  end if !lpnetcdf .or. master_task
 
  call MPI_Barrier(pom_comm,ierr)
 
  ! open netcdf
  if(lpnetcdf)then
     status=nf90_open(trim(netcdf_out_file),nf90_write,ncid, &
          & comm=mpi_comm_world,info=mpi_info_null)
     call handle_error_netcdf('nf90_open',status,nf90_noerr)     
  else if(my_task == master_task)then
     status=nf90_open(trim(netcdf_out_file),nf90_write,ncid)
     call handle_error_netcdf('nf90_open',status,nf90_noerr)
  end if
        
  ! write data
  !1D
  if(lpnetcdf)then
     call write_pnetcdf_var1d_time_sngl(ncid,1,nprint,"time",real(time))
  else if(my_task == master_task)then
     call write_netcdf_var1d_time_sngl(ncid,1,nprint,"time",real(time))
  end if
  
  !2D
  if(lpnetcdf)then

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"wubot",real(wubot))

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"wvbot",real(wvbot))

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"aam2d",real(aam2d))

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"ua",real(ua))

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"uab",real(uab))
     
     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"va",real(va))

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"vab",real(vab))
     
     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"el",real(el))

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"elb",real(elb))
     
     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"et",real(et))

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"etb",real(etb))

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"egb",real(egb))

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"utb",real(utb))

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"vtb",real(vtb))

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"adx2d",real(adx2d))

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"ady2d",real(ady2d))

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"advua",real(advua))

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"advva",real(advva))

  else
     
     call merge2d(glb2d,wubot)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"wubot",real(glb2d))

     call merge2d(glb2d,wvbot)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"wvbot",real(glb2d))

     call merge2d(glb2d,aam2d)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"aam2d",real(glb2d))

     call merge2d(glb2d,ua)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"ua",real(glb2d))

     call merge2d(glb2d,uab)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"uab",real(glb2d))
     
     call merge2d(glb2d,va)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"va",real(glb2d))

     call merge2d(glb2d,vab)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"vab",real(glb2d))
     
     call merge2d(glb2d,el)
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"el",real(glb2d))

     call merge2d(glb2d,elb)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"elb",real(glb2d))
     
     call merge2d(glb2d,et)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"et",real(glb2d))

     call merge2d(glb2d,etb)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"etb",real(glb2d))

     call merge2d(glb2d,egb)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"egb",real(glb2d))

     call merge2d(glb2d,utb)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"utb",real(glb2d))

     call merge2d(glb2d,vtb)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"vtb",real(glb2d))

     call merge2d(glb2d,adx2d)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"adx2d",real(glb2d))

     call merge2d(glb2d,ady2d)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"ady2d",real(glb2d))

     call merge2d(glb2d,advua)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"advua",real(glb2d))

     call merge2d(glb2d,advva)  
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"advva",real(glb2d))
     
  end if
  
  !3D
  if(lpnetcdf)then

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"u",real(u))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"ub",real(ub))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"v",real(v))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"vb",real(vb))
     
     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"w",real(w))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"t",real(t))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"tb",real(tb))
     
     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"s",real(s))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"sb",real(sb))
     
     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"rho",real(rho))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"km",real(km))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"kh",real(kh))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"kq",real(kq))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"l",real(l))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"q2",real(q2))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"q2b",real(q2b))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"q2l",real(q2l))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"q2lb",real(q2lb))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"aam",real(aam))

  else

     call merge3d(glb3d,u)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"u",real(glb3d))

     call merge3d(glb3d,ub)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"ub",real(glb3d))
     
     call merge3d(glb3d,v)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"v",real(glb3d))

     call merge3d(glb3d,vb)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"vb",real(glb3d))
     
     call merge3d(glb3d,w)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"w",real(glb3d))

     call merge3d(glb3d,t)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"t",real(glb3d))

     call merge3d(glb3d,tb)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"tb",real(glb3d))
     
     call merge3d(glb3d,s)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"s",real(glb3d))

     call merge3d(glb3d,sb)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"sb",real(glb3d))
     
     call merge3d(glb3d,rho)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"rho",real(glb3d))

     call merge3d(glb3d,km)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"km",real(glb3d))

     call merge3d(glb3d,kh)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"kh",real(glb3d))

     call merge3d(glb3d,kq)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"kq",real(glb3d))

     call merge3d(glb3d,l)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"l",real(glb3d))

     call merge3d(glb3d,q2)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"q2",real(glb3d))

     call merge3d(glb3d,q2b)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"q2b",real(glb3d))

     call merge3d(glb3d,q2l)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"q2l",real(glb3d))

     call merge3d(glb3d,q2lb)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"q2lb",real(glb3d))

     call merge3d(glb3d,aam)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"aam",real(glb3d))

  end if
     
  ! close file
  if(lpnetcdf .or. my_task == master_task) then
     status=nf90_close(ncid)
     call handle_error_netcdf('nf90_close',status,nf90_noerr)
  end if

  if(.not. lpnetcdf) deallocate(glb2d,glb3d)
  
end subroutine write_restart_netcdf

!_______________________________________________________________________
subroutine write_iau_netcdf

  use mpi
  use netcdf
  use common_pom_var  
  implicit none

  integer,parameter :: nprint=1
  
  integer status,ncid,varid
  integer time_dimid,x_dimid,y_dimid,z_dimid

  integer ndim
  integer,allocatable :: dim(:)

  integer ierr
  
  real(kind = r_sngl),allocatable :: glb2d(:,:),glb3d(:,:,:)

  character(120) netcdf_out_file  

  if(.not. lpnetcdf) &
       & allocate(glb2d(im_global,jm_global),glb3d(im_global,jm_global,kb))
  
  !--- create netcdf restart file ------------------------------------------------------
  netcdf_out_file="out/"//trim(write_iau_file)//".nc"

  if(my_task == master_task) &
       & write(*,'(/''writing file '',a)') trim(netcdf_out_file)

  if(lpnetcdf .or. my_task == master_task)then      

     if(lpnetcdf)then
        status=nf90_create(trim(netcdf_out_file),ior(nf90_netcdf4,nf90_mpiio),ncid, &
             & comm=mpi_comm_world,info=mpi_info_null)
        call handle_error_netcdf('nf90_create: '//trim(netcdf_out_file),status,nf90_noerr)
     else
        status=nf90_create(trim(netcdf_out_file),nf90_netcdf4,ncid)
        call handle_error_netcdf('nf90_create',status,nf90_noerr)
     end if

     status=nf90_put_att(ncid,nf90_global,"title",trim(title))
     call handle_error_netcdf('nf90_put_att: title',status,nf90_noerr)
     status=nf90_put_att(ncid,nf90_global,"description","iau file")
     call handle_error_netcdf('nf90_put_att: description',status,nf90_noerr)

     !     define dimensions
     status=nf90_def_dim(ncid,"x",im_global,x_dimid)
     call handle_error_netcdf('nf90_def_dim: x',status,nf90_noerr)

     status=nf90_def_dim(ncid,"y",jm_global,y_dimid)
     call handle_error_netcdf('nf90_def_dim: y',status,nf90_noerr)

     status=nf90_def_dim(ncid,"z",kb,z_dimid)
     call handle_error_netcdf('nf90_def_dim: z',status,nf90_noerr)

     status=nf90_def_dim(ncid,"time",nf90_unlimited,time_dimid)
     call handle_error_netcdf('nf90_def_dim: time',status,nf90_noerr)

     !     define variables and their attributes
     !1D
     ndim=1
     call define_var_netcdf(ncid,ndim,time_dimid,varid,"time","time","days since "//time_start)

     !2D+Time
     ndim=3
     allocate(dim(ndim))
     dim(1)=x_dimid
     dim(2)=y_dimid
     dim(3)=time_dimid

     call define_var_netcdf(ncid,ndim,dim,varid,"ua_fcst","forecast vertical averaged u [u]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"ua_anal","analysis vertical averaged u [u]","m/s")

     call define_var_netcdf(ncid,ndim,dim,varid,"va_fcst","forecast vertical averaged v [v]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"va_anal","analysis vertical averaged v [v]","m/s")

     call define_var_netcdf(ncid,ndim,dim,varid,"el_fcst", &
          & "forecast surface elevation in external mode [el]","meter")
     call define_var_netcdf(ncid,ndim,dim,varid,"el_anal", &
          & "analysis surface elevation in external mode [el]","meter")

     deallocate(dim)

     !3D+Time
     ndim=4
     allocate(dim(ndim))
     dim(1)=x_dimid
     dim(2)=y_dimid
     dim(3)=z_dimid
     dim(4)=time_dimid

     call define_var_netcdf(ncid,ndim,dim,varid,"u_fcst","forecast zonal velocity [u]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"u_anal","analysis zonal velocity at time -dt [u]","m/s")

     call define_var_netcdf(ncid,ndim,dim,varid,"v_fcst","forecast meridional velocity [v]","m/s")
     call define_var_netcdf(ncid,ndim,dim,varid,"v_anal","analysis meridional velocity at time -dt [v]","m/s")

     call define_var_netcdf(ncid,ndim,dim,varid,"t_fcst","forecast potential temperature [el,zz]","degree C")
     call define_var_netcdf(ncid,ndim,dim,varid,"t_anal","analysis potential temperature at time -dt [el,zz]","degree C")

     call define_var_netcdf(ncid,ndim,dim,varid,"s_fcst","forecast salinity [el,zz]","-")
     call define_var_netcdf(ncid,ndim,dim,varid,"s_anal","analysis salinity at time -dt [el,zz]","-")

     deallocate(dim)

     !     end definitions
     status=nf90_enddef(ncid)
     call handle_error_netcdf('nf90_enddef',status,nf90_noerr)
     status=nf90_close(ncid)
     call handle_error_netcdf('nf90_close',status,nf90_noerr)
     
  end if !lpnetcdf .or. master_task
 
  call MPI_Barrier(pom_comm,ierr)
 
  ! open netcdf
  if(lpnetcdf)then
     status=nf90_open(trim(netcdf_out_file),nf90_write,ncid, &
          & comm=mpi_comm_world,info=mpi_info_null)
     call handle_error_netcdf('nf90_open',status,nf90_noerr)     
  else if(my_task == master_task)then
     status=nf90_open(trim(netcdf_out_file),nf90_write,ncid)
     call handle_error_netcdf('nf90_open',status,nf90_noerr)
  end if
        
  ! write data
  !1D
  if(lpnetcdf)then
     call write_pnetcdf_var1d_time_sngl(ncid,1,nprint,"time",real(time))
  else if(my_task == master_task)then
     call write_netcdf_var1d_time_sngl(ncid,1,nprint,"time",real(time))
  end if
  
  !2D
  if(lpnetcdf)then

     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"ua_fcst",real(ua_ave))
     
     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"va_fcst",real(va_ave))
     
     call write_pnetcdf_var2d_time_sngl(ncid,im,jm,nprint,"el_fcst",real(el_ave))
     
  else
     
     call merge2d(glb2d,ua_ave)
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"ua_fcst",real(glb2d))
     
     call merge2d(glb2d,va_ave)
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"va_fcst",real(glb2d))
     
     call merge2d(glb2d,el_ave)
     if(my_task == master_task) &
          & call write_netcdf_var2d_time_sngl(ncid,im_global,jm_global,nprint,"el_fcst",real(glb2d))
          
  end if
  
  !3D
  if(lpnetcdf)then

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"u_fcst",real(u_ave))

     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"v_fcst",real(v_ave))
     
     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"t_fcst",real(t_ave))
     
     call write_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,nprint,"s_fcst",real(s_ave))     
     
  else

     call merge3d(glb3d,u_ave)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"u_fcst",real(glb3d))
     
     call merge3d(glb3d,v_ave)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"v_fcst",real(glb3d))
     
     call merge3d(glb3d,t_ave)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"t_fcst",real(glb3d))
     
     call merge3d(glb3d,s_ave)
     if(my_task == master_task) &
          & call write_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,nprint,"s_fcst",real(glb3d))
     
  end if
     
  ! close file
  if(lpnetcdf .or. my_task == master_task) then
     status=nf90_close(ncid)
     call handle_error_netcdf('nf90_close',status,nf90_noerr)
  end if

  if(.not. lpnetcdf) deallocate(glb2d,glb3d)
     
end subroutine write_iau_netcdf

