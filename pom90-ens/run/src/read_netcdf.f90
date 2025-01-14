!io_netcdf.f90

!NetCDF I/O + MPI

!_______________________________________________________________________
subroutine handle_error_netcdf(routine,status,nf90_noerr)

  use common_pom_var
  implicit none

  character(*),intent(in) :: routine
  integer,intent(in) :: status,nf90_noerr

  if(status /= nf90_noerr) then
     error_status=1
     write(*,'(/a,a)') 'Error: NetCDF routine ',routine
  end if

end subroutine handle_error_netcdf
!_______________________________________________________________________
subroutine read_netcdf_var_sngl(ncid,im,jm,km,varname,glb)

  use netcdf
  use common_pom_var, only: r_sngl
  implicit none

  integer status,varid

  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,km

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(out) :: glb(im,jm,km)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb,(/1,1,1/),(/im,jm,km/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_netcdf_var_sngl
!_______________________________________________________________________
subroutine read_pnetcdf_var_sngl(ncid,im,jm,km,varname,glb)

  use netcdf
  use common_pom_var, only: r_sngl,i_global,j_global
  implicit none

  integer status,varid

  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,km

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(out) :: glb(im,jm,km)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_var_par_access(ncid,varid,nf90_collective)
  call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb,(/i_global(1),j_global(1),1/),(/im,jm,km/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_pnetcdf_var_sngl
!_______________________________________________________________________
subroutine read_netcdf_var_dble(ncid,im,jm,km,varname,glb)

  use netcdf
  use common_pom_var, only: r_dble
  implicit none

  integer status,varid

  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,km

  character(*),intent(in) :: varname

  real(kind = r_dble),intent(out) :: glb(im,jm,km)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb,(/1,1,1/),(/im,jm,km/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_netcdf_var_dble
!_______________________________________________________________________
subroutine read_pnetcdf_var_dble(ncid,im,jm,km,varname,glb)

  use netcdf
  use common_pom_var, only: r_dble,i_global,j_global
  implicit none

  integer status,varid

  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,km

  character(*),intent(in) :: varname

  real(kind = r_dble),intent(out) :: glb(im,jm,km)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_var_par_access(ncid,varid,nf90_collective)
  call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb,(/i_global(1),j_global(1),1/),(/im,jm,km/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_pnetcdf_var_dble
!_______________________________________________________________________
subroutine read_pnetcdf_var0d_sngl(ncid,varname,glb)

  use netcdf
  use common_pom_var, only: r_sngl
  implicit none

  integer status,varid

  integer,intent(in) :: ncid

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(out) :: glb(1)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_var_par_access(ncid,varid,nf90_collective)
  call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb)
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_pnetcdf_var0d_sngl
!_______________________________________________________________________
subroutine read_netcdf_var1d_time_sngl(ncid,im,it,varname,glb)

  use netcdf
  use common_pom_var, only: r_sngl
  implicit none

  integer status,varid

  integer,intent(in) :: ncid
  integer,intent(in) :: im,it

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(out) :: glb(im)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb,(/1,it/),(/im,1/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_netcdf_var1d_time_sngl
!_______________________________________________________________________
subroutine read_pnetcdf_var1d_time_sngl(ncid,im,it,varname,glb)

  use netcdf
  use common_pom_var, only: r_sngl,im_local,jm_local,i_global,j_global
  implicit none

  integer status,varid
  integer is

  integer,intent(in) :: ncid
  integer,intent(in) :: im,it

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(out) :: glb(im)

  if(im == im_local)then
     is=i_global(1)
  else if(im == jm_local)then
     is=j_global(1)
  else
     is=1
  end if
  
  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_var_par_access(ncid,varid,nf90_collective)
  call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb,(/is,it/),(/im,1/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_pnetcdf_var1d_time_sngl
!_______________________________________________________________________
subroutine read_netcdf_var1d_time_dble(ncid,im,it,varname,glb)

  use netcdf
  use common_pom_var, only: r_dble
  implicit none

  integer status,varid

  integer,intent(in) :: ncid
  integer,intent(in) :: im,it

  character(*),intent(in) :: varname

  real(kind = r_dble),intent(out) :: glb(im)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb,(/1,it/),(/im,1/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_netcdf_var1d_time_dble
!_______________________________________________________________________
subroutine read_pnetcdf_var1d_dble(ncid,im,varname,glb)

  use netcdf
  use common_pom_var, only: r_dble
  implicit none

  integer status,varid

  integer,intent(in) :: ncid
  integer,intent(in) :: im

  character(*),intent(in) :: varname

  real(kind = r_dble),intent(out) :: glb(im)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_var_par_access(ncid,varid,nf90_collective)
  call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb,(/1/),(/im/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_pnetcdf_var1d_dble
!_______________________________________________________________________
subroutine read_netcdf_var2d_time_sngl(ncid,im,jm,it,varname,glb)

  use netcdf
  use common_pom_var, only: r_sngl
  implicit none

  integer status,varid

  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,it

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(out) :: glb(im,jm)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb,(/1,1,it/),(/im,jm,1/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_netcdf_var2d_time_sngl
!_______________________________________________________________________
subroutine read_pnetcdf_var2d_time_sngl(ncid,im,jm,it,varname,glb)

  use netcdf
  use common_pom_var, only: r_sngl,i_global,j_global,im_local,jm_local,kb
  implicit none

  integer status,varid
  integer is,js
  
  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,it

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(out) :: glb(im,jm)

  if(im == im_local .and. jm == jm_local)then
     is=i_global(1)
     js=j_global(1)
  else if(im == im_local .and. jm == kb)then
     is=i_global(1)
     js=1
  else if(im == jm_local .and. jm == kb)then
     is=j_global(1)
     js=1
  end if

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_var_par_access(ncid,varid,nf90_collective)
  call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb,(/is,js,it/),(/im,jm,1/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_pnetcdf_var2d_time_sngl
!_______________________________________________________________________
subroutine read_netcdf_var2d_time_dble(ncid,im,jm,it,varname,glb)

  use netcdf
  use common_pom_var, only: r_dble
  implicit none

  integer status,varid

  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,it

  character(*),intent(in) :: varname

  real(kind = r_dble),intent(out) :: glb(im,jm)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb,(/1,1,it/),(/im,jm,1/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_netcdf_var2d_time_dble
!_______________________________________________________________________
subroutine read_netcdf_var3d_time_sngl(ncid,im,jm,km,it,varname,glb)

  use netcdf
  use common_pom_var, only: r_sngl
  implicit none

  integer status,varid

  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,km,it

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(out) :: glb(im,jm,km)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb,(/1,1,1,it/),(/im,jm,km,1/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_netcdf_var3d_time_sngl
!_______________________________________________________________________
subroutine read_pnetcdf_var3d_time_sngl(ncid,im,jm,km,it,varname,glb)

  use netcdf
  use common_pom_var, only: r_sngl,i_global,j_global
  implicit none

  integer status,varid

  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,km,it

  character(*),intent(in) :: varname

  real(kind = r_sngl),intent(out) :: glb(im,jm,km)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_var_par_access(ncid,varid,nf90_collective)
  call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb,(/i_global(1),j_global(1),1,it/),(/im,jm,km,1/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_pnetcdf_var3d_time_sngl
!_______________________________________________________________________
subroutine read_netcdf_var3d_time_dble(ncid,im,jm,km,it,varname,glb)

  use netcdf
  use common_pom_var, only: r_dble
  implicit none

  integer status,varid

  integer,intent(in) :: ncid
  integer,intent(in) :: im,jm,km,it

  character(*),intent(in) :: varname

  real(kind = r_dble),intent(out) :: glb(im,jm,km)

  status=nf90_inq_varid(ncid,varname,varid)
  call handle_error_netcdf('nf90_inq_varid:'//trim(varname),status,nf90_noerr)
  status=nf90_get_var(ncid,varid,glb,(/1,1,1,it/),(/im,jm,km,1/))
  call handle_error_netcdf('nf90_get_var:'//trim(varname),status,nf90_noerr)

end subroutine read_netcdf_var3d_time_dble
!_______________________________________________________________________
subroutine read_netcdf_timeinfo_sngl(ncid,varname,time_julday,dayint)

  use netcdf
  use common_pom_var, only: r_sngl,r_dble,r_size,my_task,master_task,lpnetcdf
  implicit none

  integer status
  integer dimid,varid
  integer ilen                !time length
  integer iyr,imon,iday
  integer julday

  real(kind = r_sngl) time(2) !First two time

  character(30) name          !time name ex. days since yyyy-mm-dd 00:00:00
  character(4) yyyys,yyyye
  character(2) mms,mme,dds,dde

  !IN
  integer,intent(in) :: ncid
  character(*),intent(in) :: varname

  !OUT
  integer,intent(out) :: time_julday 
  real(kind = r_dble),intent(out) :: dayint  !time interval

  !Length
  status=nf90_inq_dimid(ncid,'t',dimid)
  call handle_error_netcdf('nf90_inq_dimid: t',status,nf90_noerr)
  status=nf90_inquire_dimension(ncid,dimid,len = ilen)
  call handle_error_netcdf('nf90_inquire_dimension: t',status,nf90_noerr)

  !Get name
  status=nf90_inq_varid(ncid,trim(varname),varid)
  call handle_error_netcdf('nf90_inq_varid: '//trim(varname),status,nf90_noerr)

  if(lpnetcdf)then
     status=nf90_var_par_access(ncid,varid,nf90_collective)
     call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  end if
  
  status=nf90_get_att(ncid,varid,'units',name)
  call handle_error_netcdf('nf90_get_att: '//trim(varname),status,nf90_noerr)

  !Interval
  status=nf90_get_var(ncid,varid,time,(/1/),(/2/))
  call handle_error_netcdf('nf90_get_var: '//trim(varname),status,nf90_noerr)
  dayint=dble(time(2))-dble(time(1))

  !Start time
  yyyys=name(12:15)
  read(yyyys,*) iyr
  mms=name(17:18)
  read(mms,*) imon
  dds=name(20:21)
  read(dds,*) iday
  time_julday = julday(iyr,imon,iday)

  call caldat(int(time_julday+(ilen-1)*dayint),iyr,imon,iday)

  write(yyyye,'(i4.4)') iyr
  write(mme,'(i2.2)') imon
  write(dde,'(i2.2)') iday

  if(my_task == master_task) then
     write(6,'(a,i10)')   trim(varname)//' time length: ',ilen
     write(6,'(a,f12.5)') trim(varname)//' time interval: ',dayint
     write(6,'(a)')       trim(varname)//' start time: '//yyyys//mms//dds
     write(6,'(a)')       trim(varname)//' end time  : '//yyyye//mme//dde 
  end if

end subroutine read_netcdf_timeinfo_sngl
!_______________________________________________________________________
subroutine read_netcdf_timeinfo_dble(ncid,varname,time_julday,dayint)

  use netcdf
  use common_pom_var, only: r_sngl,r_dble,r_size,my_task,master_task,lpnetcdf
  implicit none

  integer status
  integer dimid,varid
  integer ilen                !time length
  integer iyr,imon,iday
  integer julday

  real(kind = r_dble) time(2) !First two time

  character(30) name          !time name ex. days since yyyy-mm-dd 00:00:00
  character(4) yyyys,yyyye
  character(2) mms,mme,dds,dde

  !IN
  integer,intent(in) :: ncid
  character(*),intent(in) :: varname

  !OUT
  integer,intent(out) :: time_julday 
  real(kind = r_dble),intent(out) :: dayint  !time interval

  !Length
  status=nf90_inq_dimid(ncid,'t',dimid)
  call handle_error_netcdf('nf90_inq_dimid: t',status,nf90_noerr)
  status=nf90_inquire_dimension(ncid,dimid,len = ilen)
  call handle_error_netcdf('nf90_inquire_dimension: t',status,nf90_noerr)

  !Get name
  status=nf90_inq_varid(ncid,trim(varname),varid)
  call handle_error_netcdf('nf90_inq_varid: '//trim(varname),status,nf90_noerr)

  if(lpnetcdf)then
     status=nf90_var_par_access(ncid,varid,nf90_collective)
     call handle_error_netcdf('nf90_var_par_access:'//trim(varname),status,nf90_noerr)
  end if
  
  status=nf90_get_att(ncid,varid,'units',name)
  call handle_error_netcdf('nf90_get_att: '//trim(varname),status,nf90_noerr)

  !Interval
  status=nf90_get_var(ncid,varid,time,(/1/),(/2/))
  call handle_error_netcdf('nf90_get_var: '//trim(varname),status,nf90_noerr)
  dayint=dble(time(2))-dble(time(1))

  !Start time
  yyyys=name(12:15)
  read(yyyys,*) iyr
  mms=name(17:18)
  read(mms,*) imon
  dds=name(20:21)
  read(dds,*) iday
  time_julday = julday(iyr,imon,iday)

  call caldat(int(time_julday+(ilen-1)*dayint),iyr,imon,iday)

  write(yyyye,'(i4.4)') iyr
  write(mme,'(i2.2)') imon
  write(dde,'(i2.2)') iday

  if(my_task == master_task) then
     write(6,'(a,i10)')   trim(varname)//' time length: ',ilen
     write(6,'(a,f12.5)') trim(varname)//' time interval: ',dayint
     write(6,'(a)')       trim(varname)//' start time: '//yyyys//mms//dds
     write(6,'(a)')       trim(varname)//' end time  : '//yyyye//mme//dde 
  end if

end subroutine read_netcdf_timeinfo_dble
!_______________________________________________________________________
subroutine read_netcdf_timeinfo_jra55do(ncid,varname,time_julday,dayint)

  use netcdf
  use common_pom_var, only: r_sngl,r_dble,r_size,my_task,master_task,iint
  implicit none

  integer status
  integer dimid,varid
  integer ilen                !time length
  integer iyr,imon,iday
  integer julday

  real(kind = r_dble) time(2) !First two time

  character(30) name          !time name ex. days since yyyy-mm-dd 00:00:00
  character(4) yyyys,yyyye
  character(2) mms,mme,dds,dde

  !IN
  integer,intent(in) :: ncid
  character(*),intent(in) :: varname

  !OUT
  integer,intent(out) :: time_julday 
  real(kind = r_dble),intent(out) :: dayint  !time interval

  !Length
  status=nf90_inq_dimid(ncid,'time',dimid)
  call handle_error_netcdf('nf90_inq_dimid: time',status,nf90_noerr)
  status=nf90_inquire_dimension(ncid,dimid,len = ilen)
  call handle_error_netcdf('nf90_inquire_dimension: t',status,nf90_noerr)

  !Get name
  status=nf90_inq_varid(ncid,trim(varname),varid)
  call handle_error_netcdf('nf90_inq_varid: '//trim(varname),status,nf90_noerr)
  
  status=nf90_get_att(ncid,varid,'units',name)
  call handle_error_netcdf('nf90_get_att: '//trim(varname),status,nf90_noerr)

  !Interval
  status=nf90_get_var(ncid,varid,time,(/1/),(/2/))
  call handle_error_netcdf('nf90_get_var: '//trim(varname),status,nf90_noerr)
  dayint=dble(time(2))-dble(time(1))

  !Start time
  yyyys=name(12:15)
  read(yyyys,*) iyr
  mms=name(17:18)
  read(mms,*) imon
  dds=name(20:21)
  read(dds,*) iday
  time_julday = julday(iyr,imon,iday)

  call caldat(int(time_julday+time(1)),iyr,imon,iday)
  write(yyyys,'(i4.4)') iyr
  write(mms,'(i2.2)') imon
  write(dds,'(i2.2)') iday

  call caldat(int(time_julday+time(1)+(ilen-1)*dayint),iyr,imon,iday)  
  write(yyyye,'(i4.4)') iyr
  write(mme,'(i2.2)') imon
  write(dde,'(i2.2)') iday

  if(my_task == master_task .and. iint == 1)then
     write(6,'(a,i10)')   trim(varname)//' time length (ATM): ',ilen
     write(6,'(a,f12.5)') trim(varname)//' time interval (ATM): ',dayint
     write(6,'(a)')       trim(varname)//' start time (ATM): '//yyyys//mms//dds
     write(6,'(a)')       trim(varname)//' end time (ATM): '//yyyye//mme//dde 
  end if

end subroutine read_netcdf_timeinfo_jra55do
!_______________________________________________________________________
subroutine read_netcdf_timeinfo_cama(ncid,varname,time_julday,dayint)
  
  use netcdf
  use common_pom_var, only: r_sngl,r_dble,r_size,my_task,master_task,iint
  implicit none

  integer status
  integer dimid,varid
  integer ilen                !time length
  integer iyr,imon,iday
  integer julday

  integer time(2) !First two time

  character(30) name          !time name ex. days since yyyy-mm-dd 00:00:00
  character(4) yyyys,yyyye
  character(2) mms,mme,dds,dde

  !IN
  integer,intent(in) :: ncid
  character(*),intent(in) :: varname

  !OUT
  integer,intent(out) :: time_julday 
  real(kind = r_dble),intent(out) :: dayint  !time interval

  !Length
  status=nf90_inq_dimid(ncid,'time',dimid)
  call handle_error_netcdf('nf90_inq_dimid: time',status,nf90_noerr)
  status=nf90_inquire_dimension(ncid,dimid,len = ilen)
  call handle_error_netcdf('nf90_inquire_dimension: t',status,nf90_noerr)

  !Get name
  status=nf90_inq_varid(ncid,trim(varname),varid)
  call handle_error_netcdf('nf90_inq_varid: '//trim(varname),status,nf90_noerr)
  
  status=nf90_get_att(ncid,varid,'units',name)
  call handle_error_netcdf('nf90_get_att: '//trim(varname),status,nf90_noerr)

  !Interval
  status=nf90_get_var(ncid,varid,time,(/1/),(/2/))
  call handle_error_netcdf('nf90_get_var: '//trim(varname),status,nf90_noerr)
  dayint=3.d0*dble(time(2)-time(1))/24.d0 !*3.d0 due to bug in CaMa-Flood 

  !Start time
  yyyys=name(12:15)
  read(yyyys,*) iyr
  mms=name(17:18)
  read(mms,*) imon
  dds=name(20:21)
  read(dds,*) iday
  time_julday = julday(iyr,imon,iday)

  call caldat(int(time_julday+time(1)),iyr,imon,iday)
  write(yyyys,'(i4.4)') iyr
  write(mms,'(i2.2)') imon
  write(dds,'(i2.2)') iday

  call caldat(int(time_julday+time(1)+(ilen-1)*dayint),iyr,imon,iday)  
  write(yyyye,'(i4.4)') iyr
  write(mme,'(i2.2)') imon
  write(dde,'(i2.2)') iday

  if(my_task == master_task .and. iint == 1)then
     write(6,'(a,i10)')   trim(varname)//' time length (RIVER): ',ilen
     write(6,'(a,f12.5)') trim(varname)//' time interval (RIVER): ',dayint
     write(6,'(a)')       trim(varname)//' start time (RIVER): '//yyyys//mms//dds
     write(6,'(a)')       trim(varname)//' end time (RIVER): '//yyyye//mme//dde 
  end if

end subroutine read_netcdf_timeinfo_cama
!_______________________________________________________________________
subroutine read_grid_netcdf

  ! read grid data

  use mpi
  use netcdf
  use common_pom_var
  implicit none

  integer ncid,status
  real(kind = r_dble) loc2d(im,jm),loc3d(im,jm,kb)
  real(kind = r_dble),allocatable :: glb2d(:,:),glb3d(:,:,:)

  character(120) netcdf_grid_file

  if(.not. lpnetcdf) allocate(glb2d(im_global,jm_global),glb3d(im_global,jm_global,kb))
  
  ! open netcdf file
  write(netcdf_grid_file,'(''in/'',a,''.grid.nc'')') trim(netcdf_file)
  if(my_task == master_task) &
       & write(*,'(/''reading file '',a)') trim(netcdf_grid_file)

  if(lpnetcdf)then
     status=nf90_open(trim(netcdf_grid_file),nf90_nowrite,ncid, &
          & comm=mpi_comm_world,info=mpi_info_null)
  else
     status=nf90_open(trim(netcdf_grid_file),nf90_nowrite,ncid)
  end if
  call handle_error_netcdf('nf90_open: '//trim(netcdf_grid_file),status,nf90_noerr)
     
  ! get variables & data
  if(lpnetcdf)then
     
     call read_pnetcdf_var_dble(ncid,im,jm,1,'dx',loc2d)
     dx(:,:)=loc2d(:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,1,'dy',loc2d)
     dy(:,:)=loc2d(:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,1,'east_u',loc2d)
     east_u(:,:)=loc2d(:,:)     

     call read_pnetcdf_var_dble(ncid,im,jm,1,'east_v',loc2d)
     east_v(:,:)=loc2d(:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,1,'east_e',loc2d)
     east_e(:,:)=loc2d(:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,1,'east_c',loc2d)
     east_c(:,:)=loc2d(:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,1,'north_u',loc2d)
     north_u(:,:)=loc2d(:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,1,'north_v',loc2d)
     north_v(:,:)=loc2d(:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,1,'north_e',loc2d)
     north_e(:,:)=loc2d(:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,1,'north_c',loc2d)
     north_c(:,:)=loc2d(:,:)
     
     call read_pnetcdf_var_dble(ncid,im,jm,1,'h',loc2d)
     h(:,:)=loc2d(:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,1,'fsm',loc2d)
     fsm(:,:)=loc2d(:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,1,'dum',loc2d)
     dum(:,:)=loc2d(:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,1,'dvm',loc2d)
     dvm(:,:)=loc2d(:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,kb,'z_w',loc3d)
     z(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,kb,'z_e',loc3d)
     zz(:,:,:)=loc3d(:,:,:)     
     
  else
     
     call read_netcdf_var_dble(ncid,im_global,jm_global,1,'dx',glb2d)
     call scatter2d(dx,dble(glb2d))
     
     call read_netcdf_var_dble(ncid,im_global,jm_global,1,'dy',glb2d)
     call scatter2d(dy,dble(glb2d))

     call read_netcdf_var_dble(ncid,im_global,jm_global,1,'east_u',glb2d)
     call scatter2d(east_u,dble(glb2d))

     call read_netcdf_var_dble(ncid,im_global,jm_global,1,'east_v',glb2d)
     call scatter2d(east_v,dble(glb2d))
     
     call read_netcdf_var_dble(ncid,im_global,jm_global,1,'east_e',glb2d)  
     call scatter2d(east_e,dble(glb2d))

     call read_netcdf_var_dble(ncid,im_global,jm_global,1,'east_c',glb2d)  
     call scatter2d(east_c,dble(glb2d))

     call read_netcdf_var_dble(ncid,im_global,jm_global,1,'north_u',glb2d)
     call scatter2d(north_u,dble(glb2d))

     call read_netcdf_var_dble(ncid,im_global,jm_global,1,'north_v',glb2d)
     call scatter2d(north_v,dble(glb2d))

     call read_netcdf_var_dble(ncid,im_global,jm_global,1,'north_e',glb2d)
     call scatter2d(north_e,dble(glb2d))

     call read_netcdf_var_dble(ncid,im_global,jm_global,1,'north_c',glb2d)
     call scatter2d(north_c,dble(glb2d))
     
     call read_netcdf_var_dble(ncid,im_global,jm_global,1,'h',glb2d)
     call scatter2d(h,dble(glb2d))

     call read_netcdf_var_dble(ncid,im_global,jm_global,1,'fsm',glb2d)
     call scatter2d(fsm,dble(glb2d))

     call read_netcdf_var_dble(ncid,im_global,jm_global,1,'dum',glb2d)
     call scatter2d(dum,dble(glb2d))
     
     call read_netcdf_var_dble(ncid,im_global,jm_global,1,'dvm',glb2d)
     call scatter2d(dvm,dble(glb2d))

     call read_netcdf_var_dble(ncid,im_global,jm_global,kb,'z_w',glb3d)
     call scatter3d(z,dble(glb3d))

     call read_netcdf_var_dble(ncid,im_global,jm_global,kb,'z_e',glb3d)
     call scatter3d(zz,dble(glb3d))     
     
  end if
  
  ! close file:
  status=nf90_close(ncid)
  call handle_error_netcdf('nf90_close: grid',status,nf90_noerr)

  if(.not. lpnetcdf) deallocate(glb2d,glb3d)
  
end subroutine read_grid_netcdf
!_______________________________________________________________________
subroutine read_restart_netcdf
  ! read data for a seamless restart

  !$use omp_lib  
  use mpi
  use netcdf
  use common_pom_var  
  implicit none

  integer i,j
  integer ncid,status

  !sngl
  real(kind = r_sngl) loc2d(im,jm),loc3d(im,jm,kb)
  real(kind = r_sngl) glb1d
  real(kind = r_sngl),allocatable :: glb2d(:,:),glb3d(:,:,:)

  character(120) netcdf_in_file

  if(.not. lpnetcdf) &
       & allocate(glb2d(im_global,jm_global),glb3d(im_global,jm_global,kb))
  
  ! open netcdf restart file
  write(netcdf_in_file,'(''in/'',a)') trim(read_rst_file)
  if(my_task == master_task) &
       & write(*,'(/''reading file '',a)')trim(netcdf_in_file)

  if(lpnetcdf)then
     status=nf90_open(netcdf_in_file,nf90_nowrite,ncid, &
          & comm=mpi_comm_world,info=mpi_info_null)
  else
     status=nf90_open(netcdf_in_file,nf90_nowrite,ncid)
  end if
  call handle_error_netcdf('nf90_open: '//trim(read_rst_file),status,nf90_noerr)     
     
  !---get variables & data
  !---1d
  !time0
  if(lpnetcdf)then
     call read_pnetcdf_var0d_sngl(ncid,'time',glb1d)
  else
     call read_netcdf_var_sngl(ncid,1,1,1,'time',glb1d)
  end if
  time0=dble(glb1d)

  !---2d
  if(lpnetcdf)then

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'wubot',loc2d)
     wubot(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'wvbot',loc2d)
     wvbot(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'aam2d',loc2d)
     aam2d(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'ua',loc2d)
     ua(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'uab',loc2d)
     uab(:,:)=loc2d(:,:)    

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'va',loc2d)
     va(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'vab',loc2d)
     vab(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'el',loc2d)
     el(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'elb',loc2d)
     elb(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'et',loc2d)
     et(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'etb',loc2d)
     etb(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'egb',loc2d)
     egb(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'utb',loc2d)
     utb(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'vtb',loc2d)
     vtb(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'adx2d',loc2d)
     adx2d(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'ady2d',loc2d)
     ady2d(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'advua',loc2d)
     advua(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'advva',loc2d)
     advva(:,:)=loc2d(:,:)

  else
     
     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'wubot',glb2d)
     call scatter2d(wubot,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'wvbot',glb2d)
     call scatter2d(wvbot,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'aam2d',glb2d)
     call scatter2d(aam2d,dble(glb2d))
     
     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'ua',glb2d)
     call scatter2d(ua,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'uab',glb2d)
     call scatter2d(uab,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'va',glb2d)
     call scatter2d(va,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'vab',glb2d)
     call scatter2d(vab,dble(glb2d))
     
     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'el',glb2d)
     call scatter2d(el,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'elb',glb2d)
     call scatter2d(elb,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'et',glb2d)
     call scatter2d(et,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'etb',glb2d)
     call scatter2d(etb,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'egb',glb2d)
     call scatter2d(egb,dble(glb2d))
     
     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'utb',glb2d)
     call scatter2d(utb,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'vtb',glb2d)
     call scatter2d(vtb,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'adx2d',glb2d)
     call scatter2d(adx2d,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'ady2d',glb2d)
     call scatter2d(ady2d,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'advua',glb2d)
     call scatter2d(advua,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'advva',glb2d)
     call scatter2d(advva,dble(glb2d))
          
  end if
                                    
  !---3d
  if(lpnetcdf)then
     
     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'u',loc3d)
     u(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'ub',loc3d)
     ub(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'v',loc3d)
     v(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'vb',loc3d)
     vb(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'w',loc3d)
     w(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'t',loc3d)
     t(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'tb',loc3d)
     tb(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'s',loc3d)
     s(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'sb',loc3d)
     sb(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'rho',loc3d)
     rho(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'km',loc3d)
     km(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'kh',loc3d)
     kh(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'kq',loc3d)
     kq(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'l',loc3d)
     l(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'q2',loc3d)
     q2(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'q2b',loc3d)
     q2b(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'q2l',loc3d)
     q2l(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'q2lb',loc3d)
     q2lb(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'aam',loc3d)
     aam(:,:,:)=loc3d(:,:,:)
     
  else
     
     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'u',glb3d)
     call scatter3d(u,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'ub',glb3d)
     call scatter3d(ub,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'v',glb3d)
     call scatter3d(v,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'vb',glb3d)
     call scatter3d(vb,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'w',glb3d)
     call scatter3d(w,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'t',glb3d)
     call scatter3d(t,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'tb',glb3d)
     call scatter3d(tb,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'s',glb3d)
     call scatter3d(s,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'sb',glb3d)
     call scatter3d(sb,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'rho',glb3d)
     call scatter3d(rho,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'km',glb3d)
     call scatter3d(km,dble(glb3d))
     
     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'kh',glb3d)
     call scatter3d(kh,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'kq',glb3d)
     call scatter3d(kq,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'l',glb3d)
     call scatter3d(l,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'q2',glb3d)
     call scatter3d(q2,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'q2b',glb3d)
     call scatter3d(q2b,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'q2l',glb3d)
     call scatter3d(q2l,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'q2lb',glb3d)
     call scatter3d(q2lb,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'aam',glb3d)
     call scatter3d(aam,dble(glb3d))
     
  end if
       
  ! close file
  status=nf90_close(ncid)
  call handle_error_netcdf('nf90_close',status,nf90_noerr)

  ! Intermittent method
  if(assim == 3)then

     ua(1:im,1:jm)=uab(1:im,1:jm)
     va(1:im,1:jm)=vab(1:im,1:jm)
     el(1:im,1:jm)=elb(1:im,1:jm)

     u(1:im,1:jm,1:kb)=ub(1:im,1:jm,1:kb)
     v(1:im,1:jm,1:kb)=vb(1:im,1:jm,1:kb)
     t(1:im,1:jm,1:kb)=tb(1:im,1:jm,1:kb)
     s(1:im,1:jm,1:kb)=sb(1:im,1:jm,1:kb)

  end if

  ! initialize elevation
  !$omp parallel
  !$omp do private(i,j)
  do j=1,jm
     do i=1,im
        d(i,j)=h(i,j)+el(i,j)
        dt(i,j)=h(i,j)+et(i,j)
     end do
  end do
  !$omp end do
  !$omp end parallel

  ! initialize time
  time=time0

  if(.not. lpnetcdf) deallocate(glb2d,glb3d)
  
end subroutine read_restart_netcdf
!_______________________________________________________________________
subroutine read_iau_netcdf
  ! read data for a seamless restart

  !$use omp_lib  
  use mpi
  use netcdf
  use common_pom_var
  implicit none

  integer i,j,k
  integer ncid,status

  character(120) netcdf_in_file

  !r_sngl  
  real(kind = r_sngl) loc2d(im,jm),loc3d(im,jm,kb)
  real(kind = r_sngl),allocatable :: glb2d(:,:),glb3d(:,:,:)

  !r_dble
  real(kind = r_dble) ua_fcst(im,jm),ua_anal(im,jm)
  real(kind = r_dble) va_fcst(im,jm),va_anal(im,jm)
  real(kind = r_dble) el_fcst(im,jm),el_anal(im,jm)
  real(kind = r_dble) t_fcst(im,jm,kb),t_anal(im,jm,kb)
  real(kind = r_dble) s_fcst(im,jm,kb),s_anal(im,jm,kb)
  real(kind = r_dble) u_fcst(im,jm,kb),u_anal(im,jm,kb)
  real(kind = r_dble) v_fcst(im,jm,kb),v_anal(im,jm,kb)

  if(.not. lpnetcdf) &
       & allocate(glb2d(im_global,jm_global),glb3d(im_global,jm_global,kb))
  
  ! open netcdf restart file
  write(netcdf_in_file,'(''in/'',a)') trim(read_iau_file)
  if(my_task == master_task) &
       & write(*,'(/''reading file '',a)') trim(netcdf_in_file)

  if(lpnetcdf)then
     status=nf90_open(netcdf_in_file,nf90_nowrite,ncid, &
          & comm=mpi_comm_world,info=mpi_info_null)
  else
     status=nf90_open(netcdf_in_file,nf90_nowrite,ncid)
  end if
  call handle_error_netcdf('nf90_open: '//trim(read_iau_file),status,nf90_noerr)
     
  ! get variables & data
  !2d
  if(lpnetcdf)then

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'ua_fcst',loc2d)
     ua_fcst(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'ua_anal',loc2d)
     ua_anal(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'va_fcst',loc2d)
     va_fcst(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'va_anal',loc2d)
     va_anal(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'el_fcst',loc2d)
     el_fcst(:,:)=loc2d(:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,1,'el_anal',loc2d)
     el_anal(:,:)=loc2d(:,:)
     
  else
     
     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'ua_fcst',glb2d)
     call scatter2d(ua_fcst,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'ua_anal',glb2d)
     call scatter2d(ua_anal,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'va_fcst',glb2d)
     call scatter2d(va_fcst,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'va_anal',glb2d)
     call scatter2d(va_anal,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'el_fcst',glb2d)
     call scatter2d(el_fcst,dble(glb2d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,1,'el_anal',glb2d)
     call scatter2d(el_anal,dble(glb2d))
     
  end if

  !3d
  if(lpnetcdf)then
     
     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'u_fcst',loc3d)
     u_fcst(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'u_anal',loc3d)
     u_anal(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'v_fcst',loc3d)
     v_fcst(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'v_anal',loc3d)
     v_anal(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'t_fcst',loc3d)
     t_fcst(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'t_anal',loc3d)
     t_anal(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'s_fcst',loc3d)
     s_fcst(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'s_anal',loc3d)
     s_anal(:,:,:)=loc3d(:,:,:)
     
  else

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'u_fcst',glb3d)
     call scatter3d(u_fcst,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'u_anal',glb3d)
     call scatter3d(u_anal,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'v_fcst',glb3d)
     call scatter3d(v_fcst,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'v_anal',glb3d)
     call scatter3d(v_anal,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'t_fcst',glb3d)
     call scatter3d(t_fcst,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'t_anal',glb3d)
     call scatter3d(t_anal,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'s_fcst',glb3d)
     call scatter3d(s_fcst,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'s_anal',glb3d)
     call scatter3d(s_anal,dble(glb3d))     
     
  end if
  
  ! close file
  status=nf90_close(ncid)
  call handle_error_netcdf('nf90_close',status,nf90_noerr)

  ! Calculate variable_iau (analysis-forecast)

  !$omp parallel
  !$omp do private(i,j)
  do j=1,jm
     do i=1,im
        ua_iau(i,j)=ua_anal(i,j)-ua_fcst(i,j)
        va_iau(i,j)=va_anal(i,j)-va_fcst(i,j)
        el_iau(i,j)=el_anal(i,j)-el_fcst(i,j)
     end do
  end do
  !$omp end do

  !$omp do private(i,j,k)  
  do k=1,kb
     do j=1,jm
        do i=1,im
           t_iau(i,j,k)=t_anal(i,j,k)-t_fcst(i,j,k)
           s_iau(i,j,k)=s_anal(i,j,k)-s_fcst(i,j,k)
           u_iau(i,j,k)=u_anal(i,j,k)-u_fcst(i,j,k)
           v_iau(i,j,k)=v_anal(i,j,k)-v_fcst(i,j,k)
        end do
     end do
  end do
  !$omp end do
  !$omp end parallel

  if(.not. lpnetcdf) deallocate(glb2d,glb3d)
  
end subroutine read_iau_netcdf
!_______________________________________________________________________
subroutine read_initial_tsuv_netcdf(temp,salt,xvel,yvel)
  ! read initial temperature, salinity, x-velocity, and y-velocity

  use mpi
  use netcdf
  use common_pom_var
  implicit none

  integer ncid,status  
  real(kind = r_dble) loc3d(im,jm,kb)  
  real(kind = r_dble),allocatable :: glb3d(:,:,:)  
  character(120) netcdf_ic_file

  real(kind = r_size),intent(out) :: temp(im,jm,kb),salt(im,jm,kb)
  real(kind = r_size),intent(out) :: xvel(im,jm,kb),yvel(im,jm,kb)

  if(.not. lpnetcdf) allocate(glb3d(im_global,jm_global,kb))
  
  ! open netcdf ic file
  write(netcdf_ic_file,'(''in/'',a,''.ic.nc'')') trim(netcdf_file)
  if(my_task == master_task) &
       & write(*,'(/''reading file '',a)') trim(netcdf_ic_file)

  if(lpnetcdf)then
     status=nf90_open(trim(netcdf_ic_file),nf90_nowrite,ncid, &
          & comm=mpi_comm_world,info=mpi_info_null)
  else
     status=nf90_open(trim(netcdf_ic_file),nf90_nowrite,ncid)
  end if
  call handle_error_netcdf('nf90_open: '//netcdf_ic_file,status,nf90_noerr)
     
  ! get variables & data
  if(lpnetcdf)then

     call read_pnetcdf_var_dble(ncid,im,jm,kb,'t',loc3d)
     temp(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,kb,'s',loc3d)
     salt(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,kb,'u',loc3d)
     xvel(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_dble(ncid,im,jm,kb,'v',loc3d)
     yvel(:,:,:)=loc3d(:,:,:)     
     
  else
     
     call read_netcdf_var_dble(ncid,im_global,jm_global,kb,'t',glb3d)
     call scatter3d(temp,dble(glb3d))

     call read_netcdf_var_dble(ncid,im_global,jm_global,kb,'s',glb3d)
     call scatter3d(salt,dble(glb3d))

     call read_netcdf_var_dble(ncid,im_global,jm_global,kb,'u',glb3d)
     call scatter3d(xvel,dble(glb3d))

     call read_netcdf_var_dble(ncid,im_global,jm_global,kb,'v',glb3d)  
     call scatter3d(yvel,dble(glb3d))
     
  end if
     
  ! close file
  status=nf90_close(ncid)
  call handle_error_netcdf('nf90_close',status,nf90_noerr)

  if(.not. lpnetcdf) deallocate(glb3d)
  
end subroutine read_initial_tsuv_netcdf
!_______________________________________________________________________
subroutine read_tsclim_netcdf

  ! read temperature and salinity climatological/mean data

  use mpi
  use netcdf
  use common_pom_var
  implicit none

  integer ncid,status

  !single precision
  real(kind = r_sngl) :: loc3d(im,jm,kb) 
  real(kind = r_sngl),allocatable :: glb3d(:,:,:)

  character(120) netcdf_tsclim_file

  if(.not. lpnetcdf) allocate(glb3d(im_global,jm_global,kb))
  
  ! open netcdf file
  write(netcdf_tsclim_file,'(''in/'',a,''.tsclim.nc'')') trim(netcdf_file)
  if(my_task == master_task) &
       & write(*,'(/''reading file '',a)') trim(netcdf_tsclim_file)

  if(lpnetcdf)then
     status=nf90_open(netcdf_tsclim_file,nf90_nowrite,ncid, &
          & comm=mpi_comm_world,info=mpi_info_null)
  else
     status=nf90_open(netcdf_tsclim_file,nf90_nowrite,ncid)
  end if
  call handle_error_netcdf('nf90_open: '//netcdf_tsclim_file,status,nf90_noerr)
     
  ! read data
  if(lpnetcdf)then
     
     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'tclim',loc3d)
     tclim(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'sclim',loc3d)
     sclim(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var_sngl(ncid,im,jm,kb,'rmean',loc3d)
     rmean(:,:,:)=loc3d(:,:,:)
     
  else
     
     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'tclim',glb3d)  
     call scatter3d(tclim,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'sclim',glb3d)  
     call scatter3d(sclim,dble(glb3d))

     call read_netcdf_var_sngl(ncid,im_global,jm_global,kb,'rmean',glb3d)  
     call scatter3d(rmean,dble(glb3d))
     
  end if
     
  ! close file
  status=nf90_close(ncid)
  call handle_error_netcdf('nf90_close',status,nf90_noerr)

  if(.not. lpnetcdf) deallocate(glb3d)
  
end subroutine read_tsclim_netcdf
!______________________________________________________________________________
subroutine read_atm_netcdf

  ! read atmospheric condition 

  use mpi
  use netcdf
  use common_pom_var
  implicit none

  !---Parameter
  integer,parameter :: nvar=8

  !---Common
  integer ncid,status
  integer it
  integer i,j,k
  integer ivar
  integer iyr,imon,iday,ihour
  integer iyr_ens(nens)
  integer ijul
  integer ijul_ens
  
  real(kind = r_sngl) loc2d_sngl(im,jm)
  real(kind = r_size) loc2d_size(im,jm)
  real(kind = r_sngl),allocatable :: glb2d(:,:)
  
  character(120) netcdf_atm_file
  character(5) var(nvar)
  character(4) yyyy
  character(2) mm,dd,hh
  
  if(.not. lpnetcdf) allocate(glb2d(im_global,jm_global))

  !---variable name
  var(1)="windu"
  var(2)="windv"
  var(3)="airt"
  var(4)="airh"
  var(5)="swrad"
  var(6)="slp"
  var(7)="lwrad"
  var(8)="prep"
  
  do k=1,2

     !---iyr,imon,iday,ihour
     if(k == 1)then
        call julian_ymd(int(julday_start+time),iyr,imon,iday)
        ihour=(julday_start+time-int(julday_start+time))*24.d0
     else if(k == 2)then
        call julian_ymd(int(julday_start+time+atmtime_dayint),iyr,imon,iday)
        ihour=(julday_start+time+atmtime_dayint-int(julday_start+time+atmtime_dayint))*24.d0
     end if
     
     write(yyyy,'(i4.4)') iyr
     write(mm,'(i2.2)') imon
     write(dd,'(i2.2)') iday
     write(hh,'(i2.2)') ihour

     netcdf_atm_file="in/atm/atm."//yyyy//mm//dd//".nc"
     if(my_task == master_task) write(*,'(/''reading file '',a)') trim(netcdf_atm_file)

     if(lpnetcdf)then
        status=nf90_open(trim(netcdf_atm_file),nf90_nowrite,ncid, &
             & comm=mpi_comm_world,info=mpi_info_null)
     else
        status=nf90_open(trim(netcdf_atm_file),nf90_nowrite,ncid)
     end if
     call handle_error_netcdf('nf90_open: '//trim(netcdf_atm_file),status,nf90_noerr)

     if(iint == 1 .and. k == 1)then
        call read_netcdf_timeinfo_sngl(ncid,'atmtime',atmtime_julday,atmtime_dayint)
     end if
        
     !---it
     it=nint(ihour/(atmtime_dayint*24.d0))+1
     
     if(my_task == master_task)then
        write(6,'(a,i10)') 'atm step: ',it
        write(6,'(a)') 'atm date: '//yyyy//mm//dd//hh       
     end if
     
     if(it <= 0) then
        if(my_task == master_task) &
             & write(*,'(/a)') 'wrong time of atmospheric data'
        call finalize_mpi
        stop
     end if

     !---Read data
     do ivar=1,nvar
     
        ! get variables
        if(lpnetcdf)then
           call read_pnetcdf_var2d_time_sngl(ncid,im,jm,it,trim(var(ivar)),loc2d_sngl)
           loc2d_size(:,:)=REAL(loc2d_sngl(:,:),r_size)
        else
           call read_netcdf_var2d_time_sngl(ncid,im_global,jm_global,it,trim(var(ivar)),glb2d)
           call scatter2d(loc2d_size,dble(glb2d))
        end if

        if(k == 1 .and. var(ivar) == "windu")then
           windu0(:,:)=loc2d_size(:,:)
        else if(k == 1 .and. var(ivar) == "windv")then
           windv0(:,:)=loc2d_size(:,:)           
        else if(k == 1 .and. var(ivar) == "airt")then
           airt0(:,:)=loc2d_size(:,:)
        else if(k == 1 .and. var(ivar) == "airh")then
           airh0(:,:)=loc2d_size(:,:)
        else if(k == 1 .and. var(ivar) == "swrad")then
           swrad0(:,:)=loc2d_size(:,:)
        else if(k == 1 .and. var(ivar) == "slp")then
           slp0(:,:)=loc2d_size(:,:)
        else if(k == 1 .and. var(ivar) == "lwrad")then
           lwrad0(:,:)=loc2d_size(:,:)
        else if(k == 1 .and. var(ivar) == "prep")then
           prep0(:,:)=loc2d_size(:,:)
        else if(k == 2 .and. var(ivar) == "windu")then
           windu1(:,:)=loc2d_size(:,:)
        else if(k == 2 .and. var(ivar) == "windv")then
           windv1(:,:)=loc2d_size(:,:)           
        else if(k == 2 .and. var(ivar) == "airt")then
           airt1(:,:)=loc2d_size(:,:)
        else if(k == 2 .and. var(ivar) == "airh")then
           airh1(:,:)=loc2d_size(:,:)
        else if(k == 2 .and. var(ivar) == "swrad")then
           swrad1(:,:)=loc2d_size(:,:)
        else if(k == 2 .and. var(ivar) == "slp")then
           slp1(:,:)=loc2d_size(:,:)
        else if(k == 2 .and. var(ivar) == "lwrad")then
           lwrad1(:,:)=loc2d_size(:,:)
        else if(k == 2 .and. var(ivar) == "prep")then
           prep1(:,:)=loc2d_size(:,:)
        end if
           
     end do !ivar

     ! close file
     status=nf90_close(ncid)
     call handle_error_netcdf('nf90_close',status,nf90_noerr)

  end do !k

  if(.not. lpnetcdf) deallocate(glb2d)

  if(iens == 0) return

  !---Perturbed atmospheric forcing-------------------------------------
  
  do k=1,2

     !---iyr,imon,iday,ihour
     if(k == 1)then
        call julian_ymd(int(julday_start+time),iyr,imon,iday)
        ihour=(julday_start+time-int(julday_start+time))*24.d0
     else if(k == 2)then
        call julian_ymd(int(julday_start+time+atmtime_dayint),iyr,imon,iday)
        ihour=(julday_start+time+atmtime_dayint-int(julday_start+time+atmtime_dayint))*24.d0
     end if

     call set_ensemble_year(syr_atm,eyr_atm,iyr,imon,nens,iyr_ens)     

     write(yyyy,'(i4.4)') iyr_ens(iens)
     write(mm,'(i2.2)') imon
     write(dd,'(i2.2)') iday
     write(hh,'(i2.2)') ihour

     netcdf_atm_file="in/atm/atm."//yyyy//mm//dd//".nc"
     if(my_task == master_task) write(*,'(/''reading file '',a)') trim(netcdf_atm_file)

     if(lpnetcdf)then
        status=nf90_open(trim(netcdf_atm_file),nf90_nowrite,ncid, &
             & comm=mpi_comm_world,info=mpi_info_null)
     else
        status=nf90_open(trim(netcdf_atm_file),nf90_nowrite,ncid)
     end if
     call handle_error_netcdf('nf90_open: '//trim(netcdf_atm_file),status,nf90_noerr)
          
     !---it
     it=nint(ihour/(atmtime_dayint*24.d0))+1
     
     if(my_task == master_task)then
        write(6,'(a,i10)') 'atm step (ens): ',it
        write(6,'(a)') 'atm date (ens): '//yyyy//mm//dd//hh       
     end if
     
     if(it <= 0) then
        if(my_task == master_task) &
             & write(*,'(/a)') 'wrong time of atmospheric data'
        call finalize_mpi
        stop
     end if
     
     !---Read data
     do ivar=1,nvar-1 !except for precipitation
     
        ! get variables
        if(lpnetcdf)then
           call read_pnetcdf_var2d_time_sngl(ncid,im,jm,it,trim(var(ivar)),loc2d_sngl)
           loc2d_size(:,:)=REAL(loc2d_sngl(:,:),r_size)
        else
           call read_netcdf_var2d_time_sngl(ncid,im_global,jm_global,it,trim(var(ivar)),glb2d)
           call scatter2d(loc2d_size,dble(glb2d))
        end if

        !$omp parallel
        !$omp do private(i,j)
        do j=1,jm
           do i=1,im

              !x=x+
              if(k == 1 .and. var(ivar) == "windu")then
                 windu0(i,j)=windu0(i,j)+alpha_atm*loc2d_size(i,j)
              else if(k == 1 .and. var(ivar) == "windv")then
                 windv0(i,j)=windv0(i,j)+alpha_atm*loc2d_size(i,j)
              else if(k == 1 .and. var(ivar) == "airt")then
                 airt0(i,j)=airt0(i,j)+alpha_atm*loc2d_size(i,j)
              else if(k == 1 .and. var(ivar) == "airh")then
                 airh0(i,j)=airh0(i,j)+alpha_atm*loc2d_size(i,j)
              else if(k == 1 .and. var(ivar) == "swrad")then
                 swrad0(i,j)=swrad0(i,j)+alpha_atm*loc2d_size(i,j)
              else if(k == 1 .and. var(ivar) == "slp")then
                 slp0(i,j)=slp0(i,j)+alpha_atm*loc2d_size(i,j)
              else if(k == 1 .and. var(ivar) == "lwrad")then
                 lwrad0(i,j)=lwrad0(i,j)+alpha_atm*loc2d_size(i,j)
              else if(k == 2 .and. var(ivar) == "windu")then
                 windu1(i,j)=windu1(i,j)+alpha_atm*loc2d_size(i,j)
              else if(k == 2 .and. var(ivar) == "windv")then
                 windv1(i,j)=windv1(i,j)+alpha_atm*loc2d_size(i,j)
              else if(k == 2 .and. var(ivar) == "airt")then
                 airt1(i,j)=airt1(i,j)+alpha_atm*loc2d_size(i,j)
              else if(k == 2 .and. var(ivar) == "airh")then
                 airh1(i,j)=airh1(i,j)+alpha_atm*loc2d_size(i,j)
              else if(k == 2 .and. var(ivar) == "swrad")then
                 swrad1(i,j)=swrad1(i,j)+alpha_atm*loc2d_size(i,j)
              else if(k == 2 .and. var(ivar) == "slp")then
                 slp1(i,j)=slp1(i,j)+alpha_atm*loc2d_size(i,j)
              else if(k == 2 .and. var(ivar) == "lwrad")then
                 lwrad1(i,j)=lwrad1(i,j)+alpha_atm*loc2d_size(i,j)
              end if
              
           end do
        end do
        !$omp end do
        !$omp end parallel

     end do !ivar

     ! close file
     status=nf90_close(ncid)
     call handle_error_netcdf('nf90_close',status,nf90_noerr)
     
  end do !k  

  !---open netcdf atm_clim file
  write(netcdf_atm_file,'(''in/'',a,''.atm_clim.nc'')') trim(netcdf_file)
  if(my_task == master_task) write(*,'(/''reading file '',a)') trim(netcdf_atm_file)

  if(lpnetcdf)then
     status=nf90_open(trim(netcdf_atm_file),nf90_nowrite,ncid, &
          & comm=mpi_comm_world,info=mpi_info_null)
  else
     status=nf90_open(trim(netcdf_atm_file),nf90_nowrite,ncid)
  end if
  call handle_error_netcdf('nf90_open: '//trim(netcdf_atm_file),status,nf90_noerr)     
  
  do k=1,2

     !---imon,iday,ihour
     if(k == 1)then
        call julian_ymd(int(julday_start+time),iyr,imon,iday)
        ihour=(julday_start+time-int(julday_start+time))*24.d0
     else if(k == 2)then
        call julian_ymd(int(julday_start+time+atmtime_dayint),iyr,imon,iday)
        ihour=(julday_start+time+atmtime_dayint-int(julday_start+time+atmtime_dayint))*24.d0
     end if

     !Leap year
     if(imon == 2 .and. iday == 29) iday=28
     
     write(mm,'(i2.2)') imon
     write(dd,'(i2.2)') iday
     write(hh,'(i2.2)') ihour

     call ymd_julian(1,imon,iday,ijul_ens)
     call ymd_julian(1,1,1,ijul)
     
     !---it
     it=(ijul_ens-ijul)/atmtime_dayint+nint(ihour/(atmtime_dayint*24.d0))+1

     if(my_task == master_task)then
        write(6,'(a,i10)') 'atm clim step: ',it
        write(6,'(a)') 'atm clime date: '//mm//dd//hh
     end if

     if(it <= 0) then
        if(my_task == master_task) &
             & write(*,'(/a)') 'wrong time of atmospheric data'
        call finalize_mpi
        stop
     end if
     
     !---Read data
     do ivar=1,nvar-1 !except for precipitation
     
        ! get variables
        if(lpnetcdf)then
           call read_pnetcdf_var2d_time_sngl(ncid,im,jm,it,trim(var(ivar)),loc2d_sngl)
           loc2d_size(:,:)=REAL(loc2d_sngl(:,:),r_size)
        else
           call read_netcdf_var2d_time_sngl(ncid,im_global,jm_global,it,trim(var(ivar)),glb2d)
           call scatter2d(loc2d_size,dble(glb2d))
        end if

        !$omp parallel
        !$omp do private(i,j)
        do j=1,jm
           do i=1,im

              !x=x+
              if(k == 1 .and. var(ivar) == "windu")then
                 windu0(i,j)=windu0(i,j)-alpha_atm*loc2d_size(i,j)
              else if(k == 1 .and. var(ivar) == "windv")then
                 windv0(i,j)=windv0(i,j)-alpha_atm*loc2d_size(i,j)
              else if(k == 1 .and. var(ivar) == "airt")then
                 airt0(i,j)=airt0(i,j)-alpha_atm*loc2d_size(i,j)
              else if(k == 1 .and. var(ivar) == "airh")then
                 airh0(i,j)=airh0(i,j)-alpha_atm*loc2d_size(i,j)
              else if(k == 1 .and. var(ivar) == "swrad")then
                 swrad0(i,j)=swrad0(i,j)-alpha_atm*loc2d_size(i,j)
              else if(k == 1 .and. var(ivar) == "slp")then
                 slp0(i,j)=slp0(i,j)-alpha_atm*loc2d_size(i,j)
              else if(k == 1 .and. var(ivar) == "lwrad")then
                 lwrad0(i,j)=lwrad0(i,j)-alpha_atm*loc2d_size(i,j)
              else if(k == 2 .and. var(ivar) == "windu")then
                 windu1(i,j)=windu1(i,j)-alpha_atm*loc2d_size(i,j)
              else if(k == 2 .and. var(ivar) == "windv")then
                 windv1(i,j)=windv1(i,j)-alpha_atm*loc2d_size(i,j)
              else if(k == 2 .and. var(ivar) == "airt")then
                 airt1(i,j)=airt1(i,j)-alpha_atm*loc2d_size(i,j)
              else if(k == 2 .and. var(ivar) == "airh")then
                 airh1(i,j)=airh1(i,j)-alpha_atm*loc2d_size(i,j)
              else if(k == 2 .and. var(ivar) == "swrad")then
                 swrad1(i,j)=swrad1(i,j)-alpha_atm*loc2d_size(i,j)
              else if(k == 2 .and. var(ivar) == "slp")then
                 slp1(i,j)=slp1(i,j)-alpha_atm*loc2d_size(i,j)
              else if(k == 2 .and. var(ivar) == "lwrad")then
                 lwrad1(i,j)=lwrad1(i,j)-alpha_atm*loc2d_size(i,j)
              end if
              
           end do
        end do
        !$omp end do
        !$omp end parallel

     end do !ivar
  end do !k  

  ! close file
  status=nf90_close(ncid)
  call handle_error_netcdf('nf90_close',status,nf90_noerr)
  
end subroutine read_atm_netcdf
!_______________________________________________________________________
subroutine read_riv_netcdf

  ! read river condition
  !     2018.06.06 S.Ohishi
  !     2018.08.20 S.Ohishi

  use mpi
  use netcdf
  use common_pom_var  
  implicit none

  integer ncid,status
  integer k
  integer iyr,imon,iday,ihour
  integer it

  real(kind = r_sngl) loc2d_sngl(im,jm)
  real(kind = r_size) loc2d_size(im,jm)
  real(kind = r_sngl),allocatable :: glb2d(:,:)

  character(120) netcdf_riv_file
  character(4) yyyy
  character(2) mm,dd,hh
  
  if(.not. lpnetcdf) allocate(glb2d(im_global,jm_global))
  
  do k=1,2
  
     !---iyr,imon,iday,ihour
     if(k == 1)then
        call julian_ymd(int(julday_start+time),iyr,imon,iday)
        ihour=(julday_start+time-int(julday_start+time))*24.d0
     else if(k == 2)then
        call julian_ymd(int(julday_start+time+rivtime_dayint),iyr,imon,iday)
        ihour=(julday_start+time+rivtime_dayint-int(julday_start+time+rivtime_dayint))*24.d0
     end if
     
     write(yyyy,'(i4.4)') iyr
     write(mm,'(i2.2)') imon
     write(dd,'(i2.2)') iday
     write(hh,'(i2.2)') ihour

     ! open netcdf fresh water flux file
     netcdf_riv_file="in/river/river."//yyyy//mm//dd//".nc"
     if(my_task == master_task) write(*,'(/''reading file '',a)') trim(netcdf_riv_file)

     if(lpnetcdf)then
        status=nf90_open(trim(netcdf_riv_file),nf90_nowrite,ncid, &
             & comm=mpi_comm_world,info=mpi_info_null)
     else
        status=nf90_open(trim(netcdf_riv_file),nf90_nowrite,ncid)
     end if
     call handle_error_netcdf('nf90_open: '//trim(netcdf_riv_file),status,nf90_noerr)

     ! get time
     if(iint == 1 .and. k == 1)then
        call read_netcdf_timeinfo_sngl(ncid,'rivtime',rivtime_julday,rivtime_dayint)
     end if
        
     !---it
     it=nint(ihour/(rivtime_dayint*24.d0))+1
     
     if(my_task == master_task)then
        write(6,'(a,i10)') 'river step: ',it
        write(6,'(a)') 'river date: '//yyyy//mm//dd//hh       
     end if

     if(it <= 0)then
        if(my_task == master_task) &
             & write(*,'(/a)') 'wrong time of river data'
        call finalize_mpi
        stop
     end if
     
     if(lpnetcdf)then
        call read_pnetcdf_var2d_time_sngl(ncid,im,jm,it,'river',loc2d_sngl)
        loc2d_size(:,:)=REAL(loc2d_sngl(:,:),r_size)
     else
        call read_netcdf_var2d_time_sngl(ncid,im_global,jm_global,it,'river',glb2d)
        call scatter2d(loc2d_size,dble(glb2d))
     end if
          
     if(k == 1)then
        river0(:,:)=loc2d_size(:,:)
     else if(k == 2)then
        river1(:,:)=loc2d_size(:,:)
     end if

     ! close file
     status=nf90_close(ncid)
     call handle_error_netcdf('nf90_close',status,nf90_noerr)
     
  end do
     
  if(.not. lpnetcdf) deallocate(glb2d)
  
end subroutine read_riv_netcdf
!______________________________________________________________________________
subroutine read_jra55do_netcdf

  !$use omp_lib
  use mpi
  use netcdf
  use common_pom_var
  implicit none

  !---Parameter
  integer,parameter :: nvar=8
  integer,parameter :: ncount=10
  real(kind = r_sngl),parameter :: dmiss=1.e20

  !---Common
  integer ncid,status
  integer iyr,imon,iday,ihour
  integer i,j,k
  integer it
  integer ivar
  integer ijul,ijul0

  real(kind = r_dble) time0_atm  !Initial atmdata time [day]
  real(kind = r_dble),allocatable :: tmp1d_dble(:)
  real(kind = r_sngl),allocatable :: tmp2d_sngl(:,:)
  
  character(200) netcdf_atm_file !Atm filename
  character(4) var(nvar)

  character(10) yyyymmddhh
  character(4) yyyy
  character(2) mm,dd,hh

  !---JRAd55do  
  real(kind = r_size) dat_atm(im_atm,jm_atm)

  !---Model
  real(kind = r_size) loc2d(im,jm)  

  !---Variable name
  var(1)="uas"
  var(2)="vas"
  var(3)="tas"
  var(4)="huss"
  var(5)="psl"
  var(6)="rlds"
  var(7)="rsds"
  var(8)="prra"

  !---Read lon_atm,lat_atm,land_atm
  if(iint == 1)then
     netcdf_atm_file="in/atm/sftof_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr.nc"

     status=nf90_open(trim(netcdf_atm_file),nf90_nowrite,ncid)
     call handle_error_netcdf("nf90_open:"//trim(netcdf_atm_file),status,nf90_noerr)

     allocate(tmp1d_dble(im_atm))
     call read_netcdf_var_dble(ncid,im_atm,1,1,"lon",tmp1d_dble)
     lon_atm(:)=tmp1d_dble(:)
     deallocate(tmp1d_dble)
     
     allocate(tmp1d_dble(jm_atm))
     call read_netcdf_var_dble(ncid,jm_atm,1,1,"lat",tmp1d_dble)
     lat_atm(:)=tmp1d_dble(:)
     deallocate(tmp1d_dble)
     
     allocate(tmp2d_sngl(im_atm,jm_atm))
     call read_netcdf_var_sngl(ncid,im_atm,jm_atm,1,"sftof",tmp2d_sngl)
     land_atm(:,:)=tmp2d_sngl(:,:)
     deallocate(tmp2d_sngl)
     
     status=nf90_close(ncid)
     call handle_error_netcdf("nf90_close:"//trim(netcdf_atm_file),status,nf90_noerr)

     !$omp parallel
     !$omp do private(i,j)    
     do j=1,jm_atm
        do i=1,im_atm
           if(land_atm(i,j) == 100.d0)then
              land_atm(i,j)=1.d0
           else
              land_atm(i,j)=0.d0
           end if
        end do
     end do
     !$omp end do
     !$omp end parallel

     !---ID
     call set_idlon(im_atm,lon_atm,im_local,east_e(:,1),idx_atm)
     call set_idlat(jm_atm,lat_atm,jm_local,north_e(1,:),idy_atm)

  end if

  !---Read atm data
  do k=1,2

     !---iyr,imon,iday,ihour
     if(k == 1)then
        call julian_ymd(int(julday_start+time),iyr,imon,iday)
        ihour=(julday_start+time-int(julday_start+time))*24.d0
     else if(k == 2)then
        call julian_ymd(int(julday_start+time+dt_atm/24.d0),iyr,imon,iday)
        ihour=(julday_start+time+dt_atm/24.d0-int(julday_start+time+dt_atm/24.d0))*24.d0
     end if
     
     write(yyyy,'(i4.4)') iyr
     write(mm,'(i2.2)') imon
     write(dd,'(i2.2)') iday
     write(hh,'(i2.2)') ihour
     yyyymmddhh=yyyy//mm//dd//hh        
     
     if(my_task == master_task) write(6,'(a)') 'atm step: '//yyyymmddhh

     !---it
     call ymd_julian(iyr,imon,iday,ijul)
     call ymd_julian(iyr,1,1,ijul0)
     it=(ijul-ijul0)*8+ihour/3+1
     
     do ivar=1,nvar

        !filename
        if(var(ivar) == "uas" .or. var(ivar) == "vas" &
             & .or. var(ivar) == "tas" .or. var(ivar) == "huss" .or. var(ivar) == "psl")then
           if(2016 <= iyr .and. iyr <= 2023)then
              netcdf_atm_file= &
                   & "in/atm/"//trim(var(ivar))//"_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0-1_gr_"&
                   &//yyyy//"01010000-"//yyyy//"12312100.nc"
           else if(iyr == 2024)then
              netcdf_atm_file= &
                   & "in/atm/"//trim(var(ivar))//"_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0-1_gr_"&
                   &//yyyy//"01010000-"//yyyy//"02012100.nc"
           else
              netcdf_atm_file= &
                   & "in/atm/"//trim(var(ivar))//"_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_"&
                   &//yyyy//"01010000-"//yyyy//"12312100.nc"

           end if
        else if(var(ivar) == "rlds" .or. var(ivar) == "rsds" .or. var(ivar) == "prra")then
           if(2016 <= iyr .and. iyr <= 2023)then
              netcdf_atm_file= &
                   & "in/atm/"//trim(var(ivar))//"_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0-1_gr_"&
                   &//yyyy//"01010130-"//yyyy//"12312230.nc"
           else if(iyr == 2024)then
              netcdf_atm_file= &
                   & "in/atm/"//trim(var(ivar))//"_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0-1_gr_"&
                   &//yyyy//"01010130-"//yyyy//"02012230.nc"
           else
              netcdf_atm_file= &
                   & "in/atm/"//trim(var(ivar))//"_input4MIPs_atmosphericState_OMIP_MRI-JRA55-do-1-5-0_gr_"&
                   &//yyyy//"01010130-"//yyyy//"12312230.nc"
           end if
        end if

        !---Get atmtime_julday,atmtime_dayint
        if(ivar == 1)then
           
           status=nf90_open(trim(netcdf_atm_file),nf90_nowrite,ncid)
           call handle_error_netcdf('nf90_open: '//trim(netcdf_atm_file),status,nf90_noerr)

           call read_netcdf_timeinfo_jra55do(ncid,'time',atmtime_julday,atmtime_dayint)

           time0_atm=dble(atmtime_julday)-dble(julday_start)

           status=nf90_close(ncid)
           call handle_error_netcdf('nf90_close: '//trim(netcdf_atm_file),status,nf90_noerr)

        end if

        !Read data
        status=nf90_open(trim(netcdf_atm_file),nf90_nowrite,ncid)
        call handle_error_netcdf('nf90_open: '//trim(netcdf_atm_file),status,nf90_noerr)

        allocate(tmp2d_sngl(im_atm,jm_atm))
        call read_netcdf_var2d_time_sngl(ncid,im_atm,jm_atm,it,trim(var(ivar)),tmp2d_sngl)
        dat_atm(:,:)=tmp2d_sngl

        status=nf90_close(ncid)
        call handle_error_netcdf("nf90_close:"//trim(netcdf_atm_file),status,nf90_noerr)
        
        !Missing value & Unit
        !$omp parallel
        !$omp do private(i,j)
        do j=1,jm_atm
           do i=1,im_atm

              if(tmp2d_sngl(i,j) == dmiss)then !Missing value
                 dat_atm(i,j)=0.d0
              else if(var(ivar) == "tas")then
                 dat_atm(i,j)=(tmp2d_sngl(i,j)-273.15d0)*land_atm(i,j) ![K] -> [degree C]
              else if(var(ivar) == "prra")then
                 ![kg m-2 s-1] / 1000.[kg/m3] = [m s-1]
                 ![m s-1] * 1000 [mm/m] = ![mm s-1]
                 ![m s-1] * 24*60*60*[s/day] -> [mm day -1]
                 !prep(i,j)=24.d0*60.d0*60.d0*dble(tmp(i,j))
                 dat_atm(i,j)=86400.d0*tmp2d_sngl(i,j)
              else if(var(ivar) == "psl")then
                 dat_atm(i,j)=tmp2d_sngl(i,j)
              else
                 dat_atm(i,j)=tmp2d_sngl(i,j)*land_atm(i,j)
              end if
              
           end do
        end do
        !$omp end do
        !$omp end parallel   

        deallocate(tmp2d_sngl)
        
        !Fill value on the land
        if(var(ivar) /= "prra")then
           call fillvalue_2d(ncount,im_atm,jm_atm,dat_atm,land_atm)
        end if

        !Bilinear interpolation
        call bilinear_interpolation_2d(idx_atm,idy_atm, &
             & im_atm,jm_atm,lon_atm,lat_atm,dat_atm,land_atm, &
             & im,jm,east_e(:,1),north_e(1,:),loc2d,fsm)        
        
        if(k == 1 .and. var(ivar) == "uas")then
           windu0(:,:)=loc2d(:,:)
        else if(k == 1 .and. var(ivar) == "vas")then
           windv0(:,:)=loc2d(:,:)
        else if(k == 1 .and. var(ivar) == "tas")then
           airt0(:,:)=loc2d(:,:)
        else if(k == 1 .and. var(ivar) == "huss")then
           airh0(:,:)=loc2d(:,:)
        else if(k == 1 .and. var(ivar) == "psl")then
           slp0(:,:)=loc2d(:,:)
        else if(k == 1 .and. var(ivar) == "rlds")then
           lwrad0(:,:)=loc2d(:,:)         
        else if(k == 1 .and. var(ivar) == "rsds")then
           swrad0(:,:)=loc2d(:,:)         
        else if(k == 1 .and. var(ivar) == "prra")then
           prep0(:,:)=loc2d(:,:)         
        else if(k == 2 .and. var(ivar) == "uas")then
           windu1(:,:)=loc2d(:,:)
        else if(k == 2 .and. var(ivar) == "vas")then
           windv1(:,:)=loc2d(:,:)
        else if(k == 2 .and. var(ivar) == "tas")then
           airt1(:,:)=loc2d(:,:)
        else if(k == 2 .and. var(ivar) == "huss")then
           airh1(:,:)=loc2d(:,:)
        else if(k == 2 .and. var(ivar) == "psl")then
           slp1(:,:)=loc2d(:,:)         
        else if(k == 2 .and. var(ivar) == "rlds")then
           lwrad1(:,:)=loc2d(:,:)         
        else if(k == 2 .and. var(ivar) == "rsds")then
           swrad1(:,:)=loc2d(:,:)         
        else if(k == 2 .and. var(ivar) == "prra")then
           prep1(:,:)=loc2d(:,:)         
        end if

     end do !ivar
  end do !it
  
end subroutine read_jra55do_netcdf
!______________________________________________________________________________
subroutine read_cama_netcdf

  !$use omp_lib
  use mpi
  use netcdf
  use common_pom_var
  implicit none

  !---Parameter
  real(kind = r_sngl),parameter :: dmiss=9.999e20

  !---Common
  integer ncid,status
  integer iyr,imon,iday,ihour
  integer i,j,k
  integer it
  integer ijul,ijul0

  real(kind = r_dble) time0_riv  !Initial atmdata time [day]
  real(kind = r_dble),allocatable :: tmp1d_sngl(:)
  real(kind = r_sngl),allocatable :: tmp2d_sngl(:,:)
  
  character(200) netcdf_riv_file !River filename

  character(10) yyyymmddhh
  character(4) yyyy
  character(2) mm,dd,hh

  !---JRAd55do  
  real(kind = r_size) land_riv(im_riv,jm_riv)
  real(kind = r_size) dat_riv(im_riv,jm_riv)

  !---Model
  real(kind = r_size) loc2d(im,jm)  

  if(iint == 1)then
  
     !---Model coast line
     call coast_line(im,jm,fsm,cline)
  
     !---Read lon_riv,lat_riv,land_riv
     netcdf_riv_file="in/river/YEE2_JRA-55_outflw_H198101_GLB025.nc"

     status=nf90_open(trim(netcdf_riv_file),nf90_nowrite,ncid)
     call handle_error_netcdf("nf90_open:"//trim(netcdf_riv_file),status,nf90_noerr)
     
     allocate(tmp1d_sngl(im_riv))
     call read_netcdf_var_dble(ncid,im_riv,1,1,"lon",tmp1d_sngl)
     lon_riv(:)=tmp1d_sngl(:)
     deallocate(tmp1d_sngl)
     
     allocate(tmp1d_sngl(jm_riv))
     call read_netcdf_var_dble(ncid,jm_riv,1,1,"lat",tmp1d_sngl)
     lat_riv(:)=tmp1d_sngl(:)
     deallocate(tmp1d_sngl)
     
     allocate(tmp2d_sngl(im_riv,jm_riv))
     call read_netcdf_var3d_time_sngl(ncid,im_riv,jm_riv,1,1,"outflw",tmp2d_sngl)
     land_riv(:,:)=tmp2d_sngl(:,:)
     deallocate(tmp2d_sngl)
     
     status=nf90_close(ncid)
     call handle_error_netcdf("nf90_close:"//trim(netcdf_riv_file),status,nf90_noerr)

     !$omp parallel
     !$omp do private(i)    
     do i=1,im_riv-1
        dlon_riv(i)=lon_riv(i+1)-lon_riv(i)
     end do
     !$omp end do
     
     !$omp do private(j)
     do j=1,jm_riv-1
        dlat_riv(j)=lat_riv(j+1)-lat_riv(j)
     end do
     !$omp end do
  
     !$omp do private(i,j)    
     do j=1,jm_riv
        do i=1,im_riv
           if(land_riv(i,j) == dmiss)then
              land_riv(i,j)=0.d0
           else
              land_riv(i,j)=1.d0
           end if
        end do
     end do
     !$omp end do
     !$omp end parallel
     
     dlon_riv(im)=dlon_riv(im-1)  
     dlat_riv(jm)=dlat_riv(jm-1)
     
     call coast_line(im_riv,jm_riv,land_riv,cline_riv)
  
     !---ID
     call set_id_cama(im,jm,east_e(:,1),north_e(1,:),cline, &
          & im_riv,jm_riv,lon_riv,lat_riv,cline_riv, &
          & idx_riv,idy_riv,dnum_riv)     

  end if
     
  !---Read river data
  do k=1,2

     !---iyr,imon,iday,ihour
     if(k == 1)then
        call julian_ymd(int(julday_start+time),iyr,imon,iday)
        ihour=(julday_start+time-int(julday_start+time))*24.d0
     else if(k == 2)then
        call julian_ymd(int(julday_start+time+dt_riv/24.d0),iyr,imon,iday)
        ihour=(julday_start+time+dt_riv/24.d0-int(julday_start+time+dt_riv/24.d0))*24.d0
     end if
     
     write(yyyy,'(i4.4)') iyr
     write(mm,'(i2.2)') imon
     write(dd,'(i2.2)') iday
     write(hh,'(i2.2)') ihour
     yyyymmddhh=yyyy//mm//dd//hh        
     
     if(my_task == master_task) write(6,'(a)') 'river step: '//yyyymmddhh

     !---it
     call ymd_julian(iyr,imon,iday,ijul)
     call ymd_julian(iyr,imon,1,ijul0)
     it=(ijul-ijul0)*8+ihour/3+1
     
     !filename
     netcdf_riv_file="in/river/YEE2_JRA-55_outflw_H"//yyyy//mm//"_GLB025.nc"

     !---Get rivtime_julday,rivtime_dayint
           
     status=nf90_open(trim(netcdf_riv_file),nf90_nowrite,ncid)
     call handle_error_netcdf('nf90_open: '//trim(netcdf_riv_file),status,nf90_noerr)

     call read_netcdf_timeinfo_cama(ncid,'time',rivtime_julday,rivtime_dayint)
     
     time0_riv=dble(rivtime_julday)-dble(julday_start)

     !---Read data
     allocate(tmp2d_sngl(im_riv,jm_riv))
     call read_netcdf_var3d_time_sngl(ncid,im_riv,jm_riv,1,it,"outflw",tmp2d_sngl)
     dat_riv(:,:)=tmp2d_sngl(:,:)
     deallocate(tmp2d_sngl)

     status=nf90_close(ncid)
     call handle_error_netcdf("nf90_close:"//trim(netcdf_riv_file),status,nf90_noerr)
        
     !Missing value & Unit
     !$omp parallel
     !$omp do private(i,j)
     do j=1,jm_riv
        do i=1,im_riv

           if(dat_riv(i,j) == dmiss)then
              dat_riv(i,j)=0.d0
           else
              ![m^3/s] --> [mm * m^2/day]/(dx*dy) --> [mm/day]
              !dat(i,j)=dat(i,j)*1.d3*24.d0*60.d0*60.d0/(dx(j)*dy)
              !dx(j)=dlon*pi*a*cos(pi*lati(j)/180.d0)/180.d0
              !dy=dlat*pi*a/180.d0
              dat_riv(i,j)=dat_riv(i,j)*86400.d3*180.d0*180.d0&
                   & /(dlon_riv(i)*dlat_riv(j)*pi*pi*earth*earth*cos(pi*lat_riv(j)/180.d0))
              dat_riv(i,j)=dat_riv(i,j)*cline_riv(i,j)
           end if

        end do
     end do
     !$omp end do
     !$omp end parallel

     !loc2d
     call distribute_river(im_riv,jm_riv,dat_riv,im,jm,loc2d,idx_riv,idy_riv,dnum_riv)

     !loc2d -> river
     if(k == 1)then
        river0(:,:)=loc2d(:,:)
     else if(k == 2)then
        river1(:,:)=loc2d(:,:)
     end if

  end do !it
  
end subroutine read_cama_netcdf
!_______________________________________________________________________
subroutine read_lbc_mclim_netcdf(imon)
  !     read monthly climatology of lateral boundary condition
  !     Created by 2018.08 S.Ohishi

  use mpi
  use netcdf
  use common_pom_var
  implicit none

  integer,intent(in) ::  imon

  integer jmon
  integer ncid,status

  real(kind = r_sngl) loc2dew(jm),loc2dns(im)
  real(kind = r_sngl) loc3dew(jm,kb),loc3dns(im,kb)
  real(kind = r_sngl),allocatable :: glb2dew(:),glb2dns(:)
  real(kind = r_sngl),allocatable :: glb3dew(:,:),glb3dns(:,:)

  character(120) netcdf_lbc_file

  if(.not. lpnetcdf) &
       & allocate(glb2dew(jm_global),glb2dns(im_global),glb3dew(jm_global,kb),glb3dns(im_global,kb))
  
  !     open netcdf file
  write(netcdf_lbc_file,'(''in/'',a,''.lbc.nc'')') trim(netcdf_file)
  if(my_task == master_task)then 
     write(*,'(/''reading file '',a)') trim(netcdf_lbc_file)
  endif

  if(lpnetcdf)then
     status=nf90_open(trim(netcdf_lbc_file),nf90_nowrite,ncid, &
          & comm=mpi_comm_world,info=mpi_info_null)
  else
     status=nf90_open(trim(netcdf_lbc_file),nf90_nowrite,ncid)
  end if
  call handle_error_netcdf('nf90_open: '//netcdf_lbc_file,status,nf90_noerr)
     
  ! check month

  if(my_task == master_task) write(6,*) 'lbc month: ',imon
  if(imon <= 0 .or. 12 < imon)then
     if(my_task == master_task) &
          & write(*,'(/a)') 'wrong month in read_lbc_mclim_netcdf'
     call finalize_mpi
     stop
  end if
  
  !     specify 2d array shape at wester/eastern boudary  
  !     eastern boudary
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,jm,imon,'ele',loc2dew)
     ele0(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,imon,'uabe',loc2dew)
     uabe0(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,imon,'vabe',loc2dew)
     vabe0(:)=loc2dew(:)
     
  else

     call read_netcdf_var1d_time_sngl(ncid,jm_global,imon,'ele',glb2dew)
     call scatter_bnd_ew(ele0,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,imon,'uabe',glb2dew)
     call scatter_bnd_ew(uabe0,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,imon,'vabe',glb2dew)
     call scatter_bnd_ew(vabe0,dble(glb2dew),1)
     
  end if
  
  !     western boudary
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,jm,imon,'elw',loc2dew)
     elw0(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,imon,'uabw',loc2dew)
     uabw0(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,imon,'vabw',loc2dew)
     vabw0(:)=loc2dew(:)
     
  else

     call read_netcdf_var1d_time_sngl(ncid,jm_global,imon,'elw',glb2dew)
     call scatter_bnd_ew(elw0,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,imon,'uabw',glb2dew)
     call scatter_bnd_ew(uabw0,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,imon,'vabw',glb2dew)
     call scatter_bnd_ew(vabw0,dble(glb2dew),1)
     
  end if

  !     specify 3d array shape at wester/eastern boudary      
  !     eastern boudary
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,imon,'tbe',loc3dew)
     tbe0(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,imon,'sbe',loc3dew)
     sbe0(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,imon,'ube',loc3dew)
     ube0(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,imon,'vbe',loc3dew)
     vbe0(:,:)=loc3dew(:,:)
     
  else
     
     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,imon,'tbe',glb3dew)
     call scatter_bnd_ew(tbe0,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,imon,'sbe',glb3dew)
     call scatter_bnd_ew(sbe0,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,imon,'ube',glb3dew)
     call scatter_bnd_ew(ube0,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,imon,'vbe',glb3dew)
     call scatter_bnd_ew(vbe0,dble(glb3dew),kb)
     
  end if
  
  !     western boudary
  if(lpnetcdf)then
     
     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,imon,'tbw',loc3dew)
     tbw0(:,:)=loc3dew(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,imon,'sbw',loc3dew)
     sbw0(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,imon,'ubw',loc3dew)
     ubw0(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,imon,'vbw',loc3dew)
     vbw0(:,:)=loc3dew(:,:)
     
  else
     
     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,imon,'tbw',glb3dew)
     call scatter_bnd_ew(tbw0,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,imon,'sbw',glb3dew)
     call scatter_bnd_ew(sbw0,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,imon,'ubw',glb3dew)
     call scatter_bnd_ew(ubw0,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,imon,'vbw',glb3dew)
     call scatter_bnd_ew(vbw0,dble(glb3dew),kb)
     
  end if
     
  !     specify 2d array shape at northern/southern boudary
  !     northern boudary
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,im,imon,'eln',loc2dns)
     eln0(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,im,imon,'uabn',loc2dns)
     uabn0(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,im,imon,'vabn',loc2dns)
     vabn0(:)=loc2dns(:)
     
  else
     
     call read_netcdf_var1d_time_sngl(ncid,im_global,imon,'eln',glb2dns)
     call scatter_bnd_ns(eln0,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,imon,'uabn',glb2dns)
     call scatter_bnd_ns(uabn0,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,imon,'vabn',glb2dns)
     call scatter_bnd_ns(vabn0,dble(glb2dns),1)
     
  end if
     
  !     southern boudary
  if(lpnetcdf)then
     
     call read_pnetcdf_var1d_time_sngl(ncid,im,imon,'els',loc2dns)
     els0(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,im,imon,'uabs',loc2dns)
     uabs0(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,im,imon,'vabs',loc2dns)
     vabs0(:)=loc2dns(:)
     
  else
     
     call read_netcdf_var1d_time_sngl(ncid,im_global,imon,'els',glb2dns)
     call scatter_bnd_ns(els0,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,imon,'uabs',glb2dns)
     call scatter_bnd_ns(uabs0,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,imon,'vabs',glb2dns)
     call scatter_bnd_ns(vabs0,dble(glb2dns),1)
          
  end if
     
  !     specify 3d array shape at northern/southern boudary      
  !     northern boudary
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,imon,'tbn',loc3dns)
     tbn0(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,imon,'sbn',loc3dns)
     sbn0(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,imon,'ubn',loc3dns)
     ubn0(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,imon,'vbn',loc3dns)
     vbn0(:,:)=loc3dns(:,:)
     
  else

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,imon,'tbn',glb3dns)
     call scatter_bnd_ns(tbn0,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,imon,'sbn',glb3dns)
     call scatter_bnd_ns(sbn0,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,imon,'ubn',glb3dns)
     call scatter_bnd_ns(ubn0,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,imon,'vbn',glb3dns)
     call scatter_bnd_ns(vbn0,dble(glb3dns),kb)

  end if
     
  !     southern boudary
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,imon,'tbs',loc3dns)
     tbs0(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,imon,'sbs',loc3dns)
     sbs0(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,imon,'ubs',loc3dns)
     ubs0(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,imon,'vbs',loc3dns)
     vbs0(:,:)=loc3dns(:,:)
     
  else
     
     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,imon,'tbs',glb3dns)
     call scatter_bnd_ns(tbs0,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,imon,'sbs',glb3dns)
     call scatter_bnd_ns(sbs0,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,imon,'ubs',glb3dns)
     call scatter_bnd_ns(ubs0,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,imon,'vbs',glb3dns)
     call scatter_bnd_ns(vbs0,dble(glb3dns),kb)

  end if
     
  !     next month
  jmon=imon+1
  if(12 < jmon) jmon=1
  if(my_task == master_task) write(6,'(a,i10)') 'lbcmclim month: ', jmon

  !     specify 2d array shape at wester/eastern boudary
  !     eastern boudary
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,jm,jmon,'ele',loc2dew)
     ele1(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,jmon,'uabe',loc2dew)
     uabe1(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,jmon,'vabe',loc2dew)
     vabe1(:)=loc2dew(:)

  else
  
     call read_netcdf_var1d_time_sngl(ncid,jm_global,jmon,'ele',glb2dew)
     call scatter_bnd_ew(ele1,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,jmon,'uabe',glb2dew)
     call scatter_bnd_ew(uabe1,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,jmon,'vabe',glb2dew)
     call scatter_bnd_ew(vabe1,dble(glb2dew),1)

  end if

  !     western boudary
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,jm,jmon,'elw',loc2dew)
     elw1(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,jmon,'uabw',loc2dew)
     uabw1(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,jmon,'vabw',loc2dew)
     vabw1(:)=loc2dew(:)
     
  else
     
     call read_netcdf_var1d_time_sngl(ncid,jm_global,jmon,'elw',glb2dew)
     call scatter_bnd_ew(elw1,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,jmon,'uabw',glb2dew)
     call scatter_bnd_ew(uabw1,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,jmon,'vabw',glb2dew)
     call scatter_bnd_ew(vabw1,dble(glb2dew),1)

  end if
     
  !     specify 3d array shape at wester/eastern boudary  
  !     eastern boudary
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,jmon,'tbe',loc3dew)
     tbe1(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,jmon,'sbe',loc3dew)
     sbe1(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,jmon,'ube',loc3dew)
     ube1(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,jmon,'vbe',loc3dew)
     vbe1(:,:)=loc3dew(:,:)
     
  else

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,jmon,'tbe',glb3dew)
     call scatter_bnd_ew(tbe1,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,jmon,'sbe',glb3dew)
     call scatter_bnd_ew(sbe1,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,jmon,'ube',glb3dew)
     call scatter_bnd_ew(ube1,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,jmon,'vbe',glb3dew)
     call scatter_bnd_ew(vbe1,dble(glb3dew),kb)

  end if
     
  !     western boudary
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,jmon,'tbw',loc3dew)
     tbw1(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,jmon,'sbw',loc3dew)
     sbw1(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,jmon,'ubw',loc3dew)
     ubw1(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,jmon,'vbw',loc3dew)
     vbw1(:,:)=loc3dew(:,:)

  else
     
     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,jmon,'tbw',glb3dew)
     call scatter_bnd_ew(tbw1,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,jmon,'sbw',glb3dew)
     call scatter_bnd_ew(sbw1,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,jmon,'ubw',glb3dew)
     call scatter_bnd_ew(ubw1,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,jmon,'vbw',glb3dew)
     call scatter_bnd_ew(vbw1,dble(glb3dew),kb)

  end if
     
  !     specify 2d array shape at northern/southern boudary
  !     northern boudary
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,im,jmon,'eln',loc2dns)
     eln1(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,im,jmon,'uabn',loc2dns)
     uabn1(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,im,jmon,'vabn',loc2dns)
     vabn1(:)=loc2dns(:)
     
  else

     call read_netcdf_var1d_time_sngl(ncid,im_global,jmon,'eln',glb2dns)
     call scatter_bnd_ns(eln1,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,jmon,'uabn',glb2dns)
     call scatter_bnd_ns(uabn1,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,jmon,'vabn',glb2dns)
     call scatter_bnd_ns(vabn1,dble(glb2dns),1)

  end if
     
  !     southern boudary
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,im,jmon,'els',loc2dns)
     els1(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,im,jmon,'uabs',loc2dns)
     uabs1(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,im,jmon,'vabs',loc2dns)
     vabs1(:)=loc2dns(:)
     
  else

     call read_netcdf_var1d_time_sngl(ncid,im_global,jmon,'els',glb2dns)
     call scatter_bnd_ns(els1,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,jmon,'uabs',glb2dns)
     call scatter_bnd_ns(uabs1,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,jmon,'vabs',glb2dns)
     call scatter_bnd_ns(vabs1,dble(glb2dns),1)

  end if
     
  !     specify 3d array shape at northern/southern boudary          
  !     northern boudary
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,jmon,'tbn',loc3dns)
     tbn1(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,jmon,'sbn',loc3dns)
     sbn1(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,jmon,'ubn',loc3dns)
     ubn1(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,jmon,'vbn',loc3dns)
     vbn1(:,:)=loc3dns(:,:)

  else
     
     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,jmon,'tbn',glb3dns)
     call scatter_bnd_ns(tbn1,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,jmon,'sbn',glb3dns)
     call scatter_bnd_ns(sbn1,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,jmon,'ubn',glb3dns)
     call scatter_bnd_ns(ubn1,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,jmon,'vbn',glb3dns)
     call scatter_bnd_ns(vbn1,dble(glb3dns),kb)

  end if
     
  !     southern boudary
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,jmon,'tbs',loc3dns)
     tbs1(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,jmon,'sbs',loc3dns)
     sbs1(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,jmon,'ubs',loc3dns)
     ubs1(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,jmon,'vbs',loc3dns)
     vbs1(:,:)=loc3dns(:,:)
     
  else
     
     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,jmon,'tbs',glb3dns)
     call scatter_bnd_ns(tbs1,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,jmon,'sbs',glb3dns)
     call scatter_bnd_ns(sbs1,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,jmon,'ubs',glb3dns)
     call scatter_bnd_ns(ubs1,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,jmon,'vbs',glb3dns)
     call scatter_bnd_ns(vbs1,dble(glb3dns),kb)

  end if
  
  !     close file
  status=nf90_close(ncid)
  call handle_error_netcdf('nf90_close',status,nf90_noerr)

  if(.not. lpnetcdf) deallocate(glb2dew,glb2dns,glb3dew,glb3dns)
  
end subroutine read_lbc_mclim_netcdf

!_______________________________________________________________________
subroutine read_lbc_netcdf

  ! read lateral boundary condition

  use mpi
  use netcdf
  use common_pom_var
  implicit none

  integer ncid,status
  integer it

  real(kind = r_sngl) loc2dew(jm),loc2dns(im)
  real(kind = r_sngl) loc3dew(jm,kb),loc3dns(im,kb)
  real(kind = r_sngl),allocatable :: glb2dew(:),glb2dns(:)
  real(kind = r_sngl),allocatable :: glb3dew(:,:),glb3dns(:,:)

  real(kind = r_dble) time0_lbc

  character(120) netcdf_lbc_file

  if(.not. lpnetcdf) &
       & allocate(glb2dew(jm_global),glb2dns(im_global),glb3dew(jm_global,kb),glb3dns(im_global,kb))
  
  ! open netcdf lbc file
  write(netcdf_lbc_file,'(''in/'',a,''.lbc.nc'')') trim(netcdf_file)
  if(my_task == master_task) write(*,'(/''reading file '',a)') trim(netcdf_lbc_file)
  if(lpnetcdf)then
     status=nf90_open(trim(netcdf_lbc_file),nf90_nowrite,ncid, &
          & comm=mpi_comm_world,info=mpi_info_null)
     call handle_error_netcdf('nf90_open: '//trim(netcdf_lbc_file),status,nf90_noerr)
  else
     status=nf90_open(trim(netcdf_lbc_file),nf90_nowrite,ncid)
     call handle_error_netcdf('nf90_open: '//trim(netcdf_lbc_file),status,nf90_noerr)
  end if
     
  ! get time
  if(iint == 1)then
     call read_netcdf_timeinfo_sngl(ncid,'lbctime',lbctime_julday,lbctime_dayint)
  end if

  time0_lbc=dble(lbctime_julday)-dble(julday_start)
  if(my_task == master_task)then 
     write(6,'(a,i10)') 'lbctime_julday: ',lbctime_julday
     write(6,'(a,i10)') 'julday_start: ',julday_start
  end if

  ! get number on step-0

  it=(time-dti/86400.d0-time0_lbc)/lbctime_dayint+1
  if(my_task == master_task)then
     write(6,'(a,f12.5)') 'time: ',time-dti/86400.d0
     write(6,'(a,f12.5)') 'time0_lbc: ',time0_lbc
     write(6,'(a,f12.5)') 'lbctime_dayint: ',lbctime_dayint
     write(6,'(a,i10)')   'lbc step: ',it
  end if

  if(it <= 0) then
     if(my_task == master_task) &
          & write(*,'(/a)') 'wrong time of boundary data'
     call finalize_mpi
     stop
  end if

  ! get variables on step-0
  !east
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'ele',loc2dew)
     ele0(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'uabe',loc2dew)
     uabe0(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'vabe',loc2dew)
     vabe0(:)=loc2dew(:)
     
  else
     
     call read_netcdf_var1d_time_sngl(ncid,jm_global,it,'ele',glb2dew)
     call scatter_bnd_ew(ele0,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,it,'uabe',glb2dew)
     call scatter_bnd_ew(uabe0,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,it,'vabe',glb2dew)
     call scatter_bnd_ew(vabe0,dble(glb2dew),1)

  end if
     
  !west
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'elw',loc2dew)
     elw0(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'uabw',loc2dew)
     uabw0(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'vabw',loc2dew)
     vabw0(:)=loc2dew(:)  
     
  else
     
     call read_netcdf_var1d_time_sngl(ncid,jm_global,it,'elw',glb2dew)
     call scatter_bnd_ew(elw0,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,it,'uabw',glb2dew)
     call scatter_bnd_ew(uabw0,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,it,'vabw',glb2dew)
     call scatter_bnd_ew(vabw0,dble(glb2dew),1)

  end if
     
  !east
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'tbe',loc3dew)
     tbe0(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'sbe',loc3dew)
     sbe0(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'ube',loc3dew)
     ube0(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'vbe',loc3dew)
     vbe0(:,:)=loc3dew(:,:)
     
  else
     
     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'tbe',glb3dew)
     call scatter_bnd_ew(tbe0,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'sbe',glb3dew)
     call scatter_bnd_ew(sbe0,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'ube',glb3dew)
     call scatter_bnd_ew(ube0,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'vbe',glb3dew)
     call scatter_bnd_ew(vbe0,dble(glb3dew),kb)

  end if
     
  !west
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'tbw',loc3dew)
     tbw0(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'sbw',loc3dew)
     sbw0(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'ubw',loc3dew)
     ubw0(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'vbw',loc3dew)
     vbw0(:,:)=loc3dew(:,:)
     
  else
     
     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'tbw',glb3dew)
     call scatter_bnd_ew(tbw0,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'sbw',glb3dew)
     call scatter_bnd_ew(sbw0,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'ubw',glb3dew)
     call scatter_bnd_ew(ubw0,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'vbw',glb3dew)
     call scatter_bnd_ew(vbw0,dble(glb3dew),kb)

  end if
     
  !north
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'eln',loc2dns)
     eln0(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'uabn',loc2dns)
     uabn0(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'vabn',loc2dns)
     vabn0(:)=loc2dns(:)
     
  else
     
     call read_netcdf_var1d_time_sngl(ncid,im_global,it,'eln',glb2dns)
     call scatter_bnd_ns(eln0,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,it,'uabn',glb2dns)
     call scatter_bnd_ns(uabn0,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,it,'vabn',glb2dns)
     call scatter_bnd_ns(vabn0,dble(glb2dns),1)

  end if

  !south
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'els',loc2dns)
     els0(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'uabs',loc2dns)
     uabs0(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'vabs',loc2dns)
     vabs0(:)=loc2dns(:)
     
  else
     
     call read_netcdf_var1d_time_sngl(ncid,im_global,it,'els',glb2dns)
     call scatter_bnd_ns(els0,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,it,'uabs',glb2dns)
     call scatter_bnd_ns(uabs0,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,it,'vabs',glb2dns)
     call scatter_bnd_ns(vabs0,dble(glb2dns),1)

  end if
     
  !north
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'tbn',loc3dns)
     tbn0(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'sbn',loc3dns)
     sbn0(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'ubn',loc3dns)
     ubn0(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'vbn',loc3dns)
     vbn0(:,:)=loc3dns(:,:)
     
  else
     
     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'tbn',glb3dns)
     call scatter_bnd_ns(tbn0,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'sbn',glb3dns)
     call scatter_bnd_ns(sbn0,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'ubn',glb3dns)
     call scatter_bnd_ns(ubn0,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'vbn',glb3dns)
     call scatter_bnd_ns(vbn0,dble(glb3dns),kb)

  end if
     
  !south
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'tbs',loc3dns)
     tbs0(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'sbs',loc3dns)
     sbs0(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'ubs',loc3dns)
     ubs0(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'vbs',loc3dns)
     vbs0(:,:)=loc3dns(:,:)
     
  else
     
     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'tbs',glb3dns)
     call scatter_bnd_ns(tbs0,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'sbs',glb3dns)
     call scatter_bnd_ns(sbs0,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'ubs',glb3dns)
     call scatter_bnd_ns(ubs0,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'vbs',glb3dns)
     call scatter_bnd_ns(vbs0,dble(glb3dns),kb)

  end if
     
  ! get number on step-1

  it=(time-dti/86400.d0+lbctime_dayint-time0_lbc)/lbctime_dayint+1
  if(my_task == master_task) write(6,'(a,i10)') 'lbc step: ',it

  if(it <= 0) then
     if(my_task == master_task) &
          & write(*,'(/a)') 'wrong time of boundary data'
     call finalize_mpi
     stop
  end if

  ! get variables on step-1
  !east
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'ele',loc2dew)
     ele1(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'uabe',loc2dew)
     uabe1(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'vabe',loc2dew)
     vabe1(:)=loc2dew(:)
     
  else
     
     call read_netcdf_var1d_time_sngl(ncid,jm_global,it,'ele',glb2dew)
     call scatter_bnd_ew(ele1,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,it,'uabe',glb2dew)
     call scatter_bnd_ew(uabe1,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,it,'vabe',glb2dew)
     call scatter_bnd_ew(vabe1,dble(glb2dew),1)

  end if
     
  !west
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'elw',loc2dew)
     elw1(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'uabw',loc2dew)
     uabw1(:)=loc2dew(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'vabw',loc2dew)
     vabw1(:)=loc2dew(:)  
     
  else
     
     call read_netcdf_var1d_time_sngl(ncid,jm_global,it,'elw',glb2dew)
     call scatter_bnd_ew(elw1,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,it,'uabw',glb2dew)
     call scatter_bnd_ew(uabw1,dble(glb2dew),1)

     call read_netcdf_var1d_time_sngl(ncid,jm_global,it,'vabw',glb2dew)
     call scatter_bnd_ew(vabw1,dble(glb2dew),1)

  end if

  !east
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'tbe',loc3dew)
     tbe1(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'sbe',loc3dew)
     sbe1(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'ube',loc3dew)
     ube1(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'vbe',loc3dew)
     vbe1(:,:)=loc3dew(:,:)
     
  else
     
     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'tbe',glb3dew)
     call scatter_bnd_ew(tbe1,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'sbe',glb3dew)
     call scatter_bnd_ew(sbe1,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'ube',glb3dew)
     call scatter_bnd_ew(ube1,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'vbe',glb3dew)
     call scatter_bnd_ew(vbe1,dble(glb3dew),kb)

  end if
     
  !west
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'tbw',loc3dew)
     tbw1(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'sbw',loc3dew)
     sbw1(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'ubw',loc3dew)
     ubw1(:,:)=loc3dew(:,:)

     call read_pnetcdf_var2d_time_sngl(ncid,jm,kb,it,'vbw',loc3dew)
     vbw1(:,:)=loc3dew(:,:)
     
  else
     
     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'tbw',glb3dew)
     call scatter_bnd_ew(tbw1,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'sbw',glb3dew)
     call scatter_bnd_ew(sbw1,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'ubw',glb3dew)
     call scatter_bnd_ew(ubw1,dble(glb3dew),kb)

     call read_netcdf_var2d_time_sngl(ncid,jm_global,kb,it,'vbw',glb3dew)
     call scatter_bnd_ew(vbw1,dble(glb3dew),kb)
     
  end if
     
  !north
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'eln',loc2dns)
     eln1(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'uabn',loc2dns)
     uabn1(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'vabn',loc2dns)
     vabn1(:)=loc2dns(:)
     
  else
     
     call read_netcdf_var1d_time_sngl(ncid,im_global,it,'eln',glb2dns)
     call scatter_bnd_ns(eln1,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,it,'uabn',glb2dns)
     call scatter_bnd_ns(uabn1,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,it,'vabn',glb2dns)
     call scatter_bnd_ns(vabn1,dble(glb2dns),1)

  end if

  !south
  if(lpnetcdf)then

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'els',loc2dns)
     els1(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'uabs',loc2dns)
     uabs1(:)=loc2dns(:)

     call read_pnetcdf_var1d_time_sngl(ncid,jm,it,'vabs',loc2dns)
     vabs1(:)=loc2dns(:)
     
  else
     
     call read_netcdf_var1d_time_sngl(ncid,im_global,it,'els',glb2dns)
     call scatter_bnd_ns(els1,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,it,'uabs',glb2dns)
     call scatter_bnd_ns(uabs1,dble(glb2dns),1)

     call read_netcdf_var1d_time_sngl(ncid,im_global,it,'vabs',glb2dns)
     call scatter_bnd_ns(vabs1,dble(glb2dns),1)

  end if
     
  !north
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'tbn',loc3dns)
     tbn1(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'sbn',loc3dns)
     sbn1(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'ubn',loc3dns)
     ubn1(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'vbn',loc3dns)
     vbn1(:,:)=loc3dns(:,:)
     
  else
     
     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'tbn',glb3dns)
     call scatter_bnd_ns(tbn1,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'sbn',glb3dns)
     call scatter_bnd_ns(sbn1,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'ubn',glb3dns)
     call scatter_bnd_ns(ubn1,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'vbn',glb3dns)
     call scatter_bnd_ns(vbn1,dble(glb3dns),kb)

  end if
     
  !south
  if(lpnetcdf)then

     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'tbs',loc3dns)
     tbs1(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'sbs',loc3dns)
     sbs1(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'ubs',loc3dns)
     ubs1(:,:)=loc3dns(:,:)
     
     call read_pnetcdf_var2d_time_sngl(ncid,im,kb,it,'vbs',loc3dns)
     vbs1(:,:)=loc3dns(:,:)
     
  else
     
     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'tbs',glb3dns)
     call scatter_bnd_ns(tbs1,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'sbs',glb3dns)
     call scatter_bnd_ns(sbs1,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'ubs',glb3dns)
     call scatter_bnd_ns(ubs1,dble(glb3dns),kb)

     call read_netcdf_var2d_time_sngl(ncid,im_global,kb,it,'vbs',glb3dns)
     call scatter_bnd_ns(vbs1,dble(glb3dns),kb)

  end if
     
  ! close file
  status=nf90_close(ncid)
  call handle_error_netcdf('nf90_close',status,nf90_noerr)

  if(.not. lpnetcdf) deallocate(glb2dew,glb2dns,glb3dew,glb3dns)
  
end subroutine read_lbc_netcdf
!________________________________________________________________________
subroutine read_tsdata_mclim_netcdf(imon)
  !     read temperature salinity data

  use mpi
  use netcdf
  use common_pom_var
  implicit none

  integer,intent(in) :: imon

  integer jmon
  integer ncid,status

  real(kind = r_sngl) :: loc3d(im,jm,kb)
  real(kind = r_sngl),allocatable :: glb3d(:,:,:)

  character(120) netcdf_tsdata_file

  if(.not. lpnetcdf) allocate(glb3d(im_global,jm_global,kb))
  
  !     open ntecdf file
  write(netcdf_tsdata_file,'(''in/'',a,''.tsdata.nc'')') trim(netcdf_file)
  if(my_task == master_task) &
       & write(*,'(''reading file '',a)') trim(netcdf_tsdata_file)

  if(lpnetcdf)then
     status=nf90_open(trim(netcdf_tsdata_file),nf90_nowrite,ncid, &
          & comm=mpi_comm_world,info=mpi_info_null)
  else
     status=nf90_open(trim(netcdf_tsdata_file),nf90_nowrite,ncid)
  end if
  call handle_error_netcdf('nf90_open: '//trim(netcdf_tsdata_file),status,nf90_noerr)
     
  !     present month
  if(lpnetcdf)then

     call read_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,imon,'temp',loc3d)
     tref0(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,imon,'sal',loc3d)
     sref0(:,:,:)=loc3d(:,:,:)
     
  else

     call read_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,imon,'temp',glb3d)
     call scatter3d(tref0,dble(glb3d))

     call read_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,imon,'sal',glb3d)
     call scatter3d(sref0,dble(glb3d))

  end if
     
  !     next month

  jmon=imon+1
  if(jmon > 12) jmon=jmon-12
  if(my_task == master_task) write(6,'(a,i2)') 'tsdata month ',jmon

  if(lpnetcdf)then

     call read_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,jmon,'temp',loc3d)
     tref1(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,jmon,'sal',loc3d)
     sref1(:,:,:)=loc3d(:,:,:)
     
  else
     
     call read_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,jmon,'temp',glb3d)
     call scatter3d(tref1,dble(glb3d))

     call read_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,jmon,'sal',glb3d)
     call scatter3d(sref1,dble(glb3d))

  end if
     
  !     close file
  status=nf90_close(ncid)
  call handle_error_netcdf('nf90_close',status,nf90_noerr)

  if(.not. lpnetcdf) deallocate(glb3d)

end subroutine read_tsdata_mclim_netcdf
!_______________________________________________________________________
subroutine read_tsdata_netcdf

  ! read temperature salinity data

  use mpi
  use netcdf
  use common_pom_var  
  implicit none

  integer ncid,status
  integer it

  real(kind = r_sngl) loc3d(im,jm,kb)
  real(kind = r_sngl),allocatable :: glb3d(:,:,:)
  real(kind = r_dble) time0_tsdata

  character(120) netcdf_tsdata_file

  if(.not. lpnetcdf) allocate(glb3d(im_global,jm_global,kb))
  
  ! open netcdf tsdata file
  write(netcdf_tsdata_file,'(''in/'',a,''.tsdata.nc'')') trim(netcdf_file)
  if(my_task == master_task) write(*,'(/''reading file '',a)') trim(netcdf_tsdata_file)

  if(lpnetcdf)then
     status=nf90_open(trim(netcdf_tsdata_file),nf90_nowrite,ncid, &
          & comm=mpi_comm_world,info=mpi_info_null)
  else
     status=nf90_open(trim(netcdf_tsdata_file),nf90_nowrite,ncid)
  end if
  call handle_error_netcdf('nf90_open: '//trim(netcdf_tsdata_file),status,nf90_noerr)
     
  ! get time
  if(iint == 1)then
     call read_netcdf_timeinfo_sngl(ncid,'tsdatatime',tsdatatime_julday,tsdatatime_dayint)
  end if

  time0_tsdata=dble(tsdatatime_julday)-dble(julday_start)
  if(my_task == master_task)then
     write(6,'(a,i10)') 'tsdatatime_julday: ',tsdatatime_julday
     write(6,'(a,i10)') 'julday_start: ',julday_start
  end if

  ! get number on step-0

  it=(time-dti/86400.d0-time0_tsdata)/tsdatatime_dayint+1
  if(my_task == master_task)then
     write(6,'(a,f12.5)') 'time: ',time-dti/86400.d0
     write(6,'(a,f12.5)') 'time0_tsdata: ',time0_tsdata
     write(6,'(a,f12.5)') 'tsdatatime_dayint: ',tsdatatime_dayint
     write(6,'(a,i10)')   'tsdata step: ',it
  end if

  if(it <= 0)then
     if(my_task == master_task) write(*,'(/a)') 'wrong time of ts data'
     call finalize_mpi
     stop
  end if

  ! get variables on step-0
  if(lpnetcdf)then

     call read_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,it,'temp',loc3d)
     tref0(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,it,'sal',loc3d)
     sref0(:,:,:)=loc3d(:,:,:)
          
  else
     
     call read_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,it,'temp',glb3d)
     call scatter3d(tref0,dble(glb3d))

     call read_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,it,'sal',glb3d)
     call scatter3d(sref0,dble(glb3d))

  end if
     
  ! get number on step-1

  it=(time-dti/86400.d0+tsdatatime_dayint-time0_tsdata)/tsdatatime_dayint+1
  if(my_task == master_task) write(6,'(a,i10)') 'tsdata step: ',it

  if(it <= 0) then
     if(my_task == master_task) &
          & write(*,'(/a)') 'wrong time of ts data'
     call finalize_mpi
     stop
  end if

  ! get variables on step-1
  if(lpnetcdf)then

     call read_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,it,'temp',loc3d)
     tref1(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,it,'sal',loc3d)
     sref1(:,:,:)=loc3d(:,:,:)
     
  else
     
     call read_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,it,'temp',glb3d)
     call scatter3d(tref1,dble(glb3d))
     
     call read_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,it,'sal',glb3d)
     call scatter3d(sref1,dble(glb3d))

  end if
     
  ! close file
  status=nf90_close(ncid)
  call handle_error_netcdf('nf90_close',status,nf90_noerr)

  if(.not. lpnetcdf) deallocate(glb3d)
  
end subroutine read_tsdata_netcdf
!_______________________________________________________________________
subroutine read_tsclim_monthly_netcdf(imon)

  ! read temperature and salinity climatological data

  use mpi
  use netcdf
  use common_pom_var
  implicit none


  integer,intent(in) :: imon

  integer jmon

  integer ncid,status
  
  real(kind = r_sngl) loc3d(im,jm,kb)
  real(kind = r_sngl),allocatable :: glb3d(:,:,:)

  character(120) netcdf_tsclim_file

  if(.not. lpnetcdf) allocate(glb3d(im_global,jm_global,kb))
  
  ! open netcdf atm file
  write(netcdf_tsclim_file,'(''in/'',a,''.tsclim.nc'')') trim(netcdf_file)
  if(my_task == master_task) write(*,'(/''reading file '',a)') trim(netcdf_tsclim_file)

  if(lpnetcdf)then
     status=nf90_open(trim(netcdf_tsclim_file),nf90_nowrite,ncid, &
          & comm=mpi_comm_world,info=mpi_info_null)
  else
     status=nf90_open(trim(netcdf_tsclim_file),nf90_nowrite,ncid)
  end if
  call handle_error_netcdf('nf90_open: '//trim(netcdf_tsclim_file),status,nf90_noerr)
  
  ! check month
  if(my_task == master_task) write(6,'(a,i2)') 'tsclim month: ',imon
  if(imon <= 0.or. imon > 12) then
     if(my_task == master_task) &
          & write(*,'(/a)') 'wrong month in read_tsclim_netcdf'
     call finalize_mpi
     stop
  end if

  ! get variables on step-0
  if(lpnetcdf)then

     call read_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,imon,'tclimm',loc3d)
     tclim0(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,imon,'sclimm',loc3d)
     sclim0(:,:,:)=loc3d(:,:,:)
     
  else
  
     call read_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,imon,'tclimm',glb3d)
     call scatter3d(tclim0,dble(glb3d))      

     call read_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,imon,'sclimm',glb3d)
     call scatter3d(sclim0,dble(glb3d))      

  end if
     
  ! get month on step-1

  jmon=imon+1
  if(jmon > 12) jmon=jmon-12
  if(my_task == master_task) write(6,'(a,i2)') 'tsclim month: ',jmon

  ! get variables on step-1
    if(lpnetcdf)then

     call read_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,jmon,'tclimm',loc3d)
     tclim1(:,:,:)=loc3d(:,:,:)

     call read_pnetcdf_var3d_time_sngl(ncid,im,jm,kb,jmon,'sclimm',loc3d)
     sclim1(:,:,:)=loc3d(:,:,:)

  else
     
     call read_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,jmon,'tclimm',glb3d)
     call scatter3d(tclim1,dble(glb3d))      

     call read_netcdf_var3d_time_sngl(ncid,im_global,jm_global,kb,jmon,'sclimm',glb3d)
     call scatter3d(sclim1,dble(glb3d))      

  end if
     
  ! close file
  status=nf90_close(ncid)
  call handle_error_netcdf('nf90_close',status,nf90_noerr)

  if(.not. lpnetcdf) deallocate(glb3d)
  
end subroutine read_tsclim_monthly_netcdf
