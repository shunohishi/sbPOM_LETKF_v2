!-----------------------------------------------------------
! Read Model SSH |
!-----------------------------------------------------------

subroutine read_ssh(iyr,imon,iday,ssh)

  use netcdf
  use setting, only: pom_dirname,pom_filename
  use mod_rmiss
  use mod_gridinfo
  implicit none

  !---Parameter
  real(kind = 4),parameter :: dmiss=0.e0

  !---Common
  integer i,j
  integer status,access,ncid,varid

  real(kind = 4) tmp(im,jm)

  character(100) filename
  character(4) yyyy
  character(2) mm,dd

  !---IN
  integer,intent(in) :: iyr,imon,iday

  !---OUT
  real(kind = 8),intent(out) :: ssh(im,jm)
    
  write(yyyy,'(i4.4)') iyr
  write(mm,'(i2.2)') imon
  write(dd,'(i2.2)') iday

  filename="../../"//trim(pom_dirname)//"/output/"//yyyy//mm//&
       & "/mean/"//trim(pom_filename)//yyyy//mm//dd//".nc"

  status=access(trim(filename)," ")
  if(status /= 0)then
     write(*,'(a)') "***Error: Not Found "//trim(filename)
     stop
  end if

  status=nf90_open(trim(filename),nf90_nowrite,ncid)

  status=nf90_inq_varid(ncid,"el",varid)
  status=nf90_get_var(ncid,varid,tmp,(/1,1,1/),(/im,jm,1/))
  
  status=nf90_close(ncid)

  do j=1,jm
     do i=1,im
        if(tmp(i,j) == dmiss)then
           ssh(i,j)=rmiss
        else
           ssh(i,j)=dble(tmp(i,j))
        end if
     end do
  end do

end subroutine read_ssh
