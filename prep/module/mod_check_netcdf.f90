module mod_check_netcdf

contains

  !-------------------------------------------
  ! Check netcdf |
  !-------------------------------------------
  
  subroutine check_netcdf(status)
    
    use netcdf
    implicit none

    integer,intent(in) :: status

    if(status /= nf90_noerr)then

       write(*,'(A)') '  '//trim(nf90_strerror(status))
       stop 1
       
    end if

  end subroutine check_netcdf

end module mod_check_netcdf
