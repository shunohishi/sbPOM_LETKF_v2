program main

  use mod_read_gtspp
  implicit none

  integer iyr,imon

  call read_argument(iyr,imon)

  write(*,*) "Make gtspp filename",iyr,imon
  call make_filename(imon,12,iyr,imon)

end program main

!------------------------------------------------------------------------------

subroutine read_argument(iyr,imon)

  implicit none

  !---Common
  integer i,length,status

  character(:),allocatable :: arg
  
  intrinsic :: command_argument_count, get_command_argument
  
  !---Out
  integer,intent(out) :: iyr,imon

  do i=1,command_argument_count()

     call get_command_argument(i,length=length,status=status)

     if(status /= 0)then
        write(*,*) "Error: arugument ",status
     else

        allocate(character(length) :: arg)

        call get_command_argument(i,arg,status=status)

        if(i == 1)then
           read(arg,'(I4)') iyr
        else if(i == 2)then
           read(arg,'(I2)') imon
        end if
        
        deallocate(arg)

     end if

  end do
  
end subroutine read_argument
