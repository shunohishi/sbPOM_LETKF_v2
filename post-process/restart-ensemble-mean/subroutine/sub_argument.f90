subroutine read_argument(dir,letkf,iyr,imon,iday,nmem)

  implicit none

  !---Common
  integer i,length,status

  character(:),allocatable :: arg

  intrinsic :: command_argument_count, get_command_argument

  !---Out
  character(100),intent(out) :: dir,letkf

  integer,intent(out) :: iyr,imon,iday
  integer,intent(out) :: nmem

  if(command_argument_count() /= 6)then
     write(*,*) "***Error: command_argument_count => ",command_argument_count()
     stop
  end if

  do i=1,command_argument_count()

     call get_command_argument(i,length=length,status=status)

     if(status /= 0)then
        write(*,*) "Error: arugument ",status
     else

        allocate(character(length) :: arg)

        call get_command_argument(i,arg,status=status)

        if(i == 1)then
           dir=arg
        else if(i == 2)then
           letkf=arg
        else if(i == 3)then
           read(arg,'(I4)') iyr
        else if(i == 4)then
           read(arg,'(I2)') imon
        else if(i == 5)then
           read(arg,'(I2)') iday
        else if(i == 6)then
           read(arg,'(I5)') nmem
        end if

        deallocate(arg)

     end if

  end do

end subroutine read_argument
