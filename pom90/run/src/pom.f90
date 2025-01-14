! pom.f

! main program

program pom
  
  use common_pom_var
  use mpi
  implicit none
  
  integer ndd,nhh,nmm,nss
  real(kind = r_size) :: realtime,realtime0

  call sub_realtime(realtime0)
  
  !---initialize model
  call initialize

  !---initialize tide
  !if(ltide)then
  !   call tide_initialize
  !end if
  
  !---main loop

  do iint=1,iend
     
     !if(ltide) call tide_advance ! update time dependent tidal parameters
     !if(mod(iint,20) == 0) call diag_w(6)
     
     call get_time
     call advance ! advance model
     
     !if(error_status /= 0 .and. my_task == master_task)then
     if(error_status /= 0)then
        write(5,'(/a)') 'POM terminated with error'
        call finalize_mpi
        stop
     end if
     
  end do !iint

  !---computer time usage
  if(my_task == master_task)then
     call sub_realtime(realtime)
     realtime0 = realtime - realtime0
     ndd = int(realtime0/86400.d0)
     nhh = int((realtime0-ndd*86400.d0)/3600.d0)
     nmm = int((realtime0-ndd*86400.d0-nhh*3600.d0)/60.d0)
     nss = int(realtime0-ndd*86400.d0-nhh*3600.d0-nmm*60.d0)
     write (6,'(/a,i5.5,a,i5.5,a,i5.5,a,i5.5,a)')  &
          & 'Elapsed time usage = ', &
          & ndd,' days ',nhh,' hours ',nmm,' minutes ',nss,' seconds'

     write (6,'(/a)') 'POM terminated successfully '

  end if

  !---finalize mpi
  call finalize_mpi

  if (my_task == master_task) then
     write(6,'(/a)') "***normal end***"
  end if
  
end program pom

!_______________________________________________________________________
subroutine sub_realtime(realtime)

  use common_pom_var
  implicit none
  
  integer i,j,k
  real(kind = r_size),intent(out) :: realtime
  
  call system_clock(i,j,k)
  realtime=dble(i)/dble(j)
  
end subroutine sub_realtime
