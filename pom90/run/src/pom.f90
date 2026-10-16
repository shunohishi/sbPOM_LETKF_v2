program pom
  
  use common_pom_var
  use mpi
  implicit none
  
  integer ndd,nhh,nmm,nss
  real(kind = r_dble) :: realtime,realtime0

  call sub_realtime(realtime0)
  
  !---initialize model
  call initialize
  
  !---main loop
  do iint=1,iend
          
     call advance ! advance model
     
     if(error_status /= 0)then
        if(my_task == master_task) write(6,'(/a)') 'POM terminated with error'
        call finalize_mpi
        stop 1
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
  if (my_task == master_task) then
     write(6,'(/a)') "***normal end***"
  end if

  call finalize_mpi
  
end program pom
!_______________________________________________________________________
subroutine sub_realtime(realtime)

  use common_pom_var
  implicit none
  
  integer(kind = 8) icount,irate,imax
  real(kind = r_dble),intent(out) :: realtime
  
  call system_clock(icount,irate,imax)
  realtime=dble(icount)/dble(irate)
  
end subroutine sub_realtime
