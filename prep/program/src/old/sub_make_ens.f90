subroutine make_ens(nens,im,jm,iyr,imon,iday,ihour,ivar,dat,sd,ens)

  use setting, only: factor
  implicit none
  
  integer iens,nens
  integer i,j,im,jm
  integer iyr,imon,iday,ihour
  integer ivar
  
  real dat(im,jm),sd(im,jm)
  real ens(nens,im,jm)

  real noise(nens)

  do j=1,jm
     do i=1,im

        call standarized_noise(i,j,iyr,imon,iday,ihour,ivar,nens,noise)
        
        do iens=1,nens
           ens(iens,i,j)=dat(i,j)+factor*sd(i,j)*noise(iens)
        end do

     end do
  end do

end subroutine make_ens

!-----------------------------------------------------------------
! Make stanrarized noise |
!-----------------------------------------------------------------

subroutine standarized_noise(i,j,iyr,imon,iday,ihour,ivar,nens,noise)

  implicit none

  integer i,j
  integer iyr,imon,iday,ihour
  integer ivar
  integer iens,nens
  integer iseed,nseed
  integer status,system
  integer,allocatable :: seed(:)

  real noise(nens)

  real mean,sd

  !Read random number
  call random_seed(size=nseed)
  allocate(seed(nseed))
  call random_seed(get=seed)
  seed(:)=seed(:)+i+j+iyr+imon+iday+ihour+ivar
  call random_seed(put=seed)
  call random_number(noise)
  
  !Standarization
  mean=0.
  sd=0.

  do iens=1,nens
     mean=mean+noise(iens)
  end do
  
  mean=mean/real(nens)

  do iens=1,nens
     sd=sd+(noise(iens)-mean)*(noise(iens)-mean)
  end do
  
  sd=sqrt(sd/real(nens))

  do iens=1,nens
     noise(iens)=(noise(iens)-mean)/sd
  end do

  deallocate(seed)

end subroutine standarized_noise
