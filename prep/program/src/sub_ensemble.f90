!---------------------------------------------------------------------
! Set ensemble year |
!---------------------------------------------------------------------
!
! syr: dataset start year 
! nyr: dataset period
! iyr,imon,iday,ihour: date for make random seed
! nens: number of ensemble member
! ens_year: different random year
!
!---------------------------------------------------------------------

subroutine set_ens_year(syr,nyr,iyr,imon,iday,ihour,nens,ens_year)

  implicit none

  integer iens
  integer nseed
  integer,allocatable :: seed(:)
  real(kind = 8) noise(nens)
  
  integer,intent(in) :: syr,nyr
  integer,intent(in) :: iyr,imon,iday,ihour
  integer,intent(in) :: nens
  integer,intent(out) :: ens_year(nens)



  call random_seed(size=nseed)
  allocate(seed(nseed))
  call random_seed(get=seed)
  seed(:)=seed(:)+iyr+imon+iday+ihour
  call random_seed(put=seed)
  call random_number(noise)

  do iens=1,nens
     ens_year(iens)=syr+nint(noise(iens)*dble(nyr-1)) !JRA55:1979-2017 (39 years)
  end do

end subroutine set_ens_year

!----------------------------------------------------------------------
! Calculate ensemble |
!----------------------------------------------------------------------

subroutine cal_ens(nens,im,jm,km,dat,ens,factor,rmiss)

  implicit none

  integer iens,nens
  integer i,j,k,im,jm,km
  
  real dat(im,jm,km)
  real ens(nens,im,jm,km)
  real ensmean(im,jm,km)

  real factor
  real rmiss

  call cal_ensmean(nens,im,jm,km,ens,ensmean,rmiss)

  do k=1,km
     do j=1,jm
        do i=1,im
           do iens=1,nens
              ens(iens,i,j,k)=dat(i,j,k)+factor*(ens(iens,i,j,k)-ensmean(i,j,k))
           end do
        end do
     end do
  end do

end subroutine cal_ens

!-------------------------------------------------------------------
! Ensemble mean |
!-------------------------------------------------------------------

subroutine cal_ensmean(nens,im,jm,km,ens,ensmean,rmiss)

  implicit none

  integer iens,nens
  integer i,j,k,im,jm,km
  
  real ens(nens,im,jm,km)
  real ensmean(im,jm,km),pass(im,jm,km),miss(im,jm,km)

  real rmiss

  ensmean(:,:,:)=0.
  pass(:,:,:)=0.
  miss(:,:,:)=0.

  do k=1,km
     do j=1,jm
        do i=1,im
           do iens=1,nens
              if(ens(iens,i,j,k) == rmiss)then
                 miss(i,j,k)=miss(i,j,k)+1.
              else
                 ensmean(i,j,k)=ensmean(i,j,k)+ens(iens,i,j,k)
                 pass(i,j,k)=pass(i,j,k)+1.
              end if
           end do
        end do
     end do
  end do
  
  do k=1,km
     do j=1,jm
        do i=1,im
           if(pass(i,j,k) == 0.)then
              ensmean(i,j,k)=rmiss
           else
              ensmean(i,j,k)=ensmean(i,j,k)/pass(i,j,k)
           end if
        end do
     end do
  end do

end subroutine cal_ensmean

!----------------------------------------------------------------------------
! Gaussian-shaped random number |
!--------------------------------
!
! Use of Box Muller's method
! Mean:0, Sd:1
!
!----------------------------------------------------------------------------

subroutine make_gaussian_random(imon,mean,sigma,n,out)

  implicit none

  real(8),parameter :: pi=4.d0*datan(1.d0)
    
  integer imon
  integer i,n
  integer seedsize
  integer,allocatable :: seed(:)

  real(8) ave,sd,ske
  real(8) mean,sigma
  real(8) noise1(n),noise2(n)
  real(8) out(n)
  
  !Random number
  call random_seed(size=seedsize)
  allocate(seed(seedsize))
  call random_seed(get=seed)
  seed(:)=seed(:)+imon
  call random_seed(put=seed)
  call random_number(noise1)
  call random_number(noise2)
  deallocate(seed)

  do i=1,n
     out(i)=dsqrt(-2.d0*dlog(1.d0-noise1(i)))*dcos(2.d0*pi*noise2(i))
  end do
  
  !Average
  ave=0.d0
  do i=1,n
     ave=ave+out(i)
  end do
  ave=ave/dble(n)

  !Standard deviation
  sd=0.d0
  do i=1,n
     sd=sd+(out(i)-ave)**2.d0
  end do
  sd=dsqrt(sd/dble(n-1.d0))

  !Skewness
  ske=0.d0
  do i=1,n
     ske=ske+((out(i)-ave)/sd)**3.d0
  end do
  ske=dble(n)/(dble(n-1.d0)*dble(n-2.d0))*ske
  
  !Normalization
  do i=1,n
     out(i)=(out(i)-ave)/sd
  end do

  !Add/Substract mean, inflate/deflate standard deviation
  do i=1,n
     out(i)=out(i)*sigma+mean
  end do  

  write(*,*) "Mean:",ave,"--->",mean
  write(*,*) "SD:",sd,"--->",sigma
  write(*,*) "Skewness:",ske
  
end subroutine make_gaussian_random
