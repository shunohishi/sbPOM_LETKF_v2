!--------------------------------------------------------------------------
! Gauss filter |
!--------------------------------------------------------------------------
!
! sigma: exp(-r^2/sigma^2) in window function
! Apply Gaussian filter within exp(-r^2/sigma^2) > 0.1
!
!--------------------------------------------------------------------------

subroutine gauss_filter(im,jm,lon,lat,fsm,sigma,dat)

  use mod_parameter
  implicit none
  
  real(8),parameter :: rmiss=0.d0

  integer i,j,i1,j1,i2,j2
  integer itmp
  integer idlon(im,jm),idlat(im,jm)
  integer count(im,jm)
  integer pass

  real(8) ratio(im,jm)
  real(8) r
  real(8) sum1(im,jm),sum2(im,jm)
  real(8) tmp(im,jm)

  integer,intent(in) :: im,jm
  real(8),intent(in) :: lon(im),lat(jm)
  real(8),intent(in) :: sigma(im,jm),fsm(im,jm)
  real(8),intent(inout) :: dat(im,jm)

  !---Estimate calculating idlon,idlat
  do j=1,jm
     do i=1,im
        idlon(i,j)=nint(sigma(i,j)*log(10.d0)*180.d0/(pi*earth*abs(lon(2)-lon(1))*cos(lat(j)*pi/180.d0)))
        idlat(i,j)=nint(sigma(i,j)*log(10.d0)*180.d0/(pi*earth*abs(lat(2)-lat(1))))
     end do !i
  end  do !j
  
  write(*,*) "idlon:",maxval(idlon)
  write(*,*) "idlat:",maxval(idlat)
  
  !---Apply Gaussian filter
  do j1=1,jm
     
     if(mod(j1,100) == 1)then
        write(*,'(i6,a,i6,a)') j1,"/",jm," in Gaussian Filter"
     end if
     
     do i1=1,im

        if(fsm(i1,j1) == 0.d0)then
           tmp(i1,j1)=rmiss
           cycle
        end if

        if(sigma(i1,j1) == 0.d0)then
           tmp(i1,j1)=dat(i1,j1)
           cycle
        end if
        
        !---Gaussian filter
        sum1(i1,j1)=0.d0
        sum2(i1,j1)=0.d0
        pass=0
        
        do j2=j1-idlat(i1,j1),j1+idlat(i1,j1)
           if(j2 < 1 .or. jm < j2)cycle
           do i2=i1-idlon(i1,j1),i1+idlon(i1,j1)
              if(i2 < 1 .or. im < i2)cycle

              if(fsm(i2,j2) == 0.d0)cycle
        
              !Calculate distance between (lon(i1),lat(j1)) and (lon(i2),lat(j2))
              call distance(lon(i1),lat(j1),lon(i2),lat(j2),r)
              
              if(r <= sigma(i1,j1)*sqrt(log(10.d0)))then 
                 sum1(i1,j1)= &
                      & sum1(i1,j1)+dat(i2,j2)*exp(-1.d0*r*r/(sigma(i1,j1)*sigma(i1,j1)))
                 sum2(i1,j1)= &
                      & sum2(i1,j1)+exp(-1.d0*r*r/(sigma(i1,j1)*sigma(i1,j1)))
                 pass=pass+1
              end if
              
           end do !i2
        end do !j2
        
        if(pass == 0)then
           tmp(i1,j1)=rmiss
        else
           tmp(i1,j1)=sum1(i1,j1)/sum2(i1,j1)
        end if
                
     end do !i1
  end do !j1

  !dat <-- tmp
  dat(:,:)=tmp(:,:)
  
end subroutine gauss_filter
