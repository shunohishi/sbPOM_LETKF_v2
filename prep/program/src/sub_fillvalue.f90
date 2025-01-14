!-----------------------------------------------------------------------
! Fill value in horizontal direction |
!----------------------------------------
!
! ncount: the number of times in fill value
!
!-----------------------------------------------------------------------

subroutine fillvalue_2d(ncount,im,jm,km,dat,rmiss)

  !$use omp_lib
  implicit none
  
  integer i,j,k
  integer di,dj
  integer count,pass

  real(kind = 8) tmp(im,jm,km)
  
  integer,intent(in) :: ncount
  integer,intent(in) :: im,jm,km
  
  real(kind = 8),intent(inout) :: dat(im,jm,km)
  real(kind = 8),intent(in) :: rmiss

  count=1

  do while(count <= ncount)

     !     write(*,*) "The number of times filling value:",count,"/",ncount 

     !$omp parallel
     !$omp do private(i,j,k,pass,di,dj)       
     do k=1,km
        do j=1,jm
           do i=1,im
              if(dat(i,j,k) == rmiss)then

                 tmp(i,j,k)=0.d0
                 pass=0

                 do dj=-1,1
                    do di=-1,1
                       if(i+di < 1 .or. im < i+di)cycle
                       if(j+dj < 1 .or. jm < j+dj)cycle
                       if(dat(i+di,j+dj,k) /= rmiss)then
                          tmp(i,j,k)=tmp(i,j,k)+dat(i+di,j+dj,k)
                          pass=pass+1
                       end if
                    end do
                 end do
                    
                 if(pass==0)then
                    tmp(i,j,k)=rmiss
                 else
                    tmp(i,j,k)=tmp(i,j,k)/dble(pass)
                 end if

              else
                 tmp(i,j,k)=dat(i,j,k)
              end if
           end do
        end do
     end do
     !$omp end do
     !$omp end parallel
     
     dat(:,:,:)=tmp(:,:,:)
     count=count+1

  end do

end subroutine fillvalue_2d

!----------------------------------------------------------
! Fill value in vertical direction |
!-----------------------------------
!
! Use the upper layer value dat(i,j,k-1) if dat(i,j,k) is rmiss
!
!----------------------------------------------------------

subroutine fillvalue_vertical(im,jm,km,fsm,dat,rmiss)

  implicit none

  integer i,j,k

  integer,intent(in) :: im,jm,km

  real(kind = 8),intent(in) :: fsm(im,jm)
  real(kind = 8),intent(in) :: rmiss

  real(kind = 8),intent(inout) :: dat(im,jm,km)

  do j=1,jm
     do i=1,im

        if(fsm(i,j) == 0.d0)then
           dat(i,j,:)=0.d0
           cycle
        end if
        
        do k=2,km
           if(dat(i,j,k) == rmiss)then
              dat(i,j,k)=dat(i,j,k-1)
           end if
        end do

     end do
  end do

end subroutine fillvalue_vertical
