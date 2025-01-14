!------------------------------------------------------------------
! Detect nearshore grid |
!------------------------------------------------------------------

subroutine detect_nearshore(lon,lat,fsm,nshore)

  use mod_gridinfo
  use setting, only : range => nearshore_range
  use mod_parameter
  implicit none
  
  integer i,j,i1,j1,i2,j2
  integer idlon(jm),idlat
  integer itmp

  real(kind = 8) r
  real(kind = 8),intent(in) :: lon(im),lat(jm),fsm(im,jm)
  real(kind = 8),intent(out) :: nshore(im,jm) !1: Ocean, 0: Land or Nearshore
  
  !---Estimate calculating idlon,idlat
  idlat=nint(range*180.d0/(pi*earth*abs(lat(2)-lat(1))))

  do j1=1,jm

     idlon(j1)=0

     do j2=j1-idlat,j1+idlat

        if(j2 < 1 .or. jm < j2)cycle
           
        itmp=nint(range*180.d0/(pi*earth*abs(lon(2)-lon(1))*cos(lat(j2)*pi/180.d0)))
        if(idlon(j1) < itmp)then
           idlon(j1)=itmp
        end if

     end do
  end do

  write(*,'(a,i6)') "idlat: ",idlat
  write(*,'(a,i6)') "idlon: ",maxval(idlon)
  
  !---Detect nearshore grid
  nshore(:,:)=1.d0
  do j1=1,jm
     
     if(mod(j1,100) == 1)then
        write(*,'(i4.4,a1,i4.4,a)') j1,"/",jm," in Nearshore"
     end if
     
     do i1=1,im
                
        if(fsm(i1,j1) == 0.d0)then
           nshore(i1,j1)=0.d0
           cycle
        end if
        
        do j2=j1-idlat,j1+idlat
           
           if(j2 < 1 .or. jm < j2)cycle
           if(nshore(i1,j1) == 0.d0) cycle

           do i2=i1-idlon(j2),i1+idlon(j2)
              
              if(i2 < 1 .or. im < i2)cycle
              if(nshore(i1,j1) == 0.d0) cycle

              !Calculate distance between (lon(i1),lat(j1)) and (lon(i2),lat(j2))
              call distance(lon(i1),lat(j1),lon(i2),lat(j2),r)
              
              if(r <= range .and. fsm(i2,j2) == 0.d0)then
                 nshore(i1,j1)=0.d0
                 cycle
              end if
              
           end do !i2
        end do !j2
        
     end do !i1
  end do !j1

end subroutine detect_nearshore
