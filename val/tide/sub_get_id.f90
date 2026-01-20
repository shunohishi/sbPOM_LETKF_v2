subroutine get_idt(ntime,ijul,tgt,idt)

  implicit none

  !---Common
  integer i
  
  !---IN
  integer,intent(in) :: ntime
  integer,intent(in) :: ijul(ntime)
  integer,intent(in) :: tgt

  !---OUT
  integer,intent(out) :: idt

  idt=0
  
  do i=1,ntime
     if(ijul(i) == tgt)then
        idt=i
        return
     end if
  end do

end subroutine get_idt

!----------------------------------------------------

subroutine get_id(im,jm,lon,lat,mask,tlon,tlat,idx,idy,dist)

  implicit none

  !---Parameter
  integer,parameter :: ni=5,nj=5

  !---Common
  integer i,j
  integer di,dj
  integer idx_min,idy_min
  integer idx_tmp,idy_tmp
  
  real(kind = 8) dlon(im),dlat(jm)
  real(kind = 8) dlon_min,dlat_min

  real(kind = 8) dist_tmp(-ni:ni,-nj:nj)
  
  !---IN
  integer,intent(in) :: im,jm

  real(kind = 8),intent(in) :: lon(im),lat(jm)
  real(kind = 8),intent(in) :: mask(im,jm)

  !Target
  real(kind = 8),intent(in) :: tlon,tlat

  !---OUT
  integer,intent(out) :: idx,idy

  real(kind = 8),intent(out) :: dist !Distance [km]
  
  !Difference
  dlon(:)=lon(:)-tlon
  dlat(:)=lat(:)-tlat

  !Minimum value
  dlon_min=minval(abs(dlon(:)))
  dlat_min=minval(abs(dlat(:)))

  !idx(tmp)
  idx_min=0
  do i=1,im
     if(abs(dlon(i)) == dlon_min)then
        idx_min=i
        exit
     end if
  end do

  !idy(tmp)
  idy_min=0
  do j=1,jm
     if(abs(dlat(j)) == dlat_min)then
        idy_min=j
        exit
     end if
  end do

  !Check mask
  if(mask(idx_min,idy_min) == 1.d0)then
     idx=idx_min
     idy=idy_min
     return
  end if

  !If idx(tmp),idy(tmp) is land,
  idx=0
  idy=0
  do dj=-nj,nj
     do di=-ni,ni

        idx_tmp=idx_min+di
        idy_tmp=idy_min+dj
        if(idx_tmp < 1 .or. im < idx_tmp) cycle
        if(idy_tmp < 1 .or. jm < idy_tmp) cycle

        call distance_2p(lon(idx_tmp),lat(idy_tmp),tlon,tlat,dist_tmp(di,dj))
        
     end do
  end do
     
  dist=maxval(dist_tmp(:,:))
  do dj=-nj,nj
     do di=-ni,ni

        idx_tmp=idx_min+di
        idy_tmp=idy_min+dj
        if(idx_tmp < 1 .or. im < idx_tmp) cycle
        if(idy_tmp < 1 .or. jm < idy_tmp) cycle
        
        if(dist_tmp(di,dj) <= dist .and. mask(idx_tmp,idy_tmp) == 1.d0)then
           idx=idx_tmp
           idy=idy_tmp
           dist=dist_tmp(di,dj)
        end if
        
     end do
  end do
  
end subroutine get_id

!-----------------------------------------------------------------------
! DISTANCE BETWEEN TWO POINTS (LONa,LATa)-(LONb,LATb)
!-----------------------------------------------------------------------
SUBROUTINE distance_2p(lona,lata,lonb,latb,dist)

  IMPLICIT NONE

  !Parameter
  REAL(kind = 8),PARAMETER :: pi=4.d0*atan(1.d0)
  REAL(kind = 8),PARAMETER :: re=6371.3d3
  
  !Common
  REAL(kind = 8) :: cosd

  !IN
  REAL(kind = 8),INTENT(IN) :: lona,lata
  REAL(kind = 8),INTENT(IN) :: lonb,latb
  
  !OUT
  REAL(kind = 8),INTENT(OUT) :: dist
  
  cosd = SIN(pi*lata/180.d0)*SIN(pi*latb/180.d0) &
       & + COS(pi*lata/180.d0)*COS(pi*latb/180.d0)*COS(pi*(lonb-lona)/180.d0)
  cosd = MIN( 1.d0,cosd)
  cosd = MAX(-1.d0,cosd)
  
  dist = re*ACOS(cosd)

END SUBROUTINE distance_2p
