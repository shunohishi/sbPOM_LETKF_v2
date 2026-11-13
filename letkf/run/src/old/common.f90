MODULE common

  !=======================================================================
  !
  ! [PURPOSE:] General procedures
  !
  ! [HISTORY:]
  !   2004/07/20 Takemasa MIYOSHI  created
  !   2009/01/23 Takemasa MIYOSHI  modified for SFMT
  !   2023/01/04 Shun OHISHI       modified for sbPOM-LETKF v2
  !   2026/05/29 Shun OHISHI       modified
  !   2026/10/21 Shun OHISHI       moved subroutine monitor*
  !                                modified using Claude
  !=======================================================================
  IMPLICIT NONE

CONTAINS

  !-----------------------------------------------------------------------
  ! Timer |
  !-----------------------------------------------------------------------

  SUBROUTINE TIMER_INIT

    USE common_setting
    IMPLICIT NONE

    INTEGER(kind = 8) i

    CALL SYSTEM_CLOCK(rtimerxx,rate,i)
    CALL SYSTEM_CLOCK(rtimer00)
    CALL CPU_TIME(ctimer00)

  END SUBROUTINE TIMER_INIT

  !------------------------

  SUBROUTINE TIMER(section_name)

    USE common_setting
    USE MPI
    IMPLICIT NONE

    CHARACTER(*),INTENT(IN) :: section_name

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)  
    
    CALL SYSTEM_CLOCK(rtimer)
    CALL CPU_TIME(ctimer)
    IF(myrank == root)THEN
       WRITE(6,'(A)')       "=== TIMER: "//trim(section_name)//" ====================="
       WRITE(6,'(A,F10.2)') " TOTAL CPUTIME [s]:",ctimer
       WRITE(6,'(A,F10.2)') " SECTION CPUTIME [s]:",ctimer-ctimer00
       WRITE(6,'(A,F10.2)') " TOTAL REALTIME [s]:",dble(rtimer-rtimerxx)/REAL(rate, r_dble)
       WRITE(6,'(A,F10.2)') " SECTION REALTIME [s]:",dble(rtimer-rtimer00)/REAL(rate, r_dble)
       WRITE(6,'(A)')       "========================================"
    END IF
    rtimer00=rtimer
    ctimer00=ctimer    

  END SUBROUTINE TIMER

  !-----------------------

  SUBROUTINE TIMER_END

    USE common_setting
    IMPLICIT NONE

    IF(myrank == root)THEN
       CALL SYSTEM_CLOCK(rtimer)  
       CALL CPU_TIME(ctimer)
       WRITE(6,'(A)') "=== SUCCESS: Data Assimilation ===================="
       WRITE(6,'(A,F10.2)') " TOTAL CPUTIME [s]:",ctimer
       WRITE(6,'(A,F10.2)') " TOTAL REALTIME [s]:",dble(rtimer-rtimerxx)/dble(rate)
       WRITE(6,'(A)') "==================================================="
    END IF

  END SUBROUTINE TIMER_END

  !-----------------------------------------------------------------------
  ! DISTANCE BETWEEN TWO POINTS (LONa,LATa)-(LONb,LATb) |
  !------------------------------------------------------
  ! 2026.10 S.Ohishi (Claude) changed to the haversine formula.
  !-----------------------------------------------------------------------

  SUBROUTINE distance_2p(lona,lata,lonb,latb,dist)

    USE common_setting, only: r_size, r_dble, pi, re
    IMPLICIT NONE

    REAL(r_dble),PARAMETER :: d2r = pi/180.d0

    REAL(r_dble) :: phia,phib,dphi,dlam,a

    REAL(r_size),INTENT(IN)  :: lona,lata,lonb,latb
    REAL(r_size),INTENT(OUT) :: dist

    phia = REAL(lata,r_dble)*d2r
    phib = REAL(latb,r_dble)*d2r
    dphi = phib-phia
    dlam = REAL(lonb-lona,r_dble)*d2r

    a = SIN(0.5d0*dphi)**2 + COS(phia)*COS(phib)*SIN(0.5d0*dlam)**2
    dist = REAL(2.d0*REAL(re,r_dble)*ASIN(MIN(1.d0,SQRT(a))), r_size)

  END SUBROUTINE distance_2p

  !-----------------------------------------------------------------------
  ! Observation ID (for Interpolation) |
  !-----------------------------------------------------------------------
  ! idx: western edge of model grid point
  ! idy: southern edge of model grid point
  ! Created by S.Ohishi 2025.06
  ! Modified by S.Ohishi 2026.05
  ! Modified by S.Ohishi (Claude) 2026.10
  !-----------------------------------------------------------------------

  SUBROUTINE obs_id_uniform(no,lono,lato,lon,lat,idx,idy)

    !$ USE OMP_LIB
    USE common_setting, only: r_size, undef, nlon, nlat
    IMPLICIT NONE

    !---Common
    INTEGER ilon,ilat,io

    REAL(r_size) :: lon_min(nlon-1),lon_max(nlon-1)
    REAL(r_size) :: lat_min(nlat-1),lat_max(nlat-1)
    REAL(r_size) :: rdlon,rdlat !Inverse of grid interval [1/degree]

    !---IN
    INTEGER,INTENT(IN) :: no                                 !Number of obs.
    REAL(r_size),INTENT(IN) :: lono(no),lato(no)             !Longitude and Latitude (Obs)
    REAL(r_size),INTENT(IN) :: lon(nlon,nlat),lat(nlon,nlat) !Longitude and Latitude (Model)

    !---OUT
    INTEGER,INTENT(OUT) :: idx(no),idy(no)               !Model grid point

    idx(:)=NINT(undef)
    idy(:)=NINT(undef)

    !---Grid cell range
    DO ilon=1,nlon-1
       lon_min(ilon)=MINVAL(lon(ilon:ilon+1, 1))
       lon_max(ilon)=MAXVAL(lon(ilon:ilon+1, 1))
    END DO

    DO ilat=1,nlat-1
       lat_min(ilat)=MINVAL(lat(1, ilat:ilat+1))
       lat_max(ilat)=MAXVAL(lat(1, ilat:ilat+1))
    END DO

    !---Mean grid interval (Uniform grid)
    rdlon = REAL(nlon-1,r_size) / (lon_max(nlon-1)-lon_min(1))
    rdlat = REAL(nlat-1,r_size) / (lat_max(nlat-1)-lat_min(1))

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(io,ilon,ilat)
    DO io=1,no

       !---Out of domain
       IF(lono(io) < lon_min(1) .OR. lon_max(nlon-1) <= lono(io)) CYCLE
       IF(lato(io) < lat_min(1) .OR. lat_max(nlat-1) <= lato(io)) CYCLE

       !---First guess from the uniform grid interval
       ilon = MIN( MAX( INT((lono(io)-lon_min(1))*rdlon)+1, 1 ), nlon-1 )
       ilat = MIN( MAX( INT((lato(io)-lat_min(1))*rdlat)+1, 1 ), nlat-1 )

       !---Correction (round-off error near the cell boundary)
       !ilon
       DO WHILE(lono(io) < lon_min(ilon))
          ilon = ilon - 1
       END DO
       DO WHILE(lon_max(ilon) <= lono(io))
          ilon = ilon + 1
       END DO

       !ilat
       DO WHILE(lato(io) < lat_min(ilat))
          ilat = ilat - 1
       END DO
       DO WHILE(lat_max(ilat) <= lato(io))
          ilat = ilat + 1
       END DO

       idx(io)=ilon
       idy(io)=ilat

    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE obs_id_uniform

  !----------------------------

  SUBROUTINE obs_id_nonuniform(no,lono,lato,lon,lat,idx,idy)

    !$ USE OMP_LIB
    USE common_setting, only: r_size, undef, nlon, nlat
    IMPLICIT NONE

    !---Parameter
    INTEGER,PARAMETER :: nbx=nlon-1, nby=nlat-1 !Number of search bins (lon,lat)

    !---Common
    INTEGER ilon,ilat,io
    INTEGER ib,jb,ib1,ib2,jb1,jb2,kb,n,ic

    REAL(r_size),ALLOCATABLE :: lon_min(:,:),lon_max(:,:) !Bounding box of each cell
    REAL(r_size),ALLOCATABLE :: lat_min(:,:),lat_max(:,:) !(Allocatable: too large for stack)
    REAL(r_size) :: xmin,xmax,ymin,ymax !Domain range
    REAL(r_size) :: rbx,rby             !Inverse of bin interval

    INTEGER,ALLOCATABLE :: nbin(:)      !Number of cells in each bin
    INTEGER,ALLOCATABLE :: ibin(:)      !Start index of each bin (CSR)
    INTEGER,ALLOCATABLE :: cell(:)      !Cell index (ilon + (ilat-1)*(nlon-1))

    !---IN
    INTEGER,INTENT(IN) :: no                                 !Number of obs.
    REAL(r_size),INTENT(IN) :: lono(no),lato(no)             !Longitude and Latitude (Obs)
    REAL(r_size),INTENT(IN) :: lon(nlon,nlat),lat(nlon,nlat) !Longitude and Latitude (Model)

    !---OUT
    INTEGER,INTENT(OUT) :: idx(no),idy(no)               !Model grid point

    idx(:)=NINT(undef)
    idy(:)=NINT(undef)

    !---Bounding box of each grid cell
    ALLOCATE(lon_min(nlon-1,nlat-1),lon_max(nlon-1,nlat-1))
    ALLOCATE(lat_min(nlon-1,nlat-1),lat_max(nlon-1,nlat-1))

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ilon,ilat)
    DO ilat=1,nlat-1
       DO ilon=1,nlon-1

          lon_min(ilon,ilat)=MINVAL(lon(ilon:ilon+1, ilat:ilat+1))
          lon_max(ilon,ilat)=MAXVAL(lon(ilon:ilon+1, ilat:ilat+1))

          lat_min(ilon,ilat)=MINVAL(lat(ilon:ilon+1, ilat:ilat+1))
          lat_max(ilon,ilat)=MAXVAL(lat(ilon:ilon+1, ilat:ilat+1))

       END DO
    END DO
    !$OMP END PARALLEL DO

    !---Search bins covering the domain
    xmin=MINVAL(lon_min); xmax=MAXVAL(lon_max)
    ymin=MINVAL(lat_min); ymax=MAXVAL(lat_max)
    rbx=REAL(nbx,r_size)/(xmax-xmin)
    rby=REAL(nby,r_size)/(ymax-ymin)

    !---Register each cell to all bins its bounding box overlaps
    !   (cells are registered in the order ilat -> ilon, same as the original search)
    ALLOCATE(nbin(nbx*nby),ibin(nbx*nby+1))

    DO n=1,2 !1: count, 2: fill

       IF(n == 1)THEN
          nbin(:)=0
       ELSE
          ibin(1)=1
          DO kb=1,nbx*nby
             ibin(kb+1)=ibin(kb)+nbin(kb)
          END DO
          ALLOCATE(cell(ibin(nbx*nby+1)-1))
          nbin(:)=0
       END IF

       DO ilat=1,nlat-1
          DO ilon=1,nlon-1

             ib1=MIN(MAX(INT((lon_min(ilon,ilat)-xmin)*rbx)+1,1),nbx)
             ib2=MIN(MAX(INT((lon_max(ilon,ilat)-xmin)*rbx)+1,1),nbx)
             jb1=MIN(MAX(INT((lat_min(ilon,ilat)-ymin)*rby)+1,1),nby)
             jb2=MIN(MAX(INT((lat_max(ilon,ilat)-ymin)*rby)+1,1),nby)

             DO jb=jb1,jb2
                DO ib=ib1,ib2
                   kb=ib+(jb-1)*nbx
                   IF(n == 2) cell(ibin(kb)+nbin(kb))=ilon+(ilat-1)*(nlon-1)
                   nbin(kb)=nbin(kb)+1
                END DO
             END DO

          END DO
       END DO

    END DO

    !---Search the cell containing each obs. in its bin
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(io,ib,jb,kb,ic,ilon,ilat)
    DO io=1,no

       IF(lono(io) < xmin .OR. xmax <= lono(io)) CYCLE
       IF(lato(io) < ymin .OR. ymax <= lato(io)) CYCLE

       ib=MIN(MAX(INT((lono(io)-xmin)*rbx)+1,1),nbx)
       jb=MIN(MAX(INT((lato(io)-ymin)*rby)+1,1),nby)
       kb=ib+(jb-1)*nbx

       DO ic=ibin(kb),ibin(kb+1)-1

          ilon=MOD(cell(ic)-1,nlon-1)+1
          ilat=(cell(ic)-1)/(nlon-1)+1

          IF(lono(io) < lon_min(ilon,ilat) .OR. lon_max(ilon,ilat) <= lono(io)) CYCLE
          IF(lato(io) < lat_min(ilon,ilat) .OR. lat_max(ilon,ilat) <= lato(io)) CYCLE

          idx(io)=ilon
          idy(io)=ilat
          EXIT

       END DO

    END DO
    !$OMP END PARALLEL DO

    DEALLOCATE(nbin,ibin,cell)
    DEALLOCATE(lon_min,lon_max,lat_min,lat_max)

  END SUBROUTINE obs_id_nonuniform

  !-----------------------------------------------------------
  ! Bilinear Interpolation |
  !-----------------------------------------------------------
  !
  ! f01(x0,y1) ___________ f11(x1,y1)
  !  |                       |
  ! y|-----------fxy(x,y)    |
  !  |            |          |
  !  |            |          |
  ! f00(x0,y0) ___|_______ f10(x1,y0)
  !               x
  !-----------------------------------------------------------

  SUBROUTINE bilinear_interpolation(x0,x1,y0,y1,x,y,f00,f10,f01,f11,fxy)

    use common_setting, only: r_size, undef
    IMPLICIT NONE

    !---IN
    REAL(r_size),INTENT(IN)  :: x0,x1,y0,y1,x,y
    REAL(r_size),INTENT(IN)  :: f00,f10,f01,f11

    !---OUT
    REAL(r_size),INTENT(OUT) :: fxy

    IF(x0 == x1 .or. y0 == y1)THEN

       fxy=undef

    ELSE

       fxy=(y1-y)/(y1-y0)*((x1-x)/(x1-x0)*f00+(x-x0)/(x1-x0)*f10) &
            & +(y-y0)/(y1-y0)*((x1-x)/(x1-x0)*f01+(x-x0)/(x1-x0)*f11)

    END IF

  END SUBROUTINE bilinear_interpolation

  !-----------------------------------------------------------------------
  ! Akima interpolation (1D) |
  !-----------------------------------------------------------------------
  ! [Reference] Akima, H., 1970: A new method of interpolation and smooth
  !             curve fitting based on local procedures. J. ACM, 17, 589-602.
  !
  ! [IN]  ndim     : Number of data points (ndim >= 2)
  !       x(ndim)  : Coordinate (strictly monotonic; ascending or descending)
  !       y(ndim)  : Data at x
  !       n        : Number of target points
  !       x5(n)    : Target coordinate
  ! [OUT] y5(n)    : Interpolated data (undef if x5 is outside the range of x)
  !
  ! [Method]
  !   1. Segment slope s(m) = (y(m+1)-y(m))/(x(m+1)-x(m)), m=1,...,ndim-1
  !      Two slopes beyond each end are linearly extrapolated
  !        s(0)    = 2s(1)-s(2),           s(-1)     = 2s(0)-s(1)
  !        s(ndim) = 2s(ndim-1)-s(ndim-2), s(ndim+1) = 2s(ndim)-s(ndim-1)
  !   2. Node slope t(p): weighted mean of the two adjacent segment slopes
  !        t(p) = ( |s(p+1)-s(p)|*s(p-1) + |s(p-1)-s(p-2)|*s(p) )
  !             / ( |s(p+1)-s(p)| + |s(p-1)-s(p-2)| )
  !        (simple mean of s(p-1) and s(p) if both weights are zero)
  !   3. In segment k [x(k),x(k+1)], cubic polynomial with y and t at both ends
  !        y5 = y(k) + t1*d + (3s-2t1-t2)*d**2/h + (t1+t2-2s)*d**3/h**2
  !        (d = x5-x(k), h = x(k+1)-x(k), s = s(k), t1 = t(k), t2 = t(k+1))
  !
  ! [Note] Descending x is treated as ascending -x (sgn = -1)
  !        (Akima interpolation does not depend on the direction of x).
  !        If ndim = 2, linear interpolation is used.
  !-----------------------------------------------------------------------
  ! 2026.10 S.Ohishi (Claude) created
  !-----------------------------------------------------------------------

  SUBROUTINE akima_spline(ndim,x,y,n,x5,y5)

    USE common_setting, only: r_size, undef
    IMPLICIT NONE

    !---Common
    INTEGER j
    INTEGER k,klo,khi,kmid     !Segment index (binary search)
    INTEGER m,ks
    INTEGER sgn           !+1: x ascending, -1: x descending

    REAL(r_size) xt            !sgn*x5(j)
    REAL(r_size) h,d           !Segment width, distance from x(k) (in sgn*x)
    REAL(r_size) s1,s2         !Slopes of the two end segments (extrapolation)
    REAL(r_size) sk(-2:2)      !Segment slopes s(k-2),...,s(k+2)
    REAL(r_size) w1,w2         !Akima weights
    REAL(r_size) t(2)          !Node slopes at x(k) and x(k+1)

    !---IN
    INTEGER,INTENT(IN) :: ndim
    REAL(r_size),INTENT(IN) :: x(ndim),y(ndim)

    INTEGER,INTENT(IN) :: n
    REAL(r_size),INTENT(IN) :: x5(n)

    !---OUT
    REAL(r_size),INTENT(OUT) :: y5(n)

    IF(ndim < 2)THEN
       y5(:)=undef
       RETURN
    END IF

    !---Direction of x
    IF(x(ndim) < x(1))THEN
       sgn = -1
    ELSE
       sgn =  1
    END IF

    DO j=1,n

       xt = sgn*x5(j)

       !---Out of range
       IF(xt < sgn*x(1) .or. sgn*x(ndim) < xt)THEN
          y5(j) = undef
          CYCLE
       END IF

       !---Segment k: x(k) <= x5 <= x(k+1) in sgn*x (binary search)
       klo = 1
       khi = ndim
       DO WHILE(khi-klo > 1)
          kmid = (klo+khi)/2
          IF(xt < sgn*x(kmid))THEN
             khi = kmid
          ELSE
             klo = kmid
          END IF
       END DO
       k = klo

       !---On a data point
       IF(x5(j) == x(k))THEN
          y5(j) = y(k)
          CYCLE
       ELSE IF(x5(j) == x(k+1))THEN
          y5(j) = y(k+1)
          CYCLE
       END IF

       h = sgn*(x(k+1) - x(k))
       d = sgn*(x5(j) - x(k))

       !---Only 2 points: linear interpolation
       IF(ndim == 2)THEN
          y5(j) = y(k) + (y(k+1)-y(k))/h*d
          CYCLE
       END IF

       !---1. Segment slopes s(k-2),...,s(k+2)
       DO m=-2,2

          ks = k+m

          IF(ks < 1)THEN               !Beyond x(1): s(0), s(-1)
             s1 = (y(2)-y(1)) / (sgn*(x(2)-x(1)))
             s2 = (y(3)-y(2)) / (sgn*(x(3)-x(2)))
             sk(m) = s1 + REAL(1-ks,r_size)*(s1-s2)
          ELSE IF(ndim-1 < ks)THEN     !Beyond x(ndim): s(ndim), s(ndim+1)
             s1 = (y(ndim)-y(ndim-1))   / (sgn*(x(ndim)-x(ndim-1)))
             s2 = (y(ndim-1)-y(ndim-2)) / (sgn*(x(ndim-1)-x(ndim-2)))
             sk(m) = s1 + REAL(ks-(ndim-1),r_size)*(s1-s2)
          ELSE
             sk(m) = (y(ks+1)-y(ks)) / (sgn*(x(ks+1)-x(ks)))
          END IF

       END DO

       !---2. Node slopes t(1) at x(k), t(2) at x(k+1)
       !      node p=k  : s(p-2),s(p-1),s(p),s(p+1) = sk(-2),sk(-1),sk(0),sk(1)
       !      node p=k+1: s(p-2),s(p-1),s(p),s(p+1) = sk(-1),sk(0),sk(1),sk(2)
       DO m=1,2

          w1 = ABS(sk(m)   - sk(m-1)) !|s(p+1)-s(p)|   : weight for s(p-1)
          w2 = ABS(sk(m-2) - sk(m-3)) !|s(p-1)-s(p-2)| : weight for s(p)

          IF(w1+w2 == REAL(0.d0, r_size))THEN
             t(m) = 0.5d0*(sk(m-2)+sk(m-1))
          ELSE
             t(m) = (w1*sk(m-2) + w2*sk(m-1)) / (w1+w2)
          END IF

       END DO

       !---3. Cubic polynomial in segment k
       y5(j) = y(k) + t(1)*d &
            & + (3.d0*sk(0) - 2.d0*t(1) - t(2))*d*d/h &
            & + (t(1) + t(2) - 2.d0*sk(0))*d*d*d/(h*h)

    END DO !j

  END SUBROUTINE akima_spline

  !---------------------------------------------------------------------------------
  ! Transformation of a state vector in model space to observation space [x -> H(x)]
  !---------------------------------------------------------------------------------

  SUBROUTINE Trans_XtoY(islot,no,elm,idt,idxt,idxu,idxv,idyt,idyu,idyv,lono,lato,levo, &
       & v3d,v2d,hx)

    !$ USE OMP_LIB
    USE common_setting
    USE MPI
    IMPLICIT NONE

    !---Common
    INTEGER io
    INTEGER i,j
    INTEGER k_min,k_max,km
    INTEGER iv2d,iv3d
    INTEGER ierr

    REAL(r_size) lon0,lon1,lat0,lat1
    REAL(r_size) dat00(1),dat01(1),dat10(1),dat11(1) !Akima interpolated data

    !---IN
    INTEGER,INTENT(IN) :: islot             !ith Time slot
    INTEGER,INTENT(IN) :: no                !Number of observation
    INTEGER,INTENT(IN) :: elm(no)           !Obs element
    INTEGER,INTENT(IN) :: idt(no)           !Time ID
    INTEGER,INTENT(IN) :: idxt(no),idyt(no) !Longitude and Latitude ID at t point
    INTEGER,INTENT(IN) :: idxu(no),idyu(no) !Longitude and Latitude ID at u point
    INTEGER,INTENT(IN) :: idxv(no),idyv(no) !Longitude and Latitude ID at v point

    REAL(r_size),INTENT(IN) :: lono(no)  !Obs longitude (deg E)
    REAL(r_size),INTENT(IN) :: lato(no)  !Obs latitude (deg N)
    REAL(r_size),INTENT(IN) :: levo(no)  !Obs depth [m]

    REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d) !3D state variable in model space
    REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)      !2D state variable in model space

    !---OUT
    REAL(r_size),INTENT(INOUT) :: hx(no) !State variable in obs space

    k_min=MIN(k_sfc,k_btm)
    k_max=MAX(k_sfc,k_btm)
    km=k_max-k_min+1

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(io,iv2d,iv3d,i,j,lon0,lon1,lat0,lat1,dat00,dat01,dat10,dat11)
    DO io=1,no

       !---Check time id
       IF(idt(io) /= islot) CYCLE

       !---Get iv2d/iv3d
       SELECT CASE(elm(io))
       CASE(id_z_obs)
          iv2d=iv2d_z
          iv3d=0
       CASE(id_u_obs)
          iv2d=0
          iv3d=iv3d_u
       CASE(id_v_obs)
          iv2d=0
          iv3d=iv3d_v
       CASE(id_t_obs)
          iv2d=0
          iv3d=iv3d_t
       CASE(id_s_obs)
          iv2d=0
          iv3d=iv3d_s
       CASE DEFAULT
          iv2d=0
          iv3d=0
       END SELECT

       !---Get i,j,lon0,lon1,lat0,lat1
       SELECT CASE(elm(io))
       CASE(id_z_obs,id_t_obs,id_s_obs)
          i=idxt(io)
          j=idyt(io)
          lon0=lont(i,j)
          lon1=lont(i+1,j)
          lat0=latt(i,j)
          lat1=latt(i,j+1)
       CASE(id_u_obs)
          i=idxu(io)
          j=idyu(io)
          lon0=lonu(i,j)
          lon1=lonu(i+1,j)
          lat0=latu(i,j)
          lat1=latu(i,j+1)
       CASE(id_v_obs)
          i=idxv(io)
          j=idyv(io)
          lon0=lonv(i,j)
          lon1=lonv(i+1,j)
          lat0=latv(i,j)
          lat1=latv(i,j+1)
       CASE DEFAULT
          i=0
          j=0
          lon0=undef
          lon1=undef
          lat0=undef
          lat1=undef          
       END SELECT


       !---Get 4 grid data [i:i+1,j:j+1]       
       IF(iv2d /= 0 .and. iv3d == 0)THEN
          dat00(1)=v2d(i  ,j  ,iv2d)
          dat10(1)=v2d(i+1,j  ,iv2d)
          dat01(1)=v2d(i  ,j+1,iv2d)
          dat11(1)=v2d(i+1,j+1,iv2d)
       ELSE IF(iv2d == 0 .and. iv3d /= 0)THEN       
          IF(  depth(i,j  ,k_sfc) <  levo(io) .or. depth(i+1,j  ,k_sfc) <  levo(io) .or. &
               depth(i,j+1,k_sfc) <  levo(io) .or. depth(i+1,j+1,k_sfc) <  levo(io))THEN
             !Surface
             dat00(1)=v3d(i  ,j  ,k_sfc,iv3d)
             dat10(1)=v3d(i+1,j  ,k_sfc,iv3d)
             dat01(1)=v3d(i  ,j+1,k_sfc,iv3d)
             dat11(1)=v3d(i+1,j+1,k_sfc,iv3d)
          ELSE       
             !Akima interpolation in vertical direction
             CALL akima_spline(km,depth(i  ,j  ,k_min:k_max),v3d(i  ,j  ,k_min:k_max,iv3d), &
                  & 1,levo(io),dat00(1))
             CALL akima_spline(km,depth(i+1,j  ,k_min:k_max),v3d(i+1,j  ,k_min:k_max,iv3d), &
                  & 1,levo(io),dat10(1))
             CALL akima_spline(km,depth(i  ,j+1,k_min:k_max),v3d(i  ,j+1,k_min:k_max,iv3d), &
                  & 1,levo(io),dat01(1))
             CALL akima_spline(km,depth(i+1,j+1,k_min:k_max),v3d(i+1,j+1,k_min:k_max,iv3d), &
                  & 1,levo(io),dat11(1))
          END IF
       ELSE
          write(6,'(A,I0)') "***Error: Incorrect obs element ==> ", elm(io)
          call MPI_Abort(MPI_COMM_WORLD,1,ierr)
       END IF

       !---BI
       CALL bilinear_interpolation( &
            & lon0,lon1,lat0,lat1, &
            & lono(io),lato(io), &
            & dat00(1),dat10(1),dat01(1),dat11(1), &
            & hx(io))

    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE Trans_XtoY

  !-----------------------------------------------------------------------
  ! Ensemble Mean/Spread
  !-----------------------------------------------------------------------

  !---Obs space (transposed layout: dx(iobs,imem) = dxT_flat(imem+(iobs-1)*nmem))
  !   Same loop/summation order as the original (nobs,nmem) version
  SUBROUTINE ensemble_mesp_obs(dxT_flat,hxmean,hxsprd)

    !$ USE OMP_LIB
    USE common_setting
    USE MPI
    IMPLICIT NONE

    !---Parameter
    INTEGER,PARAMETER :: nblk=512 !Block size of obs. (nblk*nmem*4 byte per thread in L2)

    !---common
    INTEGER iobs,is,ie
    INTEGER imem
    INTEGER ierr
    INTEGER(8) ib

    LOGICAL lundef

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: dxT_flat(INT(nmem,8)*INT(nobs,8)) !HX^T => HdX^T

    !---OUT
    REAL(r_size),INTENT(OUT) :: hxmean(nobs)   !Ensemble mean of H(X)
    REAL(r_size),INTENT(OUT) :: hxsprd(nobs)   !Ensemble spread of H(X)

    IF(nmem <= 1)THEN
       write(6,'(A,I0)') "***Error: Ensemble size ==> ",nmem
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    END IF

    lundef=.FALSE.

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(is,ie,iobs,imem,ib) REDUCTION(.OR.:lundef)
    DO is=1,nobs,nblk

       ie=MIN(is+nblk-1,nobs)

       !Ensemble mean in obs. space
       hxmean(is:ie)=REAL(0.d0,r_size)
       DO imem=1,nmem
          DO iobs=is,ie
             ib=INT(iobs-1,8)*INT(nmem,8)+INT(imem,8)
             IF(dxT_flat(ib) == undef) lundef=.TRUE.
             hxmean(iobs)=hxmean(iobs)+dxT_flat(ib)
          END DO
       END DO
       hxmean(is:ie)=hxmean(is:ie)/REAL(nmem,r_size)

       !Ensemble perturbation dY=HX-Hxmean and spread
       hxsprd(is:ie)=REAL(0.d0,r_size)
       DO imem=1,nmem
          DO iobs=is,ie
             ib=INT(iobs-1,8)*INT(nmem,8)+INT(imem,8)
             dxT_flat(ib)=dxT_flat(ib)-hxmean(iobs)
             hxsprd(iobs)=hxsprd(iobs)+dxT_flat(ib)**2
          END DO
       END DO
       hxsprd(is:ie)=SQRT(hxsprd(is:ie)/REAL(nmem-1,r_size))

    END DO
    !$OMP END PARALLEL DO

    IF(lundef)THEN
       write(6,'(A)') "***Error: H(X) includes undef"
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    END IF

  END SUBROUTINE ensemble_mesp_obs  

  !---Model space
  SUBROUTINE ensemble_mesp(v3d,v2d,v3dm,v3ds,v2dm,v2ds)

    !$ USE OMP_LIB
    USE common_setting
    USE MPI
    IMPLICIT NONE

    INTEGER i,k
    INTEGER imem,ivd
    INTEGER ierr

    !---IN
    REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nmem,nv3d)
    REAL(r_size),INTENT(IN) :: v2d(nij1,nmem,nv2d)

    !---OUT
    REAL(r_size),INTENT(OUT) :: v3dm(nij1,nlev,nv3d),v3ds(nij1,nlev,nv3d)
    REAL(r_size),INTENT(OUT) :: v2dm(nij1,nv2d),v2ds(nij1,nv2d)

    IF(nmem <= 1)THEN
       write(6,'(A,I0)') "***Error: Ensemble size ==> ",nmem
       call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    END IF

    !---3D ensemble mean and spread
    !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC) PRIVATE(ivd,k,imem,i)
    DO ivd=1,nv3d
       DO k=1,nlev

          !Mean
          v3dm(:,k,ivd) = REAL(0.0d0, r_size)
          DO imem=1,nmem
             DO i=1,nij1
                v3dm(i,k,ivd) = v3dm(i,k,ivd) + v3d(i,k,imem,ivd)
             END DO
          END DO
          v3dm(:,k,ivd) = v3dm(:,k,ivd) / REAL(nmem,r_size)

          !Spread
          v3ds(:,k,ivd) = REAL(0.0d0, r_size)
          DO imem=1,nmem
             DO i=1,nij1
                v3ds(i,k,ivd) = v3ds(i,k,ivd) + (v3d(i,k,imem,ivd)-v3dm(i,k,ivd))**2
             END DO
          END DO

          v3ds(:,k,ivd) = SQRT(v3ds(:,k,ivd) / REAL(nmem-1,r_size))

       END DO
    END DO
    !$OMP END PARALLEL DO

    !---2D ensemble mean and spread
    !$OMP PARALLEL PRIVATE(ivd,imem,i)
    DO ivd=1,nv2d

       !Mean
       !$OMP DO SCHEDULE(STATIC)
       DO i=1,nij1
          v2dm(i,ivd) = REAL(0.d0, r_size)
       END DO
       !$OMP END DO

       DO imem=1,nmem
          !$OMP DO SCHEDULE(STATIC)
          DO i=1,nij1
             v2dm(i,ivd) = v2dm(i,ivd) + v2d(i,imem,ivd)
          END DO
          !$OMP END DO
       END DO

       !$OMP DO SCHEDULE(STATIC)       
       DO i=1,nij1
          v2dm(i,ivd) = v2dm(i,ivd) / REAL(nmem,r_size)
          v2ds(i,ivd) = REAL(0.d0, r_size)
       END DO
       !$OMP END DO

       !Spread
       DO imem=1,nmem
          !$OMP DO SCHEDULE(STATIC)
          DO i=1,nij1
             v2ds(i,ivd) = v2ds(i,ivd) + (v2d(i,imem,ivd)-v2dm(i,ivd))**2
          END DO
          !$OMP END DO
       END DO

       !$OMP DO SCHEDULE(STATIC)
       DO i=1,nij1
          v2ds(i,ivd) = SQRT(v2ds(i,ivd) / REAL(nmem-1,r_size))
       END DO
       !$OMP END DO

    END DO
    !$OMP END PARALLEL

  END SUBROUTINE ensemble_mesp

  !-----------------------------------------------------------------------
  ! Get/Put member imem of transposed matrix dxT_flat (64-bit index)
  !   dx(iobs,imem) = dxT_flat(imem+(iobs-1)*nmem)
  !-----------------------------------------------------------------------
  
  SUBROUTINE get_member_dxT_flat(imem,dxT_flat,vec)

    !$ USE OMP_LIB
    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER iobs

    !---IN
    INTEGER,INTENT(IN)       :: imem

    REAL(r_size),INTENT(IN)  :: dxT_flat(INT(nmem,8)*INT(nobs,8))

    !---OUT
    REAL(r_size),INTENT(OUT) :: vec(nobs)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iobs)
    DO iobs=1,nobs
       vec(iobs)=dxT_flat(INT(iobs-1,8)*INT(nmem,8)+INT(imem,8))
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE get_member_dxT_flat

  !----------------------------
  
  SUBROUTINE put_member_dxT_flat(imem,dxT_flat,vec)

    !$ USE OMP_LIB
    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER iobs

    !---IN
    INTEGER,INTENT(IN)         :: imem

    REAL(r_size),INTENT(IN)    :: vec(nobs)

    !---IN/OUT
    REAL(r_size),INTENT(INOUT) :: dxT_flat(INT(nmem,8)*INT(nobs,8))

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(iobs)
    DO iobs=1,nobs
       dxT_flat(INT(iobs-1,8)*INT(nmem,8)+INT(imem,8))=vec(iobs)
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE put_member_dxT_flat
  
  !=======================================================================
  !  Eigenvalue decomposition
  !    INPUT
  !      n      : Matrix dimension
  !      A(n,n) : Input Matrix
  !
  !    OUTPUT
  !      eval(n)   : Eigenvalue in descending order (i.e. eval(1) is the largest)
  !      evec(n,n) : Eigenvectors
  !      np        : Number of Positive eigenvalues
  !
  !    NOTES
  !      - Work arrays are ALLOCATABLE (heap) to reduce the OpenMP thread stack
  !        (about 400 KB for n=128 if they are automatic arrays)
  !=======================================================================

  SUBROUTINE mtx_eigen(n,A,eval,evec,np)

    USE common_setting,only: r_size,r_dble
    USE mpi
    IMPLICIT NONE

    !---Common
    INTEGER i
    INTEGER lwork,liwork
    INTEGER info
    INTEGER ierr

    INTEGER,ALLOCATABLE :: iwork(:)          !Work for DSYEVD [3+5n]
    REAL(r_dble),ALLOCATABLE :: eval8(:)     !Eigenvalue (ascending) [n]
    REAL(r_dble),ALLOCATABLE :: evec8(:,:)   !Eigenvector [n*n]
    REAL(r_dble),ALLOCATABLE :: work(:)      !Work for DSYEVD [1+6n+2n^2]

    !---IN
    INTEGER,INTENT(IN) :: n

    REAL(r_size),INTENT(IN) :: A(n,n)

    !---OUT
    REAL(r_dble),INTENT(OUT) :: eval(n),evec(n,n)
    INTEGER,INTENT(OUT) :: np

    lwork=1+6*n+2*n*n
    liwork=3+5*n
    ALLOCATE(iwork(liwork),work(lwork),eval8(n),evec8(n,n))
    evec8(:,:)=REAL(A(:,:),r_dble)

    CALL DSYEVD('V','U',n,evec8,n,eval8,work,lwork,iwork,liwork,info)
    !eval8: Ascending order (i.e., eval8(1) is the smallest)

    IF(info /= 0)THEN
       WRITE(6,'(A,I10)') "***Error: DSYEVD, INFO: ",info
       DO i=1,MIN(3,n)
          WRITE(6,'(a,3f12.5)') "A:",A(i,1:MIN(3,n))
       END DO
       DO i=1,MIN(3,n)
          WRITE(6,'(a,3f12.5)') "Evec:",evec8(i,1:MIN(3,n))
       END DO
       CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
    END IF

    np=n
    IF(eval8(n) < 0)THEN
       WRITE(6,'(A)') "***Error (mtx_eigen): All Eigenvalues are Negative"
       CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
    ELSE
       !Reconditioning
       DO i=1,n
          IF(eval8(i) < 0.d0)THEN
             WRITE(6,'(A,I5,ES12.4)') "***Error (mtx_eigen): Negative eigenvalue: ",i,eval8(i)
             CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
             !np=np-1
             !eval8(i)=0.0d0
             !evec8(:,i)=0.0d0
          END IF
       END DO
    END IF

    !---Descending order
    DO i=1,n
       eval(i)=eval8(n+1-i)
       evec(:,i)=evec8(:,n+1-i)
    END DO

    DEALLOCATE(iwork,work,eval8,evec8)

  END SUBROUTINE mtx_eigen

  !=======================================================================
  !  Random number generation with optional seeding (wrapper)
  !=======================================================================
  !    INPUT
  !      ndim      : Dimension size of "rand"
  !      iseed     : Seed number; 0 => use DATE_AND_TIME; otherwise use provided
  !
  !    OUTPUT
  !      rand(ndim): Uniform [0,1) random numbers (r_dble)
  !
  !    NOTES
  !      - Uses Mersenne Twister interface: init_gen_rand, genrand_res53 (r_dble).
  !=======================================================================

  SUBROUTINE com_rand_seed(ndim,iseed,rand)

    USE common_setting,only: r_dble
    USE mt19937,only: init_genrand,genrand_res53
    IMPLICIT NONE

    INTEGER i,jseed    
    INTEGER idate(8)

    !---IN
    INTEGER,INTENT(IN) :: ndim  !Dimension size
    INTEGER,INTENT(IN) :: iseed !iseed: 0 ==> Date-base; otherwise ==> iseed-base

    !---OUT
    REAL(r_dble),INTENT(OUT) :: rand(ndim)

    IF(iseed==0)THEN
       CALL DATE_AND_TIME(VALUES=idate) !idate(8) => msec, idate(7) => sec
       jseed=idate(8)+idate(7)*1000
    ELSE
       jseed = iseed
    ENDIF

    CALL init_genrand(jseed)

    DO i=1,ndim
       rand(i) = genrand_res53()
    END DO

  END SUBROUTINE com_rand_seed

  !=======================================================================
  !  Thread-safe uniform random numbers [0,1) (L'Ecuyer 1988, Schrage method)
  !    seed1, seed2: positive integers (e.g., from grid index, level, monte loop)
  !    The state is local, so this routine can be called from OpenMP threads.
  !=======================================================================
  SUBROUTINE rand_lecuyer(seed1,seed2,ndim,rand)

    USE common_setting,only: r_dble
    IMPLICIT NONE

    INTEGER,PARAMETER :: im1=2147483563, im2=2147483399
    INTEGER,PARAMETER :: ia1=40014, ia2=40692
    INTEGER,PARAMETER :: iq1=53668, iq2=52774
    INTEGER,PARAMETER :: ir1=12211, ir2=3791
    INTEGER,PARAMETER :: imm1=im1-1

    INTEGER,INTENT(IN) :: seed1,seed2,ndim
    REAL(r_dble),INTENT(OUT) :: rand(ndim)

    INTEGER i,k,s1,s2,z

    s1 = MAX(1,MOD(ABS(seed1),imm1))
    s2 = MAX(1,MOD(ABS(seed2),im2-1))

    !warm-up to decorrelate nearby seeds
    DO i=1,8
       k=s1/iq1
       s1=ia1*(s1-k*iq1)-k*ir1
       IF(s1 < 0) s1=s1+im1
       k=s2/iq2
       s2=ia2*(s2-k*iq2)-k*ir2
       IF(s2 < 0) s2=s2+im2
    END DO

    DO i=1,ndim
       k=s1/iq1
       s1=ia1*(s1-k*iq1)-k*ir1
       IF(s1 < 0) s1=s1+im1
       k=s2/iq2
       s2=ia2*(s2-k*iq2)-k*ir2
       IF(s2 < 0) s2=s2+im2
       z=s1-s2
       IF(z < 1) z=z+imm1
       rand(i)=REAL(z,r_dble)/REAL(im1,r_dble)
    END DO

  END SUBROUTINE rand_lecuyer

  !=======================================================================
  !  Quick sort (ascending) with index tracking
  !=======================================================================
  !
  !    INPUT
  !      n: Size of "dat"
  !      idx_s, idx_e : Sorting range [idx_s,idx_e]
  !
  !    IN/OUTPUT
  !      idx: Array Index
  !      dat: Array to be sorted (r_dble)
  !
  !=======================================================================

  RECURSIVE SUBROUTINE quick_sort_asend(n,dat,idx,idx_s,idx_e)

    USE common_setting,only: r_dble
    IMPLICIT NONE

    INTEGER i,j,it
    REAL(r_dble) x,t

    !---IN
    INTEGER,INTENT(IN) :: n
    INTEGER,INTENT(IN) :: idx_s,idx_e

    !---INOUT
    INTEGER,INTENT(INOUT) :: idx(n)

    REAL(r_dble),INTENT(INOUT) :: dat(n)

    x = dat((idx_s+idx_e)/2)
    i = idx_s
    j = idx_e

    do

       do while (dat(i) < x)
          i=i+1
       end do

       do while (x < dat(j))
          j=j-1
       end do

       if(j <= i) exit

       t=dat(i)
       dat(i)=dat(j)
       dat(j)=t

       it=idx(i)
       idx(i)=idx(j)
       idx(j)=it

       i=i+1
       j=j-1

    end do

    if(idx_s < i-1) call quick_sort_asend(n,dat,idx,idx_s,i-1)
    if(j+1 < idx_e) call quick_sort_asend(n,dat,idx,j+1,idx_e)

  END SUBROUTINE quick_sort_asend

END MODULE common
