MODULE common
  !=======================================================================
  !
  ! [PURPOSE:] General constants and procedures
  !
  ! [ATTENTION:] This module calls 'SFMT.f90'
  !
  ! [HISTORY:]
  !   07/20/2004 Takemasa MIYOSHI  created
  !   01/23/2009 Takemasa MIYOSHI  modified for SFMT
  !   01/04/2023 Shun OHISHI       modified for sbPOM-LETKF v2
  !
  !=======================================================================
  IMPLICIT NONE
  PUBLIC

CONTAINS
  
  !-----------------------------------------------------------------------
  ! DISTANCE BETWEEN TWO POINTS (LONa,LATa)-(LONb,LATb)
  !-----------------------------------------------------------------------
  SUBROUTINE distance_2p(lona,lata,lonb,latb,dist)

    USE common_setting, only: r_size, pi, re
    IMPLICIT NONE

    !Common
    REAL(r_size) :: cosd

    !IN
    REAL(r_size),INTENT(IN) :: lona,lata
    REAL(r_size),INTENT(IN) :: lonb,latb

    !OUT
    REAL(r_size),INTENT(OUT) :: dist

    cosd = SIN(pi*lata/180.d0)*SIN(pi*latb/180.d0) &
         & + COS(pi*lata/180.d0)*COS(pi*latb/180.d0)*COS(pi*(lonb-lona)/180.d0)
    cosd = MIN( 1.d0,cosd)
    cosd = MAX(-1.d0,cosd)

    dist = re*ACOS(cosd)

  END SUBROUTINE distance_2p
  
  !-----------------------------------------------------------------------
  ! (LON,LAT) --> (i,j) conversion
  !   [ORIGINAL AUTHOR:] Masaru Kunii
  !-----------------------------------------------------------------------
  SUBROUTINE obs_id(im,jm,rlon,rlat,no,olon,olat,oi,oj)

    USE common_setting, only: r_size, undef, file_unit
    IMPLICIT NONE

    ! --- Common
    INTEGER i,j,io
    INTEGER idx,idy

    REAL(r_size) :: rlon_min,rlon_max,rlat_min,rlat_max    
    REAL(r_size) :: dist(4),dist2(4),ratio(4)
    REAL(r_size) :: sum_dist

    ! --- IN
    INTEGER,INTENT(IN) :: im,jm !Grid size(Forecast)
    INTEGER,INTENT(IN) :: no    !Number of obs.

    REAL(r_size),INTENT(IN) :: rlon(im,jm),rlat(im,jm)     !Lon, Lat (Forecast)    
    REAL(r_size),INTENT(IN) :: olon(no),olat(no) !Lon, Lat (Obs.)

    ! --- OUT
    REAL(r_size),INTENT(OUT) :: oi(no),oj(no) !i,j (Obs.)

    !$omp parallel
    !$omp do private(io,idx,idy,rlon_min,rlon_max,rlat_min,rlat_max,dist,dist2,sum_dist,ratio)    
    DO io=1,no

       !idx,idy
       idx = int(undef)
       idy = int(undef)

       DO j=1,jm-1

          rlat_min = MINVAL(rlat(1, j:j+1))          
          rlat_max = MAXVAL(rlat(1, j:j+1))

          IF(olat(io) <  rlat_min .OR. rlat_max <= olat(io)) CYCLE

          DO i=1,im-1

             rlon_min = MINVAL(rlon(i:i+1, j:j+1))             
             rlon_max = MAXVAL(rlon(i:i+1, j:j+1))

             IF(olon(io) < rlon_min .OR. rlon_max <= olon(io)) CYCLE

             idx=i
             idy=j

          END DO !i
       END DO !j

       IF(idx == int(undef) .OR. idy == int(undef))THEN
          oi(io) = undef
          oj(io) = undef
          CYCLE
       END IF

       !dist
       CALL distance_2p(rlon(idx  ,idy  ),rlat(idx  ,idy  ),olon(io),olat(io),dist(1))
       CALL distance_2p(rlon(idx+1,idy  ),rlat(idx+1,idy  ),olon(io),olat(io),dist(2))
       CALL distance_2p(rlon(idx  ,idy+1),rlat(idx  ,idy+1),olon(io),olat(io),dist(3))
       CALL distance_2p(rlon(idx+1,idy+1),rlat(idx+1,idy+1),olon(io),olat(io),dist(4))
       
       !ratio
       do i=1,4
          dist2(i)=1.d-6*dist(i)*dist(i)
       end do
       sum_dist = &
            &  dist2(1)*dist2(2)*dist2(3) &
            & +dist2(2)*dist2(3)*dist2(4) &
            & +dist2(3)*dist2(4)*dist2(1) &
            & +dist2(4)*dist2(1)*dist2(2)
            
       ratio(1) = (dist2(2)*dist2(3)*dist2(4))/sum_dist
       ratio(2) = (dist2(3)*dist2(4)*dist2(1))/sum_dist
       ratio(3) = (dist2(4)*dist2(1)*dist2(2))/sum_dist
       ratio(4) = (dist2(1)*dist2(2)*dist2(3))/sum_dist

       !oi,oj
       oi(io)=ratio(1)*idx+ratio(2)*(idx+1) &
            & +ratio(3)*idx+ratio(4)*(idx+1)
       oj(io)=ratio(1)*idy+ratio(2)*idy &
            & +ratio(3)*(idy+1)+ratio(4)*(idy+1)

    END DO !io
    !$omp end do
    !$omp end parallel

  END SUBROUTINE obs_id

  !-----------------------------------------------------------------------
  ! Cubic spline interpolation
  !   [Reference:] Akima, H., 1970: J. ACM, 17, 589-602.
  !-----------------------------------------------------------------------

  SUBROUTINE akima_spline(ndim,x,y,n,x5,y5)

    USE common_setting, only: r_size, undef  
    IMPLICIT NONE

    !IN
    INTEGER,INTENT(IN) :: ndim         ! number of grid points
    REAL(r_size),INTENT(IN) :: x(ndim) ! coordinate
    REAL(r_size),INTENT(IN) :: y(ndim) ! variable

    INTEGER,INTENT(IN) :: n            ! number of targets
    REAL(r_size),INTENT(IN) :: x5(n)   ! target coordinates

    !OUT
    REAL(r_size),INTENT(OUT) :: y5(n)  ! target values

    INTEGER :: i,j,m
    INTEGER :: iflg

    REAL(r_size) :: dydx(5),ddydx(4),t(2),dx21,dx
    REAL(r_size) :: wk

    DO j=1,n

       iflg = 0

       !Get i
       DO i=1,ndim

          IF(x5(j) == x(i))THEN
             y5(j) = y(i)
             iflg = 1
             CYCLE
          END IF

          IF(x5(j) < x(i)) THEN
             iflg = 1
             EXIT
          END IF

       END DO

       !       i-3   i-2   i-1    i    i+1   i+2
       !     ---+-----+-----+---*-+-----+-----+---
       !dydx       1     2     3     4     5
       !ddydx         1     2     3     4
       !t                   1     2

       IF(i==2)THEN
          DO m=3,5
             dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
          END DO
          dydx(2) = 2.0d0*dydx(3) - dydx(4)
          dydx(1) = 2.0d0*dydx(2) - dydx(3)
       ELSE IF(i==3)THEN
          DO m=2,5
             dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
          END DO
          dydx(1) = 2.0d0*dydx(2) - dydx(3)
       ELSE IF(i==ndim)THEN
          DO m=1,3
             dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
          END DO
          dydx(4) = 2.0d0*dydx(3) - dydx(2)
          dydx(5) = 2.0d0*dydx(4) - dydx(3)
       ELSE IF(i==ndim-1)THEN
          DO m=1,4
             dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
          END DO
          dydx(5) = 2.0d0*dydx(4) - dydx(3)
       ELSE
          DO m=1,5
             dydx(m) = (y(i-3+m)-y(i-4+m)) / (x(i-3+m)-x(i-4+m))
          END DO
       END IF

       DO m=1,4
          ddydx(m) = ABS(dydx(m+1) - dydx(m))
       END DO

       DO m=1,2
          wk = ddydx(m+2) + ddydx(m)
          IF(wk == 0)THEN
             t(m) = 0.0d0
          ELSE
             t(m) = (ddydx(m+2)*dydx(m+1)+ddydx(m)*dydx(m+2))/wk
          END IF
       END DO

       dx21 = x(i)-x(i-1)
       dx = x5(j) - x(i-1)

       !y5
       IF(iflg == 0)THEN
          y5(j)=undef
       ELSE
          y5(j) = y(i-1) &
               & + dx*t(1) &
               & + dx*dx*(3.0d0*dydx(3)-2.0d0*t(1)-t(2))/dx21 &
               & + dx*dx*dx*(t(1)+t(2)-2.0d0*dydx(3))/dx21/dx21
       END IF

    END DO !j

  END SUBROUTINE akima_spline
  
  !-----------------------------------------------------------------------
  ! Transformation of variables from model space to an observation space
  !-----------------------------------------------------------------------
  
  SUBROUTINE Trans_XtoY(elm,ri,rj,rlev,v3d,v2d,yobs)

    USE common_setting
    IMPLICIT NONE

    !Common
    INTEGER :: i,j

    REAL(r_size) x1d(nlev) !1D variable in model space at same lon/lat as obs.
    REAL(r_size) y1d(nlev) !1D variable in model space at same lon/lat as obs.
        
    !IN
    REAL(r_size),INTENT(IN) :: elm   !Obs. element
    REAL(r_size),INTENT(IN) :: ri,rj !Obs. ID
    REAL(r_size),INTENT(IN) :: rlev(1)  !Obs. depth [m]
    
    REAL(r_size),INTENT(IN) :: v3d(nlon,nlat,nlev,nv3d) !3D variable in model space
    REAL(r_size),INTENT(IN) :: v2d(nlon,nlat,nv2d)      !2D variable in model space

    !OUT
    REAL(r_size),INTENT(OUT) :: yobs(1) !variable in obs. space

    !i,j
    IF(ri < 1.d0 .or. REAL(nlon,r_size) < ri &
         & .or. rj < 1.d0 .or. REAL(nlat,r_size) < rj)THEN
       yobs(1) = undef
       RETURN
    ELSE
       i=NINT(ri)
       j=NINT(rj)
    END IF
    
    !Land
    IF(fsm(i,j) == hnosea)THEN
       yobs(1) = undef
       RETURN
    END IF

    SELECT CASE(NINT(elm))
    CASE(id_z_obs)
       yobs(1)=v2d(i,j,iv2d_z)
       RETURN
    CASE(id_u_obs)
       x1d(:)=REAL(depth(i,j,:),r_size)
       y1d(:)=v3d(i,j,:,iv3d_u)
    CASE(id_v_obs)
       x1d(:)=REAL(depth(i,j,:),r_size)
       y1d(:)=v3d(i,j,:,iv3d_v)
    CASE(id_t_obs)
       x1d(:)=REAL(depth(i,j,:),r_size)
       y1d(:)=v3d(i,j,:,iv3d_t)
    CASE(id_s_obs)
       x1d(:)=REAL(depth(i,j,:),r_size)
       y1d(:)=v3d(i,j,:,iv3d_s)
    END SELECT

    IF(REAL(rlev(1),r_size) == REAL(0.d0,r_size))THEN
       yobs(1)=y1d(nlev)
    ELSE IF(rlev(1) < x1d(2) .or. REAL(0.d0,r_size) < rlev(1))THEN !*x1d(1):Deepest and dummy, x1d(2): 2nd deepest
       yobs(1)=undef
    ELSE
       CALL akima_spline(nlev,x1d(:),y1d(:),1,rlev(1),yobs(1))
    END IF
    
  END SUBROUTINE Trans_XtoY

  !-----------------------------------------------------------------------
  ! Ensemble Mean/Spread
  !-----------------------------------------------------------------------
  SUBROUTINE ensemble_mesp(v3d,v2d,v3dm,v3ds,v2dm,v2ds)

    !$USE OMP_LIB
    USE common_setting
    IMPLICIT NONE

    INTEGER i,k,ibv,n
    
    !---IN
    REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nbv,nv3d)
    REAL(r_size),INTENT(IN) :: v2d(nij1,nbv,nv2d)

    !---OUT
    REAL(r_size),INTENT(OUT) :: v3dm(nij1,nlev,nv3d),v3ds(nij1,nlev,nv3d)
    REAL(r_size),INTENT(OUT) :: v2dm(nij1,nv2d),v2ds(nij1,nv2d)

    !---3D ensemble mean and spread
    DO n=1,nv3d

       !$OMP PARALLEL DO PRIVATE(k,i,ibv)
       DO k=1,nlev
          DO i=1,nij1

             v3dm(i,k,n) = v3d(i,k,1,n)
             DO ibv=2,nbv
                v3dm(i,k,n) = v3dm(i,k,n) + v3d(i,k,ibv,n)
             END DO
             v3dm(i,k,n) = v3dm(i,k,n) / REAL(nbv,r_size)

          END DO
       END DO
       !$OMP END PARALLEL DO

       !$OMP PARALLEL DO PRIVATE(k,i,ibv)           
       DO k=1,nlev
          DO i=1,nij1

             v3ds(i,k,n) = (v3d(i,k,1,n)-v3dm(i,k,n))**2
             DO ibv=2,nbv
                v3ds(i,k,n) = v3ds(i,k,n) + (v3d(i,k,ibv,n)-v3dm(i,k,n))**2
             END DO
             v3ds(i,k,n) = SQRT(v3ds(i,k,n) / REAL(nbv-1,r_size))
          END DO
       END DO
       !$OMP END PARALLEL DO

    END DO
    
    !---2D ensemble mean and spread
    DO n=1,nv2d

       !$OMP PARALLEL DO PRIVATE(i,ibv)       
       DO i=1,nij1

          v2dm(i,n) = v2d(i,1,n)
          DO ibv=2,nbv
             v2dm(i,n) = v2dm(i,n) + v2d(i,ibv,n)
          END DO
          v2dm(i,n) = v2dm(i,n) / REAL(nbv,r_size)
          
       END DO
       !$OMP END PARALLEL DO

       !$OMP PARALLEL DO PRIVATE(i,ibv)              
       DO i=1,nij1

          v2ds(i,n) = (v2d(i,1,n)-v2dm(i,n))**2
          DO ibv=2,nbv
             v2ds(i,n) = v2ds(i,n) + (v2d(i,ibv,n)-v2dm(i,n))**2
          END DO
          v2ds(i,n) = SQRT(v2ds(i,n) / REAL(nbv-1,r_size))
       END DO
       !$OMP END PARALLEL DO
       
    END DO

  END SUBROUTINE ensemble_mesp

  !-----------------------------------------------------------------------
  ! Monitor departure
  !-----------------------------------------------------------------------
  SUBROUTINE monit_dep(nn,elm,dep,qc)

    USE common_setting
    IMPLICIT NONE

    !Common
    REAL(r_size) :: rmsd_u,rmsd_v,rmsd_t,rmsd_s,rmsd_z
    REAL(r_size) :: bias_u,bias_v,bias_t,bias_s,bias_z
    INTEGER :: n,iu,iv,it,is,iz

    !IN
    INTEGER,INTENT(IN) :: nn
    REAL(r_size),INTENT(IN) :: elm(nn)
    REAL(r_size),INTENT(IN) :: dep(nn)
    INTEGER,INTENT(IN) :: qc(nn)


    rmsd_u = 0.0d0
    rmsd_v = 0.0d0
    rmsd_t = 0.0d0
    rmsd_s = 0.0d0
    rmsd_z = 0.0d0
    bias_u = 0.0d0
    bias_v = 0.0d0
    bias_t = 0.0d0
    bias_s = 0.0d0
    bias_z = 0.0d0
    iu = 0
    iv = 0
    it = 0
    is = 0
    iz = 0

    DO n=1,nn
       IF(qc(n) == 1)THEN
          SELECT CASE(NINT(elm(n)))
          CASE(id_u_obs)
             rmsd_u = rmsd_u + dep(n)**2
             bias_u = bias_u + dep(n)
             iu = iu + 1
          CASE(id_v_obs)
             rmsd_v = rmsd_v + dep(n)**2
             bias_v = bias_v + dep(n)
             iv = iv + 1
          CASE(id_t_obs)
             rmsd_t = rmsd_t + dep(n)**2
             bias_t = bias_t + dep(n)
             it = it + 1
          CASE(id_s_obs)
             rmsd_s = rmsd_s + dep(n)**2
             bias_s = bias_s + dep(n)
             is = is + 1
          CASE(id_z_obs)
             rmsd_z = rmsd_z + dep(n)**2
             bias_z = bias_z + dep(n)
             iz = iz + 1
          END SELECT
       END IF
    END DO !n

    !U
    IF(iu == 0) THEN
       rmsd_u = undef
       bias_u = undef
    ELSE
       rmsd_u = SQRT(rmsd_u / REAL(iu,r_size))
       bias_u = bias_u / REAL(iu,r_size)
    END IF

    !V
    IF(iv == 0) THEN
       rmsd_v = undef
       bias_v = undef
    ELSE
       rmsd_v = SQRT(rmsd_v / REAL(iv,r_size))
       bias_v = bias_v / REAL(iv,r_size)
    END IF

    !T
    IF(it == 0) THEN
       rmsd_t = undef
       bias_t = undef
    ELSE
       rmsd_t = SQRT(rmsd_t / REAL(it,r_size))
       bias_t = bias_t / REAL(it,r_size)
    END IF

    !S
    IF(is == 0) THEN
       rmsd_s = undef
       bias_s = undef
    ELSE
       rmsd_s = SQRT(rmsd_s / REAL(is,r_size))
       bias_s = bias_s / REAL(is,r_size)
    END IF

    !SSH
    IF(iz == 0) THEN
       rmsd_z = undef
       bias_z = undef
    ELSE
       rmsd_z = SQRT(rmsd_z / REAL(iz,r_size))
       bias_z = bias_z / REAL(iz,r_size)
    END IF

    WRITE(file_unit,'(A)') "== OBSERVATIONAL DEPARTURE ================================="
    WRITE(file_unit,'(6A12)') "Variable: ", "U", "V", "T", "S", "SSH"
    WRITE(file_unit,'(A12,5F12.5)') "Bias: ", bias_u, bias_v, bias_t, bias_s, bias_z
    WRITE(file_unit,'(A12,5F12.5)') "RMSD: ", rmsd_u, rmsd_v, rmsd_t, rmsd_s, rmsd_z
    WRITE(file_unit,'(A)') "== NUMBER OF OBSERVATIONS TO BE ASSIMILATED ================"
    WRITE(file_unit,'(6A12)') "Variable: ","U", "V", "T", "S", "SSH"
    WRITE(file_unit,'(A12,5I12)') "#OBS: ", iu, iv, it, is, iz
    WRITE(file_unit,'(A)') "============================================================"

  END SUBROUTINE monit_dep

  !--------------------------------------------------------------------------------
  ! Monitor departure of forecast and analysis mean relative to assimilated obs.|  
  !--------------------------------------------------------------------------------

  SUBROUTINE monit_mean(file,v3dg,v2dg)

    !$USE OMP_LIB
    USE common_setting
    IMPLICIT NONE

    !---Common
    INTEGER i,ivd
    
    REAL(r_size) hdxf(nobs) !Ensenbme mean in obs space
    REAL(r_size) dep(nobs)  !Innovation

    !3D+2D variable
    INTEGER no(nv3d+nv2d)
    REAL(r_size) bias(nv3d+nv2d),rmsd(nv3d+nv2d)

    !3D variable
    INTEGER no3ds(nv3d),no3di(nv3d)
    REAL(r_size) bias3ds(nv3d),rmsd3ds(nv3d) !Sea surface
    REAL(r_size) bias3di(nv3d),rmsd3di(nv3d) !Internal

    CHARACTER(12) :: filename="file_mean.nc"

    !---IN    
    CHARACTER(4),INTENT(IN) :: file
    REAL(r_size),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d),v2dg(nlon,nlat,nv2d)
    
    !---Initialization
    no(:)=0
    no3ds(:)=0
    no3di(:)=0

    bias(:)=0.d0
    bias3ds(:)=0.d0
    bias3di(:)=0.d0

    rmsd(:)=0.d0    
    rmsd3ds(:)=0.d0
    rmsd3di(:)=0.d0

    WRITE(filename(1:4),'(A4)') file

    !---Read data
    !call read_state_vector(filename,v3dg,v2dg)

    !---Innovation
    !$OMP PARALLEL
    !$OMP DO PRIVATE(i,ivd) REDUCTION(+:no,no3ds,no3di,bias,bias3ds,bias3di,rmsd,rmsd3ds,rmsd3di)
    DO i=1,nobs

       !State vector in model space --> obs space
       CALL Trans_XtoY(obselm(i),obsi(i),obsj(i),obslev(i),v3dg,v2dg,hdxf(i))
       IF(hdxf(i) == undef) CYCLE

       !Innovation
       dep(i) = obsdat(i) - hdxf(i)
       
       !Get ID
       SELECT CASE(NINT(obselm(i)))
       CASE(id_u_obs)
          ivd=iv3d_u
       CASE(id_v_obs)
          ivd=iv3d_v
       CASE(id_t_obs)
          ivd=iv3d_t
       CASE(id_s_obs)
          ivd=iv3d_s
       CASE(id_z_obs)
          ivd=nv3d+iv2d_z
       END SELECT
       
       no(ivd)=no(ivd)+1
       bias(ivd)=bias(ivd)+dep(i)
       rmsd(ivd)=rmsd(ivd)+dep(i)**2

       IF(obselm(i) == id_z_obs) CYCLE
       
       !3D surface and internal
       IF(obslev(i) == 0.d0)THEN
          no3ds(ivd)=no3ds(ivd)+1
          bias3ds(ivd)=bias3ds(ivd)+dep(i)
          rmsd3ds(ivd)=rmsd3ds(ivd)+dep(i)**2
       ELSE
          no3di(ivd)=no3di(ivd)+1
          bias3di(ivd)=bias3di(ivd)+dep(i)
          rmsd3di(ivd)=rmsd3di(ivd)+dep(i)**2
       END IF

    END DO !i
    !$OMP END DO

    !$OMP DO PRIVATE(ivd)    
    DO ivd=1,nv3d+nv2d
       IF(no(ivd) == 0)THEN
          bias(ivd)=undef
          rmsd(ivd)=undef
       ELSE
          bias(ivd)=bias(ivd)/no(ivd)
          rmsd(ivd)=SQRT(rmsd(ivd)/no(ivd))
       END IF
    END DO !ivd
    !$OMP END DO

    !$OMP DO PRIVATE(ivd)        
    DO ivd=1,nv3d
       !Surface
       IF(no3ds(ivd) == 0)THEN
          bias3ds(ivd)=undef
          rmsd3ds(ivd)=undef
       ELSE
          bias3ds(ivd)=bias3ds(ivd)/no3ds(ivd)
          rmsd3ds(ivd)=SQRT(rmsd3ds(ivd)/no3ds(ivd))
       END IF

       !Internal
       IF(no3di(ivd) == 0)THEN
          bias3di(ivd)=undef
          rmsd3di(ivd)=undef
       ELSE
          bias3di(ivd)=bias3di(ivd)/no3di(ivd)
          rmsd3di(ivd)=SQRT(rmsd3di(ivd)/no3di(ivd))
       END IF
    END DO !ivd
    !$OMP END DO
    !$OMP END PARALLEL

    !---Write Information
    WRITE(6,'(A,I10)') "TOTAL #OBS:",nobs
    WRITE(6,'(A)') "== PARTIAL OBSERVATIONAL DEPARTURE ("//file//") =================="
    WRITE(6,'(14A11)') &
         & "Variable:", &
         & "U","SSU","Interior U", &
         & "V","SSV","Interior V", &
         & "T","SST","Interior T", &
         & "S","SSS","Interior S", &
         & "SSH"
    WRITE(6,'(A11,13F11.5)') &
         "Bias:", &
         & bias(iv3d_u),bias3ds(iv3d_u),bias3di(iv3d_u), &
         & bias(iv3d_v),bias3ds(iv3d_v),bias3di(iv3d_v), &
         & bias(iv3d_t),bias3ds(iv3d_t),bias3di(iv3d_t), &
         & bias(iv3d_s),bias3ds(iv3d_s),bias3di(iv3d_s), &
         & bias(nv3d+iv2d_z)
    WRITE(6,'(A11,13F11.5)') &
         "RMSD:", &
         & rmsd(iv3d_u),rmsd3ds(iv3d_u),rmsd3di(iv3d_u), &
         & rmsd(iv3d_v),rmsd3ds(iv3d_v),rmsd3di(iv3d_v), &
         & rmsd(iv3d_t),rmsd3ds(iv3d_t),rmsd3di(iv3d_t), &
         & rmsd(iv3d_s),rmsd3ds(iv3d_s),rmsd3di(iv3d_s), &
         & rmsd(nv3d+iv2d_z)

    WRITE(6,'(A)') "== NUMBER OF OBSERVATIONS =================================="
    WRITE(6,'(14A11)') &
         & "Variable:", &
         & "U","SSU","Interior U", &
         & "V","SSV","Interior V", &
         & "T","SST","Interior T", &
         & "S","SSS","Interior S", &
         & "SSH"
    WRITE(6,'(A11,13I11)') &
         & "#OBS:", &
         & no(iv3d_u),no3ds(iv3d_u),no3di(iv3d_u), &
         & no(iv3d_v),no3ds(iv3d_v),no3di(iv3d_v), &
         & no(iv3d_t),no3ds(iv3d_t),no3di(iv3d_t), &
         & no(iv3d_s),no3ds(iv3d_s),no3di(iv3d_s), &
         & no(nv3d+iv2d_z)
    WRITE(6,'(A)') "============================================================"

  END SUBROUTINE monit_mean
  
END MODULE common
