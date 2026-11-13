MODULE common_mpi
  
  !=======================================================================
  !
  ! [PURPOSE:] General MPI procedures
  !
  ! [HISTORY:]
  !   2005/09/06 Takemasa MIYOSHI  created
  !   2025/07/31 Shun OHISHI       added bcast
  !   2026/06/04 Shun OHISHI       updated
  !
  !=======================================================================

CONTAINS

  !----------------------------------------------------------------------

  SUBROUTINE initialize_mpi

    USE common_setting
    USE MPI    
    IMPLICIT NONE

    INTEGER ierr
    INTEGER iprov   !Thread support level provided by MPI
    INTEGER color   !Color for MPI_COMM_SPLIT (0: node leader, MPI_UNDEFINED: others)
    INTEGER nnode       !Number of nodes
    INTEGER nprocs_node !Number of processes in my node

    CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,iprov,ierr)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

    IF(myrank == root)THEN
       WRITE(6,'(A,I0)') "Number of Processor:",nprocs
       IF(iprov < MPI_THREAD_FUNNELED)THEN
          WRITE(6,'(A)') "***Warning: MPI_THREAD_FUNNELED is not supported"
       END IF
    END IF

    IF(MOD(nprocs,nmem) /= 0 .and. MOD(nmem,nprocs) /= 0)THEN
       IF(myrank == root) WRITE(6,'(A,2I8)') &
            & "***Error: nprocs must be a multiple or divisor of nmem: ", &
            & nprocs,nmem
       CALL MPI_Abort(MPI_COMM_WORLD,1,ierr)
    END IF

    root_out=nprocs-1

    !---MPI(r_size)
    IF(r_size == r_dble)THEN
       MPI_R_SIZE=MPI_DOUBLE_PRECISION
    ELSE IF(r_size == r_sngl)THEN
       MPI_R_SIZE=MPI_REAL
    END IF

    !---Node communicator (shared memory) and node-leader communicator
    CALL MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,comm_node,ierr)
    CALL MPI_COMM_RANK(comm_node,node_rank,ierr)
    color=MPI_UNDEFINED
    IF(node_rank == 0) color=0
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,color,myrank,comm_leader,ierr)

    !---Node information
    CALL MPI_COMM_SIZE(comm_node,nprocs_node,ierr)
    nnode=0
    IF(node_rank == 0) nnode=1
    CALL MPI_ALLREDUCE(MPI_IN_PLACE,nnode,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    IF(myrank == root) WRITE(6,'(A,I0,A,I0)') "Number of Nodes:",nnode,", Processes per node:",nprocs_node

  END SUBROUTINE initialize_mpi

  !-----------------------------------------------------------------------
  ! Scatter gridded data to processes (nrank -> all)
  !-----------------------------------------------------------------------

  SUBROUTINE scatter_grd_mpi(nrank,v3dg,v2dg,v3d,v2d)

    !$ USE OMP_LIB
    USE common_setting
    USE MPI
    IMPLICIT NONE

    !---Common
    INTEGER k
    INTEGER n
    INTEGER ivd
    INTEGER ierr

    REAL(r_sngl),ALLOCATABLE :: bufs(:,:,:) !Send buffer (nrank only) [nij1max*nlev*nprocs]
    REAL(r_sngl),ALLOCATABLE :: bufr(:,:)   !Recv buffer [nij1max*nlev]

    !---IN
    INTEGER,INTENT(IN) :: nrank
    REAL(r_size),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
    REAL(r_size),INTENT(IN) :: v2dg(nlon,nlat,nv2d)

    !---OUT
    REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
    REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)

    !*** Buffers are allocated for one 3D variable (nlev levels) and reused
    !*** for each variable to reduce memory (1/nv3d of the full buffer).

    ALLOCATE(bufr(nij1max,nlev))
    bufr(:,:) = 0.e0

    IF(myrank == nrank)THEN
       ALLOCATE(bufs(nij1max,nlev,nprocs))
       bufs(:,:,:)=0.e0
    ELSE
       ALLOCATE(bufs(1,1,1))
    END IF

    !---3D variables (one variable per MPI_SCATTER)
    n = nij1max * nlev
    DO ivd=1,nv3d

       !v3dg --> bufs
       IF(myrank == nrank)THEN
          DO k=1,nlev
             CALL grd_to_buf(v3dg(:,:,k,ivd),bufs(:,k,:))
          END DO
       END IF

       !bufs --> bufr
       CALL MPI_SCATTER(bufs,n,MPI_REAL,&
            & bufr,n,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

       !bufr --> v3d
       !$OMP PARALLEL DO PRIVATE(k)
       DO k=1,nlev
          v3d(:,k,ivd) = REAL(bufr(1:nij1,k),r_size)
       END DO
       !$OMP END PARALLEL DO

    END DO

    !---2D variables (bufs(:,1:nv2d,:) is packed as [nij1max*nv2d] per process)
    DEALLOCATE(bufs,bufr)
    ALLOCATE(bufr(nij1max,nv2d))
    IF(myrank == nrank)THEN
       ALLOCATE(bufs(nij1max,nv2d,nprocs))
       bufs(:,:,:)=0.e0
       DO ivd=1,nv2d
          CALL grd_to_buf(v2dg(:,:,ivd),bufs(:,ivd,:))
       END DO
    ELSE
       ALLOCATE(bufs(1,1,1))
    END IF

    n = nij1max * nv2d
    CALL MPI_SCATTER(bufs,n,MPI_REAL,&
         & bufr,n,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

    DO ivd=1,nv2d
       v2d(:,ivd) = REAL(bufr(1:nij1,ivd),r_size)
    END DO

    DEALLOCATE(bufs,bufr)

  END SUBROUTINE scatter_grd_mpi

  !-----------------------------------------------------------------------
  ! Gather gridded data (all -> nrank)
  !-----------------------------------------------------------------------

  SUBROUTINE gather_grd_mpi(nrank,v3d,v2d,v3dg,v2dg)

    !$ USE OMP_LIB
    USE common_setting
    USE MPI
    IMPLICIT NONE

    !---Common
    INTEGER k
    INTEGER n
    INTEGER ivd
    INTEGER ierr

    REAL(r_sngl),ALLOCATABLE :: bufs(:,:)   !Send buffer [nij1max*nlev]
    REAL(r_sngl),ALLOCATABLE :: bufr(:,:,:) !Recv buffer (nrank only) [nij1max*nlev*nprocs]

    !---IN
    INTEGER,INTENT(IN) :: nrank
    REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
    REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)

    !---OUT
    REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
    REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)

    !*** Buffers are allocated for one 3D variable (nlev levels) and reused
    !*** for each variable to reduce memory (1/nv3d of the full buffer).

    ALLOCATE(bufs(nij1max,nlev))
    bufs(:,:) = 0.e0
    IF(myrank == nrank)THEN
       ALLOCATE(bufr(nij1max,nlev,nprocs))
    ELSE
       ALLOCATE(bufr(1,1,1))
    END IF

    !---3D variables (one variable per MPI_GATHER)
    n = nij1max * nlev
    DO ivd=1,nv3d

       !v3d --> bufs
       !$OMP PARALLEL DO PRIVATE(k)
       DO k=1,nlev
          bufs(1:nij1,k) = REAL(v3d(:,k,ivd),r_sngl)
       END DO
       !$OMP END PARALLEL DO

       !bufs --> bufr
       CALL MPI_GATHER(bufs,n,MPI_REAL,&
            & bufr,n,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

       !bufr --> v3dg
       IF(myrank == nrank)THEN
          DO k=1,nlev
             CALL buf_to_grd(bufr(:,k,:),v3dg(:,:,k,ivd))
          END DO
       END IF

    END DO

    !---2D variables
    DEALLOCATE(bufs,bufr)
    ALLOCATE(bufs(nij1max,nv2d))
    bufs(:,:) = 0.e0
    IF(myrank == nrank)THEN
       ALLOCATE(bufr(nij1max,nv2d,nprocs))
    ELSE
       ALLOCATE(bufr(1,1,1))
    END IF

    DO ivd=1,nv2d
       bufs(1:nij1,ivd) = REAL(v2d(:,ivd),r_sngl)
    END DO

    n = nij1max * nv2d
    CALL MPI_GATHER(bufs,n,MPI_REAL,&
         & bufr,n,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

    IF(myrank == nrank)THEN
       DO ivd=1,nv2d
          CALL buf_to_grd(bufr(:,ivd,:),v2dg(:,:,ivd))
       END DO
    END IF

    DEALLOCATE(bufs,bufr)

  END SUBROUTINE gather_grd_mpi

  !-----------------------------------------------------------------------
  ! Scatter gridded data of all members in block ib at once (owner -> all)
  !
  !   Member imem in block ib is owned by rank  imem-1-(ib-1)*nprocs
  !   (= the rank that read the member file in read_ens_mpi).
  !   Each owner sends its member to all ranks by one MPI_ALLTOALLV
  !   (all owners send simultaneously).
  !
  !   send : owner only, nsize to every rank
  !   recv : all ranks , nsize from each owner (rank 0 ~ nmem_blk-1)
  !-----------------------------------------------------------------------
  SUBROUTINE scatter_grd_mpi_all(ib,v3dg,v2dg,v3d,v2d)

    !$ USE OMP_LIB
    USE common_setting
    USE MPI
    IMPLICIT NONE

    !---Common
    INTEGER k,ivd,ierr
    INTEGER iprocs
    INTEGER imem_own        !Member owned by myrank
    INTEGER imem_blk        !Member index within block (1 ~ nmem_blk)
    INTEGER imem_offset     !Member offset of block ib
    INTEGER nmem_blk        !#Members in block ib
    INTEGER nsize3d,nsize2d !#Elements per rank (3D: one variable, 2D: all variables)
    LOGICAL lowner          !.TRUE.: myrank owns a member in block ib

    INTEGER,ALLOCATABLE :: send_cnt1(:),send_dsp1(:) !Send count/displacement (per element)
    INTEGER,ALLOCATABLE :: recv_cnt1(:),recv_dsp1(:) !Recv count/displacement (per element)

    REAL(r_sngl),ALLOCATABLE :: bufs(:,:,:) !Send buffer (owner only) [nij1max,nlev or nv2d,nprocs]
    REAL(r_sngl),ALLOCATABLE :: bufr(:,:,:) !Recv buffer              [nij1max,nlev or nv2d,nmem_blk]

    !---IN
    INTEGER,INTENT(IN) :: ib
    REAL(r_size),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d) !Global data (owner only)
    REAL(r_size),INTENT(IN) :: v2dg(nlon,nlat,nv2d)

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: v3d(nij1,nlev,nmem,nv3d) !Local data
    REAL(r_size),INTENT(INOUT) :: v2d(nij1,nmem,nv2d)

    !=====================================================================
    ! (1) Owner and number of members in block ib
    !=====================================================================

    imem_offset = (ib-1)*nprocs
    imem_own    = myrank+1 + imem_offset
    lowner      = (1 <= imem_own .and. imem_own <= nmem)
    nmem_blk    = MIN(nmem-imem_offset,nprocs)

    nsize3d = nij1max*nlev
    nsize2d = nij1max*nv2d

    !=====================================================================
    ! (2) Send/Recv count and displacement (per element; multiplied by nsize)
    !=====================================================================

    ALLOCATE(send_cnt1(nprocs),send_dsp1(nprocs))
    ALLOCATE(recv_cnt1(nprocs),recv_dsp1(nprocs))

    send_cnt1(:)=0
    send_dsp1(:)=0
    recv_cnt1(:)=0
    recv_dsp1(:)=0

    DO iprocs=0,nprocs-1

       !Send: owner -> every rank
       IF(lowner)THEN
          send_cnt1(iprocs+1) = 1
          send_dsp1(iprocs+1) = iprocs
       END IF

       !Recv: from owner rank (0 ~ nmem_blk-1)
       IF(iprocs < nmem_blk)THEN
          recv_cnt1(iprocs+1) = 1
          recv_dsp1(iprocs+1) = iprocs
       END IF

    END DO

    !=====================================================================
    ! (3) 3D variables (one variable per MPI_ALLTOALLV)
    !=====================================================================

    IF(lowner)THEN
       ALLOCATE(bufs(nij1max,nlev,nprocs))
    ELSE
       ALLOCATE(bufs(1,1,1))
    END IF
    ALLOCATE(bufr(nij1max,nlev,nmem_blk))

    DO ivd=1,nv3d

       !Global --> send buffer (owner only)
       IF(lowner)THEN
          DO k=1,nlev
             CALL grd_to_buf(v3dg(:,:,k,ivd),bufs(:,k,:))
          END DO
       END IF

       !Owner --> all ranks
       CALL MPI_ALLTOALLV(bufs,send_cnt1*nsize3d,send_dsp1*nsize3d,MPI_REAL, &
            &             bufr,recv_cnt1*nsize3d,recv_dsp1*nsize3d,MPI_REAL, &
            &             MPI_COMM_WORLD,ierr)

       !Recv buffer --> local data
       !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(imem_blk,k)
       DO imem_blk=1,nmem_blk
          DO k=1,nlev
             v3d(:,k,imem_offset+imem_blk,ivd) = REAL(bufr(1:nij1,k,imem_blk),r_size)
          END DO
       END DO
       !$OMP END PARALLEL DO

    END DO

    DEALLOCATE(bufs,bufr)

    !=====================================================================
    ! (4) 2D variables (all variables in one MPI_ALLTOALLV)
    !=====================================================================

    IF(lowner)THEN
       ALLOCATE(bufs(nij1max,nv2d,nprocs))
       DO ivd=1,nv2d
          CALL grd_to_buf(v2dg(:,:,ivd),bufs(:,ivd,:))
       END DO
    ELSE
       ALLOCATE(bufs(1,1,1))
    END IF
    ALLOCATE(bufr(nij1max,nv2d,nmem_blk))

    CALL MPI_ALLTOALLV(bufs,send_cnt1*nsize2d,send_dsp1*nsize2d,MPI_REAL, &
         &             bufr,recv_cnt1*nsize2d,recv_dsp1*nsize2d,MPI_REAL, &
         &             MPI_COMM_WORLD,ierr)

    DO imem_blk=1,nmem_blk
       DO ivd=1,nv2d
          v2d(:,imem_offset+imem_blk,ivd) = REAL(bufr(1:nij1,ivd,imem_blk),r_size)
       END DO
    END DO

    DEALLOCATE(bufs,bufr)
    DEALLOCATE(send_cnt1,send_dsp1,recv_cnt1,recv_dsp1)

  END SUBROUTINE scatter_grd_mpi_all

  !-----------------------------------------------------------------------
  ! Gather gridded data of all tasks in block i at once (all -> owner)
  !
  !   Task itask in block i is owned by rank  itask-1-(i-1)*nprocs
  !     itask = [1,nmem]: analysis member
  !             nmem+1  : analysis mean,  nmem+2: analysis spread
  !             nmem+3  : forecast mean,  nmem+4: forecast spread
  !   All ranks send their local data of all tasks by one MPI_ALLTOALLV.
  !
  !   send : all ranks , nsize to each owner (rank 0 ~ ntask_blk-1)
  !   recv : owner only, nsize from every rank
  !-----------------------------------------------------------------------
  SUBROUTINE gather_grd_mpi_all(i,ntask,v3da,v2da,v3dam,v2dam,v3das,v2das, &
       &                        v3dfm,v2dfm,v3dfs,v2dfs,v3dg,v2dg)

    !$ USE OMP_LIB
    USE common_setting
    USE MPI
    IMPLICIT NONE

    !---Common
    INTEGER k,ivd,ierr
    INTEGER iprocs
    INTEGER itask           !Task index (1 ~ ntask)
    INTEGER itask_own       !Task owned by myrank
    INTEGER itask_blk       !Task index within block (1 ~ ntask_blk)
    INTEGER itask_offset    !Task offset of block i
    INTEGER ntask_blk       !#Tasks in block i
    INTEGER nsize3d,nsize2d !#Elements per rank (3D: one variable, 2D: all variables)
    LOGICAL lowner          !.TRUE.: myrank owns a task in block i

    INTEGER,ALLOCATABLE :: send_cnt1(:),send_dsp1(:) !Send count/displacement (per element)
    INTEGER,ALLOCATABLE :: recv_cnt1(:),recv_dsp1(:) !Recv count/displacement (per element)

    REAL(r_sngl),ALLOCATABLE :: bufs(:,:,:) !Send buffer              [nij1max,nlev or nv2d,ntask_blk]
    REAL(r_sngl),ALLOCATABLE :: bufr(:,:,:) !Recv buffer (owner only) [nij1max,nlev or nv2d,nprocs]

    !---IN
    INTEGER,INTENT(IN) :: i,ntask
    REAL(r_size),INTENT(IN) :: v3da(nij1,nlev,nmem,nv3d),v2da(nij1,nmem,nv2d) !Analysis ensemble
    REAL(r_size),INTENT(IN) :: v3dam(nij1,nlev,nv3d),v2dam(nij1,nv2d)         !Analysis mean
    REAL(r_size),INTENT(IN) :: v3das(nij1,nlev,nv3d),v2das(nij1,nv2d)         !Analysis spread
    REAL(r_size),INTENT(IN) :: v3dfm(nij1,nlev,nv3d),v2dfm(nij1,nv2d)         !Forecast mean
    REAL(r_size),INTENT(IN) :: v3dfs(nij1,nlev,nv3d),v2dfs(nij1,nv2d)         !Forecast spread

    !---INOUT (updated on owner only)
    REAL(r_sngl),INTENT(INOUT) :: v3dg(nlon,nlat,nlev,nv3d)
    REAL(r_sngl),INTENT(INOUT) :: v2dg(nlon,nlat,nv2d)

    !=====================================================================
    ! (1) Owner and number of tasks in block i
    !=====================================================================

    itask_offset = (i-1)*nprocs
    itask_own    = myrank+1 + itask_offset
    lowner       = (1 <= itask_own .and. itask_own <= ntask)
    ntask_blk    = MIN(ntask-itask_offset,nprocs)

    nsize3d = nij1max*nlev
    nsize2d = nij1max*nv2d

    !=====================================================================
    ! (2) Send/Recv count and displacement (per element; multiplied by nsize)
    !=====================================================================

    ALLOCATE(send_cnt1(nprocs),send_dsp1(nprocs))
    ALLOCATE(recv_cnt1(nprocs),recv_dsp1(nprocs))

    send_cnt1(:)=0
    send_dsp1(:)=0
    recv_cnt1(:)=0
    recv_dsp1(:)=0

    DO iprocs=0,nprocs-1

       !Send: every rank -> owner rank (0 ~ ntask_blk-1)
       IF(iprocs < ntask_blk)THEN
          send_cnt1(iprocs+1) = 1
          send_dsp1(iprocs+1) = iprocs
       END IF

       !Recv: owner <- every rank
       IF(lowner)THEN
          recv_cnt1(iprocs+1) = 1
          recv_dsp1(iprocs+1) = iprocs
       END IF

    END DO

    !=====================================================================
    ! (3) 3D variables (one variable per MPI_ALLTOALLV)
    !=====================================================================

    ALLOCATE(bufs(nij1max,nlev,ntask_blk))
    bufs(:,:,:)=0.e0 !Padding (nij1 < nij1max)
    IF(lowner)THEN
       ALLOCATE(bufr(nij1max,nlev,nprocs))
    ELSE
       ALLOCATE(bufr(1,1,1))
    END IF

    DO ivd=1,nv3d

       !Local data of all tasks --> send buffer
       !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(itask_blk,k,itask)
       DO itask_blk=1,ntask_blk
          DO k=1,nlev

             itask = itask_offset+itask_blk

             IF(itask <= nmem)THEN
                bufs(1:nij1,k,itask_blk) = REAL(v3da(:,k,itask,ivd),r_sngl)
             ELSE IF(itask == nmem+1)THEN
                bufs(1:nij1,k,itask_blk) = REAL(v3dam(:,k,ivd),r_sngl)
             ELSE IF(itask == nmem+2)THEN
                bufs(1:nij1,k,itask_blk) = REAL(v3das(:,k,ivd),r_sngl)
             ELSE IF(itask == nmem+3)THEN
                bufs(1:nij1,k,itask_blk) = REAL(v3dfm(:,k,ivd),r_sngl)
             ELSE IF(itask == nmem+4)THEN
                bufs(1:nij1,k,itask_blk) = REAL(v3dfs(:,k,ivd),r_sngl)
             END IF

          END DO
       END DO
       !$OMP END PARALLEL DO

       !All ranks --> owner
       CALL MPI_ALLTOALLV(bufs,send_cnt1*nsize3d,send_dsp1*nsize3d,MPI_REAL, &
            &             bufr,recv_cnt1*nsize3d,recv_dsp1*nsize3d,MPI_REAL, &
            &             MPI_COMM_WORLD,ierr)

       !Recv buffer --> global data (owner only)
       IF(lowner)THEN
          DO k=1,nlev
             CALL buf_to_grd(bufr(:,k,:),v3dg(:,:,k,ivd))
          END DO
       END IF

    END DO

    DEALLOCATE(bufs,bufr)

    !=====================================================================
    ! (4) 2D variables (all variables in one MPI_ALLTOALLV)
    !=====================================================================

    ALLOCATE(bufs(nij1max,nv2d,ntask_blk))
    bufs(:,:,:)=0.e0 !Padding (nij1 < nij1max)
    IF(lowner)THEN
       ALLOCATE(bufr(nij1max,nv2d,nprocs))
    ELSE
       ALLOCATE(bufr(1,1,1))
    END IF

    DO itask_blk=1,ntask_blk

       itask = itask_offset+itask_blk

       DO ivd=1,nv2d
          IF(itask <= nmem)THEN
             bufs(1:nij1,ivd,itask_blk) = REAL(v2da(:,itask,ivd),r_sngl)
          ELSE IF(itask == nmem+1)THEN
             bufs(1:nij1,ivd,itask_blk) = REAL(v2dam(:,ivd),r_sngl)
          ELSE IF(itask == nmem+2)THEN
             bufs(1:nij1,ivd,itask_blk) = REAL(v2das(:,ivd),r_sngl)
          ELSE IF(itask == nmem+3)THEN
             bufs(1:nij1,ivd,itask_blk) = REAL(v2dfm(:,ivd),r_sngl)
          ELSE IF(itask == nmem+4)THEN
             bufs(1:nij1,ivd,itask_blk) = REAL(v2dfs(:,ivd),r_sngl)
          END IF
       END DO

    END DO

    CALL MPI_ALLTOALLV(bufs,send_cnt1*nsize2d,send_dsp1*nsize2d,MPI_REAL, &
         &             bufr,recv_cnt1*nsize2d,recv_dsp1*nsize2d,MPI_REAL, &
         &             MPI_COMM_WORLD,ierr)

    IF(lowner)THEN
       DO ivd=1,nv2d
          CALL buf_to_grd(bufr(:,ivd,:),v2dg(:,:,ivd))
       END DO
    END IF

    DEALLOCATE(bufs,bufr)
    DEALLOCATE(send_cnt1,send_dsp1,recv_cnt1,recv_dsp1)

  END SUBROUTINE gather_grd_mpi_all

  !-----------------------------------------------------------------------
  ! Share nrank data with all processer
  !-----------------------------------------------------------------------

  SUBROUTINE bcast_mpi_1d(nrank,im,dat)

    USE common_setting
    USE MPI
    IMPLICIT NONE

    !---Common
    INTEGER ierr

    !---IN
    INTEGER,INTENT(IN) :: nrank
    INTEGER,INTENT(IN) :: im

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: dat(im)

    CALL MPI_BCAST(dat,im,MPI_R_SIZE,nrank,MPI_COMM_WORLD,ierr)

  END SUBROUTINE bcast_mpi_1d

  !--------------------------

  SUBROUTINE bcast_mpi_2d(nrank,im,jm,dat)

    USE common_setting
    USE MPI
    IMPLICIT NONE

    !---Common
    INTEGER ierr

    !---IN
    INTEGER,INTENT(IN) :: nrank
    INTEGER,INTENT(IN) :: im,jm

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: dat(im,jm)

    CALL MPI_BCAST(dat,im*jm,MPI_R_SIZE,nrank,MPI_COMM_WORLD,ierr)

  END SUBROUTINE bcast_mpi_2d

  !----------------------------

  SUBROUTINE bcast_mpi_3d(nrank,im,jm,km,dat)

    USE common_setting
    USE MPI
    IMPLICIT NONE

    !---Common
    INTEGER ierr

    !---IN
    INTEGER,INTENT(IN) :: nrank
    INTEGER,INTENT(IN) :: im,jm,km

    !---INOUT
    REAL(r_size),INTENT(INOUT) :: dat(im,jm,km)

    CALL MPI_BCAST(dat,im*jm*km,MPI_R_SIZE,nrank,MPI_COMM_WORLD,ierr)

  END SUBROUTINE bcast_mpi_3d

  !-----------------------------------------------------------------------
  ! gridded data -> buffer
  !-----------------------------------------------------------------------
  SUBROUTINE grd_to_buf(grd,buf)

    !$ USE OMP_LIB
    USE common_setting
    IMPLICIT NONE

    INTEGER i,j,m,ilon,ilat

    !IN
    REAL(r_size),INTENT(IN) :: grd(nlon,nlat)

    !OUT
    REAL(r_sngl),INTENT(OUT) :: buf(nij1max,nprocs)

    !$OMP PARALLEL DO PRIVATE(m,i,j,ilon,ilat)    
    DO m=1,nprocs
       DO i=1,nij1node(m)
          j = m-1 + nprocs * (i-1)
          ilon = MOD(j,nlon) + 1
          ilat = (j-ilon+1) / nlon + 1
          buf(i,m) = REAL(grd(ilon,ilat),r_sngl)
       END DO
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE grd_to_buf

  !-----------------------------------------------------------------------
  ! buffer -> gridded data
  !-----------------------------------------------------------------------
  SUBROUTINE buf_to_grd(buf,grd)

    !$ USE OMP_LIB
    USE common_setting
    IMPLICIT NONE

    INTEGER i,j,m,ilon,ilat

    !IN
    REAL(r_sngl),INTENT(IN) :: buf(nij1max,nprocs)

    !OUT
    REAL(r_sngl),INTENT(OUT) :: grd(nlon,nlat)

    !$OMP PARALLEL DO PRIVATE(m,i,j,ilon,ilat)    
    DO m=1,nprocs
       DO i=1,nij1node(m)
          j = m-1 + nprocs * (i-1)
          ilon = MOD(j,nlon) + 1
          ilat = (j-ilon+1) / nlon + 1
          grd(ilon,ilat) = buf(i,m)
       END DO
    END DO
    !$OMP END PARALLEL DO

  END SUBROUTINE buf_to_grd

  !-----------------------------------------------------------------------
  ! Allocate dYf in shared memory (one copy per node), transposed layout
  !
  !   hdxfT(imem,iobs) = hdxfT_flat( imem + (iobs-1)*nmem )
  !
  ! NOTE:
  !   Fujitsu compiler may overflow 32-bit address offset (> 2 GiB) in
  !   strided access of 2-D pointer (e.g., hdxfT(imem,:) over all iobs).
  !   - hdxfT(:,iobs)            : contiguous access only (e.g., set_obs_local)
  !   - hdxfT_flat (INTEGER(8))  : strided access over iobs
  !                                (read_ens_mpi, ensemble_mesp_obs, output)
  !-----------------------------------------------------------------------

  SUBROUTINE allocate_shared_hdxfT

    !$ USE OMP_LIB
    USE,INTRINSIC :: ISO_C_BINDING, only: C_PTR, C_F_POINTER
    USE common_setting
    USE MPI
    IMPLICIT NONE

    !---Parameter
    INTEGER,PARAMETER :: nblk=4096 !Block size of obs. for first touch

    !---Common
    INTEGER ierr,disp_unit,nbyte
    INTEGER nprocs_node
    INTEGER is,ie
    INTEGER(8) is8,ie8
    INTEGER(KIND=MPI_ADDRESS_KIND) winsize
    TYPE(C_PTR) baseptr

    nbyte=STORAGE_SIZE(REAL(0.d0,r_size))/8

    !---Memory is allocated only by node leader (others: size 0)
    winsize=0
    IF(node_rank == 0) winsize=INT(MAX(nobs,1),MPI_ADDRESS_KIND)*INT(nmem,MPI_ADDRESS_KIND)*nbyte

    CALL MPI_WIN_ALLOCATE_SHARED(winsize,nbyte,MPI_INFO_NULL,comm_node,baseptr,win_hdxfT,ierr)
    CALL MPI_WIN_SHARED_QUERY(win_hdxfT,0,winsize,disp_unit,baseptr,ierr) !Address of node leader

    !---2-D view (nmem,nobs) and 1-D view (64-bit index) of the same memory
    CALL C_F_POINTER(baseptr,hdxfT,(/nmem,MAX(nobs,1)/))
    CALL C_F_POINTER(baseptr,hdxfT_flat,(/INT(MAX(nobs,1),8)*INT(nmem,8)/))

    !---Initialization: first touch distributed over node ranks and threads
    !   (block of obs.; spread pages over CMGs)
    CALL MPI_COMM_SIZE(comm_node,nprocs_node,ierr)
    CALL MPI_WIN_FENCE(0,win_hdxfT,ierr)

    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(is,ie,is8,ie8)
    DO is=node_rank*nblk+1,nobs,nprocs_node*nblk
       ie  = MIN(is+nblk-1,nobs)
       is8 = INT(is-1,8)*INT(nmem,8) + 1_8 !hdxfT(1,is)
       ie8 = INT(ie,8)*INT(nmem,8)         !hdxfT(nmem,ie)
       hdxfT_flat(is8:ie8) = REAL(0.d0,r_size)
    END DO
    !$OMP END PARALLEL DO

    CALL MPI_WIN_FENCE(0,win_hdxfT,ierr)

  END SUBROUTINE allocate_shared_hdxfT

  !-----------------------------------------------------------------------
  ! Deallocate dYf in shared memory
  !-----------------------------------------------------------------------
  SUBROUTINE deallocate_shared_hdxfT

    USE common_setting
    USE MPI
    IMPLICIT NONE

    INTEGER ierr

    NULLIFY(hdxfT)
    NULLIFY(hdxfT_flat)
    CALL MPI_WIN_FREE(win_hdxfT,ierr)

  END SUBROUTINE deallocate_shared_hdxfT

END MODULE common_mpi
