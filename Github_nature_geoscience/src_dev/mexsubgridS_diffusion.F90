# include <fintrf.h>
!======================================================================
#if 0
!     Input
!     (MXN, MYN, MCXN, MCYN, MDX, MDY, MCDX, MCDY, META, MMU, MSXX, MSXY, sxx1, sxy1, dsxx, dsxy, dsubgrids, timestep)
!     Output
!     structure G
#endif
!
subroutine mexFunction(nlhs,plhs,nrhs,prhs)

    ! Declarations
    implicit none

    !     mexFunction arguments:
    mwPointer :: plhs(*), prhs(*)
    integer :: nlhs, nrhs

    !     Function declarations:
    mwPointer :: mxGetPr
    mwPointer :: mxCreateDoubleMatrix
    integer :: mexPrintf
    mwPointer :: mxGetM, mxGetN

    !     Pointers to input/output mxArrays:
    mwPointer :: inmarker_ptr, ingrid_ptr

    !     Arguments for computational routine:
    double precision, allocatable, dimension(:) :: MDX, MDY, MCDX, MCDY, MSXX, MSXY, META, MMU
    integer, allocatable, dimension (:) :: MXN, MYN, MCXN, MCYN
    double precision, allocatable, dimension(:) :: tmpD

    double precision, allocatable, dimension(:,:) :: sxx1, sxy1, dsxx, dsxy
    double precision, allocatable, dimension(:,:) :: dsxxn, dsxyn, wtetas, wtetan

    ! Local variables
    integer :: marknum, xnum, ynum, numel, numelc
    integer :: i, k, mm1
    integer :: xn, yn

    double precision :: dx, dy, timestep, dsubgrids, sdif
    double precision :: sxxm, sxym, dsxxm, dsxym, sdm, mwt

    character*120 line

    !-----------------------------------------------------------------------
    ! Number of markers
    marknum = mxGetM(prhs(1))*mxGetN(prhs(1))

    allocate(MXN(marknum),MYN(marknum),MDX(marknum),MDY(marknum))
    allocate(MCXN(marknum),MCYN(marknum),MCDX(marknum),MCDY(marknum))
    allocate(MSXX(marknum),MSXY(marknum),META(marknum),MMU(marknum))

    ! Populate marker arrays: MXN, MYN, MCXN, MCYN, MDX, MDY, MCDX, MCDY, META, MMU, MSXX, MSXY

    allocate(tmpD(marknum))

    inmarker_ptr = mxGetPr(prhs(1))
    call mxCopyPtrToReal8(inmarker_ptr,tmpD,marknum)
    MXN = int(tmpD)
    inmarker_ptr = mxGetPr(prhs(2))
    call mxCopyPtrToReal8(inmarker_ptr,tmpD,marknum)
    MYN = int(tmpD)

    inmarker_ptr = mxGetPr(prhs(3))
    call mxCopyPtrToReal8(inmarker_ptr,tmpD,marknum)
    MCXN = int(tmpD)
    inmarker_ptr = mxGetPr(prhs(4))
    call mxCopyPtrToReal8(inmarker_ptr,tmpD,marknum)
    MCYN = int(tmpD)

    inmarker_ptr = mxGetPr(prhs(5))
    call mxCopyPtrToReal8(inmarker_ptr,MDX,marknum)
    inmarker_ptr = mxGetPr(prhs(6))
    call mxCopyPtrToReal8(inmarker_ptr,MDY,marknum)

    inmarker_ptr = mxGetPr(prhs(7))
    call mxCopyPtrToReal8(inmarker_ptr,MCDX,marknum)
    inmarker_ptr = mxGetPr(prhs(8))
    call mxCopyPtrToReal8(inmarker_ptr,MCDY,marknum)

    inmarker_ptr = mxGetPr(prhs(9))
    call mxCopyPtrToReal8(inmarker_ptr,META,marknum)
    inmarker_ptr = mxGetPr(prhs(10))
    call mxCopyPtrToReal8(inmarker_ptr,MMU,marknum)
    inmarker_ptr = mxGetPr(prhs(11))
    call mxCopyPtrToReal8(inmarker_ptr,MSXX,marknum)
    inmarker_ptr = mxGetPr(prhs(12))
    call mxCopyPtrToReal8(inmarker_ptr,MSXY,marknum)

        ! Populate grid arrays: sxx1, sxy1, dsxx, dsxy
    ynum = mxGetM(prhs(14))
    xnum = mxGetN(prhs(14))

    allocate(sxy1(ynum,xnum),dsxy(ynum,xnum))
    allocate(sxx1(ynum-1,xnum-1),dsxx(ynum-1,xnum-1))

    numel = ynum*xnum
    numelc = (ynum-1)*(xnum-1)

    ingrid_ptr = mxGetPr(prhs(13))
    call mxCopyPtrToReal8(ingrid_ptr,sxx1,numelc)
    ingrid_ptr = mxGetPr(prhs(14))
    call mxCopyPtrToReal8(ingrid_ptr,sxy1,numel)
    ingrid_ptr = mxGetPr(prhs(15))
    call mxCopyPtrToReal8(ingrid_ptr,dsxx,numelc)
    ingrid_ptr = mxGetPr(prhs(16))
    call mxCopyPtrToReal8(ingrid_ptr,dsxy,numel)

    ingrid_ptr = mxGetPr(prhs(17))
    call mxCopyPtrToReal8(ingrid_ptr,dsubgrids,1)
    ingrid_ptr = mxGetPr(prhs(18))
    call mxCopyPtrToReal8(ingrid_ptr,timestep,1)

    allocate(dsxyn(ynum,xnum),dsxxn(ynum-1,xnum-1))
    allocate(wtetas(ynum,xnum),wtetan(ynum-1,xnum-1))
    dsxyn = 0.0d0
    dsxxn = 0.0d0
    wtetas = 0.0d0
    wtetan = 0.0d0

    ! Marker cycle
    do mm1=1,marknum,1

        ! Compute local stress relaxation timescale (Maxwell time) for the marker
        sdm = META(mm1)/MMU(mm1)
        ! Computing degree of subgrid stress relaxation
        sdif=-dsubgrids*timestep/sdm
        if (sdif<-30.0d0) then
            sdif=-30.0d0
        end if
        sdif=(1-dexp(sdif))

        ! yn sxy(yn,xn)--------------------sxy(yn,xn+1)
        ! ? ^ ?
        ! ? ? ?
        ! ? dy ?
        ! ? ? ?
        ! ? v ?
        ! ?<----dx--->o MSXY(mm1) ?
        ! ? ?
        ! ? ?
        ! yn+1 sxy(yn+1,xn)-------------------sxy(yn+1,xn+1)
        !
        !
        ! Interpolating old shear stress from Sxy nodes
        !
        ! Define indexes for upper left node in the cell where the marker is
        xn=MXN(mm1)
        yn=MYN(mm1)

        ! Define normalized distances from marker to the upper left node;
        dx=MDX(mm1)
        dy=MDY(mm1)

        ! Compute marker weight koefficient from cell dimensions
        ! Number of markers in a cell is in invert proportion to the cell volume
        mwt=1 !/xstp1(xn)/ystp1(yn);

        ! Interpolate old Sxy stress for the marker
        sxym=0.0d0
        sxym=sxym+(1.0-dx)*(1.0-dy)*sxy1(yn,xn)
        sxym=sxym+(1.0-dx)*dy*sxy1(yn+1,xn)
        sxym=sxym+dx*(1.0-dy)*sxy1(yn,xn+1)
        sxym=sxym+dx*dy*sxy1(yn+1,xn+1)
        ! Calculate Nodal-Marker subgrid Sxy stress difference
        dsxym=sxym-MSXY(mm1)
        ! Relaxing Nodal-Marker subgrid Sxy stress difference
        dsxym=dsxym*sdif

        ! Correcting old stress for the marker
        MSXY(mm1)=MSXY(mm1)+dsxym

        ! Interpolating subgrid Sxy stress changes to 4 nodes
        ! only using markers located at <=0.5 gridstep distances from nodes
        if (dx<=0.5 .and. dy<=0.5) then
            dsxyn(yn,xn)=dsxyn(yn,xn)+(1.0-dx)*(1.0-dy)*dsxym*mwt
            wtetas(yn,xn)=wtetas(yn,xn)+(1.0-dx)*(1.0-dy)*mwt
        end if
        if (dx<=0.5 .and. dy>=0.5) then
            dsxyn(yn+1,xn)=dsxyn(yn+1,xn)+(1.0-dx)*dy*dsxym*mwt
            wtetas(yn+1,xn)=wtetas(yn+1,xn)+(1.0-dx)*dy*mwt
        end if
        if (dx>=0.5 .and. dy<=0.5) then
            dsxyn(yn,xn+1)=dsxyn(yn,xn+1)+dx*(1.0-dy)*dsxym*mwt
            wtetas(yn,xn+1)=wtetas(yn,xn+1)+dx*(1.0-dy)*mwt
        end if
        if (dx>=0.5 .and. dy>=0.5) then
            dsxyn(yn+1,xn+1)=dsxyn(yn+1,xn+1)+dx*dy*dsxym*mwt
            wtetas(yn+1,xn+1)=wtetas(yn+1,xn+1)+dx*dy*mwt
        end if

        ! Computing marker weight for the center of current
        ! basic cell where Sxx stress is located
        mwt=mwt*(1.0d0-dabs(0.5d0-dx))*(1.0d0-dabs(0.5d0-dy))
        ! yn sxx(yn,xn)--------------------sxx(yn,xn+1)
        ! ? ^ ?
        ! ? ? ?
        ! ? dy ?
        ! ? ? ?
        ! ? v ?
        ! ?<----dx--->o MSXX(mm1) ?
        ! ? ?
        ! ? ?
        ! yn+1 sxx(yn+1,xn)-------------------sxx(yn+1,xn+1)
        !
        !
        ! Interpolating old normal stress from Sxx nodes
        !
        ! Define, check indexes for upper left node in the Sxx cell where the marker is
        xn = MCXN(mm1)
        yn = MCYN(mm1)
        ! Define normalized distances from marker to the upper left node;
        dx = MCDX(mm1)
        dy = MCDY(mm1)

        ! Interpolate old Sxx stress for the marker
        sxxm=0.d0
        sxxm=sxxm+(1.0-dx)*(1.0-dy)*sxx1(yn,xn)
        sxxm=sxxm+(1.0-dx)*dy*sxx1(yn+1,xn)
        sxxm=sxxm+dx*(1.0-dy)*sxx1(yn,xn+1)
        sxxm=sxxm+dx*dy*sxx1(yn+1,xn+1)
        ! Calculate Nodal-Marker subgrid Sxx stress difference
        dsxxm=sxxm-MSXX(mm1)
        ! Relaxing Nodal-Marker subgrid Sxx stress difference
        dsxxm=dsxxm*sdif

        ! Correcting old stress for the marker
        MSXX(mm1)=MSXX(mm1)+dsxxm

        ! Interpolating subgrid Sxx stress changes for the center of current basic cell
        xn=MXN(mm1)
        yn=MYN(mm1)
        dsxxn(yn,xn)=dsxxn(yn,xn)+dsxxm*mwt
        wtetan(yn,xn)=wtetan(yn,xn)+mwt

    end do

    ! Computing subgrid stress changes for nodes
    where (wtetas>0)
        dsxyn = dsxyn/wtetas
    elsewhere
        dsxyn = 0.0d0
    end where
    where (wtetan>0)
        dsxxn = dsxxn/wtetan
    elsewhere
        dsxxn = 0.0d0
    end where

    ! Subtracting subgrid stress change part from nodal stress changes
    dsxy=dsxy-dsxyn
    dsxx=dsxx-dsxxn

    ! Updating stress for markers
    do mm1=1,marknum,1

        ! Interpolating old shear stress changes from Sxy nodes
        !
        ! Define indexes for upper left node in the cell where the marker is
        xn=MXN(mm1)
        yn=MYN(mm1)

        ! Define normalized distances from marker to the upper left node;
        dx=MDX(mm1)
        dy=MDY(mm1)

        ! Interpolate old Sxy stress change for the marker
        dsxym=0
        dsxym=dsxym+(1.0-dx)*(1.0-dy)*dsxy(yn,xn)
        dsxym=dsxym+(1.0-dx)*dy*dsxy(yn+1,xn)
        dsxym=dsxym+dx*(1.0-dy)*dsxy(yn,xn+1)
        dsxym=dsxym+dx*dy*dsxy(yn+1,xn+1)

        ! Update stress for the marker
        MSXY(mm1)=MSXY(mm1)+dsxym

        ! Interpolating old normal stress changes from Sxx nodes
        !
        ! Define, check indexes for upper left node in the Sxx cell where the marker is
        xn = MCXN(mm1)
        yn = MCYN(mm1)

        ! Define normalized distances from marker to the upper left node;
        dx=MCDX(mm1)
        dy=MCDY(mm1)

        ! Interpolate old Sxx stress for the marker
        dsxxm=0
        dsxxm=dsxxm+(1.0-dx)*(1.0-dy)*dsxx(yn,xn)
        dsxxm=dsxxm+(1.0-dx)*dy*dsxx(yn+1,xn)
        dsxxm=dsxxm+dx*(1.0-dy)*dsxx(yn,xn+1)
        dsxxm=dsxxm+dx*dy*dsxx(yn+1,xn+1)

        ! Correcting old stress for the marker
        MSXX(mm1)=MSXX(mm1)+dsxxm

    end do

    ! Prepare for output
    plhs(1) = mxCreateDoubleMatrix(marknum,1,0)
    call mxCopyReal8ToPtr(MSXX,mxGetPr(plhs(1)),marknum)
    ! Prepare for output
    plhs(2) = mxCreateDoubleMatrix(marknum,1,0)
    call mxCopyReal8ToPtr(MSXY,mxGetPr(plhs(2)),marknum)
    return


end subroutine  mexFunction
