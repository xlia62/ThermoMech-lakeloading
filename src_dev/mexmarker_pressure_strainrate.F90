# include <fintrf.h>
!======================================================================
#if 0
!     Input
!     (MCXN, MCYN, MCDX, MCDY, MXN, MYN, MDX, MDY,sxx,sxy,exx,exy,pr1)
!
!     Output
!     structure G
#endif
!
subroutine mexFunction(nlhs,plhs,nrhs,prhs)

    ! marker_pressure_stress_strainrate(MCXN, MCYN, MCDX, MCDY, MXN, MYN, MDX, MDY,...
    !              sxx, sxy, exx, exy, pr1)

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
    double precision, allocatable, dimension(:) :: MX, MY, MDX, MDY, MCDX, MCDY
    double precision, allocatable, dimension(:) :: MEXX, MEXY, META, MGII, MMU, MPR, MRAT
    integer, allocatable, dimension (:) :: MXN, MYN, MCXN, MCYN
    double precision, allocatable, dimension(:) :: tmpD

    double precision, allocatable, dimension(:,:) :: sxx2, dsxx, sxy2, dsxy, exx, exy, pr1
    ! Local variables
    integer :: marknum, xnum, ynum, numel, numelc
    integer :: i, k, mm1
    integer :: xn, yn

    double precision :: dx, dy, timestep
    double precision :: exxm, exym, sxxm, sxym, dsxxm, dsxym, prm, eiim, eiimg

    character*120 line

    !-----------------------------------------------------------------------
    ! Number of markers
    marknum = mxGetM(prhs(1))*mxGetN(prhs(1))

    allocate(MXN(marknum),MYN(marknum),MCXN(marknum),MCYN(marknum))
    allocate(MDX(marknum),MDY(marknum),MCDX(marknum),MCDY(marknum))
    allocate(MEXX(marknum),MEXY(marknum),MGII(marknum),META(marknum),MMU(marknum),MPR(marknum),MRAT(marknum))

    ! Populate marker arrays: MCXN, MCYN, MCDX, MCDY, MXN, MYN, MDX, MDY, META, MMU
    allocate(tmpD(marknum))

    inmarker_ptr = mxGetPr(prhs(1))
    call mxCopyPtrToReal8(inmarker_ptr,tmpD,marknum)
    MCXN = int(tmpD)
    inmarker_ptr = mxGetPr(prhs(2))
    call mxCopyPtrToReal8(inmarker_ptr,tmpD,marknum)
    MCYN = int(tmpD)

    inmarker_ptr = mxGetPr(prhs(3))
    call mxCopyPtrToReal8(inmarker_ptr,MCDX,marknum)
    inmarker_ptr = mxGetPr(prhs(4))
    call mxCopyPtrToReal8(inmarker_ptr,MCDY,marknum)

    inmarker_ptr = mxGetPr(prhs(5))
    call mxCopyPtrToReal8(inmarker_ptr,tmpD,marknum)
    MXN = int(tmpD)
    inmarker_ptr = mxGetPr(prhs(6))
    call mxCopyPtrToReal8(inmarker_ptr,tmpD,marknum)
    MYN = int(tmpD)

    inmarker_ptr = mxGetPr(prhs(7))
    call mxCopyPtrToReal8(inmarker_ptr,MDX,marknum)
    inmarker_ptr = mxGetPr(prhs(8))
    call mxCopyPtrToReal8(inmarker_ptr,MDY,marknum)

    inmarker_ptr = mxGetPr(prhs(9))
    call mxCopyPtrToReal8(inmarker_ptr,META,marknum)
    inmarker_ptr = mxGetPr(prhs(10))
    call mxCopyPtrToReal8(inmarker_ptr,MMU,marknum)

    deallocate(tmpD)

    ! Populate grid arrays: sxx2, dsxx, sxy2, dsxy, exx, exy, pr1
    ynum = mxGetM(prhs(13))
    xnum = mxGetN(prhs(13))

    allocate(sxx2(ynum-1,xnum-1),dsxx(ynum-1,xnum-1),exx(ynum-1,xnum-1),pr1(ynum-1,xnum-1))
    allocate(sxy2(ynum,xnum),dsxy(ynum,xnum),exy(ynum,xnum))

    numel = ynum*xnum
    numelc = (ynum-1)*(xnum-1)

    ingrid_ptr = mxGetPr(prhs(11))
    call mxCopyPtrToReal8(ingrid_ptr,sxx2,numelc)
    ingrid_ptr = mxGetPr(prhs(12))
    call mxCopyPtrToReal8(ingrid_ptr,dsxx,numelc)

    ingrid_ptr = mxGetPr(prhs(13))
    call mxCopyPtrToReal8(ingrid_ptr,sxy2,numel)
    ingrid_ptr = mxGetPr(prhs(14))
    call mxCopyPtrToReal8(ingrid_ptr,dsxy,numel)

    ingrid_ptr = mxGetPr(prhs(15))
    call mxCopyPtrToReal8(ingrid_ptr,exx,numelc)

    ingrid_ptr = mxGetPr(prhs(16))
    call mxCopyPtrToReal8(ingrid_ptr,exy,numel)

    ingrid_ptr = mxGetPr(prhs(17))
    call mxCopyPtrToReal8(ingrid_ptr,pr1,numelc)

    call mxCopyPtrToReal8(mxGetPr(prhs(18)),timestep,1)

    ! Calculate Strain, Pressure, and strain rate marker-grid ration
    MEXX = 0.0d0
    MEXY = 0.0d0
    MPR = 0.0d0
    MRAT = 1.0d0

    do mm1=1,marknum
        ! Calculating strain rate for marker
        !
        ! Interpolating squares of EPS'xx=-EPS'yy from cell centers
        ! EPS'xx-nodes are displaced rightward and downward for 1/2 of gridsteps
        ! Horizontal EPS'xx index
        xn = MCXN(mm1)
        ! Vertical EPS'xx index
        yn = MCYN(mm1)
        ! Define and check normalized distances from marker to the upper left EPS'xx-node;
        ! dx=(MX(mm1)-gridcx(xn+1))/xstpc1(xn+1);
        ! dy=(MY(mm1)-gridcy(yn+1))/ystpc1(yn+1);
        dx = MCDX(mm1)
        dy = MCDY(mm1)
        ! Calculate and save Marker EPS'xx from four surrounding nodes
        exxm=0.0d0
        exxm=exxm+(1.0d0-dx)*(1.0d0-dy)*exx(yn,xn)
        exxm=exxm+(1.0d0-dx)*dy*exx(yn+1,xn)
        exxm=exxm+dx*(1.0d0-dy)*exx(yn,xn+1)
        exxm=exxm+dx*dy*exx(yn+1,xn+1)
        MEXX(mm1)=exxm
        ! Calculate Marker SIG'xx from four surrounding nodes
        sxxm=0.0d0
        sxxm=sxxm+(1.0d0-dx)*(1.0d0-dy)*sxx2(yn,xn)
        sxxm=sxxm+(1.0d0-dx)*dy*sxx2(yn+1,xn)
        sxxm=sxxm+dx*(1.0d0-dy)*sxx2(yn,xn+1)
        sxxm=sxxm+dx*dy*sxx2(yn+1,xn+1)
        ! Calculate Marker dSIG'xx from four surrounding nodes
        dsxxm=0.0d0
        dsxxm=dsxxm+(1.0d0-dx)*(1.0d0-dy)*dsxx(yn,xn)
        dsxxm=dsxxm+(1.0d0-dx)*dy*dsxx(yn+1,xn)
        dsxxm=dsxxm+dx*(1.0d0-dy)*dsxx(yn,xn+1)
        dsxxm=dsxxm+dx*dy*dsxx(yn+1,xn+1)

        ! Calculate and save Marker pressure from four surrounding nodes
        prm=0.d0
        prm=prm+(1.0d0-dx)*(1.0d0-dy)*pr1(yn,xn)
        prm=prm+(1.0d0-dx)*dy*pr1(yn+1,xn)
        prm=prm+dx*(1.0d0-dy)*pr1(yn,xn+1)
        prm=prm+dx*dy*pr1(yn+1,xn+1)
        MPR(mm1)=prm

        ! Interpolating EPSxy=EPSyx from basic nodes
        ! Horizontal EPSxy index
        xn=MXN(mm1)
        ! Vertical EPSxy index
        yn=MYN(mm1)

        ! Define and check normalized distances from marker to the upper left VX-node;
        ! dx=(MX(mm1)-gridx(xn))/xstp1(xn);
        ! dy=(MY(mm1)-gridy(yn))/ystp1(yn);
        dx = MDX(mm1)
        dy = MDY(mm1)

        ! Calculate and save Marker EPSxy from four surrounding nodes
        exym=0.0d0
        exym=exym+(1.0d0-dx)*(1.0d0-dy)*exy(yn,xn)
        exym=exym+(1.0d0-dx)*dy*exy(yn+1,xn)
        exym=exym+dx*(1.0d0-dy)*exy(yn,xn+1)
        exym=exym+dx*dy*exy(yn+1,xn+1)
        MEXY(mm1)=exym
        ! Calculate Marker SIGxy from four surrounding nodes
        sxym=0.0d0
        sxym=sxym+(1.0d0-dx)*(1.0d0-dy)*sxy2(yn,xn)
        sxym=sxym+(1.0d0-dx)*dy*sxy2(yn+1,xn)
        sxym=sxym+dx*(1.0d0-dy)*sxy2(yn,xn+1)
        sxym=sxym+dx*dy*sxy2(yn+1,xn+1)
        ! Calculate Marker SIGxy from four surrounding nodes
        dsxym=0.0d0
        dsxym=dsxym+(1.0d0-dx)*(1.0d0-dy)*dsxy(yn,xn)
        dsxym=dsxym+(1.0d0-dx)*dy*dsxy(yn+1,xn)
        dsxym=dsxym+dx*(1.0d0-dy)*dsxy(yn,xn+1)
        dsxym=dsxym+dx*dy*dsxy(yn+1,xn+1)

        ! Computing second strain rate invariant
        ! for the marker using grid values
        eiimg=dsqrt(MEXX(mm1)**2+MEXY(mm1)**2)
        ! Correcting strain rate for the marker using Maxwell model
        if (eiimg>0.d0) then
            ! Computing second strain rate invariant for the marker
            ! from stresses using Maxwell model
            eiim=dsqrt((sxxm/2/META(mm1)+dsxxm/2/timestep/MMU(mm1))**2+(sxym/2/META(mm1)+dsxym/2/timestep/MMU(mm1))**2)
            ! Computing EiiMarker/EiiGrid ratio
            MRAT(mm1)=(eiim/eiimg)
        end if
    end do

    ! Prepare for output
    do i = 1,4
        plhs(i) = mxCreateDoubleMatrix(marknum,1,0)
    end do

    call mxCopyReal8ToPtr(MPR,mxGetPr(plhs(1)),marknum)
    call mxCopyReal8ToPtr(MRAT,mxGetPr(plhs(2)),marknum)
    call mxCopyReal8ToPtr(MEXX,mxGetPr(plhs(3)),marknum)
    call mxCopyReal8ToPtr(MEXY,mxGetPr(plhs(4)),marknum)

    return
end subroutine
