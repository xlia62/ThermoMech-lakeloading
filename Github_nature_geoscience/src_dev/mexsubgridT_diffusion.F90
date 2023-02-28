# include <fintrf.h>
!======================================================================
#if 0
!     Input
!     (MXN,MYN,MDX,MDY,G.tk1,G.rhocp1,dtk1,G.kt1,xstp1,ystp1,dsubgridt,timestep)
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
    double precision, allocatable, dimension(:) :: MDX, MDY, MTK
    integer, allocatable, dimension (:) :: MXN, MYN
    double precision, allocatable, dimension(:) :: xstp1, ystp1, tmpD

    double precision, allocatable, dimension(:,:) :: tk1, rhocp1, dtk1, kt1, dtkn, wtnodes
    ! Local variables
    integer :: marknum, xnum, ynum, numel
    integer :: i, k, mm1
    integer :: xn, yn

    double precision :: dx, dy, timestep, dsubgridt, sdif
    double precision :: tkm, ktm, rhocpm, tdm, dtkm, mwt

    character*120 line

    !-----------------------------------------------------------------------
    ! Number of markers
    marknum = mxGetM(prhs(1))*mxGetN(prhs(1))

    allocate(MXN(marknum),MYN(marknum),MDX(marknum),MDY(marknum),MTK(marknum))

    ! Populate marker arrays: MXN,MYN,MDX,MDY,MTK

    allocate(tmpD(marknum))

    inmarker_ptr = mxGetPr(prhs(1))
    call mxCopyPtrToReal8(inmarker_ptr,tmpD,marknum)
    MXN = int(tmpD)
    inmarker_ptr = mxGetPr(prhs(2))
    call mxCopyPtrToReal8(inmarker_ptr,tmpD,marknum)
    MYN = int(tmpD)

    inmarker_ptr = mxGetPr(prhs(3))
    call mxCopyPtrToReal8(inmarker_ptr,MDX,marknum)
    inmarker_ptr = mxGetPr(prhs(4))
    call mxCopyPtrToReal8(inmarker_ptr,MDY,marknum)
    inmarker_ptr = mxGetPr(prhs(5))
    call mxCopyPtrToReal8(inmarker_ptr,MTK,marknum)

    ! Populate grids: G.tk1,G.rhocp1,dtk1,G.kt1,xstp1,ystp1
    ynum = mxGetM(prhs(6))
    xnum = mxGetN(prhs(6))

    allocate(tk1(ynum,xnum),rhocp1(ynum,xnum),dtk1(ynum,xnum),kt1(ynum,xnum))
    allocate(xstp1(xnum-1),ystp1(ynum-1))

    numel = ynum*xnum

    ingrid_ptr = mxGetPr(prhs(6))
    call mxCopyPtrToReal8(ingrid_ptr,tk1,numel)
    ingrid_ptr = mxGetPr(prhs(7))
    call mxCopyPtrToReal8(ingrid_ptr,rhocp1,numel)
    ingrid_ptr = mxGetPr(prhs(8))
    call mxCopyPtrToReal8(ingrid_ptr,dtk1,numel)
    ingrid_ptr = mxGetPr(prhs(9))
    call mxCopyPtrToReal8(ingrid_ptr,kt1,numel)

    ingrid_ptr = mxGetPr(prhs(10))
    call mxCopyPtrToReal8(ingrid_ptr,xstp1,xnum-1)
    ingrid_ptr = mxGetPr(prhs(11))
    call mxCopyPtrToReal8(ingrid_ptr,ystp1,ynum-1)


    ingrid_ptr = mxGetPr(prhs(12))
    call mxCopyPtrToReal8(ingrid_ptr,dsubgridt,1)
    ingrid_ptr = mxGetPr(prhs(13))
    call mxCopyPtrToReal8(ingrid_ptr,timestep,1)


    allocate(dtkn(ynum,xnum))
    dtkn = 0.0d0

    allocate(wtnodes(ynum,xnum))
    wtnodes = 0.0d0

    ! Marker cycle
    do mm1=1,marknum

        ! Interpolating temperature changes from basic nodes
        !
        ! Define indexes for upper left node in the cell where the marker is
        xn=MXN(mm1)
        yn=MYN(mm1)

        ! Define normalized distances from marker to the upper left node;
        dx=MDX(mm1)
        dy=MDY(mm1)

        ! Compute marker weight koefficient from cell dimensions
        ! Number of markers in a cell is in invert proportion to the cell volume
        mwt=1.0d0 !/xstp1(xn)/ystp1(yn);


        ! Interpolate old nodal temperature for the marker
        tkm=0.0d0
        tkm=tkm+(1.0-dx)*(1.0-dy)*tk1(yn,xn)
        tkm=tkm+(1.0-dx)*dy*tk1(yn+1,xn)
        tkm=tkm+dx*(1.0-dy)*tk1(yn,xn+1)
        tkm=tkm+dx*dy*tk1(yn+1,xn+1)
        ! Calculate Nodal-Marker subgrid temperature difference
        dtkm=tkm-MTK(mm1)
        ! Compute nodal k and RHO*Cp for the marker
        ! k
        ktm=0.0d0
        ktm=ktm+(1.0-dx)*(1.0-dy)*kt1(yn,xn)
        ktm=ktm+(1.0-dx)*dy*kt1(yn+1,xn)
        ktm=ktm+dx*(1.0-dy)*kt1(yn,xn+1)
        ktm=ktm+dx*dy*kt1(yn+1,xn+1)
        ! RHO*Cp
        rhocpm=0.0d0
        rhocpm=rhocpm+(1.0-dx)*(1.0-dy)*rhocp1(yn,xn)
        rhocpm=rhocpm+(1.0-dx)*dy*rhocp1(yn+1,xn)
        rhocpm=rhocpm+dx*(1.0-dy)*rhocp1(yn,xn+1)
        rhocpm=rhocpm+dx*dy*rhocp1(yn+1,xn+1)

        ! Compute local thermal diffusion timescale for the marker
        tdm=rhocpm/ktm/(2/xstp1(xn)**2+2/ystp1(yn)**2)

        ! Computing subgrid diffusion
        sdif=-dsubgridt*timestep/tdm
        if (sdif<-30.0d0) then
            sdif=-30.0d0
        end if
        dtkm=dtkm*(1-exp(sdif))

        ! Correcting old temperature for the marker
        MTK(mm1)=MTK(mm1)+dtkm

        ! Interpolating subgrid temperature changes to 4 nodes
        dtkn(yn,xn)=dtkn(yn,xn)+(1.0-dx)*(1.0-dy)*dtkm*mwt
        wtnodes(yn,xn)=wtnodes(yn,xn)+(1.0-dx)*(1.0-dy)*mwt

        dtkn(yn+1,xn)=dtkn(yn+1,xn)+(1.0-dx)*dy*dtkm*mwt
        wtnodes(yn+1,xn)=wtnodes(yn+1,xn)+(1.0-dx)*dy*mwt

        dtkn(yn,xn+1)=dtkn(yn,xn+1)+dx*(1.0-dy)*dtkm*mwt
        wtnodes(yn,xn+1)=wtnodes(yn,xn+1)+dx*(1.0-dy)*mwt

        dtkn(yn+1,xn+1)=dtkn(yn+1,xn+1)+dx*dy*dtkm*mwt
        wtnodes(yn+1,xn+1)=wtnodes(yn+1,xn+1)+dx*dy*mwt

    end do

    ! Computing subgrid diffusion for nodes
    where(wtnodes/=0) dtkn = dtkn/wtnodes

    ! Subtracting subgrid diffusion part from nodal temperature changes
    dtk1=dtk1-dtkn

    ! Updating temperature for markers
    do mm1=1,marknum,1

        ! Interpolating temperature changes from basic nodes
        !
        ! Define indexes for upper left node in the cell where the marker is
        xn=MXN(mm1)
        yn=MYN(mm1)

        ! Define normalized distances from marker to the upper left node;
        dx=MDX(mm1)
        dy=MDY(mm1)

        ! Calculate Marker temperature change from four surrounding nodes
        dtkm=0.0d0
        dtkm=dtkm+(1.0-dx)*(1.0-dy)*dtk1(yn,xn)
        dtkm=dtkm+(1.0-dx)*dy*dtk1(yn+1,xn)
        dtkm=dtkm+dx*(1.0-dy)*dtk1(yn,xn+1)
        dtkm=dtkm+dx*dy*dtk1(yn+1,xn+1)
        !
        !Computing new temperature for the marker
        MTK(mm1)=MTK(mm1)+dtkm

    end do

    ! Prepare for output
    plhs(1) = mxCreateDoubleMatrix(marknum,1,0)
    call mxCopyReal8ToPtr(MTK,mxGetPr(plhs(1)),marknum)

    return

end subroutine mexFunction
