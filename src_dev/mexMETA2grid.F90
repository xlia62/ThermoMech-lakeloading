# include <fintrf.h>
!======================================================================
#if 0
!     Input
!     (MXN,MYN,MDX,MDY,MRHOCUR,MTK,MKTCUR,MRHOCPCUR,MHR,MHE,MHACUR,METACUR,G)
!
!     Output
!     structure G
#endif
!
!     timestwo.f
!
!    Computational function that takes a scalar and doubles it.

!     This is a MEX-file for MATLAB.
!    Copyright 1984-2011 The MathWorks, Inc.
!
!======================================================================


!     Gateway routine
subroutine mexFunction(nlhs, plhs, nrhs, prhs)


  implicit none

  !     mexFunction arguments:
  mwPointer :: plhs(*), prhs(*)
  integer :: nlhs, nrhs

  !     Function declarations:
  mwPointer :: mxGetPr, mxGetField
  mwPointer :: mxCreateDoubleMatrix, mxCreateStructArray
  !integer :: mexPrintf
  integer :: mxGetNumberOfFields, mxAddField, mxGetFieldNumber
  mwPointer :: mxGetM, mxGetN, mxGetData, mxDuplicateArray
  character(10) :: mxGetFieldNameByNumber

  !     Pointers to input/output mxArrays:
  mwPointer :: inmarker_ptr, ingrid_ptr

  !  Arguments for computational routine:
  double precision, allocatable, dimension(:) :: MX, MY, MDX, MDY
  double precision, allocatable, dimension(:) :: METACUR
  integer, allocatable, dimension (:) :: MXN, MYN
  double precision, allocatable, dimension(:) :: tmpD
  double precision, allocatable, dimension(:,:) :: etas1, etan1, etas0, etan0
  double precision, allocatable, dimension(:,:) :: wtetas, wtetan, wtnodes

  ! Local variables
  integer :: nfields, marknum, xnum, ynum, numel, numelc
  integer :: i, k, mm1
  integer :: xn, yn
  integer*4 p
  
  double precision :: indble, dx, dy, mwt


  character*120 line
    
  !  call mexprintf("Here 1 ")

  !-----------------------------------------------------------------------
  ! Number of markers
  marknum = mxGetM(prhs(1))*mxGetN(prhs(1))

  allocate(MXN(marknum),MYN(marknum))
  allocate(MDX(marknum),MDY(marknum))
  allocate(METACUR(marknum))
  
  ! Populate marker arrays: ! MXN,MYN,MDX,MDY,METACUR,etas0,etan0
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
  call mxCopyPtrToReal8(inmarker_ptr,METACUR,marknum)
   
  ! size of etas0 
  ynum = mxGetM(prhs(6))
  xnum = mxGetN(prhs(6))

  ! Allocate output grid viscosity 
  allocate(etas0(ynum,xnum),etas1(ynum,xnum))
  allocate(etan0(ynum-1,xnum-1),etan1(ynum-1,xnum-1))
 
  numel = ynum*xnum
  numelc = (ynum-1)*(xnum-1)

  ! Read stress
  ingrid_ptr = mxGetPr(prhs(6))
  call mxCopyPtrToReal8(ingrid_ptr,etas0,numel)
  ingrid_ptr = mxGetPr(prhs(7))
  call mxCopyPtrToReal8(ingrid_ptr,etan0,numelc)

  
  allocate(wtetas(ynum,xnum),wtetan(ynum-1,xnum-1))
  
  wtetas = 0.0d0
  wtetan = 0.0d0
  etas1 = 0.0d0
  etan1 = 0.0d0
    
  ! Interpolating parameters from markers to nodes
  do mm1=1,marknum

    xn = MXN(mm1)
    yn = MYN(mm1)

    ! Define normalized distances from marker to the upper left node;
    dx=MDX(mm1)
    dy=MDY(mm1)

    ! Compute marker weight koefficient from cell dimensions
    ! Number of markers in a cell is in invert proportion to the cell volume
    mwt=1 !/xstp1(xn)/ystp1(yn);

    ! Add viscosity etas(), shear stress sxy(),shear modulus mus() and rock type typ() to 4 surrounding basic nodes
    ! only using markers located at <=0.5 gridstep distances from nodes
    if (dx<=0.5 .and. dy<=0.5) then
      etas1(yn,xn) = etas1(yn,xn) + (1.0-dx)*(1.0-dy)*METACUR(mm1)*mwt
      wtetas(yn,xn) = wtetas(yn,xn) + (1.0-dx)*(1.0-dy)*mwt
    end if
    if (dx<=0.5 .and. dy>=0.5) then
      etas1(yn+1,xn) = etas1(yn+1,xn) + (1.0-dx)*dy*METACUR(mm1)*mwt
      wtetas(yn+1,xn) = wtetas(yn+1,xn) + (1.0-dx)*dy*mwt
    end if
    if (dx>=0.5 .and. dy<=0.5) then
      etas1(yn,xn+1) = etas1(yn,xn+1) + dx*(1.0-dy)*METACUR(mm1)*mwt
      wtetas(yn,xn+1) = wtetas(yn,xn+1) + dx*(1.0-dy)*mwt
    end if
    if (dx>=0.5 .and. dy>=0.5) then
      etas1(yn+1,xn+1) = etas1(yn+1,xn+1) + dx*dy*METACUR(mm1)*mwt
      wtetas(yn+1,xn+1) = wtetas(yn+1,xn+1) + dx*dy*mwt
    end if

    ! Add viscosity etan(), normal stress sxx() and sGin(6)%propar modulus mun() to tGin(6)%prop center of current cell
    etan1(yn,xn) = etan1(yn,xn) + (1.0-dabs(0.5-dx))*(1.0-dabs(0.5-dy))*METACUR(mm1)*mwt
    wtetan(yn,xn) = wtetan(yn,xn) + (1.0-dabs(0.5-dx))*(1.0-dabs(0.5-dy))*mwt
  end do

  ! Set etas viscosity
  where (wtetas > 0.0d0)
    ! Compute new value interpolated from markers
    etas1=etas1/wtetas
  elsewhere
    ! If no new value is interpolated from markers old value is used
    etas1=etas0
  end where
      
  ! Normal viscosity
  where (wtetan > 0.0d0)
    ! Compute new value interpolated from markers
    etan1=etan1/wtetan
  elsewhere
    ! If no new value is interpolated from markers old value is used
    etan1=etan0
  end where
   
  ! Prepare for output
  plhs(1) = mxCreateDoubleMatrix(ynum,xnum,0.0d0)
  plhs(2) = mxCreateDoubleMatrix(ynum-1,xnum-1,0.0d0)

  ! Send to out
  call mxCopyReal8ToPtr(etas1,mxGetPr(plhs(1)),numel) 
  call mxCopyReal8ToPtr(etan1,mxGetPr(plhs(2)),numelc) 
      
end subroutine
