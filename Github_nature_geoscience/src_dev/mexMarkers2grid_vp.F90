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

module gridtype

  implicit none

  type grid1

    mwPointer :: ptr
    double precision, allocatable, dimension(:,:) :: prop

  end type grid1

  type grid2

    mwPointer :: ptr
    double precision, allocatable, dimension(:) :: prop

  end type grid2


end module gridtype

!     Gateway routine
subroutine mexFunction(nlhs, plhs, nrhs, prhs)

  use gridtype

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
  double precision, allocatable, dimension(:) :: MHR, MHE, MTK
  double precision, allocatable, dimension(:) :: MRHOCUR, MRHOCPCUR, METACUR
  double precision, allocatable, dimension(:) :: MHACUR, MKTCUR
  integer, allocatable, dimension (:) :: MXN, MYN
  double precision, allocatable, dimension(:) :: tmpD
  double precision, allocatable, dimension(:,:) :: wtetas, wtetan, wtnodes
  type(grid1), dimension(9) ::  Gin, Gout

  ! Local variables
  integer :: nfields, marknum, xnum, ynum, numel, numelc
  integer :: i, k, mm1
  integer :: xn, yn
  integer*4 p
  
  double precision :: dx, dy, mwt


  character*120 line
    
!  call mexprintf("Here 1 ")

  !-----------------------------------------------------------------------
  ! Number of markers
  marknum = mxGetM(prhs(1))*mxGetN(prhs(1))

  allocate(MXN(marknum),MYN(marknum))
  allocate(MDX(marknum),MDY(marknum))
  allocate(MHR(marknum),MHE(marknum),MTK(marknum))
  allocate(METACUR(marknum),MKTCUR(marknum),MRHOCUR(marknum))
  allocate(MHACUR(marknum),MRHOCPCUR(marknum))

  ! Populate marker arrays: ! MXN,MYN,MDX,MDY,MRHOCUR,MTK,MKTCUR,MRHOCPCUR,MHR,MHE,MHACUR,METACUR,MMUCUR,G)

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
  call mxCopyPtrToReal8(inmarker_ptr,MRHOCUR,marknum)
  inmarker_ptr = mxGetPr(prhs(6))
  call mxCopyPtrToReal8(inmarker_ptr,MTK,marknum)
  inmarker_ptr = mxGetPr(prhs(7))
  call mxCopyPtrToReal8(inmarker_ptr,MKTCUR,marknum)
  inmarker_ptr = mxGetPr(prhs(8))
  call mxCopyPtrToReal8(inmarker_ptr,MRHOCPCUR,marknum)

  inmarker_ptr = mxGetPr(prhs(9))
  call mxCopyPtrToReal8(inmarker_ptr,MHR,marknum)
  inmarker_ptr = mxGetPr(prhs(10))
  call mxCopyPtrToReal8(inmarker_ptr,MHE,marknum)
  inmarker_ptr = mxGetPr(prhs(11))
  call mxCopyPtrToReal8(inmarker_ptr,MHACUR,marknum)
  inmarker_ptr = mxGetPr(prhs(12))
  call mxCopyPtrToReal8(inmarker_ptr,METACUR,marknum)
   
!   call mexprintf("Here 2 ")
 
  ! Extract the fields from the grid structure
  nfields = mxGetNumberOfFields(prhs(13))

  Gin(1)%ptr = mxGetField(prhs(13), 1,'rho1')
  Gin(2)%ptr = mxGetField(prhs(13), 1,'tk1')
  Gin(3)%ptr = mxGetField(prhs(13), 1,'kt1')
  Gin(4)%ptr = mxGetField(prhs(13), 1,'rhocp1')
  Gin(5)%ptr = mxGetField(prhs(13), 1,'hr1')
  Gin(6)%ptr = mxGetField(prhs(13), 1,'he')
  Gin(7)%ptr = mxGetField(prhs(13), 1,'ha1')
  Gin(8)%ptr = mxGetField(prhs(13), 1,'etas1')
  Gin(9)%ptr = mxGetField(prhs(13), 1,'etan1')

!  call mexprintf("Here 3 ")
    
  ! Create structue for the return argument
  plhs(1) = mxCreateStructArray(1,1,1,mxGetFieldNameByNumber(prhs(13), 1))

  do i = 2,nfields
    k = mxAddField(plhs(1),mxGetFieldNameByNumber(prhs(13), i))
  end do

  ! Grid Arrays
  ynum = mxGetM(Gin(1)%ptr)
  xnum = mxGetN(Gin(1)%ptr)

  !     Allocate Grid Arrays
  do k = 1,8
    allocate(Gin(k)%prop(ynum,xnum))
    allocate(Gout(k)%prop(ynum,xnum))
  end do
  allocate(Gin(9)%prop(ynum-1,xnum-1))
  allocate(Gout(9)%prop(ynum-1,xnum-1))
  
!  call mexprintf("Here 4 ")

  ! Populate the Input Arrays
  numel = ynum*xnum
  numelc = (ynum-1)*(xnum-1)
  
  do k = 1,8
    call mxCopyPtrToReal8(mxGetPr(Gin(k)%ptr),Gin(k)%prop,numel)
  end do
  call mxCopyPtrToReal8(mxGetPr(Gin(9)%ptr),Gin(9)%prop,numelc)

  allocate(wtetas(ynum,xnum),wtetan(ynum-1,xnum-1),wtnodes(ynum,xnum))
  wtetas = 0.0d0
  wtetan = 0.0d0
  wtnodes = 0.0d0
  
   do i = 1,9
         Gout(i)%prop = 0.0d0
   end do
  
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
      
    ! Add properties to 4 surrounding nodes
    Gout(1)%prop(yn,xn)=Gout(1)%prop(yn,xn)+(1.0-dx)*(1.0-dy)*MRHOCUR(mm1)*mwt
    Gout(2)%prop(yn,xn)=Gout(2)%prop(yn,xn)+(1.0-dx)*(1.0-dy)*MTK(mm1)*mwt
    Gout(3)%prop(yn,xn)=Gout(3)%prop(yn,xn)+(1.0-dx)*(1.0-dy)*MKTCUR(mm1)*mwt
    Gout(4)%prop(yn,xn)=Gout(4)%prop(yn,xn)+(1.0-dx)*(1.0-dy)*MRHOCPCUR(mm1)*mwt
    Gout(5)%prop(yn,xn)=Gout(5)%prop(yn,xn)+(1.0-dx)*(1.0-dy)*MHR(mm1)*mwt
    Gout(7)%prop(yn,xn)=Gout(7)%prop(yn,xn)+(1.0-dx)*(1.0-dy)*MHACUR(mm1)*mwt
    Gout(6)%prop(yn,xn)=Gout(6)%prop(yn,xn)+(1.0-dx)*(1.0-dy)*MHE(mm1)*mwt
    wtnodes(yn,xn)=wtnodes(yn,xn)+(1.0-dx)*(1.0-dy)*mwt

    Gout(1)%prop(yn+1,xn)=Gout(1)%prop(yn+1,xn)+(1.0-dx)*dy*MRHOCUR(mm1)*mwt
    Gout(2)%prop(yn+1,xn)=Gout(2)%prop(yn+1,xn)+(1.0-dx)*dy*MTK(mm1)*mwt
    Gout(3)%prop(yn+1,xn)=Gout(3)%prop(yn+1,xn)+(1.0-dx)*dy*MKTCUR(mm1)*mwt
    Gout(4)%prop(yn+1,xn)=Gout(4)%prop(yn+1,xn)+(1.0-dx)*dy*MRHOCPCUR(mm1)*mwt
    Gout(5)%prop(yn+1,xn)=Gout(5)%prop(yn+1,xn)+(1.0-dx)*dy*MHR(mm1)*mwt
    Gout(7)%prop(yn+1,xn)=Gout(7)%prop(yn+1,xn)+(1.0-dx)*dy*MHACUR(mm1)*mwt
    Gout(6)%prop(yn+1,xn)=Gout(6)%prop(yn+1,xn)+(1.0-dx)*dy*MHE(mm1)*mwt
    wtnodes(yn+1,xn)=wtnodes(yn+1,xn)+(1.0-dx)*dy*mwt

    Gout(1)%prop(yn,xn+1)=Gout(1)%prop(yn,xn+1)+dx*(1.0-dy)*MRHOCUR(mm1)*mwt
    Gout(2)%prop(yn,xn+1)=Gout(2)%prop(yn,xn+1)+dx*(1.0-dy)*MTK(mm1)*mwt
    Gout(3)%prop(yn,xn+1)=Gout(3)%prop(yn,xn+1)+dx*(1.0-dy)*MKTCUR(mm1)*mwt
    Gout(4)%prop(yn,xn+1)=Gout(4)%prop(yn,xn+1)+dx*(1.0-dy)*MRHOCPCUR(mm1)*mwt
    Gout(5)%prop(yn,xn+1)=Gout(5)%prop(yn,xn+1)+dx*(1.0-dy)*MHR(mm1)*mwt
    Gout(7)%prop(yn,xn+1)=Gout(7)%prop(yn,xn+1)+dx*(1.0-dy)*MHACUR(mm1)*mwt
    Gout(6)%prop(yn,xn+1)=Gout(6)%prop(yn,xn+1)+dx*(1.0-dy)*MHE(mm1)*mwt
    wtnodes(yn,xn+1)=wtnodes(yn,xn+1)+dx*(1.0-dy)*mwt

    Gout(1)%prop(yn+1,xn+1)=Gout(1)%prop(yn+1,xn+1)+dx*dy*MRHOCUR(mm1)*mwt
    Gout(2)%prop(yn+1,xn+1)=Gout(2)%prop(yn+1,xn+1)+dx*dy*MTK(mm1)*mwt
    Gout(3)%prop(yn+1,xn+1)=Gout(3)%prop(yn+1,xn+1)+dx*dy*MKTCUR(mm1)*mwt
    Gout(4)%prop(yn+1,xn+1)=Gout(4)%prop(yn+1,xn+1)+dx*dy*MRHOCPCUR(mm1)*mwt
    Gout(5)%prop(yn+1,xn+1)=Gout(5)%prop(yn+1,xn+1)+dx*dy*MHR(mm1)*mwt
    Gout(7)%prop(yn+1,xn+1)=Gout(7)%prop(yn+1,xn+1)+dx*dy*MHACUR(mm1)*mwt
    Gout(6)%prop(yn+1,xn+1)=Gout(6)%prop(yn+1,xn+1)+dx*dy*MHE(mm1)*mwt
    wtnodes(yn+1,xn+1)=wtnodes(yn+1,xn+1)+dx*dy*mwt


    ! Add viscosity etas(), shear stress sxy(),shear modulus mus() and rock type typ() to 4 surrounding basic nodes
    ! only using markers located at <=0.5 gridstep distances from nodes
    if (dx<=0.5 .and. dy<=0.5) then
      Gout(8)%prop(yn,xn)=Gout(8)%prop(yn,xn)+(1.0-dx)*(1.0-dy)*METACUR(mm1)*mwt
      wtetas(yn,xn)=wtetas(yn,xn)+(1.0-dx)*(1.0-dy)*mwt
    end if
    if (dx<=0.5 .and. dy>=0.5) then
      Gout(8)%prop(yn+1,xn)=Gout(8)%prop(yn+1,xn)+(1.0-dx)*dy*METACUR(mm1)*mwt
      wtetas(yn+1,xn)=wtetas(yn+1,xn)+(1.0-dx)*dy*mwt
    end if
    if (dx>=0.5 .and. dy<=0.5) then
      Gout(8)%prop(yn,xn+1)=Gout(8)%prop(yn,xn+1)+dx*(1.0-dy)*METACUR(mm1)*mwt
      wtetas(yn,xn+1)=wtetas(yn,xn+1)+dx*(1.0-dy)*mwt
    end if
    if (dx>=0.5 .and. dy>=0.5) then
      Gout(8)%prop(yn+1,xn+1)=Gout(8)%prop(yn+1,xn+1)+dx*dy*METACUR(mm1)*mwt
      wtetas(yn+1,xn+1)=wtetas(yn+1,xn+1)+dx*dy*mwt
    end if

    ! Add viscosity etan(), normal stress sxx() and sGin(6)%propar modulus mun() to tGin(6)%prop center of current cell
    Gout(9)%prop(yn,xn)=Gout(9)%prop(yn,xn)+(1.0-dabs(0.5-dx))*(1.0-dabs(0.5-dy))*METACUR(mm1)*mwt
    wtetan(yn,xn)=wtetan(yn,xn)+(1.0-dabs(0.5-dx))*(1.0-dabs(0.5-dy))*mwt
  end do

  ! Computing Viscosity, density, rock type for nodal points
  where (wtnodes > 0.0d0) 
    Gout(1)%prop=Gout(1)%prop/wtnodes
    Gout(2)%prop=Gout(2)%prop/wtnodes
    Gout(3)%prop=Gout(3)%prop/wtnodes
    Gout(4)%prop=Gout(4)%prop/wtnodes
    Gout(5)%prop=Gout(5)%prop/wtnodes
    Gout(6)%prop=Gout(6)%prop/wtnodes
    Gout(7)%prop=Gout(7)%prop/wtnodes
  elsewhere
    ! If no new value is interpolated from markers old value is used
    Gout(1)%prop=Gin(1)%prop
    Gout(2)%prop=Gin(2)%prop
    Gout(3)%prop=Gin(3)%prop
    Gout(4)%prop=Gin(4)%prop
    Gout(5)%prop=Gin(5)%prop
    Gout(6)%prop=Gin(6)%prop
    Gout(7)%prop=Gin(7)%prop
  end where
        
  ! Set etas viscosity
  where (wtetas > 0.0d0)
    ! Compute new value interpolated from markers
    Gout(8)%prop=Gout(8)%prop/wtetas
  elsewhere
    ! If no new value is interpolated from markers old value is used
    Gout(8)%prop=Gin(8)%prop
  end where
      
  ! Normal viscosity
  where (wtetan > 0.0d0)
    ! Compute new value interpolated from markers
    Gout(9)%prop=Gout(9)%prop/wtetan
  elsewhere
    ! If no new value is interpolated from markers old value is used
    Gout(9)%prop=Gin(9)%prop
  end where
      
  ! Populate the Output Structure
  call mxSetField(plhs(1),1,'xsize',mxDuplicateArray(mxGetField(prhs(13),1,'xsize')))
  call mxSetField(plhs(1),1,'ysize',mxDuplicateArray(mxGetField(prhs(13),1,'ysize')))
  call mxSetField(plhs(1),1,'xnum',mxDuplicateArray(mxGetField(prhs(13),1,'xnum')))
  call mxSetField(plhs(1),1,'ynum',mxDuplicateArray(mxGetField(prhs(13),1,'ynum')))
  call mxSetField(plhs(1),1,'bx',mxDuplicateArray(mxGetField(prhs(13),1,'bx')))
  call mxSetField(plhs(1),1,'by',mxDuplicateArray(mxGetField(prhs(13),1,'by')))
  call mxSetField(plhs(1),1,'wx',mxDuplicateArray(mxGetField(prhs(13),1,'wx')))
  call mxSetField(plhs(1),1,'wy',mxDuplicateArray(mxGetField(prhs(13),1,'wy')))
  call mxSetField(plhs(1),1,'f',mxDuplicateArray(mxGetField(prhs(13),1,'f')))
  call mxSetField(plhs(1),1,'gridx',mxDuplicateArray(mxGetField(prhs(13),1,'gridx')))
  call mxSetField(plhs(1),1,'gridy',mxDuplicateArray(mxGetField(prhs(13),1,'gridy')))
  call mxSetField(plhs(1),1,'tk2',mxDuplicateArray(mxGetField(prhs(13),1,'tk2')))
  call mxSetField(plhs(1),1,'sxx1',mxDuplicateArray(mxGetField(prhs(13),1,'sxx1')))
  call mxSetField(plhs(1),1,'sxy1',mxDuplicateArray(mxGetField(prhs(13),1,'sxy1')))
        
  do i = 1,8
    Gout(i)%ptr = mxCreateDoubleMatrix(ynum,xnum,0.d0)
  end do
               
  call mxCopyReal8ToPtr(Gout(1)%prop,mxGetPr(Gout(1)%ptr),numel)
  call mxSetField(plhs(1),1,'rho1',Gout(1)%ptr)
  call mxCopyReal8ToPtr(Gout(2)%prop,mxGetPr(Gout(2)%ptr),numel)
  call mxSetField(plhs(1),1,'tk1',Gout(2)%ptr)
  call mxCopyReal8ToPtr(Gout(3)%prop,mxGetPr(Gout(3)%ptr),numel)
  call mxSetField(plhs(1),1,'kt1',Gout(3)%ptr)
  call mxCopyReal8ToPtr(Gout(4)%prop,mxGetPr(Gout(4)%ptr),numel)
  call mxSetField(plhs(1),1,'rhocp1',Gout(4)%ptr)
  call mxCopyReal8ToPtr(Gout(5)%prop,mxGetPr(Gout(5)%ptr),numel)
  call mxSetField(plhs(1),1,'hr1',Gout(5)%ptr)
  call mxCopyReal8ToPtr(Gout(6)%prop,mxGetPr(Gout(6)%ptr),numel)
  call mxSetField(plhs(1),1,'he',Gout(6)%ptr)
  call mxCopyReal8ToPtr(Gout(7)%prop,mxGetPr(Gout(7)%ptr),numel)
  call mxSetField(plhs(1),1,'ha1',Gout(7)%ptr)
  call mxCopyReal8ToPtr(Gout(8)%prop,mxGetPr(Gout(8)%ptr),numel)
  call mxSetField(plhs(1),1,'etas1',Gout(8)%ptr)

  Gout(9)%ptr = mxCreateDoubleMatrix(ynum-1,xnum-1,0.d0)

  call mxCopyReal8ToPtr(Gout(9)%prop,mxGetPr(Gout(9)%ptr),numelc)
  call mxSetField(plhs(1),1,'etan1',Gout(9)%ptr)
      
end subroutine
