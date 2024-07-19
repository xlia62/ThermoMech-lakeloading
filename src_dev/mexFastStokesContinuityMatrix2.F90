# include <fintrf.h>
!======================================================================
#if 0
!     Input
!     (MXN,MYN,MDX,MDY,MRHOCUR,MTK,MKTCUR,MRHOCPCUR,MHR,MHE,MHACUR,METACUR,MMUCUR,MSXX,MSXY,G)
!
!     Output
!     structure G
#endif
! Function Stokes_Continuity_solver_sandbox() 
! This function formulates and solves 
! Stokes and Continuity equations defined on 2D staggered irregularly spaced grid 
! with specified resolution (xnum, ynum) and grid lines positions (gridx, gridy) 
! given distribution of right parts for all equations (RX,RY,RC) on the grid 
! and given variable shear (etas) and normal (etan) viscosity distributions 
! pressure is normalized relative to given value (prnorm) in the first cell 
! 
! Velocity Boundary condition specified by bleft,bright,btop,bbottom,bintern 
! are implemented from ghost nodes 
! directly into Stokes and continuity equations 
! 
! Function returns solution for velocity and pressure (vx,vy,pr) 
! and distribution of residuals (resx,resy,resc) 
!

!  Gateway routine
!  function [Im,Jm,Lm,Rm]=Fast_Stokes_Continuity_solver_sandbox(prfirst,etas,etan,gridx,gridy,RX,RY,RC,bleft,bright,btop,bbottom,bintern,psale)

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

  !     Declarations
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
  mwPointer :: press_ptr, etas_ptr, etan_ptr, gridx_ptr, gridy_ptr, RX_ptr, RY_ptr, RC_ptr, bc_ptr

  !     Arguments for computational routine:
  double precision, allocatable, dimension(:) :: gridx,gridy, xstp, ystp, xstpc, ystpc
  double precision, allocatable, dimension(:,:) :: etas, etan, RX, RY, RC
  double precision, allocatable, dimension(:,:) :: btop, bbottom, bleft, bright
  double precision, allocatable, dimension(:,:) :: btopx, btopy, bbottomx, bbottomy, bleftx, blefty, brightx, brighty
  double precision, dimension(8,4) :: bintern

  !    Output sparse array structure
  integer, allocatable, dimension(:) :: Im, Jm, tmpI, tmpJ
  double precision, allocatable, dimension (:) :: Lm, Rm, tmpL

  ! Local variables
  integer :: xnum, ynum, ynum3, nzeros, nzerosnew, matsize
  integer :: bpres
  integer :: i, j, k, c, ivx, ivy, ipr, tmp
  integer :: iprleft, iprright, iprtop, iprbottom
  integer :: ivxleft, ivxright, ivxtop, ivxbottom
  integer :: ivxtopleft, ivxbottomleft, ivxtopright, ivxbottomright
  integer :: ivytop, ivybottom, ivyleft, ivyright
  integer :: ivytopleft, ivybottomleft, ivytopright, ivybottomright

  double precision, dimension (2) :: prfirst
  double precision :: prnorm, pscale, xstpavr, ystpavr

  logical :: notintern

  character*120 line

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------

  ! write(line,*) xnum, ynum
  ! k = mexPrintf(line//achar(10))

  !  Check for proper number of arguments.
  !      if(nrhs .ne. 10) then
  !         call mexErrMsgIdAndTxt ('MATLAB:mexFastStokesSolver:nInput','15 inputs required.')
  !      elseif(nlhs .gt. 6) then
  !         call mexErrMsgIdAndTxt ('MATLAB:mexFastStokesSolver:nOutput','Too many output arguments.')
  !      endif


  !  Get the size of the input arrays.
  ynum = mxGetM(prhs(2))
  xnum = mxGetN(prhs(2))

     !     write(line,*) xnum, ynum
     !     k = mexPrintf(line//achar(10))
  ! Arguments declarations

  !
  ! Staggered Grid for Multigrid
  !
  ! vx vx vx
  !
  ! vy +---vy---+---vy---+ vy
  ! | | |
  ! vx P vx P vx
  ! | | |
  ! vy +---vy---+---vy---+ vy
  ! | | |
  ! vx P vx P vx
  ! | | |
  ! vy +---vy---+---vy---+ vy
  !
  ! vx vx vx
  !
  ! Lines show basic grid
  ! Basic (density) nodes are shown with +
  ! Ghost nodes shown outside the basic grid
  ! are used for boundary conditions

  ! Boundary conditions
  ! Pressure boundary condition
  ! Pressure in first cell

  press_ptr = mxGetPr(prhs(1));
  call mxCopyPtrToReal8(press_ptr,prfirst,2)

  bpres=0
  prnorm=prfirst(2)
  ! Channel flow top->bottom
  if (prfirst(1)==1) then
    bpres=1
    prnorm=prfirst(2)
  end if

  ! Input arrays
  etas_ptr = mxGetPr(prhs(2))
  allocate(etas(ynum,xnum))
  call mxCopyPtrToReal8(etas_ptr,etas,xnum*ynum)

  etan_ptr = mxGetPr(prhs(3))
  allocate(etan(ynum-1,xnum-1))
  call mxCopyPtrToReal8(etan_ptr,etan,(xnum-1)*(ynum-1))

  gridx_ptr = mxGetPr(prhs(4))
  allocate(gridx(xnum))
  call mxCopyPtrToReal8(gridx_ptr,gridx,xnum)

  gridy_ptr = mxGetPr(prhs(5))
  allocate(gridy(ynum))
  call mxCopyPtrToReal8(gridy_ptr,gridy,ynum)

  RX_ptr = mxGetPr(prhs(6))
  RY_ptr = mxGetPr(prhs(7))
  RC_ptr = mxGetPr(prhs(8))
  allocate(RX(ynum+1,xnum),RY(ynum,xnum+1),RC(ynum-1,xnum-1))
  call mxCopyPtrToReal8(RX_ptr,RX,(ynum+1)*xnum)
  call mxCopyPtrToReal8(RY_ptr,RY,ynum*(xnum+1))
  call mxCopyPtrToReal8(RC_ptr,RC,(ynum-1)*(xnum-1))

  ! Input Velocity boundary conditions
  allocate(btop(xnum+1,4),bbottom(xnum+1,4))
  allocate(bleft(ynum+1,4),bright(ynum+1,4))
  allocate(btopx(xnum+1,2),btopy(xnum+1,2),bbottomx(xnum+1,2),bbottomy(xnum+1,2))
  allocate(bleftx(ynum+1,2),blefty(ynum+1,2),brightx(ynum+1,2),brighty(ynum+1,2))


  bc_ptr = mxGetPr(prhs(9))
  call mxCopyPtrToReal8(bc_ptr,bleft,4*(ynum+1))
  bleftx(:,1)=bleft(:,1)
  bleftx(:,2)=bleft(:,2)
  blefty(:,1)=bleft(:,3)
  blefty(:,2)=bleft(:,4)

  bc_ptr = mxGetPr(prhs(10))
  call mxCopyPtrToReal8(bc_ptr,bright,4*(ynum+1))
  brightx(:,1)=bright(:,1)
  brightx(:,2)=bright(:,2)
  brighty(:,1)=bright(:,3)
  brighty(:,2)=bright(:,4)


  bc_ptr = mxGetPr(prhs(11))
  call mxCopyPtrToReal8(bc_ptr,btop,4*(xnum+1))
  btopx(:,1) = btop(:,1);
  btopx(:,2) = btop(:,2);
  btopy(:,1) = btop(:,3);
  btopy(:,2) = btop(:,4);

  bc_ptr = mxGetPr(prhs(12))
  call mxCopyPtrToReal8(bc_ptr,bbottom,4*(xnum+1))
  bbottomx(:,1)=bbottom(:,1)
  bbottomx(:,2)=bbottom(:,2)
  bbottomy(:,1)=bbottom(:,3)
  bbottomy(:,2)=bbottom(:,4)

  deallocate(btop,bbottom,bleft,bright)

  ! Prescribed internal velocity condition
  bc_ptr = mxGetPr(prhs(13))
  call mxCopyPtrToReal8(bc_ptr,bintern,8*4)

  bc_ptr = mxGetPr(prhs(14))
  call mxCopyPtrToReal8(bc_ptr,pscale,1)

  ! Computing grid steps for basic nodes
  allocate(xstp(xnum-1),ystp(ynum-1))
  xstp = 0.0d0
  ystp = 0.0d0

  do i=1,xnum-1,1
    xstp(i)=gridx(i+1)-gridx(i)
  end do
  do i=1,ynum-1,1
    ystp(i)=gridy(i+1)-gridy(i)
  end do

  ! Computing grid steps for vx and vy nodes
  allocate(xstpc(xnum),ystpc(ynum))
  xstpc = 0.0d0
  ystpc = 0.0d0

  ! First and last steps (for ghost nodes)
  xstpc(1)=xstp(1)
  xstpc(xnum)=xstp(xnum-1)
  ystpc(1)=ystp(1)
  ystpc(ynum)=ystp(ynum-1)
  do i=2,xnum-1,1
    xstpc(i)=(gridx(i+1)-gridx(i-1))/2
  end do
  do i=2,ynum-1,1
    ystpc(i)=(gridy(i+1)-gridy(i-1))/2
  end do

  ! Average x and y steps
  xstpavr=(gridx(xnum)-gridx(1))/(xnum-1)
  ystpavr=(gridy(ynum)-gridy(1))/(ynum-1)

  ! Koefficient for scaling pressure
  !    pscale=2.0d0*etan(1,1)/(xstpavr+ystpavr)

  ! Horizontal shift index
  ynum3=(ynum-1)*3

  ! Creating output matrix indeces for lhs sparse matrix
  nzeros = xnum*ynum*26
  allocate(Im(nzeros),Jm(nzeros),Lm(nzeros))
  Im = 0
  Jm = 0
  Lm = 0.0d0

  ! Allocate rhs matrix
  matsize = 3*(xnum-1)*(ynum-1)
  allocate(Rm(matsize))
  Rm = 0.0d0

  c = 0
  ! Solving of Stokes and continuity equations on nodes

  do i=1,ynum-1,1
    do j=1,xnum-1,1
      ! Indexes for P,vx,vy
      ivx=((j-1)*(ynum-1)+(i-1))*3+1
      ivy=ivx+1
      ipr=ivx+2

      ! Increase the number of non-zero elements, grow the arrays
      if (c == nzeros) then
        nzerosnew = int(dble(nzeros)*1.25d0)
        allocate(tmpI(nzerosnew),tmpJ(nzerosnew),tmpL(nzerosnew))
        tmpI = 0
        tmpJ = 0
        tmpL = 0.0d0
        tmpI(1:nzeros) = Im
        tmpJ(1:nzeros) = Jm
        tmpL(1:nzeros) = Lm
        deallocate(Im,Jm, Lm)
        nzeros = nzerosnew
        allocate(Im(nzeros),Jm(nzeros),Lm(nzeros))
        Im = tmpI
        Jm = tmpJ
        Lm = tmpL
        deallocate(tmpI,tmpJ,tmpL)
      end if

      ! Logical test for internal boundaries
      notintern = (j/=bintern(1,1) .or. i<bintern(2,1) .or. i>bintern(3,1))
      notintern = notintern .and. (j/=bintern(1,2) .or. i<bintern(2,2) .or. i>bintern(3,2))
      notintern = notintern .and. (j/=bintern(1,3) .or. i<bintern(2,3) .or. i>bintern(3,3))
      notintern = notintern .and. (j/=bintern(1,4) .or. i<bintern(2,4) .or. i>bintern(3,4))

      ! x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy-dP/dx=RX
      if (j<xnum-1 .and. notintern) then
        ! x-Stokes equation stensil
        ! +----------------------+----------------------+
        ! | | |
        ! | | |
        ! | vx(i-1,j) ------------
        ! | | |
        ! | | |
        ! +-----vy(i-1,j)---etas(i,j+1)---vy(i-1,j+1)---- ystpc(i)--------
        ! | | |
        ! | | |
        ! vx(i,j-1) pr(i,j) vx(i,j) P(i,j+1) vx(i,j+1)----- ystp(i)
        ! | etan(i,j) | etan(i,j+1) |
        ! | | |
        ! +------vy(i,j)---etas(i+1,j+1)---vy(i,j+1)----- ystpc(i+1)------
        ! | | |
        ! | | |
        ! | vx(i+1,j) -----------
        ! | | |
        ! | | |
        ! +----------------------+----------------------+
        ! | xstp(j) | xstp(j+1) |
        ! | xstpc(j+1) |
        ! Right Part
        Rm(ivx)=RX(i+1,j+1)

        ! Computing Current x-Stokes coefficients

        ! Central Vx node
        c = c + 1
        Im(c) = ivx
        Jm(c) = ivx
        Lm(c) = -2*(etan(i,j+1)/xstp(j+1)+etan(i,j)/xstp(j))/xstpc(j+1)-(etas(i+1,j+1)/ystpc(i+1)+etas(i,j+1)/ystpc(i))/ystp(i)
        k = c

        ! Left Vx node
        if (j>1) then
          ivxleft=ivx-ynum3
          c = c + 1
          Im(c) = ivx
          Jm(c) = ivxleft
          Lm(c) = 2.0d0*etan(i,j)/xstp(j)/xstpc(j+1)
        else ! Left Vx boundary condition
          Lm(k) = Lm(k) + bleftx(i+1,2)*2*etan(i,j)/xstp(j)/xstpc(j+1)
          Rm(ivx)=Rm(ivx)-bleftx(i+1,1)*2*etan(i,j)/xstp(j)/xstpc(j+1)
        end if
        ! Right Vx node
        if (j<xnum-2) then
          ivxright=ivx+ynum3
          c = c + 1
          Im(c) = ivx
          Jm(c) = ivxright
          Lm(c) = 2.0d0*etan(i,j+1)/xstp(j+1)/xstpc(j+1)
        else
          Lm(k) = Lm(k) + brightx(i+1,2)*2*etan(i,j+1)/xstp(j+1)/xstpc(j+1)
          Rm(ivx)=Rm(ivx)-brightx(i+1,1)*2*etan(i,j+1)/xstp(j+1)/xstpc(j+1)
        end if
        ! Top Vx node
        if (i>1) then
          ivxtop=ivx-3
          c = c + 1
          Im(c) = ivx
          Jm(c) = ivxtop
          Lm(c) = etas(i,j+1)/ystpc(i)/ystp(i)
        else
          Lm(k) = Lm(k) + btopx(j+1,2)*etas(i,j+1)/ystpc(i)/ystp(i)
          Rm(ivx)=Rm(ivx)-btopx(j+1,1)*etas(i,j+1)/ystpc(i)/ystp(i)
        end if
        ! Bottom Vx node
        if (i<ynum-1) then
          ivxbottom=ivx+3
          c = c + 1
          Im(c) = ivx
          Jm(c) = ivxbottom
          Lm(c) = etas(i+1,j+1)/ystpc(i+1)/ystp(i)
        else
          Lm(k) = Lm(k) + bbottomx(j+1,2)*etas(i+1,j+1)/ystpc(i+1)/ystp(i)
          Rm(ivx)=Rm(ivx)-bbottomx(j+1,1)*etas(i+1,j+1)/ystpc(i+1)/ystp(i)
        end if
        ! Top Left Vy node
        if (i>1) then
          ivytopleft=ivx-3+1
          c = c + 1
          Im(c) = ivx
          Jm(c) = ivytopleft
          Lm(c) = etas(i,j+1)/xstpc(j+1)/ystp(i)
        else
          ivybottomleft=ivx+1
          c = c + 1
          Im(c) = ivx
          Jm(c) = ivybottomleft
          Lm(c) = btopy(j+1,2)*etas(i,j+1)/xstpc(j+1)/ystp(i)
          Rm(ivx)=Rm(ivx)-btopy(j+1,1)*etas(i,j+1)/xstpc(j+1)/ystp(i)
        end if
        ! Top Right Vy node
        if (i>1) then
          ivytopright=ivx-3+1+ynum3
          c = c + 1
          Im(c) = ivx
          Jm(c) = ivytopright
          Lm(c) = -etas(i,j+1)/xstpc(j+1)/ystp(i)
        else
          ivybottomright=ivx+1+ynum3
          c = c + 1
          Im(c) = ivx
          Jm(c) = ivybottomright
          Lm(c) = -btopy(j+2,2)*etas(i,j+1)/xstpc(j+1)/ystp(i)
          Rm(ivx)=Rm(ivx)+btopy(j+2,1)*etas(i,j+1)/xstpc(j+1)/ystp(i)
        end if
        ! Bottom Left Vy node
        if (i<ynum-1) then
          ivybottomleft=ivx+1
          if (i>1) then
            c = c + 1
            Im(c) = ivx
            Jm(c) = ivybottomleft
            Lm(c) = -etas(i+1,j+1)/xstpc(j+1)/ystp(i)
          else
            !k = find(Im ==ivx & Jm == ivybottomleft,1)
            where (Im ==ivx .and. Jm == ivybottomleft) Lm = Lm - etas(i+1,j+1)/xstpc(j+1)/ystp(i)
          end if
        else
          ivytopleft=ivx-3+1
          !k = find(Im == ivx & Jm == ivytopleft,1)
          where(Im == ivx .and. Jm == ivytopleft) Lm = Lm - bbottomy(j+1,2)*etas(i+1,j+1)/xstpc(j+1)/ystp(i)
          Rm(ivx)=Rm(ivx)+bbottomy(j+1,1)*etas(i+1,j+1)/xstpc(j+1)/ystp(i)
        end if
        ! Bottom Right Vy node
        if (i<ynum-1) then
          ivybottomright=ivx+1+ynum3
          if (i>1) then
            c = c + 1
            Im(c) = ivx
            Jm(c) = ivybottomright
            Lm(c) = etas(i+1,j+1)/xstpc(j+1)/ystp(i)
          else
            !k = find(Im ==ivx & Jm == ivybottomright,1)
            where(Im ==ivx .and. Jm == ivybottomright) Lm = Lm + etas(i+1,j+1)/xstpc(j+1)/ystp(i)
          end if
        else
          ivytopright=ivx-3+1+ynum3
          !k = find(Im ==ivx & Jm == ivytopright,1)
          where(Im ==ivx .and. Jm == ivytopright) Lm = Lm + bbottomy(j+2,2)*etas(i+1,j+1)/xstpc(j+1)/ystp(i)
          Rm(ivx) = Rm(ivx) - bbottomy(j+2,1)*etas(i+1,j+1)/xstpc(j+1)/ystp(i)
        end if
        ! Left P node
        iprleft=ivx+2
        c = c + 1
        Im(c) = ivx
        Jm(c) = iprleft
        Lm(c) = pscale/xstpc(j+1)
        ! Right P node
        iprright=ivx+2+ynum3
        c = c + 1
        Im(c) = ivx
        Jm(c) = iprright
        Lm(c) = -pscale/xstpc(j+1)

      ! Ghost Vx_parameter=0 used for numbering, internal prescribed horizontal velocity
      else
        c = c + 1
        Im(c) = ivx
        Jm(c) = ivx
        Lm(c) = 2*pscale/(xstpavr+ystpavr)
        if (notintern) then
          Rm(ivx)=0
        elseif (j==bintern(1,1)) then
          ! Internal prescribed horizontal velocity
          Rm(ivx)=2*pscale/(xstpavr+ystpavr)*bintern(4,1)
        elseif (j==bintern(1,2)) then
          ! Internal prescribed horizontal velocity
          Rm(ivx)=2*pscale/(xstpavr+ystpavr)*bintern(4,2)
        elseif (j==bintern(1,3)) then
          ! Internal prescribed horizontal velocity
          Rm(ivx)=2*pscale/(xstpavr+ystpavr)*bintern(4,3)
        elseif (j==bintern(1,4)) then
          ! Internal prescribed horizontal velocity
          Rm(ivx)=2*pscale/(xstpavr+ystpavr)*bintern(4,4)
        end if

      end if

      ! Logical test for internal boundaries
      notintern = (j/=bintern(5,1) .or. i<bintern(6,1) .or. i>bintern(7,1))
      notintern = notintern .and. (j/=bintern(5,2) .or. i<bintern(6,2) .or. i>bintern(7,2))
      notintern = notintern .and. (j/=bintern(5,3) .or. i<bintern(6,3) .or. i>bintern(7,3))
      notintern = notintern .and. (j/=bintern(5,4) .or. i<bintern(6,4) .or. i>bintern(7,4))

      ! y-Stokes equation dSIGMAyy/dy+dSIGMAyx/dx-dP/dy=RY
      if (i<ynum-1 .and. notintern) then
        ! y-Stokes equation stensil
        ! +-------------------- -+-------vy(i-1,j)------+----------------------+-----
        ! | | | |
        ! | | | |
        ! | vx(i,j-1) P(i,j) vx(i,j) |ystp(i)-------
        ! | | etan(i,j) | |
        ! | | | |
        ! +-----vy(i,j-1)---etas(i+1,j)---vy(i,j)--etas(i+1,j+1)---vy(i,j+1)---+----- ystpc(i+1)
        ! | | | |
        ! | | | |
        ! | vx(i+1,j-1) P(i+1,j) vx(i+1,j) |ystp(i+1)------
        ! | | etan(i+1,j) | |
        ! | | | |
        ! +----------------------+-------vy(i+1,j)------+----------------------+-----
        ! | xstpc(j) | xstpc(j+1) |
        ! | xstp(j) |
        !
        ! Right Part
        Rm(ivy)=RY(i+1,j+1)
        ! Computing Current y-Stokes coefficients
        ! Central Vy node
        c = c + 1
        Im(c) = ivy
        Jm(c) = ivy
        Lm(c) = -2*(etan(i+1,j)/ystp(i+1)+etan(i,j)/ystp(i))/ystpc(i+1)-(etas(i+1,j+1)/xstpc(j+1)+etas(i+1,j)/xstpc(j))/xstp(j)
        k = c
        ! Top Vy node
        if (i>1) then
          ivytop=ivy-3
          c = c + 1
          Im(c) = ivy
          Jm(c) = ivytop
          Lm(c)=2.0d0*etan(i,j)/ystp(i)/ystpc(i+1)
        else
          ! Add boundary condition for the top Vy node
          Lm(k) = Lm(k) + btopy(j+1,2)*2.0d0*etan(i,j)/ystp(i)/ystpc(i+1)
          Rm(ivy) = Rm(ivy)-btopy(j+1,1)*2.0d0*etan(i,j)/ystp(i)/ystpc(i+1)
        end if
        ! Bottom Vy node
        if (i<ynum-2) then
          ivybottom=ivy+3
          c = c + 1
          Im(c) = ivy
          Jm(c) = ivybottom
          Lm(c) = 2*etan(i+1,j)/ystp(i+1)/ystpc(i+1)
        else
          ! Add boundary condition for the bottom Vy node
          Lm(k) = Lm(k) + bbottomy(j+1,2)*2.0d0*etan(i+1,j)/ystp(i+1)/ystpc(i+1)
          Rm(ivy)=Rm(ivy)-bbottomy(j+1,1)*2.0d0*etan(i+1,j)/ystp(i+1)/ystpc(i+1)
        end if
        ! Left Vy node
        if (j>1) then
          ivyleft=ivy-ynum3
          c = c + 1
          Im(c) = ivy
          Jm(c) = ivyleft
          Lm(c) = etas(i+1,j)/xstpc(j)/xstp(j)
        else
          ! Add boundary condition for the left Vy node
          Lm(k) = Lm(k) + blefty(i+1,2)*etas(i+1,j)/xstpc(j)/xstp(j)
          Rm(ivy)=Rm(ivy) - blefty(i+1,1)*etas(i+1,j)/xstpc(j)/xstp(j)
        end if
        ! Right Vy node
        if (j<xnum-1) then
          ivyright = ivy+ynum3
          c = c + 1
          Im(c) = ivy
          Jm(c) = ivyright
          Lm(c) = etas(i+1,j+1)/xstpc(j+1)/xstp(j)
        else
          ! Add boundary condition for the right Vy node
          Lm(k) = Lm(k) + brighty(i+1,2)*etas(i+1,j+1)/xstpc(j+1)/xstp(j)
          Rm(ivy)=Rm(ivy)-brighty(i+1,1)*etas(i+1,j+1)/xstpc(j+1)/xstp(j)
        end if
        ! Top left Vx node
        if (j>1) then
          ivxtopleft=ivy-1-ynum3
          c = c + 1
          Im(c) = ivy
          Jm(c) = ivxtopleft
          Lm(c) = etas(i+1,j)/ystpc(i+1)/xstp(j)
        else
          ivxtopright=ivy-1
          c = c + 1
          Im(c) = ivy
          Jm(c) = ivxtopright
          Lm(c) = bleftx(i+1,2)*etas(i+1,j)/ystpc(i+1)/xstp(j)
          Rm(ivy)=Rm(ivy)-bleftx(i+1,1)*etas(i+1,j)/ystpc(i+1)/xstp(j)
        end if
        ! Bottom left Vx node
        if (j>1) then
          ivxbottomleft=ivy-1+3-ynum3
          c = c + 1
          Im(c) = ivy
          Jm(c) = ivxbottomleft
          Lm(c) = -etas(i+1,j)/ystpc(i+1)/xstp(j)
        else
          ivxbottomright=ivy-1+3
          c = c + 1
          Im(c) = ivy
          Jm(c) = ivxbottomright
          Lm(c) = -bleftx(i+2,2)*etas(i+1,j)/ystpc(i+1)/xstp(j)
          Rm(ivy)=Rm(ivy)+bleftx(i+2,1)*etas(i+1,j)/ystpc(i+1)/xstp(j)
        end if
        ! Top right Vx node
        if (j<xnum-1) then
          ivxtopright=ivy-1
          if (j>1) then
            c = c + 1
            Im(c) = ivy
            Jm(c) = ivxtopright
            Lm(c)=-etas(i+1,j+1)/ystpc(i+1)/xstp(j)
          else
            !k = find(Im == ivy & Jm == ivxtopright,1)
            where(Im == ivy .and. Jm == ivxtopright) Lm = Lm-etas(i+1,j+1)/ystpc(i+1)/xstp(j)
          end if
        else
          ivxtopleft=ivy-1-ynum3
          !k = find(Im == ivy & Jm == ivxtopleft,1)
          where(Im == ivy .and. Jm == ivxtopleft) Lm = Lm - brightx(i+1,2)*etas(i+1,j+1)/ystpc(i+1)/xstp(j)
          Rm(ivy)=Rm(ivy)+brightx(i+1,1)*etas(i+1,j+1)/ystpc(i+1)/xstp(j)
        end if
        ! Bottom right Vx node
        if (j<xnum-1) then
          ivxbottomright=ivy-1+3
          if (j>1) then
            c = c + 1
            Im(c) = ivy
            Jm(c) = ivxbottomright
            Lm(c) = etas(i+1,j+1)/ystpc(i+1)/xstp(j)
          else
            !k = find(Im == ivy & Jm == ivxbottomright,1)
            where(Im == ivy .and. Jm == ivxbottomright) Lm = Lm + etas(i+1,j+1)/ystpc(i+1)/xstp(j)
          end if
        else
          ivxbottomleft=ivy-1+3-ynum3
          !k = find(Im == ivy & Jm == ivxbottomleft,1)
          where(Im == ivy .and. Jm == ivxbottomleft) Lm = Lm + brightx(i+2,2)*etas(i+1,j+1)/ystpc(i+1)/xstp(j)
          Rm(ivy)=Rm(ivy)-brightx(i+2,1)*etas(i+1,j+1)/ystpc(i+1)/xstp(j)
        end if
        ! Top P node
        iprtop=ivy+1
        c = c + 1
        Im(c) = ivy
        Jm(c) = iprtop
        Lm(c)=pscale/ystpc(i+1)
        ! Bottom P node
        iprbottom=ivy+1+3
        c = c + 1
        Im(c) = ivy
        Jm(c) = iprbottom
        Lm(c)=-pscale/ystpc(i+1)

      ! Ghost Vy_parameter=0 used for numbering
      else
        c = c + 1
        Im(c) = ivy
        Jm(c) = ivy
        Lm(c) = 2*pscale/(xstpavr+ystpavr)
        if (notintern) then
          Rm(ivy)=0
        elseif (j==bintern(5,1)) then
          ! Internal prescribed vertical velocity
          Rm(ivy)=2*pscale/(xstpavr+ystpavr)*bintern(8,1)
        elseif (j==bintern(5,2)) then
          ! Internal prescribed vertical velocity
          Rm(ivy)=2*pscale/(xstpavr+ystpavr)*bintern(8,2)
        elseif (j==bintern(5,3)) then
          ! Internal prescribed vertical velocity
          Rm(ivy)=2*pscale/(xstpavr+ystpavr)*bintern(8,3)
        elseif (j==bintern(5,4)) then
          ! Internal prescribed vertical velocity
          Rm(ivy)=2*pscale/(xstpavr+ystpavr)*bintern(8,4)
        end if
      end if


      ! Continuity equation dvx/dx+dvy/dy=RC
      if ( ((j>1 .or. i>1) .and. bpres==0) .or. (i>1 .and. i<ynum-1 .and. bpres==1) .or. (j>1 .and. j<xnum-1 .and. bpres==2) ) then
        ! Continuity equation stensil
        ! +-----vy(i-1,j)--------+--------
        ! | |
        ! | |
        ! vx(i,j-1) pr(i,j) vx(i,j) ystp(i)
        ! | |
        ! | |
        ! +------vy(i,j)---------+--------
        ! | xstp(j) |
        !
        ! Right Part
        Rm(ipr)=RC(i,j)
        ! Computing Current Continuity coefficients
        ! Left Vx node
        if (j>1) then
          ivxleft=ipr-2-ynum3
          c = c + 1
          Im(c) = ipr
          Jm(c) = ivxleft
          Lm(c)=-pscale/xstp(j)
          ! Add boundary condition for the right Vx node
          if (j==xnum-1) then
            Lm(c) = Lm(c)+brightx(i+1,2)*pscale/xstp(j)
            Rm(ipr)=Rm(ipr)-brightx(i+1,1)*pscale/xstp(j)
          end if
        end if
        ! Right Vx node
        if (j<xnum-1) then
          ivxright=ipr-2
          c = c + 1
          Im(c) = ipr
          Jm(c) = ivxright
          Lm(c)=pscale/xstp(j)
          ! Add boundary condition for the left Vx node
          if (j==1) then
            Lm(c) = Lm(c)-bleftx(i+1,2)*pscale/xstp(j)
            Rm(ipr)=Rm(ipr)+bleftx(i+1,1)*pscale/xstp(j)
          end if
        end if
        ! Top Vy node
        if (i>1) then
          ivytop=ipr-1-3
          c = c + 1
          Im(c) = ipr
          Jm(c) = ivytop
          Lm(c)=-pscale/ystp(i)
          ! Add boundary condition for the bottom Vy node
          if (i==ynum-1) then
            Lm(c) = Lm(c)+bbottomy(j+1,2)*pscale/ystp(i)
            Rm(ipr)=Rm(ipr)-bbottomy(j+1,1)*pscale/ystp(i)
          end if
        end if
        ! Bottom Vy node
        if (i<ynum-1) then
          ivybottom=ipr-1
          c = c + 1
          Im(c) = ipr
          Jm(c) = ivybottom
          Lm(c)=pscale/ystp(i)
          ! Add boundary condition for the top Vy node
          if (i==1) then
            Lm(c) = Lm(c)-btopy(j+1,2)*pscale/ystp(i)
            Rm(ipr)=Rm(ipr)+btopy(j+1,1)*pscale/ystp(i)
          end if
        end if

      ! Pressure definition for the boundary condition regions
      else
        c = c + 1
        Im(c) = ipr
        Jm(c) = ipr
        ! Pressure definition in one cell
        if (bpres==0) then
          Lm(c)=2.0d0*pscale/(xstpavr+ystpavr)
          Rm(ipr)=2.0d0*prnorm/(xstpavr+ystpavr)
        end if
        ! Pressure definition at the top and bottom
        if (bpres==1) then
          Lm(c)=2.0d0*pscale/(xstpavr+ystpavr)
          if (i==1) then
            Rm(ipr)=2.0d0*prnorm/(xstpavr+ystpavr)
          else
            Rm(ipr)=0.0d0
          end if
        end if
        ! Pressure definition at the left and right
        if (bpres==2) then
          Lm(c)=2.0d0*pscale/(xstpavr+ystpavr)
          if (j==1) then
            Rm(ipr)=2.0d0*prnorm/(xstpavr+ystpavr)
          else
            Rm(ipr)=0.0d0
          end if
        end if
      end if

    end do
  end do

  ! Total final nonzero elements
  nzeros = c

  ! Prepare for output
  plhs(1) = mxCreateDoubleMatrix(nzeros,1,0)
  plhs(2) = mxCreateDoubleMatrix(nzeros,1,0)
  plhs(3) = mxCreateDoubleMatrix(nzeros,1,0)
  plhs(4) = mxCreateDoubleMatrix((xnum-1)*(ynum-1)*3,1,0)

  call mxCopyReal8ToPtr(dble(Im(1:nzeros)),mxGetPr(plhs(1)),nzeros)
  call mxCopyReal8ToPtr(dble(Jm(1:nzeros)),mxGetPr(plhs(2)),nzeros)
  call mxCopyReal8ToPtr(Lm(1:nzeros),mxGetPr(plhs(3)),nzeros)
  call mxCopyReal8ToPtr(Rm,mxGetPr(plhs(4)),matsize)

  return
end subroutine
