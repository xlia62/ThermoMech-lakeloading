# include <fintrf.h>
!> 

! function [MX,MY,MXN,MYN,MCXN,MCYN,MSXX,MSXY] = movemarkersFast(MX,MY,MXN,MYN,MCXN,MCYN,MSXX,MSXY,...
! vx1,vy1,esp,gridx,gridy,gridcx,gridcy,markmove,timestep,xstp1,ystp1,xstpc1,ystpc1)

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
    mwPointer :: inmarker_ptr, ingrid_ptr
    mwPointer :: outmarker_ptr

    !     Arguments for computational routine:
    double precision, allocatable, dimension(:) :: MX, MY, MSXX, MSXY
    double precision, allocatable, dimension(:) :: tmpD
    integer, allocatable, dimension (:) :: MXN, MYN, MCXN, MCYN
    double precision, allocatable, dimension(:,:) :: vx1, vy1, esp
    double precision, allocatable, dimension(:) :: gridx,gridy, gridcx, gridcy, xstp1, ystp1, xstpc1, ystpc1

    ! Local variables
    integer :: marknum, xnum, ynum, markmove
    integer :: i, k, mm1, rk
    integer :: xn, yn, xnmin, ynmin

    double precision, dimension (4) :: vxm, vym, espm
    double precision :: timestep, xcur, ycur, dx, dy
    double precision :: msxxold, msxyold, tmp

    character*120 line

    ! Number of markers
    marknum = mxGetM(prhs(1))*mxGetN(prhs(1))

    allocate(MX(marknum),MY(marknum),MXN(marknum),MYN(marknum),MCXN(marknum),MCYN(marknum))
    allocate(MSXX(marknum),MSXY(marknum))

    ! Populate marker arrays
    inmarker_ptr = mxGetPr(prhs(1));
    call mxCopyPtrToReal8(inmarker_ptr,MX,marknum)
    inmarker_ptr = mxGetPr(prhs(2));
    call mxCopyPtrToReal8(inmarker_ptr,MY,marknum)
    inmarker_ptr = mxGetPr(prhs(3));

    allocate(tmpD(marknum))

    call mxCopyPtrToReal8(inmarker_ptr,tmpD,marknum)
    MXN = int(tmpD)
    inmarker_ptr = mxGetPr(prhs(4));
    call mxCopyPtrToReal8(inmarker_ptr,tmpD,marknum)
    MYN = int(tmpD)
    inmarker_ptr = mxGetPr(prhs(5));
    call mxCopyPtrToReal8(inmarker_ptr,tmpD,marknum)
    MCXN = int(tmpD)
    inmarker_ptr = mxGetPr(prhs(6));
    call mxCopyPtrToReal8(inmarker_ptr,tmpD,marknum)
    MCYN = int(tmpD)

    deallocate(tmpD)

    inmarker_ptr = mxGetPr(prhs(7));
    call mxCopyPtrToReal8(inmarker_ptr,MSXX,marknum)
    inmarker_ptr = mxGetPr(prhs(8));
    call mxCopyPtrToReal8(inmarker_ptr,MSXY,marknum)
    
    ! Populate grids
    ynum = mxGetM(prhs(11))
    xnum = mxGetN(prhs(11))

    allocate(vx1(ynum+1,xnum),vy1(ynum,xnum+1),esp(ynum,xnum))

    ingrid_ptr = mxGetPr(prhs(9));
    call mxCopyPtrToReal8(ingrid_ptr,vx1,xnum*(ynum+1))
    ingrid_ptr = mxGetPr(prhs(10));
    call mxCopyPtrToReal8(ingrid_ptr,vy1,(xnum+1)*ynum)
    ingrid_ptr = mxGetPr(prhs(11));
    call mxCopyPtrToReal8(ingrid_ptr,esp,xnum*ynum)

    allocate(gridx(xnum),gridy(ynum),gridcx(xnum+1),gridcy(ynum+1))

    ingrid_ptr = mxGetPr(prhs(12));
    call mxCopyPtrToReal8(ingrid_ptr,gridx,xnum)
    ingrid_ptr = mxGetPr(prhs(13));
    call mxCopyPtrToReal8(ingrid_ptr,gridy,ynum)

    ingrid_ptr = mxGetPr(prhs(14));
    call mxCopyPtrToReal8(ingrid_ptr,gridcx,xnum+1)
    ingrid_ptr = mxGetPr(prhs(15));
    call mxCopyPtrToReal8(ingrid_ptr,gridcy,ynum+1)

    ingrid_ptr = mxGetPr(prhs(16));
    call mxCopyPtrToReal8(ingrid_ptr,tmp,1)
    markmove = int(tmp)

    ingrid_ptr = mxGetPr(prhs(17));
    call mxCopyPtrToReal8(ingrid_ptr,timestep,1)

    allocate(xstp1(xnum-1),ystp1(ynum-1),xstpc1(xnum),ystpc1(ynum))

    ingrid_ptr = mxGetPr(prhs(18));
    call mxCopyPtrToReal8(ingrid_ptr,xstp1,xnum-1)
    ingrid_ptr = mxGetPr(prhs(19));
    call mxCopyPtrToReal8(ingrid_ptr,ystp1,ynum-1)
    ingrid_ptr = mxGetPr(prhs(20));
    call mxCopyPtrToReal8(ingrid_ptr,xstpc1,xnum)
    ingrid_ptr = mxGetPr(prhs(21));
    call mxCopyPtrToReal8(ingrid_ptr,ystpc1,ynum)

    do mm1=1,marknum

        xcur = MX(mm1)
        ycur = MY(mm1)

        xnmin = MXN(mm1)
        ynmin = MYN(mm1)

        ! Defining number of Runge-Kutta cycles
        do rk=1,markmove

            ! xn V(xn,yn)--------------------V(xn+1,yn)
            ! ? ^ ?
            ! ? ? ?
            ! ? dy ?
            ! ? ? ?
            ! ? v ?
            ! ?<----dx--->o Mrho(xm,ym) ?
            ! ? ?
            ! ? ?
            ! xn+1 V(xn,yn+1)-------------------V(xn+1,yn+1)

            ! Define indexes for upper left BASIC node in the cell where the marker is
            ! using bisection
            ! Load horizontal and vertical indexes

            ! Define indexes for upper left node in the Vx-cell where the marker is
            ! Horizontal Vx index
            xn=xnmin
            ! Vertical Vx index
            yn=ynmin
            if (ycur>gridcy(yn+1)) then
                yn=yn+1
            end if
            if (yn>ynum) then
                yn=ynum
            end if

            ! Define and check normalized distances from marker to the upper left VX-node;
            dx=(xcur-gridx(xn))/xstp1(xn)
            dy=(ycur-gridcy(yn))/ystpc1(yn)

            ! Calculate Marker velocity from four surrounding Vx nodes
            vxm(rk)=0.0d0
            vxm(rk)=vxm(rk)+(1.0-dx)*(1.0-dy)*vx1(yn,xn)
            vxm(rk)=vxm(rk)+(1.0-dx)*dy*vx1(yn+1,xn)
            vxm(rk)=vxm(rk)+dx*(1.0-dy)*vx1(yn,xn+1)
            vxm(rk)=vxm(rk)+dx*dy*vx1(yn+1,xn+1)

            ! Define indexes for upper left node in the VY-cell where the marker is
            ! Vertical Vy index
            yn=ynmin
            ! Horizontal Vy index
            xn=xnmin
            if (xcur>gridcx(xn+1)) then
                xn=xn+1
            end if
            if (xn>xnum) then
                xn=xnum
            end if

            ! Define and check normalized distances from marker to the upper left VX-node;
            dx=(xcur-gridcx(xn))/xstpc1(xn)
            dy=(ycur-gridy(yn))/ystp1(yn)

            ! Calculate Marker velocity from four surrounding nodes
            vym(rk)=0.0d0
            vym(rk)=vym(rk)+(1.0-dx)*(1.0-dy)*vy1(yn,xn)
            vym(rk)=vym(rk)+(1.0-dx)*dy*vy1(yn+1,xn)
            vym(rk)=vym(rk)+dx*(1.0-dy)*vy1(yn,xn+1)
            vym(rk)=vym(rk)+dx*dy*vy1(yn+1,xn+1)

            ! Define indexes for upper left node in the Espin cell where the marker is
            xn=xnmin
            yn=ynmin

            ! Define normalized distances from marker to the upper left node;
            dx=(xcur-gridx(xn))/xstp1(xn)
            dy=(ycur-gridy(yn))/ystp1(yn)

            ! Interpolate old Sxy stress for the marker
            espm(rk)=0.0d0
            espm(rk)=espm(rk)+(1.0-dx)*(1.0-dy)*esp(yn,xn)
            espm(rk)=espm(rk)+(1.0-dx)*dy*esp(yn+1,xn)
            espm(rk)=espm(rk)+dx*(1.0-dy)*esp(yn,xn+1)
            espm(rk)=espm(rk)+dx*dy*esp(yn+1,xn+1)

            ! Update coordinates for the next cycle
            if (rk<4) then
                if (rk<3) then
                    xcur=MX(mm1)+timestep/2*vxm(rk)
                    ycur=MY(mm1)+timestep/2*vym(rk)
                else
                    xcur=MX(mm1)+timestep*vxm(rk)
                    ycur=MY(mm1)+timestep*vym(rk)
                end if
            end if

            ! Update marker grid-index to next point in advection scheme
            ! Use a while scheme since it will take only a couple if any
            ! iterations to find the new index

            do while (xn < xnum-1 .and. xcur > gridx(xn+1)) ! Search right
                xn = xn + 1
            end do

            do while (xn > 1 .and. xcur < gridx(xn)) ! Search left
                xn = xn - 1
            end do

            do while (yn < ynum-1 .and. ycur > gridy(yn+1)) ! Search down
                yn = yn + 1
            end do

            do while (yn > 1 .and. ycur < gridy(yn)) ! Search up
                yn = yn - 1
            end do

            xnmin = xn
            ynmin = yn

        end do

        ! Recompute velocity and spin using 4-th order Runge_Kutta
        if (markmove==4) then
            vxm(1)=(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4))/6
            vym(1)=(vym(1)+2*vym(2)+2*vym(3)+vym(4))/6
            espm(1)=(espm(1)+2*espm(2)+2*espm(3)+espm(4))/6
        end if

        ! Displacing Marker according to its velocity
        xcur=MX(mm1)+timestep*vxm(1)
        ycur=MY(mm1)+timestep*vym(1)

        ! Update marker grid-index to next point in advection scheme
        ! Use a while scheme since it will take only a couple if any
        ! iterations to find the new index

        do while (xn < xnum-1 .and. xcur > gridx(xn+1)) ! Search right
            xn = xn + 1
        end do

        do while (xn > 1 .and. xcur < gridx(xn)) ! Search left
            xn = xn - 1
        end do

        do while (yn < ynum-1 .and. ycur > gridy(yn+1)) ! Search down
            yn = yn + 1
        end do

        do while (yn > 1 .and. ycur < gridy(yn)) ! Search up
            yn = yn - 1
        end do

        ! Update marker position and indeces
        MXN(mm1) = xn
        MYN(mm1) = yn
        MX(mm1) = xcur
        MY(mm1) = ycur

        ! Now define the index for the central nodes (Presure, SXX, etc..)
        if (MX(mm1) < gridcx(xn+1)) then
            xn = xn-1
        end if
        if (xn<1) then
            xn=1
        end if
        if (xn>xnum-2) then
            xn=xnum-2
        end if
        MCXN(mm1) = xn

        if (MY(mm1) < gridcy(yn+1)) then
            yn = yn-1
        end if
        if (yn<1) then
            yn = 1
        end if
        if (yn>ynum-2) then
            yn = ynum-2
        end if
        MCYN(mm1) = yn

        ! Rotate stress on marker according to its spin
        ! Compute amount of rotation from spin rate:
        ! Espin=1/2(dvy/dx-dvx/dy) i.e. positive for clockwise rotation
        ! (when x axis is directed rightward and y axis is directed downward)
        espm(1)=espm(1)*timestep
        ! Save old stresses
        msxxold=MSXX(mm1)
        msxyold=MSXY(mm1)
        ! SxyNEW=0.5(Sxx-Syy)*sin(2*Espin*dt)+Sxy*cos(2*Espin*dt)
        ! where Sxx-Syy=2Sxx
        MSXY(mm1)=msxxold*dsin(2*espm(1))+msxyold*dcos(2*espm(1))
        ! SxxNEW=Sxx*(cos(Espin*dt))^2+Syy*(sin(Espin*dt))^2-Sxy*sin(2*Espin*dt)
        ! where Sxx=-Syy
        MSXX(mm1)=msxxold*((dcos(espm(1)))**2-(dsin(espm(1)))**2)-msxyold*dsin(2*espm(1))

    end do

    ! Prepare for output
    do i = 1,8
        plhs(i) = mxCreateDoubleMatrix(marknum,1,0)
    end do

    call mxCopyReal8ToPtr(MX,mxGetPr(plhs(1)),marknum)
    call mxCopyReal8ToPtr(MY,mxGetPr(plhs(2)),marknum)
    call mxCopyReal8ToPtr(dble(MXN),mxGetPr(plhs(3)),marknum)
    call mxCopyReal8ToPtr(dble(MYN),mxGetPr(plhs(4)),marknum)
    call mxCopyReal8ToPtr(dble(MCXN),mxGetPr(plhs(5)),marknum)
    call mxCopyReal8ToPtr(dble(MCYN),mxGetPr(plhs(6)),marknum)
    call mxCopyReal8ToPtr(MSXX,mxGetPr(plhs(7)),marknum)
    call mxCopyReal8ToPtr(MSXY,mxGetPr(plhs(8)),marknum)

    return
end subroutine
