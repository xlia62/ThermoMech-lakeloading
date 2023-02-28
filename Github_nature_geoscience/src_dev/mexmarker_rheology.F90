# include <fintrf.h>
!======================================================================
#if 0
!     Input
!     (MSXX, MSXY, MFLOW(MI), MTK, MEXX, MEXY, MRAT, MMU(MI), MPR, siiyeld, plastic_yield, ttop, stressmin, etamin, etamax, RGAS, timestep, ntimestep)
!     Output
!     (METACUR, MMUCUR, MSXX, MSXY, plastic_yield)
#endif
!
subroutine mexFunction(nlhs,plhs,nrhs,prhs)

    use omp_lib
    
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
    mwPointer :: inmarker_ptr, in_ptr

    mwSize :: marknum, n, m
    
    !     Arguments for computational routine:
    double precision, allocatable, dimension(:) :: MSXX, MSXY, MTK, MEXX, MEXY, MRAT, MMU, MPR
    double precision, allocatable, dimension(:) :: METACUR, MMUCUR
    double precision, allocatable, dimension(:,:) :: MFLOW

    !
    double precision, allocatable, dimension(:) :: sii0, eii0, sii1, eta0, eii1, siiyeld
    double precision, allocatable, dimension(:) :: plawexp, xelvis1, sxx1new, sxy1new
    double precision, allocatable, dimension(:) :: plastic_yield

    ! Local variables
    integer :: ntimestep
    integer :: i, j, k, plawiter

    double precision :: ttop, stressmin, etamin, etamax, RGAS, timestep, tmp
    double precision :: eiicur, eiiold, siicur, siinew, siiold, sxxnew, sxynew, xelvis, tmp_eta
    double precision :: sii0cur, sii1cur, tmp_sxx, tmp_sxy
    
    integer :: nprocs, nthreads

    character*10 line
    
   ! write(line,*) 1, 2
    !k = mexPrintf(line//achar(10))
    
    !k = mexPrintf('Hello world'//char(10)//char(0))
    
    !-----------------------------------------------------------------------
    ! Number of markers
    marknum = mxGetM(prhs(1))*mxGetN(prhs(1))
    
    !k = mexPrintf('Point 1'//char(10)//char(0))
    
    allocate(MSXX(marknum),MSXY(marknum),MFLOW(marknum,5),MTK(marknum))
    allocate(MEXX(marknum),MEXY(marknum),MRAT(marknum))
    allocate(MMU(marknum),MPR(marknum),plastic_yield(marknum))    
    allocate(siiyeld(marknum))

    !k = mexPrintf('Point 2'//char(10)//char(0))

    ! Populate marker arrays: MSXX, MSXY, MFLOW, MTK, MEXX, MEXY, MRAT, MPL, MGII, MMU, MPR

    inmarker_ptr = mxGetPr(prhs(1))
    call mxCopyPtrToReal8(inmarker_ptr,MSXX,marknum)    
    inmarker_ptr = mxGetPr(prhs(2))
    call mxCopyPtrToReal8(inmarker_ptr,MSXY,marknum)
    inmarker_ptr = mxGetPr(prhs(3))    
    call mxCopyPtrToReal8(inmarker_ptr,MFLOW,marknum*5)
    inmarker_ptr = mxGetPr(prhs(4))
    call mxCopyPtrToReal8(inmarker_ptr,MTK,marknum)
    inmarker_ptr = mxGetPr(prhs(5))
    call mxCopyPtrToReal8(inmarker_ptr,MEXX,marknum)    
    
    !k = mexPrintf('Point 3'//char(10)//char(0))

    
    inmarker_ptr = mxGetPr(prhs(6))
    call mxCopyPtrToReal8(inmarker_ptr,MEXY,marknum)
    inmarker_ptr = mxGetPr(prhs(7))
    call mxCopyPtrToReal8(inmarker_ptr,MRAT,marknum)
    inmarker_ptr = mxGetPr(prhs(8))
    call mxCopyPtrToReal8(inmarker_ptr,MMU,marknum)
    inmarker_ptr = mxGetPr(prhs(9))
    call mxCopyPtrToReal8(inmarker_ptr,MPR,marknum)
    inmarker_ptr = mxGetPr(prhs(10))
    call mxCopyPtrToReal8(inmarker_ptr,siiyeld,marknum)
    inmarker_ptr = mxGetPr(prhs(11))
    call mxCopyPtrToReal8(inmarker_ptr,plastic_yield,marknum)
    
    ! ttop, stressmin, etamin, etamax, RGAS, timestep, ntimestep
    n = 1
    in_ptr = mxGetPr(prhs(12))
    call mxCopyPtrToReal8(in_ptr,ttop,n)
    in_ptr = mxGetPr(prhs(13))
    call mxCopyPtrToReal8(in_ptr,stressmin,n)
    in_ptr = mxGetPr(prhs(14))
    call mxCopyPtrToReal8(in_ptr,etamin,n)
    in_ptr = mxGetPr(prhs(15))
    call mxCopyPtrToReal8(in_ptr,etamax,n)
    in_ptr = mxGetPr(prhs(16))
    call mxCopyPtrToReal8(in_ptr,RGAS,n)
    in_ptr = mxGetPr(prhs(17))
    call mxCopyPtrToReal8(in_ptr,timestep,n)
    in_ptr = mxGetPr(prhs(18))
    call mxCopyPtrToReal8(in_ptr,tmp,n)
    ntimestep = int(tmp)

    allocate(sii0(marknum),sii1(marknum))
    allocate(sxx1new(marknum),sxy1new(marknum))
    allocate(plawexp(marknum))
    allocate(eii0(marknum),eii1(marknum))
    allocate(eta0(marknum),xelvis1(marknum))

    ! Power-law: EPSILONii=AD*SIGMAii^n*exp[-(Ea+Va*P)/RT)
    ! Iterate for viscosity
    ! First viscosity value
    ! Compute and check old marker stress invariant in Pa
    sii0 = dsqrt(MSXX**2 + MSXY**2)

    ! Check old marker stress invariant (should be allways positive to be used in power law)
    where (sii0 < stressmin) sii0 = stressmin

    ! Check marker temperature
    plawexp = MTK
    where (plawexp < ttop) plawexp = ttop

    ! Compute exponential term:
    ! Ea is in J/mol(=1000*kJ/mol)
    ! Va is in J/Pa (=1e-6*cm^3)
    ! Cut if too big (at cold temperature);
    plawexp = (MFLOW(:,4)*1000 + MFLOW(:,5)*(1.0d-6)*MPR)/RGAS/plawexp
    where (plawexp > 150) plawexp = 150.0d0

    ! Compute AD*exp[-Ea/RT)
    plawexp = MFLOW(:,2)*dexp(-plawexp)

    ! Compute strain rate invariant from power law
    eii0 = plawexp*((1.0d-6)*sii0)**MFLOW(:,3)

    ! Compute effective viscosity
    eta0 = sii0/2/eii0

    ! Forcasting second invariant of future marker stress for given viscoelastic timestep
    xelvis1 = eta0/(MMU*timestep + eta0)
    sxx1new = MSXX*xelvis1 + 2.0d0*eta0*MEXX*MRAT*(1.0d0-xelvis1)
    sxy1new = MSXY*xelvis1 + 2.0d0*eta0*MEXY*MRAT*(1.0d0-xelvis1)
    sii1 = dsqrt(sxx1new**2 + sxy1new**2)

    ! Check new marker stress invariant (should be allways positive to be used in power law)
    where (sii1 < stressmin) sii1 = stressmin

    ! Compute strain rate invariant from power law
    eii1 = plawexp*((1.0d-6)*sii1)**MFLOW(:,3)

    ! Compute effective viscosity
    METACUR = sii1/2/eii1

!    call omp_set_num_threads(16)

    !$OMP PARALLEL DO SCHEDULE(DYNAMIC,32) PRIVATE(plawiter, tmp_eta, siicur, eiicur), &
    !$OMP& PRIVATE(xelvis, sxxnew, sxynew, siinew, sii0cur, sii1cur)
    
    do j=1,marknum
        ! Iteration counter
        plawiter = 0
        sii0cur = sii0(j)
        sii1cur = sii1(j)

        if (MFLOW(j,1) > 0.01d0 ) then
         
			tmp_eta = METACUR(j)
		
			do while (plawiter < 20 .and. dabs(sii1cur-sii0cur)>1.0d0)

				! Add iteration counter
				plawiter = plawiter + 1
				! Compute middle stress
				siicur = (sii0cur + sii1cur)*0.5d0
				! Compute strain rate invariant from power law
				eiicur = plawexp(j)*(1.0d-6*siicur)**MFLOW(j,3)
				! Compute effective viscosity
				tmp_eta = siicur/2.0d0/eiicur
				! Forcasting second invariant of future marker stress for given viscoelastic timestep
				xelvis = tmp_eta/(MMU(j)*timestep + tmp_eta)
				sxxnew = MSXX(j)*xelvis + 2.0d0*tmp_eta*MEXX(j)*MRAT(j)*(1.0d0-xelvis)
				sxynew = MSXY(j)*xelvis + 2.0d0*tmp_eta*MEXY(j)*MRAT(j)*(1.0d0-xelvis)
				siinew = dsqrt(sxxnew**2 + sxynew**2)
				! Changing bisection limits
				if ((sii0cur < sii1cur .and. siicur < siinew) .or. (sii0cur > sii1cur .and. siicur > siinew)) then
					sii0cur = siicur
				else
					sii1cur = siicur
				end if

			end do
		
			METACUR(j) = tmp_eta
			
        else
        
			METACUR(j) = MFLOW(j,2)
        
        end if
        
    end do
    !$OMP END PARALLEL DO
        
    ! Limiting viscosity for the power law
    where (METACUR < etamin) METACUR = etamin

    where (METACUR > etamax) METACUR = etamax

    ! Check if any plastic yeiding condition is present

    if (ntimestep > 1) then

        !$OMP PARALLEL DO SCHEDULE(DYNAMIC,32) PRIVATE(xelvis, sxxnew, sxynew, siinew, siiold, eiiold), &
        !$OMP& PRIVATE(tmp_eta, tmp_sxx, tmp_sxy)
        
        do j=1,marknum
              
                tmp_eta = METACUR(j)
                tmp_sxx = MSXX(j)
                tmp_sxy = MSXY(j)
                
                if (MFLOW(j,1) > 0.01d0 .and. siiyeld(j) > 0.01d0) then
					! Checking for plastic yeilding
					! Forcasting second invariant of future marker stress for given viscoelastic timestep
					xelvis = tmp_eta/(MMU(j)*timestep + tmp_eta)
					sxxnew = tmp_sxx*xelvis + 2.0d0*tmp_eta*MEXX(j)*MRAT(j)*(1.0d0-xelvis)
					sxynew = tmp_sxy*xelvis + 2.0d0*tmp_eta*MEXY(j)*MRAT(j)*(1.0d0-xelvis)
					siinew = dsqrt(sxxnew**2 + sxynew**2)
			
			        siiold = dsqrt(tmp_sxx**2+tmp_sxy**2)

					! Correcting rock properties for yeilding
					if (siiyeld(j) < siinew) then
						! Bringing marker stresses to yeilding stress
						tmp_sxx = tmp_sxx*siiyeld(j)/siiold
						tmp_sxy = tmp_sxy*siiyeld(j)/siiold
						! Bringing marker viscosity to yeilding stress
						eiiold = MRAT(j)*dsqrt(MEXX(j)**2+MEXY(j)**2)
						tmp_eta = siiyeld(j)/2.0d0/eiiold
						! Mark that plastic yeildind occur
						plastic_yield(j) = 1.0d0
					end if
					
                end if
                
                METACUR(j) = tmp_eta
                MSXX(j) = tmp_sxx
                MSXY(j) = tmp_sxy
                            
        end do
        !$OMP END PARALLEL DO
        
		! Limiting viscosity for the power law
		where (METACUR < etamin) METACUR = etamin

		where (METACUR > etamax) METACUR = etamax

		! Constant viscosity
		where (MFLOW(:,1) < 0.01d0) METACUR = MFLOW(:,2)
 
    end if
    
    ! Compute 1/MU values (MU is shear modulus)
    MMUCUR = 1.0d0/MMU
   
    ! Prepare for output
    n = 1
    plhs(1) = mxCreateDoubleMatrix(marknum,n,0)
    call mxCopyReal8ToPtr(METACUR,mxGetPr(plhs(1)),marknum)
    plhs(2) = mxCreateDoubleMatrix(marknum,n,0)
    call mxCopyReal8ToPtr(MMUCUR,mxGetPr(plhs(2)),marknum)
    plhs(3) = mxCreateDoubleMatrix(marknum,n,0)
    call mxCopyReal8ToPtr(MSXX,mxGetPr(plhs(3)),marknum)
    plhs(4) = mxCreateDoubleMatrix(marknum,n,0)
    call mxCopyReal8ToPtr(MSXY,mxGetPr(plhs(4)),marknum)
    plhs(5) = mxCreateDoubleMatrix(marknum,n,0)
    call mxCopyReal8ToPtr(plastic_yield,mxGetPr(plhs(5)),marknum)

    return

end subroutine  mexFunction
