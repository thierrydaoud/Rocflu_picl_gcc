!*********************************************************************
!* Illinois Open Source License                                      *
!*                                                                   *
!* University of Illinois/NCSA                                       * 
!* Open Source License                                               *
!*                                                                   *
!* Copyright@2008, University of Illinois.  All rights reserved.     *
!*                                                                   *
!*  Developed by:                                                    *
!*                                                                   *
!*     Center for Simulation of Advanced Rockets                     *
!*                                                                   *
!*     University of Illinois                                        *
!*                                                                   *
!*     www.csar.uiuc.edu                                             *
!*                                                                   *
!* Permission is hereby granted, free of charge, to any person       *
!* obtaining a copy of this software and associated documentation    *
!* files (the "Software"), to deal with the Software without         *
!* restriction, including without limitation the rights to use,      *
!* copy, modify, merge, publish, distribute, sublicense, and/or      *
!* sell copies of the Software, and to permit persons to whom the    *
!* Software is furnished to do so, subject to the following          *
!* conditions:                                                       *
!*                                                                   *
!*                                                                   *
!* @ Redistributions of source code must retain the above copyright  * 
!*   notice, this list of conditions and the following disclaimers.  *
!*                                                                   * 
!* @ Redistributions in binary form must reproduce the above         *
!*   copyright notice, this list of conditions and the following     *
!*   disclaimers in the documentation and/or other materials         *
!*   provided with the distribution.                                 *
!*                                                                   *
!* @ Neither the names of the Center for Simulation of Advanced      *
!*   Rockets, the University of Illinois, nor the names of its       *
!*   contributors may be used to endorse or promote products derived * 
!*   from this Software without specific prior written permission.   *
!*                                                                   *
!* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,   *
!* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES   *
!* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND          *
!* NONINFRINGEMENT.  IN NO EVENT SHALL THE CONTRIBUTORS OR           *
!* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       * 
!* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,   *
!* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE    *
!* USE OR OTHER DEALINGS WITH THE SOFTWARE.                          *
!*********************************************************************
!* Please acknowledge The University of Illinois Center for          *
!* Simulation of Advanced Rockets in works and publications          *
!* resulting from this software or its derivatives.                  *
!*********************************************************************
! ******************************************************************************
!
! Purpose: calculate solution at a new time level.
!
! Description: the governing equations are integrated in time using
!              the classical 4-stage Runge-Kutta method (4th-order in
!              time) in low-storage formulation.
!
! Input: regions = data of all regions.
!
! Output: regions%levels%mixt = new solution after one time step.
!
! Notes: none.
!
! ******************************************************************************
!
! $Id: RungeKuttaMP.F90,v 1.1.1.1 2015/01/23 22:57:50 tbanerjee Exp $
!
! Copyright: (c) 2003-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE RungeKuttaMP( regions )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal, ONLY     : t_global
  USE ModError
  USE ModParameters

  USE RFLU_ModBoundXvUtils
  USE RFLU_ModGridSpeedUtils, ONLY: RFLU_DescaleGridSpeeds, &
                                    RFLU_ScaleGridSpeeds, &
                                    RFLU_SetGridSpeedScaleFactor
  USE RFLU_ModMovingFrame, ONLY: RFLU_MVF_ComputeAcceleration, &
                                 RFLU_MVF_SetVelocity, &
                                 RFLU_MVF_UpdateBC
  USE RFLU_ModTimeZoom, ONLY: RFLU_TimeZoomDriver
  USE RFLU_ModMPI
  USE RFLU_ModNSCBC
  USE RFLU_ModGFM
  USE RFLU_ModRelatedPatches, ONLY: RFLU_RELP_TransformWrapper
  USE RFLU_ModTime, ONLY: RFLU_SetTimeRK

  USE ModInterfaces, ONLY: AfterUpdateMP, &
                           CellGradientsMP, ConvectiveFluxesMP, &
                           GlobalCommunicationMP, InitCommunicationMP, &
                           RKInitMP, RKUpdateMP, &
                           SourceTermsMP, UpdateBoundaryConditionsMP, &
                           UpdateDependentVarsMP, ViscousFluxesMP, &
                           ZeroDummyCellsMP, ZeroResidualsMP
  USE ModInterfaces, ONLY: RFLU_EquilibriumEulerian, &
                           RFLU_SetVarsContWrapper, &
                           RFLU_SetVarsDiscWrapper, & 
                           RFLU_UpdateBoundaryValues,&
                           RFLU_DecideWrite !BRAD added for picl     
#ifdef GENX
  USE RFLU_ModGENXTools, ONLY: RFLU_GENX_InitBFLAG
#endif
#ifdef SPEC
  USE SPEC_RFLU_ModPBA, ONLY: SPEC_RFLU_PBA_ProgramBurn, &
                              SPEC_RFLU_PBA_ReactionZone
#endif
#ifdef PLAG
  USE PLAG_RFLU_ModComm
#endif

#ifdef PICL
USE RFLU_ModConvertCv, ONLY: RFLU_ConvertCvCons2Prim, &
                               RFLU_ConvertCvPrim2Cons
#endif



#ifdef PICL
!DEC$ NOFREEFORM
!#include "/home/rahul.koneru/codes/Rocflu-ppiclF/ppiclF/source/PPICLF"
!#include
!"/home/rahul.koneru/codes/Rocflu-ppiclF/ppiclF_new/ppiclF/source/PPICLF"
#include "/home/t.daoud/Rocflu_picl_gcc/ppiclf/source/PPICLF_USER.h"
#include "/home/t.daoud/Rocflu_picl_gcc/ppiclf/source/PPICLF_STD.h"
!DEC$ FREEFORM
#endif



  IMPLICIT NONE

! ... parameters
  TYPE(t_region), POINTER :: regions(:)

! ... loop variables
  INTEGER :: iReg, iRegLocal, istage, icg

! ... local variables
  INTEGER :: flowModel

  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

#ifdef PICL
  LOGICAL :: doWrite      
  INTEGER(KIND=4) :: i,piclIO,nCells,lx,ly,lz
  INTEGER :: errorFlag      
  REAL(KIND=8) :: piclDtMin,piclCurrentTime,drudtMixt,drvdtMixt,drwdtMixt,energydotg
  REAL(KIND=8), DIMENSION(3) :: ug      
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: rhoF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: uxF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: uyF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: uzF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: csF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: vfP
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: dpxF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: dpyF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: dpzF
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: SDRX
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: SDRY
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: SDRZ
  REAL(KIND=8), DIMENSION(:,:,:), POINTER :: pGc 
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: rhsR        
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: pGcX 
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: pGcY
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: pGcZ
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: JFX
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFXCell
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: JFY
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFYCell
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: JFZ
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: JFZCell
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: PhiP
  REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: YTEMP


#endif




!******************************************************************************

  global => regions(1)%global

  CALL RegisterFunction( global,'RungeKuttaMP',__FILE__ )

! loop over stages and regions ================================================

  DO istage=1,global%nrkSteps

! ----- compute particle accelerations and velocities  ------------------------
    IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
      CALL RFLU_MVF_ComputeAcceleration(regions)
    
      DO iRegLocal=1,global%nRegionsLocal
        iReg = iRegLocal

        pRegion => regions(iReg)

        IF ( global%mvfAccFlag .EQV. .TRUE. ) THEN
          CALL RFLU_MVF_SetVelocity(pRegion)
        END IF ! global%mvfAccFlag

        CALL RFLU_MVF_UpdateBC(pRegion)
      END DO ! iRegLocal
    END IF ! global%mvFrameFlag

    DO iRegLocal=1,global%nRegionsLocal
      iReg = iRegLocal

! ----- set pointer and get models --------------------------------------------

        pRegion => regions(iReg)

        flowModel  = regions(iReg)%mixtInput%flowModel
        regions(iReg)%irkStep = istage
        regions(iReg)%dummyStep = .FALSE.

! ----- Set RK time -----------------------------------------------------------

        CALL RFLU_SetTimeRK(pRegion,iStage)


! ----- RFLU fill GENX incoming buffers ---------------------------------------

#ifdef GENX
	CALL RFLU_GENX_InitBFLAG(pRegion)
#endif
        CALL RFLU_UpdateBoundaryValues(regions(iReg),istage)

! ----- Scale grid speeds -----------------------------------------------------

        CALL RFLU_SetGridSpeedScaleFactor(pRegion)
        CALL RFLU_ScaleGridSpeeds(pRegion)

! ----- program burn, before RKInitMP -----------------------------------------
#ifdef SPEC        
! DEBUG: Manoj-PBA1D, Notes: To match Jianghui's implementation
!                         1. Move call to ProgramBurn to before update dependent vars
!                    changing it back to old implementation of Manoj
!                         2. Remove call to ReactionZone
! END DEBUG
        IF ( (global%specUsed .EQV. .TRUE.) .AND. &
             (global%pbaFlag .EQV. .TRUE.) ) THEN
! DEBUG: Manoj-PBA1D
          CALL SPEC_RFLU_PBA_ProgramBurn( pRegion )
! END DEBUG
        END IF ! pbaFlag
#endif

! ----- set ghost fluid solution ----------------------------------------------

! TEMPORARY: Manoj: GFM, need to put more thinking
        CALL RFLU_GFM_SetGhostFluid( pRegion )

! ----- store previous solution; set dissipation to zero ----------------------
        CALL RKInitMP( regions(iReg),istage )

! ----- compute cell-gradients for higher-order scheme ------------------------

        CALL CellGradientsMP( regions(iReg) )

! ----- compute viscous fluxes ------------------------------------------------

        IF ( flowModel == FLOW_NAVST ) THEN
          CALL ViscousFluxesMP( regions(iReg) )
        END IF ! flowModel

! ----- compute convective fluxes; form residual ------------------------------

        CALL ConvectiveFluxesMP( regions(iReg) )

! ----- zero residuals --------------------------------------------------------

        CALL ZeroResidualsMP(regions(iReg))

! ----- add source terms ------------------------------------------------------

        CALL SourceTermsMP( regions(iReg) )

! ----- zero residuals --------------------------------------------------------

        CALL ZeroResidualsMP(regions(iReg))

! ----- add Equilibrium Eulerian corrections ----------------------------------

        CALL RFLU_EquilibriumEulerian( pRegion )

! ----- zero out residuals in dummy cells -------------------------------------

        CALL ZeroDummyCellsMP( regions(iReg) )

! ----- Descale grid speeds -----------------------------------------------------
        CALL RFLU_DescaleGridSpeeds(pRegion)

    ENDDO    ! iReg

    IF(global%zoomFactor > 1) THEN
       CALL RFLU_TimeZoomDriver(regions)
    ENDIF ! global%zoomFactor    

    DO iRegLocal=1,global%nRegionsLocal
      iReg = iRegLocal

! ----- set Region pointer

        pRegion => regions(iReg)

!Most likly remove this call
! ----- update solution; sum up residuals -------------------------------------

!        CALL RKUpdateMP( regions(iReg),iReg,istage )

!PPICLF Integration
#ifdef PICL

     piclIO = 100000000
     piclDtMin = REAL(global%dtMin,8)
     piclCurrentTime = REAL(global%currentTime,8)
     doWrite = RFLU_DecideWrite(global)
!Figure out piclIO call, might need to look into timestepping
    IF ( (doWrite .EQV. .TRUE.)) piclIO = 1

!PARTICLE stuff possbile needed
!    CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)


!allocate arrays to send to picl
    nCells = pRegion%grid%nCells
    ALLOCATE(rhoF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uxF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uyF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uzF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(csF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(vfP(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpxF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error
    
    ALLOCATE(dpyF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpzF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRX(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRY(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(SDRZ(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(rhsR(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcX(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcY(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(pGcZ(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFX(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFXCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFY(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFYCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFZ(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(JFZCell(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(PhiP(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(YTEMP(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error


!Might need to update prim like plag does
pGc => pRegion%mixt%gradCell
!Fill arrays for interp field
    DO i = 1,pRegion%grid%nCells
!Zero out phip
        PhiP(i) = 0

        ug(XCOORD) = pRegion%mixt%cvOld(CV_MIXT_XMOM,i)&
                        /pRegion%mixt%cvOld(CV_MIXT_DENS,i)

        ug(YCOORD) = pRegion%mixt%cvOld(CV_MIXT_YMOM,i)&
                        /pRegion%mixt%cvOld(CV_MIXT_DENS,i)

        ug(ZCOORD) = pRegion%mixt%cvOld(CV_MIXT_ZMOM,i)&
                        /pRegion%mixt%cvOld(CV_MIXT_DENS,i)

        drudtMixt = -pRegion%mixt%rhs(CV_MIXT_XMOM,i)/pRegion%grid%vol(i) &
                  +pRegion%mixt%cvOld(CV_MIXT_DENS,i)*DOT_PRODUCT(ug,pGc(:,2,i)) &
                  +ug(XCOORD)*DOT_PRODUCT(ug,pGc(:,1,i))

        drvdtMixt = -pRegion%mixt%rhs(CV_MIXT_YMOM,i)/pRegion%grid%vol(i) &
                  +pRegion%mixt%cvOld(CV_MIXT_DENS,i)*DOT_PRODUCT(ug,pGc(:,3,i))&
                  +ug(YCOORD)*DOT_PRODUCT(ug,pGc(:,1,i))

        drwdtMixt = -pRegion%mixt%rhs(CV_MIXT_ZMOM,i)/pRegion%grid%vol(i) &
                  +pRegion%mixt%cvOld(CV_MIXT_DENS,i)*DOT_PRODUCT(ug,pGc(:,4,i))&
                  +ug(ZCOORD)*DOT_PRODUCT(ug,pGc(:,1,i))


       do lz=1,2
       do ly=1,2
       do lx=1,2 
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JPHIP,vfP(lx,ly,lz,i))
       !PhiP(i) = PhiP(i) +  (0.125*vfP(lx,ly,lz,i))*(0.0002**3)*(0.0006/0.0002) ! / pRegion%grid%vol(i)    
       ! Sam - generalize Brad's hardcode

       ! the factor of 0.0006/0.0002 was to compensate for making the grid look
       ! larger in the z direction to ppiclF so that the particles would not go
       ! through the thin layer of the wall. That's my best guess anyway.
       ! zpf_factor. This should be a global specified in the input file instead
       ! of this weird ad-hoc implementation - Sam
       PhiP(i) = PhiP(i) + (0.125*vfP(lx,ly,lz,i))*global%ppiclFInitZpfFactor*pRegion%grid%vol(i) ! / pRegion%grid%vol(i)    
       rhoF(lx,ly,lz,i) = pRegion%mixt%cvOld(CV_MIXT_DENS,i)
       uxF(lx,ly,lz,i) = pRegion%mixt%cvOld(CV_MIXT_XMOM,i) &
                        /pRegion%mixt%cvOld(CV_MIXT_DENS,i)
       uyF(lx,ly,lz,i) = pRegion%mixt%cvOld(CV_MIXT_YMOM,i) &
                        /pRegion%mixt%cvOld(CV_MIXT_DENS,i)
       uzF(lx,ly,lz,i) = pRegion%mixt%cv(CV_MIXT_ZMOM,i) &
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)
       csF(lx,ly,lz,i) = pRegion%mixt%dv(DV_MIXT_SOUN,i)
       dpxF(lx,ly,lz,i) = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_PRES,i)
       dpyF(lx,ly,lz,i) = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_PRES,i) 
       dpzF(lx,ly,lz,i) = pRegion%mixt%gradCell(ZCOORD,GRC_MIXT_PRES,i) 
!DEBUG TOP ISSUE
!        dpyF(lx,ly,lz,i) = dpyF(lx,ly,lz,i) *global%pi/6.0*(5.0E-6)**3         

 
       SDRX(lx,ly,lz,i) = drudtMixt 
       SDRY(lx,ly,lz,i) = drvdtMixt 
       SDRZ(lx,ly,lz,i) = drwdtMixt 
       rhsR(lx,ly,lz,i) = -pRegion%mixt%rhs(CV_MIXT_DENS,i)/pRegion%grid%vol(i) 
       pGcX(lx,ly,lz,i) = pGc(XCOORD,1,i)
       pGcY(lx,ly,lz,i) = pGc(YCOORD,1,i)
       pGcZ(lx,ly,lz,i) = pGc(ZCOORD,1,i)

!Debug Top Issue
!        if (dpyF(lx,ly,lz,i) .lt. -0.5 .or. dpyF(lx,ly,lz,i) .gt. 0.5) then
!                write(*,*) "TOP BAD GC i vol", &
!               pRegion%mixt%gradCell(YCOORD,GRC_MIXT_PRES,i),i,global%pi/6.0*(5.0E-6)**3
!        end if


       end do
       end do
       end do 
       
      !Dump back VolFrac
       PhiP(i) = Phip(i) / (pRegion%grid%vol(i))!*(0.0006/0.0002))!*2.0*global%pi/0.0001)
       if (PhiP(i) .gt. 0.6) PhiP(i) = 0.6  

!VOL Frac cap
       do lz=1,2
       do ly=1,2
       do lx=1,2 
             vfp(lx,ly,lz,i) = PhiP(i)      

       end do
       end do
       end do   

    END DO

!Interp field calls

DO i = 1,pRegion%grid%nCells


       do lz=1,2
       do ly=1,2
       do lx=1,2
        YTEMP(lx,ly,lz,i) = pRegion%grid%cofg(YCOORD,i)
       end do
       end do
       end do
        if ((pRegion%grid%cofg(YCOORD,i) .gt. 1.0) .or.&
                (pRegion%grid%cofg(YCOORD,i) .le. 0) )then
                ! Sam - not sure why this was here
                !write(*,*) "Bad YCOORD",i,pRegion%grid%cofg(YCOORD,i) 
        endif 
end do
 call ppiclf_solve_InterpFieldUser(PPICLF_R_JSPT,YTEMP)



!write(*,*) 'Pdtmin,Pcurr time', piclIO,piclDtMin,piclCurrentTime 

      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHOF,rhoF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUX,uxF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUY,uyF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUZ,uzF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPDX,dpxF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPDY,dpyF)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPDZ,dpzF)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JCS,csF)
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JPHIP,vfP)  
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JSDRX,SDRX)
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JSDRY,SDRY)
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JSDRZ,SDRZ)
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JRHSR,rhsR)  
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JPGCX,pGcX) 
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JPGCY,pGcY) 
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JPGCZ,pGcZ) 


!SOLVE
    CALL ppiclf_solve_IntegrateParticle(1,piclIO,piclDtMin,piclCurrentTime)


!FEED BACK TERM
!Fill arrays for interp field
if (.true.) then ! Sam - kill feedback to stick with 1 way coupling for now
    DO i = 1,pRegion%grid%nCells

        ug(XCOORD) = pRegion%mixt%cvOld(CV_MIXT_XMOM,i)&
                        /pRegion%mixt%cvOld(CV_MIXT_DENS,i)

        ug(YCOORD) = pRegion%mixt%cvOld(CV_MIXT_YMOM,i)&
                        /pRegion%mixt%cvOld(CV_MIXT_DENS,i)

        ug(ZCOORD) = pRegion%mixt%cvOld(CV_MIXT_ZMOM,i)&
                        /pRegion%mixt%cvOld(CV_MIXT_DENS,i)

       JFXCell(i) = 0.0 
       JFYCell(i) = 0.0 
       JFZCell(i) = 0.0 
       do lz=1,2
       do ly=1,2
       do lx=1,2 
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JFX,JFX(lx,ly,lz,i))  
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JFY,JFY(lx,ly,lz,i))      
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JFZ,JFZ(lx,ly,lz,i))      
       JFXCell(i) = JFXCell(i) + JFX(lx,ly,lz,i) ! / pRegion%grid%vol(i)    
       JFYCell(i) = JFYCell(i) + JFY(lx,ly,lz,i)  
       JFZCell(i) = JFZCell(i) + JFZ(lx,ly,lz,i)  
       end do
       end do
       end do 
        !JFXCell(i) = JFXCell(i) * 0.125 * (0.0002**3)*(0.0006/0.0002)
        !JFYCell(i) = JFYCell(i) * 0.125 * (0.0002**3)*(0.0006/0.0002)

        ! Sam - remove hardcode
        !JFXCell(i) = JFXCell(i) * 0.125 * pRegion%grid%vol(i) * &
        !             global%ppiclFInitZpfFactor !(0.0002**3)*(0.0006/0.0002)
        !JFYCell(i) = JFYCell(i) * 0.125 * pRegion%grid%vol(i) * &
        !             global%ppiclFInitZpfFactor !(0.0002**3)*(0.0006/0.0002)
        !JFZCell(i) = JFZCell(i) * 0.125 * pRegion%grid%vol(i) * &
        !             global%ppiclFInitZpfFactor !(0.0002**3)*(0.0006/0.0002)

        ! Sam - z factor is stupid, remove it
        JFXCell(i) = JFXCell(i) * 0.125 * pRegion%grid%vol(i)
        JFYCell(i) = JFYCell(i) * 0.125 * pRegion%grid%vol(i)
        JFZCell(i) = JFZCell(i) * 0.125 * pRegion%grid%vol(i)

        ! Sam - debugging log
!        if ((abs(JFXCell(i)) .gt. 1.0E-10) .or. &
!            (abs(JFYCell(i)) .gt. 1.0E-10) .or. &
!            (abs(JFZCell(i)) .gt. 1.0E-10)) then
!
!            WRITE(202, "(3(1x,E23.16))") JFXCell(i), JFYCell(i), JFZCell(i)
!        elseif (abs(JFZCell(i)) .gt. 1.0_RFREAL*10.0**-10) then
!            WRITE(203, *) JFZCell(i)
!        end if

        energydotg = JFXCell(i) * ug(1) + JFYCell(i) * ug(2) +&
                     JFZCell(i) * ug(3)

        pRegion%mixt%rhs(CV_MIXT_XMOM,i) &
                         = pRegion%mixt%rhs(CV_MIXT_XMOM,i) &
                         + JFXCell(i)
!        pRegion%mixt%rhs(CV_MIXT_ENER,icg) &
!                         = pRegion%mixt%rhs(CV_MIXT_ENER,icg) &
!                         + energydotg        
        
        pRegion%mixt%rhs(CV_MIXT_YMOM,i) &
                         = pRegion%mixt%rhs(CV_MIXT_YMOM,i) &
                         + JFYCell(i)

        pRegion%mixt%rhs(CV_MIXT_ZMOM,i) &
                         = pRegion%mixt%rhs(CV_MIXT_ZMOM,i) &
                         + JFZCell(i)

! Sam - I have no idea why icg is used here. It's not even defined.
!        pRegion%mixt%rhs(CV_MIXT_ENER,icg) &
!                         = pRegion%mixt%rhs(CV_MIXT_ENER,icg) &
!                         + energydotg


        pRegion%mixt%rhs(CV_MIXT_ENER,i) &
                         = pRegion%mixt%rhs(CV_MIXT_ENER,i) &
                         + energydotg


 !       if (JFXCell(i) .gt. 0) write(*,*) "Brad rhs stuff ",JFXCell(i)  


!For 3d will use later
!        region%mixt%rhs(CV_MIXT_XMOM:CV_MIXT_ZMOM,icg) &
!                         = region%mixt%rhs(CV_MIXT_XMOM:CV_MIXT_ZMOM,icg) &
!                         + contFac*forceTotal
!        region%mixt%rhs(CV_MIXT_ENER,icg) &
!                         = region%mixt%rhs(CV_MIXT_ENER,icg) &
!                         + contFac*energydotg
    END DO
end if

!
 DO i = 1,pRegion%grid%nCells
!zero out PhiP
       PhiP(i) = 0 
       do lz=1,2
       do ly=1,2
       do lx=1,2 
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JPHIP,vfP(lx,ly,lz,i))
       ! Sam - generalize Brad's hard code
       !PhiP(i) = PhiP(i) +  (0.125*vfP(lx,ly,lz,i))*(0.0002**3)*(0.0006/0.0002) ! / pRegion%grid%vol(i) 

       PhiP(i) = PhiP(i) + (0.125*vfP(lx,ly,lz,i))*global%ppiclfInitZpfFactor*pRegion%grid%vol(i) ! / pRegion%grid%vol(i)    
       end do
       end do
       end do 
       !Particles have moved  
      !Dump back VolFrac
       PhiP(i) = Phip(i) / (pRegion%grid%vol(i) * global%ppiclfInitZpfFactor)!*(0.0006/0.0002))!*2.0*global%pi/0.0001) 
!VOL Frac Cap
       if (Phip(i) .gt. 0.6) phip(i) = 0.6  
       pRegion%mixt%piclVF(i) = PhiP(i) 
end DO



!Deallocate arrays

    DEALLOCATE(YTEMP,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error


    DEALLOCATE(rhoF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uxF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uyF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uzF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(csF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(vfP,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpxF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpyF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpzF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error
        
    DEALLOCATE(SDRX,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDRY,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SDRZ,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(rhsR,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(pGcX,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(pGcY,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF !global%error    

    DEALLOCATE(pGcZ,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF !global%error    

    DEALLOCATE(JFX,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFXCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFY,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFYCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFZ,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(JFZCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(PhiP,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

#endif
!PPICLF Integration END

! ----- update solution; sum up residuals -------------------------------------

        CALL RKUpdateMP( regions(iReg),iReg,istage )

        IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_BXV_ComputeVarsCv(pRegion)
        END IF ! RFLU_NSCBC_DecideHaveNSCBC

! ----- program burn, after RKUpdateMP ----------------------------------------
#ifdef SPEC        
        IF ( (global%specUsed .EQV. .TRUE.) .AND. &
             (global%pbaFlag .EQV. .TRUE.) ) THEN
! DEBUG: Manoj-PBA1D
!          CALL SPEC_RFLU_PBA_ReactionZone( pRegion )
! END DEBUG
        END IF ! pbaFlag
#endif
        
! ----- set ghost fluid solution ----------------------------------------------

! TEMPORARY: Manoj: GFM, need to put more thinking
        CALL RFLU_GFM_SetGhostFluid( pRegion )

! ----- perform checks and enforce after-update conditions --------------------

!        CALL AfterUpdateMP( pRegion,istage )

! ----- Descale grid speeds -----------------------------------------------------

        CALL RFLU_DescaleGridSpeeds(pRegion)

! ----- program burn, before SetDependentVars ---------------------------------
#ifdef SPEC        
        IF ( (global%specUsed .EQV. .TRUE.) .AND. &
             (global%pbaFlag .EQV. .TRUE.) ) THEN
! DEBUG: Manoj-PBA1D
!          CALL SPEC_RFLU_PBA_ProgramBurn( pRegion )
! END DEBUG
        END IF ! pbaFlag
#endif

! ----- update dependent variables --------------------------------------------

        CALL RFLU_MPI_ISendWrapper(pRegion)
        CALL RFLU_SetVarsContWrapper(pRegion,1,pRegion%grid%nCells)

! Subbu - Perform check after computing pressure & temperature
! ----- perform checks and enforce after-update conditions --------------------

        CALL AfterUpdateMP( pRegion,istage )
! Subbu - End Perform check

! ----- update dependent variables on boundary faces --------------------------

        IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
          CALL RFLU_BXV_SetDependentVars(pRegion)
        END IF !
    END DO ! iReg

    CALL RFLU_MPI_CopyWrapper(regions)

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)

      CALL RFLU_MPI_RecvWrapper(pRegion)
      CALL RFLU_SetVarsContWrapper(pRegion,pRegion%grid%nCells+1, & 
                                   pRegion%grid%nCellsTot)
      CALL RFLU_RELP_TransformWrapper(pRegion)                                   
    END DO ! iReg

    DO iReg = 1,global%nRegionsLocal
      pRegion => regions(iReg)

      CALL RFLU_MPI_ClearRequestWrapper(pRegion)
    END DO ! iReg

#ifdef PLAG
    IF ( global%plagUsed .EQV. .TRUE. ) THEN
      CALL PLAG_RFLU_CommDriver(regions)

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_SetVarsDiscWrapper(pRegion)
      END DO ! iReg

! --- update particle volume fraction in virtual cells ------------------------      
! TEMPORARY: Manoj: 2012-05-29: checking effect of not communicating vFracE
!IF (1==2) THEN
      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_MPI_PLAG_ISendWrapper(pRegion)
      END DO ! iReg 

      CALL RFLU_MPI_PLAG_CopyWrapper(regions)
    
      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_MPI_PLAG_RecvWrapper(pRegion)
      END DO ! iReg 

      DO iReg = 1,global%nRegionsLocal
        pRegion => regions(iReg)

        CALL RFLU_MPI_ClearRequestWrapper(pRegion)
      END DO ! iReg 
!END IF ! 1==2
! END TEMPORARY
    END IF ! global%plagUsed
 
#endif
 
!Kept for ref, set to delete$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#ifdef PICL
if (1.eq.2) then
     piclIO = 100000000
     piclDtMin = REAL(global%dtMin,8)
     piclCurrentTime = REAL(global%currentTime,8)
     doWrite = RFLU_DecideWrite(global)
!Figure out piclIO call, might need to look into timestepping
    IF ( (doWrite .EQV. .TRUE.)) piclIO = 1

!PARTICLE stuff possbile needed
!    CALL RFLU_ConvertCvCons2Prim(pRegion,CV_MIXT_STATE_DUVWP)


!allocate arrays to send to picl
    nCells = pRegion%grid%nCells
    ALLOCATE(rhoF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uxF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(uyF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(csF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(vfP(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(dpxF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error
    
    ALLOCATE(dpyF(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error


!Might need to update prim like plag does

!Fill arrays for interp field
    DO i = 1,pRegion%grid%nCells
       do lz=1,2
       do ly=1,2
       do lx=1,2 
       call ppiclf_solve_GetProFldIJKEF(lx,ly,lz,i,PPICLF_P_JPHIP,vfP(lx,ly,lz,i))
       vfP(lx,ly,lz,i) = vfP(lx,ly,lz,i) ! / pRegion%grid%vol(i)    
       rhoF(lx,ly,lz,i) = pRegion%mixt%cv(CV_MIXT_DENS,i)
       uxF(lx,ly,lz,i) = pRegion%mixt%cv(CV_MIXT_XMOM,i) &
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)
       uyF(lx,ly,lz,i) = pRegion%mixt%cv(CV_MIXT_YMOM,i) &
                        /pRegion%mixt%cv(CV_MIXT_DENS,i)
!       uzF(lx,ly,lz,i) = pRegion%mixt%cv(CV_MIXT_ZMOM,i) &
!                        /pRegion%mixt%cv(CV_MIXT_DENS,i)
       csF(lx,ly,lz,i) = pRegion%mixt%dv(DV_MIXT_SOUN,i)
       dpxF(lx,ly,lz,i) = pRegion%mixt%gradCell(XCOORD,GRC_MIXT_PRES,i)
       dpyF(lx,ly,lz,i) = pRegion%mixt%gradCell(YCOORD,GRC_MIXT_PRES,i) 
       end do
       end do
       end do 
    END DO

!Interp field calls

!write(*,*) 'Pdtmin,Pcurr time', piclIO,piclDtMin,piclCurrentTime 

      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JRHOF,rhoF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUX,uxF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUY,uyF)
!      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JUZ,uzF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPDX,dpxF)
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JDPDY,dpyF)  
      CALL ppiclf_solve_InterpFieldUser(PPICLF_R_JCS,csF)
      call ppiclf_solve_InterpFieldUser(PPICLF_R_JPHIP,vfP)  

!SOLVE
     
  WRITE(*,*) "Calling ppiclf_solve_IntegrateParticle loc2"
  CALL ppiclf_solve_IntegrateParticle(1,piclIO,piclDtMin,piclCurrentTime)
  WRITE(*,*) "Finished Calling ppiclf_solve_IntegrateParticle loc2"

  !Deallocate arrays
    DEALLOCATE(rhoF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uxF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(uyF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(csF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(vfP,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpxF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(dpyF,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

end if
#endif
!Kept for ref set to delete END$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$




  
  END DO ! istage

! finalize ====================================================================

  CALL DeregisterFunction( global )

END SUBROUTINE RungeKuttaMP

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RungeKuttaMP.F90,v $
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.8  2009/08/28 18:29:39  mtcampbe
! RocfluMP integration with Rocstar and some makefile tweaks.  To build
! Rocstar with new Rocflu:
! make ROCFLU=RocfluMP
! To build Rocstar with the new RocfluND:
! make ROCFLU=RocfluMP HYPRE=/the/hypre/install/path
!
! Revision 1.7  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:16:48  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2008/05/29 01:35:11  mparmar
! Added Setting of reference frame velocity and update of BC here
!
! Revision 1.4  2007/12/04 13:36:59  haselbac
! Bug fix: Removed USE RFLU_ModPatchVelocity
!
! Revision 1.3  2007/12/03 16:34:03  mparmar
! Removed RFLU_SetPatchVelocity
!
! Revision 1.2  2007/06/18 17:42:13  mparmar
! Added calls for moving reference frame implementation
!
! Revision 1.1  2007/04/09 18:48:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.16  2007/03/27 00:39:50  haselbac
! Removed call to PLAG_CalcnPclsTotGlobal, now in RFLU_TimeStepping
!
! Revision 1.15  2007/03/20 22:02:29  fnajjar
! Included call to PLAG_CalcnPclsTotGlobal
!
! Revision 1.14  2006/08/21 16:10:01  haselbac
! Adapted to name change
!
! Revision 1.13  2006/08/19 15:48:25  mparmar
! Added computations of boundary Cv and Dv for NSCBC implementation
!
! Revision 1.12  2006/08/18 21:09:27  fnajjar
! Removed IF around PLAG_RFLU_CommDriver for serial periodic cases
!
! Revision 1.11  2006/03/25 21:40:03  haselbac
! Added call to transforming data on related patches, cosmetics
!
! Revision 1.10  2005/12/03 19:44:54  haselbac
! Apparent bug fix: Separated call to RFLU_MPI_ClearRequestWrapper into separate loop
!
! Revision 1.9  2005/12/01 21:52:14  fnajjar
! Added IF statement around PLAG_RFLU_CommDriver, only active for more than one nRegions
!
! Revision 1.8  2005/11/10 22:21:07  fnajjar
! ACH: Proper fix for updating PLAG dv
!
! Revision 1.7  2005/11/10 16:51:28  fnajjar
! Added plagUsed IF statement around PLAG routines
!
! Revision 1.6  2005/11/02 14:53:24  haselbac
! Fady: Temporary fix so comm particles get non-cv vars updated properly
!
! Revision 1.5  2005/05/18 22:04:41  fnajjar
! Added PLAG communication routines; only initial implementation
!
! Revision 1.4  2005/04/29 00:06:09  haselbac
! Added routines to clear send requests
!
! Revision 1.3  2005/04/15 15:06:06  haselbac
! Converted to MPI
!
! Revision 1.2  2005/03/31 16:31:02  haselbac
! Added call to RFLU_SetTimeRK
!
! Revision 1.1  2004/12/01 16:51:13  haselbac
! Initial revision after changing case
!
! Revision 1.18  2004/11/14 19:36:23  haselbac
! Replaced call to UpdateDependentVarsMP by RFLU_SetVarsWrapper
!
! Revision 1.17  2004/07/30 22:47:33  jferry
! Implemented Equilibrium Eulerian method for Rocflu
!
! Revision 1.16  2004/04/14 02:07:02  haselbac
! Added grid-speed scaling calls for RFLU
!
! Revision 1.15  2004/03/25 21:14:20  jferry
! changed AfterUpdate to call most subroutines only after final RK stage
!
! Revision 1.14  2004/03/02 21:47:28  jferry
! Added After Update interactions
!
! Revision 1.13  2004/02/26 21:11:58  wasistho
! added globalCommunication
!
! Revision 1.12  2004/02/26 21:01:46  haselbac
! Enclosed updateBoundaryConditionsMP within ifdef RFLO
!
! Revision 1.11  2004/01/29 22:52:47  haselbac
! Added calls to RFLU_EnforceBoundsWrapper and updateDependentVarsMP
!
! Revision 1.10  2003/12/04 03:23:06  haselbac
! Added call to CellGradientsMP and validity check
!
! Revision 1.9  2003/11/25 21:01:45  haselbac
! Added calls to RFLU_UpdateDummyCells and ZeroResidualsMP
!
! Revision 1.8  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.5  2003/10/03 20:42:07  haselbac
! Added Rocflu calls
!
! Revision 1.4  2003/10/01 23:52:09  jblazek
! Corrected bug in moving noslip wall BC and grid speeds.
!
! Revision 1.3  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.2  2003/04/10 01:22:41  jblazek
! Got rid of pRegion in ViscousFluxesMP.
!
! Revision 1.1  2003/03/28 19:42:55  fnajjar
! Initial import for RocfluidMP
!
! ******************************************************************************

