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
!******************************************************************************
!
! Purpose: compute interaction source for drag forces on Lagrangian particles.
!
! Description: none.
!
! Input: region  = current region.
!
! Output: region%levels(iLev)%plag%inrtSources
!
! Notes: none.
!
!******************************************************************************
!
! $Id: INRT_CalcDragUnsteady_AMExplicit.F90,v 1.3 2016/02/08 22:26:29 rahul Exp $
!
! Copyright: (c) 2003 by the University of Illinois
!
!******************************************************************************

SUBROUTINE INRT_CalcDragUnsteady_AMExplicit( region )

  USE ModDataTypes
  USE ModDataStruct, ONLY : t_region
  USE ModGlobal,     ONLY : t_global
  USE ModPartLag,    ONLY : t_plag
  USE ModError
  USE ModParameters
  USE INRT_ModParameters
  USE ModMPI

#ifdef PLAG
  USE PLAG_ModParameters
#endif

  USE RFLU_ModDifferentiationCells, ONLY: RFLU_ComputeGradCellsGGScalar, &
                                          RFLU_ComputeGradCellsGGVector

  IMPLICIT NONE

! ... parameters
  TYPE(t_region), INTENT(INOUT), TARGET :: region

! ... loop variables
  INTEGER :: iPcls

! ... local variables
  CHARACTER(CHRLEN)  :: RCSIdentString

  INTEGER :: dragUnsteady, iCont, nCont, nPcls
  INTEGER, POINTER, DIMENSION(:)   :: pCvPlagMass

  REAL(RFREAL) :: CdTotal,diamL,factor,gamma,machL,massL,mixtVolR,pi,psiL, &
                  relVelMagL,reyL,tauLR
  REAL(RFREAL),          DIMENSION(3)   :: relVel, accelL
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pCv, pDv, pTv
  REAL(RFREAL), POINTER, DIMENSION(:,:,:) :: dudtMixt,dudtPlag

  TYPE(t_plag)  , POINTER :: pPlag
  TYPE(t_global), POINTER :: global
  TYPE(t_region), POINTER :: pRegion

! New variables for unsteady force ----------
  INTEGER :: icg,iT,nUnsteadyData
  REAL(RFREAL), DIMENSION(3) :: forceIU,forcePG,forceTotal,forceVU
  REAL(RFREAL) :: A,B,CamEff,dt,fH,kernelVU,mf,mu,nu,refArea,rhoMixt, &
                  speedSound,time,vFrac,vFracCorr,volL,volMixt
  ! Subbu - Vars cyldet case
  REAL(RFREAL), DIMENSION(:,:,:), POINTER :: pGc
  REAL(RFREAL) :: ConvFlux,tBeg,tEnd
  INTEGER :: Coord
  ! Subbu - End Vars cyldet case 

! New variables for augmenting rhs ----------
  REAL(RFREAL) :: coeffIU,contFac,energydotg,energydotp
  REAL(RFREAL) :: drudtMixt,drvdtMixt,drwdtMixt,drudtPlag,drvdtPlag,drwdtPlag 
  REAL(RFREAL), DIMENSION(3) :: ug,up
! -------------------------------------------
!******************************************************************************

  RCSIdentString = '$RCSfile: INRT_CalcDragUnsteady_AMExplicit.F90,v $'

  global => region%global
  pRegion => region

  CALL RegisterFunction( global,'INRT_CalcDragUnsteady_AMExplicit',__FILE__ )

#ifdef PLAG
! Check if there are any particles

  nPcls = 0
  IF (global%plagUsed) nPcls = region%plag%nPcls

  IF (nPcls < 1) GO TO 999

! Get dimensions --------------------------------------------------------------

  pi = global%pi
  dt = global%dtMin

  nCont        = region%plagInput%nCont
  dragUnsteady = region%inrtInput%inrts(INRT_TYPE_DRAG)%switches( &
    INRT_SWI_DRAG_UNSTEADY)

  SELECT CASE (dragUnsteady)

    CASE (INRT_DRAG_UNSTEADY_NONE)

    CASE (INRT_DRAG_UNSTEADY_USE)

    pPlag     => region%plag

    pCv       => pPlag%cv
    pDv       => pPlag%dv
    pTv       => pPlag%tv

    dudtMixt  => pPlag%dudtMixt
    dudtPlag  => pPlag%dudtPlag

    pCvPlagMass => pPlag%cvPlagMass

    nUnsteadyData = region%plagInput%nUnsteadyData

! =============================================================================
!   Step 1: Pressure gradient and viscous unsteady force
! =============================================================================

    DO iPcls = 1,nPcls
      icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcls)

      diamL   = pDv(DV_PLAG_DIAM,iPcls)
      volL    = (pi*(diamL**3.0_RFREAL)/6.0_RFREAL)
      volMixt = region%grid%vol(icg) 

      vFrac = pPlag%vFracL(1,iPcls)

      massL = SUM( pCv(pCvPlagMass(:),iPcls) )

      relVel(1) = pDv(DV_PLAG_UVELMIXT,iPcls)-pDv(DV_PLAG_UVEL,iPcls)
      relVel(2) = pDv(DV_PLAG_VVELMIXT,iPcls)-pDv(DV_PLAG_VVEL,iPcls)
      relVel(3) = pDv(DV_PLAG_WVELMIXT,iPcls)-pDv(DV_PLAG_WVEL,iPcls)

      relVelMagL = SQRT( relVel(1)*relVel(1)+ &
                         relVel(2)*relVel(2)+ &
                         relVel(3)*relVel(3)  )

! TEMPORARY: Manoj, Need proper method to interpolate mixture gamma at particle location. 
      gamma = 1.4_RFREAL   ! For shock-particle-cloud, rectangle to triangle.
! END TEMPORARY

      mu         = pTv(TV_PLAG_MUELMIXT,iPcls)
      rhoMixt    = pDv(DV_PLAG_DENSMIXT,iPcls)/(1.0_RFREAL-vFrac)
      nu         = mu/rhoMixt
      mf         = rhoMixt*volL
      speedSound = SQRT(gamma*pDv(DV_PLAG_PRESMIXT,iPcls)/rhoMixt)

      reyL  = diamL*relVelMagL*rhoMixt/mu
      machL = relVelMagL/speedSound

! -----------------------------------------------------------------------------
!     Code for pressure gradient and viscous unsteady force
! -----------------------------------------------------------------------------

      forcePG    = 0.0_RFREAL
      forceVU    = 0.0_RFREAL
      forceTotal = 0.0_RFREAL

! --- Pressure gradient force
! --- Rahul, modified for AUSM+up
!      IF ((region%mixtInput%spaceDiscr .EQ. DISCR_UPW_AUSMPLUSUP)) THEN
!        forcePG(XCOORD) = 0.0_RFREAL
!        forcePG(YCOORD) = 0.0_RFREAL
!        forcePG(ZCOORD) = 0.0_RFREAL
!      ELSE
        forcePG(XCOORD) = -volL*pPlag%pgMixt(XCOORD,iPcls)
        forcePG(YCOORD) = -volL*pPlag%pgMixt(YCOORD,iPcls)
        forcePG(ZCOORD) = -volL*pPlag%pgMixt(ZCOORD,iPcls)
!     END IF ! spaceDiscr
! rahul - end

! Subbu - Copy forces to new vars to wite to output file
        pPlag%inrtPG(INRT_PG_L_XMOM_G,iPcls) = forcePG(XCOORD)
        pPlag%inrtPG(INRT_PG_L_YMOM_G,iPcls) = forcePG(YCOORD)
        pPlag%inrtPG(INRT_PG_L_ZMOM_G,iPcls) = forcePG(ZCOORD)
! Subbu - End copy forces to new vars to wite to output file

! --- Viscous unsteady force
      ! Subbu - Turning off viscous unsteady force
      !F (1 == 2) THEN
      IF ( region%inrtInput%inrts(INRT_TYPE_DRAG)%switches(INRT_SWI_DRAG_VU) &
           == INRT_DRAG_VISCUNST_USE) THEN

      fH     = (0.75_RFREAL + .105_RFREAL*reyL)
      factor = 3.0_RFREAL*pi*mu*diamL*dt

      IF ( region%plagInput%nTimeBH > 1 ) THEN
        DO iT=2,region%plagInput%nTimeBH-1
          time = pPlag%timeBH(iT)

          A  = (4.0_RFREAL*pi*time*nu/diamL**2.0_RFREAL)**(.25_RFREAL)
          B  = (0.5_RFREAL*pi*(relVelMagL**3.0_RFREAL)*(time**2.0_RFREAL)/ &
               (0.5_RFREAL*diamL*nu*fH**3.0_RFREAL))**(.5_RFREAL)

          kernelVU = factor*(A+B)**(-2.0_RFREAL)

          forceVU(XCOORD) = forceVU(XCOORD) &
                      + kernelVU*(dudtMixt(XCOORD,iT,iPcls)-dudtPlag(XCOORD,iT,iPcls)) 
          forceVU(YCOORD) = forceVU(YCOORD) &
                      + kernelVU*(dudtMixt(YCOORD,iT,iPcls)-dudtPlag(YCOORD,iT,iPcls)) 
          forceVU(ZCOORD) = forceVU(ZCOORD) &
                      + kernelVU*(dudtMixt(ZCOORD,iT,iPcls)-dudtPlag(ZCOORD,iT,iPcls)) 
        END DO ! iT

        iT   = region%plagInput%nTimeBH
        time = pPlag%timeBH(iT)

        A  = (4.0_RFREAL*pi*time*nu/diamL**2.0_RFREAL)**(.25_RFREAL)
        B  = (0.5_RFREAL*pi*(relVelMagL**3.0_RFREAL)*(time**2.0_RFREAL)/ &
             (0.5_RFREAL*diamL*nu*fH**3.0_RFREAL))**(.5_RFREAL)

        kernelVU = 0.5_RFREAL*factor*(A+B)**(-2.0_RFREAL)

        forceVU(XCOORD) = forceVU(XCOORD) &
                      + kernelVU*(dudtMixt(XCOORD,iT,iPcls)-dudtPlag(XCOORD,iT,iPcls)) 
        forceVU(YCOORD) = forceVU(YCOORD) &
                      + kernelVU*(dudtMixt(YCOORD,iT,iPcls)-dudtPlag(YCOORD,iT,iPcls)) 
        forceVU(ZCOORD) = forceVU(ZCOORD) &
                      + kernelVU*(dudtMixt(ZCOORD,iT,iPcls)-dudtPlag(ZCOORD,iT,iPcls)) 
         ! Subbu - Copy forces to new vars to wite to output file
         pPlag%inrtVU(INRT_VU_L_XMOM_G,iPcls) = forceVU(XCOORD)
         pPlag%inrtVU(INRT_VU_L_YMOM_G,iPcls) = forceVU(YCOORD)
         pPlag%inrtVU(INRT_VU_L_ZMOM_G,iPcls) = forceVU(ZCOORD)
         ! Subbu - End copy forces to new vars to wite to output file

      END IF ! region%plagInput%nTimeBH

      END IF  
      ! Subbu - End turning off viscous unsteady force

! -----------------------------------------------------------------------------
!     The forceTotal represents force on a single particle, NOT superparticle
! -----------------------------------------------------------------------------

      forceTotal = forcePG + forceVU

      contFac = pPlag%arv(ARV_PLAG_SPLOAD,iPcls)

      IF ( region%dummyStep .EQV. .FALSE. ) THEN
        ug = region%plag%dv(DV_PLAG_UVELMIXT:DV_PLAG_WVELMIXT,iPcls)
        up = region%plag%dv(DV_PLAG_UVEL:DV_PLAG_WVEL,iPcls)
        energydotg = DOT_PRODUCT(forceTotal,ug)
        energydotp = DOT_PRODUCT(forceTotal,up)

! Rahul - Consistent energy and momentum formulation for AUSM+up 
! NOTE: forceIU being inviscid should not contribute to fluid phase internal energy.
!     
        IF ((region%mixtInput%spaceDiscr .EQ. DISCR_UPW_AUSMPLUSUP)) THEN

! ----- Augment Gas Sources - momentum. 
! AUSM+up correction remove contribution of pressure gradient
         region%mixt%rhs(CV_MIXT_XMOM:CV_MIXT_ZMOM,icg) &
                         = region%mixt%rhs(CV_MIXT_XMOM:CV_MIXT_ZMOM,icg) &
                         + contFac*(forceTotal - forcePG)

! ----- Augment Gas Sources - energy. 
        region%mixt%rhs(CV_MIXT_ENER,icg) &
                        = region%mixt%rhs(CV_MIXT_ENER,icg) &
                        + contFac*(energydotg - DOT_PRODUCT(forcePG,ug))
        ELSE
         region%mixt%rhs(CV_MIXT_XMOM:CV_MIXT_ZMOM,icg) &
                         = region%mixt%rhs(CV_MIXT_XMOM:CV_MIXT_ZMOM,icg) &
                         + contFac*forceTotal
         region%mixt%rhs(CV_MIXT_ENER,icg) &
                         = region%mixt%rhs(CV_MIXT_ENER,icg) &
                         + contFac*energydotg
        END IF
! Rahul - end

! ----- Augment Gas Sources
!        region%mixt%rhs(CV_MIXT_XMOM:CV_MIXT_ZMOM,icg) &
!                        = region%mixt%rhs(CV_MIXT_XMOM:CV_MIXT_ZMOM,icg) &
!                        + contFac*forceTotal
!        region%mixt%rhs(CV_MIXT_ENER,icg) &
!                        = region%mixt%rhs(CV_MIXT_ENER,icg) &
!                        + contFac*energydotg

! ----- Augment Particle Sources
        region%plag%rhs(CV_PLAG_XMOM:CV_PLAG_ZMOM,iPcls) &
                        = region%plag%rhs(CV_PLAG_XMOM:CV_PLAG_ZMOM,iPcls) &
                        - forceTotal
        region%plag%rhs(CV_PLAG_ENER,iPcls) &
                        = region%plag%rhs(CV_PLAG_ENER,iPcls) &
                        - energydotp
      END IF ! region%dummyStep 

! -----------------------------------------------------------------------------
!     Update force total
! -----------------------------------------------------------------------------

      ! Subbu - Copy forces to new vars to wite to output file
      pPlag%inrtTOT(INRT_TOT_L_XMOM_G,iPcls) = forceTotal(XCOORD) 
      pPlag%inrtTOT(INRT_TOT_L_YMOM_G,iPcls) = forceTotal(YCOORD)
      pPlag%inrtTOT(INRT_TOT_L_ZMOM_G,iPcls) = forceTotal(ZCOORD)        
      ! Subbu - End copy forces to new vars to wite to output file

      pPlag%forceTotal(XCOORD,icg) = pPlag%forceTotal(XCOORD,icg) + forceTotal(1)
      pPlag%forceTotal(YCOORD,icg) = pPlag%forceTotal(YCOORD,icg) + forceTotal(2)
      pPlag%forceTotal(ZCOORD,icg) = pPlag%forceTotal(ZCOORD,icg) + forceTotal(3)
    END DO ! iPcls



! =============================================================================
!  Inviscid-unsteady force - Explicit scheme
! =============================================================================


! =============================================================================
!   Step 2: Compute Eulerian field gradients
! =============================================================================

! Subbu - Turn off computing cell gradients 
!IF  (1==2) THEN
!    CALL RFLU_ComputeGradCellsGGScalar(pRegion,1,1,1,1,pRegion%mixt%cv, &
!                                       pRegion%mixt%gradCellE)

!    CALL RFLU_ComputeGradCellsGGVector(pRegion,2,4,2,4,pRegion%mixt%cv, &
!                                       pRegion%mixt%gradCellE)
!END IF
! Subbu - End turn off computing cell gradients 

! =============================================================================
!   Step 3: Fluid momentum equation source terms due to Inviscid unsteady force
!           Accumulate coeffIU
! =============================================================================

    pGc  => pRegion%mixt%gradCell
    DO iPcls = 1,nPcls
      icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcls)

      diamL   = pDv(DV_PLAG_DIAM,iPcls)
      volL    = (pi*(diamL**3.0_RFREAL)/6.0_RFREAL)
      volMixt = region%grid%vol(icg)

      vFrac = pPlag%vFracL(1,iPcls)

      massL = SUM( pCv(pCvPlagMass(:),iPcls) )

      relVel(1) = pDv(DV_PLAG_UVELMIXT,iPcls)-pDv(DV_PLAG_UVEL,iPcls)
      relVel(2) = pDv(DV_PLAG_VVELMIXT,iPcls)-pDv(DV_PLAG_VVEL,iPcls)
      relVel(3) = pDv(DV_PLAG_WVELMIXT,iPcls)-pDv(DV_PLAG_WVEL,iPcls)

      relVelMagL = SQRT( relVel(1)*relVel(1)+ &
                         relVel(2)*relVel(2)+ &
                         relVel(3)*relVel(3)  )

! TEMPORARY: Manoj, Need proper method to interpolate mixture gamma at particle
! location. 
      gamma = 1.4_RFREAL   ! For shock-particle-cloud, rectangle to triangle.
! END TEMPORARY

      rhoMixt    = pDv(DV_PLAG_DENSMIXT,iPcls)/(1.0_RFREAL-vFrac)
      mf         = rhoMixt*volL
      speedSound = SQRT(gamma*pDv(DV_PLAG_PRESMIXT,iPcls)/rhoMixt)

      machL = relVelMagL/speedSound

      if ( machL < 0.6_RFREAL ) THEN
        CamEff = 0.5_RFREAL*(1.0_RFREAL + 1.8_RFREAL*machL*machL &
                                        + 7.6_RFREAL*machL**4.0_RFREAL)
      ELSE
        CamEff = 0.5_RFREAL*(1.0_RFREAL + 1.8_RFREAL*0.6_RFREAL*0.6_RFREAL &
                                        + 7.6_RFREAL*0.6_RFREAL**4.0_RFREAL)
      END IF ! machL

! -----------------------------------------------------------------------------
!     Code for inviscid unsteady force
! -----------------------------------------------------------------------------
      ! Subbu - Use already existing convective flux without computing gradient
      forceIU    = 0.0_RFREAL
      forceTotal = 0.0_RFREAL

      vFracCorr = (1.0_RFREAL + 2.0_RFREAL*vFrac)/(1.0_RFREAL-vFrac)
      coeffIU   = vFracCorr*CamEff*volL

      ug = region%plag%dv(DV_PLAG_UVELMIXT:DV_PLAG_WVELMIXT,iPcls)
      up = region%plag%dv(DV_PLAG_UVEL:DV_PLAG_WVEL,iPcls)
      
      drudtMixt = -region%mixt%rhs(CV_MIXT_XMOM,icg)/volMixt &
                  +rhoMixt*DOT_PRODUCT(ug,pGc(:,2,icg)) &
                  +ug(XCOORD)*DOT_PRODUCT(ug,pGc(:,1,icg))

      drvdtMixt = -region%mixt%rhs(CV_MIXT_YMOM,icg)/volMixt &
                  +rhoMixt*DOT_PRODUCT(ug,pGc(:,3,icg)) &
                  +ug(YCOORD)*DOT_PRODUCT(ug,pGc(:,1,icg))

      drwdtMixt = -region%mixt%rhs(CV_MIXT_ZMOM,icg)/volMixt &
                  +rhoMixt*DOT_PRODUCT(ug,pGc(:,4,icg)) &
                  +ug(ZCOORD)*DOT_PRODUCT(ug,pGc(:,1,icg))

      drudtPlag = rhoMixt*dudtPlag(XCOORD,1,iPcls) &
                - up(XCOORD)*region%mixt%rhs(CV_MIXT_DENS,icg)/volMixt &
                + up(XCOORD)*DOT_PRODUCT(up,pGc(:,1,icg))

      drvdtPlag = rhoMixt*dudtPlag(YCOORD,1,iPcls) &
                - up(YCOORD)*region%mixt%rhs(CV_MIXT_DENS,icg)/volMixt &
                + up(YCOORD)*DOT_PRODUCT(up,pGc(:,1,icg))

      drwdtPlag = rhoMixt*dudtPlag(ZCOORD,1,iPcls) &
                - up(ZCOORD)*region%mixt%rhs(CV_MIXT_DENS,icg)/volMixt &
                + up(ZCOORD)*DOT_PRODUCT(up,pGc(:,1,icg))

      forceIU(XCOORD) = coeffIU*( drudtMixt - drudtPlag )
      forceIU(YCOORD) = coeffIU*( drvdtMixt - drvdtPlag )
      forceIU(ZCOORD) = coeffIU*( drwdtMixt - drwdtPlag )

      ! Subbu - Copy forces to new vars to wite to output file
      pPlag%inrtIU(INRT_IU_L_XMOM_G,iPcls) = forceIU(XCOORD)
      pPlag%inrtIU(INRT_IU_L_YMOM_G,iPcls) = forceIU(YCOORD)
      pPlag%inrtIU(INRT_IU_L_ZMOM_G,iPcls) = forceIU(ZCOORD)
      ! Subbu - End copy forces to new vars to wite to output file

      ! Subbu - End use of already existing convective flux without computing gradient
! -----------------------------------------------------------------------------
!     The forceIU represents force on a single particle, NOT superparticle
! -----------------------------------------------------------------------------

      forceTotal = forceIU

      contFac = pPlag%arv(ARV_PLAG_SPLOAD,iPcls)

      IF ( region%dummyStep .EQV. .FALSE. ) THEN
        ug = region%plag%dv(DV_PLAG_UVELMIXT:DV_PLAG_WVELMIXT,iPcls)
        up = region%plag%dv(DV_PLAG_UVEL:DV_PLAG_WVEL,iPcls)
        energydotg = DOT_PRODUCT(forceTotal,ug)
        energydotp = DOT_PRODUCT(forceTotal,up)


! ----- Augment Gas Sources - momentum. 
        region%mixt%rhs(CV_MIXT_XMOM:CV_MIXT_ZMOM,icg) &
                        = region%mixt%rhs(CV_MIXT_XMOM:CV_MIXT_ZMOM,icg) &
                        + contFac*forceTotal 

! ----- Augment Particle Sources - momentum
        region%plag%rhs(CV_PLAG_XMOM:CV_PLAG_ZMOM,iPcls) &
                        = (region%plag%rhs(CV_PLAG_XMOM:CV_PLAG_ZMOM,iPcls) &
                           -forceTotal)
! ----- Augment Gas Sources - energy
        region%mixt%rhs(CV_MIXT_ENER,icg) &
                        = region%mixt%rhs(CV_MIXT_ENER,icg) &
                        + contFac*energydotg !&

! ----- Augment Particle Sources - energy
        region%plag%rhs(CV_PLAG_ENER,iPcls) &
                        = region%plag%rhs(CV_PLAG_ENER,iPcls) &
                        - energydotp 
      END IF ! region%dummyStep 
! -----------------------------------------------------------------------------
!     Update force total
! -----------------------------------------------------------------------------

      pPlag%forceTotal(XCOORD,icg) = pPlag%forceTotal(XCOORD,icg) + forceTotal(1)
      pPlag%forceTotal(YCOORD,icg) = pPlag%forceTotal(YCOORD,icg) + forceTotal(2)
      pPlag%forceTotal(ZCOORD,icg) = pPlag%forceTotal(ZCOORD,icg) + forceTotal(3)

      ! Subbu - Copy forces to new vars to wite to output file
      pPlag%inrtTOT(INRT_TOT_L_XMOM_G,iPcls) = pPlag%inrtTOT(INRT_TOT_L_XMOM_G,iPcls) &
                                               + pPlag%inrtSources(INRT_DRAG_L_XMOM_G,iPcls) &
                                               + forceTotal(XCOORD) 
      pPlag%inrtTOT(INRT_TOT_L_YMOM_G,iPcls) = pPlag%inrtTOT(INRT_TOT_L_YMOM_G,iPcls) &
                                               + pPlag%inrtSources(INRT_DRAG_L_YMOM_G,iPcls) &
                                               + forceTotal(YCOORD) 
      pPlag%inrtTOT(INRT_TOT_L_ZMOM_G,iPcls) = pPlag%inrtTOT(INRT_TOT_L_ZMOM_G,iPcls) &
                                               + pPlag%inrtSources(INRT_DRAG_L_ZMOM_G,iPcls) &
                                               + forceTotal(ZCOORD) 
      ! Subbu - End copy forces to new vars to wite to output file

! -----------------------------------------------------------------------------
!     Update particle acceleration
!     inrtSources now consist of quasi-steady and unsteady forces
! -----------------------------------------------------------------------------

      IF ( region%dummyStep .EQV. .FALSE. ) THEN
        dudtPlag(XCOORD,1,iPcls) = -region%plag%rhs(CV_PLAG_XMOM,iPcls)/massL
        dudtPlag(YCOORD,1,iPcls) = -region%plag%rhs(CV_PLAG_YMOM,iPcls)/massL
        dudtPlag(ZCOORD,1,iPcls) = -region%plag%rhs(CV_PLAG_ZMOM,iPcls)/massL
      END IF ! region%dummyStep


    END DO ! iPcls

    
    CASE DEFAULT
      CALL ErrorStop( global,ERR_REACHED_DEFAULT,__LINE__ )

  END SELECT ! dragUnsteady
! finalize --------------------------------------------------------------------

  ! Subbu - Print forces to output file forcesPart.dat
  ! IF (1==2) THEN
  IF (region%irkStep == global%nrkSteps) THEN
   IF(global%currentTime .EQ. 0.00000E+00)THEN
     OPEN(1113,file='forcesPart.dat',form='formatted',status='unknown')
   ELSE
     OPEN(1113,file='forcesPart.dat',form='formatted',status='old' &
         ,position='append')
   END IF !currentTime
   DO iPcls = 1,nPcls
    icg = pPlag%aiv(AIV_PLAG_ICELLS,iPcls)
    WRITE(1113,'(E23.16,1X,I8,1X,4(E23.16,1X))') &
                                global%currentTime, iPcls,&
                                pPlag%inrtSources(1,iPcls),& 
                                !pPlag%inrtPG(1,iPcls),&
                                !pPlag%inrtIU(1,iPcls),&
                                pPlag%inrtTOT(1,iPcls),&
                                pPlag%inrtTOT(2,iPcls),&
                                pPlag%inrtTOT(3,iPcls)
!                region%plag%dv(DV_PLAG_UVELMIXT,iPcls),&
!                region%plag%dv(DV_PLAG_VVELMIXT,iPcls),&
!                     pRegion%mixt%cv(CV_MIXT_XVEL,icg),&
!                     pRegion%mixt%cv(CV_MIXT_YVEL,icg)
   END DO
   CLOSE(1113)
  END IF !irkStep
  ! END IF
  ! Subbu - End print forces to output file forcesPart.dat

999  CONTINUE
#endif
  CALL DeregisterFunction( global )

END SUBROUTINE INRT_CalcDragUnsteady_AMExplicit

!******************************************************************************
!
! RCS Revision history:
!
! $Log: INRT_CalcDragUnsteady_AMExplicit.F90,v $
! Revision 1.3  2016/02/08 22:26:29  rahul
! Subtracted work done contribution from F_pg to fluid phase energy equation.
!
! Revision 1.2  2015/12/18 23:34:41  rahul
! Suppressed the computation of pressure gradient force. This is a
! consequence of governing equations' formulation of multiphase AUSM+up
! scheme.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!
!******************************************************************************

