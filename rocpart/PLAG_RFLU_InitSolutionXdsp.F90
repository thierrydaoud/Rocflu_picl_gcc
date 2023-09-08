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
! Garno
!
! Purpose: Initialize particle solution for Barrel Ejection Explosive Dispersal (xdsp) cases
!
! Description: None.
!
! Input: 
!   pRegion     Pointer to region
!   nPclsSumReg Sum of particles in regions 1 to n-1, n being the current region
!
! Output: None.
!
! ******************************************************************************
!
! $Id: PLAG_RFLU_InitSolutionXdsp.F90,(previously: v 1.1 2015/08/12 03:55:40 brollin Exp) $
!
! Copyright: (c) 2004-2005 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE PLAG_RFLU_InitSolutionXdsp(pRegion,nPclsSumReg)

  USE ModDataTypes
  USE ModError
  USE ModDataStruct, ONLY: t_region
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModGrid, ONLY: t_grid
  USE ModMixture, ONLY: t_mixt_input
  USE ModMPI
  USE ModPartLag, ONLY: t_plag
  USE ModParameters  
  USE ModRandom, ONLY: Rand1Uniform
  
  USE PLAG_ModParameters    

  USE RFLU_ModFaceList, ONLY: RFLU_CreateCell2FaceList, &
                              RFLU_BuildCell2FaceList, &
                              RFLU_DestroyCell2FaceList
  USE RFLU_ModInCellTest, ONLY: RFLU_ICT_TestInCell
     
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_region), POINTER :: pRegion
  INTEGER :: nPclsSumReg


! ==============================================================================
! Locals
! ==============================================================================

  CHARACTER(CHRLEN) :: RCSIdentString
  INTEGER :: errorFlag,icg,icl,iCont,iPcl,loopCounter,m,n,nPclsBeg,nPclsEnd, &
             iPclJ,iPclK,partNum,iYpos ! Josh
  INTEGER, POINTER, DIMENSION(:,:) :: pAiv
  INTEGER, POINTER, DIMENSION(:) :: pCvPlagMass
  LOGICAL :: foundFlag
  REAL(RFREAL) :: dia,gaussAmp,heatCapSum,massRatio,massFluxRatioSum, &
                  massFluxRatioSumR,massFluxRatioLimit,massSum,massSumR, &
                  meanVfrac,perturb,rad,rand,rMinCell,rMaxCell,radtmp,the,spLoad, &
                  tMinCell,thetmp,tMaxCell,tol,T,u,v,vFrac, &
                  volPcl,volPclsSum,w,xLoc,yLoc,z,zLoc,zMinCell,&
                  zMaxCell 
  REAL(RFREAL), POINTER, DIMENSION(:,:) :: pArv,pCv
  REAL(RFREAL), POINTER, DIMENSION(:) :: pDens,pIniComp,pSpcHeat
  TYPE(t_global), POINTER :: global
  TYPE(t_grid),   POINTER :: pGrid
  TYPE(t_plag),   POINTER :: pPlag
  REAL(RFREAL) :: x,y,xMinCell,xMaxCell,yMinCell ,yMaxCell
! Josh
  REAL(RFREAL), DIMENSION(5) :: meanDia0 !  Josh,
  REAL(RFREAL), DIMENSION(5) :: partDens !  Josh,
  REAL(RFREAL), DIMENSION(8) :: Aypos !  Josh,
  REAL(RFREAL) :: meanDia
 
! Josh - Read how to place particles ( UQ or Rocflu )
  TYPE(t_mixt_input), POINTER :: pMixtInput
!=============================================================================  
! ******************************************************************************
! Start
! ******************************************************************************

  RCSIdentString = & 
    '$RCSfile: PLAG_RFLU_InitSolutionXdsp.F90,v $'

  global => pRegion%global

  CALL RegisterFunction(global,'PLAG_RFLU_InitSolutionXdsp', &
                        __FILE__)
 
  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel > VERBOSE_LOW ) THEN
       WRITE(STDOUT,'(A,3X,A,1X,I5.5)') SOLVER_NAME,'Global region:', & 
                                         pRegion%iRegionGlobal          
       WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
                          'Initializing particle solution...' 
                           
  END IF ! global%verbLevel
  
! ******************************************************************************
! Set pointers and values
! ******************************************************************************

  pAiv  => pRegion%plag%aiv
  pArv  => pRegion%plag%arv
  pCv   => pRegion%plag%cv
  pCvPlagMass => pRegion%plag%cvPlagMass
  pDens => pRegion%plagInput%dens
  pGrid => pRegion%grid
  pIniComp => pRegion%plagInput%iniComp
  pPlag => pRegion%plag
  pSpcHeat => pRegion%plagInput%spht

 ! Josh

  pMixtInput => pRegion%mixtInput

! ******************************************************************************
! Build cell-to-face list (needed for in-cell test and getting number of cells)
! ******************************************************************************

  CALL RFLU_CreateCell2FaceList(pRegion)
  CALL RFLU_BuildCell2FaceList(pRegion)

! ==============================================================================  
! Set inverse of massFluxRatioSum to avoid division by zero 
! ==============================================================================  
 write(*,*) "NULLSPL:RCP Is this even running?"!BRAD DEBUG ASVF issue
  massFluxRatioLimit = 1.0E-10_RFREAL
  massFluxRatioSum = SUM(pIniComp)
  IF ( massFluxRatioSum > massFluxRatioLimit ) THEN
    massFluxRatioSumR = 1.0_RFREAL/massFluxRatioSum
  ELSE
    massFluxRatioSumR = 1.0_RFREAL
  END IF ! massFluxRatioSum 

! ******************************************************************************  
! Read initial pcl velocities,temperature,diameter and vFrac from user input
! ******************************************************************************

  u = pRegion%plagInput%iniRandUMin
  v = pRegion%plagInput%iniRandVMin
  w = pRegion%plagInput%iniRandWMin
  T = pRegion%plagInput%iniRandTempMin

   IF ( pMixtInput%prepRealVal27 .eq. 1.00 ) THEN ! UQ particles


        WRITE(*,*) 'Particle Init Type: UQ Particles'
 
!  meanDia = pRegion%plagInput%iniRandDiamMin
        ! Josh - updated August 7, 2018
  partDens (1:5) = (/ 15590.000_RFREAL, 15695.000_RFREAL, 15800.000_RFREAL, 15905.000_RFREAL, 16010.000_RFREAL   /)
  meanDia0 (1:5) =  (/ 1.9920E-03_RFREAL, 2.0040E-03_RFREAL, 2.0160E-03_RFREAL, 2.0280E-03_RFREAL, 2.040E-03_RFREAL /)
 
  meanVfrac = pRegion%plagInput%iniRandSpLoadMin

  ! Temporarily set pPlag%nPcls = 0 
  ! It will be recomputed as we place particles one by one in this routine

  pPlag%nPcls = 0
  tol = 1.0E-14_RFREAL

! ******************************************************************************
! Sum up the number of particles in the all other regions except the current one
! to assign proper global initial PCL_ID.
! ******************************************************************************
        
        ! Josh - looping through cells with particles
        WRITE(*,*) 'rocpart/PLAG_RFLU_InitSolutionXdsp.F90, look through cells'

  DO icl = 1, pGrid%nCells
    IF (INT(pGrid%nPclsPerCell(icl)) .GT. 0) THEN

        !Josh
        WRITE(*,*) 'icl with nPclsPerCell > 0',icl,INT(pGrid%nPclsPerCell(icl))

      x = pGrid%cofg(XCOORD,icl)
      y = pGrid%cofg(YCOORD,icl)
      z = pGrid%cofg(ZCOORD,icl)

! ******************************************************************************
! Compute the min and max r,the,z coordiantes for each cell
! ******************************************************************************
      ! pGrid%dx,pGrid%dy,pGrid%dz computed in 
      ! libflu/RFLU_ComputeGridSpacingXdsp.F90       

      xMinCell = x - 0.5_RFREAL*pGrid%dx
      xMaxCell = x + 0.5_RFREAL*pGrid%dx
      yMinCell = y - 0.5_RFREAL*pGrid%dy
      yMaxCell = y + 0.5_RFREAL*pGrid%dy
      zMinCell = z - 0.5_RFREAL*pGrid%dz
      zMaxCell = z + 0.5_RFREAL*pGrid%dz

      volPclsSum = 0.0_RFREAL
      nPclsBeg = pPlag%nPcls + 1
      nPclsEnd = nPclsBeg + INT(pGrid%nPclsPerCell(icl)) - 1

        ! Josh
              WRITE(*,*) 'pGrid%nPclsPerCell(icl)',pGrid%nPclsPerCell(icl)

! ******************************************************************************
! Use random number generator to generate a random location for each particle
! within a cell
! ******************************************************************************

      loopCounter = 0
      iPcl = nPclsBeg
!      DO WHILE (iPcl .LE. nPclsEnd)
      DO iPcl = nPclsBeg,nPclsBeg
        xLoc = xMinCell &
             + Rand1Uniform(pRegion%randData) * (xMaxCell - xMinCell)
        yLoc = yMinCell &
             + Rand1Uniform(pRegion%randData) * (yMaxCell - yMinCell)
        IF (pRegion%mixtInput%dimens == 2) THEN
          zLoc = z
        ELSE
          zLoc = zMinCell &
               + Rand1Uniform(pRegion%randData) * (zMaxCell - zMinCell)
        END IF

       ! pCv(CV_PLAG_XPOS,iPcl) = xLoc
       ! pCv(CV_PLAG_YPOS,iPcl) = yLoc
       ! pCv(CV_PLAG_ZPOS,iPcl) = zLoc

        Aypos (1:8) =  (/ 1.000E-08_RFREAL, 0.125E-03_RFREAL, 0.254E-03_RFREAL, 1.619E-03_RFREAL, 3.0463E-03_RFREAL, &
                         3.622E-03_RFREAL, 4.4972E-03_RFREAL, 4.86E-03_RFREAL /)

        !pCv(CV_PLAG_XPOS,nPclsBeg:nPclsEnd) = xLoc
        pCv(CV_PLAG_XPOS,nPclsBeg:nPclsEnd) = xMaxCell - 1.00E-06_RFREAL
        !pCv(CV_PLAG_YPOS,nPclsBeg:nPclsEnd) = yLoc
    
        iYpos = minloc(abs(y-Aypos(:)),1)
        pCv(CV_PLAG_YPOS,nPclsBeg:nPclsEnd) = Aypos(iYpos)
        
      
    !        IF ( nPclsEnd .le. 30 ) THEN
    !   ! pCv(CV_PLAG_YPOS,nPclsBeg:nPclsEnd) = yMinCell + 1.000E-06_RFREAL
    !    pCv(CV_PLAG_YPOS,nPclsBeg:nPclsEnd) = 1.000E-08_RFREAL
    !        ELSEIF ( nPclsEnd .le. 60 ) THEN
    !    !pCv(CV_PLAG_YPOS,nPclsBeg:nPclsEnd) = yMinCell + 10.000E-06_RFREAL
    !    pCv(CV_PLAG_YPOS,nPclsBeg:nPclsEnd) = 0.125E-03_RFREAL
    !        ELSEIF ( nPclsEnd .le. 90 ) THEN
    !    pCv(CV_PLAG_YPOS,nPclsBeg:nPclsEnd) = 0.254E-03_RFREAL

    !        ELSEIF ( nPclsEnd .le. 120 ) THEN ! For Tracer Particles - Contact Tracking
    !    pCv(CV_PLAG_YPOS,nPclsBeg:nPclsEnd) = 1.619E-03_RFREAL
    !        ELSEIF ( nPclsEnd .le. 150 ) THEN ! For Tracer Particles - Contact Tracking
    !    pCv(CV_PLAG_YPOS,nPclsBeg:nPclsEnd) = 3.0463E-03_RFREAL
    !        ELSEIF ( nPclsEnd .le. 180 ) THEN ! For Tracer Particles - Contact Tracking
    !    pCv(CV_PLAG_YPOS,nPclsBeg:nPclsEnd) = 3.622E-03_RFREAL
    !        ELSEIF ( nPclsEnd .le. 210 ) THEN ! For Tracer Particles - Contact Tracking
    !    pCv(CV_PLAG_YPOS,nPclsBeg:nPclsEnd) = 4.4972E-03_RFREAL
    !        ELSEIF ( nPclsEnd .le. 240 ) THEN ! For Tracer Particles - Contact Tracking
    !    pCv(CV_PLAG_YPOS,nPclsBeg:nPclsEnd) = 4.86E-03_RFREAL
    !        ENDIF

        pCv(CV_PLAG_ZPOS,nPclsBeg:nPclsEnd) = zLoc

        xLoc = pCv(CV_PLAG_XPOS,nPclsBeg)
        yLoc = pCv(CV_PLAG_YPOS,nPclsBeg)

! ******************************************************************************
! Check if the particle is placed within the cell of interest, if not generate
! a random number until particle (x,y,z) falls inside cell
! The nth particle's initial cell location and regID is also initialized
! ******************************************************************************

        icg = pGrid%hex2CellGlob(icl)
        foundFlag = .FALSE.
        IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,icg) .EQV. .TRUE. ) THEN
         ! pPlag%aiv(AIV_PLAG_PIDINI,iPcl) = iPcl + nPclsSumReg
         ! pAiv(AIV_PLAG_REGINI,iPcl) = pRegion%iRegionGlobal
         ! pAiv(AIV_PLAG_ICELLS,iPcl) = icg

        ! Josh
          DO iPclJ = nPclsBeg,nPclsEnd
                pPlag%aiv(AIV_PLAG_PIDINI,iPclJ) = iPclJ + nPclsSumReg
          ENDDO

          pAiv(AIV_PLAG_REGINI,nPclsBeg:nPclsEnd) = pRegion%iRegionGlobal
          pAiv(AIV_PLAG_ICELLS,nPclsBeg:nPclsEnd) = icg
                

          foundFlag = .TRUE.

! Josh


        DO iPclJ = 1,SIZE(partDens)
              !  IF ( iPclJ .eq. size(partDens) ) THEN
              !      meanDia0 (1:5) =  (/ 5.00E-06_RFREAL, 10.00E-06_RFREAL, 50.00E-06_RFREAL, 100.00E-06_RFREAL, 200.00E-06_RFREAL /)
              !  ELSE
              !      meanDia0 (1:5) =  (/ 1.99E-03_RFREAL, 2.0028E-03_RFREAL, 2.0156E-03_RFREAL, 2.0228E-03_RFREAL, 2.03E-03_RFREAL /)
              !  ENDIF

           DO iPclK = 1,SIZE(meanDia0)
               
               IF ( nPclsEnd .le. size(partDens)*size(meanDia0) ) THEN 
           partNum = (iPclJ-1)*SIZE(meanDia0) + iPclK
               ELSE
           partNum = (iPclJ-1)*SIZE(meanDia0) + iPclK + (nPclsBeg-1)
               ENDIF             
 
              IF ( partNum .LE. nPclsEnd ) THEN 


          dia = meanDia0(iPclK) !+ Rand1Normal(meanDia,stdDev) ! Add stdDev later
          volPcl = global%pi*dia**3.0_RFREAL/6.0_RFREAL
! ******************************************************************************
! Compute mass of each particle and initialize each particle with x,y and z
! momentum in addition to total energy
! ******************************************************************************

          DO iCont = 1,pRegion%plagInput%nCont
            massRatio = pIniComp(iCont)*massFluxRatioSumR
!            pCv(pCvPlagMass(iCont),iPclJ) = pDens(iCont)*massRatio*volPcl
             pCv(pCvPlagMass(iCont),partNum) = partDens(iPclJ)*massRatio*volPcl
          END DO ! iCont
          heatCapSum = SUM(pCv(pCvPlagMass(:),partNum)*pSpcHeat(:))
          massSum    = SUM(pCv(pCvPlagMass(:),partNum))
          massSumR = 1.0_RFREAL/massSum

          pCv(CV_PLAG_XMOM,partNum) = massSum*u
          pCv(CV_PLAG_YMOM,partNum) = massSum*v
          pCv(CV_PLAG_ZMOM,partNum) = massSum*w
          pCv(CV_PLAG_ENER,partNum) = heatCapSum*T + &
                                   massSum*0.5_RFREAL* &
                                      ((massSumR*pCv(CV_PLAG_XMOM,partNum))**2.0_RFREAL+ &
                                       (massSumR*pCv(CV_PLAG_YMOM,partNum))**2.0_RFREAL+ &
                                       (massSumR*pCv(CV_PLAG_ZMOM,partNum))**2.0_RFREAL)
          volPclsSum = volPclsSum + volPcl

              ENDIF ! (partNum .LE. nPclsEnd)
                
           ENDDO ! iPclK
        ENDDO ! iPclJ

!          iPcl = iPcl + 1
!          pPlag%nPcls = pPlag%nPcls + 1
          pPlag%nPcls = nPclsEnd  ! Josh
            WRITE(802,*) 'pRegion%iRegionGlobal, icg, nPclsBeg, nPclsEnd,pIDBeg, pIDEnd, Xi, Yi',pRegion%iRegionGlobal,icg,nPclsBeg,nPclsEnd,&
                                                                 pPlag%aiv(AIV_PLAG_PIDINI,nPclsBeg),pPlag%aiv(AIV_PLAG_PIDINI,nPclsEnd), & 
                                                                                           xLoc,yLoc 
 
      ELSE   ! Increment loop counter if pcl falls outside cell
           ! WRITE(802,*) 'OUTSIDE CELL, pRegion%iRegionGlobal, icg, nPclsBeg, nPclsEnd, Xi, Yi',pRegion%iRegionGlobal,icg,nPclsBeg,nPclsEnd,xLoc,yLoc 

          loopCounter = loopCounter + 1
          IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN ! Prevent infinite loop
            CALL ErrorStop(global,ERR_INFINITE_LOOP ,__LINE__)
          END IF ! loopCounter

      END IF ! RFLU_ICT_TestInCell
    END DO ! iPcl

! ******************************************************************************
! Compute superparticle loading : spLoad for all particles within a given cell
! is taken to be the same
! ******************************************************************************

      vFrac = meanVfrac !+ Rand1Normal(meanVolFrac,stdDev) ! Add stdDev later
!      gaussAmp = pRegion%mixtInput%prepRealVal5
!      m = pRegion%mixtInput%prepIntVal1 ! Wave no in theta
!      n = pRegion%mixtInput%prepIntVal2 ! Wave number in z
!      perturb = gaussAmp*DCOS(m*the)
!      vFrac = meanVfrac*(1.0_RFREAL+perturb)

      spLoad = vFrac * pGrid%vol(icl)/volPclsSum
      DO iPcl = nPclsBeg,nPclsEnd
        pArv(ARV_PLAG_SPLOAD,iPcl) = spLoad
      END DO
    END IF ! nPclsPerCell > 0
  END DO ! icl

      ELSE ! Rocflu Particles !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        WRITE(*,*) 'Particle Init Type: Rocflu Particles'

  meanDia = pRegion%plagInput%iniRandDiamMin
  meanVfrac = pRegion%plagInput%iniRandSpLoadMin

  ! Temporarily set pPlag%nPcls = 0 
  ! It will be recomputed as we place particles one by one in this routine

  pPlag%nPcls = 0
  tol = 1.0E-14_RFREAL

! ******************************************************************************
! Sum up the number of particles in the all other regions except the current one
! to assign proper global initial PCL_ID.
! ******************************************************************************

  DO icl = 1, pGrid%nCells
    IF (INT(pGrid%nPclsPerCell(icl)) .GT. 0) THEN

      x = pGrid%cofg(XCOORD,icl)
      y = pGrid%cofg(YCOORD,icl)
      z = pGrid%cofg(ZCOORD,icl)

! ******************************************************************************
! Compute the min and max r,the,z coordiantes for each cell
! ******************************************************************************
      ! pGrid%dx,pGrid%dy,pGrid%dz computed in 
      ! libflu/RFLU_ComputeGridSpacingXdsp.F90       

      xMinCell = x - 0.5_RFREAL*pGrid%dx
      xMaxCell = x + 0.5_RFREAL*pGrid%dx
      yMinCell = y - 0.5_RFREAL*pGrid%dy
      yMaxCell = y + 0.5_RFREAL*pGrid%dy
      zMinCell = z   - 0.5_RFREAL*pGrid%dz
      zMaxCell = z   + 0.5_RFREAL*pGrid%dz

      volPclsSum = 0.0_RFREAL
      nPclsBeg = pPlag%nPcls + 1
      nPclsEnd = nPclsBeg + INT(pGrid%nPclsPerCell(icl)) - 1
! ******************************************************************************
! Use random number generator to generate a random location for each particle
! within a cell
! ******************************************************************************

      loopCounter = 0
      iPcl = nPclsBeg
      DO WHILE (iPcl .LE. nPclsEnd)
        xLoc = xMinCell &
             + Rand1Uniform(pRegion%randData) * (xMaxCell - xMinCell)
        yLoc = yMinCell &
             + Rand1Uniform(pRegion%randData) * (yMaxCell - yMinCell)
        IF (pRegion%mixtInput%dimens == 2) THEN
          zLoc = z
        ELSE
          zLoc = zMinCell &
               + Rand1Uniform(pRegion%randData) * (zMaxCell - zMinCell)
        END IF

        pCv(CV_PLAG_XPOS,iPcl) = xLoc
        pCv(CV_PLAG_YPOS,iPcl) = yLoc
        pCv(CV_PLAG_ZPOS,iPcl) = zLoc

! ******************************************************************************
! Check if the particle is placed within the cell of interest, if not generate
! a random number until particle (x,y,z) falls inside cell
! The nth particle's initial cell location and regID is also initialized
! ******************************************************************************

        icg = pGrid%hex2CellGlob(icl)
        foundFlag = .FALSE.
        IF ( RFLU_ICT_TestInCell(pRegion,xLoc,yLoc,zLoc,icg) .EQV. .TRUE. ) THEN
          pPlag%aiv(AIV_PLAG_PIDINI,iPcl) = iPcl + nPclsSumReg
          pAiv(AIV_PLAG_REGINI,iPcl) = pRegion%iRegionGlobal
          pAiv(AIV_PLAG_ICELLS,iPcl) = icg
          foundFlag = .TRUE.

        !BRAD debug yloc particle issue
       !if (yLoc .lt. 0.0_RFREAL) write(*,*) "DEBERG STOP BAD",yLoc
        !BRAD debug yloc particle issue 


          dia = meanDia !+ Rand1Normal(meanDia,stdDev) ! Add stdDev later
          volPcl = global%pi*dia**3.0_RFREAL/6.0_RFREAL
! ******************************************************************************
! Compute mass of each particle and initialize each particle with x,y and z
! momentum in addition to total energy
! ******************************************************************************

          DO iCont = 1,pRegion%plagInput%nCont
            massRatio = pIniComp(iCont)*massFluxRatioSumR
            pCv(pCvPlagMass(iCont),iPcl) = pDens(iCont)*massRatio*volPcl
          END DO ! iCont
          heatCapSum = SUM(pCv(pCvPlagMass(:),iPcl)*pSpcHeat(:))
          massSum    = SUM(pCv(pCvPlagMass(:),iPcl))
          massSumR = 1.0_RFREAL/massSum

          pCv(CV_PLAG_XMOM,iPcl) = massSum*u
          pCv(CV_PLAG_YMOM,iPcl) = massSum*v
          pCv(CV_PLAG_ZMOM,iPcl) = massSum*w
          pCv(CV_PLAG_ENER,iPcl) = heatCapSum*T + &
                                   massSum*0.5_RFREAL* &
                                      ((massSumR*pCv(CV_PLAG_XMOM,iPcl))**2.0_RFREAL+ &
                                       (massSumR*pCv(CV_PLAG_YMOM,iPcl))**2.0_RFREAL+ &
                                       (massSumR*pCv(CV_PLAG_ZMOM,iPcl))**2.0_RFREAL)
          volPclsSum = volPclsSum + volPcl

          iPcl = iPcl + 1
          pPlag%nPcls = pPlag%nPcls + 1

        ELSE   ! Increment loop counter if pcl falls outside cell

          loopCounter = loopCounter + 1
          IF ( loopCounter >= LIMIT_INFINITE_LOOP ) THEN ! Prevent infinite loop
            CALL ErrorStop(global,ERR_INFINITE_LOOP ,__LINE__)
          END IF ! loopCounter

        END IF ! RFLU_ICT_TestInCell
      END DO ! iPcl


! ******************************************************************************
! Compute superparticle loading : spLoad for all particles within a given cell
! is taken to be the same
! ******************************************************************************

      vFrac = meanVfrac !+ Rand1Normal(meanVolFrac,stdDev) ! Add stdDev later
!      gaussAmp = pRegion%mixtInput%prepRealVal5
!      m = pRegion%mixtInput%prepIntVal1 ! Wave no in theta
!      n = pRegion%mixtInput%prepIntVal2 ! Wave number in z
!      perturb = gaussAmp*DCOS(m*the)
!      vFrac = meanVfrac*(1.0_RFREAL+perturb)

!BRAD TANH PART CURTAIN
     ! if (x .lt. 0.0432) then !curt mid
      !  vFrac = 0.5_RFREAL * meanVfrac * (1.0_RFREAL + &
       !       DTANH(1.0E2_RFREAL*(x-(0.042 + 0.5_RFREAL*meanDia)))) !curtstart
       !else
        !vFrac = 0.5_RFREAL * meanVfrac * & (1.0_RFREAL - &
         !       DTANH(1.0E5_RFREAL*(x-(0.0435 - 0.5_RFREAL*meanDia))))!curt end

       ! vFrac = meanVfrac
      ! end if  
!BRAD TANH H

      spLoad = vFrac * pGrid%vol(icl)/volPclsSum

       !write(*,*) "NULLSPL: ",spLoad ,vFrac, pGrid%vol(icl), volPclsSum !BRAD DEBUG ASVF issue 

      DO iPcl = nPclsBeg,nPclsEnd
        pArv(ARV_PLAG_SPLOAD,iPcl) = spLoad
      END DO
    END IF ! nPclsPerCell > 0
  END DO ! icl


      ENDIF ! pMixtInput%prepRealVal27

! ******************************************************************************
! Destroy cell-to-face list
! ******************************************************************************

  CALL RFLU_DestroyCell2FaceList(pRegion)

! ******************************************************************************
! Destroy memory for number of particles per cell (only if a region contains
! pcls)
! NOTE memory deallocated in PLAG_RFLU_ComputeCellsContiningPcls if no particles
! present
! ******************************************************************************

  IF (pGrid%initPclPresent .EQV. .TRUE.) THEN
    DEALLOCATE(pGrid%nPclsPerCell,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'nPclsPerCell')
    END IF ! global%error
  END IF

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_LOW ) THEN
       WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME, &
            'Initializing particle solution done.'
  END IF ! global%verbLevel

!*******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE PLAG_RFLU_InitSolutionXdsp 

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_RFLU_InitSolutionShktb.F90,v $
! Revision 1.1  2015/08/12 03:55:40  brollin
! New Subroutine to Initialize Shock Tube related problems.
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
!******************************************************************************
