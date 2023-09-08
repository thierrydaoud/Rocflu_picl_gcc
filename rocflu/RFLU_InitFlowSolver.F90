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
! Purpose: Main initialization routine of Rocflu-MP.
!
! Description: None.
!
! Input: None.
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: RFLU_InitFlowSolver.F90,v 1.2 2017/01/23 16:48:31 rahul Exp $
!
! Copyright: (c) 2001-2006 by the University of Illinois
!
! ******************************************************************************

#ifdef GENX
SUBROUTINE RFLU_InitFlowSolver(globalGenx,initialTime,communicator,genxHandle, &
                               surfWinNameInput,volWinNameInput,handleObtain)
#else
SUBROUTINE RFLU_InitFlowSolver(casename,verbLevel,global,levels)
#endif

#ifdef PICL

USE ModRandom, ONLY: Rand1Uniform,Rand1Normal
#endif

USE ModDataTypes
#ifdef GENX
  USE ModGenx, ONLY: t_globalGenx
#endif  
  USE ModBndPatch, ONLY: t_patch
  USE ModDataStruct, ONLY: t_level,t_region
  USE ModMixture, ONLY: t_mixt_input  
  USE ModGlobal, ONLY: t_global
  USE ModGrid, ONLY: t_grid
  USE ModBndPatch, ONLY: t_patch
  USE ModError
  USE ModMPI
  USE ModParameters  

  USE RFLU_ModABC
  USE RFLU_ModAxisymmetry     
  USE RFLU_ModBFaceGradAccessList
  USE RFLU_ModBoundLists   
  USE RFLU_ModBoundXvUtils
  USE RFLU_ModNSCBC, ONLY: RFLU_NSCBC_DecideHaveNSCBC
  USE RFLU_ModCellMapping
  USE RFLU_ModCommLists
  USE RFLU_ModDifferentiationCells, ONLY: RFLU_ComputeGradCellsGGScalar, &
                                          RFLU_ComputeGradCellsGGVector, &
                                          RFLU_AXI_ComputeGradCellsGGScalar, &
                                          RFLU_AXI_ComputeGradCellsGGVector
  USE RFLU_ModDimensionality  
  USE RFLU_ModDimensions
  USE RFLU_ModEdgeList   
  USE RFLU_ModFaceList
  USE RFLU_ModForcesMoments
  USE RFLU_ModGeometry
  USE RFLU_ModGFM
  USE RFLU_ModGlobalIds
  USE RFLU_ModGridSpeedUtils
  USE RFLU_ModHouMahesh
  USE RFLU_ModHypre
  USE RFLU_ModInCellTest
  USE RFLU_ModInterpolation
  USE RFLU_ModMPI
  USE RFLU_ModMovingFrame, ONLY: RFLU_MVF_CreatePatchVelAccel, &
                                 RFLU_MVF_InitPatchVelAccel, &
                                 RFLU_MVF_ReadPatchVelAccel
  USE RFLU_ModOLES
  USE RFLU_ModPatchCoeffs
  USE RFLU_ModPatchUtils
  USE RFLU_ModProbes
  USE RFLU_ModReadBcInputFile
  USE RFLU_ModReadWriteAuxVars
  USE RFLU_ModReadWriteBcDataFile
  USE RFLU_ModReadWriteFlow   
  USE RFLU_ModReadWriteGrid 
  USE RFLU_ModReadWriteGridSpeeds
  USE RFLU_ModRegionMapping
  USE RFLU_ModRenumberings
  USE RFLU_ModStencilsBFaces
  USE RFLU_ModStencilsCells
  USE RFLU_ModStencilsFaces
  USE RFLU_ModStencilsUtils
  USE RFLU_ModSymmetryPeriodic, ONLY: RFLU_SYPE_HaveSyPePatches, & 
                                      RFLU_SYPE_ReadTransforms
  USE RFLU_ModVertexLists
  USE RFLU_ModWeights
    
#ifdef GENX
  USE RFLU_ModGENXAdmin
  USE RFLU_ModGENXIO, ONLY: RFLU_GENX_GetDataInterf, & 
                            RFLU_GENX_GetDimensionsDerived, &
                            RFLU_GENX_GetGlobalData
#endif  
    
#ifdef PETSC
  USE RFLU_ModPETScAdmin
!  USE RFLU_ModPETScPoisson
  USE RFLU_ModPETScNewtonKrylov
#endif
  
#ifdef PLAG
  USE ModInterfacesLagrangian, ONLY: PLAG_RFLU_ComputeVolFrac
  USE PLAG_ModSurfStats, ONLY: PLAG_CreateSurfStats, & 
                               PLAG_DecideHaveSurfStats, & 
                               PLAG_ReadSurfStatsWrapper
  ! Subbu - Check for hardcode pcl tracking
  USE PLAG_ModParameters
  ! Subbu - End Check for hardcode pcl tracking
#endif    
    
  USE RFLU_ModRelatedPatches, ONLY: RFLU_RELP_TransformWrapper
    
  USE ModInterfaces, ONLY: RFLU_AllocateMemoryTbc, &
                           RFLU_AllocateMemoryWrapper, &
                           RFLU_BuildDataStruct, & 
                           RFLU_ComputeIntegrals1245OLES, &
                           RFLU_ComputeIntegralValues, &
                           RFLU_CreateGrid, &
                           RFLU_DecideNeedBGradFace, &
                           RFLU_DecideNeedStencils, &
                           RFLU_DecideNeedWeights, &
                           RFLU_GetUserInput, & 
                           RFLU_InitGlobal, &
                           RFLU_MoveGridWrapper, &
                           RFLU_OpenConverFile, &
                           !begin BBR
                           RFLU_OpenPMFile, &
                           RFLU_OpenIntegFile, &
                           !end BBR
                           RFLU_OpenStatsFileOLES, &
                           RFLU_OpenTotalMassFile, &
                           RFLU_PrintFlowInfo, &
                           RFLU_PrintGridInfo, & 
                           RFLU_PrintHeader, &
                           RFLU_RandomInit, &
                           RFLU_ReadIntegrals1245OLES, &
                           RFLU_ReadRestartInfo, &
                           RFLU_ReadTbcInputFile, &
                           RFLU_SetModuleType, & 
                           RFLU_SetMoveGridOptions, & 
                           RFLU_SetRestartTimeFlag, &
                           RFLU_SetVarsContWrapper, &
                           RFLU_SetVarsDiscWrapper, & 
                           RFLU_SetVarInfoWrapper, & 
                           RFLU_SetVarsWrapper, & 
                           RFLU_WriteIntegrals1245OLES, &
                           RFLU_WriteVersionString, &
                           WriteTotalMass
                                                
#ifdef GENX
  USE ModInterfaces, ONLY: RFLU_CheckCouplingInput, &
                           RFLU_GetBoundaryValues, &
                           RFLU_PutBoundaryValues
  USE RFLU_ModGENXTools, ONLY: RFLU_GENX_InitBFLAG
#endif
#ifdef STATS
  USE ModStatsRoutines, ONLY: StatMapping, InitStatistics
#endif
  
  IMPLICIT NONE

#ifdef GENX
  INCLUDE 'roccomf90.h'
#endif

#ifdef PICL
!DEC$ NOFREEFORM
!#include "/ascldap/users/shbrine/ppiclF/ppiclF_new/ppiclF/source/PPICLF_USER.h"
!#include "/ascldap/users/shbrine/ppiclF/ppiclF_new/ppiclF/source/PPICLF_STD.h"
#include "/home/t.daoud/Rocflu_picl_gcc/ppiclf/source/PPICLF_USER.h"
#include "/home/t.daoud/Rocflu_picl_gcc/ppiclf/source/PPICLF_STD.h"
!#include "/home/rahul.koneru/codes/Rocflu-ppiclF/ppiclF/source/PPICLF"
!DEC$ FREEFORM
#endif

!
! ******************************************************************************
! Arguments
! ******************************************************************************

#ifdef GENX
  CHARACTER(*), INTENT(IN) :: surfWinNameInput,volWinNameInput
  INTEGER, INTENT(IN) :: communicator,genxHandle,handleObtain
  DOUBLE PRECISION, INTENT(IN) :: initialTime
  TYPE(t_globalGenx), POINTER :: globalGenx
#else
  CHARACTER(CHRLEN), INTENT(IN) :: casename
  INTEGER, INTENT(IN) :: verbLevel
  TYPE(t_global), POINTER :: global
#endif
  TYPE(t_level), POINTER :: levels(:)  
  
! ******************************************************************************
! Locals
! ******************************************************************************

  CHARACTER(CHRLEN) :: RCSIdentString,msg
  LOGICAL :: fileExists,moveGrid
  INTEGER :: errorFlag,flag,icg,iPatch,iReg,iVar
  REAL(RFREAL) :: currentTime
  TYPE(t_grid) :: grid
  TYPE(t_grid), POINTER :: pGrid
  TYPE(t_level), POINTER :: pLevel
  TYPE(t_mixt_input), POINTER :: pMixtInput  
  TYPE(t_patch), POINTER :: pPatch
  TYPE(t_region), POINTER :: pRegion

#ifdef GENX
  INTEGER :: verbLevel
  CHARACTER(CHRLEN) :: casename
  TYPE(t_global), POINTER :: global  
#endif 

  ! Subbu - Vars cyldet case
   REAL(RFREAL) :: nr,nTol,nx,ny,nz,rad,the,theLow,theUpp,Tol,z
   INTEGER :: ifl,bcSlipFlag,bcPeriLeftFlag,bcPeriRightFlag,bcOutfFlag
   TYPE(t_region), POINTER :: pRegionSerial
  ! Subbu - End vars cyldet case

  ! Rahul - TEMP Debug
   INTEGER :: iPcl
  ! Rahul - end
#ifdef PICL
   CHARACTER(CHRLEN) :: endString
   LOGICAL :: notfoundFlag,pf_fluidInit,wall_exists
   INTEGER :: i,npart,nCells,lx,ly,lz,vi,vii,ii,i3,i4,i5,i6,i7,jj,kk,&
              loopCounter,&
              ipart,icl,npart_local,i_global,i_global_min,i_global_max
   INTEGER :: PPC,numPclCells      
   REAL(RFREAL) :: dp_min,dp_max,rhop,tester,ratio,total_vol,filter,xMinCurt,&
                   xMaxCurt,yMinCurt,yMaxCurt,xMinCell,xMaxCell,yMinCell,&
                   yMaxCell,zMinCell,zMaxCell,x,vFrac,volpclsum,xLoc,yLoc,zLoc,yL, &
                   zpf_factor,xpf_factor,vol,diam,x0,y0,z0,ymax,dy,ddiam,dcut, &
                   dRho, dSigma, dMu, uid
   REAL(KIND=8) :: y(PPICLF_LRS, PPICLF_LPART), &
                   rprop(PPICLF_LRP, PPICLF_LPART), sigma, mutab
   REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: xGrid, yGrid, zGrid,vfP
   REAL(RFREAL),ALLOCATABLE,DIMENSION(:) :: xData,yData,zData,rData,dumData     
   REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: volp,SPL 
   REAL(KIND=8), DIMENSION(3) :: tpw1,tpw2,tpw3         
   ! Sam
   REAL(RFREAL) :: da_axi, y_axi1, y_axi0
   INTEGER :: iFile
   CHARACTER(CHRLEN) :: iFileName
#endif

! ******************************************************************************
! Start, initialize some variables
! ******************************************************************************

  RCSIdentString = '$RCSfile: RFLU_InitFlowSolver.F90,v $ $Revision: 1.2 $'
 
  moveGrid = .FALSE.
 
! ******************************************************************************
! Set global pointer and initialize global type, register function
! ****************************************************************************** 
 
#ifdef GENX
  global => globalGenx%global 
     
  CALL RFLU_InitGlobal(" ",VERBOSE_LOW,communicator,global)
  CALL RFLU_GENX_ReadCtrlFile(global)
  CALL RFLU_GENX_StoreCommunicator(global,communicator)

! ******************************************************************************
! Initialize MPI - ROCSTAR
! ****************************************************************************** 

  CALL MPI_Comm_size(global%mpiComm,global%nProcAlloc,errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
  END IF ! global%error
    
  CALL MPI_Comm_rank(global%mpiComm,global%myProcid,errorFlag)  
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
  END IF ! global%error

  CALL RFLU_GENX_StoreNamesHandles(global,surfWinNameInput,volWinNameInput, & 
                                   handleObtain)

  CALL COM_set_verbose(global%verbLevelCOM)

  global%timeStamp   = initialTime
  global%currentTime = initialTime
#else
  CALL RFLU_InitGlobal(casename,verbLevel,MPI_COMM_WORLD,global)

! ******************************************************************************
! Initialize MPI Native
! ****************************************************************************** 

  CALL MPI_Init(errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
  END IF ! global%error
  
  CALL MPI_Comm_size(global%mpiComm,global%nProcAlloc,errorFlag)
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
  END IF ! global%error
    
  CALL MPI_Comm_rank(global%mpiComm,global%myProcid,errorFlag)  
  global%error = errorFlag
  IF ( global%error /= ERR_NONE ) THEN 
    CALL ErrorStop(global,ERR_MPI_OUTPUT,__LINE__)
  END IF ! global%error

#ifdef PICL
  call ppiclf_comm_InitMPI(global%mpiComm,global%myProcid,global%nProcAlloc)
#endif
#endif

  CALL RFLU_SetModuleType(global,MODULE_TYPE_SOLVER)

  CALL RegisterFunction(global,'RFLU_InitFlowSolver',__FILE__)

#ifndef GENX
#ifdef ROCPROF
! ******************************************************************************
! Initialize Rocprof. NOTE done here so everything can be profiled and also bcos
! RFLU_SetVars is called here.
! ****************************************************************************** 

  CALL Rocprof_Init("Rocflu",global%myProcid)
#endif
#endif

! ******************************************************************************
! Print header and check for stop file 
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC ) THEN 
    CALL RFLU_WriteVersionString(global)
  
    IF ( global%verbLevel /= VERBOSE_NONE ) THEN
      CALL RFLU_PrintHeader(global)
    END IF ! global%verbLevel
  END IF ! global

  IF ( global%myProcid == MASTERPROC ) THEN 
    INQUIRE(FILE="STOP",EXIST=fileExists)
    IF ( fileExists .EQV. .TRUE. ) THEN
      CALL ErrorStop(global,ERR_STOPFILE_FOUND,__LINE__)
    END IF ! fileExists
  END IF ! global%myProcid  

! ******************************************************************************
! Read processor mapping file, prepare data structure, set level pointer
! ******************************************************************************
 
  CALL RFLU_ReadRegionMappingFile(global,MAPFILE_READMODE_PEEK,global%myProcId)
  CALL RFLU_CreateRegionMapping(global,MAPTYPE_REG)
  CALL RFLU_ReadRegionMappingFile(global,MAPFILE_READMODE_ALL,global%myProcId)
  
  IF ( global%nRegions == 1 ) THEN 
    CALL RFLU_ImposeRegionMappingSerial(global)
  END IF ! global
  
  write(*,*) "global%nProcs" , global%nProcs ,"global%nProcAlloc" , global%nProcAlloc
  IF ( global%nProcs /= global%nProcAlloc ) THEN 
    CALL ErrorStop(global,ERR_PROC_MISMATCH,__LINE__)
  END IF ! global%nProcs
  
  CALL RFLU_BuildDataStruct(global,levels)
  CALL RFLU_ApplyRegionMapping(global,levels)
  CALL RFLU_DestroyRegionMapping(global,MAPTYPE_REG)
   
 
#ifdef GENX    
  globalGenx%levels => levels
  pLevel => globalGenx%levels(1) ! single-level grids for now      
#else
  pLevel => levels(1) ! single-level grids for now   
#endif

#ifdef GENX
! ******************************************************************************
! Create windows and new attributes. NOTE must only be done once for every 
! process.
! ******************************************************************************
  
  pRegion => pLevel%regions(1)
  CALL RFLU_GENX_CreateWindows(pRegion,communicator)
  CALL RFLU_GENX_CreateAttrWrapper(pRegion)
  
! ******************************************************************************
! Load modules
! ******************************************************************************

  CALL Rocmop_load_module(TRIM(global%winName)//'-MOP')  
#endif           
           
! ******************************************************************************
! Read user input
! ******************************************************************************

  CALL RFLU_GetUserInput(pLevel%regions)

! ******************************************************************************
! Initialize random number generator
! ******************************************************************************

  CALL RFLU_RandomInit(pLevel%regions)

! ******************************************************************************
! Read restart info. NOTE must not read restart info file for GENx runs because 
! will overwrite global%currentTime.
! ******************************************************************************
  
#ifndef GENX
  CALL RFLU_ReadRestartInfo(global)
#endif
  CALL RFLU_SetRestartTimeFlag(global)
  
! ******************************************************************************
! Read dimensions file
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_ReadDimensionsWrapper(pRegion)     
  END DO ! iReg


! ******************************************************************************
! Determine whether have moving grid
! ******************************************************************************

  moveGridLoop: DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    
    IF ( pRegion%mixtInput%movegrid .EQV. .TRUE. ) THEN
      moveGrid = .TRUE.       
      EXIT moveGridLoop
    END IF ! pRegion
  END DO moveGridLoop

! ******************************************************************************
! Allocate memory for grid and borders
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)

    CALL RFLU_CreateGrid(pRegion)
    CALL RFLU_COMM_CreateBorders(pRegion,CREATE_BORDERS_MODE_DIM_KNOWN)
  END DO ! iReg
 
! ******************************************************************************
! Read boundary condition file, must be done after creation of grid, because 
! grid and patches are created there.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    
    IF ( pRegion%grid%nPatches > 0 ) THEN      
      CALL RFLU_ReadBcInputFileWrapper(pRegion)
      CALL RFLU_AllocateMemoryTbc(pRegion,iReg)
      CALL RFLU_ReadTbcInputFile(pRegion)
    END IF ! pRegion%grid%nPatches    
  END DO ! iReg

#ifdef GENX
! ******************************************************************************
! Create surface and and register volume and surface grid. NOTE must be done 
! after having read comm lists because of pconn.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)

    CALL RFLU_GENX_GetDimensionsDerived(pRegion)
    CALL RFLU_GENX_CreateGridSurf(pRegion)    
    CALL RFLU_GENX_CreatePConn(pRegion)
    CALL RFLU_GENX_RegisterGridVol(pRegion)
    CALL RFLU_GENX_RegisterGridSurf(pRegion)
  END DO ! iReg  
#endif
  
! ******************************************************************************
! Read grid file
! ******************************************************************************
  
  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_ReadGridWrapper(pRegion)
  END DO ! iReg
 
! ==============================================================================
! Non-dimensionalize mesh for SOLV_IMPLICIT_HM
! ==============================================================================

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    CALL RFLU_HM_ConvGridCoordD2ND(pRegion)
  END IF
 
! ******************************************************************************
! Print grid information
! ******************************************************************************

  IF ( global%myProcid == MASTERPROC .AND. &
       global%verbLevel > VERBOSE_NONE ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)
      CALL RFLU_PrintGridInfo(pRegion)
    END DO ! iReg
  END IF ! global%verbLevel

!*******************************************************************************
! Populate Time Zooming Parameter
!*******************************************************************************
   
#ifdef GENX
  CALL RFLU_GENX_GetGlobalData(global)
#endif

! ******************************************************************************
! Build data stucture, part 1
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg) 
    CALL RFLU_CreateCellMapping(pRegion)
    CALL RFLU_ReadLoc2GlobCellMapping(pRegion)
    CALL RFLU_BuildGlob2LocCellMapping(pRegion)      
  END DO ! iReg 

! ******************************************************************************
! Create and read comm lists, set proc ids for borders. Must be done before 
! face list is constructed because reorientation of actual-virtual faces
! requires knowledge of communication lists.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    pGrid   => pRegion%grid   
    
    IF ( pGrid%nBorders > 0 ) THEN 
      CALL RFLU_COMM_CreateCommLists(pRegion)
      CALL RFLU_COMM_ReadCommLists(pRegion)
      CALL RFLU_COMM_GetProcLocRegIds(pRegion)
    END IF ! pRegion%grid%nBorders
  END DO ! iReg  

! ******************************************************************************
! Build data stucture, part 2. Building the face list must be done after the 
! communication lists have been updated to take into account renumbering of 
! cells. Build cell-to-face list anyway, although only needed when have probes 
! or when running with higher-order scheme.
! ******************************************************************************
 
  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg) 
    CALL RFLU_CreateBVertexLists(pRegion)
    CALL RFLU_BuildBVertexLists(pRegion)        
  END DO ! iReg

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg) 
    CALL RFLU_CreateFaceList(pRegion)
    CALL RFLU_BuildFaceList(pRegion)
    CALL RFLU_RenumberBFaceLists(pRegion)

! TEMPORARY
!    IF ( pRegion%grid%nCellsTot > pRegion%grid%nCells ) THEN 
!      CALL RFLU_RNMB_CreatePBF2SBFMap(pRegion)
!      CALL RFLU_RNMB_CreatePC2SCMap(pRegion)
!      CALL RFLU_RNMB_CreatePV2SVMap(pRegion)
!      CALL RFLU_RNMB_ReadPxx2SxxMaps(pRegion)
!      CALL RFLU_ReorientFaces(pRegion) 
!! TEMPORARY
!!      CALL RFLU_RNMB_DestroyPBF2SBFMap(pRegion)
!!      CALL RFLU_RNMB_DestroyPC2SCMap(pRegion)
!!      CALL RFLU_RNMB_DestroyPV2SVMap(pRegion)
!! END TEMPORARY 
!    END IF ! pGrid%nCellsTot
! END TEMPORARY
  END DO ! iReg 

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_CreateCell2FaceList(pRegion)
    CALL RFLU_BuildCell2FaceList(pRegion)
  END DO ! iReg

#ifdef PLAG
  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)

    IF ( (pRegion%global%plagUsed .EQV. .TRUE.) .AND. & 
         (pRegion%grid%nFacesAV > 0) ) THEN
      CALL RFLU_CreateAVFace2BorderList(pRegion)
      CALL RFLU_BuildAVFace2BorderList(pRegion)
      
      CALL RFLU_CreateAVFace2PatchList(pRegion)
      CALL RFLU_BuildAVFace2PatchList(pRegion)      
    END IF ! pRegion%global%plagUsed
  END DO ! iReg
#endif
 
! ******************************************************************************
! Build data stucture, part 3
! ******************************************************************************

  IF ( moveGrid .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)
      CALL RFLU_CreateEdgeList(pRegion)       
      CALL RFLU_BuildEdgeList(pRegion)      
      CALL RFLU_CreateEdge2CellList(pRegion)              
      CALL RFLU_BuildEdge2CellList(pRegion)
      CALL RFLU_DestroyEdge2CellList(pRegion)      
    END DO ! iReg 
  END IF ! moveGrid

! ******************************************************************************
! Check topology
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal 
    pRegion => pLevel%regions(iReg)              
    CALL RFLU_123D_CheckTopology(pRegion)
  END DO ! iReg  

! ******************************************************************************
! Build boundary-face gradient access list for viscous flows
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal 
    pRegion => pLevel%regions(iReg) 
    
    IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN          
      CALL RFLU_CreateBFaceGradAccessList(pRegion)
      CALL RFLU_BuildBFaceGradAccessList(pRegion)
    END IF ! pRegion%mixtInput
  END DO ! iReg  

! ******************************************************************************
! Allocate memory. NOTE must be done after having built face lists and after 
! reading boundary-condition file because of particle module.
! ******************************************************************************
 
  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_AllocateMemoryWrapper(pRegion)
    
    CALL RFLU_CreatePatchCoeffs(pRegion)

    IF ( global%forceFlag .EQV. .TRUE. ) THEN 
      CALL RFLU_CreateForcesMoments(pRegion)
      CALL RFLU_CreateGlobalThrustFlags(pRegion)
    END IF ! global%forceFlag

   IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
     CALL RFLU_MVF_CreatePatchVelAccel(pRegion)
   END IF ! global%mvFrameFlag
  END DO ! iReg
    
! ------------------------------------------------------------------------------
! Set Global thrust flags from local thrust flags 
! ------------------------------------------------------------------------------

  CALL RFLU_SetGlobalThrustFlags(pLevel%regions)

! ------------------------------------------------------------------------------
! Allocate memory for boundary variables if required, i.e., if NSCBC==1 
! ------------------------------------------------------------------------------

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
      CALL RFLU_BXV_CreateVarsCv(pRegion)
      CALL RFLU_BXV_CreateVarsDv(pRegion)
      CALL RFLU_BXV_CreateVarsTStep(pRegion)
    END IF ! RFLU_NSCBC_DecideHaveNSCBC(pRegion)    

#ifdef PLAG
    IF ( global%plagUsed .EQV. .TRUE. ) THEN 
      CALL PLAG_CreateSurfStats(pRegion)
    END IF ! global%plagUsed
#endif    
  END DO ! iReg

! ******************************************************************************
! Set variable info. NOTE must be done after allocation of memory for solution.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_SetVarInfoWrapper(pRegion)
  END DO ! iReg  

#ifdef STATS
! ******************************************************************************
! Statistics mapping (must be done before RFLU_GENX_RegisterData)
! ******************************************************************************

  IF (global%flowType == FLOW_UNSTEADY .AND. global%doStat==ACTIVE) THEN
    CALL StatMapping( global )
  ENDIF
#endif

#ifdef GENX
! ******************************************************************************
! Register data with Roccom for interaction in GENX. NOTE must be done after 
! memory allocation.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    
    CALL RFLU_GENX_RegisterDataWrapper(pRegion)
    CALL RFLU_GENX_CreateWindowsDone(pRegion)
  END DO ! iReg  
#endif  

! ******************************************************************************
! Compute geometry
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)   
    CALL RFLU_CreateGeometry(pRegion)    
    CALL RFLU_BuildGeometry(pRegion)
    CALL RFLU_123D_CheckGeometryWrapper(pRegion)    
  END DO ! iReg 

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)   
    CALL RFLU_ComputePatchNormalsLocal(pRegion)       
  END DO ! iReg 

  CALL RFLU_ComputePatchNormalsGlobal(pLevel%regions)
  
#ifdef PLAG
 ! Subbu - Compute cell-size, nrad,nthe,nz in each region
  IF ( (pRegion%global%plagUsed .EQV. .TRUE.) .AND. &
       (pRegion%plagInput%findPclMethod == FIND_PCL_METHOD_HARDCODE) ) THEN
  SELECT CASE ( global%casename )
    CASE("cyldet")
     Tol = 1.0E-8_RFREAL

     DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)
      pMixtInput => pRegion%mixtInput
      DO icg = 1,pGrid%nCells
        IF (pGrid%hex2CellGlob(icg) == 1) THEN
          pGrid%radMin = SQRT(pGrid%cofg(XCOORD,icg)**2.0_RFREAL + &
                               pGrid%cofg(YCOORD,icg)**2.0_RFREAL)
          pGrid%theMin = ATAN2(pGrid%cofg(YCOORD,icg),pGrid%cofg(XCOORD,icg))
          pGrid%zMin   = pGrid%cofg(ZCOORD,icg)
          pGrid%SgntheMin = INT(DSIGN(1.0_RFREAL,pGrid%theMin))
        END IF
        IF (pGrid%hex2CellGlob(icg) == 2) THEN
          the = ATAN2(pGrid%cofg(YCOORD,icg),pGrid%cofg(XCOORD,icg))
          pGrid%dthe = the -pGrid%theMin
        END IF
      END DO
          
      DO icg = 2,pGrid%nCells
        rad = SQRT(pGrid%cofg(XCOORD,icg)**2.0_RFREAL + &
                   pGrid%cofg(YCOORD,icg)**2.0_RFREAL)
        the = ATAN2(pGrid%cofg(YCOORD,icg),pGrid%cofg(XCOORD,icg))
        z   = pGrid%cofg(ZCOORD,icg)
        IF (ABS(the-pGrid%theMin) .LE. Tol .AND. &
            ABS(z-pGrid%zMin) .LE. Tol) THEN
          pGrid%Imax = icg - 1
          pGrid%drad = rad -pGrid%radMin
          EXIT
        END IF
      END DO
    
      IF (pRegion%mixtInput%dimens == 2) THEN
        pGrid%dz = 2.0_RFREAL*pGrid%zMin
        pGrid%Jmax = 1
        pGrid%Kmax = pGrid%nCells/(pGrid%Imax*pGrid%Jmax)
      ELSE
        DO icg = 2,pGrid%nCells
          rad = SQRT(pGrid%cofg(XCOORD,icg)**2.0_RFREAL + &
                     pGrid%cofg(YCOORD,icg)**2.0_RFREAL)
          the = ATAN2(pGrid%cofg(YCOORD,icg),pGrid%cofg(XCOORD,icg))
          z   = pGrid%cofg(ZCOORD,icg)
          IF (ABS(the-pGrid%theMin) .LE. Tol .AND. &
              ABS(rad-pGrid%radMin) .LE. Tol) THEN
             pGrid%dz = z - pGrid%zMin
             pGrid%Kmax = (icg-1)/pGrid%Imax
             pGrid%Jmax = pGrid%nCells/(pGrid%Imax*pGrid%Kmax)
             EXIT
          END IF
        END DO
      END IF 

      pRegionSerial => levels(1)%regions(0)
      CALL RFLU_ReadDimensions(pRegionSerial)
      CALL RFLU_CreateGrid(pRegionSerial)
      global%nz   = pRegionSerial%grid%nCells/pRegionSerial%patches(3)%nBQuads
      global%nrad = pRegionSerial%grid%nCells/pRegionSerial%patches(1)%nBQuads
      global%nthe = pRegionSerial%grid%nCells/(global%nz * global%nrad)
     
      IF (pRegion%grid%Imax .LT. global%nthe) THEN
        pGrid%nPartThetaFlag = 1
        pGrid%ImaxVirt = pGrid%Imax + 6
      ELSE
        pGrid%nPartThetaFlag = 0
        pGrid%ImaxVirt = pGrid%Imax 
      END IF

      theLow = pGrid%theMin
      icg    = pGrid%nCells 
      theUpp = ATAN2(pGrid%cofg(YCOORD,icg),pGrid%cofg(XCOORD,icg)) 
      IF (theLow .LT. 0.0_RFREAL) THEN
        theLow = theLow + 2.0_RFREAL*global%pi      
      END IF
      IF (theUpp .LT. 0.0_RFREAL) THEN
        theUpp = theUpp + 2.0_RFREAL*global%pi      
      END IF
      pGrid%theLowLim = ABS(MIN(theLow,theUpp) + 0.5_RFREAL*pGrid%dthe) ! dthe is always -ve
      pGrid%theUppLim = ABS(MAX(theLow,theUpp) - 0.5_RFREAL*pGrid%dthe)

      bcSlipFlag = 0
      bcOutfFlag = 0
       
      IF (pRegion%grid%nPatches > 0) THEN
        DO iPatch = 1,pRegion%grid%nPatches
          pPatch => pRegion%patches(iPatch)
          IF (pPatch%bcType == BC_SLIPWALL .AND. pPatch%iPatchGlobal == 1) THEN
            bcSlipFlag = 1
          END IF
          IF (pPatch%bcType == BC_OUTFLOW) THEN
            bcOutfFlag = 1
          END IF
        END DO
      END IF

      IF ( bcOutfFlag == 1 .OR. bcSlipFlag == 1 ) THEN
        pGrid%KmaxVirt = pGrid%Kmax+3
      ELSE
        pGrid%KmaxVirt = pGrid%Kmax+6
      END IF
      
      icg  = pGrid%nCells + 1
      pGrid%radMinVirtExt = SQRT(pGrid%cofg(XCOORD,icg)**2.0_RFREAL + &
                                 pGrid%cofg(YCOORD,icg)**2.0_RFREAL)
      pGrid%theMinVirtExt = ATAN2(pGrid%cofg(YCOORD,icg),pGrid%cofg(XCOORD,icg))
      pGrid%zMinVirtExt   = pGrid%cofg(ZCOORD,icg)

      ! Add 3 layer of cells if there is sharing in radial direction (r < rMin)
      IF (ABS(pGrid%radMinVirtExt - pGrid%radMin) .LE. Tol) THEN
        pGrid%nCellsShare = 0 ! r-dirn - one layer only
      ELSE
        pGrid%nCellsShare = 3 * pGrid%ImaxVirt ! r-dirn - one layer only
      END IF

      ! Add  virtual cells if there is sharing in z-direction (z < zMin)
      IF (ABS(pGrid%zMinVirtExt - pGrid%zMin) .LE. Tol) THEN
        pGrid%nCellsVirt = 0
      ELSE
        pGrid%nCellsVirt = 3 * pGrid%ImaxVirt * pGrid%KmaxVirt
      END IF

      IF (pGrid%nPartThetaFlag == 1) THEN
        icg = pGrid%nCells + pGrid%nCellsVirt + pGrid%nCellsShare + 1
        pGrid%theMinVirtInt = ATAN2(pGrid%cofg(YCOORD,icg),pGrid%cofg(XCOORD,icg))
      END IF

      pGrid%SgntheMinVirtExt = INT(DSIGN(1.0_RFREAL,pGrid%theMinVirtExt))
      pGrid%SgntheMinVirtInt = INT(DSIGN(1.0_RFREAL,pGrid%theMinVirtInt))

     END DO

     nTol = 1.0E-14_RFREAL
     pGrid%radMaxflag = 0
     pGrid%zMaxflag = 0
     pGrid%zMaxflag = 0
     DO iPatch=1,pGrid%nPatches
       pPatch => pRegion%patches(iPatch)
       DO ifl = 1,pPatch%nBFaces
         nx = pPatch%fn(XCOORD,ifl)
         ny = pPatch%fn(YCOORD,ifl)
         nz = pPatch%fn(ZCOORD,ifl)
         nr = DSQRT(nx**2.0_RFREAL + ny**2.0_RFREAL)
         IF ( DABS(nr - 1.0_RFREAL) .LE. nTol ) THEN
           pGrid%radMaxflag = 1
           pGrid%radMaxGlob = pGrid%radMin+(pGrid%Kmax-1)*pGrid%drad
           EXIT
         END IF
       END DO
       
       IF (pRegion%mixtInput%dimens == 3) THEN
         DO ifl = 1,pPatch%nBFaces
           nz = pPatch%fn(ZCOORD,ifl)
           IF ( DABS(DABS(nz) - 1.0_RFREAL) .LE. nTol ) THEN
             IF ( INT(DSIGN(1.0_RFREAL,nz)) == 1) THEN
               pGrid%zMaxflag = 1
               pGrid%zMaxGlob = pGrid%zMin+(pGrid%Jmax-1)*pGrid%dz
               EXIT
             ELSE
               pGrid%zMinflag = 1
               pGrid%zMinGlob = pGrid%zMin+(pGrid%Jmax-1)*pGrid%dz
               EXIT
             END IF
           END IF
         END DO
       END IF
     END DO

   IF (1==2) THEN
    IF (pRegion%iRegionGlobal == 2) THEN
     DO icg = pRegion%grid%nCells+1,pGrid%nCellsTot
       rad = SQRT(pRegion%grid%cofg(XCOORD,icg)**2.0_RFREAL + &
                  pRegion%grid%cofg(YCOORD,icg)**2.0_RFREAL)
       the = ATAN2(pRegion%grid%cofg(YCOORD,icg),pRegion%grid%cofg(XCOORD,icg))
       z   = pRegion%grid%cofg(ZCOORD,icg)
       !iBorder = pGrid%vc2border(1,icg-pGrid%nCells)
       WRITE(*,*) 'r,t,z=',pRegion%iRegionGlobal,icg,rad,the*180/global%pi,z
     END DO
    END IF
   END IF 

  END SELECT
  END IF ! plag find pcl hardcode
  ! Subbu - End compute cell-size,nrad,nthe,nz in each region for cyldet case
#endif
  
! TEMPORARY: commented by Manoj as there is error in subroutine  
!  DO iReg = 1,global%nRegionsLocal
!    pRegion => pLevel%regions(iReg)   
!  CALL RFLU_CheckPatchBcConsistency(pRegion)
!  END DO ! iReg 
! END TEMPORARY

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    pMixtInput => pRegion%mixtInput

    IF ( pMixtInput%fluidModel == FLUID_MODEL_INCOMP ) THEN
      CALL RFLU_CreateFaceDist(pRegion)
      CALL RFLU_ComputeFaceDist(pRegion)
    END IF ! pMixtInput%fluidModel
  END DO ! iReg

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg) 
    
    IF ( RFLU_SYPE_HaveSyPePatches(pRegion) .EQV. .TRUE. ) THEN   
      CALL RFLU_SYPE_ReadTransforms(pRegion)       
    END IF ! RFLU_SYPE_HaveSyPePatches
  END DO ! iReg


! ******************************************************************************
! Build stencils. NOTE This must be done after the geometry is computed 
! because the stencils are now tested for singularity, which requires the
! computation of the weights and hence the geometry.
! ******************************************************************************

  IF ( RFLU_DecideNeedStencils(pRegion) .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal 
      pRegion => pLevel%regions(iReg)
      pMixtInput => pRegion%mixtInput 

      CALL RFLU_CreateVert2CellList(pRegion)
      CALL RFLU_BuildVert2CellList(pRegion)
    
      IF ( pMixtInput%spaceOrder > 1 ) THEN 
        CALL RFLU_SetInfoC2CStencilWrapper(pRegion,pMixtInput%spaceOrder-1) 
        CALL RFLU_CreateC2CStencilWrapper(pRegion)
        CALL RFLU_BuildC2CStencilWrapper(pRegion)
        CALL RFLU_BuildListCC2CStencil(pRegion) 
      END IF ! pMixtInput%spaceOrder      
    
      IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN
        CALL RFLU_SetInfoF2CStencilWrapper(pRegion,pMixtInput%spaceOrder-1)    
        CALL RFLU_CreateF2CStencilWrapper(pRegion)      
        CALL RFLU_BuildF2CStencilWrapper(pRegion)   
        CALL RFLU_BuildListCF2CStencil(pRegion) 
      END IF ! pMixtInput%flowModel
      
      DO iPatch = 1,pRegion%grid%nPatches 
        pPatch => pRegion%patches(iPatch)

        IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
          CALL RFLU_SetInfoBF2CStencilWrapper(pRegion,pPatch, &
                                              pPatch%spaceOrder)    
          CALL RFLU_CreateBF2CStencilWrapper(pRegion,pPatch)      
          CALL RFLU_BuildBF2CStencilWrapper(pRegion,pPatch) 
        END IF ! RFLU_DecideNeedBGradFace
      END DO ! iPatch
    
      CALL RFLU_DestroyVert2CellList(pRegion)
    END DO ! iReg
  END IF ! RFLU_DecideNeedStencils

! ******************************************************************************
! Compute weights. NOTE must be done after stencils are built.
! ******************************************************************************

! ==============================================================================
! Weights for cell and face gradients
! ==============================================================================

  IF ( RFLU_DecideNeedWeights(pRegion) .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal 
      pRegion => pLevel%regions(iReg)
      pMixtInput => pRegion%mixtInput   

      IF ( pMixtInput%spaceOrder > 1 ) THEN 
        CALL RFLU_CreateWtsC2CWrapper(pRegion,pMixtInput%spaceOrder-1)
        CALL RFLU_ComputeWtsC2CWrapper(pRegion,pMixtInput%spaceOrder-1)
      END IF ! pMixtInput%spaceOrder      

      IF ( pMixtInput%flowModel == FLOW_NAVST ) THEN
        CALL RFLU_CreateWtsF2CWrapper(pRegion,pMixtInput%spaceOrder-1)    
        CALL RFLU_ComputeWtsF2CWrapper(pRegion,pMixtInput%spaceOrder-1)
      END IF ! pMixtInput%flowModel
      
      DO iPatch = 1,pRegion%grid%nPatches
        pPatch => pRegion%patches(iPatch)
      
        IF ( RFLU_DecideNeedBGradFace(pRegion,pPatch) .EQV. .TRUE. ) THEN
          CALL RFLU_CreateWtsBF2CWrapper(pRegion,pPatch,pPatch%spaceOrder)    
          CALL RFLU_ComputeWtsBF2CWrapper(pRegion,pPatch,pPatch%spaceOrder)    
        END IF ! RFLU_DecideNeedBGradFace
      END DO ! iPatch
    END DO ! iReg   
  END IF ! RFLU_DecideNeedWeights

! ==============================================================================
! Weights for optimal LES approach. NOTE that at present, this will work only
! for non-moving grids.
! ==============================================================================

  IF ( RFLU_DecideNeedWeights(pRegion) .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal   
      pRegion => pLevel%regions(iReg)
      pMixtInput => pRegion%mixtInput 

      IF ( pMixtInput%spaceDiscr == DISCR_OPT_LES ) THEN            
        CALL RFLU_FindPrototypeFacesOLES(pRegion)
        CALL RFLU_BuildStencilsOLES(pRegion)
        CALL RFLU_ComputeGeometricTermsOLES(pRegion)
        CALL RFLU_BuildSymmetryMapsOLES(pRegion)
      
!        CALL RFLU_ComputeIntegrals1245OLES(pRegion)
!        CALL RFLU_EnforceSymmetryOLES(pRegion)
!        CALL RFLU_WriteIntegrals1245OLES(pRegion)
      
        CALL RFLU_ReadIntegrals1245OLES(pRegion)
        CALL RFLU_EnforceSymmetryOLES(pRegion) 
!        CALL RFLU_WriteIntegrals1245OLES(pRegion)           
      END IF ! pMixtInput 
    END DO ! iReg 
  END IF ! RFLU_DecideNeedWeights

#ifdef PLAG
! ******************************************************************************
! Build stencils and weights.
! Needed for particle volume fraction computation even for spaceOrder < 2 
! ******************************************************************************

  IF ( RFLU_DecideNeedStencils(pRegion) .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal 
      pRegion => pLevel%regions(iReg)
      pMixtInput => pRegion%mixtInput 

      CALL RFLU_CreateVert2CellList(pRegion)
      CALL RFLU_BuildVert2CellList(pRegion)
    
      IF ( (pRegion%global%plagUsed .EQV. .TRUE.) .AND. &
           (pMixtInput%spaceOrder < 2) ) THEN 
        CALL RFLU_SetInfoC2CStencilWrapper(pRegion,1) ! Forcing to create 
        CALL RFLU_CreateC2CStencilWrapper(pRegion)    ! stencil irrespective of
        CALL RFLU_BuildC2CStencilWrapper(pRegion)     ! spaceOrder
        CALL RFLU_BuildListCC2CStencil(pRegion) 
      END IF ! pMixtInput%spaceOrder      
    
      CALL RFLU_DestroyVert2CellList(pRegion)
    END DO ! iReg
  END IF ! RFLU_DecideNeedStencils

  IF ( RFLU_DecideNeedWeights(pRegion) .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal 
      pRegion => pLevel%regions(iReg)
      pMixtInput => pRegion%mixtInput   

      IF ( (pRegion%global%plagUsed .EQV. .TRUE.) .AND. &
           (pMixtInput%spaceOrder < 2) ) THEN 
        CALL RFLU_CreateWtsC2CWrapper(pRegion,1)  ! Forcing to create Wts
        CALL RFLU_ComputeWtsC2CWrapper(pRegion,1) ! irrespective of spaceOrder
      END IF ! pMixtInput%spaceOrder      
    END DO ! iReg   
  END IF ! RFLU_DecideNeedWeights

! ******************************************************************************
! Compute in-cell test tolerance
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)

    IF ( pRegion%global%plagUsed .EQV. .TRUE. ) THEN 
      CALL RFLU_ICT_ComputeTolerance(pRegion)
    END IF ! pRegion%global%plagUsed
  END DO ! iReg
#endif

! ******************************************************************************
! Allocate memory for absorbing boundary condition
! ******************************************************************************

  IF ( global%abcFlag .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)
      pMixtInput => pRegion%mixtInput

      IF ( global%abcKind == 0 ) THEN
        CALL RFLU_ABC_CreateSigma(pRegion)
        CALL RFLU_ABC_InitSigma(pRegion)
        CALL RFLU_ABC_SetSigma(pRegion)
        CALL RFLU_ABC_SetRefSoln(pRegion)
      END IF ! global%abcKind
    END DO ! iReg
  END IF ! global%abcFlag

! ******************************************************************************
! Modify geometry for axisymmetric flow
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)

    IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
      CALL RFLU_AXI_ScaleGeometry(pRegion)
    END IF ! pRegion%mixtInput%axiFlag
  END DO ! iReg

! ******************************************************************************
! Read solution data files
! ******************************************************************************
  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)     
    CALL RFLU_ReadFlowWrapper(pRegion)

    IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
      CALL RFLU_ReadAuxVarsWrapper(pRegion)
    END IF ! global%solverType

    CALL RFLU_ReadPatchCoeffsWrapper(pRegion)   
 
#ifdef PLAG
    IF ( global%plagUsed .EQV. .TRUE. ) THEN
      IF ( PLAG_DecideHaveSurfStats(pRegion) .EQV. .TRUE. ) THEN 
        CALL PLAG_ReadSurfStatsWrapper(pRegion)
      END IF ! PLAG_DecideHaveSurfStats
    END IF ! plagUsed
#endif          
  END DO ! iReg

  IF ( global%mvFrameFlag .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)
      CALL RFLU_MVF_InitPatchVelAccel(pRegion)
      CALL RFLU_MVF_ReadPatchVelAccel(pRegion)
    END DO ! iReg
  END IF ! global%mvFrameFlag

! ******************************************************************************
! Setup Ghost Fluid Method (GFM) related stuff
! ******************************************************************************

  IF ( global%gfmFlag .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)     

!     Renumbering cells to access partitioned cell from serial cell
      IF ( global%nRegionsLocal > 1 ) THEN
        CALL RFLU_RNMB_CreatePC2SCMap(pRegion)    
        CALL RFLU_RNMB_CreatePV2SVMap(pRegion)
        CALL RFLU_RNMB_CreatePBF2SBFMap(pRegion)
      
        CALL RFLU_RNMB_ReadPxx2SxxMaps(pRegion)
        CALL RFLU_RNMB_BuildSC2PCMap(pRegion)
            
        CALL RFLU_RNMB_DestroyPV2SVMap(pRegion)
        CALL RFLU_RNMB_DestroyPBF2SBFMap(pRegion)
      END IF ! global%nRegionsLocal > 1

!     Setting up levelset function
      CALL RFLU_GFM_CreateLevelSet(pRegion)
      CALL RFLU_GFM_InitLevelSet(pRegion)
      CALL RFLU_GFM_SetLevelSet(pRegion)

!      CALL RFLU_GFM_SetGhostFluid(pRegion)

! TEMPORARY: Manoj: Testing LevelSet
      CALL RFLU_WriteFlowWrapper(pRegion)
      WRITE(*,*) "Calling MPI_Barrier from RFLU_TimeStepping"
      CALL MPI_Barrier(global%mpiComm,errorFlag)
      WRITE(*,*) "Done Calling MPI_Barrier from RFLU_TimeStepping"
! END TEMPORARY
    END DO ! iReg
! TEMPORARY: Manoj: Testing LevelSet
      WRITE(*,*) "Stopping here..."
      STOP
! END TEMPORARY
  END IF ! global%gfmFlag

! ******************************************************************************
! Initialize variables neeeded for SOLV_IMPLICIT_HM
! ******************************************************************************

! ==============================================================================
! Compute cell number offsets in each region
! ==============================================================================

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)

      CALL RFLU_GID_CreatenCellsOffset(pRegion)
    END DO ! iReg

    CALL RFLU_GID_ComputenCellsOffset(pLevel%regions)
  END IF ! global%solverType

! ==============================================================================
! Create and initialize global cell numbers of virtual cells.
! ==============================================================================

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)
      pGrid   => pRegion%grid

      IF ( pGrid%nBorders > 0 ) THEN
        CALL RFLU_GID_CreateGlobalIds(pRegion)
      END IF ! pRegion%grid%nBorders
    END DO ! iReg 

    CALL RFLU_GID_ComputeGlobalIds(pLevel%regions)
  END IF ! global%solverType

! ******************************************************************************
! Initialize Hypre solver neeeded for SOLV_IMPLICIT_HM
! ******************************************************************************

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    CALL RFLU_HYPRE_CreateObjects(pLevel%regions)
  END IF ! global%solverType

#ifdef GENX
! ******************************************************************************
! Get interface data
! ******************************************************************************
  
  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)     
    CALL RFLU_GENX_GetDataInterf(pRegion)
  END DO ! iReg  
#endif

! ******************************************************************************
! Read boundary condition data (only relevant if distribution specified)
! Also Read boundary variables if required, i.e., if NSCBC==1 
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)

    CALL RFLU_BXV_ReadVarsWrapper(pRegion)

    IF ( RFLU_DecideReadWriteBcDataFile(pRegion) .EQV. .TRUE. ) THEN 
      CALL RFLU_ReadBcDataFile(pRegion)
    END IF ! RFLU_DecideReadWriteBcDataFile
  END DO ! iReg
  
! ******************************************************************************
! Read grid speeds. NOTE must be done after face list is constructed, because 
! need to know number of faces to be able to read grid speeds. 
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    
    IF ( RFLU_DecideNeedGridSpeeds(pRegion) .EQV. .TRUE. ) THEN 
      CALL RFLU_ReadGridSpeedsWrapper(pRegion)    
    END IF ! RFLU_DecideNeedGridSpeeds
  END DO ! iReg

! ******************************************************************************
! Initialize grid speed scaling factor. NOTE needs to be done for steady flows
! to work properly because grid speed is used together with scaling routines.
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    
    CALL RFLU_InitGridSpeedScaleFactor(pRegion)
  END DO ! iReg

#ifdef STATS
! ******************************************************************************
! Statistics initialization
! ******************************************************************************

  CALL InitStatistics(pLevel%regions)
#endif

! ******************************************************************************
! Create buffers and tags for communication
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    pGrid   => pRegion%grid   
    
    IF ( pGrid%nBorders > 0 ) THEN 
      CALL RFLU_MPI_CreateBuffersWrapper(pRegion)
      CALL RFLU_MPI_SetTagsWrapper(pRegion)
      
#ifdef PLAG
      IF ( pRegion%global%plagUsed .EQV. .TRUE. ) THEN 
        CALL RFLU_MPI_CreateBufferIPclSend(pRegion)
      END IF ! pRegion%global%plagUsed
#endif      
    END IF ! pRegion%grid%nBorders
  END DO ! iReg  

! *****************************************************************************
! Initialize virtual-cell data via communication (done for GENX runs after 
! remeshing because remeshing tool does not initialize ghosts)
!
! Routines to communicate variables at region boundaries needs information
! of which variable to communicate. Not all variables needed to be 
! communicated in this predictor-corrector solver unlike in dissipative 
! solver. Information about which variable to communicate is passed through
! variable 'iVar' which can have values from 1 to 9 which has following
! meaning,
!      1          density from cv
!      2,3,4      x,y,z velocities from cv
!      5          pressure from cv
!      6          delp
!      7          density from cvOld
!      8          pressure from cvold
!      9          gradients
! *****************************************************************************

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    DO iVar = 1,8
      DO iReg = 1,global%nRegionsLocal
        pRegion => pLevel%regions(iReg)

        CALL RFLU_MPI_HM_ISendWrapper(pRegion,iVar)
      END DO ! iReg 

      CALL RFLU_MPI_HM_CopyWrapper(pLevel%regions,iVar)
    
      DO iReg = 1,global%nRegionsLocal
        pRegion => pLevel%regions(iReg)

        CALL RFLU_MPI_HM_RecvWrapper(pRegion,iVar)
      END DO ! iReg 

      DO iReg = 1,global%nRegionsLocal
        pRegion => pLevel%regions(iReg)

        CALL RFLU_MPI_ClearRequestWrapper(pRegion)
      END DO ! iReg 
    END DO ! iVar
      
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)

      CALL RFLU_RELP_TransformWrapper(pRegion)
    END DO ! iReg 

  ELSE

    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)

      CALL RFLU_MPI_ISendWrapper(pRegion)
    END DO ! iReg 

    CALL RFLU_MPI_CopyWrapper(pLevel%regions)
  
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)
        
      CALL RFLU_MPI_RecvWrapper(pRegion)
      CALL RFLU_RELP_TransformWrapper(pRegion)

    END DO ! iReg 

    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)

      CALL RFLU_MPI_ClearRequestWrapper(pRegion)
    END DO ! iReg 
  END IF ! global%solverType

! ******************************************************************************
! Initialize remaining variables neeeded for SOLV_IMPLICIT_HM (done after 
! virtual cells are initialized)
! ******************************************************************************

  IF ( global%solverType == SOLV_IMPLICIT_HM ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)

! ==============================================================================
! Face-normal velocity
! ==============================================================================

      CALL RFLU_HM_PredictFaceNormalVelocity(pRegion,pRegion%mixt%cv, &
                                             pRegion%mixt%vfMixt)

! ==============================================================================
! Velocity and pressure gradients
! ==============================================================================

      IF ( pRegion%mixtInput%flowModel == FLOW_NAVST ) THEN
        IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
          CALL RFLU_AXI_ComputeGradCellsGGVector(pRegion,CV_MIXT_XVEL, &
                                                 CV_MIXT_ZVEL,GRC_MIXT_XVEL, &
                                                 GRC_MIXT_ZVEL,pRegion%mixt%cv, &
                                                 pRegion%mixt%gradCell)
        ELSE
          CALL RFLU_ComputeGradCellsGGVector(pRegion,CV_MIXT_XVEL,CV_MIXT_ZVEL, &
                                             GRC_MIXT_XVEL,GRC_MIXT_ZVEL, &
                                             pRegion%mixt%cv, &
                                             pRegion%mixt%gradCell)
        END IF ! pRegion%mixtInput%axiFlag
      END IF ! pRegion%mixtInput%flowModel

      flag = 1
      IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
        CALL RFLU_AXI_ComputeGradCellsGGScalar(pRegion,CV_MIXT_PRES, &
                                               CV_MIXT_PRES,GRC_MIXT_PRES, &
                                               GRC_MIXT_PRES,pRegion%mixt%cv, &
                                               pRegion%mixt%gradCell,flag)
      ELSE
        CALL RFLU_ComputeGradCellsGGScalar(pRegion,CV_MIXT_PRES,CV_MIXT_PRES, &
                                           GRC_MIXT_PRES,GRC_MIXT_PRES, &
                                           pRegion%mixt%cv, &
                                           pRegion%mixt%gradCell,flag)
      END IF ! pRegion%mixtInput%axiFlag

      flag = 2
      IF ( pRegion%mixtInput%axiFlag .EQV. .TRUE. ) THEN
        CALL RFLU_AXI_ComputeGradCellsGGScalar(pRegion,CV_MIXT_PRES, &
                                               CV_MIXT_PRES,GRC_MIXT_PRES, &
                                               GRC_MIXT_PRES, &
                                               pRegion%mixt%cvOld, &
                                               pRegion%mixt%gradCellOld, &
                                               flag)
      ELSE
        CALL RFLU_ComputeGradCellsGGScalar(pRegion,CV_MIXT_PRES,CV_MIXT_PRES, &
                                           GRC_MIXT_PRES,GRC_MIXT_PRES, &
                                           pRegion%mixt%cvOld, &
                                           pRegion%mixt%gradCellOld,flag)
      END IF ! pRegion%mixtInput%axiFlag
    END DO ! iReg
  END IF ! global%solverType

! *****************************************************************************
! Initialize virtual-cell gradient data via communication
! *****************************************************************************

  IF ( (global%solverType == SOLV_IMPLICIT_HM) .AND. &
       (pLevel%regions(1)%mixtInput%flowModel == FLOW_NAVST) ) THEN
    iVar = 9

    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)

      CALL RFLU_MPI_HM_ISendWrapper(pRegion,iVar)
    END DO ! iReg 

    CALL RFLU_MPI_HM_CopyWrapper(pLevel%regions,iVar)
    
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)

      CALL RFLU_MPI_HM_RecvWrapper(pRegion,iVar)
    END DO ! iReg 

    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)

      CALL RFLU_MPI_ClearRequestWrapper(pRegion)
    END DO ! iReg 
      
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)

      CALL RFLU_RELP_TransformWrapper(pRegion)
    END DO ! iReg 
  END IF ! global%solverType

! ******************************************************************************
! Initialize dependent variables 
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)

#ifdef PLAG
    IF ( global%plagUsed .AND. (pRegion%plag%nPcls > 0) ) THEN
      CALL PLAG_RFLU_ComputeVolFrac(pRegion)
    END IF ! global%plagUsed
#endif

! Rahul BUG FIX -  Dependent variables are updated after volume fraction is
!                   communicated to virtual cells  
!    CALL RFLU_SetVarsWrapper(pRegion,1,pRegion%grid%nCellsTot)   

! TEMPORARY : move it to proper place ....
!    IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
!      CALL RFLU_BXV_SetDependentVars(pRegion)
!    END IF !
! END TEMPORARY
! Rahul - End BUG FIX

  END DO ! iReg

 !BBR - begin
 !CALL MPI_Barrier(global%mpiComm,errorFlag)
 !WRITE(777,*) "Line 1576...so far so good"
 !BBR - end

! *****************************************************************************
! Initialize virtual-cell particle volume fraction for computation of gradient 
! *****************************************************************************

! TEMPORARY: Manoj: 2012-05-29: checking effect of not communicating vFracE
!IF (1==2) THEN
#ifdef PLAG
  IF ( global%plagUsed .EQV. .TRUE. ) THEN
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)

      CALL RFLU_MPI_PLAG_ISendWrapper(pRegion)
    END DO ! iReg 

    CALL RFLU_MPI_PLAG_CopyWrapper(pLevel%regions)
    
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)

      CALL RFLU_MPI_PLAG_RecvWrapper(pRegion)
    END DO ! iReg 

    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)

      CALL RFLU_MPI_ClearRequestWrapper(pRegion)
    END DO ! iReg 
  END IF ! global%plagUsed
#endif 

!MOVING PICL HERE TO AVOID SEG FAULT OF ROCFLU STORED VF
#ifdef PICL

!Num of cells with particles
      !npart = 3600*Pval ! 40000 = 11.11  
        !10800 = 3
      !per_cell = npart / 3600 !hardcoded number of cells
 !Bounds of particle curtain       
  !    xMinCell = 0.034
  !    xMaxCell = 0.04
  !    yMinCell = 0.0
  !    yMaxCell = 0.006
      !zMinCell = z   - 0.5_RFREAL*pGrid%dz
      !zMaxCell = z   + 0.5_RFREAL*pGrid%dz

!      numPclCells = 0  
!      DO icl = 1, pGrid%nCells
!        x = pGrid%cofg(XCOORD,icl)
!        y = pGrid%cofg(YCOORD,icl)
        !z = pGrid%cofg(ZCOORD,icl)

!        if ((x .ge. xMinCell) .and. (x .le. xMaxCell)) then
!        if ((y .ge. yMinCell) .and. (y .le. yMaxCell)) then
        !cell within particle curtain        
!                numPclCells = numPclCells + 1
!        endif
!        endif


      !enddo     
 
      PPC = 1 !Hardcoded 
      !npart = numPclCells * PPC !Hardcodded for 40,000  
  
      !npart   = 0    ! particles/rank to distribute
       !IF ( global%myProcid == MASTERPROC) npart = 5000
       !IF ( global%myProcid == MASTERPROC) write(*,*) "BADPQ NP =", npart 
      !dp_min  = 115.0E-6   ! particle diameter min
      !dp_max  = 200E-7  ! particle diameter max
      rhop    = 200.0 !8000.0 !2500 ! particle density
      filter = 0.0002*sqrt(3.0)+1.0E-10 !0.0002  
!      rdum    = ran2(-1-ppiclf_nid) ! initialize random number generator
!Init particle props such as loc, vel,dens...
!Place Particles
 !Loop through Cells
  !Loop through particles per cell
  !Place then confirm
  !Apply properties including spload

!Restart from file
pf_fluidInit = .true.
npart = 0
        
!default values
!zpf_factor = 0.0006/0.0002!1.0!0.0012/0.0002!0.0002
zpf_factor = 1.0 !global%ppiclFInitZpfFactor
!xpf_factor = 0.003
xpf_factor = 0.0

! Sam - hard coded initialization for testing purposes
!npart = 2
!do i=1,npart
!  y(PPICLF_JX,i) = 0.52_RFREAL            ! x position
!  y(PPICLF_JY,i) = 0.4_RFREAL*i          ! y position
!  y(PPICLF_JZ,i) = 0.52_RFREAL
!  ! zero velocity in all directions
!  y(PPICLF_JVX,i) = 0.0
!  y(PPICLF_JVY,i) = 0.0
!  y(PPICLF_JVZ,i) = 0.0
!
!  ! initially zero out all properties
!  do ii=1,PPICLF_LRP
!    rprop(i7, i) = 0.0
!  end do
!
!  ! now set properties that are not interpolated from Rocflu onto the particles
!  rprop(PPICLF_R_JRHOP,i) = 2710.0 ! particle density
!  dp_min = 1.0*10.0**(-3)
!  rprop(PPICLF_R_JDP,i)   = dp_min ! particle diameter
!  rprop(PPICLF_R_JVOLP,i) = (4.0_RFREAL/3.0_RFREAL)*global%pi*&
!                            (0.5_RFREAL*dp_min)**3 ! particle volume
!
!  rprop(PPICLF_R_JSPL,i) = 1.0_RFREAL
!
!end do

! Sam - hard code initialization reading .dat file with points set all other
! properties manually except location
IF ( global%verbLevel > VERBOSE_NONE ) THEN 
  WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Reading points.dat file...'
END IF ! global%verbLevel

iFileName = 'points.dat'
OPEN(iFile,FILE=iFileName,FORM="FORMATTED",STATUS="OLD",IOSTAT=errorFlag)
global%error = errorFlag   
IF ( global%error /= ERR_NONE ) THEN 
  CALL ErrorStop(global,ERR_FILE_OPEN,__LINE__,iFileName)
END IF

READ(iFile,'(I16)') npart ! global number of particles

if (npart .gt. PPICLF_LPART*global%nProcs) then
  CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__,'PPICLF:too &
    many particles to initialize')
end if

npart_local = npart/global%nProcs+1
i = 1
i_global = 1
i_global_min = npart_local*global%myProcid
i_global_max = npart_local*(global%myProcid+1)

do i_global=1,npart
  ! ASSUME: PPICLF_JX=1, PPICLF_JY=2, ...
  READ(iFile,'(3E23.16)') (y(ii,i),ii=PPICLF_JX,PPICLF_JZ)

  ! if in range for this processor set all the other properties and increment i
  if ((i_global .gt. i_global_min) .and. (i_global .le. i_global_max)) then
    y(PPICLF_JVX,i) = 0.0
    y(PPICLF_JVY,i) = 0.0
    y(PPICLF_JVZ,i) = 0.0

    ! initially zero out all properties
    do ii=1,PPICLF_LRP
      rprop(i7, i) = 0.0
    end do
  
    ! now set properties that are not interpolated from Rocflu onto the particles
    rprop(PPICLF_R_JRHOP,i) = 200.0 ! particle density
    dp_min = 9.5*10**(-5) ! particle diameter
    rprop(PPICLF_R_JVOLP,i) = (4.0_RFREAL/3.0_RFREAL)*global%pi*&
                              (0.5_RFREAL*dp_min)**3 ! particle volume
  
    rprop(PPICLF_R_JSPL,i) = 1.0_RFREAL

    i = i + 1
  end if
end do

! Close points.dat file
CLOSE(iFile, IOSTAT=errorFlag)
global%error = errorFlag   
IF ( global%error /= ERR_NONE ) THEN 
  CALL ErrorStop(global,ERR_FILE_CLOSE,__LINE__,iFileName)
END IF ! global%error  
 
! Sam - Particle initialization from input file

! The restart routine must be called from each process, so it must be outside
! the switch
!if (global%ppiclFInitType .eq. 0) then
!  call ppiclf_io_ReadParticleVTU('par_restart.vtu')
!  npart = 0
!end if

!PPICLF_SETTINGS_R_HARD_CUTOFF = global%ppiclFInitRHardCutoff
!call ppiclf_solve_InitUserSetting('breakup_model', &
!     global%ppiclFInitBreakupModel, 0.0)
!
!call ppiclf_solve_InitUserSetting('CF_TAB', &
!    0, global%ppiclFInitCFTAB)
!
!call ppiclf_solve_InitUserSetting('n_child', &
!    global%ppiclFInitNChild, 0.0) 

! inintialize with no particles
!call ppiclf_solve_InitParticle(-1,3,0,0,y,rprop) 

!select case (global%ppiclFInitType)
!  ! Particle initialization for a line in the y direction
!  ! Still a line, but now we vary the properties linearly
!  case (2)
!    y0 = global%ppiclFInitYMin!0.0
!    ymax = global%ppiclFInitYMax !0.04 ! max y value for a particle
!    npart = global%ppiclFInitNPart*global%ppiclFInitND*&
!            global%ppiclFInitNRho*global%ppiclFInitNSigma*&
!            global%ppiclFInitNMu*global%ppiclFInitNDuplicates
!
!    if (npart .gt. PPICLF_LPART*global%nProcs) then
!        CALL ErrorStop(global,ERR_ILLEGAL_VALUE,__LINE__,'PPICLF:too &
!          many particles to initialize')
!    end if
!
!    npart_local = npart/global%nProcs+1
!
!    if (global%ppiclFInitNPart .gt. 1) then
!        dy = (ymax - y0) / (global%ppiclFInitNPart - 1)
!    else
!        dy = 0.0
!    end if
!
!    if (pRegion%mixtInput%axiFlag) then
!        da_axi = global%pi*(ymax**2 - y0**2) / global%ppiclFInitNPart
!        y_axi0 = 0.0
!    end if
!
!    if (global%ppiclFInitND .gt. 1) then
!        ddiam = (global%ppiclFInitDMax - global%ppiclFInitD) / &
!                (global%ppiclFInitND - 1)
!    else
!        ddiam = 0.0
!    end if
!
!    if (global%ppiclFInitNRho .gt. 1) then
!        dRho = (global%ppiclFInitRhoMax - global%ppiclFInitRho) /&
!             (global%ppiclFInitNRho - 1)
!    else
!        dRho = 0.0
!    end if
!
!    if (global%ppiclFInitNSigma .gt. 1) then
!        dSigma = (global%ppiclFInitSigmaMax-global%ppiclFInitSigma)/&
!                 (global%ppiclFInitNSigma - 1)
!    else
!        dSigma = 0.0
!    end if
!
!    if (global%ppiclFInitNMu .gt. 1) then
!        dMu = (global%ppiclFInitMuMax-global%ppiclFInitMu)/&
!                 (global%ppiclFInitNMu - 1)
!    else
!        dMu = 0.0
!    end if
!
!    x0 = global%ppiclFInitX !-0.09
!    z0 = global%ppiclFInitZ !0.0005
!
!    rhop = global%ppiclFInitRho !1000.0 ! density
!    sigma = global%ppiclFInitSigma !71.99*10.0**(-3.0)
!    mutab = global%ppiclFInitMu !8.891*10.0**(-4.0)
!    
!    i = 1
!    i_global = 1
!    i_global_min = npart_local*global%myProcid
!    i_global_max = npart_local*(global%myProcid+1)
!    
!    do i5=1,global%ppiclFInitNRho
!    rhop = global%ppiclFInitRho + dRho*(i5-1)
!
!    do i4=1,global%ppiclFInitNSigma
!    sigma = global%ppiclFInitSigma + dSigma*(i4-1)
!    
!    do i3=1,global%ppiclFInitNMu
!    mutab = global%ppiclFInitMu + dMu*(i3-1)
!    
!    do ii=1,global%ppiclFInitND
!      diam = global%ppiclFInitD + ddiam*(ii-1) !1.0*10.0**(-3.0) ! diameter
!      vol = 4.0/3.0*global%pi*(diam/2.0)**3.0
!      
!      dcut = global%ppiclFInitDCutoff
!
!      if (pRegion%mixtInput%axiFlag) then
!        y_axi0 = 0.0
!      end if
!
!      do jj=1,global%ppiclFInitNPart
!        if (pRegion%mixtInput%axiFlag) then
!          y_axi1 = SQRT(da_axi/global%pi + y_axi0**2)
!        end if
!        do i6=1,global%ppiclFInitNDuplicates
!          if ((i_global .gt. i_global_min) .and. &
!              (i_global .le. i_global_max)) then
!            ! initial position
!            y(PPICLF_JX, i) = x0
!            if (pRegion%mixtInput%axiFlag) then
!              y(PPICLF_JY,i) = SQRT(0.5*da_axi/global%pi + y_axi0**2)
!            else
!              y(PPICLF_JY, i) = y0 + dy*(jj-1)
!            end if
!            y(PPICLF_JZ, i) = z0*zpf_factor
!            ! initial velocity
!            y(PPICLF_JVX, i) = global%ppiclFInitVX
!            y(PPICLF_JVY, i) = global%ppiclFInitVY
!            y(PPICLF_JVZ, i) = global%ppiclFInitVZ
!
!            do i7=1,PPICLF_LRP
!              rprop(i7, i) = 0.0
!            end do
!
!            rprop(PPICLF_R_JRHOP, i) = rhop
!            rprop(PPICLF_R_JDP, i) = diam
!            rprop(PPICLF_R_JVOLP, i) = vol
!            rprop(PPICLF_R_JPHIP, i) = vol ! WRONG?? irrelevant for 1-way coupling anyway
!
!            select case(global%ppiclFInitBreakupModel)
!              case (1)
!                ! TAB
!                rprop(PPICLF_R_JYTAB, i) = 0.0
!                rprop(PPICLF_R_JYDOTTAB, i) = 0.0
!                rprop(PPICLF_R_JR0TAB, i) = 0.5*diam
!                rprop(PPICLF_R_JSIGMATAB, i) = sigma !71.99*10**(-3.0)
!                rprop(PPICLF_R_JMUPTAB, i) =  mutab !8.891*10.0**(-4.0)
!                rprop(PPICLF_R_JFTAB,i) = 0.0
!              case (2)
!                ! Ashgriz internal flow
!                rprop(PPICLF_R_JD0ASH,i) = diam
!                rprop(PPICLF_R_JtIASH,i) = 0.0
!                rprop(PPICLF_R_Jt0ASH,i) = 0.0
!                rprop(PPICLF_R_JRDOTASH,i) = 0.0
!                rprop(PPICLF_R_JU0ASH,i) = 0.0
!                rprop(PPICLF_R_JSIGMAASH,i) = sigma
!              case (3) 
!                ! KHRT model
!                rprop(PPICLF_R_JR0KHRT,i) = 0.5*diam
!                rprop(PPICLF_R_JRPREVKHRT,i) = 0.5*diam
!                rprop(PPICLF_R_JtIKHRT,i) = 0.0
!                rprop(PPICLF_R_JSIGMAKHRT,i) = sigma
!                rprop(PPICLF_R_JMUPKHRT,i) = mutab
!
!                y(PPICLF_JRKHRT,i) = 0.5*diam
!            end select
!
!            ! super particle loading
!            rprop(PPICLF_R_JSPL,i) = 1.0
!            ! other stuff
!            rprop(PPICLF_R_JDPi,i) = (rprop(PPICLF_R_JDP,i)/2.0) &
!                *(rprop(PPICLF_R_JSPL,i))**(1.0/3.0) 
!            rprop(PPICLF_R_JDPe,i) = diam
!            rprop(PPICLF_R_JSPT,i) = 1!Maybe need to init
!
!            rprop(PPICLF_R_JPARENT,i) = i_global
!            rprop(PPICLF_R_JRCUTOFF,i) = dcut/2.0
!            
!            rprop(PPICLF_R_JVXOLD,i) = 0.0
!            rprop(PPICLF_R_JVYOLD,i) = 0.0
!
!            ! generate unique ID based on Hopcroft and Ullman pairing function
!            ! assume tindex=0 which will guarantee uniqueness from breakup
!            ! model generated IDs
!            ! do this the same way its done in ppiclf: tindex, pid, i, j
!            uid = 0.5*(-1+global%myProcid-2)*(-1+global%myProcid-1) - 1
!            uid = 0.5*(uid+i-2)*(uid+i-1)+uid
!            ! assume j=1 since we're only adding one particle each i
!            !uid = 0.5*(uid+1-2)*(uid+1-1)+uid
!            uid = 0.5*(uid-1)*(uid)+uid ! simplified from previous line
!            rprop(PPICLF_R_JUID,i) = uid
!
!            i = i + 1
!          end if
!          i_global = i_global + 1
!        end do
!        if (pRegion%mixtInput%axiFlag) then
!          y_axi0 = y_axi1
!        end if
!      end do
!    end do
!    end do
!    end do
!    end do
!  case default
!      npart = 0
!end select

!Init particles
!write(*,*) 'PID',global%myProcid 
    call ppiclf_solve_InitParticle(-1,3,0,i-1,y,rprop) 
    !call ppiclf_user_trimPuck    
! User sets up any filter parameters (i.e., Box or Gaussian):
!    call ppiclf_solve_InitGaussianFilter(3.0*dp_min,0,0)!3.0*dp_min,1E-3,0
!    call ppiclf_solve_InitPeriodicY(0.0,5E-3)
     !call ppiclf_solve_InitPeriodicZ(0.0,1E-4)
IF ( global%myProcid == MASTERPROC) then
      write(*,*) "BOXF = ", filter 
      write(*,*) "PHIPF = ",  global%pi/6.0*(filter)**3
      write(*,*) "XLOC = ", y(PPICLF_JX ,1) 
      write(*,*) "YLOC = ",   y(PPICLF_JY ,1) 
      write(*,*) "ZLOC = ",  y(PPICLF_JZ ,1)   

end if  
     !call ppiclf_solve_InitBoxFilter(dp_min*2.0d0,0) !10*dp_min,0

!SINGLE TEST
     !call ppiclf_solve_InitBoxFilter(0.0,0,1) !TEST
     !call ppiclf_solve_InitBoxFilter(0.1d0,0,0) !TEST
     call ppiclf_solve_InitGaussianFilter(dp_min*3.0d0,0.01d0,0)
!Single TEST   
     !tester = tester * 0.5          
     !call ppiclf_solve_InitGaussianFilter(tester,0.01,0)
     

!BRAD STARTS HERE
!Taking what was already done and building 1-way coupled first
!NEED TO BUILD overlap mesh 
!Loop through cells
 !Pull face id's
  !grab face center coords 
   !dump into approperate array
!How cooord is grab needs to match how props are grabbed


     ! User sets up overlap mesh:
    nCells = pRegion%grid%nCells
    lx = 2
    ly = 2
    lz = 2
     
    ALLOCATE(xGrid(lx,ly,lz,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    ALLOCATE(yGrid(lx,ly,lz,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:yGrid')
    END IF ! global%error

    ALLOCATE(zGrid(lx,ly,lz,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:zGrid')
    END IF ! global%error
!Loop cells
    DO i = 1, nCells 

vi = pRegion%grid%hex2v(1,i) 
         xGrid(1,1,1,i) = pRegion%grid%xyz(XCOORD,vi) 
         yGrid(1,1,1,i) = pRegion%grid%xyz(YCOORD,vi) 
         zGrid(1,1,1,i) = zpf_factor*pRegion%grid%xyz(ZCOORD,vi) 
vi = pRegion%grid%hex2v(4,i) 
         xGrid(2,1,1,i) = pRegion%grid%xyz(XCOORD,vi)  
         yGrid(2,1,1,i) = pRegion%grid%xyz(YCOORD,vi) 
         zGrid(2,1,1,i) = zpf_factor*pRegion%grid%xyz(ZCOORD,vi) 
vi = pRegion%grid%hex2v(5,i) 
         xGrid(1,2,1,i) = pRegion%grid%xyz(XCOORD,vi) 
         yGrid(1,2,1,i) = pRegion%grid%xyz(YCOORD,vi) 
         zGrid(1,2,1,i) = zpf_factor*pRegion%grid%xyz(ZCOORD,vi) 
vi = pRegion%grid%hex2v(8,i) 
         xGrid(2,2,1,i) = pRegion%grid%xyz(XCOORD,vi) 
         yGrid(2,2,1,i) = pRegion%grid%xyz(YCOORD,vi) 
         zGrid(2,2,1,i) = zpf_factor*pRegion%grid%xyz(ZCOORD,vi) 
vi = pRegion%grid%hex2v(2,i) 
         xGrid(1,1,2,i) = pRegion%grid%xyz(XCOORD,vi) 
         yGrid(1,1,2,i) = pRegion%grid%xyz(YCOORD,vi)
         zGrid(1,1,2,i) = zpf_factor*pRegion%grid%xyz(ZCOORD,vi)
vi = pRegion%grid%hex2v(3,i) 
         xGrid(2,1,2,i) = pRegion%grid%xyz(XCOORD,vi) 
         yGrid(2,1,2,i) = pRegion%grid%xyz(YCOORD,vi) 
         zGrid(2,1,2,i) = zpf_factor*pRegion%grid%xyz(ZCOORD,vi) 
vi = pRegion%grid%hex2v(6,i) 
         xGrid(1,2,2,i) = pRegion%grid%xyz(XCOORD,vi) 
         yGrid(1,2,2,i) = pRegion%grid%xyz(YCOORD,vi) 
         zGrid(1,2,2,i) = zpf_factor*pRegion%grid%xyz(ZCOORD,vi) 
vi = pRegion%grid%hex2v(7,i) 
         xGrid(2,2,2,i) = pRegion%grid%xyz(XCOORD,vi) 
         yGrid(2,2,2,i) = pRegion%grid%xyz(YCOORD,vi)
         zGrid(2,2,2,i) = zpf_factor*pRegion%grid%xyz(ZCOORD,vi) 
IF ( global%myProcid == MASTERPROC .and. (1 .eq.1 )) then
if ((i.eq.1)) then !.or.(i.eq.2)) then
        vii = 0
write(*,*) "Cell ",i
        do kk = 1,lz
        do jj = 1,ly
        do ii = 1,lx
         vii = vii + 1
         vi = pRegion%grid%hex2v(vii,i) 
         write(*,*) "X ", xGrid(ii,jj,kk,i)      
         write(*,*) "Y ", yGrid(ii,jj,kk,i)
         write(*,*) "Z ", zGrid(ii,jj,kk,i)
!         xGrid(ii,jj,kk,i) = pRegion%grid%xyz(XCOORD,vi) 
!         yGrid(ii,jj,kk,i) = pRegion%grid%xyz(YCOORD,vi) 
!         zGrid(ii,jj,kk,i) = pRegion%grid%xyz(ZCOORD,vi) 
        end do !lx
        end do !ly
        end do !lz
end if
end if
    END DO !nCells
        
 !IF (( global%myProcid == MASTERPROC) .or. npart .gt. 0) write(*,*) "AG,BW1"
 call ppiclf_comm_InitOverlapMesh(nCells,lx,ly,lz,xGrid,yGrid,zGrid)
     !call ppiclf_solve_InitTargetBins('z', 3, 0)
     call ppiclf_solve_InitNeighborBin(4*dp_min)
IF (( global%myProcid == MASTERPROC) .or. npart .gt. 0) write(*,*) "AZ,BW2"

!call ppiclf_solve_InitPeriodicZ(0.0,0.001)!0.002
!call ppiclf_solve_InitPeriodicZ(0.0d0,0.012d0)!0.002
!call ppiclf_solve_InitPeriodicZ(minval(zGrid), maxval(zGrid))!0.002

call ppiclf_solve_InitPeriodicZ(0.0d0, 0.0025d0)
!Init walls
IF (( global%myProcid == MASTERPROC) .or. npart .gt. 0) write(*,*) "AZ,BW3"

! check for existence of wall vtk
INQUIRE(FILE='filein.vtk', EXIST=wall_exists)
if (wall_exists) then
  call ppiclf_io_ReadWallVTK('filein.vtk')
else if (global%myProcid == MASTERPROC) then
  WRITE(*,*) 'Could not find filein.vtk'
end if
 
IF ( global%myProcid == MASTERPROC) write(*,*) "PFINIT: Sending yWall"
    !write(*,*) "Back to Rocflu from Per Z"
      
     !Might need this call first for zero order   

!***Giving ppicl the y-coord for axi-sym collision correction
DO i = 1,pRegion%grid%nCells
       
        
       do lz=1,2
       do ly=1,2
       do lx=1,2
        xGrid(lx,ly,lz,i) = pGrid%cofg(YCOORD,i)
     !   yGrid(lx,ly,lz,i) = pRegion%mixt%cv(CV_MIXT_XMOM,i) &
     !                   /pRegion%mixt%cv(CV_MIXT_DENS,i)
     !   zGrid(lx,ly,lz,i) = pRegion%mixt%cv(CV_MIXT_YMOM,i) &
     !                   /pRegion%mixt%cv(CV_MIXT_DENS,i)
end do
end do
end do 
end do

    call ppiclf_solve_InterpFieldUser(PPICLF_R_JSPT,xGrid)
!!!!!!!!!!!!!!!!!!!!!!!!WARNING!!!!!!!!!!!!!!!!!!!!
!Make sure this call is done specifically when starting the simulation
!for the First Time. As in t = 0.000 else the collision factor will be messed up 
if (pf_fluidInit .eqv. .true. ) call ppiclf_solve_InterpFieldUser(PPICLF_R_JDPi,xGrid)    

!*** axi-sym collision correction

    ! call ppiclf_solve_InterpFieldUser(PPICLF_R_JORF,xGrid)
    ! call ppiclf_solve_InterpFieldUser(PPICLF_R_JOUX,yGrid)
    ! call ppiclf_solve_InterpFieldUser(PPICLF_R_JOUY,zGrid)

!    DEALLOCATE(xGrid,STAT=errorFlag)
!    global%error = errorFlag
!    IF ( global%error /= ERR_NONE ) THEN
!      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
!    END IF ! global%error

    DEALLOCATE(yGrid,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:yGrid')
    END IF ! global%error

    DEALLOCATE(zGrid,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:zGrid')
    END IF ! global%error


!Checking for seg faults on init turning off this off for now
!Usinging this for init vol frac 
if (1.eq.1) then!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ALLOCATE(vfP(2,2,2,nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

ALLOCATE(volp(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

ALLOCATE(SPL(nCells),STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

!    ALLOCATE(pRegion%mixt%piclVF(pGrid%nCellsTot),STAT=errorFlag)
!    global%error = errorFlag
!    IF ( global%error /= ERR_NONE ) THEN
!      CALL ErrorStop(global,ERR_ALLOCATE,__LINE__,'PPICLF:xGrid')
!    END IF ! global%error



!Testing particle count
!if Box
ratio = 1!TESTH0.0001**3 !(global%pi/6.0*((filter)+1.0E-10)**3)
!write(*,*) "PICLVF INIT START"
!Init rocflu vol frac array
do i=1,pGrid%nCellsTot
        pRegion%mixt%piclVF(i) = 0.0        
end do

!DO i = 1,pRegion%grid%nCells
!       do lz=1,2
!       do ly=1,2
!       do lx=1,2
!        xGrid(lx,ly,lz,i) = SPL(i)
!       enddo
!        enddo
!        enddo
!enddo 
!call ppiclf_solve_InterpFieldUser(PPICLF_R_JSPL,xGrid)

!write(*,*) "PICLVF INIT DONE"
!if gauss
!ratio = 1 
IF ( global%myProcid == MASTERPROC) write(*,*) "PFINIT: Calc Init VolP"
total_vol = 0
DO i = 1,pRegion%grid%nCells
       volp(i) = 0 
       tester = 0
       SPL(i) = 1   
       do lz=1,2
       do ly=1,2
       do lx=1,2
       call ppiclf_solve_GetProFldIJKEF(lx, ly, lz, i, PPICLF_P_JPHIP,&
                        vfP(lx,ly,lz,i))
!AXIS SYM
!       tester = tester + (0.125*vfP(lx,ly,lz,i)/pRegion%grid%vol(i)*0.0001)*&
!                                (global%pi/6.0*(0.0001*sqrt(3.0)+1.0E-10)**3)
!       vfp(lx,ly,lz,i) = (0.125*vfP(lx,ly,lz,i) /pRegion%grid%vol(i)*0.0001)* &
!                        (global%pi/6.0*(0.0001*sqrt(3.0)+1.0E-10)**3) 

       !if (vfP(lx,ly,lz,i) .gt. 0.0) write(*,*) "BADPQ V:", vfP(lz,ly,lz,i) 
       tester = tester + (0.125*vfP(lx,ly,lz,i))*pRegion%grid%vol(i)*global%ppiclfInitZpfFactor
       !TESTHvfp(lx,ly,lz,i) = (0.125*vfP(lx,ly,lz,i) /pRegion%grid%vol(i))* &
       !                 (ratio)  


!       if (vfP(lx,ly,lz,i) .gt. 0) write(*,*) "Val VFP,i = ", vfP(lx,ly,lz,i), lx,ly,lz,i
       !write(*,*) "Val VFP,i = ", vfP(lx,ly,lz,i) 
       !vfP(lx,ly,lz,i) = vfP(lx,ly,lz,i) !/ pRegion%grid%vol(i)
       end do 
       end do 
       end do
       total_vol = total_vol + tester         
       volp(i) = (tester/(pRegion%grid%vol(i)))!*(0.0006/0.0002)))!*(0.0002/0.0012)
       !volp(i) = tester/(pRegion%grid%vol(i)*2.0*global%pi/0.0001)
       if (tester .gt. 0 .and. pf_fluidInit .eqv. .true. ) then  
        SPL(i) = (0.56 * pRegion%grid%vol(i)) / tester
        if (pGrid%cofg(YCOORD,i) .gt. 0.0064) SPL(i) = (0.18 * pRegion%grid%vol(i)) / tester
       endif 
       !SPL(i) = (0.6 * 2.400000190734853E-015)/tester  
       !if (volp(i) .gt. 0) volp(i) = 0.6  
!*** VOL FRAC CAP
       if (volp(i) .gt. 0.6) then
        !write(*,*) "OVFP: vf, vol",volp(i),pRegion%grid%vol(i)
        volp(i) = 0.6 
        endif
       pRegion%mixt%piclVF(i) = volp(i) 
       !tester = tester / pRegion%grid%vol(i) 
       if (tester .gt. 0) write(*,*) "VFP,S,#,Y=",&
         volp(i),(pRegion%grid%vol(i)*pRegion%grid%vol(i)*global%ppiclfInitZpfFactor),&
         SPL(i),tester/(global%pi/6.0*(dp_min)**3),pGrid%cofg(YCOORD,i)!pRegion%grid%vol(i)  
!       if (vfP(1,1,1,i) .gt. 1E-20) write(*,*) "Val VFP,i = ", vfP(1,1,1,i), i
       !if (vfP(1,1,1,i) .eq. 0.0) write(*,*) "Zero P V i = ",vfP(1,1,1,i),pRegion%grid%vol(i), i  
end do


!Super Particle Load Factor
IF ( global%myProcid == MASTERPROC) write(*,*) "PFINIT: Sending SPL to particles"
DO i = 1,pRegion%grid%nCells
       
        
       do lz=1,2
       do ly=1,2
       do lx=1,2
        xGrid(lx,ly,lz,i) = SPL(i)
end do
end do
end do 
end do
!if (pf_fluidInit .eq. .true. ) call ppiclf_solve_InterpFieldUser(PPICLF_R_JSPL,xGrid)

    DEALLOCATE(xGrid,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error

    DEALLOCATE(SPL,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:xGrid')
    END IF ! global%error


!Turn off for settling
if (1.eq.1 .and. pf_fluidInit .eqv. .true.) then!!!!!!!!!!!!!!!!!!
DO icg = 1,pGrid%nCellsTot
     vFrac = 1.0_RFREAL - pRegion%mixt%piclVF(i)
    pRegion%mixt%cv(CV_MIXT_DENS,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_DENS,icg)
    pRegion%mixt%cv(CV_MIXT_XMOM,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_XMOM,icg)
    pRegion%mixt%cv(CV_MIXT_YMOM,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_YMOM,icg)
    pRegion%mixt%cv(CV_MIXT_ZMOM,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_ZMOM,icg)
    pRegion%mixt%cv(CV_MIXT_ENER,icg) = vFrac*pRegion%mixt%cv(CV_MIXT_ENER,icg)
    if (pRegion%mixt%cv(CV_MIXT_DENS,icg) .le. 0.0) then
         write(*,*) "FAILED",pRegion%mixt%cv(CV_MIXT_DENS,icg)      
    end if    

END DO ! icg
endif


!total_vol = 0
!do i=1,pRegion%grid%nCells
!        if (volp(i) .gt. 0.0) then
!                ratio = pRegion%grid%vol(i)/ &
!                        (global%pi/6.0*(dp_min)**3)
!                tester = volp(i) !TESTH* ratio
!                total_vol = total_vol + tester
       !         write(*,*) "BADPQ = ", tester
!        end if
!end do

if (total_vol .gt. 0.0) then
        write(*,*) "BADPQ = ", total_vol/(global%pi/6.0*(dp_min)**3)
end if

call ppiclf_solve_InterpFieldUser(PPICLF_R_JPHIP,vfP)

DEALLOCATE(vfP,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:zGrid')
    END IF ! global%error
DEALLOCATE(volp,STAT=errorFlag)
    global%error = errorFlag
    IF ( global%error /= ERR_NONE ) THEN
      CALL ErrorStop(global,ERR_DEALLOCATE,__LINE__,'PPICLF:zGrid')
    END IF ! global%error



end if !Check seg fualt turn off for now



!    call ppiclf_exittr('Exit ppiclf $',0.0,0)
#endif
!PICL MOVE HERE END
  
! Rahul BUG FIX - Dependent vars are updated after phi communication. This is 
!                 necessary to update virtual cells
  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    CALL RFLU_SetVarsWrapper(pRegion,1,pRegion%grid%nCellsTot)   
      
    IF ( RFLU_NSCBC_DecideHaveNSCBC(pRegion) .EQV. .TRUE. ) THEN
      CALL RFLU_BXV_SetDependentVars(pRegion)
    END IF ! RFLU_NSCBC
  END DO ! iReg 
! Rahul - End BUG FIX

#ifdef GENX
! ******************************************************************************
! Initialize Rocman
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)        
    CALL RFLU_GENX_InitRocman(pRegion,genxHandle)     
  END DO ! iReg

! ******************************************************************************
! Check that have coupled boundaries defined, must be done after field flags
! are created
! ******************************************************************************

! TEMPORARY - Breaks intermittently, without this you need to make sure that
!             your input files are correct; there is no checking...
!  CALL RFLU_CheckCouplingInput(pLevel%regions)
! END TEMPORARY

! ******************************************************************************
! Put initial data into GENX buffers. NOTE no longer get data to fill incoming
! buffers to make sure restarts are always the same regardless of the time. 
! Incoming buffers are filled at beginning of time-stepping routine
! ******************************************************************************

  DO iReg = 1,global%nRegionsLocal
     pRegion => pLevel%regions(iReg)
     IF ( global%currentTime /= 0.0_RFREAL ) THEN
        CALL RFLU_GENX_InitBFLAG(pRegion)
        CALL RFLU_GetBoundaryValues(pLevel%regions(iReg))
     END IF
     CALL RFLU_PutBoundaryValues(pLevel%regions(iReg))
  ENDDO ! iReg
#endif

! ******************************************************************************
! Set options for grid motion, must be done after boundary normals are built 
! and data fields for communication are created
! ******************************************************************************

! TEMPORARY - Disabled until merged with changes for periodic and symmetry
!             boundaries
!  IF ( moveGrid .EQV. .TRUE. ) THEN 
!    DO iReg = 1,global%nRegionsLocal
!      pRegion => pLevel%regions(iReg)
!      CALL RFLU_SetMoveGridOptions(pRegion) 
!    END DO ! iReg
!  END IF ! moveGrid
! END TEMPORARY

! ******************************************************************************
! Find probe cells and print information on probe locations
! ******************************************************************************

  IF ( global%nProbes > 0 ) THEN 
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)
      CALL RFLU_FindProbeCells(pRegion)
    END DO ! iReg
    
    CALL RFLU_PrintProbeInfo(global)
  END IF ! global 

! ******************************************************************************
! Open files for convergence data, probe data, and OLES output
! ******************************************************************************

  CALL RFLU_OpenConverFile(global)
  
#ifndef GENX  
  IF ( moveGrid .EQV. .TRUE. ) THEN 
    CALL RFLU_OpenTotalMassFile(global)
  END IF ! moveGrid
#endif

  IF ( global%nProbes > 0 ) THEN 
    DO iReg = 1,global%nRegionsLocal
      pRegion => pLevel%regions(iReg)
      CALL RFLU_OpenProbeFiles(pRegion)
    END DO ! iReg
  END IF ! global      

! BEGIN TEMPORARY - only works for single regions at the moment
  DO iReg = 1,global%nRegionsLocal
    pRegion => pLevel%regions(iReg)
    
    IF ( pRegion%mixtInput%spaceDiscr == DISCR_OPT_LES ) THEN 
      CALL RFLU_OpenStatsFileOLES(global)
    END IF ! pRegion%mixtInput
  END DO ! iReg  
! END TEMPORARY

!begin BBR
 CALL RFLU_OpenPMFile(global)
! CALL RFLU_OpenIntegFile(global)
!end BBR

#ifndef GENX
! ******************************************************************************
! Compute integral values and write to file
! ******************************************************************************
      
  global%massIn  = 0.0_RFREAL
  global%massOut = 0.0_RFREAL    
      
  IF ( moveGrid .EQV. .TRUE. ) THEN
    CALL RFLU_ComputeIntegralValues(pLevel%regions) 
    CALL WriteTotalMass(pLevel%regions)
  END IF ! moveGrid
#endif  

#ifdef PETSC
! ******************************************************************************
! Initialize PETSc-related arrays and quantities
! ******************************************************************************
      
!  pRegion => pLevel%regions(1)
!  CALL RFLU_PETSC_Init(global)
!  CALL RFLU_PETSC_CreatePoisson(pRegion)
!  CALL RFLU_PETSC_BuildPoisson(pRegion)
!  CALL RFLU_PETSC_SetSolverPoisson(pRegion)
!  CALL RFLU_PETSC_SolvePressurePoisson(pRegion)
!  CALL RFLU_PETSC_Finalize(global)

  IF ( global%solverType == SOLV_IMPLICIT_NK ) THEN
    pRegion => pLevel%regions(1)
    CALL RFLU_PETSC_Init(global)
    CALL RFLU_PETSC_CreateVectors(pRegion)
    CALL RFLU_PETSC_CreateJacobian(pRegion)
  END IF ! global%solverType
#endif
 

! ******************************************************************************
! End
! ******************************************************************************

 !BBR - begin
 !CALL MPI_Barrier(global%mpiComm,errorFlag)
 !WRITE(777,*) "Line 1756...so far so good"
 !BBR - end

  IF ( global%myProcid == MASTERPROC .AND. & 
       global%verbLevel /= VERBOSE_NONE ) THEN
    WRITE(STDOUT,'(A)') SOLVER_NAME    
    WRITE(STDOUT,'(A,1X,A)') SOLVER_NAME,'Initialization done.'
    WRITE(STDOUT,'(A)') SOLVER_NAME 
  END IF ! global%myProcid

  CALL DeregisterFunction(global)

END SUBROUTINE RFLU_InitFlowSolver

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: RFLU_InitFlowSolver.F90,v $
! Revision 1.2  2017/01/23 16:48:31  rahul
! Fixed a bug related to updating virtual cells with dependent variables
! after the volume fraction is communicated.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.12  2009/09/28 14:21:55  mparmar
! Modified gradient computation for axisymm computation, updating virtual cells gradient data
!
! Revision 1.11  2009/08/28 18:29:48  mtcampbe
! RocfluMP integration with Rocstar and some makefile tweaks.  To build
! Rocstar with new Rocflu:
! make ROCFLU=RocfluMP
! To build Rocstar with the new RocfluND:
! make ROCFLU=RocfluMP HYPRE=/the/hypre/install/path
!
! Revision 1.10  2009/07/08 20:54:00  mparmar
! Removed RFLU_ModHouMaheshBoundCond
!
! Revision 1.9  2009/07/08 19:12:15  mparmar
! Added allocation for absorbing layer and call for RFLU_MVF_InitPatchVelAccel
!
! Revision 1.8  2008/12/06 08:43:48  mtcampbe
! Updated license.
!
! Revision 1.7  2008/11/19 22:17:00  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.6  2008/03/27 12:13:25  haselbac
! Added axisymmetry capability
!
! Revision 1.5  2008/01/19 20:19:45  haselbac
! Added calls to PLAG_DecideHaveSurfStats
!
! Revision 1.4  2007/12/03 16:34:35  mparmar
! Removed RFLU_CreatePatchVelocity
!
! Revision 1.3  2007/11/28 23:05:32  mparmar
! Allocating SOLV_IMPLICIT_HM related arrays
!
! Revision 1.2  2007/06/18 18:09:44  mparmar
! Added initialization of moving reference frame
!
! Revision 1.1  2007/04/09 18:49:57  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:01  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.119  2007/02/27 13:12:29  haselbac
! Adapted to changes in RFLU_ModDimensionality
!
! Revision 1.118  2006/12/21 12:22:12  haselbac
! Added call to check patch-bc consistency
!
! Revision 1.117  2006/10/20 21:32:23  mparmar
! Added calls to create and compute global thrustFlags
!
! Revision 1.116  2006/08/19 15:46:27  mparmar
! Changed because of NSCBC implementation
!
! Revision 1.115  2006/08/18 14:04:50  haselbac
! Added calls to build AVFace2Patch list, read transforms, cosmetics
!
! Revision 1.114  2006/06/23 21:37:47  mtcampbe
! Added a GetBoundary call for genx restarts - makes surface data correct.
!
! Revision 1.113  2006/06/16 19:39:30  mtcampbe
! Fixed MPI_Init sequencing for Rocstar/Native initialization
!
! Revision 1.112  2006/06/06 21:36:04  mtcampbe
! Remeshing/Rocstar mods DummyCell population, time 0 inits
!
! Revision 1.111  2006/04/07 16:04:03  haselbac
! Adapted to changes in bf2c wts computation
!
! Revision 1.110  2006/04/07 15:19:22  haselbac
! Removed tabs
!
! Revision 1.109  2006/04/07 14:53:27  haselbac
! Adapted to changes in bface stencil routines
!
! Revision 1.108  2006/03/25 22:03:26  haselbac
! Fix comment
!
! Revision 1.107  2006/03/25 22:02:51  haselbac
! Changes bcos of sype patches: Comment out Rxx2Sxx maps and face 
! reorientation, compute patch normals
!
! Revision 1.106  2006/03/09 14:10:03  haselbac
! Now call wrapper routines for stencils
!
! Revision 1.105  2006/02/06 23:55:54  haselbac
! Added comm argument to RFLU_InitGlobal
!
! Revision 1.104  2006/01/10 05:04:58  wasistho
! moved InitStatistics to before InitRocman
!
! Revision 1.103  2006/01/06 22:15:36  haselbac
! Adapted to name changes
!
! Revision 1.102  2006/01/03 06:32:03  wasistho
! moved Genx new-attr and registr. routines to wrappers
!
! Revision 1.101  2005/12/24 21:34:49  haselbac
! Added computation of ICT tolerance
!
! Revision 1.100  2005/12/24 02:25:31  wasistho
! moved statistics mapping befor Genx registration
!
! Revision 1.99  2005/12/24 02:15:22  wasistho
! activated statistics mapping
!
! Revision 1.98  2005/12/01 17:14:36  fnajjar
! Moved call to RFLU_RandomInit after reading user input
!
! Revision 1.97  2005/11/10 16:51:29  fnajjar
! Added plagUsed IF statement around PLAG routines
!
! Revision 1.96  2005/11/04 14:07:24  haselbac
! Renamed dim check routine, added dim geom check routine
!
! Revision 1.95  2005/10/28 19:18:15  haselbac
! Added check for nProcs not being equal to nProcAlloc
!
! Revision 1.94  2005/10/27 19:20:19  haselbac
! Adapted to changes in stencil routine names
!
! Revision 1.93  2005/10/25 19:39:23  haselbac
! Added IF on forceFlag
!
! Revision 1.92  2005/10/05 14:19:23  haselbac
! Adapted to changes in stencil mods, added call to create and comp bface wts
!
! Revision 1.91  2005/09/14 15:59:33  haselbac
! Minor clean-up
!
! Revision 1.90  2005/09/13 20:39:39  mtcampbe
! Moved profiling call into this file
!
! Revision 1.89  2005/08/03 18:42:39  hdewey2
! Enclosed PETSc init calls inside IF
!
! Revision 1.88  2005/08/02 18:57:58  hdewey2
! Temporarily commented out USE RFLU_ModPETScPoisson
!
! Revision 1.87  2005/08/02 18:29:17  hdewey2
! Added init of PETSc data and quantities
!
! Revision 1.86  2005/07/01 15:15:17  haselbac
! Added setting of Roccom verbosity level
!
! Revision 1.85  2005/06/09 20:30:14  haselbac
! Disabled checking of move grid options
!
! Revision 1.84  2005/05/18 22:13:27  fnajjar
! ACH: Added creation of iPclSend buffers, now use nFacesAV
!
! Revision 1.83  2005/05/12 18:02:53  haselbac
! Removed call to COM_set_verbose
!
! Revision 1.82  2005/04/29 23:02:54  haselbac
! Added building of avf2b list
!
! Revision 1.81  2005/04/29 13:00:56  haselbac
! Fixed bug in name of routine for opening probe files
!
! Revision 1.80  2005/04/29 12:49:31  haselbac
! Added USE RFLU_ModProbes, removed interfaces for probe routines
!
! Revision 1.79  2005/04/15 16:31:18  haselbac
! Removed calls to XyzEdge2RegionDegrList routines
!
! Revision 1.78  2005/04/15 15:07:18  haselbac
! Converted to MPI, integrated MPI code with GENx
!
! Revision 1.77  2005/03/09 15:08:25  haselbac
! Added dimensionality check
!
! Revision 1.76  2005/01/18 15:18:18  haselbac
! Commented out COMM calls for now
!
! Revision 1.75  2005/01/14 21:35:29  haselbac
! Added calls to create and read comm lists
!
! Revision 1.74  2005/01/13 21:40:29  haselbac
! Bug fix in setting pRegion for PETSc testing
!
! Revision 1.73  2005/01/07 19:25:56  fnajjar
! Added call to PLAG_ReadSurfStatsWrapper
!
! Revision 1.72  2005/01/03 15:58:34  haselbac
! Adapted to changes in RFLU_ModStencils
!
! Revision 1.71  2004/12/28 20:28:13  wasistho
! moved statistics routines into module ModStatsRoutines
!
! Revision 1.70  2004/12/21 23:34:32  fnajjar
! Added definition of pMixtInput pointer in incompressible part
!
! Revision 1.69  2004/12/21 15:05:11  fnajjar
! Included calls for PLAG surface statistics
!
! Revision 1.68  2004/12/19 15:49:50  haselbac
! Added incompressible stuff
!
! Revision 1.67  2004/12/04 03:34:04  haselbac
! Adapted to changes in RFLU_ModCellMapping
!
! Revision 1.66  2004/11/29 17:17:20  wasistho
! use ModInterfacesStatistics
!
! Revision 1.65  2004/11/14 19:47:15  haselbac
! Now call RFLU_SetVarsWrapper instead of UpdateDependentVarsMP
!
! Revision 1.64  2004/11/03 17:05:08  haselbac
! Removed HACK_PERIODIC ifdef and call to RFLU_GENX_CreateAttrGridVol
!
! Revision 1.63  2004/11/02 02:32:51  haselbac
! Added call to RFLU_SetVarInfoWrapper
!
! Revision 1.62  2004/10/21 15:54:13  haselbac
! Added ifdef for static linking
!
! Revision 1.61  2004/10/20 15:01:13  haselbac
! Bug fix: Changes made to allow compilation within GENx without CHARM=1
!
! Revision 1.60  2004/10/19 19:29:19  haselbac
! Substantial changes because of new GENX logic
!
! Revision 1.59  2004/07/06 15:14:52  haselbac
! Adapted to changes in libflu and modflu, cosmetics
!
! Revision 1.58  2004/06/25 20:08:25  haselbac
! Added call to RFLU_SetRestartTimeFlag
!
! Revision 1.57  2004/06/22 15:55:05  haselbac
! Bug fix for running serial jobs with Charm code
!
! Revision 1.56  2004/06/16 20:01:09  haselbac
! Added allocation of patch and force and moment coeffs, cosmetics
!
! Revision 1.55  2004/06/07 23:09:57  wasistho
! moved statistics mapping from initStatistics to before initGenxInterfaces
!
! Revision 1.54  2004/04/14 02:09:24  haselbac
! Added initialization of grid-speed factor
!
! Revision 1.53  2004/02/26 21:02:09  haselbac
! Added PLAG support, changed alloc logic, added updateDependentVarsMP
!
! Revision 1.52  2004/01/29 22:59:22  haselbac
! Added call to reading of bc data file, clean-up
!
! Revision 1.51  2004/01/22 16:04:33  haselbac
! Changed declaration to eliminate warning on ALC
!
! Revision 1.50  2003/12/07 04:59:22  jiao
! When GENX is defined, added a local variable "levels" to be
! used when calling RFLU_BuildDataStruct and RFLU_AssignRegionMapping
! to work around a complaint of PGI compiler.
!
! Revision 1.49  2003/12/04 03:30:04  haselbac
! Added calls for stencils and weights, cleaned up
!
! Revision 1.48  2003/11/25 21:04:41  haselbac
! Added call to RFLU_ReadBcInputFileWrapper, commented out communication
!
! Revision 1.47  2003/11/03 03:51:31  haselbac
! Added call to build boundary-face gradient access list
!
! Revision 1.46  2003/09/12 21:35:51  haselbac
! Fixed bug: Needed fem files for 1 proc cases
!
! Revision 1.45  2003/08/29 22:48:33  haselbac
! mpif.h should not be included directly here...
!
! Revision 1.44  2003/08/28 20:27:19  olawlor
! Minor tweaks for non-genx compilation--
!    - We need mpif.h
!    - Read input files from current directory, not Rocflu/inDir/.
!
! Revision 1.43  2003/08/27 15:38:00  haselbac
! Removed single and double quotes from last comments
!
! Revision 1.42  2003/08/26 22:48:56  olawlor
! Changes to startup sequence for latest Charm++/FEM framework:
!   - Call FEM_Init with an MPI communicator to set up FEM.
!     This replaces the old FEM_Attach call, and corresponding
!     -fem parameter in genx.C.
!
!   - Get partition number and number of partitions from MPI.
!     The FEM partition numbers are not needed any more.
!
!   - Read in the FEM input files yourself, using FEM_Mesh_read.
!     This allows the FEM input files to reside in the Rocflu
!     directory.
!
! Revision 1.41  2003/08/13 20:29:01  haselbac
! Fixed bug with writing probe data within GENx
!
! Revision 1.40  2003/07/22 15:39:50  haselbac
! Added Nullify routines
!
! Revision 1.39  2003/06/20 22:35:43  haselbac
! Added call to RFLU_ReadRestartInfo
!
! Revision 1.38  2003/06/04 20:05:53  jferry
! re-worked implementation of TBCs in unstructured code
!
! Revision 1.37  2003/05/16 21:51:29  mtcampbe
! ACH: Temporarily disabled checking of coupling input
!
! Revision 1.36  2003/05/13 23:49:46  haselbac
! Changed format for writing out number of procs
!
! Revision 1.35  2003/05/01 14:11:19  haselbac
! Added call to RFLU_CheckCouplingInput
!
! Revision 1.34  2003/04/24 15:43:13  haselbac
! Adapted interface to RFLU_PutBoundaryValues
!
! Revision 1.33  2003/04/12 21:38:08  haselbac
! Added setting of FEMRocfluGrid, verb and check level now read
!
! Revision 1.32  2003/04/07 14:26:50  haselbac
! Added cell-to-face list and probe info calls
!
! Revision 1.31  2003/03/31 16:17:21  haselbac
! Added writing of version string
!
! Revision 1.30  2003/03/25 19:16:58  haselbac
! Changed calling sequence for reorientation of AV faces
!
! Revision 1.29  2003/03/18 21:34:28  haselbac
! Modified call to RFLU_AllocateMemoryWrapper
!
! Revision 1.28  2003/03/15 18:38:12  haselbac
! Some adaptations for || gm, added serial mapping imposition
!
! Revision 1.27  2003/02/17 19:31:12  jferry
! Implemented portable random number generator ModRandom
!
! Revision 1.26  2003/01/31 13:53:57  haselbac
! Removed call to RFLU_GetBoundaryValues if not restarting
!
! Revision 1.25  2003/01/28 14:40:33  haselbac
! Extensive reorganization to get GENx restart to work with gm
!
! Revision 1.24  2002/11/27 20:25:14  haselbac
! Moved RFLU_SetMoveGridOptions to after having read bc file
!
! Revision 1.23  2002/11/26 15:27:55  haselbac
! Added RFLU_SetMoveGridOptions
!
! Revision 1.22  2002/11/26 15:18:39  haselbac
! Added ifdef to opening of total mass file
!
! Revision 1.21  2002/11/08 21:32:01  haselbac
! Added opening of total-mass file
!
! Revision 1.20  2002/10/27 19:13:03  haselbac
! New calls for grid motion, cosmetic redesign
!
! Revision 1.19  2002/10/19 16:15:19  haselbac
! Cosmetic changes to output
!
! Revision 1.18  2002/10/17 22:35:03  jiao
! ACH: Deactivate RFLU_PrintFlowInfo with GENX, gives small discrepancy for restart
!
! Revision 1.17  2002/10/17 14:14:00  haselbac
! Removed RFLU_GetBValues: Moved to RFLU_FlowSolver (discussion with Jim J.)
!
! Revision 1.16  2002/10/16 21:17:37  haselbac
! Added writing of header when running
!
! Revision 1.15  2002/10/12 14:59:13  haselbac
! Changed order of BC and GENX, added call to RFLU_GetBValues
!
! Revision 1.14  2002/10/05 19:22:53  haselbac
! GENX integration, added interface for RFLU_CreateFields
!
! Revision 1.13  2002/09/09 16:27:12  haselbac
! CVS appears to have got confused...
!
! Revision 1.12  2002/09/09 15:51:56  haselbac
! RFLU_HACK_PeriodicCellMapping.F90
!
! Revision 1.11  2002/08/16 21:33:48  jblazek
! Changed interface to MixtureProperties.
!
! Revision 1.10  2002/07/25 14:20:51  haselbac
! Added call to OLES routines, completed CHARM code segment
!
! Revision 1.9  2002/06/27 15:28:22  haselbac
! Added CHARM stuff, deleted MPI stuff
!
! Revision 1.8  2002/06/14 21:54:35  wasistho
! Added time avg statistics
!
! Revision 1.7  2002/06/14 20:21:26  haselbac
! Deleted ModLocal, changed local%nRegions to global%nRegionsLocal
!
! Revision 1.6  2002/06/10 21:31:14  haselbac
! Added call to RFLU_PrintHeader, changed flag to CHECK_UNIFLOW
!
! Revision 1.5  2002/06/05 18:59:53  haselbac
! Added RFLU_PrintGridInfo and changed version number
!
! Revision 1.4  2002/05/04 17:10:20  haselbac
! Added memory allocation, mixture properties, and checking
!
! Revision 1.3  2002/04/11 19:03:45  haselbac
! Added calls and cosmetic changes
!
! Revision 1.2  2002/03/26 19:23:39  haselbac
! Some cleaning and added reading of input files
!
! Revision 1.1  2002/03/14 19:12:48  haselbac
! Initial revision
!
! ******************************************************************************

