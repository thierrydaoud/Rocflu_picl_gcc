################################################################################
#
# $Id: Filelist.txt,v 1.4 2016/08/11 16:36:39 rahul Exp $
#
# Purpose: Filelist for Rocpart module.
#
# Description: None.
#
# Notes: None.
#
# Copyright: (c) 2003 by the University of Illinois
#
################################################################################

SRCF90+=	PLAG_BuildVersionString.F90\
		PLAG_ModBcData.F90\
		PLAG_ModCheckVars.F90\
		PLAG_ModDataStruct.F90\
		PLAG_ModDimensions.F90\
		PLAG_ModEulerian.F90\
		PLAG_ModInflow.F90\
		PLAG_ModInterfaces.F90\
		PLAG_ModParameters.F90\
		PLAG_ModPlotting.F90\
		PLAG_ModReallocateMemory.F90\
		PLAG_ModRkInit.F90\
		PLAG_ModSurfStats.F90\
		PLAG_CalcBreakup.F90\
		PLAG_CalcDerivedVariables.F90\
		PLAG_CalcRhsPosition.F90\
		PLAG_CheckUserInput.F90\
		PLAG_ContinuousRandomWalk.F90\
		PLAG_DerivedInputValues.F90\
		PLAG_InitPatchData.F90\
		PLAG_InitInputValues.F90\
		PLAG_InjcSetInjection.F90\
		PLAG_InjcTileRKUpdate.F90\
		PLAG_InjcTileUpdate.F90\
		PLAG_InjcTileZeroRhs.F90\
		PLAG_IntrpMixtProperties.F90\
		PLAG_NonCvUpdate.F90\
		PLAG_PrintUserInput.F90\
		PLAG_ReadDisPartInitSection.F90\
		PLAG_ReadDisPartnContSection.F90\
		PLAG_ReadDisPartSection.F90\
		PLAG_ReadInputFile.F90\
		PLAG_ReadPdfFromFile.F90\
		PLAG_ReflectParticleData.F90\
		PLAG_RkInit.F90\
		PLAG_RkUpdateWrapper.F90\
		PLAG_SetDependentVarsOld.F90\
		PLAG_UserInput.F90\
		PLAG_UpdateDataStruct.F90\
		PLAG_ZeroRhs.F90\
		PLAG_INRT_AllocMemTStep.F90\
		PLAG_INRT_DeallocMemTStep.F90\
		PLAG_RFLU_AllocMemSol.F90\
		PLAG_RFLU_AllocMemSolTile.F90\
		PLAG_RFLU_AllocMemTStep.F90\
		PLAG_RFLU_AllocMemTStepTile.F90\
		PLAG_RFLU_CalcForceInterface.F90\
                PLAG_RFLU_ComputeCellsContainingPcls.F90\
                PLAG_RFLU_ComputeLagPclsPerReg.F90\
                PLAG_RFLU_ComputeMaxImpulse.F90\
		PLAG_RFLU_ComputeSourcePint.F90\
                PLAG_RFLU_ComputeVolFrac.F90\
                PLAG_RFLU_ComputeVolFracGradL.F90\
		PLAG_RFLU_CorrectMixtProperties.F90\
		PLAG_RFLU_DeallocMemSol.F90\
		PLAG_RFLU_DeallocMemSolTile.F90\
		PLAG_RFLU_DeallocMemTStep.F90\
		PLAG_RFLU_DeallocMemTStepTile.F90\
		PLAG_RFLU_GetMixtPG.F90\
		PLAG_RFLU_GetMixtSD.F90\
		PLAG_RFLU_InitSolutionCyldet.F90\
		PLAG_RFLU_InitSolutionFile.F90\
		PLAG_RFLU_InitSolutionHardcode.F90\
		PLAG_RFLU_InitSolutionRandom.F90\
		PLAG_RFLU_InitSolutionScratch.F90\
		PLAG_RFLU_InitSolutionShktb.F90\
		PLAG_RFLU_InitSolutionXdsp.F90\
		PLAG_RFLU_InitSolFromSerial.F90\
		PLAG_RFLU_InitSolFromSerialCopy.F90\
		PLAG_RFLU_InitSolSerialWrapper.F90\
		PLAG_RFLU_InitSolSerial_1D.F90\
		PLAG_RFLU_InitSolSerial_2D.F90\
		PLAG_RFLU_InitSolSerial_3D.F90\
		PLAG_RFLU_InjcTileCalcRhs.F90\
		PLAG_RFLU_InjectionDriver.F90\
		PLAG_RFLU_InterParticleForce.F90\
		PLAG_RFLU_ModComm.F90\
		PLAG_RFLU_ModFindCells.F90\
                PLAG_RFLU_ModifyFlowFieldVolFrac.F90\
		PLAG_RFLU_ReadSolutionASCII.F90\
		PLAG_RFLU_ReadSolutionBinary.F90\
		PLAG_RFLU_ReadUnsteadyDataASCII.F90\
		PLAG_RFLU_ReadUnsteadyDataBinary.F90\
		PLAG_RFLU_ShiftUnsteadyData.F90\
		PLAG_RFLU_Update.F90\
		PLAG_RFLU_WriteSolutionASCII.F90\
		PLAG_RFLU_WriteSolutionBinary.F90\
		PLAG_RFLU_WriteUnsteadyDataASCII.F90\
		PLAG_RFLU_WriteUnsteadyDataBinary.F90

ifdef STATS
  SRCF90+=	PLAG_StatMapping.F90
endif

################################################################################
#
# RCS Revision history:
#
#   $Log: Filelist.txt,v $
#   Revision 1.4  2016/08/11 16:36:39  rahul
#   Added PLAG_RFLU_ComputeSourcePint.F90.
#
#   Revision 1.3  2016/05/06 00:03:17  rahul
#   Added new subroutine PLAG_RFLU_CalcForceInterface.F90.
#
#   Revision 1.2  2015/07/27 04:45:42  brollin
#   1) Corrected bug in RFLUCONV where global%gridFormat was used instead of global%gridSrcFormat
#   2) Implemented new subroutine for shock tube problems (Shktb)
#
#   Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
#   merged rocflu micro and macro
#
#   Revision 1.2  2014/07/21 16:41:20  subbu
#   Added PLAG_RFLU_InitSolutionFile
#
#   Revision 1.1.1.1  2014/05/05 21:47:47  tmish
#   Initial checkin for rocflu macro.
#
#   Revision 1.5  2007/05/16 22:21:47  fnajjar
#   Added new calls and deleted obsolete ones
#
#   Revision 1.4  2007/04/16 23:13:25  fnajjar
#   Deleted irrelevant interface routine
#
#   Revision 1.3  2007/04/15 02:34:50  haselbac
#   Added entry for PLAG_RFLU_InitSolSerial_2D.F90
#
#   Revision 1.2  2007/04/12 17:55:55  haselbac
#   Added entry for PLAG_ModPlotting.F90
#
#   Revision 1.1  2007/04/09 18:50:25  haselbac
#   Initial revision after split from RocfloMP
#
#   Revision 1.1  2007/04/09 18:01:32  haselbac
#   Initial revision after split from RocfloMP
#
#   Revision 1.51  2007/03/20 17:35:56  fnajjar
#   Added entry to PLAG_ModDimensions
#
#   Revision 1.50  2007/03/15 21:58:52  haselbac
#   Added/deleted entries
#
#   Revision 1.49  2007/03/12 23:33:53  haselbac
#   Added entry for PLAG_ModDataStruct.F90, some clean-up
#
#   Revision 1.48  2006/05/05 17:27:29  haselbac
#   Added entry for PLAG_RFLU_InitSolSerial.F90
#
#   Revision 1.47  2005/12/01 21:53:28  fnajjar
#   Added PLAG_ModCheckVars
#
#   Revision 1.46  2005/11/30 22:17:02  fnajjar
#   Added entries for routines for RFLU
#
#   Revision 1.45  2005/05/19 16:01:57  fnajjar
#   Added call to PLAG_ModRkInit
#
#   Revision 1.44  2005/05/18 22:14:22  fnajjar
#   Added entries for PLAG_RFLU_ModComm and PLAG_RFLU_InitSolFromSerial
#
#   Revision 1.43  2005/04/27 14:56:12  fnajjar
#   Included call to PLAG_RFLU_ModFindCells
#
#   Revision 1.42  2005/04/25 18:39:10  luca1
#   Imposed PDF from file option for random particle ejection
#
#   Revision 1.41  2005/03/11 02:21:42  haselbac
#   Added and removed PLAG_RFLU_FindCellsTrajXYZ routines
#
#   Revision 1.40  2005/01/08 20:51:24  fnajjar
#   Included calls to PLAG statistics
#
#   Revision 1.39  2005/01/01 21:33:28  haselbac
#   Added entry for PLAG_RFLU_FindCellsApte
#
#   Revision 1.38  2004/12/29 23:30:16  wasistho
#   prepared statistics for PLAG
#
#   Revision 1.37  2004/12/21 15:06:32  fnajjar
#   Added PLAG_ModSurfStats call
#
#   Revision 1.36  2004/12/01 21:11:54  fnajjar
#   Changed to upper case
#
#   Revision 1.35  2004/12/01 00:01:43  wasistho
#   added BuildVersionString
#
#   Revision 1.34  2004/11/17 16:46:42  haselbac
#   Removed PLAG_rkUpdate and PLAG_rkUpdateGeneric
#
#   Revision 1.33  2004/11/05 21:51:21  fnajjar
#   Added entries for particle-cell search routines
#
#   Revision 1.32  2004/11/04 16:42:16  fnajjar
#   Added entry to PLAG_SetDimensions
#
#   Revision 1.31  2004/10/11 22:09:11  haselbac
#   Renamed procedures
#
#   Revision 1.30  2004/10/10 20:06:03  fnajjar
#   Added call to PLAG_RFLU_InitSolutionRandom
#
#   Revision 1.29  2004/10/08 22:09:13  haselbac
#   Added entry for PLAG_RFLU_FindParticleCellsBrut
#
#   Revision 1.28  2004/08/23 23:07:19  fnajjar
#   Added binary IO
#
#   Revision 1.27  2004/08/20 23:27:13  fnajjar
#   Added Infrastructure for Plag prep tool
#
#   Revision 1.26  2004/07/28 18:59:20  fnajjar
#   Included calls to routines for dynamic memory reallocation
#
#   Revision 1.25  2004/07/26 19:01:01  fnajjar
#   Included call to PLAG_INRT_DeallocMemTStep routine
#
#   Revision 1.24  2004/07/26 17:05:51  fnajjar
#   moved allocation of inrtSources into Rocpart
#
#   Revision 1.23  2004/04/09 22:57:07  fnajjar
#   Added RFLO specific routines
#
#   Revision 1.22  2004/04/08 01:33:21  haselbac
#   Added entry for PLAG_ReflectParticleData
#
#   Revision 1.21  2004/03/26 21:26:26  fnajjar
#   Added new routines for RFLU-specific routines
#
#   Revision 1.20  2004/03/18 21:40:31  fnajjar
#   Added routines of sending-receiving MPI-based buffer data
#
#   Revision 1.19  2004/03/15 21:05:34  haselbac
#   Deleted/added file
#
#   Revision 1.18  2004/03/10 23:12:54  fnajjar
#   Included interfaces for MPI-based routines for corner-edge cells
#
#   Revision 1.17  2004/03/08 22:15:58  fnajjar
#   Added calls to injection routines in RFLU section
#
#   Revision 1.16  2004/03/05 23:16:41  haselbac
#   Added routines
#
#   Revision 1.15  2004/02/27 16:08:29  haselbac
#   Added RFLO entry for PLAG_derivedInputValues.F90
#
#   Revision 1.14  2004/02/26 21:02:11  haselbac
#   Added RFLU routines
#
#   Revision 1.13  2004/02/25 21:56:45  fnajjar
#   Included generic RKUpdate for PLAG
#
#   Revision 1.12  2004/02/10 21:23:08  fnajjar
#   Added call and interfaces for index mapping between corner-edge regions
#
#   Revision 1.11  2004/02/06 21:19:03  fnajjar
#   Added more routines to RFLU file list
#
#   Revision 1.10  2004/02/02 22:52:37  haselbac
#   Added routines to RFLU file list
#
#   Revision 1.9  2004/01/26 22:55:31  fnajjar
#   Included calls to routines for corner-edge data loading
#
#   Revision 1.8  2004/01/15 21:10:15  fnajjar
#   Added calls to corner-edge cell metrics routines
#
#   Revision 1.7  2003/11/12 21:34:00  fnajjar
#   Added Corner-Edge cells subroutine calls
#
#   Revision 1.6  2003/11/03 21:22:20  fnajjar
#   Added PLAG_copyFaceVectors
#
#   Revision 1.5  2003/09/13 20:14:21  fnajjar
#   Added infrastructure for Breakup model
#
#   Revision 1.4  2003/05/28 15:16:11  fnajjar
#   Removed obsolete PLAG_mixt calls as embedded in Rocinteract
#
#   Revision 1.3  2003/04/14 14:32:20  fnajjar
#   Added PLAG_initInputValues for proper initialization
#
#   Revision 1.2  2003/03/28 19:49:35  fnajjar
#   Added PLAG wrapper routines
#
#   Revision 1.1  2003/03/20 19:26:21  haselbac
#   Initial revision
#
################################################################################
