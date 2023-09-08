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
! Purpose: Read in user input related to flow initialization.
!
! Description: None.
!
! Input: 
!   regions     Data for regions
!
! Output: None.
!
! Notes: None.
!
! ******************************************************************************
!
! $Id: ReadInitFlowSection.F90,v 1.4 2016/03/22 19:32:14 fred Exp $
!
! Copyright: (c) 2002-2006 by the University of Illinois
!
! ******************************************************************************

SUBROUTINE ReadInitPPICLFSection(global)

  USE ModDataTypes
  USE ModGlobal, ONLY: t_global
  USE ModError
  USE ModParameters
  
  USE ModInterfaces, ONLY: ReadSection
  
  IMPLICIT NONE

! ******************************************************************************
! Declarations and definitions
! ******************************************************************************

! ==============================================================================
! Arguments
! ==============================================================================

  TYPE(t_global), POINTER :: global

! ==============================================================================
! Local variables
! ==============================================================================

  INTEGER :: i,nVals
  INTEGER, PARAMETER :: NVALS_MAX = 30

  CHARACTER(50) :: keys(NVALS_MAX)
  LOGICAL :: defined(NVALS_MAX)
  REAL(RFREAL) :: vals(NVALS_MAX)

! ******************************************************************************
! Start
! ******************************************************************************

  CALL RegisterFunction(global,'ReadInitPPICLFSection',__FILE__)

! specify keywords and search for them

#ifdef PICL
  nVals = NVALS_MAX

  ! these properties apply to ppiclF particles
  keys( 1) = 'INITTYPE' ! type of initialization routine. 1=line in y
  keys( 2) = 'X' ! x location of y line for inittype=1
  keys( 3) = 'YMIN' ! y minimum location line for inittype=1
  keys( 4) = 'YMAX' ! y maximum location line for inittype=1
  keys( 5) = 'Z' ! z location of y line for inittype=1
  keys( 6) = 'VX' ! initial x velocity
  keys( 7) = 'VY' ! initial y velocity
  keys( 8) = 'VZ' ! initial z velocity
  keys( 9) = 'NPART' ! number of particles
  keys(10) = 'MU' ! viscosity for TAB breakup model
  keys(11) = 'SIGMA' ! surface tension for TAB breakup model
  keys(12) = 'RHO' ! density
  keys(13) = 'D' ! diameter
  keys(14) = 'DCUTOFF' ! TAB diameter cutoff 
  keys(15) = 'DMAX' ! max diameter
  keys(16) = 'ND' ! number of diameters to vary
  keys(18) = 'RHOMAX'
  keys(19) = 'SIGMAMAX'
  keys(20) = 'NRHO'
  keys(21) = 'NSIGMA'
  keys(22) = 'MUMAX'
  keys(23) = 'NMU'
  keys(24) = 'BREAKUPMODEL' ! 0=None, 1=TAB, 2=Ashgriz
  keys(27) = 'CFTAB' ! C_F in TAB model, C2=2CF, C2 from Schmehl
  keys(28) = 'NCHILD'
  keys(29) = 'ZPFFACTOR'
  keys(30) = 'NDUPLICATES'

  CALL ReadSection(global,IF_INPUT,nVals,keys,vals,defined )

  if (defined(1) .EQV. .FALSE.) then
    vals(1) = 0
  else
    if (defined(2) .EQV. .TRUE.) then 
        global%ppiclFInitX = vals(2)
    else
        global%ppiclFInitX = 0.0
    end if

    if (defined(3) .EQV. .TRUE.) then
        global%ppiclFInitYmin = vals(3)
    else
        global%ppiclFInitYmin = 0.0
    end if

    if (defined(4) .EQV. .TRUE.) then 
        global%ppiclFInitYmax = vals(4)
    else
        global%ppiclFInitYmax = 1.0
    end if

    if (defined(5) .EQV. .TRUE.) then
        global%ppiclFInitZ = vals(5)
    else
        global%ppiclFInitZ = 0.0005
    end if

    if (defined(6) .EQV. .TRUE.) then
        global%ppiclFInitVX = vals(6)
    else
        global%ppiclFInitVX = 0.0
    end if

    if (defined(7) .EQV. .TRUE.) then 
        global%ppiclFInitVY = vals(7)
    else
        global%ppiclFInitVY = 0.0
    end if

    if (defined(8) .EQV. .TRUE.) then 
        global%ppiclFInitVZ = vals(8)
    else
        global%ppiclFInitVz = 0.0
    end if

    if (defined(9) .EQV. .TRUE.) then 
        global%ppiclFInitNPart = vals(9)
    else
        global%ppiclFInitNPart = 0
    end if

    if (defined(10) .EQV. .TRUE.) then 
        global%ppiclFInitMu = vals(10)
    else
        global%ppiclFInitMu = 8.891*10.0**(-4.0) ! default = water
    end if

    if (defined(11) .EQV. .TRUE.) then
        global%ppiclFInitSigma = vals(11)
    else
        global%ppiclFInitSigma = 71.99*10.0**(-3.0) ! default = water
    end if

    if (defined(12) .EQV. .TRUE.) then 
        global%ppiclFInitRho = vals(12)
    else
        global%ppiclFInitRho = 1000.0 ! default = water
    end if

    if (defined(13) .EQV. .TRUE.) then 
        global%ppiclFInitD = vals(13)
    else
        global%ppiclFInitD = 1.0*10**(-3.0) ! default = 1mm
    end if

    if (defined(14) .EQV. .TRUE.) then
        global%ppiclFInitDCutoff = vals(14)
    else
        global%ppiclFInitDCutoff = 0.0
    end if

    if (defined(15) .EQV. .TRUE.) then
        global%ppiclFInitDMAX = vals(15)
    else
        global%ppiclFInitDMAX = 0.0
    end if

    if (defined(16) .EQV. .TRUE.) then
        global%ppiclFInitND = vals(16)
    else
        global%ppiclFInitND = 1
    end if

    if (defined(18) .EQV. .TRUE.) then
        global%ppiclFInitRHOMAX = vals(18)
    else
        global%ppiclFInitRHOMAX = 0.0
    end if

    if (defined(19) .EQV. .TRUE.) then
        global%ppiclFInitSigmaMAX = vals(19)
    else
        global%ppiclFInitSigmaMAX = 0.0
    end if

    if (defined(20) .EQV. .TRUE.) then
        global%ppiclFInitNRho = vals(20)
    else
        global%ppiclFInitNRho = 1
    end if

    if (defined(21) .EQV. .TRUE.) then
        global%ppiclFInitNSigma = vals(21)
    else
        global%ppiclFInitNSigma = 1
    end if

    if (defined(22) .EQV. .TRUE.) then
        global%ppiclFInitMuMAX = vals(22)
    else
        global%ppiclFInitMuMAX = 0.0
    end if

    if (defined(23) .EQV. .TRUE.) then
        global%ppiclFInitNMu = vals(23)
    else
        global%ppiclFInitNMu = 1
    end if

    if (defined(24) .EQV. .TRUE.) then
        global%ppiclFInitBreakupModel = vals(24)
    else
        global%ppiclFInitBreakupModel = 0
    end if

    if (defined(27) .EQV. .TRUE.) then
        global%ppiclFInitCFTAB = vals(27)
    else
        global%ppiclFInitCFTAB = 1.0/3.0 ! default from TAB paper
    end if

    if (defined(28) .EQV. .TRUE.) then
        global%ppiclFInitNChild = vals(28)
    else
        global%ppiclFInitNChild = 100 ! sane default for number of child droplets for each breakup
    end if

    if (defined(29) .EQV. .TRUE.) then
        global%ppiclFInitZpfFactor = vals(29)
    else
        global%ppiclFInitZpfFactor = 1
    end if

    if (defined(30) .EQV. .TRUE.) then
        global%ppiclFInitNDuplicates = vals(30)
    else
        global%ppiclFInitNDuplicates = 1
    end if
  end if

  global%ppiclFInitType = vals(1)
#endif

! ******************************************************************************
! End
! ******************************************************************************

  CALL DeregisterFunction(global)

END SUBROUTINE ReadInitPPICLFSection

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: ReadInitFlowSection.F90,v $
! Revision 1.4  2016/03/22 19:32:14  fred
! Adding Glasser's model as a choice for collision modelling
!
! Revision 1.3  2016/02/06 17:24:25  fred
! Adding Variable JWL density RVAL into read list
!
! Revision 1.2  2016/02/03 20:38:34  rahul
! Added 3 additional keys to read dx, dy and dz of the shktb grid from the
! .inp file. This is a temporary fix.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:37  brollin
! New Stable version
!
! Revision 1.4  2009/07/08 19:11:22  mparmar
! Adapted to new init option
!
! Revision 1.3  2008/12/06 08:43:32  mtcampbe
! Updated license.
!
! Revision 1.2  2008/11/19 22:16:47  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.1  2007/04/09 18:48:32  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 17:59:25  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.9  2007/04/05 00:56:57  haselbac
! Added additional RVALxy params for 2p shocktube problems
!
! Revision 1.8  2006/04/07 15:19:15  haselbac
! Removed tabs
!
! Revision 1.7  2006/03/26 20:21:20  haselbac
! Added TEMP input argument
!
! Revision 1.6  2005/11/17 14:37:24  haselbac
! Added more RVAL variables
!
! Revision 1.5  2005/09/13 21:36:45  haselbac
! Adapted to new init option
!
! Revision 1.4  2005/04/20 14:38:36  haselbac
! Added more int and real vals
!
! Revision 1.3  2005/03/29 22:28:31  haselbac
! Added setting of initFlowFlag for combo option
!
! Revision 1.2  2005/03/22 03:32:39  haselbac
! Added initialization of integer and real helper variables
!
! Revision 1.1  2004/12/01 16:50:26  haselbac
! Initial revision after changing case
!
! Revision 1.13  2004/11/14 19:35:42  haselbac
! Added initialization for incompressible fluid model
!
! Revision 1.12  2004/07/28 16:41:32  haselbac
! Bug fix: Initial values not assigned to all regions
!
! Revision 1.11  2004/04/08 03:15:24  wasistho
! nDummyCells in Rocflo read from INITFLOW section
!
! Revision 1.10  2003/11/20 16:40:35  mdbrandy
! Backing out RocfluidMP changes from 11-17-03
!
! Revision 1.7  2003/09/15 00:36:19  haselbac
! Added hard-code option as input
!
! Revision 1.6  2003/05/16 02:27:43  haselbac
! Removed KIND=RFREAL from NINT statements
!
! Revision 1.5  2003/05/15 02:57:02  jblazek
! Inlined index function.
!
! Revision 1.4  2003/03/25 19:15:17  haselbac
! Fixed bug in RCSIdentString
!
! Revision 1.3  2003/03/15 16:27:39  haselbac
! Added KIND qualifyer
!
! Revision 1.2  2003/02/13 22:30:54  jferry
! removed RFLU_ prefix in RegisterFunction
!
! Revision 1.1  2003/01/28 16:12:56  haselbac
! Moved here from rfluprep
!
! Revision 1.4  2002/10/27 19:23:46  haselbac
! Removed tabs
!
! Revision 1.3  2002/10/07 14:11:29  haselbac
! Removed tabs
!
! Revision 1.2  2002/09/09 16:40:13  haselbac
! global and mixtInput now under regions
!
! Revision 1.1  2002/04/11 19:13:22  haselbac
! Initial revision
!
! ******************************************************************************

