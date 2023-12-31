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
! Purpose: define parameters for discrete Lagrangian particles
!
! Description: none
!
! Notes: 
!   1. error codes are defined in the module ModError.F90
!   2. the parameters for cv and dv are extended by nCont since 
!      the datastructure has components based on the number of constituents
!      nCv=CV_PLAG_LAST+nCont
!   3. Similar extension applied to nEv=EV_PLAG_LAST+2*nCont
!
! ******************************************************************************
!
! $Id: PLAG_ModParameters.F90,v 1.2 2016/02/24 06:04:09 rahul Exp $
!
! Copyright: (c) 2002-2005 by the University of Illinois
!
! ******************************************************************************

MODULE PLAG_ModParameters

  IMPLICIT NONE

! Lagrangian Particles: PLAG ---------------------------------------------------
 
  INTEGER, PARAMETER :: CV_PLAG_XMOM        = 1, &       ! momentum components         
                        CV_PLAG_YMOM        = 2, &               
                        CV_PLAG_ZMOM        = 3, &               
                        CV_PLAG_ENER        = 4, &
                        CV_PLAG_XPOS        = 5, &       ! position components  
                        CV_PLAG_YPOS        = 6, &               
                        CV_PLAG_ZPOS        = 7, &
                        CV_PLAG_ENERVAPOR   = 8, &       ! vapor energy
                        CV_PLAG_LAST        = 8          ! final index
                                                         ! see note 2 for index above 9

  INTEGER, PARAMETER :: DV_PLAG_UVEL        = 1,  &      ! particle velocity components
                        DV_PLAG_VVEL        = 2,  &
                        DV_PLAG_WVEL        = 3,  &
                        DV_PLAG_TEMP        = 4,  &      ! static temperature
                        DV_PLAG_DIAM        = 5,  &      ! diameter
                        DV_PLAG_UVELMIXT    = 6,  &      ! mixture velocity components
                        DV_PLAG_VVELMIXT    = 7,  &
                        DV_PLAG_WVELMIXT    = 8,  &
                        DV_PLAG_DENSMIXT    = 9,  &      ! density
                        DV_PLAG_TEMPMIXT    = 10, &      ! static temperature
                        DV_PLAG_PRESMIXT    = 11, &      ! static pressure
                        DV_PLAG_LAST        = 12, &      ! final index
                        DV_PLAG_SOUNMIXT   = 12
                       
  INTEGER, PARAMETER :: TV_PLAG_MUELMIXT    = 1, &       ! laminar kinematic viscosity 
                                                         ! at particle location
                        TV_PLAG_TCOLMIXT    = 2, &       ! laminar thermal conductivity 
                        TV_PLAG_LAST        = 2          ! final index

  INTEGER, PARAMETER :: AIV_PLAG_PIDINI     = 1, &       ! particle initial id at creation
                        AIV_PLAG_REGINI     = 2, &       ! particle initial region at creation
                        AIV_PLAG_ICELLS     = 3, &       ! particle cell index
                        AIV_PLAG_BURNSTAT   = 4, &       ! particle burning status
                        AIV_PLAG_STATUS     = 5, &       ! particle search status
                        AIV_PLAG_LAST       = 5          ! final index

  INTEGER, PARAMETER :: ARV_PLAG_SPLOAD     = 1, &       ! superparticle loading
                        ARV_PLAG_DISTOT     = 2, &       ! total distance travelled
                        ARV_PLAG_UPRIME     = 3, &       ! particle flutuating x-velocity used in CRW
                        ARV_PLAG_VPRIME     = 4, &       ! particle flutuating y-velocity used in CRW
                        ARV_PLAG_WPRIME     = 5, &       ! particle flutuating z-velocity used in CRW
                        ARV_PLAG_LAST       = 5          ! final index

  INTEGER, PARAMETER :: EV_PLAG_DIA3        = 1,  &      ! eulerian-based particle variables
                        EV_PLAG_DIA4        = 2,  &
                        EV_PLAG_NDNS        = 3,  &      ! number density
                        EV_PLAG_UVEL        = 4,  &      ! velocity components
                        EV_PLAG_VVEL        = 5,  &      !
                        EV_PLAG_WVEL        = 6,  &      ! 
                        EV_PLAG_TEMP        = 7,  &      ! 
                        EV_PLAG_MFRC        = 8,  &      ! mass fraction
                        EV_PLAG_VFRC        = 9,  &      ! volume fraction
                        EV_PLAG_REYN        = 10, &      ! reynolds number
                        EV_PLAG_LAST        = 10         ! final index
                                                         ! see note 3 for index above 9

  INTEGER, PARAMETER :: PLAG_BC_INFLOW_LOGNORM   = 1, &     ! Inflow Diameter Distribution Models
                        PLAG_BC_INFLOW_LOGSKWD   = 2, &
                        PLAG_BC_INFLOW_PDF       = 3

  INTEGER, PARAMETER :: ZEROTH_ORDER        = 0, &       ! Order of accuracy for interpolation
                        FIRST_ORDER         = 1, &
                        SECOND_ORDER        = 2                                                

  INTEGER, PARAMETER :: PLAG_BREAKUP_NOMODEL = 0, &      ! Breakup Models
                        PLAG_BREAKUP_MODEL1  = 1

  INTEGER, PARAMETER :: PLAG_BREAKUP_NOWEBSWI = 0, &     ! Weber Switch for Breakup Model
                        PLAG_BREAKUP_WEBSWI1  = 1

  INTEGER, PARAMETER :: CV_TILE_MOMNRM      = 1, &       ! tile infrastructure         
                        CV_TILE_ENER        = 2, &
                        CV_TILE_LAST        = 2          ! final index

  INTEGER, PARAMETER :: DV_TILE_DIAM        = 1, &       ! 
                        DV_TILE_SPLOAD      = 2, & 
                        DV_TILE_POOLVOLD    = 3, & 
                        DV_TILE_COUNTDOWN   = 4, &          
                        DV_TILE_LAST        = 4          ! final index  

  INTEGER, PARAMETER :: PLAG_STATUS_KEEP    = 0, &       ! Status parameters
                        PLAG_STATUS_COMM    = 1, &
                        PLAG_STATUS_DELETE  = 2, &
                        PLAG_STATUS_LOST    = 3
  
  INTEGER, PARAMETER :: PLAG_BC_NOINFLOW      = 0, &
                        PLAG_BC_INFLOW_MODEL1 = 1, &     ! Inflow Models
                        PLAG_BC_INFLOW_CRE    = 2

  INTEGER, PARAMETER :: NPCLS_TOT_MIN     = 1000         ! Minimum Size of Particle DataStructure

  INTEGER, PARAMETER :: FIND_PCL_METHOD_TRAJ_FAST = 0, & 
                        FIND_PCL_METHOD_TRAJ_SAFE = 1, &        
                        FIND_PCL_METHOD_BRUTE     = 2, & 
                        ! Subbu - Add hardcoded pcl tracking flag
                        FIND_PCL_METHOD_HARDCODE  = 5, &
                        ! Subbu - End add hardcoded pcl tracking flag
                        FIND_PCL_METHOD_OCT       = 3, & 
                        FIND_PCL_METHOD_LOHNER    = 4

  INTEGER, PARAMETER :: PLAG_SURF_STATS_DIAM3 = 1, &     ! Surface Statistics DataStructure
                        PLAG_SURF_STATS_DIAM4 = 2, &
                        PLAG_SURF_STATS_THETA = 3, &
                        PLAG_SURF_STATS_MOME1 = 4, &
                        PLAG_SURF_STATS_MOME2 = 5, &
                        PLAG_SURF_STATS_MASS  = 6, &
                        PLAG_SURF_STATS_ENER  = 7, &
                        PLAG_SURF_STATS_LAST  = 7

  INTEGER, PARAMETER :: BIN_METHOD_LINEAR  = 1, &        ! Binning Method
                        BIN_METHOD_LOGNORM = 2 

END MODULE PLAG_ModParameters

! ******************************************************************************
!
! RCS Revision history:
!
! $Log: PLAG_ModParameters.F90,v $
! Revision 1.2  2016/02/24 06:04:09  rahul
! Added PLAG_BC_NOINFLOW parameter as a 3rd inflow model. This is added to
! bypass particle injection issue in shktb case.
!
! Revision 1.1.1.1  2015/01/23 22:57:50  tbanerjee
! merged rocflu micro and macro
!
! Revision 1.1.1.1  2014/07/15 14:31:38  brollin
! New Stable version
!
! Revision 1.7  2008/12/06 08:43:51  mtcampbe
! Updated license.
!
! Revision 1.6  2008/11/19 22:17:03  mtcampbe
! Added Illinois Open Source License/Copyright
!
! Revision 1.5  2007/05/16 22:30:35  fnajjar
! Modified parameters to be aligned with new bc datastructure
!
! Revision 1.4  2007/04/26 20:43:38  fnajjar
! Cleaned up DV_PLAG names removing irrelevant parameters
!
! Revision 1.3  2007/04/24 14:05:00  fnajjar
! Modified and added to EV_PLAG variables
!
! Revision 1.2  2007/04/16 23:21:41  fnajjar
! Deleted all RFLO-pertinent calls and statements
!
! Revision 1.1  2007/04/09 18:50:26  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.1  2007/04/09 18:01:33  haselbac
! Initial revision after split from RocfloMP
!
! Revision 1.27  2006/05/05 17:28:39  haselbac
! Cosmetics only
!
! Revision 1.26  2006/04/07 15:19:23  haselbac
! Removed tabs
!
! Revision 1.25  2005/11/30 22:20:47  fnajjar
! Added EV_PLAG_TEMP
!
! Revision 1.24  2005/05/18 22:15:06  fnajjar
! Changed values of COMM and DELETE status parameters
!
! Revision 1.23  2005/05/02 21:21:23  fnajjar
! Modified parameter definitions in preparation for MPI with RFLU
!
! Revision 1.22  2005/04/27 14:56:56  fnajjar
! Remove parameter for Apte search as obsolete
!
! Revision 1.21  2005/04/25 18:39:10  luca1
! Imposed PDF from file option for random particle ejection
!
! Revision 1.20  2005/03/11 02:24:18  haselbac
! Changed parameters for particle tracking methods
!
! Revision 1.19  2005/01/08 20:41:02  fnajjar
! Added infrastructure for PLAG statistics
!
! Revision 1.18  2005/01/01 21:34:04  haselbac
! Added parameter
!
! Revision 1.17  2004/12/21 15:06:08  fnajjar
! Included definitions for surface statistics datastructure
!
! Revision 1.16  2004/11/06 21:16:08  fnajjar
! Redefined FIND_PCL_METHOD parameters so trajectory be set to 1 as default
!
! Revision 1.15  2004/11/05 21:49:02  fnajjar
! Added parameter entry for trajectory-based search
!
! Revision 1.14  2004/10/10 20:05:42  fnajjar
! Moved PLAG initialization flags to main ModParameters for prep consistency
!
! Revision 1.13  2004/10/09 16:37:10  fnajjar
! Added initialization parameters
!
! Revision 1.12  2004/10/08 22:09:58  haselbac
! Added parameters for FIND_PCL_METHOD_xyz
!
! Revision 1.11  2004/07/28 18:56:18  fnajjar
! Added minimum size for particle datastructure
!
! Revision 1.10  2004/06/17 15:19:03  fnajjar
! Added infrastructure for ejection model
!
! Revision 1.9  2004/06/16 23:00:39  fnajjar
! Renamed PLAG_INJC_MODEL1-2 to PLAG_INJC_LOGNORM-SKWD and TIMEFCTR to COUNTDOWN for CRE kernel
!
! Revision 1.8  2004/03/26 21:25:49  fnajjar
! Split aiv parameters for RFLO and RFLU, added aiv status flag and included 
! note
!
! Revision 1.7  2004/03/02 21:50:03  jferry
! Changed name of DV_PLAG_HTCP to DV_PLAG_SPHT
!
! Revision 1.6  2004/02/13 23:22:07  fnajjar
! Included new cv and aiv definitions for particle burning module
!
! Revision 1.5  2003/09/13 20:14:21  fnajjar
! Added infrastructure for Breakup model
!
! Revision 1.4  2003/09/10 23:35:50  fnajjar
! Removed flags that are subsumed with Rocinteract
!
! Revision 1.3  2003/03/24 23:29:27  jferry
! deleted temporary parameters for material handling
!
! Revision 1.2  2003/03/04 22:12:35  jferry
! Initial import of Rocinteract
!
! Revision 1.1  2002/10/25 14:13:15  f-najjar
! Initial Import of Rocpart
!
! ******************************************************************************

