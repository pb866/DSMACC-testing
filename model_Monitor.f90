! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Utility Data Module File
! 
! Generated by KPP-2.2.4_gc symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : model_Monitor.f90
! Time                 : Sun Feb  5 00:06:35 2017
! Working directory    : /work/home/dp626/DSMACC2
! Equation file        : model.kpp
! Output root filename : model
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE model_Monitor


  CHARACTER(LEN=15), PARAMETER, DIMENSION(5) :: SPC_NAMES = (/ &
     'M              ','NO2            ','O              ', & ! index 1 - 3
     'NO             ','O3             ' /)

  INTEGER, DIMENSION(1) :: LOOKAT
  INTEGER, DIMENSION(1) :: MONITOR
  CHARACTER(LEN=15), DIMENSION(1) :: SMASS
  CHARACTER(LEN=100), PARAMETER, DIMENSION(2) :: EQN_NAMES = (/ &
     'M + O --> M + O3                                                                                    ', & ! index 1
     '  NO2 --> O + NO                                                                                    ' /)

  CHARACTER(LEN=15), PARAMETER, DIMENSION(2) :: EQN_TAGS = (/ &
     '               ','               ' /)

  CHARACTER(LEN=15), DIMENSION(1) :: FAM_NAMES
! INLINED global variables

! End INLINED global variables


END MODULE model_Monitor
