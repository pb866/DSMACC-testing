! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Sparse Jacobian Data Structures File
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
! File                 : model_JacobianSP.f90
! Time                 : Wed Jul  5 16:07:46 2017
! Working directory    : /work/home/dp626/DSMACC-testing
! Equation file        : model.kpp
! Output root filename : model
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE model_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(161) :: LU_IROW = (/ &
       1,  2,  2,  3,  3,  3,  4,  4,  4,  5,  5,  6, & ! index 1 - 12
       6,  7,  7,  8,  8,  8,  9,  9,  9, 10, 10, 10, & ! index 13 - 24
      11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 14, 14, & ! index 25 - 36
      14, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 17, & ! index 37 - 48
      17, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 19, & ! index 49 - 60
      20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21, & ! index 61 - 72
      21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 22, 22, & ! index 73 - 84
      22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, & ! index 85 - 96
      24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, & ! index 97 - 108
      24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, & ! index 109 - 120
      25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 26, & ! index 121 - 132
      26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, & ! index 133 - 144
      27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 28, 28, & ! index 145 - 156
      28, 28, 28, 28, 28 /)

  INTEGER, PARAMETER, DIMENSION(161) :: LU_ICOL = (/ &
       1,  1,  2,  3,  9, 17,  4, 12, 24,  2,  5,  6, & ! index 1 - 12
      28,  7,  8,  7,  8, 24,  9, 23, 27, 10, 22, 27, & ! index 13 - 24
      11, 21, 24, 12, 21, 23, 24, 13, 24, 25, 14, 22, & ! index 25 - 36
      24, 15, 24, 26, 27, 16, 18, 19, 22, 23, 25, 17, & ! index 37 - 48
      21, 23, 24, 27, 18, 22, 24, 25, 19, 22, 24, 26, & ! index 49 - 60
      20, 23, 25, 27, 28, 14, 16, 18, 19, 21, 22, 23, & ! index 61 - 72
      24, 25, 26,  7,  8, 10, 19, 22, 23, 24, 25, 26, & ! index 73 - 84
      27,  9, 17, 20, 21, 22, 23, 24, 25, 26, 27, 28, & ! index 85 - 96
       6,  8, 11, 12, 13, 14, 15, 17, 18, 19, 21, 22, & ! index 97 - 108
      23, 24, 25, 26, 27, 28, 13, 20, 22, 23, 24, 25, & ! index 109 - 120
      26, 27, 28,  1, 11, 12, 14, 15, 16, 18, 19, 21, & ! index 121 - 132
      22, 23, 24, 25, 26, 27, 28,  9, 10, 13, 15, 17, & ! index 133 - 144
      18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 20, 23, & ! index 145 - 156
      24, 25, 26, 27, 28 /)

  INTEGER, PARAMETER, DIMENSION(29) :: LU_CROW = (/ &
       1,  2,  4,  7, 10, 12, 14, 16, 19, 22, 25, 28, & ! index 1 - 12
      32, 35, 38, 42, 48, 53, 57, 61, 66, 76, 86, 97, & ! index 13 - 24
     115,124,140,155,162 /)

  INTEGER, PARAMETER, DIMENSION(29) :: LU_DIAG = (/ &
       1,  3,  4,  7, 11, 12, 14, 17, 19, 22, 25, 28, & ! index 1 - 12
      32, 35, 38, 42, 48, 53, 57, 61, 70, 80, 91,110, & ! index 13 - 24
     120,137,153,161,162 /)


END MODULE model_JacobianSP

