  ! ! FOR COMPILATION
  ! >(GFORTRAN -O3 SET_RK.F90 FS_FC_Z_MOD_RK.F90 FS_FC_Z_DRIVER_RK.F90 -o FS_FC_Z_DRIVER_RK)
  ! ! -WALL
  ! ! OR
  ! >(IFORT -O3 SET_RK.F90 FS_FC_Z_MOD_RK.F90 FS_FC_Z_DRIVER_RK.F90 -O FS_FC_Z_DRIVER_RK)
  ! ! /FPE:0 /TRACEBACK
  ! ! /WARN:ALL
  ! ! /WARN:DECLARATIONS
  ! ! /STAND:F90
  ! ! /CU

  Program FS_FC_Z_DRIVER_RK

    ! .. USE STATEMENTS ..
    Use SET_RK, Only: RK, SP, DP, QP
    Use FS_FC_Z_MOD_RK, Only: FRESNELS_Z, FRESNELS_X, FRESNELC_Z, FRESNELC_X

    Implicit None
    ! .. PARAMETERS ..

    ! For FS real argument
    Integer, Parameter :: MMIN = 1, MMAX = 400001, NREP = 2 ! 338127
    Real (RK), Dimension (MMAX) :: XREF, YREF
    Real (RK), Dimension (MMAX) :: XP, YP, XREF_RK, YREF_RK, REL_ERR, YP_723, REL_ERR_723
    Real (DP), Dimension (MMAX) :: XP723 ! for code 732
    Real (RK) :: MAX_REL_ERR, XMAX, YMAX, MAX_REL_ERR_723, XMAX_723, YMAX_723
    Real (RK), Parameter :: REALMIN = TINY(0.0_RK), EPS_ = EPSILON(1.0_RK)
    Character (Len=45) :: FILENAME
    Integer :: I, N, INDEX_MAX = 0, INDEX_MAX_723 = 0

    ! For FC real argument
    Integer, Parameter :: MMIN2 = 1, MMAX2 = 400001 ! 338127
    Real (RK), Dimension (MMAX2) :: XREF2, YREF2
    Real (RK), Dimension (MMAX2) :: XP2, YP2, XREF2_RK, YREF2_RK, REL_ERR2, YP2_723, REL_ERR2_723
    Real (DP), Dimension (MMAX2) :: XP732_2 ! for code 732
    Real (RK) :: MAX_REL_ERR2, XMAX2, YMAX2, MAX_REL_ERR2_723, XMAX2_723, YMAX2_723
    Character (Len=45) :: FILENAME2
    Integer :: INDEX_MAX2 = 0, INDEX_MAX2_723 = 0

    ! For FS Complex argument
    Integer, Parameter :: NMAX = 39139, NREPEATS = 1, NMIN = 1
    Complex (RK), Parameter :: J1 = (0.0E0_RK, 1.0E0_RK)

    ! .. LOCAL SCALARS ..
    Real (RK) :: TIME_BEGIN, TIME_END ! , value2
    Integer :: JN, JN_ERRMAX
    Real (RK) :: ERR_RE_FS_MAX, ERR_IM_FS_MAX, ERR_FS_MAX_OLD
    ! ..
    ! .. LOCAL ARRAYS ..
    Complex (RK), Dimension (NMAX) :: Z1, FS1
    Real (RK), Dimension (NMAX) :: XTEST, YTEST, RE_FS, IM_FS

    ! For FC Complex argument
    Integer, Parameter :: NMAX2 = 39139, NREPEATS2 = 1, NMIN2 = 1

    ! .. LOCAL SCALARS ..
    Real (RK) :: TIME_BEGIN2, TIME_END2
    Integer :: JN2, JN_ERRMAX2
    Real (RK) :: ERR_RE_FC_MAX, ERR_IM_FC_MAX, ERR_FC_MAX_OLD

    ! .. LOCAL ARRAYS ..
    Complex (RK), Dimension (NMAX2) :: Z1_2, FC1
    Real (RK), Dimension (NMAX2) :: XTEST2, YTEST2, RE_FC, IM_FC

    ! .. INTRINSIC FUNCTIONS ..
    Intrinsic CMPLX, CPU_TIME, REAL




    ! -------------------------------
    ! For FresnelS_x
    ! -------------------------------

   
    FILENAME = 'FresnelS_ref_values2.txt'
    Open (Unit=13, File=FILENAME, Status='old')
    Do I = 1, MMAX
      Read (13, *) XREF(I), YREF(I)
      XREF_RK(I) = XREF(I)
      YREF_RK(I) = YREF(I)
    End Do

     Write (*, *) ' '
    Write (*, *) '*************************'
    Write (*, '(A48,ES15.7E3,A3,ES15.7E3,A8)') ' !!**** Summary of Results From FRESNELS(x) in &
    [', xref(1),       ' , ', xref(mmax), '] ****!!'
    Write (*, *) '*************************'

    XP = XREF

    Call CPU_TIME(TIME_BEGIN)
    Do N = 1, NREP
      YP = 0.0_RK
      Do I = MMIN, MMAX
        YP(I) = FRESNELS_X(XP(I))
      End Do
    End Do
    Call CPU_TIME(TIME_END)
    Write (*, 100) 'Avg. Evaluation time: FresnelS              ', (TIME_END-TIME_BEGIN)/(NREP*(MMAX-MMIN)), &
      '  sec '

    MAX_REL_ERR = 0.0_RK
    Do I = MMIN, MMAX
      If (YP(I)>REALMIN) Then
        REL_ERR(I) = ABS(YP(I)-YREF(I))/ABS(YREF(I))
      End If
      If (REL_ERR(I)>MAX_REL_ERR) Then
        MAX_REL_ERR = REL_ERR(I)
        XMAX = XP(I)
        YMAX = YP(I)
        INDEX_MAX = I
      End If
    End Do

    Write (*, *) '  '
    Write (*, *) 'index of max error                = ', INDEX_MAX
    Write (*, '(A20,ES23.15E3, A24, 2X,ES23.15E3,A24, 2X,ES23.15E3 )') 'x(indx_max_err)  = ', XMAX, &
      'FS_x(indx_max_err)  =', YMAX, 'rel_err(indx_max_err)=', MAX_REL_ERR
    Write (*, *) '  '



    ! -------------------------------
    ! For FresnelS_z
    ! -------------------------------

    Write (*, *) ' '
    Write (*, *) ' '
    Write (*, *) '*************************'
    Write (*, *) 'Summary of Results From FRESNELS(z)'
    Write (*, *) '*************************'

    Open (Unit=12, File='FresnelS308_zref.TXT', Status='OLD')

    Do JN = 1, NMAX
      Read (12, *) XTEST(JN), YTEST(JN), RE_FS(JN), IM_FS(JN)
      Z1(JN) = CMPLX(XTEST(JN), YTEST(JN), KIND=RK)
    End Do

    Call CPU_TIME(TIME_BEGIN)
    Do N = 1, 20
      FS1 = FRESNELS_Z(Z1)
    End Do

    Call CPU_TIME(TIME_END)
    Write (*, 100) 'Avg. Evaluation time:: FresnelS_z              ', (TIME_END-TIME_BEGIN)/(20*NMAX), &
      '  sec '


    ERR_RE_FS_MAX = 0.0_RK
    ERR_IM_FS_MAX = 0.0_RK
    JN_ERRMAX = 1
    Do JN = 1, NMAX
      If (ABS(RE_FS(JN))>0.0_RK .And. ABS(IM_FS(JN))>0.0_RK .And. ABS(RE_FS( &
        JN))<1.0E-2_RK*HUGE(1.0_RK) .And. ABS(IM_FS(JN))<1.0E-2_RK*HUGE(1.0_RK) .And. (ABS( &
        YTEST(JN))/ABS(XTEST(JN)))>2.7E-4_RK+(SP/RK)*21.99973_RK .And. (ABS(YTEST(JN))/ABS( &
        XTEST(JN)))<3.7E3_RK-(SP/RK)*3672.0_RK) Then

        ERR_FS_MAX_OLD = MAX(ERR_RE_FS_MAX, ERR_IM_FS_MAX)
        ERR_RE_FS_MAX = MAX(ERR_RE_FS_MAX, ABS((REAL(FS1(JN),KIND=RK)-RE_FS(JN))/RE_FS(JN)))
        ERR_IM_FS_MAX = MAX(ERR_IM_FS_MAX, ABS((AIMAG(FS1(JN))-IM_FS(JN))/IM_FS(JN)))
        If (MAX(ERR_RE_FS_MAX,ERR_IM_FS_MAX)>ERR_FS_MAX_OLD) Then
          JN_ERRMAX = JN

        End If
      End If
    End Do
    
    Write(*,*)' '
    Write (*, *) 'index of max error  = ', JN_ERRMAX
    Write (*, '(A8,ES23.15E3, A8, 2X,ES23.15E3, A12, 2X, ES23.15E3, A12, 2X, ES23.15E3 )') 'x = ', &
      REAL(Z1(JN_ERRMAX)), 'y =', AIMAG(Z1(JN_ERRMAX)), 'Real(S(z))=', REAL(FS1(JN_ERRMAX)), &
      'imag(S(z))=', AIMAG(FS1(JN_ERRMAX))
    Write (*, '(A24,ES23.15E3, 5X,A24,ES23.15E3)') 'MAX_ERR_RE_FS PRESENT=', ERR_RE_FS_MAX, &
      'MAX_ERR_IM_FS PRESENT=', ERR_IM_FS_MAX

100 Format (A, E14.7, A)
    Close (Unit=12)

    ! -------------------------------
    ! For FresnelC_x
    ! -------------------------------

    FILENAME2 = 'FresnelC_ref_values2.txt'
    Open (Unit=14, File=FILENAME2, Status='old')

    Do I = 1, MMAX2
      Read (14, *) XREF2(I), YREF2(I)
      XREF2_RK(I) = XREF2(I)
      YREF2_RK(I) = YREF2(I)
    End Do

    Write (*, *) ' '
    Write (*, *) '*************************'
    Write (*, '(A48,ES15.7E3,A3,ES15.7E3,A8)') ' !!**** Summary of Results From FRESNELC(x) in &
    [', xref(1),       ' , ', xref(mmax), '] ****!!'
    Write (*, *) '*************************'

    XP2 = XREF2


    Call CPU_TIME(TIME_BEGIN)
    Do N = 1, NREP
      YP2 = 0.0_RK
      Do I = MMIN2, MMAX2
        YP2(I) = FRESNELC_X(XP2(I))
      End Do
    End Do
    Call CPU_TIME(TIME_END)
    Write (*, 100) 'Avg. Evaluation time:: FresnelC_x              ', (TIME_END-TIME_BEGIN)/(NREP*(MMAX2- &
      MMIN2)), '  sec'

    MAX_REL_ERR2 = 0.0_RK
    Do I = MMIN2, MMAX2
      If (YP2(I)>REALMIN) Then
        REL_ERR2(I) = ABS(YP2(I)-YREF2_RK(I))/ABS(YREF2_RK(I))
      End If
      If (REL_ERR2(I)>MAX_REL_ERR2) Then
        MAX_REL_ERR2 = REL_ERR2(I)
        XMAX2 = XP2(I)
        YMAX2 = YP2(I)
        INDEX_MAX2 = I
      End If
    End Do

    Write (*, *) ''
    Write (*, *) 'index of max error                = ', INDEX_MAX2
    Write (*, '(A20,ES23.15E3, A24, 2X,ES23.15E3,A24, 2X,ES23.15E3 )') 'x(indx_max_err)  = ', XMAX2, &
      'FC_x(indx_max_err)  =', YMAX2, 'rel_err(indx_max_err)=', MAX_REL_ERR2
    Write (*, *) '  '



    Close (Unit=14)

    ! -------------------------------
    ! For FresnelC_z
    ! -------------------------------

    Write (*, *) ' '
    Write (*, *) ' '
    Write (*, *) '*************************'
    Write (*, *) 'Summary of Results From FRESNELC(z)'
    Write (*, *) '*************************'

    Open (Unit=15, File='FresnelC308_zref.TXT', Status='OLD')

    Do JN2 = 1, NMAX2
      Read (15, *) XTEST2(JN2), YTEST2(JN2), RE_FC(JN2), IM_FC(JN2)
      Z1_2(JN2) = CMPLX(XTEST2(JN2), YTEST2(JN2), KIND=RK)
    End Do

    Call CPU_TIME(TIME_BEGIN)
    Do N = 1, 20
      FC1 = FRESNELC_Z(Z1_2)
    End Do
    Call CPU_TIME(TIME_END)

    Write (*, 100) 'Avg. Evaluation time:: FresnelC_Z              ', (TIME_END-TIME_BEGIN)/(20*NMAX2), &
      '  sec '

    ERR_RE_FC_MAX = 0.0_RK
    ERR_IM_FC_MAX = 0.0_RK
    JN_ERRMAX2 = 1
    Do JN2 = 1, NMAX2
      If (ABS(RE_FC(JN2))>0.0_RK .And. ABS(IM_FC(JN2))>0.0_RK .And. ABS(RE_FC( &
        JN2))<1.0E-2_RK*HUGE(1.0_RK) .And. ABS(IM_FC(JN2))<1.0E-2_RK*HUGE(1.0_RK) .And. (ABS( &
        YTEST2(JN2))/ABS(XTEST2(JN2)))>3.9E-4_RK+(SP/RK)*0.11961_RK .And. (ABS(YTEST2(JN2))/ABS( &
        XTEST2(JN2)))<2.6E3_RK+(SP/RK)*3.4E38_RK) Then

        ERR_FC_MAX_OLD = MAX(ERR_RE_FC_MAX, ERR_IM_FC_MAX)
        ERR_RE_FC_MAX = MAX(ERR_RE_FC_MAX, ABS((REAL(FC1(JN2),KIND=RK)-RE_FC(JN2))/RE_FC(JN2)))
        ERR_IM_FC_MAX = MAX(ERR_IM_FC_MAX, ABS((AIMAG(FC1(JN2))-IM_FC(JN2))/IM_FC(JN2)))
        If (MAX(ERR_RE_FC_MAX,ERR_IM_FC_MAX)>ERR_FC_MAX_OLD) Then
          JN_ERRMAX2 = JN2

        End If
      End If
    End Do

    Write(*,*)' '
    Write (*, *) 'index of max error  = ', JN_ERRMAX2
    Write (*, '(A8,ES23.15E3, A8, 2X,ES23.15E3, A12, 2X, ES23.15E3, A12, 2X, ES23.15E3 )') 'x = ', &
      REAL(Z1_2(JN_ERRMAX2)), 'y =', AIMAG(Z1_2(JN_ERRMAX2)), 'Real(C(z))=', REAL(FC1(JN_ERRMAX2)), &
      'imag(C(z))=', AIMAG(FC1(JN_ERRMAX2))
    Write (*, '(A24,ES23.15E3, 5X,A24,ES23.15E3)') 'MAX_ERR_RE_FC PRESENT=', ERR_RE_FC_MAX, &
      'MAX_ERR_IM_FC PRESENT=', ERR_IM_FC_MAX

    Close (Unit=15)

  End Program
