MODULE FS_FC_Z_MOD_RK

  !--------------------
  ! This module provides generic interfaces for the Fresnel Sine and Cosine
  ! Integrals for a complex (z=x+iy)  or a real (z=x) argument, in any of the
  ! standard precision arithmetic (single, double or quad) based on
  ! the choice of an integer "rk" in the subsidary module "set_rk" in the
  ! file "set_rk.f90.

  ! For FresnelS function:
  ! S(z)=Int(sin(0.5*pi*t**2) dt, [0, z])
  ! The function is also defined in the form;
  ! S1(z)=sqrt(2/pi)*Int(sin(t**2) dt, [0, z]), or the form
  ! S2(z)=(1/sqrt(2*pi))*Int(sin(t)/sqrt(t) dt, [0, z])
  !       = 0.5*Int(J_(1/2)(t) dt, [0, z]) where J_(1/2) is the ordinary
  ! Bessel function of the first kind and order 1/2
  ! The relation between the above forms of the function can written as
  ! S(z) = S1(sqrt(pi/2)*z) = S2((pi/2)*z**2)
  ! The function is an odd function with the symmetry relation S(-z) =- S(z)
  ! In addition, it has the simple limiting values
  ! S(0)=0, and lim S(z) as z --> +/- Inf = +/- (1/2)

  ! For FresnelC function:
  ! C(z)=Int(cos(0.5*pi*t**2) dt, [0, z])
  ! The function is also defined in the form;
  ! C1(z)=sqrt(2/pi)*Int(cos(t**2) dt, [0, z]), or the form
  ! C2(z)=(1/sqrt(2*pi))*Int(cos(t)/sqrt(t) dt, [0, z])
  !       = 0.5*Int(J_(-1/2)(t) dt, [0, z]) where J_(-1/2) is the ordinary
  ! Bessel function of the first kind and order -1/2
  ! The relation between the above forms of the function can written as
  ! C(z) = C1(sqrt(pi/2)*z) = C2((pi/2)*z**2)
  ! The function is an odd function with the symmetry relation C(-z) = -C(z)
  ! In addition, it has the simple limiting values
  ! C(0)=0, and lim C(z) as z --> +/- Inf = +/- (1/2)

  ! Developed by 
  !Mofreh R. Zaghloul and Leen AlRawas
  ! United Arab Emirates University
  ! March 2023
  !--------------------------------

  USE SET_RK, ONLY: rk, SP, DP, QP
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN,&
    ieee_POSITIVE_INF,ieee_Negative_inf,ieee_POSITIVE_zero,ieee_Negative_zero

  IMPLICIT NONE

  INCLUDE "parameters_FS.f90"
  INCLUDE "parameters_FC.f90"


  PRIVATE
  PUBLIC :: FresnelS_z, FresnelS_x, FresnelC_z, FresnelC_x

CONTAINS

  !-----------------------------------------------
  ELEMENTAL FUNCTION FresnelS_z(z)

    !------------------------------------------------------------------------
    ! FresnelS_z(z) : Fresnel Integral S(z) for a complex (z) input
    ! Method 1: Series for small |z|, sum of half-integer Bessel functions
    ! with the use of Miller's algorithm for intermediate range of |z|
    ! and asymptotic expressions in terms of the auxiliary functions
    ! f(z) and g(z) for large z
    !------------------------------------------------------------------------

    COMPLEX(rk), INTENT(IN) :: z
    COMPLEX(rk):: FresnelS_z, T_sk, u, u_sqr, inv_u, F_k, F_K_max, F_K_max_1
    INTEGER :: k, K_max
    REAL(rk) :: xsqr_plus_ysqr, ax, x_sqr, ay
   
    ax = ABS(REAL(z, KIND = rk))
    ay = ABS(AIMAG(z))

    !--------------------------------------------
    ! pure real or pure imaginary input
    !--------------------------------------------
    IF (ax == zero) THEN
      FresnelS_z = -j1 * FresnelS_x(ay)
    ELSE IF (ay == zero) THEN
      FresnelS_z = FresnelS_x(ax)
    END IF

    x_sqr = ax * ax
    xsqr_plus_ysqr = x_sqr + ay * ay
    u = half_pi*z*z

    !-----------------------------------------------------------------------
    ! For |z|^2 <=6.9 (Series approx.  Abramowitz and Stegun 1964 7.3.13)
    !-----------------------------------------------------------------------
    IF (xsqr_plus_ysqr <= 6.9e0_rk) THEN
      u_sqr = u * u
      FresnelS_z = one_third * z * u
      T_sk = FresnelS_z

      DO k = 1, 60
        T_sk = -(cff_all(k) * u_sqr) * (T_sk)
        FresnelS_z = FresnelS_z + T_sk
        IF (ABS(T_sk) < ABS(FresnelS_z) * eps_trunc) EXIT
      END DO

      !----------------------------------------------------------
      ! For |z|^2> 6.9 and |z|^2 <= (17.0_rk+53.0_rk*(rk/qp))
      ! Using sum of half-integer order Bessel functions together
      ! with Miller's algorithm for back recurrence
      !----------------------------------------------------------
    ELSE IF (xsqr_plus_ysqr <= Miller_border1) THEN
      K_max =  K_max_real_1+K_max_real_1_1*ABS(u)
      inv_u = one / u
      FresnelS_z = (zero, zero)
      F_K_max_1 = (one, zero)
      F_K_max = (one, zero)

      DO  k = K_max, 0, -1
        F_k = (two * k + three) * F_K_max * inv_u - F_K_max_1
        IF (MOD(k,2) /= 0) THEN         !k odd
          FresnelS_z = FresnelS_z + F_k
        END IF
        F_K_max_1 = F_K_max
        F_K_max = F_k
      END DO

      FresnelS_z = FresnelS_z * SIN(u) / (F_k * (half_pi * z))

      !------------------------------------------------------
      ! For rest of the domain
      ! Asymptotic series expansion for S(z) in terms of the
      ! Auxiliary functions f & g  as z-->Inf
      !-----------------------------------------------------

    ELSE IF (xsqr_plus_ysqr <= big_x_border1) THEN
      FresnelS_z = S_Aux_large_z (13, z)
    ELSE IF (xsqr_plus_ysqr <= big_x_border2) THEN
      FresnelS_z = S_Aux_large_z (12, z)
    ELSE IF (xsqr_plus_ysqr <= big_x_border3) THEN
      FresnelS_z = S_Aux_large_z (11, z)
    ELSE IF (xsqr_plus_ysqr <= big_x_border4) THEN
      FresnelS_z = S_Aux_large_z (10, z)
    ELSE IF (xsqr_plus_ysqr <= big_x_border5) THEN
      FresnelS_z = S_Aux_large_z (9, z)
    ELSE IF (xsqr_plus_ysqr <= big_x_border6) THEN
      FresnelS_z = S_Aux_large_z (8, z)
    ELSE IF (xsqr_plus_ysqr <= big_x_border7) THEN
      FresnelS_z = S_Aux_large_z (7, z)
    ELSE
      FresnelS_z = S_Aux_large_z (6, z)
    END IF

  END FUNCTION FresnelS_z

  !-----------------------------------------------
  ! S(z) in terms of Auxiliary functions f & g
  ! as z-->Inf
  !-----------------------------------------------
  ELEMENTAL FUNCTION  S_Aux_large_z (m, z)
    IMPLICIT NONE
    COMPLEX(rk), INTENT(IN):: z
    INTEGER, INTENT (IN) :: m
    COMPLEX(rk):: S_Aux_large_z, u, theta, half_exp_1
    REAL(rk):: ax, ay

    ax = ABS(REAL(z, KIND = rk))
    ay = ABS(AIMAG(z))

    u = z * z
    theta = half_pi * u

    IF (ax < ay) THEN
      half_exp_1=half*EXP(j1*theta)
      S_Aux_large_z = half - (( (one+j1) * half_exp_1) + large_z_f(m,z)) * COS(theta) - &
        (((one-j1) * half_exp_1) + large_z_g(m,z)) * SIN(theta)
    ELSE
      S_Aux_large_z = half - large_z_f(m, z) * COS(theta) - large_z_g(m, z) * SIN(theta)
    ENDIF

    RETURN
  END FUNCTION  S_Aux_large_z

  !--------------------------------------------
  ! Auxiliary function f Expansion for large z
  !--------------------------------------------
  ELEMENTAL FUNCTION  large_z_f (m, z)
    IMPLICIT NONE
    COMPLEX(rk), INTENT(IN):: z
    INTEGER, INTENT (IN):: m
    INTEGER:: k
    COMPLEX(rk):: large_z_f, inv_z, u

    inv_z = one / z
    u = inv_z * inv_z
    u = u * u
    large_z_f = (zero, zero)

    DO k = m, 1, -1
      large_z_f = u * (large_z_f + T_fk(k))
    END DO

    large_z_f = inv_z * (inv_pi + large_z_f)

    RETURN
  END FUNCTION  large_z_f

  !--------------------------------------------
  ! Auxiliary function g Expansion for large z
  !--------------------------------------------
  ELEMENTAL FUNCTION  large_z_g (m, z)
    IMPLICIT NONE
    COMPLEX(rk), INTENT(IN):: z
    INTEGER, INTENT(IN):: m
    INTEGER :: k
    COMPLEX(rk):: large_z_g, inv_z, u

    inv_z = one / z
    u = inv_z * inv_z
    u = u * u
    large_z_g = (zero, zero)

    DO k = m, 1, -1
      large_z_g = u * (large_z_g + T_gk(k))
    END DO

    large_z_g = z * u * (inv_pi_sqr + large_z_g)

    RETURN
  END FUNCTION  large_z_g
  !-----------------------------

  ELEMENTAL FUNCTION FresnelS_x (x)

    !------------------------------------------------------------------------------
    !  FresnelS_x is an elemetal function that evaluates the
    !  Fresnel Sine Integral for the real argument x.
    !  The function FresnelS_x(x) recives "x" as an input and returns "FresnelS_x'
    !  for the Fresnel sine integral using the precision arithmetic
    !  determined by the integer "rk" set in the subsidary module "set_rk".
    !------------------------------------------------------------------------------

    IMPLICIT NONE
    REAL(rk), INTENT(IN)  :: x
    REAL(rk) :: FresnelS_x, ax ,y200, t, T_sk, u, u_sqr
    INTEGER  :: jmin, j, k, ycase

    ax = ABS ( x )
    u = half_pi * ax * ax

    !--------------------------------------------------
    ! Simple converging series for small x ( x<=1.576_rk+1.124 *(rk/qp) )
    !--------------------------------------------------

    u_sqr = u * u
    FresnelS_x = one_third  *ax * u
    T_sk = FresnelS_x
    IF (ax <= small_x_border_S) THEN

      DO k = 1, 40
        T_sk = -(cff_all(k) * u_sqr) * T_sk
        FresnelS_x = FresnelS_x + T_sk
        IF ( ABS(T_sk / FresnelS_x) <  eps_trunc ) EXIT
      END DO

    ELSE IF (ax <= cheb_ax(4)) THEN

      !--------------------------------------------------
      ! Chebyshev subinterval polynomial approximation
      ! for |x|<(12.17407407407407e0_rk - 4.320092440099588e0_rk * (rk/qp)
      !--------------------------------------------------
      DO k = 1, 4
        IF (ax <= cheb_ax(k)) THEN
          y200 = 380.0_rk / (ax + 1.9_rk)
          ycase = y200
          t = two * y200 - (two * ycase + one)
          jmin = N_step(k) + (ycase - N_shift(k)) *shift_factor(k)! (Np_pls_1 + 8 * (k-1)) + 1
          FresnelS_x = zero
          DO j = 1, loop_Nmax(k)
            FresnelS_x = t * (t * (t * (t * (FresnelS_x + &
              cffs(jmin)) + cffs(jmin+1)) + cffs(jmin+2)) + cffs(jmin+3))
            jmin = jmin + 4
          END DO
          FresnelS_x = FresnelS_x + cffs(jmin)
          EXIT
        END IF

      END DO

      !------------------------------------------------------
      ! Asymptotic series expansion for S(z) in terms of the
      ! Auxiliary functions f & g  as z-->Inf
      !------------------------------------------------------
    ELSE

      DO j = 1, 13
        IF (ax <= Aux_large_z_ax(j)) THEN
          FresnelS_x = S_Aux_large_zx ((15-j), ax)
          EXIT
        ELSE IF (ax > Aux_large_z_ax(13)) THEN
          FresnelS_x = S_Aux_large_zx (1, ax)
          EXIT
        END IF
      END DO

    END IF

    IF (x < zero)  THEN
      FresnelS_x = -FresnelS_x
    END IF

    RETURN
  END FUNCTION FresnelS_x

  !--------------------------------------------------------
  ! S(x) in terms of Auxiliary functions f & g  as x-->Inf
  !--------------------------------------------------------
  ELEMENTAL FUNCTION  S_Aux_large_zx (m, z)
    IMPLICIT NONE
    REAL(rk), INTENT(IN) :: z
    INTEGER, INTENT(IN) :: m
    REAL(rk) :: y2, u, t, cos_half_pi_z2, sin_half_pi_z2, S_Aux_large_zx

    u = z * z
    y2 = MOD(u, four)

    IF (y2 == zero) THEN
      cos_half_pi_z2 = one
      sin_half_pi_z2 = zero
    ELSE IF (y2 == two) THEN
      cos_half_pi_z2 = -one
      sin_half_pi_z2 = zero
    ELSE IF (y2 == one) THEN
      cos_half_pi_z2 = zero
      sin_half_pi_z2 = one
    ELSE IF (y2 == three) THEN
      cos_half_pi_z2 = zero
      sin_half_pi_z2 = -one

    ELSE

      t = y2 * half_pi
      sin_half_pi_z2 = SIN(t)
      cos_half_pi_z2 = SQRT(one-sin_half_pi_z2*sin_half_pi_z2)
      IF (t > half_pi .and. t < three*half_pi) THEN
        cos_half_pi_z2 = -cos_half_pi_z2

      END IF

    END IF

    S_Aux_large_zx = half - large_z_f_real(m, z) * cos_half_pi_z2 - &
      large_z_g_real(m, z) * sin_half_pi_z2

    RETURN
  END FUNCTION  S_Aux_large_zx

  !--------------------------------------------
  ! Auxiliary function f Expansion for large x
  !--------------------------------------------
  ELEMENTAL FUNCTION  large_z_f_real (m, z)
    IMPLICIT NONE
    REAL(rk), INTENT(IN):: z
    INTEGER, INTENT(IN):: m
    INTEGER:: k
    REAL(rk):: large_z_f_real, u, inv_z

    inv_z = one / z
    u = inv_z * inv_z
    u = u * u
    large_z_f_real = zero

    DO k = m, 1, -1
      large_z_f_real = u * (large_z_f_real + T_fk(k))
    END DO

    large_z_f_real = inv_z * (inv_pi + large_z_f_real)

    RETURN
  END FUNCTION  large_z_f_real

  !--------------------------------------------
  ! Auxiliary function g Expansion for large x
  !--------------------------------------------
  ELEMENTAL FUNCTION  large_z_g_real (m, z)
    IMPLICIT NONE
    REAL(rk), INTENT(IN):: z
    INTEGER, INTENT (IN):: m
    INTEGER :: k
    REAL(rk):: large_z_g_real, u, inv_z

    inv_z = one / z
    u = inv_z * inv_z
    u = u * u
    large_z_g_real = zero

    DO k = m, 1, -1
      large_z_g_real = u * (large_z_g_real + T_gk(k))
    END DO

    large_z_g_real = z * u * (inv_pi_sqr + large_z_g_real)

    RETURN
  END FUNCTION  large_z_g_real


  !------------
  ! FresnelC
  !------------


  ELEMENTAL FUNCTION FresnelC_z(z)
    !---------------------------------------------------------------------
    ! FresnelC_z(z) : Fresnel Integral C(z) for a complex (z) input
    ! Method 1: Series for small z, sum of half-integer Bessel functions
    ! with the use of Miller's algorithm and asymptotic expressions for
    ! large z
    !---------------------------------------------------------------------
    COMPLEX(rk), INTENT(IN):: z
    COMPLEX(rk):: FresnelC_z, T_ck, u, inv_u, u_sqr, F_k, F_K_max, F_K_max_1
    INTEGER:: k, K_max
    REAL(rk):: xsqr_plus_ysqr, ax, x_sqr, ay

    ax=ABS(REAL(z, KIND = rk))
    ay=ABS(AIMAG(z))

    !-----------------------------------
    ! pure real or pure imaginary input
    !-----------------------------------
    IF (ax == zero) THEN
      FresnelC_z = j1 * FresnelC_x(ay)
    ELSE IF (ay == zero) THEN
      FresnelC_z = FresnelC_x(ax)
    END IF

    x_sqr = ax * ax
    xsqr_plus_ysqr = x_sqr + ay * ay


    u = half_pi * z * z

    !-----------------------------------------------------------------------
    ! For |z|^2 <=12.75 (Series approx.  Abramowitz and Stegun 1964 7.3.13)
    !-----------------------------------------------------------------------
    IF (xsqr_plus_ysqr <= 12.75_rk) THEN
      u_sqr = u * u
      FresnelC_z = z
      T_ck = FresnelC_z

      DO k = 1, 60
        T_ck = -(cff_all2(k) * u_sqr) * (T_ck)
        FresnelC_z = FresnelC_z + T_ck
        IF (ABS(T_ck) < ABS(FresnelC_z) *eps_trunc) EXIT
      END DO

      !-------------------------------------------------------------
      ! For |z|^2> 12.75 and |z|^2 <= (15.65_rk+46.51_rk*(rk/qp))
      ! Using sum of half-integer order Bessel functions together
      ! with Miller's algorithm for back recurrence
      !-------------------------------------------------------------
    ELSE IF (xsqr_plus_ysqr <= Miller_border2) THEN
      K_max =  K_max_real_2+K_max_real_2_2*ABS(u)
      inv_u = one / u
      FresnelC_z = (zero, zero)
      F_K_max_1 = (zero, zero)
      F_K_max = (one, zero)

      DO  k = K_max, 0, -1
        F_k = (two * k + three) * F_K_max * inv_u - F_K_max_1
        IF (MOD(k,2) == 0) THEN              !k even
          FresnelC_z = FresnelC_z + F_k
        END IF
        F_K_max_1 = F_K_max
        F_K_max = F_k
      END DO

      FresnelC_z = FresnelC_z * SIN(u) / F_k / (half * pi * z)

      !------------------------------------------------------
      ! For rest of the domain
      ! Asymptotic series expansion for C(z) in terms of the
      ! Auxiliary functions f & g  as z-->Inf
      !------------------------------------------------------
    ELSE IF (xsqr_plus_ysqr <= big_x_border1_C) THEN
      FresnelC_z = C_Aux_large_z (14, z)
    ELSE IF (xsqr_plus_ysqr <= big_x_border2_C) THEN
      FresnelC_z = C_Aux_large_z (13, z)
    ELSE IF (xsqr_plus_ysqr <= big_x_border3_C) THEN
      FresnelC_z = C_Aux_large_z (11, z)
    ELSE IF (xsqr_plus_ysqr <= big_x_border4_C) THEN
      FresnelC_z = C_Aux_large_z (9, z)
    ELSE IF (xsqr_plus_ysqr <= big_x_border5_C) THEN
      FresnelC_z = C_Aux_large_z (9, z)
    ELSE IF (xsqr_plus_ysqr <= big_x_border6_C) THEN
      FresnelC_z = C_Aux_large_z (8, z)
    ELSE
      FresnelC_z = C_Aux_large_z (7, z)
    END IF

    RETURN
  END FUNCTION FresnelC_z

  !--------------------------------------------
  ! C(z) in terms of Auxiliary functions f & g
  ! as z-->Inf
  !--------------------------------------------
  ELEMENTAL FUNCTION  C_Aux_large_z (m, z)
    IMPLICIT NONE
    COMPLEX(rk), INTENT(IN):: z
    INTEGER, INTENT(IN):: m
    COMPLEX(rk):: C_Aux_large_z, u, theta, half_exp_1
    REAL(rk):: xx, yy, x_sqr, y_sqr
    INTEGER:: mm

    xx = ABS(REAL(z, KIND = rk))
    yy = ABS(AIMAG(z))
    x_sqr = xx * xx
    y_sqr = yy * yy

    u = z * z
    theta = half_pi * u

    IF (xx < yy) THEN
      half_exp_1 =half*EXP(j1*theta)
      C_Aux_large_z = half + (( (one+j1) * half_exp_1 ) + large_z_f(m,z)) * SIN(theta) - &
        (( (one-j1) * half_exp_1) + large_z_g(m, z)) * COS(theta)
    ELSE
      C_Aux_large_z = half + large_z_f(m, z) * SIN(theta) - large_z_g(m, z) * COS(theta)
    ENDIF

    RETURN
  END FUNCTION  C_Aux_large_z
  !--------------------------------

  ELEMENTAL FUNCTION FresnelC_x (x)

    !-------------------------------------------------------------------------------
    !  FresnelC_x is an elemetal function that evaluates the
    !  Fresnel Cosine Integral for a real argument x
    !  The function FresnelC_x(x) recives "x" as an input and returns "FresnelC_x'
    !  for the Fresnel Cosine integral using the precision arithmetic
    !  determined by the integer "rk" set in the subsidary module "set_rk".
    !-------------------------------------------------------------------------------
    IMPLICIT NONE
    REAL(rk), INTENT(IN)  :: x
    REAL(rk):: FresnelC_x, ax,y200, t, T_ck, u
    INTEGER:: jmin, j, i, ycase, k
    REAL(rk):: theta, cos_half_pi_x2, sin_half_pi_x2

    ax = ABS(x)

    !--------------------------------------------------------------------
    ! Simple converging series for small x (1.5760_rk + 1.424_rk * (rk/qp) )
    !--------------------------------------------------------------------
    IF (ax <= small_x_border_C ) THEN

      u = half_pi * ax * ax
      u = u * u
      FresnelC_x = ax
      T_ck = FresnelC_x

      DO k = 1, 40
        T_ck = -(cff_all2(k) * u) * (T_ck)
        FresnelC_x = FresnelC_x + T_ck
        IF (ABS(T_ck/FresnelC_x) <=  eps_trunc) EXIT
      END DO

    ELSE IF (ax <= cheb_ax2(4)) THEN

      !------------------------------------------------
      ! Chebyshev subinterval polynomial approximation
      !------------------------------------------------
      DO i = 1, 4
        IF (ax <= cheb_ax2(i)) THEN
          y200 = 380.0_rk / (ax + 1.9_rk)
          ycase = y200
          t = two * y200 - (two * ycase + one)
          jmin = N_step2(i) + (ycase-N_shift2(i) ) *shift_factor(i)  ! (Np_pls_1 + 8 * (i-1)) + 1
          FresnelC_x = zero
          DO j = 1, loop_Nmax(i)
            FresnelC_x = t * (t * (t * (t * (FresnelC_x + cffs2(jmin)) + cffs2(jmin+1)) + cffs2(jmin+2)) + cffs2(jmin+3))
            jmin = jmin + 4
          END DO
          FresnelC_x = FresnelC_x + cffs2(jmin)
          EXIT
        END IF
      END DO


      !-------------------------------------------------------
      ! Asymptotic series expansion for C(z) in terms of the
      ! Auxiliary functions f & g  as z-->Inf
      !-------------------------------------------------------
    ELSE

      DO j = 1, 15
        IF (ax <= Aux_large_z_ax2(j)) THEN
          FresnelC_x = C_Aux_large_zx ((17-j), ax)
          EXIT
        ELSE IF (ax > Aux_large_z_ax2(15)) THEN
          FresnelC_x = C_Aux_large_zx (1, ax)
          EXIT
        END IF
      END DO
    END IF


    IF (x < zero)  THEN
      FresnelC_x = -FresnelC_x
    ENDIF

    RETURN
  END FUNCTION FresnelC_x

  !---------------------------------------------
  ! C(z) in terms of Auxiliary functions f & g
  ! as z-->Inf
  !---------------------------------------------
  ELEMENTAL FUNCTION  C_Aux_large_zx (m, z)

    IMPLICIT NONE
    REAL(rk), INTENT(IN):: z
    INTEGER, INTENT(IN):: m
    REAL(rk):: C_Aux_large_zx, y2, u, t,cos_half_pi_z2, sin_half_pi_z2

    u = z * z
    y2 = MOD(u, four)
    IF (y2 == zero) THEN
      cos_half_pi_z2 = one
      sin_half_pi_z2 = zero
    ELSE IF (y2 == two) THEN
      cos_half_pi_z2 = -one
      sin_half_pi_z2 = zero
    ELSE IF (y2 == one) THEN
      cos_half_pi_z2 = zero
      sin_half_pi_z2 = one
    ELSE IF (y2 == three) THEN
      cos_half_pi_z2 = zero
      sin_half_pi_z2 = -one
    ELSE

      t = y2 * half_pi
      sin_half_pi_z2 = SIN(t)

      cos_half_pi_z2 = SQRT(one-sin_half_pi_z2*sin_half_pi_z2)
      IF (t > half_pi .and. t < three*half_pi) THEN
        cos_half_pi_z2 = -cos_half_pi_z2

      END IF

    ENDIF

    C_Aux_large_zx = half + (large_z_f_real(m, z) * sin_half_pi_z2 - &
      large_z_g_real(m, z) * cos_half_pi_z2)

    RETURN
  END FUNCTION  C_Aux_large_zx




END MODULE FS_FC_Z_MOD_RK


