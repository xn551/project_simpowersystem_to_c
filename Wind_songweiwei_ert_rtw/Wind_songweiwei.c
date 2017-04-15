/*
 * File: Wind_songweiwei.c
 *
 * Code generated for Simulink model 'Wind_songweiwei'.
 *
 * Model version                  : 1.320
 * Simulink Coder version         : 8.7 (R2014b) 08-Sep-2014
 * C/C++ source code generated on : Thu Mar 30 15:34:34 2017
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "Wind_songweiwei.h"
#include "Wind_songweiwei_private.h"

/* Block signals (auto storage) */
B_Wind_songweiwei_T Wind_songweiwei_B;

/* Block states (auto storage) */
DW_Wind_songweiwei_T Wind_songweiwei_DW;

/* External outputs (root outports fed by signals with auto storage) */
ExtY_Wind_songweiwei_T Wind_songweiwei_Y;

/* Real-time model */
RT_MODEL_Wind_songweiwei_T Wind_songweiwei_M_;
RT_MODEL_Wind_songweiwei_T *const Wind_songweiwei_M = &Wind_songweiwei_M_;

/* Lookup Binary Search Utility BINARYSEARCH_real_T */
void BINARYSEARCH_real_T(uint32_T *piLeft, uint32_T *piRght, real_T u, const
  real_T *pData, uint32_T iHi)
{
  /* Find the location of current input value in the data table. */
  *piLeft = 0U;
  *piRght = iHi;
  if (u <= pData[0] ) {
    /* Less than or equal to the smallest point in the table. */
    *piRght = 0U;
  } else if (u >= pData[iHi] ) {
    /* Greater than or equal to the largest point in the table. */
    *piLeft = iHi;
  } else {
    uint32_T i;

    /* Do a binary search. */
    while (( *piRght - *piLeft ) > 1U ) {
      /* Get the average of the left and right indices using to Floor rounding. */
      i = (*piLeft + *piRght) >> 1;

      /* Move either the right index or the left index so that */
      /*  LeftDataPoint <= CurrentValue < RightDataPoint */
      if (u < pData[i] ) {
        *piRght = i;
      } else {
        *piLeft = i;
      }
    }
  }
}

/* Lookup Utility LookUp_real_T_real_T */
void LookUp_real_T_real_T(real_T *pY, const real_T *pYData, real_T u, const
  real_T *pUData, uint32_T iHi)
{
  uint32_T iLeft;
  uint32_T iRght;
  BINARYSEARCH_real_T( &(iLeft), &(iRght), u, pUData, iHi);

  {
    real_T lambda;
    if (pUData[iRght] > pUData[iLeft] ) {
      real_T num;
      real_T den;
      den = pUData[iRght];
      den = den - pUData[iLeft];
      num = u;
      num = num - pUData[iLeft];
      lambda = num / den;
    } else {
      lambda = 0.0;
    }

    {
      real_T yLeftCast;
      real_T yRghtCast;
      yLeftCast = pYData[iLeft];
      yRghtCast = pYData[iRght];
      yLeftCast += lambda * ( yRghtCast - yLeftCast );
      (*pY) = yLeftCast;
    }
  }
}

/*
 * Start for enable system:
 *    '<S160>/Subsystem - pi//2 delay'
 *    '<S174>/Subsystem - pi//2 delay'
 *    '<S241>/Subsystem - pi//2 delay'
 *    '<S254>/Subsystem - pi//2 delay'
 *    '<S262>/Subsystem - pi//2 delay'
 *    '<S268>/Subsystem - pi//2 delay'
 *    '<S274>/Subsystem - pi//2 delay'
 */
void Wind_so_Subsystempi2delay_Start(B_Subsystempi2delay_Wind_song_T *localB,
  P_Subsystempi2delay_Wind_song_T *localP)
{
  /* VirtualOutportStart for Outport: '<S164>/dq' */
  localB->Fcn = localP->dq_Y0[0];
  localB->Fcn1 = localP->dq_Y0[1];
}

/*
 * Output and update for enable system:
 *    '<S160>/Subsystem - pi//2 delay'
 *    '<S174>/Subsystem - pi//2 delay'
 *    '<S241>/Subsystem - pi//2 delay'
 *    '<S254>/Subsystem - pi//2 delay'
 *    '<S262>/Subsystem - pi//2 delay'
 *    '<S268>/Subsystem - pi//2 delay'
 *    '<S274>/Subsystem - pi//2 delay'
 */
void Wind_songweiw_Subsystempi2delay(uint8_T rtu_Enable, const real_T
  rtu_alpha_beta[2], real_T rtu_wt, B_Subsystempi2delay_Wind_song_T *localB)
{
  /* Outputs for Enabled SubSystem: '<S160>/Subsystem - pi//2 delay' incorporates:
   *  EnablePort: '<S164>/Enable'
   */
  if (rtu_Enable > 0) {
    /* Fcn: '<S164>/Fcn' */
    localB->Fcn = rtu_alpha_beta[0] * sin(rtu_wt) - rtu_alpha_beta[1] * cos
      (rtu_wt);

    /* Fcn: '<S164>/Fcn1' */
    localB->Fcn1 = rtu_alpha_beta[0] * cos(rtu_wt) + rtu_alpha_beta[1] * sin
      (rtu_wt);
  }

  /* End of Outputs for SubSystem: '<S160>/Subsystem - pi//2 delay' */
}

/*
 * Start for enable system:
 *    '<S160>/Subsystem1'
 *    '<S174>/Subsystem1'
 *    '<S241>/Subsystem1'
 *    '<S254>/Subsystem1'
 *    '<S262>/Subsystem1'
 *    '<S268>/Subsystem1'
 *    '<S274>/Subsystem1'
 */
void Wind_songweiwe_Subsystem1_Start(B_Subsystem1_Wind_songweiwei_T *localB,
  P_Subsystem1_Wind_songweiwei_T *localP)
{
  /* VirtualOutportStart for Outport: '<S165>/dq' */
  localB->Fcn = localP->dq_Y0[0];
  localB->Fcn1 = localP->dq_Y0[1];
}

/*
 * Output and update for enable system:
 *    '<S160>/Subsystem1'
 *    '<S174>/Subsystem1'
 *    '<S241>/Subsystem1'
 *    '<S254>/Subsystem1'
 *    '<S262>/Subsystem1'
 *    '<S268>/Subsystem1'
 *    '<S274>/Subsystem1'
 */
void Wind_songweiwei_Subsystem1(uint8_T rtu_Enable, const real_T rtu_alpha_beta
  [2], real_T rtu_wt, B_Subsystem1_Wind_songweiwei_T *localB)
{
  /* Outputs for Enabled SubSystem: '<S160>/Subsystem1' incorporates:
   *  EnablePort: '<S165>/Enable'
   */
  if (rtu_Enable > 0) {
    /* Fcn: '<S165>/Fcn' */
    localB->Fcn = rtu_alpha_beta[0] * cos(rtu_wt) + rtu_alpha_beta[1] * sin
      (rtu_wt);

    /* Fcn: '<S165>/Fcn1' */
    localB->Fcn1 = -rtu_alpha_beta[0] * sin(rtu_wt) + rtu_alpha_beta[1] * cos
      (rtu_wt);
  }

  /* End of Outputs for SubSystem: '<S160>/Subsystem1' */
}

real_T rt_roundd_snf(real_T u)
{
  real_T y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

real_T rt_remd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T u1_0;
  if (!((!rtIsNaN(u0)) && (!rtIsInf(u0)) && ((!rtIsNaN(u1)) && (!rtIsInf(u1)))))
  {
    y = (rtNaN);
  } else {
    if (u1 < 0.0) {
      u1_0 = ceil(u1);
    } else {
      u1_0 = floor(u1);
    }

    if ((u1 != 0.0) && (u1 != u1_0)) {
      u1_0 = u0 / u1;
      if (fabs(u1_0 - rt_roundd_snf(u1_0)) <= DBL_EPSILON * fabs(u1_0)) {
        y = 0.0;
      } else {
        y = fmod(u0, u1);
      }
    } else {
      y = fmod(u0, u1);
    }
  }

  return y;
}

void rt_invd5x5_snf(const real_T u[25], real_T y[25])
{
  int8_T p[5];
  int32_T c;
  real_T A[25];
  int8_T ipiv[5];
  int32_T j;
  int32_T ix;
  real_T smax;
  real_T s;
  int32_T jA;
  int32_T jBcol;
  int32_T kAcol;
  int32_T i;
  int32_T i_0;
  for (i_0 = 0; i_0 < 25; i_0++) {
    y[i_0] = 0.0;
    A[i_0] = u[i_0];
  }

  for (i_0 = 0; i_0 < 5; i_0++) {
    ipiv[i_0] = (int8_T)(1 + i_0);
  }

  for (j = 0; j < 4; j++) {
    c = j * 6;
    kAcol = j * 6;
    jBcol = 1;
    ix = c;
    smax = fabs(A[kAcol]);
    for (jA = 2; jA <= 5 - j; jA++) {
      ix++;
      s = fabs(A[ix]);
      if (s > smax) {
        jBcol = jA;
        smax = s;
      }
    }

    jBcol--;
    if (A[c + jBcol] != 0.0) {
      if (jBcol != 0) {
        ipiv[j] = (int8_T)((j + jBcol) + 1);
        ix = j;
        jBcol += j;
        for (jA = 0; jA < 5; jA++) {
          smax = A[ix];
          A[ix] = A[jBcol];
          A[jBcol] = smax;
          ix += 5;
          jBcol += 5;
        }
      }

      i_0 = c - j;
      for (i = c + 1; i + 1 <= i_0 + 5; i++) {
        A[i] /= A[c];
      }
    }

    jA = kAcol;
    jBcol = kAcol + 5;
    for (kAcol = 1; kAcol <= 4 - j; kAcol++) {
      if (A[jBcol] != 0.0) {
        smax = -A[jBcol];
        ix = c;
        i_0 = jA - j;
        for (i = jA + 6; i + 1 <= i_0 + 10; i++) {
          A[i] += A[ix + 1] * smax;
          ix++;
        }
      }

      jBcol += 5;
      jA += 5;
    }
  }

  for (i_0 = 0; i_0 < 5; i_0++) {
    p[i_0] = (int8_T)(1 + i_0);
  }

  if (ipiv[0] > 1) {
    jBcol = p[ipiv[0] - 1];
    p[ipiv[0] - 1] = p[0];
    p[0] = (int8_T)jBcol;
  }

  if (ipiv[1] > 2) {
    jBcol = p[ipiv[1] - 1];
    p[ipiv[1] - 1] = p[1];
    p[1] = (int8_T)jBcol;
  }

  if (ipiv[2] > 3) {
    jBcol = p[ipiv[2] - 1];
    p[ipiv[2] - 1] = p[2];
    p[2] = (int8_T)jBcol;
  }

  if (ipiv[3] > 4) {
    jBcol = p[ipiv[3] - 1];
    p[ipiv[3] - 1] = p[3];
    p[3] = (int8_T)jBcol;
  }

  for (jA = 0; jA < 5; jA++) {
    c = p[jA] - 1;
    y[jA + 5 * (p[jA] - 1)] = 1.0;
    for (j = jA; j + 1 < 6; j++) {
      if (y[5 * c + j] != 0.0) {
        for (i = j + 1; i + 1 < 6; i++) {
          y[i + 5 * c] -= y[5 * c + j] * A[5 * j + i];
        }
      }
    }
  }

  for (j = 0; j < 5; j++) {
    jBcol = 5 * j;
    for (jA = 4; jA >= 0; jA += -1) {
      kAcol = 5 * jA;
      if (y[jA + jBcol] != 0.0) {
        y[jA + jBcol] /= A[jA + kAcol];
        for (i = 0; i + 1 <= jA; i++) {
          y[i + jBcol] -= y[jA + jBcol] * A[i + kAcol];
        }
      }
    }
  }
}

real_T rt_modd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  if (u1 == 0.0) {
    y = u0;
  } else if (!((!rtIsNaN(u0)) && (!rtIsInf(u0)) && ((!rtIsNaN(u1)) && (!rtIsInf
                (u1))))) {
    y = (rtNaN);
  } else {
    tmp = u0 / u1;
    if (u1 <= floor(u1)) {
      y = u0 - floor(tmp) * u1;
    } else if (fabs(tmp - rt_roundd_snf(tmp)) <= DBL_EPSILON * fabs(tmp)) {
      y = 0.0;
    } else {
      y = (tmp - floor(tmp)) * u1;
    }
  }

  return y;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  real_T tmp_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = (rtNaN);
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T a;
  a = fabs(u0);
  y = fabs(u1);
  if (a < y) {
    a /= y;
    y *= sqrt(a * a + 1.0);
  } else if (a > y) {
    y /= a;
    y = sqrt(y * y + 1.0) * a;
  } else {
    if (!rtIsNaN(y)) {
      y = a * 1.4142135623730951;
    }
  }

  return y;
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  int32_T u0_0;
  int32_T u1_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      u0_0 = 1;
    } else {
      u0_0 = -1;
    }

    if (u1 > 0.0) {
      u1_0 = 1;
    } else {
      u1_0 = -1;
    }

    y = atan2(u0_0, u1_0);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

/* Model step function */
void Wind_songweiwei_step(void)
{
  /* local block i/o variables */
  real_T rtb_ib1;
  real_T rtb_LookupTable;
  real_T rtb_dw_delay;
  real_T rtb_Saturation[6];
  real_T rtb_MathFunction;
  real_T rtb_Gain1[3];
  real_T rtb_Kp5;
  real_T rtb_Saturation2;
  real_T rtb_Gain1_k[3];
  real_T rtb_MathFunction_f;
  real_T rtb_Gain1_c[3];
  real_T rtb_Switch[2];
  real_T rtb_Switch_e;
  real_T rtb_Switch_d;
  real_T rtb_Vpu;
  real_T rtb_Sum;
  real_T rtb_MinMax;
  real_T rtb_Gain1_h[3];
  real_T rtb_Kp5_p;
  real_T rtb_Kp5_a[2];
  real_T rtb_xk1_k;
  real_T rtb_Gain1_i[3];
  real_T rtb_Switch_o[2];
  real_T rtb_Switch_n;
  real_T rtb_Switch_h;
  real_T rtb_varpu;
  real_T rtb_Sum2_g;
  real_T rtb_Gain3_h[3];
  real_T rtb_x1k1[3];
  real_T rtb_x2k1[3];
  real_T rtb_x1k1_n[3];
  real_T rtb_x2k1_n[3];
  real_T rtb_x1k1_a[3];
  real_T rtb_x2k1_a[3];
  real_T rtb_x1k1_p[3];
  real_T rtb_x2k1_i[3];
  real_T rtb_x1k1_ag;
  real_T rtb_x2k1_il;
  real_T rtb_x1k1_al;
  real_T rtb_x2k1_k;
  real_T rtb_x1k1_o;
  real_T rtb_x2k1_ac;
  real_T rtb_x1k1_f;
  real_T rtb_x2k1_l;
  real_T rtb_x1k1_d;
  real_T rtb_x2k1_h;
  real_T rtb_Gain1_m[3];
  real_T rtb_Switch_di[2];
  real_T rtb_Switch_g;
  real_T rtb_Divide;
  real_T rtb_Saturation1;
  real_T rtb_x1k1_oi;
  real_T rtb_x2k1_hi;
  real_T rtb_yk_c;
  real_T rtb_xk1_g;
  real_T rtb_Wpu;
  real_T rtb_Kp5_i;
  real_T rtb_Kp5_l;
  real_T rtb_Kp5_aq;
  real_T rtb_Saturation2_e;
  real_T rtb_xk1_a;
  real_T rtb_yk_a3;
  real_T rtb_xk1_go;
  real_T rtb__2H_WT;
  real_T rtb_wbase;
  real_T rtb_Vfd;
  real_T rtb_Fcn2;
  real_T rtb_Fcn3;
  real_T rtb_IC[5];
  real_T rtb_uH;
  real_T rtb_webase;
  real_T rtb_web3;
  real_T rtb_Gain1_dv[3];
  real_T rtb_Switch_m[2];
  real_T rtb_ElementaryMath;
  real_T rtb_ComplextoMagnitudeAngle_o1;
  real_T rtb_Switch4;
  real_T rtb_ElementaryMath1;
  real_T rtb_Sum1_jb;
  real_T rtb_donotdeletethisgain_l;
  real_T rtb_donotdeletethisgain_o;
  uint8_T rtb_Compare;
  uint8_T rtb_Compare_a;
  uint8_T rtb_Compare_m;
  uint8_T rtb_Compare_g;
  uint8_T rtb_Compare_b;
  uint8_T rtb_Compare_p;
  uint8_T rtb_Compare_o;
  uint8_T rtb_Compare_j;
  uint8_T rtb_Compare_f;
  uint8_T rtb_Compare_k;
  uint8_T rtb_Compare_b0;
  uint8_T rtb_Compare_c;
  uint8_T rtb_Compare_kz;
  uint8_T rtb_Compare_od;

  {
    real_T rtb_Add4;
    real_T rtb_changeIqIdcurrentsigns[5];
    real_T rtb_pu_a;
    real_T rtb_units;
    real_T rtb_n;
    real_T rtb_degrd;
    real_T rtb_yk;
    real_T rtb_Cp;
    real_T rtb_donotdeletethisgain;
    uint8_T rtb_Compare_iw;
    real_T rtb_phimd;
    real_T rtb_Switch1_p;
    real_T rtb_Ifdsat;
    real_T rtb_Lmsatd;
    real_T rtb_phimq;
    real_T rtb_Switch_e3;
    real_T rtb_Sum2_d[25];
    real_T rtb_inversion[25];
    int32_T i;
    real_T rtb_Sum2_d_0[25];
    real_T tmp;
    real_T tmp_0[5];
    real_T tmp_1[3];
    int32_T i_0;
    real_T tmp_2[5];
    real_T rtb_inversion_0[5];
    int32_T i_1;
    real_T rtb_Ll_q_idx_1;
    real_T rtb_Ll_q_idx_0;
    real_T rtb_Ll_d_idx_2;
    real_T rtb_Ll_d_idx_1;
    real_T rtb_Ll_d_idx_0;
    real_T rtb_Kv_o_idx_2;
    real_T rtb_Kv_o_idx_1;
    real_T rtb_Kv_o_idx_0;
    real_T rtb_Kv1_idx_2;
    real_T rtb_Kv1_idx_1;
    real_T rtb_Kv1_idx_0;
    real_T rtb_Kv_idx_2;
    real_T rtb_Kv_idx_1;
    real_T rtb_Kv_idx_0;
    real_T rtb_Kv1_j_idx_2;
    real_T rtb_Kv1_j_idx_1;
    real_T rtb_Kv1_j_idx_0;
    real_T rtb_Kv_ow_idx_2;
    real_T rtb_Kv_ow_idx_1;
    real_T unnamed_idx_1;
    real_T rtb_Gain2_m_idx_0;
    boolean_T rtb_RelationalOperator2_idx_2;
    boolean_T rtb_RelationalOperator2_idx_1;
    boolean_T rtb_RelationalOperator2_idx_0;

    /* DigitalClock: '<S302>/Digital Clock' */
    rtb_ElementaryMath1 = (((Wind_songweiwei_M->Timing.clockTick1+
      Wind_songweiwei_M->Timing.clockTickH1* 4294967296.0)) * 2.0E-6);

    /* Sum: '<S302>/Add1' incorporates:
     *  Constant: '<S302>/Constant3'
     */
    rtb_ElementaryMath1 += Wind_songweiwei_P.Constant3_Value_e;

    /* Math: '<S302>/Math Function' incorporates:
     *  Constant: '<S302>/Constant1'
     */
    rtb_ElementaryMath1 = rt_remd_snf(rtb_ElementaryMath1,
      Wind_songweiwei_P.Constant1_Value_m);

    /* Gain: '<S302>/1\ib1' */
    rtb_ib1 = Wind_songweiwei_P.ib1_Gain * rtb_ElementaryMath1;

    /* Lookup: '<S302>/Lookup Table'
     * About '<S302>/Lookup Table':
     * Input0  Data Type:  Floating Point real_T
     * Output0 Data Type:  Floating Point real_T
     * Lookup Method: Linear_Endpoint
     *
     * XData parameter uses the same data type and scaling as Input0
     * YData parameter uses the same data type and scaling as Output0
     */
    LookUp_real_T_real_T( &(rtb_LookupTable),
                         Wind_songweiwei_P.LookupTable_YData_k, rtb_ib1,
                         Wind_songweiwei_P.LookupTable_XData_f, 2U);

    /* Gain: '<S280>/Gain1' incorporates:
     *  Constant: '<S123>/Constant10'
     *  Sum: '<S280>/Add3'
     */
    rtb_ElementaryMath1 = (Wind_songweiwei_P.PWMGenerator2Level_MinMax[1] -
      Wind_songweiwei_P.PWMGenerator2Level_MinMax[0]) *
      Wind_songweiwei_P.Gain1_Gain_d;

    /* Sum: '<S280>/Add4' incorporates:
     *  Constant: '<S123>/Constant10'
     *  Constant: '<S302>/Constant2'
     *  Gain: '<S302>/1\ib1'
     *  Lookup: '<S302>/Lookup Table'
     *  Product: '<S280>/MUL1'
     *  Sum: '<S302>/Add3'
     */
    rtb_Add4 = ((rtb_LookupTable - Wind_songweiwei_P.Constant2_Value_a) *
                rtb_ElementaryMath1 +
                Wind_songweiwei_P.PWMGenerator2Level_MinMax[0]) +
      rtb_ElementaryMath1;

    /* RelationalOperator: '<S286>/Relational Operator2' incorporates:
     *  UnitDelay: '<S67>/Unit Delay6'
     */
    rtb_RelationalOperator2_idx_0 = (Wind_songweiwei_DW.UnitDelay6_DSTATE[0] >=
      rtb_Add4);
    rtb_RelationalOperator2_idx_1 = (Wind_songweiwei_DW.UnitDelay6_DSTATE[1] >=
      rtb_Add4);
    rtb_RelationalOperator2_idx_2 = (Wind_songweiwei_DW.UnitDelay6_DSTATE[2] >=
      rtb_Add4);

    /* DataTypeConversion: '<S123>/Data Type Conversion' incorporates:
     *  Logic: '<S123>/Logical Operator4'
     */
    Wind_songweiwei_B.DataTypeConversion[0] = rtb_RelationalOperator2_idx_0;
    Wind_songweiwei_B.DataTypeConversion[1] = !rtb_RelationalOperator2_idx_0;
    Wind_songweiwei_B.DataTypeConversion[2] = rtb_RelationalOperator2_idx_1;
    Wind_songweiwei_B.DataTypeConversion[3] = !rtb_RelationalOperator2_idx_1;
    Wind_songweiwei_B.DataTypeConversion[4] = rtb_RelationalOperator2_idx_2;
    Wind_songweiwei_B.DataTypeConversion[5] = !rtb_RelationalOperator2_idx_2;

    /* Outputs for Enabled SubSystem: '<S353>/Tail' incorporates:
     *  EnablePort: '<S354>/Enable'
     */
    /* Constant: '<S353>/2' */
    if (Wind_songweiwei_P._Value_e) {
      if (!Wind_songweiwei_DW.Tail_MODE) {
        for (i = 0; i < 6; i++) {
          /* InitializeConditions for DiscreteIntegrator: '<S354>/Discrete-Time Integrator' */
          Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_b4[i] =
            Wind_songweiwei_P.DiscreteTimeIntegrator_IC_j;
          Wind_songweiwei_DW.DiscreteTimeIntegrator_PrevRese[i] = 2;

          /* InitializeConditions for UnitDelay: '<S354>/Unit Delay' */
          Wind_songweiwei_DW.UnitDelay_DSTATE_fv[i] =
            Wind_songweiwei_P.UnitDelay_InitialCondition_nz;
        }

        Wind_songweiwei_DW.Tail_MODE = true;
      }

      /* Constant: '<S354>/2' */
      rtb_RelationalOperator2_idx_0 = (Wind_songweiwei_P.Tail_Tf +
        Wind_songweiwei_P.Tail_Tt > 0.0);

      /* Gain: '<S356>/-0.9//Tf' */
      rtb_Ifdsat = -0.9 / (Wind_songweiwei_P.Tail_Tf + 2.2204460492503131E-16);

      /* Sum: '<S356>/Add1' incorporates:
       *  Constant: '<S356>/Constant2'
       */
      rtb_Add4 = Wind_songweiwei_P.Tail_Tf + Wind_songweiwei_P.Tail_Tt;

      /* Gain: '<S356>/0.1//Tt' */
      rtb_Lmsatd = 0.1 / (Wind_songweiwei_P.Tail_Tt + 2.2204460492503131E-16);

      /* Update for DiscreteIntegrator: '<S354>/Discrete-Time Integrator' incorporates:
       *  Constant: '<S354>/1'
       */
      tmp = Wind_songweiwei_P.DiscreteTimeIntegrator_gainva_c *
        Wind_songweiwei_P._Value;
      for (i = 0; i < 6; i++) {
        /* DiscreteIntegrator: '<S354>/Discrete-Time Integrator' */
        if ((Wind_songweiwei_B.DataTypeConversion[i] <= 0.0) &&
            (Wind_songweiwei_DW.DiscreteTimeIntegrator_PrevRese[i] == 1)) {
          Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_b4[i] =
            Wind_songweiwei_P.DiscreteTimeIntegrator_IC_j;
        }

        /* Switch: '<S354>/Switch' incorporates:
         *  UnitDelay: '<S353>/Unit Delay'
         *  UnitDelay: '<S354>/Unit Delay'
         */
        if (Wind_songweiwei_B.DataTypeConversion[i] >=
            Wind_songweiwei_P.Switch_Threshold_b) {
          rtb_Kv_ow_idx_2 = Wind_songweiwei_DW.UnitDelay_DSTATE[i];
        } else {
          rtb_Kv_ow_idx_2 = Wind_songweiwei_DW.UnitDelay_DSTATE_fv[i];
        }

        /* End of Switch: '<S354>/Switch' */

        /* Gain: '<S356>/0.1//Tt' incorporates:
         *  Constant: '<S356>/Constant2'
         *  DiscreteIntegrator: '<S354>/Discrete-Time Integrator'
         *  Sum: '<S356>/Add1'
         */
        unnamed_idx_1 = (rtb_Add4 -
                         Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_b4[i]) *
          rtb_Lmsatd;

        /* Sum: '<S356>/Add' incorporates:
         *  Constant: '<S356>/Constant'
         *  DiscreteIntegrator: '<S354>/Discrete-Time Integrator'
         *  Gain: '<S356>/-0.9//Tf'
         */
        rtb_Switch_e3 = rtb_Ifdsat *
          Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_b4[i] +
          Wind_songweiwei_P.Constant_Value_c;

        /* Saturate: '<S356>/Saturation1' */
        if (rtb_Switch_e3 > Wind_songweiwei_P.Saturation1_UpperSat) {
          rtb_Switch_e3 = Wind_songweiwei_P.Saturation1_UpperSat;
        } else {
          if (rtb_Switch_e3 < Wind_songweiwei_P.Saturation1_LowerSat) {
            rtb_Switch_e3 = Wind_songweiwei_P.Saturation1_LowerSat;
          }
        }

        /* Saturate: '<S356>/Saturation2' */
        if (unnamed_idx_1 > Wind_songweiwei_P.Saturation2_UpperSat) {
          unnamed_idx_1 = Wind_songweiwei_P.Saturation2_UpperSat;
        } else {
          if (unnamed_idx_1 < Wind_songweiwei_P.Saturation2_LowerSat) {
            unnamed_idx_1 = Wind_songweiwei_P.Saturation2_LowerSat;
          }
        }

        /* Product: '<S354>/Product' incorporates:
         *  Constant: '<S354>/2'
         *  Saturate: '<S356>/Saturation1'
         *  Saturate: '<S356>/Saturation2'
         *  Sum: '<S356>/Add2'
         */
        Wind_songweiwei_B.Product[i] = (rtb_Switch_e3 + unnamed_idx_1) *
          rtb_Kv_ow_idx_2 * (real_T)rtb_RelationalOperator2_idx_0;

        /* Update for DiscreteIntegrator: '<S354>/Discrete-Time Integrator' */
        Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_b4[i] += tmp;
        if (Wind_songweiwei_B.DataTypeConversion[i] > 0.0) {
          Wind_songweiwei_DW.DiscreteTimeIntegrator_PrevRese[i] = 1;
        } else if (Wind_songweiwei_B.DataTypeConversion[i] < 0.0) {
          Wind_songweiwei_DW.DiscreteTimeIntegrator_PrevRese[i] = -1;
        } else if (Wind_songweiwei_B.DataTypeConversion[i] == 0.0) {
          Wind_songweiwei_DW.DiscreteTimeIntegrator_PrevRese[i] = 0;
        } else {
          Wind_songweiwei_DW.DiscreteTimeIntegrator_PrevRese[i] = 2;
        }

        /* Update for UnitDelay: '<S354>/Unit Delay' */
        Wind_songweiwei_DW.UnitDelay_DSTATE_fv[i] = rtb_Kv_ow_idx_2;
      }
    } else {
      if (Wind_songweiwei_DW.Tail_MODE) {
        Wind_songweiwei_DW.Tail_MODE = false;
      }
    }

    /* End of Constant: '<S353>/2' */
    /* End of Outputs for SubSystem: '<S353>/Tail' */

    /* DiscreteIntegrator: '<S329>/Rotor angle dthetae' */
    rtb_ElementaryMath1 = Wind_songweiwei_DW.Rotorangledthetae_DSTATE;

    /* DigitalClock: '<S329>/Digital Clock' */
    rtb_ElementaryMath = (((Wind_songweiwei_M->Timing.clockTick1+
      Wind_songweiwei_M->Timing.clockTickH1* 4294967296.0)) * 2.0E-6);

    /* Sum: '<S329>/Sum3' incorporates:
     *  Gain: '<S329>/web2'
     */
    rtb_ElementaryMath1 += Wind_songweiwei_P.web2_Gain * rtb_ElementaryMath;

    /* Trigonometry: '<S332>/Elementary Math' */
    rtb_ElementaryMath = sin(rtb_ElementaryMath1);

    /* Trigonometry: '<S332>/Elementary Math1' */
    rtb_ElementaryMath1 = cos(rtb_ElementaryMath1);

    /* Outputs for Enabled SubSystem: '<S331>/Saturation' incorporates:
     *  EnablePort: '<S341>/Enable'
     */
    /* Constant: '<S331>/Constant1' */
    if (Wind_songweiwei_P.Constant1_Value_f > 0.0) {
      /* Abs: '<S343>/Abs' incorporates:
       *  Constant: '<S346>/u1'
       *  Gain: '<S347>/1//Ll_d'
       *  Math: '<S346>/Math Function2'
       *  Math: '<S346>/Math Function3'
       *  Product: '<S347>/Product4'
       *  Sum: '<S346>/Sum1'
       *  Sum: '<S347>/sum2'
       *  UnitDelay: '<S342>/fluxes'
       *  UnitDelay: '<S343>/Lmd_sat'
       *
       * About '<S346>/Math Function2':
       *  Operator: reciprocal
       *
       * About '<S346>/Math Function3':
       *  Operator: reciprocal
       */
      rtb_Add4 = fabs(1.0 / (((Wind_songweiwei_P.u1_Value[0] +
        Wind_songweiwei_P.u1_Value[1]) + Wind_songweiwei_P.u1_Value[2]) + 1.0 /
        Wind_songweiwei_DW.Lmd_sat_DSTATE) * ((Wind_songweiwei_P.Ll_d_Gain[0] *
        Wind_songweiwei_DW.fluxes_DSTATE[1] + Wind_songweiwei_P.Ll_d_Gain[1] *
        Wind_songweiwei_DW.fluxes_DSTATE[2]) + Wind_songweiwei_P.Ll_d_Gain[2] *
        Wind_songweiwei_DW.fluxes_DSTATE[3]));

      /* Lookup: '<S343>/Lookup Table' */
      rtb_Ifdsat = rt_Lookup(Wind_songweiwei_P.LookupTable_XData_k, 2, rtb_Add4,
        Wind_songweiwei_P.LookupTable_YData_j);

      /* Switch: '<S343>/Switch' incorporates:
       *  Constant: '<S343>/Constant1'
       *  Product: '<S343>/Product'
       */
      if (rtb_Ifdsat != 0.0) {
        rtb_Switch1_p = rtb_Add4 / rtb_Ifdsat;
      } else {
        rtb_Switch1_p = Wind_songweiwei_P.Constant1_Value_n;
      }

      /* End of Switch: '<S343>/Switch' */

      /* Gain: '<S343>/Lmd' */
      rtb_Lmsatd = Wind_songweiwei_P.Lmd_Gain * rtb_Switch1_p;

      /* Outputs for Enabled SubSystem: '<S341>/Lmq_sat' incorporates:
       *  EnablePort: '<S344>/Enable'
       */
      /* Constant: '<S341>/Constant1' */
      if (Wind_songweiwei_P.Constant1_Value_k) {
        /* Abs: '<S344>/Abs' incorporates:
         *  Constant: '<S348>/u2'
         *  Gain: '<S349>/1//Ll_q'
         *  Math: '<S348>/Math Function'
         *  Math: '<S348>/Math Function1'
         *  Product: '<S349>/Product3'
         *  Sum: '<S348>/Sum2'
         *  Sum: '<S349>/sum1'
         *  UnitDelay: '<S342>/fluxes'
         *  UnitDelay: '<S344>/Lmq_sat'
         *
         * About '<S348>/Math Function':
         *  Operator: reciprocal
         *
         * About '<S348>/Math Function1':
         *  Operator: reciprocal
         */
        rtb_Add4 = fabs(1.0 / (((Wind_songweiwei_P.u2_Value[0] +
          Wind_songweiwei_P.u2_Value[1]) + Wind_songweiwei_P.u2_Value[2]) + 1.0 /
          Wind_songweiwei_DW.Lmq_sat_DSTATE) * (Wind_songweiwei_P.Ll_q_Gain[0] *
          Wind_songweiwei_DW.fluxes_DSTATE[0] + Wind_songweiwei_P.Ll_q_Gain[1] *
          Wind_songweiwei_DW.fluxes_DSTATE[4]));

        /* Lookup: '<S344>/Lookup Table' */
        rtb_Ifdsat = rt_Lookup(Wind_songweiwei_P.LookupTable_XData, 2, rtb_Add4,
          Wind_songweiwei_P.LookupTable_YData);

        /* Switch: '<S344>/Switch' incorporates:
         *  Constant: '<S344>/Constant1'
         *  Product: '<S344>/Product'
         */
        if (rtb_Ifdsat != 0.0) {
          rtb_Switch_e3 = rtb_Add4 / rtb_Ifdsat;
        } else {
          rtb_Switch_e3 = Wind_songweiwei_P.Constant1_Value_e;
        }

        /* End of Switch: '<S344>/Switch' */

        /* Gain: '<S344>/Lmq' */
        Wind_songweiwei_B.Lmsatq = Wind_songweiwei_P.Lmq_Gain * rtb_Switch_e3;

        /* Update for UnitDelay: '<S344>/Lmq_sat' */
        Wind_songweiwei_DW.Lmq_sat_DSTATE = Wind_songweiwei_B.Lmsatq;
      }

      /* End of Constant: '<S341>/Constant1' */
      /* End of Outputs for SubSystem: '<S341>/Lmq_sat' */

      /* Switch: '<S341>/Switch1' incorporates:
       *  Constant: '<S341>/Constant2'
       *  Constant: '<S341>/u3'
       */
      if (Wind_songweiwei_P.Constant2_Value_j) {
        rtb_Switch1_p = Wind_songweiwei_B.Lmsatq;
      } else {
        rtb_Switch1_p = Wind_songweiwei_P.u3_Value;
      }

      /* End of Switch: '<S341>/Switch1' */

      /* Assignment: '<S345>/Update Lmq rows[1,5, (6)] col[1,5, (6)] ' incorporates:
       *  Constant: '<S345>/u1'
       */
      memcpy(&rtb_Sum2_d[0], &Wind_songweiwei_P.u1_Value_m[0], 25U * sizeof
             (real_T));
      rtb_Sum2_d[0] = rtb_Switch1_p;
      rtb_Sum2_d[4] = rtb_Switch1_p;
      rtb_Sum2_d[20] = rtb_Switch1_p;
      rtb_Sum2_d[24] = rtb_Switch1_p;

      /* Assignment: '<S345>/Update Lmd rows[2,3,4] col[2,3,4]' */
      for (i = 0; i < 3; i++) {
        rtb_Sum2_d[1 + 5 * (1 + i)] = rtb_Lmsatd;
        rtb_Sum2_d[2 + 5 * (1 + i)] = rtb_Lmsatd;
        rtb_Sum2_d[3 + 5 * (1 + i)] = rtb_Lmsatd;
      }

      /* End of Assignment: '<S345>/Update Lmd rows[2,3,4] col[2,3,4]' */

      /* Sum: '<S345>/Sum2' incorporates:
       *  Constant: '<S345>/u5'
       */
      for (i = 0; i < 25; i++) {
        rtb_Sum2_d_0[i] = rtb_Sum2_d[i] + Wind_songweiwei_P.u5_Value[i];
      }

      /* End of Sum: '<S345>/Sum2' */

      /* Product: '<S341>/inversion' */
      rt_invd5x5_snf(rtb_Sum2_d_0, Wind_songweiwei_B.Linv);

      /* Product: '<S341>/Product1' incorporates:
       *  Constant: '<S341>/u1'
       */
      for (i = 0; i < 5; i++) {
        for (i_1 = 0; i_1 < 5; i_1++) {
          Wind_songweiwei_B.RLinv[i + 5 * i_1] = 0.0;
          for (i_0 = 0; i_0 < 5; i_0++) {
            Wind_songweiwei_B.RLinv[i + 5 * i_1] +=
              Wind_songweiwei_P.u1_Value_p[5 * i_0 + i] *
              Wind_songweiwei_B.Linv[5 * i_1 + i_0];
          }
        }
      }

      /* End of Product: '<S341>/Product1' */

      /* Update for UnitDelay: '<S343>/Lmd_sat' */
      Wind_songweiwei_DW.Lmd_sat_DSTATE = rtb_Lmsatd;
    }

    /* End of Constant: '<S331>/Constant1' */
    /* End of Outputs for SubSystem: '<S331>/Saturation' */

    /* Switch: '<S331>/Switch' incorporates:
     *  Constant: '<S331>/Constant3'
     *  Constant: '<S331>/Constant4'
     *  Product: '<S331>/Product3'
     */
    rtb_RelationalOperator2_idx_0 = (Wind_songweiwei_P.Constant3_Value_er >=
      Wind_songweiwei_P.Switch_Threshold_j);
    for (i = 0; i < 5; i++) {
      for (i_1 = 0; i_1 < 5; i_1++) {
        if (rtb_RelationalOperator2_idx_0) {
          rtb_Sum2_d_0[i_1 + 5 * i] = Wind_songweiwei_B.Linv[5 * i + i_1];
        } else {
          rtb_Sum2_d_0[i_1 + 5 * i] = Wind_songweiwei_P.Constant4_Value_k[5 * i
            + i_1];
        }
      }
    }

    for (i = 0; i < 5; i++) {
      /* Product: '<S331>/Product3' incorporates:
       *  Gain: '<S331>/change Iq Id  current signs'
       *  UnitDelay: '<S342>/fluxes'
       */
      tmp_0[i] = 0.0;
      for (i_1 = 0; i_1 < 5; i_1++) {
        tmp_0[i] += rtb_Sum2_d_0[5 * i_1 + i] *
          Wind_songweiwei_DW.fluxes_DSTATE[i_1];
      }

      /* Gain: '<S331>/change Iq Id  current signs' */
      rtb_changeIqIdcurrentsigns[i] =
        Wind_songweiwei_P.changeIqIdcurrentsigns_Gain[i] * tmp_0[i];
    }

    /* End of Switch: '<S331>/Switch' */

    /* Fcn: '<S333>/Fcn' */
    rtb_pu_a = rtb_changeIqIdcurrentsigns[0] * rtb_ElementaryMath1 +
      rtb_changeIqIdcurrentsigns[1] * rtb_ElementaryMath;

    /* Fcn: '<S333>/Fcn1' */
    rtb_units = ((-rtb_changeIqIdcurrentsigns[0] - 1.7320508075688772 *
                  rtb_changeIqIdcurrentsigns[1]) * rtb_ElementaryMath1 +
                 (1.7320508075688772 * rtb_changeIqIdcurrentsigns[0] -
                  rtb_changeIqIdcurrentsigns[1]) * rtb_ElementaryMath) * 0.5;

    /* Sum: '<S333>/Sum' */
    rtb_Sum1_jb = (0.0 - rtb_pu_a) - rtb_units;

    /* Gain: '<S333>/ib' */
    Wind_songweiwei_B.ib[0] = Wind_songweiwei_P.ib_Gain * rtb_pu_a;
    Wind_songweiwei_B.ib[1] = Wind_songweiwei_P.ib_Gain * rtb_units;
    Wind_songweiwei_B.ib[2] = Wind_songweiwei_P.ib_Gain * rtb_Sum1_jb;

    /* Outputs for Enabled SubSystem: '<S10>/Signal generator' incorporates:
     *  EnablePort: '<S12>/Enable'
     */
    /* Constant: '<S10>/valp_nom3' */
    if (Wind_songweiwei_P.u0kV_VariationEntity - 1.0 > 0.0) {
      if (!Wind_songweiwei_DW.Signalgenerator_MODE) {
        Wind_songweiwei_DW.Signalgenerator_MODE = true;
      }

      /* DigitalClock: '<S15>/t' */
      rtb_Switch4 = (((Wind_songweiwei_M->Timing.clockTick1+
                       Wind_songweiwei_M->Timing.clockTickH1* 4294967296.0)) *
                     2.0E-6);

      /* Lookup: '<S15>/Look-Up Table' */
      Wind_songweiwei_B.LookUpTable = rt_Lookup
        (Wind_songweiwei_P.LookUpTable_XData, 7, rtb_Switch4,
         Wind_songweiwei_P.LookUpTable_YData);

      /* DiscreteIntegrator: '<S12>/Discrete-Time Integrator' */
      Wind_songweiwei_B.DiscreteTimeIntegrator =
        Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_l;

      /* Step: '<S16>/Step1' */
      rtb_Kv_ow_idx_2 = Wind_songweiwei_M->Timing.t[1];
      if (rtb_Kv_ow_idx_2 < Wind_songweiwei_P.VariationSubSystem_Toff_Variati) {
        rtb_Switch4 = Wind_songweiwei_P.Step1_Y0_f;
      } else {
        rtb_Switch4 = Wind_songweiwei_P.Step1_YFinal_m;
      }

      /* End of Step: '<S16>/Step1' */

      /* Step: '<S16>/Step' */
      rtb_Kv_ow_idx_2 = Wind_songweiwei_M->Timing.t[1];
      if (rtb_Kv_ow_idx_2 < Wind_songweiwei_P.VariationSubSystem_Ton_Variatio) {
        rtb_Add4 = Wind_songweiwei_P.Step_Y0_i;
      } else {
        rtb_Add4 = Wind_songweiwei_P.Step_YFinal_b;
      }

      /* End of Step: '<S16>/Step' */

      /* Switch: '<S16>/Switch2' incorporates:
       *  UnitDelay: '<S16>/Unit Delay'
       */
      if (rtb_Switch4 >= Wind_songweiwei_P.Switch2_Threshold) {
        /* MultiPortSwitch: '<S16>/Multiport Switch1' incorporates:
         *  Constant: '<S12>/valp_nom5'
         *  Constant: '<S16>/Constant5'
         *  Constant: '<S16>/valp_nom6'
         *  Constant: '<S16>/valp_nom8'
         *  Constant: '<S16>/valp_nom9'
         *  DiscreteIntegrator: '<S16>/Discrete-Time Integrator1'
         *  Gain: '<S16>/Gain1'
         *  Product: '<S16>/Product'
         *  Product: '<S16>/Product1'
         *  Product: '<S16>/Product2'
         *  Trigonometry: '<S16>/Trigonometric Function1'
         */
        switch ((int32_T)Wind_songweiwei_P.u0kV_VariationType) {
         case 1:
          rtb_Lmsatd = Wind_songweiwei_P.VariationSubSystem_VariationSte *
            rtb_Add4;
          break;

         case 2:
          rtb_Lmsatd = Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTAT_a;
          break;

         case 3:
          rtb_Lmsatd = sin(Wind_songweiwei_P.Gain1_Gain *
                           Wind_songweiwei_P.VariationSubSystem_VariationFre *
                           Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTAT_a) *
            Wind_songweiwei_P.VariationSubSystem_VariationMag;
          break;

         default:
          rtb_Lmsatd = Wind_songweiwei_P.Constant5_Value_h;
          break;
        }

        /* End of MultiPortSwitch: '<S16>/Multiport Switch1' */
      } else {
        rtb_Lmsatd = Wind_songweiwei_DW.UnitDelay_DSTATE_gd;
      }

      /* End of Switch: '<S16>/Switch2' */

      /* Switch: '<S16>/Switch3' incorporates:
       *  Constant: '<S12>/valp_nom5'
       *  Constant: '<S16>/Constant1'
       *  Constant: '<S16>/Constant3'
       *  DataTypeConversion: '<S16>/Data Type  Conversion2'
       *  Logic: '<S16>/Logical Operator'
       *  Logic: '<S16>/Logical Operator1'
       *  RelationalOperator: '<S16>/Relational Operator1'
       */
      if ((!(rtb_Switch4 != 0.0)) && (Wind_songweiwei_P.u0kV_VariationType ==
           Wind_songweiwei_P.Constant3_Value_f)) {
        rtb_Ifdsat = Wind_songweiwei_P.Constant1_Value_h;
      } else {
        rtb_Ifdsat = rtb_Lmsatd;
      }

      /* End of Switch: '<S16>/Switch3' */

      /* Switch: '<S12>/Switch4' incorporates:
       *  Constant: '<S12>/Constant3'
       *  Constant: '<S12>/Constant5'
       *  Constant: '<S12>/valp_nom3'
       *  RelationalOperator: '<S12>/Relational Operator2'
       */
      if (Wind_songweiwei_P.u0kV_VariationEntity ==
          Wind_songweiwei_P.Constant3_Value) {
        rtb_Switch4 = rtb_Ifdsat;
      } else {
        rtb_Switch4 = Wind_songweiwei_P.Constant5_Value;
      }

      /* End of Switch: '<S12>/Switch4' */

      /* RelationalOperator: '<S12>/Relational Operator' incorporates:
       *  Constant: '<S12>/Constant'
       *  Constant: '<S12>/valp_nom3'
       */
      rtb_RelationalOperator2_idx_0 = (Wind_songweiwei_P.u0kV_VariationEntity ==
        Wind_songweiwei_P.Constant_Value_f);

      /* Logic: '<S12>/Logical Operator1' incorporates:
       *  Constant: '<S12>/Constant6'
       *  Constant: '<S12>/valp_nom5'
       *  RelationalOperator: '<S12>/Relational Operator3'
       */
      Wind_songweiwei_B.LogicalOperator1 =
        ((Wind_songweiwei_P.u0kV_VariationType ==
          Wind_songweiwei_P.Constant6_Value) && rtb_RelationalOperator2_idx_0);

      /* Switch: '<S12>/Switch2' incorporates:
       *  Constant: '<S12>/Constant1'
       */
      if (rtb_RelationalOperator2_idx_0) {
        Wind_songweiwei_B.Switch2 = rtb_Ifdsat;
      } else {
        Wind_songweiwei_B.Switch2 = Wind_songweiwei_P.Constant1_Value;
      }

      /* End of Switch: '<S12>/Switch2' */

      /* Switch: '<S12>/Switch3' incorporates:
       *  Constant: '<S12>/Constant2'
       *  Constant: '<S12>/Constant4'
       *  Constant: '<S12>/valp_nom3'
       *  Gain: '<S12>/Gain4'
       *  RelationalOperator: '<S12>/Relational Operator1'
       */
      if (Wind_songweiwei_P.u0kV_VariationEntity ==
          Wind_songweiwei_P.Constant2_Value_m) {
        Wind_songweiwei_B.Switch3 = Wind_songweiwei_P.Gain4_Gain * rtb_Ifdsat;
      } else {
        Wind_songweiwei_B.Switch3 = Wind_songweiwei_P.Constant4_Value;
      }

      /* End of Switch: '<S12>/Switch3' */

      /* Switch: '<S16>/Switch' */
      if (rtb_Add4 >= Wind_songweiwei_P.Switch_Threshold) {
        /* Switch: '<S16>/Switch1' incorporates:
         *  Constant: '<S12>/valp_nom5'
         *  Constant: '<S16>/Constant'
         *  Constant: '<S16>/Constant2'
         *  Constant: '<S16>/valp_nom7'
         *  RelationalOperator: '<S16>/Relational Operator'
         */
        if (Wind_songweiwei_P.u0kV_VariationType ==
            Wind_songweiwei_P.Constant_Value) {
          rtb_Add4 = Wind_songweiwei_P.VariationSubSystem_VariationRat;
        } else {
          rtb_Add4 = Wind_songweiwei_P.Constant2_Value;
        }
      } else {
        /* Switch: '<S16>/Switch1' incorporates:
         *  Constant: '<S16>/Constant4'
         */
        rtb_Add4 = Wind_songweiwei_P.Constant4_Value_d;
      }

      /* End of Switch: '<S16>/Switch' */

      /* Update for DiscreteIntegrator: '<S12>/Discrete-Time Integrator' incorporates:
       *  Gain: '<S12>/Gain2'
       */
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_l +=
        Wind_songweiwei_P.Gain2_Gain * rtb_Switch4 *
        Wind_songweiwei_P.DiscreteTimeIntegrator_gainval;

      /* Update for DiscreteIntegrator: '<S16>/Discrete-Time Integrator1' */
      Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTAT_a +=
        Wind_songweiwei_P.DiscreteTimeIntegrator1_gainval * rtb_Add4;

      /* Update for UnitDelay: '<S16>/Unit Delay' */
      Wind_songweiwei_DW.UnitDelay_DSTATE_gd = rtb_Lmsatd;
    } else {
      if (Wind_songweiwei_DW.Signalgenerator_MODE) {
        /* Disable for Outport: '<S12>/timer' */
        Wind_songweiwei_B.LookUpTable = Wind_songweiwei_P.timer_Y0;

        /* Disable for Outport: '<S12>/selector' */
        Wind_songweiwei_B.LogicalOperator1 = Wind_songweiwei_P.selector_Y0;

        /* Disable for Outport: '<S12>/magnitude' */
        Wind_songweiwei_B.Switch2 = Wind_songweiwei_P.magnitude_Y0;

        /* Disable for Outport: '<S12>/frequency' */
        Wind_songweiwei_B.DiscreteTimeIntegrator =
          Wind_songweiwei_P.frequency_Y0;

        /* Disable for Outport: '<S12>/phase' */
        Wind_songweiwei_B.Switch3 = Wind_songweiwei_P.phase_Y0;
        Wind_songweiwei_DW.Signalgenerator_MODE = false;
      }
    }

    /* End of Constant: '<S10>/valp_nom3' */
    /* End of Outputs for SubSystem: '<S10>/Signal generator' */

    /* Sum: '<S10>/Sum3' incorporates:
     *  Constant: '<S10>/valp_nom2'
     */
    rtb_Sum1_jb = Wind_songweiwei_B.Switch2 +
      Wind_songweiwei_P.valp_nom2_Value_o;

    /* Switch: '<S10>/Switch1' */
    if (Wind_songweiwei_B.LogicalOperator1) {
      rtb_units = Wind_songweiwei_B.LookUpTable;
    } else {
      rtb_units = rtb_Sum1_jb;
    }

    /* End of Switch: '<S10>/Switch1' */

    /* Switch: '<S10>/Switch5' incorporates:
     *  Constant: '<S10>/SinglePhase'
     */
    if (Wind_songweiwei_P.u0kV_VariationPhaseA >=
        Wind_songweiwei_P.Switch5_Threshold) {
      rtb_Kv_o_idx_1 = rtb_Sum1_jb;
      rtb_Kv_o_idx_2 = rtb_Sum1_jb;
    } else {
      rtb_Kv_o_idx_1 = rtb_units;
      rtb_Kv_o_idx_2 = rtb_units;
    }

    /* DigitalClock: '<S10>/t' */
    rtb_Sum1_jb = (((Wind_songweiwei_M->Timing.clockTick1+
                     Wind_songweiwei_M->Timing.clockTickH1* 4294967296.0)) *
                   2.0E-6);

    /* Sum: '<S10>/Sum' incorporates:
     *  Constant: '<S10>/valp_nom1'
     *  Gain: '<S10>/Gain'
     *  Product: '<S10>/Product'
     */
    rtb_Switch1_p = Wind_songweiwei_P.Gain_Gain *
      Wind_songweiwei_P.valp_nom1_Value * rtb_Sum1_jb +
      Wind_songweiwei_B.DiscreteTimeIntegrator;

    /* Gain: '<S10>/Gain3' incorporates:
     *  Constant: '<S10>/valp_nom'
     */
    rtb_phimq = Wind_songweiwei_P.Gain3_Gain_n *
      Wind_songweiwei_P.valp_nom_Value;

    /* Outputs for Enabled SubSystem: '<S10>/Harmonic Generator' incorporates:
     *  EnablePort: '<S11>/Enable'
     */
    /* Constant: '<S10>/valp_nom7' */
    if (Wind_songweiwei_P.u0kV_HarmonicGeneration > 0.0) {
      if (!Wind_songweiwei_DW.HarmonicGenerator_MODE) {
        Wind_songweiwei_DW.HarmonicGenerator_MODE = true;
      }

      /* Gain: '<S13>/Gain1' */
      rtb_Add4 = Wind_songweiwei_P.HarmonicAgeneration_n_Harmo * rtb_Switch1_p;

      /* Gain: '<S13>/Gain3' incorporates:
       *  Constant: '<S13>/Phase_Harmo'
       */
      rtb_Lmsatd = Wind_songweiwei_P.Gain3_Gain *
        Wind_songweiwei_P.HarmonicAgeneration_Phase_Harmo;

      /* MultiPortSwitch: '<S13>/Multiport Switch' incorporates:
       *  Constant: '<S13>/Negative-sequence'
       *  Constant: '<S13>/Phase_Harmo2'
       *  Constant: '<S13>/Positive-sequence'
       *  Constant: '<S13>/Zero-sequence'
       *  Constant: '<S13>/valp_nom2'
       *  Sum: '<S13>/Sum1'
       */
      switch ((int32_T)(Wind_songweiwei_P.HarmonicAgeneration_Seq_Harmo +
                        Wind_songweiwei_P.valp_nom2_Value)) {
       case 1:
        rtb_Ll_d_idx_0 = Wind_songweiwei_P.Zerosequence_Value[0];
        rtb_Ll_d_idx_1 = Wind_songweiwei_P.Zerosequence_Value[1];
        rtb_Ll_d_idx_2 = Wind_songweiwei_P.Zerosequence_Value[2];
        break;

       case 2:
        rtb_Ll_d_idx_0 = Wind_songweiwei_P.Positivesequence_Value[0];
        rtb_Ll_d_idx_1 = Wind_songweiwei_P.Positivesequence_Value[1];
        rtb_Ll_d_idx_2 = Wind_songweiwei_P.Positivesequence_Value[2];
        break;

       default:
        rtb_Ll_d_idx_0 = Wind_songweiwei_P.Negativesequence_Value[0];
        rtb_Ll_d_idx_1 = Wind_songweiwei_P.Negativesequence_Value[1];
        rtb_Ll_d_idx_2 = Wind_songweiwei_P.Negativesequence_Value[2];
        break;
      }

      /* End of MultiPortSwitch: '<S13>/Multiport Switch' */

      /* Step: '<S11>/Step' */
      rtb_Kv_ow_idx_2 = Wind_songweiwei_M->Timing.t[1];
      if (rtb_Kv_ow_idx_2 < Wind_songweiwei_P.Step_Time) {
        rtb_Kv_ow_idx_2 = Wind_songweiwei_P.Step_Y0;
      } else {
        rtb_Kv_ow_idx_2 = Wind_songweiwei_P.Step_YFinal;
      }

      /* End of Step: '<S11>/Step' */

      /* Step: '<S11>/Step1' */
      if (Wind_songweiwei_M->Timing.t[1] < Wind_songweiwei_P.Step1_Time) {
        tmp = Wind_songweiwei_P.Step1_Y0;
      } else {
        tmp = Wind_songweiwei_P.Step1_YFinal;
      }

      /* End of Step: '<S11>/Step1' */

      /* Sum: '<S11>/Sum4' */
      rtb_phimd = rtb_Kv_ow_idx_2 + tmp;

      /* Product: '<S13>/Product1' incorporates:
       *  Constant: '<S13>/Phase_Harmo1'
       *  Sum: '<S13>/Sum'
       *  Trigonometry: '<S13>/Trigonometric Function1'
       */
      Wind_songweiwei_B.Product1[0] = sin((rtb_Add4 + rtb_Lmsatd) +
        rtb_Ll_d_idx_0) * (rtb_phimd *
                           Wind_songweiwei_P.HarmonicAgeneration_Mag_Harmo);
      Wind_songweiwei_B.Product1[1] = sin((rtb_Add4 + rtb_Lmsatd) +
        rtb_Ll_d_idx_1) * (rtb_phimd *
                           Wind_songweiwei_P.HarmonicAgeneration_Mag_Harmo);
      Wind_songweiwei_B.Product1[2] = sin((rtb_Add4 + rtb_Lmsatd) +
        rtb_Ll_d_idx_2) * (rtb_phimd *
                           Wind_songweiwei_P.HarmonicAgeneration_Mag_Harmo);

      /* Gain: '<S14>/Gain1' */
      rtb_Ifdsat = Wind_songweiwei_P.HarmonicBgeneration_n_Harmo * rtb_Switch1_p;

      /* Gain: '<S14>/Gain3' incorporates:
       *  Constant: '<S14>/Phase_Harmo'
       */
      rtb_Lmsatd = Wind_songweiwei_P.Gain3_Gain_p *
        Wind_songweiwei_P.HarmonicBgeneration_Phase_Harmo;

      /* MultiPortSwitch: '<S14>/Multiport Switch' incorporates:
       *  Constant: '<S14>/Negative-sequence'
       *  Constant: '<S14>/Phase_Harmo2'
       *  Constant: '<S14>/Positive-sequence'
       *  Constant: '<S14>/Zero-sequence'
       *  Constant: '<S14>/valp_nom2'
       *  Sum: '<S14>/Sum1'
       */
      switch ((int32_T)(Wind_songweiwei_P.HarmonicBgeneration_Seq_Harmo +
                        Wind_songweiwei_P.valp_nom2_Value_a)) {
       case 1:
        rtb_Ll_d_idx_0 = Wind_songweiwei_P.Zerosequence_Value_b[0];
        rtb_Ll_d_idx_1 = Wind_songweiwei_P.Zerosequence_Value_b[1];
        rtb_Ll_d_idx_2 = Wind_songweiwei_P.Zerosequence_Value_b[2];
        break;

       case 2:
        rtb_Ll_d_idx_0 = Wind_songweiwei_P.Positivesequence_Value_k[0];
        rtb_Ll_d_idx_1 = Wind_songweiwei_P.Positivesequence_Value_k[1];
        rtb_Ll_d_idx_2 = Wind_songweiwei_P.Positivesequence_Value_k[2];
        break;

       default:
        rtb_Ll_d_idx_0 = Wind_songweiwei_P.Negativesequence_Value_m[0];
        rtb_Ll_d_idx_1 = Wind_songweiwei_P.Negativesequence_Value_m[1];
        rtb_Ll_d_idx_2 = Wind_songweiwei_P.Negativesequence_Value_m[2];
        break;
      }

      /* End of MultiPortSwitch: '<S14>/Multiport Switch' */

      /* Product: '<S14>/Product1' incorporates:
       *  Constant: '<S14>/Phase_Harmo1'
       *  Sum: '<S14>/Sum'
       *  Trigonometry: '<S14>/Trigonometric Function1'
       */
      Wind_songweiwei_B.Product1_k[0] = sin((rtb_Ifdsat + rtb_Lmsatd) +
        rtb_Ll_d_idx_0) * (rtb_phimd *
                           Wind_songweiwei_P.HarmonicBgeneration_Mag_Harmo);
      Wind_songweiwei_B.Product1_k[1] = sin((rtb_Ifdsat + rtb_Lmsatd) +
        rtb_Ll_d_idx_1) * (rtb_phimd *
                           Wind_songweiwei_P.HarmonicBgeneration_Mag_Harmo);
      Wind_songweiwei_B.Product1_k[2] = sin((rtb_Ifdsat + rtb_Lmsatd) +
        rtb_Ll_d_idx_2) * (rtb_phimd *
                           Wind_songweiwei_P.HarmonicBgeneration_Mag_Harmo);
    } else {
      if (Wind_songweiwei_DW.HarmonicGenerator_MODE) {
        /* Disable for Outport: '<S11>/Out1' */
        Wind_songweiwei_B.Product1[0] = Wind_songweiwei_P.Out1_Y0;
        Wind_songweiwei_B.Product1[1] = Wind_songweiwei_P.Out1_Y0;
        Wind_songweiwei_B.Product1[2] = Wind_songweiwei_P.Out1_Y0;

        /* Disable for Outport: '<S11>/Out2' */
        Wind_songweiwei_B.Product1_k[0] = Wind_songweiwei_P.Out2_Y0;
        Wind_songweiwei_B.Product1_k[1] = Wind_songweiwei_P.Out2_Y0;
        Wind_songweiwei_B.Product1_k[2] = Wind_songweiwei_P.Out2_Y0;
        Wind_songweiwei_DW.HarmonicGenerator_MODE = false;
      }
    }

    /* End of Constant: '<S10>/valp_nom7' */
    /* End of Outputs for SubSystem: '<S10>/Harmonic Generator' */

    /* Sum: '<S10>/Sum5' incorporates:
     *  Constant: '<S10>/P1'
     *  Product: '<S10>/Product1'
     *  Sum: '<S10>/Sum1'
     *  Sum: '<S10>/Sum2'
     *  Switch: '<S10>/Switch5'
     *  Trigonometry: '<S10>/Trigonometric Function'
     */
    Wind_songweiwei_B.Sum5[0] = (sin(((rtb_phimq + Wind_songweiwei_P.P1_Value[0])
      + Wind_songweiwei_B.Switch3) + rtb_Switch1_p) * rtb_units +
      Wind_songweiwei_B.Product1[0]) + Wind_songweiwei_B.Product1_k[0];
    Wind_songweiwei_B.Sum5[1] = (sin(((rtb_phimq + Wind_songweiwei_P.P1_Value[1])
      + Wind_songweiwei_B.Switch3) + rtb_Switch1_p) * rtb_Kv_o_idx_1 +
      Wind_songweiwei_B.Product1[1]) + Wind_songweiwei_B.Product1_k[1];
    Wind_songweiwei_B.Sum5[2] = (sin(((rtb_phimq + Wind_songweiwei_P.P1_Value[2])
      + Wind_songweiwei_B.Switch3) + rtb_Switch1_p) * rtb_Kv_o_idx_2 +
      Wind_songweiwei_B.Product1[2]) + Wind_songweiwei_B.Product1_k[2];

    /* S-Function block: <S357>/State-Space */
    {
      real_T accum;

      /* Circuit has switches */
      int_T *switch_status = (int_T*)
        Wind_songweiwei_DW.StateSpace_PWORK.SWITCH_STATUS;
      int_T *switch_status_init = (int_T*)
        Wind_songweiwei_DW.StateSpace_PWORK.SWITCH_STATUS_INIT;
      int_T *SwitchChange = (int_T*) Wind_songweiwei_DW.StateSpace_PWORK.SW_CHG;
      int_T *Chopper = (int_T*) Wind_songweiwei_DW.StateSpace_PWORK.CHOPPER;
      int_T *gState = (int_T*) Wind_songweiwei_DW.StateSpace_PWORK.G_STATE;
      real_T *yswitch = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.Y_SWITCH;
      int_T *switchTypes = (int_T*)
        Wind_songweiwei_DW.StateSpace_PWORK.SWITCH_TYPES;
      int_T *idxOutSw = (int_T*) Wind_songweiwei_DW.StateSpace_PWORK.IDX_OUT_SW;
      real_T *DxCol = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.DX_COL;
      real_T *tmp2 = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.TMP2;
      real_T *BDcol = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.BD_COL;
      real_T *tmp1 = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.TMP1;
      int_T newState;
      int_T swChanged = 0;
      int loopsToDo = 20;
      real_T temp;

      /* keep an initial copy of switch_status*/
      memcpy(switch_status_init, switch_status, 19 * sizeof(int_T));
      do {
        if (loopsToDo == 1) {          /* Need to reset some variables: */
          swChanged = 0;

          /* return to the original switch status*/
          {
            int_T i1;
            for (i1=0; i1 < 19; i1++) {
              swChanged = ((SwitchChange[i1] = switch_status_init[i1] -
                            switch_status[i1]) != 0) ? 1 : swChanged;
              switch_status[i1] = switch_status_init[i1];
            }
          }
        } else {
          /*
           * Compute outputs:
           * ---------------
           */

          /*
           * Chopper parameter will force zero current (y[i])
           * for an open switch.
           */
          real_T *Cs = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.CS;
          real_T *Ds = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.DS;

          {
            int_T i1;
            real_T *y0 = &Wind_songweiwei_B.StateSpace_o1[0];
            for (i1=0; i1 < 50; i1++) {
              accum = 0.0;

              {
                int_T i2;
                real_T *xd = &Wind_songweiwei_DW.StateSpace_DSTATE[0];
                for (i2=0; i2 < 45; i2++) {
                  accum += *(Cs++) * xd[i2];
                }
              }

              {
                int_T i2;
                const real_T *u0 = &Wind_songweiwei_P.SwitchCurrents_Value[0];
                for (i2=0; i2 < 13; i2++) {
                  accum += *(Ds++) * u0[i2];
                }

                u0 = &Wind_songweiwei_B.Product[0];
                for (i2=0; i2 < 6; i2++) {
                  accum += *(Ds++) * u0[i2];
                }

                accum += *(Ds++) * Wind_songweiwei_B.ib[0];
                accum += *(Ds++) * Wind_songweiwei_B.ib[1];
                accum += *(Ds++) * Wind_songweiwei_B.Sum5[0];
                accum += *(Ds++) * Wind_songweiwei_B.Sum5[1];
                accum += *(Ds++) * Wind_songweiwei_B.Sum5[2];
              }

              y0[i1] = accum * Chopper[i1];
            }
          }

          swChanged = 0;

          {
            int_T i1;
            real_T *y0 = &Wind_songweiwei_B.StateSpace_o1[0];
            for (i1=0; i1 < 19; i1++) {
              switch (switchTypes[i1]) {
               case 1 :                /* Ideal switch */
                newState = gState[i1] > 0 ? 1 : 0;
                break;

               case 3 :                /* Diodes */
                newState = y0[i1] > 0.0 ? 1 : ((y0[i1] < 0.0) ? 0 :
                  switch_status[i1]);
                break;

               case 7 :                /* MOSFETs or IGBT/Diode pairs */
                newState = ((y0[i1] > 0.0) && (gState[i1] > 0)) || (y0[i1] < 0.0)
                  ? 1 : (((y0[i1] > 0.0) && gState[i1] == 0) ? 0 :
                         switch_status[i1]);
                break;
              }

              swChanged = ((SwitchChange[i1] = newState - switch_status[i1]) !=
                           0) ? 1 : swChanged;
              switch_status[i1] = newState;/* Keep new state */
            }
          }
        }

        /*
         * Compute new As, Bs, Cs and Ds matrixes:
         * --------------------------------------
         */
        if (swChanged) {
          real_T *As = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.AS;
          real_T *Cs = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.CS;
          real_T *Bs = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.BS;
          real_T *Ds = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.DS;
          real_T a1;

          {
            int_T i1;
            for (i1=0; i1 < 19; i1++) {
              if (SwitchChange[i1] != 0) {
                if (idxOutSw[i1] > -1 ) {/* A positive index points to a switch actual measure output */
                  Chopper[idxOutSw[i1]] = switch_status[i1];
                }

                a1 = yswitch[i1]*SwitchChange[i1];
                temp = 1/(1-Ds[i1*25]*a1);

                {
                  int_T i2;
                  for (i2=0; i2 < 50; i2++) {
                    DxCol[i2]= Ds[i2 * 24 + i1]*temp*a1;
                  }
                }

                DxCol[i1] = temp;

                {
                  int_T i2;
                  for (i2=0; i2 < 45; i2++) {
                    BDcol[i2]= Bs[i2 * 24 + i1]*a1;
                  }
                }

                /* Copy row nSw of Cs into tmp1 and zero it out in Cs */
                memcpy(tmp1, &Cs[i1 * 45], 45 * sizeof(real_T));
                memset(&Cs[i1 * 45], '\0', 45 * sizeof(real_T));

                /* Copy row nSw of Ds into tmp2 and zero it out in Ds */
                memcpy(tmp2, &Ds[i1 * 24], 24 * sizeof(real_T));
                memset(&Ds[i1 * 24], '\0', 24 * sizeof(real_T));

                /* Cs = Cs + DxCol * tmp1, Ds = Ds + DxCol * tmp2 *******************/
                {
                  int_T i2;
                  for (i2=0; i2 < 50; i2++) {
                    a1 = DxCol[i2];

                    {
                      int_T i3;
                      for (i3=0; i3 < 45; i3++) {
                        Cs[i2 * 45 + i3] += a1 * tmp1[i3];
                      }
                    }

                    {
                      int_T i3;
                      for (i3=0; i3 < 24; i3++) {
                        Ds[i2 * 24 + i3] += a1 * tmp2[i3];
                      }
                    }
                  }
                }

                /* As = As + BdCol*Cs(nSw,:), Bs = Bs + BdCol*Ds(nSw,:) *************/
                {
                  int_T i2;
                  for (i2=0; i2 < 45; i2++) {
                    a1 = BDcol[i2];

                    {
                      int_T i3;
                      for (i3=0; i3 < 45; i3++) {
                        As[i2 * 45 + i3] += a1 * Cs[i1 * 45 + i3];
                      }
                    }

                    {
                      int_T i3;
                      for (i3=0; i3 < 24; i3++) {
                        Bs[i2 * 24 + i3] += a1 * Ds[i1 * 24 + i3];
                      }
                    }
                  }
                }
              }
            }
          }
        }                              /* if (swChanged) */
      } while (swChanged > 0 && --loopsToDo > 0);

      if (loopsToDo == 0) {
        real_T *Cs = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.CS;
        real_T *Ds = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.DS;

        {
          int_T i1;
          real_T *y0 = &Wind_songweiwei_B.StateSpace_o1[0];
          for (i1=0; i1 < 50; i1++) {
            accum = 0.0;

            {
              int_T i2;
              real_T *xd = &Wind_songweiwei_DW.StateSpace_DSTATE[0];
              for (i2=0; i2 < 45; i2++) {
                accum += *(Cs++) * xd[i2];
              }
            }

            {
              int_T i2;
              const real_T *u0 = &Wind_songweiwei_P.SwitchCurrents_Value[0];
              for (i2=0; i2 < 13; i2++) {
                accum += *(Ds++) * u0[i2];
              }

              u0 = &Wind_songweiwei_B.Product[0];
              for (i2=0; i2 < 6; i2++) {
                accum += *(Ds++) * u0[i2];
              }

              accum += *(Ds++) * Wind_songweiwei_B.ib[0];
              accum += *(Ds++) * Wind_songweiwei_B.ib[1];
              accum += *(Ds++) * Wind_songweiwei_B.Sum5[0];
              accum += *(Ds++) * Wind_songweiwei_B.Sum5[1];
              accum += *(Ds++) * Wind_songweiwei_B.Sum5[2];
            }

            y0[i1] = accum * Chopper[i1];
          }
        }
      }

      /* Output new switches states */
      {
        int_T i1;
        real_T *y1 = &Wind_songweiwei_B.StateSpace_o2[0];
        for (i1=0; i1 < 19; i1++) {
          y1[i1] = (real_T)switch_status[i1];
        }
      }
    }

    /* Gain: '<S49>/do not delete this gain' */
    rtb_Sum1_jb = Wind_songweiwei_P.donotdeletethisgain_Gain *
      Wind_songweiwei_B.StateSpace_o1[21];

    /* Outport: '<Root>/Out_Vabc_B575' incorporates:
     *  Gain: '<S50>/do not delete this gain'
     *  Gain: '<S51>/do not delete this gain'
     *  Gain: '<S6>/Kv1'
     */
    Wind_songweiwei_Y.Out_Vabc_B575[0] = Wind_songweiwei_P.Kv1_Gain *
      rtb_Sum1_jb;
    Wind_songweiwei_Y.Out_Vabc_B575[1] =
      Wind_songweiwei_P.donotdeletethisgain_Gain_k *
      Wind_songweiwei_B.StateSpace_o1[22] * Wind_songweiwei_P.Kv1_Gain;
    Wind_songweiwei_Y.Out_Vabc_B575[2] =
      Wind_songweiwei_P.donotdeletethisgain_Gain_g *
      Wind_songweiwei_B.StateSpace_o1[23] * Wind_songweiwei_P.Kv1_Gain;

    /* Gain: '<S46>/do not delete this gain' */
    rtb_Sum1_jb = Wind_songweiwei_P.donotdeletethisgain_Gain_p *
      Wind_songweiwei_B.StateSpace_o1[36];

    /* Outport: '<Root>/Out_labc_B575' incorporates:
     *  Gain: '<S47>/do not delete this gain'
     *  Gain: '<S48>/do not delete this gain'
     *  Gain: '<S6>/Kv'
     */
    Wind_songweiwei_Y.Out_labc_B575[0] = Wind_songweiwei_P.Kv_Gain * rtb_Sum1_jb;
    Wind_songweiwei_Y.Out_labc_B575[1] =
      Wind_songweiwei_P.donotdeletethisgain_Gain_c *
      Wind_songweiwei_B.StateSpace_o1[37] * Wind_songweiwei_P.Kv_Gain;
    Wind_songweiwei_Y.Out_labc_B575[2] =
      Wind_songweiwei_P.donotdeletethisgain_Gain_e *
      Wind_songweiwei_B.StateSpace_o1[38] * Wind_songweiwei_P.Kv_Gain;

    /* UnitDelay: '<S67>/Unit Delay3' */
    rtb_Sum1_jb = Wind_songweiwei_DW.UnitDelay3_DSTATE;

    /* Outport: '<Root>/Out_P' incorporates:
     *  Gain: '<Root>/MW'
     */
    Wind_songweiwei_Y.Out_P = Wind_songweiwei_P.MW_Gain * rtb_Sum1_jb;

    /* Outport: '<Root>/Out_Q' incorporates:
     *  Gain: '<Root>/MW'
     *  UnitDelay: '<S67>/Unit Delay7'
     */
    Wind_songweiwei_Y.Out_Q = Wind_songweiwei_P.MW_Gain *
      Wind_songweiwei_DW.UnitDelay7_DSTATE;

    /* Gain: '<S76>/Gain2' */
    rtb_Gain2_m_idx_0 = Wind_songweiwei_P.Gain2_Gain_i[0] *
      Wind_songweiwei_B.StateSpace_o1[48];

    /* Outport: '<Root>/Out_Vdc' */
    Wind_songweiwei_Y.Out_Vdc = rtb_Gain2_m_idx_0;

    /* UnitDelay: '<S351>/dw_delay' */
    rtb_dw_delay = Wind_songweiwei_DW.dw_delay_DSTATE;

    /* Sum: '<S351>/Sum1' incorporates:
     *  Gain: '<S351>/F2'
     *  UnitDelay: '<S351>/dw_predict'
     */
    rtb_Sum1_jb = Wind_songweiwei_P.F2_Gain * rtb_dw_delay -
      Wind_songweiwei_DW.dw_predict_DSTATE;

    /* Sum: '<S329>/Sum' incorporates:
     *  Constant: '<S329>/nominal speed'
     */
    rtb_n = Wind_songweiwei_P.nominalspeed_Value + rtb_Sum1_jb;

    /* Gain: '<S329>/units' */
    rtb_units = Wind_songweiwei_P.units_Gain * rtb_n;

    /* Outport: '<Root>/Out_wr' */
    Wind_songweiwei_Y.Out_wr = rtb_units;

    /* Saturate: '<S353>/Saturation' incorporates:
     *  Constant: '<S353>/0 4'
     *  Gain: '<S353>/1//Ron'
     *  Switch: '<S353>/Switch'
     */
    for (i = 0; i < 6; i++) {
      if (Wind_songweiwei_B.StateSpace_o2[i + 13] >=
          Wind_songweiwei_P.Switch_Threshold_l) {
        unnamed_idx_1 = Wind_songweiwei_B.StateSpace_o1[i + 13] *
          Wind_songweiwei_P.Ron_Gain;
      } else {
        unnamed_idx_1 = Wind_songweiwei_P.u_Value;
      }

      if (unnamed_idx_1 > Wind_songweiwei_P.Saturation_UpperSat_c) {
        rtb_Saturation[i] = Wind_songweiwei_P.Saturation_UpperSat_c;
      } else if (unnamed_idx_1 < Wind_songweiwei_P.Saturation_LowerSat_i) {
        rtb_Saturation[i] = Wind_songweiwei_P.Saturation_LowerSat_i;
      } else {
        rtb_Saturation[i] = unnamed_idx_1;
      }
    }

    /* End of Saturate: '<S353>/Saturation' */

    /* Gain: '<S8>/-> pu' incorporates:
     *  Gain: '<S76>/Gain2'
     */
    rtb_pu_a = Wind_songweiwei_P.Gain2_Gain_i[1] *
      Wind_songweiwei_B.StateSpace_o1[49] * Wind_songweiwei_P.pu_Gain;

    /* Gain: '<S64>/Kv' incorporates:
     *  Gain: '<S84>/do not delete this gain'
     *  Gain: '<S85>/do not delete this gain'
     *  Gain: '<S86>/do not delete this gain'
     */
    rtb_Kv_o_idx_0 = Wind_songweiwei_P.donotdeletethisgain_Gain_ga *
      Wind_songweiwei_B.StateSpace_o1[42] * Wind_songweiwei_P.Kv_Gain_n;
    rtb_Kv_o_idx_1 = Wind_songweiwei_P.donotdeletethisgain_Gain_i *
      Wind_songweiwei_B.StateSpace_o1[43] * Wind_songweiwei_P.Kv_Gain_n;
    rtb_Kv_o_idx_2 = Wind_songweiwei_P.donotdeletethisgain_Gain_d *
      Wind_songweiwei_B.StateSpace_o1[44] * Wind_songweiwei_P.Kv_Gain_n;

    /* Gain: '<S64>/Kv1' incorporates:
     *  Gain: '<S87>/do not delete this gain'
     *  Gain: '<S88>/do not delete this gain'
     *  Gain: '<S89>/do not delete this gain'
     */
    rtb_Kv1_idx_0 = Wind_songweiwei_P.donotdeletethisgain_Gain_m *
      Wind_songweiwei_B.StateSpace_o1[27] * Wind_songweiwei_P.Kv1_Gain_c;
    rtb_Kv1_idx_1 = Wind_songweiwei_P.donotdeletethisgain_Gain_eg *
      Wind_songweiwei_B.StateSpace_o1[28] * Wind_songweiwei_P.Kv1_Gain_c;
    rtb_Kv1_idx_2 = Wind_songweiwei_P.donotdeletethisgain_Gain_cj *
      Wind_songweiwei_B.StateSpace_o1[29] * Wind_songweiwei_P.Kv1_Gain_c;

    /* Gain: '<S65>/Kv' incorporates:
     *  Gain: '<S100>/do not delete this gain'
     *  Gain: '<S101>/do not delete this gain'
     *  Gain: '<S99>/do not delete this gain'
     */
    rtb_Kv_idx_0 = Wind_songweiwei_P.donotdeletethisgain_Gain_i1 *
      Wind_songweiwei_B.StateSpace_o1[45] * Wind_songweiwei_P.Kv_Gain_h;
    rtb_Kv_idx_1 = Wind_songweiwei_P.donotdeletethisgain_Gain_mu *
      Wind_songweiwei_B.StateSpace_o1[46] * Wind_songweiwei_P.Kv_Gain_h;
    rtb_Kv_idx_2 = Wind_songweiwei_P.donotdeletethisgain_Gain_d5 *
      Wind_songweiwei_B.StateSpace_o1[47] * Wind_songweiwei_P.Kv_Gain_h;

    /* Gain: '<S65>/Kv1' incorporates:
     *  Gain: '<S102>/do not delete this gain'
     *  Gain: '<S103>/do not delete this gain'
     *  Gain: '<S104>/do not delete this gain'
     */
    rtb_Kv1_j_idx_0 = Wind_songweiwei_P.donotdeletethisgain_Gain_m4 *
      Wind_songweiwei_B.StateSpace_o1[33] * Wind_songweiwei_P.Kv1_Gain_cr;
    rtb_Kv1_j_idx_1 = Wind_songweiwei_P.donotdeletethisgain_Gain_h *
      Wind_songweiwei_B.StateSpace_o1[34] * Wind_songweiwei_P.Kv1_Gain_cr;
    rtb_Kv1_j_idx_2 = Wind_songweiwei_P.donotdeletethisgain_Gain_kc *
      Wind_songweiwei_B.StateSpace_o1[35] * Wind_songweiwei_P.Kv1_Gain_cr;

    /* Sum: '<S221>/C*X(k)+D*u(k)' incorporates:
     *  Gain: '<S221>/D*u(k)'
     *  Gain: '<S224>/C11'
     *  Gain: '<S224>/C12'
     *  Sum: '<S224>/sum2'
     *  UnitDelay: '<S221>/Delay_x1'
     *  UnitDelay: '<S221>/Delay_x2'
     */
    rtb_yk = (Wind_songweiwei_P.C11_Gain * Wind_songweiwei_DW.Delay_x1_DSTATE +
              Wind_songweiwei_P.C12_Gain * Wind_songweiwei_DW.Delay_x2_DSTATE) +
      Wind_songweiwei_P.Duk_Gain * rtb_units;

    /* Switch: '<S120>/Switch' incorporates:
     *  Constant: '<S120>/Flux_ref'
     *  Constant: '<S120>/Flux_ref '
     *  Product: '<S120>/Divide'
     *  Saturate: '<S120>/avoid division by 0'
     */
    if (rtb_yk >= Wind_songweiwei_P.Switch_Threshold_p) {
      /* Saturate: '<S120>/avoid division by 0' */
      if (rtb_yk > Wind_songweiwei_P.avoiddivisionby0_UpperSat) {
        rtb_Kv_ow_idx_2 = Wind_songweiwei_P.avoiddivisionby0_UpperSat;
      } else if (rtb_yk < Wind_songweiwei_P.avoiddivisionby0_LowerSat) {
        rtb_Kv_ow_idx_2 = Wind_songweiwei_P.avoiddivisionby0_LowerSat;
      } else {
        rtb_Kv_ow_idx_2 = rtb_yk;
      }

      rtb_phimd = Wind_songweiwei_P.Flux_ref_Value_c / rtb_Kv_ow_idx_2;
    } else {
      rtb_phimd = Wind_songweiwei_P.Flux_ref_Value;
    }

    /* End of Switch: '<S120>/Switch' */

    /* Sum: '<S197>/C*X(k)+D*u(k)' incorporates:
     *  Gain: '<S197>/D*u(k)'
     *  Gain: '<S200>/C11'
     *  Gain: '<S200>/C12'
     *  Sum: '<S200>/sum2'
     *  UnitDelay: '<S197>/Delay_x1'
     *  UnitDelay: '<S197>/Delay_x2'
     */
    rtb_Ll_d_idx_0 = (Wind_songweiwei_P.C11_Gain_m *
                      Wind_songweiwei_DW.Delay_x1_DSTATE_c[0] +
                      Wind_songweiwei_P.C12_Gain_b *
                      Wind_songweiwei_DW.Delay_x2_DSTATE_e[0]) +
      Wind_songweiwei_P.Duk_Gain_l * rtb_Kv1_idx_0;
    rtb_Ll_d_idx_1 = (Wind_songweiwei_P.C11_Gain_m *
                      Wind_songweiwei_DW.Delay_x1_DSTATE_c[1] +
                      Wind_songweiwei_P.C12_Gain_b *
                      Wind_songweiwei_DW.Delay_x2_DSTATE_e[1]) +
      Wind_songweiwei_P.Duk_Gain_l * rtb_Kv1_idx_1;

    /* Gain: '<S8>/deg->rd' incorporates:
     *  DiscreteIntegrator: '<S329>/theta'
     *  Fcn: '<S329>/Fcn'
     *  Gain: '<S329>/t'
     */
    rtb_degrd = Wind_songweiwei_P.t_Gain * rt_remd_snf
      (Wind_songweiwei_DW.theta_DSTATE, 6.2831853071795862) *
      Wind_songweiwei_P.degrd_Gain;

    /* Math: '<S122>/Math Function' incorporates:
     *  Constant: '<S122>/Constant4'
     *  Gain: '<S122>/# pairs of poles'
     *  Gain: '<S213>/D*u(k)'
     *  Gain: '<S216>/C11'
     *  Gain: '<S216>/C12'
     *  Sum: '<S213>/C*X(k)+D*u(k)'
     *  Sum: '<S216>/sum2'
     *  UnitDelay: '<S213>/Delay_x1'
     *  UnitDelay: '<S213>/Delay_x2'
     */
    rtb_MathFunction = rt_modd_snf(((Wind_songweiwei_P.C11_Gain_k *
      Wind_songweiwei_DW.Delay_x1_DSTATE_a + Wind_songweiwei_P.C12_Gain_h *
      Wind_songweiwei_DW.Delay_x2_DSTATE_p) + Wind_songweiwei_P.Duk_Gain_g *
      rtb_degrd) * Wind_songweiwei_P.pairsofpoles_Gain,
      Wind_songweiwei_P.Constant4_Value_m);

    /* Trigonometry: '<S122>/Trigonometric Function4' */
    rtb_Switch_e3 = sin(rtb_MathFunction);

    /* Trigonometry: '<S122>/Trigonometric Function5' */
    rtb_Switch1_p = cos(rtb_MathFunction);

    /* Sum: '<S201>/C*X(k)+D*u(k)' incorporates:
     *  Gain: '<S201>/D*u(k)'
     *  Gain: '<S204>/C11'
     *  Gain: '<S204>/C12'
     *  Gain: '<S275>/Gain3'
     *  Sum: '<S204>/sum2'
     *  UnitDelay: '<S201>/Delay_x1'
     *  UnitDelay: '<S201>/Delay_x2'
     */
    rtb_Kv_ow_idx_2 = (Wind_songweiwei_P.C11_Gain_o *
                       Wind_songweiwei_DW.Delay_x1_DSTATE_i[0] +
                       Wind_songweiwei_P.C12_Gain_j *
                       Wind_songweiwei_DW.Delay_x2_DSTATE_m[0]) +
      Wind_songweiwei_P.Duk_Gain_m * rtb_Kv_o_idx_0;
    unnamed_idx_1 = (Wind_songweiwei_P.C11_Gain_o *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_i[1] +
                     Wind_songweiwei_P.C12_Gain_j *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_m[1]) +
      Wind_songweiwei_P.Duk_Gain_m * rtb_Kv_o_idx_1;
    rtb_Ifdsat = (Wind_songweiwei_P.C11_Gain_o *
                  Wind_songweiwei_DW.Delay_x1_DSTATE_i[2] +
                  Wind_songweiwei_P.C12_Gain_j *
                  Wind_songweiwei_DW.Delay_x2_DSTATE_m[2]) +
      Wind_songweiwei_P.Duk_Gain_m * rtb_Kv_o_idx_2;

    /* Gain: '<S275>/Gain3' incorporates:
     *  Gain: '<S275>/Gain1'
     */
    for (i = 0; i < 3; i++) {
      tmp_1[i] = Wind_songweiwei_P.Gain3_Gain_m[i + 6] * rtb_Ifdsat +
        (Wind_songweiwei_P.Gain3_Gain_m[i + 3] * unnamed_idx_1 +
         Wind_songweiwei_P.Gain3_Gain_m[i] * rtb_Kv_ow_idx_2);
    }

    /* Gain: '<S275>/Gain1' */
    rtb_Gain1[0] = Wind_songweiwei_P.Gain1_Gain_h * tmp_1[0];
    rtb_Gain1[1] = Wind_songweiwei_P.Gain1_Gain_h * tmp_1[1];
    rtb_Gain1[2] = Wind_songweiwei_P.Gain1_Gain_h * tmp_1[2];

    /* RelationalOperator: '<S276>/Compare' incorporates:
     *  Constant: '<S274>/Constant'
     *  Constant: '<S276>/Constant'
     */
    rtb_Compare = (uint8_T)(Wind_songweiwei_P.AlphaBetaZerotodq0_Alignment_i ==
      Wind_songweiwei_P.CompareToConstant_const_o);

    /* Outputs for Enabled SubSystem: '<S274>/Subsystem1' */
    Wind_songweiwei_Subsystem1(rtb_Compare, &rtb_Gain1[0], rtb_MathFunction,
      &Wind_songweiwei_B.Subsystem1);

    /* End of Outputs for SubSystem: '<S274>/Subsystem1' */

    /* RelationalOperator: '<S277>/Compare' incorporates:
     *  Constant: '<S274>/Constant'
     *  Constant: '<S277>/Constant'
     */
    rtb_Compare_a = (uint8_T)(Wind_songweiwei_P.AlphaBetaZerotodq0_Alignment_i ==
      Wind_songweiwei_P.CompareToConstant1_const_o);

    /* Outputs for Enabled SubSystem: '<S274>/Subsystem - pi//2 delay' */
    Wind_songweiw_Subsystempi2delay(rtb_Compare_a, &rtb_Gain1[0],
      rtb_MathFunction, &Wind_songweiwei_B.Subsystempi2delay);

    /* End of Outputs for SubSystem: '<S274>/Subsystem - pi//2 delay' */

    /* Saturate: '<S145>/avoid division by 0' */
    if (rtb_yk > Wind_songweiwei_P.avoiddivisionby0_UpperSat_d) {
      rtb_phimq = Wind_songweiwei_P.avoiddivisionby0_UpperSat_d;
    } else if (rtb_yk < Wind_songweiwei_P.avoiddivisionby0_LowerSat_n) {
      rtb_phimq = Wind_songweiwei_P.avoiddivisionby0_LowerSat_n;
    } else {
      rtb_phimq = rtb_yk;
    }

    /* End of Saturate: '<S145>/avoid division by 0' */

    /* Switch: '<S274>/Switch' */
    if (rtb_Compare != 0) {
      tmp = Wind_songweiwei_B.Subsystem1.Fcn1;
    } else {
      tmp = Wind_songweiwei_B.Subsystempi2delay.Fcn1;
    }

    if (rtb_Compare != 0) {
      unnamed_idx_1 = Wind_songweiwei_B.Subsystem1.Fcn;
    } else {
      unnamed_idx_1 = Wind_songweiwei_B.Subsystempi2delay.Fcn;
    }

    /* End of Switch: '<S274>/Switch' */

    /* Fcn: '<S145>/Magnitude of flux' incorporates:
     *  Constant: '<S145>/Rs2'
     *  Constant: '<S145>/Rs4'
     *  Fcn: '<S149>/Fcn2'
     *  Fcn: '<S149>/Fcn3'
     *  Gain: '<S145>/Gain'
     *  Product: '<S145>/Divide2'
     *  Product: '<S145>/Divide4'
     *  Product: '<S145>/Divide7'
     *  Product: '<S145>/Divide8'
     *  Sum: '<S145>/Sum2'
     *  Sum: '<S145>/Sum4'
     */
    tmp = rt_powd_snf((((2.0 * rtb_Ll_d_idx_0 + rtb_Ll_d_idx_1) * rtb_Switch1_p
                        + 1.7320508075688772 * rtb_Ll_d_idx_1 * rtb_Switch_e3) /
                       3.0 - tmp * Wind_songweiwei_P.WindTurbineType4_Resistance)
                      / rtb_phimq, 2.0) + rt_powd_snf((((2.0 * rtb_Ll_d_idx_0 +
      rtb_Ll_d_idx_1) * rtb_Switch_e3 + -1.7320508075688772 * rtb_Ll_d_idx_1 *
      rtb_Switch1_p) / 3.0 - unnamed_idx_1 *
      Wind_songweiwei_P.WindTurbineType4_Resistance) /
      (Wind_songweiwei_P.Gain_Gain_b * rtb_phimq), 2.0);
    if (tmp < 0.0) {
      tmp = -sqrt(-tmp);
    } else {
      tmp = sqrt(tmp);
    }

    /* Sum: '<S120>/Sum4' incorporates:
     *  Fcn: '<S145>/Magnitude of flux'
     *  Gain: '<S261>/C'
     *  Gain: '<S261>/D'
     *  Sum: '<S261>/sum1'
     *  UnitDelay: '<S261>/Delay_x'
     */
    rtb_phimd -= Wind_songweiwei_P.D_Gain * tmp + Wind_songweiwei_P.C_Gain *
      Wind_songweiwei_DW.Delay_x_DSTATE;

    /* Gain: '<S126>/Kp5' */
    rtb_Kp5 = Wind_songweiwei_P.DiscretePIController_Ki * rtb_phimd;

    /* Sum: '<S126>/Sum6' incorporates:
     *  DiscreteIntegrator: '<S126>/Discrete-Time Integrator'
     *  Gain: '<S126>/Kp4'
     */
    unnamed_idx_1 = Wind_songweiwei_P.DiscretePIController_Kp * rtb_phimd +
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE;

    /* Saturate: '<S126>/Saturation2' */
    if (unnamed_idx_1 > Wind_songweiwei_P.Saturation2_UpperSat_o) {
      rtb_Saturation2 = Wind_songweiwei_P.Saturation2_UpperSat_o;
    } else if (unnamed_idx_1 < Wind_songweiwei_P.Saturation2_LowerSat_f) {
      rtb_Saturation2 = Wind_songweiwei_P.Saturation2_LowerSat_f;
    } else {
      rtb_Saturation2 = unnamed_idx_1;
    }

    /* End of Saturate: '<S126>/Saturation2' */

    /* Sum: '<S217>/C*X(k)+D*u(k)' incorporates:
     *  Gain: '<S217>/D*u(k)'
     *  Gain: '<S220>/C11'
     *  Gain: '<S220>/C12'
     *  Sum: '<S220>/sum2'
     *  UnitDelay: '<S217>/Delay_x1'
     *  UnitDelay: '<S217>/Delay_x2'
     */
    rtb_Cp = (Wind_songweiwei_P.C11_Gain_e *
              Wind_songweiwei_DW.Delay_x1_DSTATE_ck +
              Wind_songweiwei_P.C12_Gain_a *
              Wind_songweiwei_DW.Delay_x2_DSTATE_ee) +
      Wind_songweiwei_P.Duk_Gain_a * rtb_Gain2_m_idx_0;

    /* Gain: '<S121>/->pu' incorporates:
     *  Constant: '<S121>/Vdc_ref (V)'
     *  Sum: '<S121>/Sum7'
     */
    rtb_phimq = (Wind_songweiwei_P.Vdc_refV_Value - rtb_Cp) *
      Wind_songweiwei_P.pu_Gain_b;

    /* Sum: '<S189>/C*X(k)+D*u(k)' incorporates:
     *  Gain: '<S189>/D*u(k)'
     *  Gain: '<S192>/C11'
     *  Gain: '<S192>/C12'
     *  Sum: '<S192>/sum2'
     *  UnitDelay: '<S189>/Delay_x1'
     *  UnitDelay: '<S189>/Delay_x2'
     */
    rtb_Ll_d_idx_0 = (Wind_songweiwei_P.C11_Gain_a *
                      Wind_songweiwei_DW.Delay_x1_DSTATE_m[0] +
                      Wind_songweiwei_P.C12_Gain_p *
                      Wind_songweiwei_DW.Delay_x2_DSTATE_o[0]) +
      Wind_songweiwei_P.Duk_Gain_j * rtb_Kv1_j_idx_0;
    rtb_Ll_d_idx_1 = (Wind_songweiwei_P.C11_Gain_a *
                      Wind_songweiwei_DW.Delay_x1_DSTATE_m[1] +
                      Wind_songweiwei_P.C12_Gain_p *
                      Wind_songweiwei_DW.Delay_x2_DSTATE_o[1]) +
      Wind_songweiwei_P.Duk_Gain_j * rtb_Kv1_j_idx_1;
    rtb_Ll_d_idx_2 = (Wind_songweiwei_P.C11_Gain_a *
                      Wind_songweiwei_DW.Delay_x1_DSTATE_m[2] +
                      Wind_songweiwei_P.C12_Gain_p *
                      Wind_songweiwei_DW.Delay_x2_DSTATE_o[2]) +
      Wind_songweiwei_P.Duk_Gain_j * rtb_Kv1_j_idx_2;

    /* Gain: '<S263>/Gain3' incorporates:
     *  Gain: '<S263>/Gain1'
     */
    for (i = 0; i < 3; i++) {
      tmp_1[i] = Wind_songweiwei_P.Gain3_Gain_a[i + 6] * rtb_Ll_d_idx_2 +
        (Wind_songweiwei_P.Gain3_Gain_a[i + 3] * rtb_Ll_d_idx_1 +
         Wind_songweiwei_P.Gain3_Gain_a[i] * rtb_Ll_d_idx_0);
    }

    /* End of Gain: '<S263>/Gain3' */

    /* Gain: '<S263>/Gain1' */
    rtb_Gain1_k[0] = Wind_songweiwei_P.Gain1_Gain_k0 * tmp_1[0];
    rtb_Gain1_k[1] = Wind_songweiwei_P.Gain1_Gain_k0 * tmp_1[1];
    rtb_Gain1_k[2] = Wind_songweiwei_P.Gain1_Gain_k0 * tmp_1[2];

    /* Math: '<S225>/Math Function' incorporates:
     *  Constant: '<S225>/Constant4'
     *  DiscreteIntegrator: '<S225>/Discrete-Time Integrator'
     */
    rtb_MathFunction_f = rt_modd_snf
      (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_h,
       Wind_songweiwei_P.Constant4_Value_mc);

    /* RelationalOperator: '<S264>/Compare' incorporates:
     *  Constant: '<S262>/Constant'
     *  Constant: '<S264>/Constant'
     */
    rtb_Compare_m = (uint8_T)(Wind_songweiwei_P.AlphaBetaZerotodq0_Alignment_n ==
      Wind_songweiwei_P.CompareToConstant_const_c);

    /* Outputs for Enabled SubSystem: '<S262>/Subsystem1' */
    Wind_songweiwei_Subsystem1(rtb_Compare_m, &rtb_Gain1_k[0],
      rtb_MathFunction_f, &Wind_songweiwei_B.Subsystem1_l);

    /* End of Outputs for SubSystem: '<S262>/Subsystem1' */

    /* RelationalOperator: '<S265>/Compare' incorporates:
     *  Constant: '<S262>/Constant'
     *  Constant: '<S265>/Constant'
     */
    rtb_Compare_g = (uint8_T)(Wind_songweiwei_P.AlphaBetaZerotodq0_Alignment_n ==
      Wind_songweiwei_P.CompareToConstant1_const_b);

    /* Outputs for Enabled SubSystem: '<S262>/Subsystem - pi//2 delay' */
    Wind_songweiw_Subsystempi2delay(rtb_Compare_g, &rtb_Gain1_k[0],
      rtb_MathFunction_f, &Wind_songweiwei_B.Subsystempi2delay_j);

    /* End of Outputs for SubSystem: '<S262>/Subsystem - pi//2 delay' */

    /* Switch: '<S262>/Switch' */
    if (rtb_Compare_m != 0) {
      rtb_Ll_q_idx_0 = Wind_songweiwei_B.Subsystem1_l.Fcn;
    } else {
      rtb_Ll_q_idx_0 = Wind_songweiwei_B.Subsystempi2delay_j.Fcn;
    }

    if (rtb_Compare_m != 0) {
      rtb_Ll_q_idx_1 = Wind_songweiwei_B.Subsystem1_l.Fcn1;
    } else {
      rtb_Ll_q_idx_1 = Wind_songweiwei_B.Subsystempi2delay_j.Fcn1;
    }

    /* End of Switch: '<S262>/Switch' */

    /* Product: '<S122>/Divide3' incorporates:
     *  Constant: '<S122>/Fnom'
     *  UnitDelay: '<S225>/Unit Delay'
     */
    rtb_Switch1_p = Wind_songweiwei_DW.UnitDelay_DSTATE_p /
      Wind_songweiwei_P.Fnom_Value;

    /* Switch: '<S131>/Switch7' incorporates:
     *  Constant: '<S121>/Constant2'
     *  DataTypeConversion: '<S121>/Data Type Conversion2'
     *  UnitDelay: '<S131>/IC = i_ic'
     */
    if (Wind_songweiwei_P.Constant2_Value_i != 0.0) {
      rtb_Switch_e3 = Wind_songweiwei_P.Constant2_Value_i;
    } else {
      rtb_Switch_e3 = Wind_songweiwei_DW.ICi_ic_DSTATE;
    }

    /* End of Switch: '<S131>/Switch7' */

    /* Sum: '<S128>/Sum6' incorporates:
     *  DiscreteIntegrator: '<S128>/Discrete-Time Integrator'
     *  Gain: '<S128>/Kp4'
     */
    rtb_Lmsatd = Wind_songweiwei_P.DiscretePIController_Kp_h * rtb_phimq +
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_g;

    /* Saturate: '<S128>/Saturation2' */
    if (rtb_Lmsatd > Wind_songweiwei_P.Saturation2_UpperSat_g) {
      rtb_Lmsatd = Wind_songweiwei_P.Saturation2_UpperSat_g;
    } else {
      if (rtb_Lmsatd < Wind_songweiwei_P.Saturation2_LowerSat_m) {
        rtb_Lmsatd = Wind_songweiwei_P.Saturation2_LowerSat_m;
      }
    }

    /* End of Saturate: '<S128>/Saturation2' */

    /* Sum: '<S121>/Sum1' incorporates:
     *  Constant: '<S121>/Imax^2'
     *  Math: '<S121>/Math Function'
     */
    rtb_phimd = Wind_songweiwei_P.Imax2_Value - rtb_Lmsatd * rtb_Lmsatd;

    /* Math: '<S121>/Math Function1'
     *
     * About '<S121>/Math Function1':
     *  Operator: sqrt
     */
    if (rtb_phimd < 0.0) {
      rtb_phimd = -sqrt(fabs(rtb_phimd));
    } else {
      rtb_phimd = sqrt(rtb_phimd);
    }

    /* End of Math: '<S121>/Math Function1' */

    /* Sum: '<S134>/sum1' incorporates:
     *  Gain: '<S134>/C'
     *  Gain: '<S134>/D'
     *  UnitDelay: '<S134>/Delay_x'
     */
    rtb_Add4 = Wind_songweiwei_P.D_Gain_b * rtb_phimd +
      Wind_songweiwei_P.C_Gain_d * Wind_songweiwei_DW.Delay_x_DSTATE_m;

    /* Gain: '<S121>/Gain' */
    rtb_donotdeletethisgain = Wind_songweiwei_P.Gain_Gain_j * rtb_Add4;

    /* Switch: '<S131>/Switch3' incorporates:
     *  Constant: '<S121>/Constant2'
     *  DataTypeConversion: '<S121>/Data Type Conversion2'
     */
    if (Wind_songweiwei_P.Constant2_Value_i != 0.0) {
      rtb_donotdeletethisgain_l = rtb_Switch_e3;
    } else {
      rtb_donotdeletethisgain_l = rtb_donotdeletethisgain;
    }

    /* End of Switch: '<S131>/Switch3' */

    /* MinMax: '<S131>/MinMax1' */
    if (!((rtb_donotdeletethisgain_l >= rtb_donotdeletethisgain) || rtIsNaN
          (rtb_donotdeletethisgain))) {
      rtb_donotdeletethisgain_l = rtb_donotdeletethisgain;
    }

    /* End of MinMax: '<S131>/MinMax1' */

    /* Switch: '<S131>/Switch2' incorporates:
     *  Constant: '<S121>/Constant2'
     *  DataTypeConversion: '<S121>/Data Type Conversion2'
     */
    if (Wind_songweiwei_P.Constant2_Value_i != 0.0) {
      rtb_donotdeletethisgain = rtb_Switch_e3;
    } else {
      rtb_donotdeletethisgain = rtb_Add4;
    }

    /* End of Switch: '<S131>/Switch2' */

    /* MinMax: '<S131>/MinMax' */
    if (!((rtb_donotdeletethisgain <= rtb_Add4) || rtIsNaN(rtb_Add4))) {
      rtb_donotdeletethisgain = rtb_Add4;
    }

    /* End of MinMax: '<S131>/MinMax' */

    /* DiscreteIntegrator: '<S121>/Discrete-Time Integrator1' incorporates:
     *  UnitDelay: '<S121>/Unit Delay1'
     */
    if (Wind_songweiwei_DW.DiscreteTimeIntegrator1_IC_LOAD != 0) {
      Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE =
        Wind_songweiwei_DW.UnitDelay1_DSTATE;
      if (Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE >=
          Wind_songweiwei_P.DiscreteTimeIntegrator1_UpperSa) {
        Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE =
          Wind_songweiwei_P.DiscreteTimeIntegrator1_UpperSa;
      } else {
        if (Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE <=
            Wind_songweiwei_P.DiscreteTimeIntegrator1_LowerSa) {
          Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE =
            Wind_songweiwei_P.DiscreteTimeIntegrator1_LowerSa;
        }
      }
    }

    if (Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE >=
        Wind_songweiwei_P.DiscreteTimeIntegrator1_UpperSa) {
      Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE =
        Wind_songweiwei_P.DiscreteTimeIntegrator1_UpperSa;
    } else {
      if (Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE <=
          Wind_songweiwei_P.DiscreteTimeIntegrator1_LowerSa) {
        Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE =
          Wind_songweiwei_P.DiscreteTimeIntegrator1_LowerSa;
      }
    }

    /* Gain: '<S161>/Gain3' incorporates:
     *  Gain: '<S122>/pu->V'
     *  Gain: '<S161>/Gain1'
     */
    for (i = 0; i < 3; i++) {
      tmp_1[i] = Wind_songweiwei_P.Gain3_Gain_nj[i + 6] *
        (Wind_songweiwei_P.puV_Gain * rtb_Ll_d_idx_2) +
        (Wind_songweiwei_P.Gain3_Gain_nj[i + 3] * (Wind_songweiwei_P.puV_Gain *
          rtb_Ll_d_idx_1) + Wind_songweiwei_P.puV_Gain * rtb_Ll_d_idx_0 *
         Wind_songweiwei_P.Gain3_Gain_nj[i]);
    }

    /* End of Gain: '<S161>/Gain3' */

    /* Gain: '<S161>/Gain1' */
    rtb_Gain1_c[0] = Wind_songweiwei_P.Gain1_Gain_dl * tmp_1[0];
    rtb_Gain1_c[1] = Wind_songweiwei_P.Gain1_Gain_dl * tmp_1[1];
    rtb_Gain1_c[2] = Wind_songweiwei_P.Gain1_Gain_dl * tmp_1[2];

    /* RelationalOperator: '<S162>/Compare' incorporates:
     *  Constant: '<S160>/Constant'
     *  Constant: '<S162>/Constant'
     */
    rtb_Compare_b = (uint8_T)(Wind_songweiwei_P.AlphaBetaZerotodq0_Alignment_iw ==
      Wind_songweiwei_P.CompareToConstant_const_g);

    /* Outputs for Enabled SubSystem: '<S160>/Subsystem1' */
    Wind_songweiwei_Subsystem1(rtb_Compare_b, &rtb_Gain1_c[0],
      rtb_MathFunction_f, &Wind_songweiwei_B.Subsystem1_o);

    /* End of Outputs for SubSystem: '<S160>/Subsystem1' */

    /* RelationalOperator: '<S163>/Compare' incorporates:
     *  Constant: '<S160>/Constant'
     *  Constant: '<S163>/Constant'
     */
    rtb_Compare_p = (uint8_T)(Wind_songweiwei_P.AlphaBetaZerotodq0_Alignment_iw ==
      Wind_songweiwei_P.CompareToConstant1_const_g);

    /* Outputs for Enabled SubSystem: '<S160>/Subsystem - pi//2 delay' */
    Wind_songweiw_Subsystempi2delay(rtb_Compare_p, &rtb_Gain1_c[0],
      rtb_MathFunction_f, &Wind_songweiwei_B.Subsystempi2delay_k);

    /* End of Outputs for SubSystem: '<S160>/Subsystem - pi//2 delay' */

    /* Switch: '<S160>/Switch' */
    if (rtb_Compare_b != 0) {
      rtb_Switch[0] = Wind_songweiwei_B.Subsystem1_o.Fcn;
    } else {
      rtb_Switch[0] = Wind_songweiwei_B.Subsystempi2delay_k.Fcn;
    }

    if (rtb_Compare_b != 0) {
      rtb_Switch[1] = Wind_songweiwei_B.Subsystem1_o.Fcn1;
    } else {
      rtb_Switch[1] = Wind_songweiwei_B.Subsystempi2delay_k.Fcn1;
    }

    /* End of Switch: '<S160>/Switch' */

    /* DiscreteIntegrator: '<S154>/Integ4' */
    if (Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE != 0) {
      Wind_songweiwei_B.Integ4 = Wind_songweiwei_DW.Integ4_DSTATE;
    } else {
      Wind_songweiwei_B.Integ4 = Wind_songweiwei_P.Integ4_gainval_g *
        rtb_Switch[0] + Wind_songweiwei_DW.Integ4_DSTATE;
    }

    /* End of DiscreteIntegrator: '<S154>/Integ4' */

    /* Saturate: '<S154>/To avoid division by zero' incorporates:
     *  UnitDelay: '<S225>/Unit Delay'
     */
    if (Wind_songweiwei_DW.UnitDelay_DSTATE_p >
        Wind_songweiwei_P.Toavoiddivisionbyzero_UpperS_ft) {
      unnamed_idx_1 = Wind_songweiwei_P.Toavoiddivisionbyzero_UpperS_ft;
    } else if (Wind_songweiwei_DW.UnitDelay_DSTATE_p <
               Wind_songweiwei_P.Toavoiddivisionbyzero_LowerSa_b) {
      unnamed_idx_1 = Wind_songweiwei_P.Toavoiddivisionbyzero_LowerSa_b;
    } else {
      unnamed_idx_1 = Wind_songweiwei_DW.UnitDelay_DSTATE_p;
    }

    /* Fcn: '<S154>/Number of samples per cycle' incorporates:
     *  Saturate: '<S154>/To avoid division by zero'
     */
    rtb_Add4 = 1.0 / unnamed_idx_1 / 2.0e-6;

    /* Rounding: '<S154>/Rounding Function' */
    rtb_Ifdsat = ceil(rtb_Add4);

    /* Gain: '<S154>/Gain' */
    Wind_songweiwei_B.Delay = Wind_songweiwei_P.Ts * rtb_Ifdsat;

    /* S-Function block: <S159>/S-Function  */
    {
      int_T indDelayed;
      int_T discreteDelay;

      /* Input present value(s) */
      ((real_T *)Wind_songweiwei_DW.SFunction_PWORK.uBuffers)
        [Wind_songweiwei_DW.SFunction_IWORK.indEnd] = Wind_songweiwei_B.Integ4;

      /* Calculate delayed index */
      discreteDelay =
        (int_T)floor((Wind_songweiwei_B.Delay/2.0E-6) + 0.5);
      if (discreteDelay > Wind_songweiwei_DW.SFunction_IWORK.maxDiscrDelay)
        discreteDelay = Wind_songweiwei_DW.SFunction_IWORK.maxDiscrDelay;
      indDelayed = Wind_songweiwei_DW.SFunction_IWORK.indEnd - ((discreteDelay >
        0) ? discreteDelay : 0);
      if (indDelayed < 0) {
        if (Wind_songweiwei_DW.SFunction_IWORK.indBeg == 0 )
          indDelayed = 0;
        else
          indDelayed += Wind_songweiwei_DW.SFunction_IWORK.bufSz;
      }

      /* Output past value(s) */
      Wind_songweiwei_B.SFunction = ((real_T *)
        Wind_songweiwei_DW.SFunction_PWORK.uBuffers)[indDelayed];
    }

    /* Step: '<S154>/Step' */
    rtb_Kv_ow_idx_2 = Wind_songweiwei_M->Timing.t[1];
    if (rtb_Kv_ow_idx_2 < 1.0 /
        Wind_songweiwei_P.Discrete3phasePLLDrivenPositive) {
      rtb_donotdeletethisgain_o = Wind_songweiwei_P.Step_Y0_c;
    } else {
      rtb_donotdeletethisgain_o = Wind_songweiwei_P.Step_YFinal_g;
    }

    /* End of Step: '<S154>/Step' */

    /* Switch: '<S154>/Switch' incorporates:
     *  Gain: '<S158>/Gain1'
     *  Product: '<S154>/Product'
     *  Product: '<S158>/Product2'
     *  Product: '<S158>/Product4'
     *  Product: '<S158>/Product5'
     *  Sum: '<S154>/Sum5'
     *  Sum: '<S154>/Sum7'
     *  Sum: '<S158>/Sum1'
     *  Sum: '<S158>/Sum4'
     *  UnitDelay: '<S154>/Unit Delay'
     *  UnitDelay: '<S158>/Unit Delay'
     *  UnitDelay: '<S225>/Unit Delay'
     */
    if (rtb_donotdeletethisgain_o >= Wind_songweiwei_P.Switch_Threshold_n) {
      /* Sum: '<S158>/Sum5' */
      rtb_Ifdsat = rtb_Add4 - rtb_Ifdsat;
      rtb_Switch_e = ((rtb_Switch[0] - Wind_songweiwei_DW.UnitDelay_DSTATE_f) *
                      rtb_Ifdsat * Wind_songweiwei_P.Gain1_Gain_g + rtb_Switch[0])
        * (rtb_Ifdsat / rtb_Add4) + (Wind_songweiwei_B.Integ4 -
        Wind_songweiwei_B.SFunction) * Wind_songweiwei_DW.UnitDelay_DSTATE_p;
    } else {
      rtb_Switch_e = Wind_songweiwei_DW.UnitDelay_DSTATE_m;
    }

    /* End of Switch: '<S154>/Switch' */

    /* DiscreteIntegrator: '<S153>/Integ4' */
    if (Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_f != 0) {
      Wind_songweiwei_B.Integ4_n = Wind_songweiwei_DW.Integ4_DSTATE_n;
    } else {
      Wind_songweiwei_B.Integ4_n = Wind_songweiwei_P.Integ4_gainval_e *
        rtb_Switch[1] + Wind_songweiwei_DW.Integ4_DSTATE_n;
    }

    /* End of DiscreteIntegrator: '<S153>/Integ4' */

    /* Saturate: '<S153>/To avoid division by zero' incorporates:
     *  UnitDelay: '<S225>/Unit Delay'
     */
    if (Wind_songweiwei_DW.UnitDelay_DSTATE_p >
        Wind_songweiwei_P.Toavoiddivisionbyzero_UpperSa_p) {
      unnamed_idx_1 = Wind_songweiwei_P.Toavoiddivisionbyzero_UpperSa_p;
    } else if (Wind_songweiwei_DW.UnitDelay_DSTATE_p <
               Wind_songweiwei_P.Toavoiddivisionbyzero_LowerSa_d) {
      unnamed_idx_1 = Wind_songweiwei_P.Toavoiddivisionbyzero_LowerSa_d;
    } else {
      unnamed_idx_1 = Wind_songweiwei_DW.UnitDelay_DSTATE_p;
    }

    /* Fcn: '<S153>/Number of samples per cycle' incorporates:
     *  Saturate: '<S153>/To avoid division by zero'
     */
    rtb_Add4 = 1.0 / unnamed_idx_1 / 2.0e-6;

    /* Rounding: '<S153>/Rounding Function' */
    rtb_Ifdsat = ceil(rtb_Add4);

    /* Gain: '<S153>/Gain' */
    Wind_songweiwei_B.Delay_f = Wind_songweiwei_P.Ts * rtb_Ifdsat;

    /* S-Function block: <S157>/S-Function  */
    {
      int_T indDelayed;
      int_T discreteDelay;

      /* Input present value(s) */
      ((real_T *)Wind_songweiwei_DW.SFunction_PWORK_c.uBuffers)
        [Wind_songweiwei_DW.SFunction_IWORK_f.indEnd] =
        Wind_songweiwei_B.Integ4_n;

      /* Calculate delayed index */
      discreteDelay =
        (int_T)floor((Wind_songweiwei_B.Delay_f/2.0E-6) + 0.5);
      if (discreteDelay > Wind_songweiwei_DW.SFunction_IWORK_f.maxDiscrDelay)
        discreteDelay = Wind_songweiwei_DW.SFunction_IWORK_f.maxDiscrDelay;
      indDelayed = Wind_songweiwei_DW.SFunction_IWORK_f.indEnd - ((discreteDelay
        > 0) ? discreteDelay : 0);
      if (indDelayed < 0) {
        if (Wind_songweiwei_DW.SFunction_IWORK_f.indBeg == 0 )
          indDelayed = 0;
        else
          indDelayed += Wind_songweiwei_DW.SFunction_IWORK_f.bufSz;
      }

      /* Output past value(s) */
      Wind_songweiwei_B.SFunction_j = ((real_T *)
        Wind_songweiwei_DW.SFunction_PWORK_c.uBuffers)[indDelayed];
    }

    /* Step: '<S153>/Step' */
    rtb_Kv_ow_idx_2 = Wind_songweiwei_M->Timing.t[1];
    if (rtb_Kv_ow_idx_2 < 1.0 /
        Wind_songweiwei_P.Discrete3phasePLLDrivenPositive) {
      rtb_donotdeletethisgain_o = Wind_songweiwei_P.Step_Y0_a;
    } else {
      rtb_donotdeletethisgain_o = Wind_songweiwei_P.Step_YFinal_f;
    }

    /* End of Step: '<S153>/Step' */

    /* Switch: '<S153>/Switch' incorporates:
     *  Gain: '<S156>/Gain1'
     *  Product: '<S153>/Product'
     *  Product: '<S156>/Product2'
     *  Product: '<S156>/Product4'
     *  Product: '<S156>/Product5'
     *  Sum: '<S153>/Sum5'
     *  Sum: '<S153>/Sum7'
     *  Sum: '<S156>/Sum1'
     *  Sum: '<S156>/Sum4'
     *  UnitDelay: '<S153>/Unit Delay'
     *  UnitDelay: '<S156>/Unit Delay'
     *  UnitDelay: '<S225>/Unit Delay'
     */
    if (rtb_donotdeletethisgain_o >= Wind_songweiwei_P.Switch_Threshold_m) {
      /* Sum: '<S156>/Sum5' */
      rtb_Ifdsat = rtb_Add4 - rtb_Ifdsat;
      rtb_Switch_d = ((rtb_Switch[1] - Wind_songweiwei_DW.UnitDelay_DSTATE_fd) *
                      rtb_Ifdsat * Wind_songweiwei_P.Gain1_Gain_k + rtb_Switch[1])
        * (rtb_Ifdsat / rtb_Add4) + (Wind_songweiwei_B.Integ4_n -
        Wind_songweiwei_B.SFunction_j) * Wind_songweiwei_DW.UnitDelay_DSTATE_p;
    } else {
      rtb_Switch_d = Wind_songweiwei_DW.UnitDelay_DSTATE_l;
    }

    /* End of Switch: '<S153>/Switch' */

    /* Fcn: '<S152>/x->r' */
    rtb_donotdeletethisgain_o = rt_hypotd_snf(rtb_Switch_e, rtb_Switch_d);

    /* Gain: '<S122>/V->pu' */
    rtb_Vpu = Wind_songweiwei_P.Vpu_Gain * rtb_donotdeletethisgain_o;

    /* Sum: '<S121>/Sum' incorporates:
     *  DiscreteIntegrator: '<S121>/Discrete-Time Integrator1'
     */
    rtb_Sum = Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE - rtb_Vpu;

    /* Sum: '<S131>/Sum' incorporates:
     *  Constant: '<S121>/Constant3'
     *  Gain: '<S131>/Gain'
     *  Product: '<S131>/Product'
     *  Sum: '<S131>/Sum1'
     *  UnitDelay: '<S131>/IC = 0'
     */
    rtb_Switch_e3 += Wind_songweiwei_P.Ts * 0.5 *
      Wind_songweiwei_P.Constant3_Value_i * (rtb_Sum +
      Wind_songweiwei_DW.IC0_DSTATE);

    /* MinMax: '<S135>/MinMax1' */
    if (!((rtb_donotdeletethisgain <= rtb_Switch_e3) || rtIsNaN(rtb_Switch_e3)))
    {
      rtb_donotdeletethisgain = rtb_Switch_e3;
    }

    /* End of MinMax: '<S135>/MinMax1' */

    /* MinMax: '<S135>/MinMax' */
    if ((rtb_donotdeletethisgain_l >= rtb_donotdeletethisgain) || rtIsNaN
        (rtb_donotdeletethisgain)) {
      rtb_MinMax = rtb_donotdeletethisgain_l;
    } else {
      rtb_MinMax = rtb_donotdeletethisgain;
    }

    /* End of MinMax: '<S135>/MinMax' */

    /* Sum: '<S193>/C*X(k)+D*u(k)' incorporates:
     *  Gain: '<S193>/D*u(k)'
     *  Gain: '<S196>/C11'
     *  Gain: '<S196>/C12'
     *  Sum: '<S196>/sum2'
     *  UnitDelay: '<S193>/Delay_x1'
     *  UnitDelay: '<S193>/Delay_x2'
     */
    rtb_Ifdsat = (Wind_songweiwei_P.C11_Gain_h *
                  Wind_songweiwei_DW.Delay_x1_DSTATE_b[0] +
                  Wind_songweiwei_P.C12_Gain_h3 *
                  Wind_songweiwei_DW.Delay_x2_DSTATE_pn[0]) +
      Wind_songweiwei_P.Duk_Gain_f * rtb_Kv_idx_0;
    rtb_Kv_ow_idx_1 = (Wind_songweiwei_P.C11_Gain_h *
                       Wind_songweiwei_DW.Delay_x1_DSTATE_b[1] +
                       Wind_songweiwei_P.C12_Gain_h3 *
                       Wind_songweiwei_DW.Delay_x2_DSTATE_pn[1]) +
      Wind_songweiwei_P.Duk_Gain_f * rtb_Kv_idx_1;
    rtb_Kv_ow_idx_2 = (Wind_songweiwei_P.C11_Gain_h *
                       Wind_songweiwei_DW.Delay_x1_DSTATE_b[2] +
                       Wind_songweiwei_P.C12_Gain_h3 *
                       Wind_songweiwei_DW.Delay_x2_DSTATE_pn[2]) +
      Wind_songweiwei_P.Duk_Gain_f * rtb_Kv_idx_2;

    /* Gain: '<S269>/Gain3' incorporates:
     *  Gain: '<S269>/Gain1'
     */
    for (i = 0; i < 3; i++) {
      tmp_1[i] = Wind_songweiwei_P.Gain3_Gain_a3[i + 6] * rtb_Kv_ow_idx_2 +
        (Wind_songweiwei_P.Gain3_Gain_a3[i + 3] * rtb_Kv_ow_idx_1 +
         Wind_songweiwei_P.Gain3_Gain_a3[i] * rtb_Ifdsat);
    }

    /* End of Gain: '<S269>/Gain3' */

    /* Gain: '<S269>/Gain1' */
    rtb_Gain1_h[0] = Wind_songweiwei_P.Gain1_Gain_db * tmp_1[0];
    rtb_Gain1_h[1] = Wind_songweiwei_P.Gain1_Gain_db * tmp_1[1];
    rtb_Gain1_h[2] = Wind_songweiwei_P.Gain1_Gain_db * tmp_1[2];

    /* RelationalOperator: '<S270>/Compare' incorporates:
     *  Constant: '<S268>/Constant'
     *  Constant: '<S270>/Constant'
     */
    rtb_Compare_o = (uint8_T)(Wind_songweiwei_P.AlphaBetaZerotodq0_Alignment_m ==
      Wind_songweiwei_P.CompareToConstant_const_i);

    /* Outputs for Enabled SubSystem: '<S268>/Subsystem1' */
    Wind_songweiwei_Subsystem1(rtb_Compare_o, &rtb_Gain1_h[0],
      rtb_MathFunction_f, &Wind_songweiwei_B.Subsystem1_p);

    /* End of Outputs for SubSystem: '<S268>/Subsystem1' */

    /* RelationalOperator: '<S271>/Compare' incorporates:
     *  Constant: '<S268>/Constant'
     *  Constant: '<S271>/Constant'
     */
    rtb_Compare_j = (uint8_T)(Wind_songweiwei_P.AlphaBetaZerotodq0_Alignment_m ==
      Wind_songweiwei_P.CompareToConstant1_const_k);

    /* Outputs for Enabled SubSystem: '<S268>/Subsystem - pi//2 delay' */
    Wind_songweiw_Subsystempi2delay(rtb_Compare_j, &rtb_Gain1_h[0],
      rtb_MathFunction_f, &Wind_songweiwei_B.Subsystempi2delay_jr);

    /* End of Outputs for SubSystem: '<S268>/Subsystem - pi//2 delay' */

    /* Switch: '<S268>/Switch' */
    if (rtb_Compare_o != 0) {
      unnamed_idx_1 = Wind_songweiwei_B.Subsystem1_p.Fcn;
    } else {
      unnamed_idx_1 = Wind_songweiwei_B.Subsystempi2delay_jr.Fcn;
    }

    /* Sum: '<S121>/Sum6' */
    rtb_donotdeletethisgain_l = rtb_Lmsatd - unnamed_idx_1;

    /* Switch: '<S268>/Switch' */
    if (rtb_Compare_o != 0) {
      unnamed_idx_1 = Wind_songweiwei_B.Subsystem1_p.Fcn1;
    } else {
      unnamed_idx_1 = Wind_songweiwei_B.Subsystempi2delay_jr.Fcn1;
    }

    /* Sum: '<S121>/Sum4' */
    rtb_donotdeletethisgain = rtb_MinMax - unnamed_idx_1;

    /* Sum: '<S129>/Sum6' incorporates:
     *  DiscreteIntegrator: '<S129>/Discrete-Time Integrator'
     *  Gain: '<S129>/Kp4'
     */
    unnamed_idx_1 = Wind_songweiwei_P.DiscretePIController1_Kp *
      rtb_donotdeletethisgain_l +
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[0];
    rtb_Switch_e3 = Wind_songweiwei_P.DiscretePIController1_Kp *
      rtb_donotdeletethisgain +
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[1];

    /* Saturate: '<S129>/Saturation2' */
    if (unnamed_idx_1 > Wind_songweiwei_P.Saturation2_UpperSat_gz) {
      unnamed_idx_1 = Wind_songweiwei_P.Saturation2_UpperSat_gz;
    } else {
      if (unnamed_idx_1 < Wind_songweiwei_P.Saturation2_LowerSat_a) {
        unnamed_idx_1 = Wind_songweiwei_P.Saturation2_LowerSat_a;
      }
    }

    /* Sum: '<S121>/Sum3' incorporates:
     *  Constant: '<S121>/Constant1'
     *  Constant: '<S121>/Constant5'
     *  Product: '<S121>/Product'
     *  Product: '<S121>/Product2'
     *  Sum: '<S121>/Sum5'
     */
    rtb_Add4 = ((Wind_songweiwei_P.Constant1_Value_f5 * rtb_Switch1_p *
                 rtb_MinMax + rtb_Ll_q_idx_0) -
                Wind_songweiwei_P.Constant5_Value_c * rtb_Lmsatd) -
      unnamed_idx_1;

    /* Saturate: '<S129>/Saturation2' */
    if (rtb_Switch_e3 > Wind_songweiwei_P.Saturation2_UpperSat_gz) {
      rtb_Switch_e3 = Wind_songweiwei_P.Saturation2_UpperSat_gz;
    } else {
      if (rtb_Switch_e3 < Wind_songweiwei_P.Saturation2_LowerSat_a) {
        rtb_Switch_e3 = Wind_songweiwei_P.Saturation2_LowerSat_a;
      }
    }

    /* Sum: '<S121>/Sum9' incorporates:
     *  Constant: '<S121>/Constant4'
     *  Constant: '<S121>/Constant6'
     *  Product: '<S121>/Product1'
     *  Product: '<S121>/Product3'
     *  Sum: '<S121>/Sum8'
     */
    rtb_Switch1_p = ((rtb_Ll_q_idx_1 - rtb_MinMax *
                      Wind_songweiwei_P.Constant6_Value_g) - rtb_Lmsatd *
                     Wind_songweiwei_P.Constant4_Value_de * rtb_Switch1_p) -
      rtb_Switch_e3;

    /* Saturate: '<S121>/Avoid division by zero' */
    if (rtb_Cp > Wind_songweiwei_P.Avoiddivisionbyzero_UpperSat) {
      rtb_Cp = Wind_songweiwei_P.Avoiddivisionbyzero_UpperSat;
    } else {
      if (rtb_Cp < Wind_songweiwei_P.Avoiddivisionbyzero_LowerSat) {
        rtb_Cp = Wind_songweiwei_P.Avoiddivisionbyzero_LowerSat;
      }
    }

    /* Product: '<S121>/Product4' incorporates:
     *  Constant: '<S121>/K'
     *  Fcn: '<S127>/x->r'
     *  Saturate: '<S121>/Avoid division by zero'
     */
    rtb_Switch_e3 = rt_hypotd_snf(rtb_Add4, rtb_Switch1_p) *
      Wind_songweiwei_P.K_Value / rtb_Cp;

    /* Saturate: '<S121>/0-Mod_index_max' */
    if (rtb_Switch_e3 > Wind_songweiwei_P.Mod_index_max_UpperSat) {
      rtb_Switch_e3 = Wind_songweiwei_P.Mod_index_max_UpperSat;
    } else {
      if (rtb_Switch_e3 < Wind_songweiwei_P.Mod_index_max_LowerSat) {
        rtb_Switch_e3 = Wind_songweiwei_P.Mod_index_max_LowerSat;
      }
    }

    /* End of Saturate: '<S121>/0-Mod_index_max' */

    /* Fcn: '<S127>/x->theta' */
    rtb_Add4 = rt_atan2d_snf(rtb_Switch1_p, rtb_Add4);

    /* Gain: '<S128>/Kp5' */
    rtb_Kp5_p = Wind_songweiwei_P.DiscretePIController_Ki_f * rtb_phimq;

    /* Gain: '<S129>/Kp5' */
    rtb_Kp5_a[0] = Wind_songweiwei_P.DiscretePIController1_Ki *
      rtb_donotdeletethisgain_l;
    rtb_Kp5_a[1] = Wind_songweiwei_P.DiscretePIController1_Ki *
      rtb_donotdeletethisgain;

    /* Sum: '<S134>/A*x(k) + B*u(k)' incorporates:
     *  Gain: '<S134>/A'
     *  Gain: '<S134>/B'
     *  UnitDelay: '<S134>/Delay_x'
     */
    rtb_xk1_k = Wind_songweiwei_P.A_Gain * Wind_songweiwei_DW.Delay_x_DSTATE_m +
      Wind_songweiwei_P.B_Gain * rtb_phimd;

    /* Fcn: '<S132>/r->x' */
    rtb_phimd = rtb_Switch_e3 * cos(rtb_Add4);

    /* Fcn: '<S132>/theta->y' */
    rtb_Lmsatd = rtb_Switch_e3 * sin(rtb_Add4);

    /* Product: '<S122>/Divide' incorporates:
     *  Constant: '<S122>/Constant1'
     */
    rtb_donotdeletethisgain_l = rtb_donotdeletethisgain_o /
      Wind_songweiwei_P.Constant1_Value_g;

    /* Math: '<S122>/Math Function1' */
    rtb_donotdeletethisgain_l *= rtb_donotdeletethisgain_l;

    /* Gain: '<S122>/pu->A' */
    rtb_Ifdsat *= Wind_songweiwei_P.puA_Gain;
    rtb_Kv_ow_idx_1 *= Wind_songweiwei_P.puA_Gain;
    rtb_Kv_ow_idx_2 *= Wind_songweiwei_P.puA_Gain;

    /* Gain: '<S175>/Gain3' incorporates:
     *  Gain: '<S175>/Gain1'
     */
    for (i = 0; i < 3; i++) {
      tmp_1[i] = Wind_songweiwei_P.Gain3_Gain_o[i + 6] * rtb_Kv_ow_idx_2 +
        (Wind_songweiwei_P.Gain3_Gain_o[i + 3] * rtb_Kv_ow_idx_1 +
         Wind_songweiwei_P.Gain3_Gain_o[i] * rtb_Ifdsat);
    }

    /* End of Gain: '<S175>/Gain3' */

    /* Gain: '<S175>/Gain1' */
    rtb_Gain1_i[0] = Wind_songweiwei_P.Gain1_Gain_ep * tmp_1[0];
    rtb_Gain1_i[1] = Wind_songweiwei_P.Gain1_Gain_ep * tmp_1[1];
    rtb_Gain1_i[2] = Wind_songweiwei_P.Gain1_Gain_ep * tmp_1[2];

    /* RelationalOperator: '<S176>/Compare' incorporates:
     *  Constant: '<S174>/Constant'
     *  Constant: '<S176>/Constant'
     */
    rtb_Compare_f = (uint8_T)(Wind_songweiwei_P.AlphaBetaZerotodq0_Alignment_l ==
      Wind_songweiwei_P.CompareToConstant_const_gf);

    /* Outputs for Enabled SubSystem: '<S174>/Subsystem1' */
    Wind_songweiwei_Subsystem1(rtb_Compare_f, &rtb_Gain1_i[0],
      rtb_MathFunction_f, &Wind_songweiwei_B.Subsystem1_oh);

    /* End of Outputs for SubSystem: '<S174>/Subsystem1' */

    /* RelationalOperator: '<S177>/Compare' incorporates:
     *  Constant: '<S174>/Constant'
     *  Constant: '<S177>/Constant'
     */
    rtb_Compare_k = (uint8_T)(Wind_songweiwei_P.AlphaBetaZerotodq0_Alignment_l ==
      Wind_songweiwei_P.CompareToConstant1_const_gp);

    /* Outputs for Enabled SubSystem: '<S174>/Subsystem - pi//2 delay' */
    Wind_songweiw_Subsystempi2delay(rtb_Compare_k, &rtb_Gain1_i[0],
      rtb_MathFunction_f, &Wind_songweiwei_B.Subsystempi2delay_h);

    /* End of Outputs for SubSystem: '<S174>/Subsystem - pi//2 delay' */

    /* Switch: '<S174>/Switch' */
    if (rtb_Compare_f != 0) {
      rtb_Switch_o[0] = Wind_songweiwei_B.Subsystem1_oh.Fcn;
    } else {
      rtb_Switch_o[0] = Wind_songweiwei_B.Subsystempi2delay_h.Fcn;
    }

    if (rtb_Compare_f != 0) {
      rtb_Switch_o[1] = Wind_songweiwei_B.Subsystem1_oh.Fcn1;
    } else {
      rtb_Switch_o[1] = Wind_songweiwei_B.Subsystempi2delay_h.Fcn1;
    }

    /* End of Switch: '<S174>/Switch' */

    /* DiscreteIntegrator: '<S168>/Integ4' */
    if (Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_l != 0) {
      Wind_songweiwei_B.Integ4_j = Wind_songweiwei_DW.Integ4_DSTATE_nw;
    } else {
      Wind_songweiwei_B.Integ4_j = Wind_songweiwei_P.Integ4_gainval_l *
        rtb_Switch_o[0] + Wind_songweiwei_DW.Integ4_DSTATE_nw;
    }

    /* End of DiscreteIntegrator: '<S168>/Integ4' */

    /* Saturate: '<S168>/To avoid division by zero' incorporates:
     *  UnitDelay: '<S225>/Unit Delay'
     */
    if (Wind_songweiwei_DW.UnitDelay_DSTATE_p >
        Wind_songweiwei_P.Toavoiddivisionbyzero_UpperSa_e) {
      unnamed_idx_1 = Wind_songweiwei_P.Toavoiddivisionbyzero_UpperSa_e;
    } else if (Wind_songweiwei_DW.UnitDelay_DSTATE_p <
               Wind_songweiwei_P.Toavoiddivisionbyzero_LowerSa_a) {
      unnamed_idx_1 = Wind_songweiwei_P.Toavoiddivisionbyzero_LowerSa_a;
    } else {
      unnamed_idx_1 = Wind_songweiwei_DW.UnitDelay_DSTATE_p;
    }

    /* Fcn: '<S168>/Number of samples per cycle' incorporates:
     *  Saturate: '<S168>/To avoid division by zero'
     */
    rtb_Add4 = 1.0 / unnamed_idx_1 / 2.0e-6;

    /* Rounding: '<S168>/Rounding Function' */
    rtb_Ifdsat = ceil(rtb_Add4);

    /* Gain: '<S168>/Gain' */
    Wind_songweiwei_B.Delay_n = Wind_songweiwei_P.Ts * rtb_Ifdsat;

    /* S-Function block: <S173>/S-Function  */
    {
      int_T indDelayed;
      int_T discreteDelay;

      /* Input present value(s) */
      ((real_T *)Wind_songweiwei_DW.SFunction_PWORK_l.uBuffers)
        [Wind_songweiwei_DW.SFunction_IWORK_c.indEnd] =
        Wind_songweiwei_B.Integ4_j;

      /* Calculate delayed index */
      discreteDelay =
        (int_T)floor((Wind_songweiwei_B.Delay_n/2.0E-6) + 0.5);
      if (discreteDelay > Wind_songweiwei_DW.SFunction_IWORK_c.maxDiscrDelay)
        discreteDelay = Wind_songweiwei_DW.SFunction_IWORK_c.maxDiscrDelay;
      indDelayed = Wind_songweiwei_DW.SFunction_IWORK_c.indEnd - ((discreteDelay
        > 0) ? discreteDelay : 0);
      if (indDelayed < 0) {
        if (Wind_songweiwei_DW.SFunction_IWORK_c.indBeg == 0 )
          indDelayed = 0;
        else
          indDelayed += Wind_songweiwei_DW.SFunction_IWORK_c.bufSz;
      }

      /* Output past value(s) */
      Wind_songweiwei_B.SFunction_f = ((real_T *)
        Wind_songweiwei_DW.SFunction_PWORK_l.uBuffers)[indDelayed];
    }

    /* Step: '<S168>/Step' */
    rtb_Kv_ow_idx_2 = Wind_songweiwei_M->Timing.t[1];
    if (rtb_Kv_ow_idx_2 < 1.0 /
        Wind_songweiwei_P.Discrete3phasePLLDrivenPositive) {
      rtb_Kv_ow_idx_2 = Wind_songweiwei_P.Step_Y0_l;
    } else {
      rtb_Kv_ow_idx_2 = Wind_songweiwei_P.Step_YFinal_j;
    }

    /* End of Step: '<S168>/Step' */

    /* Switch: '<S168>/Switch' incorporates:
     *  Gain: '<S172>/Gain1'
     *  Product: '<S168>/Product'
     *  Product: '<S172>/Product2'
     *  Product: '<S172>/Product4'
     *  Product: '<S172>/Product5'
     *  Sum: '<S168>/Sum5'
     *  Sum: '<S168>/Sum7'
     *  Sum: '<S172>/Sum1'
     *  Sum: '<S172>/Sum4'
     *  UnitDelay: '<S168>/Unit Delay'
     *  UnitDelay: '<S172>/Unit Delay'
     *  UnitDelay: '<S225>/Unit Delay'
     */
    if (rtb_Kv_ow_idx_2 >= Wind_songweiwei_P.Switch_Threshold_nb) {
      /* Sum: '<S172>/Sum5' */
      rtb_Ifdsat = rtb_Add4 - rtb_Ifdsat;
      rtb_Switch_n = ((rtb_Switch_o[0] - Wind_songweiwei_DW.UnitDelay_DSTATE_fz)
                      * rtb_Ifdsat * Wind_songweiwei_P.Gain1_Gain_e +
                      rtb_Switch_o[0]) * (rtb_Ifdsat / rtb_Add4) +
        (Wind_songweiwei_B.Integ4_j - Wind_songweiwei_B.SFunction_f) *
        Wind_songweiwei_DW.UnitDelay_DSTATE_p;
    } else {
      rtb_Switch_n = Wind_songweiwei_DW.UnitDelay_DSTATE_b;
    }

    /* End of Switch: '<S168>/Switch' */

    /* DiscreteIntegrator: '<S167>/Integ4' */
    if (Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_h != 0) {
      Wind_songweiwei_B.Integ4_m = Wind_songweiwei_DW.Integ4_DSTATE_l;
    } else {
      Wind_songweiwei_B.Integ4_m = Wind_songweiwei_P.Integ4_gainval_b *
        rtb_Switch_o[1] + Wind_songweiwei_DW.Integ4_DSTATE_l;
    }

    /* End of DiscreteIntegrator: '<S167>/Integ4' */

    /* Saturate: '<S167>/To avoid division by zero' incorporates:
     *  UnitDelay: '<S225>/Unit Delay'
     */
    if (Wind_songweiwei_DW.UnitDelay_DSTATE_p >
        Wind_songweiwei_P.Toavoiddivisionbyzero_UpperSa_a) {
      unnamed_idx_1 = Wind_songweiwei_P.Toavoiddivisionbyzero_UpperSa_a;
    } else if (Wind_songweiwei_DW.UnitDelay_DSTATE_p <
               Wind_songweiwei_P.Toavoiddivisionbyzero_LowerS_dg) {
      unnamed_idx_1 = Wind_songweiwei_P.Toavoiddivisionbyzero_LowerS_dg;
    } else {
      unnamed_idx_1 = Wind_songweiwei_DW.UnitDelay_DSTATE_p;
    }

    /* Fcn: '<S167>/Number of samples per cycle' incorporates:
     *  Saturate: '<S167>/To avoid division by zero'
     */
    rtb_Add4 = 1.0 / unnamed_idx_1 / 2.0e-6;

    /* Rounding: '<S167>/Rounding Function' */
    rtb_Ifdsat = ceil(rtb_Add4);

    /* Gain: '<S167>/Gain' */
    Wind_songweiwei_B.Delay_h = Wind_songweiwei_P.Ts * rtb_Ifdsat;

    /* S-Function block: <S171>/S-Function  */
    {
      int_T indDelayed;
      int_T discreteDelay;

      /* Input present value(s) */
      ((real_T *)Wind_songweiwei_DW.SFunction_PWORK_p.uBuffers)
        [Wind_songweiwei_DW.SFunction_IWORK_k.indEnd] =
        Wind_songweiwei_B.Integ4_m;

      /* Calculate delayed index */
      discreteDelay =
        (int_T)floor((Wind_songweiwei_B.Delay_h/2.0E-6) + 0.5);
      if (discreteDelay > Wind_songweiwei_DW.SFunction_IWORK_k.maxDiscrDelay)
        discreteDelay = Wind_songweiwei_DW.SFunction_IWORK_k.maxDiscrDelay;
      indDelayed = Wind_songweiwei_DW.SFunction_IWORK_k.indEnd - ((discreteDelay
        > 0) ? discreteDelay : 0);
      if (indDelayed < 0) {
        if (Wind_songweiwei_DW.SFunction_IWORK_k.indBeg == 0 )
          indDelayed = 0;
        else
          indDelayed += Wind_songweiwei_DW.SFunction_IWORK_k.bufSz;
      }

      /* Output past value(s) */
      Wind_songweiwei_B.SFunction_m = ((real_T *)
        Wind_songweiwei_DW.SFunction_PWORK_p.uBuffers)[indDelayed];
    }

    /* Step: '<S167>/Step' */
    rtb_Kv_ow_idx_2 = Wind_songweiwei_M->Timing.t[1];
    if (rtb_Kv_ow_idx_2 < 1.0 /
        Wind_songweiwei_P.Discrete3phasePLLDrivenPositive) {
      rtb_Kv_ow_idx_2 = Wind_songweiwei_P.Step_Y0_o;
    } else {
      rtb_Kv_ow_idx_2 = Wind_songweiwei_P.Step_YFinal_gq;
    }

    /* End of Step: '<S167>/Step' */

    /* Switch: '<S167>/Switch' incorporates:
     *  Gain: '<S170>/Gain1'
     *  Product: '<S167>/Product'
     *  Product: '<S170>/Product2'
     *  Product: '<S170>/Product4'
     *  Product: '<S170>/Product5'
     *  Sum: '<S167>/Sum5'
     *  Sum: '<S167>/Sum7'
     *  Sum: '<S170>/Sum1'
     *  Sum: '<S170>/Sum4'
     *  UnitDelay: '<S167>/Unit Delay'
     *  UnitDelay: '<S170>/Unit Delay'
     *  UnitDelay: '<S225>/Unit Delay'
     */
    if (rtb_Kv_ow_idx_2 >= Wind_songweiwei_P.Switch_Threshold_jy) {
      /* Sum: '<S170>/Sum5' */
      rtb_Ifdsat = rtb_Add4 - rtb_Ifdsat;
      rtb_Switch_h = ((rtb_Switch_o[1] - Wind_songweiwei_DW.UnitDelay_DSTATE_e) *
                      rtb_Ifdsat * Wind_songweiwei_P.Gain1_Gain_m +
                      rtb_Switch_o[1]) * (rtb_Ifdsat / rtb_Add4) +
        (Wind_songweiwei_B.Integ4_m - Wind_songweiwei_B.SFunction_m) *
        Wind_songweiwei_DW.UnitDelay_DSTATE_p;
    } else {
      rtb_Switch_h = Wind_songweiwei_DW.UnitDelay_DSTATE_d;
    }

    /* End of Switch: '<S167>/Switch' */

    /* Product: '<S142>/Product2' incorporates:
     *  Fcn: '<S166>/x->r'
     */
    rtb_donotdeletethisgain_o *= rt_hypotd_snf(rtb_Switch_n, rtb_Switch_h);

    /* Gain: '<S142>/Deg->Rad' incorporates:
     *  Fcn: '<S152>/x->theta'
     *  Fcn: '<S166>/x->theta'
     *  Gain: '<S150>/Rad->Deg.'
     *  Gain: '<S151>/Rad->Deg.'
     *  Sum: '<S142>/Sum'
     */
    rtb_donotdeletethisgain = (Wind_songweiwei_P.RadDeg_Gain * rt_atan2d_snf
      (rtb_Switch_d, rtb_Switch_e) - Wind_songweiwei_P.RadDeg_Gain_j *
      rt_atan2d_snf(rtb_Switch_h, rtb_Switch_n)) * Wind_songweiwei_P.DegRad_Gain;

    /* Product: '<S142>/Product' incorporates:
     *  Trigonometry: '<S142>/Trigonometric Function1'
     */
    rtb_Add4 = rtb_donotdeletethisgain_o * cos(rtb_donotdeletethisgain);

    /* Product: '<S142>/Product1' incorporates:
     *  Trigonometry: '<S142>/Trigonometric Function'
     */
    rtb_donotdeletethisgain_o *= sin(rtb_donotdeletethisgain);

    /* Gain: '<S142>/Gain1' */
    rtb_Ll_q_idx_0 = Wind_songweiwei_P.Gain1_Gain_ki * rtb_Add4;

    /* Gain: '<S122>/var->pu ' incorporates:
     *  Gain: '<S122>/C_var_filter'
     *  Gain: '<S142>/Gain1'
     *  Sum: '<S122>/Sum'
     */
    rtb_varpu = (Wind_songweiwei_P.Gain1_Gain_ki * rtb_donotdeletethisgain_o -
                 Wind_songweiwei_P.C_var_filter_Gain * rtb_donotdeletethisgain_l)
      * Wind_songweiwei_P.varpu_Gain;

    /* Sum: '<S121>/Sum2' incorporates:
     *  Constant: '<Root>/Constant1'
     *  Gain: '<S209>/D*u(k)'
     *  Gain: '<S212>/C11'
     *  Gain: '<S212>/C12'
     *  Sum: '<S209>/C*X(k)+D*u(k)'
     *  Sum: '<S212>/sum2'
     *  UnitDelay: '<S209>/Delay_x1'
     *  UnitDelay: '<S209>/Delay_x2'
     */
    rtb_Sum2_g = ((Wind_songweiwei_P.C11_Gain_j *
                   Wind_songweiwei_DW.Delay_x1_DSTATE_f +
                   Wind_songweiwei_P.C12_Gain_ax *
                   Wind_songweiwei_DW.Delay_x2_DSTATE_d) +
                  Wind_songweiwei_P.Duk_Gain_c *
                  Wind_songweiwei_P.Constant1_Value_fy) - rtb_varpu;

    /* RelationalOperator: '<S138>/Compare' incorporates:
     *  Constant: '<S137>/Constant'
     *  Constant: '<S138>/Constant'
     */
    rtb_Compare_iw = (uint8_T)(Wind_songweiwei_P.dq0toAlphaBetaZero_Alignment ==
      Wind_songweiwei_P.CompareToConstant_const_ia);

    /* Outputs for Enabled SubSystem: '<S137>/Subsystem1' incorporates:
     *  EnablePort: '<S141>/Enable'
     */
    if (rtb_Compare_iw > 0) {
      /* Fcn: '<S141>/Fcn' */
      Wind_songweiwei_B.Fcn = rtb_phimd * cos(rtb_MathFunction_f) - rtb_Lmsatd *
        sin(rtb_MathFunction_f);

      /* Fcn: '<S141>/Fcn1' */
      Wind_songweiwei_B.Fcn1 = rtb_phimd * sin(rtb_MathFunction_f) + rtb_Lmsatd *
        cos(rtb_MathFunction_f);
    }

    /* End of Outputs for SubSystem: '<S137>/Subsystem1' */

    /* Outputs for Enabled SubSystem: '<S137>/Subsystem - pi//2 delay' incorporates:
     *  EnablePort: '<S140>/Enable'
     */
    /* RelationalOperator: '<S139>/Compare' incorporates:
     *  Constant: '<S137>/Constant'
     *  Constant: '<S139>/Constant'
     */
    if ((Wind_songweiwei_P.dq0toAlphaBetaZero_Alignment ==
         Wind_songweiwei_P.CompareToConstant1_const_o0) > 0) {
      /* Fcn: '<S140>/Fcn' */
      Wind_songweiwei_B.Fcn_c = rtb_phimd * sin(rtb_MathFunction_f) + rtb_Lmsatd
        * cos(rtb_MathFunction_f);

      /* Fcn: '<S140>/Fcn1' */
      Wind_songweiwei_B.Fcn1_b = -rtb_phimd * cos(rtb_MathFunction_f) +
        rtb_Lmsatd * sin(rtb_MathFunction_f);
    }

    /* End of RelationalOperator: '<S139>/Compare' */
    /* End of Outputs for SubSystem: '<S137>/Subsystem - pi//2 delay' */

    /* SignalConversion: '<S136>/TmpSignal ConversionAtGain3Inport1' incorporates:
     *  Switch: '<S137>/Switch'
     */
    if (rtb_Compare_iw != 0) {
      rtb_Kv_ow_idx_2 = Wind_songweiwei_B.Fcn;
      unnamed_idx_1 = Wind_songweiwei_B.Fcn1;
    } else {
      rtb_Kv_ow_idx_2 = Wind_songweiwei_B.Fcn_c;
      unnamed_idx_1 = Wind_songweiwei_B.Fcn1_b;
    }

    /* Gain: '<S136>/Gain3' incorporates:
     *  Constant: '<S121>/V0'
     *  SignalConversion: '<S136>/TmpSignal ConversionAtGain3Inport1'
     */
    for (i = 0; i < 3; i++) {
      rtb_Gain3_h[i] = 0.0;
      rtb_Gain3_h[i] += Wind_songweiwei_P.Gain3_Gain_g[i] * rtb_Kv_ow_idx_2;
      rtb_Gain3_h[i] += Wind_songweiwei_P.Gain3_Gain_g[i + 3] * unnamed_idx_1;
      rtb_Gain3_h[i] += Wind_songweiwei_P.Gain3_Gain_g[i + 6] *
        Wind_songweiwei_P.V0_Value;
    }

    /* End of Gain: '<S136>/Gain3' */

    /* Gain: '<S122>/C_var_filter//Q' */
    rtb_Switch1_p = Wind_songweiwei_P.C_var_filterQ_Gain *
      rtb_donotdeletethisgain_l;

    /* Sum: '<S189>/A*x1(k) + B*u1(k) ' incorporates:
     *  Gain: '<S190>/A11'
     *  Gain: '<S190>/A12'
     *  Gain: '<S191>/B11'
     *  Sum: '<S190>/sum2'
     *  UnitDelay: '<S189>/Delay_x1'
     *  UnitDelay: '<S189>/Delay_x2'
     */
    rtb_x1k1[0] = (Wind_songweiwei_P.A11_Gain *
                   Wind_songweiwei_DW.Delay_x1_DSTATE_m[0] +
                   Wind_songweiwei_P.A12_Gain *
                   Wind_songweiwei_DW.Delay_x2_DSTATE_o[0]) +
      Wind_songweiwei_P.B11_Gain * rtb_Kv1_j_idx_0;
    rtb_x1k1[1] = (Wind_songweiwei_P.A11_Gain *
                   Wind_songweiwei_DW.Delay_x1_DSTATE_m[1] +
                   Wind_songweiwei_P.A12_Gain *
                   Wind_songweiwei_DW.Delay_x2_DSTATE_o[1]) +
      Wind_songweiwei_P.B11_Gain * rtb_Kv1_j_idx_1;
    rtb_x1k1[2] = (Wind_songweiwei_P.A11_Gain *
                   Wind_songweiwei_DW.Delay_x1_DSTATE_m[2] +
                   Wind_songweiwei_P.A12_Gain *
                   Wind_songweiwei_DW.Delay_x2_DSTATE_o[2]) +
      Wind_songweiwei_P.B11_Gain * rtb_Kv1_j_idx_2;

    /* Sum: '<S189>/A*x2(k) + B*u2(k)' incorporates:
     *  Gain: '<S190>/A21'
     *  Gain: '<S190>/A22'
     *  Gain: '<S191>/B21'
     *  Sum: '<S190>/sum3'
     *  UnitDelay: '<S189>/Delay_x1'
     *  UnitDelay: '<S189>/Delay_x2'
     */
    rtb_x2k1[0] = (Wind_songweiwei_P.A21_Gain *
                   Wind_songweiwei_DW.Delay_x1_DSTATE_m[0] +
                   Wind_songweiwei_P.A22_Gain *
                   Wind_songweiwei_DW.Delay_x2_DSTATE_o[0]) +
      Wind_songweiwei_P.B21_Gain * rtb_Kv1_j_idx_0;
    rtb_x2k1[1] = (Wind_songweiwei_P.A21_Gain *
                   Wind_songweiwei_DW.Delay_x1_DSTATE_m[1] +
                   Wind_songweiwei_P.A22_Gain *
                   Wind_songweiwei_DW.Delay_x2_DSTATE_o[1]) +
      Wind_songweiwei_P.B21_Gain * rtb_Kv1_j_idx_1;
    rtb_x2k1[2] = (Wind_songweiwei_P.A21_Gain *
                   Wind_songweiwei_DW.Delay_x1_DSTATE_m[2] +
                   Wind_songweiwei_P.A22_Gain *
                   Wind_songweiwei_DW.Delay_x2_DSTATE_o[2]) +
      Wind_songweiwei_P.B21_Gain * rtb_Kv1_j_idx_2;

    /* Sum: '<S193>/A*x1(k) + B*u1(k) ' incorporates:
     *  Gain: '<S194>/A11'
     *  Gain: '<S194>/A12'
     *  Gain: '<S195>/B11'
     *  Sum: '<S194>/sum2'
     *  UnitDelay: '<S193>/Delay_x1'
     *  UnitDelay: '<S193>/Delay_x2'
     */
    rtb_x1k1_n[0] = (Wind_songweiwei_P.A11_Gain_k *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_b[0] +
                     Wind_songweiwei_P.A12_Gain_g *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_pn[0]) +
      Wind_songweiwei_P.B11_Gain_k * rtb_Kv_idx_0;
    rtb_x1k1_n[1] = (Wind_songweiwei_P.A11_Gain_k *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_b[1] +
                     Wind_songweiwei_P.A12_Gain_g *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_pn[1]) +
      Wind_songweiwei_P.B11_Gain_k * rtb_Kv_idx_1;
    rtb_x1k1_n[2] = (Wind_songweiwei_P.A11_Gain_k *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_b[2] +
                     Wind_songweiwei_P.A12_Gain_g *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_pn[2]) +
      Wind_songweiwei_P.B11_Gain_k * rtb_Kv_idx_2;

    /* Sum: '<S193>/A*x2(k) + B*u2(k)' incorporates:
     *  Gain: '<S194>/A21'
     *  Gain: '<S194>/A22'
     *  Gain: '<S195>/B21'
     *  Sum: '<S194>/sum3'
     *  UnitDelay: '<S193>/Delay_x1'
     *  UnitDelay: '<S193>/Delay_x2'
     */
    rtb_x2k1_n[0] = (Wind_songweiwei_P.A21_Gain_l *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_b[0] +
                     Wind_songweiwei_P.A22_Gain_k *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_pn[0]) +
      Wind_songweiwei_P.B21_Gain_m * rtb_Kv_idx_0;
    rtb_x2k1_n[1] = (Wind_songweiwei_P.A21_Gain_l *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_b[1] +
                     Wind_songweiwei_P.A22_Gain_k *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_pn[1]) +
      Wind_songweiwei_P.B21_Gain_m * rtb_Kv_idx_1;
    rtb_x2k1_n[2] = (Wind_songweiwei_P.A21_Gain_l *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_b[2] +
                     Wind_songweiwei_P.A22_Gain_k *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_pn[2]) +
      Wind_songweiwei_P.B21_Gain_m * rtb_Kv_idx_2;

    /* Sum: '<S197>/A*x1(k) + B*u1(k) ' incorporates:
     *  Gain: '<S198>/A11'
     *  Gain: '<S198>/A12'
     *  Gain: '<S199>/B11'
     *  Sum: '<S198>/sum2'
     *  UnitDelay: '<S197>/Delay_x1'
     *  UnitDelay: '<S197>/Delay_x2'
     */
    rtb_x1k1_a[0] = (Wind_songweiwei_P.A11_Gain_j *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_c[0] +
                     Wind_songweiwei_P.A12_Gain_o *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_e[0]) +
      Wind_songweiwei_P.B11_Gain_n * rtb_Kv1_idx_0;
    rtb_x1k1_a[1] = (Wind_songweiwei_P.A11_Gain_j *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_c[1] +
                     Wind_songweiwei_P.A12_Gain_o *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_e[1]) +
      Wind_songweiwei_P.B11_Gain_n * rtb_Kv1_idx_1;
    rtb_x1k1_a[2] = (Wind_songweiwei_P.A11_Gain_j *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_c[2] +
                     Wind_songweiwei_P.A12_Gain_o *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_e[2]) +
      Wind_songweiwei_P.B11_Gain_n * rtb_Kv1_idx_2;

    /* Sum: '<S197>/A*x2(k) + B*u2(k)' incorporates:
     *  Gain: '<S198>/A21'
     *  Gain: '<S198>/A22'
     *  Gain: '<S199>/B21'
     *  Sum: '<S198>/sum3'
     *  UnitDelay: '<S197>/Delay_x1'
     *  UnitDelay: '<S197>/Delay_x2'
     */
    rtb_x2k1_a[0] = (Wind_songweiwei_P.A21_Gain_j *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_c[0] +
                     Wind_songweiwei_P.A22_Gain_i *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_e[0]) +
      Wind_songweiwei_P.B21_Gain_n * rtb_Kv1_idx_0;
    rtb_x2k1_a[1] = (Wind_songweiwei_P.A21_Gain_j *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_c[1] +
                     Wind_songweiwei_P.A22_Gain_i *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_e[1]) +
      Wind_songweiwei_P.B21_Gain_n * rtb_Kv1_idx_1;
    rtb_x2k1_a[2] = (Wind_songweiwei_P.A21_Gain_j *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_c[2] +
                     Wind_songweiwei_P.A22_Gain_i *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_e[2]) +
      Wind_songweiwei_P.B21_Gain_n * rtb_Kv1_idx_2;

    /* Sum: '<S201>/A*x1(k) + B*u1(k) ' incorporates:
     *  Gain: '<S202>/A11'
     *  Gain: '<S202>/A12'
     *  Gain: '<S203>/B11'
     *  Sum: '<S202>/sum2'
     *  UnitDelay: '<S201>/Delay_x1'
     *  UnitDelay: '<S201>/Delay_x2'
     */
    rtb_x1k1_p[0] = (Wind_songweiwei_P.A11_Gain_d *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_i[0] +
                     Wind_songweiwei_P.A12_Gain_b *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_m[0]) +
      Wind_songweiwei_P.B11_Gain_m * rtb_Kv_o_idx_0;
    rtb_x1k1_p[1] = (Wind_songweiwei_P.A11_Gain_d *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_i[1] +
                     Wind_songweiwei_P.A12_Gain_b *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_m[1]) +
      Wind_songweiwei_P.B11_Gain_m * rtb_Kv_o_idx_1;
    rtb_x1k1_p[2] = (Wind_songweiwei_P.A11_Gain_d *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_i[2] +
                     Wind_songweiwei_P.A12_Gain_b *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_m[2]) +
      Wind_songweiwei_P.B11_Gain_m * rtb_Kv_o_idx_2;

    /* Sum: '<S201>/A*x2(k) + B*u2(k)' incorporates:
     *  Gain: '<S202>/A21'
     *  Gain: '<S202>/A22'
     *  Gain: '<S203>/B21'
     *  Sum: '<S202>/sum3'
     *  UnitDelay: '<S201>/Delay_x1'
     *  UnitDelay: '<S201>/Delay_x2'
     */
    rtb_x2k1_i[0] = (Wind_songweiwei_P.A21_Gain_a *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_i[0] +
                     Wind_songweiwei_P.A22_Gain_k5 *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_m[0]) +
      Wind_songweiwei_P.B21_Gain_o * rtb_Kv_o_idx_0;
    rtb_x2k1_i[1] = (Wind_songweiwei_P.A21_Gain_a *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_i[1] +
                     Wind_songweiwei_P.A22_Gain_k5 *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_m[1]) +
      Wind_songweiwei_P.B21_Gain_o * rtb_Kv_o_idx_1;
    rtb_x2k1_i[2] = (Wind_songweiwei_P.A21_Gain_a *
                     Wind_songweiwei_DW.Delay_x1_DSTATE_i[2] +
                     Wind_songweiwei_P.A22_Gain_k5 *
                     Wind_songweiwei_DW.Delay_x2_DSTATE_m[2]) +
      Wind_songweiwei_P.B21_Gain_o * rtb_Kv_o_idx_2;

    /* Sum: '<S205>/A*x1(k) + B*u1(k) ' incorporates:
     *  Gain: '<S206>/A11'
     *  Gain: '<S206>/A12'
     *  Gain: '<S207>/B11'
     *  Sum: '<S206>/sum2'
     *  UnitDelay: '<S205>/Delay_x1'
     *  UnitDelay: '<S205>/Delay_x2'
     */
    rtb_x1k1_ag = (Wind_songweiwei_P.A11_Gain_b *
                   Wind_songweiwei_DW.Delay_x1_DSTATE_in +
                   Wind_songweiwei_P.A12_Gain_gm *
                   Wind_songweiwei_DW.Delay_x2_DSTATE_a) +
      Wind_songweiwei_P.B11_Gain_me * rtb_pu_a;

    /* Sum: '<S205>/A*x2(k) + B*u2(k)' incorporates:
     *  Gain: '<S206>/A21'
     *  Gain: '<S206>/A22'
     *  Gain: '<S207>/B21'
     *  Sum: '<S206>/sum3'
     *  UnitDelay: '<S205>/Delay_x1'
     *  UnitDelay: '<S205>/Delay_x2'
     */
    rtb_x2k1_il = (Wind_songweiwei_P.A21_Gain_lt *
                   Wind_songweiwei_DW.Delay_x1_DSTATE_in +
                   Wind_songweiwei_P.A22_Gain_o *
                   Wind_songweiwei_DW.Delay_x2_DSTATE_a) +
      Wind_songweiwei_P.B21_Gain_f * rtb_pu_a;

    /* Sum: '<S217>/A*x1(k) + B*u1(k) ' incorporates:
     *  Gain: '<S218>/A11'
     *  Gain: '<S218>/A12'
     *  Gain: '<S219>/B11'
     *  Sum: '<S218>/sum2'
     *  UnitDelay: '<S217>/Delay_x1'
     *  UnitDelay: '<S217>/Delay_x2'
     */
    rtb_x1k1_al = (Wind_songweiwei_P.A11_Gain_jo *
                   Wind_songweiwei_DW.Delay_x1_DSTATE_ck +
                   Wind_songweiwei_P.A12_Gain_c *
                   Wind_songweiwei_DW.Delay_x2_DSTATE_ee) +
      Wind_songweiwei_P.B11_Gain_h * rtb_Gain2_m_idx_0;

    /* Sum: '<S217>/A*x2(k) + B*u2(k)' incorporates:
     *  Gain: '<S218>/A21'
     *  Gain: '<S218>/A22'
     *  Gain: '<S219>/B21'
     *  Sum: '<S218>/sum3'
     *  UnitDelay: '<S217>/Delay_x1'
     *  UnitDelay: '<S217>/Delay_x2'
     */
    rtb_x2k1_k = (Wind_songweiwei_P.A21_Gain_e *
                  Wind_songweiwei_DW.Delay_x1_DSTATE_ck +
                  Wind_songweiwei_P.A22_Gain_kz *
                  Wind_songweiwei_DW.Delay_x2_DSTATE_ee) +
      Wind_songweiwei_P.B21_Gain_np * rtb_Gain2_m_idx_0;

    /* Sum: '<S209>/A*x1(k) + B*u1(k) ' incorporates:
     *  Constant: '<Root>/Constant1'
     *  Gain: '<S210>/A11'
     *  Gain: '<S210>/A12'
     *  Gain: '<S211>/B11'
     *  Sum: '<S210>/sum2'
     *  UnitDelay: '<S209>/Delay_x1'
     *  UnitDelay: '<S209>/Delay_x2'
     */
    rtb_x1k1_o = (Wind_songweiwei_P.A11_Gain_h *
                  Wind_songweiwei_DW.Delay_x1_DSTATE_f +
                  Wind_songweiwei_P.A12_Gain_j *
                  Wind_songweiwei_DW.Delay_x2_DSTATE_d) +
      Wind_songweiwei_P.B11_Gain_f * Wind_songweiwei_P.Constant1_Value_fy;

    /* Sum: '<S209>/A*x2(k) + B*u2(k)' incorporates:
     *  Constant: '<Root>/Constant1'
     *  Gain: '<S210>/A21'
     *  Gain: '<S210>/A22'
     *  Gain: '<S211>/B21'
     *  Sum: '<S210>/sum3'
     *  UnitDelay: '<S209>/Delay_x1'
     *  UnitDelay: '<S209>/Delay_x2'
     */
    rtb_x2k1_ac = (Wind_songweiwei_P.A21_Gain_h *
                   Wind_songweiwei_DW.Delay_x1_DSTATE_f +
                   Wind_songweiwei_P.A22_Gain_l *
                   Wind_songweiwei_DW.Delay_x2_DSTATE_d) +
      Wind_songweiwei_P.B21_Gain_nt * Wind_songweiwei_P.Constant1_Value_fy;

    /* Sum: '<S221>/A*x1(k) + B*u1(k) ' incorporates:
     *  Gain: '<S222>/A11'
     *  Gain: '<S222>/A12'
     *  Gain: '<S223>/B11'
     *  Sum: '<S222>/sum2'
     *  UnitDelay: '<S221>/Delay_x1'
     *  UnitDelay: '<S221>/Delay_x2'
     */
    rtb_x1k1_f = (Wind_songweiwei_P.A11_Gain_jk *
                  Wind_songweiwei_DW.Delay_x1_DSTATE +
                  Wind_songweiwei_P.A12_Gain_e *
                  Wind_songweiwei_DW.Delay_x2_DSTATE) +
      Wind_songweiwei_P.B11_Gain_g * rtb_units;

    /* Sum: '<S221>/A*x2(k) + B*u2(k)' incorporates:
     *  Gain: '<S222>/A21'
     *  Gain: '<S222>/A22'
     *  Gain: '<S223>/B21'
     *  Sum: '<S222>/sum3'
     *  UnitDelay: '<S221>/Delay_x1'
     *  UnitDelay: '<S221>/Delay_x2'
     */
    rtb_x2k1_l = (Wind_songweiwei_P.A21_Gain_i *
                  Wind_songweiwei_DW.Delay_x1_DSTATE +
                  Wind_songweiwei_P.A22_Gain_ir *
                  Wind_songweiwei_DW.Delay_x2_DSTATE) +
      Wind_songweiwei_P.B21_Gain_fv * rtb_units;

    /* Sum: '<S213>/A*x1(k) + B*u1(k) ' incorporates:
     *  Gain: '<S214>/A11'
     *  Gain: '<S214>/A12'
     *  Gain: '<S215>/B11'
     *  Sum: '<S214>/sum2'
     *  UnitDelay: '<S213>/Delay_x1'
     *  UnitDelay: '<S213>/Delay_x2'
     */
    rtb_x1k1_d = (Wind_songweiwei_P.A11_Gain_f *
                  Wind_songweiwei_DW.Delay_x1_DSTATE_a +
                  Wind_songweiwei_P.A12_Gain_i *
                  Wind_songweiwei_DW.Delay_x2_DSTATE_p) +
      Wind_songweiwei_P.B11_Gain_kp * rtb_degrd;

    /* Sum: '<S213>/A*x2(k) + B*u2(k)' incorporates:
     *  Gain: '<S214>/A21'
     *  Gain: '<S214>/A22'
     *  Gain: '<S215>/B21'
     *  Sum: '<S214>/sum3'
     *  UnitDelay: '<S213>/Delay_x1'
     *  UnitDelay: '<S213>/Delay_x2'
     */
    rtb_x2k1_h = (Wind_songweiwei_P.A21_Gain_p *
                  Wind_songweiwei_DW.Delay_x1_DSTATE_a +
                  Wind_songweiwei_P.A22_Gain_c *
                  Wind_songweiwei_DW.Delay_x2_DSTATE_p) +
      Wind_songweiwei_P.B21_Gain_l * rtb_degrd;

    /* Outputs for Enabled SubSystem: '<S225>/Automatic Gain Control' incorporates:
     *  EnablePort: '<S226>/Enable'
     */
    /* Constant: '<S225>/Constant1' */
    if (Wind_songweiwei_P.Constant1_Value_o > 0.0) {
      if (!Wind_songweiwei_DW.AutomaticGainControl_MODE) {
        /* Enable for DiscreteIntegrator: '<S235>/Integ4' */
        Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_e = 1U;

        /* Enable for DiscreteIntegrator: '<S238>/Integ4' */
        Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_m = 1U;
        Wind_songweiwei_DW.AutomaticGainControl_MODE = true;
      }

      /* Gain: '<S242>/Gain3' incorporates:
       *  Gain: '<S242>/Gain1'
       */
      for (i = 0; i < 3; i++) {
        tmp_1[i] = Wind_songweiwei_P.Gain3_Gain_l[i + 6] * rtb_Ll_d_idx_2 +
          (Wind_songweiwei_P.Gain3_Gain_l[i + 3] * rtb_Ll_d_idx_1 +
           Wind_songweiwei_P.Gain3_Gain_l[i] * rtb_Ll_d_idx_0);
      }

      /* End of Gain: '<S242>/Gain3' */

      /* Gain: '<S242>/Gain1' */
      rtb_Gain1_dv[0] = Wind_songweiwei_P.Gain1_Gain_ms * tmp_1[0];
      rtb_Gain1_dv[1] = Wind_songweiwei_P.Gain1_Gain_ms * tmp_1[1];
      rtb_Gain1_dv[2] = Wind_songweiwei_P.Gain1_Gain_ms * tmp_1[2];

      /* RelationalOperator: '<S243>/Compare' incorporates:
       *  Constant: '<S241>/Constant'
       *  Constant: '<S243>/Constant'
       */
      rtb_Compare_kz = (uint8_T)(Wind_songweiwei_P.AlphaBetaZerotodq0_Alignment ==
        Wind_songweiwei_P.CompareToConstant_const);

      /* Outputs for Enabled SubSystem: '<S241>/Subsystem1' */
      Wind_songweiwei_Subsystem1(rtb_Compare_kz, &rtb_Gain1_dv[0],
        rtb_MathFunction_f, &Wind_songweiwei_B.Subsystem1_lp);

      /* End of Outputs for SubSystem: '<S241>/Subsystem1' */

      /* RelationalOperator: '<S244>/Compare' incorporates:
       *  Constant: '<S241>/Constant'
       *  Constant: '<S244>/Constant'
       */
      rtb_Compare_od = (uint8_T)(Wind_songweiwei_P.AlphaBetaZerotodq0_Alignment ==
        Wind_songweiwei_P.CompareToConstant1_const);

      /* Outputs for Enabled SubSystem: '<S241>/Subsystem - pi//2 delay' */
      Wind_songweiw_Subsystempi2delay(rtb_Compare_od, &rtb_Gain1_dv[0],
        rtb_MathFunction_f, &Wind_songweiwei_B.Subsystempi2delay_dg);

      /* End of Outputs for SubSystem: '<S241>/Subsystem - pi//2 delay' */

      /* Switch: '<S241>/Switch' */
      if (rtb_Compare_kz != 0) {
        rtb_Switch_m[0] = Wind_songweiwei_B.Subsystem1_lp.Fcn;
      } else {
        rtb_Switch_m[0] = Wind_songweiwei_B.Subsystempi2delay_dg.Fcn;
      }

      if (rtb_Compare_kz != 0) {
        rtb_Switch_m[1] = Wind_songweiwei_B.Subsystem1_lp.Fcn1;
      } else {
        rtb_Switch_m[1] = Wind_songweiwei_B.Subsystempi2delay_dg.Fcn1;
      }

      /* End of Switch: '<S241>/Switch' */

      /* DiscreteIntegrator: '<S235>/Integ4' */
      if (Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_e != 0) {
        Wind_songweiwei_B.Integ4_e = Wind_songweiwei_DW.Integ4_DSTATE_nn;
      } else {
        Wind_songweiwei_B.Integ4_e = Wind_songweiwei_P.Integ4_gainval *
          rtb_Switch_m[0] + Wind_songweiwei_DW.Integ4_DSTATE_nn;
      }

      /* End of DiscreteIntegrator: '<S235>/Integ4' */

      /* Saturate: '<S235>/To avoid division  by zero' incorporates:
       *  UnitDelay: '<S225>/Unit Delay'
       */
      if (Wind_songweiwei_DW.UnitDelay_DSTATE_p >
          Wind_songweiwei_P.Toavoiddivisionbyzero_UpperSat) {
        unnamed_idx_1 = Wind_songweiwei_P.Toavoiddivisionbyzero_UpperSat;
      } else if (Wind_songweiwei_DW.UnitDelay_DSTATE_p <
                 Wind_songweiwei_P.Toavoiddivisionbyzero_LowerSat) {
        unnamed_idx_1 = Wind_songweiwei_P.Toavoiddivisionbyzero_LowerSat;
      } else {
        unnamed_idx_1 = Wind_songweiwei_DW.UnitDelay_DSTATE_p;
      }

      /* Fcn: '<S235>/Number of samples per cycle' incorporates:
       *  Saturate: '<S235>/To avoid division  by zero'
       */
      rtb_Add4 = 1.0 / unnamed_idx_1 / 2.0e-6;

      /* Rounding: '<S235>/Rounding Function' */
      rtb_Ifdsat = ceil(rtb_Add4);

      /* Gain: '<S235>/Gain' */
      Wind_songweiwei_B.Delay_e = Wind_songweiwei_P.Ts * rtb_Ifdsat;

      /* Level2 S-Function Block: '<S237>/S-Function' (sfun_discreteVariableDelay) */
      {
        SimStruct *rts = Wind_songweiwei_M->childSfunctions[0];
        sfcnOutputs(rts, 1);
      }

      /* DigitalClock: '<S235>/Digital  Clock' */
      rtb_ComplextoMagnitudeAngle_o1 = Wind_songweiwei_M->Timing.t[1];

      /* Switch: '<S235>/Switch' incorporates:
       *  Constant: '<S235>/Constant'
       *  Gain: '<S236>/Gain1'
       *  Product: '<S235>/Product'
       *  Product: '<S236>/Product2'
       *  Product: '<S236>/Product4'
       *  Product: '<S236>/Product5'
       *  RelationalOperator: '<S235>/Relational Operator'
       *  Sum: '<S235>/Sum5'
       *  Sum: '<S235>/Sum7'
       *  Sum: '<S236>/Sum1'
       *  Sum: '<S236>/Sum4'
       *  UnitDelay: '<S225>/Unit Delay'
       *  UnitDelay: '<S235>/Unit Delay1'
       *  UnitDelay: '<S236>/Unit Delay'
       */
      if (rtb_ComplextoMagnitudeAngle_o1 >= Wind_songweiwei_P.Constant_Value_g)
      {
        /* Sum: '<S236>/Sum5' */
        rtb_Ifdsat = rtb_Add4 - rtb_Ifdsat;
        rtb_Lmsatd = ((rtb_Switch_m[0] - Wind_songweiwei_DW.UnitDelay_DSTATE_i) *
                      rtb_Ifdsat * Wind_songweiwei_P.Gain1_Gain_me +
                      rtb_Switch_m[0]) * (rtb_Ifdsat / rtb_Add4) +
          (Wind_songweiwei_B.Integ4_e - Wind_songweiwei_B.SFunction_h) *
          Wind_songweiwei_DW.UnitDelay_DSTATE_p;
      } else {
        rtb_Lmsatd = Wind_songweiwei_DW.UnitDelay1_DSTATE_f;
      }

      /* End of Switch: '<S235>/Switch' */

      /* DiscreteIntegrator: '<S238>/Integ4' */
      if (Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_m != 0) {
        Wind_songweiwei_B.Integ4_h = Wind_songweiwei_DW.Integ4_DSTATE_c;
      } else {
        Wind_songweiwei_B.Integ4_h = Wind_songweiwei_P.Integ4_gainval_j *
          rtb_Switch_m[1] + Wind_songweiwei_DW.Integ4_DSTATE_c;
      }

      /* End of DiscreteIntegrator: '<S238>/Integ4' */

      /* Saturate: '<S238>/To avoid division  by zero' incorporates:
       *  UnitDelay: '<S225>/Unit Delay'
       */
      if (Wind_songweiwei_DW.UnitDelay_DSTATE_p >
          Wind_songweiwei_P.Toavoiddivisionbyzero_UpperSa_f) {
        unnamed_idx_1 = Wind_songweiwei_P.Toavoiddivisionbyzero_UpperSa_f;
      } else if (Wind_songweiwei_DW.UnitDelay_DSTATE_p <
                 Wind_songweiwei_P.Toavoiddivisionbyzero_LowerSa_m) {
        unnamed_idx_1 = Wind_songweiwei_P.Toavoiddivisionbyzero_LowerSa_m;
      } else {
        unnamed_idx_1 = Wind_songweiwei_DW.UnitDelay_DSTATE_p;
      }

      /* Fcn: '<S238>/Number of samples per cycle' incorporates:
       *  Saturate: '<S238>/To avoid division  by zero'
       */
      rtb_Add4 = 1.0 / unnamed_idx_1 / 2.0e-6;

      /* Rounding: '<S238>/Rounding Function' */
      rtb_Ifdsat = ceil(rtb_Add4);

      /* Gain: '<S238>/Gain' */
      Wind_songweiwei_B.Delay_l = Wind_songweiwei_P.Ts * rtb_Ifdsat;

      /* Level2 S-Function Block: '<S240>/S-Function' (sfun_discreteVariableDelay) */
      {
        SimStruct *rts = Wind_songweiwei_M->childSfunctions[1];
        sfcnOutputs(rts, 1);
      }

      /* DigitalClock: '<S238>/Digital  Clock' */
      rtb_ComplextoMagnitudeAngle_o1 = Wind_songweiwei_M->Timing.t[1];

      /* Switch: '<S238>/Switch' incorporates:
       *  Constant: '<S238>/Constant'
       *  Gain: '<S239>/Gain1'
       *  Product: '<S238>/Product'
       *  Product: '<S239>/Product2'
       *  Product: '<S239>/Product4'
       *  Product: '<S239>/Product5'
       *  RelationalOperator: '<S238>/Relational Operator'
       *  Sum: '<S238>/Sum5'
       *  Sum: '<S238>/Sum7'
       *  Sum: '<S239>/Sum1'
       *  Sum: '<S239>/Sum4'
       *  UnitDelay: '<S225>/Unit Delay'
       *  UnitDelay: '<S238>/Unit Delay1'
       *  UnitDelay: '<S239>/Unit Delay'
       */
      if (rtb_ComplextoMagnitudeAngle_o1 >= Wind_songweiwei_P.Constant_Value_d)
      {
        /* Sum: '<S239>/Sum5' */
        rtb_Ifdsat = rtb_Add4 - rtb_Ifdsat;
        rtb_Ifdsat = ((rtb_Switch_m[1] - Wind_songweiwei_DW.UnitDelay_DSTATE_g) *
                      rtb_Ifdsat * Wind_songweiwei_P.Gain1_Gain_o +
                      rtb_Switch_m[1]) * (rtb_Ifdsat / rtb_Add4) +
          (Wind_songweiwei_B.Integ4_h - Wind_songweiwei_B.SFunction_i) *
          Wind_songweiwei_DW.UnitDelay_DSTATE_p;
      } else {
        rtb_Ifdsat = Wind_songweiwei_DW.UnitDelay1_DSTATE_c;
      }

      /* End of Switch: '<S238>/Switch' */

      /* ComplexToMagnitudeAngle: '<S231>/Complex to Magnitude-Angle' incorporates:
       *  RealImagToComplex: '<S231>/Real-Imag to Complex'
       */
      rtb_ComplextoMagnitudeAngle_o1 = rt_hypotd_snf(rtb_Lmsatd, rtb_Ifdsat);

      /* Saturate: '<S226>/Saturation' */
      if (rtb_ComplextoMagnitudeAngle_o1 > Wind_songweiwei_P.Saturation_UpperSat)
      {
        unnamed_idx_1 = Wind_songweiwei_P.Saturation_UpperSat;
      } else if (rtb_ComplextoMagnitudeAngle_o1 <
                 Wind_songweiwei_P.Saturation_LowerSat) {
        unnamed_idx_1 = Wind_songweiwei_P.Saturation_LowerSat;
      } else {
        unnamed_idx_1 = rtb_ComplextoMagnitudeAngle_o1;
      }

      /* Math: '<S226>/Math Function' incorporates:
       *  Saturate: '<S226>/Saturation'
       *
       * About '<S226>/Math Function':
       *  Operator: reciprocal
       */
      Wind_songweiwei_B.MathFunction = 1.0 / unnamed_idx_1;

      /* Update for DiscreteIntegrator: '<S235>/Integ4' */
      Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_e = 0U;
      Wind_songweiwei_DW.Integ4_DSTATE_nn = Wind_songweiwei_P.Integ4_gainval *
        rtb_Switch_m[0] + Wind_songweiwei_B.Integ4_e;

      /* Level2 S-Function Block: '<S237>/S-Function' (sfun_discreteVariableDelay) */
      {
        SimStruct *rts = Wind_songweiwei_M->childSfunctions[0];
        sfcnUpdate(rts, 1);
        if (ssGetErrorStatus(rts) != (NULL))
          return;
      }

      /* Update for UnitDelay: '<S236>/Unit Delay' */
      Wind_songweiwei_DW.UnitDelay_DSTATE_i = rtb_Switch_m[0];

      /* Update for UnitDelay: '<S235>/Unit Delay1' */
      Wind_songweiwei_DW.UnitDelay1_DSTATE_f = rtb_Lmsatd;

      /* Update for DiscreteIntegrator: '<S238>/Integ4' */
      Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_m = 0U;
      Wind_songweiwei_DW.Integ4_DSTATE_c = Wind_songweiwei_P.Integ4_gainval_j *
        rtb_Switch_m[1] + Wind_songweiwei_B.Integ4_h;

      /* Level2 S-Function Block: '<S240>/S-Function' (sfun_discreteVariableDelay) */
      {
        SimStruct *rts = Wind_songweiwei_M->childSfunctions[1];
        sfcnUpdate(rts, 1);
        if (ssGetErrorStatus(rts) != (NULL))
          return;
      }

      /* Update for UnitDelay: '<S239>/Unit Delay' */
      Wind_songweiwei_DW.UnitDelay_DSTATE_g = rtb_Switch_m[1];

      /* Update for UnitDelay: '<S238>/Unit Delay1' */
      Wind_songweiwei_DW.UnitDelay1_DSTATE_c = rtb_Ifdsat;
    } else {
      if (Wind_songweiwei_DW.AutomaticGainControl_MODE) {
        Wind_songweiwei_DW.AutomaticGainControl_MODE = false;
      }
    }

    /* End of Constant: '<S225>/Constant1' */
    /* End of Outputs for SubSystem: '<S225>/Automatic Gain Control' */

    /* Gain: '<S255>/Gain3' incorporates:
     *  Gain: '<S255>/Gain1'
     */
    for (i = 0; i < 3; i++) {
      tmp_1[i] = Wind_songweiwei_P.Gain3_Gain_c[i + 6] * rtb_Ll_d_idx_2 +
        (Wind_songweiwei_P.Gain3_Gain_c[i + 3] * rtb_Ll_d_idx_1 +
         Wind_songweiwei_P.Gain3_Gain_c[i] * rtb_Ll_d_idx_0);
    }

    /* End of Gain: '<S255>/Gain3' */

    /* Gain: '<S255>/Gain1' */
    rtb_Gain1_m[0] = Wind_songweiwei_P.Gain1_Gain_b * tmp_1[0];
    rtb_Gain1_m[1] = Wind_songweiwei_P.Gain1_Gain_b * tmp_1[1];
    rtb_Gain1_m[2] = Wind_songweiwei_P.Gain1_Gain_b * tmp_1[2];

    /* RelationalOperator: '<S256>/Compare' incorporates:
     *  Constant: '<S254>/Constant'
     *  Constant: '<S256>/Constant'
     */
    rtb_Compare_b0 = (uint8_T)(Wind_songweiwei_P.AlphaBetaZerotodq0_Alignment_f ==
      Wind_songweiwei_P.CompareToConstant_const_oc);

    /* Outputs for Enabled SubSystem: '<S254>/Subsystem1' */
    Wind_songweiwei_Subsystem1(rtb_Compare_b0, &rtb_Gain1_m[0],
      rtb_MathFunction_f, &Wind_songweiwei_B.Subsystem1_l1);

    /* End of Outputs for SubSystem: '<S254>/Subsystem1' */

    /* RelationalOperator: '<S257>/Compare' incorporates:
     *  Constant: '<S254>/Constant'
     *  Constant: '<S257>/Constant'
     */
    rtb_Compare_c = (uint8_T)(Wind_songweiwei_P.AlphaBetaZerotodq0_Alignment_f ==
      Wind_songweiwei_P.CompareToConstant1_const_j);

    /* Outputs for Enabled SubSystem: '<S254>/Subsystem - pi//2 delay' */
    Wind_songweiw_Subsystempi2delay(rtb_Compare_c, &rtb_Gain1_m[0],
      rtb_MathFunction_f, &Wind_songweiwei_B.Subsystempi2delay_p);

    /* End of Outputs for SubSystem: '<S254>/Subsystem - pi//2 delay' */

    /* Switch: '<S254>/Switch' */
    if (rtb_Compare_b0 != 0) {
      rtb_Switch_di[0] = Wind_songweiwei_B.Subsystem1_l1.Fcn;
    } else {
      rtb_Switch_di[0] = Wind_songweiwei_B.Subsystempi2delay_p.Fcn;
    }

    if (rtb_Compare_b0 != 0) {
      rtb_Switch_di[1] = Wind_songweiwei_B.Subsystem1_l1.Fcn1;
    } else {
      rtb_Switch_di[1] = Wind_songweiwei_B.Subsystempi2delay_p.Fcn1;
    }

    /* End of Switch: '<S254>/Switch' */

    /* DiscreteIntegrator: '<S251>/Integ4' */
    if (Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_a != 0) {
      Wind_songweiwei_B.Integ4_o = Wind_songweiwei_DW.Integ4_DSTATE_a;
    } else {
      Wind_songweiwei_B.Integ4_o = Wind_songweiwei_P.Integ4_gainval_d *
        rtb_Switch_di[1] + Wind_songweiwei_DW.Integ4_DSTATE_a;
    }

    /* End of DiscreteIntegrator: '<S251>/Integ4' */

    /* Saturate: '<S251>/To avoid division  by zero' incorporates:
     *  UnitDelay: '<S225>/Unit Delay'
     */
    if (Wind_songweiwei_DW.UnitDelay_DSTATE_p >
        Wind_songweiwei_P.Toavoiddivisionbyzero_UpperS_en) {
      unnamed_idx_1 = Wind_songweiwei_P.Toavoiddivisionbyzero_UpperS_en;
    } else if (Wind_songweiwei_DW.UnitDelay_DSTATE_p <
               Wind_songweiwei_P.Toavoiddivisionbyzero_LowerSa_n) {
      unnamed_idx_1 = Wind_songweiwei_P.Toavoiddivisionbyzero_LowerSa_n;
    } else {
      unnamed_idx_1 = Wind_songweiwei_DW.UnitDelay_DSTATE_p;
    }

    /* Fcn: '<S251>/Number of samples per cycle' incorporates:
     *  Saturate: '<S251>/To avoid division  by zero'
     */
    rtb_Add4 = 1.0 / unnamed_idx_1 / 2.0e-6;

    /* Rounding: '<S251>/Rounding Function' */
    rtb_Ifdsat = ceil(rtb_Add4);

    /* Gain: '<S251>/Gain' */
    Wind_songweiwei_B.Delay_o = Wind_songweiwei_P.Ts * rtb_Ifdsat;

    /* Level2 S-Function Block: '<S253>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Wind_songweiwei_M->childSfunctions[2];
      sfcnOutputs(rts, 1);
    }

    /* DigitalClock: '<S251>/Digital  Clock' */
    rtb_donotdeletethisgain_o = Wind_songweiwei_M->Timing.t[1];

    /* Switch: '<S251>/Switch' incorporates:
     *  Constant: '<S251>/Constant'
     *  Gain: '<S252>/Gain1'
     *  Product: '<S251>/Product'
     *  Product: '<S252>/Product2'
     *  Product: '<S252>/Product4'
     *  Product: '<S252>/Product5'
     *  RelationalOperator: '<S251>/Relational Operator'
     *  Sum: '<S251>/Sum5'
     *  Sum: '<S251>/Sum7'
     *  Sum: '<S252>/Sum1'
     *  Sum: '<S252>/Sum4'
     *  UnitDelay: '<S225>/Unit Delay'
     *  UnitDelay: '<S251>/Unit Delay1'
     *  UnitDelay: '<S252>/Unit Delay'
     */
    if (rtb_donotdeletethisgain_o >= Wind_songweiwei_P.Constant_Value_e) {
      /* Sum: '<S252>/Sum5' */
      rtb_Ifdsat = rtb_Add4 - rtb_Ifdsat;
      rtb_Switch_g = ((rtb_Switch_di[1] - Wind_songweiwei_DW.UnitDelay_DSTATE_n)
                      * rtb_Ifdsat * Wind_songweiwei_P.Gain1_Gain_et +
                      rtb_Switch_di[1]) * (rtb_Ifdsat / rtb_Add4) +
        (Wind_songweiwei_B.Integ4_o - Wind_songweiwei_B.SFunction_k) *
        Wind_songweiwei_DW.UnitDelay_DSTATE_p;
    } else {
      rtb_Switch_g = Wind_songweiwei_DW.UnitDelay1_DSTATE_p;
    }

    /* End of Switch: '<S251>/Switch' */

    /* Product: '<S225>/Divide' */
    rtb_Divide = rtb_Switch_g * Wind_songweiwei_B.MathFunction;

    /* DiscreteTransferFcn: '<S227>/Discrete Derivative ' */
    Wind_songweiwei_DW.DiscreteDerivative_tmp = (rtb_Divide -
      Wind_songweiwei_P.DiscreteDerivative_DenCoef[1] *
      Wind_songweiwei_DW.DiscreteDerivative_states) /
      Wind_songweiwei_P.DiscreteDerivative_DenCoef[0];
    rtb_donotdeletethisgain_o = Wind_songweiwei_P.Discrete_Kd *
      Wind_songweiwei_DW.DiscreteDerivative_tmp + -Wind_songweiwei_P.Discrete_Kd
      * Wind_songweiwei_DW.DiscreteDerivative_states;

    /* DiscreteIntegrator: '<S227>/Discrete-Time Integrator' */
    rtb_donotdeletethisgain_l =
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_k;

    /* Sum: '<S227>/Sum6' incorporates:
     *  Gain: '<S227>/Kp4'
     */
    rtb_donotdeletethisgain_l = (Wind_songweiwei_P.Discrete_Kp * rtb_Divide +
      rtb_donotdeletethisgain_l) + rtb_donotdeletethisgain_o;

    /* Saturate: '<S227>/Saturation1' */
    if (rtb_donotdeletethisgain_l > Wind_songweiwei_P.Saturation1_UpperSat_h) {
      rtb_Saturation1 = Wind_songweiwei_P.Saturation1_UpperSat_h;
    } else if (rtb_donotdeletethisgain_l <
               Wind_songweiwei_P.Saturation1_LowerSat_a) {
      rtb_Saturation1 = Wind_songweiwei_P.Saturation1_LowerSat_a;
    } else {
      rtb_Saturation1 = rtb_donotdeletethisgain_l;
    }

    /* End of Saturate: '<S227>/Saturation1' */

    /* Gain: '<S225>/Gain10' */
    rtb_donotdeletethisgain_o = Wind_songweiwei_P.Gain10_Gain * rtb_Saturation1;

    /* RateLimiter: '<S225>/Rate Limiter' */
    rtb_phimd = rtb_donotdeletethisgain_o - Wind_songweiwei_DW.PrevY;
    if (rtb_phimd > Wind_songweiwei_P.RateLimiter_RisingLim) {
      rtb_donotdeletethisgain_o = Wind_songweiwei_DW.PrevY +
        Wind_songweiwei_P.RateLimiter_RisingLim;
    } else {
      if (rtb_phimd < Wind_songweiwei_P.RateLimiter_FallingLim) {
        rtb_donotdeletethisgain_o = Wind_songweiwei_DW.PrevY +
          Wind_songweiwei_P.RateLimiter_FallingLim;
      }
    }

    Wind_songweiwei_DW.PrevY = rtb_donotdeletethisgain_o;

    /* End of RateLimiter: '<S225>/Rate Limiter' */

    /* Sum: '<S247>/A*x1(k) + B*u1(k) ' incorporates:
     *  Gain: '<S248>/A11'
     *  Gain: '<S248>/A12'
     *  Gain: '<S249>/B11'
     *  Sum: '<S248>/sum2'
     *  UnitDelay: '<S247>/Delay_x1'
     *  UnitDelay: '<S247>/Delay_x2'
     */
    rtb_x1k1_oi = (Wind_songweiwei_P.A11_Gain_n *
                   Wind_songweiwei_DW.Delay_x1_DSTATE_bs +
                   Wind_songweiwei_P.A12_Gain_a *
                   Wind_songweiwei_DW.Delay_x2_DSTATE_oc) +
      Wind_songweiwei_P.B11_Gain_e * rtb_donotdeletethisgain_o;

    /* Sum: '<S247>/A*x2(k) + B*u2(k)' incorporates:
     *  Gain: '<S248>/A21'
     *  Gain: '<S248>/A22'
     *  Gain: '<S249>/B21'
     *  Sum: '<S248>/sum3'
     *  UnitDelay: '<S247>/Delay_x1'
     *  UnitDelay: '<S247>/Delay_x2'
     */
    rtb_x2k1_hi = (Wind_songweiwei_P.A21_Gain_c *
                   Wind_songweiwei_DW.Delay_x1_DSTATE_bs +
                   Wind_songweiwei_P.A22_Gain_g *
                   Wind_songweiwei_DW.Delay_x2_DSTATE_oc) +
      Wind_songweiwei_P.B21_Gain_a * rtb_donotdeletethisgain_o;

    /* Sum: '<S247>/C*X(k)+D*u(k)' incorporates:
     *  Gain: '<S247>/D*u(k)'
     *  Gain: '<S250>/C11'
     *  Gain: '<S250>/C12'
     *  Sum: '<S250>/sum2'
     *  UnitDelay: '<S247>/Delay_x1'
     *  UnitDelay: '<S247>/Delay_x2'
     */
    rtb_yk_c = (Wind_songweiwei_P.C11_Gain_g *
                Wind_songweiwei_DW.Delay_x1_DSTATE_bs +
                Wind_songweiwei_P.C12_Gain_jf *
                Wind_songweiwei_DW.Delay_x2_DSTATE_oc) +
      Wind_songweiwei_P.Duk_Gain_p * rtb_donotdeletethisgain_o;

    /* Sum: '<S261>/A*x(k) + B*u(k)' incorporates:
     *  Fcn: '<S145>/Magnitude of flux'
     *  Gain: '<S261>/A'
     *  Gain: '<S261>/B'
     *  UnitDelay: '<S261>/Delay_x'
     */
    rtb_xk1_g = Wind_songweiwei_P.A_Gain_h * Wind_songweiwei_DW.Delay_x_DSTATE +
      Wind_songweiwei_P.B_Gain_k * tmp;

    /* Gain: '<S122>/W->pu' incorporates:
     *  Sum: '<S122>/Sum1'
     */
    rtb_Wpu = (rtb_Ll_q_idx_0 + rtb_Switch1_p) * Wind_songweiwei_P.Wpu_Gain;

    /* Rounding: '<S124>/Rounding Function1' incorporates:
     *  Constant: '<S124>/Nb_of_samples_by_cycle'
     */
    rtb_donotdeletethisgain_o = rt_roundd_snf(1.0 /
      Wind_songweiwei_P.PulseGenerator_Freq_sawtooth / Wind_songweiwei_P.Ts);

    /* Product: '<S124>/Divide3' incorporates:
     *  Constant: '<S124>/Ts'
     */
    rtb_donotdeletethisgain_o *= Wind_songweiwei_P.Ts;

    /* DigitalClock: '<S124>/Digital Clock' */
    rtb_donotdeletethisgain_l = Wind_songweiwei_M->Timing.t[1];

    /* Math: '<S124>/Math Function' */
    rtb_donotdeletethisgain_l = rt_modd_snf(rtb_donotdeletethisgain_l,
      rtb_donotdeletethisgain_o);

    /* Switch: '<S124>/Switch' incorporates:
     *  Constant: '<S124>/Constant3'
     *  Constant: '<S124>/Constant4'
     *  Constant: '<S124>/Ts'
     *  Product: '<S124>/Divide2'
     *  RelationalOperator: '<S124>/Relational Operator1'
     *  Sum: '<S124>/Add'
     *  UnitDelay: '<S67>/Unit Delay2'
     */
    if ((rtb_donotdeletethisgain_o - Wind_songweiwei_P.Ts) *
        Wind_songweiwei_DW.UnitDelay2_DSTATE > rtb_donotdeletethisgain_l) {
      Wind_songweiwei_B.Switch = Wind_songweiwei_P.Constant3_Value_m;
    } else {
      Wind_songweiwei_B.Switch = Wind_songweiwei_P.Constant4_Value_c;
    }

    /* End of Switch: '<S124>/Switch' */

    /* Saturate: '<S125>/0-inf' */
    if (rtb_Wpu > Wind_songweiwei_P.inf_UpperSat) {
      tmp = Wind_songweiwei_P.inf_UpperSat;
    } else if (rtb_Wpu < Wind_songweiwei_P.inf_LowerSat) {
      tmp = Wind_songweiwei_P.inf_LowerSat;
    } else {
      tmp = rtb_Wpu;
    }

    /* Gain: '<S125>/pu_elec->pu_mec' incorporates:
     *  Saturate: '<S125>/0-inf'
     */
    rtb_Ifdsat = Wind_songweiwei_P.pu_elecpu_mec_Gain * tmp;

    /* Switch: '<S125>/Switch' incorporates:
     *  Constant: '<S125>/wref'
     *  Fcn: '<S125>/Tracking Curve'
     */
    if (rtb_Ifdsat >= Wind_songweiwei_P.Switch_Threshold_ny) {
      rtb_donotdeletethisgain_o = Wind_songweiwei_P.wref_Value;
    } else {
      rtb_donotdeletethisgain_o = (-0.5551 * rt_powd_snf(rtb_Ifdsat, 2.0) +
        1.183 * rtb_Ifdsat) + 0.425;
    }

    /* End of Switch: '<S125>/Switch' */

    /* Sum: '<S310>/sum1' incorporates:
     *  Gain: '<S310>/C'
     *  Gain: '<S310>/D'
     *  UnitDelay: '<S310>/Delay_x'
     */
    rtb_donotdeletethisgain_l = Wind_songweiwei_P.D_Gain_i *
      rtb_donotdeletethisgain_o + Wind_songweiwei_P.C_Gain_g *
      Wind_songweiwei_DW.Delay_x_DSTATE_o;

    /* Sum: '<S125>/Sum2' incorporates:
     *  Constant: '<S125>/Constant2'
     */
    rtb_donotdeletethisgain = rtb_Ifdsat - Wind_songweiwei_P.Constant2_Value_e;

    /* Gain: '<S125>/Gain2' incorporates:
     *  Sum: '<S125>/Sum'
     */
    unnamed_idx_1 = (rtb_yk - rtb_donotdeletethisgain_l) *
      Wind_songweiwei_P.Gain2_Gain_n;

    /* Sum: '<S305>/Sum6' incorporates:
     *  DiscreteIntegrator: '<S305>/Discrete-Time Integrator'
     *  Gain: '<S305>/Kp4'
     */
    rtb_Switch_e3 = Wind_songweiwei_P.DiscretePIController1_Kp_p *
      rtb_donotdeletethisgain +
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_o;

    /* Saturate: '<S125>/0-pitch_max' */
    if (unnamed_idx_1 > Wind_songweiwei_P.pitch_max_UpperSat) {
      unnamed_idx_1 = Wind_songweiwei_P.pitch_max_UpperSat;
    } else {
      if (unnamed_idx_1 < Wind_songweiwei_P.pitch_max_LowerSat) {
        unnamed_idx_1 = Wind_songweiwei_P.pitch_max_LowerSat;
      }
    }

    /* Saturate: '<S305>/Saturation2' */
    if (rtb_Switch_e3 > Wind_songweiwei_P.Saturation2_UpperSat_f) {
      rtb_Switch_e3 = Wind_songweiwei_P.Saturation2_UpperSat_f;
    } else {
      if (rtb_Switch_e3 < Wind_songweiwei_P.Saturation2_LowerSat_m2) {
        rtb_Switch_e3 = Wind_songweiwei_P.Saturation2_LowerSat_m2;
      }
    }

    /* Sum: '<S125>/Sum3' incorporates:
     *  Saturate: '<S125>/0-pitch_max'
     *  Saturate: '<S305>/Saturation2'
     */
    rtb_Lmsatd = unnamed_idx_1 + rtb_Switch_e3;

    /* RateLimiter: '<S125>/Rate Limiter   1' */
    if (Wind_songweiwei_DW.LastMajorTime == (rtInf)) {
      Wind_songweiwei_B.RateLimiter1 = rtb_Lmsatd;
    } else {
      rtb_Ifdsat = Wind_songweiwei_M->Timing.t[0] -
        Wind_songweiwei_DW.LastMajorTime;
      rtb_Add4 = rtb_Ifdsat * Wind_songweiwei_P.RateLimiter1_RisingLim;
      rtb_phimd = rtb_Lmsatd - Wind_songweiwei_DW.PrevY_m;
      if (rtb_phimd > rtb_Add4) {
        Wind_songweiwei_B.RateLimiter1 = Wind_songweiwei_DW.PrevY_m + rtb_Add4;
      } else {
        rtb_Ifdsat *= Wind_songweiwei_P.RateLimiter1_FallingLim;
        if (rtb_phimd < rtb_Ifdsat) {
          Wind_songweiwei_B.RateLimiter1 = Wind_songweiwei_DW.PrevY_m +
            rtb_Ifdsat;
        } else {
          Wind_songweiwei_B.RateLimiter1 = rtb_Lmsatd;
        }
      }
    }

    /* End of RateLimiter: '<S125>/Rate Limiter   1' */

    /* Saturate: '<S125>/0-pitch_max ' */
    if (Wind_songweiwei_B.RateLimiter1 > Wind_songweiwei_P.pitch_max_UpperSat_a)
    {
      rtb_Ifdsat = Wind_songweiwei_P.pitch_max_UpperSat_a;
    } else if (Wind_songweiwei_B.RateLimiter1 <
               Wind_songweiwei_P.pitch_max_LowerSat_n) {
      rtb_Ifdsat = Wind_songweiwei_P.pitch_max_LowerSat_n;
    } else {
      rtb_Ifdsat = Wind_songweiwei_B.RateLimiter1;
    }

    /* End of Saturate: '<S125>/0-pitch_max ' */

    /* Sum: '<S125>/Sum1' */
    rtb_donotdeletethisgain_l = rtb_yk - rtb_donotdeletethisgain_l;

    /* Gain: '<S304>/Kp5' */
    rtb_Kp5_i = Wind_songweiwei_P.DiscretePIController_Ki_p *
      rtb_donotdeletethisgain_l;

    /* Sum: '<S304>/Sum6' incorporates:
     *  DiscreteIntegrator: '<S304>/Discrete-Time Integrator'
     *  Gain: '<S304>/Kp4'
     */
    rtb_Add4 = Wind_songweiwei_P.DiscretePIController_Kp_c *
      rtb_donotdeletethisgain_l +
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_kk;

    /* Gain: '<S305>/Kp5' */
    rtb_Kp5_l = Wind_songweiwei_P.DiscretePIController1_Ki_j *
      rtb_donotdeletethisgain;

    /* DiscreteIntegrator: '<S306>/Discrete-Time Integrator' */
    rtb_donotdeletethisgain_l =
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_b;

    /* Saturate: '<S304>/Saturation2' */
    if (rtb_Add4 > Wind_songweiwei_P.Saturation2_UpperSat_b) {
      rtb_Add4 = Wind_songweiwei_P.Saturation2_UpperSat_b;
    } else {
      if (rtb_Add4 < Wind_songweiwei_P.Saturation2_LowerSat_d) {
        rtb_Add4 = Wind_songweiwei_P.Saturation2_LowerSat_d;
      }
    }

    /* Sum: '<S125>/Sum4' incorporates:
     *  Gain: '<S205>/D*u(k)'
     *  Gain: '<S208>/C11'
     *  Gain: '<S208>/C12'
     *  Saturate: '<S304>/Saturation2'
     *  Sum: '<S205>/C*X(k)+D*u(k)'
     *  Sum: '<S208>/sum2'
     *  UnitDelay: '<S205>/Delay_x1'
     *  UnitDelay: '<S205>/Delay_x2'
     */
    rtb_donotdeletethisgain = rtb_Add4 - ((Wind_songweiwei_P.C11_Gain_ai *
      Wind_songweiwei_DW.Delay_x1_DSTATE_in + Wind_songweiwei_P.C12_Gain_e *
      Wind_songweiwei_DW.Delay_x2_DSTATE_a) + Wind_songweiwei_P.Duk_Gain_k *
      rtb_pu_a);

    /* Gain: '<S306>/Kp5' */
    rtb_Kp5_aq = Wind_songweiwei_P.DiscretePIController2_Ki *
      rtb_donotdeletethisgain;

    /* Sum: '<S306>/Sum6' incorporates:
     *  Gain: '<S306>/Kp4'
     */
    rtb_donotdeletethisgain_l += Wind_songweiwei_P.DiscretePIController2_Kp *
      rtb_donotdeletethisgain;

    /* Saturate: '<S306>/Saturation2' */
    if (rtb_donotdeletethisgain_l > Wind_songweiwei_P.Saturation2_UpperSat_j) {
      rtb_Saturation2_e = Wind_songweiwei_P.Saturation2_UpperSat_j;
    } else if (rtb_donotdeletethisgain_l <
               Wind_songweiwei_P.Saturation2_LowerSat_aj) {
      rtb_Saturation2_e = Wind_songweiwei_P.Saturation2_LowerSat_aj;
    } else {
      rtb_Saturation2_e = rtb_donotdeletethisgain_l;
    }

    /* End of Saturate: '<S306>/Saturation2' */

    /* Sum: '<S309>/A*x(k) + B*u(k)' incorporates:
     *  Gain: '<S309>/A'
     *  Gain: '<S309>/B'
     *  UnitDelay: '<S309>/Delay_x'
     */
    rtb_xk1_a = Wind_songweiwei_P.A_Gain_m * Wind_songweiwei_DW.Delay_x_DSTATE_i
      + Wind_songweiwei_P.B_Gain_b * rtb_Ifdsat;

    /* Sum: '<S309>/sum1' incorporates:
     *  Gain: '<S309>/C'
     *  Gain: '<S309>/D'
     *  UnitDelay: '<S309>/Delay_x'
     */
    rtb_yk_a3 = Wind_songweiwei_P.D_Gain_p * rtb_Ifdsat +
      Wind_songweiwei_P.C_Gain_a * Wind_songweiwei_DW.Delay_x_DSTATE_i;

    /* Sum: '<S310>/A*x(k) + B*u(k)' incorporates:
     *  Gain: '<S310>/A'
     *  Gain: '<S310>/B'
     *  UnitDelay: '<S310>/Delay_x'
     */
    rtb_xk1_go = Wind_songweiwei_P.A_Gain_e *
      Wind_songweiwei_DW.Delay_x_DSTATE_o + Wind_songweiwei_P.B_Gain_h *
      rtb_donotdeletethisgain_o;

    /* UnitDelay: '<S67>/Unit Delay1' */
    rtb_donotdeletethisgain_o = Wind_songweiwei_DW.UnitDelay1_DSTATE_n;

    /* Saturate: '<S80>/Avoid div. by zero' incorporates:
     *  Constant: '<Root>/Constant'
     */
    if (Wind_songweiwei_P.Constant_Value_c5 >
        Wind_songweiwei_P.Avoiddivbyzero_UpperSat) {
      rtb_donotdeletethisgain_l = Wind_songweiwei_P.Avoiddivbyzero_UpperSat;
    } else if (Wind_songweiwei_P.Constant_Value_c5 <
               Wind_songweiwei_P.Avoiddivbyzero_LowerSat) {
      rtb_donotdeletethisgain_l = Wind_songweiwei_P.Avoiddivbyzero_LowerSat;
    } else {
      rtb_donotdeletethisgain_l = Wind_songweiwei_P.Constant_Value_c5;
    }

    /* End of Saturate: '<S80>/Avoid div. by zero' */

    /* Saturate: '<S80>/Avoid div. by zero ' incorporates:
     *  DiscreteIntegrator: '<S69>/Discrete-Time Integrator'
     */
    if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_bh >
        Wind_songweiwei_P.Avoiddivbyzero_UpperSat_l) {
      rtb_Add4 = Wind_songweiwei_P.Avoiddivbyzero_UpperSat_l;
    } else if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_bh <
               Wind_songweiwei_P.Avoiddivbyzero_LowerSat_j) {
      rtb_Add4 = Wind_songweiwei_P.Avoiddivbyzero_LowerSat_j;
    } else {
      rtb_Add4 = Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_bh;
    }

    /* End of Saturate: '<S80>/Avoid div. by zero ' */

    /* Fcn: '<S80>/Fcn1' */
    rtb_Cp = 109.4445 * rtb_Add4 / rtb_donotdeletethisgain_l;

    /* Fcn: '<S80>/Fcn2' */
    rtb_Cp = (((1.0 / ((2.5 + rtb_donotdeletethisgain_o) * 0.08 + rtb_Cp) -
                0.035 / (rt_powd_snf(2.5 + rtb_donotdeletethisgain_o, 3.0) + 1.0))
               * 116.0 + (-5.0 - (2.5 + rtb_donotdeletethisgain_o) * 0.4)) / exp
              ((1.0 / ((2.5 + rtb_donotdeletethisgain_o) * 0.08 + rtb_Cp) -
                0.035 / (rt_powd_snf(2.5 + rtb_donotdeletethisgain_o, 3.0) + 1.0))
               * 21.0) + 0.00912 * rtb_Cp) * 0.645;

    /* Fcn: '<S80>/K2*Cp*wind^3 ' */
    rtb_donotdeletethisgain_l = 11269.722013523666 * rt_powd_snf
      (rtb_donotdeletethisgain_l, 3.0) * rtb_Cp;

    /* DiscreteIntegrator: '<S69>/Discrete-Time Integrator1' */
    rtb_donotdeletethisgain_o =
      Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTAT_l;

    /* Sum: '<S69>/Sum1' incorporates:
     *  DiscreteIntegrator: '<S69>/Discrete-Time Integrator'
     */
    rtb_donotdeletethisgain = Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_bh
      - rtb_units;

    /* Sum: '<S69>/Sum2' incorporates:
     *  Gain: '<S69>/Mutual damping'
     *  Gain: '<S69>/Stiffness'
     */
    rtb_donotdeletethisgain_o = Wind_songweiwei_P.DriveTrain_Ksh *
      rtb_donotdeletethisgain_o + Wind_songweiwei_P.DriveTrain_D_mutual *
      rtb_donotdeletethisgain;

    /* Gain: '<S69>/1_2H_WT' incorporates:
     *  Gain: '<S80>/->pu'
     *  Product: '<S80>/Product2'
     *  Sum: '<S69>/Sum'
     */
    rtb__2H_WT = (1.0 / Wind_songweiwei_P.WindTurbine_Prated *
                  rtb_donotdeletethisgain_l / rtb_Add4 -
                  rtb_donotdeletethisgain_o) * (1.0 / (2.0 *
      Wind_songweiwei_P.DriveTrain_H_WT));

    /* Gain: '<S69>/wbase' */
    rtb_wbase = Wind_songweiwei_P.DriveTrain_wbase * rtb_donotdeletethisgain;

    /* Gain: '<S8>/Power base for the Generator' */
    rtb_Add4 = Wind_songweiwei_P.PowerbasefortheGenerator_Gain *
      rtb_donotdeletethisgain_o;

    /* Gain: '<S327>/N' incorporates:
     *  UnitDelay: '<S67>/Unit Delay4'
     */
    rtb_Vfd = Wind_songweiwei_P.N_Gain * Wind_songweiwei_DW.UnitDelay4_DSTATE;

    /* Gain: '<S332>/1_Vb' */
    rtb_Ll_q_idx_0 = Wind_songweiwei_P._Vb_Gain *
      Wind_songweiwei_B.StateSpace_o1[19];
    rtb_Ll_q_idx_1 = Wind_songweiwei_P._Vb_Gain *
      Wind_songweiwei_B.StateSpace_o1[20];

    /* Fcn: '<S332>/Fcn2' */
    rtb_Fcn2 = ((2.0 * rtb_Ll_q_idx_0 + rtb_Ll_q_idx_1) * rtb_ElementaryMath1 +
                1.7320508075688772 * rtb_Ll_q_idx_1 * rtb_ElementaryMath) *
      0.33333333333333331;

    /* Fcn: '<S332>/Fcn3' */
    rtb_Fcn3 = ((2.0 * rtb_Ll_q_idx_0 + rtb_Ll_q_idx_1) * rtb_ElementaryMath +
                -1.7320508075688772 * rtb_Ll_q_idx_1 * rtb_ElementaryMath1) *
      0.33333333333333331;

    /* DigitalClock: '<S342>/Digital Clock' */
    rtb_donotdeletethisgain_o = Wind_songweiwei_M->Timing.t[1];

    /* Switch: '<S342>/IC' incorporates:
     *  UnitDelay: '<S342>/fluxes'
     */
    if (rtb_donotdeletethisgain_o >= Wind_songweiwei_P.IC_Threshold) {
      /* Assignment: '<S339>/W(1,2)=wr' incorporates:
       *  Constant: '<S339>/u1'
       */
      memcpy(&rtb_Sum2_d[0], &Wind_songweiwei_P.u1_Value_f[0], 25U * sizeof
             (real_T));
      rtb_Sum2_d[5] = rtb_n;

      /* Assignment: '<S339>/W(2,1)=-wr' incorporates:
       *  Gain: '<S339>/Gain1'
       */
      rtb_Sum2_d[1] = Wind_songweiwei_P.Gain1_Gain_l * rtb_n;

      /* Switch: '<S331>/Switch1' incorporates:
       *  Constant: '<S331>/Constant2'
       */
      rtb_RelationalOperator2_idx_0 = (Wind_songweiwei_P.Constant2_Value_l >=
        Wind_songweiwei_P.Switch1_Threshold);
      for (i = 0; i < 25; i++) {
        /* Gain: '<S350>/wbase*Ts//2' incorporates:
         *  Constant: '<S331>/Constant6'
         *  Sum: '<S331>/Sum1'
         *  Switch: '<S331>/Switch1'
         */
        if (rtb_RelationalOperator2_idx_0) {
          tmp = Wind_songweiwei_B.RLinv[i];
        } else {
          tmp = Wind_songweiwei_P.Constant6_Value_d[i];
        }

        rtb_Kv_ow_idx_2 = ((0.0 - rtb_Sum2_d[i]) - tmp) *
          Wind_songweiwei_P.wbaseTs2_Gain;

        /* Sum: '<S350>/Sum1' incorporates:
         *  Constant: '<S350>/u5'
         */
        rtb_Sum2_d_0[i] = Wind_songweiwei_P.u5_Value_j[i] - rtb_Kv_ow_idx_2;

        /* Sum: '<S350>/Sum5' incorporates:
         *  Constant: '<S350>/u5'
         */
        rtb_Kv_ow_idx_2 += Wind_songweiwei_P.u5_Value_j[i];
        rtb_Sum2_d[i] = rtb_Kv_ow_idx_2;
      }

      /* Product: '<S350>/inversion' incorporates:
       *  Gain: '<S350>/wbase*Ts//2'
       *  Sum: '<S331>/Sum1'
       */
      rt_invd5x5_snf(rtb_Sum2_d_0, rtb_inversion);

      /* Product: '<S350>/Product4' incorporates:
       *  Product: '<S342>/Product2'
       */
      for (i = 0; i < 5; i++) {
        for (i_1 = 0; i_1 < 5; i_1++) {
          rtb_Sum2_d_0[i + 5 * i_1] = 0.0;
          for (i_0 = 0; i_0 < 5; i_0++) {
            rtb_Sum2_d_0[i + 5 * i_1] += rtb_inversion[5 * i_0 + i] *
              rtb_Sum2_d[5 * i_1 + i_0];
          }
        }
      }

      /* End of Product: '<S350>/Product4' */

      /* Sum: '<S342>/sum' incorporates:
       *  Constant: '<S327>/[ Vkd =0 Vkq1=0  Vkq2=0 ]'
       */
      tmp_2[0] = rtb_Fcn2;
      tmp_2[1] = rtb_Fcn3;
      tmp_2[2] = rtb_Vfd;
      tmp_2[3] = Wind_songweiwei_P.Vkd0Vkq10Vkq20_Value[0];
      tmp_2[4] = Wind_songweiwei_P.Vkd0Vkq10Vkq20_Value[1];
      for (i = 0; i < 5; i++) {
        /* Sum: '<S342>/sum' incorporates:
         *  Product: '<S342>/Product1'
         *  UnitDelay: '<S342>/voltages'
         */
        tmp_0[i] = tmp_2[i] + Wind_songweiwei_DW.voltages_DSTATE[i];

        /* Product: '<S342>/Product2' incorporates:
         *  Sum: '<S342>/Ad*x(k-1) + Bd*( u(k-1) + u(k))'
         *  UnitDelay: '<S342>/fluxes'
         */
        rtb_inversion_0[i] = 0.0;
        for (i_1 = 0; i_1 < 5; i_1++) {
          rtb_inversion_0[i] += rtb_Sum2_d_0[5 * i_1 + i] *
            Wind_songweiwei_DW.fluxes_DSTATE[i_1];
        }
      }

      /* Product: '<S342>/Product1' incorporates:
       *  Gain: '<S350>/wbase*Ts//2 '
       *  Sum: '<S342>/Ad*x(k-1) + Bd*( u(k-1) + u(k))'
       */
      for (i = 0; i < 5; i++) {
        tmp_2[i] = 0.0;
        for (i_1 = 0; i_1 < 5; i_1++) {
          tmp_2[i] += rtb_inversion[5 * i_1 + i] *
            Wind_songweiwei_P.wbaseTs2_Gain_f * tmp_0[i_1];
        }

        rtb_IC[i] = rtb_inversion_0[i] + tmp_2[i];
      }
    } else {
      for (i = 0; i < 5; i++) {
        rtb_IC[i] = Wind_songweiwei_DW.fluxes_DSTATE[i];
      }
    }

    /* End of Switch: '<S342>/IC' */

    /* Gain: '<S329>/1_Pb' incorporates:
     *  Product: '<S8>/Product'
     */
    rtb_donotdeletethisgain_o = rtb_units * rtb_Add4 *
      Wind_songweiwei_P._Pb_Gain;

    /* Fcn: '<S329>/div' */
    rtb_donotdeletethisgain_l = rtb_donotdeletethisgain_o / rtb_n;

    /* Gain: '<S329>/1 ----- 2H' incorporates:
     *  Gain: '<S329>/F'
     *  Gain: '<S338>/1-1'
     *  Product: '<S338>/Mult1'
     *  Sum: '<S329>/Sum2'
     *  Sum: '<S338>/Sum2'
     *  UnitDelay: '<S342>/fluxes'
     */
    rtb_uH = ((rtb_donotdeletethisgain_l - (Wind_songweiwei_P.u_Gain[0] *
                rtb_changeIqIdcurrentsigns[0] *
                Wind_songweiwei_DW.fluxes_DSTATE[1] + Wind_songweiwei_P.u_Gain[1]
                * rtb_changeIqIdcurrentsigns[1] *
                Wind_songweiwei_DW.fluxes_DSTATE[0])) - Wind_songweiwei_P.F_Gain
              * rtb_n) * Wind_songweiwei_P.uH_Gain;

    /* DiscreteIntegrator: '<S329>/Rotor speed deviation (dw)' */
    if (Wind_songweiwei_DW.Rotorspeeddeviationdw_SYSTEM_EN != 0) {
      Wind_songweiwei_B.dw = Wind_songweiwei_DW.Rotorspeeddeviationdw_DSTATE;
    } else {
      Wind_songweiwei_B.dw = Wind_songweiwei_P.Rotorspeeddeviationdw_gainval *
        rtb_uH + Wind_songweiwei_DW.Rotorspeeddeviationdw_DSTATE;
    }

    /* End of DiscreteIntegrator: '<S329>/Rotor speed deviation (dw)' */

    /* Gain: '<S329>/we base' */
    rtb_webase = Wind_songweiwei_P.webase_Gain * rtb_Sum1_jb;

    /* Gain: '<S329>/web3' */
    rtb_web3 = Wind_songweiwei_P.web3_Gain * rtb_n;

    /* Gain: '<S34>/do not delete this gain' */
    rtb_donotdeletethisgain_o = Wind_songweiwei_P.donotdeletethisgain_Gain_a *
      Wind_songweiwei_B.StateSpace_o1[24];

    /* Gain: '<S35>/do not delete this gain' */
    rtb_donotdeletethisgain_l = Wind_songweiwei_P.donotdeletethisgain_Gain_l *
      Wind_songweiwei_B.StateSpace_o1[25];

    /* Outport: '<Root>/Out_Vabc_B25' incorporates:
     *  Gain: '<S36>/do not delete this gain'
     *  Gain: '<S5>/Kv1'
     */
    Wind_songweiwei_Y.Out_Vabc_B25[0] = Wind_songweiwei_P.Kv1_Gain_p *
      rtb_donotdeletethisgain_o;
    Wind_songweiwei_Y.Out_Vabc_B25[1] = Wind_songweiwei_P.Kv1_Gain_p *
      rtb_donotdeletethisgain_l;
    Wind_songweiwei_Y.Out_Vabc_B25[2] =
      Wind_songweiwei_P.donotdeletethisgain_Gain_h0 *
      Wind_songweiwei_B.StateSpace_o1[26] * Wind_songweiwei_P.Kv1_Gain_p;

    /* Gain: '<S31>/do not delete this gain' */
    rtb_donotdeletethisgain_o = Wind_songweiwei_P.donotdeletethisgain_Gain_mh *
      Wind_songweiwei_B.StateSpace_o1[39];

    /* Gain: '<S32>/do not delete this gain' */
    rtb_donotdeletethisgain_l = Wind_songweiwei_P.donotdeletethisgain_Gain_gl *
      Wind_songweiwei_B.StateSpace_o1[40];

    /* Outport: '<Root>/Out_labc_B25' incorporates:
     *  Gain: '<S33>/do not delete this gain'
     *  Gain: '<S5>/Kv'
     */
    Wind_songweiwei_Y.Out_labc_B25[0] = Wind_songweiwei_P.Kv_Gain_b *
      rtb_donotdeletethisgain_o;
    Wind_songweiwei_Y.Out_labc_B25[1] = Wind_songweiwei_P.Kv_Gain_b *
      rtb_donotdeletethisgain_l;
    Wind_songweiwei_Y.Out_labc_B25[2] =
      Wind_songweiwei_P.donotdeletethisgain_Gain_b *
      Wind_songweiwei_B.StateSpace_o1[41] * Wind_songweiwei_P.Kv_Gain_b;
  }

  {
    int32_T i;

    /* Update for UnitDelay: '<S353>/Unit Delay' */
    for (i = 0; i < 6; i++) {
      Wind_songweiwei_DW.UnitDelay_DSTATE[i] = rtb_Saturation[i];
    }

    /* End of Update for UnitDelay: '<S353>/Unit Delay' */

    /* Update for UnitDelay: '<S67>/Unit Delay6' */
    Wind_songweiwei_DW.UnitDelay6_DSTATE[0] = rtb_Gain3_h[0];
    Wind_songweiwei_DW.UnitDelay6_DSTATE[1] = rtb_Gain3_h[1];
    Wind_songweiwei_DW.UnitDelay6_DSTATE[2] = rtb_Gain3_h[2];

    /* Update for DiscreteIntegrator: '<S329>/Rotor angle dthetae' */
    Wind_songweiwei_DW.Rotorangledthetae_DSTATE +=
      Wind_songweiwei_P.Rotorangledthetae_gainval * rtb_webase;

    /* Update for UnitDelay: '<S342>/fluxes' */
    for (i = 0; i < 5; i++) {
      Wind_songweiwei_DW.fluxes_DSTATE[i] = rtb_IC[i];
    }

    /* End of Update for UnitDelay: '<S342>/fluxes' */

    /* S-Function block: <S357>/State-Space */
    {
      const real_T *As = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.AS;
      const real_T *Bs = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.BS;
      real_T *xtmp = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.XTMP;
      real_T accum;

      /* Calculate new states... */
      {
        int_T i1;
        real_T *xd = &Wind_songweiwei_DW.StateSpace_DSTATE[0];
        for (i1=0; i1 < 45; i1++) {
          accum = 0.0;

          {
            int_T i2;
            real_T *xd = &Wind_songweiwei_DW.StateSpace_DSTATE[0];
            for (i2=0; i2 < 45; i2++) {
              accum += *(As++) * xd[i2];
            }
          }

          {
            int_T i2;
            const real_T *u0 = &Wind_songweiwei_P.SwitchCurrents_Value[0];
            for (i2=0; i2 < 13; i2++) {
              accum += *(Bs++) * u0[i2];
            }

            u0 = &Wind_songweiwei_B.Product[0];
            for (i2=0; i2 < 6; i2++) {
              accum += *(Bs++) * u0[i2];
            }

            accum += *(Bs++) * Wind_songweiwei_B.ib[0];
            accum += *(Bs++) * Wind_songweiwei_B.ib[1];
            accum += *(Bs++) * Wind_songweiwei_B.Sum5[0];
            accum += *(Bs++) * Wind_songweiwei_B.Sum5[1];
            accum += *(Bs++) * Wind_songweiwei_B.Sum5[2];
          }

          xtmp[i1] = accum;
        }
      }

      {
        int_T i1;
        real_T *xd = &Wind_songweiwei_DW.StateSpace_DSTATE[0];
        for (i1=0; i1 < 45; i1++) {
          xd[i1] = xtmp[i1];
        }
      }

      {
        int_T *gState = (int_T*)Wind_songweiwei_DW.StateSpace_PWORK.G_STATE;

        /* Store switch gates values for next step */
        *(gState++) = (int_T) Wind_songweiwei_P.In_ac_switch_Value;
        *(gState++) = (int_T) Wind_songweiwei_P.In_ac_switch_Value;
        *(gState++) = (int_T) Wind_songweiwei_P.In_ac_switch_Value;
        *(gState++) = (int_T) Wind_songweiwei_P.In_dc_switch1_Value;
        *(gState++) = (int_T) Wind_songweiwei_B.Switch;
        *(gState++) = (int_T) Wind_songweiwei_P.In_dc_switch1_Value;
        *(gState++) = (int_T) 0.0;

        {
          int_T i1;
          const real_T *u1 = &Wind_songweiwei_P.g_Value[0];
          for (i1=0; i1 < 6; i1++) {
            *(gState++) = (int_T) u1[i1];
          }

          u1 = &Wind_songweiwei_B.DataTypeConversion[0];
          for (i1=0; i1 < 6; i1++) {
            *(gState++) = (int_T) u1[i1];
          }
        }
      }
    }

    /* Update for UnitDelay: '<S67>/Unit Delay3' */
    Wind_songweiwei_DW.UnitDelay3_DSTATE = rtb_Wpu;

    /* Update for UnitDelay: '<S67>/Unit Delay7' */
    Wind_songweiwei_DW.UnitDelay7_DSTATE = rtb_varpu;

    /* Update for UnitDelay: '<S351>/dw_delay' */
    Wind_songweiwei_DW.dw_delay_DSTATE = Wind_songweiwei_B.dw;

    /* Update for UnitDelay: '<S351>/dw_predict' */
    Wind_songweiwei_DW.dw_predict_DSTATE = rtb_dw_delay;

    /* Update for DiscreteIntegrator: '<S126>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE +=
      Wind_songweiwei_P.DiscreteTimeIntegrator_gainva_a * rtb_Kp5;
    if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE >=
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperSat) {
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE =
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperSat;
    } else {
      if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE <=
          Wind_songweiwei_P.DiscreteTimeIntegrator_LowerSat) {
        Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE =
          Wind_songweiwei_P.DiscreteTimeIntegrator_LowerSat;
      }
    }

    /* End of Update for DiscreteIntegrator: '<S126>/Discrete-Time Integrator' */

    /* Update for UnitDelay: '<S221>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE = rtb_x1k1_f;

    /* Update for UnitDelay: '<S221>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE = rtb_x2k1_l;

    /* Update for UnitDelay: '<S197>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_c[0] = rtb_x1k1_a[0];
    Wind_songweiwei_DW.Delay_x1_DSTATE_c[1] = rtb_x1k1_a[1];
    Wind_songweiwei_DW.Delay_x1_DSTATE_c[2] = rtb_x1k1_a[2];

    /* Update for UnitDelay: '<S197>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_e[0] = rtb_x2k1_a[0];
    Wind_songweiwei_DW.Delay_x2_DSTATE_e[1] = rtb_x2k1_a[1];
    Wind_songweiwei_DW.Delay_x2_DSTATE_e[2] = rtb_x2k1_a[2];

    /* Update for DiscreteIntegrator: '<S329>/theta' */
    Wind_songweiwei_DW.theta_DSTATE += Wind_songweiwei_P.theta_gainval *
      rtb_web3;

    /* Update for UnitDelay: '<S213>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_a = rtb_x1k1_d;

    /* Update for UnitDelay: '<S213>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_p = rtb_x2k1_h;

    /* Update for UnitDelay: '<S201>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_i[0] = rtb_x1k1_p[0];
    Wind_songweiwei_DW.Delay_x1_DSTATE_i[1] = rtb_x1k1_p[1];
    Wind_songweiwei_DW.Delay_x1_DSTATE_i[2] = rtb_x1k1_p[2];

    /* Update for UnitDelay: '<S201>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_m[0] = rtb_x2k1_i[0];
    Wind_songweiwei_DW.Delay_x2_DSTATE_m[1] = rtb_x2k1_i[1];
    Wind_songweiwei_DW.Delay_x2_DSTATE_m[2] = rtb_x2k1_i[2];

    /* Update for UnitDelay: '<S261>/Delay_x' */
    Wind_songweiwei_DW.Delay_x_DSTATE = rtb_xk1_g;

    /* Update for UnitDelay: '<S217>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_ck = rtb_x1k1_al;

    /* Update for UnitDelay: '<S217>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_ee = rtb_x2k1_k;

    /* Update for UnitDelay: '<S189>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_m[0] = rtb_x1k1[0];
    Wind_songweiwei_DW.Delay_x1_DSTATE_m[1] = rtb_x1k1[1];
    Wind_songweiwei_DW.Delay_x1_DSTATE_m[2] = rtb_x1k1[2];

    /* Update for UnitDelay: '<S189>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_o[0] = rtb_x2k1[0];
    Wind_songweiwei_DW.Delay_x2_DSTATE_o[1] = rtb_x2k1[1];
    Wind_songweiwei_DW.Delay_x2_DSTATE_o[2] = rtb_x2k1[2];

    /* Update for DiscreteIntegrator: '<S225>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_h +=
      Wind_songweiwei_P.DiscreteTimeIntegrator_gainva_l * rtb_Saturation1;

    /* Update for UnitDelay: '<S225>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_p = rtb_yk_c;

    /* Update for UnitDelay: '<S131>/IC = i_ic' */
    Wind_songweiwei_DW.ICi_ic_DSTATE = rtb_MinMax;

    /* Update for DiscreteIntegrator: '<S128>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_g +=
      Wind_songweiwei_P.DiscreteTimeIntegrator_gainva_p * rtb_Kp5_p;
    if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_g >=
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_b) {
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_g =
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_b;
    } else {
      if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_g <=
          Wind_songweiwei_P.DiscreteTimeIntegrator_LowerS_p) {
        Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_g =
          Wind_songweiwei_P.DiscreteTimeIntegrator_LowerS_p;
      }
    }

    /* End of Update for DiscreteIntegrator: '<S128>/Discrete-Time Integrator' */

    /* Update for UnitDelay: '<S134>/Delay_x' */
    Wind_songweiwei_DW.Delay_x_DSTATE_m = rtb_xk1_k;

    /* Update for UnitDelay: '<S121>/Unit Delay1' */
    Wind_songweiwei_DW.UnitDelay1_DSTATE = rtb_Vpu;

    /* Update for DiscreteIntegrator: '<S121>/Discrete-Time Integrator1' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator1_IC_LOAD = 0U;
    Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE +=
      Wind_songweiwei_P.DiscreteTimeIntegrator1_gainv_d * rtb_Sum2_g;
    if (Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE >=
        Wind_songweiwei_P.DiscreteTimeIntegrator1_UpperSa) {
      Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE =
        Wind_songweiwei_P.DiscreteTimeIntegrator1_UpperSa;
    } else {
      if (Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE <=
          Wind_songweiwei_P.DiscreteTimeIntegrator1_LowerSa) {
        Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE =
          Wind_songweiwei_P.DiscreteTimeIntegrator1_LowerSa;
      }
    }

    /* End of Update for DiscreteIntegrator: '<S121>/Discrete-Time Integrator1' */

    /* Update for DiscreteIntegrator: '<S154>/Integ4' */
    Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE = 0U;
    Wind_songweiwei_DW.Integ4_DSTATE = Wind_songweiwei_P.Integ4_gainval_g *
      rtb_Switch[0] + Wind_songweiwei_B.Integ4;

    /* S-Function block: <S159>/S-Function  */
    {
      int_T indBeg = Wind_songweiwei_DW.SFunction_IWORK.indBeg;
      int_T indEnd = Wind_songweiwei_DW.SFunction_IWORK.indEnd;
      int_T bufSz = Wind_songweiwei_DW.SFunction_IWORK.bufSz;
      indEnd = indEnd < bufSz-1 ? indEnd+1 : 0;
      if (indEnd == indBeg) {
        indBeg = indBeg < bufSz-1 ? indBeg+1 : 0;
      }

      Wind_songweiwei_DW.SFunction_IWORK.indBeg = indBeg;
      Wind_songweiwei_DW.SFunction_IWORK.indEnd = indEnd;
    }

    /* Update for UnitDelay: '<S158>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_f = rtb_Switch[0];

    /* Update for UnitDelay: '<S154>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_m = rtb_Switch_e;

    /* Update for DiscreteIntegrator: '<S153>/Integ4' */
    Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_f = 0U;
    Wind_songweiwei_DW.Integ4_DSTATE_n = Wind_songweiwei_P.Integ4_gainval_e *
      rtb_Switch[1] + Wind_songweiwei_B.Integ4_n;

    /* S-Function block: <S157>/S-Function  */
    {
      int_T indBeg = Wind_songweiwei_DW.SFunction_IWORK_f.indBeg;
      int_T indEnd = Wind_songweiwei_DW.SFunction_IWORK_f.indEnd;
      int_T bufSz = Wind_songweiwei_DW.SFunction_IWORK_f.bufSz;
      indEnd = indEnd < bufSz-1 ? indEnd+1 : 0;
      if (indEnd == indBeg) {
        indBeg = indBeg < bufSz-1 ? indBeg+1 : 0;
      }

      Wind_songweiwei_DW.SFunction_IWORK_f.indBeg = indBeg;
      Wind_songweiwei_DW.SFunction_IWORK_f.indEnd = indEnd;
    }

    /* Update for UnitDelay: '<S156>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_fd = rtb_Switch[1];

    /* Update for UnitDelay: '<S153>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_l = rtb_Switch_d;

    /* Update for UnitDelay: '<S131>/IC = 0' */
    Wind_songweiwei_DW.IC0_DSTATE = rtb_Sum;

    /* Update for UnitDelay: '<S193>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_b[0] = rtb_x1k1_n[0];
    Wind_songweiwei_DW.Delay_x1_DSTATE_b[1] = rtb_x1k1_n[1];
    Wind_songweiwei_DW.Delay_x1_DSTATE_b[2] = rtb_x1k1_n[2];

    /* Update for UnitDelay: '<S193>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_pn[0] = rtb_x2k1_n[0];
    Wind_songweiwei_DW.Delay_x2_DSTATE_pn[1] = rtb_x2k1_n[1];
    Wind_songweiwei_DW.Delay_x2_DSTATE_pn[2] = rtb_x2k1_n[2];

    /* Update for DiscreteIntegrator: '<S129>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[0] +=
      Wind_songweiwei_P.DiscreteTimeIntegrator_gainv_lq * rtb_Kp5_a[0];
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[1] +=
      Wind_songweiwei_P.DiscreteTimeIntegrator_gainv_lq * rtb_Kp5_a[1];
    if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[0] >=
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_d) {
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[0] =
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_d;
    } else {
      if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[0] <=
          Wind_songweiwei_P.DiscreteTimeIntegrator_LowerS_k) {
        Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[0] =
          Wind_songweiwei_P.DiscreteTimeIntegrator_LowerS_k;
      }
    }

    if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[1] >=
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_d) {
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[1] =
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_d;
    } else {
      if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[1] <=
          Wind_songweiwei_P.DiscreteTimeIntegrator_LowerS_k) {
        Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[1] =
          Wind_songweiwei_P.DiscreteTimeIntegrator_LowerS_k;
      }
    }

    /* End of Update for DiscreteIntegrator: '<S129>/Discrete-Time Integrator' */

    /* Update for UnitDelay: '<S209>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_f = rtb_x1k1_o;

    /* Update for UnitDelay: '<S209>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_d = rtb_x2k1_ac;

    /* Update for DiscreteIntegrator: '<S168>/Integ4' */
    Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_l = 0U;
    Wind_songweiwei_DW.Integ4_DSTATE_nw = Wind_songweiwei_P.Integ4_gainval_l *
      rtb_Switch_o[0] + Wind_songweiwei_B.Integ4_j;

    /* S-Function block: <S173>/S-Function  */
    {
      int_T indBeg = Wind_songweiwei_DW.SFunction_IWORK_c.indBeg;
      int_T indEnd = Wind_songweiwei_DW.SFunction_IWORK_c.indEnd;
      int_T bufSz = Wind_songweiwei_DW.SFunction_IWORK_c.bufSz;
      indEnd = indEnd < bufSz-1 ? indEnd+1 : 0;
      if (indEnd == indBeg) {
        indBeg = indBeg < bufSz-1 ? indBeg+1 : 0;
      }

      Wind_songweiwei_DW.SFunction_IWORK_c.indBeg = indBeg;
      Wind_songweiwei_DW.SFunction_IWORK_c.indEnd = indEnd;
    }

    /* Update for UnitDelay: '<S172>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_fz = rtb_Switch_o[0];

    /* Update for UnitDelay: '<S168>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_b = rtb_Switch_n;

    /* Update for DiscreteIntegrator: '<S167>/Integ4' */
    Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_h = 0U;
    Wind_songweiwei_DW.Integ4_DSTATE_l = Wind_songweiwei_P.Integ4_gainval_b *
      rtb_Switch_o[1] + Wind_songweiwei_B.Integ4_m;

    /* S-Function block: <S171>/S-Function  */
    {
      int_T indBeg = Wind_songweiwei_DW.SFunction_IWORK_k.indBeg;
      int_T indEnd = Wind_songweiwei_DW.SFunction_IWORK_k.indEnd;
      int_T bufSz = Wind_songweiwei_DW.SFunction_IWORK_k.bufSz;
      indEnd = indEnd < bufSz-1 ? indEnd+1 : 0;
      if (indEnd == indBeg) {
        indBeg = indBeg < bufSz-1 ? indBeg+1 : 0;
      }

      Wind_songweiwei_DW.SFunction_IWORK_k.indBeg = indBeg;
      Wind_songweiwei_DW.SFunction_IWORK_k.indEnd = indEnd;
    }

    /* Update for UnitDelay: '<S170>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_e = rtb_Switch_o[1];

    /* Update for UnitDelay: '<S167>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_d = rtb_Switch_h;

    /* Update for UnitDelay: '<S205>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_in = rtb_x1k1_ag;

    /* Update for UnitDelay: '<S205>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_a = rtb_x2k1_il;

    /* Update for DiscreteIntegrator: '<S251>/Integ4' */
    Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_a = 0U;
    Wind_songweiwei_DW.Integ4_DSTATE_a = Wind_songweiwei_P.Integ4_gainval_d *
      rtb_Switch_di[1] + Wind_songweiwei_B.Integ4_o;

    /* Level2 S-Function Block: '<S253>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Wind_songweiwei_M->childSfunctions[2];
      sfcnUpdate(rts, 1);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for UnitDelay: '<S252>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_n = rtb_Switch_di[1];

    /* Update for UnitDelay: '<S251>/Unit Delay1' */
    Wind_songweiwei_DW.UnitDelay1_DSTATE_p = rtb_Switch_g;

    /* Update for DiscreteTransferFcn: '<S227>/Discrete Derivative ' */
    Wind_songweiwei_DW.DiscreteDerivative_states =
      Wind_songweiwei_DW.DiscreteDerivative_tmp;

    /* Update for DiscreteIntegrator: '<S227>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_k +=
      Wind_songweiwei_P.DiscreteTimeIntegrator_gainva_o * rtb_Divide;
    if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_k >=
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_e) {
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_k =
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_e;
    } else {
      if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_k <=
          Wind_songweiwei_P.DiscreteTimeIntegrator_LowerS_b) {
        Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_k =
          Wind_songweiwei_P.DiscreteTimeIntegrator_LowerS_b;
      }
    }

    /* End of Update for DiscreteIntegrator: '<S227>/Discrete-Time Integrator' */

    /* Update for UnitDelay: '<S247>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_bs = rtb_x1k1_oi;

    /* Update for UnitDelay: '<S247>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_oc = rtb_x2k1_hi;

    /* Update for UnitDelay: '<S67>/Unit Delay2' */
    Wind_songweiwei_DW.UnitDelay2_DSTATE = rtb_Saturation2_e;

    /* Update for UnitDelay: '<S310>/Delay_x' */
    Wind_songweiwei_DW.Delay_x_DSTATE_o = rtb_xk1_go;

    /* Update for DiscreteIntegrator: '<S305>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_o +=
      Wind_songweiwei_P.DiscreteTimeIntegrator_gainva_f * rtb_Kp5_l;
    if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_o >=
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_f) {
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_o =
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_f;
    } else {
      if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_o <=
          Wind_songweiwei_P.DiscreteTimeIntegrator_Lower_pa) {
        Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_o =
          Wind_songweiwei_P.DiscreteTimeIntegrator_Lower_pa;
      }
    }

    /* End of Update for DiscreteIntegrator: '<S305>/Discrete-Time Integrator' */

    /* Update for RateLimiter: '<S125>/Rate Limiter   1' */
    Wind_songweiwei_DW.PrevY_m = Wind_songweiwei_B.RateLimiter1;
    Wind_songweiwei_DW.LastMajorTime = Wind_songweiwei_M->Timing.t[0];

    /* Update for DiscreteIntegrator: '<S304>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_kk +=
      Wind_songweiwei_P.DiscreteTimeIntegrator_gainva_e * rtb_Kp5_i;
    if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_kk >=
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_h) {
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_kk =
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_h;
    } else {
      if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_kk <=
          Wind_songweiwei_P.DiscreteTimeIntegrator_LowerS_l) {
        Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_kk =
          Wind_songweiwei_P.DiscreteTimeIntegrator_LowerS_l;
      }
    }

    /* End of Update for DiscreteIntegrator: '<S304>/Discrete-Time Integrator' */

    /* Update for DiscreteIntegrator: '<S306>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_b +=
      Wind_songweiwei_P.DiscreteTimeIntegrator_gainv_cb * rtb_Kp5_aq;
    if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_b >=
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_i) {
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_b =
        Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_i;
    } else {
      if (Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_b <=
          Wind_songweiwei_P.DiscreteTimeIntegrator_Lower_lu) {
        Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_b =
          Wind_songweiwei_P.DiscreteTimeIntegrator_Lower_lu;
      }
    }

    /* End of Update for DiscreteIntegrator: '<S306>/Discrete-Time Integrator' */

    /* Update for UnitDelay: '<S309>/Delay_x' */
    Wind_songweiwei_DW.Delay_x_DSTATE_i = rtb_xk1_a;

    /* Update for UnitDelay: '<S67>/Unit Delay1' */
    Wind_songweiwei_DW.UnitDelay1_DSTATE_n = rtb_yk_a3;

    /* Update for UnitDelay: '<S67>/Unit Delay4' */
    Wind_songweiwei_DW.UnitDelay4_DSTATE = rtb_Saturation2;

    /* Update for DiscreteIntegrator: '<S69>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_bh +=
      Wind_songweiwei_P.DiscreteTimeIntegrator_gainva_m * rtb__2H_WT;

    /* Update for DiscreteIntegrator: '<S69>/Discrete-Time Integrator1' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTAT_l +=
      Wind_songweiwei_P.DiscreteTimeIntegrator1_gainv_i * rtb_wbase;

    /* Update for UnitDelay: '<S342>/voltages' incorporates:
     *  Constant: '<S327>/[ Vkd =0 Vkq1=0  Vkq2=0 ]'
     */
    Wind_songweiwei_DW.voltages_DSTATE[0] = rtb_Fcn2;
    Wind_songweiwei_DW.voltages_DSTATE[1] = rtb_Fcn3;
    Wind_songweiwei_DW.voltages_DSTATE[2] = rtb_Vfd;
    Wind_songweiwei_DW.voltages_DSTATE[3] =
      Wind_songweiwei_P.Vkd0Vkq10Vkq20_Value[0];
    Wind_songweiwei_DW.voltages_DSTATE[4] =
      Wind_songweiwei_P.Vkd0Vkq10Vkq20_Value[1];

    /* Update for DiscreteIntegrator: '<S329>/Rotor speed deviation (dw)' */
    Wind_songweiwei_DW.Rotorspeeddeviationdw_SYSTEM_EN = 0U;
    Wind_songweiwei_DW.Rotorspeeddeviationdw_DSTATE =
      Wind_songweiwei_P.Rotorspeeddeviationdw_gainval * rtb_uH +
      Wind_songweiwei_B.dw;
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++Wind_songweiwei_M->Timing.clockTick0)) {
    ++Wind_songweiwei_M->Timing.clockTickH0;
  }

  Wind_songweiwei_M->Timing.t[0] = Wind_songweiwei_M->Timing.clockTick0 *
    Wind_songweiwei_M->Timing.stepSize0 + Wind_songweiwei_M->Timing.clockTickH0 *
    Wind_songweiwei_M->Timing.stepSize0 * 4294967296.0;

  {
    /* Update absolute timer for sample time: [2.0E-6s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    if (!(++Wind_songweiwei_M->Timing.clockTick1)) {
      ++Wind_songweiwei_M->Timing.clockTickH1;
    }

    Wind_songweiwei_M->Timing.t[1] = Wind_songweiwei_M->Timing.clockTick1 *
      Wind_songweiwei_M->Timing.stepSize1 +
      Wind_songweiwei_M->Timing.clockTickH1 *
      Wind_songweiwei_M->Timing.stepSize1 * 4294967296.0;
  }
}

/* Model initialize function */
void Wind_songweiwei_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  Wind_songweiwei_P.avoiddivisionby0_UpperSat = rtInf;
  Wind_songweiwei_P.Saturation_UpperSat = rtInf;
  Wind_songweiwei_P.Saturation_UpperSat_c = rtInf;
  Wind_songweiwei_P.avoiddivisionby0_UpperSat_d = rtInf;
  Wind_songweiwei_P.DiscreteTimeIntegrator_UpperS_e = rtInf;
  Wind_songweiwei_P.Saturation1_UpperSat_h = rtInf;
  Wind_songweiwei_P.inf_UpperSat = rtInf;
  Wind_songweiwei_P.Avoiddivbyzero_UpperSat = rtInf;
  Wind_songweiwei_P.Avoiddivbyzero_UpperSat_l = rtInf;

  /* initialize real-time model */
  (void) memset((void *)Wind_songweiwei_M, 0,
                sizeof(RT_MODEL_Wind_songweiwei_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&Wind_songweiwei_M->solverInfo,
                          &Wind_songweiwei_M->Timing.simTimeStep);
    rtsiSetTPtr(&Wind_songweiwei_M->solverInfo, &rtmGetTPtr(Wind_songweiwei_M));
    rtsiSetStepSizePtr(&Wind_songweiwei_M->solverInfo,
                       &Wind_songweiwei_M->Timing.stepSize0);
    rtsiSetErrorStatusPtr(&Wind_songweiwei_M->solverInfo, ((const char_T **)
      (&rtmGetErrorStatus(Wind_songweiwei_M))));
    rtsiSetRTModelPtr(&Wind_songweiwei_M->solverInfo, Wind_songweiwei_M);
  }

  rtsiSetSimTimeStep(&Wind_songweiwei_M->solverInfo, MAJOR_TIME_STEP);
  rtsiSetSolverName(&Wind_songweiwei_M->solverInfo,"FixedStepDiscrete");
  Wind_songweiwei_M->solverInfoPtr = (&Wind_songweiwei_M->solverInfo);

  /* Initialize timing info */
  {
    int_T *mdlTsMap = Wind_songweiwei_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    Wind_songweiwei_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    Wind_songweiwei_M->Timing.sampleTimes =
      (&Wind_songweiwei_M->Timing.sampleTimesArray[0]);
    Wind_songweiwei_M->Timing.offsetTimes =
      (&Wind_songweiwei_M->Timing.offsetTimesArray[0]);

    /* task periods */
    Wind_songweiwei_M->Timing.sampleTimes[0] = (0.0);
    Wind_songweiwei_M->Timing.sampleTimes[1] = (2.0E-6);

    /* task offsets */
    Wind_songweiwei_M->Timing.offsetTimes[0] = (0.0);
    Wind_songweiwei_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(Wind_songweiwei_M, &Wind_songweiwei_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = Wind_songweiwei_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    Wind_songweiwei_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(Wind_songweiwei_M, -1);
  Wind_songweiwei_M->Timing.stepSize0 = 2.0E-6;
  Wind_songweiwei_M->Timing.stepSize1 = 2.0E-6;
  Wind_songweiwei_M->solverInfoPtr = (&Wind_songweiwei_M->solverInfo);
  Wind_songweiwei_M->Timing.stepSize = (2.0E-6);
  rtsiSetFixedStepSize(&Wind_songweiwei_M->solverInfo, 2.0E-6);
  rtsiSetSolverMode(&Wind_songweiwei_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  (void) memset(((void *) &Wind_songweiwei_B), 0,
                sizeof(B_Wind_songweiwei_T));

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      Wind_songweiwei_B.DataTypeConversion[i] = 0.0;
    }

    for (i = 0; i < 50; i++) {
      Wind_songweiwei_B.StateSpace_o1[i] = 0.0;
    }

    for (i = 0; i < 19; i++) {
      Wind_songweiwei_B.StateSpace_o2[i] = 0.0;
    }

    for (i = 0; i < 6; i++) {
      Wind_songweiwei_B.Product[i] = 0.0;
    }

    for (i = 0; i < 25; i++) {
      Wind_songweiwei_B.Linv[i] = 0.0;
    }

    for (i = 0; i < 25; i++) {
      Wind_songweiwei_B.RLinv[i] = 0.0;
    }

    Wind_songweiwei_B.ib[0] = 0.0;
    Wind_songweiwei_B.ib[1] = 0.0;
    Wind_songweiwei_B.ib[2] = 0.0;
    Wind_songweiwei_B.Sum5[0] = 0.0;
    Wind_songweiwei_B.Sum5[1] = 0.0;
    Wind_songweiwei_B.Sum5[2] = 0.0;
    Wind_songweiwei_B.Integ4 = 0.0;
    Wind_songweiwei_B.Delay = 0.0;
    Wind_songweiwei_B.SFunction = 0.0;
    Wind_songweiwei_B.Integ4_n = 0.0;
    Wind_songweiwei_B.Delay_f = 0.0;
    Wind_songweiwei_B.SFunction_j = 0.0;
    Wind_songweiwei_B.Integ4_j = 0.0;
    Wind_songweiwei_B.Delay_n = 0.0;
    Wind_songweiwei_B.SFunction_f = 0.0;
    Wind_songweiwei_B.Integ4_m = 0.0;
    Wind_songweiwei_B.Delay_h = 0.0;
    Wind_songweiwei_B.SFunction_m = 0.0;
    Wind_songweiwei_B.Integ4_o = 0.0;
    Wind_songweiwei_B.Delay_o = 0.0;
    Wind_songweiwei_B.SFunction_k = 0.0;
    Wind_songweiwei_B.Switch = 0.0;
    Wind_songweiwei_B.RateLimiter1 = 0.0;
    Wind_songweiwei_B.dw = 0.0;
    Wind_songweiwei_B.Lmsatq = 0.0;
    Wind_songweiwei_B.Integ4_e = 0.0;
    Wind_songweiwei_B.Delay_e = 0.0;
    Wind_songweiwei_B.SFunction_h = 0.0;
    Wind_songweiwei_B.Integ4_h = 0.0;
    Wind_songweiwei_B.Delay_l = 0.0;
    Wind_songweiwei_B.SFunction_i = 0.0;
    Wind_songweiwei_B.MathFunction = 0.0;
    Wind_songweiwei_B.Fcn = 0.0;
    Wind_songweiwei_B.Fcn1 = 0.0;
    Wind_songweiwei_B.Fcn_c = 0.0;
    Wind_songweiwei_B.Fcn1_b = 0.0;
    Wind_songweiwei_B.LookUpTable = 0.0;
    Wind_songweiwei_B.DiscreteTimeIntegrator = 0.0;
    Wind_songweiwei_B.Switch2 = 0.0;
    Wind_songweiwei_B.Switch3 = 0.0;
    Wind_songweiwei_B.Product1[0] = 0.0;
    Wind_songweiwei_B.Product1[1] = 0.0;
    Wind_songweiwei_B.Product1[2] = 0.0;
    Wind_songweiwei_B.Product1_k[0] = 0.0;
    Wind_songweiwei_B.Product1_k[1] = 0.0;
    Wind_songweiwei_B.Product1_k[2] = 0.0;
    Wind_songweiwei_B.Subsystem1.Fcn = 0.0;
    Wind_songweiwei_B.Subsystem1.Fcn1 = 0.0;
    Wind_songweiwei_B.Subsystempi2delay.Fcn = 0.0;
    Wind_songweiwei_B.Subsystempi2delay.Fcn1 = 0.0;
    Wind_songweiwei_B.Subsystem1_p.Fcn = 0.0;
    Wind_songweiwei_B.Subsystem1_p.Fcn1 = 0.0;
    Wind_songweiwei_B.Subsystempi2delay_jr.Fcn = 0.0;
    Wind_songweiwei_B.Subsystempi2delay_jr.Fcn1 = 0.0;
    Wind_songweiwei_B.Subsystem1_l.Fcn = 0.0;
    Wind_songweiwei_B.Subsystem1_l.Fcn1 = 0.0;
    Wind_songweiwei_B.Subsystempi2delay_j.Fcn = 0.0;
    Wind_songweiwei_B.Subsystempi2delay_j.Fcn1 = 0.0;
    Wind_songweiwei_B.Subsystem1_l1.Fcn = 0.0;
    Wind_songweiwei_B.Subsystem1_l1.Fcn1 = 0.0;
    Wind_songweiwei_B.Subsystempi2delay_p.Fcn = 0.0;
    Wind_songweiwei_B.Subsystempi2delay_p.Fcn1 = 0.0;
    Wind_songweiwei_B.Subsystem1_lp.Fcn = 0.0;
    Wind_songweiwei_B.Subsystem1_lp.Fcn1 = 0.0;
    Wind_songweiwei_B.Subsystempi2delay_dg.Fcn = 0.0;
    Wind_songweiwei_B.Subsystempi2delay_dg.Fcn1 = 0.0;
    Wind_songweiwei_B.Subsystem1_oh.Fcn = 0.0;
    Wind_songweiwei_B.Subsystem1_oh.Fcn1 = 0.0;
    Wind_songweiwei_B.Subsystempi2delay_h.Fcn = 0.0;
    Wind_songweiwei_B.Subsystempi2delay_h.Fcn1 = 0.0;
    Wind_songweiwei_B.Subsystem1_o.Fcn = 0.0;
    Wind_songweiwei_B.Subsystem1_o.Fcn1 = 0.0;
    Wind_songweiwei_B.Subsystempi2delay_k.Fcn = 0.0;
    Wind_songweiwei_B.Subsystempi2delay_k.Fcn1 = 0.0;
  }

  /* states (dwork) */
  (void) memset((void *)&Wind_songweiwei_DW, 0,
                sizeof(DW_Wind_songweiwei_T));

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      Wind_songweiwei_DW.UnitDelay_DSTATE[i] = 0.0;
    }
  }

  Wind_songweiwei_DW.UnitDelay6_DSTATE[0] = 0.0;
  Wind_songweiwei_DW.UnitDelay6_DSTATE[1] = 0.0;
  Wind_songweiwei_DW.UnitDelay6_DSTATE[2] = 0.0;
  Wind_songweiwei_DW.Rotorangledthetae_DSTATE = 0.0;

  {
    int_T i;
    for (i = 0; i < 5; i++) {
      Wind_songweiwei_DW.fluxes_DSTATE[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 45; i++) {
      Wind_songweiwei_DW.StateSpace_DSTATE[i] = 0.0;
    }
  }

  Wind_songweiwei_DW.UnitDelay3_DSTATE = 0.0;
  Wind_songweiwei_DW.UnitDelay7_DSTATE = 0.0;
  Wind_songweiwei_DW.dw_delay_DSTATE = 0.0;
  Wind_songweiwei_DW.dw_predict_DSTATE = 0.0;
  Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_c[0] = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_c[1] = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_c[2] = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_e[0] = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_e[1] = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_e[2] = 0.0;
  Wind_songweiwei_DW.theta_DSTATE = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_a = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_p = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_i[0] = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_i[1] = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_i[2] = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_m[0] = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_m[1] = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_m[2] = 0.0;
  Wind_songweiwei_DW.Delay_x_DSTATE = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_ck = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_ee = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_m[0] = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_m[1] = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_m[2] = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_o[0] = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_o[1] = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_o[2] = 0.0;
  Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_h = 0.0;
  Wind_songweiwei_DW.UnitDelay_DSTATE_p = 0.0;
  Wind_songweiwei_DW.ICi_ic_DSTATE = 0.0;
  Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_g = 0.0;
  Wind_songweiwei_DW.Delay_x_DSTATE_m = 0.0;
  Wind_songweiwei_DW.UnitDelay1_DSTATE = 0.0;
  Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTATE = 0.0;
  Wind_songweiwei_DW.Integ4_DSTATE = 0.0;
  Wind_songweiwei_DW.UnitDelay_DSTATE_f = 0.0;
  Wind_songweiwei_DW.UnitDelay_DSTATE_m = 0.0;
  Wind_songweiwei_DW.Integ4_DSTATE_n = 0.0;
  Wind_songweiwei_DW.UnitDelay_DSTATE_fd = 0.0;
  Wind_songweiwei_DW.UnitDelay_DSTATE_l = 0.0;
  Wind_songweiwei_DW.IC0_DSTATE = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_b[0] = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_b[1] = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_b[2] = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_pn[0] = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_pn[1] = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_pn[2] = 0.0;
  Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[0] = 0.0;
  Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[1] = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_f = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_d = 0.0;
  Wind_songweiwei_DW.Integ4_DSTATE_nw = 0.0;
  Wind_songweiwei_DW.UnitDelay_DSTATE_fz = 0.0;
  Wind_songweiwei_DW.UnitDelay_DSTATE_b = 0.0;
  Wind_songweiwei_DW.Integ4_DSTATE_l = 0.0;
  Wind_songweiwei_DW.UnitDelay_DSTATE_e = 0.0;
  Wind_songweiwei_DW.UnitDelay_DSTATE_d = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_in = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_a = 0.0;
  Wind_songweiwei_DW.Integ4_DSTATE_a = 0.0;
  Wind_songweiwei_DW.UnitDelay_DSTATE_n = 0.0;
  Wind_songweiwei_DW.UnitDelay1_DSTATE_p = 0.0;
  Wind_songweiwei_DW.DiscreteDerivative_states = 0.0;
  Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_k = 0.0;
  Wind_songweiwei_DW.Delay_x1_DSTATE_bs = 0.0;
  Wind_songweiwei_DW.Delay_x2_DSTATE_oc = 0.0;
  Wind_songweiwei_DW.UnitDelay2_DSTATE = 0.0;
  Wind_songweiwei_DW.Delay_x_DSTATE_o = 0.0;
  Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_o = 0.0;
  Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_kk = 0.0;
  Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_b = 0.0;
  Wind_songweiwei_DW.Delay_x_DSTATE_i = 0.0;
  Wind_songweiwei_DW.UnitDelay1_DSTATE_n = 0.0;
  Wind_songweiwei_DW.UnitDelay4_DSTATE = 0.0;
  Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_bh = 0.0;
  Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTAT_l = 0.0;

  {
    int_T i;
    for (i = 0; i < 5; i++) {
      Wind_songweiwei_DW.voltages_DSTATE[i] = 0.0;
    }
  }

  Wind_songweiwei_DW.Rotorspeeddeviationdw_DSTATE = 0.0;

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_b4[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 6; i++) {
      Wind_songweiwei_DW.UnitDelay_DSTATE_fv[i] = 0.0;
    }
  }

  Wind_songweiwei_DW.Lmd_sat_DSTATE = 0.0;
  Wind_songweiwei_DW.Lmq_sat_DSTATE = 0.0;
  Wind_songweiwei_DW.Integ4_DSTATE_nn = 0.0;
  Wind_songweiwei_DW.UnitDelay_DSTATE_i = 0.0;
  Wind_songweiwei_DW.UnitDelay1_DSTATE_f = 0.0;
  Wind_songweiwei_DW.Integ4_DSTATE_c = 0.0;
  Wind_songweiwei_DW.UnitDelay_DSTATE_g = 0.0;
  Wind_songweiwei_DW.UnitDelay1_DSTATE_c = 0.0;
  Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_l = 0.0;
  Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTAT_a = 0.0;
  Wind_songweiwei_DW.UnitDelay_DSTATE_gd = 0.0;
  Wind_songweiwei_DW.DiscreteDerivative_tmp = 0.0;
  Wind_songweiwei_DW.PrevY = 0.0;
  Wind_songweiwei_DW.PrevY_m = 0.0;
  Wind_songweiwei_DW.LastMajorTime = 0.0;

  {
    int_T i;
    for (i = 0; i < 25; i++) {
      Wind_songweiwei_DW.inversion_DWORK4[i] = 0.0;
    }
  }

  {
    int_T i;
    for (i = 0; i < 25; i++) {
      Wind_songweiwei_DW.inversion_DWORK4_h[i] = 0.0;
    }
  }

  Wind_songweiwei_DW.SFunction_RWORK = 0.0;
  Wind_songweiwei_DW.SFunction_RWORK_e = 0.0;
  Wind_songweiwei_DW.SFunction_RWORK_g = 0.0;

  /* external outputs */
  Wind_songweiwei_Y.Out_Vabc_B575[0] = 0.0;
  Wind_songweiwei_Y.Out_Vabc_B575[1] = 0.0;
  Wind_songweiwei_Y.Out_Vabc_B575[2] = 0.0;
  Wind_songweiwei_Y.Out_labc_B575[0] = 0.0;
  Wind_songweiwei_Y.Out_labc_B575[1] = 0.0;
  Wind_songweiwei_Y.Out_labc_B575[2] = 0.0;
  Wind_songweiwei_Y.Out_P = 0.0;
  Wind_songweiwei_Y.Out_Q = 0.0;
  Wind_songweiwei_Y.Out_Vdc = 0.0;
  Wind_songweiwei_Y.Out_wr = 0.0;
  Wind_songweiwei_Y.Out_Vabc_B25[0] = 0.0;
  Wind_songweiwei_Y.Out_Vabc_B25[1] = 0.0;
  Wind_songweiwei_Y.Out_Vabc_B25[2] = 0.0;
  Wind_songweiwei_Y.Out_labc_B25[0] = 0.0;
  Wind_songweiwei_Y.Out_labc_B25[1] = 0.0;
  Wind_songweiwei_Y.Out_labc_B25[2] = 0.0;

  /* child S-Function registration */
  {
    RTWSfcnInfo *sfcnInfo = &Wind_songweiwei_M->NonInlinedSFcns.sfcnInfo;
    Wind_songweiwei_M->sfcnInfo = (sfcnInfo);
    rtssSetErrorStatusPtr(sfcnInfo, ((const char_T **)(&rtmGetErrorStatus
      (Wind_songweiwei_M))));
    rtssSetNumRootSampTimesPtr(sfcnInfo, &Wind_songweiwei_M->Sizes.numSampTimes);
    Wind_songweiwei_M->NonInlinedSFcns.taskTimePtrs[0] = &(rtmGetTPtr
      (Wind_songweiwei_M)[0]);
    Wind_songweiwei_M->NonInlinedSFcns.taskTimePtrs[1] = &(rtmGetTPtr
      (Wind_songweiwei_M)[1]);
    rtssSetTPtrPtr(sfcnInfo,Wind_songweiwei_M->NonInlinedSFcns.taskTimePtrs);
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(Wind_songweiwei_M));
    rtssSetTFinalPtr(sfcnInfo, &rtmGetTFinal(Wind_songweiwei_M));
    rtssSetTimeOfLastOutputPtr(sfcnInfo, &rtmGetTimeOfLastOutput
      (Wind_songweiwei_M));
    rtssSetStepSizePtr(sfcnInfo, &Wind_songweiwei_M->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo, &rtmGetStopRequested(Wind_songweiwei_M));
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo,
      &Wind_songweiwei_M->ModelData.derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo,
      &Wind_songweiwei_M->ModelData.zCCacheNeedsReset);
    rtssSetBlkStateChangePtr(sfcnInfo,
      &Wind_songweiwei_M->ModelData.blkStateChange);
    rtssSetSampleHitsPtr(sfcnInfo, &Wind_songweiwei_M->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
      &Wind_songweiwei_M->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &Wind_songweiwei_M->simMode);
    rtssSetSolverInfoPtr(sfcnInfo, &Wind_songweiwei_M->solverInfoPtr);
  }

  Wind_songweiwei_M->Sizes.numSFcns = (3);

  /* register each child */
  {
    (void) memset((void *)&Wind_songweiwei_M->NonInlinedSFcns.childSFunctions[0],
                  0,
                  3*sizeof(SimStruct));
    Wind_songweiwei_M->childSfunctions =
      (&Wind_songweiwei_M->NonInlinedSFcns.childSFunctionPtrs[0]);
    Wind_songweiwei_M->childSfunctions[0] =
      (&Wind_songweiwei_M->NonInlinedSFcns.childSFunctions[0]);
    Wind_songweiwei_M->childSfunctions[1] =
      (&Wind_songweiwei_M->NonInlinedSFcns.childSFunctions[1]);
    Wind_songweiwei_M->childSfunctions[2] =
      (&Wind_songweiwei_M->NonInlinedSFcns.childSFunctions[2]);

    /* Level2 S-Function Block: Wind_songweiwei/<S237>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Wind_songweiwei_M->childSfunctions[0];

      /* timing info */
      time_T *sfcnPeriod = Wind_songweiwei_M->NonInlinedSFcns.Sfcn0.sfcnPeriod;
      time_T *sfcnOffset = Wind_songweiwei_M->NonInlinedSFcns.Sfcn0.sfcnOffset;
      int_T *sfcnTsMap = Wind_songweiwei_M->NonInlinedSFcns.Sfcn0.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts, &Wind_songweiwei_M->NonInlinedSFcns.blkInfo2[0]);
      }

      ssSetRTWSfcnInfo(rts, Wind_songweiwei_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &Wind_songweiwei_M->NonInlinedSFcns.methods2[0]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &Wind_songweiwei_M->NonInlinedSFcns.methods3[0]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &Wind_songweiwei_M->NonInlinedSFcns.statesInfo2[0]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn0.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Wind_songweiwei_M->NonInlinedSFcns.Sfcn0.UPtrs0;
          sfcnUPtrs[0] = &Wind_songweiwei_B.Integ4_e;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Wind_songweiwei_M->NonInlinedSFcns.Sfcn0.UPtrs1;
          sfcnUPtrs[0] = &Wind_songweiwei_B.Delay_e;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn0.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &Wind_songweiwei_B.SFunction_h));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/Mean (Variable Frequency)1/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,Wind_songweiwei_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn0.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)Wind_songweiwei_P.SFunction_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)Wind_songweiwei_P.SFunction_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)Wind_songweiwei_P.SFunction_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)Wind_songweiwei_P.SFunction_P4_Size);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *) &Wind_songweiwei_DW.SFunction_RWORK_e);
      ssSetIWork(rts, (int_T *) &Wind_songweiwei_DW.SFunction_IWORK_h);
      ssSetPWork(rts, (void **) &Wind_songweiwei_DW.SFunction_PWORK_ge);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn0.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn0.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &Wind_songweiwei_DW.SFunction_RWORK_e);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &Wind_songweiwei_DW.SFunction_IWORK_h);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &Wind_songweiwei_DW.SFunction_PWORK_ge);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-6);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 1);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetInputPortWidth(rts, 1, 1);
      ssSetInputPortDataType(rts, 1, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 1, 0);
      ssSetInputPortFrameData(rts, 1, 0);
      ssSetOutputPortWidth(rts, 0, 1);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumIWork(rts, 1);
      ssSetNumPWork(rts, 1);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: Wind_songweiwei/<S240>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Wind_songweiwei_M->childSfunctions[1];

      /* timing info */
      time_T *sfcnPeriod = Wind_songweiwei_M->NonInlinedSFcns.Sfcn1.sfcnPeriod;
      time_T *sfcnOffset = Wind_songweiwei_M->NonInlinedSFcns.Sfcn1.sfcnOffset;
      int_T *sfcnTsMap = Wind_songweiwei_M->NonInlinedSFcns.Sfcn1.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts, &Wind_songweiwei_M->NonInlinedSFcns.blkInfo2[1]);
      }

      ssSetRTWSfcnInfo(rts, Wind_songweiwei_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &Wind_songweiwei_M->NonInlinedSFcns.methods2[1]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &Wind_songweiwei_M->NonInlinedSFcns.methods3[1]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &Wind_songweiwei_M->NonInlinedSFcns.statesInfo2[1]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn1.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Wind_songweiwei_M->NonInlinedSFcns.Sfcn1.UPtrs0;
          sfcnUPtrs[0] = &Wind_songweiwei_B.Integ4_h;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Wind_songweiwei_M->NonInlinedSFcns.Sfcn1.UPtrs1;
          sfcnUPtrs[0] = &Wind_songweiwei_B.Delay_l;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn1.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &Wind_songweiwei_B.SFunction_i));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/Mean (Variable Frequency)2/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,Wind_songweiwei_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn1.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)Wind_songweiwei_P.SFunction_P1_Size_j);
        ssSetSFcnParam(rts, 1, (mxArray*)Wind_songweiwei_P.SFunction_P2_Size_j);
        ssSetSFcnParam(rts, 2, (mxArray*)Wind_songweiwei_P.SFunction_P3_Size_e);
        ssSetSFcnParam(rts, 3, (mxArray*)Wind_songweiwei_P.SFunction_P4_Size_c);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *) &Wind_songweiwei_DW.SFunction_RWORK_g);
      ssSetIWork(rts, (int_T *) &Wind_songweiwei_DW.SFunction_IWORK_n);
      ssSetPWork(rts, (void **) &Wind_songweiwei_DW.SFunction_PWORK_cs);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn1.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn1.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &Wind_songweiwei_DW.SFunction_RWORK_g);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &Wind_songweiwei_DW.SFunction_IWORK_n);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &Wind_songweiwei_DW.SFunction_PWORK_cs);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-6);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 1);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetInputPortWidth(rts, 1, 1);
      ssSetInputPortDataType(rts, 1, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 1, 0);
      ssSetInputPortFrameData(rts, 1, 0);
      ssSetOutputPortWidth(rts, 0, 1);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumIWork(rts, 1);
      ssSetNumPWork(rts, 1);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }

    /* Level2 S-Function Block: Wind_songweiwei/<S253>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Wind_songweiwei_M->childSfunctions[2];

      /* timing info */
      time_T *sfcnPeriod = Wind_songweiwei_M->NonInlinedSFcns.Sfcn2.sfcnPeriod;
      time_T *sfcnOffset = Wind_songweiwei_M->NonInlinedSFcns.Sfcn2.sfcnOffset;
      int_T *sfcnTsMap = Wind_songweiwei_M->NonInlinedSFcns.Sfcn2.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts, &Wind_songweiwei_M->NonInlinedSFcns.blkInfo2[2]);
      }

      ssSetRTWSfcnInfo(rts, Wind_songweiwei_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &Wind_songweiwei_M->NonInlinedSFcns.methods2[2]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &Wind_songweiwei_M->NonInlinedSFcns.methods3[2]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &Wind_songweiwei_M->NonInlinedSFcns.statesInfo2[2]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn2.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Wind_songweiwei_M->NonInlinedSFcns.Sfcn2.UPtrs0;
          sfcnUPtrs[0] = &Wind_songweiwei_B.Integ4_o;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Wind_songweiwei_M->NonInlinedSFcns.Sfcn2.UPtrs1;
          sfcnUPtrs[0] = &Wind_songweiwei_B.Delay_o;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn2.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *)
            &Wind_songweiwei_B.SFunction_k));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Variable Frequency Mean value/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,Wind_songweiwei_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn2.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)Wind_songweiwei_P.SFunction_P1_Size_c);
        ssSetSFcnParam(rts, 1, (mxArray*)Wind_songweiwei_P.SFunction_P2_Size_e);
        ssSetSFcnParam(rts, 2, (mxArray*)Wind_songweiwei_P.SFunction_P3_Size_m);
        ssSetSFcnParam(rts, 3, (mxArray*)Wind_songweiwei_P.SFunction_P4_Size_n);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *) &Wind_songweiwei_DW.SFunction_RWORK);
      ssSetIWork(rts, (int_T *) &Wind_songweiwei_DW.SFunction_IWORK_fk);
      ssSetPWork(rts, (void **) &Wind_songweiwei_DW.SFunction_PWORK_g);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn2.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &Wind_songweiwei_M->NonInlinedSFcns.Sfcn2.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &Wind_songweiwei_DW.SFunction_RWORK);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &Wind_songweiwei_DW.SFunction_IWORK_fk);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &Wind_songweiwei_DW.SFunction_PWORK_g);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 2.0E-6);
      ssSetOffsetTime(rts, 0, 0.0);
      sfcnTsMap[0] = 1;

      /* set compiled values of dynamic vector attributes */
      ssSetInputPortWidth(rts, 0, 1);
      ssSetInputPortDataType(rts, 0, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 0, 0);
      ssSetInputPortFrameData(rts, 0, 0);
      ssSetInputPortWidth(rts, 1, 1);
      ssSetInputPortDataType(rts, 1, SS_DOUBLE);
      ssSetInputPortComplexSignal(rts, 1, 0);
      ssSetInputPortFrameData(rts, 1, 0);
      ssSetOutputPortWidth(rts, 0, 1);
      ssSetOutputPortDataType(rts, 0, SS_DOUBLE);
      ssSetOutputPortComplexSignal(rts, 0, 0);
      ssSetOutputPortFrameData(rts, 0, 0);
      ssSetNumIWork(rts, 1);
      ssSetNumPWork(rts, 1);
      ssSetNumNonsampledZCs(rts, 0);

      /* Update connectivity flags for each port */
      _ssSetInputPortConnected(rts, 0, 1);
      _ssSetInputPortConnected(rts, 1, 1);
      _ssSetOutputPortConnected(rts, 0, 1);
      _ssSetOutputPortBeingMerged(rts, 0, 0);

      /* Update the BufferDstPort flags for each input port */
      ssSetInputPortBufferDstPort(rts, 0, -1);
      ssSetInputPortBufferDstPort(rts, 1, -1);
    }
  }

  {
    int32_T i;

    /* InitializeConditions for Enabled SubSystem: '<S353>/Tail' */
    /* Start for Enabled SubSystem: '<S353>/Tail' */
    for (i = 0; i < 6; i++) {
      /* InitializeConditions for DiscreteIntegrator: '<S354>/Discrete-Time Integrator' */
      Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_b4[i] =
        Wind_songweiwei_P.DiscreteTimeIntegrator_IC_j;
      Wind_songweiwei_DW.DiscreteTimeIntegrator_PrevRese[i] = 2;

      /* InitializeConditions for UnitDelay: '<S354>/Unit Delay' */
      Wind_songweiwei_DW.UnitDelay_DSTATE_fv[i] =
        Wind_songweiwei_P.UnitDelay_InitialCondition_nz;

      /* VirtualOutportStart for Outport: '<S354>/itail' */
      Wind_songweiwei_B.Product[i] = Wind_songweiwei_P.itail_Y0;
    }

    /* End of Start for SubSystem: '<S353>/Tail' */
    /* End of InitializeConditions for SubSystem: '<S353>/Tail' */

    /* Start for Enabled SubSystem: '<S331>/Saturation' */
    /* InitializeConditions for Enabled SubSystem: '<S341>/Lmq_sat' */
    /* InitializeConditions for UnitDelay: '<S344>/Lmq_sat' */
    Wind_songweiwei_DW.Lmq_sat_DSTATE =
      Wind_songweiwei_P.Lmq_sat_InitialCondition;

    /* End of InitializeConditions for SubSystem: '<S341>/Lmq_sat' */
    /* End of Start for SubSystem: '<S331>/Saturation' */

    /* InitializeConditions for Enabled SubSystem: '<S331>/Saturation' */
    /* InitializeConditions for UnitDelay: '<S343>/Lmd_sat' */
    Wind_songweiwei_DW.Lmd_sat_DSTATE =
      Wind_songweiwei_P.Lmd_sat_InitialCondition;

    /* End of InitializeConditions for SubSystem: '<S331>/Saturation' */

    /* InitializeConditions for Enabled SubSystem: '<S10>/Signal generator' */
    /* InitializeConditions for DiscreteIntegrator: '<S12>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_l =
      Wind_songweiwei_P.DiscreteTimeIntegrator_IC;

    /* InitializeConditions for DiscreteIntegrator: '<S16>/Discrete-Time Integrator1' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTAT_a =
      Wind_songweiwei_P.DiscreteTimeIntegrator1_IC;

    /* InitializeConditions for UnitDelay: '<S16>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_gd =
      Wind_songweiwei_P.UnitDelay_InitialCondition;

    /* End of InitializeConditions for SubSystem: '<S10>/Signal generator' */

    /* Start for Enabled SubSystem: '<S10>/Signal generator' */
    /* VirtualOutportStart for Outport: '<S12>/timer' */
    Wind_songweiwei_B.LookUpTable = Wind_songweiwei_P.timer_Y0;

    /* VirtualOutportStart for Outport: '<S12>/selector' */
    Wind_songweiwei_B.LogicalOperator1 = Wind_songweiwei_P.selector_Y0;

    /* VirtualOutportStart for Outport: '<S12>/magnitude' */
    Wind_songweiwei_B.Switch2 = Wind_songweiwei_P.magnitude_Y0;

    /* VirtualOutportStart for Outport: '<S12>/frequency' */
    Wind_songweiwei_B.DiscreteTimeIntegrator = Wind_songweiwei_P.frequency_Y0;

    /* VirtualOutportStart for Outport: '<S12>/phase' */
    Wind_songweiwei_B.Switch3 = Wind_songweiwei_P.phase_Y0;

    /* End of Start for SubSystem: '<S10>/Signal generator' */

    /* Start for Enabled SubSystem: '<S10>/Harmonic Generator' */
    /* VirtualOutportStart for Outport: '<S11>/Out1' */
    Wind_songweiwei_B.Product1[0] = Wind_songweiwei_P.Out1_Y0;
    Wind_songweiwei_B.Product1[1] = Wind_songweiwei_P.Out1_Y0;
    Wind_songweiwei_B.Product1[2] = Wind_songweiwei_P.Out1_Y0;

    /* VirtualOutportStart for Outport: '<S11>/Out2' */
    Wind_songweiwei_B.Product1_k[0] = Wind_songweiwei_P.Out2_Y0;
    Wind_songweiwei_B.Product1_k[1] = Wind_songweiwei_P.Out2_Y0;
    Wind_songweiwei_B.Product1_k[2] = Wind_songweiwei_P.Out2_Y0;

    /* End of Start for SubSystem: '<S10>/Harmonic Generator' */

    /* S-Function block: <S357>/State-Space */
    {
      Wind_songweiwei_DW.StateSpace_PWORK.AS = (real_T*)calloc(45 * 45, sizeof
        (real_T));
      Wind_songweiwei_DW.StateSpace_PWORK.BS = (real_T*)calloc(45 * 24, sizeof
        (real_T));
      Wind_songweiwei_DW.StateSpace_PWORK.CS = (real_T*)calloc(50 * 45, sizeof
        (real_T));
      Wind_songweiwei_DW.StateSpace_PWORK.DS = (real_T*)calloc(50 * 24, sizeof
        (real_T));
      Wind_songweiwei_DW.StateSpace_PWORK.DX_COL = (real_T*)calloc(50, sizeof
        (real_T));
      Wind_songweiwei_DW.StateSpace_PWORK.TMP2 = (real_T*)calloc(24, sizeof
        (real_T));
      Wind_songweiwei_DW.StateSpace_PWORK.BD_COL = (real_T*)calloc(45, sizeof
        (real_T));
      Wind_songweiwei_DW.StateSpace_PWORK.TMP1 = (real_T*)calloc(45, sizeof
        (real_T));
      Wind_songweiwei_DW.StateSpace_PWORK.XTMP = (real_T*)calloc(45, sizeof
        (real_T));
      Wind_songweiwei_DW.StateSpace_PWORK.CHOPPER = (int_T*)calloc(50, sizeof
        (int_T));
      Wind_songweiwei_DW.StateSpace_PWORK.SWITCH_STATUS = (int_T*)calloc(19,
        sizeof(int_T));
      Wind_songweiwei_DW.StateSpace_PWORK.SW_CHG = (int_T*)calloc(19, sizeof
        (int_T));
      Wind_songweiwei_DW.StateSpace_PWORK.G_STATE = (int_T*)calloc(19, sizeof
        (int_T));
      Wind_songweiwei_DW.StateSpace_PWORK.Y_SWITCH = (real_T*)calloc(19, sizeof
        (real_T));
      Wind_songweiwei_DW.StateSpace_PWORK.SWITCH_TYPES = (int_T*)calloc(19,
        sizeof(int_T));
      Wind_songweiwei_DW.StateSpace_PWORK.IDX_OUT_SW = (int_T*)calloc(19, sizeof
        (int_T));
      Wind_songweiwei_DW.StateSpace_PWORK.SWITCH_STATUS_INIT = (int_T*)calloc(19,
        sizeof(int_T));
    }

    /* Start for Enabled SubSystem: '<S274>/Subsystem1' */
    Wind_songweiwe_Subsystem1_Start(&Wind_songweiwei_B.Subsystem1,
      (P_Subsystem1_Wind_songweiwei_T *)&Wind_songweiwei_P.Subsystem1);

    /* End of Start for SubSystem: '<S274>/Subsystem1' */

    /* Start for Enabled SubSystem: '<S274>/Subsystem - pi//2 delay' */
    Wind_so_Subsystempi2delay_Start(&Wind_songweiwei_B.Subsystempi2delay,
      (P_Subsystempi2delay_Wind_song_T *)&Wind_songweiwei_P.Subsystempi2delay);

    /* End of Start for SubSystem: '<S274>/Subsystem - pi//2 delay' */

    /* Start for Enabled SubSystem: '<S262>/Subsystem1' */
    Wind_songweiwe_Subsystem1_Start(&Wind_songweiwei_B.Subsystem1_l,
      (P_Subsystem1_Wind_songweiwei_T *)&Wind_songweiwei_P.Subsystem1_l);

    /* End of Start for SubSystem: '<S262>/Subsystem1' */

    /* Start for Enabled SubSystem: '<S262>/Subsystem - pi//2 delay' */
    Wind_so_Subsystempi2delay_Start(&Wind_songweiwei_B.Subsystempi2delay_j,
      (P_Subsystempi2delay_Wind_song_T *)&Wind_songweiwei_P.Subsystempi2delay_j);

    /* End of Start for SubSystem: '<S262>/Subsystem - pi//2 delay' */

    /* Start for Enabled SubSystem: '<S160>/Subsystem1' */
    Wind_songweiwe_Subsystem1_Start(&Wind_songweiwei_B.Subsystem1_o,
      (P_Subsystem1_Wind_songweiwei_T *)&Wind_songweiwei_P.Subsystem1_o);

    /* End of Start for SubSystem: '<S160>/Subsystem1' */

    /* Start for Enabled SubSystem: '<S160>/Subsystem - pi//2 delay' */
    Wind_so_Subsystempi2delay_Start(&Wind_songweiwei_B.Subsystempi2delay_k,
      (P_Subsystempi2delay_Wind_song_T *)&Wind_songweiwei_P.Subsystempi2delay_k);

    /* End of Start for SubSystem: '<S160>/Subsystem - pi//2 delay' */

    /* S-Function Block: <S159>/S-Function  */
    {
      static real_T dvtd_buffer[1 * 11113];
      Wind_songweiwei_DW.SFunction_IWORK.bufSz = 11113;
      Wind_songweiwei_DW.SFunction_PWORK.uBuffers = (void *)&dvtd_buffer[0];
    }

    {
      Wind_songweiwei_DW.SFunction_IWORK.maxDiscrDelay =
        Wind_songweiwei_DW.SFunction_IWORK.bufSz - 1;

      /* Assign default sample(s) */
      /* Single initial value */
      if (Wind_songweiwei_DW.SFunction_PWORK.uBuffers != NULL)
        *((real_T *)Wind_songweiwei_DW.SFunction_PWORK.uBuffers) = (real_T)0.0;

      /* Set work values */
      Wind_songweiwei_DW.SFunction_IWORK.indBeg = 0;
      Wind_songweiwei_DW.SFunction_IWORK.indEnd = 1;
    }

    /* S-Function Block: <S157>/S-Function  */
    {
      static real_T dvtd_buffer[1 * 11113];
      Wind_songweiwei_DW.SFunction_IWORK_f.bufSz = 11113;
      Wind_songweiwei_DW.SFunction_PWORK_c.uBuffers = (void *)&dvtd_buffer[0];
    }

    {
      Wind_songweiwei_DW.SFunction_IWORK_f.maxDiscrDelay =
        Wind_songweiwei_DW.SFunction_IWORK_f.bufSz - 1;

      /* Assign default sample(s) */
      /* Single initial value */
      if (Wind_songweiwei_DW.SFunction_PWORK_c.uBuffers != NULL)
        *((real_T *)Wind_songweiwei_DW.SFunction_PWORK_c.uBuffers) = (real_T)0.0;

      /* Set work values */
      Wind_songweiwei_DW.SFunction_IWORK_f.indBeg = 0;
      Wind_songweiwei_DW.SFunction_IWORK_f.indEnd = 1;
    }

    /* Start for Enabled SubSystem: '<S268>/Subsystem1' */
    Wind_songweiwe_Subsystem1_Start(&Wind_songweiwei_B.Subsystem1_p,
      (P_Subsystem1_Wind_songweiwei_T *)&Wind_songweiwei_P.Subsystem1_p);

    /* End of Start for SubSystem: '<S268>/Subsystem1' */

    /* Start for Enabled SubSystem: '<S268>/Subsystem - pi//2 delay' */
    Wind_so_Subsystempi2delay_Start(&Wind_songweiwei_B.Subsystempi2delay_jr,
      (P_Subsystempi2delay_Wind_song_T *)&Wind_songweiwei_P.Subsystempi2delay_jr);

    /* End of Start for SubSystem: '<S268>/Subsystem - pi//2 delay' */

    /* Start for Enabled SubSystem: '<S174>/Subsystem1' */
    Wind_songweiwe_Subsystem1_Start(&Wind_songweiwei_B.Subsystem1_oh,
      (P_Subsystem1_Wind_songweiwei_T *)&Wind_songweiwei_P.Subsystem1_oh);

    /* End of Start for SubSystem: '<S174>/Subsystem1' */

    /* Start for Enabled SubSystem: '<S174>/Subsystem - pi//2 delay' */
    Wind_so_Subsystempi2delay_Start(&Wind_songweiwei_B.Subsystempi2delay_h,
      (P_Subsystempi2delay_Wind_song_T *)&Wind_songweiwei_P.Subsystempi2delay_h);

    /* End of Start for SubSystem: '<S174>/Subsystem - pi//2 delay' */

    /* S-Function Block: <S173>/S-Function  */
    {
      static real_T dvtd_buffer[1 * 11113];
      Wind_songweiwei_DW.SFunction_IWORK_c.bufSz = 11113;
      Wind_songweiwei_DW.SFunction_PWORK_l.uBuffers = (void *)&dvtd_buffer[0];
    }

    {
      Wind_songweiwei_DW.SFunction_IWORK_c.maxDiscrDelay =
        Wind_songweiwei_DW.SFunction_IWORK_c.bufSz - 1;

      /* Assign default sample(s) */
      /* Single initial value */
      if (Wind_songweiwei_DW.SFunction_PWORK_l.uBuffers != NULL)
        *((real_T *)Wind_songweiwei_DW.SFunction_PWORK_l.uBuffers) = (real_T)0.0;

      /* Set work values */
      Wind_songweiwei_DW.SFunction_IWORK_c.indBeg = 0;
      Wind_songweiwei_DW.SFunction_IWORK_c.indEnd = 1;
    }

    /* S-Function Block: <S171>/S-Function  */
    {
      static real_T dvtd_buffer[1 * 11113];
      Wind_songweiwei_DW.SFunction_IWORK_k.bufSz = 11113;
      Wind_songweiwei_DW.SFunction_PWORK_p.uBuffers = (void *)&dvtd_buffer[0];
    }

    {
      Wind_songweiwei_DW.SFunction_IWORK_k.maxDiscrDelay =
        Wind_songweiwei_DW.SFunction_IWORK_k.bufSz - 1;

      /* Assign default sample(s) */
      /* Single initial value */
      if (Wind_songweiwei_DW.SFunction_PWORK_p.uBuffers != NULL)
        *((real_T *)Wind_songweiwei_DW.SFunction_PWORK_p.uBuffers) = (real_T)0.0;

      /* Set work values */
      Wind_songweiwei_DW.SFunction_IWORK_k.indBeg = 0;
      Wind_songweiwei_DW.SFunction_IWORK_k.indEnd = 1;
    }

    /* Start for Enabled SubSystem: '<S137>/Subsystem1' */
    /* VirtualOutportStart for Outport: '<S141>/alpha_beta' */
    Wind_songweiwei_B.Fcn = Wind_songweiwei_P.alpha_beta_Y0_o[0];
    Wind_songweiwei_B.Fcn1 = Wind_songweiwei_P.alpha_beta_Y0_o[1];

    /* End of Start for SubSystem: '<S137>/Subsystem1' */

    /* Start for Enabled SubSystem: '<S137>/Subsystem - pi//2 delay' */
    /* VirtualOutportStart for Outport: '<S140>/alpha_beta' */
    Wind_songweiwei_B.Fcn_c = Wind_songweiwei_P.alpha_beta_Y0[0];
    Wind_songweiwei_B.Fcn1_b = Wind_songweiwei_P.alpha_beta_Y0[1];

    /* End of Start for SubSystem: '<S137>/Subsystem - pi//2 delay' */

    /* Start for Enabled SubSystem: '<S225>/Automatic Gain Control' */

    /* Start for Enabled SubSystem: '<S241>/Subsystem1' */
    Wind_songweiwe_Subsystem1_Start(&Wind_songweiwei_B.Subsystem1_lp,
      (P_Subsystem1_Wind_songweiwei_T *)&Wind_songweiwei_P.Subsystem1_lp);

    /* End of Start for SubSystem: '<S241>/Subsystem1' */

    /* Start for Enabled SubSystem: '<S241>/Subsystem - pi//2 delay' */
    Wind_so_Subsystempi2delay_Start(&Wind_songweiwei_B.Subsystempi2delay_dg,
      (P_Subsystempi2delay_Wind_song_T *)&Wind_songweiwei_P.Subsystempi2delay_dg);

    /* End of Start for SubSystem: '<S241>/Subsystem - pi//2 delay' */

    /* Level2 S-Function Block: '<S237>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Wind_songweiwei_M->childSfunctions[0];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Level2 S-Function Block: '<S240>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Wind_songweiwei_M->childSfunctions[1];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* End of Start for SubSystem: '<S225>/Automatic Gain Control' */

    /* InitializeConditions for Enabled SubSystem: '<S225>/Automatic Gain Control' */
    /* InitializeConditions for DiscreteIntegrator: '<S235>/Integ4' */
    Wind_songweiwei_DW.Integ4_DSTATE_nn = Wind_songweiwei_P.Integ4_IC;

    /* Level2 S-Function Block: '<S237>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Wind_songweiwei_M->childSfunctions[0];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S236>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_i =
      Wind_songweiwei_P.UnitDelay_InitialCondition_b;

    /* InitializeConditions for UnitDelay: '<S235>/Unit Delay1' */
    Wind_songweiwei_DW.UnitDelay1_DSTATE_f =
      Wind_songweiwei_P.UnitDelay1_InitialCondition;

    /* InitializeConditions for DiscreteIntegrator: '<S238>/Integ4' */
    Wind_songweiwei_DW.Integ4_DSTATE_c = Wind_songweiwei_P.Integ4_IC_k;

    /* Level2 S-Function Block: '<S240>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Wind_songweiwei_M->childSfunctions[1];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S239>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_g =
      Wind_songweiwei_P.UnitDelay_InitialCondition_n;

    /* InitializeConditions for UnitDelay: '<S238>/Unit Delay1' */
    Wind_songweiwei_DW.UnitDelay1_DSTATE_c =
      Wind_songweiwei_P.UnitDelay1_InitialCondition_f;

    /* End of InitializeConditions for SubSystem: '<S225>/Automatic Gain Control' */

    /* Start for Enabled SubSystem: '<S225>/Automatic Gain Control' */
    /* VirtualOutportStart for Outport: '<S226>/Gain' */
    Wind_songweiwei_B.MathFunction = Wind_songweiwei_P.Gain_Y0;

    /* End of Start for SubSystem: '<S225>/Automatic Gain Control' */

    /* Start for Enabled SubSystem: '<S254>/Subsystem1' */
    Wind_songweiwe_Subsystem1_Start(&Wind_songweiwei_B.Subsystem1_l1,
      (P_Subsystem1_Wind_songweiwei_T *)&Wind_songweiwei_P.Subsystem1_l1);

    /* End of Start for SubSystem: '<S254>/Subsystem1' */

    /* Start for Enabled SubSystem: '<S254>/Subsystem - pi//2 delay' */
    Wind_so_Subsystempi2delay_Start(&Wind_songweiwei_B.Subsystempi2delay_p,
      (P_Subsystempi2delay_Wind_song_T *)&Wind_songweiwei_P.Subsystempi2delay_p);

    /* End of Start for SubSystem: '<S254>/Subsystem - pi//2 delay' */

    /* Level2 S-Function Block: '<S253>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Wind_songweiwei_M->childSfunctions[2];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S353>/Unit Delay' */
    for (i = 0; i < 6; i++) {
      Wind_songweiwei_DW.UnitDelay_DSTATE[i] =
        Wind_songweiwei_P.UnitDelay_InitialCondition_o;
    }

    /* End of InitializeConditions for UnitDelay: '<S353>/Unit Delay' */

    /* InitializeConditions for UnitDelay: '<S67>/Unit Delay6' */
    Wind_songweiwei_DW.UnitDelay6_DSTATE[0] =
      Wind_songweiwei_P.UnitDelay6_InitialCondition;
    Wind_songweiwei_DW.UnitDelay6_DSTATE[1] =
      Wind_songweiwei_P.UnitDelay6_InitialCondition;
    Wind_songweiwei_DW.UnitDelay6_DSTATE[2] =
      Wind_songweiwei_P.UnitDelay6_InitialCondition;

    /* InitializeConditions for DiscreteIntegrator: '<S329>/Rotor angle dthetae' */
    Wind_songweiwei_DW.Rotorangledthetae_DSTATE =
      Wind_songweiwei_P.Rotorangledthetae_IC;

    /* InitializeConditions for S-Function (sfun_spssw_discc): '<S357>/State-Space' */
    {
      real_T *As = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.AS;
      real_T *Bs = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.BS;
      real_T *Cs = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.CS;
      real_T *Ds = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.DS;
      real_T *X0 = (real_T*)&Wind_songweiwei_DW.StateSpace_DSTATE[0];
      int_T *Chopper = (int_T*) Wind_songweiwei_DW.StateSpace_PWORK.CHOPPER;
      X0[0] = 5.5000000000916672E+8;
      X0[1] = 0.0;
      X0[2] = 0.010343729157215265;
      X0[3] = -0.005011879576829914;
      X0[4] = -515.24882632416256;
      X0[5] = 1034.3287682306886;
      X0[6] = -519.07994180062735;
      X0[7] = -2.8081452229035073E+6;
      X0[8] = 4.4547848273858326E+6;
      X0[9] = 5.6441334903563783E+6;
      X0[10] = 3.2493609751496492E+6;
      X0[11] = -1.2054238522361817E+6;
      X0[12] = -2.8359882674528738E+6;
      X0[13] = -728577.7162419646;
      X0[14] = -572573.82765827142;
      X0[15] = 1.3011515439005033E+6;
      X0[16] = 360749.08178433275;
      X0[17] = -770390.62551393092;
      X0[18] = 409641.5437295984;
      X0[19] = -5.1119100856303921E+9;
      X0[20] = -5.08717223485545E+9;
      X0[21] = 1.0199082320485846E+10;
      X0[22] = 1.9879962815711711E-8;
      X0[23] = -5.121245897425724E+9;
      X0[24] = -5.0841736683591938E+9;
      X0[25] = 1.0205419565784918E+10;
      X0[26] = -1.8450142426524452E-6;
      X0[27] = 0.0;
      X0[28] = 0.0;
      X0[29] = 0.0;
      X0[30] = 0.0;
      X0[31] = 0.0;
      X0[32] = 0.0;
      X0[33] = -204388.09649707848;
      X0[34] = 2035.286272495985;
      X0[35] = 506.14807709193718;
      X0[36] = -5.029883153969422;
      X0[37] = 203881.9484199874;
      X0[38] = -2030.256389349699;
      X0[39] = 1.1323713886544227E+6;
      X0[40] = -319640.69010972045;
      X0[41] = 1.2637731574098954E+6;
      X0[42] = 160104.13927908384;
      X0[43] = -2.396144546064266E+6;
      X0[44] = 159536.55083063786;

      /* Copy and transpose A and B */
      As[0] = 0.99999999996666666;
      As[1] = 0.0;
      As[2] = 0.0;
      As[3] = 0.0;
      As[4] = 0.0;
      As[5] = 0.0;
      As[6] = 0.0;
      As[7] = 0.0;
      As[8] = 0.0;
      As[9] = 0.0;
      As[10] = 0.0;
      As[11] = 0.0;
      As[12] = 0.0;
      As[13] = 0.0;
      As[14] = 0.0;
      As[15] = 0.0;
      As[16] = 0.0;
      As[17] = 0.0;
      As[18] = 0.0;
      As[19] = 0.0;
      As[20] = 0.0;
      As[21] = 0.0;
      As[22] = 0.0;
      As[23] = 0.0;
      As[24] = 0.0;
      As[25] = 0.0;
      As[26] = 0.0;
      As[27] = 0.0;
      As[28] = 0.0;
      As[29] = 0.0;
      As[30] = 0.0;
      As[31] = 0.0;
      As[32] = 0.0;
      As[33] = 0.0;
      As[34] = 0.0;
      As[35] = 0.0;
      As[36] = 0.0;
      As[37] = 0.0;
      As[38] = 0.0;
      As[39] = 0.0;
      As[40] = 0.0;
      As[41] = 0.0;
      As[42] = 0.0;
      As[43] = 0.0;
      As[44] = 0.0;
      As[45] = 0.0;
      As[46] = -0.99782057144712388;
      As[47] = 0.0;
      As[48] = 0.0;
      As[49] = 0.0;
      As[50] = 0.0;
      As[51] = 0.0;
      As[52] = 0.0;
      As[53] = 0.0;
      As[54] = 0.0;
      As[55] = 0.0;
      As[56] = 0.0;
      As[57] = 0.0;
      As[58] = 0.0;
      As[59] = 0.0;
      As[60] = 0.0;
      As[61] = 0.0;
      As[62] = 0.0;
      As[63] = 0.0;
      As[64] = 0.0;
      As[65] = 0.0;
      As[66] = 0.0;
      As[67] = 0.0;
      As[68] = 0.0;
      As[69] = 0.0;
      As[70] = 0.0;
      As[71] = 0.0;
      As[72] = -3.0269841012171553E-6;
      As[73] = -3.0269841012171575E-6;
      As[74] = -3.0269841012171575E-6;
      As[75] = -3.0269841012171558E-6;
      As[76] = -3.0269841012171575E-6;
      As[77] = -3.0269841012171579E-6;
      As[78] = 0.0;
      As[79] = 0.0;
      As[80] = 0.0;
      As[81] = 0.0;
      As[82] = 0.0;
      As[83] = 0.0;
      As[84] = 0.0;
      As[85] = 0.0;
      As[86] = 0.0;
      As[87] = 0.0;
      As[88] = 0.0;
      As[89] = 0.0;
      As[90] = 0.0;
      As[91] = 0.0;
      As[92] = -0.99976323533977529;
      As[93] = 2.673025536520679E-20;
      As[94] = -6.6658762481084077E-6;
      As[95] = 1.3331752496216817E-5;
      As[96] = -6.6658762481084077E-6;
      As[97] = 1.545148130846668E-18;
      As[98] = -1.8450800379097555E-17;
      As[99] = -3.0902962616933368E-18;
      As[100] = 2.9443723774288979E-33;
      As[101] = 1.8450800379097552E-17;
      As[102] = 1.5451481308466686E-18;
      As[103] = -5.7269177444198917E-17;
      As[104] = 2.863458872209949E-17;
      As[105] = 2.8634588722099465E-17;
      As[106] = -1.2376777139489712E-13;
      As[107] = 1.237677713948972E-13;
      As[108] = 2.614001792366622E-29;
      As[109] = -3.9240443053650868E-18;
      As[110] = 3.9240443053650883E-18;
      As[111] = 1.0513461515531762E-34;
      As[112] = -7.7144006102350382E-34;
      As[113] = -2.1043233889183243E-14;
      As[114] = 2.1043233889183243E-14;
      As[115] = -3.0746639168199234E-31;
      As[116] = 5.6407808750203227E-30;
      As[117] = 0.0;
      As[118] = 0.0;
      As[119] = 0.0;
      As[120] = 0.0;
      As[121] = 0.0;
      As[122] = 0.0;
      As[123] = -6.9786670387901659E-11;
      As[124] = 1.0498282581022514E-12;
      As[125] = 1.3957334077580327E-10;
      As[126] = -2.0996565162027827E-12;
      As[127] = -6.9786670387901659E-11;
      As[128] = 1.0498282581022514E-12;
      As[129] = 2.5507039832759886E-17;
      As[130] = 5.7270170990784434E-17;
      As[131] = -1.2753519916379945E-17;
      As[132] = -2.8635085495392235E-17;
      As[133] = -1.2753519916379943E-17;
      As[134] = -2.8635085495392211E-17;
      As[135] = 0.0;
      As[136] = 0.0;
      As[137] = -5.3087497934978641E-26;
      As[138] = -0.99976323533977529;
      As[139] = -6.6658762481084069E-6;
      As[140] = -6.6658762481084077E-6;
      As[141] = 1.3331752496216815E-5;
      As[142] = 1.5451481308466674E-18;
      As[143] = -1.61121355861877E-34;
      As[144] = 1.5451481308466688E-18;
      As[145] = -1.8450800379097545E-17;
      As[146] = -1.8450800379097564E-17;
      As[147] = -3.0902962616933364E-18;
      As[148] = 2.8634588722099471E-17;
      As[149] = -5.7269177444198942E-17;
      As[150] = 2.8634588722099459E-17;
      As[151] = -8.9629628349702546E-31;
      As[152] = -1.237677713948972E-13;
      As[153] = 1.237677713948971E-13;
      As[154] = 5.9264491192899615E-34;
      As[155] = -3.9240443053650883E-18;
      As[156] = 3.924044305365086E-18;
      As[157] = -9.8496945298457112E-34;
      As[158] = -8.6964060928513564E-32;
      As[159] = -2.1043233889183246E-14;
      As[160] = 2.1043233889183237E-14;
      As[161] = 5.9973635249759482E-34;
      As[162] = 0.0;
      As[163] = 0.0;
      As[164] = 0.0;
      As[165] = 0.0;
      As[166] = 0.0;
      As[167] = 0.0;
      As[168] = -6.9786670387901659E-11;
      As[169] = 1.0498282581022514E-12;
      As[170] = -6.9786670387901607E-11;
      As[171] = 1.0498282581080088E-12;
      As[172] = 1.3957334077580337E-10;
      As[173] = -2.0996565162045029E-12;
      As[174] = -1.2753519916379949E-17;
      As[175] = -2.8635085495392229E-17;
      As[176] = 2.5507039832759892E-17;
      As[177] = 5.7270170990784446E-17;
      As[178] = -1.275351991637994E-17;
      As[179] = -2.8635085495392205E-17;
      As[180] = 0.0;
      As[181] = 0.0;
      As[182] = 3.9332241969608592E-8;
      As[183] = 3.9332241969581487E-8;
      As[184] = 0.99999999446291277;
      As[185] = 1.1073605787394254E-9;
      As[186] = 1.1073605787394258E-9;
      As[187] = -5.7179795176263562E-15;
      As[188] = 3.4139541881231387E-14;
      As[189] = 2.8589897588131781E-15;
      As[190] = 3.4139541881231368E-14;
      As[191] = 7.0728870751135686E-30;
      As[192] = 2.8589897588131781E-15;
      As[193] = 5.2982619769569255E-14;
      As[194] = 5.29826197695693E-14;
      As[195] = -1.0596523953913849E-13;
      As[196] = 2.2900768141579287E-10;
      As[197] = -1.7635893023419171E-26;
      As[198] = -2.290076814157929E-10;
      As[199] = 7.2606646949898488E-15;
      As[200] = 7.9366019371554937E-31;
      As[201] = -7.2606646949898472E-15;
      As[202] = 5.1447739500351989E-31;
      As[203] = 3.8936325249618069E-11;
      As[204] = -3.2224536411585141E-27;
      As[205] = -3.8936325249618069E-11;
      As[206] = 4.1013357394524021E-27;
      As[207] = 0.0;
      As[208] = 0.0;
      As[209] = 0.0;
      As[210] = 0.0;
      As[211] = 0.0;
      As[212] = 0.0;
      As[213] = 3.8740596207292961E-7;
      As[214] = -5.8274963272620847E-9;
      As[215] = 2.6839079677817859E-11;
      As[216] = -2.2150761360648511E-16;
      As[217] = 2.6839079677817837E-11;
      As[218] = -2.2150761776822213E-16;
      As[219] = -2.3597855831317992E-14;
      As[220] = -5.2983538949891245E-14;
      As[221] = -2.3597855831317998E-14;
      As[222] = -5.2983538949891277E-14;
      As[223] = 4.7195711662636009E-14;
      As[224] = 1.0596707789978244E-13;
      As[225] = 0.0;
      As[226] = 0.0;
      As[227] = -3.9332241969599348E-8;
      As[228] = -4.4405529224011993E-24;
      As[229] = 1.1073605787394258E-9;
      As[230] = 0.99999999446291277;
      As[231] = 1.1073605787394254E-9;
      As[232] = 2.8589897588131789E-15;
      As[233] = -3.4139541881231375E-14;
      As[234] = -5.7179795176263569E-15;
      As[235] = -1.2293648476959042E-29;
      As[236] = 3.4139541881231394E-14;
      As[237] = 2.8589897588131777E-15;
      As[238] = -1.0596523953913852E-13;
      As[239] = 5.2982619769569293E-14;
      As[240] = 5.2982619769569287E-14;
      As[241] = -2.290076814157928E-10;
      As[242] = 2.29007681415793E-10;
      As[243] = 1.4110080650329069E-25;
      As[244] = -7.2606646949898488E-15;
      As[245] = 7.2606646949898488E-15;
      As[246] = 1.2495153588186624E-30;
      As[247] = 2.7255267380771876E-30;
      As[248] = -3.8936325249618075E-11;
      As[249] = 3.8936325249618075E-11;
      As[250] = 3.2480916205781966E-27;
      As[251] = 2.0796455056499183E-26;
      As[252] = 0.0;
      As[253] = 0.0;
      As[254] = 0.0;
      As[255] = 0.0;
      As[256] = 0.0;
      As[257] = 0.0;
      As[258] = 2.6839079677817837E-11;
      As[259] = -2.2150761855526038E-16;
      As[260] = 3.8740596207292982E-7;
      As[261] = -5.8274963272385026E-9;
      As[262] = 2.683907967781785E-11;
      As[263] = -2.2150761776741434E-16;
      As[264] = 4.7195711662635991E-14;
      As[265] = 1.0596707789978252E-13;
      As[266] = -2.3597855831318008E-14;
      As[267] = -5.298353894989127E-14;
      As[268] = -2.3597855831318E-14;
      As[269] = -5.2983538949891251E-14;
      As[270] = 0.0;
      As[271] = 0.0;
      As[272] = -2.81702131951354E-30;
      As[273] = -3.9332241969572243E-8;
      As[274] = 1.1073605787394256E-9;
      As[275] = 1.1073605787394256E-9;
      As[276] = 0.99999999446291277;
      As[277] = 2.8589897588131777E-15;
      As[278] = -1.2186710019484182E-29;
      As[279] = 2.8589897588131785E-15;
      As[280] = -3.4139541881231356E-14;
      As[281] = -3.4139541881231394E-14;
      As[282] = -5.7179795176263554E-15;
      As[283] = 5.2982619769569281E-14;
      As[284] = -1.0596523953913855E-13;
      As[285] = 5.2982619769569255E-14;
      As[286] = -8.8518715002485485E-26;
      As[287] = -2.2900768141579295E-10;
      As[288] = 2.2900768141579269E-10;
      As[289] = 1.6031886738208659E-32;
      As[290] = -7.2606646949898472E-15;
      As[291] = 7.2606646949898441E-15;
      As[292] = -4.341519186187369E-30;
      As[293] = -1.1981551434936453E-30;
      As[294] = -3.8936325249618062E-11;
      As[295] = 3.8936325249618062E-11;
      As[296] = -1.789333469790389E-26;
      As[297] = 0.0;
      As[298] = 0.0;
      As[299] = 0.0;
      As[300] = 0.0;
      As[301] = 0.0;
      As[302] = 0.0;
      As[303] = 2.6839079677817843E-11;
      As[304] = -2.2150761318511619E-16;
      As[305] = 2.6839079677817817E-11;
      As[306] = -2.2150761436242202E-16;
      As[307] = 3.8740596207292982E-7;
      As[308] = -5.8274963272620847E-9;
      As[309] = -2.3597855831317998E-14;
      As[310] = -5.298353894989127E-14;
      As[311] = 4.7195711662635997E-14;
      As[312] = 1.0596707789978253E-13;
      As[313] = -2.3597855831317995E-14;
      As[314] = -5.2983538949891226E-14;
      As[315] = 0.0;
      As[316] = 0.0;
      As[317] = 1.9862415655548358E-23;
      As[318] = 1.9862415651806357E-23;
      As[319] = 1.2456989112038275E-17;
      As[320] = -6.2284945607888487E-18;
      As[321] = -6.228494557262179E-18;
      As[322] = 0.99934715189965651;
      As[323] = -0.00060213100710638232;
      As[324] = -5.0424999514343813E-5;
      As[325] = -0.00060213100710638265;
      As[326] = 5.3976351761363179E-20;
      As[327] = -5.0424999514343827E-5;
      As[328] = 3.5290430969130674E-7;
      As[329] = 3.5290430969064648E-7;
      As[330] = -7.0580861938185739E-7;
      As[331] = -2.8444276795837083E-7;
      As[332] = 1.34295633419945E-20;
      As[333] = 2.8444276795800972E-7;
      As[334] = 4.836151691307132E-8;
      As[335] = -2.5246418620254425E-21;
      As[336] = -4.8361516913037478E-8;
      As[337] = 1.2990012131473517E-20;
      As[338] = 9.0301488505048827E-12;
      As[339] = 4.2979483956175405E-21;
      As[340] = -9.0301488418405761E-12;
      As[341] = 1.380032924924162E-24;
      As[342] = 0.0;
      As[343] = 0.0;
      As[344] = 0.0;
      As[345] = 0.0;
      As[346] = 0.0;
      As[347] = 0.0;
      As[348] = -1.2443710294656139E-11;
      As[349] = 2.1849811005502063E-17;
      As[350] = 6.22185513193483E-12;
      As[351] = -1.0924905566616042E-17;
      As[352] = 6.221855143094405E-12;
      As[353] = -1.09249052909838E-17;
      As[354] = -1.5717956300703978E-7;
      As[355] = -3.5291043212727584E-7;
      As[356] = -1.5717956300678263E-7;
      As[357] = -3.5291043212661547E-7;
      As[358] = 3.1435912601382233E-7;
      As[359] = 7.0582086425379528E-7;
      As[360] = 0.0;
      As[361] = 0.0;
      As[362] = -3.5137704233377777E-20;
      As[363] = -1.7568852132517454E-20;
      As[364] = -1.6527823456548841E-14;
      As[365] = 1.6527823456660875E-14;
      As[366] = 1.4951299471843863E-23;
      As[367] = -0.13380689046808503;
      As[368] = 0.73061941840154687;
      As[369] = 0.13380689046808503;
      As[370] = -0.13281518456307478;
      As[371] = 0.13281518456307478;
      As[372] = -2.9592866890638963E-20;
      As[373] = -0.00093646065973929517;
      As[374] = 1.5926225643696599E-16;
      As[375] = 0.00093646065973930341;
      As[376] = 0.0012149120489449048;
      As[377] = 0.00046011973038863596;
      As[378] = 0.00046011973038863579;
      As[379] = -0.00013906300953769448;
      As[380] = -1.0731740151889774E-5;
      As[381] = -1.0731740151889767E-5;
      As[382] = -0.00016052648984148392;
      As[383] = -2.31228424361381E-8;
      As[384] = 8.3940144125038319E-10;
      As[385] = 8.3940144123779154E-10;
      As[386] = -2.1444039562939883E-8;
      As[387] = 0.0;
      As[388] = 0.0;
      As[389] = 0.0;
      As[390] = 0.0;
      As[391] = 0.0;
      As[392] = 0.0;
      As[393] = 1.6510204962404134E-8;
      As[394] = -2.8990138276410012E-14;
      As[395] = -1.6510204962244611E-8;
      As[396] = 2.8990136116551566E-14;
      As[397] = 4.8895236044954785E-17;
      As[398] = 2.62223055002803E-23;
      As[399] = 0.00041708891965579089;
      As[400] = 0.0009364769061275352;
      As[401] = -4.4817898235658935E-18;
      As[402] = -1.592760450571751E-16;
      As[403] = -0.00041708891965578688;
      As[404] = -0.00093647690612754311;
      As[405] = 0.0;
      As[406] = 0.0;
      As[407] = -1.9862415676150243E-23;
      As[408] = -6.2104431501447947E-33;
      As[409] = -6.2284945670423571E-18;
      As[410] = 1.2456989125164839E-17;
      As[411] = -6.2284945611953416E-18;
      As[412] = -5.0424999514343813E-5;
      As[413] = 0.00060213100710638243;
      As[414] = 0.99934715189965651;
      As[415] = 1.6577859696242124E-19;
      As[416] = -0.00060213100710638265;
      As[417] = -5.0424999514343813E-5;
      As[418] = -7.0580861938218562E-7;
      As[419] = 3.5290430969132834E-7;
      As[420] = 3.5290430969119255E-7;
      As[421] = 2.8444276795826706E-7;
      As[422] = -2.8444276795819279E-7;
      As[423] = 9.6948404539946087E-20;
      As[424] = -4.8361516913044208E-8;
      As[425] = 4.8361516913059461E-8;
      As[426] = 2.5522069623122283E-21;
      As[427] = -4.0101842223397682E-20;
      As[428] = -9.03014884806517E-12;
      As[429] = 9.0301488392352335E-12;
      As[430] = -6.1111078601760614E-21;
      As[431] = -3.4718867520948457E-23;
      As[432] = 0.0;
      As[433] = 0.0;
      As[434] = 0.0;
      As[435] = 0.0;
      As[436] = 0.0;
      As[437] = 0.0;
      As[438] = 6.2218550864824021E-12;
      As[439] = -1.0924906022580809E-17;
      As[440] = -1.2443710222144362E-11;
      As[441] = 2.1849810345833185E-17;
      As[442] = 6.2218551255562833E-12;
      As[443] = -1.092490539905598E-17;
      As[444] = 3.1435912601380941E-7;
      As[445] = 7.0582086425412361E-7;
      As[446] = -1.5717956300693785E-7;
      As[447] = -3.5291043212729739E-7;
      As[448] = -1.5717956300687231E-7;
      As[449] = -3.5291043212716154E-7;
      As[450] = 0.0;
      As[451] = 0.0;
      As[452] = -1.7568852081803993E-20;
      As[453] = -3.513770416512947E-20;
      As[454] = -1.6527823425695406E-14;
      As[455] = -1.9008537174933664E-24;
      As[456] = 1.6527823423342928E-14;
      As[457] = -0.133806890468085;
      As[458] = -0.13281518456307503;
      As[459] = 2.9230180799207658E-18;
      As[460] = 0.73061941840154687;
      As[461] = -0.13281518456307503;
      As[462] = 0.133806890468085;
      As[463] = 5.6060833260794421E-17;
      As[464] = -0.000936460659738933;
      As[465] = 0.00093646065973842618;
      As[466] = -0.00046011973038846568;
      As[467] = -0.00046011973038844568;
      As[468] = -0.0012149120489442255;
      As[469] = 1.073174015188282E-5;
      As[470] = 1.0731740151880341E-5;
      As[471] = 0.00013906300953769898;
      As[472] = 0.00016052648984145025;
      As[473] = -8.3940144777186191E-10;
      As[474] = -8.3940144856234486E-10;
      As[475] = 2.3122842491194982E-8;
      As[476] = 2.1444039563004522E-8;
      As[477] = 0.0;
      As[478] = 0.0;
      As[479] = 0.0;
      As[480] = 0.0;
      As[481] = 0.0;
      As[482] = 0.0;
      As[483] = 1.6510205425643505E-8;
      As[484] = -2.8990136586749389E-14;
      As[485] = -1.558842527683578E-18;
      As[486] = -3.3405607095585638E-24;
      As[487] = -1.6510205438058368E-8;
      As[488] = 2.8990135317989395E-14;
      As[489] = -2.30935161350465E-17;
      As[490] = -5.60651008568952E-17;
      As[491] = 0.000417088919655849;
      As[492] = 0.00093647690612717275;
      As[493] = -0.0004170889196558251;
      As[494] = -0.00093647690612666578;
      As[495] = 0.0;
      As[496] = 0.0;
      As[497] = 1.7568852081811949E-20;
      As[498] = -1.7568852136931327E-20;
      As[499] = -1.8494573329666993E-24;
      As[500] = -1.6527823425660605E-14;
      As[501] = 1.6527823473790154E-14;
      As[502] = 2.5671771274414686E-20;
      As[503] = 0.132815184563075;
      As[504] = -0.133806890468085;
      As[505] = -0.132815184563075;
      As[506] = 0.73061941840154687;
      As[507] = 0.133806890468085;
      As[508] = 0.00093646065973842314;
      As[509] = -0.00093646065973895353;
      As[510] = 6.8581121820192994E-17;
      As[511] = 0.00046011973038847669;
      As[512] = 0.001214912048944279;
      As[513] = 0.00046011973038845663;
      As[514] = -1.0731740151884553E-5;
      As[515] = -0.00013906300953771242;
      As[516] = -1.0731740151883811E-5;
      As[517] = -0.00016052648984145025;
      As[518] = 8.3940144776508151E-10;
      As[519] = -2.3122842495660388E-8;
      As[520] = 8.394014468117421E-10;
      As[521] = -2.1444039562827373E-8;
      As[522] = 0.0;
      As[523] = 0.0;
      As[524] = 0.0;
      As[525] = 0.0;
      As[526] = 0.0;
      As[527] = 0.0;
      As[528] = -1.1635573955627294E-17;
      As[529] = -2.944020688902497E-24;
      As[530] = 1.651020545208791E-8;
      As[531] = -2.8990135297034617E-14;
      As[532] = -1.6510205289277242E-8;
      As[533] = 2.8990135800457526E-14;
      As[534] = -0.00041708891965589373;
      As[535] = -0.00093647690612666307;
      As[536] = 0.00041708891965590837;
      As[537] = 0.00093647690612719346;
      As[538] = -1.3580362725654386E-17;
      As[539] = -6.8588393466863081E-17;
      As[540] = 0.0;
      As[541] = 0.0;
      As[542] = 3.8635508707230738E-33;
      As[543] = -1.9862415639913302E-23;
      As[544] = -6.2284945484580884E-18;
      As[545] = -6.2284945520964844E-18;
      As[546] = 1.2456989109654892E-17;
      As[547] = -5.0424999514343827E-5;
      As[548] = -2.1290515325334625E-21;
      As[549] = -5.0424999514343827E-5;
      As[550] = 0.00060213100710638276;
      As[551] = 0.00060213100710638276;
      As[552] = 0.99934715189965651;
      As[553] = 3.5290430969088396E-7;
      As[554] = -7.0580861938203749E-7;
      As[555] = 3.5290430969087195E-7;
      As[556] = 3.146433306785459E-21;
      As[557] = 2.8444276795810496E-7;
      As[558] = -2.8444276795813228E-7;
      As[559] = -1.2268184430162786E-22;
      As[560] = -4.83615169130222E-8;
      As[561] = 4.8361516913023912E-8;
      As[562] = -2.9223126434903779E-22;
      As[563] = 8.7022111770949277E-23;
      As[564] = -9.0301488284850214E-12;
      As[565] = 9.0301488273336538E-12;
      As[566] = -8.511098673834391E-24;
      As[567] = 0.0;
      As[568] = 0.0;
      As[569] = 0.0;
      As[570] = 0.0;
      As[571] = 0.0;
      As[572] = 0.0;
      As[573] = 6.221855060679894E-12;
      As[574] = -1.0924905278444844E-17;
      As[575] = 6.2218550565758785E-12;
      As[576] = -1.0924905667766797E-17;
      As[577] = -1.2443710087543905E-11;
      As[578] = 2.1849810821755909E-17;
      As[579] = -1.5717956300672225E-7;
      As[580] = -3.5291043212685296E-7;
      As[581] = 3.1435912601345594E-7;
      As[582] = 7.0582086425397538E-7;
      As[583] = -1.5717956300673329E-7;
      As[584] = -3.5291043212684084E-7;
      As[585] = 0.0;
      As[586] = 0.0;
      As[587] = -1.2666600800617312E-20;
      As[588] = -4.8100772027957E-33;
      As[589] = -3.9720170757817437E-15;
      As[590] = 7.94403415155551E-15;
      As[591] = -3.9720170757746022E-15;
      As[592] = 1.2144066877969172E-5;
      As[593] = -0.00014501376876600237;
      As[594] = -2.4288133755938348E-5;
      As[595] = -6.85802507015598E-20;
      As[596] = 0.00014501376876600256;
      As[597] = 1.2144066877969175E-5;
      As[598] = 0.32189695418714037;
      As[599] = 0.04739923107025492;
      As[600] = 0.047399231070254892;
      As[601] = 0.00018139399813298112;
      As[602] = -0.00018139399813298112;
      As[603] = -4.9272472251541394E-20;
      As[604] = -3.084096309287552E-5;
      As[605] = 3.0840963092875533E-5;
      As[606] = 4.8319945176671063E-21;
      As[607] = 2.3633674021777163E-21;
      As[608] = -5.7586797310493709E-9;
      As[609] = 5.7586797310658533E-9;
      As[610] = 1.5893587927836035E-21;
      As[611] = -7.1324547771099508E-22;
      As[612] = 0.0;
      As[613] = 0.0;
      As[614] = 0.0;
      As[615] = 0.0;
      As[616] = 0.0;
      As[617] = 0.0;
      As[618] = 3.9677829791557755E-9;
      As[619] = -6.9669973717165656E-15;
      As[620] = -7.9355659583289086E-9;
      As[621] = 1.3933996075423774E-14;
      As[622] = 3.9677829791358388E-9;
      As[623] = -6.9669992329050552E-15;
      As[624] = -0.24431416306927842;
      As[625] = 0.67807746348067477;
      As[626] = 0.016923316451174508;
      As[627] = -0.04740270250583644;
      As[628] = 0.016923316451174512;
      As[629] = -0.047402702505836412;
      As[630] = 0.0;
      As[631] = 0.0;
      As[632] = 1.8016746823245469E-32;
      As[633] = -1.2666600800619978E-20;
      As[634] = -3.972017075778059E-15;
      As[635] = -3.9720170757952861E-15;
      As[636] = 7.9440341515586944E-15;
      As[637] = 1.2144066877969184E-5;
      As[638] = -1.978034288614478E-20;
      As[639] = 1.2144066877969182E-5;
      As[640] = -0.00014501376876600256;
      As[641] = -0.00014501376876600259;
      As[642] = -2.4288133755938364E-5;
      As[643] = 0.047399231070254913;
      As[644] = 0.3218969541871407;
      As[645] = 0.047399231070254878;
      As[646] = -4.9901325735620512E-20;
      As[647] = 0.00018139399813298139;
      As[648] = -0.00018139399813298137;
      As[649] = -3.4998314337855439E-21;
      As[650] = -3.0840963092875554E-5;
      As[651] = 3.0840963092875554E-5;
      As[652] = -1.6580295738562972E-22;
      As[653] = -3.3772517730359226E-23;
      As[654] = -5.758679731042143E-9;
      As[655] = 5.7586797310527615E-9;
      As[656] = -1.4265758695885622E-20;
      As[657] = 0.0;
      As[658] = 0.0;
      As[659] = 0.0;
      As[660] = 0.0;
      As[661] = 0.0;
      As[662] = 0.0;
      As[663] = 3.9677829791154149E-9;
      As[664] = -6.9669983391863453E-15;
      As[665] = 3.9677829790633538E-9;
      As[666] = -6.9669973121163571E-15;
      As[667] = -7.9355659582269817E-9;
      As[668] = 1.3933995645797841E-14;
      As[669] = 0.016923316451174529;
      As[670] = -0.047402702505836405;
      As[671] = -0.24431416306927858;
      As[672] = 0.67807746348067466;
      As[673] = 0.016923316451174536;
      As[674] = -0.047402702505836392;
      As[675] = 0.0;
      As[676] = 0.0;
      As[677] = 1.2666600800562406E-20;
      As[678] = 1.2666600800575483E-20;
      As[679] = 7.94403415153449E-15;
      As[680] = -3.9720170757503092E-15;
      As[681] = -3.9720170757561294E-15;
      As[682] = -2.4288133755938341E-5;
      As[683] = 0.00014501376876600237;
      As[684] = 1.2144066877969174E-5;
      As[685] = 0.00014501376876600251;
      As[686] = -1.0661043538301263E-19;
      As[687] = 1.2144066877969169E-5;
      As[688] = 0.047399231070254878;
      As[689] = 0.047399231070254913;
      As[690] = 0.32189695418714021;
      As[691] = -0.00018139399813298112;
      As[692] = -3.3305481533983743E-20;
      As[693] = 0.00018139399813298112;
      As[694] = 3.084096309287552E-5;
      As[695] = -9.9705441716164454E-21;
      As[696] = -3.0840963092875527E-5;
      As[697] = -6.9424232659766407E-22;
      As[698] = 5.7586797310459836E-9;
      As[699] = -1.7646935508805398E-20;
      As[700] = -5.7586797310658558E-9;
      As[701] = 1.2477223357262405E-20;
      As[702] = 0.0;
      As[703] = 0.0;
      As[704] = 0.0;
      As[705] = 0.0;
      As[706] = 0.0;
      As[707] = 0.0;
      As[708] = -7.9355659583378156E-9;
      As[709] = 1.3933996684934597E-14;
      As[710] = 3.96778297920102E-9;
      As[711] = -6.966999392106111E-15;
      As[712] = 3.9677829791357726E-9;
      As[713] = -6.9669987109431053E-15;
      As[714] = 0.016923316451174519;
      As[715] = -0.0474027025058364;
      As[716] = 0.016923316451174512;
      As[717] = -0.047402702505836447;
      As[718] = -0.24431416306927864;
      As[719] = 0.6780774634806751;
      As[720] = 0.0;
      As[721] = 0.0;
      As[722] = -9.3060665839515753E-17;
      As[723] = -4.6530332919788944E-17;
      As[724] = -4.3773214252566949E-11;
      As[725] = 4.37732142525671E-11;
      As[726] = 1.4653556203484588E-27;
      As[727] = -2.4956449207745748E-5;
      As[728] = 0.0002875502345273766;
      As[729] = 2.4956449207745748E-5;
      As[730] = 1.0457738501649563E-5;
      As[731] = -1.0457738501649563E-5;
      As[732] = 1.0200829413479651E-23;
      As[733] = 0.00046249135908813335;
      As[734] = 3.9312090717181776E-20;
      As[735] = -0.00046249135908813324;
      As[736] = 0.99909533287543206;
      As[737] = 6.1040493303412454E-5;
      As[738] = 6.1040493303412454E-5;
      As[739] = 4.8938570736512966E-5;
      As[740] = -1.4440603363854119E-5;
      As[741] = -1.4440603363854119E-5;
      As[742] = 2.0057364008804732E-5;
      As[743] = -4.9001569450193306E-5;
      As[744] = 1.4461380700905895E-5;
      As[745] = 1.4461380700905895E-5;
      As[746] = -2.0078808048342793E-5;
      As[747] = 0.0;
      As[748] = 0.0;
      As[749] = 0.0;
      As[750] = 0.0;
      As[751] = 0.0;
      As[752] = 0.0;
      As[753] = 4.3726552816725049E-5;
      As[754] = -7.6779103877763351E-11;
      As[755] = -4.3726552816725029E-5;
      As[756] = 7.6779103870817412E-11;
      As[757] = 4.84452509869281E-21;
      As[758] = -5.2420093598159333E-27;
      As[759] = -0.00020598838755913067;
      As[760] = -0.0004624993827184901;
      As[761] = 4.0595442053165938E-21;
      As[762] = -4.21810834222611E-20;
      As[763] = 0.00020598838755913072;
      As[764] = 0.00046249938271848977;
      As[765] = 0.0;
      As[766] = 0.0;
      As[767] = 4.6530332919549771E-17;
      As[768] = -4.6530332919745705E-17;
      As[769] = -3.3117646094463866E-23;
      As[770] = -4.3773214252378636E-11;
      As[771] = 4.377321425247013E-11;
      As[772] = 3.6213975066207138E-21;
      As[773] = -1.0457738501649715E-5;
      As[774] = -2.4956449207745765E-5;
      As[775] = 1.0457738501649675E-5;
      As[776] = 0.00028755023452737665;
      As[777] = 2.4956449207745762E-5;
      As[778] = -0.00046249135908813351;
      As[779] = 0.00046249135908813362;
      As[780] = 6.41798540259714E-20;
      As[781] = 6.1040493303497753E-5;
      As[782] = 0.999095332875433;
      As[783] = 6.1040493303719757E-5;
      As[784] = -1.4440603363854133E-5;
      As[785] = 4.893857073651298E-5;
      As[786] = -1.4440603363854122E-5;
      As[787] = 2.0057364008804718E-5;
      As[788] = 1.4461380700892026E-5;
      As[789] = -4.90015694501492E-5;
      As[790] = 1.4461380700910891E-5;
      As[791] = -2.0078808048292581E-5;
      As[792] = 0.0;
      As[793] = 0.0;
      As[794] = 0.0;
      As[795] = 0.0;
      As[796] = 0.0;
      As[797] = 0.0;
      As[798] = 2.8683311350576514E-18;
      As[799] = -2.499892107465431E-22;
      As[800] = 4.3726552816723863E-5;
      As[801] = -7.67791080881712E-11;
      As[802] = -4.3726552816726675E-5;
      As[803] = 7.67791051704855E-11;
      As[804] = 0.00020598838755913075;
      As[805] = 0.00046249938271849026;
      As[806] = -0.00020598838755913078;
      As[807] = -0.00046249938271849026;
      As[808] = -1.8648624479819824E-20;
      As[809] = -5.9255351479594968E-20;
      As[810] = 0.0;
      As[811] = 0.0;
      As[812] = 4.6530332919679952E-17;
      As[813] = 9.3060665839460755E-17;
      As[814] = 4.3773214252524788E-11;
      As[815] = 1.3966882109412924E-27;
      As[816] = -4.3773214252470181E-11;
      As[817] = 2.4956449207745762E-5;
      As[818] = -1.0457738501649765E-5;
      As[819] = -1.938815336602461E-22;
      As[820] = -0.00028755023452737655;
      As[821] = -1.0457738501649771E-5;
      As[822] = -2.4956449207745762E-5;
      As[823] = -1.7976904006754611E-21;
      As[824] = -0.00046249135908813373;
      As[825] = 0.00046249135908813351;
      As[826] = 6.10404933031262E-5;
      As[827] = 6.1040493303126211E-5;
      As[828] = 0.99909533287543251;
      As[829] = -1.4440603363854129E-5;
      As[830] = -1.4440603363854131E-5;
      As[831] = 4.8938570736512966E-5;
      As[832] = 2.0057364008804694E-5;
      As[833] = 1.4461380700892019E-5;
      As[834] = 1.4461380700892021E-5;
      As[835] = -4.9001569450171249E-5;
      As[836] = -2.0078808048367685E-5;
      As[837] = 0.0;
      As[838] = 0.0;
      As[839] = 0.0;
      As[840] = 0.0;
      As[841] = 0.0;
      As[842] = 0.0;
      As[843] = -4.3726552816745358E-5;
      As[844] = 7.6779102437576338E-11;
      As[845] = -3.0664617928207952E-21;
      As[846] = 6.0701308275715029E-23;
      As[847] = 4.37265528167454E-5;
      As[848] = -7.677910516748267E-11;
      As[849] = -6.3792916522708422E-21;
      As[850] = -7.5741298043947327E-22;
      As[851] = 0.0002059883875591308;
      As[852] = 0.00046249938271849042;
      As[853] = -0.0002059883875591308;
      As[854] = -0.0004624993827184901;
      As[855] = 0.0;
      As[856] = 0.0;
      As[857] = 5.4667094528802031E-16;
      As[858] = 2.733354726440707E-16;
      As[859] = 2.5713919192220604E-10;
      As[860] = -2.5713919192292476E-10;
      As[861] = -3.213468560966218E-22;
      As[862] = -0.78617878122916129;
      As[863] = 10.172917941560177;
      As[864] = 0.78617878122916129;
      As[865] = -0.78506219805093458;
      As[866] = 0.78506219805093458;
      As[867] = -3.0503349911308396E-18;
      As[868] = 14.56941610523968;
      As[869] = 1.4318669302719348E-15;
      As[870] = -14.569416105239682;
      As[871] = -11.745187981837121;
      As[872] = -0.0021603941695526939;
      As[873] = -0.0021603941695553077;
      As[874] = 0.997360899195679;
      As[875] = 0.00078833253132878749;
      As[876] = 0.00078833253132878759;
      As[877] = -0.0010624357416640687;
      As[878] = 0.00028787047785449726;
      As[879] = -8.4933118991603322E-5;
      As[880] = -8.4933118991825434E-5;
      As[881] = 0.00011800423987103967;
      As[882] = 0.0;
      As[883] = 0.0;
      As[884] = 0.0;
      As[885] = 0.0;
      As[886] = 0.0;
      As[887] = 0.0;
      As[888] = -0.00025686508630557944;
      As[889] = 4.5102737134512163E-10;
      As[890] = 0.00025686508630459276;
      As[891] = -4.510273305018592E-10;
      As[892] = 2.5736324790007915E-16;
      As[893] = -5.2863226715149694E-21;
      As[894] = -6.4890521135648935;
      As[895] = -14.569668865874116;
      As[896] = -4.5637908167696623E-17;
      As[897] = -1.5746063431082566E-15;
      As[898] = 6.4890521135648962;
      As[899] = 14.569668865874108;
      As[900] = 0.0;
      As[901] = 0.0;
      As[902] = -2.7333547264452853E-16;
      As[903] = 2.733354726442372E-16;
      As[904] = 1.1474220818077869E-22;
      As[905] = 2.5713919192325238E-10;
      As[906] = -2.5713919192227904E-10;
      As[907] = 1.6001563395473364E-18;
      As[908] = 0.78506219805093447;
      As[909] = -0.78617878122916152;
      As[910] = -0.78506219805093458;
      As[911] = 10.172917941560176;
      As[912] = 0.78617878122916152;
      As[913] = -14.569416105239688;
      As[914] = 14.569416105239695;
      As[915] = 2.3791479309636305E-16;
      As[916] = -0.0021603941695594047;
      As[917] = -11.745187981837127;
      As[918] = -0.002160394169559409;
      As[919] = 0.00078833253132863494;
      As[920] = 0.99736089919567883;
      As[921] = 0.000788332531328643;
      As[922] = -0.0010624357416637391;
      As[923] = -8.4933118991762442E-5;
      As[924] = 0.00028787047785418165;
      As[925] = -8.4933118991755029E-5;
      As[926] = 0.00011800423987002866;
      As[927] = 0.0;
      As[928] = 0.0;
      As[929] = 0.0;
      As[930] = 0.0;
      As[931] = 0.0;
      As[932] = 0.0;
      As[933] = 4.184435731912125E-16;
      As[934] = 6.4624634419530563E-22;
      As[935] = -0.00025686508630367249;
      As[936] = 4.51027314122938E-10;
      As[937] = 0.0002568650863042527;
      As[938] = -4.5102735323396489E-10;
      As[939] = 6.4890521135648935;
      As[940] = 14.569668865874119;
      As[941] = -6.4890521135648962;
      As[942] = -14.569668865874123;
      As[943] = 2.2237004533864521E-16;
      As[944] = -2.2391257503834919E-16;
      As[945] = 0.0;
      As[946] = 0.0;
      As[947] = -2.7333547264349355E-16;
      As[948] = -5.4667094528750646E-16;
      As[949] = -2.5713919192187476E-10;
      As[950] = 1.3928558523922156E-23;
      As[951] = 2.5713919192214116E-10;
      As[952] = 0.78617878122916163;
      As[953] = 0.78506219805093447;
      As[954] = -1.6001563395473364E-18;
      As[955] = -10.172917941560176;
      As[956] = 0.78506219805093458;
      As[957] = -0.78617878122916152;
      As[958] = -5.72164307593589E-16;
      As[959] = -14.569416105239691;
      As[960] = 14.569416105239686;
      As[961] = -0.0021603941695561724;
      As[962] = -0.0021603941695562687;
      As[963] = -11.745187981837125;
      As[964] = 0.000788332531328635;
      As[965] = 0.000788332531328643;
      As[966] = 0.99736089919567883;
      As[967] = -0.0010624357416637393;
      As[968] = -8.4933118991754243E-5;
      As[969] = -8.4933118991754622E-5;
      As[970] = 0.00028787047785422848;
      As[971] = 0.00011800423987066232;
      As[972] = 0.0;
      As[973] = 0.0;
      As[974] = 0.0;
      As[975] = 0.0;
      As[976] = 0.0;
      As[977] = 0.0;
      As[978] = 0.00025686508630452272;
      As[979] = -4.5102737345620323E-10;
      As[980] = 4.7962599303692577E-17;
      As[981] = 1.7862933179178464E-21;
      As[982] = -0.00025686508630452505;
      As[983] = 4.5102736013391379E-10;
      As[984] = 8.7251415816606426E-18;
      As[985] = 4.5070040695794807E-16;
      As[986] = 6.4890521135648953;
      As[987] = 14.569668865874121;
      As[988] = -6.4890521135648962;
      As[989] = -14.569668865874114;
      As[990] = 0.0;
      As[991] = 0.0;
      As[992] = 1.1772337165010519E-28;
      As[993] = -2.5186682538247821E-29;
      As[994] = 5.65190388504358E-23;
      As[995] = -5.42287260250875E-23;
      As[996] = 8.021330411512742E-23;
      As[997] = 3.0598552448928167E-18;
      As[998] = 4.9374124489386046;
      As[999] = -3.0598552448928167E-18;
      As[1000] = -4.9374124489386046;
      As[1001] = 4.9374124489386046;
      As[1002] = 3.0598552448928167E-18;
      As[1003] = -1.4607439142113302E-15;
      As[1004] = -3.8816092975772889E-18;
      As[1005] = 1.3541817113055472E-15;
      As[1006] = -4.940131288186369;
      As[1007] = -4.94013128818637;
      As[1008] = -4.94013128818637;
      As[1009] = -0.00044670565823168876;
      As[1010] = -0.000446705658231675;
      As[1011] = -0.00044670565823166805;
      As[1012] = 0.998659883025305;
      As[1013] = 4.9615388092228606E-5;
      As[1014] = 4.9615388092103428E-5;
      As[1015] = 4.9615388092110496E-5;
      As[1016] = 0.00014884616427595412;
      As[1017] = 0.0;
      As[1018] = 0.0;
      As[1019] = 0.0;
      As[1020] = 0.0;
      As[1021] = 0.0;
      As[1022] = 0.0;
      As[1023] = -5.1028847454498812E-16;
      As[1024] = 9.913554065778106E-23;
      As[1025] = 5.5944946406261478E-16;
      As[1026] = -9.5118284304130369E-23;
      As[1027] = 2.2029449727382057E-16;
      As[1028] = 1.4069576532927425E-22;
      As[1029] = 1.4189414231852873E-16;
      As[1030] = 1.4607694627619426E-15;
      As[1031] = 3.9535935670072396E-17;
      As[1032] = 3.88188316059688E-18;
      As[1033] = -1.8116617317915591E-16;
      As[1034] = -1.3542049980947888E-15;
      As[1035] = 0.0;
      As[1036] = 0.0;
      As[1037] = 2.9315990511099932E-12;
      As[1038] = 1.4657995255544194E-12;
      As[1039] = 1.3789447153539517E-6;
      As[1040] = -1.3789447153552737E-6;
      As[1041] = -2.06538058424869E-22;
      As[1042] = -0.00014679670666609812;
      As[1043] = 0.0016915122121546554;
      As[1044] = 0.00014679670666609812;
      As[1045] = 6.1404985189379577E-5;
      As[1046] = -6.14049851893793E-5;
      As[1047] = 2.2800898501401506E-20;
      As[1048] = 0.0027204274057783043;
      As[1049] = 5.8898329274624671E-19;
      As[1050] = -0.0027204274057783047;
      As[1051] = 11.759723390130716;
      As[1052] = 0.0011735983914471894;
      As[1053] = 0.0011735983914454285;
      As[1054] = 0.00028787047785419487;
      As[1055] = -8.4933118991766156E-5;
      As[1056] = -8.493311899176621E-5;
      As[1057] = 0.00011800423987066248;
      As[1058] = 0.99943474949493882;
      As[1059] = 0.00022356005133550579;
      As[1060] = 0.00022356005133535611;
      As[1061] = -0.00011813040239100823;
      As[1062] = 0.0;
      As[1063] = 0.0;
      As[1064] = 0.0;
      As[1065] = 0.0;
      As[1066] = 0.0;
      As[1067] = 0.0;
      As[1068] = -1.3774747858202308;
      As[1069] = 2.41869689799831E-6;
      As[1070] = 1.3774747858202305;
      As[1071] = -2.4186968319045603E-6;
      As[1072] = 2.0637273318222924E-16;
      As[1073] = -3.6293251722657294E-22;
      As[1074] = -0.0012116474043813539;
      As[1075] = -0.0027204746016955628;
      As[1076] = -1.7327726781057152E-19;
      As[1077] = -5.68557507328622E-19;
      As[1078] = 0.0012116474043813543;
      As[1079] = 0.0027204746016955619;
      As[1080] = 0.0;
      As[1081] = 0.0;
      As[1082] = -1.4657995255551048E-12;
      As[1083] = 1.4657995255558391E-12;
      As[1084] = -3.6193168141090024E-22;
      As[1085] = 1.378944715354833E-6;
      As[1086] = -1.3789447153552716E-6;
      As[1087] = 5.2968634989258707E-21;
      As[1088] = -6.1404985189379414E-5;
      As[1089] = -0.00014679670666609818;
      As[1090] = 6.14049851893793E-5;
      As[1091] = 0.0016915122121546559;
      As[1092] = 0.00014679670666609815;
      As[1093] = -0.0027204274057783051;
      As[1094] = 0.0027204274057783056;
      As[1095] = 3.1883395495577606E-19;
      As[1096] = 0.001173598391450815;
      As[1097] = 11.759723390130723;
      As[1098] = 0.0011735983914508146;
      As[1099] = -8.49331189917662E-5;
      As[1100] = 0.00028787047785419487;
      As[1101] = -8.4933118991766183E-5;
      As[1102] = 0.00011800423987066258;
      As[1103] = 0.00022356005133528317;
      As[1104] = 0.99943474949493882;
      As[1105] = 0.00022356005133528312;
      As[1106] = -0.00011813040238974583;
      As[1107] = 0.0;
      As[1108] = 0.0;
      As[1109] = 0.0;
      As[1110] = 0.0;
      As[1111] = 0.0;
      As[1112] = 0.0;
      As[1113] = -1.9261403614703751E-16;
      As[1114] = 1.5148099615500428E-17;
      As[1115] = -1.3774747858202299;
      As[1116] = 2.4186968539693942E-6;
      As[1117] = 1.3774747858202308;
      As[1118] = -2.4186968320289086E-6;
      As[1119] = 0.0012116474043813543;
      As[1120] = 0.0027204746016955632;
      As[1121] = -0.0012116474043813545;
      As[1122] = -0.0027204746016955637;
      As[1123] = -5.499746633445371E-20;
      As[1124] = -3.218451355268742E-19;
      As[1125] = 0.0;
      As[1126] = 0.0;
      As[1127] = -1.4657995255535389E-12;
      As[1128] = -2.9315990511079159E-12;
      As[1129] = -1.3789447153539513E-6;
      As[1130] = 1.1420865825448816E-22;
      As[1131] = 1.3789447153532936E-6;
      As[1132] = 0.00014679670666609815;
      As[1133] = -6.1404985189380051E-5;
      As[1134] = -1.043961099417464E-21;
      As[1135] = -0.0016915122121546552;
      As[1136] = -6.1404985189380051E-5;
      As[1137] = -0.00014679670666609815;
      As[1138] = 1.5075398066653314E-20;
      As[1139] = -0.0027204274057783051;
      As[1140] = 0.0027204274057783051;
      As[1141] = 0.0011735983914461633;
      As[1142] = 0.0011735983914461636;
      As[1143] = 11.759723390130716;
      As[1144] = -8.4933118991766237E-5;
      As[1145] = -8.4933118991766237E-5;
      As[1146] = 0.00028787047785419476;
      As[1147] = 0.00011800423987066241;
      As[1148] = 0.00022356005133528317;
      As[1149] = 0.0002235600513352832;
      As[1150] = 0.99943474949493871;
      As[1151] = -0.00011813040239068626;
      As[1152] = 0.0;
      As[1153] = 0.0;
      As[1154] = 0.0;
      As[1155] = 0.0;
      As[1156] = 0.0;
      As[1157] = 0.0;
      As[1158] = 1.3774747858202292;
      As[1159] = -2.4186968980502136E-6;
      As[1160] = -8.2449113539926813E-20;
      As[1161] = -5.7092632739833575E-18;
      As[1162] = -1.3774747858202288;
      As[1163] = 2.4186969309401252E-6;
      As[1164] = -1.942552056579236E-20;
      As[1165] = -3.7292990748650848E-20;
      As[1166] = 0.0012116474043813543;
      As[1167] = 0.0027204746016955628;
      As[1168] = -0.0012116474043813547;
      As[1169] = -0.0027204746016955624;
      As[1170] = 0.0;
      As[1171] = 0.0;
      As[1172] = -2.736801705705222E-28;
      As[1173] = -1.3681846547009532E-28;
      As[1174] = -1.2871139095107806E-22;
      As[1175] = 1.2875207409043582E-22;
      As[1176] = 2.8199112080141376E-30;
      As[1177] = 6.6052364936919682E-21;
      As[1178] = 0.00065956757665613952;
      As[1179] = -6.6052364936919682E-21;
      As[1180] = -0.00065956757665613962;
      As[1181] = 0.00065956757665613962;
      As[1182] = -1.7102708434243452E-22;
      As[1183] = -3.9570594845647374E-19;
      As[1184] = 1.7575373396718565E-19;
      As[1185] = 5.0783740297279446E-19;
      As[1186] = 4.9454129578386583;
      As[1187] = 4.9454129578386592;
      As[1188] = 4.9454129578386592;
      As[1189] = 4.9615388092089971E-5;
      As[1190] = 4.9615388092089991E-5;
      As[1191] = 4.9615388092089991E-5;
      As[1192] = 0.00014884616427626995;
      As[1193] = -4.9668433664074092E-5;
      As[1194] = -4.966843366397852E-5;
      As[1195] = -4.9668433663970815E-5;
      As[1196] = 0.99985099469900818;
      As[1197] = 0.0;
      As[1198] = 0.0;
      As[1199] = 0.0;
      As[1200] = 0.0;
      As[1201] = 0.0;
      As[1202] = 0.0;
      As[1203] = 1.8743488665561843E-16;
      As[1204] = -2.2576238802575796E-22;
      As[1205] = -1.41933321060806E-16;
      As[1206] = 2.258337276383148E-22;
      As[1207] = -4.536788666701399E-17;
      As[1208] = 4.94404588853109E-30;
      As[1209] = 1.2974696882884903E-19;
      As[1210] = 3.9571227512096741E-19;
      As[1211] = -1.5788999868343103E-20;
      As[1212] = -1.7575732139333486E-19;
      As[1213] = -8.643720568206758E-20;
      As[1214] = -5.0784675162357012E-19;
      As[1215] = 0.0;
      As[1216] = 1.3135322195321801E-5;
      As[1217] = 0.0;
      As[1218] = 0.0;
      As[1219] = 0.0;
      As[1220] = 0.0;
      As[1221] = 0.0;
      As[1222] = 0.0;
      As[1223] = 0.0;
      As[1224] = 0.0;
      As[1225] = 0.0;
      As[1226] = 0.0;
      As[1227] = 0.0;
      As[1228] = 0.0;
      As[1229] = 0.0;
      As[1230] = 0.0;
      As[1231] = 0.0;
      As[1232] = 0.0;
      As[1233] = 0.0;
      As[1234] = 0.0;
      As[1235] = 0.0;
      As[1236] = 0.0;
      As[1237] = 0.0;
      As[1238] = 0.0;
      As[1239] = 0.0;
      As[1240] = 0.0;
      As[1241] = 0.0;
      As[1242] = 0.9933993216904905;
      As[1243] = 0.0033003117895002527;
      As[1244] = 0.0033003117895002523;
      As[1245] = -0.0066006783095096519;
      As[1246] = 0.0033003117895002518;
      As[1247] = 0.0033003117895002518;
      As[1248] = 0.0;
      As[1249] = 0.0;
      As[1250] = 0.0;
      As[1251] = 0.0;
      As[1252] = 0.0;
      As[1253] = 0.0;
      As[1254] = 0.0;
      As[1255] = 0.0;
      As[1256] = 0.0;
      As[1257] = 0.0;
      As[1258] = 0.0;
      As[1259] = 0.0;
      As[1260] = 0.0;
      As[1261] = 1.3135322195322191E-5;
      As[1262] = 0.0;
      As[1263] = 0.0;
      As[1264] = 0.0;
      As[1265] = 0.0;
      As[1266] = 0.0;
      As[1267] = 0.0;
      As[1268] = 0.0;
      As[1269] = 0.0;
      As[1270] = 0.0;
      As[1271] = 0.0;
      As[1272] = 0.0;
      As[1273] = 0.0;
      As[1274] = 0.0;
      As[1275] = 0.0;
      As[1276] = 0.0;
      As[1277] = 0.0;
      As[1278] = 0.0;
      As[1279] = 0.0;
      As[1280] = 0.0;
      As[1281] = 0.0;
      As[1282] = 0.0;
      As[1283] = 0.0;
      As[1284] = 0.0;
      As[1285] = 0.0;
      As[1286] = 0.0;
      As[1287] = 0.0033003117895002518;
      As[1288] = 0.99339932169049017;
      As[1289] = 0.0033003117895002518;
      As[1290] = 0.0033003117895002514;
      As[1291] = -0.006600678309509651;
      As[1292] = 0.0033003117895002514;
      As[1293] = 0.0;
      As[1294] = 0.0;
      As[1295] = 0.0;
      As[1296] = 0.0;
      As[1297] = 0.0;
      As[1298] = 0.0;
      As[1299] = 0.0;
      As[1300] = 0.0;
      As[1301] = 0.0;
      As[1302] = 0.0;
      As[1303] = 0.0;
      As[1304] = 0.0;
      As[1305] = 0.0;
      As[1306] = 1.3135322195321135E-5;
      As[1307] = 0.0;
      As[1308] = 0.0;
      As[1309] = 0.0;
      As[1310] = 0.0;
      As[1311] = 0.0;
      As[1312] = 0.0;
      As[1313] = 0.0;
      As[1314] = 0.0;
      As[1315] = 0.0;
      As[1316] = 0.0;
      As[1317] = 0.0;
      As[1318] = 0.0;
      As[1319] = 0.0;
      As[1320] = 0.0;
      As[1321] = 0.0;
      As[1322] = 0.0;
      As[1323] = 0.0;
      As[1324] = 0.0;
      As[1325] = 0.0;
      As[1326] = 0.0;
      As[1327] = 0.0;
      As[1328] = 0.0;
      As[1329] = 0.0;
      As[1330] = 0.0;
      As[1331] = 0.0;
      As[1332] = 0.0033003117895002518;
      As[1333] = 0.0033003117895002518;
      As[1334] = 0.99339932169049028;
      As[1335] = 0.0033003117895002518;
      As[1336] = 0.0033003117895002514;
      As[1337] = -0.006600678309509651;
      As[1338] = 0.0;
      As[1339] = 0.0;
      As[1340] = 0.0;
      As[1341] = 0.0;
      As[1342] = 0.0;
      As[1343] = 0.0;
      As[1344] = 0.0;
      As[1345] = 0.0;
      As[1346] = 0.0;
      As[1347] = 0.0;
      As[1348] = 0.0;
      As[1349] = 0.0;
      As[1350] = 0.0;
      As[1351] = 1.3135322195320705E-5;
      As[1352] = 0.0;
      As[1353] = 0.0;
      As[1354] = 0.0;
      As[1355] = 0.0;
      As[1356] = 0.0;
      As[1357] = 0.0;
      As[1358] = 0.0;
      As[1359] = 0.0;
      As[1360] = 0.0;
      As[1361] = 0.0;
      As[1362] = 0.0;
      As[1363] = 0.0;
      As[1364] = 0.0;
      As[1365] = 0.0;
      As[1366] = 0.0;
      As[1367] = 0.0;
      As[1368] = 0.0;
      As[1369] = 0.0;
      As[1370] = 0.0;
      As[1371] = 0.0;
      As[1372] = 0.0;
      As[1373] = 0.0;
      As[1374] = 0.0;
      As[1375] = 0.0;
      As[1376] = 0.0;
      As[1377] = -0.006600678309509651;
      As[1378] = 0.0033003117895002523;
      As[1379] = 0.0033003117895002518;
      As[1380] = 0.99339932169049028;
      As[1381] = 0.0033003117895002514;
      As[1382] = 0.0033003117895002518;
      As[1383] = 0.0;
      As[1384] = 0.0;
      As[1385] = 0.0;
      As[1386] = 0.0;
      As[1387] = 0.0;
      As[1388] = 0.0;
      As[1389] = 0.0;
      As[1390] = 0.0;
      As[1391] = 0.0;
      As[1392] = 0.0;
      As[1393] = 0.0;
      As[1394] = 0.0;
      As[1395] = 0.0;
      As[1396] = 1.3135322195320707E-5;
      As[1397] = 0.0;
      As[1398] = 0.0;
      As[1399] = 0.0;
      As[1400] = 0.0;
      As[1401] = 0.0;
      As[1402] = 0.0;
      As[1403] = 0.0;
      As[1404] = 0.0;
      As[1405] = 0.0;
      As[1406] = 0.0;
      As[1407] = 0.0;
      As[1408] = 0.0;
      As[1409] = 0.0;
      As[1410] = 0.0;
      As[1411] = 0.0;
      As[1412] = 0.0;
      As[1413] = 0.0;
      As[1414] = 0.0;
      As[1415] = 0.0;
      As[1416] = 0.0;
      As[1417] = 0.0;
      As[1418] = 0.0;
      As[1419] = 0.0;
      As[1420] = 0.0;
      As[1421] = 0.0;
      As[1422] = 0.0033003117895002514;
      As[1423] = -0.00660067830950965;
      As[1424] = 0.0033003117895002514;
      As[1425] = 0.0033003117895002514;
      As[1426] = 0.99339932169049028;
      As[1427] = 0.0033003117895002514;
      As[1428] = 0.0;
      As[1429] = 0.0;
      As[1430] = 0.0;
      As[1431] = 0.0;
      As[1432] = 0.0;
      As[1433] = 0.0;
      As[1434] = 0.0;
      As[1435] = 0.0;
      As[1436] = 0.0;
      As[1437] = 0.0;
      As[1438] = 0.0;
      As[1439] = 0.0;
      As[1440] = 0.0;
      As[1441] = 1.3135322195319651E-5;
      As[1442] = 0.0;
      As[1443] = 0.0;
      As[1444] = 0.0;
      As[1445] = 0.0;
      As[1446] = 0.0;
      As[1447] = 0.0;
      As[1448] = 0.0;
      As[1449] = 0.0;
      As[1450] = 0.0;
      As[1451] = 0.0;
      As[1452] = 0.0;
      As[1453] = 0.0;
      As[1454] = 0.0;
      As[1455] = 0.0;
      As[1456] = 0.0;
      As[1457] = 0.0;
      As[1458] = 0.0;
      As[1459] = 0.0;
      As[1460] = 0.0;
      As[1461] = 0.0;
      As[1462] = 0.0;
      As[1463] = 0.0;
      As[1464] = 0.0;
      As[1465] = 0.0;
      As[1466] = 0.0;
      As[1467] = 0.0033003117895002514;
      As[1468] = 0.0033003117895002514;
      As[1469] = -0.006600678309509651;
      As[1470] = 0.0033003117895002514;
      As[1471] = 0.0033003117895002514;
      As[1472] = 0.99339932169049028;
      As[1473] = 0.0;
      As[1474] = 0.0;
      As[1475] = 0.0;
      As[1476] = 0.0;
      As[1477] = 0.0;
      As[1478] = 0.0;
      As[1479] = 0.0;
      As[1480] = 0.0;
      As[1481] = 0.0;
      As[1482] = 0.0;
      As[1483] = 0.0;
      As[1484] = 0.0;
      As[1485] = 0.0;
      As[1486] = 0.0;
      As[1487] = -2.4919075328810973E-13;
      As[1488] = -2.4919075328809206E-13;
      As[1489] = -2.3444139097650854E-7;
      As[1490] = -1.6241854251584936E-11;
      As[1491] = -1.6241854251588939E-11;
      As[1492] = -3.4565861692519442E-9;
      As[1493] = 2.0637756383613037E-8;
      As[1494] = 1.7282930846259725E-9;
      As[1495] = 2.0637756383613017E-8;
      As[1496] = 4.7498967454352377E-24;
      As[1497] = 1.7282930846259719E-9;
      As[1498] = 3.2028619574743118E-8;
      As[1499] = 3.2028619574743138E-8;
      As[1500] = -6.4057239149486263E-8;
      As[1501] = 0.00013843784885799788;
      As[1502] = -1.0665794540782114E-20;
      As[1503] = -0.00013843784885799788;
      As[1504] = 4.3891575838830726E-9;
      As[1505] = 5.1069458127138E-26;
      As[1506] = -4.3891575838830709E-9;
      As[1507] = 1.1320755998159064E-24;
      As[1508] = 2.3537468597857916E-5;
      As[1509] = 1.8048737011165425E-21;
      As[1510] = -2.353746859785791E-5;
      As[1511] = -6.4479147035602763E-21;
      As[1512] = 0.0;
      As[1513] = 0.0;
      As[1514] = 0.0;
      As[1515] = 0.0;
      As[1516] = 0.0;
      As[1517] = 0.0;
      As[1518] = -0.76580852621493622;
      As[1519] = -4.1121492860861024E-7;
      As[1520] = 1.6217531037632857E-5;
      As[1521] = -2.8476234272369759E-11;
      As[1522] = 1.6217531037632844E-5;
      As[1523] = -2.8476234072104687E-11;
      As[1524] = -1.4265182629474617E-8;
      As[1525] = -3.2029175230106752E-8;
      As[1526] = -1.4265182629474612E-8;
      As[1527] = -3.2029175230106752E-8;
      As[1528] = 2.8530365258949233E-8;
      As[1529] = 6.4058350460213477E-8;
      As[1530] = 0.0;
      As[1531] = 0.0;
      As[1532] = 2.1259022226482947E-11;
      As[1533] = 2.1259022226482947E-11;
      As[1534] = 1.9999336168423295E-5;
      As[1535] = 7.6019014100191173E-13;
      As[1536] = 7.6019014100175118E-13;
      As[1537] = 3.4419992602974417E-11;
      As[1538] = -2.0550664363147615E-10;
      As[1539] = -1.7209996301487202E-11;
      As[1540] = -2.0550664363147607E-10;
      As[1541] = -5.9464305849873823E-26;
      As[1542] = -1.7209996301487208E-11;
      As[1543] = -3.1893457731583767E-10;
      As[1544] = -3.1893457731583783E-10;
      As[1545] = 6.3786915463167565E-10;
      As[1546] = -1.3785363651718824E-6;
      As[1547] = 1.061888163207754E-22;
      As[1548] = 1.3785363651718828E-6;
      As[1549] = -4.3706351924459389E-11;
      As[1550] = -5.5953500481896835E-27;
      As[1551] = 4.3706351924459383E-11;
      As[1552] = -3.5290773904606413E-27;
      As[1553] = -2.3438139695106807E-7;
      As[1554] = 1.5527522761427934E-23;
      As[1555] = 2.3438139695106809E-7;
      As[1556] = -1.5481941954535261E-23;
      As[1557] = 0.0;
      As[1558] = 0.0;
      As[1559] = 0.0;
      As[1560] = 0.0;
      As[1561] = 0.0;
      As[1562] = 0.0;
      As[1563] = -0.002332031779524834;
      As[1564] = -0.99996492076023435;
      As[1565] = -1.6149092514970468E-7;
      As[1566] = 2.8356051783370643E-13;
      As[1567] = -1.614909251497047E-7;
      As[1568] = 2.8356052584633971E-13;
      As[1569] = 1.420498308285647E-10;
      As[1570] = 3.1894011042063211E-10;
      As[1571] = 1.4204983082856478E-10;
      As[1572] = 3.1894011042063226E-10;
      As[1573] = -2.8409966165712956E-10;
      As[1574] = -6.3788022084126432E-10;
      As[1575] = 0.0;
      As[1576] = 0.0;
      As[1577] = 2.4919075328802082E-13;
      As[1578] = 6.3888388527414035E-29;
      As[1579] = -1.62418542515622E-11;
      As[1580] = -2.3444139097658419E-7;
      As[1581] = -1.6241854251588942E-11;
      As[1582] = 1.7282930846259736E-9;
      As[1583] = -2.0637756383613044E-8;
      As[1584] = -3.4565861692519463E-9;
      As[1585] = -2.5853006756134556E-24;
      As[1586] = 2.0637756383613037E-8;
      As[1587] = 1.7282930846259727E-9;
      As[1588] = -6.4057239149486276E-8;
      As[1589] = 3.2028619574743145E-8;
      As[1590] = 3.2028619574743138E-8;
      As[1591] = -0.00013843784885799793;
      As[1592] = 0.000138437848857998;
      As[1593] = 5.25706222751371E-20;
      As[1594] = -4.3891575838830759E-9;
      As[1595] = 4.3891575838830742E-9;
      As[1596] = 4.6586819034474116E-25;
      As[1597] = 7.3372656269504864E-25;
      As[1598] = -2.3537468597857916E-5;
      As[1599] = 2.3537468597857923E-5;
      As[1600] = 4.6680213872640374E-21;
      As[1601] = -4.7609667748466689E-22;
      As[1602] = 0.0;
      As[1603] = 0.0;
      As[1604] = 0.0;
      As[1605] = 0.0;
      As[1606] = 0.0;
      As[1607] = 0.0;
      As[1608] = 1.6217531037632874E-5;
      As[1609] = -2.8476235409597062E-11;
      As[1610] = -0.76580852621493634;
      As[1611] = -4.1121492482639119E-7;
      As[1612] = 1.6217531037632871E-5;
      As[1613] = -2.8476234072104687E-11;
      As[1614] = 2.853036525894924E-8;
      As[1615] = 6.4058350460213517E-8;
      As[1616] = -1.4265182629474618E-8;
      As[1617] = -3.2029175230106758E-8;
      As[1618] = -1.426518262947463E-8;
      As[1619] = -3.2029175230106745E-8;
      As[1620] = 0.0;
      As[1621] = 0.0;
      As[1622] = -2.1259022226454264E-11;
      As[1623] = -2.4000133979514369E-27;
      As[1624] = 7.6019014100192506E-13;
      As[1625] = 1.9999336168423298E-5;
      As[1626] = 7.6019014100175158E-13;
      As[1627] = -1.7209996301487208E-11;
      As[1628] = 2.0550664363147607E-10;
      As[1629] = 3.4419992602974417E-11;
      As[1630] = 4.2638767308831431E-26;
      As[1631] = -2.0550664363147617E-10;
      As[1632] = -1.7209996301487208E-11;
      As[1633] = 6.3786915463167565E-10;
      As[1634] = -3.1893457731583788E-10;
      As[1635] = -3.1893457731583788E-10;
      As[1636] = 1.3785363651718817E-6;
      As[1637] = -1.3785363651718834E-6;
      As[1638] = -7.41874707805748E-22;
      As[1639] = 4.3706351924459383E-11;
      As[1640] = -4.3706351924459389E-11;
      As[1641] = -4.2006439677589347E-27;
      As[1642] = -1.2518717197707276E-26;
      As[1643] = 2.3438139695106812E-7;
      As[1644] = -2.3438139695106809E-7;
      As[1645] = 7.6674617198685858E-24;
      As[1646] = -1.6818688941418582E-22;
      As[1647] = 0.0;
      As[1648] = 0.0;
      As[1649] = 0.0;
      As[1650] = 0.0;
      As[1651] = 0.0;
      As[1652] = 0.0;
      As[1653] = -1.6149092514970454E-7;
      As[1654] = 2.83560517182327E-13;
      As[1655] = -0.0023320317795248345;
      As[1656] = -0.99996492076023458;
      As[1657] = -1.614909251497047E-7;
      As[1658] = 2.8356052584633971E-13;
      As[1659] = -2.8409966165712945E-10;
      As[1660] = -6.3788022084126473E-10;
      As[1661] = 1.4204983082856478E-10;
      As[1662] = 3.1894011042063231E-10;
      As[1663] = 1.4204983082856473E-10;
      As[1664] = 3.1894011042063221E-10;
      As[1665] = 0.0;
      As[1666] = 0.0;
      As[1667] = -1.2311113365659006E-29;
      As[1668] = 2.4919075328791874E-13;
      As[1669] = -1.6241854251582571E-11;
      As[1670] = -1.6241854251576182E-11;
      As[1671] = -2.3444139097655446E-7;
      As[1672] = 1.7282930846259731E-9;
      As[1673] = -6.5896309532083185E-24;
      As[1674] = 1.7282930846259736E-9;
      As[1675] = -2.0637756383613034E-8;
      As[1676] = -2.0637756383613053E-8;
      As[1677] = -3.4565861692519467E-9;
      As[1678] = 3.2028619574743138E-8;
      As[1679] = -6.4057239149486263E-8;
      As[1680] = 3.2028619574743138E-8;
      As[1681] = -4.6198095946151665E-20;
      As[1682] = -0.000138437848857998;
      As[1683] = 0.00013843784885799783;
      As[1684] = -1.6072417043077026E-25;
      As[1685] = -4.389157583883075E-9;
      As[1686] = 4.3891575838830734E-9;
      As[1687] = -2.0333418182508569E-24;
      As[1688] = -2.0470155316832677E-25;
      As[1689] = -2.3537468597857916E-5;
      As[1690] = 2.3537468597857913E-5;
      As[1691] = -9.3391103833774581E-21;
      As[1692] = 0.0;
      As[1693] = 0.0;
      As[1694] = 0.0;
      As[1695] = 0.0;
      As[1696] = 0.0;
      As[1697] = 0.0;
      As[1698] = 1.6217531037632854E-5;
      As[1699] = -2.8476234390649591E-11;
      As[1700] = 1.6217531037632867E-5;
      As[1701] = -2.8476234710419554E-11;
      As[1702] = -0.76580852621493622;
      As[1703] = -4.112149263130176E-7;
      As[1704] = -1.4265182629474627E-8;
      As[1705] = -3.2029175230106758E-8;
      As[1706] = 2.853036525894926E-8;
      As[1707] = 6.40583504602135E-8;
      As[1708] = -1.4265182629474627E-8;
      As[1709] = -3.2029175230106745E-8;
      As[1710] = 0.0;
      As[1711] = 0.0;
      As[1712] = -3.2165389969871059E-31;
      As[1713] = -2.1259022226476944E-11;
      As[1714] = 7.60190141002058E-13;
      As[1715] = 7.601901410023646E-13;
      As[1716] = 1.9999336168423295E-5;
      As[1717] = -1.72099963014872E-11;
      As[1718] = 6.2815251268518591E-26;
      As[1719] = -1.7209996301487208E-11;
      As[1720] = 2.0550664363147594E-10;
      As[1721] = 2.055066436314762E-10;
      As[1722] = 3.4419992602974404E-11;
      As[1723] = -3.1893457731583783E-10;
      As[1724] = 6.3786915463167576E-10;
      As[1725] = -3.1893457731583767E-10;
      As[1726] = 4.63398550055088E-22;
      As[1727] = 1.3785363651718832E-6;
      As[1728] = -1.3785363651718813E-6;
      As[1729] = -2.7509513482359723E-27;
      As[1730] = 4.3706351924459383E-11;
      As[1731] = -4.3706351924459357E-11;
      As[1732] = 2.5569611199912214E-26;
      As[1733] = 1.0283643956410022E-26;
      As[1734] = 2.3438139695106809E-7;
      As[1735] = -2.3438139695106802E-7;
      As[1736] = 9.3664751130837253E-23;
      As[1737] = 0.0;
      As[1738] = 0.0;
      As[1739] = 0.0;
      As[1740] = 0.0;
      As[1741] = 0.0;
      As[1742] = 0.0;
      As[1743] = -1.6149092514970441E-7;
      As[1744] = 2.8356051051030569E-13;
      As[1745] = -1.6149092514970475E-7;
      As[1746] = 2.835604951876558E-13;
      As[1747] = -0.002332031779524834;
      As[1748] = -0.99996492076023435;
      As[1749] = 1.4204983082856473E-10;
      As[1750] = 3.1894011042063231E-10;
      As[1751] = -2.8409966165712945E-10;
      As[1752] = -6.3788022084126463E-10;
      As[1753] = 1.4204983082856473E-10;
      As[1754] = 3.1894011042063195E-10;
      As[1755] = 0.0;
      As[1756] = 0.0;
      As[1757] = 5.3509066708801543E-20;
      As[1758] = 2.1390599793897255E-32;
      As[1759] = 1.6779476200609223E-14;
      As[1760] = -3.355895240117439E-14;
      As[1761] = 1.6779476200592966E-14;
      As[1762] = -5.1301662925854235E-5;
      As[1763] = 0.00061259935074446073;
      As[1764] = 0.00010260332585170848;
      As[1765] = 1.6001216714383811E-19;
      As[1766] = -0.00061259935074446106;
      As[1767] = -5.1301662925854242E-5;
      As[1768] = -2.5705055061763424;
      As[1769] = -0.092720260420749517;
      As[1770] = -0.092720260420749476;
      As[1771] = -0.00076628479096019491;
      As[1772] = 0.000766284790960195;
      As[1773] = 1.6365558087616359E-19;
      As[1774] = 0.00013028524206908833;
      As[1775] = -0.0001302852420690883;
      As[1776] = 6.1707411731195124E-21;
      As[1777] = -4.7471824522588141E-20;
      As[1778] = 2.4327093174710415E-8;
      As[1779] = -2.432709317471649E-8;
      As[1780] = 9.66320620466285E-21;
      As[1781] = -2.171190755012199E-20;
      As[1782] = 0.0;
      As[1783] = 0.0;
      As[1784] = 0.0;
      As[1785] = 0.0;
      As[1786] = 0.0;
      As[1787] = 0.0;
      As[1788] = -1.676158958961707E-8;
      As[1789] = 2.9431538280839229E-14;
      As[1790] = 3.3523179179307428E-8;
      As[1791] = -5.8863079930924678E-14;
      As[1792] = -1.6761589589639652E-8;
      As[1793] = 2.9431541036654669E-14;
      As[1794] = 0.07088428867714483;
      As[1795] = 2.5705552545648871;
      As[1796] = -0.032693059226855475;
      As[1797] = 0.092727022421495456;
      As[1798] = -0.032693059226855475;
      As[1799] = 0.092727022421495414;
      As[1800] = 0.0;
      As[1801] = 0.0;
      As[1802] = 2.5029637406485289E-24;
      As[1803] = 5.4788737851044109E-36;
      As[1804] = 7.8488419066843534E-19;
      As[1805] = -1.5697683813348681E-18;
      As[1806] = 7.8488419066326412E-19;
      As[1807] = -2.3997092462340671E-9;
      As[1808] = 2.8655217830718891E-8;
      As[1809] = 4.7994184924681341E-9;
      As[1810] = 6.0673451337619828E-24;
      As[1811] = -2.8655217830718908E-8;
      As[1812] = -2.3997092462340671E-9;
      As[1813] = 0.0001486087930494043;
      As[1814] = 5.2539122504697077E-6;
      As[1815] = 5.2539122504697094E-6;
      As[1816] = -3.5844075868912934E-8;
      As[1817] = 3.5844075868913027E-8;
      As[1818] = 2.5617768894101513E-23;
      As[1819] = 6.0942800332398573E-9;
      As[1820] = -6.0942800332398556E-9;
      As[1821] = 1.1907764167491587E-24;
      As[1822] = -4.0892933858894761E-24;
      As[1823] = 1.137934856218618E-12;
      As[1824] = -1.1379348562175939E-12;
      As[1825] = -2.5093297521367987E-25;
      As[1826] = 4.5909521236248637E-24;
      As[1827] = 0.0;
      As[1828] = 0.0;
      As[1829] = 0.0;
      As[1830] = 0.0;
      As[1831] = 0.0;
      As[1832] = 0.0;
      As[1833] = -7.840475186456503E-13;
      As[1834] = 1.3767027965755217E-18;
      As[1835] = 1.5680950372954185E-12;
      As[1836] = -2.7534054927173647E-18;
      As[1837] = -7.8404751865000774E-13;
      As[1838] = 1.3767028508872655E-18;
      As[1839] = 5.35532344701018E-5;
      As[1840] = 0.99985138833124221;
      As[1841] = 1.931812967114485E-6;
      As[1842] = -5.2543009358721546E-6;
      As[1843] = 1.9318129671144846E-6;
      As[1844] = -5.2543009358721546E-6;
      As[1845] = 0.0;
      As[1846] = 0.0;
      As[1847] = -8.9845605089761351E-33;
      As[1848] = 5.3509066708699619E-20;
      As[1849] = 1.6779476200614745E-14;
      As[1850] = 1.6779476200630324E-14;
      As[1851] = -3.3558952401042272E-14;
      As[1852] = -5.1301662925854276E-5;
      As[1853] = 2.4905260676471863E-20;
      As[1854] = -5.1301662925854269E-5;
      As[1855] = 0.00061259935074446138;
      As[1856] = 0.00061259935074446127;
      As[1857] = 0.00010260332585170855;
      As[1858] = -0.092720260420749462;
      As[1859] = -2.5705055061763438;
      As[1860] = -0.092720260420749434;
      As[1861] = -1.1901946201734029E-19;
      As[1862] = -0.00076628479096019545;
      As[1863] = 0.00076628479096019534;
      As[1864] = 5.127284363702182E-21;
      As[1865] = 0.00013028524206908841;
      As[1866] = -0.00013028524206908844;
      As[1867] = 1.638274468432435E-21;
      As[1868] = -4.5325085825850447E-22;
      As[1869] = 2.4327093174702951E-8;
      As[1870] = -2.4327093174723974E-8;
      As[1871] = -1.5782235762476655E-20;
      As[1872] = 0.0;
      As[1873] = 0.0;
      As[1874] = 0.0;
      As[1875] = 0.0;
      As[1876] = 0.0;
      As[1877] = 0.0;
      As[1878] = -1.6761589589659458E-8;
      As[1879] = 2.9431539406321951E-14;
      As[1880] = -1.6761589589572749E-8;
      As[1881] = 2.9431538142547205E-14;
      As[1882] = 3.3523179179346967E-8;
      As[1883] = -5.8863085931413745E-14;
      As[1884] = -0.032693059226855517;
      As[1885] = 0.09272702242149547;
      As[1886] = 0.070884288677145149;
      As[1887] = 2.5705552545648875;
      As[1888] = -0.032693059226855538;
      As[1889] = 0.0927270224214954;
      As[1890] = 0.0;
      As[1891] = 0.0;
      As[1892] = 2.2152558954377212E-36;
      As[1893] = 2.5029637406496353E-24;
      As[1894] = 7.848841906695598E-19;
      As[1895] = 7.848841906678723E-19;
      As[1896] = -1.5697683813359347E-18;
      As[1897] = -2.3997092462340679E-9;
      As[1898] = -1.7263132965400081E-24;
      As[1899] = -2.3997092462340675E-9;
      As[1900] = 2.8655217830718911E-8;
      As[1901] = 2.8655217830718905E-8;
      As[1902] = 4.7994184924681349E-9;
      As[1903] = 5.2539122504697077E-6;
      As[1904] = 0.00014860879304940427;
      As[1905] = 5.2539122504697035E-6;
      As[1906] = 8.1490896648221863E-25;
      As[1907] = -3.5844075868912993E-8;
      As[1908] = 3.5844075868912954E-8;
      As[1909] = -3.6899930217057871E-25;
      As[1910] = 6.0942800332398573E-9;
      As[1911] = -6.0942800332398581E-9;
      As[1912] = 1.4897361418492407E-25;
      As[1913] = 3.0718339373660696E-26;
      As[1914] = 1.1379348562180071E-12;
      As[1915] = -1.1379348562219263E-12;
      As[1916] = -2.1650641855908445E-25;
      As[1917] = 0.0;
      As[1918] = 0.0;
      As[1919] = 0.0;
      As[1920] = 0.0;
      As[1921] = 0.0;
      As[1922] = 0.0;
      As[1923] = -7.8404751866457639E-13;
      As[1924] = 1.3767027287776575E-18;
      As[1925] = -7.84047518646543E-13;
      As[1926] = 1.3767027306174441E-18;
      As[1927] = 1.5680950373183441E-12;
      As[1928] = -2.7534054236542118E-18;
      As[1929] = 1.9318129671144871E-6;
      As[1930] = -5.2543009358721529E-6;
      As[1931] = 5.3553234470101795E-5;
      As[1932] = 0.99985138833124232;
      As[1933] = 1.9318129671144871E-6;
      As[1934] = -5.2543009358721487E-6;
      As[1935] = 0.0;
      As[1936] = 0.0;
      As[1937] = -5.3509066708647421E-20;
      As[1938] = -5.3509066708648522E-20;
      As[1939] = -3.355895240106874E-14;
      As[1940] = 1.6779476200617945E-14;
      As[1941] = 1.6779476200557834E-14;
      As[1942] = 0.0001026033258517085;
      As[1943] = -0.000612599350744461;
      As[1944] = -5.1301662925854255E-5;
      As[1945] = -0.00061259935074446127;
      As[1946] = 1.5188073024845491E-19;
      As[1947] = -5.1301662925854255E-5;
      As[1948] = -0.092720260420749476;
      As[1949] = -0.0927202604207495;
      As[1950] = -2.5705055061763424;
      As[1951] = 0.000766284790960195;
      As[1952] = -1.1739778028905536E-20;
      As[1953] = -0.00076628479096019534;
      As[1954] = -0.00013028524206908836;
      As[1955] = -3.9799622688209553E-21;
      As[1956] = 0.00013028524206908836;
      As[1957] = 3.96331784118649E-20;
      As[1958] = -2.4327093174710422E-8;
      As[1959] = 9.3845296549403263E-21;
      As[1960] = 2.4327093174716497E-8;
      As[1961] = 4.1941591636447214E-21;
      As[1962] = 0.0;
      As[1963] = 0.0;
      As[1964] = 0.0;
      As[1965] = 0.0;
      As[1966] = 0.0;
      As[1967] = 0.0;
      As[1968] = 3.352317917934929E-8;
      As[1969] = -5.8863084572364822E-14;
      As[1970] = -1.676158958965557E-8;
      As[1971] = 2.9431538971543687E-14;
      As[1972] = -1.6761589589630394E-8;
      As[1973] = 2.9431541218813286E-14;
      As[1974] = -0.0326930592268555;
      As[1975] = 0.092727022421495414;
      As[1976] = -0.032693059226855489;
      As[1977] = 0.092727022421495484;
      As[1978] = 0.070884288677145343;
      As[1979] = 2.5705552545648862;
      As[1980] = 0.0;
      As[1981] = 0.0;
      As[1982] = -2.5029637406487188E-24;
      As[1983] = -2.5029637406481076E-24;
      As[1984] = -1.5697683813371469E-18;
      As[1985] = 7.8488419066857334E-19;
      As[1986] = 7.8488419066677442E-19;
      As[1987] = 4.7994184924681358E-9;
      As[1988] = -2.8655217830718921E-8;
      As[1989] = -2.3997092462340675E-9;
      As[1990] = -2.86552178307189E-8;
      As[1991] = -1.3354998606833773E-23;
      As[1992] = -2.3997092462340679E-9;
      As[1993] = 5.2539122504697068E-6;
      As[1994] = 5.2539122504697094E-6;
      As[1995] = 0.00014860879304940427;
      As[1996] = 3.5844075868912987E-8;
      As[1997] = 2.4346553027288697E-25;
      As[1998] = -3.5844075868913007E-8;
      As[1999] = -6.09428003323986E-9;
      As[2000] = -1.2131351874816978E-24;
      As[2001] = 6.0942800332398589E-9;
      As[2002] = 1.4355247143307806E-27;
      As[2003] = -1.1379348562198592E-12;
      As[2004] = -1.1556253309930041E-24;
      As[2005] = 1.1379348562163533E-12;
      As[2006] = -4.5037552561795514E-29;
      As[2007] = 0.0;
      As[2008] = 0.0;
      As[2009] = 0.0;
      As[2010] = 0.0;
      As[2011] = 0.0;
      As[2012] = 0.0;
      As[2013] = 1.5680950372885569E-12;
      As[2014] = -2.7534054838104379E-18;
      As[2015] = -7.8404751864710755E-13;
      As[2016] = 1.376702752425468E-18;
      As[2017] = -7.8404751864875626E-13;
      As[2018] = 1.3767027094672585E-18;
      As[2019] = 1.9318129671144862E-6;
      As[2020] = -5.2543009358721529E-6;
      As[2021] = 1.9318129671144858E-6;
      As[2022] = -5.2543009358721563E-6;
      As[2023] = 5.3553234470101795E-5;
      As[2024] = 0.99985138833124243;
      Bs[0] = 0.0;
      Bs[1] = 0.0;
      Bs[2] = 0.0;
      Bs[3] = 0.0;
      Bs[4] = 0.0;
      Bs[5] = 0.0;
      Bs[6] = 2.2222222221851853;
      Bs[7] = 0.0;
      Bs[8] = 0.0;
      Bs[9] = 0.0;
      Bs[10] = 0.0;
      Bs[11] = 0.0;
      Bs[12] = 0.0;
      Bs[13] = -1.1111111110925924;
      Bs[14] = -1.1111111110925926;
      Bs[15] = -1.1111111110925924;
      Bs[16] = -1.1111111110925926;
      Bs[17] = -1.1111111110925924;
      Bs[18] = -1.1111111110925926;
      Bs[19] = 0.0;
      Bs[20] = 0.0;
      Bs[21] = 0.0;
      Bs[22] = 0.0;
      Bs[23] = 0.0;
      Bs[24] = 0.0;
      Bs[25] = 0.0;
      Bs[26] = 0.0;
      Bs[27] = -454047.61518257356;
      Bs[28] = 90809.523036514714;
      Bs[29] = 454047.61518257356;
      Bs[30] = 90809.523036514714;
      Bs[31] = 2.7638907119321683;
      Bs[32] = 2.7638907119321687;
      Bs[33] = 2.7638907119321687;
      Bs[34] = 2.7638907119321692;
      Bs[35] = 2.7638907119321687;
      Bs[36] = 2.7638907119321687;
      Bs[37] = 0.0;
      Bs[38] = 0.0;
      Bs[39] = 0.0;
      Bs[40] = 0.0;
      Bs[41] = 0.0;
      Bs[42] = 0.0;
      Bs[43] = 4.9946609044459236E-16;
      Bs[44] = -5.1587778558400328E-18;
      Bs[45] = 0.0;
      Bs[46] = 0.0;
      Bs[47] = 0.0;
      Bs[48] = -0.029926346251814642;
      Bs[49] = 0.059852692503629269;
      Bs[50] = -0.029926346251814642;
      Bs[51] = 0.0;
      Bs[52] = 0.0;
      Bs[53] = 0.0;
      Bs[54] = 0.0;
      Bs[55] = 0.0;
      Bs[56] = 0.0;
      Bs[57] = 0.0;
      Bs[58] = 0.0;
      Bs[59] = 0.0;
      Bs[60] = 0.0;
      Bs[61] = 333293.84233176667;
      Bs[62] = -333293.84233176667;
      Bs[63] = -666587.68466353335;
      Bs[64] = 666587.68466353335;
      Bs[65] = 333293.84233176667;
      Bs[66] = -333293.84233176667;
      Bs[67] = 0.0;
      Bs[68] = 0.0;
      Bs[69] = -3.5661572517644239E-16;
      Bs[70] = 1.7830786258822134E-16;
      Bs[71] = 1.783078625882213E-16;
      Bs[72] = -0.029926346251814642;
      Bs[73] = -0.029926346251814635;
      Bs[74] = 0.059852692503629283;
      Bs[75] = 0.0;
      Bs[76] = 0.0;
      Bs[77] = 0.0;
      Bs[78] = 0.0;
      Bs[79] = 0.0;
      Bs[80] = 0.0;
      Bs[81] = 0.0;
      Bs[82] = 0.0;
      Bs[83] = 0.0;
      Bs[84] = 0.0;
      Bs[85] = 333293.84233176662;
      Bs[86] = -333293.84233176662;
      Bs[87] = 333293.84233176662;
      Bs[88] = -333293.84233176662;
      Bs[89] = -666587.68466353323;
      Bs[90] = 666587.68466353323;
      Bs[91] = 0.0;
      Bs[92] = 0.0;
      Bs[93] = 1.7830786258822127E-16;
      Bs[94] = -3.5661572517644259E-16;
      Bs[95] = 1.7830786258822125E-16;
      Bs[96] = 166.11828794289275;
      Bs[97] = 6.3142838177635162E-6;
      Bs[98] = 6.3142838177593582E-6;
      Bs[99] = 0.0;
      Bs[100] = 0.0;
      Bs[101] = 0.0;
      Bs[102] = 0.0;
      Bs[103] = 0.0;
      Bs[104] = 0.0;
      Bs[105] = 0.0;
      Bs[106] = 0.0;
      Bs[107] = 0.0;
      Bs[108] = 0.0;
      Bs[109] = 110.73607050251023;
      Bs[110] = -110.73607050251023;
      Bs[111] = -55.368035251255115;
      Bs[112] = 55.368035251255115;
      Bs[113] = -55.368035251255108;
      Bs[114] = 55.368035251255108;
      Bs[115] = 0.0;
      Bs[116] = 0.0;
      Bs[117] = 3.2992328882814369E-13;
      Bs[118] = 3.2992328882814394E-13;
      Bs[119] = -6.5984657765628758E-13;
      Bs[120] = 6.31428381775857E-6;
      Bs[121] = 166.11828794289278;
      Bs[122] = 6.3142838177593582E-6;
      Bs[123] = 0.0;
      Bs[124] = 0.0;
      Bs[125] = 0.0;
      Bs[126] = 0.0;
      Bs[127] = 0.0;
      Bs[128] = 0.0;
      Bs[129] = 0.0;
      Bs[130] = 0.0;
      Bs[131] = 0.0;
      Bs[132] = 0.0;
      Bs[133] = -55.368035251255115;
      Bs[134] = 55.368035251255115;
      Bs[135] = 110.73607050251023;
      Bs[136] = -110.73607050251023;
      Bs[137] = -55.368035251255108;
      Bs[138] = 55.368035251255108;
      Bs[139] = 0.0;
      Bs[140] = 0.0;
      Bs[141] = -6.5984657765628738E-13;
      Bs[142] = 3.2992328882814394E-13;
      Bs[143] = 3.2992328882814389E-13;
      Bs[144] = 6.3142838177639406E-6;
      Bs[145] = 6.3142838177627641E-6;
      Bs[146] = 166.11828794289275;
      Bs[147] = 0.0;
      Bs[148] = 0.0;
      Bs[149] = 0.0;
      Bs[150] = 0.0;
      Bs[151] = 0.0;
      Bs[152] = 0.0;
      Bs[153] = 0.0;
      Bs[154] = 0.0;
      Bs[155] = 0.0;
      Bs[156] = 0.0;
      Bs[157] = -55.368035251255108;
      Bs[158] = 55.368035251255108;
      Bs[159] = -55.368035251255108;
      Bs[160] = 55.368035251255108;
      Bs[161] = 110.73607050251022;
      Bs[162] = -110.73607050251022;
      Bs[163] = 0.0;
      Bs[164] = 0.0;
      Bs[165] = 3.2992328882814379E-13;
      Bs[166] = -6.5984657765628768E-13;
      Bs[167] = 3.2992328882814379E-13;
      Bs[168] = -6.2284951152259743E-7;
      Bs[169] = 3.1142475599978428E-7;
      Bs[170] = 3.1142475582345082E-7;
      Bs[171] = 0.0;
      Bs[172] = 0.0;
      Bs[173] = 0.0;
      Bs[174] = 0.0;
      Bs[175] = 0.0;
      Bs[176] = 0.0;
      Bs[177] = 0.0;
      Bs[178] = 0.0;
      Bs[179] = 0.0;
      Bs[180] = 0.0;
      Bs[181] = 5.5920683643585638E-14;
      Bs[182] = -5.5920683643585638E-14;
      Bs[183] = -2.7960341829708615E-14;
      Bs[184] = 2.7960341829708615E-14;
      Bs[185] = -2.7960341813877023E-14;
      Bs[186] = 2.7960341813877023E-14;
      Bs[187] = 0.0;
      Bs[188] = 0.0;
      Bs[189] = 2.1975385702209368E-6;
      Bs[190] = 2.1975385702168253E-6;
      Bs[191] = -4.3950771404373317E-6;
      Bs[192] = 0.00082639124702251279;
      Bs[193] = -0.00082639124702811454;
      Bs[194] = -7.4756501816979557E-13;
      Bs[195] = 0.0;
      Bs[196] = 0.0;
      Bs[197] = 0.0;
      Bs[198] = 0.0;
      Bs[199] = 0.0;
      Bs[200] = 0.0;
      Bs[201] = 0.0;
      Bs[202] = 0.0;
      Bs[203] = 0.0;
      Bs[204] = 0.0;
      Bs[205] = -7.4195070626625164E-11;
      Bs[206] = 7.4195070626625164E-11;
      Bs[207] = 7.4195070582047561E-11;
      Bs[208] = -7.4195070582047561E-11;
      Bs[209] = 4.4577598223599725E-20;
      Bs[210] = -4.4577598223599725E-20;
      Bs[211] = 0.0;
      Bs[212] = 0.0;
      Bs[213] = -0.0058313496399968294;
      Bs[214] = 8.9582432182710918E-16;
      Bs[215] = 0.0058313496399966915;
      Bs[216] = 3.1142475631245977E-7;
      Bs[217] = -6.2284951217892568E-7;
      Bs[218] = 3.1142475602010895E-7;
      Bs[219] = 0.0;
      Bs[220] = 0.0;
      Bs[221] = 0.0;
      Bs[222] = 0.0;
      Bs[223] = 0.0;
      Bs[224] = 0.0;
      Bs[225] = 0.0;
      Bs[226] = 0.0;
      Bs[227] = 0.0;
      Bs[228] = 0.0;
      Bs[229] = -2.7960341862180405E-14;
      Bs[230] = 2.7960341862180405E-14;
      Bs[231] = 5.5920683698112966E-14;
      Bs[232] = -5.5920683698112966E-14;
      Bs[233] = -2.7960341835932561E-14;
      Bs[234] = 2.7960341835932561E-14;
      Bs[235] = 0.0;
      Bs[236] = 0.0;
      Bs[237] = -4.39507714043951E-6;
      Bs[238] = 2.1975385702209346E-6;
      Bs[239] = 2.1975385702200897E-6;
      Bs[240] = 0.00082639124547984088;
      Bs[241] = 9.5042688043235131E-14;
      Bs[242] = -0.00082639124536221687;
      Bs[243] = 0.0;
      Bs[244] = 0.0;
      Bs[245] = 0.0;
      Bs[246] = 0.0;
      Bs[247] = 0.0;
      Bs[248] = 0.0;
      Bs[249] = 0.0;
      Bs[250] = 0.0;
      Bs[251] = 0.0;
      Bs[252] = 0.0;
      Bs[253] = -7.4195070459216138E-11;
      Bs[254] = 7.4195070459216138E-11;
      Bs[255] = -2.1685702980814023E-21;
      Bs[256] = 2.1685702980814023E-21;
      Bs[257] = 7.4195070461384722E-11;
      Bs[258] = -7.4195070461384722E-11;
      Bs[259] = 0.0;
      Bs[260] = 0.0;
      Bs[261] = 6.0715321659188248E-16;
      Bs[262] = -0.0058313496399942221;
      Bs[263] = 0.0058313496399915819;
      Bs[264] = 9.2472944202730534E-14;
      Bs[265] = 0.00082639124547810073;
      Bs[266] = -0.00082639124788457831;
      Bs[267] = 0.0;
      Bs[268] = 0.0;
      Bs[269] = 0.0;
      Bs[270] = 0.0;
      Bs[271] = 0.0;
      Bs[272] = 0.0;
      Bs[273] = 0.0;
      Bs[274] = 0.0;
      Bs[275] = 0.0;
      Bs[276] = 0.0;
      Bs[277] = -7.755439507685594E-20;
      Bs[278] = 7.755439507685594E-20;
      Bs[279] = -7.4195070534676465E-11;
      Bs[280] = 7.4195070534676465E-11;
      Bs[281] = 7.4195070612230865E-11;
      Bs[282] = -7.4195070612230865E-11;
      Bs[283] = 0.0;
      Bs[284] = 0.0;
      Bs[285] = 0.0058313496399915662;
      Bs[286] = -0.0058313496399943427;
      Bs[287] = 6.9161776021935342E-16;
      Bs[288] = 3.1142475538324621E-7;
      Bs[289] = 3.11424755565166E-7;
      Bs[290] = -6.2284951140342826E-7;
      Bs[291] = 0.0;
      Bs[292] = 0.0;
      Bs[293] = 0.0;
      Bs[294] = 0.0;
      Bs[295] = 0.0;
      Bs[296] = 0.0;
      Bs[297] = 0.0;
      Bs[298] = 0.0;
      Bs[299] = 0.0;
      Bs[300] = 0.0;
      Bs[301] = -2.7960341796969275E-14;
      Bs[302] = 2.7960341796969275E-14;
      Bs[303] = -2.7960341813302407E-14;
      Bs[304] = 2.7960341813302407E-14;
      Bs[305] = 5.5920683610271682E-14;
      Bs[306] = -5.5920683610271682E-14;
      Bs[307] = 0.0;
      Bs[308] = 0.0;
      Bs[309] = 2.19753857021852E-6;
      Bs[310] = -4.3950771404382355E-6;
      Bs[311] = 2.1975385702184456E-6;
      Bs[312] = 0.00019860087161987347;
      Bs[313] = -0.00039720174323934813;
      Bs[314] = 0.00019860087161951642;
      Bs[315] = 0.0;
      Bs[316] = 0.0;
      Bs[317] = 0.0;
      Bs[318] = 0.0;
      Bs[319] = 0.0;
      Bs[320] = 0.0;
      Bs[321] = 0.0;
      Bs[322] = 0.0;
      Bs[323] = 0.0;
      Bs[324] = 0.0;
      Bs[325] = -1.7830786258880667E-11;
      Bs[326] = 1.7830786258880667E-11;
      Bs[327] = 3.5661572517729275E-11;
      Bs[328] = -3.5661572517729275E-11;
      Bs[329] = -1.7830786258848608E-11;
      Bs[330] = 1.7830786258848608E-11;
      Bs[331] = 0.0;
      Bs[332] = 0.0;
      Bs[333] = 7.4208621094747382;
      Bs[334] = -0.51544761328856536;
      Bs[335] = -0.51544761328856536;
      Bs[336] = 0.00019860087161968924;
      Bs[337] = 0.00019860087162055056;
      Bs[338] = -0.00039720174323950723;
      Bs[339] = 0.0;
      Bs[340] = 0.0;
      Bs[341] = 0.0;
      Bs[342] = 0.0;
      Bs[343] = 0.0;
      Bs[344] = 0.0;
      Bs[345] = 0.0;
      Bs[346] = 0.0;
      Bs[347] = 0.0;
      Bs[348] = 0.0;
      Bs[349] = -1.7830786258843447E-11;
      Bs[350] = 1.7830786258843447E-11;
      Bs[351] = -1.7830786258920776E-11;
      Bs[352] = 1.7830786258920776E-11;
      Bs[353] = 3.5661572517764223E-11;
      Bs[354] = -3.5661572517764223E-11;
      Bs[355] = 0.0;
      Bs[356] = 0.0;
      Bs[357] = -0.51544761328856559;
      Bs[358] = 7.4208621094747427;
      Bs[359] = -0.515447613288566;
      Bs[360] = -0.0003972017432382971;
      Bs[361] = 0.00019860087161830176;
      Bs[362] = 0.00019860087161859273;
      Bs[363] = 0.0;
      Bs[364] = 0.0;
      Bs[365] = 0.0;
      Bs[366] = 0.0;
      Bs[367] = 0.0;
      Bs[368] = 0.0;
      Bs[369] = 0.0;
      Bs[370] = 0.0;
      Bs[371] = 0.0;
      Bs[372] = 0.0;
      Bs[373] = 3.5661572517591678E-11;
      Bs[374] = -3.5661572517591678E-11;
      Bs[375] = -1.7830786258782789E-11;
      Bs[376] = 1.7830786258782789E-11;
      Bs[377] = -1.7830786258808893E-11;
      Bs[378] = 1.7830786258808893E-11;
      Bs[379] = 0.0;
      Bs[380] = 0.0;
      Bs[381] = -0.51544761328856536;
      Bs[382] = -0.51544761328856525;
      Bs[383] = 7.4208621094747391;
      Bs[384] = 2.1886609091307321;
      Bs[385] = -2.1886609091307396;
      Bs[386] = -7.3267530588588332E-17;
      Bs[387] = 0.0;
      Bs[388] = 0.0;
      Bs[389] = 0.0;
      Bs[390] = 0.0;
      Bs[391] = 0.0;
      Bs[392] = 0.0;
      Bs[393] = 0.0;
      Bs[394] = 0.0;
      Bs[395] = 0.0;
      Bs[396] = 0.0;
      Bs[397] = -1.9650238463492784E-7;
      Bs[398] = 1.9650238463492784E-7;
      Bs[399] = 1.965023846349281E-7;
      Bs[400] = -1.965023846349281E-7;
      Bs[401] = -2.722140438625182E-22;
      Bs[402] = 2.722140438625182E-22;
      Bs[403] = 0.0;
      Bs[404] = 0.0;
      Bs[405] = 0.0028799381931014592;
      Bs[406] = 1.6124354859399419E-19;
      Bs[407] = -0.0028799381931014596;
      Bs[408] = 1.6558825407306325E-12;
      Bs[409] = 2.1886609091213165;
      Bs[410] = -2.1886609091258915;
      Bs[411] = 0.0;
      Bs[412] = 0.0;
      Bs[413] = 0.0;
      Bs[414] = 0.0;
      Bs[415] = 0.0;
      Bs[416] = 0.0;
      Bs[417] = 0.0;
      Bs[418] = 0.0;
      Bs[419] = 0.0;
      Bs[420] = 0.0;
      Bs[421] = -2.3599980159508313E-19;
      Bs[422] = 2.3599980159508313E-19;
      Bs[423] = -1.9650238463416961E-7;
      Bs[424] = 1.9650238463416961E-7;
      Bs[425] = 1.9650238463440561E-7;
      Bs[426] = -1.9650238463440561E-7;
      Bs[427] = 0.0;
      Bs[428] = 0.0;
      Bs[429] = -0.00287993819310146;
      Bs[430] = 0.0028799381931014613;
      Bs[431] = 1.6783958471455009E-19;
      Bs[432] = -2.1886609091286244;
      Bs[433] = -6.9752649900975711E-17;
      Bs[434] = 2.1886609091258937;
      Bs[435] = 0.0;
      Bs[436] = 0.0;
      Bs[437] = 0.0;
      Bs[438] = 0.0;
      Bs[439] = 0.0;
      Bs[440] = 0.0;
      Bs[441] = 0.0;
      Bs[442] = 0.0;
      Bs[443] = 0.0;
      Bs[444] = 0.0;
      Bs[445] = 1.9650238463465663E-7;
      Bs[446] = -1.9650238463465663E-7;
      Bs[447] = -8.1781743193057082E-20;
      Bs[448] = 8.1781743193057082E-20;
      Bs[449] = -1.9650238463457486E-7;
      Bs[450] = 1.9650238463457486E-7;
      Bs[451] = 0.0;
      Bs[452] = 0.0;
      Bs[453] = 1.144846922391657E-19;
      Bs[454] = -0.0028799381931014618;
      Bs[455] = 0.0028799381931014613;
      Bs[456] = -12.856960750434135;
      Bs[457] = 12.856960750470073;
      Bs[458] = 1.6067342690836219E-11;
      Bs[459] = 0.0;
      Bs[460] = 0.0;
      Bs[461] = 0.0;
      Bs[462] = 0.0;
      Bs[463] = 0.0;
      Bs[464] = 0.0;
      Bs[465] = 0.0;
      Bs[466] = 0.0;
      Bs[467] = 0.0;
      Bs[468] = 0.0;
      Bs[469] = 1.1543238315636564E-6;
      Bs[470] = -1.1543238315636564E-6;
      Bs[471] = -1.1543238315637706E-6;
      Bs[472] = 1.1543238315637706E-6;
      Bs[473] = 1.1411271306639956E-19;
      Bs[474] = -1.1411271306639956E-19;
      Bs[475] = 0.0;
      Bs[476] = 0.0;
      Bs[477] = 90.7238958483362;
      Bs[478] = 9.8322002140358478E-15;
      Bs[479] = -90.723895848336241;
      Bs[480] = -5.7371092957635561E-12;
      Bs[481] = -12.856960750486452;
      Bs[482] = 12.856960750437786;
      Bs[483] = 0.0;
      Bs[484] = 0.0;
      Bs[485] = 0.0;
      Bs[486] = 0.0;
      Bs[487] = 0.0;
      Bs[488] = 0.0;
      Bs[489] = 0.0;
      Bs[490] = 0.0;
      Bs[491] = 0.0;
      Bs[492] = 0.0;
      Bs[493] = -1.1132895703682332E-18;
      Bs[494] = 1.1132895703682332E-18;
      Bs[495] = 1.154323831565169E-6;
      Bs[496] = -1.154323831565169E-6;
      Bs[497] = -1.1543238315640558E-6;
      Bs[498] = 1.1543238315640558E-6;
      Bs[499] = 0.0;
      Bs[500] = 0.0;
      Bs[501] = -90.723895848336227;
      Bs[502] = 90.7238958483363;
      Bs[503] = -4.6576426504441113E-15;
      Bs[504] = 12.856960750417571;
      Bs[505] = -6.9642756938981512E-13;
      Bs[506] = -12.856960750430892;
      Bs[507] = 0.0;
      Bs[508] = 0.0;
      Bs[509] = 0.0;
      Bs[510] = 0.0;
      Bs[511] = 0.0;
      Bs[512] = 0.0;
      Bs[513] = 0.0;
      Bs[514] = 0.0;
      Bs[515] = 0.0;
      Bs[516] = 0.0;
      Bs[517] = -1.1543238315610321E-6;
      Bs[518] = 1.1543238315610321E-6;
      Bs[519] = -3.5670674991246723E-19;
      Bs[520] = 3.5670674991246723E-19;
      Bs[521] = 1.1543238315613889E-6;
      Bs[522] = -1.1543238315613889E-6;
      Bs[523] = 0.0;
      Bs[524] = 0.0;
      Bs[525] = -1.7763568394002505E-15;
      Bs[526] = -90.723895848336284;
      Bs[527] = 90.723895848336284;
      Bs[528] = -2.8259520727857769E-12;
      Bs[529] = 2.7114366681483019E-12;
      Bs[530] = -4.0106654423863111E-12;
      Bs[531] = 0.0;
      Bs[532] = 0.0;
      Bs[533] = 0.0;
      Bs[534] = 0.0;
      Bs[535] = 0.0;
      Bs[536] = 0.0;
      Bs[537] = 0.0;
      Bs[538] = 0.0;
      Bs[539] = 0.0;
      Bs[540] = 0.0;
      Bs[541] = 1.3026398720796875E-19;
      Bs[542] = -1.3026398720796875E-19;
      Bs[543] = -3.6689392688677759E-19;
      Bs[544] = 3.6689392688677759E-19;
      Bs[545] = 2.3662993967880886E-19;
      Bs[546] = -2.3662993967880886E-19;
      Bs[547] = 0.0;
      Bs[548] = 0.0;
      Bs[549] = -9.0328729587173649E-15;
      Bs[550] = 3.9022789309798385E-17;
      Bs[551] = 8.4956963007165621E-15;
      Bs[552] = -68947.241957919847;
      Bs[553] = 68947.241957985942;
      Bs[554] = 1.0324924783564313E-11;
      Bs[555] = 0.0;
      Bs[556] = 0.0;
      Bs[557] = 0.0;
      Bs[558] = 0.0;
      Bs[559] = 0.0;
      Bs[560] = 0.0;
      Bs[561] = 0.0;
      Bs[562] = 0.0;
      Bs[563] = 0.0;
      Bs[564] = 0.0;
      Bs[565] = 0.0061902222506079645;
      Bs[566] = -0.0061902222506079645;
      Bs[567] = -0.006190222250609943;
      Bs[568] = 0.006190222250609943;
      Bs[569] = 1.9780104602924559E-15;
      Bs[570] = -1.9780104602924559E-15;
      Bs[571] = 0.0;
      Bs[572] = 0.0;
      Bs[573] = 0.01694012792565032;
      Bs[574] = 3.93186976668021E-18;
      Bs[575] = -0.01694012792565033;
      Bs[576] = 1.8095929552650256E-11;
      Bs[577] = -68947.241957963925;
      Bs[578] = 68947.24195798584;
      Bs[579] = 0.0;
      Bs[580] = 0.0;
      Bs[581] = 0.0;
      Bs[582] = 0.0;
      Bs[583] = 0.0;
      Bs[584] = 0.0;
      Bs[585] = 0.0;
      Bs[586] = 0.0;
      Bs[587] = 0.0;
      Bs[588] = 0.0;
      Bs[589] = 6.5460524690502431E-16;
      Bs[590] = -6.5460524690502431E-16;
      Bs[591] = 0.0061902222506105995;
      Bs[592] = -0.0061902222506105995;
      Bs[593] = -0.0061902222506112544;
      Bs[594] = 0.0061902222506112544;
      Bs[595] = 0.0;
      Bs[596] = 0.0;
      Bs[597] = -0.016940127925650326;
      Bs[598] = 0.016940127925650337;
      Bs[599] = 8.8333656522023788E-19;
      Bs[600] = 68947.241957919832;
      Bs[601] = -5.7114189872985041E-12;
      Bs[602] = -68947.241957886945;
      Bs[603] = 0.0;
      Bs[604] = 0.0;
      Bs[605] = 0.0;
      Bs[606] = 0.0;
      Bs[607] = 0.0;
      Bs[608] = 0.0;
      Bs[609] = 0.0;
      Bs[610] = 0.0;
      Bs[611] = 0.0;
      Bs[612] = 0.0;
      Bs[613] = -0.006190222250605;
      Bs[614] = 0.006190222250605;
      Bs[615] = 9.8628993326203186E-16;
      Bs[616] = -9.8628993326203186E-16;
      Bs[617] = 0.0061902222506040137;
      Bs[618] = -0.0061902222506040137;
      Bs[619] = 0.0;
      Bs[620] = 0.0;
      Bs[621] = 2.4198360288954365E-19;
      Bs[622] = -0.01694012792565033;
      Bs[623] = 0.016940127925650333;
      Bs[624] = 6.4355701254132335E-12;
      Bs[625] = -6.4376042824419889E-12;
      Bs[626] = -1.4093469193128503E-19;
      Bs[627] = 0.0;
      Bs[628] = 0.0;
      Bs[629] = 0.0;
      Bs[630] = 0.0;
      Bs[631] = 0.0;
      Bs[632] = 0.0;
      Bs[633] = 0.0;
      Bs[634] = 0.0;
      Bs[635] = 0.0;
      Bs[636] = 0.0;
      Bs[637] = -5.7785932895605107E-19;
      Bs[638] = 5.7785932895605107E-19;
      Bs[639] = 5.7792019742547294E-19;
      Bs[640] = -5.7792019742547294E-19;
      Bs[641] = -6.086846942182349E-23;
      Bs[642] = 6.086846942182349E-23;
      Bs[643] = 0.0;
      Bs[644] = 0.0;
      Bs[645] = -2.6287865935834384E-18;
      Bs[646] = 9.2969840919454E-19;
      Bs[647] = 2.9975866472800142E-18;
      Bs[648] = 0.0;
      Bs[649] = 0.0;
      Bs[650] = 0.0;
      Bs[651] = -2736.5254573585835;
      Bs[652] = 547.30509147171676;
      Bs[653] = 2736.5254573585835;
      Bs[654] = 547.30509147171676;
      Bs[655] = -12053.897362855107;
      Bs[656] = -3013.4618473259;
      Bs[657] = -3013.4618473259;
      Bs[658] = 6026.9736682033081;
      Bs[659] = -3013.4618473258993;
      Bs[660] = -3013.4618473258993;
      Bs[661] = 0.0;
      Bs[662] = 0.0;
      Bs[663] = 0.0;
      Bs[664] = 0.0;
      Bs[665] = 0.0;
      Bs[666] = 0.0;
      Bs[667] = 9040.4355155292069;
      Bs[668] = 6.6178456227328588E-15;
      Bs[669] = 0.0;
      Bs[670] = 0.0;
      Bs[671] = 0.0;
      Bs[672] = 0.0;
      Bs[673] = 0.0;
      Bs[674] = 0.0;
      Bs[675] = -2736.5254573585835;
      Bs[676] = 547.30509147171665;
      Bs[677] = 2736.5254573585835;
      Bs[678] = 547.30509147171665;
      Bs[679] = -3013.4618473258993;
      Bs[680] = -12053.897362855103;
      Bs[681] = -3013.4618473258993;
      Bs[682] = -3013.4618473258993;
      Bs[683] = 6026.9736682033044;
      Bs[684] = -3013.4618473258993;
      Bs[685] = 0.0;
      Bs[686] = 0.0;
      Bs[687] = 0.0;
      Bs[688] = 0.0;
      Bs[689] = 0.0;
      Bs[690] = 0.0;
      Bs[691] = 9040.4355155292051;
      Bs[692] = 9040.4355155292051;
      Bs[693] = 0.0;
      Bs[694] = 0.0;
      Bs[695] = 0.0;
      Bs[696] = 0.0;
      Bs[697] = 0.0;
      Bs[698] = 0.0;
      Bs[699] = -2736.5254573585835;
      Bs[700] = 547.30509147171676;
      Bs[701] = 2736.5254573585835;
      Bs[702] = 547.30509147171676;
      Bs[703] = -3013.4618473258993;
      Bs[704] = -3013.4618473258993;
      Bs[705] = -12053.897362855103;
      Bs[706] = -3013.4618473258993;
      Bs[707] = -3013.4618473258993;
      Bs[708] = 6026.9736682033044;
      Bs[709] = 0.0;
      Bs[710] = 0.0;
      Bs[711] = 0.0;
      Bs[712] = 0.0;
      Bs[713] = 0.0;
      Bs[714] = 0.0;
      Bs[715] = 4.8414887833326083E-15;
      Bs[716] = 9040.4355155292051;
      Bs[717] = 0.0;
      Bs[718] = 0.0;
      Bs[719] = 0.0;
      Bs[720] = 0.0;
      Bs[721] = 0.0;
      Bs[722] = 0.0;
      Bs[723] = -2736.5254573585835;
      Bs[724] = 547.30509147171665;
      Bs[725] = 2736.5254573585835;
      Bs[726] = 547.30509147171665;
      Bs[727] = 6026.9736682033063;
      Bs[728] = -3013.4618473258997;
      Bs[729] = -3013.4618473258997;
      Bs[730] = -12053.897362855105;
      Bs[731] = -3013.4618473258988;
      Bs[732] = -3013.4618473258988;
      Bs[733] = 0.0;
      Bs[734] = 0.0;
      Bs[735] = 0.0;
      Bs[736] = 0.0;
      Bs[737] = 0.0;
      Bs[738] = 0.0;
      Bs[739] = -9040.4355155292051;
      Bs[740] = 8.210227336302985E-15;
      Bs[741] = 0.0;
      Bs[742] = 0.0;
      Bs[743] = 0.0;
      Bs[744] = 0.0;
      Bs[745] = 0.0;
      Bs[746] = 0.0;
      Bs[747] = -2736.5254573585835;
      Bs[748] = 547.30509147171676;
      Bs[749] = 2736.5254573585835;
      Bs[750] = 547.30509147171676;
      Bs[751] = -3013.4618473258993;
      Bs[752] = 6026.9736682033054;
      Bs[753] = -3013.4618473258993;
      Bs[754] = -3013.4618473258993;
      Bs[755] = -12053.897362855105;
      Bs[756] = -3013.4618473258993;
      Bs[757] = 0.0;
      Bs[758] = 0.0;
      Bs[759] = 0.0;
      Bs[760] = 0.0;
      Bs[761] = 0.0;
      Bs[762] = 0.0;
      Bs[763] = -9040.4355155292051;
      Bs[764] = -9040.4355155292033;
      Bs[765] = 0.0;
      Bs[766] = 0.0;
      Bs[767] = 0.0;
      Bs[768] = 0.0;
      Bs[769] = 0.0;
      Bs[770] = 0.0;
      Bs[771] = -2736.5254573585839;
      Bs[772] = 547.30509147171688;
      Bs[773] = 2736.5254573585839;
      Bs[774] = 547.30509147171688;
      Bs[775] = -3013.4618473258997;
      Bs[776] = -3013.4618473258997;
      Bs[777] = 6026.9736682033054;
      Bs[778] = -3013.4618473258997;
      Bs[779] = -3013.4618473258997;
      Bs[780] = -12053.897362855105;
      Bs[781] = 0.0;
      Bs[782] = 0.0;
      Bs[783] = 0.0;
      Bs[784] = 0.0;
      Bs[785] = 0.0;
      Bs[786] = 0.0;
      Bs[787] = 6.1721824383383966E-16;
      Bs[788] = -9040.4355155292033;
      Bs[789] = 0.0;
      Bs[790] = 0.0;
      Bs[791] = 0.0;
      Bs[792] = 11722.070250397561;
      Bs[793] = 0.811741926513632;
      Bs[794] = 0.81174192651383215;
      Bs[795] = 0.0;
      Bs[796] = 0.0;
      Bs[797] = 0.0;
      Bs[798] = 0.0;
      Bs[799] = 0.0;
      Bs[800] = 0.0;
      Bs[801] = 0.0;
      Bs[802] = 0.0;
      Bs[803] = 0.0;
      Bs[804] = 0.0;
      Bs[805] = -0.00070157213122970736;
      Bs[806] = 0.00070157213122970736;
      Bs[807] = 0.00035078606561485373;
      Bs[808] = -0.00035078606561485373;
      Bs[809] = 0.00035078606561485368;
      Bs[810] = -0.00035078606561485368;
      Bs[811] = 0.0;
      Bs[812] = 0.0;
      Bs[813] = 1.994425257316909E-7;
      Bs[814] = 1.9944252573169109E-7;
      Bs[815] = -3.9888505146338217E-7;
      Bs[816] = -999966.86827385728;
      Bs[817] = -0.0080831607982809472;
      Bs[818] = -0.0080831607982729328;
      Bs[819] = 0.0;
      Bs[820] = 0.0;
      Bs[821] = 0.0;
      Bs[822] = 0.0;
      Bs[823] = 0.0;
      Bs[824] = 0.0;
      Bs[825] = 0.0;
      Bs[826] = 0.0;
      Bs[827] = 0.0;
      Bs[828] = 0.0;
      Bs[829] = 0.059852692503629269;
      Bs[830] = -0.059852692503629269;
      Bs[831] = -0.029926346251814648;
      Bs[832] = 0.029926346251814648;
      Bs[833] = -0.029926346251814628;
      Bs[834] = 0.029926346251814628;
      Bs[835] = 0.0;
      Bs[836] = 0.0;
      Bs[837] = -1.9860087161920744E-9;
      Bs[838] = -1.9860087161920761E-9;
      Bs[839] = 3.9720174323841513E-9;
      Bs[840] = 0.8117419265124951;
      Bs[841] = 11722.070250401342;
      Bs[842] = 0.81174192651383215;
      Bs[843] = 0.0;
      Bs[844] = 0.0;
      Bs[845] = 0.0;
      Bs[846] = 0.0;
      Bs[847] = 0.0;
      Bs[848] = 0.0;
      Bs[849] = 0.0;
      Bs[850] = 0.0;
      Bs[851] = 0.0;
      Bs[852] = 0.0;
      Bs[853] = 0.00035078606561496703;
      Bs[854] = -0.00035078606561496703;
      Bs[855] = -0.000701572131229934;
      Bs[856] = 0.000701572131229934;
      Bs[857] = 0.00035078606561496692;
      Bs[858] = -0.00035078606561496692;
      Bs[859] = 0.0;
      Bs[860] = 0.0;
      Bs[861] = -3.98885051463382E-7;
      Bs[862] = 1.9944252573169106E-7;
      Bs[863] = 1.9944252573169111E-7;
      Bs[864] = -0.008083160798281603;
      Bs[865] = -999966.86827385752;
      Bs[866] = -0.0080831607982729328;
      Bs[867] = 0.0;
      Bs[868] = 0.0;
      Bs[869] = 0.0;
      Bs[870] = 0.0;
      Bs[871] = 0.0;
      Bs[872] = 0.0;
      Bs[873] = 0.0;
      Bs[874] = 0.0;
      Bs[875] = 0.0;
      Bs[876] = 0.0;
      Bs[877] = -0.029926346251814652;
      Bs[878] = 0.029926346251814652;
      Bs[879] = 0.059852692503629311;
      Bs[880] = -0.059852692503629311;
      Bs[881] = -0.029926346251814652;
      Bs[882] = 0.029926346251814652;
      Bs[883] = 0.0;
      Bs[884] = 0.0;
      Bs[885] = 3.9720174323841504E-9;
      Bs[886] = -1.9860087161920756E-9;
      Bs[887] = -1.9860087161920756E-9;
      Bs[888] = 0.81174192651351373;
      Bs[889] = 0.81174192651319421;
      Bs[890] = 11722.070250399855;
      Bs[891] = 0.0;
      Bs[892] = 0.0;
      Bs[893] = 0.0;
      Bs[894] = 0.0;
      Bs[895] = 0.0;
      Bs[896] = 0.0;
      Bs[897] = 0.0;
      Bs[898] = 0.0;
      Bs[899] = 0.0;
      Bs[900] = 0.0;
      Bs[901] = 0.00035078606561492236;
      Bs[902] = -0.00035078606561492236;
      Bs[903] = 0.00035078606561492242;
      Bs[904] = -0.00035078606561492242;
      Bs[905] = -0.00070157213122984472;
      Bs[906] = 0.00070157213122984472;
      Bs[907] = 0.0;
      Bs[908] = 0.0;
      Bs[909] = 1.9944252573169103E-7;
      Bs[910] = -3.9888505146338217E-7;
      Bs[911] = 1.9944252573169109E-7;
      Bs[912] = -0.00808316079828827;
      Bs[913] = -0.0080831607983035958;
      Bs[914] = -999966.86827385728;
      Bs[915] = 0.0;
      Bs[916] = 0.0;
      Bs[917] = 0.0;
      Bs[918] = 0.0;
      Bs[919] = 0.0;
      Bs[920] = 0.0;
      Bs[921] = 0.0;
      Bs[922] = 0.0;
      Bs[923] = 0.0;
      Bs[924] = 0.0;
      Bs[925] = -0.029926346251814631;
      Bs[926] = 0.029926346251814631;
      Bs[927] = -0.029926346251814631;
      Bs[928] = 0.029926346251814631;
      Bs[929] = 0.059852692503629262;
      Bs[930] = -0.059852692503629262;
      Bs[931] = 0.0;
      Bs[932] = 0.0;
      Bs[933] = -1.9860087161920752E-9;
      Bs[934] = 3.9720174323841513E-9;
      Bs[935] = -1.9860087161920752E-9;
      Bs[936] = -0.00083897388535522643;
      Bs[937] = 0.0016779477707082502;
      Bs[938] = -0.00083897388535441361;
      Bs[939] = 0.0;
      Bs[940] = 0.0;
      Bs[941] = 0.0;
      Bs[942] = 0.0;
      Bs[943] = 0.0;
      Bs[944] = 0.0;
      Bs[945] = 0.0;
      Bs[946] = 0.0;
      Bs[947] = 0.0;
      Bs[948] = 0.0;
      Bs[949] = 7.5324765216425754E-11;
      Bs[950] = -7.5324765216425754E-11;
      Bs[951] = -1.5064953043277854E-10;
      Bs[952] = 1.5064953043277854E-10;
      Bs[953] = 7.5324765216352782E-11;
      Bs[954] = -7.5324765216352782E-11;
      Bs[955] = 0.0;
      Bs[956] = 0.0;
      Bs[957] = -14.429668901488018;
      Bs[958] = 0.99952393278284035;
      Bs[959] = 0.99952393278284057;
      Bs[960] = -3.9244213056846263E-8;
      Bs[961] = 7.84884261135924E-8;
      Bs[962] = -3.9244213056587706E-8;
      Bs[963] = 0.0;
      Bs[964] = 0.0;
      Bs[965] = 0.0;
      Bs[966] = 0.0;
      Bs[967] = 0.0;
      Bs[968] = 0.0;
      Bs[969] = 0.0;
      Bs[970] = 0.0;
      Bs[971] = 0.0;
      Bs[972] = 0.0;
      Bs[973] = 3.5234244905779534E-15;
      Bs[974] = -3.5234244905779534E-15;
      Bs[975] = -7.0468489811326921E-15;
      Bs[976] = 7.0468489811326921E-15;
      Bs[977] = 3.5234244905547396E-15;
      Bs[978] = -3.5234244905547396E-15;
      Bs[979] = 0.0;
      Bs[980] = 0.0;
      Bs[981] = 0.0008343452446523109;
      Bs[982] = -5.8326992931449861E-5;
      Bs[983] = -5.8326992931449847E-5;
      Bs[984] = -0.00083897388535550258;
      Bs[985] = -0.00083897388535628158;
      Bs[986] = 0.0016779477707016444;
      Bs[987] = 0.0;
      Bs[988] = 0.0;
      Bs[989] = 0.0;
      Bs[990] = 0.0;
      Bs[991] = 0.0;
      Bs[992] = 0.0;
      Bs[993] = 0.0;
      Bs[994] = 0.0;
      Bs[995] = 0.0;
      Bs[996] = 0.0;
      Bs[997] = 7.5324765216188651E-11;
      Bs[998] = -7.5324765216188651E-11;
      Bs[999] = 7.53247652162586E-11;
      Bs[1000] = -7.53247652162586E-11;
      Bs[1001] = -1.5064953043244725E-10;
      Bs[1002] = 1.5064953043244725E-10;
      Bs[1003] = 0.0;
      Bs[1004] = 0.0;
      Bs[1005] = 0.99952393278284113;
      Bs[1006] = -14.429668901488029;
      Bs[1007] = 0.99952393278284191;
      Bs[1008] = -3.9244213056902485E-8;
      Bs[1009] = -3.9244213056818112E-8;
      Bs[1010] = 7.8488426113645727E-8;
      Bs[1011] = 0.0;
      Bs[1012] = 0.0;
      Bs[1013] = 0.0;
      Bs[1014] = 0.0;
      Bs[1015] = 0.0;
      Bs[1016] = 0.0;
      Bs[1017] = 0.0;
      Bs[1018] = 0.0;
      Bs[1019] = 0.0;
      Bs[1020] = 0.0;
      Bs[1021] = 3.5234244905760183E-15;
      Bs[1022] = -3.5234244905760183E-15;
      Bs[1023] = 3.5234244905684429E-15;
      Bs[1024] = -3.5234244905684429E-15;
      Bs[1025] = -7.0468489811444612E-15;
      Bs[1026] = 7.0468489811444612E-15;
      Bs[1027] = 0.0;
      Bs[1028] = 0.0;
      Bs[1029] = -5.8326992931449867E-5;
      Bs[1030] = 0.00083434524465231111;
      Bs[1031] = -5.8326992931449915E-5;
      Bs[1032] = 0.0016779477707029678;
      Bs[1033] = -0.0008389738853556625;
      Bs[1034] = -0.000838973885352657;
      Bs[1035] = 0.0;
      Bs[1036] = 0.0;
      Bs[1037] = 0.0;
      Bs[1038] = 0.0;
      Bs[1039] = 0.0;
      Bs[1040] = 0.0;
      Bs[1041] = 0.0;
      Bs[1042] = 0.0;
      Bs[1043] = 0.0;
      Bs[1044] = 0.0;
      Bs[1045] = -1.5064953043242282E-10;
      Bs[1046] = 1.5064953043242282E-10;
      Bs[1047] = 7.5324765216346384E-11;
      Bs[1048] = -7.5324765216346384E-11;
      Bs[1049] = 7.5324765216076439E-11;
      Bs[1050] = -7.5324765216076439E-11;
      Bs[1051] = 0.0;
      Bs[1052] = 0.0;
      Bs[1053] = 0.9995239327828408;
      Bs[1054] = 0.9995239327828408;
      Bs[1055] = -14.429668901488025;
      Bs[1056] = 7.8488426113706329E-8;
      Bs[1057] = -3.9244213056853165E-8;
      Bs[1058] = -3.9244213056763214E-8;
      Bs[1059] = 0.0;
      Bs[1060] = 0.0;
      Bs[1061] = 0.0;
      Bs[1062] = 0.0;
      Bs[1063] = 0.0;
      Bs[1064] = 0.0;
      Bs[1065] = 0.0;
      Bs[1066] = 0.0;
      Bs[1067] = 0.0;
      Bs[1068] = 0.0;
      Bs[1069] = -7.04684898114497E-15;
      Bs[1070] = 7.04684898114497E-15;
      Bs[1071] = 3.5234244905765243E-15;
      Bs[1072] = -3.5234244905765243E-15;
      Bs[1073] = 3.5234244905684452E-15;
      Bs[1074] = -3.5234244905684452E-15;
      Bs[1075] = 0.0;
      Bs[1076] = 0.0;
      Bs[1077] = -5.8326992931449867E-5;
      Bs[1078] = -5.8326992931449861E-5;
      Bs[1079] = 0.00083434524465231111;
      Cs[0] = 0.0;
      Cs[1] = 0.0;
      Cs[2] = -2.1259022226465764E-12;
      Cs[3] = -2.125902222646576E-12;
      Cs[4] = -1.999933616842329E-6;
      Cs[5] = -7.601901410019118E-14;
      Cs[6] = -7.60190141001751E-14;
      Cs[7] = -3.4419992602974414E-12;
      Cs[8] = 2.0550664363147613E-11;
      Cs[9] = 1.7209996301487205E-12;
      Cs[10] = 2.0550664363147606E-11;
      Cs[11] = 6.0218088549010634E-27;
      Cs[12] = 1.7209996301487211E-12;
      Cs[13] = 3.1893457731583761E-11;
      Cs[14] = 3.1893457731583787E-11;
      Cs[15] = -6.378691546316756E-11;
      Cs[16] = 1.3785363651718822E-7;
      Cs[17] = -1.0618921845918996E-23;
      Cs[18] = -1.3785363651718825E-7;
      Cs[19] = 4.3706351924459388E-12;
      Cs[20] = 1.0097419586828951E-27;
      Cs[21] = -4.3706351924459371E-12;
      Cs[22] = 1.4901813611218238E-28;
      Cs[23] = 2.3438139695106804E-8;
      Cs[24] = -3.105770170849696E-24;
      Cs[25] = -2.3438139695106811E-8;
      Cs[26] = 3.0969451882490321E-24;
      Cs[27] = 0.0;
      Cs[28] = 0.0;
      Cs[29] = 0.0;
      Cs[30] = 0.0;
      Cs[31] = 0.0;
      Cs[32] = 0.0;
      Cs[33] = 0.00023320317795248353;
      Cs[34] = -3.5079239765565253E-6;
      Cs[35] = 1.614909251497046E-8;
      Cs[36] = -2.8356053055635841E-14;
      Cs[37] = 1.6149092514970456E-8;
      Cs[38] = -2.8356053055607732E-14;
      Cs[39] = -1.4204983082856465E-11;
      Cs[40] = -3.1894011042063221E-11;
      Cs[41] = -1.4204983082856474E-11;
      Cs[42] = -3.1894011042063234E-11;
      Cs[43] = 2.8409966165712942E-11;
      Cs[44] = 6.3788022084126429E-11;
      Cs[45] = 0.0;
      Cs[46] = 0.0;
      Cs[47] = 2.1259022226465772E-12;
      Cs[48] = 2.3999207555933726E-28;
      Cs[49] = -7.6019014100192493E-14;
      Cs[50] = -1.99993361684233E-6;
      Cs[51] = -7.6019014100175138E-14;
      Cs[52] = 1.7209996301487207E-12;
      Cs[53] = -2.0550664363147606E-11;
      Cs[54] = -3.4419992602974422E-12;
      Cs[55] = -4.17153560170518E-27;
      Cs[56] = 2.0550664363147616E-11;
      Cs[57] = 1.7209996301487211E-12;
      Cs[58] = -6.378691546316756E-11;
      Cs[59] = 3.1893457731583787E-11;
      Cs[60] = 3.1893457731583787E-11;
      Cs[61] = -1.3785363651718817E-7;
      Cs[62] = 1.3785363651718835E-7;
      Cs[63] = 7.4187502738565969E-23;
      Cs[64] = -4.3706351924459371E-12;
      Cs[65] = 4.3706351924459388E-12;
      Cs[66] = 2.8272774843121061E-28;
      Cs[67] = 9.5921403103392935E-28;
      Cs[68] = -2.3438139695106807E-8;
      Cs[69] = 2.3438139695106804E-8;
      Cs[70] = -1.531940110035341E-24;
      Cs[71] = 3.3638922412245672E-23;
      Cs[72] = 0.0;
      Cs[73] = 0.0;
      Cs[74] = 0.0;
      Cs[75] = 0.0;
      Cs[76] = 0.0;
      Cs[77] = 0.0;
      Cs[78] = 1.6149092514970466E-8;
      Cs[79] = -2.8356053055638144E-14;
      Cs[80] = 0.00023320317795248356;
      Cs[81] = -3.5079239765565257E-6;
      Cs[82] = 1.6149092514970456E-8;
      Cs[83] = -2.8356053055607732E-14;
      Cs[84] = 2.8409966165712942E-11;
      Cs[85] = 6.3788022084126481E-11;
      Cs[86] = -1.4204983082856474E-11;
      Cs[87] = -3.1894011042063227E-11;
      Cs[88] = -1.4204983082856474E-11;
      Cs[89] = -3.1894011042063227E-11;
      Cs[90] = 0.0;
      Cs[91] = 0.0;
      Cs[92] = 3.2581919834171059E-32;
      Cs[93] = 2.1259022226465756E-12;
      Cs[94] = -7.60190141002058E-14;
      Cs[95] = -7.6019014100236442E-14;
      Cs[96] = -1.9999336168423295E-6;
      Cs[97] = 1.7209996301487197E-12;
      Cs[98] = -6.1411445587103495E-27;
      Cs[99] = 1.7209996301487207E-12;
      Cs[100] = -2.0550664363147593E-11;
      Cs[101] = -2.0550664363147619E-11;
      Cs[102] = -3.441999260297441E-12;
      Cs[103] = 3.189345773158378E-11;
      Cs[104] = -6.3786915463167573E-11;
      Cs[105] = 3.1893457731583767E-11;
      Cs[106] = -4.6339945379143055E-23;
      Cs[107] = -1.3785363651718833E-7;
      Cs[108] = 1.3785363651718811E-7;
      Cs[109] = 2.3505764316289037E-28;
      Cs[110] = -4.3706351924459379E-12;
      Cs[111] = 4.3706351924459355E-12;
      Cs[112] = -2.58371768130869E-27;
      Cs[113] = -9.4915744116192144E-28;
      Cs[114] = -2.3438139695106811E-8;
      Cs[115] = 2.3438139695106794E-8;
      Cs[116] = -1.8735641989439513E-23;
      Cs[117] = 0.0;
      Cs[118] = 0.0;
      Cs[119] = 0.0;
      Cs[120] = 0.0;
      Cs[121] = 0.0;
      Cs[122] = 0.0;
      Cs[123] = 1.614909251497046E-8;
      Cs[124] = -2.8356053055661529E-14;
      Cs[125] = 1.6149092514970466E-8;
      Cs[126] = -2.8356053055715295E-14;
      Cs[127] = 0.00023320317795248353;
      Cs[128] = -3.5079239765565253E-6;
      Cs[129] = -1.4204983082856471E-11;
      Cs[130] = -3.1894011042063234E-11;
      Cs[131] = 2.8409966165712949E-11;
      Cs[132] = 6.3788022084126455E-11;
      Cs[133] = -1.4204983082856466E-11;
      Cs[134] = -3.18940110420632E-11;
      Cs[135] = 0.0;
      Cs[136] = -0.00021794285528763531;
      Cs[137] = 0.0;
      Cs[138] = 0.0;
      Cs[139] = 0.0;
      Cs[140] = 0.0;
      Cs[141] = 0.0;
      Cs[142] = 0.0;
      Cs[143] = 0.0;
      Cs[144] = 0.0;
      Cs[145] = 0.0;
      Cs[146] = 0.0;
      Cs[147] = 0.0;
      Cs[148] = 0.0;
      Cs[149] = 0.0;
      Cs[150] = 0.0;
      Cs[151] = 0.0;
      Cs[152] = 0.0;
      Cs[153] = 0.0;
      Cs[154] = 0.0;
      Cs[155] = 0.0;
      Cs[156] = 0.0;
      Cs[157] = 0.0;
      Cs[158] = 0.0;
      Cs[159] = 0.0;
      Cs[160] = 0.0;
      Cs[161] = 0.0;
      Cs[162] = 3.026984101217156E-7;
      Cs[163] = 3.0269841012171581E-7;
      Cs[164] = 3.0269841012171576E-7;
      Cs[165] = 3.0269841012171554E-7;
      Cs[166] = 3.026984101217157E-7;
      Cs[167] = 3.0269841012171576E-7;
      Cs[168] = 0.0;
      Cs[169] = 0.0;
      Cs[170] = 0.0;
      Cs[171] = 0.0;
      Cs[172] = 0.0;
      Cs[173] = 0.0;
      Cs[174] = 0.0;
      Cs[175] = 0.0;
      Cs[176] = 0.0;
      Cs[177] = 0.0;
      Cs[178] = 0.0;
      Cs[179] = 0.0;
      Cs[180] = 0.0;
      Cs[181] = 4.3588571057527068E-5;
      Cs[182] = 0.0;
      Cs[183] = 0.0;
      Cs[184] = 0.0;
      Cs[185] = 0.0;
      Cs[186] = 0.0;
      Cs[187] = 0.0;
      Cs[188] = 0.0;
      Cs[189] = 0.0;
      Cs[190] = 0.0;
      Cs[191] = 0.0;
      Cs[192] = 0.0;
      Cs[193] = 0.0;
      Cs[194] = 0.0;
      Cs[195] = 0.0;
      Cs[196] = 0.0;
      Cs[197] = 0.0;
      Cs[198] = 0.0;
      Cs[199] = 0.0;
      Cs[200] = 0.0;
      Cs[201] = 0.0;
      Cs[202] = 0.0;
      Cs[203] = 0.0;
      Cs[204] = 0.0;
      Cs[205] = 0.0;
      Cs[206] = 0.0;
      Cs[207] = -6.0539682024343114E-8;
      Cs[208] = -6.0539682024343167E-8;
      Cs[209] = -6.0539682024343154E-8;
      Cs[210] = -6.0539682024343114E-8;
      Cs[211] = -6.0539682024343141E-8;
      Cs[212] = -6.0539682024343154E-8;
      Cs[213] = 0.0;
      Cs[214] = 0.0;
      Cs[215] = 0.0;
      Cs[216] = 0.0;
      Cs[217] = 0.0;
      Cs[218] = 0.0;
      Cs[219] = 0.0;
      Cs[220] = 0.0;
      Cs[221] = 0.0;
      Cs[222] = 0.0;
      Cs[223] = 0.0;
      Cs[224] = 0.0;
      Cs[225] = 0.0;
      Cs[226] = 0.00021794285528763531;
      Cs[227] = 0.0;
      Cs[228] = 0.0;
      Cs[229] = 0.0;
      Cs[230] = 0.0;
      Cs[231] = 0.0;
      Cs[232] = 0.0;
      Cs[233] = 0.0;
      Cs[234] = 0.0;
      Cs[235] = 0.0;
      Cs[236] = 0.0;
      Cs[237] = 0.0;
      Cs[238] = 0.0;
      Cs[239] = 0.0;
      Cs[240] = 0.0;
      Cs[241] = 0.0;
      Cs[242] = 0.0;
      Cs[243] = 0.0;
      Cs[244] = 0.0;
      Cs[245] = 0.0;
      Cs[246] = 0.0;
      Cs[247] = 0.0;
      Cs[248] = 0.0;
      Cs[249] = 0.0;
      Cs[250] = 0.0;
      Cs[251] = 0.0;
      Cs[252] = -3.026984101217156E-7;
      Cs[253] = -3.0269841012171581E-7;
      Cs[254] = -3.0269841012171576E-7;
      Cs[255] = -3.0269841012171554E-7;
      Cs[256] = -3.026984101217157E-7;
      Cs[257] = -3.0269841012171576E-7;
      Cs[258] = 0.0;
      Cs[259] = 0.0;
      Cs[260] = 0.0;
      Cs[261] = 0.0;
      Cs[262] = 0.0;
      Cs[263] = 0.0;
      Cs[264] = 0.0;
      Cs[265] = 0.0;
      Cs[266] = 0.0;
      Cs[267] = 0.0;
      Cs[268] = 0.0;
      Cs[269] = 0.0;
      Cs[270] = -1.9999999999666666E-6;
      Cs[271] = 4.3588571057527068E-5;
      Cs[272] = 0.0;
      Cs[273] = 0.0;
      Cs[274] = 0.0;
      Cs[275] = 0.0;
      Cs[276] = 0.0;
      Cs[277] = 0.0;
      Cs[278] = 0.0;
      Cs[279] = 0.0;
      Cs[280] = 0.0;
      Cs[281] = 0.0;
      Cs[282] = 0.0;
      Cs[283] = 0.0;
      Cs[284] = 0.0;
      Cs[285] = 0.0;
      Cs[286] = 0.0;
      Cs[287] = 0.0;
      Cs[288] = 0.0;
      Cs[289] = 0.0;
      Cs[290] = 0.0;
      Cs[291] = 0.0;
      Cs[292] = 0.0;
      Cs[293] = 0.0;
      Cs[294] = 0.0;
      Cs[295] = 0.0;
      Cs[296] = 0.0;
      Cs[297] = -6.0539682024343114E-8;
      Cs[298] = -6.0539682024343167E-8;
      Cs[299] = -6.0539682024343154E-8;
      Cs[300] = -6.0539682024343114E-8;
      Cs[301] = -6.0539682024343141E-8;
      Cs[302] = -6.0539682024343154E-8;
      Cs[303] = 0.0;
      Cs[304] = 0.0;
      Cs[305] = 0.0;
      Cs[306] = 0.0;
      Cs[307] = 0.0;
      Cs[308] = 0.0;
      Cs[309] = 0.0;
      Cs[310] = 0.0;
      Cs[311] = 0.0;
      Cs[312] = 0.0;
      Cs[313] = 0.0;
      Cs[314] = 0.0;
      Cs[315] = 0.0;
      Cs[316] = 1.3266675417274417E-9;
      Cs[317] = 0.0;
      Cs[318] = 0.0;
      Cs[319] = 0.0;
      Cs[320] = 0.0;
      Cs[321] = 0.0;
      Cs[322] = 0.0;
      Cs[323] = 0.0;
      Cs[324] = 0.0;
      Cs[325] = 0.0;
      Cs[326] = 0.0;
      Cs[327] = 0.0;
      Cs[328] = 0.0;
      Cs[329] = 0.0;
      Cs[330] = 0.0;
      Cs[331] = 0.0;
      Cs[332] = 0.0;
      Cs[333] = 0.0;
      Cs[334] = 0.0;
      Cs[335] = 0.0;
      Cs[336] = 0.0;
      Cs[337] = 0.0;
      Cs[338] = 0.0;
      Cs[339] = 0.0;
      Cs[340] = 0.0;
      Cs[341] = 0.0;
      Cs[342] = 1.3333314907395251E-6;
      Cs[343] = 3.3333149073952561E-7;
      Cs[344] = 3.3333149073952561E-7;
      Cs[345] = -6.66668509260475E-7;
      Cs[346] = 3.3333149073952561E-7;
      Cs[347] = 3.3333149073952561E-7;
      Cs[348] = 0.0;
      Cs[349] = 0.0;
      Cs[350] = 0.0;
      Cs[351] = 0.0;
      Cs[352] = 0.0;
      Cs[353] = 0.0;
      Cs[354] = 0.0;
      Cs[355] = 0.0;
      Cs[356] = 0.0;
      Cs[357] = 0.0;
      Cs[358] = 0.0;
      Cs[359] = 0.0;
      Cs[360] = 0.0;
      Cs[361] = 1.3266675417274413E-9;
      Cs[362] = 0.0;
      Cs[363] = 0.0;
      Cs[364] = 0.0;
      Cs[365] = 0.0;
      Cs[366] = 0.0;
      Cs[367] = 0.0;
      Cs[368] = 0.0;
      Cs[369] = 0.0;
      Cs[370] = 0.0;
      Cs[371] = 0.0;
      Cs[372] = 0.0;
      Cs[373] = 0.0;
      Cs[374] = 0.0;
      Cs[375] = 0.0;
      Cs[376] = 0.0;
      Cs[377] = 0.0;
      Cs[378] = 0.0;
      Cs[379] = 0.0;
      Cs[380] = 0.0;
      Cs[381] = 0.0;
      Cs[382] = 0.0;
      Cs[383] = 0.0;
      Cs[384] = 0.0;
      Cs[385] = 0.0;
      Cs[386] = 0.0;
      Cs[387] = 3.333314907395255E-7;
      Cs[388] = 1.3333314907395253E-6;
      Cs[389] = 3.3333149073952534E-7;
      Cs[390] = 3.3333149073952539E-7;
      Cs[391] = -6.6666850926047451E-7;
      Cs[392] = 3.3333149073952534E-7;
      Cs[393] = 0.0;
      Cs[394] = 0.0;
      Cs[395] = 0.0;
      Cs[396] = 0.0;
      Cs[397] = 0.0;
      Cs[398] = 0.0;
      Cs[399] = 0.0;
      Cs[400] = 0.0;
      Cs[401] = 0.0;
      Cs[402] = 0.0;
      Cs[403] = 0.0;
      Cs[404] = 0.0;
      Cs[405] = 0.0;
      Cs[406] = 1.3266675417274415E-9;
      Cs[407] = 0.0;
      Cs[408] = 0.0;
      Cs[409] = 0.0;
      Cs[410] = 0.0;
      Cs[411] = 0.0;
      Cs[412] = 0.0;
      Cs[413] = 0.0;
      Cs[414] = 0.0;
      Cs[415] = 0.0;
      Cs[416] = 0.0;
      Cs[417] = 0.0;
      Cs[418] = 0.0;
      Cs[419] = 0.0;
      Cs[420] = 0.0;
      Cs[421] = 0.0;
      Cs[422] = 0.0;
      Cs[423] = 0.0;
      Cs[424] = 0.0;
      Cs[425] = 0.0;
      Cs[426] = 0.0;
      Cs[427] = 0.0;
      Cs[428] = 0.0;
      Cs[429] = 0.0;
      Cs[430] = 0.0;
      Cs[431] = 0.0;
      Cs[432] = 3.3333149073952555E-7;
      Cs[433] = 3.3333149073952545E-7;
      Cs[434] = 1.3333314907395251E-6;
      Cs[435] = 3.3333149073952545E-7;
      Cs[436] = 3.3333149073952539E-7;
      Cs[437] = -6.6666850926047461E-7;
      Cs[438] = 0.0;
      Cs[439] = 0.0;
      Cs[440] = 0.0;
      Cs[441] = 0.0;
      Cs[442] = 0.0;
      Cs[443] = 0.0;
      Cs[444] = 0.0;
      Cs[445] = 0.0;
      Cs[446] = 0.0;
      Cs[447] = 0.0;
      Cs[448] = 0.0;
      Cs[449] = 0.0;
      Cs[450] = 0.0;
      Cs[451] = 1.3266675417274413E-9;
      Cs[452] = 0.0;
      Cs[453] = 0.0;
      Cs[454] = 0.0;
      Cs[455] = 0.0;
      Cs[456] = 0.0;
      Cs[457] = 0.0;
      Cs[458] = 0.0;
      Cs[459] = 0.0;
      Cs[460] = 0.0;
      Cs[461] = 0.0;
      Cs[462] = 0.0;
      Cs[463] = 0.0;
      Cs[464] = 0.0;
      Cs[465] = 0.0;
      Cs[466] = 0.0;
      Cs[467] = 0.0;
      Cs[468] = 0.0;
      Cs[469] = 0.0;
      Cs[470] = 0.0;
      Cs[471] = 0.0;
      Cs[472] = 0.0;
      Cs[473] = 0.0;
      Cs[474] = 0.0;
      Cs[475] = 0.0;
      Cs[476] = 0.0;
      Cs[477] = -6.6666850926047451E-7;
      Cs[478] = 3.3333149073952529E-7;
      Cs[479] = 3.3333149073952529E-7;
      Cs[480] = 1.3333314907395255E-6;
      Cs[481] = 3.3333149073952529E-7;
      Cs[482] = 3.3333149073952529E-7;
      Cs[483] = 0.0;
      Cs[484] = 0.0;
      Cs[485] = 0.0;
      Cs[486] = 0.0;
      Cs[487] = 0.0;
      Cs[488] = 0.0;
      Cs[489] = 0.0;
      Cs[490] = 0.0;
      Cs[491] = 0.0;
      Cs[492] = 0.0;
      Cs[493] = 0.0;
      Cs[494] = 0.0;
      Cs[495] = 0.0;
      Cs[496] = 1.3266675417274417E-9;
      Cs[497] = 0.0;
      Cs[498] = 0.0;
      Cs[499] = 0.0;
      Cs[500] = 0.0;
      Cs[501] = 0.0;
      Cs[502] = 0.0;
      Cs[503] = 0.0;
      Cs[504] = 0.0;
      Cs[505] = 0.0;
      Cs[506] = 0.0;
      Cs[507] = 0.0;
      Cs[508] = 0.0;
      Cs[509] = 0.0;
      Cs[510] = 0.0;
      Cs[511] = 0.0;
      Cs[512] = 0.0;
      Cs[513] = 0.0;
      Cs[514] = 0.0;
      Cs[515] = 0.0;
      Cs[516] = 0.0;
      Cs[517] = 0.0;
      Cs[518] = 0.0;
      Cs[519] = 0.0;
      Cs[520] = 0.0;
      Cs[521] = 0.0;
      Cs[522] = 3.3333149073952524E-7;
      Cs[523] = -6.666685092604744E-7;
      Cs[524] = 3.3333149073952555E-7;
      Cs[525] = 3.3333149073952508E-7;
      Cs[526] = 1.3333314907395253E-6;
      Cs[527] = 3.3333149073952555E-7;
      Cs[528] = 0.0;
      Cs[529] = 0.0;
      Cs[530] = 0.0;
      Cs[531] = 0.0;
      Cs[532] = 0.0;
      Cs[533] = 0.0;
      Cs[534] = 0.0;
      Cs[535] = 0.0;
      Cs[536] = 0.0;
      Cs[537] = 0.0;
      Cs[538] = 0.0;
      Cs[539] = 0.0;
      Cs[540] = 0.0;
      Cs[541] = 1.3266675417274415E-9;
      Cs[542] = 0.0;
      Cs[543] = 0.0;
      Cs[544] = 0.0;
      Cs[545] = 0.0;
      Cs[546] = 0.0;
      Cs[547] = 0.0;
      Cs[548] = 0.0;
      Cs[549] = 0.0;
      Cs[550] = 0.0;
      Cs[551] = 0.0;
      Cs[552] = 0.0;
      Cs[553] = 0.0;
      Cs[554] = 0.0;
      Cs[555] = 0.0;
      Cs[556] = 0.0;
      Cs[557] = 0.0;
      Cs[558] = 0.0;
      Cs[559] = 0.0;
      Cs[560] = 0.0;
      Cs[561] = 0.0;
      Cs[562] = 0.0;
      Cs[563] = 0.0;
      Cs[564] = 0.0;
      Cs[565] = 0.0;
      Cs[566] = 0.0;
      Cs[567] = 3.3333149073952508E-7;
      Cs[568] = 3.3333149073952555E-7;
      Cs[569] = -6.666685092604744E-7;
      Cs[570] = 3.3333149073952497E-7;
      Cs[571] = 3.333314907395255E-7;
      Cs[572] = 1.3333314907395257E-6;
      Cs[573] = 0.0;
      Cs[574] = 0.0;
      Cs[575] = 0.0;
      Cs[576] = 0.0;
      Cs[577] = 0.0;
      Cs[578] = 0.0;
      Cs[579] = 0.0;
      Cs[580] = 0.0;
      Cs[581] = 0.0;
      Cs[582] = 0.0;
      Cs[583] = 0.0;
      Cs[584] = 0.0;
      Cs[585] = 9.9999999998333331E-7;
      Cs[586] = 0.0;
      Cs[587] = 2.3676466022461935E-5;
      Cs[588] = 2.3676466022461924E-5;
      Cs[589] = -1.3331752496216812E-6;
      Cs[590] = 6.66587624810841E-7;
      Cs[591] = 6.6658762481084071E-7;
      Cs[592] = 3.0902962616933355E-19;
      Cs[593] = -1.8450800379097555E-18;
      Cs[594] = -1.5451481308466678E-19;
      Cs[595] = -1.8450800379097539E-18;
      Cs[596] = -7.087915233410791E-34;
      Cs[597] = -1.5451481308466675E-19;
      Cs[598] = -2.8634588722099446E-18;
      Cs[599] = -2.8634588722099458E-18;
      Cs[600] = 5.7269177444198931E-18;
      Cs[601] = -1.2376777139489712E-14;
      Cs[602] = -6.9793364184161707E-31;
      Cs[603] = 1.2376777139489712E-14;
      Cs[604] = -3.9240443053650858E-19;
      Cs[605] = -1.4381920378310571E-34;
      Cs[606] = 3.9240443053650858E-19;
      Cs[607] = -2.40500670782669E-34;
      Cs[608] = -2.1043233889183238E-15;
      Cs[609] = -1.1672617042374267E-31;
      Cs[610] = 2.1043233889183226E-15;
      Cs[611] = 1.128165341042545E-30;
      Cs[612] = 0.0;
      Cs[613] = 0.0;
      Cs[614] = 0.0;
      Cs[615] = 0.0;
      Cs[616] = 0.0;
      Cs[617] = 0.0;
      Cs[618] = -1.3957334077580342E-11;
      Cs[619] = 2.09965651619415E-13;
      Cs[620] = 6.978667038790168E-12;
      Cs[621] = -1.0498282580970747E-13;
      Cs[622] = 6.978667038790168E-12;
      Cs[623] = -1.0498282580970748E-13;
      Cs[624] = 1.2753519916379938E-18;
      Cs[625] = 2.8635085495392202E-18;
      Cs[626] = 1.275351991637994E-18;
      Cs[627] = 2.8635085495392225E-18;
      Cs[628] = -2.5507039832759884E-18;
      Cs[629] = -5.727017099078442E-18;
      Cs[630] = 9.9999999998333331E-7;
      Cs[631] = 0.0;
      Cs[632] = -2.3676466022461935E-5;
      Cs[633] = -2.3676466022461924E-5;
      Cs[634] = 1.3331752496216812E-6;
      Cs[635] = -6.66587624810841E-7;
      Cs[636] = -6.6658762481084071E-7;
      Cs[637] = -3.0902962616933355E-19;
      Cs[638] = 1.8450800379097555E-18;
      Cs[639] = 1.5451481308466678E-19;
      Cs[640] = 1.8450800379097539E-18;
      Cs[641] = 7.087915233410791E-34;
      Cs[642] = 1.5451481308466675E-19;
      Cs[643] = 2.8634588722099446E-18;
      Cs[644] = 2.8634588722099458E-18;
      Cs[645] = -5.7269177444198931E-18;
      Cs[646] = 1.2376777139489712E-14;
      Cs[647] = 6.9793364184161707E-31;
      Cs[648] = -1.2376777139489712E-14;
      Cs[649] = 3.9240443053650858E-19;
      Cs[650] = 1.4381920378310571E-34;
      Cs[651] = -3.9240443053650858E-19;
      Cs[652] = 2.40500670782669E-34;
      Cs[653] = 2.1043233889183238E-15;
      Cs[654] = 1.1672617042374267E-31;
      Cs[655] = -2.1043233889183226E-15;
      Cs[656] = -1.128165341042545E-30;
      Cs[657] = 0.0;
      Cs[658] = 0.0;
      Cs[659] = 0.0;
      Cs[660] = 0.0;
      Cs[661] = 0.0;
      Cs[662] = 0.0;
      Cs[663] = 1.3957334077580342E-11;
      Cs[664] = -2.09965651619415E-13;
      Cs[665] = -6.978667038790168E-12;
      Cs[666] = 1.0498282580970747E-13;
      Cs[667] = -6.978667038790168E-12;
      Cs[668] = 1.0498282580970748E-13;
      Cs[669] = -1.2753519916379938E-18;
      Cs[670] = -2.8635085495392202E-18;
      Cs[671] = -1.275351991637994E-18;
      Cs[672] = -2.8635085495392225E-18;
      Cs[673] = 2.5507039832759884E-18;
      Cs[674] = 5.727017099078442E-18;
      Cs[675] = 9.9999999998333331E-7;
      Cs[676] = 0.0;
      Cs[677] = -2.3676466022461935E-5;
      Cs[678] = -2.6730301104972934E-21;
      Cs[679] = 6.6658762481084071E-7;
      Cs[680] = -1.3331752496216816E-6;
      Cs[681] = 6.6658762481084071E-7;
      Cs[682] = -1.545148130846668E-19;
      Cs[683] = 1.8450800379097555E-18;
      Cs[684] = 3.0902962616933365E-19;
      Cs[685] = -2.9129833672230949E-34;
      Cs[686] = -1.8450800379097555E-18;
      Cs[687] = -1.5451481308466685E-19;
      Cs[688] = 5.7269177444198916E-18;
      Cs[689] = -2.8634588722099481E-18;
      Cs[690] = -2.863458872209947E-18;
      Cs[691] = 1.2376777139489712E-14;
      Cs[692] = -1.2376777139489719E-14;
      Cs[693] = -2.6139936544608445E-30;
      Cs[694] = 3.9240443053650863E-19;
      Cs[695] = -3.9240443053650873E-19;
      Cs[696] = -9.6296497219361785E-36;
      Cs[697] = 1.6358376832359693E-34;
      Cs[698] = 2.1043233889183238E-15;
      Cs[699] = -2.1043233889183242E-15;
      Cs[700] = 3.0067377550530479E-32;
      Cs[701] = -1.1281654707601834E-30;
      Cs[702] = 0.0;
      Cs[703] = 0.0;
      Cs[704] = 0.0;
      Cs[705] = 0.0;
      Cs[706] = 0.0;
      Cs[707] = 0.0;
      Cs[708] = 6.9786670387901712E-12;
      Cs[709] = -1.049828258097075E-13;
      Cs[710] = -1.3957334077580339E-11;
      Cs[711] = 2.0996565161941494E-13;
      Cs[712] = 6.9786670387901712E-12;
      Cs[713] = -1.049828258097075E-13;
      Cs[714] = -2.5507039832759884E-18;
      Cs[715] = -5.7270170990784435E-18;
      Cs[716] = 1.2753519916379952E-18;
      Cs[717] = 2.8635085495392221E-18;
      Cs[718] = 1.2753519916379944E-18;
      Cs[719] = 2.8635085495392218E-18;
      Cs[720] = 9.9999999998333331E-7;
      Cs[721] = 0.0;
      Cs[722] = 2.3676466022461935E-5;
      Cs[723] = 2.6730301104972934E-21;
      Cs[724] = -6.6658762481084071E-7;
      Cs[725] = 1.3331752496216816E-6;
      Cs[726] = -6.6658762481084071E-7;
      Cs[727] = 1.545148130846668E-19;
      Cs[728] = -1.8450800379097555E-18;
      Cs[729] = -3.0902962616933365E-19;
      Cs[730] = 2.9129833672230949E-34;
      Cs[731] = 1.8450800379097555E-18;
      Cs[732] = 1.5451481308466685E-19;
      Cs[733] = -5.7269177444198916E-18;
      Cs[734] = 2.8634588722099481E-18;
      Cs[735] = 2.863458872209947E-18;
      Cs[736] = -1.2376777139489712E-14;
      Cs[737] = 1.2376777139489719E-14;
      Cs[738] = 2.6139936544608445E-30;
      Cs[739] = -3.9240443053650863E-19;
      Cs[740] = 3.9240443053650873E-19;
      Cs[741] = 9.6296497219361785E-36;
      Cs[742] = -1.6358376832359693E-34;
      Cs[743] = -2.1043233889183238E-15;
      Cs[744] = 2.1043233889183242E-15;
      Cs[745] = -3.0067377550530479E-32;
      Cs[746] = 1.1281654707601834E-30;
      Cs[747] = 0.0;
      Cs[748] = 0.0;
      Cs[749] = 0.0;
      Cs[750] = 0.0;
      Cs[751] = 0.0;
      Cs[752] = 0.0;
      Cs[753] = -6.9786670387901712E-12;
      Cs[754] = 1.049828258097075E-13;
      Cs[755] = 1.3957334077580339E-11;
      Cs[756] = -2.0996565161941494E-13;
      Cs[757] = -6.9786670387901712E-12;
      Cs[758] = 1.049828258097075E-13;
      Cs[759] = 2.5507039832759884E-18;
      Cs[760] = 5.7270170990784435E-18;
      Cs[761] = -1.2753519916379952E-18;
      Cs[762] = -2.8635085495392221E-18;
      Cs[763] = -1.2753519916379944E-18;
      Cs[764] = -2.8635085495392218E-18;
      Cs[765] = 9.9999999998333331E-7;
      Cs[766] = 0.0;
      Cs[767] = -2.5243645263569597E-30;
      Cs[768] = -2.3676466022461921E-5;
      Cs[769] = 6.665876248108405E-7;
      Cs[770] = 6.665876248108406E-7;
      Cs[771] = -1.3331752496216814E-6;
      Cs[772] = -1.5451481308466675E-19;
      Cs[773] = 9.9298966930751157E-36;
      Cs[774] = -1.5451481308466687E-19;
      Cs[775] = 1.8450800379097543E-18;
      Cs[776] = 1.8450800379097562E-18;
      Cs[777] = 3.090296261693336E-19;
      Cs[778] = -2.863458872209947E-18;
      Cs[779] = 5.7269177444198939E-18;
      Cs[780] = -2.8634588722099458E-18;
      Cs[781] = 8.9638005363754027E-32;
      Cs[782] = 1.2376777139489719E-14;
      Cs[783] = -1.2376777139489709E-14;
      Cs[784] = -5.099206403174884E-35;
      Cs[785] = 3.9240443053650892E-19;
      Cs[786] = -3.9240443053650858E-19;
      Cs[787] = 7.6916902459072073E-35;
      Cs[788] = 1.523169844767255E-32;
      Cs[789] = 2.1043233889183242E-15;
      Cs[790] = -2.1043233889183226E-15;
      Cs[791] = 1.2971763832316259E-37;
      Cs[792] = 0.0;
      Cs[793] = 0.0;
      Cs[794] = 0.0;
      Cs[795] = 0.0;
      Cs[796] = 0.0;
      Cs[797] = 0.0;
      Cs[798] = 6.9786670387901712E-12;
      Cs[799] = -1.049828258097075E-13;
      Cs[800] = 6.9786670387901712E-12;
      Cs[801] = -1.0498282580970747E-13;
      Cs[802] = -1.3957334077580339E-11;
      Cs[803] = 2.09965651619415E-13;
      Cs[804] = 1.2753519916379946E-18;
      Cs[805] = 2.8635085495392233E-18;
      Cs[806] = -2.5507039832759892E-18;
      Cs[807] = -5.7270170990784443E-18;
      Cs[808] = 1.2753519916379942E-18;
      Cs[809] = 2.8635085495392206E-18;
      Cs[810] = 9.9999999998333331E-7;
      Cs[811] = 0.0;
      Cs[812] = 2.5243645263569597E-30;
      Cs[813] = 2.3676466022461921E-5;
      Cs[814] = -6.665876248108405E-7;
      Cs[815] = -6.665876248108406E-7;
      Cs[816] = 1.3331752496216814E-6;
      Cs[817] = 1.5451481308466675E-19;
      Cs[818] = -9.9298966930751157E-36;
      Cs[819] = 1.5451481308466687E-19;
      Cs[820] = -1.8450800379097543E-18;
      Cs[821] = -1.8450800379097562E-18;
      Cs[822] = -3.090296261693336E-19;
      Cs[823] = 2.863458872209947E-18;
      Cs[824] = -5.7269177444198939E-18;
      Cs[825] = 2.8634588722099458E-18;
      Cs[826] = -8.9638005363754027E-32;
      Cs[827] = -1.2376777139489719E-14;
      Cs[828] = 1.2376777139489709E-14;
      Cs[829] = 5.099206403174884E-35;
      Cs[830] = -3.9240443053650892E-19;
      Cs[831] = 3.9240443053650858E-19;
      Cs[832] = -7.6916902459072073E-35;
      Cs[833] = -1.523169844767255E-32;
      Cs[834] = -2.1043233889183242E-15;
      Cs[835] = 2.1043233889183226E-15;
      Cs[836] = -1.2971763832316259E-37;
      Cs[837] = 0.0;
      Cs[838] = 0.0;
      Cs[839] = 0.0;
      Cs[840] = 0.0;
      Cs[841] = 0.0;
      Cs[842] = 0.0;
      Cs[843] = -6.9786670387901712E-12;
      Cs[844] = 1.049828258097075E-13;
      Cs[845] = -6.9786670387901712E-12;
      Cs[846] = 1.0498282580970747E-13;
      Cs[847] = 1.3957334077580339E-11;
      Cs[848] = -2.09965651619415E-13;
      Cs[849] = -1.2753519916379946E-18;
      Cs[850] = -2.8635085495392233E-18;
      Cs[851] = 2.5507039832759892E-18;
      Cs[852] = 5.7270170990784443E-18;
      Cs[853] = -1.2753519916379942E-18;
      Cs[854] = -2.8635085495392206E-18;
      Cs[855] = 0.0;
      Cs[856] = 2.4432101739489869E-25;
      Cs[857] = 0.0;
      Cs[858] = 0.0;
      Cs[859] = 0.0;
      Cs[860] = 0.0;
      Cs[861] = 0.0;
      Cs[862] = 0.0;
      Cs[863] = 0.0;
      Cs[864] = 0.0;
      Cs[865] = 0.0;
      Cs[866] = 0.0;
      Cs[867] = 0.0;
      Cs[868] = 0.0;
      Cs[869] = 0.0;
      Cs[870] = 0.0;
      Cs[871] = 0.0;
      Cs[872] = 0.0;
      Cs[873] = 0.0;
      Cs[874] = 0.0;
      Cs[875] = 0.0;
      Cs[876] = 0.0;
      Cs[877] = 0.0;
      Cs[878] = 0.0;
      Cs[879] = 0.0;
      Cs[880] = 0.0;
      Cs[881] = 0.0;
      Cs[882] = 9.9999999999999953E-7;
      Cs[883] = 2.19958495703998E-22;
      Cs[884] = -9.9999999999999953E-7;
      Cs[885] = -1.0000000000000004E-6;
      Cs[886] = 2.201753361384951E-22;
      Cs[887] = 1.0000000000000004E-6;
      Cs[888] = 0.0;
      Cs[889] = 0.0;
      Cs[890] = 0.0;
      Cs[891] = 0.0;
      Cs[892] = 0.0;
      Cs[893] = 0.0;
      Cs[894] = 0.0;
      Cs[895] = 0.0;
      Cs[896] = 0.0;
      Cs[897] = 0.0;
      Cs[898] = 0.0;
      Cs[899] = 0.0;
      Cs[900] = 0.0;
      Cs[901] = -2.4286270368118714E-25;
      Cs[902] = 0.0;
      Cs[903] = 0.0;
      Cs[904] = 0.0;
      Cs[905] = 0.0;
      Cs[906] = 0.0;
      Cs[907] = 0.0;
      Cs[908] = 0.0;
      Cs[909] = 0.0;
      Cs[910] = 0.0;
      Cs[911] = 0.0;
      Cs[912] = 0.0;
      Cs[913] = 0.0;
      Cs[914] = 0.0;
      Cs[915] = 0.0;
      Cs[916] = 0.0;
      Cs[917] = 0.0;
      Cs[918] = 0.0;
      Cs[919] = 0.0;
      Cs[920] = 0.0;
      Cs[921] = 0.0;
      Cs[922] = 0.0;
      Cs[923] = 0.0;
      Cs[924] = 0.0;
      Cs[925] = 0.0;
      Cs[926] = 0.0;
      Cs[927] = 3.3053209482886144E-22;
      Cs[928] = 1.0E-6;
      Cs[929] = 9.9999999999999953E-7;
      Cs[930] = 3.3031525439436431E-22;
      Cs[931] = -9.9999999999999974E-7;
      Cs[932] = -1.0000000000000002E-6;
      Cs[933] = 0.0;
      Cs[934] = 0.0;
      Cs[935] = 0.0;
      Cs[936] = 0.0;
      Cs[937] = 0.0;
      Cs[938] = 0.0;
      Cs[939] = 0.0;
      Cs[940] = 0.0;
      Cs[941] = 0.0;
      Cs[942] = 0.0;
      Cs[943] = 0.0;
      Cs[944] = 0.0;
      Cs[945] = 0.0;
      Cs[946] = 0.0;
      Cs[947] = 7.0263438718684876E-17;
      Cs[948] = 7.0263438718060567E-17;
      Cs[949] = 6.60838689858856E-11;
      Cs[950] = -1.6164350224291195E-14;
      Cs[951] = -1.6164350224275131E-14;
      Cs[952] = -3.4419995693636558E-12;
      Cs[953] = 2.0550666208446102E-11;
      Cs[954] = 1.7209997846818275E-12;
      Cs[955] = 2.0550666208446095E-11;
      Cs[956] = 6.0218093927089814E-27;
      Cs[957] = 1.7209997846818281E-12;
      Cs[958] = 3.1893460595381663E-11;
      Cs[959] = 3.1893460595381689E-11;
      Cs[960] = -6.3786921190763352E-11;
      Cs[961] = 1.3785364889543075E-7;
      Cs[962] = -1.0618922799415422E-23;
      Cs[963] = -1.3785364889543077E-7;
      Cs[964] = 4.3706355848968286E-12;
      Cs[965] = 1.0097420490533719E-27;
      Cs[966] = -4.370635584896827E-12;
      Cs[967] = 1.4901814865895853E-28;
      Cs[968] = 2.3438141799679339E-8;
      Cs[969] = -3.1057704510113295E-24;
      Cs[970] = -2.3438141799679345E-8;
      Cs[971] = 3.0969454693915771E-24;
      Cs[972] = 0.0;
      Cs[973] = 0.0;
      Cs[974] = 0.0;
      Cs[975] = 0.0;
      Cs[976] = 0.0;
      Cs[977] = 0.0;
      Cs[978] = 0.00023320319889241412;
      Cs[979] = -3.5079242915423036E-6;
      Cs[980] = 1.6149093965666937E-8;
      Cs[981] = -2.8356065028487405E-14;
      Cs[982] = 1.6149093965666934E-8;
      Cs[983] = -2.8356065028459293E-14;
      Cs[984] = -1.4204984358359455E-11;
      Cs[985] = -3.18940139059108E-11;
      Cs[986] = -1.4204984358359465E-11;
      Cs[987] = -3.1894013905910812E-11;
      Cs[988] = 2.8409968716718923E-11;
      Cs[989] = 6.3788027811821586E-11;
      Cs[990] = 0.0;
      Cs[991] = 0.0;
      Cs[992] = -7.0263438717715526E-17;
      Cs[993] = -2.636569695429025E-32;
      Cs[994] = -1.616435022429252E-14;
      Cs[995] = 6.6083868985207966E-11;
      Cs[996] = -1.6164350224275175E-14;
      Cs[997] = 1.7209997846818279E-12;
      Cs[998] = -2.0550666208446095E-11;
      Cs[999] = -3.4419995693636562E-12;
      Cs[1000] = -4.1715359818585954E-27;
      Cs[1001] = 2.0550666208446108E-11;
      Cs[1002] = 1.7209997846818281E-12;
      Cs[1003] = -6.3786921190763352E-11;
      Cs[1004] = 3.1893460595381689E-11;
      Cs[1005] = 3.1893460595381689E-11;
      Cs[1006] = -1.3785364889543069E-7;
      Cs[1007] = 1.3785364889543085E-7;
      Cs[1008] = 7.418750941791597E-23;
      Cs[1009] = -4.3706355848968278E-12;
      Cs[1010] = 4.3706355848968286E-12;
      Cs[1011] = 2.8272777413726698E-28;
      Cs[1012] = 9.592141190048298E-28;
      Cs[1013] = -2.3438141799679342E-8;
      Cs[1014] = 2.3438141799679342E-8;
      Cs[1015] = -1.5319402385484054E-24;
      Cs[1016] = 3.3638925418486292E-23;
      Cs[1017] = 0.0;
      Cs[1018] = 0.0;
      Cs[1019] = 0.0;
      Cs[1020] = 0.0;
      Cs[1021] = 0.0;
      Cs[1022] = 0.0;
      Cs[1023] = 1.6149093965666944E-8;
      Cs[1024] = -2.8356065028489709E-14;
      Cs[1025] = 0.00023320319889241414;
      Cs[1026] = -3.507924291542304E-6;
      Cs[1027] = 1.6149093965666934E-8;
      Cs[1028] = -2.8356065028459293E-14;
      Cs[1029] = 2.8409968716718923E-11;
      Cs[1030] = 6.3788027811821638E-11;
      Cs[1031] = -1.4204984358359465E-11;
      Cs[1032] = -3.1894013905910806E-11;
      Cs[1033] = -1.4204984358359465E-11;
      Cs[1034] = -3.1894013905910806E-11;
      Cs[1035] = 0.0;
      Cs[1036] = 0.0;
      Cs[1037] = 3.2581700205301994E-32;
      Cs[1038] = -7.0263438717960159E-17;
      Cs[1039] = -1.6164350224305833E-14;
      Cs[1040] = -1.6164350224336488E-14;
      Cs[1041] = 6.6083868985546774E-11;
      Cs[1042] = 1.7209997846818269E-12;
      Cs[1043] = -6.1411451121696791E-27;
      Cs[1044] = 1.7209997846818279E-12;
      Cs[1045] = -2.0550666208446079E-11;
      Cs[1046] = -2.0550666208446111E-11;
      Cs[1047] = -3.441999569363655E-12;
      Cs[1048] = 3.1893460595381676E-11;
      Cs[1049] = -6.3786921190763378E-11;
      Cs[1050] = 3.1893460595381663E-11;
      Cs[1051] = -4.6339949551666365E-23;
      Cs[1052] = -1.3785364889543083E-7;
      Cs[1053] = 1.3785364889543064E-7;
      Cs[1054] = 2.3505766401477092E-28;
      Cs[1055] = -4.3706355848968278E-12;
      Cs[1056] = 4.3706355848968254E-12;
      Cs[1057] = -2.5837179137465508E-27;
      Cs[1058] = -9.4915752614236936E-28;
      Cs[1059] = -2.3438141799679345E-8;
      Cs[1060] = 2.3438141799679332E-8;
      Cs[1061] = -1.8735643676428758E-23;
      Cs[1062] = 0.0;
      Cs[1063] = 0.0;
      Cs[1064] = 0.0;
      Cs[1065] = 0.0;
      Cs[1066] = 0.0;
      Cs[1067] = 0.0;
      Cs[1068] = 1.6149093965666937E-8;
      Cs[1069] = -2.8356065028513094E-14;
      Cs[1070] = 1.6149093965666944E-8;
      Cs[1071] = -2.8356065028566856E-14;
      Cs[1072] = 0.00023320319889241412;
      Cs[1073] = -3.5079242915423036E-6;
      Cs[1074] = -1.4204984358359461E-11;
      Cs[1075] = -3.1894013905910812E-11;
      Cs[1076] = 2.8409968716718929E-11;
      Cs[1077] = 6.3788027811821612E-11;
      Cs[1078] = -1.4204984358359457E-11;
      Cs[1079] = -3.1894013905910787E-11;
      Cs[1080] = 0.0;
      Cs[1081] = 0.0;
      Cs[1082] = 5.4667094528831493E-22;
      Cs[1083] = 2.7333547264403286E-22;
      Cs[1084] = 2.5713919192226255E-16;
      Cs[1085] = -2.57139191922979E-16;
      Cs[1086] = -2.4113355198243367E-28;
      Cs[1087] = -7.8617878122916136E-7;
      Cs[1088] = 1.5110330390498779E-5;
      Cs[1089] = 7.8617878122916136E-7;
      Cs[1090] = -5.7224746469895382E-6;
      Cs[1091] = 5.7224746469895382E-6;
      Cs[1092] = 0.0;
      Cs[1093] = 1.4569416105239678E-5;
      Cs[1094] = 1.6018554777131082E-21;
      Cs[1095] = -1.456941610523968E-5;
      Cs[1096] = -1.6685319270023485E-5;
      Cs[1097] = -4.9422916823559218E-6;
      Cs[1098] = -4.9422916823559243E-6;
      Cs[1099] = 1.9969141935374475E-6;
      Cs[1100] = 3.4162687309727265E-10;
      Cs[1101] = 3.4162687309728656E-10;
      Cs[1102] = 1.9975974472836407E-6;
      Cs[1103] = 3.3748586594716696E-10;
      Cs[1104] = -3.5317730899323527E-11;
      Cs[1105] = -3.5317730899753719E-11;
      Cs[1106] = 2.668504041470552E-10;
      Cs[1107] = 0.0;
      Cs[1108] = 0.0;
      Cs[1109] = 0.0;
      Cs[1110] = 0.0;
      Cs[1111] = 0.0;
      Cs[1112] = 0.0;
      Cs[1113] = -2.5686508630556174E-10;
      Cs[1114] = 4.5102735213575377E-16;
      Cs[1115] = 2.5686508630546155E-10;
      Cs[1116] = -4.5102735213701042E-16;
      Cs[1117] = 2.4097081809160765E-22;
      Cs[1118] = -4.2295308514019789E-28;
      Cs[1119] = -6.4890521135648919E-6;
      Cs[1120] = -1.4569668865874114E-5;
      Cs[1121] = -7.2288125548153679E-22;
      Cs[1122] = -1.585809541092895E-21;
      Cs[1123] = 6.4890521135648927E-6;
      Cs[1124] = 1.4569668865874105E-5;
      Cs[1125] = 0.0;
      Cs[1126] = 0.0;
      Cs[1127] = -2.7333547264458911E-22;
      Cs[1128] = 2.7333547264389072E-22;
      Cs[1129] = 1.7126124703152139E-28;
      Cs[1130] = 2.5713919192319814E-16;
      Cs[1131] = -2.571391919221988E-16;
      Cs[1132] = 0.0;
      Cs[1133] = 5.7224746469895382E-6;
      Cs[1134] = -7.8617878122916158E-7;
      Cs[1135] = -5.7224746469895382E-6;
      Cs[1136] = 1.5110330390498781E-5;
      Cs[1137] = 7.8617878122916158E-7;
      Cs[1138] = -1.4569416105239687E-5;
      Cs[1139] = 1.4569416105239693E-5;
      Cs[1140] = 1.6145512835029284E-21;
      Cs[1141] = -4.9422916823559268E-6;
      Cs[1142] = -1.6685319270023495E-5;
      Cs[1143] = -4.9422916823559277E-6;
      Cs[1144] = 3.4162687309693961E-10;
      Cs[1145] = 1.9969141935374471E-6;
      Cs[1146] = 3.4162687309699684E-10;
      Cs[1147] = 1.9975974472836412E-6;
      Cs[1148] = -3.5317730899391479E-11;
      Cs[1149] = 3.3748586594628534E-10;
      Cs[1150] = -3.531773089961284E-11;
      Cs[1151] = 2.6685040414503311E-10;
      Cs[1152] = 0.0;
      Cs[1153] = 0.0;
      Cs[1154] = 0.0;
      Cs[1155] = 0.0;
      Cs[1156] = 0.0;
      Cs[1157] = 0.0;
      Cs[1158] = -1.7108753961908385E-22;
      Cs[1159] = 3.0039563550052118E-28;
      Cs[1160] = -2.5686508630535329E-10;
      Cs[1161] = 4.510273521373948E-16;
      Cs[1162] = 2.5686508630513636E-10;
      Cs[1163] = -4.51027352135642E-16;
      Cs[1164] = 6.4890521135648961E-6;
      Cs[1165] = 1.4569668865874122E-5;
      Cs[1166] = -6.4890521135648986E-6;
      Cs[1167] = -1.4569668865874121E-5;
      Cs[1168] = -4.3752191441774572E-22;
      Cs[1169] = -1.6145969095671509E-21;
      Cs[1170] = 0.0;
      Cs[1171] = 0.0;
      Cs[1172] = -2.7333547264326145E-22;
      Cs[1173] = -5.4667094528703684E-22;
      Cs[1174] = -2.5713919192181823E-16;
      Cs[1175] = -4.0300167501873974E-29;
      Cs[1176] = 2.5713919192222138E-16;
      Cs[1177] = 7.8617878122916158E-7;
      Cs[1178] = 5.7224746469895382E-6;
      Cs[1179] = 0.0;
      Cs[1180] = -1.5110330390498781E-5;
      Cs[1181] = 5.7224746469895382E-6;
      Cs[1182] = -7.8617878122916158E-7;
      Cs[1183] = -2.0028786700137359E-21;
      Cs[1184] = -1.4569416105239692E-5;
      Cs[1185] = 1.456941610523969E-5;
      Cs[1186] = -4.9422916823559243E-6;
      Cs[1187] = -4.9422916823559251E-6;
      Cs[1188] = -1.6685319270023495E-5;
      Cs[1189] = 3.4162687309693961E-10;
      Cs[1190] = 3.41626873096983E-10;
      Cs[1191] = 1.9969141935374471E-6;
      Cs[1192] = 1.9975974472836412E-6;
      Cs[1193] = -3.5317730899375E-11;
      Cs[1194] = -3.5317730899626178E-11;
      Cs[1195] = 3.3748586594639308E-10;
      Cs[1196] = 2.668504041463005E-10;
      Cs[1197] = 0.0;
      Cs[1198] = 0.0;
      Cs[1199] = 0.0;
      Cs[1200] = 0.0;
      Cs[1201] = 0.0;
      Cs[1202] = 0.0;
      Cs[1203] = 2.5686508630516009E-10;
      Cs[1204] = -4.5102735213497428E-16;
      Cs[1205] = 4.0222024706156658E-23;
      Cs[1206] = -7.0687329298380579E-29;
      Cs[1207] = -2.5686508630529673E-10;
      Cs[1208] = 4.5102735213568149E-16;
      Cs[1209] = 7.5865590252763527E-22;
      Cs[1210] = 2.0779783368662731E-21;
      Cs[1211] = 6.4890521135648978E-6;
      Cs[1212] = 1.4569668865874121E-5;
      Cs[1213] = -6.4890521135648969E-6;
      Cs[1214] = -1.4569668865874117E-5;
      Cs[1215] = 0.0;
      Cs[1216] = 2.4432101739489869E-25;
      Cs[1217] = 0.0;
      Cs[1218] = 0.0;
      Cs[1219] = 0.0;
      Cs[1220] = 0.0;
      Cs[1221] = 0.0;
      Cs[1222] = 0.0;
      Cs[1223] = 0.0;
      Cs[1224] = 0.0;
      Cs[1225] = 0.0;
      Cs[1226] = 0.0;
      Cs[1227] = 0.0;
      Cs[1228] = 0.0;
      Cs[1229] = 0.0;
      Cs[1230] = 0.0;
      Cs[1231] = 0.0;
      Cs[1232] = 0.0;
      Cs[1233] = 0.0;
      Cs[1234] = 0.0;
      Cs[1235] = 0.0;
      Cs[1236] = 0.0;
      Cs[1237] = 0.0;
      Cs[1238] = 0.0;
      Cs[1239] = 0.0;
      Cs[1240] = 0.0;
      Cs[1241] = 0.0;
      Cs[1242] = 9.9999999999999953E-7;
      Cs[1243] = 2.19958495703998E-22;
      Cs[1244] = -9.9999999999999953E-7;
      Cs[1245] = -1.0000000000000004E-6;
      Cs[1246] = 2.201753361384951E-22;
      Cs[1247] = 1.0000000000000004E-6;
      Cs[1248] = 0.0;
      Cs[1249] = 0.0;
      Cs[1250] = 0.0;
      Cs[1251] = 0.0;
      Cs[1252] = 0.0;
      Cs[1253] = 0.0;
      Cs[1254] = 0.0;
      Cs[1255] = 0.0;
      Cs[1256] = 0.0;
      Cs[1257] = 0.0;
      Cs[1258] = 0.0;
      Cs[1259] = 0.0;
      Cs[1260] = 0.0;
      Cs[1261] = -2.4286270368118714E-25;
      Cs[1262] = 0.0;
      Cs[1263] = 0.0;
      Cs[1264] = 0.0;
      Cs[1265] = 0.0;
      Cs[1266] = 0.0;
      Cs[1267] = 0.0;
      Cs[1268] = 0.0;
      Cs[1269] = 0.0;
      Cs[1270] = 0.0;
      Cs[1271] = 0.0;
      Cs[1272] = 0.0;
      Cs[1273] = 0.0;
      Cs[1274] = 0.0;
      Cs[1275] = 0.0;
      Cs[1276] = 0.0;
      Cs[1277] = 0.0;
      Cs[1278] = 0.0;
      Cs[1279] = 0.0;
      Cs[1280] = 0.0;
      Cs[1281] = 0.0;
      Cs[1282] = 0.0;
      Cs[1283] = 0.0;
      Cs[1284] = 0.0;
      Cs[1285] = 0.0;
      Cs[1286] = 0.0;
      Cs[1287] = 3.3053209482886144E-22;
      Cs[1288] = 1.0E-6;
      Cs[1289] = 9.9999999999999953E-7;
      Cs[1290] = 3.3031525439436431E-22;
      Cs[1291] = -9.9999999999999974E-7;
      Cs[1292] = -1.0000000000000002E-6;
      Cs[1293] = 0.0;
      Cs[1294] = 0.0;
      Cs[1295] = 0.0;
      Cs[1296] = 0.0;
      Cs[1297] = 0.0;
      Cs[1298] = 0.0;
      Cs[1299] = 0.0;
      Cs[1300] = 0.0;
      Cs[1301] = 0.0;
      Cs[1302] = 0.0;
      Cs[1303] = 0.0;
      Cs[1304] = 0.0;
      Cs[1305] = 0.0;
      Cs[1306] = 1.4473635564792414E-27;
      Cs[1307] = 0.0;
      Cs[1308] = 0.0;
      Cs[1309] = 0.0;
      Cs[1310] = 0.0;
      Cs[1311] = 0.0;
      Cs[1312] = 0.0;
      Cs[1313] = 0.0;
      Cs[1314] = 0.0;
      Cs[1315] = 0.0;
      Cs[1316] = 0.0;
      Cs[1317] = 0.0;
      Cs[1318] = 0.0;
      Cs[1319] = 0.0;
      Cs[1320] = 0.0;
      Cs[1321] = 0.0;
      Cs[1322] = 0.0;
      Cs[1323] = 0.0;
      Cs[1324] = 0.0;
      Cs[1325] = 0.0;
      Cs[1326] = 0.0;
      Cs[1327] = 0.0;
      Cs[1328] = 0.0;
      Cs[1329] = 0.0;
      Cs[1330] = 0.0;
      Cs[1331] = 0.0;
      Cs[1332] = -1.0E-6;
      Cs[1333] = -1.0000000000000002E-6;
      Cs[1334] = -5.5294309277683945E-23;
      Cs[1335] = 1.0000000000000002E-6;
      Cs[1336] = 9.9999999999999974E-7;
      Cs[1337] = -5.4860628408689733E-23;
      Cs[1338] = 0.0;
      Cs[1339] = 0.0;
      Cs[1340] = 0.0;
      Cs[1341] = 0.0;
      Cs[1342] = 0.0;
      Cs[1343] = 0.0;
      Cs[1344] = 0.0;
      Cs[1345] = 0.0;
      Cs[1346] = 0.0;
      Cs[1347] = 0.0;
      Cs[1348] = 0.0;
      Cs[1349] = 0.0;
      Cs[1350] = 0.0;
      Cs[1351] = 0.0;
      Cs[1352] = -2.3676466022461935E-5;
      Cs[1353] = -2.3676466022461931E-5;
      Cs[1354] = 1.9998418564285624E-6;
      Cs[1355] = 7.8981996039653213E-11;
      Cs[1356] = 7.898199604019428E-11;
      Cs[1357] = -3.090296261693335E-19;
      Cs[1358] = 1.8450800379097555E-18;
      Cs[1359] = 1.5451481308466675E-19;
      Cs[1360] = 1.8450800379097543E-18;
      Cs[1361] = 6.5258724807712384E-34;
      Cs[1362] = 1.5451481308466675E-19;
      Cs[1363] = 2.863458872209945E-18;
      Cs[1364] = 2.8634588722099466E-18;
      Cs[1365] = -5.7269177444198939E-18;
      Cs[1366] = 1.2376777139489712E-14;
      Cs[1367] = 1.2232164498678737E-31;
      Cs[1368] = -1.2376777139489711E-14;
      Cs[1369] = 3.9240443053650858E-19;
      Cs[1370] = 1.6046837578707044E-34;
      Cs[1371] = -3.9240443053650858E-19;
      Cs[1372] = 1.9652727594186308E-34;
      Cs[1373] = 2.1043233889183238E-15;
      Cs[1374] = -3.0443330096005687E-31;
      Cs[1375] = -2.104323388918323E-15;
      Cs[1376] = -5.9470069772072334E-31;
      Cs[1377] = 0.0;
      Cs[1378] = 0.0;
      Cs[1379] = 0.0;
      Cs[1380] = 0.0;
      Cs[1381] = 0.0;
      Cs[1382] = 0.0;
      Cs[1383] = 2.0938278067832848E-11;
      Cs[1384] = -3.1496091907980215E-13;
      Cs[1385] = 2.2769514623383949E-15;
      Cs[1386] = -1.2441650679759871E-17;
      Cs[1387] = 2.2769514623392769E-15;
      Cs[1388] = -1.2441650679752264E-17;
      Cs[1389] = -1.2753519916379932E-18;
      Cs[1390] = -2.8635085495392202E-18;
      Cs[1391] = -1.2753519916379938E-18;
      Cs[1392] = -2.8635085495392229E-18;
      Cs[1393] = 2.550703983275988E-18;
      Cs[1394] = 5.7270170990784428E-18;
      Cs[1395] = 0.0;
      Cs[1396] = 0.0;
      Cs[1397] = 2.3676466022461935E-5;
      Cs[1398] = -7.7234916200485644E-22;
      Cs[1399] = 7.8981996040111342E-11;
      Cs[1400] = 1.9998418564285628E-6;
      Cs[1401] = 7.898199604019428E-11;
      Cs[1402] = 1.5451481308466675E-19;
      Cs[1403] = -1.8450800379097547E-18;
      Cs[1404] = -3.090296261693336E-19;
      Cs[1405] = 5.791388270476916E-34;
      Cs[1406] = 1.8450800379097547E-18;
      Cs[1407] = 1.5451481308466685E-19;
      Cs[1408] = -5.7269177444198924E-18;
      Cs[1409] = 2.8634588722099481E-18;
      Cs[1410] = 2.863458872209947E-18;
      Cs[1411] = -1.2376777139489711E-14;
      Cs[1412] = 1.2376777139489719E-14;
      Cs[1413] = 1.7235380106143467E-30;
      Cs[1414] = -3.9240443053650858E-19;
      Cs[1415] = 3.9240443053650878E-19;
      Cs[1416] = -9.2145570638893566E-36;
      Cs[1417] = -2.0755716316440288E-34;
      Cs[1418] = -2.1043233889183238E-15;
      Cs[1419] = 2.1043233889183238E-15;
      Cs[1420] = -3.6892567362620958E-31;
      Cs[1421] = 1.661630114082005E-30;
      Cs[1422] = 0.0;
      Cs[1423] = 0.0;
      Cs[1424] = 0.0;
      Cs[1425] = 0.0;
      Cs[1426] = 0.0;
      Cs[1427] = 0.0;
      Cs[1428] = 2.2769514623351034E-15;
      Cs[1429] = -1.2441650679726298E-17;
      Cs[1430] = 2.0938278067832845E-11;
      Cs[1431] = -3.1496091907980215E-13;
      Cs[1432] = 2.2769514623358887E-15;
      Cs[1433] = -1.2441650679752264E-17;
      Cs[1434] = 2.550703983275988E-18;
      Cs[1435] = 5.7270170990784443E-18;
      Cs[1436] = -1.275351991637995E-18;
      Cs[1437] = -2.8635085495392229E-18;
      Cs[1438] = -1.2753519916379942E-18;
      Cs[1439] = -2.8635085495392221E-18;
      Cs[1440] = 0.0;
      Cs[1441] = 0.0;
      Cs[1442] = 5.0945674233941068E-30;
      Cs[1443] = 2.3676466022461918E-5;
      Cs[1444] = 7.8981996040222365E-11;
      Cs[1445] = 7.89819960400973E-11;
      Cs[1446] = 1.9998418564285628E-6;
      Cs[1447] = 1.545148130846667E-19;
      Cs[1448] = 8.5799755144913863E-35;
      Cs[1449] = 1.5451481308466685E-19;
      Cs[1450] = -1.8450800379097539E-18;
      Cs[1451] = -1.8450800379097562E-18;
      Cs[1452] = -3.090296261693335E-19;
      Cs[1453] = 2.8634588722099473E-18;
      Cs[1454] = -5.7269177444198947E-18;
      Cs[1455] = 2.8634588722099458E-18;
      Cs[1456] = 1.1417229085798943E-30;
      Cs[1457] = -1.237677713948972E-14;
      Cs[1458] = 1.2376777139489706E-14;
      Cs[1459] = 1.0442642775757775E-34;
      Cs[1460] = -3.9240443053650892E-19;
      Cs[1461] = 3.9240443053650853E-19;
      Cs[1462] = -1.20890297299878E-34;
      Cs[1463] = -2.8659637759534537E-32;
      Cs[1464] = -2.1043233889183254E-15;
      Cs[1465] = 2.1043233889183226E-15;
      Cs[1466] = 5.3346451360418329E-31;
      Cs[1467] = 0.0;
      Cs[1468] = 0.0;
      Cs[1469] = 0.0;
      Cs[1470] = 0.0;
      Cs[1471] = 0.0;
      Cs[1472] = 0.0;
      Cs[1473] = 2.2769514623359506E-15;
      Cs[1474] = -1.2441650679739533E-17;
      Cs[1475] = 2.2769514623358536E-15;
      Cs[1476] = -1.2441650679759871E-17;
      Cs[1477] = 2.0938278067832845E-11;
      Cs[1478] = -3.1496091907980215E-13;
      Cs[1479] = -1.2753519916379942E-18;
      Cs[1480] = -2.8635085495392233E-18;
      Cs[1481] = 2.5507039832759888E-18;
      Cs[1482] = 5.7270170990784451E-18;
      Cs[1483] = -1.275351991637994E-18;
      Cs[1484] = -2.863508549539221E-18;
      Cs[1485] = 0.0;
      Cs[1486] = 0.0;
      Cs[1487] = 2.1259724860852949E-12;
      Cs[1488] = 2.1259724860852937E-12;
      Cs[1489] = 1.9999997007113151E-6;
      Cs[1490] = 5.9854663875899985E-14;
      Cs[1491] = 5.9854663875899973E-14;
      Cs[1492] = -3.0906621442710843E-19;
      Cs[1493] = 1.8452984903114843E-18;
      Cs[1494] = 1.5453310721355419E-19;
      Cs[1495] = 1.8452984903114839E-18;
      Cs[1496] = 5.3780791803324339E-34;
      Cs[1497] = 1.5453310721355424E-19;
      Cs[1498] = 2.8637978978646771E-18;
      Cs[1499] = 2.8637978978646794E-18;
      Cs[1500] = -5.7275957957293565E-18;
      Cs[1501] = 1.2378242515861597E-14;
      Cs[1502] = -9.5349642629419209E-31;
      Cs[1503] = -1.23782425158616E-14;
      Cs[1504] = 3.9245089014180426E-19;
      Cs[1505] = 9.0370476880492864E-35;
      Cs[1506] = -3.9245089014180411E-19;
      Cs[1507] = 1.2546776145428408E-35;
      Cs[1508] = 2.1045725350197816E-15;
      Cs[1509] = -2.8016163367421805E-31;
      Cs[1510] = -2.1045725350197824E-15;
      Cs[1511] = 2.8114254507597906E-31;
      Cs[1512] = 0.0;
      Cs[1513] = 0.0;
      Cs[1514] = 0.0;
      Cs[1515] = 0.0;
      Cs[1516] = 0.0;
      Cs[1517] = 0.0;
      Cs[1518] = 2.0939930577798992E-11;
      Cs[1519] = -3.1498577843545848E-13;
      Cs[1520] = 1.4506964792651322E-15;
      Cs[1521] = -1.1972851563968945E-20;
      Cs[1522] = 1.450696479265132E-15;
      Cs[1523] = -1.1972851563966325E-20;
      Cs[1524] = -1.2755029898060397E-18;
      Cs[1525] = -2.8638475810756121E-18;
      Cs[1526] = -1.2755029898060404E-18;
      Cs[1527] = -2.8638475810756137E-18;
      Cs[1528] = 2.5510059796120805E-18;
      Cs[1529] = 5.7276951621512227E-18;
      Cs[1530] = 0.0;
      Cs[1531] = 0.0;
      Cs[1532] = -2.1259724860852949E-12;
      Cs[1533] = -2.4001844125629154E-28;
      Cs[1534] = 5.9854663875899973E-14;
      Cs[1535] = 1.9999997007113151E-6;
      Cs[1536] = 5.9854663875899973E-14;
      Cs[1537] = 1.5453310721355421E-19;
      Cs[1538] = -1.8452984903114839E-18;
      Cs[1539] = -3.0906621442710847E-19;
      Cs[1540] = -3.8015341498097066E-34;
      Cs[1541] = 1.8452984903114847E-18;
      Cs[1542] = 1.5453310721355424E-19;
      Cs[1543] = -5.7275957957293565E-18;
      Cs[1544] = 2.8637978978646786E-18;
      Cs[1545] = 2.8637978978646786E-18;
      Cs[1546] = -1.2378242515861591E-14;
      Cs[1547] = 1.2378242515861607E-14;
      Cs[1548] = 6.6793500054945339E-30;
      Cs[1549] = -3.9245089014180416E-19;
      Cs[1550] = 3.9245089014180426E-19;
      Cs[1551] = 2.570605633138817E-35;
      Cs[1552] = 8.7970900435169852E-35;
      Cs[1553] = -2.1045725350197824E-15;
      Cs[1554] = 2.104572535019782E-15;
      Cs[1555] = -1.2851306430175919E-31;
      Cs[1556] = 3.0062406287685356E-30;
      Cs[1557] = 0.0;
      Cs[1558] = 0.0;
      Cs[1559] = 0.0;
      Cs[1560] = 0.0;
      Cs[1561] = 0.0;
      Cs[1562] = 0.0;
      Cs[1563] = 1.4506964792651328E-15;
      Cs[1564] = -1.1972851563968975E-20;
      Cs[1565] = 2.0939930577798995E-11;
      Cs[1566] = -3.1498577843545858E-13;
      Cs[1567] = 1.450696479265132E-15;
      Cs[1568] = -1.1972851563966324E-20;
      Cs[1569] = 2.5510059796120805E-18;
      Cs[1570] = 5.7276951621512281E-18;
      Cs[1571] = -1.2755029898060404E-18;
      Cs[1572] = -2.8638475810756133E-18;
      Cs[1573] = -1.2755029898060404E-18;
      Cs[1574] = -2.8638475810756125E-18;
      Cs[1575] = 0.0;
      Cs[1576] = 0.0;
      Cs[1577] = -2.1962886906723967E-37;
      Cs[1578] = -2.1259724860852933E-12;
      Cs[1579] = 5.985466387589996E-14;
      Cs[1580] = 5.985466387589996E-14;
      Cs[1581] = 1.9999997007113151E-6;
      Cs[1582] = 1.5453310721355412E-19;
      Cs[1583] = -5.53459329093419E-34;
      Cs[1584] = 1.5453310721355421E-19;
      Cs[1585] = -1.8452984903114824E-18;
      Cs[1586] = -1.8452984903114851E-18;
      Cs[1587] = -3.0906621442710838E-19;
      Cs[1588] = 2.8637978978646782E-18;
      Cs[1589] = -5.7275957957293573E-18;
      Cs[1590] = 2.8637978978646775E-18;
      Cs[1591] = -4.17252330920538E-30;
      Cs[1592] = -1.2378242515861603E-14;
      Cs[1593] = 1.2378242515861588E-14;
      Cs[1594] = 2.0851880588576509E-35;
      Cs[1595] = -3.9245089014180421E-19;
      Cs[1596] = 3.9245089014180392E-19;
      Cs[1597] = -2.3243786110301596E-34;
      Cs[1598] = -8.49804478805489E-35;
      Cs[1599] = -2.1045725350197824E-15;
      Cs[1600] = 2.1045725350197812E-15;
      Cs[1601] = -1.6869892438790483E-30;
      Cs[1602] = 0.0;
      Cs[1603] = 0.0;
      Cs[1604] = 0.0;
      Cs[1605] = 0.0;
      Cs[1606] = 0.0;
      Cs[1607] = 0.0;
      Cs[1608] = 1.4506964792651322E-15;
      Cs[1609] = -1.1972851563971226E-20;
      Cs[1610] = 1.4506964792651328E-15;
      Cs[1611] = -1.1972851563975919E-20;
      Cs[1612] = 2.0939930577798992E-11;
      Cs[1613] = -3.1498577843545848E-13;
      Cs[1614] = -1.2755029898060402E-18;
      Cs[1615] = -2.8638475810756137E-18;
      Cs[1616] = 2.5510059796120809E-18;
      Cs[1617] = 5.7276951621512266E-18;
      Cs[1618] = -1.2755029898060399E-18;
      Cs[1619] = -2.8638475810756106E-18;
      Cs[1620] = 0.0;
      Cs[1621] = 0.0;
      Cs[1622] = -2.1259022226465763E-17;
      Cs[1623] = -2.1259022226465756E-17;
      Cs[1624] = -1.9999336168423291E-11;
      Cs[1625] = -7.6019014100191184E-19;
      Cs[1626] = -7.60190141001751E-19;
      Cs[1627] = -3.4419992602974416E-17;
      Cs[1628] = 2.0550664363147614E-16;
      Cs[1629] = 1.7209996301487205E-17;
      Cs[1630] = 2.0550664363147606E-16;
      Cs[1631] = 6.0218088549010633E-32;
      Cs[1632] = 1.7209996301487211E-17;
      Cs[1633] = 3.1893457731583766E-16;
      Cs[1634] = 3.1893457731583791E-16;
      Cs[1635] = -6.3786915463167562E-16;
      Cs[1636] = 1.3785363651718823E-12;
      Cs[1637] = -1.0618921845918994E-28;
      Cs[1638] = -1.3785363651718827E-12;
      Cs[1639] = 4.3706351924459388E-17;
      Cs[1640] = 1.009741958682895E-32;
      Cs[1641] = -4.3706351924459369E-17;
      Cs[1642] = 1.4901813611218238E-33;
      Cs[1643] = 2.34381396951068E-13;
      Cs[1644] = -3.1057701708496958E-29;
      Cs[1645] = -2.3438139695106811E-13;
      Cs[1646] = 3.0969451882490318E-29;
      Cs[1647] = 0.0;
      Cs[1648] = 0.0;
      Cs[1649] = 0.0;
      Cs[1650] = 0.0;
      Cs[1651] = 0.0;
      Cs[1652] = 0.0;
      Cs[1653] = 2.3320317795248351E-9;
      Cs[1654] = -3.507923976556525E-11;
      Cs[1655] = 1.614909251497046E-13;
      Cs[1656] = -2.8356053055635841E-19;
      Cs[1657] = 1.6149092514970455E-13;
      Cs[1658] = -2.8356053055607732E-19;
      Cs[1659] = -1.4204983082856464E-16;
      Cs[1660] = -3.1894011042063221E-16;
      Cs[1661] = -1.4204983082856473E-16;
      Cs[1662] = -3.1894011042063236E-16;
      Cs[1663] = 2.8409966165712942E-16;
      Cs[1664] = 6.3788022084126433E-16;
      Cs[1665] = 0.0;
      Cs[1666] = 0.0;
      Cs[1667] = 2.1259022226465772E-17;
      Cs[1668] = 2.3999207555933728E-33;
      Cs[1669] = -7.6019014100192493E-19;
      Cs[1670] = -1.9999336168423298E-11;
      Cs[1671] = -7.6019014100175141E-19;
      Cs[1672] = 1.7209996301487208E-17;
      Cs[1673] = -2.0550664363147606E-16;
      Cs[1674] = -3.4419992602974422E-17;
      Cs[1675] = -4.1715356017051804E-32;
      Cs[1676] = 2.0550664363147616E-16;
      Cs[1677] = 1.7209996301487211E-17;
      Cs[1678] = -6.3786915463167562E-16;
      Cs[1679] = 3.1893457731583786E-16;
      Cs[1680] = 3.1893457731583786E-16;
      Cs[1681] = -1.3785363651718817E-12;
      Cs[1682] = 1.3785363651718833E-12;
      Cs[1683] = 7.4187502738565965E-28;
      Cs[1684] = -4.3706351924459375E-17;
      Cs[1685] = 4.3706351924459388E-17;
      Cs[1686] = 2.8272774843121063E-33;
      Cs[1687] = 9.592140310339293E-33;
      Cs[1688] = -2.3438139695106811E-13;
      Cs[1689] = 2.3438139695106806E-13;
      Cs[1690] = -1.531940110035341E-29;
      Cs[1691] = 3.363892241224567E-28;
      Cs[1692] = 0.0;
      Cs[1693] = 0.0;
      Cs[1694] = 0.0;
      Cs[1695] = 0.0;
      Cs[1696] = 0.0;
      Cs[1697] = 0.0;
      Cs[1698] = 1.6149092514970468E-13;
      Cs[1699] = -2.8356053055638142E-19;
      Cs[1700] = 2.3320317795248355E-9;
      Cs[1701] = -3.5079239765565257E-11;
      Cs[1702] = 1.6149092514970455E-13;
      Cs[1703] = -2.8356053055607732E-19;
      Cs[1704] = 2.8409966165712942E-16;
      Cs[1705] = 6.3788022084126482E-16;
      Cs[1706] = -1.4204983082856473E-16;
      Cs[1707] = -3.1894011042063231E-16;
      Cs[1708] = -1.4204983082856473E-16;
      Cs[1709] = -3.1894011042063226E-16;
      Cs[1710] = 0.0;
      Cs[1711] = 0.0;
      Cs[1712] = 3.2581919834171061E-37;
      Cs[1713] = 2.1259022226465753E-17;
      Cs[1714] = -7.6019014100205792E-19;
      Cs[1715] = -7.6019014100236443E-19;
      Cs[1716] = -1.9999336168423295E-11;
      Cs[1717] = 1.72099963014872E-17;
      Cs[1718] = -6.1411445587103493E-32;
      Cs[1719] = 1.7209996301487208E-17;
      Cs[1720] = -2.0550664363147592E-16;
      Cs[1721] = -2.0550664363147619E-16;
      Cs[1722] = -3.441999260297441E-17;
      Cs[1723] = 3.1893457731583781E-16;
      Cs[1724] = -6.3786915463167572E-16;
      Cs[1725] = 3.1893457731583771E-16;
      Cs[1726] = -4.6339945379143054E-28;
      Cs[1727] = -1.3785363651718831E-12;
      Cs[1728] = 1.3785363651718812E-12;
      Cs[1729] = 2.3505764316289037E-33;
      Cs[1730] = -4.3706351924459381E-17;
      Cs[1731] = 4.3706351924459357E-17;
      Cs[1732] = -2.58371768130869E-32;
      Cs[1733] = -9.4915744116192135E-33;
      Cs[1734] = -2.3438139695106811E-13;
      Cs[1735] = 2.3438139695106796E-13;
      Cs[1736] = -1.8735641989439513E-28;
      Cs[1737] = 0.0;
      Cs[1738] = 0.0;
      Cs[1739] = 0.0;
      Cs[1740] = 0.0;
      Cs[1741] = 0.0;
      Cs[1742] = 0.0;
      Cs[1743] = 1.614909251497046E-13;
      Cs[1744] = -2.8356053055661528E-19;
      Cs[1745] = 1.6149092514970468E-13;
      Cs[1746] = -2.8356053055715295E-19;
      Cs[1747] = 2.3320317795248351E-9;
      Cs[1748] = -3.507923976556525E-11;
      Cs[1749] = -1.4204983082856471E-16;
      Cs[1750] = -3.1894011042063236E-16;
      Cs[1751] = 2.8409966165712947E-16;
      Cs[1752] = 6.3788022084126462E-16;
      Cs[1753] = -1.4204983082856466E-16;
      Cs[1754] = -3.18940110420632E-16;
      Cs[1755] = 0.0;
      Cs[1756] = 0.0;
      Cs[1757] = -1.4215583765083687E-25;
      Cs[1758] = -7.1077918841153746E-26;
      Cs[1759] = -6.6866252058226807E-20;
      Cs[1760] = 6.68662520584532E-20;
      Cs[1761] = 1.4916166960040638E-29;
      Cs[1762] = -1.3365298547930745E-7;
      Cs[1763] = 1.7293942197000577E-6;
      Cs[1764] = 1.3365298547930745E-7;
      Cs[1765] = -1.3342778391381919E-7;
      Cs[1766] = 1.3342778391381919E-7;
      Cs[1767] = -1.3552527156068805E-26;
      Cs[1768] = 2.4768487850958532E-6;
      Cs[1769] = 1.3877787807814456E-22;
      Cs[1770] = -2.476848785095854E-6;
      Cs[1771] = 2.7474816308652691E-9;
      Cs[1772] = -3.0616506057155405E-10;
      Cs[1773] = -3.0616506057155482E-10;
      Cs[1774] = -3.9963349367587549E-10;
      Cs[1775] = 1.1955350191719026E-10;
      Cs[1776] = 1.1955350191719029E-10;
      Cs[1777] = -1.6052648984151441E-10;
      Cs[1778] = -7.1777028764573168E-14;
      Cs[1779] = 2.5166494610178103E-14;
      Cs[1780] = 2.5166494610117802E-14;
      Cs[1781] = -2.1444039562829117E-14;
      Cs[1782] = 0.0;
      Cs[1783] = 0.0;
      Cs[1784] = 0.0;
      Cs[1785] = 0.0;
      Cs[1786] = 0.0;
      Cs[1787] = 0.0;
      Cs[1788] = 6.6794973871717143E-14;
      Cs[1789] = -1.1728476078505118E-19;
      Cs[1790] = -6.6794973871821247E-14;
      Cs[1791] = 1.1728476078544831E-19;
      Cs[1792] = -1.4903905357686488E-23;
      Cs[1793] = 2.6163253497476967E-29;
      Cs[1794] = -1.1031602589843443E-6;
      Cs[1795] = -2.4768917552372647E-6;
      Cs[1796] = -9.0205620750793963E-23;
      Cs[1797] = -1.2490009027033011E-22;
      Cs[1798] = 1.1031602589843445E-6;
      Cs[1799] = 2.4768917552372635E-6;
      Cs[1800] = 0.0;
      Cs[1801] = 0.0;
      Cs[1802] = 7.1077918790627233E-26;
      Cs[1803] = -7.1077918845551724E-26;
      Cs[1804] = -1.8549780196346213E-30;
      Cs[1805] = -6.6866252027465324E-20;
      Cs[1806] = 6.6866252075425385E-20;
      Cs[1807] = 4.0657581468206417E-26;
      Cs[1808] = 1.3342778391381946E-7;
      Cs[1809] = -1.3365298547930742E-7;
      Cs[1810] = -1.3342778391381946E-7;
      Cs[1811] = 1.7293942197000579E-6;
      Cs[1812] = 1.3365298547930742E-7;
      Cs[1813] = -2.476848785095854E-6;
      Cs[1814] = 2.4768487850958553E-6;
      Cs[1815] = 2.7755575615628912E-23;
      Cs[1816] = -3.0616506057172165E-10;
      Cs[1817] = 2.7474816308646739E-9;
      Cs[1818] = -3.0616506057174213E-10;
      Cs[1819] = 1.1955350191720073E-10;
      Cs[1820] = -3.9963349367591132E-10;
      Cs[1821] = 1.1955350191720233E-10;
      Cs[1822] = -1.6052648984144728E-10;
      Cs[1823] = 2.5166494623177237E-14;
      Cs[1824] = -7.1777028883625608E-14;
      Cs[1825] = 2.5166494621320656E-14;
      Cs[1826] = -2.1444039562677324E-14;
      Cs[1827] = 0.0;
      Cs[1828] = 0.0;
      Cs[1829] = 0.0;
      Cs[1830] = 0.0;
      Cs[1831] = 0.0;
      Cs[1832] = 0.0;
      Cs[1833] = 1.8754798377312806E-24;
      Cs[1834] = -3.2536725808021672E-30;
      Cs[1835] = 6.6794973840834538E-14;
      Cs[1836] = -1.1728476073109491E-19;
      Cs[1837] = -6.6794973888929545E-14;
      Cs[1838] = 1.1728476081521784E-19;
      Cs[1839] = 1.1031602589843447E-6;
      Cs[1840] = 2.4768917552372656E-6;
      Cs[1841] = -1.1031602589843452E-6;
      Cs[1842] = -2.4768917552372652E-6;
      Cs[1843] = 2.7755575615628912E-23;
      Cs[1844] = -4.1633363423443371E-23;
      Cs[1845] = 0.0;
      Cs[1846] = 0.0;
      Cs[1847] = 7.1077918790473443E-26;
      Cs[1848] = 1.4215583758242502E-25;
      Cs[1849] = 6.6866252027378886E-20;
      Cs[1850] = 1.9132338197633825E-30;
      Cs[1851] = -6.6866252024943042E-20;
      Cs[1852] = 1.3365298547930742E-7;
      Cs[1853] = 1.3342778391381946E-7;
      Cs[1854] = -3.4829994791096832E-24;
      Cs[1855] = -1.7293942197000579E-6;
      Cs[1856] = 1.3342778391381949E-7;
      Cs[1857] = -1.3365298547930742E-7;
      Cs[1858] = -8.3266726846886743E-23;
      Cs[1859] = -2.4768487850958553E-6;
      Cs[1860] = 2.4768487850958549E-6;
      Cs[1861] = -3.0616506057172279E-10;
      Cs[1862] = -3.0616506057174326E-10;
      Cs[1863] = 2.7474816308646735E-9;
      Cs[1864] = 1.1955350191720422E-10;
      Cs[1865] = 1.1955350191720923E-10;
      Cs[1866] = -3.9963349367588448E-10;
      Cs[1867] = -1.6052648984144723E-10;
      Cs[1868] = 2.5166494623188972E-14;
      Cs[1869] = 2.5166494624736629E-14;
      Cs[1870] = -7.1777028874731788E-14;
      Cs[1871] = -2.1444039563059755E-14;
      Cs[1872] = 0.0;
      Cs[1873] = 0.0;
      Cs[1874] = 0.0;
      Cs[1875] = 0.0;
      Cs[1876] = 0.0;
      Cs[1877] = 0.0;
      Cs[1878] = -6.6794973840891917E-14;
      Cs[1879] = 1.1728476073094331E-19;
      Cs[1880] = -1.8946340319955583E-24;
      Cs[1881] = 3.3558521561336349E-30;
      Cs[1882] = 6.6794973838496746E-14;
      Cs[1883] = -1.1728476072667078E-19;
      Cs[1884] = 2.0816681711721686E-23;
      Cs[1885] = 8.3266726846886743E-23;
      Cs[1886] = 1.1031602589843452E-6;
      Cs[1887] = 2.4768917552372652E-6;
      Cs[1888] = -1.103160258984345E-6;
      Cs[1889] = -2.4768917552372643E-6;
      Cs[1890] = 0.0;
      Cs[1891] = 0.0;
      Cs[1892] = 0.0;
      Cs[1893] = 0.0;
      Cs[1894] = 0.0;
      Cs[1895] = 0.0;
      Cs[1896] = 0.0;
      Cs[1897] = 0.0;
      Cs[1898] = 0.0;
      Cs[1899] = 0.0;
      Cs[1900] = 0.0;
      Cs[1901] = 0.0;
      Cs[1902] = 0.0;
      Cs[1903] = 0.0;
      Cs[1904] = 0.0;
      Cs[1905] = 0.0;
      Cs[1906] = 0.0;
      Cs[1907] = 0.0;
      Cs[1908] = 0.0;
      Cs[1909] = 0.0;
      Cs[1910] = 0.0;
      Cs[1911] = 0.0;
      Cs[1912] = 0.0;
      Cs[1913] = 0.0;
      Cs[1914] = 0.0;
      Cs[1915] = 0.0;
      Cs[1916] = 0.0;
      Cs[1917] = 0.0;
      Cs[1918] = 0.0;
      Cs[1919] = 0.0;
      Cs[1920] = 0.0;
      Cs[1921] = 0.0;
      Cs[1922] = 0.0;
      Cs[1923] = 0.0;
      Cs[1924] = 0.0;
      Cs[1925] = 0.0;
      Cs[1926] = 0.0;
      Cs[1927] = 0.0;
      Cs[1928] = 0.0;
      Cs[1929] = 0.0;
      Cs[1930] = 0.0;
      Cs[1931] = 0.0;
      Cs[1932] = 0.0;
      Cs[1933] = 0.0;
      Cs[1934] = 0.0;
      Cs[1935] = 0.0;
      Cs[1936] = 0.0;
      Cs[1937] = 0.0;
      Cs[1938] = 0.0;
      Cs[1939] = 0.0;
      Cs[1940] = 0.0;
      Cs[1941] = 0.0;
      Cs[1942] = 0.0;
      Cs[1943] = 0.0;
      Cs[1944] = 0.0;
      Cs[1945] = 0.0;
      Cs[1946] = 0.0;
      Cs[1947] = 0.0;
      Cs[1948] = 0.0;
      Cs[1949] = 0.0;
      Cs[1950] = 0.0;
      Cs[1951] = 0.0;
      Cs[1952] = 0.0;
      Cs[1953] = 0.0;
      Cs[1954] = 0.0;
      Cs[1955] = 0.0;
      Cs[1956] = 0.0;
      Cs[1957] = 0.0;
      Cs[1958] = 0.0;
      Cs[1959] = 0.0;
      Cs[1960] = 0.0;
      Cs[1961] = 0.0;
      Cs[1962] = 0.0;
      Cs[1963] = 0.0;
      Cs[1964] = 0.0;
      Cs[1965] = 0.0;
      Cs[1966] = 0.0;
      Cs[1967] = 0.0;
      Cs[1968] = 0.0;
      Cs[1969] = 0.0;
      Cs[1970] = 0.0;
      Cs[1971] = 0.0;
      Cs[1972] = 0.0;
      Cs[1973] = 0.0;
      Cs[1974] = 0.0;
      Cs[1975] = 0.0;
      Cs[1976] = 0.0;
      Cs[1977] = 0.0;
      Cs[1978] = 0.0;
      Cs[1979] = 0.0;
      Cs[1980] = 0.0;
      Cs[1981] = 0.0;
      Cs[1982] = 0.0;
      Cs[1983] = 0.0;
      Cs[1984] = 0.0;
      Cs[1985] = 0.0;
      Cs[1986] = 0.0;
      Cs[1987] = 0.0;
      Cs[1988] = 0.0;
      Cs[1989] = 0.0;
      Cs[1990] = 0.0;
      Cs[1991] = 0.0;
      Cs[1992] = 0.0;
      Cs[1993] = 0.0;
      Cs[1994] = 0.0;
      Cs[1995] = 0.0;
      Cs[1996] = 0.0;
      Cs[1997] = 0.0;
      Cs[1998] = 0.0;
      Cs[1999] = 0.0;
      Cs[2000] = 0.0;
      Cs[2001] = 0.0;
      Cs[2002] = 0.0;
      Cs[2003] = 0.0;
      Cs[2004] = 0.0;
      Cs[2005] = 0.0;
      Cs[2006] = 0.0;
      Cs[2007] = 0.0;
      Cs[2008] = 0.0;
      Cs[2009] = 0.0;
      Cs[2010] = 0.0;
      Cs[2011] = 0.0;
      Cs[2012] = 0.0;
      Cs[2013] = 0.0;
      Cs[2014] = 0.0;
      Cs[2015] = 0.0;
      Cs[2016] = 0.0;
      Cs[2017] = 0.0;
      Cs[2018] = 0.0;
      Cs[2019] = 0.0;
      Cs[2020] = 0.0;
      Cs[2021] = 0.0;
      Cs[2022] = 0.0;
      Cs[2023] = 0.0;
      Cs[2024] = 0.0;
      Cs[2025] = 0.0;
      Cs[2026] = 0.0;
      Cs[2027] = -2.3676466022461934E-10;
      Cs[2028] = -2.3676466022461923E-10;
      Cs[2029] = 1.3331752496216813E-11;
      Cs[2030] = -6.66587624810841E-12;
      Cs[2031] = -6.6658762481084073E-12;
      Cs[2032] = -3.0902962616933353E-24;
      Cs[2033] = 1.8450800379097553E-23;
      Cs[2034] = 1.5451481308466677E-24;
      Cs[2035] = 1.8450800379097541E-23;
      Cs[2036] = 6.1629758220391545E-39;
      Cs[2037] = 1.5451481308466677E-24;
      Cs[2038] = 2.8634588722099446E-23;
      Cs[2039] = 2.8634588722099457E-23;
      Cs[2040] = -5.7269177444198927E-23;
      Cs[2041] = 1.2376777139489712E-19;
      Cs[2042] = 0.0;
      Cs[2043] = -1.2376777139489712E-19;
      Cs[2044] = 3.9240443053650859E-24;
      Cs[2045] = 1.5407439555097886E-39;
      Cs[2046] = -3.9240443053650859E-24;
      Cs[2047] = 2.40500670782669E-39;
      Cs[2048] = 2.1043233889183238E-20;
      Cs[2049] = 3.1554436208840471E-36;
      Cs[2050] = -2.1043233889183229E-20;
      Cs[2051] = -1.128165341042545E-35;
      Cs[2052] = 0.0;
      Cs[2053] = 0.0;
      Cs[2054] = 0.0;
      Cs[2055] = 0.0;
      Cs[2056] = 0.0;
      Cs[2057] = 0.0;
      Cs[2058] = 1.3957334077580342E-16;
      Cs[2059] = -2.0996565161941497E-18;
      Cs[2060] = -6.9786670387901684E-17;
      Cs[2061] = 1.0498282580970747E-18;
      Cs[2062] = -6.9786670387901684E-17;
      Cs[2063] = 1.0498282580970749E-18;
      Cs[2064] = -1.2753519916379937E-23;
      Cs[2065] = -2.8635085495392195E-23;
      Cs[2066] = -1.2753519916379939E-23;
      Cs[2067] = -2.8635085495392219E-23;
      Cs[2068] = 2.5507039832759885E-23;
      Cs[2069] = 5.7270170990784414E-23;
      Cs[2070] = 0.0;
      Cs[2071] = 0.0;
      Cs[2072] = 2.3676466022461934E-10;
      Cs[2073] = 2.6730301104972933E-26;
      Cs[2074] = -6.6658762481084073E-12;
      Cs[2075] = 1.3331752496216816E-11;
      Cs[2076] = -6.6658762481084073E-12;
      Cs[2077] = 1.545148130846668E-24;
      Cs[2078] = -1.8450800379097553E-23;
      Cs[2079] = -3.0902962616933368E-24;
      Cs[2080] = 2.9129833672230946E-39;
      Cs[2081] = 1.8450800379097553E-23;
      Cs[2082] = 1.5451481308466686E-24;
      Cs[2083] = -5.7269177444198915E-23;
      Cs[2084] = 2.8634588722099481E-23;
      Cs[2085] = 2.8634588722099469E-23;
      Cs[2086] = -1.2376777139489712E-19;
      Cs[2087] = 1.2376777139489719E-19;
      Cs[2088] = 2.6139936544608446E-35;
      Cs[2089] = -3.9240443053650866E-24;
      Cs[2090] = 3.9240443053650873E-24;
      Cs[2091] = 9.6296497219361788E-41;
      Cs[2092] = -1.6358376832359691E-39;
      Cs[2093] = -2.1043233889183238E-20;
      Cs[2094] = 2.1043233889183241E-20;
      Cs[2095] = -3.0067377550530475E-37;
      Cs[2096] = 1.1281654707601834E-35;
      Cs[2097] = 0.0;
      Cs[2098] = 0.0;
      Cs[2099] = 0.0;
      Cs[2100] = 0.0;
      Cs[2101] = 0.0;
      Cs[2102] = 0.0;
      Cs[2103] = -6.9786670387901708E-17;
      Cs[2104] = 1.0498282580970749E-18;
      Cs[2105] = 1.3957334077580339E-16;
      Cs[2106] = -2.0996565161941493E-18;
      Cs[2107] = -6.9786670387901708E-17;
      Cs[2108] = 1.0498282580970749E-18;
      Cs[2109] = 2.5507039832759882E-23;
      Cs[2110] = 5.7270170990784426E-23;
      Cs[2111] = -1.2753519916379951E-23;
      Cs[2112] = -2.8635085495392219E-23;
      Cs[2113] = -1.2753519916379944E-23;
      Cs[2114] = -2.8635085495392213E-23;
      Cs[2115] = 0.0;
      Cs[2116] = 0.0;
      Cs[2117] = 2.5243645263569593E-35;
      Cs[2118] = 2.3676466022461923E-10;
      Cs[2119] = -6.6658762481084057E-12;
      Cs[2120] = -6.6658762481084065E-12;
      Cs[2121] = 1.3331752496216815E-11;
      Cs[2122] = 1.5451481308466675E-24;
      Cs[2123] = -9.9298966930751153E-41;
      Cs[2124] = 1.5451481308466689E-24;
      Cs[2125] = -1.8450800379097544E-23;
      Cs[2126] = -1.8450800379097559E-23;
      Cs[2127] = -3.0902962616933361E-24;
      Cs[2128] = 2.8634588722099469E-23;
      Cs[2129] = -5.7269177444198938E-23;
      Cs[2130] = 2.8634588722099457E-23;
      Cs[2131] = -8.9638005363754037E-37;
      Cs[2132] = -1.2376777139489719E-19;
      Cs[2133] = 1.2376777139489709E-19;
      Cs[2134] = 5.0992064031748831E-40;
      Cs[2135] = -3.9240443053650888E-24;
      Cs[2136] = 3.9240443053650859E-24;
      Cs[2137] = -7.6916902459072074E-40;
      Cs[2138] = -1.523169844767255E-37;
      Cs[2139] = -2.1043233889183244E-20;
      Cs[2140] = 2.1043233889183229E-20;
      Cs[2141] = -1.2971763832316257E-42;
      Cs[2142] = 0.0;
      Cs[2143] = 0.0;
      Cs[2144] = 0.0;
      Cs[2145] = 0.0;
      Cs[2146] = 0.0;
      Cs[2147] = 0.0;
      Cs[2148] = -6.9786670387901708E-17;
      Cs[2149] = 1.0498282580970749E-18;
      Cs[2150] = -6.9786670387901708E-17;
      Cs[2151] = 1.0498282580970747E-18;
      Cs[2152] = 1.3957334077580339E-16;
      Cs[2153] = -2.0996565161941497E-18;
      Cs[2154] = -1.2753519916379945E-23;
      Cs[2155] = -2.8635085495392237E-23;
      Cs[2156] = 2.5507039832759891E-23;
      Cs[2157] = 5.7270170990784438E-23;
      Cs[2158] = -1.2753519916379941E-23;
      Cs[2159] = -2.86350854953922E-23;
      Cs[2160] = 1.9999999999666666E-6;
      Cs[2161] = 0.0;
      Cs[2162] = 0.0;
      Cs[2163] = 0.0;
      Cs[2164] = 0.0;
      Cs[2165] = 0.0;
      Cs[2166] = 0.0;
      Cs[2167] = 0.0;
      Cs[2168] = 0.0;
      Cs[2169] = 0.0;
      Cs[2170] = 0.0;
      Cs[2171] = 0.0;
      Cs[2172] = 0.0;
      Cs[2173] = 0.0;
      Cs[2174] = 0.0;
      Cs[2175] = 0.0;
      Cs[2176] = 0.0;
      Cs[2177] = 0.0;
      Cs[2178] = 0.0;
      Cs[2179] = 0.0;
      Cs[2180] = 0.0;
      Cs[2181] = 0.0;
      Cs[2182] = 0.0;
      Cs[2183] = 0.0;
      Cs[2184] = 0.0;
      Cs[2185] = 0.0;
      Cs[2186] = 0.0;
      Cs[2187] = 0.0;
      Cs[2188] = 0.0;
      Cs[2189] = 0.0;
      Cs[2190] = 0.0;
      Cs[2191] = 0.0;
      Cs[2192] = 0.0;
      Cs[2193] = 0.0;
      Cs[2194] = 0.0;
      Cs[2195] = 0.0;
      Cs[2196] = 0.0;
      Cs[2197] = 0.0;
      Cs[2198] = 0.0;
      Cs[2199] = 0.0;
      Cs[2200] = 0.0;
      Cs[2201] = 0.0;
      Cs[2202] = 0.0;
      Cs[2203] = 0.0;
      Cs[2204] = 0.0;
      Cs[2205] = 0.0;
      Cs[2206] = 2.1794285528763532E-9;
      Cs[2207] = 0.0;
      Cs[2208] = 0.0;
      Cs[2209] = 0.0;
      Cs[2210] = 0.0;
      Cs[2211] = 0.0;
      Cs[2212] = 0.0;
      Cs[2213] = 0.0;
      Cs[2214] = 0.0;
      Cs[2215] = 0.0;
      Cs[2216] = 0.0;
      Cs[2217] = 0.0;
      Cs[2218] = 0.0;
      Cs[2219] = 0.0;
      Cs[2220] = 0.0;
      Cs[2221] = 0.0;
      Cs[2222] = 0.0;
      Cs[2223] = 0.0;
      Cs[2224] = 0.0;
      Cs[2225] = 0.0;
      Cs[2226] = 0.0;
      Cs[2227] = 0.0;
      Cs[2228] = 0.0;
      Cs[2229] = 0.0;
      Cs[2230] = 0.0;
      Cs[2231] = 0.0;
      Cs[2232] = -3.0269841012171559E-12;
      Cs[2233] = -3.0269841012171583E-12;
      Cs[2234] = -3.0269841012171579E-12;
      Cs[2235] = -3.0269841012171555E-12;
      Cs[2236] = -3.0269841012171571E-12;
      Cs[2237] = -3.0269841012171579E-12;
      Cs[2238] = 0.0;
      Cs[2239] = 0.0;
      Cs[2240] = 0.0;
      Cs[2241] = 0.0;
      Cs[2242] = 0.0;
      Cs[2243] = 0.0;
      Cs[2244] = 0.0;
      Cs[2245] = 0.0;
      Cs[2246] = 0.0;
      Cs[2247] = 0.0;
      Cs[2248] = 0.0;
      Cs[2249] = 0.0;
      Ds[0] = -3.3131726142746629;
      Ds[1] = 0.00080831607982809466;
      Ds[2] = 0.00080831607982729333;
      Ds[3] = 0.0;
      Ds[4] = 0.0;
      Ds[5] = 0.0;
      Ds[6] = 0.0;
      Ds[7] = 0.0;
      Ds[8] = 0.0;
      Ds[9] = 0.0;
      Ds[10] = 0.0;
      Ds[11] = 0.0;
      Ds[12] = 0.0;
      Ds[13] = -0.0059852692503629262;
      Ds[14] = 0.0059852692503629262;
      Ds[15] = 0.0029926346251814635;
      Ds[16] = -0.0029926346251814635;
      Ds[17] = 0.0029926346251814631;
      Ds[18] = -0.0029926346251814631;
      Ds[19] = 0.0;
      Ds[20] = 0.0;
      Ds[21] = 1.9860087161920738E-10;
      Ds[22] = 1.9860087161920762E-10;
      Ds[23] = -3.9720174323841508E-10;
      Ds[24] = 0.00080831607982816036;
      Ds[25] = -3.313172614260111;
      Ds[26] = 0.00080831607982729333;
      Ds[27] = 0.0;
      Ds[28] = 0.0;
      Ds[29] = 0.0;
      Ds[30] = 0.0;
      Ds[31] = 0.0;
      Ds[32] = 0.0;
      Ds[33] = 0.0;
      Ds[34] = 0.0;
      Ds[35] = 0.0;
      Ds[36] = 0.0;
      Ds[37] = 0.0029926346251814653;
      Ds[38] = -0.0029926346251814653;
      Ds[39] = -0.0059852692503629306;
      Ds[40] = 0.0059852692503629306;
      Ds[41] = 0.0029926346251814653;
      Ds[42] = -0.0029926346251814653;
      Ds[43] = 0.0;
      Ds[44] = 0.0;
      Ds[45] = -3.9720174323841497E-10;
      Ds[46] = 1.9860087161920756E-10;
      Ds[47] = 1.9860087161920762E-10;
      Ds[48] = 0.00080831607982882682;
      Ds[49] = 0.00080831607983035945;
      Ds[50] = -3.3131726142746629;
      Ds[51] = 0.0;
      Ds[52] = 0.0;
      Ds[53] = 0.0;
      Ds[54] = 0.0;
      Ds[55] = 0.0;
      Ds[56] = 0.0;
      Ds[57] = 0.0;
      Ds[58] = 0.0;
      Ds[59] = 0.0;
      Ds[60] = 0.0;
      Ds[61] = 0.0029926346251814631;
      Ds[62] = -0.0029926346251814631;
      Ds[63] = 0.0029926346251814631;
      Ds[64] = -0.0029926346251814631;
      Ds[65] = -0.0059852692503629262;
      Ds[66] = 0.0059852692503629262;
      Ds[67] = 0.0;
      Ds[68] = 0.0;
      Ds[69] = 1.9860087161920751E-10;
      Ds[70] = -3.9720174323841513E-10;
      Ds[71] = 1.9860087161920749E-10;
      Ds[72] = 0.0;
      Ds[73] = 0.0;
      Ds[74] = 0.0;
      Ds[75] = -54595.23848174265;
      Ds[76] = -9080.9523036514711;
      Ds[77] = -45404.76151825735;
      Ds[78] = -9080.9523036514711;
      Ds[79] = -0.276389071193217;
      Ds[80] = -0.27638907119321704;
      Ds[81] = -0.27638907119321704;
      Ds[82] = -0.2763890711932171;
      Ds[83] = -0.27638907119321704;
      Ds[84] = -0.27638907119321704;
      Ds[85] = 0.0;
      Ds[86] = 0.0;
      Ds[87] = 0.0;
      Ds[88] = 0.0;
      Ds[89] = 0.0;
      Ds[90] = 0.0;
      Ds[91] = -4.9510097955490009E-17;
      Ds[92] = 2.6118815047155281E-19;
      Ds[93] = 0.0;
      Ds[94] = 0.0;
      Ds[95] = 0.0;
      Ds[96] = 0.0;
      Ds[97] = 0.0;
      Ds[98] = 0.0;
      Ds[99] = -9080.9523036514711;
      Ds[100] = -18183.809539269707;
      Ds[101] = 9080.9523036514711;
      Ds[102] = -18183.809539269707;
      Ds[103] = 0.055277814238643364;
      Ds[104] = 0.055277814238643377;
      Ds[105] = 0.055277814238643377;
      Ds[106] = 0.055277814238643384;
      Ds[107] = 0.055277814238643377;
      Ds[108] = 0.055277814238643377;
      Ds[109] = 0.0;
      Ds[110] = 0.0;
      Ds[111] = 0.0;
      Ds[112] = 0.0;
      Ds[113] = 0.0;
      Ds[114] = 0.0;
      Ds[115] = 1.0029429994202406E-17;
      Ds[116] = -5.7333251577734961E-20;
      Ds[117] = 0.0;
      Ds[118] = 0.0;
      Ds[119] = 0.0;
      Ds[120] = 0.0;
      Ds[121] = 0.0;
      Ds[122] = 0.0;
      Ds[123] = -45404.76151825735;
      Ds[124] = 9080.9523036514711;
      Ds[125] = -54595.23848174265;
      Ds[126] = 9080.9523036514711;
      Ds[127] = 0.276389071193217;
      Ds[128] = 0.27638907119321704;
      Ds[129] = 0.27638907119321704;
      Ds[130] = 0.2763890711932171;
      Ds[131] = 0.27638907119321704;
      Ds[132] = 0.27638907119321704;
      Ds[133] = 0.0;
      Ds[134] = 0.0;
      Ds[135] = 0.0;
      Ds[136] = 0.0;
      Ds[137] = 0.0;
      Ds[138] = 0.0;
      Ds[139] = 4.9510097955490009E-17;
      Ds[140] = -2.6118815047155281E-19;
      Ds[141] = 0.0;
      Ds[142] = 0.0;
      Ds[143] = 0.0;
      Ds[144] = 0.0;
      Ds[145] = 0.0;
      Ds[146] = 0.0;
      Ds[147] = -9080.9523036514711;
      Ds[148] = -18183.809539269707;
      Ds[149] = 9080.9523036514711;
      Ds[150] = -18183.80954149193;
      Ds[151] = 0.055277814238643364;
      Ds[152] = 0.055277814238643377;
      Ds[153] = 0.055277814238643377;
      Ds[154] = 0.055277814238643384;
      Ds[155] = 0.055277814238643377;
      Ds[156] = 0.055277814238643377;
      Ds[157] = 1.1111111110925924E-6;
      Ds[158] = 1.1111111110925926E-6;
      Ds[159] = 1.1111111110925924E-6;
      Ds[160] = 1.1111111110925926E-6;
      Ds[161] = 1.1111111110925924E-6;
      Ds[162] = 1.1111111110925926E-6;
      Ds[163] = 1.0029429994202406E-17;
      Ds[164] = -5.7333251577734961E-20;
      Ds[165] = 0.0;
      Ds[166] = 0.0;
      Ds[167] = 0.0;
      Ds[168] = 0.0;
      Ds[169] = 0.0;
      Ds[170] = 0.0;
      Ds[171] = -0.276389071193217;
      Ds[172] = 0.055277814238643412;
      Ds[173] = 0.276389071193217;
      Ds[174] = 0.055277814238643412;
      Ds[175] = -1.2174436336483658;
      Ds[176] = -0.30435964657991588;
      Ds[177] = -0.30435964657991588;
      Ds[178] = 0.608724340488534;
      Ds[179] = -0.30435964657991604;
      Ds[180] = -0.30435964657991604;
      Ds[181] = 0.0;
      Ds[182] = 0.0;
      Ds[183] = 0.0;
      Ds[184] = 0.0;
      Ds[185] = 0.0;
      Ds[186] = 0.0;
      Ds[187] = 0.91308398706844984;
      Ds[188] = -2.221250989929003E-16;
      Ds[189] = 0.0;
      Ds[190] = 0.0;
      Ds[191] = 0.0;
      Ds[192] = 0.0;
      Ds[193] = 0.0;
      Ds[194] = 0.0;
      Ds[195] = -0.27638907119321693;
      Ds[196] = 0.055277814238643384;
      Ds[197] = 0.27638907119321693;
      Ds[198] = 0.055277814238643384;
      Ds[199] = -0.30435964657991577;
      Ds[200] = -1.2174436336483654;
      Ds[201] = -0.30435964657991582;
      Ds[202] = -0.30435964657991588;
      Ds[203] = 0.60872434048853374;
      Ds[204] = -0.30435964657991582;
      Ds[205] = 0.0;
      Ds[206] = 0.0;
      Ds[207] = 0.0;
      Ds[208] = 0.0;
      Ds[209] = 0.0;
      Ds[210] = 0.0;
      Ds[211] = 0.91308398706844962;
      Ds[212] = 0.91308398706844962;
      Ds[213] = 0.0;
      Ds[214] = 0.0;
      Ds[215] = 0.0;
      Ds[216] = 0.0;
      Ds[217] = 0.0;
      Ds[218] = 0.0;
      Ds[219] = -0.276389071193217;
      Ds[220] = 0.0552778142386434;
      Ds[221] = 0.276389071193217;
      Ds[222] = 0.0552778142386434;
      Ds[223] = -0.30435964657991582;
      Ds[224] = -0.30435964657991588;
      Ds[225] = -1.2174436336483656;
      Ds[226] = -0.30435964657991593;
      Ds[227] = -0.30435964657991588;
      Ds[228] = 0.60872434048853385;
      Ds[229] = 0.0;
      Ds[230] = 0.0;
      Ds[231] = 0.0;
      Ds[232] = 0.0;
      Ds[233] = 0.0;
      Ds[234] = 0.0;
      Ds[235] = -5.4815329008372074E-17;
      Ds[236] = 0.91308398706844973;
      Ds[237] = 0.0;
      Ds[238] = 0.0;
      Ds[239] = 0.0;
      Ds[240] = 0.0;
      Ds[241] = 0.0;
      Ds[242] = 0.0;
      Ds[243] = -0.27638907119321693;
      Ds[244] = 0.055277814238643384;
      Ds[245] = 0.27638907119321693;
      Ds[246] = 0.055277814238643384;
      Ds[247] = 0.60872434048853419;
      Ds[248] = -0.30435964657991582;
      Ds[249] = -0.30435964657991582;
      Ds[250] = -1.2174436336483658;
      Ds[251] = -0.30435964657991582;
      Ds[252] = -0.30435964657991582;
      Ds[253] = 0.0;
      Ds[254] = 0.0;
      Ds[255] = 0.0;
      Ds[256] = 0.0;
      Ds[257] = 0.0;
      Ds[258] = 0.0;
      Ds[259] = -0.91308398706845006;
      Ds[260] = 6.0798828768568259E-20;
      Ds[261] = 0.0;
      Ds[262] = 0.0;
      Ds[263] = 0.0;
      Ds[264] = 0.0;
      Ds[265] = 0.0;
      Ds[266] = 0.0;
      Ds[267] = -0.276389071193217;
      Ds[268] = 0.055277814238643412;
      Ds[269] = 0.276389071193217;
      Ds[270] = 0.055277814238643412;
      Ds[271] = -0.30435964657991582;
      Ds[272] = 0.60872434048853374;
      Ds[273] = -0.30435964657991588;
      Ds[274] = -0.30435964657991604;
      Ds[275] = -1.2174436336483656;
      Ds[276] = -0.30435964657991604;
      Ds[277] = 0.0;
      Ds[278] = 0.0;
      Ds[279] = 0.0;
      Ds[280] = 0.0;
      Ds[281] = 0.0;
      Ds[282] = 0.0;
      Ds[283] = -0.91308398706844984;
      Ds[284] = -0.91308398706844984;
      Ds[285] = 0.0;
      Ds[286] = 0.0;
      Ds[287] = 0.0;
      Ds[288] = 0.0;
      Ds[289] = 0.0;
      Ds[290] = 0.0;
      Ds[291] = -0.276389071193217;
      Ds[292] = 0.0552778142386434;
      Ds[293] = 0.276389071193217;
      Ds[294] = 0.0552778142386434;
      Ds[295] = -0.30435964657991582;
      Ds[296] = -0.30435964657991582;
      Ds[297] = 0.608724340488534;
      Ds[298] = -0.30435964657991593;
      Ds[299] = -0.304359646579916;
      Ds[300] = -1.2174436336483658;
      Ds[301] = 0.0;
      Ds[302] = 0.0;
      Ds[303] = 0.0;
      Ds[304] = 0.0;
      Ds[305] = 0.0;
      Ds[306] = 0.0;
      Ds[307] = -2.2116850665891483E-16;
      Ds[308] = -0.91308398706845;
      Ds[309] = 0.0;
      Ds[310] = 0.0;
      Ds[311] = 0.0;
      Ds[312] = -0.0059852692503629288;
      Ds[313] = 0.0029926346251814631;
      Ds[314] = 0.0029926346251814635;
      Ds[315] = 0.0;
      Ds[316] = 0.0;
      Ds[317] = 0.0;
      Ds[318] = 1.1111111110925926E-6;
      Ds[319] = 0.0;
      Ds[320] = 0.0;
      Ds[321] = 0.0;
      Ds[322] = 0.0;
      Ds[323] = 0.0;
      Ds[324] = 0.0;
      Ds[325] = -33341.231534202234;
      Ds[326] = 33341.231533091122;
      Ds[327] = -33329.384233732228;
      Ds[328] = 33329.384232621116;
      Ds[329] = -33329.384233732206;
      Ds[330] = 33329.384232621094;
      Ds[331] = 0.0;
      Ds[332] = 0.0;
      Ds[333] = -1.783078625882211E-17;
      Ds[334] = -1.783078625882212E-17;
      Ds[335] = 3.5661572517644252E-17;
      Ds[336] = 0.0059852692503629288;
      Ds[337] = -0.0029926346251814631;
      Ds[338] = -0.0029926346251814635;
      Ds[339] = 0.0;
      Ds[340] = 0.0;
      Ds[341] = 0.0;
      Ds[342] = 1.1111111110925926E-6;
      Ds[343] = 0.0;
      Ds[344] = 0.0;
      Ds[345] = 0.0;
      Ds[346] = 0.0;
      Ds[347] = 0.0;
      Ds[348] = 0.0;
      Ds[349] = 33341.231533091122;
      Ds[350] = -33341.231534202234;
      Ds[351] = 33329.384232621116;
      Ds[352] = -33329.384233732228;
      Ds[353] = 33329.384232621094;
      Ds[354] = -33329.384233732206;
      Ds[355] = 0.0;
      Ds[356] = 0.0;
      Ds[357] = 1.783078625882211E-17;
      Ds[358] = 1.783078625882212E-17;
      Ds[359] = -3.5661572517644252E-17;
      Ds[360] = 0.0029926346251814644;
      Ds[361] = -0.0059852692503629262;
      Ds[362] = 0.0029926346251814644;
      Ds[363] = 0.0;
      Ds[364] = 0.0;
      Ds[365] = 0.0;
      Ds[366] = 1.1111111110925926E-6;
      Ds[367] = 0.0;
      Ds[368] = 0.0;
      Ds[369] = 0.0;
      Ds[370] = 0.0;
      Ds[371] = 0.0;
      Ds[372] = 0.0;
      Ds[373] = -33329.38423373222;
      Ds[374] = 33329.384232621109;
      Ds[375] = -33341.231534202219;
      Ds[376] = 33341.231533091108;
      Ds[377] = -33329.38423373222;
      Ds[378] = 33329.384232621109;
      Ds[379] = 0.0;
      Ds[380] = 0.0;
      Ds[381] = 3.5661572517644233E-17;
      Ds[382] = -1.7830786258822135E-17;
      Ds[383] = -1.7830786258822126E-17;
      Ds[384] = -0.0029926346251814644;
      Ds[385] = 0.0059852692503629262;
      Ds[386] = -0.0029926346251814644;
      Ds[387] = 0.0;
      Ds[388] = 0.0;
      Ds[389] = 0.0;
      Ds[390] = 1.1111111110925926E-6;
      Ds[391] = 0.0;
      Ds[392] = 0.0;
      Ds[393] = 0.0;
      Ds[394] = 0.0;
      Ds[395] = 0.0;
      Ds[396] = 0.0;
      Ds[397] = 33329.384232621109;
      Ds[398] = -33329.38423373222;
      Ds[399] = 33341.231533091108;
      Ds[400] = -33341.231534202219;
      Ds[401] = 33329.384232621109;
      Ds[402] = -33329.38423373222;
      Ds[403] = 0.0;
      Ds[404] = 0.0;
      Ds[405] = -3.5661572517644233E-17;
      Ds[406] = 1.7830786258822135E-17;
      Ds[407] = 1.7830786258822126E-17;
      Ds[408] = 0.0029926346251814644;
      Ds[409] = 0.0029926346251814631;
      Ds[410] = -0.0059852692503629288;
      Ds[411] = 0.0;
      Ds[412] = 0.0;
      Ds[413] = 0.0;
      Ds[414] = 1.1111111110925926E-6;
      Ds[415] = 0.0;
      Ds[416] = 0.0;
      Ds[417] = 0.0;
      Ds[418] = 0.0;
      Ds[419] = 0.0;
      Ds[420] = 0.0;
      Ds[421] = -33329.384233732213;
      Ds[422] = 33329.3842326211;
      Ds[423] = -33329.384233732213;
      Ds[424] = 33329.3842326211;
      Ds[425] = -33341.231534202248;
      Ds[426] = 33341.231533091122;
      Ds[427] = 0.0;
      Ds[428] = 0.0;
      Ds[429] = -1.7830786258822123E-17;
      Ds[430] = 3.5661572517644252E-17;
      Ds[431] = -1.7830786258822123E-17;
      Ds[432] = -0.0029926346251814644;
      Ds[433] = -0.0029926346251814631;
      Ds[434] = 0.0059852692503629288;
      Ds[435] = 0.0;
      Ds[436] = 0.0;
      Ds[437] = 0.0;
      Ds[438] = 1.1111111110925926E-6;
      Ds[439] = 0.0;
      Ds[440] = 0.0;
      Ds[441] = 0.0;
      Ds[442] = 0.0;
      Ds[443] = 0.0;
      Ds[444] = 0.0;
      Ds[445] = 33329.3842326211;
      Ds[446] = -33329.384233732213;
      Ds[447] = 33329.3842326211;
      Ds[448] = -33329.384233732213;
      Ds[449] = 33341.231533091122;
      Ds[450] = -33341.231534202248;
      Ds[451] = 0.0;
      Ds[452] = 0.0;
      Ds[453] = 1.7830786258822123E-17;
      Ds[454] = -3.5661572517644252E-17;
      Ds[455] = 1.7830786258822123E-17;
      Ds[456] = 0.0;
      Ds[457] = 0.0;
      Ds[458] = 0.0;
      Ds[459] = -5.0900211957270554E-17;
      Ds[460] = 1.0180042391454112E-17;
      Ds[461] = 5.0900211957270554E-17;
      Ds[462] = 1.0180042391454112E-17;
      Ds[463] = -0.91308398706845006;
      Ds[464] = -6.691266441319681E-19;
      Ds[465] = 0.91308398706844984;
      Ds[466] = 0.91308398706845;
      Ds[467] = -1.6720258033790545E-16;
      Ds[468] = -0.91308398706845;
      Ds[469] = 0.0;
      Ds[470] = 0.0;
      Ds[471] = 0.0;
      Ds[472] = 0.0;
      Ds[473] = 0.0;
      Ds[474] = 0.0;
      Ds[475] = 0.91308398706844984;
      Ds[476] = -0.91308398706845;
      Ds[477] = 0.0;
      Ds[478] = 0.0;
      Ds[479] = 0.0;
      Ds[480] = 0.0;
      Ds[481] = 0.0;
      Ds[482] = 0.0;
      Ds[483] = 5.0596396600247315E-17;
      Ds[484] = -1.0119279320049464E-17;
      Ds[485] = -5.0596396600247315E-17;
      Ds[486] = -1.0119279320049464E-17;
      Ds[487] = -3.475577042637121E-19;
      Ds[488] = -0.91308398706844962;
      Ds[489] = -0.91308398706844984;
      Ds[490] = -3.475577042637121E-19;
      Ds[491] = 0.91308398706844984;
      Ds[492] = 0.91308398706845;
      Ds[493] = 0.0;
      Ds[494] = 0.0;
      Ds[495] = 0.0;
      Ds[496] = 0.0;
      Ds[497] = 0.0;
      Ds[498] = 0.0;
      Ds[499] = 0.91308398706844984;
      Ds[500] = 1.8261679741368997;
      Ds[501] = 0.0;
      Ds[502] = 0.0;
      Ds[503] = 0.0;
      Ds[504] = -3.3041936470981454;
      Ds[505] = 0.000808316421125539;
      Ds[506] = 0.00080831642112473762;
      Ds[507] = 0.0;
      Ds[508] = 0.0;
      Ds[509] = 0.0;
      Ds[510] = 0.0;
      Ds[511] = 0.0;
      Ds[512] = 0.0;
      Ds[513] = 0.0;
      Ds[514] = 0.0;
      Ds[515] = 0.0;
      Ds[516] = 0.0;
      Ds[517] = 1.9781982195883449E-7;
      Ds[518] = -1.9781982195883449E-7;
      Ds[519] = -9.89099109807355E-8;
      Ds[520] = 9.89099109807355E-8;
      Ds[521] = -9.8909910978098982E-8;
      Ds[522] = 9.8909910978098982E-8;
      Ds[523] = 0.0;
      Ds[524] = 0.0;
      Ds[525] = 1.986008894521048E-10;
      Ds[526] = 1.9860088945210493E-10;
      Ds[527] = -3.9720177890420986E-10;
      Ds[528] = 0.00080831642112560455;
      Ds[529] = -3.3041936470835935;
      Ds[530] = 0.00080831642112473762;
      Ds[531] = 0.0;
      Ds[532] = 0.0;
      Ds[533] = 0.0;
      Ds[534] = 0.0;
      Ds[535] = 0.0;
      Ds[536] = 0.0;
      Ds[537] = 0.0;
      Ds[538] = 0.0;
      Ds[539] = 0.0;
      Ds[540] = 0.0;
      Ds[541] = -9.8909910978492138E-8;
      Ds[542] = 9.8909910978492138E-8;
      Ds[543] = 1.978198219569842E-7;
      Ds[544] = -1.978198219569842E-7;
      Ds[545] = -9.8909910978492045E-8;
      Ds[546] = 9.8909910978492045E-8;
      Ds[547] = 0.0;
      Ds[548] = 0.0;
      Ds[549] = -3.9720177890420975E-10;
      Ds[550] = 1.9860088945210493E-10;
      Ds[551] = 1.9860088945210495E-10;
      Ds[552] = 0.00080831642112627123;
      Ds[553] = 0.00080831642112780386;
      Ds[554] = -3.3041936470981454;
      Ds[555] = 0.0;
      Ds[556] = 0.0;
      Ds[557] = 0.0;
      Ds[558] = 0.0;
      Ds[559] = 0.0;
      Ds[560] = 0.0;
      Ds[561] = 0.0;
      Ds[562] = 0.0;
      Ds[563] = 0.0;
      Ds[564] = 0.0;
      Ds[565] = -9.8909910978836457E-8;
      Ds[566] = 9.8909910978836457E-8;
      Ds[567] = -9.8909910978836589E-8;
      Ds[568] = 9.8909910978836589E-8;
      Ds[569] = 1.9781982195767305E-7;
      Ds[570] = -1.9781982195767305E-7;
      Ds[571] = 0.0;
      Ds[572] = 0.0;
      Ds[573] = 1.9860088945210488E-10;
      Ds[574] = -3.9720177890420986E-10;
      Ds[575] = 1.9860088945210488E-10;
      Ds[576] = -1.2856960750436961E-5;
      Ds[577] = 1.2856960750472785E-5;
      Ds[578] = 1.2056677248449906E-17;
      Ds[579] = 0.0;
      Ds[580] = 0.0;
      Ds[581] = 0.0;
      Ds[582] = 0.0;
      Ds[583] = 0.0;
      Ds[584] = 0.0;
      Ds[585] = 0.0;
      Ds[586] = 0.0;
      Ds[587] = 0.0;
      Ds[588] = 0.0;
      Ds[589] = 1.1543238315637869E-12;
      Ds[590] = -1.1543238315637869E-12;
      Ds[591] = -1.1543238315641375E-12;
      Ds[592] = 1.1543238315641375E-12;
      Ds[593] = 3.5074546017513777E-25;
      Ds[594] = -3.5074546017513777E-25;
      Ds[595] = 0.0;
      Ds[596] = 0.0;
      Ds[597] = 9.0723895848336192E-5;
      Ds[598] = 1.0592163324834809E-20;
      Ds[599] = -9.0723895848336219E-5;
      Ds[600] = -8.5630613685493308E-18;
      Ds[601] = -1.2856960750483741E-5;
      Ds[602] = 1.2856960750433776E-5;
      Ds[603] = 0.0;
      Ds[604] = 0.0;
      Ds[605] = 0.0;
      Ds[606] = 0.0;
      Ds[607] = 0.0;
      Ds[608] = 0.0;
      Ds[609] = 0.0;
      Ds[610] = 0.0;
      Ds[611] = 0.0;
      Ds[612] = 0.0;
      Ds[613] = -9.83001773903465E-25;
      Ds[614] = 9.83001773903465E-25;
      Ds[615] = 1.1543238315648021E-12;
      Ds[616] = -1.1543238315648021E-12;
      Ds[617] = -1.154323831563819E-12;
      Ds[618] = 1.154323831563819E-12;
      Ds[619] = 0.0;
      Ds[620] = 0.0;
      Ds[621] = -9.0723895848336233E-5;
      Ds[622] = 9.0723895848336287E-5;
      Ds[623] = 3.0000259725979364E-21;
      Ds[624] = 1.2856960750414743E-5;
      Ds[625] = 2.0150090987584866E-18;
      Ds[626] = -1.2856960750434902E-5;
      Ds[627] = 0.0;
      Ds[628] = 0.0;
      Ds[629] = 0.0;
      Ds[630] = 0.0;
      Ds[631] = 0.0;
      Ds[632] = 0.0;
      Ds[633] = 0.0;
      Ds[634] = 0.0;
      Ds[635] = 0.0;
      Ds[636] = 0.0;
      Ds[637] = -1.1543238315609018E-12;
      Ds[638] = 1.1543238315609018E-12;
      Ds[639] = -7.2358653190616468E-25;
      Ds[640] = 7.2358653190616468E-25;
      Ds[641] = 1.1543238315616254E-12;
      Ds[642] = -1.1543238315616254E-12;
      Ds[643] = 0.0;
      Ds[644] = 0.0;
      Ds[645] = -9.74392505670362E-21;
      Ds[646] = -9.0723895848336274E-5;
      Ds[647] = 9.0723895848336287E-5;
      Ds[648] = 0.0;
      Ds[649] = 0.0;
      Ds[650] = 0.0;
      Ds[651] = -5.0900211957270554E-17;
      Ds[652] = 1.0180042391454112E-17;
      Ds[653] = 5.0900211957270554E-17;
      Ds[654] = 1.0180042391454112E-17;
      Ds[655] = -0.91308398706845006;
      Ds[656] = -6.691266441319681E-19;
      Ds[657] = 0.91308398706844984;
      Ds[658] = 0.91308398706845;
      Ds[659] = -1.6720258033790545E-16;
      Ds[660] = -0.91308398706845;
      Ds[661] = 0.0;
      Ds[662] = 0.0;
      Ds[663] = 0.0;
      Ds[664] = 0.0;
      Ds[665] = 0.0;
      Ds[666] = 0.0;
      Ds[667] = 0.91308398706844984;
      Ds[668] = -0.91308398706845;
      Ds[669] = 0.0;
      Ds[670] = 0.0;
      Ds[671] = 0.0;
      Ds[672] = 0.0;
      Ds[673] = 0.0;
      Ds[674] = 0.0;
      Ds[675] = 5.0596396600247315E-17;
      Ds[676] = -1.0119279320049464E-17;
      Ds[677] = -5.0596396600247315E-17;
      Ds[678] = -1.0119279320049464E-17;
      Ds[679] = -3.475577042637121E-19;
      Ds[680] = -0.91308398706844962;
      Ds[681] = -0.91308398706844984;
      Ds[682] = -3.475577042637121E-19;
      Ds[683] = 0.91308398706844984;
      Ds[684] = 0.91308398706845;
      Ds[685] = 0.0;
      Ds[686] = 0.0;
      Ds[687] = 0.0;
      Ds[688] = 0.0;
      Ds[689] = 0.0;
      Ds[690] = 0.0;
      Ds[691] = 0.91308398706844984;
      Ds[692] = 1.8261679741368997;
      Ds[693] = 0.0;
      Ds[694] = 0.0;
      Ds[695] = 0.0;
      Ds[696] = 0.0;
      Ds[697] = 0.0;
      Ds[698] = 0.0;
      Ds[699] = -3.015340742665086E-19;
      Ds[700] = 6.030681485330172E-20;
      Ds[701] = 3.015340742665086E-19;
      Ds[702] = 6.030681485330172E-20;
      Ds[703] = 0.91308398706845006;
      Ds[704] = 0.91308398706844962;
      Ds[705] = 8.5136302390562374E-19;
      Ds[706] = -0.91308398706845006;
      Ds[707] = -0.91308398706844962;
      Ds[708] = 8.9257964648932882E-19;
      Ds[709] = 0.0;
      Ds[710] = 0.0;
      Ds[711] = 0.0;
      Ds[712] = 0.0;
      Ds[713] = 0.0;
      Ds[714] = 0.0;
      Ds[715] = -1.8261679741368997;
      Ds[716] = -0.91308398706844962;
      Ds[717] = 0.0;
      Ds[718] = 0.0;
      Ds[719] = 0.0;
      Ds[720] = 0.00897825853655043;
      Ds[721] = 3.5466100604065681E-7;
      Ds[722] = 3.5466100604043992E-7;
      Ds[723] = 0.0;
      Ds[724] = 0.0;
      Ds[725] = 0.0;
      Ds[726] = 0.0;
      Ds[727] = 0.0;
      Ds[728] = 0.0;
      Ds[729] = 0.0;
      Ds[730] = 0.0;
      Ds[731] = 0.0;
      Ds[732] = 0.0;
      Ds[733] = 7.8982003133423859;
      Ds[734] = -7.8982003133423859;
      Ds[735] = -3.9491001566493651;
      Ds[736] = 3.9491001566493651;
      Ds[737] = -3.9491001566857449;
      Ds[738] = 3.9491001566857449;
      Ds[739] = 0.0;
      Ds[740] = 0.0;
      Ds[741] = 1.7830786258822114E-17;
      Ds[742] = 1.7830786258822126E-17;
      Ds[743] = -3.5661572517644264E-17;
      Ds[744] = 3.5466100603969977E-7;
      Ds[745] = 0.00897825853655043;
      Ds[746] = 3.5466100604043992E-7;
      Ds[747] = 0.0;
      Ds[748] = 0.0;
      Ds[749] = 0.0;
      Ds[750] = 0.0;
      Ds[751] = 0.0;
      Ds[752] = 0.0;
      Ds[753] = 0.0;
      Ds[754] = 0.0;
      Ds[755] = 0.0;
      Ds[756] = 0.0;
      Ds[757] = -3.949100156663917;
      Ds[758] = 3.949100156663917;
      Ds[759] = 7.8982003133423859;
      Ds[760] = -7.8982003133423859;
      Ds[761] = -3.949100156671193;
      Ds[762] = 3.949100156671193;
      Ds[763] = 0.0;
      Ds[764] = 0.0;
      Ds[765] = -3.5661572517644239E-17;
      Ds[766] = 1.7830786258822135E-17;
      Ds[767] = 1.7830786258822126E-17;
      Ds[768] = 3.54661006040077E-7;
      Ds[769] = 3.5466100604065681E-7;
      Ds[770] = 0.00897825853655043;
      Ds[771] = 0.0;
      Ds[772] = 0.0;
      Ds[773] = 0.0;
      Ds[774] = 0.0;
      Ds[775] = 0.0;
      Ds[776] = 0.0;
      Ds[777] = 0.0;
      Ds[778] = 0.0;
      Ds[779] = 0.0;
      Ds[780] = 0.0;
      Ds[781] = -3.9491001566784689;
      Ds[782] = 3.9491001566784689;
      Ds[783] = -3.949100156671193;
      Ds[784] = 3.949100156671193;
      Ds[785] = 7.8982003133569378;
      Ds[786] = -7.8982003133569378;
      Ds[787] = 0.0;
      Ds[788] = 0.0;
      Ds[789] = 1.7830786258822126E-17;
      Ds[790] = -3.5661572517644264E-17;
      Ds[791] = 1.7830786258822123E-17;
      Ds[792] = 0.0089789671759676216;
      Ds[793] = 3.4129744437854866E-10;
      Ds[794] = 3.41297444378474E-10;
      Ds[795] = 0.0;
      Ds[796] = 0.0;
      Ds[797] = 0.0;
      Ds[798] = 0.0;
      Ds[799] = 0.0;
      Ds[800] = 0.0;
      Ds[801] = 0.0;
      Ds[802] = 0.0;
      Ds[803] = 0.0;
      Ds[804] = 0.0;
      Ds[805] = 0.0059854670701848852;
      Ds[806] = -0.0059854670701848852;
      Ds[807] = -0.0029927335350924448;
      Ds[808] = 0.0029927335350924448;
      Ds[809] = -0.0029927335350924409;
      Ds[810] = 0.0029927335350924409;
      Ds[811] = 0.0;
      Ds[812] = 0.0;
      Ds[813] = 1.7832897374872862E-17;
      Ds[814] = 1.7832897374872881E-17;
      Ds[815] = -3.5665794749745755E-17;
      Ds[816] = 3.4129744437854954E-10;
      Ds[817] = 0.0089789671759676233;
      Ds[818] = 3.4129744437847395E-10;
      Ds[819] = 0.0;
      Ds[820] = 0.0;
      Ds[821] = 0.0;
      Ds[822] = 0.0;
      Ds[823] = 0.0;
      Ds[824] = 0.0;
      Ds[825] = 0.0;
      Ds[826] = 0.0;
      Ds[827] = 0.0;
      Ds[828] = 0.0;
      Ds[829] = -0.0029927335350924435;
      Ds[830] = 0.0029927335350924435;
      Ds[831] = 0.0059854670701848869;
      Ds[832] = -0.0059854670701848869;
      Ds[833] = -0.0029927335350924435;
      Ds[834] = 0.0029927335350924435;
      Ds[835] = 0.0;
      Ds[836] = 0.0;
      Ds[837] = -3.5665794749745736E-17;
      Ds[838] = 1.7832897374872877E-17;
      Ds[839] = 1.7832897374872877E-17;
      Ds[840] = 3.412974443786137E-10;
      Ds[841] = 3.4129744437874749E-10;
      Ds[842] = 0.0089789671759676216;
      Ds[843] = 0.0;
      Ds[844] = 0.0;
      Ds[845] = 0.0;
      Ds[846] = 0.0;
      Ds[847] = 0.0;
      Ds[848] = 0.0;
      Ds[849] = 0.0;
      Ds[850] = 0.0;
      Ds[851] = 0.0;
      Ds[852] = 0.0;
      Ds[853] = -0.0029927335350924417;
      Ds[854] = 0.0029927335350924417;
      Ds[855] = -0.0029927335350924417;
      Ds[856] = 0.0029927335350924417;
      Ds[857] = 0.0059854670701848835;
      Ds[858] = -0.0059854670701848835;
      Ds[859] = 0.0;
      Ds[860] = 0.0;
      Ds[861] = 1.7832897374872874E-17;
      Ds[862] = -3.5665794749745749E-17;
      Ds[863] = 1.7832897374872874E-17;
      Ds[864] = 0.99996686827385728;
      Ds[865] = 8.0831607982809469E-9;
      Ds[866] = 8.0831607982729332E-9;
      Ds[867] = 0.0;
      Ds[868] = 0.0;
      Ds[869] = 0.0;
      Ds[870] = 0.0;
      Ds[871] = 0.0;
      Ds[872] = 0.0;
      Ds[873] = 0.0;
      Ds[874] = 0.0;
      Ds[875] = 0.0;
      Ds[876] = 0.0;
      Ds[877] = -5.9852692503629268E-8;
      Ds[878] = 5.9852692503629268E-8;
      Ds[879] = 2.9926346251814647E-8;
      Ds[880] = -2.9926346251814647E-8;
      Ds[881] = 2.9926346251814627E-8;
      Ds[882] = -2.9926346251814627E-8;
      Ds[883] = 0.0;
      Ds[884] = 0.0;
      Ds[885] = 1.9860087161920745E-15;
      Ds[886] = 1.986008716192076E-15;
      Ds[887] = -3.9720174323841513E-15;
      Ds[888] = 8.083160798281602E-9;
      Ds[889] = 0.9999668682738575;
      Ds[890] = 8.0831607982729332E-9;
      Ds[891] = 0.0;
      Ds[892] = 0.0;
      Ds[893] = 0.0;
      Ds[894] = 0.0;
      Ds[895] = 0.0;
      Ds[896] = 0.0;
      Ds[897] = 0.0;
      Ds[898] = 0.0;
      Ds[899] = 0.0;
      Ds[900] = 0.0;
      Ds[901] = 2.9926346251814654E-8;
      Ds[902] = -2.9926346251814654E-8;
      Ds[903] = -5.9852692503629307E-8;
      Ds[904] = 5.9852692503629307E-8;
      Ds[905] = 2.9926346251814654E-8;
      Ds[906] = -2.9926346251814654E-8;
      Ds[907] = 0.0;
      Ds[908] = 0.0;
      Ds[909] = -3.9720174323841505E-15;
      Ds[910] = 1.9860087161920757E-15;
      Ds[911] = 1.9860087161920757E-15;
      Ds[912] = 8.0831607982882691E-9;
      Ds[913] = 8.0831607983035951E-9;
      Ds[914] = 0.99996686827385728;
      Ds[915] = 0.0;
      Ds[916] = 0.0;
      Ds[917] = 0.0;
      Ds[918] = 0.0;
      Ds[919] = 0.0;
      Ds[920] = 0.0;
      Ds[921] = 0.0;
      Ds[922] = 0.0;
      Ds[923] = 0.0;
      Ds[924] = 0.0;
      Ds[925] = 2.9926346251814627E-8;
      Ds[926] = -2.9926346251814627E-8;
      Ds[927] = 2.9926346251814627E-8;
      Ds[928] = -2.9926346251814627E-8;
      Ds[929] = -5.9852692503629255E-8;
      Ds[930] = 5.9852692503629255E-8;
      Ds[931] = 0.0;
      Ds[932] = 0.0;
      Ds[933] = 1.9860087161920753E-15;
      Ds[934] = -3.9720174323841513E-15;
      Ds[935] = 1.9860087161920753E-15;
      Ds[936] = 3.3433129030807066E-9;
      Ds[937] = -3.3433129030920274E-9;
      Ds[938] = -7.45808256922626E-19;
      Ds[939] = 0.0;
      Ds[940] = 0.0;
      Ds[941] = 0.0;
      Ds[942] = 0.0;
      Ds[943] = 0.0;
      Ds[944] = 0.0;
      Ds[945] = 0.0;
      Ds[946] = 0.0;
      Ds[947] = 0.0;
      Ds[948] = 0.0;
      Ds[949] = -3.0016936627547367E-16;
      Ds[950] = 3.0016936627547367E-16;
      Ds[951] = 3.0016936623117247E-16;
      Ds[952] = -3.0016936623117247E-16;
      Ds[953] = 4.4301234274844369E-26;
      Ds[954] = -4.4301234274844369E-26;
      Ds[955] = 0.0;
      Ds[956] = 0.0;
      Ds[957] = 1.5423361484630858E-5;
      Ds[958] = 1.4444491929825759E-21;
      Ds[959] = -1.5423361484630868E-5;
      Ds[960] = 9.2749010604473E-20;
      Ds[961] = 3.3433129015426321E-9;
      Ds[962] = -3.3433129039406365E-9;
      Ds[963] = 0.0;
      Ds[964] = 0.0;
      Ds[965] = 0.0;
      Ds[966] = 0.0;
      Ds[967] = 0.0;
      Ds[968] = 0.0;
      Ds[969] = 0.0;
      Ds[970] = 0.0;
      Ds[971] = 0.0;
      Ds[972] = 0.0;
      Ds[973] = -7.7317304911634345E-26;
      Ds[974] = 7.7317304911634345E-26;
      Ds[975] = -3.0016936618371357E-16;
      Ds[976] = 3.0016936618371357E-16;
      Ds[977] = 3.0016936626103086E-16;
      Ds[978] = -3.0016936626103086E-16;
      Ds[979] = 0.0;
      Ds[980] = 0.0;
      Ds[981] = -1.5423361484630865E-5;
      Ds[982] = 1.5423361484630872E-5;
      Ds[983] = -7.7046462251914867E-22;
      Ds[984] = -3.3433129015383109E-9;
      Ds[985] = -9.5661735926589552E-20;
      Ds[986] = 3.3433129014165185E-9;
      Ds[987] = 0.0;
      Ds[988] = 0.0;
      Ds[989] = 0.0;
      Ds[990] = 0.0;
      Ds[991] = 0.0;
      Ds[992] = 0.0;
      Ds[993] = 0.0;
      Ds[994] = 0.0;
      Ds[995] = 0.0;
      Ds[996] = 0.0;
      Ds[997] = 3.0016936610782761E-16;
      Ds[998] = -3.0016936610782761E-16;
      Ds[999] = 2.0807824258537781E-27;
      Ds[1000] = -2.0807824258537781E-27;
      Ds[1001] = -3.0016936610990843E-16;
      Ds[1002] = 3.0016936610990843E-16;
      Ds[1003] = 0.0;
      Ds[1004] = 0.0;
      Ds[1005] = -1.2334836783517453E-22;
      Ds[1006] = -1.5423361484630875E-5;
      Ds[1007] = 1.5423361484630875E-5;
      Ds[1008] = 0.0;
      Ds[1009] = 0.0;
      Ds[1010] = 0.0;
      Ds[1011] = 0.0;
      Ds[1012] = 0.0;
      Ds[1013] = 0.0;
      Ds[1014] = 0.0;
      Ds[1015] = 0.0;
      Ds[1016] = 0.0;
      Ds[1017] = 0.0;
      Ds[1018] = 0.0;
      Ds[1019] = 0.0;
      Ds[1020] = 0.0;
      Ds[1021] = 0.0;
      Ds[1022] = 0.0;
      Ds[1023] = 0.0;
      Ds[1024] = 0.0;
      Ds[1025] = 0.0;
      Ds[1026] = 0.0;
      Ds[1027] = -1.0;
      Ds[1028] = 0.0;
      Ds[1029] = 0.0;
      Ds[1030] = 0.0;
      Ds[1031] = 0.0;
      Ds[1032] = 0.0;
      Ds[1033] = 0.0;
      Ds[1034] = 0.0;
      Ds[1035] = 0.0;
      Ds[1036] = 0.0;
      Ds[1037] = 0.0;
      Ds[1038] = 0.0;
      Ds[1039] = 0.0;
      Ds[1040] = 0.0;
      Ds[1041] = 0.0;
      Ds[1042] = 0.0;
      Ds[1043] = 0.0;
      Ds[1044] = 0.0;
      Ds[1045] = 0.0;
      Ds[1046] = 0.0;
      Ds[1047] = 0.0;
      Ds[1048] = 0.0;
      Ds[1049] = 0.0;
      Ds[1050] = 0.0;
      Ds[1051] = 0.0;
      Ds[1052] = -1.0;
      Ds[1053] = 0.0;
      Ds[1054] = 0.0;
      Ds[1055] = 0.0;
      Ds[1056] = 0.0;
      Ds[1057] = 0.0;
      Ds[1058] = 0.0;
      Ds[1059] = 0.0;
      Ds[1060] = 0.0;
      Ds[1061] = 0.0;
      Ds[1062] = 0.0;
      Ds[1063] = 0.0;
      Ds[1064] = 0.0;
      Ds[1065] = 0.0;
      Ds[1066] = 0.0;
      Ds[1067] = 0.0;
      Ds[1068] = 0.0;
      Ds[1069] = 0.0;
      Ds[1070] = 0.0;
      Ds[1071] = 0.0;
      Ds[1072] = 0.0;
      Ds[1073] = 0.0;
      Ds[1074] = 0.0;
      Ds[1075] = 1.0;
      Ds[1076] = 1.0;
      Ds[1077] = 0.0;
      Ds[1078] = 0.0;
      Ds[1079] = 0.0;
      Ds[1080] = 5.9852692503629281E-8;
      Ds[1081] = -2.9926346251814634E-8;
      Ds[1082] = -2.992634625181464E-8;
      Ds[1083] = 0.0;
      Ds[1084] = 0.0;
      Ds[1085] = 0.0;
      Ds[1086] = 0.0;
      Ds[1087] = 0.0;
      Ds[1088] = 0.0;
      Ds[1089] = 0.0;
      Ds[1090] = 0.0;
      Ds[1091] = 0.0;
      Ds[1092] = 0.0;
      Ds[1093] = -0.66658768466353324;
      Ds[1094] = 0.66658768466353324;
      Ds[1095] = 0.33329384233176673;
      Ds[1096] = -0.33329384233176673;
      Ds[1097] = 0.33329384233176657;
      Ds[1098] = -0.33329384233176657;
      Ds[1099] = 0.0;
      Ds[1100] = 0.0;
      Ds[1101] = 1.7830786258822112E-22;
      Ds[1102] = 1.7830786258822122E-22;
      Ds[1103] = -3.5661572517644253E-22;
      Ds[1104] = -2.992634625181464E-8;
      Ds[1105] = 5.9852692503629268E-8;
      Ds[1106] = -2.992634625181464E-8;
      Ds[1107] = 0.0;
      Ds[1108] = 0.0;
      Ds[1109] = 0.0;
      Ds[1110] = 0.0;
      Ds[1111] = 0.0;
      Ds[1112] = 0.0;
      Ds[1113] = 0.0;
      Ds[1114] = 0.0;
      Ds[1115] = 0.0;
      Ds[1116] = 0.0;
      Ds[1117] = 0.33329384233176668;
      Ds[1118] = -0.33329384233176668;
      Ds[1119] = -0.66658768466353335;
      Ds[1120] = 0.66658768466353335;
      Ds[1121] = 0.33329384233176668;
      Ds[1122] = -0.33329384233176668;
      Ds[1123] = 0.0;
      Ds[1124] = 0.0;
      Ds[1125] = -3.5661572517644238E-22;
      Ds[1126] = 1.7830786258822133E-22;
      Ds[1127] = 1.7830786258822129E-22;
      Ds[1128] = -2.992634625181464E-8;
      Ds[1129] = -2.9926346251814634E-8;
      Ds[1130] = 5.9852692503629281E-8;
      Ds[1131] = 0.0;
      Ds[1132] = 0.0;
      Ds[1133] = 0.0;
      Ds[1134] = 0.0;
      Ds[1135] = 0.0;
      Ds[1136] = 0.0;
      Ds[1137] = 0.0;
      Ds[1138] = 0.0;
      Ds[1139] = 0.0;
      Ds[1140] = 0.0;
      Ds[1141] = 0.33329384233176662;
      Ds[1142] = -0.33329384233176662;
      Ds[1143] = 0.33329384233176662;
      Ds[1144] = -0.33329384233176662;
      Ds[1145] = -0.66658768466353324;
      Ds[1146] = 0.66658768466353324;
      Ds[1147] = 0.0;
      Ds[1148] = 0.0;
      Ds[1149] = 1.7830786258822126E-22;
      Ds[1150] = -3.5661572517644257E-22;
      Ds[1151] = 1.7830786258822124E-22;
      Ds[1152] = 0.0;
      Ds[1153] = 0.0;
      Ds[1154] = 0.0;
      Ds[1155] = 0.0;
      Ds[1156] = 0.0;
      Ds[1157] = 0.0;
      Ds[1158] = 2.2222222221851851E-6;
      Ds[1159] = 0.0;
      Ds[1160] = 0.0;
      Ds[1161] = 0.0;
      Ds[1162] = 0.0;
      Ds[1163] = 0.0;
      Ds[1164] = 0.0;
      Ds[1165] = -1.1111111110925924E-6;
      Ds[1166] = -1.1111111110925926E-6;
      Ds[1167] = -1.1111111110925924E-6;
      Ds[1168] = -1.1111111110925926E-6;
      Ds[1169] = -1.1111111110925924E-6;
      Ds[1170] = -1.1111111110925926E-6;
      Ds[1171] = 0.0;
      Ds[1172] = 0.0;
      Ds[1173] = 0.0;
      Ds[1174] = 0.0;
      Ds[1175] = 0.0;
      Ds[1176] = 0.0;
      Ds[1177] = 0.0;
      Ds[1178] = 0.0;
      Ds[1179] = -0.45404761518257353;
      Ds[1180] = 0.090809523036514711;
      Ds[1181] = 0.45404761518257353;
      Ds[1182] = 0.090809523036514711;
      Ds[1183] = 2.763890711932168E-6;
      Ds[1184] = 2.7638907119321684E-6;
      Ds[1185] = 2.7638907119321684E-6;
      Ds[1186] = 2.7638907119321689E-6;
      Ds[1187] = 2.7638907119321684E-6;
      Ds[1188] = 2.7638907119321684E-6;
      Ds[1189] = 0.0;
      Ds[1190] = 0.0;
      Ds[1191] = 0.0;
      Ds[1192] = 0.0;
      Ds[1193] = 0.0;
      Ds[1194] = 0.0;
      Ds[1195] = 4.9946609044459234E-22;
      Ds[1196] = -5.1587778558400329E-24;
      Ds[1197] = 0.0;
      Ds[1198] = 0.0;
      Ds[1199] = 0.0;

      {
        int_T i1;
        for (i1=0; i1 < 50; i1++) {
          Chopper[i1] = 1;
        }
      }

      {
        /* Switches work vectors */
        int_T *switch_status = (int_T*)
          Wind_songweiwei_DW.StateSpace_PWORK.SWITCH_STATUS;
        int_T *gState = (int_T*)Wind_songweiwei_DW.StateSpace_PWORK.G_STATE;
        real_T *yswitch = (real_T*)Wind_songweiwei_DW.StateSpace_PWORK.Y_SWITCH;
        int_T *switchTypes = (int_T*)
          Wind_songweiwei_DW.StateSpace_PWORK.SWITCH_TYPES;
        int_T *idxOutSw = (int_T*)Wind_songweiwei_DW.StateSpace_PWORK.IDX_OUT_SW;
        int_T *switch_status_init = (int_T*)
          Wind_songweiwei_DW.StateSpace_PWORK.SWITCH_STATUS_INIT;

        /* Initialize work vectors */
        switch_status[0] = 0;
        switch_status_init[0] = 0;
        gState[0] = (int_T) 0.0;
        yswitch[0] = 1/0.001;
        switchTypes[0] = (int_T)1.0;
        idxOutSw[0] = ((int_T)0.0) - 1;
        switch_status[1] = 0;
        switch_status_init[1] = 0;
        gState[1] = (int_T) 0.0;
        yswitch[1] = 1/0.001;
        switchTypes[1] = (int_T)1.0;
        idxOutSw[1] = ((int_T)0.0) - 1;
        switch_status[2] = 0;
        switch_status_init[2] = 0;
        gState[2] = (int_T) 0.0;
        yswitch[2] = 1/0.001;
        switchTypes[2] = (int_T)1.0;
        idxOutSw[2] = ((int_T)0.0) - 1;
        switch_status[3] = 0;
        switch_status_init[3] = 0;
        gState[3] = (int_T) 0.0;
        yswitch[3] = 1/0.001;
        switchTypes[3] = (int_T)1.0;
        idxOutSw[3] = ((int_T)0.0) - 1;
        switch_status[4] = 0;
        switch_status_init[4] = 0;
        gState[4] = (int_T) 0.0;
        yswitch[4] = 1/2.0E-7;
        switchTypes[4] = (int_T)1.0;
        idxOutSw[4] = ((int_T)0.0) - 1;
        switch_status[5] = 0;
        switch_status_init[5] = 0;
        gState[5] = (int_T) 0.0;
        yswitch[5] = 1/0.001;
        switchTypes[5] = (int_T)1.0;
        idxOutSw[5] = ((int_T)0.0) - 1;
        switch_status[6] = 0;
        switch_status_init[6] = 0;
        gState[6] = (int_T) 0.0;
        yswitch[6] = 1/2.0E-7;
        switchTypes[6] = (int_T)3.0;
        idxOutSw[6] = ((int_T)0.0) - 1;
        switch_status[7] = 0;
        switch_status_init[7] = 0;
        gState[7] = (int_T) 0.0;
        yswitch[7] = 1/2.0E-7;
        switchTypes[7] = (int_T)3.0;
        idxOutSw[7] = ((int_T)0.0) - 1;
        switch_status[8] = 0;
        switch_status_init[8] = 0;
        gState[8] = (int_T) 0.0;
        yswitch[8] = 1/2.0E-7;
        switchTypes[8] = (int_T)3.0;
        idxOutSw[8] = ((int_T)0.0) - 1;
        switch_status[9] = 0;
        switch_status_init[9] = 0;
        gState[9] = (int_T) 0.0;
        yswitch[9] = 1/2.0E-7;
        switchTypes[9] = (int_T)3.0;
        idxOutSw[9] = ((int_T)0.0) - 1;
        switch_status[10] = 0;
        switch_status_init[10] = 0;
        gState[10] = (int_T) 0.0;
        yswitch[10] = 1/2.0E-7;
        switchTypes[10] = (int_T)3.0;
        idxOutSw[10] = ((int_T)0.0) - 1;
        switch_status[11] = 0;
        switch_status_init[11] = 0;
        gState[11] = (int_T) 0.0;
        yswitch[11] = 1/2.0E-7;
        switchTypes[11] = (int_T)3.0;
        idxOutSw[11] = ((int_T)0.0) - 1;
        switch_status[12] = 0;
        switch_status_init[12] = 0;
        gState[12] = (int_T) 0.0;
        yswitch[12] = 1/2.0E-7;
        switchTypes[12] = (int_T)3.0;
        idxOutSw[12] = ((int_T)0.0) - 1;
        switch_status[13] = 0;
        switch_status_init[13] = 0;
        gState[13] = (int_T) 0.0;
        yswitch[13] = 1/2.0E-7;
        switchTypes[13] = (int_T)7.0;
        idxOutSw[13] = ((int_T)0.0) - 1;
        switch_status[14] = 0;
        switch_status_init[14] = 0;
        gState[14] = (int_T) 0.0;
        yswitch[14] = 1/2.0E-7;
        switchTypes[14] = (int_T)7.0;
        idxOutSw[14] = ((int_T)0.0) - 1;
        switch_status[15] = 0;
        switch_status_init[15] = 0;
        gState[15] = (int_T) 0.0;
        yswitch[15] = 1/2.0E-7;
        switchTypes[15] = (int_T)7.0;
        idxOutSw[15] = ((int_T)0.0) - 1;
        switch_status[16] = 0;
        switch_status_init[16] = 0;
        gState[16] = (int_T) 0.0;
        yswitch[16] = 1/2.0E-7;
        switchTypes[16] = (int_T)7.0;
        idxOutSw[16] = ((int_T)0.0) - 1;
        switch_status[17] = 0;
        switch_status_init[17] = 0;
        gState[17] = (int_T) 0.0;
        yswitch[17] = 1/2.0E-7;
        switchTypes[17] = (int_T)7.0;
        idxOutSw[17] = ((int_T)0.0) - 1;
        switch_status[18] = 0;
        switch_status_init[18] = 0;
        gState[18] = (int_T) 0.0;
        yswitch[18] = 1/2.0E-7;
        switchTypes[18] = (int_T)7.0;
        idxOutSw[18] = ((int_T)0.0) - 1;
      }
    }

    /* InitializeConditions for UnitDelay: '<S67>/Unit Delay3' */
    Wind_songweiwei_DW.UnitDelay3_DSTATE =
      Wind_songweiwei_P.UnitDelay3_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S67>/Unit Delay7' */
    Wind_songweiwei_DW.UnitDelay7_DSTATE =
      Wind_songweiwei_P.UnitDelay7_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S351>/dw_delay' */
    Wind_songweiwei_DW.dw_delay_DSTATE =
      Wind_songweiwei_P.dw_delay_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S351>/dw_predict' */
    Wind_songweiwei_DW.dw_predict_DSTATE =
      Wind_songweiwei_P.dw_predict_InitialCondition;

    /* InitializeConditions for DiscreteIntegrator: '<S126>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE =
      Wind_songweiwei_P.DiscretePIController_Init;

    /* InitializeConditions for UnitDelay: '<S221>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE =
      Wind_songweiwei_P.Delay_x1_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S221>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE =
      Wind_songweiwei_P.Delay_x2_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S197>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_c[0] =
      Wind_songweiwei_P.Delay_x1_InitialCondition_d;
    Wind_songweiwei_DW.Delay_x1_DSTATE_c[1] =
      Wind_songweiwei_P.Delay_x1_InitialCondition_d;
    Wind_songweiwei_DW.Delay_x1_DSTATE_c[2] =
      Wind_songweiwei_P.Delay_x1_InitialCondition_d;

    /* InitializeConditions for UnitDelay: '<S197>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_e[0] =
      Wind_songweiwei_P.Delay_x2_InitialCondition_d;
    Wind_songweiwei_DW.Delay_x2_DSTATE_e[1] =
      Wind_songweiwei_P.Delay_x2_InitialCondition_d;
    Wind_songweiwei_DW.Delay_x2_DSTATE_e[2] =
      Wind_songweiwei_P.Delay_x2_InitialCondition_d;

    /* InitializeConditions for DiscreteIntegrator: '<S329>/theta' */
    Wind_songweiwei_DW.theta_DSTATE = Wind_songweiwei_P.theta_IC;

    /* InitializeConditions for UnitDelay: '<S213>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_a =
      Wind_songweiwei_P.Delay_x1_InitialCondition_h;

    /* InitializeConditions for UnitDelay: '<S213>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_p =
      Wind_songweiwei_P.Delay_x2_InitialCondition_n;

    /* InitializeConditions for UnitDelay: '<S201>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_i[0] =
      Wind_songweiwei_P.Delay_x1_InitialCondition_e;
    Wind_songweiwei_DW.Delay_x1_DSTATE_i[1] =
      Wind_songweiwei_P.Delay_x1_InitialCondition_e;
    Wind_songweiwei_DW.Delay_x1_DSTATE_i[2] =
      Wind_songweiwei_P.Delay_x1_InitialCondition_e;

    /* InitializeConditions for UnitDelay: '<S201>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_m[0] =
      Wind_songweiwei_P.Delay_x2_InitialCondition_m;
    Wind_songweiwei_DW.Delay_x2_DSTATE_m[1] =
      Wind_songweiwei_P.Delay_x2_InitialCondition_m;
    Wind_songweiwei_DW.Delay_x2_DSTATE_m[2] =
      Wind_songweiwei_P.Delay_x2_InitialCondition_m;

    /* InitializeConditions for UnitDelay: '<S261>/Delay_x' */
    Wind_songweiwei_DW.Delay_x_DSTATE =
      Wind_songweiwei_P.Delay_x_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S217>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_ck =
      Wind_songweiwei_P.Delay_x1_InitialCondition_f;

    /* InitializeConditions for UnitDelay: '<S217>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_ee =
      Wind_songweiwei_P.Delay_x2_InitialCondition_e;

    /* InitializeConditions for UnitDelay: '<S189>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_m[0] =
      Wind_songweiwei_P.Delay_x1_InitialCondition_l;
    Wind_songweiwei_DW.Delay_x1_DSTATE_m[1] =
      Wind_songweiwei_P.Delay_x1_InitialCondition_l;
    Wind_songweiwei_DW.Delay_x1_DSTATE_m[2] =
      Wind_songweiwei_P.Delay_x1_InitialCondition_l;

    /* InitializeConditions for UnitDelay: '<S189>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_o[0] =
      Wind_songweiwei_P.Delay_x2_InitialCondition_di;
    Wind_songweiwei_DW.Delay_x2_DSTATE_o[1] =
      Wind_songweiwei_P.Delay_x2_InitialCondition_di;
    Wind_songweiwei_DW.Delay_x2_DSTATE_o[2] =
      Wind_songweiwei_P.Delay_x2_InitialCondition_di;

    /* InitializeConditions for DiscreteIntegrator: '<S225>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_h =
      Wind_songweiwei_P.DiscreteTimeIntegrator_IC_n;

    /* InitializeConditions for UnitDelay: '<S225>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_p =
      Wind_songweiwei_P.UnitDelay_InitialCondition_d;

    /* InitializeConditions for UnitDelay: '<S131>/IC = i_ic' */
    Wind_songweiwei_DW.ICi_ic_DSTATE = Wind_songweiwei_P.ICi_ic_InitialCondition;

    /* InitializeConditions for DiscreteIntegrator: '<S128>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_g =
      Wind_songweiwei_P.DiscretePIController_Init_n;

    /* InitializeConditions for UnitDelay: '<S134>/Delay_x' */
    Wind_songweiwei_DW.Delay_x_DSTATE_m =
      Wind_songweiwei_P.Delay_x_InitialCondition_c;

    /* InitializeConditions for UnitDelay: '<S121>/Unit Delay1' */
    Wind_songweiwei_DW.UnitDelay1_DSTATE =
      Wind_songweiwei_P.UnitDelay1_InitialCondition_a;

    /* InitializeConditions for DiscreteIntegrator: '<S121>/Discrete-Time Integrator1' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator1_IC_LOAD = 1U;

    /* InitializeConditions for DiscreteIntegrator: '<S154>/Integ4' */
    Wind_songweiwei_DW.Integ4_DSTATE = Wind_songweiwei_P.Integ4_IC_c;

    /* InitializeConditions for UnitDelay: '<S158>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_f =
      Wind_songweiwei_P.UnitDelay_InitialCondition_g;

    /* InitializeConditions for UnitDelay: '<S154>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_m =
      Wind_songweiwei_P.DiscreteVariableFrequencyMeanva;

    /* InitializeConditions for DiscreteIntegrator: '<S153>/Integ4' */
    Wind_songweiwei_DW.Integ4_DSTATE_n = Wind_songweiwei_P.Integ4_IC_m;

    /* InitializeConditions for UnitDelay: '<S156>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_fd =
      Wind_songweiwei_P.UnitDelay_InitialCondition_k;

    /* InitializeConditions for UnitDelay: '<S153>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_l =
      Wind_songweiwei_P.DiscreteVariableFrequencyMean_m;

    /* InitializeConditions for UnitDelay: '<S131>/IC = 0' */
    Wind_songweiwei_DW.IC0_DSTATE = Wind_songweiwei_P.IC0_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S193>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_b[0] =
      Wind_songweiwei_P.Delay_x1_InitialCondition_i;
    Wind_songweiwei_DW.Delay_x1_DSTATE_b[1] =
      Wind_songweiwei_P.Delay_x1_InitialCondition_i;
    Wind_songweiwei_DW.Delay_x1_DSTATE_b[2] =
      Wind_songweiwei_P.Delay_x1_InitialCondition_i;

    /* InitializeConditions for UnitDelay: '<S193>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_pn[0] =
      Wind_songweiwei_P.Delay_x2_InitialCondition_nr;
    Wind_songweiwei_DW.Delay_x2_DSTATE_pn[1] =
      Wind_songweiwei_P.Delay_x2_InitialCondition_nr;
    Wind_songweiwei_DW.Delay_x2_DSTATE_pn[2] =
      Wind_songweiwei_P.Delay_x2_InitialCondition_nr;

    /* InitializeConditions for DiscreteIntegrator: '<S129>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[0] =
      Wind_songweiwei_P.DiscretePIController1_Init;
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_n[1] =
      Wind_songweiwei_P.DiscretePIController1_Init;

    /* InitializeConditions for UnitDelay: '<S209>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_f =
      Wind_songweiwei_P.Delay_x1_InitialCondition_p;

    /* InitializeConditions for UnitDelay: '<S209>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_d =
      Wind_songweiwei_P.Delay_x2_InitialCondition_g;

    /* InitializeConditions for DiscreteIntegrator: '<S168>/Integ4' */
    Wind_songweiwei_DW.Integ4_DSTATE_nw = Wind_songweiwei_P.Integ4_IC_n;

    /* InitializeConditions for UnitDelay: '<S172>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_fz =
      Wind_songweiwei_P.UnitDelay_InitialCondition_f;

    /* InitializeConditions for UnitDelay: '<S168>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_b =
      Wind_songweiwei_P.DiscreteVariableFrequencyMean_j;

    /* InitializeConditions for DiscreteIntegrator: '<S167>/Integ4' */
    Wind_songweiwei_DW.Integ4_DSTATE_l = Wind_songweiwei_P.Integ4_IC_f;

    /* InitializeConditions for UnitDelay: '<S170>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_e =
      Wind_songweiwei_P.UnitDelay_InitialCondition_h;

    /* InitializeConditions for UnitDelay: '<S167>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_d =
      Wind_songweiwei_P.DiscreteVariableFrequencyMean_l;

    /* InitializeConditions for UnitDelay: '<S205>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_in =
      Wind_songweiwei_P.Delay_x1_InitialCondition_j;

    /* InitializeConditions for UnitDelay: '<S205>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_a =
      Wind_songweiwei_P.Delay_x2_InitialCondition_b;

    /* InitializeConditions for DiscreteIntegrator: '<S251>/Integ4' */
    Wind_songweiwei_DW.Integ4_DSTATE_a = Wind_songweiwei_P.Integ4_IC_ka;

    /* Level2 S-Function Block: '<S253>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Wind_songweiwei_M->childSfunctions[2];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S252>/Unit Delay' */
    Wind_songweiwei_DW.UnitDelay_DSTATE_n =
      Wind_songweiwei_P.UnitDelay_InitialCondition_i;

    /* InitializeConditions for UnitDelay: '<S251>/Unit Delay1' */
    Wind_songweiwei_DW.UnitDelay1_DSTATE_p =
      Wind_songweiwei_P.UnitDelay1_InitialCondition_c;

    /* InitializeConditions for DiscreteTransferFcn: '<S227>/Discrete Derivative ' */
    Wind_songweiwei_DW.DiscreteDerivative_states =
      Wind_songweiwei_P.DiscreteDerivative_InitialState;

    /* InitializeConditions for DiscreteIntegrator: '<S227>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_k =
      Wind_songweiwei_P.Discrete_Init;

    /* InitializeConditions for RateLimiter: '<S225>/Rate Limiter' */
    Wind_songweiwei_DW.PrevY = Wind_songweiwei_P.RateLimiter_IC;

    /* InitializeConditions for UnitDelay: '<S247>/Delay_x1' */
    Wind_songweiwei_DW.Delay_x1_DSTATE_bs =
      Wind_songweiwei_P.Delay_x1_InitialCondition_m;

    /* InitializeConditions for UnitDelay: '<S247>/Delay_x2' */
    Wind_songweiwei_DW.Delay_x2_DSTATE_oc =
      Wind_songweiwei_P.Delay_x2_InitialCondition_o;

    /* InitializeConditions for UnitDelay: '<S67>/Unit Delay2' */
    Wind_songweiwei_DW.UnitDelay2_DSTATE =
      Wind_songweiwei_P.UnitDelay2_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S310>/Delay_x' */
    Wind_songweiwei_DW.Delay_x_DSTATE_o =
      Wind_songweiwei_P.Delay_x_InitialCondition_b;

    /* InitializeConditions for DiscreteIntegrator: '<S305>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_o =
      Wind_songweiwei_P.DiscretePIController1_Init_i;

    /* InitializeConditions for RateLimiter: '<S125>/Rate Limiter   1' */
    Wind_songweiwei_DW.LastMajorTime = (rtInf);

    /* InitializeConditions for DiscreteIntegrator: '<S304>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_kk =
      Wind_songweiwei_P.DiscretePIController_Init_a;

    /* InitializeConditions for DiscreteIntegrator: '<S306>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTATE_b =
      Wind_songweiwei_P.DiscretePIController2_Init;

    /* InitializeConditions for UnitDelay: '<S309>/Delay_x' */
    Wind_songweiwei_DW.Delay_x_DSTATE_i =
      Wind_songweiwei_P.Delay_x_InitialCondition_d;

    /* InitializeConditions for UnitDelay: '<S67>/Unit Delay1' */
    Wind_songweiwei_DW.UnitDelay1_DSTATE_n =
      Wind_songweiwei_P.UnitDelay1_InitialCondition_i;

    /* InitializeConditions for UnitDelay: '<S67>/Unit Delay4' */
    Wind_songweiwei_DW.UnitDelay4_DSTATE =
      Wind_songweiwei_P.UnitDelay4_InitialCondition;

    /* InitializeConditions for DiscreteIntegrator: '<S69>/Discrete-Time Integrator' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator_DSTAT_bh =
      Wind_songweiwei_P.DriveTrain_w_wt0;

    /* InitializeConditions for DiscreteIntegrator: '<S69>/Discrete-Time Integrator1' */
    Wind_songweiwei_DW.DiscreteTimeIntegrator1_DSTAT_l =
      Wind_songweiwei_P.DriveTrain_torque0 / Wind_songweiwei_P.DriveTrain_Ksh;
    for (i = 0; i < 5; i++) {
      /* InitializeConditions for UnitDelay: '<S342>/fluxes' */
      Wind_songweiwei_DW.fluxes_DSTATE[i] =
        Wind_songweiwei_P.fluxes_InitialCondition[i];

      /* InitializeConditions for UnitDelay: '<S342>/voltages' */
      Wind_songweiwei_DW.voltages_DSTATE[i] =
        Wind_songweiwei_P.voltages_InitialCondition;
    }

    /* InitializeConditions for DiscreteIntegrator: '<S329>/Rotor speed deviation (dw)' */
    Wind_songweiwei_DW.Rotorspeeddeviationdw_DSTATE =
      Wind_songweiwei_P.Rotorspeeddeviationdw_IC;

    /* Enable for DiscreteIntegrator: '<S154>/Integ4' */
    Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE = 1U;

    /* Enable for DiscreteIntegrator: '<S153>/Integ4' */
    Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_f = 1U;

    /* Enable for DiscreteIntegrator: '<S168>/Integ4' */
    Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_l = 1U;

    /* Enable for DiscreteIntegrator: '<S167>/Integ4' */
    Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_h = 1U;

    /* Enable for DiscreteIntegrator: '<S251>/Integ4' */
    Wind_songweiwei_DW.Integ4_SYSTEM_ENABLE_a = 1U;

    /* Enable for DiscreteIntegrator: '<S329>/Rotor speed deviation (dw)' */
    Wind_songweiwei_DW.Rotorspeeddeviationdw_SYSTEM_EN = 1U;
  }
}

/* Model terminate function */
void Wind_songweiwei_terminate(void)
{
  /* S-Function block: <S357>/State-Space */
  {
    /* Free memory */
    free(Wind_songweiwei_DW.StateSpace_PWORK.AS);
    free(Wind_songweiwei_DW.StateSpace_PWORK.BS);
    free(Wind_songweiwei_DW.StateSpace_PWORK.CS);
    free(Wind_songweiwei_DW.StateSpace_PWORK.DS);
    free(Wind_songweiwei_DW.StateSpace_PWORK.DX_COL);
    free(Wind_songweiwei_DW.StateSpace_PWORK.TMP2);
    free(Wind_songweiwei_DW.StateSpace_PWORK.BD_COL);
    free(Wind_songweiwei_DW.StateSpace_PWORK.TMP1);
    free(Wind_songweiwei_DW.StateSpace_PWORK.XTMP);

    /*
     * Circuit has switches*/
    free(Wind_songweiwei_DW.StateSpace_PWORK.CHOPPER);
    free(Wind_songweiwei_DW.StateSpace_PWORK.G_STATE);
    free(Wind_songweiwei_DW.StateSpace_PWORK.SWITCH_STATUS);
    free(Wind_songweiwei_DW.StateSpace_PWORK.SW_CHG);
    free(Wind_songweiwei_DW.StateSpace_PWORK.SWITCH_STATUS_INIT);
  }

  /* S-Function block: <S159>/S-Function  */
  {
    /* Nothing to do! */
  }

  /* S-Function block: <S157>/S-Function  */
  {
    /* Nothing to do! */
  }

  /* S-Function block: <S173>/S-Function  */
  {
    /* Nothing to do! */
  }

  /* S-Function block: <S171>/S-Function  */
  {
    /* Nothing to do! */
  }

  /* Terminate for Enabled SubSystem: '<S225>/Automatic Gain Control' */

  /* Level2 S-Function Block: '<S237>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = Wind_songweiwei_M->childSfunctions[0];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S240>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = Wind_songweiwei_M->childSfunctions[1];
    sfcnTerminate(rts);
  }

  /* End of Terminate for SubSystem: '<S225>/Automatic Gain Control' */

  /* Level2 S-Function Block: '<S253>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = Wind_songweiwei_M->childSfunctions[2];
    sfcnTerminate(rts);
  }
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
