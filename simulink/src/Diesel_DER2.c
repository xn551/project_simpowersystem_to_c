/*
 * File: Diesel_DER2.c
 *
 * Code generated for Simulink model 'Diesel_DER2'.
 *
 * Model version                  : 1.99
 * Simulink Coder version         : 8.7 (R2014b) 08-Sep-2014
 * C/C++ source code generated on : Thu Mar 30 10:58:06 2017
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: 32-bit Generic
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "Diesel_DER2.h"
#include "Diesel_DER2_private.h"
#include "sfcn_bridge.h"

/* Block signals (auto storage) */
B_Diesel_DER2_T Diesel_DER2_B;

/* Block states (auto storage) */
DW_Diesel_DER2_T Diesel_DER2_DW;

/* External outputs (root outports fed by signals with auto storage) */
ExtY_Diesel_DER2_T Diesel_DER2_Y;

/* Real-time model */
RT_MODEL_Diesel_DER2_T Diesel_DER2_M_;
RT_MODEL_Diesel_DER2_T *const Diesel_DER2_M = &Diesel_DER2_M_;

/*
 * Time delay interpolation routine
 *
 * The linear interpolation is performed using the formula:
 *
 *          (t2 - tMinusDelay)         (tMinusDelay - t1)
 * u(t)  =  ----------------- * u1  +  ------------------- * u2
 *              (t2 - t1)                  (t2 - t1)
 */
real_T rt_TDelayInterpolate(
  real_T tMinusDelay,                  /* tMinusDelay = currentSimTime - delay */
  real_T tStart,
  real_T *tBuf,
  real_T *uBuf,
  int_T bufSz,
  int_T *lastIdx,
  int_T oldestIdx,
  int_T newIdx,
  real_T initOutput,
  boolean_T discrete,
  boolean_T minorStepAndTAtLastMajorOutput)
{
  int_T i;
  real_T yout, t1, t2, u1, u2;

  /*
   * If there is only one data point in the buffer, this data point must be
   * the t= 0 and tMinusDelay > t0, it ask for something unknown. The best
   * guess if initial output as well
   */
  if ((newIdx == 0) && (oldestIdx ==0 ) && (tMinusDelay > tStart))
    return initOutput;

  /*
   * If tMinusDelay is less than zero, should output initial value
   */
  if (tMinusDelay <= tStart)
    return initOutput;

  /* For fixed buffer extrapolation:
   * if tMinusDelay is small than the time at oldestIdx, if discrete, output
   * tailptr value,  else use tailptr and tailptr+1 value to extrapolate
   * It is also for fixed buffer. Note: The same condition can happen for transport delay block where
   * use tStart and and t[tail] other than using t[tail] and t[tail+1].
   * See below
   */
  if ((tMinusDelay <= tBuf[oldestIdx] ) ) {
    if (discrete) {
      return(uBuf[oldestIdx]);
    } else {
      int_T tempIdx= oldestIdx + 1;
      if (oldestIdx == bufSz-1)
        tempIdx = 0;
      t1= tBuf[oldestIdx];
      t2= tBuf[tempIdx];
      u1= uBuf[oldestIdx];
      u2= uBuf[tempIdx];
      if (t2 == t1) {
        if (tMinusDelay >= t2) {
          yout = u2;
        } else {
          yout = u1;
        }
      } else {
        real_T f1 = (t2-tMinusDelay) / (t2-t1);
        real_T f2 = 1.0 - f1;

        /*
         * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
         */
        yout = f1*u1 + f2*u2;
      }

      return yout;
    }
  }

  /*
   * When block does not have direct feedthrough, we use the table of
   * values to extrapolate off the end of the table for delays that are less
   * than 0 (less then step size).  This is not completely accurate.  The
   * chain of events is as follows for a given time t.  Major output - look
   * in table.  Update - add entry to table.  Now, if we call the output at
   * time t again, there is a new entry in the table. For very small delays,
   * this means that we will have a different answer from the previous call
   * to the output fcn at the same time t.  The following code prevents this
   * from happening.
   */
  if (minorStepAndTAtLastMajorOutput) {
    /* pretend that the new entry has not been added to table */
    if (newIdx != 0) {
      if (*lastIdx == newIdx) {
        (*lastIdx)--;
      }

      newIdx--;
    } else {
      if (*lastIdx == newIdx) {
        *lastIdx = bufSz-1;
      }

      newIdx = bufSz - 1;
    }
  }

  i = *lastIdx;
  if (tBuf[i] < tMinusDelay) {
    /* Look forward starting at last index */
    while (tBuf[i] < tMinusDelay) {
      /* May occur if the delay is less than step-size - extrapolate */
      if (i == newIdx)
        break;
      i = ( i < (bufSz-1) ) ? (i+1) : 0;/* move through buffer */
    }
  } else {
    /*
     * Look backwards starting at last index which can happen when the
     * delay time increases.
     */
    while (tBuf[i] >= tMinusDelay) {
      /*
       * Due to the entry condition at top of function, we
       * should never hit the end.
       */
      i = (i > 0) ? i-1 : (bufSz-1);   /* move through buffer */
    }

    i = ( i < (bufSz-1) ) ? (i+1) : 0;
  }

  *lastIdx = i;
  if (discrete) {
    /*
     * tempEps = 128 * eps;
     * localEps = max(tempEps, tempEps*fabs(tBuf[i]))/2;
     */
    double tempEps = (DBL_EPSILON) * 128.0;
    double localEps = tempEps * fabs(tBuf[i]);
    if (tempEps > localEps) {
      localEps = tempEps;
    }

    localEps = localEps / 2.0;
    if (tMinusDelay >= (tBuf[i] - localEps)) {
      yout = uBuf[i];
    } else {
      if (i == 0) {
        yout = uBuf[bufSz-1];
      } else {
        yout = uBuf[i-1];
      }
    }
  } else {
    if (i == 0) {
      t1 = tBuf[bufSz-1];
      u1 = uBuf[bufSz-1];
    } else {
      t1 = tBuf[i-1];
      u1 = uBuf[i-1];
    }

    t2 = tBuf[i];
    u2 = uBuf[i];
    if (t2 == t1) {
      if (tMinusDelay >= t2) {
        yout = u2;
      } else {
        yout = u1;
      }
    } else {
      real_T f1 = (t2-tMinusDelay) / (t2-t1);
      real_T f2 = 1.0 - f1;

      /*
       * Use Lagrange's interpolation formula.  Exact outputs at t1, t2.
       */
      yout = f1*u1 + f2*u2;
    }
  }

  return(yout);
}

/*
 * Initial conditions for enable system:
 *    '<S66>/Lmq_sat'
 *    '<S131>/Lmq_sat'
 */
void Diesel_DER2_Lmq_sat_Init(DW_Lmq_sat_Diesel_DER2_T *localDW,
  P_Lmq_sat_Diesel_DER2_T *localP)
{
  /* InitializeConditions for UnitDelay: '<S69>/Lmq_sat' */
  localDW->Lmq_sat_DSTATE = localP->Lmq_sat_InitialCondition;
}

/*
 * Start for enable system:
 *    '<S66>/Lmq_sat'
 *    '<S131>/Lmq_sat'
 */
void Diesel_DER2_Lmq_sat_Start(DW_Lmq_sat_Diesel_DER2_T *localDW,
  P_Lmq_sat_Diesel_DER2_T *localP)
{
  /* InitializeConditions for Enabled SubSystem: '<S66>/Lmq_sat' */
  Diesel_DER2_Lmq_sat_Init(localDW, localP);

  /* End of InitializeConditions for SubSystem: '<S66>/Lmq_sat' */
}

/*
 * Output and update for enable system:
 *    '<S66>/Lmq_sat'
 *    '<S131>/Lmq_sat'
 */
void Diesel_DER2_Lmq_sat(boolean_T rtu_Enable, const real_T rtu_phi[5],
  B_Lmq_sat_Diesel_DER2_T *localB, DW_Lmq_sat_Diesel_DER2_T *localDW,
  P_Lmq_sat_Diesel_DER2_T *localP)
{
  real_T rtb_phimq_c;
  real_T rtb_Switch_gd;

  /* Outputs for Enabled SubSystem: '<S66>/Lmq_sat' incorporates:
   *  EnablePort: '<S69>/Enable'
   */
  if (rtu_Enable) {
    /* Math: '<S73>/Math Function1' incorporates:
     *  Constant: '<S73>/u2'
     *  Math: '<S73>/Math Function'
     *  Sum: '<S73>/Sum2'
     *  UnitDelay: '<S69>/Lmq_sat'
     *
     * About '<S73>/Math Function1':
     *  Operator: reciprocal
     *
     * About '<S73>/Math Function':
     *  Operator: reciprocal
     */
    localB->MathFunction1 = 1.0 / (((localP->u2_Value[0] + localP->u2_Value[1])
      + localP->u2_Value[2]) + 1.0 / localDW->Lmq_sat_DSTATE);

    /* Abs: '<S69>/Abs' incorporates:
     *  Gain: '<S74>/1//Ll_q'
     *  Product: '<S74>/Product3'
     *  Sum: '<S74>/sum1'
     */
    rtb_Switch_gd = fabs((localP->Ll_q_Gain[0] * rtu_phi[0] + localP->Ll_q_Gain
                          [1] * rtu_phi[4]) * localB->MathFunction1);

    /* Lookup: '<S69>/Lookup Table' */
    rtb_phimq_c = rt_Lookup(localP->LookupTable_XData, 2, rtb_Switch_gd,
      localP->LookupTable_YData);

    /* Switch: '<S69>/Switch' incorporates:
     *  Constant: '<S69>/Constant1'
     *  Product: '<S69>/Product'
     */
    if (rtb_phimq_c != 0.0) {
      rtb_Switch_gd /= rtb_phimq_c;
    } else {
      rtb_Switch_gd = localP->Constant1_Value;
    }

    /* End of Switch: '<S69>/Switch' */

    /* Gain: '<S69>/Lmq' */
    localB->Lmsatq = localP->Lmq_Gain * rtb_Switch_gd;

    /* Update for UnitDelay: '<S69>/Lmq_sat' */
    localDW->Lmq_sat_DSTATE = localB->Lmsatq;
  }

  /* End of Outputs for SubSystem: '<S66>/Lmq_sat' */
}

/*
 * Output and update for enable system:
 *    '<S86>/Neg. Seq. Computation'
 *    '<S86>/Pos. Seq. Computation'
 */
void Diesel_DER2_NegSeqComputation(real_T rtu_Enable, creal_T rtu_In, creal_T
  rtu_In_i, creal_T rtu_In_o, B_NegSeqComputation_Diesel_DE_T *localB,
  P_NegSeqComputation_Diesel_DE_T *localP)
{
  /* Outputs for Enabled SubSystem: '<S86>/Neg. Seq. Computation' incorporates:
   *  EnablePort: '<S90>/Enable'
   */
  if (rtu_Enable > 0.0) {
    /* Gain: '<S90>/Gain3' incorporates:
     *  Gain: '<S90>/Gain1'
     *  Sum: '<S90>/Sum'
     */
    localB->Gain3.re = (((localP->Gain1_Gain[0].re * rtu_In.re -
                          localP->Gain1_Gain[0].im * rtu_In.im) +
                         (localP->Gain1_Gain[1].re * rtu_In_i.re -
                          localP->Gain1_Gain[1].im * rtu_In_i.im)) +
                        (localP->Gain1_Gain[2].re * rtu_In_o.re -
                         localP->Gain1_Gain[2].im * rtu_In_o.im)) *
      localP->Gain3_Gain;
    localB->Gain3.im = (((localP->Gain1_Gain[0].re * rtu_In.im +
                          localP->Gain1_Gain[0].im * rtu_In.re) +
                         (localP->Gain1_Gain[1].re * rtu_In_i.im +
                          localP->Gain1_Gain[1].im * rtu_In_i.re)) +
                        (localP->Gain1_Gain[2].re * rtu_In_o.im +
                         localP->Gain1_Gain[2].im * rtu_In_o.re)) *
      localP->Gain3_Gain;
  }

  /* End of Outputs for SubSystem: '<S86>/Neg. Seq. Computation' */
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

/* Model step function */
void Diesel_DER2_step(void)
{
  /* local block i/o variables */
  real_T rtb_pu2Watts;
  creal_T rtb_MagnitudeAngletoComplex;
  creal_T rtb_MagnitudeAngletoComplex1;
  creal_T rtb_MagnitudeAngletoComplex2;
  real_T rtb_xk1[5];
  real_T rtb_xk1_p[5];
  real_T rtb_Fcn2;
  real_T rtb_Fcn3;
  real_T rtb_outputformatting[18];
  real_T rtb_dw_delay;
  real_T rtb_xk1_c;
  real_T rtb_yk_m;
  real_T rtb_xk1_h;
  real_T rtb_xk1_i;
  real_T rtb_UnitDelay4;
  real_T rtb_donotdeletethisgain;
  real_T rtb_Vfd;
  real_T rtb_Fcn2_j;
  real_T rtb_Fcn3_m;
  real_T rtb_IC[5];
  real_T rtb_web1;
  real_T rtb_Kp5;
  real_T rtb_Product_n;
  real_T rtb_Switch;
  real_T rtb_Product1_f;
  real_T rtb_Switch_k;
  real_T rtb_Product_o;
  real_T rtb_Switch_i;
  real_T rtb_Product1_a;
  real_T rtb_Switch_o;
  real_T rtb_Product_h;
  real_T rtb_Switch_h;
  real_T rtb_Product1_e;
  real_T rtb_Switch_of;
  real_T rtb_ComplextoMagnitudeAngle_o1[3];
  real_T rtb_Vfd_o;
  real_T rtb_IC_o[5];
  real_T rtb_uH;
  real_T rtb_webase;
  real_T rtb_Switch_a;
  real_T rtb_Sum_k;
  real_T rtb_Sum_o;
  real_T rtb_Sum_b;
  real_T rtb_ElementaryMath;
  real_T rtb_donotdeletethisgain_i;
  real_T rtb_donotdeletethisgain_k;
  real_T rtb_DataTypeConversion;
  real_T rtb_Switch3_e;

  {
    real_T rtb_ElementaryMath1;
    real_T rtb_changeIqIdcurrentsigns[5];
    real_T rtb_Fcn1;
    real_T rtb_yk_n;
    real_T rtb_ElementaryMath1_d;
    real_T rtb_changeIqIdcurrentsigns_b[5];
    real_T rtb_phimd;
    real_T rtb_Switch_g;
    real_T rtb_Phisat;
    real_T rtb_Sum_nq;
    real_T rtb_Switch1_h;
    real_T rtb_Sum2_m3[25];
    real_T rtb_inversion_g[25];
    int32_T i;
    real_T rtb_Sum2_f[25];
    boolean_T tmp;
    real_T tmp_0[5];
    real_T tmp_1[5];
    int32_T i_0;
    real_T tmp_2[5];
    int32_T i_1;
    real_T rtb_Mult1_idx_0;
    real_T rtb_Mult1_idx_1;
    real_T rtb_Product_idx_0;
    real_T rtb_Product_idx_1;

    /* DigitalClock: '<S119>/Digital Clock' */
    rtb_ElementaryMath = ((Diesel_DER2_M->Timing.clockTick1) * 5.0E-5);

    /* Sum: '<S119>/Sum3' incorporates:
     *  DiscreteIntegrator: '<S119>/Rotor angle dthetae'
     *  Gain: '<S119>/web2'
     */
    rtb_ElementaryMath1 = Diesel_DER2_P.web2_Gain * rtb_ElementaryMath +
      Diesel_DER2_DW.Rotorangledthetae_DSTATE;

    /* Trigonometry: '<S122>/Elementary Math' */
    rtb_ElementaryMath = sin(rtb_ElementaryMath1);

    /* Trigonometry: '<S122>/Elementary Math1' */
    rtb_ElementaryMath1 = cos(rtb_ElementaryMath1);

    /* UnitDelay: '<S132>/fluxes' */
    for (i = 0; i < 5; i++) {
      rtb_xk1[i] = Diesel_DER2_DW.fluxes_DSTATE[i];
    }

    /* End of UnitDelay: '<S132>/fluxes' */

    /* Outputs for Enabled SubSystem: '<S121>/Saturation' incorporates:
     *  EnablePort: '<S131>/Enable'
     */
    /* Constant: '<S121>/Constant1' */
    if (Diesel_DER2_P.Constant1_Value_d > 0.0) {
      /* Math: '<S136>/Math Function3' incorporates:
       *  Constant: '<S136>/u1'
       *  Math: '<S136>/Math Function2'
       *  Sum: '<S136>/Sum1'
       *  UnitDelay: '<S133>/Lmd_sat'
       *
       * About '<S136>/Math Function3':
       *  Operator: reciprocal
       *
       * About '<S136>/Math Function2':
       *  Operator: reciprocal
       */
      Diesel_DER2_B.MathFunction3 = 1.0 / (((Diesel_DER2_P.u1_Value_ek[0] +
        Diesel_DER2_P.u1_Value_ek[1]) + Diesel_DER2_P.u1_Value_ek[2]) + 1.0 /
        Diesel_DER2_DW.Lmd_sat_DSTATE);

      /* Abs: '<S133>/Abs' incorporates:
       *  Gain: '<S137>/1//Ll_d'
       *  Product: '<S137>/Product4'
       *  Sum: '<S137>/sum2'
       */
      rtb_Phisat = fabs(((Diesel_DER2_P.Ll_d_Gain_o[0] * rtb_xk1[1] +
                          Diesel_DER2_P.Ll_d_Gain_o[1] * rtb_xk1[2]) +
                         Diesel_DER2_P.Ll_d_Gain_o[2] * rtb_xk1[3]) *
                        Diesel_DER2_B.MathFunction3);

      /* Lookup: '<S133>/Lookup Table' */
      rtb_Switch_g = rt_Lookup(Diesel_DER2_P.LookupTable_XData_f, 2, rtb_Phisat,
        Diesel_DER2_P.LookupTable_YData_i);

      /* Switch: '<S133>/Switch' incorporates:
       *  Constant: '<S133>/Constant1'
       *  Product: '<S133>/Product'
       */
      if (rtb_Switch_g != 0.0) {
        rtb_Switch_g = rtb_Phisat / rtb_Switch_g;
      } else {
        rtb_Switch_g = Diesel_DER2_P.Constant1_Value_e;
      }

      /* End of Switch: '<S133>/Switch' */

      /* Gain: '<S133>/Lmd' */
      Diesel_DER2_B.Lmsatd = Diesel_DER2_P.Lmd_Gain_j * rtb_Switch_g;

      /* Outputs for Enabled SubSystem: '<S131>/Lmq_sat' */

      /* Constant: '<S131>/Constant1' */
      Diesel_DER2_Lmq_sat(Diesel_DER2_P.Constant1_Value_f, rtb_xk1,
                          &Diesel_DER2_B.Lmq_sat_l, &Diesel_DER2_DW.Lmq_sat_l,
                          (P_Lmq_sat_Diesel_DER2_T *)&Diesel_DER2_P.Lmq_sat_l);

      /* End of Outputs for SubSystem: '<S131>/Lmq_sat' */

      /* Switch: '<S131>/Switch1' incorporates:
       *  Constant: '<S131>/Constant2'
       *  Constant: '<S131>/u3'
       */
      if (Diesel_DER2_P.Constant2_Value_o) {
        Diesel_DER2_B.Switch1 = Diesel_DER2_B.Lmq_sat_l.Lmsatq;
      } else {
        Diesel_DER2_B.Switch1 = Diesel_DER2_P.u3_Value_m;
      }

      /* End of Switch: '<S131>/Switch1' */

      /* Assignment: '<S135>/Update Lmq rows[1,5, (6)] col[1,5, (6)] ' incorporates:
       *  Constant: '<S135>/u1'
       */
      memcpy(&rtb_Sum2_m3[0], &Diesel_DER2_P.u1_Value_kt[0], 25U * sizeof(real_T));
      rtb_Sum2_m3[0] = Diesel_DER2_B.Switch1;
      rtb_Sum2_m3[4] = Diesel_DER2_B.Switch1;
      rtb_Sum2_m3[20] = Diesel_DER2_B.Switch1;
      rtb_Sum2_m3[24] = Diesel_DER2_B.Switch1;

      /* Assignment: '<S135>/Update Lmd rows[2,3,4] col[2,3,4]' */
      for (i = 0; i < 3; i++) {
        rtb_Sum2_m3[1 + 5 * (1 + i)] = Diesel_DER2_B.Lmsatd;
        rtb_Sum2_m3[2 + 5 * (1 + i)] = Diesel_DER2_B.Lmsatd;
        rtb_Sum2_m3[3 + 5 * (1 + i)] = Diesel_DER2_B.Lmsatd;
      }

      /* End of Assignment: '<S135>/Update Lmd rows[2,3,4] col[2,3,4]' */

      /* Sum: '<S135>/Sum2' incorporates:
       *  Constant: '<S135>/u5'
       */
      for (i = 0; i < 25; i++) {
        rtb_Sum2_f[i] = rtb_Sum2_m3[i] + Diesel_DER2_P.u5_Value_l[i];
      }

      /* End of Sum: '<S135>/Sum2' */

      /* Product: '<S131>/inversion' */
      rt_invd5x5_snf(rtb_Sum2_f, Diesel_DER2_B.Linv);

      /* Product: '<S131>/Product1' incorporates:
       *  Constant: '<S131>/u1'
       */
      for (i = 0; i < 5; i++) {
        for (i_1 = 0; i_1 < 5; i_1++) {
          Diesel_DER2_B.RLinv[i + 5 * i_1] = 0.0;
          for (i_0 = 0; i_0 < 5; i_0++) {
            Diesel_DER2_B.RLinv[i + 5 * i_1] += Diesel_DER2_P.u1_Value_k[5 * i_0
              + i] * Diesel_DER2_B.Linv[5 * i_1 + i_0];
          }
        }
      }

      /* End of Product: '<S131>/Product1' */

      /* Switch: '<S131>/Switch2' incorporates:
       *  Constant: '<S131>/Constant3'
       *  Constant: '<S131>/u2'
       */
      if (Diesel_DER2_P.Constant3_Value_p) {
        Diesel_DER2_B.Switch2 = Diesel_DER2_B.Lmq_sat_l.MathFunction1;
      } else {
        Diesel_DER2_B.Switch2 = Diesel_DER2_P.u2_Value;
      }

      /* End of Switch: '<S131>/Switch2' */

      /* Update for UnitDelay: '<S133>/Lmd_sat' */
      Diesel_DER2_DW.Lmd_sat_DSTATE = Diesel_DER2_B.Lmsatd;
    }

    /* End of Constant: '<S121>/Constant1' */
    /* End of Outputs for SubSystem: '<S121>/Saturation' */

    /* Switch: '<S121>/Switch' incorporates:
     *  Constant: '<S121>/Constant3'
     *  Constant: '<S121>/Constant4'
     *  Product: '<S121>/Product3'
     */
    tmp = (Diesel_DER2_P.Constant3_Value >= Diesel_DER2_P.Switch_Threshold_o);
    for (i = 0; i < 5; i++) {
      for (i_1 = 0; i_1 < 5; i_1++) {
        if (tmp) {
          rtb_Sum2_f[i_1 + 5 * i] = Diesel_DER2_B.Linv[5 * i + i_1];
        } else {
          rtb_Sum2_f[i_1 + 5 * i] = Diesel_DER2_P.Constant4_Value_m[5 * i + i_1];
        }
      }
    }

    for (i = 0; i < 5; i++) {
      /* Product: '<S121>/Product3' incorporates:
       *  Gain: '<S121>/change Iq Id  current signs'
       */
      tmp_0[i] = 0.0;
      for (i_1 = 0; i_1 < 5; i_1++) {
        tmp_0[i] += rtb_Sum2_f[5 * i_1 + i] * rtb_xk1[i_1];
      }

      /* Gain: '<S121>/change Iq Id  current signs' */
      rtb_changeIqIdcurrentsigns[i] =
        Diesel_DER2_P.changeIqIdcurrentsigns_Gain[i] * tmp_0[i];
    }

    /* End of Switch: '<S121>/Switch' */

    /* Fcn: '<S123>/Fcn' */
    rtb_phimd = rtb_changeIqIdcurrentsigns[0] * rtb_ElementaryMath1 +
      rtb_changeIqIdcurrentsigns[1] * rtb_ElementaryMath;

    /* Fcn: '<S123>/Fcn1' */
    rtb_Fcn1 = ((-rtb_changeIqIdcurrentsigns[0] - 1.7320508075688772 *
                 rtb_changeIqIdcurrentsigns[1]) * rtb_ElementaryMath1 +
                (1.7320508075688772 * rtb_changeIqIdcurrentsigns[0] -
                 rtb_changeIqIdcurrentsigns[1]) * rtb_ElementaryMath) * 0.5;

    /* Sum: '<S123>/Sum' */
    rtb_Sum_nq = (0.0 - rtb_phimd) - rtb_Fcn1;

    /* Switch: '<S121>/Switch3' incorporates:
     *  Constant: '<S121>/Constant8'
     *  Constant: '<S121>/Laqd_nosat'
     */
    if (Diesel_DER2_P.Constant8_Value >= Diesel_DER2_P.Switch3_Threshold) {
      rtb_Mult1_idx_0 = Diesel_DER2_B.Switch2;
      rtb_Mult1_idx_1 = Diesel_DER2_B.MathFunction3;
    } else {
      rtb_Mult1_idx_0 = Diesel_DER2_P.Laqd_nosat_Value[0];
      rtb_Mult1_idx_1 = Diesel_DER2_P.Laqd_nosat_Value[1];
    }

    /* End of Switch: '<S121>/Switch3' */

    /* Product: '<S130>/Product' incorporates:
     *  Gain: '<S130>/1//Ll_q'
     */
    rtb_Product_idx_0 = Diesel_DER2_P.Ll_q_Gain[0] * rtb_xk1[0] *
      rtb_Mult1_idx_0;
    rtb_Product_idx_1 = Diesel_DER2_P.Ll_q_Gain[1] * rtb_xk1[4] *
      rtb_Mult1_idx_0;

    /* Sum: '<S130>/sum1' incorporates:
     *  Gain: '<S130>/1//Ll_d'
     *  Product: '<S130>/Product1'
     */
    rtb_yk_n = (Diesel_DER2_P.Ll_d_Gain_ov[0] * rtb_xk1[1] * rtb_Mult1_idx_1 +
                Diesel_DER2_P.Ll_d_Gain_ov[1] * rtb_xk1[2] * rtb_Mult1_idx_1) +
      Diesel_DER2_P.Ll_d_Gain_ov[2] * rtb_xk1[3] * rtb_Mult1_idx_1;

    /* Gain: '<S123>/ib' */
    Diesel_DER2_B.ib[0] = Diesel_DER2_P.ib_Gain * rtb_phimd;
    Diesel_DER2_B.ib[1] = Diesel_DER2_P.ib_Gain * rtb_Fcn1;
    Diesel_DER2_B.ib[2] = Diesel_DER2_P.ib_Gain * rtb_Sum_nq;

    /* DigitalClock: '<S54>/Digital Clock' */
    rtb_donotdeletethisgain_i = ((Diesel_DER2_M->Timing.clockTick1) * 5.0E-5);

    /* Sum: '<S54>/Sum3' incorporates:
     *  DiscreteIntegrator: '<S54>/Rotor angle dtheta'
     *  Gain: '<S54>/web2'
     */
    rtb_ElementaryMath1_d = Diesel_DER2_P.web2_Gain_e *
      rtb_donotdeletethisgain_i + Diesel_DER2_DW.Rotorangledtheta_DSTATE;

    /* Trigonometry: '<S57>/Elementary Math' */
    rtb_donotdeletethisgain_i = sin(rtb_ElementaryMath1_d);

    /* Trigonometry: '<S57>/Elementary Math1' */
    rtb_ElementaryMath1_d = cos(rtb_ElementaryMath1_d);

    /* UnitDelay: '<S67>/fluxes' */
    for (i = 0; i < 5; i++) {
      rtb_xk1_p[i] = Diesel_DER2_DW.fluxes_DSTATE_o[i];
    }

    /* End of UnitDelay: '<S67>/fluxes' */

    /* Outputs for Enabled SubSystem: '<S56>/Saturation' incorporates:
     *  EnablePort: '<S66>/Enable'
     */
    /* Constant: '<S56>/Constant1' */
    if (Diesel_DER2_P.Constant1_Value_b > 0.0) {
      /* Abs: '<S68>/Abs' incorporates:
       *  Constant: '<S71>/u1'
       *  Gain: '<S72>/1//Ll_d'
       *  Math: '<S71>/Math Function2'
       *  Math: '<S71>/Math Function3'
       *  Product: '<S72>/Product4'
       *  Sum: '<S71>/Sum1'
       *  Sum: '<S72>/sum2'
       *  UnitDelay: '<S68>/Lmd_sat'
       *
       * About '<S71>/Math Function2':
       *  Operator: reciprocal
       *
       * About '<S71>/Math Function3':
       *  Operator: reciprocal
       */
      rtb_Switch_g = fabs(1.0 / (((Diesel_DER2_P.u1_Value[0] +
        Diesel_DER2_P.u1_Value[1]) + Diesel_DER2_P.u1_Value[2]) + 1.0 /
        Diesel_DER2_DW.Lmd_sat_DSTATE_h) * ((Diesel_DER2_P.Ll_d_Gain[0] *
        rtb_xk1_p[1] + Diesel_DER2_P.Ll_d_Gain[1] * rtb_xk1_p[2]) +
        Diesel_DER2_P.Ll_d_Gain[2] * rtb_xk1_p[3]));

      /* Lookup: '<S68>/Lookup Table' */
      rtb_Phisat = rt_Lookup(Diesel_DER2_P.LookupTable_XData, 2, rtb_Switch_g,
        Diesel_DER2_P.LookupTable_YData);

      /* Switch: '<S68>/Switch' incorporates:
       *  Constant: '<S68>/Constant1'
       *  Product: '<S68>/Product'
       */
      if (rtb_Phisat != 0.0) {
        rtb_Switch1_h = rtb_Switch_g / rtb_Phisat;
      } else {
        rtb_Switch1_h = Diesel_DER2_P.Constant1_Value;
      }

      /* End of Switch: '<S68>/Switch' */

      /* Gain: '<S68>/Lmd' */
      rtb_Switch_g = Diesel_DER2_P.Lmd_Gain * rtb_Switch1_h;

      /* Outputs for Enabled SubSystem: '<S66>/Lmq_sat' */

      /* Constant: '<S66>/Constant1' */
      Diesel_DER2_Lmq_sat(Diesel_DER2_P.Constant1_Value_n, rtb_xk1_p,
                          &Diesel_DER2_B.Lmq_sat, &Diesel_DER2_DW.Lmq_sat,
                          (P_Lmq_sat_Diesel_DER2_T *)&Diesel_DER2_P.Lmq_sat);

      /* End of Outputs for SubSystem: '<S66>/Lmq_sat' */

      /* Switch: '<S66>/Switch1' incorporates:
       *  Constant: '<S66>/Constant2'
       *  Constant: '<S66>/u3'
       */
      if (Diesel_DER2_P.Constant2_Value_ea) {
        rtb_Switch1_h = Diesel_DER2_B.Lmq_sat.Lmsatq;
      } else {
        rtb_Switch1_h = Diesel_DER2_P.u3_Value;
      }

      /* End of Switch: '<S66>/Switch1' */

      /* Assignment: '<S70>/Update Lmq rows[1,5, (6)] col[1,5, (6)] ' incorporates:
       *  Constant: '<S70>/u1'
       */
      memcpy(&rtb_Sum2_m3[0], &Diesel_DER2_P.u1_Value_p[0], 25U * sizeof(real_T));
      rtb_Sum2_m3[0] = rtb_Switch1_h;
      rtb_Sum2_m3[4] = rtb_Switch1_h;
      rtb_Sum2_m3[20] = rtb_Switch1_h;
      rtb_Sum2_m3[24] = rtb_Switch1_h;

      /* Assignment: '<S70>/Update Lmd rows[2,3,4] col[2,3,4]' */
      for (i = 0; i < 3; i++) {
        rtb_Sum2_m3[1 + 5 * (1 + i)] = rtb_Switch_g;
        rtb_Sum2_m3[2 + 5 * (1 + i)] = rtb_Switch_g;
        rtb_Sum2_m3[3 + 5 * (1 + i)] = rtb_Switch_g;
      }

      /* End of Assignment: '<S70>/Update Lmd rows[2,3,4] col[2,3,4]' */

      /* Sum: '<S70>/Sum2' incorporates:
       *  Constant: '<S70>/u5'
       */
      for (i = 0; i < 25; i++) {
        rtb_Sum2_f[i] = rtb_Sum2_m3[i] + Diesel_DER2_P.u5_Value[i];
      }

      /* End of Sum: '<S70>/Sum2' */

      /* Product: '<S66>/inversion' */
      rt_invd5x5_snf(rtb_Sum2_f, Diesel_DER2_B.Linv_a);

      /* Product: '<S66>/Product1' incorporates:
       *  Constant: '<S66>/u1'
       */
      for (i = 0; i < 5; i++) {
        for (i_1 = 0; i_1 < 5; i_1++) {
          Diesel_DER2_B.RLinv_e[i + 5 * i_1] = 0.0;
          for (i_0 = 0; i_0 < 5; i_0++) {
            Diesel_DER2_B.RLinv_e[i + 5 * i_1] += Diesel_DER2_P.u1_Value_c[5 *
              i_0 + i] * Diesel_DER2_B.Linv_a[5 * i_1 + i_0];
          }
        }
      }

      /* End of Product: '<S66>/Product1' */

      /* Update for UnitDelay: '<S68>/Lmd_sat' */
      Diesel_DER2_DW.Lmd_sat_DSTATE_h = rtb_Switch_g;
    }

    /* End of Constant: '<S56>/Constant1' */
    /* End of Outputs for SubSystem: '<S56>/Saturation' */

    /* Switch: '<S56>/Switch' incorporates:
     *  Constant: '<S56>/Constant3'
     *  Constant: '<S56>/Constant4'
     *  Product: '<S56>/Product3'
     */
    tmp = (Diesel_DER2_P.Constant3_Value_b >= Diesel_DER2_P.Switch_Threshold_a);
    for (i = 0; i < 5; i++) {
      for (i_1 = 0; i_1 < 5; i_1++) {
        if (tmp) {
          rtb_Sum2_f[i_1 + 5 * i] = Diesel_DER2_B.Linv_a[5 * i + i_1];
        } else {
          rtb_Sum2_f[i_1 + 5 * i] = Diesel_DER2_P.Constant4_Value[5 * i + i_1];
        }
      }
    }

    for (i = 0; i < 5; i++) {
      /* Product: '<S56>/Product3' incorporates:
       *  Gain: '<S56>/change Iq Id  current signs'
       */
      tmp_0[i] = 0.0;
      for (i_1 = 0; i_1 < 5; i_1++) {
        tmp_0[i] += rtb_Sum2_f[5 * i_1 + i] * rtb_xk1_p[i_1];
      }

      /* Gain: '<S56>/change Iq Id  current signs' */
      rtb_changeIqIdcurrentsigns_b[i] =
        Diesel_DER2_P.changeIqIdcurrentsigns_Gain_k[i] * tmp_0[i];
    }

    /* End of Switch: '<S56>/Switch' */

    /* Fcn: '<S58>/Fcn' */
    rtb_donotdeletethisgain_k = rtb_changeIqIdcurrentsigns_b[0] *
      rtb_ElementaryMath1_d + rtb_changeIqIdcurrentsigns_b[1] *
      rtb_donotdeletethisgain_i;

    /* Fcn: '<S58>/Fcn1' */
    rtb_Switch3_e = ((-rtb_changeIqIdcurrentsigns_b[0] - 1.7320508075688772 *
                      rtb_changeIqIdcurrentsigns_b[1]) * rtb_ElementaryMath1_d +
                     (1.7320508075688772 * rtb_changeIqIdcurrentsigns_b[0] -
                      rtb_changeIqIdcurrentsigns_b[1]) *
                     rtb_donotdeletethisgain_i) * 0.5;

    /* Sum: '<S58>/Sum' */
    rtb_DataTypeConversion = (0.0 - rtb_donotdeletethisgain_k) - rtb_Switch3_e;

    /* Gain: '<S58>/ib' */
    Diesel_DER2_B.ib_g[0] = Diesel_DER2_P.ib_Gain_f * rtb_donotdeletethisgain_k;
    Diesel_DER2_B.ib_g[1] = Diesel_DER2_P.ib_Gain_f * rtb_Switch3_e;
    Diesel_DER2_B.ib_g[2] = Diesel_DER2_P.ib_Gain_f * rtb_DataTypeConversion;

    /* UnitDelay: '<S33>/Unit Delay2' */
    Diesel_DER2_B.UnitDelay2 = Diesel_DER2_DW.UnitDelay2_DSTATE;

    /* S-Function block: <S156>/State-Space */
    {
      real_T accum;

      /* Circuit has switches */
      int_T *switch_status = (int_T*)
        Diesel_DER2_DW.StateSpace_PWORK.SWITCH_STATUS;
      int_T *switch_status_init = (int_T*)
        Diesel_DER2_DW.StateSpace_PWORK.SWITCH_STATUS_INIT;
      int_T *SwitchChange = (int_T*) Diesel_DER2_DW.StateSpace_PWORK.SW_CHG;
      int_T *Chopper = (int_T*) Diesel_DER2_DW.StateSpace_PWORK.CHOPPER;
      int_T *gState = (int_T*) Diesel_DER2_DW.StateSpace_PWORK.G_STATE;
      real_T *yswitch = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.Y_SWITCH;
      int_T *switchTypes = (int_T*) Diesel_DER2_DW.StateSpace_PWORK.SWITCH_TYPES;
      int_T *idxOutSw = (int_T*) Diesel_DER2_DW.StateSpace_PWORK.IDX_OUT_SW;
      real_T *DxCol = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.DX_COL;
      real_T *tmp2 = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.TMP2;
      real_T *BDcol = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.BD_COL;
      real_T *tmp1 = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.TMP1;
      int_T newState;
      int_T swChanged = 0;
      int loopsToDo = 20;
      real_T temp;

      /* keep an initial copy of switch_status*/
      memcpy(switch_status_init, switch_status, 9 * sizeof(int_T));
      do {
        if (loopsToDo == 1) {          /* Need to reset some variables: */
          swChanged = 0;

          /* return to the original switch status*/
          {
            int_T i1;
            for (i1=0; i1 < 9; i1++) {
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
          real_T *Cs = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.CS;
          real_T *Ds = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.DS;

          {
            int_T i1;
            real_T *y0 = &Diesel_DER2_B.StateSpace_o1[0];
            for (i1=0; i1 < 42; i1++) {
              accum = 0.0;

              {
                int_T i2;
                real_T *xd = &Diesel_DER2_DW.StateSpace_DSTATE[0];
                for (i2=0; i2 < 6; i2++) {
                  accum += *(Cs++) * xd[i2];
                }
              }

              {
                int_T i2;
                const real_T *u0 = &Diesel_DER2_P.SwitchCurrents_Value[0];
                for (i2=0; i2 < 9; i2++) {
                  accum += *(Ds++) * u0[i2];
                }

                accum += *(Ds++) * Diesel_DER2_B.ib[0];
                accum += *(Ds++) * Diesel_DER2_B.ib[1];
                accum += *(Ds++) * Diesel_DER2_B.ib_g[0];
                accum += *(Ds++) * Diesel_DER2_B.ib_g[1];
                accum += *(Ds++) * Diesel_DER2_B.UnitDelay2;
              }

              y0[i1] = accum * Chopper[i1];
            }
          }

          swChanged = 0;

          {
            int_T i1;
            real_T *y0 = &Diesel_DER2_B.StateSpace_o1[0];
            for (i1=0; i1 < 9; i1++) {
              switch (switchTypes[i1]) {
               case 2 :                /* Breaker */
                newState = gState[i1] > 0 ? 1 : 0;
                break;

               case 3 :                /* Diodes */
                newState = y0[i1] > 0.0 ? 1 : ((y0[i1] < 0.0) ? 0 :
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
          real_T *As = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.AS;
          real_T *Cs = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.CS;
          real_T *Bs = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.BS;
          real_T *Ds = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.DS;
          real_T a1;

          {
            int_T i1;
            for (i1=0; i1 < 9; i1++) {
              if (SwitchChange[i1] != 0) {
                if (idxOutSw[i1] > -1 ) {/* A positive index points to a switch actual measure output */
                  Chopper[idxOutSw[i1]] = switch_status[i1];
                }

                a1 = yswitch[i1]*SwitchChange[i1];
                temp = 1/(1-Ds[i1*15]*a1);

                {
                  int_T i2;
                  for (i2=0; i2 < 42; i2++) {
                    DxCol[i2]= Ds[i2 * 14 + i1]*temp*a1;
                  }
                }

                DxCol[i1] = temp;

                {
                  int_T i2;
                  for (i2=0; i2 < 6; i2++) {
                    BDcol[i2]= Bs[i2 * 14 + i1]*a1;
                  }
                }

                /* Copy row nSw of Cs into tmp1 and zero it out in Cs */
                memcpy(tmp1, &Cs[i1 * 6], 6 * sizeof(real_T));
                memset(&Cs[i1 * 6], '\0', 6 * sizeof(real_T));

                /* Copy row nSw of Ds into tmp2 and zero it out in Ds */
                memcpy(tmp2, &Ds[i1 * 14], 14 * sizeof(real_T));
                memset(&Ds[i1 * 14], '\0', 14 * sizeof(real_T));

                /* Cs = Cs + DxCol * tmp1, Ds = Ds + DxCol * tmp2 *******************/
                {
                  int_T i2;
                  for (i2=0; i2 < 42; i2++) {
                    a1 = DxCol[i2];

                    {
                      int_T i3;
                      for (i3=0; i3 < 6; i3++) {
                        Cs[i2 * 6 + i3] += a1 * tmp1[i3];
                      }
                    }

                    {
                      int_T i3;
                      for (i3=0; i3 < 14; i3++) {
                        Ds[i2 * 14 + i3] += a1 * tmp2[i3];
                      }
                    }
                  }
                }

                /* As = As + BdCol*Cs(nSw,:), Bs = Bs + BdCol*Ds(nSw,:) *************/
                {
                  int_T i2;
                  for (i2=0; i2 < 6; i2++) {
                    a1 = BDcol[i2];

                    {
                      int_T i3;
                      for (i3=0; i3 < 6; i3++) {
                        As[i2 * 6 + i3] += a1 * Cs[i1 * 6 + i3];
                      }
                    }

                    {
                      int_T i3;
                      for (i3=0; i3 < 14; i3++) {
                        Bs[i2 * 14 + i3] += a1 * Ds[i1 * 14 + i3];
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
        real_T *Cs = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.CS;
        real_T *Ds = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.DS;

        {
          int_T i1;
          real_T *y0 = &Diesel_DER2_B.StateSpace_o1[0];
          for (i1=0; i1 < 42; i1++) {
            accum = 0.0;

            {
              int_T i2;
              real_T *xd = &Diesel_DER2_DW.StateSpace_DSTATE[0];
              for (i2=0; i2 < 6; i2++) {
                accum += *(Cs++) * xd[i2];
              }
            }

            {
              int_T i2;
              const real_T *u0 = &Diesel_DER2_P.SwitchCurrents_Value[0];
              for (i2=0; i2 < 9; i2++) {
                accum += *(Ds++) * u0[i2];
              }

              accum += *(Ds++) * Diesel_DER2_B.ib[0];
              accum += *(Ds++) * Diesel_DER2_B.ib[1];
              accum += *(Ds++) * Diesel_DER2_B.ib_g[0];
              accum += *(Ds++) * Diesel_DER2_B.ib_g[1];
              accum += *(Ds++) * Diesel_DER2_B.UnitDelay2;
            }

            y0[i1] = accum * Chopper[i1];
          }
        }
      }

      /* Output new switches states */
      {
        int_T i1;
        real_T *y1 = &Diesel_DER2_B.StateSpace_o2[0];
        for (i1=0; i1 < 9; i1++) {
          y1[i1] = (real_T)switch_status[i1];
        }
      }
    }

    /* Gain: '<S122>/1_Vb' */
    rtb_Mult1_idx_0 = Diesel_DER2_P._Vb_Gain * Diesel_DER2_B.StateSpace_o1[9];
    rtb_Mult1_idx_1 = Diesel_DER2_P._Vb_Gain * Diesel_DER2_B.StateSpace_o1[10];

    /* Fcn: '<S122>/Fcn2' */
    rtb_Fcn2 = ((2.0 * rtb_Mult1_idx_0 + rtb_Mult1_idx_1) * rtb_ElementaryMath1
                + 1.7320508075688772 * rtb_Mult1_idx_1 * rtb_ElementaryMath) *
      0.33333333333333331;

    /* Fcn: '<S122>/Fcn3' */
    rtb_Fcn3 = ((2.0 * rtb_Mult1_idx_0 + rtb_Mult1_idx_1) * rtb_ElementaryMath +
                -1.7320508075688772 * rtb_Mult1_idx_1 * rtb_ElementaryMath1) *
      0.33333333333333331;

    /* Fcn: '<S126>/x->theta' */
    rtb_DataTypeConversion = rt_atan2d_snf(rtb_Fcn3, rtb_Fcn2);

    /* Gain: '<S127>/Gain' */
    rtb_DataTypeConversion *= Diesel_DER2_P.Gain_Gain_k;

    /* Gain: '<S125>/Gain' incorporates:
     *  Product: '<S125>/Product2'
     *  Product: '<S125>/Product3'
     *  Sum: '<S125>/Sum'
     */
    rtb_Switch3_e = (rtb_Fcn2 * rtb_changeIqIdcurrentsigns[0] + rtb_Fcn3 *
                     rtb_changeIqIdcurrentsigns[1]) * Diesel_DER2_P.Gain_Gain_g;

    /* Gain: '<S125>/Gain1' incorporates:
     *  Product: '<S125>/Product4'
     *  Product: '<S125>/Product5'
     *  Sum: '<S125>/Sum1'
     */
    rtb_donotdeletethisgain_k = (rtb_Fcn2 * rtb_changeIqIdcurrentsigns[1] -
      rtb_Fcn3 * rtb_changeIqIdcurrentsigns[0]) * Diesel_DER2_P.Gain1_Gain_p;

    /* Gain: '<S120>/output formatting' */
    rtb_outputformatting[0] = Diesel_DER2_P.outputformatting_Gain[0] * rtb_phimd;
    rtb_outputformatting[1] = Diesel_DER2_P.outputformatting_Gain[1] * rtb_Fcn1;
    rtb_outputformatting[2] = Diesel_DER2_P.outputformatting_Gain[2] *
      rtb_Sum_nq;
    rtb_outputformatting[3] = Diesel_DER2_P.outputformatting_Gain[3] *
      rtb_changeIqIdcurrentsigns[0];
    rtb_outputformatting[4] = Diesel_DER2_P.outputformatting_Gain[4] *
      rtb_changeIqIdcurrentsigns[1];
    rtb_outputformatting[5] = Diesel_DER2_P.outputformatting_Gain[5] *
      rtb_changeIqIdcurrentsigns[2];
    rtb_outputformatting[6] = Diesel_DER2_P.outputformatting_Gain[6] *
      rtb_changeIqIdcurrentsigns[4];

    /* Switch: '<S120>/Switch' incorporates:
     *  Constant: '<S120>/Constant1'
     *  Constant: '<S120>/Constant2'
     */
    if (Diesel_DER2_P.Constant1_Value_a) {
      rtb_ElementaryMath1 = rtb_changeIqIdcurrentsigns[0];
    } else {
      rtb_ElementaryMath1 = Diesel_DER2_P.Constant2_Value_l;
    }

    /* End of Switch: '<S120>/Switch' */

    /* Gain: '<S120>/output formatting' incorporates:
     *  Sum: '<S130>/sum2'
     */
    rtb_outputformatting[7] = Diesel_DER2_P.outputformatting_Gain[7] *
      rtb_ElementaryMath1;
    rtb_outputformatting[8] = Diesel_DER2_P.outputformatting_Gain[8] *
      rtb_changeIqIdcurrentsigns[3];
    rtb_outputformatting[9] = (rtb_Product_idx_0 + rtb_Product_idx_1) *
      Diesel_DER2_P.outputformatting_Gain[9];
    rtb_outputformatting[10] = Diesel_DER2_P.outputformatting_Gain[10] *
      rtb_yk_n;
    rtb_outputformatting[11] = Diesel_DER2_P.outputformatting_Gain[11] *
      rtb_Fcn2;
    rtb_outputformatting[12] = Diesel_DER2_P.outputformatting_Gain[12] *
      rtb_Fcn3;

    /* Switch: '<S121>/Switch2' incorporates:
     *  Constant: '<S121>/Constant5'
     *  Constant: '<S121>/Lmqd_nosat'
     */
    if (Diesel_DER2_P.Constant5_Value >= Diesel_DER2_P.Switch2_Threshold_j) {
      rtb_ElementaryMath1 = Diesel_DER2_B.Switch1;
    } else {
      rtb_ElementaryMath1 = Diesel_DER2_P.Lmqd_nosat_Value[0];
    }

    /* Gain: '<S120>/output formatting' */
    rtb_outputformatting[13] = Diesel_DER2_P.outputformatting_Gain[13] *
      rtb_ElementaryMath1;

    /* Switch: '<S121>/Switch2' incorporates:
     *  Constant: '<S121>/Constant5'
     *  Constant: '<S121>/Lmqd_nosat'
     */
    if (Diesel_DER2_P.Constant5_Value >= Diesel_DER2_P.Switch2_Threshold_j) {
      rtb_ElementaryMath1 = Diesel_DER2_B.Lmsatd;
    } else {
      rtb_ElementaryMath1 = Diesel_DER2_P.Lmqd_nosat_Value[1];
    }

    /* Gain: '<S120>/output formatting' */
    rtb_outputformatting[14] = Diesel_DER2_P.outputformatting_Gain[14] *
      rtb_ElementaryMath1;
    rtb_outputformatting[15] = Diesel_DER2_P.outputformatting_Gain[15] *
      rtb_DataTypeConversion;
    rtb_outputformatting[16] = Diesel_DER2_P.outputformatting_Gain[16] *
      rtb_Switch3_e;
    rtb_outputformatting[17] = Diesel_DER2_P.outputformatting_Gain[17] *
      rtb_donotdeletethisgain_k;

    /* Outport: '<Root>/out_If' */
    Diesel_DER2_Y.out_If = rtb_outputformatting[5];

    /* TransportDelay: '<S2>/Delay Td' */
    {
      real_T **uBuffer = (real_T**)&Diesel_DER2_DW.DelayTd_PWORK.TUbufferPtrs[0];
      real_T **tBuffer = (real_T**)&Diesel_DER2_DW.DelayTd_PWORK.TUbufferPtrs[1];
      real_T simTime = Diesel_DER2_M->Timing.t[0];
      real_T tMinusDelay = simTime - Diesel_DER2_P.DieselEngineSpeedRegulator_Td;
      rtb_pu2Watts = rt_TDelayInterpolate(
        tMinusDelay,
        0.0,
        *tBuffer,
        *uBuffer,
        Diesel_DER2_DW.DelayTd_IWORK.CircularBufSize,
        &Diesel_DER2_DW.DelayTd_IWORK.Last,
        Diesel_DER2_DW.DelayTd_IWORK.Tail,
        Diesel_DER2_DW.DelayTd_IWORK.Head,
        Diesel_DER2_P.DieselEngineSpeedRegulator_Pm0,
        1,
        0);
    }

    /* UnitDelay: '<S141>/dw_delay' */
    rtb_dw_delay = Diesel_DER2_DW.dw_delay_DSTATE;

    /* Sum: '<S141>/Sum1' incorporates:
     *  Gain: '<S141>/F2'
     *  UnitDelay: '<S141>/dw_predict'
     */
    rtb_DataTypeConversion = Diesel_DER2_P.F2_Gain * rtb_dw_delay -
      Diesel_DER2_DW.dw_predict_DSTATE;

    /* Sum: '<S119>/Sum' incorporates:
     *  Constant: '<S119>/nominal speed'
     */
    rtb_Phisat = Diesel_DER2_P.nominalspeed_Value + rtb_DataTypeConversion;

    /* Gain: '<S119>/units' */
    rtb_Switch3_e = Diesel_DER2_P.units_Gain * rtb_Phisat;

    /* Outport: '<Root>/out_w_pu' incorporates:
     *  Gain: '<S2>/rad//s to pu'
     */
    Diesel_DER2_Y.out_w_pu = Diesel_DER2_P.radstopu_Gain * rtb_Switch3_e;

    /* Gain: '<S2>/pu2Watts' incorporates:
     *  Product: '<S2>/Product'
     */
    rtb_pu2Watts = rtb_pu2Watts * Diesel_DER2_Y.out_w_pu *
      Diesel_DER2_P.pu2Watts_Gain;

    /* Gain: '<S54>/1_Pb' */
    rtb_phimd = Diesel_DER2_P._Pb_Gain * rtb_Switch3_e;

    /* Gain: '<S76>/units1' incorporates:
     *  Gain: '<S63>/1-1'
     *  Product: '<S63>/Mult1'
     *  Product: '<S76>/Product'
     *  Sum: '<S63>/Sum2'
     */
    rtb_Switch_g = (Diesel_DER2_P.u_Gain[0] * rtb_changeIqIdcurrentsigns_b[0] *
                    rtb_xk1_p[1] + Diesel_DER2_P.u_Gain[1] *
                    rtb_changeIqIdcurrentsigns_b[1] * rtb_xk1_p[0]) * rtb_phimd *
      Diesel_DER2_P.units1_Gain;

    /* Outport: '<Root>/out_Pmec' incorporates:
     *  Sum: '<Root>/Sum2'
     */
    Diesel_DER2_Y.out_Pmec[0] = rtb_pu2Watts;
    Diesel_DER2_Y.out_Pmec[1] = rtb_pu2Watts - rtb_Switch_g;

    /* Sum: '<S2>/Sum1' incorporates:
     *  Constant: '<Root>/Reference speed (pu)'
     */
    rtb_Switch3_e = Diesel_DER2_P.Referencespeedpu_Value -
      Diesel_DER2_Y.out_w_pu;

    /* Gain: '<S2>/Gain K' incorporates:
     *  Gain: '<S29>/C'
     *  Gain: '<S29>/D'
     *  Sum: '<S29>/sum1'
     *  UnitDelay: '<S29>/Delay_x'
     */
    rtb_donotdeletethisgain_k = (Diesel_DER2_P.D_Gain * rtb_Switch3_e +
      Diesel_DER2_P.C_Gain * Diesel_DER2_DW.Delay_x_DSTATE_i) *
      Diesel_DER2_P.DieselEngineSpeedRegulator_K;

    /* Sum: '<S30>/sum1' incorporates:
     *  Gain: '<S30>/C'
     *  Gain: '<S30>/D'
     *  UnitDelay: '<S30>/Delay_x'
     */
    rtb_yk_n = Diesel_DER2_P.D_Gain_i * rtb_donotdeletethisgain_k +
      Diesel_DER2_P.C_Gain_g * Diesel_DER2_DW.Delay_x_DSTATE_g;

    /* Sum: '<S28>/A*x(k) + B*u(k)' incorporates:
     *  Gain: '<S28>/A'
     *  Gain: '<S28>/B'
     *  UnitDelay: '<S28>/Delay_x'
     */
    rtb_xk1_c = Diesel_DER2_P.A_Gain * Diesel_DER2_DW.Delay_x_DSTATE +
      Diesel_DER2_P.B_Gain * rtb_yk_n;

    /* Sum: '<S28>/sum1' incorporates:
     *  Gain: '<S28>/C'
     *  Gain: '<S28>/D'
     *  UnitDelay: '<S28>/Delay_x'
     */
    rtb_yk_m = Diesel_DER2_P.D_Gain_j * rtb_yk_n + Diesel_DER2_P.C_Gain_a *
      Diesel_DER2_DW.Delay_x_DSTATE;

    /* DiscreteIntegrator: '<S2>/Discrete-Time Integrator' */
    Diesel_DER2_B.Gazflow = Diesel_DER2_DW.DiscreteTimeIntegrator_DSTATE;

    /* Sum: '<S29>/A*x(k) + B*u(k)' incorporates:
     *  Gain: '<S29>/A'
     *  Gain: '<S29>/B'
     *  UnitDelay: '<S29>/Delay_x'
     */
    rtb_xk1_h = Diesel_DER2_P.A_Gain_o * Diesel_DER2_DW.Delay_x_DSTATE_i +
      Diesel_DER2_P.B_Gain_o * rtb_Switch3_e;

    /* Sum: '<S30>/A*x(k) + B*u(k)' incorporates:
     *  Gain: '<S30>/A'
     *  Gain: '<S30>/B'
     *  UnitDelay: '<S30>/Delay_x'
     */
    rtb_xk1_i = Diesel_DER2_P.A_Gain_e * Diesel_DER2_DW.Delay_x_DSTATE_g +
      Diesel_DER2_P.B_Gain_n * rtb_donotdeletethisgain_k;

    /* UnitDelay: '<S33>/Unit Delay4' */
    rtb_UnitDelay4 = Diesel_DER2_DW.UnitDelay4_DSTATE;

    /* Gain: '<S78>/do not delete this gain' */
    rtb_donotdeletethisgain = Diesel_DER2_P.donotdeletethisgain_Gain *
      Diesel_DER2_B.StateSpace_o1[16];

    /* Sum: '<S36>/Sum' incorporates:
     *  Constant: '<S3>/Reference voltage (pu)'
     *  UnitDelay: '<S36>/Unit Delay1'
     */
    rtb_Switch3_e = Diesel_DER2_P.Referencevoltagepu_Value -
      Diesel_DER2_DW.UnitDelay1_DSTATE;

    /* DiscreteIntegrator: '<S85>/Discrete-Time Integrator' */
    rtb_donotdeletethisgain_k = Diesel_DER2_DW.DiscreteTimeIntegrator_DSTATE_f;

    /* Sum: '<S85>/Sum6' incorporates:
     *  Gain: '<S85>/Kp4'
     */
    rtb_donotdeletethisgain_k += Diesel_DER2_P.DiscretePIController_Kp *
      rtb_Switch3_e;

    /* Saturate: '<S85>/Saturation2' */
    if (rtb_donotdeletethisgain_k > Diesel_DER2_P.Saturation2_UpperSat) {
      rtb_ElementaryMath1 = Diesel_DER2_P.Saturation2_UpperSat;
    } else if (rtb_donotdeletethisgain_k < Diesel_DER2_P.Saturation2_LowerSat) {
      rtb_ElementaryMath1 = Diesel_DER2_P.Saturation2_LowerSat;
    } else {
      rtb_ElementaryMath1 = rtb_donotdeletethisgain_k;
    }

    /* Gain: '<S52>/N' incorporates:
     *  Gain: '<S36>/pu to volts'
     *  Saturate: '<S85>/Saturation2'
     */
    rtb_Vfd = Diesel_DER2_P.putovolts_Gain * rtb_ElementaryMath1 *
      Diesel_DER2_P.N_Gain;

    /* Gain: '<S57>/1_Vb' */
    rtb_Mult1_idx_0 = Diesel_DER2_P._Vb_Gain_p * Diesel_DER2_B.StateSpace_o1[11];
    rtb_Mult1_idx_1 = Diesel_DER2_P._Vb_Gain_p * Diesel_DER2_B.StateSpace_o1[12];

    /* Fcn: '<S57>/Fcn2' */
    rtb_Fcn2_j = ((2.0 * rtb_Mult1_idx_0 + rtb_Mult1_idx_1) *
                  rtb_ElementaryMath1_d + 1.7320508075688772 * rtb_Mult1_idx_1 *
                  rtb_donotdeletethisgain_i) * 0.33333333333333331;

    /* Fcn: '<S57>/Fcn3' */
    rtb_Fcn3_m = ((2.0 * rtb_Mult1_idx_0 + rtb_Mult1_idx_1) *
                  rtb_donotdeletethisgain_i + -1.7320508075688772 *
                  rtb_Mult1_idx_1 * rtb_ElementaryMath1_d) * 0.33333333333333331;

    /* DigitalClock: '<S67>/Digital Clock' */
    rtb_donotdeletethisgain_k = ((Diesel_DER2_M->Timing.clockTick1) * 5.0E-5);

    /* Switch: '<S67>/IC' */
    if (rtb_donotdeletethisgain_k >= Diesel_DER2_P.IC_Threshold) {
      /* Assignment: '<S64>/W(1,2)=wr' incorporates:
       *  Constant: '<S64>/u1'
       */
      memcpy(&rtb_Sum2_m3[0], &Diesel_DER2_P.u1_Value_e[0], 25U * sizeof(real_T));
      rtb_Sum2_m3[5] = rtb_phimd;

      /* Assignment: '<S64>/W(2,1)=-wr' incorporates:
       *  Gain: '<S64>/Gain1'
       */
      rtb_Sum2_m3[1] = Diesel_DER2_P.Gain1_Gain * rtb_phimd;

      /* Switch: '<S56>/Switch1' incorporates:
       *  Constant: '<S56>/Constant2'
       */
      tmp = (Diesel_DER2_P.Constant2_Value >= Diesel_DER2_P.Switch1_Threshold);
      for (i = 0; i < 25; i++) {
        /* Gain: '<S75>/wbase*Ts//2' incorporates:
         *  Constant: '<S56>/Constant6'
         *  Sum: '<S56>/Sum1'
         *  Switch: '<S56>/Switch1'
         */
        if (tmp) {
          rtb_ElementaryMath1 = Diesel_DER2_B.RLinv_e[i];
        } else {
          rtb_ElementaryMath1 = Diesel_DER2_P.Constant6_Value[i];
        }

        rtb_ElementaryMath1 = ((0.0 - rtb_Sum2_m3[i]) - rtb_ElementaryMath1) *
          Diesel_DER2_P.wbaseTs2_Gain;

        /* Sum: '<S75>/Sum1' incorporates:
         *  Constant: '<S75>/u5'
         */
        rtb_Sum2_f[i] = Diesel_DER2_P.u5_Value_g[i] - rtb_ElementaryMath1;

        /* Sum: '<S75>/Sum5' incorporates:
         *  Constant: '<S75>/u5'
         */
        rtb_ElementaryMath1 += Diesel_DER2_P.u5_Value_g[i];
        rtb_Sum2_m3[i] = rtb_ElementaryMath1;
      }

      /* Product: '<S75>/inversion' incorporates:
       *  Gain: '<S75>/wbase*Ts//2'
       *  Sum: '<S56>/Sum1'
       */
      rt_invd5x5_snf(rtb_Sum2_f, rtb_inversion_g);

      /* Product: '<S75>/Product4' incorporates:
       *  Product: '<S67>/Product2'
       */
      for (i = 0; i < 5; i++) {
        for (i_1 = 0; i_1 < 5; i_1++) {
          rtb_Sum2_f[i + 5 * i_1] = 0.0;
          for (i_0 = 0; i_0 < 5; i_0++) {
            rtb_Sum2_f[i + 5 * i_1] += rtb_inversion_g[5 * i_0 + i] *
              rtb_Sum2_m3[5 * i_1 + i_0];
          }
        }
      }

      /* End of Product: '<S75>/Product4' */

      /* Sum: '<S67>/sum' incorporates:
       *  Constant: '<S52>/[ Vkd =0 Vkq1=0  Vkq2=0 ]'
       */
      tmp_1[0] = rtb_Fcn2_j;
      tmp_1[1] = rtb_Fcn3_m;
      tmp_1[2] = rtb_Vfd;
      tmp_1[3] = Diesel_DER2_P.Vkd0Vkq10Vkq20_Value[0];
      tmp_1[4] = Diesel_DER2_P.Vkd0Vkq10Vkq20_Value[1];
      for (i = 0; i < 5; i++) {
        /* Sum: '<S67>/sum' incorporates:
         *  Product: '<S67>/Product1'
         *  UnitDelay: '<S67>/voltages'
         */
        tmp_0[i] = tmp_1[i] + Diesel_DER2_DW.voltages_DSTATE[i];

        /* Product: '<S67>/Product2' incorporates:
         *  Sum: '<S67>/Ad*x(k-1) + Bd*( u(k-1) + u(k))'
         */
        rtb_changeIqIdcurrentsigns_b[i] = 0.0;
        for (i_1 = 0; i_1 < 5; i_1++) {
          rtb_changeIqIdcurrentsigns_b[i] += rtb_Sum2_f[5 * i_1 + i] *
            rtb_xk1_p[i_1];
        }
      }

      /* Product: '<S67>/Product1' incorporates:
       *  Gain: '<S75>/wbase*Ts//2 '
       *  Sum: '<S67>/Ad*x(k-1) + Bd*( u(k-1) + u(k))'
       */
      for (i = 0; i < 5; i++) {
        tmp_1[i] = 0.0;
        for (i_1 = 0; i_1 < 5; i_1++) {
          tmp_1[i] += rtb_inversion_g[5 * i_1 + i] *
            Diesel_DER2_P.wbaseTs2_Gain_d * tmp_0[i_1];
        }

        rtb_IC[i] = rtb_changeIqIdcurrentsigns_b[i] + tmp_1[i];
      }
    } else {
      for (i = 0; i < 5; i++) {
        rtb_IC[i] = rtb_xk1_p[i];
      }
    }

    /* End of Switch: '<S67>/IC' */

    /* Gain: '<S54>/web1' incorporates:
     *  Constant: '<S54>/Nominal speed'
     *  Sum: '<S54>/Add'
     */
    rtb_web1 = (rtb_phimd - Diesel_DER2_P.Nominalspeed_Value) *
      Diesel_DER2_P.web1_Gain;

    /* Gain: '<S85>/Kp5' */
    rtb_Kp5 = Diesel_DER2_P.DiscretePIController_Ki * rtb_Switch3_e;

    /* Gain: '<S16>/do not delete this gain' */
    rtb_Switch3_e = Diesel_DER2_P.donotdeletethisgain_Gain_h *
      Diesel_DER2_B.StateSpace_o1[13];

    /* Gain: '<S17>/do not delete this gain' */
    rtb_donotdeletethisgain_k = Diesel_DER2_P.donotdeletethisgain_Gain_l *
      Diesel_DER2_B.StateSpace_o1[14];

    /* Gain: '<S18>/do not delete this gain' */
    rtb_donotdeletethisgain_i = Diesel_DER2_P.donotdeletethisgain_Gain_ho *
      Diesel_DER2_B.StateSpace_o1[15];

    /* Gain: '<S1>/Kv1' */
    rtb_Switch1_h = Diesel_DER2_P.Kv1_Gain * rtb_Switch3_e;
    rtb_Mult1_idx_0 = Diesel_DER2_P.Kv1_Gain * rtb_donotdeletethisgain_k;
    rtb_Mult1_idx_1 = Diesel_DER2_P.Kv1_Gain * rtb_donotdeletethisgain_i;

    /* Sin: '<S87>/sin(wt)' */
    if (Diesel_DER2_DW.systemEnable != 0) {
      Diesel_DER2_DW.lastSin = sin(6.2831853071795862 *
        Diesel_DER2_P.Fourier_A_Freq * Diesel_DER2_M->Timing.t[1]);
      Diesel_DER2_DW.lastCos = cos(6.2831853071795862 *
        Diesel_DER2_P.Fourier_A_Freq * Diesel_DER2_M->Timing.t[1]);
      Diesel_DER2_DW.systemEnable = 0;
    }

    rtb_Switch3_e = ((Diesel_DER2_DW.lastSin * Diesel_DER2_P.sinwt_PCos +
                      Diesel_DER2_DW.lastCos * Diesel_DER2_P.sinwt_PSin) *
                     Diesel_DER2_P.sinwt_HCos + (Diesel_DER2_DW.lastCos *
      Diesel_DER2_P.sinwt_PCos - Diesel_DER2_DW.lastSin *
      Diesel_DER2_P.sinwt_PSin) * Diesel_DER2_P.sinwt_Hsin) *
      Diesel_DER2_P.sinwt_Amp + Diesel_DER2_P.sinwt_Bias;

    /* End of Sin: '<S87>/sin(wt)' */

    /* Product: '<S87>/Product' */
    rtb_Product_n = rtb_Switch1_h * rtb_Switch3_e;

    /* DiscreteIntegrator: '<S97>/Integ4' */
    if (Diesel_DER2_DW.Integ4_SYSTEM_ENABLE != 0) {
      Diesel_DER2_B.Integ4 = Diesel_DER2_DW.Integ4_DSTATE;
    } else {
      Diesel_DER2_B.Integ4 = Diesel_DER2_P.Integ4_gainval * rtb_Product_n +
        Diesel_DER2_DW.Integ4_DSTATE;
    }

    /* End of DiscreteIntegrator: '<S97>/Integ4' */

    /* Constant: '<S97>/K1' */
    Diesel_DER2_B.K1 = Diesel_DER2_P.K1_Value;

    /* Level2 S-Function Block: '<S98>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[0];
      sfcnOutputs(rts, 1);
    }

    /* DigitalClock: '<S97>/Digital  Clock' */
    rtb_Switch3_e = Diesel_DER2_M->Timing.t[1];

    /* Switch: '<S97>/Switch' incorporates:
     *  Constant: '<S97>/K2'
     *  Gain: '<S97>/Gain'
     *  Gain: '<S97>/Gain1'
     *  Product: '<S97>/Product'
     *  RelationalOperator: '<S97>/Relational Operator'
     *  Sum: '<S97>/Sum1'
     *  Sum: '<S97>/Sum5'
     *  Sum: '<S97>/Sum7'
     *  UnitDelay: '<S97>/Unit Delay'
     *  UnitDelay: '<S97>/Unit Delay1'
     */
    if (rtb_Switch3_e >= Diesel_DER2_B.K1) {
      rtb_Switch = (Diesel_DER2_B.Integ4 - Diesel_DER2_B.SFunction) *
        Diesel_DER2_P.K2_Value + (Diesel_DER2_P.Gain1_Gain_m * rtb_Product_n -
        Diesel_DER2_P.Gain_Gain_n * Diesel_DER2_DW.UnitDelay_DSTATE);
    } else {
      rtb_Switch = Diesel_DER2_DW.UnitDelay1_DSTATE_c;
    }

    /* End of Switch: '<S97>/Switch' */

    /* Sin: '<S87>/cos(wt)' */
    if (Diesel_DER2_DW.systemEnable_b != 0) {
      Diesel_DER2_DW.lastSin_d = sin(6.2831853071795862 *
        Diesel_DER2_P.Fourier_A_Freq * Diesel_DER2_M->Timing.t[1]);
      Diesel_DER2_DW.lastCos_e = cos(6.2831853071795862 *
        Diesel_DER2_P.Fourier_A_Freq * Diesel_DER2_M->Timing.t[1]);
      Diesel_DER2_DW.systemEnable_b = 0;
    }

    rtb_Switch3_e = ((Diesel_DER2_DW.lastSin_d * Diesel_DER2_P.coswt_PCos +
                      Diesel_DER2_DW.lastCos_e * Diesel_DER2_P.coswt_PSin) *
                     Diesel_DER2_P.coswt_HCos + (Diesel_DER2_DW.lastCos_e *
      Diesel_DER2_P.coswt_PCos - Diesel_DER2_DW.lastSin_d *
      Diesel_DER2_P.coswt_PSin) * Diesel_DER2_P.coswt_Hsin) *
      Diesel_DER2_P.coswt_Amp + Diesel_DER2_P.coswt_Bias;

    /* End of Sin: '<S87>/cos(wt)' */

    /* Product: '<S87>/Product1' */
    rtb_Product1_f = rtb_Switch1_h * rtb_Switch3_e;

    /* DiscreteIntegrator: '<S95>/Integ4' */
    if (Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_f != 0) {
      Diesel_DER2_B.Integ4_h = Diesel_DER2_DW.Integ4_DSTATE_e;
    } else {
      Diesel_DER2_B.Integ4_h = Diesel_DER2_P.Integ4_gainval_o * rtb_Product1_f +
        Diesel_DER2_DW.Integ4_DSTATE_e;
    }

    /* End of DiscreteIntegrator: '<S95>/Integ4' */

    /* Constant: '<S95>/K1' */
    Diesel_DER2_B.K1_e = Diesel_DER2_P.K1_Value_f;

    /* Level2 S-Function Block: '<S96>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[1];
      sfcnOutputs(rts, 1);
    }

    /* DigitalClock: '<S95>/Digital  Clock' */
    rtb_Switch3_e = Diesel_DER2_M->Timing.t[1];

    /* Switch: '<S95>/Switch' incorporates:
     *  Constant: '<S95>/K2'
     *  Gain: '<S95>/Gain'
     *  Gain: '<S95>/Gain1'
     *  Product: '<S95>/Product'
     *  RelationalOperator: '<S95>/Relational Operator'
     *  Sum: '<S95>/Sum1'
     *  Sum: '<S95>/Sum5'
     *  Sum: '<S95>/Sum7'
     *  UnitDelay: '<S95>/Unit Delay'
     *  UnitDelay: '<S95>/Unit Delay1'
     */
    if (rtb_Switch3_e >= Diesel_DER2_B.K1_e) {
      rtb_Switch_k = (Diesel_DER2_B.Integ4_h - Diesel_DER2_B.SFunction_a) *
        Diesel_DER2_P.K2_Value_p + (Diesel_DER2_P.Gain1_Gain_n * rtb_Product1_f
        - Diesel_DER2_P.Gain_Gain * Diesel_DER2_DW.UnitDelay_DSTATE_c);
    } else {
      rtb_Switch_k = Diesel_DER2_DW.UnitDelay1_DSTATE_l;
    }

    /* End of Switch: '<S95>/Switch' */

    /* ComplexToMagnitudeAngle: '<S87>/Complex to Magnitude-Angle' incorporates:
     *  RealImagToComplex: '<S87>/Real-Imag to Complex'
     */
    rtb_Switch3_e = rt_hypotd_snf(rtb_Switch, rtb_Switch_k);
    rtb_donotdeletethisgain_k = rt_atan2d_snf(rtb_Switch_k, rtb_Switch);

    /* Gain: '<S86>/deg->rad' incorporates:
     *  Gain: '<S87>/Rad->Deg.'
     */
    rtb_donotdeletethisgain_k = Diesel_DER2_P.RadDeg_Gain *
      rtb_donotdeletethisgain_k * Diesel_DER2_P.degrad_Gain;

    /* MagnitudeAngleToComplex: '<S86>/Magnitude-Angle to Complex' */
    rtb_MagnitudeAngletoComplex.re = rtb_Switch3_e * cos
      (rtb_donotdeletethisgain_k);
    rtb_MagnitudeAngletoComplex.im = rtb_Switch3_e * sin
      (rtb_donotdeletethisgain_k);

    /* Sin: '<S88>/sin(wt)' */
    if (Diesel_DER2_DW.systemEnable_a != 0) {
      Diesel_DER2_DW.lastSin_a = sin(6.2831853071795862 *
        Diesel_DER2_P.Fourier_B_Freq * Diesel_DER2_M->Timing.t[1]);
      Diesel_DER2_DW.lastCos_c = cos(6.2831853071795862 *
        Diesel_DER2_P.Fourier_B_Freq * Diesel_DER2_M->Timing.t[1]);
      Diesel_DER2_DW.systemEnable_a = 0;
    }

    rtb_Switch3_e = ((Diesel_DER2_DW.lastSin_a * Diesel_DER2_P.sinwt_PCos_d +
                      Diesel_DER2_DW.lastCos_c * Diesel_DER2_P.sinwt_PSin_l) *
                     Diesel_DER2_P.sinwt_HCos_p + (Diesel_DER2_DW.lastCos_c *
      Diesel_DER2_P.sinwt_PCos_d - Diesel_DER2_DW.lastSin_a *
      Diesel_DER2_P.sinwt_PSin_l) * Diesel_DER2_P.sinwt_Hsin_j) *
      Diesel_DER2_P.sinwt_Amp_l + Diesel_DER2_P.sinwt_Bias_i;

    /* End of Sin: '<S88>/sin(wt)' */

    /* Product: '<S88>/Product' */
    rtb_Product_o = rtb_Mult1_idx_0 * rtb_Switch3_e;

    /* DiscreteIntegrator: '<S103>/Integ4' */
    if (Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_j != 0) {
      Diesel_DER2_B.Integ4_m = Diesel_DER2_DW.Integ4_DSTATE_b;
    } else {
      Diesel_DER2_B.Integ4_m = Diesel_DER2_P.Integ4_gainval_d * rtb_Product_o +
        Diesel_DER2_DW.Integ4_DSTATE_b;
    }

    /* End of DiscreteIntegrator: '<S103>/Integ4' */

    /* Constant: '<S103>/K1' */
    Diesel_DER2_B.K1_c = Diesel_DER2_P.K1_Value_i;

    /* Level2 S-Function Block: '<S104>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[2];
      sfcnOutputs(rts, 1);
    }

    /* DigitalClock: '<S103>/Digital  Clock' */
    rtb_Switch3_e = Diesel_DER2_M->Timing.t[1];

    /* Switch: '<S103>/Switch' incorporates:
     *  Constant: '<S103>/K2'
     *  Gain: '<S103>/Gain'
     *  Gain: '<S103>/Gain1'
     *  Product: '<S103>/Product'
     *  RelationalOperator: '<S103>/Relational Operator'
     *  Sum: '<S103>/Sum1'
     *  Sum: '<S103>/Sum5'
     *  Sum: '<S103>/Sum7'
     *  UnitDelay: '<S103>/Unit Delay'
     *  UnitDelay: '<S103>/Unit Delay1'
     */
    if (rtb_Switch3_e >= Diesel_DER2_B.K1_c) {
      rtb_Switch_i = (Diesel_DER2_B.Integ4_m - Diesel_DER2_B.SFunction_e) *
        Diesel_DER2_P.K2_Value_m + (Diesel_DER2_P.Gain1_Gain_a * rtb_Product_o -
        Diesel_DER2_P.Gain_Gain_ne * Diesel_DER2_DW.UnitDelay_DSTATE_j);
    } else {
      rtb_Switch_i = Diesel_DER2_DW.UnitDelay1_DSTATE_e;
    }

    /* End of Switch: '<S103>/Switch' */

    /* Sin: '<S88>/cos(wt)' */
    if (Diesel_DER2_DW.systemEnable_o != 0) {
      Diesel_DER2_DW.lastSin_f = sin(6.2831853071795862 *
        Diesel_DER2_P.Fourier_B_Freq * Diesel_DER2_M->Timing.t[1]);
      Diesel_DER2_DW.lastCos_eg = cos(6.2831853071795862 *
        Diesel_DER2_P.Fourier_B_Freq * Diesel_DER2_M->Timing.t[1]);
      Diesel_DER2_DW.systemEnable_o = 0;
    }

    rtb_Switch3_e = ((Diesel_DER2_DW.lastSin_f * Diesel_DER2_P.coswt_PCos_n +
                      Diesel_DER2_DW.lastCos_eg * Diesel_DER2_P.coswt_PSin_k) *
                     Diesel_DER2_P.coswt_HCos_m + (Diesel_DER2_DW.lastCos_eg *
      Diesel_DER2_P.coswt_PCos_n - Diesel_DER2_DW.lastSin_f *
      Diesel_DER2_P.coswt_PSin_k) * Diesel_DER2_P.coswt_Hsin_p) *
      Diesel_DER2_P.coswt_Amp_h + Diesel_DER2_P.coswt_Bias_h;

    /* End of Sin: '<S88>/cos(wt)' */

    /* Product: '<S88>/Product1' */
    rtb_Product1_a = rtb_Mult1_idx_0 * rtb_Switch3_e;

    /* DiscreteIntegrator: '<S101>/Integ4' */
    if (Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_fj != 0) {
      Diesel_DER2_B.Integ4_o = Diesel_DER2_DW.Integ4_DSTATE_f;
    } else {
      Diesel_DER2_B.Integ4_o = Diesel_DER2_P.Integ4_gainval_e * rtb_Product1_a +
        Diesel_DER2_DW.Integ4_DSTATE_f;
    }

    /* End of DiscreteIntegrator: '<S101>/Integ4' */

    /* Constant: '<S101>/K1' */
    Diesel_DER2_B.K1_j = Diesel_DER2_P.K1_Value_b;

    /* Level2 S-Function Block: '<S102>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[3];
      sfcnOutputs(rts, 1);
    }

    /* DigitalClock: '<S101>/Digital  Clock' */
    rtb_Switch3_e = Diesel_DER2_M->Timing.t[1];

    /* Switch: '<S101>/Switch' incorporates:
     *  Constant: '<S101>/K2'
     *  Gain: '<S101>/Gain'
     *  Gain: '<S101>/Gain1'
     *  Product: '<S101>/Product'
     *  RelationalOperator: '<S101>/Relational Operator'
     *  Sum: '<S101>/Sum1'
     *  Sum: '<S101>/Sum5'
     *  Sum: '<S101>/Sum7'
     *  UnitDelay: '<S101>/Unit Delay'
     *  UnitDelay: '<S101>/Unit Delay1'
     */
    if (rtb_Switch3_e >= Diesel_DER2_B.K1_j) {
      rtb_Switch_o = (Diesel_DER2_B.Integ4_o - Diesel_DER2_B.SFunction_i) *
        Diesel_DER2_P.K2_Value_o + (Diesel_DER2_P.Gain1_Gain_i * rtb_Product1_a
        - Diesel_DER2_P.Gain_Gain_p * Diesel_DER2_DW.UnitDelay_DSTATE_l);
    } else {
      rtb_Switch_o = Diesel_DER2_DW.UnitDelay1_DSTATE_ep;
    }

    /* End of Switch: '<S101>/Switch' */

    /* ComplexToMagnitudeAngle: '<S88>/Complex to Magnitude-Angle' incorporates:
     *  RealImagToComplex: '<S88>/Real-Imag to Complex'
     */
    rtb_Switch3_e = rt_hypotd_snf(rtb_Switch_i, rtb_Switch_o);
    rtb_donotdeletethisgain_k = rt_atan2d_snf(rtb_Switch_o, rtb_Switch_i);

    /* Gain: '<S86>/deg->rad1' incorporates:
     *  Gain: '<S88>/Rad->Deg.'
     */
    rtb_donotdeletethisgain_k = Diesel_DER2_P.RadDeg_Gain_e *
      rtb_donotdeletethisgain_k * Diesel_DER2_P.degrad1_Gain;

    /* MagnitudeAngleToComplex: '<S86>/Magnitude-Angle to Complex1' */
    rtb_MagnitudeAngletoComplex1.re = rtb_Switch3_e * cos
      (rtb_donotdeletethisgain_k);
    rtb_MagnitudeAngletoComplex1.im = rtb_Switch3_e * sin
      (rtb_donotdeletethisgain_k);

    /* Sin: '<S89>/sin(wt)' */
    if (Diesel_DER2_DW.systemEnable_j != 0) {
      Diesel_DER2_DW.lastSin_o = sin(6.2831853071795862 *
        Diesel_DER2_P.Fourier_C_Freq * Diesel_DER2_M->Timing.t[1]);
      Diesel_DER2_DW.lastCos_h = cos(6.2831853071795862 *
        Diesel_DER2_P.Fourier_C_Freq * Diesel_DER2_M->Timing.t[1]);
      Diesel_DER2_DW.systemEnable_j = 0;
    }

    rtb_Switch3_e = ((Diesel_DER2_DW.lastSin_o * Diesel_DER2_P.sinwt_PCos_c +
                      Diesel_DER2_DW.lastCos_h * Diesel_DER2_P.sinwt_PSin_k) *
                     Diesel_DER2_P.sinwt_HCos_k + (Diesel_DER2_DW.lastCos_h *
      Diesel_DER2_P.sinwt_PCos_c - Diesel_DER2_DW.lastSin_o *
      Diesel_DER2_P.sinwt_PSin_k) * Diesel_DER2_P.sinwt_Hsin_d) *
      Diesel_DER2_P.sinwt_Amp_e + Diesel_DER2_P.sinwt_Bias_n;

    /* End of Sin: '<S89>/sin(wt)' */

    /* Product: '<S89>/Product' */
    rtb_Product_h = rtb_Mult1_idx_1 * rtb_Switch3_e;

    /* DiscreteIntegrator: '<S109>/Integ4' */
    if (Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_h != 0) {
      Diesel_DER2_B.Integ4_i = Diesel_DER2_DW.Integ4_DSTATE_o;
    } else {
      Diesel_DER2_B.Integ4_i = Diesel_DER2_P.Integ4_gainval_l * rtb_Product_h +
        Diesel_DER2_DW.Integ4_DSTATE_o;
    }

    /* End of DiscreteIntegrator: '<S109>/Integ4' */

    /* Constant: '<S109>/K1' */
    Diesel_DER2_B.K1_k = Diesel_DER2_P.K1_Value_g;

    /* Level2 S-Function Block: '<S110>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[4];
      sfcnOutputs(rts, 1);
    }

    /* DigitalClock: '<S109>/Digital  Clock' */
    rtb_Switch3_e = Diesel_DER2_M->Timing.t[1];

    /* Switch: '<S109>/Switch' incorporates:
     *  Constant: '<S109>/K2'
     *  Gain: '<S109>/Gain'
     *  Gain: '<S109>/Gain1'
     *  Product: '<S109>/Product'
     *  RelationalOperator: '<S109>/Relational Operator'
     *  Sum: '<S109>/Sum1'
     *  Sum: '<S109>/Sum5'
     *  Sum: '<S109>/Sum7'
     *  UnitDelay: '<S109>/Unit Delay'
     *  UnitDelay: '<S109>/Unit Delay1'
     */
    if (rtb_Switch3_e >= Diesel_DER2_B.K1_k) {
      rtb_Switch_h = (Diesel_DER2_B.Integ4_i - Diesel_DER2_B.SFunction_g) *
        Diesel_DER2_P.K2_Value_mj + (Diesel_DER2_P.Gain1_Gain_l * rtb_Product_h
        - Diesel_DER2_P.Gain_Gain_b * Diesel_DER2_DW.UnitDelay_DSTATE_h);
    } else {
      rtb_Switch_h = Diesel_DER2_DW.UnitDelay1_DSTATE_k;
    }

    /* End of Switch: '<S109>/Switch' */

    /* Sin: '<S89>/cos(wt)' */
    if (Diesel_DER2_DW.systemEnable_p != 0) {
      Diesel_DER2_DW.lastSin_i = sin(6.2831853071795862 *
        Diesel_DER2_P.Fourier_C_Freq * Diesel_DER2_M->Timing.t[1]);
      Diesel_DER2_DW.lastCos_g = cos(6.2831853071795862 *
        Diesel_DER2_P.Fourier_C_Freq * Diesel_DER2_M->Timing.t[1]);
      Diesel_DER2_DW.systemEnable_p = 0;
    }

    rtb_Switch3_e = ((Diesel_DER2_DW.lastSin_i * Diesel_DER2_P.coswt_PCos_d +
                      Diesel_DER2_DW.lastCos_g * Diesel_DER2_P.coswt_PSin_j) *
                     Diesel_DER2_P.coswt_HCos_b + (Diesel_DER2_DW.lastCos_g *
      Diesel_DER2_P.coswt_PCos_d - Diesel_DER2_DW.lastSin_i *
      Diesel_DER2_P.coswt_PSin_j) * Diesel_DER2_P.coswt_Hsin_ph) *
      Diesel_DER2_P.coswt_Amp_a + Diesel_DER2_P.coswt_Bias_a;

    /* End of Sin: '<S89>/cos(wt)' */

    /* Product: '<S89>/Product1' */
    rtb_Product1_e = rtb_Mult1_idx_1 * rtb_Switch3_e;

    /* DiscreteIntegrator: '<S107>/Integ4' */
    if (Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_c != 0) {
      Diesel_DER2_B.Integ4_k = Diesel_DER2_DW.Integ4_DSTATE_l;
    } else {
      Diesel_DER2_B.Integ4_k = Diesel_DER2_P.Integ4_gainval_a * rtb_Product1_e +
        Diesel_DER2_DW.Integ4_DSTATE_l;
    }

    /* End of DiscreteIntegrator: '<S107>/Integ4' */

    /* Constant: '<S107>/K1' */
    Diesel_DER2_B.K1_o = Diesel_DER2_P.K1_Value_l;

    /* Level2 S-Function Block: '<S108>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[5];
      sfcnOutputs(rts, 1);
    }

    /* DigitalClock: '<S107>/Digital  Clock' */
    rtb_Switch3_e = Diesel_DER2_M->Timing.t[1];

    /* Switch: '<S107>/Switch' incorporates:
     *  Constant: '<S107>/K2'
     *  Gain: '<S107>/Gain'
     *  Gain: '<S107>/Gain1'
     *  Product: '<S107>/Product'
     *  RelationalOperator: '<S107>/Relational Operator'
     *  Sum: '<S107>/Sum1'
     *  Sum: '<S107>/Sum5'
     *  Sum: '<S107>/Sum7'
     *  UnitDelay: '<S107>/Unit Delay'
     *  UnitDelay: '<S107>/Unit Delay1'
     */
    if (rtb_Switch3_e >= Diesel_DER2_B.K1_o) {
      rtb_Switch_of = (Diesel_DER2_B.Integ4_k - Diesel_DER2_B.SFunction_p) *
        Diesel_DER2_P.K2_Value_g + (Diesel_DER2_P.Gain1_Gain_md * rtb_Product1_e
        - Diesel_DER2_P.Gain_Gain_m * Diesel_DER2_DW.UnitDelay_DSTATE_ch);
    } else {
      rtb_Switch_of = Diesel_DER2_DW.UnitDelay1_DSTATE_j;
    }

    /* End of Switch: '<S107>/Switch' */

    /* ComplexToMagnitudeAngle: '<S89>/Complex to Magnitude-Angle' incorporates:
     *  RealImagToComplex: '<S89>/Real-Imag to Complex'
     */
    rtb_Switch3_e = rt_hypotd_snf(rtb_Switch_h, rtb_Switch_of);
    rtb_donotdeletethisgain_k = rt_atan2d_snf(rtb_Switch_of, rtb_Switch_h);

    /* Gain: '<S86>/deg->rad2' incorporates:
     *  Gain: '<S89>/Rad->Deg.'
     */
    rtb_donotdeletethisgain_k = Diesel_DER2_P.RadDeg_Gain_b *
      rtb_donotdeletethisgain_k * Diesel_DER2_P.degrad2_Gain;

    /* MagnitudeAngleToComplex: '<S86>/Magnitude-Angle to Complex2' */
    rtb_MagnitudeAngletoComplex2.re = rtb_Switch3_e * cos
      (rtb_donotdeletethisgain_k);
    rtb_MagnitudeAngletoComplex2.im = rtb_Switch3_e * sin
      (rtb_donotdeletethisgain_k);

    /* Outputs for Enabled SubSystem: '<S86>/Pos. Seq. Computation' */

    /* Constant: '<S86>/Constant' */
    Diesel_DER2_NegSeqComputation(Diesel_DER2_P.Constant_Value,
      rtb_MagnitudeAngletoComplex, rtb_MagnitudeAngletoComplex1,
      rtb_MagnitudeAngletoComplex2, &Diesel_DER2_B.PosSeqComputation,
      (P_NegSeqComputation_Diesel_DE_T *)&Diesel_DER2_P.PosSeqComputation);

    /* End of Outputs for SubSystem: '<S86>/Pos. Seq. Computation' */

    /* Outputs for Enabled SubSystem: '<S86>/Neg. Seq. Computation' */

    /* Constant: '<S86>/Constant1' */
    Diesel_DER2_NegSeqComputation(Diesel_DER2_P.Constant1_Value_j,
      rtb_MagnitudeAngletoComplex, rtb_MagnitudeAngletoComplex1,
      rtb_MagnitudeAngletoComplex2, &Diesel_DER2_B.NegSeqComputation,
      (P_NegSeqComputation_Diesel_DE_T *)&Diesel_DER2_P.NegSeqComputation);

    /* End of Outputs for SubSystem: '<S86>/Neg. Seq. Computation' */

    /* Outputs for Enabled SubSystem: '<S86>/Zero Seq. Computation' incorporates:
     *  EnablePort: '<S92>/Enable'
     */
    /* Constant: '<S86>/Constant2' */
    if (Diesel_DER2_P.Constant2_Value_m > 0.0) {
      /* Gain: '<S92>/Gain3' incorporates:
       *  Sum: '<S92>/Sum'
       */
      Diesel_DER2_B.Gain3.re = ((rtb_MagnitudeAngletoComplex.re +
        rtb_MagnitudeAngletoComplex1.re) + rtb_MagnitudeAngletoComplex2.re) *
        Diesel_DER2_P.Gain3_Gain;
      Diesel_DER2_B.Gain3.im = ((rtb_MagnitudeAngletoComplex.im +
        rtb_MagnitudeAngletoComplex1.im) + rtb_MagnitudeAngletoComplex2.im) *
        Diesel_DER2_P.Gain3_Gain;
    }

    /* End of Constant: '<S86>/Constant2' */
    /* End of Outputs for SubSystem: '<S86>/Zero Seq. Computation' */

    /* ComplexToMagnitudeAngle: '<S86>/Complex to Magnitude-Angle' */
    rtb_ComplextoMagnitudeAngle_o1[0] = rt_hypotd_snf
      (Diesel_DER2_B.PosSeqComputation.Gain3.re,
       Diesel_DER2_B.PosSeqComputation.Gain3.im);
    rtb_ComplextoMagnitudeAngle_o1[1] = rt_hypotd_snf
      (Diesel_DER2_B.NegSeqComputation.Gain3.re,
       Diesel_DER2_B.NegSeqComputation.Gain3.im);
    rtb_ComplextoMagnitudeAngle_o1[2] = rt_hypotd_snf(Diesel_DER2_B.Gain3.re,
      Diesel_DER2_B.Gain3.im);

    /* Gain: '<S40>/do not delete this gain' */
    rtb_Switch3_e = Diesel_DER2_P.donotdeletethisgain_Gain_p *
      Diesel_DER2_B.StateSpace_o1[23];

    /* Gain: '<S41>/do not delete this gain' */
    rtb_donotdeletethisgain_k = Diesel_DER2_P.donotdeletethisgain_Gain_b *
      Diesel_DER2_B.StateSpace_o1[24];

    /* Gain: '<S42>/do not delete this gain' */
    rtb_donotdeletethisgain_i = Diesel_DER2_P.donotdeletethisgain_Gain_i *
      Diesel_DER2_B.StateSpace_o1[25];

    /* Outport: '<Root>/out_Iabc_Exc' incorporates:
     *  Gain: '<S31>/Kv'
     */
    Diesel_DER2_Y.out_Iabc_Exc[0] = Diesel_DER2_P.Kv_Gain * rtb_Switch3_e;
    Diesel_DER2_Y.out_Iabc_Exc[1] = Diesel_DER2_P.Kv_Gain *
      rtb_donotdeletethisgain_k;
    Diesel_DER2_Y.out_Iabc_Exc[2] = Diesel_DER2_P.Kv_Gain *
      rtb_donotdeletethisgain_i;

    /* Gain: '<S43>/do not delete this gain' */
    rtb_Switch3_e = Diesel_DER2_P.donotdeletethisgain_Gain_lw *
      Diesel_DER2_B.StateSpace_o1[17];

    /* Gain: '<S44>/do not delete this gain' */
    rtb_donotdeletethisgain_k = Diesel_DER2_P.donotdeletethisgain_Gain_a *
      Diesel_DER2_B.StateSpace_o1[18];

    /* Gain: '<S45>/do not delete this gain' */
    rtb_donotdeletethisgain_i = Diesel_DER2_P.donotdeletethisgain_Gain_f *
      Diesel_DER2_B.StateSpace_o1[19];

    /* Outport: '<Root>/out_Vabc_Exc' incorporates:
     *  Gain: '<S31>/Kv1'
     */
    Diesel_DER2_Y.out_Vabc_Exc[0] = Diesel_DER2_P.Kv1_Gain_m * rtb_Switch3_e;
    Diesel_DER2_Y.out_Vabc_Exc[1] = Diesel_DER2_P.Kv1_Gain_m *
      rtb_donotdeletethisgain_k;
    Diesel_DER2_Y.out_Vabc_Exc[2] = Diesel_DER2_P.Kv1_Gain_m *
      rtb_donotdeletethisgain_i;

    /* Gain: '<S117>/N' */
    rtb_Vfd_o = Diesel_DER2_P.N_Gain_j * rtb_UnitDelay4;

    /* DigitalClock: '<S132>/Digital Clock' */
    rtb_Switch3_e = Diesel_DER2_M->Timing.t[1];

    /* Switch: '<S132>/IC' */
    if (rtb_Switch3_e >= Diesel_DER2_P.IC_Threshold_e) {
      /* Assignment: '<S129>/W(1,2)=wr' incorporates:
       *  Constant: '<S129>/u1'
       */
      memcpy(&rtb_Sum2_m3[0], &Diesel_DER2_P.u1_Value_l[0], 25U * sizeof(real_T));
      rtb_Sum2_m3[5] = rtb_Phisat;

      /* Assignment: '<S129>/W(2,1)=-wr' incorporates:
       *  Gain: '<S129>/Gain1'
       */
      rtb_Sum2_m3[1] = Diesel_DER2_P.Gain1_Gain_l3 * rtb_Phisat;

      /* Switch: '<S121>/Switch1' incorporates:
       *  Constant: '<S121>/Constant2'
       */
      tmp = (Diesel_DER2_P.Constant2_Value_e >=
             Diesel_DER2_P.Switch1_Threshold_a);
      for (i = 0; i < 25; i++) {
        /* Gain: '<S140>/wbase*Ts//2' incorporates:
         *  Constant: '<S121>/Constant6'
         *  Sum: '<S121>/Sum1'
         *  Switch: '<S121>/Switch1'
         */
        if (tmp) {
          rtb_ElementaryMath1 = Diesel_DER2_B.RLinv[i];
        } else {
          rtb_ElementaryMath1 = Diesel_DER2_P.Constant6_Value_o[i];
        }

        rtb_ElementaryMath1 = ((0.0 - rtb_Sum2_m3[i]) - rtb_ElementaryMath1) *
          Diesel_DER2_P.wbaseTs2_Gain_p;

        /* Sum: '<S140>/Sum1' incorporates:
         *  Constant: '<S140>/u5'
         */
        rtb_Sum2_f[i] = Diesel_DER2_P.u5_Value_e[i] - rtb_ElementaryMath1;

        /* Sum: '<S140>/Sum5' incorporates:
         *  Constant: '<S140>/u5'
         */
        rtb_ElementaryMath1 += Diesel_DER2_P.u5_Value_e[i];
        rtb_Sum2_m3[i] = rtb_ElementaryMath1;
      }

      /* Product: '<S140>/inversion' incorporates:
       *  Gain: '<S140>/wbase*Ts//2'
       *  Sum: '<S121>/Sum1'
       */
      rt_invd5x5_snf(rtb_Sum2_f, rtb_inversion_g);

      /* Product: '<S140>/Product4' incorporates:
       *  Product: '<S132>/Product2'
       */
      for (i = 0; i < 5; i++) {
        for (i_1 = 0; i_1 < 5; i_1++) {
          rtb_Sum2_f[i + 5 * i_1] = 0.0;
          for (i_0 = 0; i_0 < 5; i_0++) {
            rtb_Sum2_f[i + 5 * i_1] += rtb_inversion_g[5 * i_0 + i] *
              rtb_Sum2_m3[5 * i_1 + i_0];
          }
        }
      }

      /* End of Product: '<S140>/Product4' */

      /* Sum: '<S132>/sum' incorporates:
       *  Constant: '<S117>/[ Vkd =0 Vkq1=0  Vkq2=0 ]'
       */
      tmp_2[0] = rtb_Fcn2;
      tmp_2[1] = rtb_Fcn3;
      tmp_2[2] = rtb_Vfd_o;
      tmp_2[3] = Diesel_DER2_P.Vkd0Vkq10Vkq20_Value_a[0];
      tmp_2[4] = Diesel_DER2_P.Vkd0Vkq10Vkq20_Value_a[1];
      for (i = 0; i < 5; i++) {
        /* Sum: '<S132>/sum' incorporates:
         *  Product: '<S132>/Product1'
         *  UnitDelay: '<S132>/voltages'
         */
        tmp_0[i] = tmp_2[i] + Diesel_DER2_DW.voltages_DSTATE_h[i];

        /* Product: '<S132>/Product2' incorporates:
         *  Sum: '<S132>/Ad*x(k-1) + Bd*( u(k-1) + u(k))'
         */
        rtb_changeIqIdcurrentsigns_b[i] = 0.0;
        for (i_1 = 0; i_1 < 5; i_1++) {
          rtb_changeIqIdcurrentsigns_b[i] += rtb_Sum2_f[5 * i_1 + i] *
            rtb_xk1[i_1];
        }
      }

      /* Product: '<S132>/Product1' incorporates:
       *  Gain: '<S140>/wbase*Ts//2 '
       *  Sum: '<S132>/Ad*x(k-1) + Bd*( u(k-1) + u(k))'
       */
      for (i = 0; i < 5; i++) {
        tmp_1[i] = 0.0;
        for (i_1 = 0; i_1 < 5; i_1++) {
          tmp_1[i] += rtb_inversion_g[5 * i_1 + i] *
            Diesel_DER2_P.wbaseTs2_Gain_j * tmp_0[i_1];
        }

        rtb_IC_o[i] = rtb_changeIqIdcurrentsigns_b[i] + tmp_1[i];
      }
    } else {
      for (i = 0; i < 5; i++) {
        rtb_IC_o[i] = rtb_xk1[i];
      }
    }

    /* End of Switch: '<S132>/IC' */

    /* Gain: '<S119>/1_Pb' incorporates:
     *  Sum: '<Root>/Sum1'
     */
    rtb_Switch3_e = (rtb_pu2Watts - rtb_Switch_g) * Diesel_DER2_P._Pb_Gain_f;

    /* Fcn: '<S119>/div' */
    rtb_donotdeletethisgain_k = rtb_Switch3_e / rtb_Phisat;

    /* Gain: '<S119>/1 ----- 2H' incorporates:
     *  Gain: '<S119>/F'
     *  Gain: '<S128>/1-1'
     *  Product: '<S128>/Mult1'
     *  Sum: '<S119>/Sum2'
     *  Sum: '<S128>/Sum2'
     */
    rtb_uH = ((rtb_donotdeletethisgain_k - (Diesel_DER2_P.u_Gain_k[0] *
                rtb_changeIqIdcurrentsigns[0] * rtb_xk1[1] +
                Diesel_DER2_P.u_Gain_k[1] * rtb_changeIqIdcurrentsigns[1] *
                rtb_xk1[0])) - Diesel_DER2_P.F_Gain * rtb_Phisat) *
      Diesel_DER2_P.uH_Gain;

    /* DiscreteIntegrator: '<S119>/Rotor speed deviation (dw)' */
    if (Diesel_DER2_DW.Rotorspeeddeviationdw_SYSTEM_EN != 0) {
      Diesel_DER2_B.dw = Diesel_DER2_DW.Rotorspeeddeviationdw_DSTATE;
    } else {
      Diesel_DER2_B.dw = Diesel_DER2_P.Rotorspeeddeviationdw_gainval * rtb_uH +
        Diesel_DER2_DW.Rotorspeeddeviationdw_DSTATE;
    }

    /* End of DiscreteIntegrator: '<S119>/Rotor speed deviation (dw)' */

    /* Gain: '<S119>/we base' */
    rtb_webase = Diesel_DER2_P.webase_Gain * rtb_DataTypeConversion;

    /* DiscreteIntegrator: '<S111>/Integ4' */
    if (Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_d != 0) {
      Diesel_DER2_B.Integ4_a = Diesel_DER2_DW.Integ4_DSTATE_m;
    } else {
      Diesel_DER2_B.Integ4_a = Diesel_DER2_P.Integ4_gainval_j * rtb_UnitDelay4 +
        Diesel_DER2_DW.Integ4_DSTATE_m;
    }

    /* End of DiscreteIntegrator: '<S111>/Integ4' */

    /* Constant: '<S111>/K1' */
    Diesel_DER2_B.K1_m = Diesel_DER2_P.K1_Value_fz;

    /* Level2 S-Function Block: '<S112>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[6];
      sfcnOutputs(rts, 1);
    }

    /* DigitalClock: '<S111>/Digital  Clock' */
    rtb_DataTypeConversion = Diesel_DER2_M->Timing.t[1];

    /* Switch: '<S111>/Switch' incorporates:
     *  Constant: '<S111>/K2'
     *  Gain: '<S111>/Gain'
     *  Gain: '<S111>/Gain1'
     *  Product: '<S111>/Product'
     *  RelationalOperator: '<S111>/Relational Operator'
     *  Sum: '<S111>/Sum1'
     *  Sum: '<S111>/Sum5'
     *  Sum: '<S111>/Sum7'
     *  UnitDelay: '<S111>/Unit Delay'
     *  UnitDelay: '<S111>/Unit Delay1'
     */
    if (rtb_DataTypeConversion >= Diesel_DER2_B.K1_m) {
      rtb_Switch_a = (Diesel_DER2_B.Integ4_a - Diesel_DER2_B.SFunction_iq) *
        Diesel_DER2_P.K2_Value_h + (Diesel_DER2_P.Gain1_Gain_g * rtb_UnitDelay4
        - Diesel_DER2_P.Gain_Gain_mq * Diesel_DER2_DW.UnitDelay_DSTATE_a);
    } else {
      rtb_Switch_a = Diesel_DER2_DW.UnitDelay1_DSTATE_m;
    }

    /* End of Switch: '<S111>/Switch' */

    /* Outport: '<Root>/out_Vf' */
    Diesel_DER2_Y.out_Vf[0] = rtb_UnitDelay4;
    Diesel_DER2_Y.out_Vf[1] = rtb_Switch_a;

    /* Outport: '<Root>/out_Vt_SM' incorporates:
     *  UnitDelay: '<S36>/Unit Delay1'
     */
    Diesel_DER2_Y.out_Vt_SM = Diesel_DER2_DW.UnitDelay1_DSTATE;

    /* Outport: '<Root>/out_Vabc_SM' */
    Diesel_DER2_Y.out_Vabc_SM[0] = rtb_Switch1_h;
    Diesel_DER2_Y.out_Vabc_SM[1] = rtb_Mult1_idx_0;
    Diesel_DER2_Y.out_Vabc_SM[2] = rtb_Mult1_idx_1;

    /* Gain: '<S13>/do not delete this gain' */
    rtb_DataTypeConversion = Diesel_DER2_P.donotdeletethisgain_Gain_lw2 *
      Diesel_DER2_B.StateSpace_o1[20];

    /* Gain: '<S14>/do not delete this gain' */
    rtb_Switch3_e = Diesel_DER2_P.donotdeletethisgain_Gain_lq *
      Diesel_DER2_B.StateSpace_o1[21];

    /* Gain: '<S15>/do not delete this gain' */
    rtb_donotdeletethisgain_k = Diesel_DER2_P.donotdeletethisgain_Gain_ad *
      Diesel_DER2_B.StateSpace_o1[22];

    /* Outport: '<Root>/out_Iabc_SM' incorporates:
     *  Gain: '<S1>/Kv'
     */
    Diesel_DER2_Y.out_Iabc_SM[0] = Diesel_DER2_P.Kv_Gain_d *
      rtb_DataTypeConversion;
    Diesel_DER2_Y.out_Iabc_SM[1] = Diesel_DER2_P.Kv_Gain_d * rtb_Switch3_e;
    Diesel_DER2_Y.out_Iabc_SM[2] = Diesel_DER2_P.Kv_Gain_d *
      rtb_donotdeletethisgain_k;

    /* Outport: '<Root>/out_VD1' incorporates:
     *  Gain: '<S5>/Gain2'
     */
    Diesel_DER2_Y.out_VD1 = Diesel_DER2_P.Gain2_Gain *
      Diesel_DER2_B.StateSpace_o1[26];

    /* Outport: '<Root>/out_I' incorporates:
     *  Gain: '<S6>/Gain2'
     */
    Diesel_DER2_Y.out_I = Diesel_DER2_P.Gain2_Gain_g *
      Diesel_DER2_B.StateSpace_o1[26];

    /* Sum: '<S146>/Sum' incorporates:
     *  Gain: '<S146>/Gain'
     *  UnitDelay: '<S146>/Delay Ts'
     */
    rtb_Sum_k = Diesel_DER2_P.Gain_Gain_f * Diesel_DER2_B.StateSpace_o1[0] -
      Diesel_DER2_DW.DelayTs_DSTATE;

    /* DataTypeConversion: '<S146>/Data Type  Conversion2' incorporates:
     *  Constant: '<S146>/0 1'
     *  Product: '<S146>/Product'
     *  RelationalOperator: '<S146>/b5'
     *  UnitDelay: '<S146>/Delay Ts '
     */
    rtb_DataTypeConversion = (rtb_Sum_k * Diesel_DER2_DW.DelayTs_DSTATE_f >
      Diesel_DER2_P.u_Value);

    /* DigitalClock: '<S155>/Digital Clock' */
    rtb_Switch3_e = Diesel_DER2_M->Timing.t[1];

    /* Switch: '<S8>/Switch3' incorporates:
     *  Constant: '<S8>/C4'
     *  Constant: '<S8>/com'
     *  Lookup: '<S155>/Look-Up Table'
     */
    if (Diesel_DER2_P.C4_Value >= Diesel_DER2_P.Switch3_Threshold_f) {
      rtb_Switch_g = Diesel_DER2_P.com_Value;
    } else {
      rtb_Switch_g = rt_Lookup(Diesel_DER2_P.LookUpTable_XData, 4, rtb_Switch3_e,
        Diesel_DER2_P.LookUpTable_YData);
    }

    /* End of Switch: '<S8>/Switch3' */

    /* DigitalClock: '<S148>/Digital Clock' */
    rtb_Switch3_e = Diesel_DER2_M->Timing.t[1];

    /* Switch: '<S146>/Switch3' incorporates:
     *  Constant: '<S146>/C4'
     */
    if (Diesel_DER2_P.C4_Value_n >= Diesel_DER2_P.Switch3_Threshold_p) {
      /* Switch: '<S8>/Switch' incorporates:
       *  Constant: '<S8>/Constant1'
       *  Constant: '<S8>/Constant5'
       */
      if (Diesel_DER2_P.ThreePhaseBreaker_SwitchA >=
          Diesel_DER2_P.Switch_Threshold) {
        rtb_Switch3_e = rtb_Switch_g;
      } else {
        rtb_Switch3_e = Diesel_DER2_P.Constant5_Value_p;
      }
    } else {
      /* Switch: '<S8>/Switch' incorporates:
       *  Lookup: '<S148>/Look-Up Table'
       */
      rtb_Switch3_e = rt_Lookup(Diesel_DER2_P.LookUpTable_XData_d, 4,
        rtb_Switch3_e, Diesel_DER2_P.LookUpTable_YData_o);
    }

    /* End of Switch: '<S146>/Switch3' */

    /* Sum: '<S146>/sum' */
    rtb_DataTypeConversion += rtb_Switch3_e;

    /* Relay: '<S146>/>1.5' */
    if (rtb_DataTypeConversion >= Diesel_DER2_P.u5_OnVal) {
      Diesel_DER2_DW.u5_Mode = true;
    } else {
      if (rtb_DataTypeConversion <= Diesel_DER2_P.u5_OffVal) {
        Diesel_DER2_DW.u5_Mode = false;
      }
    }

    if (Diesel_DER2_DW.u5_Mode) {
      rtb_DataTypeConversion = Diesel_DER2_P.u5_YOn;
    } else {
      rtb_DataTypeConversion = Diesel_DER2_P.u5_YOff;
    }

    /* End of Relay: '<S146>/>1.5' */

    /* DataTypeConversion: '<S146>/Data Type  Conversion' incorporates:
     *  DataTypeConversion: '<S146>/Data Type  Conversion1'
     *  DataTypeConversion: '<S146>/Data Type  Conversion3'
     *  Logic: '<S146>/or'
     */
    rtb_DataTypeConversion = ((rtb_DataTypeConversion != 0.0) || (rtb_Switch3_e
      != 0.0));

    /* InitialCondition: '<S146>/IC' */
    if (Diesel_DER2_DW.IC_FirstOutputTime) {
      Diesel_DER2_DW.IC_FirstOutputTime = false;
      Diesel_DER2_B.IC = Diesel_DER2_P.IC_Value;
    } else {
      Diesel_DER2_B.IC = rtb_DataTypeConversion;
    }

    /* End of InitialCondition: '<S146>/IC' */

    /* Sum: '<S149>/Sum' incorporates:
     *  Gain: '<S149>/Gain'
     *  UnitDelay: '<S149>/Delay Ts'
     */
    rtb_Sum_o = Diesel_DER2_P.Gain_Gain_j * Diesel_DER2_B.StateSpace_o1[1] -
      Diesel_DER2_DW.DelayTs_DSTATE_a;

    /* DataTypeConversion: '<S149>/Data Type  Conversion2' incorporates:
     *  Constant: '<S149>/0 1'
     *  Product: '<S149>/Product'
     *  RelationalOperator: '<S149>/b5'
     *  UnitDelay: '<S149>/Delay Ts '
     */
    rtb_DataTypeConversion = (rtb_Sum_o * Diesel_DER2_DW.DelayTs_DSTATE_n >
      Diesel_DER2_P.u_Value_p);

    /* DigitalClock: '<S151>/Digital Clock' */
    rtb_Switch3_e = Diesel_DER2_M->Timing.t[1];

    /* Switch: '<S149>/Switch3' incorporates:
     *  Constant: '<S149>/C4'
     */
    if (Diesel_DER2_P.C4_Value_j >= Diesel_DER2_P.Switch3_Threshold_i) {
      /* Switch: '<S8>/Switch1' incorporates:
       *  Constant: '<S8>/Constant2'
       *  Constant: '<S8>/Constant5'
       */
      if (Diesel_DER2_P.ThreePhaseBreaker_SwitchB >=
          Diesel_DER2_P.Switch1_Threshold_al) {
        rtb_Switch3_e = rtb_Switch_g;
      } else {
        rtb_Switch3_e = Diesel_DER2_P.Constant5_Value_p;
      }
    } else {
      /* Switch: '<S8>/Switch1' incorporates:
       *  Lookup: '<S151>/Look-Up Table'
       */
      rtb_Switch3_e = rt_Lookup(Diesel_DER2_P.LookUpTable_XData_a, 4,
        rtb_Switch3_e, Diesel_DER2_P.LookUpTable_YData_c);
    }

    /* End of Switch: '<S149>/Switch3' */

    /* Sum: '<S149>/sum' */
    rtb_DataTypeConversion += rtb_Switch3_e;

    /* Relay: '<S149>/>1.5' */
    if (rtb_DataTypeConversion >= Diesel_DER2_P.u5_OnVal_p) {
      Diesel_DER2_DW.u5_Mode_n = true;
    } else {
      if (rtb_DataTypeConversion <= Diesel_DER2_P.u5_OffVal_m) {
        Diesel_DER2_DW.u5_Mode_n = false;
      }
    }

    if (Diesel_DER2_DW.u5_Mode_n) {
      rtb_DataTypeConversion = Diesel_DER2_P.u5_YOn_f;
    } else {
      rtb_DataTypeConversion = Diesel_DER2_P.u5_YOff_g;
    }

    /* End of Relay: '<S149>/>1.5' */

    /* DataTypeConversion: '<S149>/Data Type  Conversion' incorporates:
     *  DataTypeConversion: '<S149>/Data Type  Conversion1'
     *  DataTypeConversion: '<S149>/Data Type  Conversion3'
     *  Logic: '<S149>/or'
     */
    rtb_DataTypeConversion = ((rtb_DataTypeConversion != 0.0) || (rtb_Switch3_e
      != 0.0));

    /* InitialCondition: '<S149>/IC' */
    if (Diesel_DER2_DW.IC_FirstOutputTime_l) {
      Diesel_DER2_DW.IC_FirstOutputTime_l = false;
      Diesel_DER2_B.IC_a = Diesel_DER2_P.IC_Value_j;
    } else {
      Diesel_DER2_B.IC_a = rtb_DataTypeConversion;
    }

    /* End of InitialCondition: '<S149>/IC' */

    /* Sum: '<S152>/Sum' incorporates:
     *  Gain: '<S152>/Gain'
     *  UnitDelay: '<S152>/Delay Ts'
     */
    rtb_Sum_b = Diesel_DER2_P.Gain_Gain_bi * Diesel_DER2_B.StateSpace_o1[2] -
      Diesel_DER2_DW.DelayTs_DSTATE_l;

    /* DataTypeConversion: '<S152>/Data Type  Conversion2' incorporates:
     *  Constant: '<S152>/0 1'
     *  Product: '<S152>/Product'
     *  RelationalOperator: '<S152>/b5'
     *  UnitDelay: '<S152>/Delay Ts '
     */
    rtb_DataTypeConversion = (rtb_Sum_b * Diesel_DER2_DW.DelayTs_DSTATE_d >
      Diesel_DER2_P.u_Value_d);

    /* DigitalClock: '<S154>/Digital Clock' */
    rtb_Switch3_e = Diesel_DER2_M->Timing.t[1];

    /* Switch: '<S152>/Switch3' incorporates:
     *  Constant: '<S152>/C4'
     */
    if (Diesel_DER2_P.C4_Value_m >= Diesel_DER2_P.Switch3_Threshold_c) {
      /* Switch: '<S8>/Switch2' incorporates:
       *  Constant: '<S8>/Constant3'
       *  Constant: '<S8>/Constant5'
       */
      if (Diesel_DER2_P.ThreePhaseBreaker_SwitchC >=
          Diesel_DER2_P.Switch2_Threshold) {
        rtb_Switch3_e = rtb_Switch_g;
      } else {
        rtb_Switch3_e = Diesel_DER2_P.Constant5_Value_p;
      }
    } else {
      /* Switch: '<S8>/Switch2' incorporates:
       *  Lookup: '<S154>/Look-Up Table'
       */
      rtb_Switch3_e = rt_Lookup(Diesel_DER2_P.LookUpTable_XData_m, 4,
        rtb_Switch3_e, Diesel_DER2_P.LookUpTable_YData_n);
    }

    /* End of Switch: '<S152>/Switch3' */

    /* Sum: '<S152>/sum' */
    rtb_DataTypeConversion += rtb_Switch3_e;

    /* Relay: '<S152>/>1.5' */
    if (rtb_DataTypeConversion >= Diesel_DER2_P.u5_OnVal_o) {
      Diesel_DER2_DW.u5_Mode_i = true;
    } else {
      if (rtb_DataTypeConversion <= Diesel_DER2_P.u5_OffVal_j) {
        Diesel_DER2_DW.u5_Mode_i = false;
      }
    }

    if (Diesel_DER2_DW.u5_Mode_i) {
      rtb_DataTypeConversion = Diesel_DER2_P.u5_YOn_g;
    } else {
      rtb_DataTypeConversion = Diesel_DER2_P.u5_YOff_b;
    }

    /* End of Relay: '<S152>/>1.5' */

    /* DataTypeConversion: '<S152>/Data Type  Conversion' incorporates:
     *  DataTypeConversion: '<S152>/Data Type  Conversion1'
     *  DataTypeConversion: '<S152>/Data Type  Conversion3'
     *  Logic: '<S152>/or'
     */
    rtb_DataTypeConversion = ((rtb_DataTypeConversion != 0.0) || (rtb_Switch3_e
      != 0.0));

    /* InitialCondition: '<S152>/IC' */
    if (Diesel_DER2_DW.IC_FirstOutputTime_n) {
      Diesel_DER2_DW.IC_FirstOutputTime_n = false;
      Diesel_DER2_B.IC_c = Diesel_DER2_P.IC_Value_n;
    } else {
      Diesel_DER2_B.IC_c = rtb_DataTypeConversion;
    }

    /* End of InitialCondition: '<S152>/IC' */
  }

  {
    real_T HoldSine;
    int32_T i;

    /* Update for DiscreteIntegrator: '<S119>/Rotor angle dthetae' */
    Diesel_DER2_DW.Rotorangledthetae_DSTATE +=
      Diesel_DER2_P.Rotorangledthetae_gainval * rtb_webase;

    /* Update for DiscreteIntegrator: '<S54>/Rotor angle dtheta' */
    Diesel_DER2_DW.Rotorangledtheta_DSTATE +=
      Diesel_DER2_P.Rotorangledtheta_gainval * rtb_web1;
    for (i = 0; i < 5; i++) {
      /* Update for UnitDelay: '<S132>/fluxes' */
      Diesel_DER2_DW.fluxes_DSTATE[i] = rtb_IC_o[i];

      /* Update for UnitDelay: '<S67>/fluxes' */
      Diesel_DER2_DW.fluxes_DSTATE_o[i] = rtb_IC[i];
    }

    /* Update for UnitDelay: '<S33>/Unit Delay2' */
    Diesel_DER2_DW.UnitDelay2_DSTATE = rtb_outputformatting[5];

    /* S-Function block: <S156>/State-Space */
    {
      const real_T *As = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.AS;
      const real_T *Bs = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.BS;
      real_T *xtmp = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.XTMP;
      real_T accum;

      /* Calculate new states... */
      {
        int_T i1;
        real_T *xd = &Diesel_DER2_DW.StateSpace_DSTATE[0];
        for (i1=0; i1 < 6; i1++) {
          accum = 0.0;

          {
            int_T i2;
            real_T *xd = &Diesel_DER2_DW.StateSpace_DSTATE[0];
            for (i2=0; i2 < 6; i2++) {
              accum += *(As++) * xd[i2];
            }
          }

          {
            int_T i2;
            const real_T *u0 = &Diesel_DER2_P.SwitchCurrents_Value[0];
            for (i2=0; i2 < 9; i2++) {
              accum += *(Bs++) * u0[i2];
            }

            accum += *(Bs++) * Diesel_DER2_B.ib[0];
            accum += *(Bs++) * Diesel_DER2_B.ib[1];
            accum += *(Bs++) * Diesel_DER2_B.ib_g[0];
            accum += *(Bs++) * Diesel_DER2_B.ib_g[1];
            accum += *(Bs++) * Diesel_DER2_B.UnitDelay2;
          }

          xtmp[i1] = accum;
        }
      }

      {
        int_T i1;
        real_T *xd = &Diesel_DER2_DW.StateSpace_DSTATE[0];
        for (i1=0; i1 < 6; i1++) {
          xd[i1] = xtmp[i1];
        }
      }

      {
        int_T *gState = (int_T*)Diesel_DER2_DW.StateSpace_PWORK.G_STATE;

        /* Store switch gates values for next step */
        *(gState++) = (int_T) Diesel_DER2_B.IC;
        *(gState++) = (int_T) Diesel_DER2_B.IC_a;
        *(gState++) = (int_T) Diesel_DER2_B.IC_c;

        {
          int_T i1;
          const real_T *u1 = &Diesel_DER2_P.g_Value[0];
          for (i1=0; i1 < 6; i1++) {
            *(gState++) = (int_T) u1[i1];
          }
        }
      }
    }

    /* Update for TransportDelay: '<S2>/Delay Td' */
    {
      real_T **uBuffer = (real_T**)&Diesel_DER2_DW.DelayTd_PWORK.TUbufferPtrs[0];
      real_T **tBuffer = (real_T**)&Diesel_DER2_DW.DelayTd_PWORK.TUbufferPtrs[1];
      real_T simTime = Diesel_DER2_M->Timing.t[0];
      Diesel_DER2_DW.DelayTd_IWORK.Head = ((Diesel_DER2_DW.DelayTd_IWORK.Head <
        (Diesel_DER2_DW.DelayTd_IWORK.CircularBufSize-1)) ?
        (Diesel_DER2_DW.DelayTd_IWORK.Head+1) : 0);
      if (Diesel_DER2_DW.DelayTd_IWORK.Head == Diesel_DER2_DW.DelayTd_IWORK.Tail)
      {
        Diesel_DER2_DW.DelayTd_IWORK.Tail = ((Diesel_DER2_DW.DelayTd_IWORK.Tail <
          (Diesel_DER2_DW.DelayTd_IWORK.CircularBufSize-1)) ?
          (Diesel_DER2_DW.DelayTd_IWORK.Tail+1) : 0);
      }

      (*tBuffer)[Diesel_DER2_DW.DelayTd_IWORK.Head] = simTime;
      (*uBuffer)[Diesel_DER2_DW.DelayTd_IWORK.Head] = Diesel_DER2_B.Gazflow;
    }

    /* Update for UnitDelay: '<S141>/dw_delay' */
    Diesel_DER2_DW.dw_delay_DSTATE = Diesel_DER2_B.dw;

    /* Update for UnitDelay: '<S141>/dw_predict' */
    Diesel_DER2_DW.dw_predict_DSTATE = rtb_dw_delay;

    /* Update for UnitDelay: '<S28>/Delay_x' */
    Diesel_DER2_DW.Delay_x_DSTATE = rtb_xk1_c;

    /* Update for UnitDelay: '<S29>/Delay_x' */
    Diesel_DER2_DW.Delay_x_DSTATE_i = rtb_xk1_h;

    /* Update for UnitDelay: '<S30>/Delay_x' */
    Diesel_DER2_DW.Delay_x_DSTATE_g = rtb_xk1_i;

    /* Update for DiscreteIntegrator: '<S2>/Discrete-Time Integrator' */
    Diesel_DER2_DW.DiscreteTimeIntegrator_DSTATE +=
      Diesel_DER2_P.DiscreteTimeIntegrator_gainval * rtb_yk_m;
    if (Diesel_DER2_DW.DiscreteTimeIntegrator_DSTATE >=
        Diesel_DER2_P.DiscreteTimeIntegrator_UpperSat) {
      Diesel_DER2_DW.DiscreteTimeIntegrator_DSTATE =
        Diesel_DER2_P.DiscreteTimeIntegrator_UpperSat;
    } else {
      if (Diesel_DER2_DW.DiscreteTimeIntegrator_DSTATE <=
          Diesel_DER2_P.DiscreteTimeIntegrator_LowerSat) {
        Diesel_DER2_DW.DiscreteTimeIntegrator_DSTATE =
          Diesel_DER2_P.DiscreteTimeIntegrator_LowerSat;
      }
    }

    /* End of Update for DiscreteIntegrator: '<S2>/Discrete-Time Integrator' */

    /* Update for UnitDelay: '<S33>/Unit Delay4' */
    Diesel_DER2_DW.UnitDelay4_DSTATE = rtb_donotdeletethisgain;

    /* Update for UnitDelay: '<S36>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE = rtb_ComplextoMagnitudeAngle_o1[0];

    /* Update for DiscreteIntegrator: '<S85>/Discrete-Time Integrator' */
    Diesel_DER2_DW.DiscreteTimeIntegrator_DSTATE_f +=
      Diesel_DER2_P.DiscreteTimeIntegrator_gainva_n * rtb_Kp5;
    if (Diesel_DER2_DW.DiscreteTimeIntegrator_DSTATE_f >=
        Diesel_DER2_P.DiscreteTimeIntegrator_UpperS_p) {
      Diesel_DER2_DW.DiscreteTimeIntegrator_DSTATE_f =
        Diesel_DER2_P.DiscreteTimeIntegrator_UpperS_p;
    } else {
      if (Diesel_DER2_DW.DiscreteTimeIntegrator_DSTATE_f <=
          Diesel_DER2_P.DiscreteTimeIntegrator_LowerS_b) {
        Diesel_DER2_DW.DiscreteTimeIntegrator_DSTATE_f =
          Diesel_DER2_P.DiscreteTimeIntegrator_LowerS_b;
      }
    }

    /* End of Update for DiscreteIntegrator: '<S85>/Discrete-Time Integrator' */

    /* Update for UnitDelay: '<S67>/voltages' incorporates:
     *  Constant: '<S52>/[ Vkd =0 Vkq1=0  Vkq2=0 ]'
     */
    Diesel_DER2_DW.voltages_DSTATE[0] = rtb_Fcn2_j;
    Diesel_DER2_DW.voltages_DSTATE[1] = rtb_Fcn3_m;
    Diesel_DER2_DW.voltages_DSTATE[2] = rtb_Vfd;
    Diesel_DER2_DW.voltages_DSTATE[3] = Diesel_DER2_P.Vkd0Vkq10Vkq20_Value[0];
    Diesel_DER2_DW.voltages_DSTATE[4] = Diesel_DER2_P.Vkd0Vkq10Vkq20_Value[1];

    /* Update for Sin: '<S87>/sin(wt)' */
    HoldSine = Diesel_DER2_DW.lastSin;
    Diesel_DER2_DW.lastSin = Diesel_DER2_DW.lastSin * Diesel_DER2_P.sinwt_HCos +
      Diesel_DER2_DW.lastCos * Diesel_DER2_P.sinwt_Hsin;
    Diesel_DER2_DW.lastCos = Diesel_DER2_DW.lastCos * Diesel_DER2_P.sinwt_HCos -
      HoldSine * Diesel_DER2_P.sinwt_Hsin;

    /* Update for DiscreteIntegrator: '<S97>/Integ4' */
    Diesel_DER2_DW.Integ4_SYSTEM_ENABLE = 0U;
    Diesel_DER2_DW.Integ4_DSTATE = Diesel_DER2_P.Integ4_gainval * rtb_Product_n
      + Diesel_DER2_B.Integ4;

    /* Level2 S-Function Block: '<S98>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[0];
      sfcnUpdate(rts, 1);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for UnitDelay: '<S97>/Unit Delay' */
    Diesel_DER2_DW.UnitDelay_DSTATE = rtb_Product_n;

    /* Update for UnitDelay: '<S97>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE_c = rtb_Switch;

    /* Update for Sin: '<S87>/cos(wt)' */
    HoldSine = Diesel_DER2_DW.lastSin_d;
    Diesel_DER2_DW.lastSin_d = Diesel_DER2_DW.lastSin_d *
      Diesel_DER2_P.coswt_HCos + Diesel_DER2_DW.lastCos_e *
      Diesel_DER2_P.coswt_Hsin;
    Diesel_DER2_DW.lastCos_e = Diesel_DER2_DW.lastCos_e *
      Diesel_DER2_P.coswt_HCos - HoldSine * Diesel_DER2_P.coswt_Hsin;

    /* Update for DiscreteIntegrator: '<S95>/Integ4' */
    Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_f = 0U;
    Diesel_DER2_DW.Integ4_DSTATE_e = Diesel_DER2_P.Integ4_gainval_o *
      rtb_Product1_f + Diesel_DER2_B.Integ4_h;

    /* Level2 S-Function Block: '<S96>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[1];
      sfcnUpdate(rts, 1);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for UnitDelay: '<S95>/Unit Delay' */
    Diesel_DER2_DW.UnitDelay_DSTATE_c = rtb_Product1_f;

    /* Update for UnitDelay: '<S95>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE_l = rtb_Switch_k;

    /* Update for Sin: '<S88>/sin(wt)' */
    HoldSine = Diesel_DER2_DW.lastSin_a;
    Diesel_DER2_DW.lastSin_a = Diesel_DER2_DW.lastSin_a *
      Diesel_DER2_P.sinwt_HCos_p + Diesel_DER2_DW.lastCos_c *
      Diesel_DER2_P.sinwt_Hsin_j;
    Diesel_DER2_DW.lastCos_c = Diesel_DER2_DW.lastCos_c *
      Diesel_DER2_P.sinwt_HCos_p - HoldSine * Diesel_DER2_P.sinwt_Hsin_j;

    /* Update for DiscreteIntegrator: '<S103>/Integ4' */
    Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_j = 0U;
    Diesel_DER2_DW.Integ4_DSTATE_b = Diesel_DER2_P.Integ4_gainval_d *
      rtb_Product_o + Diesel_DER2_B.Integ4_m;

    /* Level2 S-Function Block: '<S104>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[2];
      sfcnUpdate(rts, 1);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for UnitDelay: '<S103>/Unit Delay' */
    Diesel_DER2_DW.UnitDelay_DSTATE_j = rtb_Product_o;

    /* Update for UnitDelay: '<S103>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE_e = rtb_Switch_i;

    /* Update for Sin: '<S88>/cos(wt)' */
    HoldSine = Diesel_DER2_DW.lastSin_f;
    Diesel_DER2_DW.lastSin_f = Diesel_DER2_DW.lastSin_f *
      Diesel_DER2_P.coswt_HCos_m + Diesel_DER2_DW.lastCos_eg *
      Diesel_DER2_P.coswt_Hsin_p;
    Diesel_DER2_DW.lastCos_eg = Diesel_DER2_DW.lastCos_eg *
      Diesel_DER2_P.coswt_HCos_m - HoldSine * Diesel_DER2_P.coswt_Hsin_p;

    /* Update for DiscreteIntegrator: '<S101>/Integ4' */
    Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_fj = 0U;
    Diesel_DER2_DW.Integ4_DSTATE_f = Diesel_DER2_P.Integ4_gainval_e *
      rtb_Product1_a + Diesel_DER2_B.Integ4_o;

    /* Level2 S-Function Block: '<S102>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[3];
      sfcnUpdate(rts, 1);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for UnitDelay: '<S101>/Unit Delay' */
    Diesel_DER2_DW.UnitDelay_DSTATE_l = rtb_Product1_a;

    /* Update for UnitDelay: '<S101>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE_ep = rtb_Switch_o;

    /* Update for Sin: '<S89>/sin(wt)' */
    HoldSine = Diesel_DER2_DW.lastSin_o;
    Diesel_DER2_DW.lastSin_o = Diesel_DER2_DW.lastSin_o *
      Diesel_DER2_P.sinwt_HCos_k + Diesel_DER2_DW.lastCos_h *
      Diesel_DER2_P.sinwt_Hsin_d;
    Diesel_DER2_DW.lastCos_h = Diesel_DER2_DW.lastCos_h *
      Diesel_DER2_P.sinwt_HCos_k - HoldSine * Diesel_DER2_P.sinwt_Hsin_d;

    /* Update for DiscreteIntegrator: '<S109>/Integ4' */
    Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_h = 0U;
    Diesel_DER2_DW.Integ4_DSTATE_o = Diesel_DER2_P.Integ4_gainval_l *
      rtb_Product_h + Diesel_DER2_B.Integ4_i;

    /* Level2 S-Function Block: '<S110>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[4];
      sfcnUpdate(rts, 1);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for UnitDelay: '<S109>/Unit Delay' */
    Diesel_DER2_DW.UnitDelay_DSTATE_h = rtb_Product_h;

    /* Update for UnitDelay: '<S109>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE_k = rtb_Switch_h;

    /* Update for Sin: '<S89>/cos(wt)' */
    HoldSine = Diesel_DER2_DW.lastSin_i;
    Diesel_DER2_DW.lastSin_i = Diesel_DER2_DW.lastSin_i *
      Diesel_DER2_P.coswt_HCos_b + Diesel_DER2_DW.lastCos_g *
      Diesel_DER2_P.coswt_Hsin_ph;
    Diesel_DER2_DW.lastCos_g = Diesel_DER2_DW.lastCos_g *
      Diesel_DER2_P.coswt_HCos_b - HoldSine * Diesel_DER2_P.coswt_Hsin_ph;

    /* Update for DiscreteIntegrator: '<S107>/Integ4' */
    Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_c = 0U;
    Diesel_DER2_DW.Integ4_DSTATE_l = Diesel_DER2_P.Integ4_gainval_a *
      rtb_Product1_e + Diesel_DER2_B.Integ4_k;

    /* Level2 S-Function Block: '<S108>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[5];
      sfcnUpdate(rts, 1);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for UnitDelay: '<S107>/Unit Delay' */
    Diesel_DER2_DW.UnitDelay_DSTATE_ch = rtb_Product1_e;

    /* Update for UnitDelay: '<S107>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE_j = rtb_Switch_of;

    /* Update for UnitDelay: '<S132>/voltages' incorporates:
     *  Constant: '<S117>/[ Vkd =0 Vkq1=0  Vkq2=0 ]'
     */
    Diesel_DER2_DW.voltages_DSTATE_h[0] = rtb_Fcn2;
    Diesel_DER2_DW.voltages_DSTATE_h[1] = rtb_Fcn3;
    Diesel_DER2_DW.voltages_DSTATE_h[2] = rtb_Vfd_o;
    Diesel_DER2_DW.voltages_DSTATE_h[3] = Diesel_DER2_P.Vkd0Vkq10Vkq20_Value_a[0];
    Diesel_DER2_DW.voltages_DSTATE_h[4] = Diesel_DER2_P.Vkd0Vkq10Vkq20_Value_a[1];

    /* Update for DiscreteIntegrator: '<S119>/Rotor speed deviation (dw)' */
    Diesel_DER2_DW.Rotorspeeddeviationdw_SYSTEM_EN = 0U;
    Diesel_DER2_DW.Rotorspeeddeviationdw_DSTATE =
      Diesel_DER2_P.Rotorspeeddeviationdw_gainval * rtb_uH + Diesel_DER2_B.dw;

    /* Update for DiscreteIntegrator: '<S111>/Integ4' */
    Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_d = 0U;
    Diesel_DER2_DW.Integ4_DSTATE_m = Diesel_DER2_P.Integ4_gainval_j *
      rtb_UnitDelay4 + Diesel_DER2_B.Integ4_a;

    /* Level2 S-Function Block: '<S112>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[6];
      sfcnUpdate(rts, 1);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Update for UnitDelay: '<S111>/Unit Delay' */
    Diesel_DER2_DW.UnitDelay_DSTATE_a = rtb_UnitDelay4;

    /* Update for UnitDelay: '<S111>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE_m = rtb_Switch_a;

    /* Update for UnitDelay: '<S146>/Delay Ts' */
    Diesel_DER2_DW.DelayTs_DSTATE = Diesel_DER2_B.StateSpace_o1[0];

    /* Update for UnitDelay: '<S146>/Delay Ts ' */
    Diesel_DER2_DW.DelayTs_DSTATE_f = rtb_Sum_k;

    /* Update for UnitDelay: '<S149>/Delay Ts' */
    Diesel_DER2_DW.DelayTs_DSTATE_a = Diesel_DER2_B.StateSpace_o1[1];

    /* Update for UnitDelay: '<S149>/Delay Ts ' */
    Diesel_DER2_DW.DelayTs_DSTATE_n = rtb_Sum_o;

    /* Update for UnitDelay: '<S152>/Delay Ts' */
    Diesel_DER2_DW.DelayTs_DSTATE_l = Diesel_DER2_B.StateSpace_o1[2];

    /* Update for UnitDelay: '<S152>/Delay Ts ' */
    Diesel_DER2_DW.DelayTs_DSTATE_d = rtb_Sum_b;
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   */
  Diesel_DER2_M->Timing.t[0] =
    (++Diesel_DER2_M->Timing.clockTick0) * Diesel_DER2_M->Timing.stepSize0;

  {
    /* Update absolute timer for sample time: [5.0E-5s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick1"
     * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
     * overflow during the application lifespan selected.
     */
    Diesel_DER2_M->Timing.t[1] =
      (++Diesel_DER2_M->Timing.clockTick1) * Diesel_DER2_M->Timing.stepSize1;
  }
}

/* Model initialize function */
void Diesel_DER2_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)Diesel_DER2_M, 0,
                sizeof(RT_MODEL_Diesel_DER2_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&Diesel_DER2_M->solverInfo,
                          &Diesel_DER2_M->Timing.simTimeStep);
    rtsiSetTPtr(&Diesel_DER2_M->solverInfo, &rtmGetTPtr(Diesel_DER2_M));
    rtsiSetStepSizePtr(&Diesel_DER2_M->solverInfo,
                       &Diesel_DER2_M->Timing.stepSize0);
    rtsiSetErrorStatusPtr(&Diesel_DER2_M->solverInfo, ((const char_T **)
      (&rtmGetErrorStatus(Diesel_DER2_M))));
    rtsiSetRTModelPtr(&Diesel_DER2_M->solverInfo, Diesel_DER2_M);
  }

  rtsiSetSimTimeStep(&Diesel_DER2_M->solverInfo, MAJOR_TIME_STEP);
  rtsiSetSolverName(&Diesel_DER2_M->solverInfo,"FixedStepDiscrete");
  Diesel_DER2_M->solverInfoPtr = (&Diesel_DER2_M->solverInfo);

  /* Initialize timing info */
  {
    int_T *mdlTsMap = Diesel_DER2_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    Diesel_DER2_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    Diesel_DER2_M->Timing.sampleTimes = (&Diesel_DER2_M->
      Timing.sampleTimesArray[0]);
    Diesel_DER2_M->Timing.offsetTimes = (&Diesel_DER2_M->
      Timing.offsetTimesArray[0]);

    /* task periods */
    Diesel_DER2_M->Timing.sampleTimes[0] = (0.0);
    Diesel_DER2_M->Timing.sampleTimes[1] = (5.0E-5);

    /* task offsets */
    Diesel_DER2_M->Timing.offsetTimes[0] = (0.0);
    Diesel_DER2_M->Timing.offsetTimes[1] = (0.0);
  }

  rtmSetTPtr(Diesel_DER2_M, &Diesel_DER2_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = Diesel_DER2_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    mdlSampleHits[1] = 1;
    Diesel_DER2_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(Diesel_DER2_M, -1);
  Diesel_DER2_M->Timing.stepSize0 = 5.0E-5;
  Diesel_DER2_M->Timing.stepSize1 = 5.0E-5;
  Diesel_DER2_M->solverInfoPtr = (&Diesel_DER2_M->solverInfo);
  Diesel_DER2_M->Timing.stepSize = (5.0E-5);
  rtsiSetFixedStepSize(&Diesel_DER2_M->solverInfo, 5.0E-5);
  rtsiSetSolverMode(&Diesel_DER2_M->solverInfo, SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  (void) memset(((void *) &Diesel_DER2_B), 0,
                sizeof(B_Diesel_DER2_T));

  /* states (dwork) */
  (void) memset((void *)&Diesel_DER2_DW, 0,
                sizeof(DW_Diesel_DER2_T));

  /* external outputs */
  (void) memset((void *)&Diesel_DER2_Y, 0,
                sizeof(ExtY_Diesel_DER2_T));

  /* child S-Function registration */
  {
    RTWSfcnInfo *sfcnInfo = &Diesel_DER2_M->NonInlinedSFcns.sfcnInfo;
    Diesel_DER2_M->sfcnInfo = (sfcnInfo);
    rtssSetErrorStatusPtr(sfcnInfo, ((const char_T **)(&rtmGetErrorStatus
      (Diesel_DER2_M))));
    rtssSetNumRootSampTimesPtr(sfcnInfo, &Diesel_DER2_M->Sizes.numSampTimes);
    Diesel_DER2_M->NonInlinedSFcns.taskTimePtrs[0] = &(rtmGetTPtr(Diesel_DER2_M)
      [0]);
    Diesel_DER2_M->NonInlinedSFcns.taskTimePtrs[1] = &(rtmGetTPtr(Diesel_DER2_M)
      [1]);
    rtssSetTPtrPtr(sfcnInfo,Diesel_DER2_M->NonInlinedSFcns.taskTimePtrs);
    rtssSetTStartPtr(sfcnInfo, &rtmGetTStart(Diesel_DER2_M));
    rtssSetTFinalPtr(sfcnInfo, &rtmGetTFinal(Diesel_DER2_M));
    rtssSetTimeOfLastOutputPtr(sfcnInfo, &rtmGetTimeOfLastOutput(Diesel_DER2_M));
    rtssSetStepSizePtr(sfcnInfo, &Diesel_DER2_M->Timing.stepSize);
    rtssSetStopRequestedPtr(sfcnInfo, &rtmGetStopRequested(Diesel_DER2_M));
    rtssSetDerivCacheNeedsResetPtr(sfcnInfo,
      &Diesel_DER2_M->ModelData.derivCacheNeedsReset);
    rtssSetZCCacheNeedsResetPtr(sfcnInfo,
      &Diesel_DER2_M->ModelData.zCCacheNeedsReset);
    rtssSetBlkStateChangePtr(sfcnInfo, &Diesel_DER2_M->ModelData.blkStateChange);
    rtssSetSampleHitsPtr(sfcnInfo, &Diesel_DER2_M->Timing.sampleHits);
    rtssSetPerTaskSampleHitsPtr(sfcnInfo,
      &Diesel_DER2_M->Timing.perTaskSampleHits);
    rtssSetSimModePtr(sfcnInfo, &Diesel_DER2_M->simMode);
    rtssSetSolverInfoPtr(sfcnInfo, &Diesel_DER2_M->solverInfoPtr);
  }

  Diesel_DER2_M->Sizes.numSFcns = (7);

  /* register each child */
  {
    (void) memset((void *)&Diesel_DER2_M->NonInlinedSFcns.childSFunctions[0], 0,
                  7*sizeof(SimStruct));
    Diesel_DER2_M->childSfunctions =
      (&Diesel_DER2_M->NonInlinedSFcns.childSFunctionPtrs[0]);

    {
      int_T i;
      for (i = 0; i < 7; i++) {
        Diesel_DER2_M->childSfunctions[i] =
          (&Diesel_DER2_M->NonInlinedSFcns.childSFunctions[i]);
      }
    }

    /* Level2 S-Function Block: Diesel_DER2/<S98>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[0];

      /* timing info */
      time_T *sfcnPeriod = Diesel_DER2_M->NonInlinedSFcns.Sfcn0.sfcnPeriod;
      time_T *sfcnOffset = Diesel_DER2_M->NonInlinedSFcns.Sfcn0.sfcnOffset;
      int_T *sfcnTsMap = Diesel_DER2_M->NonInlinedSFcns.Sfcn0.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts, &Diesel_DER2_M->NonInlinedSFcns.blkInfo2[0]);
      }

      ssSetRTWSfcnInfo(rts, Diesel_DER2_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &Diesel_DER2_M->NonInlinedSFcns.methods2[0]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &Diesel_DER2_M->NonInlinedSFcns.methods3[0]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &Diesel_DER2_M->NonInlinedSFcns.statesInfo2[0]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn0.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Diesel_DER2_M->NonInlinedSFcns.Sfcn0.UPtrs0;
          sfcnUPtrs[0] = &Diesel_DER2_B.Integ4;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Diesel_DER2_M->NonInlinedSFcns.Sfcn0.UPtrs1;
          sfcnUPtrs[0] = &Diesel_DER2_B.K1;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn0.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *) &Diesel_DER2_B.SFunction));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_A/Mean value1/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,Diesel_DER2_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn0.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)Diesel_DER2_P.SFunction_P1_Size);
        ssSetSFcnParam(rts, 1, (mxArray*)Diesel_DER2_P.SFunction_P2_Size);
        ssSetSFcnParam(rts, 2, (mxArray*)Diesel_DER2_P.SFunction_P3_Size);
        ssSetSFcnParam(rts, 3, (mxArray*)Diesel_DER2_P.SFunction_P4_Size);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *) &Diesel_DER2_DW.SFunction_RWORK);
      ssSetIWork(rts, (int_T *) &Diesel_DER2_DW.SFunction_IWORK);
      ssSetPWork(rts, (void **) &Diesel_DER2_DW.SFunction_PWORK);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn0.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn0.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &Diesel_DER2_DW.SFunction_RWORK);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &Diesel_DER2_DW.SFunction_IWORK);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &Diesel_DER2_DW.SFunction_PWORK);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
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

    /* Level2 S-Function Block: Diesel_DER2/<S96>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[1];

      /* timing info */
      time_T *sfcnPeriod = Diesel_DER2_M->NonInlinedSFcns.Sfcn1.sfcnPeriod;
      time_T *sfcnOffset = Diesel_DER2_M->NonInlinedSFcns.Sfcn1.sfcnOffset;
      int_T *sfcnTsMap = Diesel_DER2_M->NonInlinedSFcns.Sfcn1.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts, &Diesel_DER2_M->NonInlinedSFcns.blkInfo2[1]);
      }

      ssSetRTWSfcnInfo(rts, Diesel_DER2_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &Diesel_DER2_M->NonInlinedSFcns.methods2[1]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &Diesel_DER2_M->NonInlinedSFcns.methods3[1]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &Diesel_DER2_M->NonInlinedSFcns.statesInfo2[1]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn1.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Diesel_DER2_M->NonInlinedSFcns.Sfcn1.UPtrs0;
          sfcnUPtrs[0] = &Diesel_DER2_B.Integ4_h;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Diesel_DER2_M->NonInlinedSFcns.Sfcn1.UPtrs1;
          sfcnUPtrs[0] = &Diesel_DER2_B.K1_e;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn1.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *) &Diesel_DER2_B.SFunction_a));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_A/Mean/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,Diesel_DER2_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn1.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)Diesel_DER2_P.SFunction_P1_Size_j);
        ssSetSFcnParam(rts, 1, (mxArray*)Diesel_DER2_P.SFunction_P2_Size_l);
        ssSetSFcnParam(rts, 2, (mxArray*)Diesel_DER2_P.SFunction_P3_Size_k);
        ssSetSFcnParam(rts, 3, (mxArray*)Diesel_DER2_P.SFunction_P4_Size_f);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *) &Diesel_DER2_DW.SFunction_RWORK_k);
      ssSetIWork(rts, (int_T *) &Diesel_DER2_DW.SFunction_IWORK_e);
      ssSetPWork(rts, (void **) &Diesel_DER2_DW.SFunction_PWORK_c);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn1.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn1.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &Diesel_DER2_DW.SFunction_RWORK_k);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &Diesel_DER2_DW.SFunction_IWORK_e);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &Diesel_DER2_DW.SFunction_PWORK_c);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
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

    /* Level2 S-Function Block: Diesel_DER2/<S104>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[2];

      /* timing info */
      time_T *sfcnPeriod = Diesel_DER2_M->NonInlinedSFcns.Sfcn2.sfcnPeriod;
      time_T *sfcnOffset = Diesel_DER2_M->NonInlinedSFcns.Sfcn2.sfcnOffset;
      int_T *sfcnTsMap = Diesel_DER2_M->NonInlinedSFcns.Sfcn2.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts, &Diesel_DER2_M->NonInlinedSFcns.blkInfo2[2]);
      }

      ssSetRTWSfcnInfo(rts, Diesel_DER2_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &Diesel_DER2_M->NonInlinedSFcns.methods2[2]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &Diesel_DER2_M->NonInlinedSFcns.methods3[2]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &Diesel_DER2_M->NonInlinedSFcns.statesInfo2[2]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn2.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Diesel_DER2_M->NonInlinedSFcns.Sfcn2.UPtrs0;
          sfcnUPtrs[0] = &Diesel_DER2_B.Integ4_m;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Diesel_DER2_M->NonInlinedSFcns.Sfcn2.UPtrs1;
          sfcnUPtrs[0] = &Diesel_DER2_B.K1_c;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn2.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *) &Diesel_DER2_B.SFunction_e));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_B/Mean value1/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,Diesel_DER2_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn2.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)Diesel_DER2_P.SFunction_P1_Size_jj);
        ssSetSFcnParam(rts, 1, (mxArray*)Diesel_DER2_P.SFunction_P2_Size_d);
        ssSetSFcnParam(rts, 2, (mxArray*)Diesel_DER2_P.SFunction_P3_Size_a);
        ssSetSFcnParam(rts, 3, (mxArray*)Diesel_DER2_P.SFunction_P4_Size_c);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *) &Diesel_DER2_DW.SFunction_RWORK_m);
      ssSetIWork(rts, (int_T *) &Diesel_DER2_DW.SFunction_IWORK_m);
      ssSetPWork(rts, (void **) &Diesel_DER2_DW.SFunction_PWORK_m);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn2.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn2.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &Diesel_DER2_DW.SFunction_RWORK_m);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &Diesel_DER2_DW.SFunction_IWORK_m);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &Diesel_DER2_DW.SFunction_PWORK_m);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
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

    /* Level2 S-Function Block: Diesel_DER2/<S102>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[3];

      /* timing info */
      time_T *sfcnPeriod = Diesel_DER2_M->NonInlinedSFcns.Sfcn3.sfcnPeriod;
      time_T *sfcnOffset = Diesel_DER2_M->NonInlinedSFcns.Sfcn3.sfcnOffset;
      int_T *sfcnTsMap = Diesel_DER2_M->NonInlinedSFcns.Sfcn3.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts, &Diesel_DER2_M->NonInlinedSFcns.blkInfo2[3]);
      }

      ssSetRTWSfcnInfo(rts, Diesel_DER2_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &Diesel_DER2_M->NonInlinedSFcns.methods2[3]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &Diesel_DER2_M->NonInlinedSFcns.methods3[3]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &Diesel_DER2_M->NonInlinedSFcns.statesInfo2[3]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn3.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Diesel_DER2_M->NonInlinedSFcns.Sfcn3.UPtrs0;
          sfcnUPtrs[0] = &Diesel_DER2_B.Integ4_o;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Diesel_DER2_M->NonInlinedSFcns.Sfcn3.UPtrs1;
          sfcnUPtrs[0] = &Diesel_DER2_B.K1_j;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn3.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *) &Diesel_DER2_B.SFunction_i));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_B/Mean/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,Diesel_DER2_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn3.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)Diesel_DER2_P.SFunction_P1_Size_a);
        ssSetSFcnParam(rts, 1, (mxArray*)Diesel_DER2_P.SFunction_P2_Size_i);
        ssSetSFcnParam(rts, 2, (mxArray*)Diesel_DER2_P.SFunction_P3_Size_b);
        ssSetSFcnParam(rts, 3, (mxArray*)Diesel_DER2_P.SFunction_P4_Size_cp);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *) &Diesel_DER2_DW.SFunction_RWORK_j);
      ssSetIWork(rts, (int_T *) &Diesel_DER2_DW.SFunction_IWORK_k);
      ssSetPWork(rts, (void **) &Diesel_DER2_DW.SFunction_PWORK_a);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn3.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn3.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &Diesel_DER2_DW.SFunction_RWORK_j);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &Diesel_DER2_DW.SFunction_IWORK_k);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &Diesel_DER2_DW.SFunction_PWORK_a);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
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

    /* Level2 S-Function Block: Diesel_DER2/<S110>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[4];

      /* timing info */
      time_T *sfcnPeriod = Diesel_DER2_M->NonInlinedSFcns.Sfcn4.sfcnPeriod;
      time_T *sfcnOffset = Diesel_DER2_M->NonInlinedSFcns.Sfcn4.sfcnOffset;
      int_T *sfcnTsMap = Diesel_DER2_M->NonInlinedSFcns.Sfcn4.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts, &Diesel_DER2_M->NonInlinedSFcns.blkInfo2[4]);
      }

      ssSetRTWSfcnInfo(rts, Diesel_DER2_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &Diesel_DER2_M->NonInlinedSFcns.methods2[4]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &Diesel_DER2_M->NonInlinedSFcns.methods3[4]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &Diesel_DER2_M->NonInlinedSFcns.statesInfo2[4]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn4.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Diesel_DER2_M->NonInlinedSFcns.Sfcn4.UPtrs0;
          sfcnUPtrs[0] = &Diesel_DER2_B.Integ4_i;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Diesel_DER2_M->NonInlinedSFcns.Sfcn4.UPtrs1;
          sfcnUPtrs[0] = &Diesel_DER2_B.K1_k;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn4.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *) &Diesel_DER2_B.SFunction_g));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_C/Mean value1/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,Diesel_DER2_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn4.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)Diesel_DER2_P.SFunction_P1_Size_l);
        ssSetSFcnParam(rts, 1, (mxArray*)Diesel_DER2_P.SFunction_P2_Size_g);
        ssSetSFcnParam(rts, 2, (mxArray*)Diesel_DER2_P.SFunction_P3_Size_d);
        ssSetSFcnParam(rts, 3, (mxArray*)Diesel_DER2_P.SFunction_P4_Size_a);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *) &Diesel_DER2_DW.SFunction_RWORK_b);
      ssSetIWork(rts, (int_T *) &Diesel_DER2_DW.SFunction_IWORK_n);
      ssSetPWork(rts, (void **) &Diesel_DER2_DW.SFunction_PWORK_ao);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn4.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn4.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &Diesel_DER2_DW.SFunction_RWORK_b);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &Diesel_DER2_DW.SFunction_IWORK_n);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &Diesel_DER2_DW.SFunction_PWORK_ao);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
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

    /* Level2 S-Function Block: Diesel_DER2/<S108>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[5];

      /* timing info */
      time_T *sfcnPeriod = Diesel_DER2_M->NonInlinedSFcns.Sfcn5.sfcnPeriod;
      time_T *sfcnOffset = Diesel_DER2_M->NonInlinedSFcns.Sfcn5.sfcnOffset;
      int_T *sfcnTsMap = Diesel_DER2_M->NonInlinedSFcns.Sfcn5.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts, &Diesel_DER2_M->NonInlinedSFcns.blkInfo2[5]);
      }

      ssSetRTWSfcnInfo(rts, Diesel_DER2_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &Diesel_DER2_M->NonInlinedSFcns.methods2[5]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &Diesel_DER2_M->NonInlinedSFcns.methods3[5]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &Diesel_DER2_M->NonInlinedSFcns.statesInfo2[5]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn5.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Diesel_DER2_M->NonInlinedSFcns.Sfcn5.UPtrs0;
          sfcnUPtrs[0] = &Diesel_DER2_B.Integ4_k;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Diesel_DER2_M->NonInlinedSFcns.Sfcn5.UPtrs1;
          sfcnUPtrs[0] = &Diesel_DER2_B.K1_o;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn5.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *) &Diesel_DER2_B.SFunction_p));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_C/Mean/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,Diesel_DER2_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn5.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)Diesel_DER2_P.SFunction_P1_Size_b);
        ssSetSFcnParam(rts, 1, (mxArray*)Diesel_DER2_P.SFunction_P2_Size_la);
        ssSetSFcnParam(rts, 2, (mxArray*)Diesel_DER2_P.SFunction_P3_Size_a3);
        ssSetSFcnParam(rts, 3, (mxArray*)Diesel_DER2_P.SFunction_P4_Size_n);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *) &Diesel_DER2_DW.SFunction_RWORK_g);
      ssSetIWork(rts, (int_T *) &Diesel_DER2_DW.SFunction_IWORK_d);
      ssSetPWork(rts, (void **) &Diesel_DER2_DW.SFunction_PWORK_i);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn5.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn5.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &Diesel_DER2_DW.SFunction_RWORK_g);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &Diesel_DER2_DW.SFunction_IWORK_d);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &Diesel_DER2_DW.SFunction_PWORK_i);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
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

    /* Level2 S-Function Block: Diesel_DER2/<S112>/S-Function (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[6];

      /* timing info */
      time_T *sfcnPeriod = Diesel_DER2_M->NonInlinedSFcns.Sfcn6.sfcnPeriod;
      time_T *sfcnOffset = Diesel_DER2_M->NonInlinedSFcns.Sfcn6.sfcnOffset;
      int_T *sfcnTsMap = Diesel_DER2_M->NonInlinedSFcns.Sfcn6.sfcnTsMap;
      (void) memset((void*)sfcnPeriod, 0,
                    sizeof(time_T)*1);
      (void) memset((void*)sfcnOffset, 0,
                    sizeof(time_T)*1);
      ssSetSampleTimePtr(rts, &sfcnPeriod[0]);
      ssSetOffsetTimePtr(rts, &sfcnOffset[0]);
      ssSetSampleTimeTaskIDPtr(rts, sfcnTsMap);

      /* Set up the mdlInfo pointer */
      {
        ssSetBlkInfo2Ptr(rts, &Diesel_DER2_M->NonInlinedSFcns.blkInfo2[6]);
      }

      ssSetRTWSfcnInfo(rts, Diesel_DER2_M->sfcnInfo);

      /* Allocate memory of model methods 2 */
      {
        ssSetModelMethods2(rts, &Diesel_DER2_M->NonInlinedSFcns.methods2[6]);
      }

      /* Allocate memory of model methods 3 */
      {
        ssSetModelMethods3(rts, &Diesel_DER2_M->NonInlinedSFcns.methods3[6]);
      }

      /* Allocate memory for states auxilliary information */
      {
        ssSetStatesInfo2(rts, &Diesel_DER2_M->NonInlinedSFcns.statesInfo2[6]);
      }

      /* inputs */
      {
        _ssSetNumInputPorts(rts, 2);
        ssSetPortInfoForInputs(rts,
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn6.inputPortInfo[0]);

        /* port 0 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Diesel_DER2_M->NonInlinedSFcns.Sfcn6.UPtrs0;
          sfcnUPtrs[0] = &Diesel_DER2_B.Integ4_a;
          ssSetInputPortSignalPtrs(rts, 0, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 0, 1);
          ssSetInputPortWidth(rts, 0, 1);
        }

        /* port 1 */
        {
          real_T const **sfcnUPtrs = (real_T const **)
            &Diesel_DER2_M->NonInlinedSFcns.Sfcn6.UPtrs1;
          sfcnUPtrs[0] = &Diesel_DER2_B.K1_m;
          ssSetInputPortSignalPtrs(rts, 1, (InputPtrsType)&sfcnUPtrs[0]);
          _ssSetInputPortNumDimensions(rts, 1, 1);
          ssSetInputPortWidth(rts, 1, 1);
        }
      }

      /* outputs */
      {
        ssSetPortInfoForOutputs(rts,
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn6.outputPortInfo[0]);
        _ssSetNumOutputPorts(rts, 1);

        /* port 0 */
        {
          _ssSetOutputPortNumDimensions(rts, 0, 1);
          ssSetOutputPortWidth(rts, 0, 1);
          ssSetOutputPortSignal(rts, 0, ((real_T *) &Diesel_DER2_B.SFunction_iq));
        }
      }

      /* path info */
      ssSetModelName(rts, "S-Function");
      ssSetPath(rts,
                "Diesel_DER2/Mean/Model/Discrete Variable Time Delay/S-Function");
      ssSetRTModel(rts,Diesel_DER2_M);
      ssSetParentSS(rts, (NULL));
      ssSetRootSS(rts, rts);
      ssSetVersion(rts, SIMSTRUCT_VERSION_LEVEL2);

      /* parameters */
      {
        mxArray **sfcnParams = (mxArray **)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn6.params;
        ssSetSFcnParamsCount(rts, 4);
        ssSetSFcnParamsPtr(rts, &sfcnParams[0]);
        ssSetSFcnParam(rts, 0, (mxArray*)Diesel_DER2_P.SFunction_P1_Size_jv);
        ssSetSFcnParam(rts, 1, (mxArray*)Diesel_DER2_P.SFunction_P2_Size_k);
        ssSetSFcnParam(rts, 2, (mxArray*)Diesel_DER2_P.SFunction_P3_Size_i);
        ssSetSFcnParam(rts, 3, (mxArray*)Diesel_DER2_P.SFunction_P4_Size_p);
      }

      /* work vectors */
      ssSetRWork(rts, (real_T *) &Diesel_DER2_DW.SFunction_RWORK_g5);
      ssSetIWork(rts, (int_T *) &Diesel_DER2_DW.SFunction_IWORK_i);
      ssSetPWork(rts, (void **) &Diesel_DER2_DW.SFunction_PWORK_h);

      {
        struct _ssDWorkRecord *dWorkRecord = (struct _ssDWorkRecord *)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn6.dWork;
        struct _ssDWorkAuxRecord *dWorkAuxRecord = (struct _ssDWorkAuxRecord *)
          &Diesel_DER2_M->NonInlinedSFcns.Sfcn6.dWorkAux;
        ssSetSFcnDWork(rts, dWorkRecord);
        ssSetSFcnDWorkAux(rts, dWorkAuxRecord);
        _ssSetNumDWork(rts, 3);

        /* RWORK */
        ssSetDWorkWidth(rts, 0, 1);
        ssSetDWorkDataType(rts, 0,SS_DOUBLE);
        ssSetDWorkComplexSignal(rts, 0, 0);
        ssSetDWork(rts, 0, &Diesel_DER2_DW.SFunction_RWORK_g5);

        /* IWORK */
        ssSetDWorkWidth(rts, 1, 1);
        ssSetDWorkDataType(rts, 1,SS_INTEGER);
        ssSetDWorkComplexSignal(rts, 1, 0);
        ssSetDWork(rts, 1, &Diesel_DER2_DW.SFunction_IWORK_i);

        /* PWORK */
        ssSetDWorkWidth(rts, 2, 1);
        ssSetDWorkDataType(rts, 2,SS_POINTER);
        ssSetDWorkComplexSignal(rts, 2, 0);
        ssSetDWork(rts, 2, &Diesel_DER2_DW.SFunction_PWORK_h);
      }

      /* registration */
      sfun_discreteVariableDelay(rts);
      sfcnInitializeSizes(rts);
      sfcnInitializeSampleTimes(rts);

      /* adjust sample time */
      ssSetSampleTime(rts, 0, 5.0E-5);
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

    /* Start for Enabled SubSystem: '<S121>/Saturation' */

    /* Start for Enabled SubSystem: '<S131>/Lmq_sat' */
    Diesel_DER2_Lmq_sat_Start(&Diesel_DER2_DW.Lmq_sat_l,
      (P_Lmq_sat_Diesel_DER2_T *)&Diesel_DER2_P.Lmq_sat_l);

    /* End of Start for SubSystem: '<S131>/Lmq_sat' */

    /* End of Start for SubSystem: '<S121>/Saturation' */

    /* InitializeConditions for Enabled SubSystem: '<S121>/Saturation' */
    /* InitializeConditions for UnitDelay: '<S133>/Lmd_sat' */
    Diesel_DER2_DW.Lmd_sat_DSTATE = Diesel_DER2_P.Lmd_sat_InitialCondition_i;

    /* End of InitializeConditions for SubSystem: '<S121>/Saturation' */

    /* Start for Enabled SubSystem: '<S56>/Saturation' */

    /* Start for Enabled SubSystem: '<S66>/Lmq_sat' */
    Diesel_DER2_Lmq_sat_Start(&Diesel_DER2_DW.Lmq_sat, (P_Lmq_sat_Diesel_DER2_T *)
      &Diesel_DER2_P.Lmq_sat);

    /* End of Start for SubSystem: '<S66>/Lmq_sat' */

    /* End of Start for SubSystem: '<S56>/Saturation' */

    /* InitializeConditions for Enabled SubSystem: '<S56>/Saturation' */
    /* InitializeConditions for UnitDelay: '<S68>/Lmd_sat' */
    Diesel_DER2_DW.Lmd_sat_DSTATE_h = Diesel_DER2_P.Lmd_sat_InitialCondition;

    /* End of InitializeConditions for SubSystem: '<S56>/Saturation' */

    /* S-Function block: <S156>/State-Space */
    {
      Diesel_DER2_DW.StateSpace_PWORK.AS = (real_T*)calloc(6 * 6, sizeof(real_T));
      Diesel_DER2_DW.StateSpace_PWORK.BS = (real_T*)calloc(6 * 14, sizeof(real_T));
      Diesel_DER2_DW.StateSpace_PWORK.CS = (real_T*)calloc(42 * 6, sizeof(real_T));
      Diesel_DER2_DW.StateSpace_PWORK.DS = (real_T*)calloc(42 * 14, sizeof
        (real_T));
      Diesel_DER2_DW.StateSpace_PWORK.DX_COL = (real_T*)calloc(42, sizeof(real_T));
      Diesel_DER2_DW.StateSpace_PWORK.TMP2 = (real_T*)calloc(14, sizeof(real_T));
      Diesel_DER2_DW.StateSpace_PWORK.BD_COL = (real_T*)calloc(6, sizeof(real_T));
      Diesel_DER2_DW.StateSpace_PWORK.TMP1 = (real_T*)calloc(6, sizeof(real_T));
      Diesel_DER2_DW.StateSpace_PWORK.XTMP = (real_T*)calloc(6, sizeof(real_T));
      Diesel_DER2_DW.StateSpace_PWORK.CHOPPER = (int_T*)calloc(42, sizeof(int_T));
      Diesel_DER2_DW.StateSpace_PWORK.SWITCH_STATUS = (int_T*)calloc(9, sizeof
        (int_T));
      Diesel_DER2_DW.StateSpace_PWORK.SW_CHG = (int_T*)calloc(9, sizeof(int_T));
      Diesel_DER2_DW.StateSpace_PWORK.G_STATE = (int_T*)calloc(9, sizeof(int_T));
      Diesel_DER2_DW.StateSpace_PWORK.Y_SWITCH = (real_T*)calloc(9, sizeof
        (real_T));
      Diesel_DER2_DW.StateSpace_PWORK.SWITCH_TYPES = (int_T*)calloc(9, sizeof
        (int_T));
      Diesel_DER2_DW.StateSpace_PWORK.IDX_OUT_SW = (int_T*)calloc(9, sizeof
        (int_T));
      Diesel_DER2_DW.StateSpace_PWORK.SWITCH_STATUS_INIT = (int_T*)calloc(9,
        sizeof(int_T));
    }

    /* Start for TransportDelay: '<S2>/Delay Td' */
    {
      real_T *pBuffer = &Diesel_DER2_DW.DelayTd_RWORK.TUbufferArea[0];
      Diesel_DER2_DW.DelayTd_IWORK.Tail = 0;
      Diesel_DER2_DW.DelayTd_IWORK.Head = 0;
      Diesel_DER2_DW.DelayTd_IWORK.Last = 0;
      Diesel_DER2_DW.DelayTd_IWORK.CircularBufSize = 1024;
      pBuffer[0] = Diesel_DER2_P.DieselEngineSpeedRegulator_Pm0;
      pBuffer[1024] = Diesel_DER2_M->Timing.t[0];
      Diesel_DER2_DW.DelayTd_PWORK.TUbufferPtrs[0] = (void *) &pBuffer[0];
      Diesel_DER2_DW.DelayTd_PWORK.TUbufferPtrs[1] = (void *) &pBuffer[1024];
    }

    /* Start for Constant: '<S97>/K1' */
    Diesel_DER2_B.K1 = Diesel_DER2_P.K1_Value;

    /* Level2 S-Function Block: '<S98>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[0];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Start for Constant: '<S95>/K1' */
    Diesel_DER2_B.K1_e = Diesel_DER2_P.K1_Value_f;

    /* Level2 S-Function Block: '<S96>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[1];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Start for Constant: '<S103>/K1' */
    Diesel_DER2_B.K1_c = Diesel_DER2_P.K1_Value_i;

    /* Level2 S-Function Block: '<S104>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[2];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Start for Constant: '<S101>/K1' */
    Diesel_DER2_B.K1_j = Diesel_DER2_P.K1_Value_b;

    /* Level2 S-Function Block: '<S102>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[3];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Start for Constant: '<S109>/K1' */
    Diesel_DER2_B.K1_k = Diesel_DER2_P.K1_Value_g;

    /* Level2 S-Function Block: '<S110>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[4];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Start for Constant: '<S107>/K1' */
    Diesel_DER2_B.K1_o = Diesel_DER2_P.K1_Value_l;

    /* Level2 S-Function Block: '<S108>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[5];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Start for Constant: '<S111>/K1' */
    Diesel_DER2_B.K1_m = Diesel_DER2_P.K1_Value_fz;

    /* Level2 S-Function Block: '<S112>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[6];
      sfcnStart(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* Start for InitialCondition: '<S146>/IC' */
    Diesel_DER2_B.IC = Diesel_DER2_P.IC_Value;
    Diesel_DER2_DW.IC_FirstOutputTime = true;

    /* Start for InitialCondition: '<S149>/IC' */
    Diesel_DER2_B.IC_a = Diesel_DER2_P.IC_Value_j;
    Diesel_DER2_DW.IC_FirstOutputTime_l = true;

    /* Start for InitialCondition: '<S152>/IC' */
    Diesel_DER2_B.IC_c = Diesel_DER2_P.IC_Value_n;
    Diesel_DER2_DW.IC_FirstOutputTime_n = true;

    /* InitializeConditions for DiscreteIntegrator: '<S119>/Rotor angle dthetae' */
    Diesel_DER2_DW.Rotorangledthetae_DSTATE = Diesel_DER2_P.Rotorangledthetae_IC;

    /* InitializeConditions for DiscreteIntegrator: '<S54>/Rotor angle dtheta' */
    Diesel_DER2_DW.Rotorangledtheta_DSTATE = Diesel_DER2_P.Rotorangledtheta_IC;

    /* InitializeConditions for UnitDelay: '<S33>/Unit Delay2' */
    Diesel_DER2_DW.UnitDelay2_DSTATE = Diesel_DER2_P.UnitDelay2_InitialCondition;

    /* InitializeConditions for S-Function (sfun_spssw_discc): '<S156>/State-Space' */
    {
      real_T *As = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.AS;
      real_T *Bs = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.BS;
      real_T *Cs = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.CS;
      real_T *Ds = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.DS;
      real_T *X0 = (real_T*)&Diesel_DER2_DW.StateSpace_DSTATE[0];
      int_T *Chopper = (int_T*) Diesel_DER2_DW.StateSpace_PWORK.CHOPPER;
      X0[0] = 2.6353173715415909E+6;
      X0[1] = 1.8486890262445915E+6;
      X0[2] = 2.6486998503715885E+6;
      X0[3] = 2.1504273972634533E+6;
      X0[4] = 2.9370557425604556E+6;
      X0[5] = 2.13704491843346E+6;

      /* Copy and transpose A and B */
      As[0] = 0.75829636054664229;
      As[1] = 0.012481513606654542;
      As[2] = 0.12085181972667888;
      As[3] = -0.02496302721330914;
      As[4] = 0.12085181972667876;
      As[5] = 0.012481513606654594;
      As[6] = 0.012481513606654511;
      As[7] = 0.75829636054664018;
      As[8] = 0.01248151360665207;
      As[9] = 0.12085181972667877;
      As[10] = -0.024963027213306597;
      As[11] = 0.12085181972668135;
      As[12] = 0.12085181972667884;
      As[13] = 0.012481513606652056;
      As[14] = 0.75829636054663951;
      As[15] = 0.012481513606654546;
      As[16] = 0.12085181972668126;
      As[17] = -0.0249630272133066;
      As[18] = -0.024963027213309143;
      As[19] = 0.12085181972667884;
      As[20] = 0.012481513606654563;
      As[21] = 0.75829636054664251;
      As[22] = 0.012481513606654584;
      As[23] = 0.12085181972667877;
      As[24] = 0.12085181972667876;
      As[25] = -0.024963027213306552;
      As[26] = 0.12085181972668127;
      As[27] = 0.012481513606654539;
      As[28] = 0.75829636054663985;
      As[29] = 0.012481513606652013;
      As[30] = 0.012481513606654553;
      As[31] = 0.12085181972668133;
      As[32] = -0.024963027213306566;
      As[33] = 0.12085181972667881;
      As[34] = 0.012481513606652016;
      As[35] = 0.75829636054664007;
      Bs[0] = 0.0;
      Bs[1] = 0.0;
      Bs[2] = 0.0;
      Bs[3] = -7914.8180273321132;
      Bs[4] = -624.07568033272707;
      Bs[5] = -6042.5909863339439;
      Bs[6] = 1248.1513606654569;
      Bs[7] = -6042.59098633393;
      Bs[8] = -624.07568033271355;
      Bs[9] = 0.0;
      Bs[10] = 0.0;
      Bs[11] = 72054.845769439169;
      Bs[12] = 36027.422884719606;
      Bs[13] = 6666.6666666666561;
      Bs[14] = 0.0;
      Bs[15] = 0.0;
      Bs[16] = 0.0;
      Bs[17] = -624.0756803327256;
      Bs[18] = -7914.8180273319967;
      Bs[19] = -624.07568033260372;
      Bs[20] = -6042.5909863339366;
      Bs[21] = 1248.1513606653339;
      Bs[22] = -6042.5909863340585;
      Bs[23] = 0.0;
      Bs[24] = 0.0;
      Bs[25] = 36027.422884719628;
      Bs[26] = 72054.84576944093;
      Bs[27] = 6666.6666666666661;
      Bs[28] = 0.0;
      Bs[29] = 0.0;
      Bs[30] = 0.0;
      Bs[31] = -6042.59098633394;
      Bs[32] = -624.07568033260293;
      Bs[33] = -7914.8180273319886;
      Bs[34] = -624.07568033272707;
      Bs[35] = -6042.5909863340657;
      Bs[36] = 1248.1513606653205;
      Bs[37] = 0.0;
      Bs[38] = 0.0;
      Bs[39] = -36027.422884719592;
      Bs[40] = 36027.422884721265;
      Bs[41] = 6666.66666666667;
      Bs[42] = 0.0;
      Bs[43] = 0.0;
      Bs[44] = 0.0;
      Bs[45] = 1248.1513606654573;
      Bs[46] = -6042.5909863339421;
      Bs[47] = -624.0756803327281;
      Bs[48] = -7914.8180273321213;
      Bs[49] = -624.0756803327231;
      Bs[50] = -6042.5909863339366;
      Bs[51] = 0.0;
      Bs[52] = 0.0;
      Bs[53] = -72054.845769439227;
      Bs[54] = -36027.422884719614;
      Bs[55] = 6666.6666666666661;
      Bs[56] = 0.0;
      Bs[57] = 0.0;
      Bs[58] = 0.0;
      Bs[59] = -6042.5909863339375;
      Bs[60] = 1248.1513606653273;
      Bs[61] = -6042.5909863340639;
      Bs[62] = -624.075680332727;
      Bs[63] = -7914.8180273319913;
      Bs[64] = -624.0756803326002;
      Bs[65] = 0.0;
      Bs[66] = 0.0;
      Bs[67] = -36027.4228847196;
      Bs[68] = -72054.845769440915;
      Bs[69] = 6666.6666666666642;
      Bs[70] = 0.0;
      Bs[71] = 0.0;
      Bs[72] = 0.0;
      Bs[73] = -624.07568033272776;
      Bs[74] = -6042.5909863340657;
      Bs[75] = 1248.1513606653284;
      Bs[76] = -6042.5909863339393;
      Bs[77] = -624.07568033260077;
      Bs[78] = -7914.8180273319949;
      Bs[79] = 0.0;
      Bs[80] = 0.0;
      Bs[81] = 36027.422884719606;
      Bs[82] = -36027.422884721309;
      Bs[83] = 6666.666666666667;
      Cs[0] = 0.0;
      Cs[1] = 0.0;
      Cs[2] = 0.0;
      Cs[3] = 0.0;
      Cs[4] = 0.0;
      Cs[5] = 0.0;
      Cs[6] = 0.0;
      Cs[7] = 0.0;
      Cs[8] = 0.0;
      Cs[9] = 0.0;
      Cs[10] = 0.0;
      Cs[11] = 0.0;
      Cs[12] = 0.0;
      Cs[13] = 0.0;
      Cs[14] = 0.0;
      Cs[15] = 0.0;
      Cs[16] = 0.0;
      Cs[17] = 0.0;
      Cs[18] = 1.978704506814405E-5;
      Cs[19] = 1.5601892009516626E-6;
      Cs[20] = 1.5106477465932938E-5;
      Cs[21] = -3.1203784017886056E-6;
      Cs[22] = 1.5106477465923008E-5;
      Cs[23] = 1.5601892008369424E-6;
      Cs[24] = 1.5601892015672139E-6;
      Cs[25] = 1.97870450675784E-5;
      Cs[26] = 1.5601892000634461E-6;
      Cs[27] = 1.510647746648361E-5;
      Cs[28] = -3.120378401630661E-6;
      Cs[29] = 1.5106477465938002E-5;
      Cs[30] = 1.5106477465788948E-5;
      Cs[31] = 1.5601892008651157E-6;
      Cs[32] = 1.9787045068288554E-5;
      Cs[33] = 1.5601892007789547E-6;
      Cs[34] = 1.5106477465922491E-5;
      Cs[35] = -3.1203784016440704E-6;
      Cs[36] = -3.1203784009083232E-6;
      Cs[37] = 1.5106477465194613E-5;
      Cs[38] = 1.5601892000644066E-6;
      Cs[39] = 1.9787045068872007E-5;
      Cs[40] = 1.560189200843917E-6;
      Cs[41] = 1.5106477465933389E-5;
      Cs[42] = 1.5106477465668508E-5;
      Cs[43] = -3.1203784014321246E-6;
      Cs[44] = 1.5106477465933897E-5;
      Cs[45] = 1.5601892005997911E-6;
      Cs[46] = 1.9787045068397584E-5;
      Cs[47] = 1.5601892008323326E-6;
      Cs[48] = 1.5601892014467764E-6;
      Cs[49] = 1.5106477465281163E-5;
      Cs[50] = -3.120378402291207E-6;
      Cs[51] = 1.510647746630445E-5;
      Cs[52] = 1.5601892008444314E-6;
      Cs[53] = 1.9787045068414406E-5;
      Cs[54] = 0.0;
      Cs[55] = 0.0;
      Cs[56] = 0.0;
      Cs[57] = 0.0;
      Cs[58] = 0.0;
      Cs[59] = 0.0;
      Cs[60] = 0.0;
      Cs[61] = 0.0;
      Cs[62] = 0.0;
      Cs[63] = 0.0;
      Cs[64] = 0.0;
      Cs[65] = 0.0;
      Cs[66] = 9.006855720660424E-5;
      Cs[67] = -9.0068557208439849E-5;
      Cs[68] = -0.00018013711441808005;
      Cs[69] = -9.00685572145316E-5;
      Cs[70] = 9.0068557211475764E-5;
      Cs[71] = 0.00018013711442297142;
      Cs[72] = 9.0068557213114658E-5;
      Cs[73] = 0.00018013711442000036;
      Cs[74] = 9.0068557209764392E-5;
      Cs[75] = -9.006855720863056E-5;
      Cs[76] = -0.00018013711442287905;
      Cs[77] = -9.0068557211369824E-5;
      Cs[78] = 0.0;
      Cs[79] = 0.0;
      Cs[80] = 0.0;
      Cs[81] = 0.0;
      Cs[82] = 0.0;
      Cs[83] = 0.0;
      Cs[84] = 0.0;
      Cs[85] = 0.0;
      Cs[86] = 0.0;
      Cs[87] = 0.0;
      Cs[88] = 0.0;
      Cs[89] = 0.0;
      Cs[90] = 0.0;
      Cs[91] = 0.0;
      Cs[92] = 0.0;
      Cs[93] = 0.0;
      Cs[94] = 0.0;
      Cs[95] = 0.0;
      Cs[96] = -1.6666666667235721E-5;
      Cs[97] = -1.6666666666146281E-5;
      Cs[98] = -1.6666666665997345E-5;
      Cs[99] = -1.66666666670834E-5;
      Cs[100] = -1.6666666666766926E-5;
      Cs[101] = -1.6666666666770331E-5;
      Cs[102] = 9.00685572089775E-5;
      Cs[103] = 1.2387701975313803E-15;
      Cs[104] = -9.0068557208997129E-5;
      Cs[105] = -9.0068557212767645E-5;
      Cs[106] = 1.9611656049134397E-17;
      Cs[107] = 9.0068557211528876E-5;
      Cs[108] = 2.3732651260877091E-15;
      Cs[109] = 9.0068557209678618E-5;
      Cs[110] = 9.0068557209082876E-5;
      Cs[111] = 1.7639431582061605E-15;
      Cs[112] = -9.0068557211456167E-5;
      Cs[113] = -9.006855721144256E-5;
      Cs[114] = -9.006855721074138E-5;
      Cs[115] = -9.0068557210321726E-5;
      Cs[116] = -6.81516278634712E-16;
      Cs[117] = 9.00685572103945E-5;
      Cs[118] = 9.00685572114229E-5;
      Cs[119] = -7.273471173362432E-17;
      Cs[120] = 0.0;
      Cs[121] = 0.0;
      Cs[122] = 0.0;
      Cs[123] = 0.0;
      Cs[124] = 0.0;
      Cs[125] = 0.0;
      Cs[126] = 0.0;
      Cs[127] = 0.0;
      Cs[128] = 0.0;
      Cs[129] = 0.0;
      Cs[130] = 0.0;
      Cs[131] = 0.0;
      Cs[132] = 0.0;
      Cs[133] = 0.0;
      Cs[134] = 0.0;
      Cs[135] = 0.0;
      Cs[136] = 0.0;
      Cs[137] = 0.0;
      Cs[138] = -2.8146424128687244E-7;
      Cs[139] = -1.9138878849800369E-20;
      Cs[140] = 2.8146424128685333E-7;
      Cs[141] = 2.8146424128687181E-7;
      Cs[142] = 1.923171366081944E-20;
      Cs[143] = -2.814642412868527E-7;
      Cs[144] = 1.3992984288641043E-22;
      Cs[145] = -2.8146424128690431E-7;
      Cs[146] = -2.8146424128690426E-7;
      Cs[147] = 6.5052130349130272E-23;
      Cs[148] = 2.8146424128690421E-7;
      Cs[149] = 2.8146424128690415E-7;
      Cs[150] = 2.814642412868717E-7;
      Cs[151] = 2.8146424128685296E-7;
      Cs[152] = -1.9126097442337397E-20;
      Cs[153] = -2.8146424128687191E-7;
      Cs[154] = -2.8146424128685259E-7;
      Cs[155] = 1.9107751022110992E-20;
      Cs[156] = 1.978704506814405E-5;
      Cs[157] = 1.5601892009516626E-6;
      Cs[158] = 1.5106477465932938E-5;
      Cs[159] = -3.1203784017886056E-6;
      Cs[160] = 1.5106477465923008E-5;
      Cs[161] = 1.5601892008369424E-6;
      Cs[162] = 1.5601892015672139E-6;
      Cs[163] = 1.97870450675784E-5;
      Cs[164] = 1.5601892000634461E-6;
      Cs[165] = 1.510647746648361E-5;
      Cs[166] = -3.120378401630661E-6;
      Cs[167] = 1.5106477465938002E-5;
      Cs[168] = 1.5106477465788948E-5;
      Cs[169] = 1.5601892008651157E-6;
      Cs[170] = 1.9787045068288554E-5;
      Cs[171] = 1.5601892007789547E-6;
      Cs[172] = 1.5106477465922491E-5;
      Cs[173] = -3.1203784016440704E-6;
      Cs[174] = -3.1203784009083232E-6;
      Cs[175] = 1.5106477465194613E-5;
      Cs[176] = 1.5601892000644066E-6;
      Cs[177] = 1.9787045068872007E-5;
      Cs[178] = 1.560189200843917E-6;
      Cs[179] = 1.5106477465933389E-5;
      Cs[180] = 1.5106477465668508E-5;
      Cs[181] = -3.1203784014321246E-6;
      Cs[182] = 1.5106477465933897E-5;
      Cs[183] = 1.5601892005997911E-6;
      Cs[184] = 1.9787045068397584E-5;
      Cs[185] = 1.5601892008323326E-6;
      Cs[186] = 1.5601892014467764E-6;
      Cs[187] = 1.5106477465281163E-5;
      Cs[188] = -3.120378402291207E-6;
      Cs[189] = 1.510647746630445E-5;
      Cs[190] = 1.5601892008444314E-6;
      Cs[191] = 1.9787045068414406E-5;
      Cs[192] = 4.6805676023551E-6;
      Cs[193] = 8.65463516824816E-17;
      Cs[194] = -4.6805676023556129E-6;
      Cs[195] = -4.6805676025675609E-6;
      Cs[196] = 5.1460571914851987E-19;
      Cs[197] = 4.6805676024810136E-6;
      Cs[198] = 1.2043775265285773E-16;
      Cs[199] = 4.68056760229724E-6;
      Cs[200] = 4.6805676023546549E-6;
      Cs[201] = 1.7916154670949425E-16;
      Cs[202] = -4.680567602475092E-6;
      Cs[203] = -4.6805676024764041E-6;
      Cs[204] = -4.6805676024755367E-6;
      Cs[205] = -4.6805676023837869E-6;
      Cs[206] = 9.6020488339287821E-19;
      Cs[207] = 4.6805676023883965E-6;
      Cs[208] = 4.6805676024745779E-6;
      Cs[209] = -4.6102632928492672E-18;
      Cs[210] = -1.6666666667235721E-5;
      Cs[211] = -1.6666666666146281E-5;
      Cs[212] = -1.6666666665997345E-5;
      Cs[213] = -1.66666666670834E-5;
      Cs[214] = -1.6666666666766926E-5;
      Cs[215] = -1.6666666666770331E-5;
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
      Cs[226] = 0.0;
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
      Ds[0] = -0.47999976960011076;
      Ds[1] = 0.0;
      Ds[2] = 0.0;
      Ds[3] = 0.0;
      Ds[4] = 0.0;
      Ds[5] = 0.0;
      Ds[6] = 0.0;
      Ds[7] = 0.0;
      Ds[8] = 0.0;
      Ds[9] = 0.31999984640007384;
      Ds[10] = 0.0;
      Ds[11] = 0.0;
      Ds[12] = 0.0;
      Ds[13] = 0.0;
      Ds[14] = 0.0;
      Ds[15] = -0.47999976960011076;
      Ds[16] = 0.0;
      Ds[17] = 0.0;
      Ds[18] = 0.0;
      Ds[19] = 0.0;
      Ds[20] = 0.0;
      Ds[21] = 0.0;
      Ds[22] = 0.0;
      Ds[23] = 0.0;
      Ds[24] = 0.31999984640007384;
      Ds[25] = 0.0;
      Ds[26] = 0.0;
      Ds[27] = 0.0;
      Ds[28] = 0.0;
      Ds[29] = 0.0;
      Ds[30] = -0.47999976960011076;
      Ds[31] = 0.0;
      Ds[32] = 0.0;
      Ds[33] = 0.0;
      Ds[34] = 0.0;
      Ds[35] = 0.0;
      Ds[36] = 0.0;
      Ds[37] = -0.31999984640007384;
      Ds[38] = -0.31999984640007384;
      Ds[39] = 0.0;
      Ds[40] = 0.0;
      Ds[41] = 0.0;
      Ds[42] = 0.0;
      Ds[43] = 0.0;
      Ds[44] = 0.0;
      Ds[45] = -0.98935225341968913;
      Ds[46] = -0.078009460041332163;
      Ds[47] = -0.75532387329563544;
      Ds[48] = 0.15601892008220672;
      Ds[49] = -0.75532387329615036;
      Ds[50] = -0.078009460041847;
      Ds[51] = 0.0;
      Ds[52] = 0.0;
      Ds[53] = 9.0068557211433316;
      Ds[54] = 4.5034278605568012;
      Ds[55] = 0.8333333333374825;
      Ds[56] = 0.0;
      Ds[57] = 0.0;
      Ds[58] = 0.0;
      Ds[59] = -0.078009460041794418;
      Ds[60] = -0.9893522534193987;
      Ds[61] = -0.078009460040965789;
      Ds[62] = -0.7553238732960712;
      Ds[63] = 0.15601892008153312;
      Ds[64] = -0.75532387329689976;
      Ds[65] = 0.0;
      Ds[66] = 0.0;
      Ds[67] = 4.5034278605785847;
      Ds[68] = 9.0068557211093765;
      Ds[69] = 0.83333333333786552;
      Ds[70] = 0.0;
      Ds[71] = 0.0;
      Ds[72] = 0.0;
      Ds[73] = -0.75532387329617556;
      Ds[74] = -0.07800946004174078;
      Ds[75] = -0.98935225342006894;
      Ds[76] = -0.078009460041689779;
      Ds[77] = -0.7553238732961246;
      Ds[78] = 0.15601892008220355;
      Ds[79] = 0.0;
      Ds[80] = 0.0;
      Ds[81] = -4.5034278605725477;
      Ds[82] = 4.5034278605755409;
      Ds[83] = 0.83333333333786541;
      Ds[84] = 0.0;
      Ds[85] = 0.0;
      Ds[86] = 0.0;
      Ds[87] = 0.15601892008203644;
      Ds[88] = -0.75532387329528672;
      Ds[89] = -0.078009460040813258;
      Ds[90] = -0.98935225341951893;
      Ds[91] = -0.078009460042195777;
      Ds[92] = -0.75532387329666928;
      Ds[93] = 0.0;
      Ds[94] = 0.0;
      Ds[95] = -9.0068557211335047;
      Ds[96] = -4.50342786060667;
      Ds[97] = 0.8333333333374825;
      Ds[98] = 0.0;
      Ds[99] = 0.0;
      Ds[100] = 0.0;
      Ds[101] = -0.75532387329585826;
      Ds[102] = 0.15601892008277979;
      Ds[103] = -0.75532387329548278;
      Ds[104] = -0.0780094600412411;
      Ds[105] = -0.9893522534198792;
      Ds[106] = -0.078009460041616574;
      Ds[107] = 0.0;
      Ds[108] = 0.0;
      Ds[109] = -4.503427860568757;
      Ds[110] = -9.0068557211592442;
      Ds[111] = 0.83333333333709947;
      Ds[112] = 0.0;
      Ds[113] = 0.0;
      Ds[114] = 0.0;
      Ds[115] = -0.078009460041477088;
      Ds[116] = -0.75532387329487816;
      Ds[117] = 0.15601892008362025;
      Ds[118] = -0.75532387329562245;
      Ds[119] = -0.078009460042221493;
      Ds[120] = -0.98935225342071986;
      Ds[121] = 0.0;
      Ds[122] = 0.0;
      Ds[123] = 4.5034278605823763;
      Ds[124] = -4.5034278606254095;
      Ds[125] = 0.83333333333709958;
      Ds[126] = -0.31999984640007384;
      Ds[127] = 0.31999984640007384;
      Ds[128] = 0.0;
      Ds[129] = 0.0;
      Ds[130] = 0.0;
      Ds[131] = 0.0;
      Ds[132] = 0.0;
      Ds[133] = 0.0;
      Ds[134] = 0.0;
      Ds[135] = 0.31999989760004927;
      Ds[136] = -0.31999989760004927;
      Ds[137] = 0.0;
      Ds[138] = 0.0;
      Ds[139] = 0.0;
      Ds[140] = 0.0;
      Ds[141] = -0.31999984640007384;
      Ds[142] = 0.31999984640007384;
      Ds[143] = 0.0;
      Ds[144] = 0.0;
      Ds[145] = 0.0;
      Ds[146] = 0.0;
      Ds[147] = 0.0;
      Ds[148] = 0.0;
      Ds[149] = 0.31999989760004927;
      Ds[150] = 0.63999979520009842;
      Ds[151] = 0.0;
      Ds[152] = 0.0;
      Ds[153] = 0.0;
      Ds[154] = 0.0;
      Ds[155] = 0.0;
      Ds[156] = 0.0;
      Ds[157] = -4.5034278605575695;
      Ds[158] = 4.5034278605959006;
      Ds[159] = 9.0068557211706821;
      Ds[160] = 4.5034278605796816;
      Ds[161] = -4.5034278605737885;
      Ds[162] = -9.00685572114857;
      Ds[163] = 0.0;
      Ds[164] = 0.0;
      Ds[165] = 260.01152344517175;
      Ds[166] = -260.01152344611097;
      Ds[167] = -2.2112280914135589E-11;
      Ds[168] = 0.0;
      Ds[169] = 0.0;
      Ds[170] = 0.0;
      Ds[171] = -4.5034278605803619;
      Ds[172] = -9.0068557211660671;
      Ds[173] = -4.5034278605906044;
      Ds[174] = 4.503427860558249;
      Ds[175] = 9.0068557211439533;
      Ds[176] = 4.5034278605684914;
      Ds[177] = 0.0;
      Ds[178] = 0.0;
      Ds[179] = 260.01152344521097;
      Ds[180] = 520.02304689101618;
      Ds[181] = 2.2112602079452153E-11;
      Ds[182] = -0.31999984640007384;
      Ds[183] = 0.0;
      Ds[184] = 0.0;
      Ds[185] = 0.0;
      Ds[186] = 0.0;
      Ds[187] = 0.0;
      Ds[188] = 0.0;
      Ds[189] = 0.0;
      Ds[190] = 0.0;
      Ds[191] = 0.31999989760004927;
      Ds[192] = 0.0;
      Ds[193] = 0.0;
      Ds[194] = 0.0;
      Ds[195] = 0.0;
      Ds[196] = 0.0;
      Ds[197] = -0.31999984640007384;
      Ds[198] = 0.0;
      Ds[199] = 0.0;
      Ds[200] = 0.0;
      Ds[201] = 0.0;
      Ds[202] = 0.0;
      Ds[203] = 0.0;
      Ds[204] = 0.0;
      Ds[205] = 0.0;
      Ds[206] = 0.31999989760004927;
      Ds[207] = 0.0;
      Ds[208] = 0.0;
      Ds[209] = 0.0;
      Ds[210] = 0.0;
      Ds[211] = 0.0;
      Ds[212] = -0.31999984640007384;
      Ds[213] = 0.0;
      Ds[214] = 0.0;
      Ds[215] = 0.0;
      Ds[216] = 0.0;
      Ds[217] = 0.0;
      Ds[218] = 0.0;
      Ds[219] = -0.31999989760004927;
      Ds[220] = -0.31999989760004927;
      Ds[221] = 0.0;
      Ds[222] = 0.0;
      Ds[223] = 0.0;
      Ds[224] = 0.0;
      Ds[225] = 0.0;
      Ds[226] = 0.0;
      Ds[227] = 0.8333333333376528;
      Ds[228] = 0.833333333336619;
      Ds[229] = 0.83333333333644866;
      Ds[230] = 0.83333333333731219;
      Ds[231] = 0.83333333333834614;
      Ds[232] = 0.83333333333851634;
      Ds[233] = 0.0;
      Ds[234] = 0.0;
      Ds[235] = -9.8280878546574965E-12;
      Ds[236] = 4.9867838924910758E-11;
      Ds[237] = -1.666666666674965;
      Ds[238] = 0.0;
      Ds[239] = 0.0;
      Ds[240] = 0.0;
      Ds[241] = -4.5034278605656208;
      Ds[242] = 8.3506691853970285E-12;
      Ds[243] = 4.5034278605838125;
      Ds[244] = 4.5034278605729918;
      Ds[245] = -9.8019228289558656E-13;
      Ds[246] = -4.5034278605764415;
      Ds[247] = 0.0;
      Ds[248] = 0.0;
      Ds[249] = 260.011523445211;
      Ds[250] = -3.8875143673067216E-10;
      Ds[251] = -7.3708321224330574E-12;
      Ds[252] = 0.0;
      Ds[253] = 0.0;
      Ds[254] = 0.0;
      Ds[255] = -8.0513144934477636E-12;
      Ds[256] = -4.50342786058755;
      Ds[257] = -4.5034278605868678;
      Ds[258] = -6.689901240237437E-12;
      Ds[259] = 4.503427860572808;
      Ds[260] = 4.5034278605721267;
      Ds[261] = 0.0;
      Ds[262] = 0.0;
      Ds[263] = 3.9306883081735579E-11;
      Ds[264] = 260.01152344572222;
      Ds[265] = 1.4741357842232353E-11;
      Ds[266] = 0.0;
      Ds[267] = 0.0;
      Ds[268] = 0.0;
      Ds[269] = 4.5034278605723106;
      Ds[270] = 4.5034278605785163;
      Ds[271] = 3.7355593588557224E-12;
      Ds[272] = -4.5034278605649387;
      Ds[273] = -4.5034278605711453;
      Ds[274] = 3.6355633588837213E-12;
      Ds[275] = 0.0;
      Ds[276] = 0.0;
      Ds[277] = -260.01152344517169;
      Ds[278] = -260.011523445294;
      Ds[279] = -7.371221988839216E-12;
      Ds[280] = 0.0;
      Ds[281] = 0.0;
      Ds[282] = 0.0;
      Ds[283] = 0.0;
      Ds[284] = 0.0;
      Ds[285] = 0.0;
      Ds[286] = 0.0;
      Ds[287] = 0.0;
      Ds[288] = 0.0;
      Ds[289] = -1.0;
      Ds[290] = 0.0;
      Ds[291] = 0.0;
      Ds[292] = 0.0;
      Ds[293] = 0.0;
      Ds[294] = 0.0;
      Ds[295] = 0.0;
      Ds[296] = 0.0;
      Ds[297] = 0.0;
      Ds[298] = 0.0;
      Ds[299] = 0.0;
      Ds[300] = 0.0;
      Ds[301] = 0.0;
      Ds[302] = 0.0;
      Ds[303] = 0.0;
      Ds[304] = -1.0;
      Ds[305] = 0.0;
      Ds[306] = 0.0;
      Ds[307] = 0.0;
      Ds[308] = 0.0;
      Ds[309] = 0.0;
      Ds[310] = 0.0;
      Ds[311] = 0.0;
      Ds[312] = 0.0;
      Ds[313] = 0.0;
      Ds[314] = 0.0;
      Ds[315] = 0.0;
      Ds[316] = 0.0;
      Ds[317] = 1.0;
      Ds[318] = 1.0;
      Ds[319] = 0.0;
      Ds[320] = 0.0;
      Ds[321] = 0.0;
      Ds[322] = 0.0;
      Ds[323] = 0.0;
      Ds[324] = 0.0;
      Ds[325] = 0.014073212064343624;
      Ds[326] = 9.570788606083624E-16;
      Ds[327] = -0.014073212064342666;
      Ds[328] = -0.014073212064343588;
      Ds[329] = -9.21002163423168E-16;
      Ds[330] = 0.014073212064342703;
      Ds[331] = 0.0;
      Ds[332] = 0.0;
      Ds[333] = 0.18746398923033453;
      Ds[334] = 1.2797496395933194E-14;
      Ds[335] = -3.5321745528449354E-17;
      Ds[336] = 0.0;
      Ds[337] = 0.0;
      Ds[338] = 0.0;
      Ds[339] = -6.7834626804597066E-18;
      Ds[340] = 0.014073212064345213;
      Ds[341] = 0.014073212064345215;
      Ds[342] = -2.9365399001335392E-18;
      Ds[343] = -0.014073212064345223;
      Ds[344] = -0.014073212064345225;
      Ds[345] = 0.0;
      Ds[346] = 0.0;
      Ds[347] = -2.1316282072803008E-18;
      Ds[348] = 0.18746398923035612;
      Ds[349] = 1.0053069487980793E-17;
      Ds[350] = 0.0;
      Ds[351] = 0.0;
      Ds[352] = 0.0;
      Ds[353] = -0.014073212064343586;
      Ds[354] = -0.014073212064342644;
      Ds[355] = 9.56041561846283E-16;
      Ds[356] = 0.014073212064343595;
      Ds[357] = 0.014073212064342651;
      Ds[358] = -9.4818019131622815E-16;
      Ds[359] = 0.0;
      Ds[360] = 0.0;
      Ds[361] = -0.18746398923033458;
      Ds[362] = -0.18746398923032193;
      Ds[363] = -8.0834151349797339E-18;
      Ds[364] = 0.0;
      Ds[365] = 0.0;
      Ds[366] = 0.0;
      Ds[367] = -0.98935225341968913;
      Ds[368] = -0.078009460041332163;
      Ds[369] = -0.75532387329563544;
      Ds[370] = 0.15601892008220672;
      Ds[371] = -0.75532387329615036;
      Ds[372] = -0.078009460041847;
      Ds[373] = 0.0;
      Ds[374] = 0.0;
      Ds[375] = 9.0068557211433316;
      Ds[376] = 4.5034278605568012;
      Ds[377] = 0.8333333333374825;
      Ds[378] = 0.0;
      Ds[379] = 0.0;
      Ds[380] = 0.0;
      Ds[381] = -0.078009460041794418;
      Ds[382] = -0.9893522534193987;
      Ds[383] = -0.078009460040965789;
      Ds[384] = -0.7553238732960712;
      Ds[385] = 0.15601892008153312;
      Ds[386] = -0.75532387329689976;
      Ds[387] = 0.0;
      Ds[388] = 0.0;
      Ds[389] = 4.5034278605785847;
      Ds[390] = 9.0068557211093765;
      Ds[391] = 0.83333333333786552;
      Ds[392] = 0.0;
      Ds[393] = 0.0;
      Ds[394] = 0.0;
      Ds[395] = -0.75532387329617556;
      Ds[396] = -0.07800946004174078;
      Ds[397] = -0.98935225342006894;
      Ds[398] = -0.078009460041689779;
      Ds[399] = -0.7553238732961246;
      Ds[400] = 0.15601892008220355;
      Ds[401] = 0.0;
      Ds[402] = 0.0;
      Ds[403] = -4.5034278605725477;
      Ds[404] = 4.5034278605755409;
      Ds[405] = 0.83333333333786541;
      Ds[406] = 0.0;
      Ds[407] = 0.0;
      Ds[408] = 0.0;
      Ds[409] = 0.15601892008203644;
      Ds[410] = -0.75532387329528672;
      Ds[411] = -0.078009460040813258;
      Ds[412] = -0.98935225341951893;
      Ds[413] = -0.078009460042195777;
      Ds[414] = -0.75532387329666928;
      Ds[415] = 0.0;
      Ds[416] = 0.0;
      Ds[417] = -9.0068557211335047;
      Ds[418] = -4.50342786060667;
      Ds[419] = 0.8333333333374825;
      Ds[420] = 0.0;
      Ds[421] = 0.0;
      Ds[422] = 0.0;
      Ds[423] = -0.75532387329585826;
      Ds[424] = 0.15601892008277979;
      Ds[425] = -0.75532387329548278;
      Ds[426] = -0.0780094600412411;
      Ds[427] = -0.9893522534198792;
      Ds[428] = -0.078009460041616574;
      Ds[429] = 0.0;
      Ds[430] = 0.0;
      Ds[431] = -4.503427860568757;
      Ds[432] = -9.0068557211592442;
      Ds[433] = 0.83333333333709947;
      Ds[434] = 0.0;
      Ds[435] = 0.0;
      Ds[436] = 0.0;
      Ds[437] = -0.078009460041477088;
      Ds[438] = -0.75532387329487816;
      Ds[439] = 0.15601892008362025;
      Ds[440] = -0.75532387329562245;
      Ds[441] = -0.078009460042221493;
      Ds[442] = -0.98935225342071986;
      Ds[443] = 0.0;
      Ds[444] = 0.0;
      Ds[445] = 4.5034278605823763;
      Ds[446] = -4.5034278606254095;
      Ds[447] = 0.83333333333709958;
      Ds[448] = 0.0;
      Ds[449] = 0.0;
      Ds[450] = 0.0;
      Ds[451] = -0.23402838012351351;
      Ds[452] = 4.0861583006090994E-13;
      Ds[453] = 0.23402838012443353;
      Ds[454] = 0.23402838012389654;
      Ds[455] = -2.5712765250318615E-14;
      Ds[456] = -0.23402838012405056;
      Ds[457] = 0.0;
      Ds[458] = 0.0;
      Ds[459] = 13.51028358171588;
      Ds[460] = -1.873952584219296E-11;
      Ds[461] = -3.82950027244533E-13;
      Ds[462] = 0.0;
      Ds[463] = 0.0;
      Ds[464] = 0.0;
      Ds[465] = -3.1731297589487895E-13;
      Ds[466] = -0.23402838012452054;
      Ds[467] = -0.23402838012458602;
      Ds[468] = -4.4865551274142492E-13;
      Ds[469] = 0.2340283801237546;
      Ds[470] = 0.23402838012382013;
      Ds[471] = 0.0;
      Ds[472] = 0.0;
      Ds[473] = -3.7895574678259437E-12;
      Ds[474] = 13.510283581734784;
      Ds[475] = 7.6595423337266763E-13;
      Ds[476] = 0.0;
      Ds[477] = 0.0;
      Ds[478] = 0.0;
      Ds[479] = 0.23402838012383084;
      Ds[480] = 0.2340283801241119;
      Ds[481] = 1.5253380305082106E-13;
      Ds[482] = -0.23402838012344779;
      Ds[483] = -0.2340283801237289;
      Ds[484] = 2.3045751712408759E-13;
      Ds[485] = 0.0;
      Ds[486] = 0.0;
      Ds[487] = -13.510283581712089;
      Ds[488] = -13.510283581716045;
      Ds[489] = -3.8300128997766975E-13;
      Ds[490] = 0.0;
      Ds[491] = 0.0;
      Ds[492] = 0.0;
      Ds[493] = 0.8333333333376528;
      Ds[494] = 0.833333333336619;
      Ds[495] = 0.83333333333644866;
      Ds[496] = 0.83333333333731219;
      Ds[497] = 0.83333333333834614;
      Ds[498] = 0.83333333333851634;
      Ds[499] = 0.0;
      Ds[500] = 0.0;
      Ds[501] = -9.8280878546574965E-12;
      Ds[502] = 4.9867838924910758E-11;
      Ds[503] = -1.666666666674965;
      Ds[504] = 0.0;
      Ds[505] = 0.0;
      Ds[506] = 0.0;
      Ds[507] = 1.0;
      Ds[508] = 0.0;
      Ds[509] = 0.0;
      Ds[510] = 0.0;
      Ds[511] = 0.0;
      Ds[512] = 0.0;
      Ds[513] = 0.0;
      Ds[514] = 0.0;
      Ds[515] = 0.0;
      Ds[516] = 0.0;
      Ds[517] = 0.0;
      Ds[518] = 0.0;
      Ds[519] = 0.0;
      Ds[520] = 0.0;
      Ds[521] = 0.0;
      Ds[522] = 1.0;
      Ds[523] = 0.0;
      Ds[524] = 0.0;
      Ds[525] = 0.0;
      Ds[526] = 0.0;
      Ds[527] = 0.0;
      Ds[528] = 0.0;
      Ds[529] = 0.0;
      Ds[530] = 0.0;
      Ds[531] = 0.0;
      Ds[532] = 0.0;
      Ds[533] = 0.0;
      Ds[534] = 0.0;
      Ds[535] = 0.0;
      Ds[536] = 0.0;
      Ds[537] = 1.0;
      Ds[538] = 0.0;
      Ds[539] = 0.0;
      Ds[540] = 0.0;
      Ds[541] = 0.0;
      Ds[542] = 0.0;
      Ds[543] = 0.0;
      Ds[544] = 0.0;
      Ds[545] = 0.0;
      Ds[546] = 0.0;
      Ds[547] = 0.0;
      Ds[548] = 0.0;
      Ds[549] = 0.0;
      Ds[550] = 0.0;
      Ds[551] = 0.0;
      Ds[552] = 1.0;
      Ds[553] = 0.0;
      Ds[554] = 0.0;
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
      Ds[565] = 0.0;
      Ds[566] = 0.0;
      Ds[567] = 1.0;
      Ds[568] = 0.0;
      Ds[569] = 0.0;
      Ds[570] = 0.0;
      Ds[571] = 0.0;
      Ds[572] = 0.0;
      Ds[573] = 0.0;
      Ds[574] = 0.0;
      Ds[575] = 0.0;
      Ds[576] = 0.0;
      Ds[577] = 0.0;
      Ds[578] = 0.0;
      Ds[579] = 0.0;
      Ds[580] = 0.0;
      Ds[581] = 0.0;
      Ds[582] = 1.0;
      Ds[583] = 0.0;
      Ds[584] = 0.0;
      Ds[585] = 0.0;
      Ds[586] = 0.0;
      Ds[587] = 0.0;

      {
        int_T i1;
        for (i1=0; i1 < 42; i1++) {
          Chopper[i1] = 1;
        }
      }

      {
        /* Switches work vectors */
        int_T *switch_status = (int_T*)
          Diesel_DER2_DW.StateSpace_PWORK.SWITCH_STATUS;
        int_T *gState = (int_T*)Diesel_DER2_DW.StateSpace_PWORK.G_STATE;
        real_T *yswitch = (real_T*)Diesel_DER2_DW.StateSpace_PWORK.Y_SWITCH;
        int_T *switchTypes = (int_T*)
          Diesel_DER2_DW.StateSpace_PWORK.SWITCH_TYPES;
        int_T *idxOutSw = (int_T*)Diesel_DER2_DW.StateSpace_PWORK.IDX_OUT_SW;
        int_T *switch_status_init = (int_T*)
          Diesel_DER2_DW.StateSpace_PWORK.SWITCH_STATUS_INIT;

        /* Initialize work vectors */
        switch_status[0] = 0;
        switch_status_init[0] = 0;
        gState[0] = (int_T) 0.0;
        yswitch[0] = 1/0.001;
        switchTypes[0] = (int_T)2.0;
        idxOutSw[0] = ((int_T)0.0) - 1;
        switch_status[1] = 0;
        switch_status_init[1] = 0;
        gState[1] = (int_T) 0.0;
        yswitch[1] = 1/0.001;
        switchTypes[1] = (int_T)2.0;
        idxOutSw[1] = ((int_T)0.0) - 1;
        switch_status[2] = 0;
        switch_status_init[2] = 0;
        gState[2] = (int_T) 0.0;
        yswitch[2] = 1/0.001;
        switchTypes[2] = (int_T)2.0;
        idxOutSw[2] = ((int_T)0.0) - 1;
        switch_status[3] = 0;
        switch_status_init[3] = 0;
        gState[3] = (int_T) 0.0;
        yswitch[3] = 1/1.0E-5;
        switchTypes[3] = (int_T)3.0;
        idxOutSw[3] = ((int_T)37.0) - 1;
        switch_status[4] = 0;
        switch_status_init[4] = 0;
        gState[4] = (int_T) 0.0;
        yswitch[4] = 1/1.0E-5;
        switchTypes[4] = (int_T)3.0;
        idxOutSw[4] = ((int_T)38.0) - 1;
        switch_status[5] = 0;
        switch_status_init[5] = 0;
        gState[5] = (int_T) 0.0;
        yswitch[5] = 1/1.0E-5;
        switchTypes[5] = (int_T)3.0;
        idxOutSw[5] = ((int_T)39.0) - 1;
        switch_status[6] = 0;
        switch_status_init[6] = 0;
        gState[6] = (int_T) 0.0;
        yswitch[6] = 1/1.0E-5;
        switchTypes[6] = (int_T)3.0;
        idxOutSw[6] = ((int_T)40.0) - 1;
        switch_status[7] = 0;
        switch_status_init[7] = 0;
        gState[7] = (int_T) 0.0;
        yswitch[7] = 1/1.0E-5;
        switchTypes[7] = (int_T)3.0;
        idxOutSw[7] = ((int_T)41.0) - 1;
        switch_status[8] = 0;
        switch_status_init[8] = 0;
        gState[8] = (int_T) 0.0;
        yswitch[8] = 1/1.0E-5;
        switchTypes[8] = (int_T)3.0;
        idxOutSw[8] = ((int_T)42.0) - 1;
      }
    }

    /* InitializeConditions for UnitDelay: '<S141>/dw_delay' */
    Diesel_DER2_DW.dw_delay_DSTATE = Diesel_DER2_P.dw_delay_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S141>/dw_predict' */
    Diesel_DER2_DW.dw_predict_DSTATE = Diesel_DER2_P.dw_predict_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S28>/Delay_x' */
    Diesel_DER2_DW.Delay_x_DSTATE = Diesel_DER2_P.Delay_x_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S29>/Delay_x' */
    Diesel_DER2_DW.Delay_x_DSTATE_i = Diesel_DER2_P.Delay_x_InitialCondition_i;

    /* InitializeConditions for UnitDelay: '<S30>/Delay_x' */
    Diesel_DER2_DW.Delay_x_DSTATE_g = Diesel_DER2_P.Delay_x_InitialCondition_f;

    /* InitializeConditions for DiscreteIntegrator: '<S2>/Discrete-Time Integrator' */
    Diesel_DER2_DW.DiscreteTimeIntegrator_DSTATE =
      Diesel_DER2_P.DieselEngineSpeedRegulator_Pm0;

    /* InitializeConditions for UnitDelay: '<S33>/Unit Delay4' */
    Diesel_DER2_DW.UnitDelay4_DSTATE = Diesel_DER2_P.UnitDelay4_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S36>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE = Diesel_DER2_P.UnitDelay1_InitialCondition;

    /* InitializeConditions for DiscreteIntegrator: '<S85>/Discrete-Time Integrator' */
    Diesel_DER2_DW.DiscreteTimeIntegrator_DSTATE_f =
      Diesel_DER2_P.DiscretePIController_Init;

    /* InitializeConditions for DiscreteIntegrator: '<S97>/Integ4' */
    Diesel_DER2_DW.Integ4_DSTATE = Diesel_DER2_P.Integ4_IC;

    /* Level2 S-Function Block: '<S98>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[0];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S97>/Unit Delay' */
    Diesel_DER2_DW.UnitDelay_DSTATE = Diesel_DER2_P.UnitDelay_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S97>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE_c =
      Diesel_DER2_P.UnitDelay1_InitialCondition_l;

    /* InitializeConditions for DiscreteIntegrator: '<S95>/Integ4' */
    Diesel_DER2_DW.Integ4_DSTATE_e = Diesel_DER2_P.Integ4_IC_k;

    /* Level2 S-Function Block: '<S96>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[1];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S95>/Unit Delay' */
    Diesel_DER2_DW.UnitDelay_DSTATE_c =
      Diesel_DER2_P.UnitDelay_InitialCondition_b;

    /* InitializeConditions for UnitDelay: '<S95>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE_l =
      Diesel_DER2_P.UnitDelay1_InitialCondition_f;

    /* InitializeConditions for DiscreteIntegrator: '<S103>/Integ4' */
    Diesel_DER2_DW.Integ4_DSTATE_b = Diesel_DER2_P.Integ4_IC_l;

    /* Level2 S-Function Block: '<S104>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[2];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S103>/Unit Delay' */
    Diesel_DER2_DW.UnitDelay_DSTATE_j =
      Diesel_DER2_P.UnitDelay_InitialCondition_f;

    /* InitializeConditions for UnitDelay: '<S103>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE_e =
      Diesel_DER2_P.UnitDelay1_InitialCondition_h;

    /* InitializeConditions for DiscreteIntegrator: '<S101>/Integ4' */
    Diesel_DER2_DW.Integ4_DSTATE_f = Diesel_DER2_P.Integ4_IC_ly;

    /* Level2 S-Function Block: '<S102>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[3];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S101>/Unit Delay' */
    Diesel_DER2_DW.UnitDelay_DSTATE_l =
      Diesel_DER2_P.UnitDelay_InitialCondition_i;

    /* InitializeConditions for UnitDelay: '<S101>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE_ep =
      Diesel_DER2_P.UnitDelay1_InitialCondition_k;

    /* InitializeConditions for DiscreteIntegrator: '<S109>/Integ4' */
    Diesel_DER2_DW.Integ4_DSTATE_o = Diesel_DER2_P.Integ4_IC_o;

    /* Level2 S-Function Block: '<S110>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[4];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S109>/Unit Delay' */
    Diesel_DER2_DW.UnitDelay_DSTATE_h =
      Diesel_DER2_P.UnitDelay_InitialCondition_d;

    /* InitializeConditions for UnitDelay: '<S109>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE_k =
      Diesel_DER2_P.UnitDelay1_InitialCondition_f5;

    /* InitializeConditions for DiscreteIntegrator: '<S107>/Integ4' */
    Diesel_DER2_DW.Integ4_DSTATE_l = Diesel_DER2_P.Integ4_IC_a;

    /* Level2 S-Function Block: '<S108>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[5];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S107>/Unit Delay' */
    Diesel_DER2_DW.UnitDelay_DSTATE_ch =
      Diesel_DER2_P.UnitDelay_InitialCondition_k;

    /* InitializeConditions for UnitDelay: '<S107>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE_j =
      Diesel_DER2_P.UnitDelay1_InitialCondition_e;
    for (i = 0; i < 5; i++) {
      /* InitializeConditions for UnitDelay: '<S132>/fluxes' */
      Diesel_DER2_DW.fluxes_DSTATE[i] = Diesel_DER2_P.fluxes_InitialCondition[i];

      /* InitializeConditions for UnitDelay: '<S67>/fluxes' */
      Diesel_DER2_DW.fluxes_DSTATE_o[i] =
        Diesel_DER2_P.fluxes_InitialCondition_h[i];

      /* InitializeConditions for UnitDelay: '<S67>/voltages' */
      Diesel_DER2_DW.voltages_DSTATE[i] =
        Diesel_DER2_P.voltages_InitialCondition;

      /* InitializeConditions for UnitDelay: '<S132>/voltages' */
      Diesel_DER2_DW.voltages_DSTATE_h[i] =
        Diesel_DER2_P.voltages_InitialCondition_c;
    }

    /* InitializeConditions for DiscreteIntegrator: '<S119>/Rotor speed deviation (dw)' */
    Diesel_DER2_DW.Rotorspeeddeviationdw_DSTATE =
      Diesel_DER2_P.Rotorspeeddeviationdw_IC;

    /* InitializeConditions for DiscreteIntegrator: '<S111>/Integ4' */
    Diesel_DER2_DW.Integ4_DSTATE_m = Diesel_DER2_P.Integ4_IC_oa;

    /* Level2 S-Function Block: '<S112>/S-Function' (sfun_discreteVariableDelay) */
    {
      SimStruct *rts = Diesel_DER2_M->childSfunctions[6];
      sfcnInitializeConditions(rts);
      if (ssGetErrorStatus(rts) != (NULL))
        return;
    }

    /* InitializeConditions for UnitDelay: '<S111>/Unit Delay' */
    Diesel_DER2_DW.UnitDelay_DSTATE_a =
      Diesel_DER2_P.UnitDelay_InitialCondition_h;

    /* InitializeConditions for UnitDelay: '<S111>/Unit Delay1' */
    Diesel_DER2_DW.UnitDelay1_DSTATE_m =
      Diesel_DER2_P.UnitDelay1_InitialCondition_lb;

    /* InitializeConditions for UnitDelay: '<S146>/Delay Ts' */
    Diesel_DER2_DW.DelayTs_DSTATE = Diesel_DER2_P.DelayTs_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S146>/Delay Ts ' */
    Diesel_DER2_DW.DelayTs_DSTATE_f = Diesel_DER2_P.DelayTs_InitialCondition_c;

    /* InitializeConditions for UnitDelay: '<S149>/Delay Ts' */
    Diesel_DER2_DW.DelayTs_DSTATE_a = Diesel_DER2_P.DelayTs_InitialCondition_p;

    /* InitializeConditions for UnitDelay: '<S149>/Delay Ts ' */
    Diesel_DER2_DW.DelayTs_DSTATE_n = Diesel_DER2_P.DelayTs_InitialCondition_e;

    /* InitializeConditions for UnitDelay: '<S152>/Delay Ts' */
    Diesel_DER2_DW.DelayTs_DSTATE_l = Diesel_DER2_P.DelayTs_InitialCondition_d;

    /* InitializeConditions for UnitDelay: '<S152>/Delay Ts ' */
    Diesel_DER2_DW.DelayTs_DSTATE_d = Diesel_DER2_P.DelayTs_InitialCondition_p4;

    /* Enable for Sin: '<S87>/sin(wt)' */
    Diesel_DER2_DW.systemEnable = 1;

    /* Enable for DiscreteIntegrator: '<S97>/Integ4' */
    Diesel_DER2_DW.Integ4_SYSTEM_ENABLE = 1U;

    /* Enable for Sin: '<S87>/cos(wt)' */
    Diesel_DER2_DW.systemEnable_b = 1;

    /* Enable for DiscreteIntegrator: '<S95>/Integ4' */
    Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_f = 1U;

    /* Enable for Sin: '<S88>/sin(wt)' */
    Diesel_DER2_DW.systemEnable_a = 1;

    /* Enable for DiscreteIntegrator: '<S103>/Integ4' */
    Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_j = 1U;

    /* Enable for Sin: '<S88>/cos(wt)' */
    Diesel_DER2_DW.systemEnable_o = 1;

    /* Enable for DiscreteIntegrator: '<S101>/Integ4' */
    Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_fj = 1U;

    /* Enable for Sin: '<S89>/sin(wt)' */
    Diesel_DER2_DW.systemEnable_j = 1;

    /* Enable for DiscreteIntegrator: '<S109>/Integ4' */
    Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_h = 1U;

    /* Enable for Sin: '<S89>/cos(wt)' */
    Diesel_DER2_DW.systemEnable_p = 1;

    /* Enable for DiscreteIntegrator: '<S107>/Integ4' */
    Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_c = 1U;

    /* Enable for DiscreteIntegrator: '<S119>/Rotor speed deviation (dw)' */
    Diesel_DER2_DW.Rotorspeeddeviationdw_SYSTEM_EN = 1U;

    /* Enable for DiscreteIntegrator: '<S111>/Integ4' */
    Diesel_DER2_DW.Integ4_SYSTEM_ENABLE_d = 1U;
  }
}

/* Model terminate function */
void Diesel_DER2_terminate(void)
{
  /* S-Function block: <S156>/State-Space */
  {
    /* Free memory */
    free(Diesel_DER2_DW.StateSpace_PWORK.AS);
    free(Diesel_DER2_DW.StateSpace_PWORK.BS);
    free(Diesel_DER2_DW.StateSpace_PWORK.CS);
    free(Diesel_DER2_DW.StateSpace_PWORK.DS);
    free(Diesel_DER2_DW.StateSpace_PWORK.DX_COL);
    free(Diesel_DER2_DW.StateSpace_PWORK.TMP2);
    free(Diesel_DER2_DW.StateSpace_PWORK.BD_COL);
    free(Diesel_DER2_DW.StateSpace_PWORK.TMP1);
    free(Diesel_DER2_DW.StateSpace_PWORK.XTMP);

    /*
     * Circuit has switches*/
    free(Diesel_DER2_DW.StateSpace_PWORK.CHOPPER);
    free(Diesel_DER2_DW.StateSpace_PWORK.G_STATE);
    free(Diesel_DER2_DW.StateSpace_PWORK.SWITCH_STATUS);
    free(Diesel_DER2_DW.StateSpace_PWORK.SW_CHG);
    free(Diesel_DER2_DW.StateSpace_PWORK.SWITCH_STATUS_INIT);
  }

  /* Level2 S-Function Block: '<S98>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = Diesel_DER2_M->childSfunctions[0];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S96>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = Diesel_DER2_M->childSfunctions[1];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S104>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = Diesel_DER2_M->childSfunctions[2];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S102>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = Diesel_DER2_M->childSfunctions[3];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S110>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = Diesel_DER2_M->childSfunctions[4];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S108>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = Diesel_DER2_M->childSfunctions[5];
    sfcnTerminate(rts);
  }

  /* Level2 S-Function Block: '<S112>/S-Function' (sfun_discreteVariableDelay) */
  {
    SimStruct *rts = Diesel_DER2_M->childSfunctions[6];
    sfcnTerminate(rts);
  }
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
