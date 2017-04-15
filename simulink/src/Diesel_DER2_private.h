/*
 * File: Diesel_DER2_private.h
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

#ifndef RTW_HEADER_Diesel_DER2_private_h_
#define RTW_HEADER_Diesel_DER2_private_h_
#include "rtwtypes.h"
#include <math.h>
#include <stdlib.h>

/* Private macros used by the generated code to access rtModel */
#ifndef rtmIsMajorTimeStep
# define rtmIsMajorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
# define rtmIsMinorTimeStep(rtm)       (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTFinal
# define rtmSetTFinal(rtm, val)        ((rtm)->Timing.tFinal = (val))
#endif

#ifndef rtmGetTPtr
# define rtmGetTPtr(rtm)               ((rtm)->Timing.t)
#endif

#ifndef rtmSetTPtr
# define rtmSetTPtr(rtm, val)          ((rtm)->Timing.t = (val))
#endif

#ifndef CodeFormat
#define CodeFormat                     S-Function
#else
#undef CodeFormat
#define CodeFormat                     S-Function
#endif

#ifndef S_FUNCTION_NAME
#define S_FUNCTION_NAME                simulink_only_sfcn
#else
#undef S_FUNCTION_NAME
#define S_FUNCTION_NAME                simulink_only_sfcn
#endif

#ifndef S_FUNCTION_LEVEL
#define S_FUNCTION_LEVEL               2
#else
#undef S_FUNCTION_LEVEL
#define S_FUNCTION_LEVEL               2
#endif

#ifndef RTW_GENERATED_S_FUNCTION
#define RTW_GENERATED_S_FUNCTION
#endif

#ifndef rtmGetDataMapInfo
# define rtmGetDataMapInfo(rtm)        NULL
#endif

#ifndef rtmSetDataMapInfo
# define rtmSetDataMapInfo(rtm, val)
#endif

#if !defined(RTW_SFUNCTION_DEFINES)
#define RTW_SFUNCTION_DEFINES
#ifndef _RTW_COMMON_DEFINES_
#define _RTW_COMMON_DEFINES_
#endif
#endif

extern void rt_invd5x5_snf(const real_T u[25], real_T y[25]);
extern real_T rt_atan2d_snf(real_T u0, real_T u1);
extern real_T rt_hypotd_snf(real_T u0, real_T u1);
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
  ;
extern void sfun_discreteVariableDelay(SimStruct *rts);
extern void Diesel_DER2_Lmq_sat_Init(DW_Lmq_sat_Diesel_DER2_T *localDW,
  P_Lmq_sat_Diesel_DER2_T *localP);
extern void Diesel_DER2_Lmq_sat_Start(DW_Lmq_sat_Diesel_DER2_T *localDW,
  P_Lmq_sat_Diesel_DER2_T *localP);
extern void Diesel_DER2_Lmq_sat(boolean_T rtu_Enable, const real_T rtu_phi[5],
  B_Lmq_sat_Diesel_DER2_T *localB, DW_Lmq_sat_Diesel_DER2_T *localDW,
  P_Lmq_sat_Diesel_DER2_T *localP);
extern void Diesel_DER2_NegSeqComputation(real_T rtu_Enable, creal_T rtu_In,
  creal_T rtu_In_i, creal_T rtu_In_o, B_NegSeqComputation_Diesel_DE_T *localB,
  P_NegSeqComputation_Diesel_DE_T *localP);

#endif                                 /* RTW_HEADER_Diesel_DER2_private_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
