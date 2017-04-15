/*
 * File: Wind_songweiwei_private.h
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

#ifndef RTW_HEADER_Wind_songweiwei_private_h_
#define RTW_HEADER_Wind_songweiwei_private_h_
#include "rtwtypes.h"
#include <math.h>
#include <stdlib.h>
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

#define MAX_DELAY_BUFFER_SIZE          32768
#ifndef ssGetFixedStepSize
#define ssGetFixedStepSize(S)          (S).stepSize
#endif                                 /* ssGetFixedStepSize */

extern real_T rt_roundd_snf(real_T u);
extern real_T rt_remd_snf(real_T u0, real_T u1);
extern void rt_invd5x5_snf(const real_T u[25], real_T y[25]);
extern real_T rt_modd_snf(real_T u0, real_T u1);
extern real_T rt_powd_snf(real_T u0, real_T u1);
extern real_T rt_hypotd_snf(real_T u0, real_T u1);
extern real_T rt_atan2d_snf(real_T u0, real_T u1);
void BINARYSEARCH_real_T(uint32_T *piLeft, uint32_T *piRght, real_T u, const
  real_T *pData, uint32_T iHi);
void LookUp_real_T_real_T(real_T *pY, const real_T *pYData, real_T u, const
  real_T *pUData, uint32_T iHi);
extern void sfun_discreteVariableDelay(SimStruct *rts);
extern void Wind_so_Subsystempi2delay_Start(B_Subsystempi2delay_Wind_song_T
  *localB, P_Subsystempi2delay_Wind_song_T *localP);
extern void Wind_songweiw_Subsystempi2delay(uint8_T rtu_Enable, const real_T
  rtu_alpha_beta[2], real_T rtu_wt, B_Subsystempi2delay_Wind_song_T *localB);
extern void Wind_songweiwe_Subsystem1_Start(B_Subsystem1_Wind_songweiwei_T
  *localB, P_Subsystem1_Wind_songweiwei_T *localP);
extern void Wind_songweiwei_Subsystem1(uint8_T rtu_Enable, const real_T
  rtu_alpha_beta[2], real_T rtu_wt, B_Subsystem1_Wind_songweiwei_T *localB);

#endif                                 /* RTW_HEADER_Wind_songweiwei_private_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
