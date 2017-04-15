/*
 * File: Diesel_DER2.h
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

#ifndef RTW_HEADER_Diesel_DER2_h_
#define RTW_HEADER_Diesel_DER2_h_
#include <stddef.h>
#include <float.h>
#include <math.h>
#include <string.h>
#ifndef Diesel_DER2_COMMON_INCLUDES_
# define Diesel_DER2_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#endif                                 /* Diesel_DER2_COMMON_INCLUDES_ */

#include "Diesel_DER2_types.h"
#include "rt_defines.h"
#include "rt_look.h"
#include "rt_look1d.h"
#include "rt_nonfinite.h"
#include "rtGetInf.h"
#include "sfcn_bridge.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetFinalTime
# define rtmGetFinalTime(rtm)          ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetSampleHitArray
# define rtmGetSampleHitArray(rtm)     ((rtm)->Timing.sampleHitArray)
#endif

#ifndef rtmGetStepSize
# define rtmGetStepSize(rtm)           ((rtm)->Timing.stepSize)
#endif

#ifndef rtmGet_TimeOfLastOutput
# define rtmGet_TimeOfLastOutput(rtm)  ((rtm)->Timing.timeOfLastOutput)
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTFinal
# define rtmGetTFinal(rtm)             ((rtm)->Timing.tFinal)
#endif

#ifndef rtmGetTStart
# define rtmGetTStart(rtm)             ((rtm)->Timing.tStart)
#endif

#ifndef rtmGetTimeOfLastOutput
# define rtmGetTimeOfLastOutput(rtm)   ((rtm)->Timing.timeOfLastOutput)
#endif

/* Block signals for system '<S66>/Lmq_sat' */
typedef struct {
  real_T MathFunction1;                /* '<S73>/Math Function1' */
  real_T Lmsatq;                       /* '<S69>/Lmq' */
} B_Lmq_sat_Diesel_DER2_T;

/* Block states (auto storage) for system '<S66>/Lmq_sat' */
typedef struct {
  real_T Lmq_sat_DSTATE;               /* '<S69>/Lmq_sat' */
} DW_Lmq_sat_Diesel_DER2_T;

/* Block signals for system '<S86>/Neg. Seq. Computation' */
typedef struct {
  creal_T Gain3;                       /* '<S90>/Gain3' */
} B_NegSeqComputation_Diesel_DE_T;

/* Block signals (auto storage) */
typedef struct {
  creal_T Gain3;                       /* '<S92>/Gain3' */
  real_T ib[3];                        /* '<S123>/ib' */
  real_T ib_g[3];                      /* '<S58>/ib' */
  real_T UnitDelay2;                   /* '<S33>/Unit Delay2' */
  real_T StateSpace_o1[42];            /* '<S156>/State-Space' */
  real_T StateSpace_o2[9];             /* '<S156>/State-Space' */
  real_T Gazflow;                      /* '<S2>/Discrete-Time Integrator' */
  real_T Integ4;                       /* '<S97>/Integ4' */
  real_T K1;                           /* '<S97>/K1' */
  real_T SFunction;                    /* '<S98>/S-Function' */
  real_T Integ4_h;                     /* '<S95>/Integ4' */
  real_T K1_e;                         /* '<S95>/K1' */
  real_T SFunction_a;                  /* '<S96>/S-Function' */
  real_T Integ4_m;                     /* '<S103>/Integ4' */
  real_T K1_c;                         /* '<S103>/K1' */
  real_T SFunction_e;                  /* '<S104>/S-Function' */
  real_T Integ4_o;                     /* '<S101>/Integ4' */
  real_T K1_j;                         /* '<S101>/K1' */
  real_T SFunction_i;                  /* '<S102>/S-Function' */
  real_T Integ4_i;                     /* '<S109>/Integ4' */
  real_T K1_k;                         /* '<S109>/K1' */
  real_T SFunction_g;                  /* '<S110>/S-Function' */
  real_T Integ4_k;                     /* '<S107>/Integ4' */
  real_T K1_o;                         /* '<S107>/K1' */
  real_T SFunction_p;                  /* '<S108>/S-Function' */
  real_T dw;                           /* '<S119>/Rotor speed deviation (dw)' */
  real_T Integ4_a;                     /* '<S111>/Integ4' */
  real_T K1_m;                         /* '<S111>/K1' */
  real_T SFunction_iq;                 /* '<S112>/S-Function' */
  real_T IC;                           /* '<S146>/IC' */
  real_T IC_a;                         /* '<S149>/IC' */
  real_T IC_c;                         /* '<S152>/IC' */
  real_T MathFunction3;                /* '<S136>/Math Function3' */
  real_T Lmsatd;                       /* '<S133>/Lmd' */
  real_T Switch1;                      /* '<S131>/Switch1' */
  real_T Linv[25];                     /* '<S131>/inversion' */
  real_T RLinv[25];                    /* '<S131>/Product1' */
  real_T Switch2;                      /* '<S131>/Switch2' */
  real_T Linv_a[25];                   /* '<S66>/inversion' */
  real_T RLinv_e[25];                  /* '<S66>/Product1' */
  B_Lmq_sat_Diesel_DER2_T Lmq_sat_l;   /* '<S131>/Lmq_sat' */
  B_NegSeqComputation_Diesel_DE_T PosSeqComputation;/* '<S86>/Pos. Seq. Computation' */
  B_NegSeqComputation_Diesel_DE_T NegSeqComputation;/* '<S86>/Neg. Seq. Computation' */
  B_Lmq_sat_Diesel_DER2_T Lmq_sat;     /* '<S66>/Lmq_sat' */
} B_Diesel_DER2_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T Rotorangledthetae_DSTATE;     /* '<S119>/Rotor angle dthetae' */
  real_T fluxes_DSTATE[5];             /* '<S132>/fluxes' */
  real_T Rotorangledtheta_DSTATE;      /* '<S54>/Rotor angle dtheta' */
  real_T fluxes_DSTATE_o[5];           /* '<S67>/fluxes' */
  real_T UnitDelay2_DSTATE;            /* '<S33>/Unit Delay2' */
  real_T StateSpace_DSTATE[6];         /* '<S156>/State-Space' */
  real_T dw_delay_DSTATE;              /* '<S141>/dw_delay' */
  real_T dw_predict_DSTATE;            /* '<S141>/dw_predict' */
  real_T Delay_x_DSTATE;               /* '<S28>/Delay_x' */
  real_T Delay_x_DSTATE_i;             /* '<S29>/Delay_x' */
  real_T Delay_x_DSTATE_g;             /* '<S30>/Delay_x' */
  real_T DiscreteTimeIntegrator_DSTATE;/* '<S2>/Discrete-Time Integrator' */
  real_T UnitDelay4_DSTATE;            /* '<S33>/Unit Delay4' */
  real_T UnitDelay1_DSTATE;            /* '<S36>/Unit Delay1' */
  real_T DiscreteTimeIntegrator_DSTATE_f;/* '<S85>/Discrete-Time Integrator' */
  real_T voltages_DSTATE[5];           /* '<S67>/voltages' */
  real_T Integ4_DSTATE;                /* '<S97>/Integ4' */
  real_T UnitDelay_DSTATE;             /* '<S97>/Unit Delay' */
  real_T UnitDelay1_DSTATE_c;          /* '<S97>/Unit Delay1' */
  real_T Integ4_DSTATE_e;              /* '<S95>/Integ4' */
  real_T UnitDelay_DSTATE_c;           /* '<S95>/Unit Delay' */
  real_T UnitDelay1_DSTATE_l;          /* '<S95>/Unit Delay1' */
  real_T Integ4_DSTATE_b;              /* '<S103>/Integ4' */
  real_T UnitDelay_DSTATE_j;           /* '<S103>/Unit Delay' */
  real_T UnitDelay1_DSTATE_e;          /* '<S103>/Unit Delay1' */
  real_T Integ4_DSTATE_f;              /* '<S101>/Integ4' */
  real_T UnitDelay_DSTATE_l;           /* '<S101>/Unit Delay' */
  real_T UnitDelay1_DSTATE_ep;         /* '<S101>/Unit Delay1' */
  real_T Integ4_DSTATE_o;              /* '<S109>/Integ4' */
  real_T UnitDelay_DSTATE_h;           /* '<S109>/Unit Delay' */
  real_T UnitDelay1_DSTATE_k;          /* '<S109>/Unit Delay1' */
  real_T Integ4_DSTATE_l;              /* '<S107>/Integ4' */
  real_T UnitDelay_DSTATE_ch;          /* '<S107>/Unit Delay' */
  real_T UnitDelay1_DSTATE_j;          /* '<S107>/Unit Delay1' */
  real_T voltages_DSTATE_h[5];         /* '<S132>/voltages' */
  real_T Rotorspeeddeviationdw_DSTATE; /* '<S119>/Rotor speed deviation (dw)' */
  real_T Integ4_DSTATE_m;              /* '<S111>/Integ4' */
  real_T UnitDelay_DSTATE_a;           /* '<S111>/Unit Delay' */
  real_T UnitDelay1_DSTATE_m;          /* '<S111>/Unit Delay1' */
  real_T DelayTs_DSTATE;               /* '<S146>/Delay Ts' */
  real_T DelayTs_DSTATE_f;             /* '<S146>/Delay Ts ' */
  real_T DelayTs_DSTATE_a;             /* '<S149>/Delay Ts' */
  real_T DelayTs_DSTATE_n;             /* '<S149>/Delay Ts ' */
  real_T DelayTs_DSTATE_l;             /* '<S152>/Delay Ts' */
  real_T DelayTs_DSTATE_d;             /* '<S152>/Delay Ts ' */
  real_T Lmd_sat_DSTATE;               /* '<S133>/Lmd_sat' */
  real_T Lmd_sat_DSTATE_h;             /* '<S68>/Lmd_sat' */
  real_T lastSin;                      /* '<S87>/sin(wt)' */
  real_T lastCos;                      /* '<S87>/sin(wt)' */
  real_T lastSin_d;                    /* '<S87>/cos(wt)' */
  real_T lastCos_e;                    /* '<S87>/cos(wt)' */
  real_T lastSin_a;                    /* '<S88>/sin(wt)' */
  real_T lastCos_c;                    /* '<S88>/sin(wt)' */
  real_T lastSin_f;                    /* '<S88>/cos(wt)' */
  real_T lastCos_eg;                   /* '<S88>/cos(wt)' */
  real_T lastSin_o;                    /* '<S89>/sin(wt)' */
  real_T lastCos_h;                    /* '<S89>/sin(wt)' */
  real_T lastSin_i;                    /* '<S89>/cos(wt)' */
  real_T lastCos_g;                    /* '<S89>/cos(wt)' */
  real_T inversion_DWORK4[25];         /* '<S140>/inversion' */
  real_T inversion_DWORK4_i[25];       /* '<S131>/inversion' */
  real_T inversion_DWORK4_k[25];       /* '<S75>/inversion' */
  real_T inversion_DWORK4_c[25];       /* '<S66>/inversion' */
  struct {
    real_T modelTStart;
    real_T TUbufferArea[2048];
  } DelayTd_RWORK;                     /* '<S2>/Delay Td' */

  real_T SFunction_RWORK;              /* '<S98>/S-Function' */
  real_T SFunction_RWORK_k;            /* '<S96>/S-Function' */
  real_T SFunction_RWORK_m;            /* '<S104>/S-Function' */
  real_T SFunction_RWORK_j;            /* '<S102>/S-Function' */
  real_T SFunction_RWORK_b;            /* '<S110>/S-Function' */
  real_T SFunction_RWORK_g;            /* '<S108>/S-Function' */
  real_T SFunction_RWORK_g5;           /* '<S112>/S-Function' */
  struct {
    void *AS;
    void *BS;
    void *CS;
    void *DS;
    void *DX_COL;
    void *BD_COL;
    void *TMP1;
    void *TMP2;
    void *XTMP;
    void *SWITCH_STATUS;
    void *SWITCH_STATUS_INIT;
    void *SW_CHG;
    void *CHOPPER;
    void *G_STATE;
    void *XKM12;
    void *XKP12;
    void *XLAST;
    void *ULAST;
    void *IDX_SW_CHG;
    void *Y_SWITCH;
    void *SWITCH_TYPES;
    void *IDX_OUT_SW;
  } StateSpace_PWORK;                  /* '<S156>/State-Space' */

  struct {
    void *TUbufferPtrs[2];
  } DelayTd_PWORK;                     /* '<S2>/Delay Td' */

  void *SFunction_PWORK;               /* '<S98>/S-Function' */
  void *SFunction_PWORK_c;             /* '<S96>/S-Function' */
  void *SFunction_PWORK_m;             /* '<S104>/S-Function' */
  void *SFunction_PWORK_a;             /* '<S102>/S-Function' */
  void *SFunction_PWORK_ao;            /* '<S110>/S-Function' */
  void *SFunction_PWORK_i;             /* '<S108>/S-Function' */
  void *SFunction_PWORK_h;             /* '<S112>/S-Function' */
  int32_T systemEnable;                /* '<S87>/sin(wt)' */
  int32_T systemEnable_b;              /* '<S87>/cos(wt)' */
  int32_T systemEnable_a;              /* '<S88>/sin(wt)' */
  int32_T systemEnable_o;              /* '<S88>/cos(wt)' */
  int32_T systemEnable_j;              /* '<S89>/sin(wt)' */
  int32_T systemEnable_p;              /* '<S89>/cos(wt)' */
  int_T StateSpace_IWORK[5];           /* '<S156>/State-Space' */
  struct {
    int_T Tail;
    int_T Head;
    int_T Last;
    int_T CircularBufSize;
  } DelayTd_IWORK;                     /* '<S2>/Delay Td' */

  int_T SFunction_IWORK;               /* '<S98>/S-Function' */
  int_T SFunction_IWORK_e;             /* '<S96>/S-Function' */
  int_T SFunction_IWORK_m;             /* '<S104>/S-Function' */
  int_T SFunction_IWORK_k;             /* '<S102>/S-Function' */
  int_T SFunction_IWORK_n;             /* '<S110>/S-Function' */
  int_T SFunction_IWORK_d;             /* '<S108>/S-Function' */
  int_T SFunction_IWORK_i;             /* '<S112>/S-Function' */
  uint8_T Integ4_SYSTEM_ENABLE;        /* '<S97>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_f;      /* '<S95>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_j;      /* '<S103>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_fj;     /* '<S101>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_h;      /* '<S109>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_c;      /* '<S107>/Integ4' */
  uint8_T Rotorspeeddeviationdw_SYSTEM_EN;/* '<S119>/Rotor speed deviation (dw)' */
  uint8_T Integ4_SYSTEM_ENABLE_d;      /* '<S111>/Integ4' */
  boolean_T u5_Mode;                   /* '<S146>/>1.5' */
  boolean_T IC_FirstOutputTime;        /* '<S146>/IC' */
  boolean_T u5_Mode_n;                 /* '<S149>/>1.5' */
  boolean_T IC_FirstOutputTime_l;      /* '<S149>/IC' */
  boolean_T u5_Mode_i;                 /* '<S152>/>1.5' */
  boolean_T IC_FirstOutputTime_n;      /* '<S152>/IC' */
  DW_Lmq_sat_Diesel_DER2_T Lmq_sat_l;  /* '<S131>/Lmq_sat' */
  DW_Lmq_sat_Diesel_DER2_T Lmq_sat;    /* '<S66>/Lmq_sat' */
} DW_Diesel_DER2_T;

/* External outputs (root outports fed by signals with auto storage) */
typedef struct {
  real_T out_If;                       /* '<Root>/out_If' */
  real_T out_Vf[2];                    /* '<Root>/out_Vf' */
  real_T out_Vt_SM;                    /* '<Root>/out_Vt_SM' */
  real_T out_Pmec[2];                  /* '<Root>/out_Pmec' */
  real_T out_w_pu;                     /* '<Root>/out_w_pu' */
  real_T out_Vabc_SM[3];               /* '<Root>/out_Vabc_SM' */
  real_T out_Iabc_SM[3];               /* '<Root>/out_Iabc_SM' */
  real_T out_Vabc_Exc[3];              /* '<Root>/out_Vabc_Exc' */
  real_T out_Iabc_Exc[3];              /* '<Root>/out_Iabc_Exc' */
  real_T out_VD1;                      /* '<Root>/out_VD1' */
  real_T out_I;                        /* '<Root>/out_I' */
} ExtY_Diesel_DER2_T;

/* Parameters for system: '<S66>/Lmq_sat' */
struct P_Lmq_sat_Diesel_DER2_T_ {
  real_T Constant1_Value;              /* Expression: SM.Lmsatq(1)
                                        * Referenced by: '<S69>/Constant1'
                                        */
  real_T Ll_q_Gain[2];                 /* Expression: SM.One_Llq
                                        * Referenced by: '<S74>/1//Ll_q'
                                        */
  real_T u2_Value[3];                  /* Expression: [ 1/SM.Ll 1/SM.Llkq1 1/SM.Llkq2]
                                        * Referenced by: '<S73>/u2'
                                        */
  real_T Lmq_sat_InitialCondition;     /* Expression: SM.Lmsatq(1)
                                        * Referenced by: '<S69>/Lmq_sat'
                                        */
  real_T LookupTable_XData[2];         /* Expression: SM.Phisat
                                        * Referenced by: '<S69>/Lookup Table'
                                        */
  real_T LookupTable_YData[2];         /* Expression: [ 0 SM.Phisat(2:end)./SM.Lmsatq(2:end)*SM.Lmq ]
                                        * Referenced by: '<S69>/Lookup Table'
                                        */
  real_T Lmq_Gain;                     /* Expression: SM.Lmq
                                        * Referenced by: '<S69>/Lmq'
                                        */
};

/* Parameters for system: '<S86>/Neg. Seq. Computation' */
struct P_NegSeqComputation_Diesel_DE_T_ {
  real_T Gain3_Gain;                   /* Expression: 1/3
                                        * Referenced by: '<S90>/Gain3'
                                        */
  creal_T Gain1_Gain[3];               /* Expression: [1 sps.a2 sps.a]
                                        * Referenced by: '<S90>/Gain1'
                                        */
};

/* Parameters (auto storage) */
struct P_Diesel_DER2_T_ {
  real_T Fourier_A_Freq;               /* Mask Parameter: Fourier_A_Freq
                                        * Referenced by:
                                        *   '<S87>/cos(wt)'
                                        *   '<S87>/sin(wt)'
                                        */
  real_T Fourier_B_Freq;               /* Mask Parameter: Fourier_B_Freq
                                        * Referenced by:
                                        *   '<S88>/cos(wt)'
                                        *   '<S88>/sin(wt)'
                                        */
  real_T Fourier_C_Freq;               /* Mask Parameter: Fourier_C_Freq
                                        * Referenced by:
                                        *   '<S89>/cos(wt)'
                                        *   '<S89>/sin(wt)'
                                        */
  real_T DiscretePIController_Init;    /* Mask Parameter: DiscretePIController_Init
                                        * Referenced by: '<S85>/Discrete-Time Integrator'
                                        */
  real_T DieselEngineSpeedRegulator_K; /* Mask Parameter: DieselEngineSpeedRegulator_K
                                        * Referenced by: '<S2>/Gain K'
                                        */
  real_T DiscretePIController_Ki;      /* Mask Parameter: DiscretePIController_Ki
                                        * Referenced by: '<S85>/Kp5'
                                        */
  real_T DiscretePIController_Kp;      /* Mask Parameter: DiscretePIController_Kp
                                        * Referenced by: '<S85>/Kp4'
                                        */
  real_T DieselEngineSpeedRegulator_Pm0;/* Mask Parameter: DieselEngineSpeedRegulator_Pm0
                                         * Referenced by:
                                         *   '<S2>/Discrete-Time Integrator'
                                         *   '<S2>/Delay Td'
                                         */
  real_T ThreePhaseBreaker_SwitchA;    /* Mask Parameter: ThreePhaseBreaker_SwitchA
                                        * Referenced by: '<S8>/Constant1'
                                        */
  real_T ThreePhaseBreaker_SwitchB;    /* Mask Parameter: ThreePhaseBreaker_SwitchB
                                        * Referenced by: '<S8>/Constant2'
                                        */
  real_T ThreePhaseBreaker_SwitchC;    /* Mask Parameter: ThreePhaseBreaker_SwitchC
                                        * Referenced by: '<S8>/Constant3'
                                        */
  real_T DieselEngineSpeedRegulator_Td;/* Mask Parameter: DieselEngineSpeedRegulator_Td
                                        * Referenced by: '<S2>/Delay Td'
                                        */
  real_T Constant4_Value[25];          /* Expression: SM.Linv
                                        * Referenced by: '<S56>/Constant4'
                                        */
  real_T u3_Value;                     /* Expression: SM.Lmq
                                        * Referenced by: '<S66>/u3'
                                        */
  real_T Constant1_Value;              /* Expression: SM.Lmsatd(1)
                                        * Referenced by: '<S68>/Constant1'
                                        */
  real_T Ll_d_Gain[3];                 /* Expression: [ 1/SM.Ll   1/SM.Llfd   1/SM.Llkd ]
                                        * Referenced by: '<S72>/1//Ll_d'
                                        */
  real_T u1_Value[3];                  /* Expression: [1/SM.Ll 1/SM.Llkd 1/SM.Llfd]
                                        * Referenced by: '<S71>/u1'
                                        */
  real_T Lmd_sat_InitialCondition;     /* Expression: SM.Lmsatd(1)
                                        * Referenced by: '<S68>/Lmd_sat'
                                        */
  real_T LookupTable_XData[2];         /* Expression: SM.Phisat
                                        * Referenced by: '<S68>/Lookup Table'
                                        */
  real_T LookupTable_YData[2];         /* Expression: [ 0 SM.Phisat(2:end)./SM.Lmsatd(2:end)*SM.Lmd ]
                                        * Referenced by: '<S68>/Lookup Table'
                                        */
  real_T Lmd_Gain;                     /* Expression: SM.Lmd
                                        * Referenced by: '<S68>/Lmd'
                                        */
  real_T u1_Value_c[25];               /* Expression: SM.R
                                        * Referenced by: '<S66>/u1'
                                        */
  real_T u1_Value_p[25];               /* Expression: zeros(SM.nState,SM.nState)
                                        * Referenced by: '<S70>/u1'
                                        */
  real_T u5_Value[25];                 /* Expression: SM.Llqd
                                        * Referenced by: '<S70>/u5'
                                        */
  real_T Constant6_Value[25];          /* Expression: SM.RLinv
                                        * Referenced by: '<S56>/Constant6'
                                        */
  real_T Constant2_Value;              /* Expression: SM.Sat
                                        * Referenced by: '<S56>/Constant2'
                                        */
  real_T Switch1_Threshold;            /* Expression: 0.5
                                        * Referenced by: '<S56>/Switch1'
                                        */
  real_T u1_Value_e[25];               /* Expression: zeros(SM.nState,SM.nState)
                                        * Referenced by: '<S64>/u1'
                                        */
  real_T Gain1_Gain;                   /* Expression: -1
                                        * Referenced by: '<S64>/Gain1'
                                        */
  real_T wbaseTs2_Gain;                /* Expression: SM.web*Ts/2
                                        * Referenced by: '<S75>/wbase*Ts//2'
                                        */
  real_T u5_Value_g[25];               /* Expression: eye(SM.nState,SM.nState)
                                        * Referenced by: '<S75>/u5'
                                        */
  real_T wbaseTs2_Gain_d;              /* Expression: SM.web*Ts/2
                                        * Referenced by: '<S75>/wbase*Ts//2 '
                                        */
  real_T Gain_Gain;                    /* Expression: sps.K1
                                        * Referenced by: '<S95>/Gain'
                                        */
  real_T Gain1_Gain_n;                 /* Expression: sps.K2
                                        * Referenced by: '<S95>/Gain1'
                                        */
  real_T Gain_Gain_n;                  /* Expression: sps.K1
                                        * Referenced by: '<S97>/Gain'
                                        */
  real_T Gain1_Gain_m;                 /* Expression: sps.K2
                                        * Referenced by: '<S97>/Gain1'
                                        */
  real_T Gain_Gain_p;                  /* Expression: sps.K1
                                        * Referenced by: '<S101>/Gain'
                                        */
  real_T Gain1_Gain_i;                 /* Expression: sps.K2
                                        * Referenced by: '<S101>/Gain1'
                                        */
  real_T Gain_Gain_ne;                 /* Expression: sps.K1
                                        * Referenced by: '<S103>/Gain'
                                        */
  real_T Gain1_Gain_a;                 /* Expression: sps.K2
                                        * Referenced by: '<S103>/Gain1'
                                        */
  real_T Gain_Gain_m;                  /* Expression: sps.K1
                                        * Referenced by: '<S107>/Gain'
                                        */
  real_T Gain1_Gain_md;                /* Expression: sps.K2
                                        * Referenced by: '<S107>/Gain1'
                                        */
  real_T Gain_Gain_b;                  /* Expression: sps.K1
                                        * Referenced by: '<S109>/Gain'
                                        */
  real_T Gain1_Gain_l;                 /* Expression: sps.K2
                                        * Referenced by: '<S109>/Gain1'
                                        */
  real_T Gain3_Gain;                   /* Expression: 1/3
                                        * Referenced by: '<S92>/Gain3'
                                        */
  real_T Gain_Gain_mq;                 /* Expression: sps.K1
                                        * Referenced by: '<S111>/Gain'
                                        */
  real_T Gain1_Gain_g;                 /* Expression: sps.K2
                                        * Referenced by: '<S111>/Gain1'
                                        */
  real_T Constant2_Value_l;            /* Expression: 0
                                        * Referenced by: '<S120>/Constant2'
                                        */
  real_T Laqd_nosat_Value[2];          /* Expression: [SM.Laq  SM.Lad]
                                        * Referenced by: '<S121>/Laqd_nosat'
                                        */
  real_T Lmqd_nosat_Value[2];          /* Expression: [SM.Lmq SM.Lmd]
                                        * Referenced by: '<S121>/Lmqd_nosat'
                                        */
  real_T Constant4_Value_m[25];        /* Expression: SM.Linv
                                        * Referenced by: '<S121>/Constant4'
                                        */
  real_T u2_Value;                     /* Expression: SM.Laq
                                        * Referenced by: '<S131>/u2'
                                        */
  real_T u3_Value_m;                   /* Expression: SM.Lmq
                                        * Referenced by: '<S131>/u3'
                                        */
  real_T Constant1_Value_e;            /* Expression: SM.Lmsatd(1)
                                        * Referenced by: '<S133>/Constant1'
                                        */
  real_T Ll_d_Gain_o[3];               /* Expression: [ 1/SM.Ll   1/SM.Llfd   1/SM.Llkd ]
                                        * Referenced by: '<S137>/1//Ll_d'
                                        */
  real_T u1_Value_ek[3];               /* Expression: [1/SM.Ll 1/SM.Llkd 1/SM.Llfd]
                                        * Referenced by: '<S136>/u1'
                                        */
  real_T Lmd_sat_InitialCondition_i;   /* Expression: SM.Lmsatd(1)
                                        * Referenced by: '<S133>/Lmd_sat'
                                        */
  real_T LookupTable_XData_f[2];       /* Expression: SM.Phisat
                                        * Referenced by: '<S133>/Lookup Table'
                                        */
  real_T LookupTable_YData_i[2];       /* Expression: [ 0 SM.Phisat(2:end)./SM.Lmsatd(2:end)*SM.Lmd ]
                                        * Referenced by: '<S133>/Lookup Table'
                                        */
  real_T Lmd_Gain_j;                   /* Expression: SM.Lmd
                                        * Referenced by: '<S133>/Lmd'
                                        */
  real_T u1_Value_k[25];               /* Expression: SM.R
                                        * Referenced by: '<S131>/u1'
                                        */
  real_T u1_Value_kt[25];              /* Expression: zeros(SM.nState,SM.nState)
                                        * Referenced by: '<S135>/u1'
                                        */
  real_T u5_Value_l[25];               /* Expression: SM.Llqd
                                        * Referenced by: '<S135>/u5'
                                        */
  real_T Constant6_Value_o[25];        /* Expression: SM.RLinv
                                        * Referenced by: '<S121>/Constant6'
                                        */
  real_T Constant2_Value_e;            /* Expression: SM.Sat
                                        * Referenced by: '<S121>/Constant2'
                                        */
  real_T Switch1_Threshold_a;          /* Expression: 0.5
                                        * Referenced by: '<S121>/Switch1'
                                        */
  real_T u1_Value_l[25];               /* Expression: zeros(SM.nState,SM.nState)
                                        * Referenced by: '<S129>/u1'
                                        */
  real_T Gain1_Gain_l3;                /* Expression: -1
                                        * Referenced by: '<S129>/Gain1'
                                        */
  real_T wbaseTs2_Gain_p;              /* Expression: SM.web*Ts/2
                                        * Referenced by: '<S140>/wbase*Ts//2'
                                        */
  real_T u5_Value_e[25];               /* Expression: eye(SM.nState,SM.nState)
                                        * Referenced by: '<S140>/u5'
                                        */
  real_T wbaseTs2_Gain_j;              /* Expression: SM.web*Ts/2
                                        * Referenced by: '<S140>/wbase*Ts//2 '
                                        */
  real_T com_Value;                    /* Expression: 1
                                        * Referenced by: '<S8>/com'
                                        */
  real_T Switch_Threshold;             /* Expression: 0.5
                                        * Referenced by: '<S8>/Switch'
                                        */
  real_T Switch1_Threshold_al;         /* Expression: 0.5
                                        * Referenced by: '<S8>/Switch1'
                                        */
  real_T Switch2_Threshold;            /* Expression: 0.5
                                        * Referenced by: '<S8>/Switch2'
                                        */
  real_T Rotorangledthetae_gainval;    /* Computed Parameter: Rotorangledthetae_gainval
                                        * Referenced by: '<S119>/Rotor angle dthetae'
                                        */
  real_T Rotorangledthetae_IC;         /* Expression: SM.tho
                                        * Referenced by: '<S119>/Rotor angle dthetae'
                                        */
  real_T web2_Gain;                    /* Expression: SM.web
                                        * Referenced by: '<S119>/web2'
                                        */
  real_T fluxes_InitialCondition[5];   /* Expression: SM.phiqd0
                                        * Referenced by: '<S132>/fluxes'
                                        */
  real_T Constant1_Value_d;            /* Expression: SM.Sat
                                        * Referenced by: '<S121>/Constant1'
                                        */
  real_T Constant3_Value;              /* Expression: SM.Sat
                                        * Referenced by: '<S121>/Constant3'
                                        */
  real_T Switch_Threshold_o;           /* Expression: 0.5
                                        * Referenced by: '<S121>/Switch'
                                        */
  real_T changeIqIdcurrentsigns_Gain[5];/* Expression: SM.IqdSign
                                         * Referenced by: '<S121>/change Iq Id  current signs'
                                         */
  real_T Ll_q_Gain[2];                 /* Expression: SM.One_Llq
                                        * Referenced by: '<S130>/1//Ll_q'
                                        */
  real_T Constant8_Value;              /* Expression: SM.Sat
                                        * Referenced by: '<S121>/Constant8'
                                        */
  real_T Switch3_Threshold;            /* Expression: 0.5
                                        * Referenced by: '<S121>/Switch3'
                                        */
  real_T Ll_d_Gain_ov[3];              /* Expression: [ 1/SM.Ll   1/SM.Llfd   1/SM.Llkd ]
                                        * Referenced by: '<S130>/1//Ll_d'
                                        */
  real_T SwitchCurrents_Value[9];      /* Expression: zeros(9,1)
                                        * Referenced by: '<S158>/SwitchCurrents'
                                        */
  real_T ib_Gain;                      /* Expression: SM.ib
                                        * Referenced by: '<S123>/ib'
                                        */
  real_T Rotorangledtheta_gainval;     /* Computed Parameter: Rotorangledtheta_gainval
                                        * Referenced by: '<S54>/Rotor angle dtheta'
                                        */
  real_T Rotorangledtheta_IC;          /* Expression: SM.tho
                                        * Referenced by: '<S54>/Rotor angle dtheta'
                                        */
  real_T web2_Gain_e;                  /* Expression: SM.web
                                        * Referenced by: '<S54>/web2'
                                        */
  real_T fluxes_InitialCondition_h[5]; /* Expression: SM.phiqd0
                                        * Referenced by: '<S67>/fluxes'
                                        */
  real_T Constant1_Value_b;            /* Expression: SM.Sat
                                        * Referenced by: '<S56>/Constant1'
                                        */
  real_T Constant3_Value_b;            /* Expression: SM.Sat
                                        * Referenced by: '<S56>/Constant3'
                                        */
  real_T Switch_Threshold_a;           /* Expression: 0.5
                                        * Referenced by: '<S56>/Switch'
                                        */
  real_T changeIqIdcurrentsigns_Gain_k[5];/* Expression: SM.IqdSign
                                           * Referenced by: '<S56>/change Iq Id  current signs'
                                           */
  real_T ib_Gain_f;                    /* Expression: SM.ib
                                        * Referenced by: '<S58>/ib'
                                        */
  real_T UnitDelay2_InitialCondition;  /* Expression: 1.1254*100
                                        * Referenced by: '<S33>/Unit Delay2'
                                        */
  real_T _Vb_Gain;                     /* Expression: 1/SM.Vb
                                        * Referenced by: '<S122>/1_Vb'
                                        */
  real_T Constant5_Value;              /* Expression: SM.Sat
                                        * Referenced by: '<S121>/Constant5'
                                        */
  real_T Switch2_Threshold_j;          /* Expression: 0.5
                                        * Referenced by: '<S121>/Switch2'
                                        */
  real_T Gain_Gain_k;                  /* Expression: 180/pi
                                        * Referenced by: '<S127>/Gain'
                                        */
  real_T Gain_Gain_g;                  /* Expression: SM.Gain1
                                        * Referenced by: '<S125>/Gain'
                                        */
  real_T Gain1_Gain_p;                 /* Expression: SM.Gain1
                                        * Referenced by: '<S125>/Gain1'
                                        */
  real_T outputformatting_Gain[18];    /* Expression: [SM.Ib2*ones(1,5),SM.N2,SM.Ib2,SM.Ib2,SM.Ib2,SM.phib,SM.phib,SM.Vb2,SM.Vb2,(SM.phib/SM.Ib2)*ones(1,2),1,1,1]
                                        * Referenced by: '<S120>/output formatting'
                                        */
  real_T nominalspeed_Value;           /* Expression: 1
                                        * Referenced by: '<S119>/nominal speed'
                                        */
  real_T dw_delay_InitialCondition;    /* Expression: SM.dwo
                                        * Referenced by: '<S141>/dw_delay'
                                        */
  real_T F2_Gain;                      /* Expression: 2
                                        * Referenced by: '<S141>/F2'
                                        */
  real_T dw_predict_InitialCondition;  /* Expression: SM.dwo
                                        * Referenced by: '<S141>/dw_predict'
                                        */
  real_T units_Gain;                   /* Expression: SM.Nb
                                        * Referenced by: '<S119>/units'
                                        */
  real_T radstopu_Gain;                /* Expression: 1/(pi*50)
                                        * Referenced by: '<S2>/rad//s to pu'
                                        */
  real_T pu2Watts_Gain;                /* Expression: 2e6
                                        * Referenced by: '<S2>/pu2Watts'
                                        */
  real_T _Pb_Gain;                     /* Expression: 1/SM.Nb
                                        * Referenced by: '<S54>/1_Pb'
                                        */
  real_T u_Gain[2];                    /* Expression: [1 -1]
                                        * Referenced by: '<S63>/1-1'
                                        */
  real_T units1_Gain;                  /* Expression: SM.Pb
                                        * Referenced by: '<S76>/units1'
                                        */
  real_T Delay_x_InitialCondition;     /* Expression: sps.x0
                                        * Referenced by: '<S28>/Delay_x'
                                        */
  real_T A_Gain;                       /* Expression: sps.A
                                        * Referenced by: '<S28>/A'
                                        */
  real_T Referencespeedpu_Value;       /* Expression: 1
                                        * Referenced by: '<Root>/Reference speed (pu)'
                                        */
  real_T D_Gain;                       /* Expression: sps.D
                                        * Referenced by: '<S29>/D'
                                        */
  real_T Delay_x_InitialCondition_i;   /* Expression: sps.x0
                                        * Referenced by: '<S29>/Delay_x'
                                        */
  real_T C_Gain;                       /* Expression: sps.C
                                        * Referenced by: '<S29>/C'
                                        */
  real_T D_Gain_i;                     /* Expression: sps.D
                                        * Referenced by: '<S30>/D'
                                        */
  real_T Delay_x_InitialCondition_f;   /* Expression: sps.x0
                                        * Referenced by: '<S30>/Delay_x'
                                        */
  real_T C_Gain_g;                     /* Expression: sps.C
                                        * Referenced by: '<S30>/C'
                                        */
  real_T B_Gain;                       /* Expression: sps.B
                                        * Referenced by: '<S28>/B'
                                        */
  real_T C_Gain_a;                     /* Expression: sps.C
                                        * Referenced by: '<S28>/C'
                                        */
  real_T D_Gain_j;                     /* Expression: sps.D
                                        * Referenced by: '<S28>/D'
                                        */
  real_T DiscreteTimeIntegrator_gainval;/* Computed Parameter: DiscreteTimeIntegrator_gainval
                                         * Referenced by: '<S2>/Discrete-Time Integrator'
                                         */
  real_T DiscreteTimeIntegrator_UpperSat;/* Expression: Tlim(2)
                                          * Referenced by: '<S2>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_LowerSat;/* Expression: Tlim(1)
                                          * Referenced by: '<S2>/Discrete-Time Integrator'
                                          */
  real_T A_Gain_o;                     /* Expression: sps.A
                                        * Referenced by: '<S29>/A'
                                        */
  real_T B_Gain_o;                     /* Expression: sps.B
                                        * Referenced by: '<S29>/B'
                                        */
  real_T A_Gain_e;                     /* Expression: sps.A
                                        * Referenced by: '<S30>/A'
                                        */
  real_T B_Gain_n;                     /* Expression: sps.B
                                        * Referenced by: '<S30>/B'
                                        */
  real_T UnitDelay4_InitialCondition;  /* Expression: 10.4482
                                        * Referenced by: '<S33>/Unit Delay4'
                                        */
  real_T donotdeletethisgain_Gain;     /* Expression: 1
                                        * Referenced by: '<S78>/do not delete this gain'
                                        */
  real_T Referencevoltagepu_Value;     /* Expression: 1
                                        * Referenced by: '<S3>/Reference voltage (pu)'
                                        */
  real_T UnitDelay1_InitialCondition;  /* Expression: 1
                                        * Referenced by: '<S36>/Unit Delay1'
                                        */
  real_T DiscreteTimeIntegrator_gainva_n;/* Computed Parameter: DiscreteTimeIntegrator_gainva_n
                                          * Referenced by: '<S85>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_UpperS_p;/* Expression: UpperLimit
                                          * Referenced by: '<S85>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_LowerS_b;/* Expression: LowerLimit
                                          * Referenced by: '<S85>/Discrete-Time Integrator'
                                          */
  real_T Saturation2_UpperSat;         /* Expression: UpperLimit
                                        * Referenced by: '<S85>/Saturation2'
                                        */
  real_T Saturation2_LowerSat;         /* Expression: LowerLimit
                                        * Referenced by: '<S85>/Saturation2'
                                        */
  real_T putovolts_Gain;               /* Expression: 11.56
                                        * Referenced by: '<S36>/pu to volts'
                                        */
  real_T N_Gain;                       /* Expression: SM.N
                                        * Referenced by: '<S52>/N'
                                        */
  real_T _Vb_Gain_p;                   /* Expression: 1/SM.Vb
                                        * Referenced by: '<S57>/1_Vb'
                                        */
  real_T Vkd0Vkq10Vkq20_Value[2];      /* Expression: zeros(1, SM.nState-3)
                                        * Referenced by: '<S52>/[ Vkd =0 Vkq1=0  Vkq2=0 ]'
                                        */
  real_T voltages_InitialCondition;    /* Expression: 0
                                        * Referenced by: '<S67>/voltages'
                                        */
  real_T IC_Threshold;                 /* Expression: Ts
                                        * Referenced by: '<S67>/IC'
                                        */
  real_T Nominalspeed_Value;           /* Expression: 1
                                        * Referenced by: '<S54>/Nominal speed'
                                        */
  real_T web1_Gain;                    /* Expression: SM.web
                                        * Referenced by: '<S54>/web1'
                                        */
  real_T donotdeletethisgain_Gain_h;   /* Expression: 1
                                        * Referenced by: '<S16>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_l;   /* Expression: 1
                                        * Referenced by: '<S17>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_ho;  /* Expression: 1
                                        * Referenced by: '<S18>/do not delete this gain'
                                        */
  real_T Kv1_Gain;                     /* Expression: Kv
                                        * Referenced by: '<S1>/Kv1'
                                        */
  real_T sinwt_Amp;                    /* Expression: sps.k
                                        * Referenced by: '<S87>/sin(wt)'
                                        */
  real_T sinwt_Bias;                   /* Expression: 0
                                        * Referenced by: '<S87>/sin(wt)'
                                        */
  real_T sinwt_Hsin;                   /* Computed Parameter: sinwt_Hsin
                                        * Referenced by: '<S87>/sin(wt)'
                                        */
  real_T sinwt_HCos;                   /* Computed Parameter: sinwt_HCos
                                        * Referenced by: '<S87>/sin(wt)'
                                        */
  real_T sinwt_PSin;                   /* Computed Parameter: sinwt_PSin
                                        * Referenced by: '<S87>/sin(wt)'
                                        */
  real_T sinwt_PCos;                   /* Computed Parameter: sinwt_PCos
                                        * Referenced by: '<S87>/sin(wt)'
                                        */
  real_T Integ4_gainval;               /* Computed Parameter: Integ4_gainval
                                        * Referenced by: '<S97>/Integ4'
                                        */
  real_T Integ4_IC;                    /* Expression: 0
                                        * Referenced by: '<S97>/Integ4'
                                        */
  real_T K1_Value;                     /* Expression: sps.Delay
                                        * Referenced by: '<S97>/K1'
                                        */
  real_T SFunction_P1_Size[2];         /* Computed Parameter: SFunction_P1_Size
                                        * Referenced by: '<S98>/S-Function'
                                        */
  real_T SFunction_P1;                 /* Expression: MaxDelay
                                        * Referenced by: '<S98>/S-Function'
                                        */
  real_T SFunction_P2_Size[2];         /* Computed Parameter: SFunction_P2_Size
                                        * Referenced by: '<S98>/S-Function'
                                        */
  real_T SFunction_P2;                 /* Expression: Ts
                                        * Referenced by: '<S98>/S-Function'
                                        */
  real_T SFunction_P3_Size[2];         /* Computed Parameter: SFunction_P3_Size
                                        * Referenced by: '<S98>/S-Function'
                                        */
  real_T SFunction_P3;                 /* Expression: InitialValue
                                        * Referenced by: '<S98>/S-Function'
                                        */
  real_T SFunction_P4_Size[2];         /* Computed Parameter: SFunction_P4_Size
                                        * Referenced by: '<S98>/S-Function'
                                        */
  real_T SFunction_P4;                 /* Expression: DFT
                                        * Referenced by: '<S98>/S-Function'
                                        */
  real_T K2_Value;                     /* Expression: sps.Freq
                                        * Referenced by: '<S97>/K2'
                                        */
  real_T UnitDelay_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S97>/Unit Delay'
                                        */
  real_T UnitDelay1_InitialCondition_l;/* Expression: sps.Vinit
                                        * Referenced by: '<S97>/Unit Delay1'
                                        */
  real_T coswt_Amp;                    /* Expression: sps.k
                                        * Referenced by: '<S87>/cos(wt)'
                                        */
  real_T coswt_Bias;                   /* Expression: 0
                                        * Referenced by: '<S87>/cos(wt)'
                                        */
  real_T coswt_Hsin;                   /* Computed Parameter: coswt_Hsin
                                        * Referenced by: '<S87>/cos(wt)'
                                        */
  real_T coswt_HCos;                   /* Computed Parameter: coswt_HCos
                                        * Referenced by: '<S87>/cos(wt)'
                                        */
  real_T coswt_PSin;                   /* Computed Parameter: coswt_PSin
                                        * Referenced by: '<S87>/cos(wt)'
                                        */
  real_T coswt_PCos;                   /* Computed Parameter: coswt_PCos
                                        * Referenced by: '<S87>/cos(wt)'
                                        */
  real_T Integ4_gainval_o;             /* Computed Parameter: Integ4_gainval_o
                                        * Referenced by: '<S95>/Integ4'
                                        */
  real_T Integ4_IC_k;                  /* Expression: 0
                                        * Referenced by: '<S95>/Integ4'
                                        */
  real_T K1_Value_f;                   /* Expression: sps.Delay
                                        * Referenced by: '<S95>/K1'
                                        */
  real_T SFunction_P1_Size_j[2];       /* Computed Parameter: SFunction_P1_Size_j
                                        * Referenced by: '<S96>/S-Function'
                                        */
  real_T SFunction_P1_l;               /* Expression: MaxDelay
                                        * Referenced by: '<S96>/S-Function'
                                        */
  real_T SFunction_P2_Size_l[2];       /* Computed Parameter: SFunction_P2_Size_l
                                        * Referenced by: '<S96>/S-Function'
                                        */
  real_T SFunction_P2_e;               /* Expression: Ts
                                        * Referenced by: '<S96>/S-Function'
                                        */
  real_T SFunction_P3_Size_k[2];       /* Computed Parameter: SFunction_P3_Size_k
                                        * Referenced by: '<S96>/S-Function'
                                        */
  real_T SFunction_P3_m;               /* Expression: InitialValue
                                        * Referenced by: '<S96>/S-Function'
                                        */
  real_T SFunction_P4_Size_f[2];       /* Computed Parameter: SFunction_P4_Size_f
                                        * Referenced by: '<S96>/S-Function'
                                        */
  real_T SFunction_P4_n;               /* Expression: DFT
                                        * Referenced by: '<S96>/S-Function'
                                        */
  real_T K2_Value_p;                   /* Expression: sps.Freq
                                        * Referenced by: '<S95>/K2'
                                        */
  real_T UnitDelay_InitialCondition_b; /* Expression: 0
                                        * Referenced by: '<S95>/Unit Delay'
                                        */
  real_T UnitDelay1_InitialCondition_f;/* Expression: sps.Vinit
                                        * Referenced by: '<S95>/Unit Delay1'
                                        */
  real_T RadDeg_Gain;                  /* Expression: 180/pi
                                        * Referenced by: '<S87>/Rad->Deg.'
                                        */
  real_T degrad_Gain;                  /* Expression: pi/180
                                        * Referenced by: '<S86>/deg->rad'
                                        */
  real_T sinwt_Amp_l;                  /* Expression: sps.k
                                        * Referenced by: '<S88>/sin(wt)'
                                        */
  real_T sinwt_Bias_i;                 /* Expression: 0
                                        * Referenced by: '<S88>/sin(wt)'
                                        */
  real_T sinwt_Hsin_j;                 /* Computed Parameter: sinwt_Hsin_j
                                        * Referenced by: '<S88>/sin(wt)'
                                        */
  real_T sinwt_HCos_p;                 /* Computed Parameter: sinwt_HCos_p
                                        * Referenced by: '<S88>/sin(wt)'
                                        */
  real_T sinwt_PSin_l;                 /* Computed Parameter: sinwt_PSin_l
                                        * Referenced by: '<S88>/sin(wt)'
                                        */
  real_T sinwt_PCos_d;                 /* Computed Parameter: sinwt_PCos_d
                                        * Referenced by: '<S88>/sin(wt)'
                                        */
  real_T Integ4_gainval_d;             /* Computed Parameter: Integ4_gainval_d
                                        * Referenced by: '<S103>/Integ4'
                                        */
  real_T Integ4_IC_l;                  /* Expression: 0
                                        * Referenced by: '<S103>/Integ4'
                                        */
  real_T K1_Value_i;                   /* Expression: sps.Delay
                                        * Referenced by: '<S103>/K1'
                                        */
  real_T SFunction_P1_Size_jj[2];      /* Computed Parameter: SFunction_P1_Size_jj
                                        * Referenced by: '<S104>/S-Function'
                                        */
  real_T SFunction_P1_b;               /* Expression: MaxDelay
                                        * Referenced by: '<S104>/S-Function'
                                        */
  real_T SFunction_P2_Size_d[2];       /* Computed Parameter: SFunction_P2_Size_d
                                        * Referenced by: '<S104>/S-Function'
                                        */
  real_T SFunction_P2_j;               /* Expression: Ts
                                        * Referenced by: '<S104>/S-Function'
                                        */
  real_T SFunction_P3_Size_a[2];       /* Computed Parameter: SFunction_P3_Size_a
                                        * Referenced by: '<S104>/S-Function'
                                        */
  real_T SFunction_P3_c;               /* Expression: InitialValue
                                        * Referenced by: '<S104>/S-Function'
                                        */
  real_T SFunction_P4_Size_c[2];       /* Computed Parameter: SFunction_P4_Size_c
                                        * Referenced by: '<S104>/S-Function'
                                        */
  real_T SFunction_P4_a;               /* Expression: DFT
                                        * Referenced by: '<S104>/S-Function'
                                        */
  real_T K2_Value_m;                   /* Expression: sps.Freq
                                        * Referenced by: '<S103>/K2'
                                        */
  real_T UnitDelay_InitialCondition_f; /* Expression: 0
                                        * Referenced by: '<S103>/Unit Delay'
                                        */
  real_T UnitDelay1_InitialCondition_h;/* Expression: sps.Vinit
                                        * Referenced by: '<S103>/Unit Delay1'
                                        */
  real_T coswt_Amp_h;                  /* Expression: sps.k
                                        * Referenced by: '<S88>/cos(wt)'
                                        */
  real_T coswt_Bias_h;                 /* Expression: 0
                                        * Referenced by: '<S88>/cos(wt)'
                                        */
  real_T coswt_Hsin_p;                 /* Computed Parameter: coswt_Hsin_p
                                        * Referenced by: '<S88>/cos(wt)'
                                        */
  real_T coswt_HCos_m;                 /* Computed Parameter: coswt_HCos_m
                                        * Referenced by: '<S88>/cos(wt)'
                                        */
  real_T coswt_PSin_k;                 /* Computed Parameter: coswt_PSin_k
                                        * Referenced by: '<S88>/cos(wt)'
                                        */
  real_T coswt_PCos_n;                 /* Computed Parameter: coswt_PCos_n
                                        * Referenced by: '<S88>/cos(wt)'
                                        */
  real_T Integ4_gainval_e;             /* Computed Parameter: Integ4_gainval_e
                                        * Referenced by: '<S101>/Integ4'
                                        */
  real_T Integ4_IC_ly;                 /* Expression: 0
                                        * Referenced by: '<S101>/Integ4'
                                        */
  real_T K1_Value_b;                   /* Expression: sps.Delay
                                        * Referenced by: '<S101>/K1'
                                        */
  real_T SFunction_P1_Size_a[2];       /* Computed Parameter: SFunction_P1_Size_a
                                        * Referenced by: '<S102>/S-Function'
                                        */
  real_T SFunction_P1_p;               /* Expression: MaxDelay
                                        * Referenced by: '<S102>/S-Function'
                                        */
  real_T SFunction_P2_Size_i[2];       /* Computed Parameter: SFunction_P2_Size_i
                                        * Referenced by: '<S102>/S-Function'
                                        */
  real_T SFunction_P2_d;               /* Expression: Ts
                                        * Referenced by: '<S102>/S-Function'
                                        */
  real_T SFunction_P3_Size_b[2];       /* Computed Parameter: SFunction_P3_Size_b
                                        * Referenced by: '<S102>/S-Function'
                                        */
  real_T SFunction_P3_n;               /* Expression: InitialValue
                                        * Referenced by: '<S102>/S-Function'
                                        */
  real_T SFunction_P4_Size_cp[2];      /* Computed Parameter: SFunction_P4_Size_cp
                                        * Referenced by: '<S102>/S-Function'
                                        */
  real_T SFunction_P4_l;               /* Expression: DFT
                                        * Referenced by: '<S102>/S-Function'
                                        */
  real_T K2_Value_o;                   /* Expression: sps.Freq
                                        * Referenced by: '<S101>/K2'
                                        */
  real_T UnitDelay_InitialCondition_i; /* Expression: 0
                                        * Referenced by: '<S101>/Unit Delay'
                                        */
  real_T UnitDelay1_InitialCondition_k;/* Expression: sps.Vinit
                                        * Referenced by: '<S101>/Unit Delay1'
                                        */
  real_T RadDeg_Gain_e;                /* Expression: 180/pi
                                        * Referenced by: '<S88>/Rad->Deg.'
                                        */
  real_T degrad1_Gain;                 /* Expression: pi/180
                                        * Referenced by: '<S86>/deg->rad1'
                                        */
  real_T sinwt_Amp_e;                  /* Expression: sps.k
                                        * Referenced by: '<S89>/sin(wt)'
                                        */
  real_T sinwt_Bias_n;                 /* Expression: 0
                                        * Referenced by: '<S89>/sin(wt)'
                                        */
  real_T sinwt_Hsin_d;                 /* Computed Parameter: sinwt_Hsin_d
                                        * Referenced by: '<S89>/sin(wt)'
                                        */
  real_T sinwt_HCos_k;                 /* Computed Parameter: sinwt_HCos_k
                                        * Referenced by: '<S89>/sin(wt)'
                                        */
  real_T sinwt_PSin_k;                 /* Computed Parameter: sinwt_PSin_k
                                        * Referenced by: '<S89>/sin(wt)'
                                        */
  real_T sinwt_PCos_c;                 /* Computed Parameter: sinwt_PCos_c
                                        * Referenced by: '<S89>/sin(wt)'
                                        */
  real_T Integ4_gainval_l;             /* Computed Parameter: Integ4_gainval_l
                                        * Referenced by: '<S109>/Integ4'
                                        */
  real_T Integ4_IC_o;                  /* Expression: 0
                                        * Referenced by: '<S109>/Integ4'
                                        */
  real_T K1_Value_g;                   /* Expression: sps.Delay
                                        * Referenced by: '<S109>/K1'
                                        */
  real_T SFunction_P1_Size_l[2];       /* Computed Parameter: SFunction_P1_Size_l
                                        * Referenced by: '<S110>/S-Function'
                                        */
  real_T SFunction_P1_a;               /* Expression: MaxDelay
                                        * Referenced by: '<S110>/S-Function'
                                        */
  real_T SFunction_P2_Size_g[2];       /* Computed Parameter: SFunction_P2_Size_g
                                        * Referenced by: '<S110>/S-Function'
                                        */
  real_T SFunction_P2_k;               /* Expression: Ts
                                        * Referenced by: '<S110>/S-Function'
                                        */
  real_T SFunction_P3_Size_d[2];       /* Computed Parameter: SFunction_P3_Size_d
                                        * Referenced by: '<S110>/S-Function'
                                        */
  real_T SFunction_P3_i;               /* Expression: InitialValue
                                        * Referenced by: '<S110>/S-Function'
                                        */
  real_T SFunction_P4_Size_a[2];       /* Computed Parameter: SFunction_P4_Size_a
                                        * Referenced by: '<S110>/S-Function'
                                        */
  real_T SFunction_P4_g;               /* Expression: DFT
                                        * Referenced by: '<S110>/S-Function'
                                        */
  real_T K2_Value_mj;                  /* Expression: sps.Freq
                                        * Referenced by: '<S109>/K2'
                                        */
  real_T UnitDelay_InitialCondition_d; /* Expression: 0
                                        * Referenced by: '<S109>/Unit Delay'
                                        */
  real_T UnitDelay1_InitialCondition_f5;/* Expression: sps.Vinit
                                         * Referenced by: '<S109>/Unit Delay1'
                                         */
  real_T coswt_Amp_a;                  /* Expression: sps.k
                                        * Referenced by: '<S89>/cos(wt)'
                                        */
  real_T coswt_Bias_a;                 /* Expression: 0
                                        * Referenced by: '<S89>/cos(wt)'
                                        */
  real_T coswt_Hsin_ph;                /* Computed Parameter: coswt_Hsin_ph
                                        * Referenced by: '<S89>/cos(wt)'
                                        */
  real_T coswt_HCos_b;                 /* Computed Parameter: coswt_HCos_b
                                        * Referenced by: '<S89>/cos(wt)'
                                        */
  real_T coswt_PSin_j;                 /* Computed Parameter: coswt_PSin_j
                                        * Referenced by: '<S89>/cos(wt)'
                                        */
  real_T coswt_PCos_d;                 /* Computed Parameter: coswt_PCos_d
                                        * Referenced by: '<S89>/cos(wt)'
                                        */
  real_T Integ4_gainval_a;             /* Computed Parameter: Integ4_gainval_a
                                        * Referenced by: '<S107>/Integ4'
                                        */
  real_T Integ4_IC_a;                  /* Expression: 0
                                        * Referenced by: '<S107>/Integ4'
                                        */
  real_T K1_Value_l;                   /* Expression: sps.Delay
                                        * Referenced by: '<S107>/K1'
                                        */
  real_T SFunction_P1_Size_b[2];       /* Computed Parameter: SFunction_P1_Size_b
                                        * Referenced by: '<S108>/S-Function'
                                        */
  real_T SFunction_P1_e;               /* Expression: MaxDelay
                                        * Referenced by: '<S108>/S-Function'
                                        */
  real_T SFunction_P2_Size_la[2];      /* Computed Parameter: SFunction_P2_Size_la
                                        * Referenced by: '<S108>/S-Function'
                                        */
  real_T SFunction_P2_c;               /* Expression: Ts
                                        * Referenced by: '<S108>/S-Function'
                                        */
  real_T SFunction_P3_Size_a3[2];      /* Computed Parameter: SFunction_P3_Size_a3
                                        * Referenced by: '<S108>/S-Function'
                                        */
  real_T SFunction_P3_j;               /* Expression: InitialValue
                                        * Referenced by: '<S108>/S-Function'
                                        */
  real_T SFunction_P4_Size_n[2];       /* Computed Parameter: SFunction_P4_Size_n
                                        * Referenced by: '<S108>/S-Function'
                                        */
  real_T SFunction_P4_k;               /* Expression: DFT
                                        * Referenced by: '<S108>/S-Function'
                                        */
  real_T K2_Value_g;                   /* Expression: sps.Freq
                                        * Referenced by: '<S107>/K2'
                                        */
  real_T UnitDelay_InitialCondition_k; /* Expression: 0
                                        * Referenced by: '<S107>/Unit Delay'
                                        */
  real_T UnitDelay1_InitialCondition_e;/* Expression: sps.Vinit
                                        * Referenced by: '<S107>/Unit Delay1'
                                        */
  real_T RadDeg_Gain_b;                /* Expression: 180/pi
                                        * Referenced by: '<S89>/Rad->Deg.'
                                        */
  real_T degrad2_Gain;                 /* Expression: pi/180
                                        * Referenced by: '<S86>/deg->rad2'
                                        */
  real_T Constant_Value;               /* Expression: sps.PosOn
                                        * Referenced by: '<S86>/Constant'
                                        */
  real_T Constant1_Value_j;            /* Expression: sps.NegOn
                                        * Referenced by: '<S86>/Constant1'
                                        */
  real_T Constant2_Value_m;            /* Expression: sps.ZeroOn
                                        * Referenced by: '<S86>/Constant2'
                                        */
  real_T donotdeletethisgain_Gain_p;   /* Expression: 1
                                        * Referenced by: '<S40>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_b;   /* Expression: 1
                                        * Referenced by: '<S41>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_i;   /* Expression: 1
                                        * Referenced by: '<S42>/do not delete this gain'
                                        */
  real_T Kv_Gain;                      /* Expression: Ki
                                        * Referenced by: '<S31>/Kv'
                                        */
  real_T donotdeletethisgain_Gain_lw;  /* Expression: 1
                                        * Referenced by: '<S43>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_a;   /* Expression: 1
                                        * Referenced by: '<S44>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_f;   /* Expression: 1
                                        * Referenced by: '<S45>/do not delete this gain'
                                        */
  real_T Kv1_Gain_m;                   /* Expression: Kv
                                        * Referenced by: '<S31>/Kv1'
                                        */
  real_T g_Value[6];                   /* Expression: zeros(1,6)
                                        * Referenced by: '<S34>/g'
                                        */
  real_T N_Gain_j;                     /* Expression: SM.N
                                        * Referenced by: '<S117>/N'
                                        */
  real_T u_Gain_k[2];                  /* Expression: [1 -1]
                                        * Referenced by: '<S128>/1-1'
                                        */
  real_T Vkd0Vkq10Vkq20_Value_a[2];    /* Expression: zeros(1, SM.nState-3)
                                        * Referenced by: '<S117>/[ Vkd =0 Vkq1=0  Vkq2=0 ]'
                                        */
  real_T voltages_InitialCondition_c;  /* Expression: 0
                                        * Referenced by: '<S132>/voltages'
                                        */
  real_T IC_Threshold_e;               /* Expression: Ts
                                        * Referenced by: '<S132>/IC'
                                        */
  real_T _Pb_Gain_f;                   /* Expression: 1/SM.Pb
                                        * Referenced by: '<S119>/1_Pb'
                                        */
  real_T F_Gain;                       /* Expression: SM.F
                                        * Referenced by: '<S119>/F'
                                        */
  real_T uH_Gain;                      /* Expression: 1/(2*SM.H)
                                        * Referenced by: '<S119>/1 ----- 2H'
                                        */
  real_T Rotorspeeddeviationdw_gainval;/* Computed Parameter: Rotorspeeddeviationdw_gainval
                                        * Referenced by: '<S119>/Rotor speed deviation (dw)'
                                        */
  real_T Rotorspeeddeviationdw_IC;     /* Expression: SM.dwo
                                        * Referenced by: '<S119>/Rotor speed deviation (dw)'
                                        */
  real_T webase_Gain;                  /* Expression: SM.web
                                        * Referenced by: '<S119>/we base'
                                        */
  real_T Integ4_gainval_j;             /* Computed Parameter: Integ4_gainval_j
                                        * Referenced by: '<S111>/Integ4'
                                        */
  real_T Integ4_IC_oa;                 /* Expression: 0
                                        * Referenced by: '<S111>/Integ4'
                                        */
  real_T K1_Value_fz;                  /* Expression: sps.Delay
                                        * Referenced by: '<S111>/K1'
                                        */
  real_T SFunction_P1_Size_jv[2];      /* Computed Parameter: SFunction_P1_Size_jv
                                        * Referenced by: '<S112>/S-Function'
                                        */
  real_T SFunction_P1_pd;              /* Expression: MaxDelay
                                        * Referenced by: '<S112>/S-Function'
                                        */
  real_T SFunction_P2_Size_k[2];       /* Computed Parameter: SFunction_P2_Size_k
                                        * Referenced by: '<S112>/S-Function'
                                        */
  real_T SFunction_P2_jc;              /* Expression: Ts
                                        * Referenced by: '<S112>/S-Function'
                                        */
  real_T SFunction_P3_Size_i[2];       /* Computed Parameter: SFunction_P3_Size_i
                                        * Referenced by: '<S112>/S-Function'
                                        */
  real_T SFunction_P3_k;               /* Expression: InitialValue
                                        * Referenced by: '<S112>/S-Function'
                                        */
  real_T SFunction_P4_Size_p[2];       /* Computed Parameter: SFunction_P4_Size_p
                                        * Referenced by: '<S112>/S-Function'
                                        */
  real_T SFunction_P4_o;               /* Expression: DFT
                                        * Referenced by: '<S112>/S-Function'
                                        */
  real_T K2_Value_h;                   /* Expression: sps.Freq
                                        * Referenced by: '<S111>/K2'
                                        */
  real_T UnitDelay_InitialCondition_h; /* Expression: 0
                                        * Referenced by: '<S111>/Unit Delay'
                                        */
  real_T UnitDelay1_InitialCondition_lb;/* Expression: sps.Vinit
                                         * Referenced by: '<S111>/Unit Delay1'
                                         */
  real_T donotdeletethisgain_Gain_lw2; /* Expression: 1
                                        * Referenced by: '<S13>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_lq;  /* Expression: 1
                                        * Referenced by: '<S14>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_ad;  /* Expression: 1
                                        * Referenced by: '<S15>/do not delete this gain'
                                        */
  real_T Kv_Gain_d;                    /* Expression: Ki
                                        * Referenced by: '<S1>/Kv'
                                        */
  real_T Gain2_Gain;                   /* Expression: Gain
                                        * Referenced by: '<S5>/Gain2'
                                        */
  real_T Mode_Value;                   /* Expression: Display
                                        * Referenced by: '<S5>/Mode'
                                        */
  real_T Gain2_Gain_g;                 /* Expression: Gain
                                        * Referenced by: '<S6>/Gain2'
                                        */
  real_T Mode_Value_c;                 /* Expression: Display
                                        * Referenced by: '<S6>/Mode'
                                        */
  real_T u_Value;                      /* Expression: 0
                                        * Referenced by: '<S146>/0 1'
                                        */
  real_T Gain_Gain_f;                  /* Expression: 2
                                        * Referenced by: '<S146>/Gain'
                                        */
  real_T DelayTs_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S146>/Delay Ts'
                                        */
  real_T DelayTs_InitialCondition_c;   /* Expression: 0
                                        * Referenced by: '<S146>/Delay Ts '
                                        */
  real_T C4_Value;                     /* Expression: External
                                        * Referenced by: '<S8>/C4'
                                        */
  real_T LookUpTable_XData[4];         /* Expression: sps.tv
                                        * Referenced by: '<S155>/Look-Up Table'
                                        */
  real_T LookUpTable_YData[4];         /* Expression: sps.opv
                                        * Referenced by: '<S155>/Look-Up Table'
                                        */
  real_T Switch3_Threshold_f;          /* Expression: 0.5
                                        * Referenced by: '<S8>/Switch3'
                                        */
  real_T Constant5_Value_p;            /* Expression: InitialState
                                        * Referenced by: '<S8>/Constant5'
                                        */
  real_T C4_Value_n;                   /* Expression: BR.com
                                        * Referenced by: '<S146>/C4'
                                        */
  real_T LookUpTable_XData_d[4];       /* Expression: sps.tv
                                        * Referenced by: '<S148>/Look-Up Table'
                                        */
  real_T LookUpTable_YData_o[4];       /* Expression: sps.opv
                                        * Referenced by: '<S148>/Look-Up Table'
                                        */
  real_T Switch3_Threshold_p;          /* Expression: 0.5
                                        * Referenced by: '<S146>/Switch3'
                                        */
  real_T u5_OnVal;                     /* Expression: 1.5
                                        * Referenced by: '<S146>/>1.5'
                                        */
  real_T u5_OffVal;                    /* Expression: 0
                                        * Referenced by: '<S146>/>1.5'
                                        */
  real_T u5_YOn;                       /* Expression: 1
                                        * Referenced by: '<S146>/>1.5'
                                        */
  real_T u5_YOff;                      /* Expression: 0
                                        * Referenced by: '<S146>/>1.5'
                                        */
  real_T IC_Value;                     /* Expression: double(InitialState)
                                        * Referenced by: '<S146>/IC'
                                        */
  real_T u_Value_p;                    /* Expression: 0
                                        * Referenced by: '<S149>/0 1'
                                        */
  real_T Gain_Gain_j;                  /* Expression: 2
                                        * Referenced by: '<S149>/Gain'
                                        */
  real_T DelayTs_InitialCondition_p;   /* Expression: 0
                                        * Referenced by: '<S149>/Delay Ts'
                                        */
  real_T DelayTs_InitialCondition_e;   /* Expression: 0
                                        * Referenced by: '<S149>/Delay Ts '
                                        */
  real_T C4_Value_j;                   /* Expression: BR.com
                                        * Referenced by: '<S149>/C4'
                                        */
  real_T LookUpTable_XData_a[4];       /* Expression: sps.tv
                                        * Referenced by: '<S151>/Look-Up Table'
                                        */
  real_T LookUpTable_YData_c[4];       /* Expression: sps.opv
                                        * Referenced by: '<S151>/Look-Up Table'
                                        */
  real_T Switch3_Threshold_i;          /* Expression: 0.5
                                        * Referenced by: '<S149>/Switch3'
                                        */
  real_T u5_OnVal_p;                   /* Expression: 1.5
                                        * Referenced by: '<S149>/>1.5'
                                        */
  real_T u5_OffVal_m;                  /* Expression: 0
                                        * Referenced by: '<S149>/>1.5'
                                        */
  real_T u5_YOn_f;                     /* Expression: 1
                                        * Referenced by: '<S149>/>1.5'
                                        */
  real_T u5_YOff_g;                    /* Expression: 0
                                        * Referenced by: '<S149>/>1.5'
                                        */
  real_T IC_Value_j;                   /* Expression: double(InitialState)
                                        * Referenced by: '<S149>/IC'
                                        */
  real_T u_Value_d;                    /* Expression: 0
                                        * Referenced by: '<S152>/0 1'
                                        */
  real_T Gain_Gain_bi;                 /* Expression: 2
                                        * Referenced by: '<S152>/Gain'
                                        */
  real_T DelayTs_InitialCondition_d;   /* Expression: 0
                                        * Referenced by: '<S152>/Delay Ts'
                                        */
  real_T DelayTs_InitialCondition_p4;  /* Expression: 0
                                        * Referenced by: '<S152>/Delay Ts '
                                        */
  real_T C4_Value_m;                   /* Expression: BR.com
                                        * Referenced by: '<S152>/C4'
                                        */
  real_T LookUpTable_XData_m[4];       /* Expression: sps.tv
                                        * Referenced by: '<S154>/Look-Up Table'
                                        */
  real_T LookUpTable_YData_n[4];       /* Expression: sps.opv
                                        * Referenced by: '<S154>/Look-Up Table'
                                        */
  real_T Switch3_Threshold_c;          /* Expression: 0.5
                                        * Referenced by: '<S152>/Switch3'
                                        */
  real_T u5_OnVal_o;                   /* Expression: 1.5
                                        * Referenced by: '<S152>/>1.5'
                                        */
  real_T u5_OffVal_j;                  /* Expression: 0
                                        * Referenced by: '<S152>/>1.5'
                                        */
  real_T u5_YOn_g;                     /* Expression: 1
                                        * Referenced by: '<S152>/>1.5'
                                        */
  real_T u5_YOff_b;                    /* Expression: 0
                                        * Referenced by: '<S152>/>1.5'
                                        */
  real_T IC_Value_n;                   /* Expression: double(InitialState)
                                        * Referenced by: '<S152>/IC'
                                        */
  boolean_T Constant1_Value_n;         /* Expression: SM.nState==6
                                        * Referenced by: '<S66>/Constant1'
                                        */
  boolean_T Constant2_Value_ea;        /* Expression: SM.nState==6
                                        * Referenced by: '<S66>/Constant2'
                                        */
  boolean_T Constant1_Value_f;         /* Expression: SM.nState==6
                                        * Referenced by: '<S131>/Constant1'
                                        */
  boolean_T Constant2_Value_o;         /* Expression: SM.nState==6
                                        * Referenced by: '<S131>/Constant2'
                                        */
  boolean_T Constant3_Value_p;         /* Expression: SM.nState==6
                                        * Referenced by: '<S131>/Constant3'
                                        */
  boolean_T Constant1_Value_a;         /* Expression: SM.nState==6
                                        * Referenced by: '<S120>/Constant1'
                                        */
  P_Lmq_sat_Diesel_DER2_T Lmq_sat_l;   /* '<S131>/Lmq_sat' */
  P_NegSeqComputation_Diesel_DE_T PosSeqComputation;/* '<S86>/Pos. Seq. Computation' */
  P_NegSeqComputation_Diesel_DE_T NegSeqComputation;/* '<S86>/Neg. Seq. Computation' */
  P_Lmq_sat_Diesel_DER2_T Lmq_sat;     /* '<S66>/Lmq_sat' */
};

/* Real-time Model Data Structure */
struct tag_RTM_Diesel_DER2_T {
  struct SimStruct_tag * *childSfunctions;
  const char_T * volatile errorStatus;
  SS_SimMode simMode;
  RTWSolverInfo solverInfo;
  RTWSolverInfo *solverInfoPtr;
  void *sfcnInfo;

  /*
   * NonInlinedSFcns:
   * The following substructure contains information regarding
   * non-inlined s-functions used in the model.
   */
  struct {
    RTWSfcnInfo sfcnInfo;
    time_T *taskTimePtrs[2];
    SimStruct childSFunctions[7];
    SimStruct *childSFunctionPtrs[7];
    struct _ssBlkInfo2 blkInfo2[7];
    struct _ssSFcnModelMethods2 methods2[7];
    struct _ssSFcnModelMethods3 methods3[7];
    struct _ssStatesInfo2 statesInfo2[7];
    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn0;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn1;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn2;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn3;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn4;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn5;

    struct {
      time_T sfcnPeriod[1];
      time_T sfcnOffset[1];
      int_T sfcnTsMap[1];
      struct _ssPortInputs inputPortInfo[2];
      real_T const *UPtrs0[1];
      real_T const *UPtrs1[1];
      struct _ssPortOutputs outputPortInfo[1];
      uint_T attribs[4];
      mxArray *params[4];
      struct _ssDWorkRecord dWork[3];
      struct _ssDWorkAuxRecord dWorkAux[3];
    } Sfcn6;
  } NonInlinedSFcns;

  /*
   * ModelData:
   * The following substructure contains information regarding
   * the data used in the model.
   */
  struct {
    boolean_T zCCacheNeedsReset;
    boolean_T derivCacheNeedsReset;
    boolean_T blkStateChange;
  } ModelData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    uint32_T options;
    int_T numContStates;
    int_T numU;
    int_T numY;
    int_T numSampTimes;
    int_T numBlocks;
    int_T numBlockIO;
    int_T numBlockPrms;
    int_T numDwork;
    int_T numSFcnPrms;
    int_T numSFcns;
    int_T numIports;
    int_T numOports;
    int_T numNonSampZCs;
    int_T sysDirFeedThru;
    int_T rtwGenSfcn;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    time_T stepSize;
    uint32_T clockTick0;
    time_T stepSize0;
    uint32_T clockTick1;
    time_T stepSize1;
    time_T tStart;
    time_T tFinal;
    time_T timeOfLastOutput;
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *sampleTimes;
    time_T *offsetTimes;
    int_T *sampleTimeTaskIDPtr;
    int_T *sampleHits;
    int_T *perTaskSampleHits;
    time_T *t;
    time_T sampleTimesArray[2];
    time_T offsetTimesArray[2];
    int_T sampleTimeTaskIDArray[2];
    int_T sampleHitArray[2];
    int_T perTaskSampleHitsArray[4];
    time_T tArray[2];
  } Timing;
};

/* Block parameters (auto storage) */
extern P_Diesel_DER2_T Diesel_DER2_P;

/* Block signals (auto storage) */
extern B_Diesel_DER2_T Diesel_DER2_B;

/* Block states (auto storage) */
extern DW_Diesel_DER2_T Diesel_DER2_DW;

/* External outputs (root outports fed by signals with auto storage) */
extern ExtY_Diesel_DER2_T Diesel_DER2_Y;

/* Model entry point functions */
extern void Diesel_DER2_initialize(void);
extern void Diesel_DER2_step(void);
extern void Diesel_DER2_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Diesel_DER2_T *const Diesel_DER2_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'Diesel_DER2'
 * '<S1>'   : 'Diesel_DER2/BSM'
 * '<S2>'   : 'Diesel_DER2/Diesel Engine & Speed Regulator'
 * '<S3>'   : 'Diesel_DER2/Exciter'
 * '<S4>'   : 'Diesel_DER2/Mean'
 * '<S5>'   : 'Diesel_DER2/Multimeter'
 * '<S6>'   : 'Diesel_DER2/Multimeter1'
 * '<S7>'   : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V'
 * '<S8>'   : 'Diesel_DER2/Three-Phase Breaker'
 * '<S9>'   : 'Diesel_DER2/powergui'
 * '<S10>'  : 'Diesel_DER2/BSM/Mode I'
 * '<S11>'  : 'Diesel_DER2/BSM/Mode V'
 * '<S12>'  : 'Diesel_DER2/BSM/Model'
 * '<S13>'  : 'Diesel_DER2/BSM/Model/I A:'
 * '<S14>'  : 'Diesel_DER2/BSM/Model/I B:'
 * '<S15>'  : 'Diesel_DER2/BSM/Model/I C:'
 * '<S16>'  : 'Diesel_DER2/BSM/Model/U A:'
 * '<S17>'  : 'Diesel_DER2/BSM/Model/U B:'
 * '<S18>'  : 'Diesel_DER2/BSM/Model/U C:'
 * '<S19>'  : 'Diesel_DER2/BSM/Model/I A:/Model'
 * '<S20>'  : 'Diesel_DER2/BSM/Model/I B:/Model'
 * '<S21>'  : 'Diesel_DER2/BSM/Model/I C:/Model'
 * '<S22>'  : 'Diesel_DER2/BSM/Model/U A:/Model'
 * '<S23>'  : 'Diesel_DER2/BSM/Model/U B:/Model'
 * '<S24>'  : 'Diesel_DER2/BSM/Model/U C:/Model'
 * '<S25>'  : 'Diesel_DER2/Diesel Engine & Speed Regulator/1//(1+T*s)'
 * '<S26>'  : 'Diesel_DER2/Diesel Engine & Speed Regulator/Lead-Lag Filter'
 * '<S27>'  : 'Diesel_DER2/Diesel Engine & Speed Regulator/Lead-Lag Filter2'
 * '<S28>'  : 'Diesel_DER2/Diesel Engine & Speed Regulator/1//(1+T*s)/Model'
 * '<S29>'  : 'Diesel_DER2/Diesel Engine & Speed Regulator/Lead-Lag Filter/Model'
 * '<S30>'  : 'Diesel_DER2/Diesel Engine & Speed Regulator/Lead-Lag Filter2/Model'
 * '<S31>'  : 'Diesel_DER2/Exciter/BExc'
 * '<S32>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V'
 * '<S33>'  : 'Diesel_DER2/Exciter/Field Connection'
 * '<S34>'  : 'Diesel_DER2/Exciter/Rectifier'
 * '<S35>'  : 'Diesel_DER2/Exciter/Transformer 10 kVA 400V// 12V'
 * '<S36>'  : 'Diesel_DER2/Exciter/Voltage Regulator'
 * '<S37>'  : 'Diesel_DER2/Exciter/BExc/Mode I'
 * '<S38>'  : 'Diesel_DER2/Exciter/BExc/Mode V'
 * '<S39>'  : 'Diesel_DER2/Exciter/BExc/Model'
 * '<S40>'  : 'Diesel_DER2/Exciter/BExc/Model/I A:'
 * '<S41>'  : 'Diesel_DER2/Exciter/BExc/Model/I B:'
 * '<S42>'  : 'Diesel_DER2/Exciter/BExc/Model/I C:'
 * '<S43>'  : 'Diesel_DER2/Exciter/BExc/Model/U A:'
 * '<S44>'  : 'Diesel_DER2/Exciter/BExc/Model/U B:'
 * '<S45>'  : 'Diesel_DER2/Exciter/BExc/Model/U C:'
 * '<S46>'  : 'Diesel_DER2/Exciter/BExc/Model/I A:/Model'
 * '<S47>'  : 'Diesel_DER2/Exciter/BExc/Model/I B:/Model'
 * '<S48>'  : 'Diesel_DER2/Exciter/BExc/Model/I C:/Model'
 * '<S49>'  : 'Diesel_DER2/Exciter/BExc/Model/U A:/Model'
 * '<S50>'  : 'Diesel_DER2/Exciter/BExc/Model/U B:/Model'
 * '<S51>'  : 'Diesel_DER2/Exciter/BExc/Model/U C:/Model'
 * '<S52>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model'
 * '<S53>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Measurements'
 * '<S54>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Mechanical model'
 * '<S55>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Machine measurements'
 * '<S56>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Synchronous Machine Discrete Model'
 * '<S57>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/abc to qd transformation'
 * '<S58>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/qd to abc transformation'
 * '<S59>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Machine measurements/Delta angle'
 * '<S60>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Machine measurements/PQ'
 * '<S61>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Machine measurements/Delta angle/Cartesian to Polar'
 * '<S62>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Machine measurements/Delta angle/Radians to Degrees'
 * '<S63>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Synchronous Machine Discrete Model/Electromagnetic Torque'
 * '<S64>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Synchronous Machine Discrete Model/Matrix W'
 * '<S65>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Synchronous Machine Discrete Model/Phimqd'
 * '<S66>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation'
 * '<S67>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Synchronous Machine Discrete Model/phi'
 * '<S68>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmd_sat'
 * '<S69>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmq_sat'
 * '<S70>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation/Update Matrix L'
 * '<S71>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmd_sat/Lad'
 * '<S72>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmd_sat/phimd'
 * '<S73>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmq_sat/Laq'
 * '<S74>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmq_sat/phimq'
 * '<S75>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Electrical model/Synchronous Machine Discrete Model/phi/Subsystem'
 * '<S76>'  : 'Diesel_DER2/Exciter/Exciter System 8.1 kVA 400V/Mechanical model/measurements'
 * '<S77>'  : 'Diesel_DER2/Exciter/Field Connection/Ifd'
 * '<S78>'  : 'Diesel_DER2/Exciter/Field Connection/Vf '
 * '<S79>'  : 'Diesel_DER2/Exciter/Field Connection/Vf /Model'
 * '<S80>'  : 'Diesel_DER2/Exciter/Rectifier/Model'
 * '<S81>'  : 'Diesel_DER2/Exciter/Transformer 10 kVA 400V// 12V/Model'
 * '<S82>'  : 'Diesel_DER2/Exciter/Transformer 10 kVA 400V// 12V/Model/Linear'
 * '<S83>'  : 'Diesel_DER2/Exciter/Transformer 10 kVA 400V// 12V/Model/Linear1'
 * '<S84>'  : 'Diesel_DER2/Exciter/Transformer 10 kVA 400V// 12V/Model/Linear2'
 * '<S85>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Discrete PI Controller'
 * '<S86>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer'
 * '<S87>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_A'
 * '<S88>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_B'
 * '<S89>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_C'
 * '<S90>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Neg. Seq. Computation'
 * '<S91>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Pos. Seq. Computation'
 * '<S92>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Zero Seq. Computation'
 * '<S93>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_A/Mean'
 * '<S94>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_A/Mean value1'
 * '<S95>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_A/Mean/Model'
 * '<S96>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_A/Mean/Model/Discrete Variable Time Delay'
 * '<S97>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_A/Mean value1/Model'
 * '<S98>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_A/Mean value1/Model/Discrete Variable Time Delay'
 * '<S99>'  : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_B/Mean'
 * '<S100>' : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_B/Mean value1'
 * '<S101>' : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_B/Mean/Model'
 * '<S102>' : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_B/Mean/Model/Discrete Variable Time Delay'
 * '<S103>' : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_B/Mean value1/Model'
 * '<S104>' : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_B/Mean value1/Model/Discrete Variable Time Delay'
 * '<S105>' : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_C/Mean'
 * '<S106>' : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_C/Mean value1'
 * '<S107>' : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_C/Mean/Model'
 * '<S108>' : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_C/Mean/Model/Discrete Variable Time Delay'
 * '<S109>' : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_C/Mean value1/Model'
 * '<S110>' : 'Diesel_DER2/Exciter/Voltage Regulator/Sequence Analyzer/Fourier_C/Mean value1/Model/Discrete Variable Time Delay'
 * '<S111>' : 'Diesel_DER2/Mean/Model'
 * '<S112>' : 'Diesel_DER2/Mean/Model/Discrete Variable Time Delay'
 * '<S113>' : 'Diesel_DER2/Multimeter/Model'
 * '<S114>' : 'Diesel_DER2/Multimeter/StoreData'
 * '<S115>' : 'Diesel_DER2/Multimeter1/Model'
 * '<S116>' : 'Diesel_DER2/Multimeter1/StoreData'
 * '<S117>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model'
 * '<S118>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Measurements'
 * '<S119>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Mechanical model'
 * '<S120>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Machine measurements'
 * '<S121>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Synchronous Machine Discrete Model'
 * '<S122>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/abc to qd transformation'
 * '<S123>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/qd to abc transformation'
 * '<S124>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Machine measurements/Delta angle'
 * '<S125>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Machine measurements/PQ'
 * '<S126>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Machine measurements/Delta angle/Cartesian to Polar'
 * '<S127>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Machine measurements/Delta angle/Radians to Degrees'
 * '<S128>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Synchronous Machine Discrete Model/Electromagnetic Torque'
 * '<S129>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Synchronous Machine Discrete Model/Matrix W'
 * '<S130>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Synchronous Machine Discrete Model/Phimqd'
 * '<S131>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation'
 * '<S132>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Synchronous Machine Discrete Model/phi'
 * '<S133>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmd_sat'
 * '<S134>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmq_sat'
 * '<S135>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation/Update Matrix L'
 * '<S136>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmd_sat/Lad'
 * '<S137>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmd_sat/phimd'
 * '<S138>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmq_sat/Laq'
 * '<S139>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmq_sat/phimq'
 * '<S140>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Electrical model/Synchronous Machine Discrete Model/phi/Subsystem'
 * '<S141>' : 'Diesel_DER2/Synchronous  Generator 2 MVA 400V/Mechanical model/Delay Prediction'
 * '<S142>' : 'Diesel_DER2/Three-Phase Breaker/Breaker A'
 * '<S143>' : 'Diesel_DER2/Three-Phase Breaker/Breaker B'
 * '<S144>' : 'Diesel_DER2/Three-Phase Breaker/Breaker C'
 * '<S145>' : 'Diesel_DER2/Three-Phase Breaker/Stair Generator'
 * '<S146>' : 'Diesel_DER2/Three-Phase Breaker/Breaker A/Model'
 * '<S147>' : 'Diesel_DER2/Three-Phase Breaker/Breaker A/Model/Stair Generator'
 * '<S148>' : 'Diesel_DER2/Three-Phase Breaker/Breaker A/Model/Stair Generator/Model'
 * '<S149>' : 'Diesel_DER2/Three-Phase Breaker/Breaker B/Model'
 * '<S150>' : 'Diesel_DER2/Three-Phase Breaker/Breaker B/Model/Stair Generator'
 * '<S151>' : 'Diesel_DER2/Three-Phase Breaker/Breaker B/Model/Stair Generator/Model'
 * '<S152>' : 'Diesel_DER2/Three-Phase Breaker/Breaker C/Model'
 * '<S153>' : 'Diesel_DER2/Three-Phase Breaker/Breaker C/Model/Stair Generator'
 * '<S154>' : 'Diesel_DER2/Three-Phase Breaker/Breaker C/Model/Stair Generator/Model'
 * '<S155>' : 'Diesel_DER2/Three-Phase Breaker/Stair Generator/Model'
 * '<S156>' : 'Diesel_DER2/powergui/EquivalentModel1'
 * '<S157>' : 'Diesel_DER2/powergui/EquivalentModel1/Gates'
 * '<S158>' : 'Diesel_DER2/powergui/EquivalentModel1/Sources'
 * '<S159>' : 'Diesel_DER2/powergui/EquivalentModel1/Status'
 * '<S160>' : 'Diesel_DER2/powergui/EquivalentModel1/Yout'
 */
#endif                                 /* RTW_HEADER_Diesel_DER2_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
