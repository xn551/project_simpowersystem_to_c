/*
 * File: Wind_songweiwei.h
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

#ifndef RTW_HEADER_Wind_songweiwei_h_
#define RTW_HEADER_Wind_songweiwei_h_
#include <stddef.h>
#include <string.h>
#include <float.h>
#include <math.h>
#ifndef Wind_songweiwei_COMMON_INCLUDES_
# define Wind_songweiwei_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "simstruc.h"
#include "fixedpoint.h"
#endif                                 /* Wind_songweiwei_COMMON_INCLUDES_ */

#include "Wind_songweiwei_types.h"
#include "sfcn_bridge.h" 
#include "rtGetInf.h"
#include "rt_nonfinite.h"
#include "rt_defines.h"
#include "rt_look.h"
#include "rt_look1d.h"

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

/* Block signals for system '<S160>/Subsystem - pi//2 delay' */
typedef struct {
  real_T Fcn;                          /* '<S164>/Fcn' */
  real_T Fcn1;                         /* '<S164>/Fcn1' */
} B_Subsystempi2delay_Wind_song_T;

/* Block signals for system '<S160>/Subsystem1' */
typedef struct {
  real_T Fcn;                          /* '<S165>/Fcn' */
  real_T Fcn1;                         /* '<S165>/Fcn1' */
} B_Subsystem1_Wind_songweiwei_T;

/* Block signals (auto storage) */
typedef struct {
  real_T DataTypeConversion[6];        /* '<S123>/Data Type Conversion' */
  real_T ib[3];                        /* '<S333>/ib' */
  real_T Sum5[3];                      /* '<S10>/Sum5' */
  real_T StateSpace_o1[50];            /* '<S357>/State-Space' */
  real_T StateSpace_o2[19];            /* '<S357>/State-Space' */
  real_T Integ4;                       /* '<S154>/Integ4' */
  real_T Delay;                        /* '<S154>/Gain' */
  real_T SFunction;                    /* '<S159>/S-Function ' */
  real_T Integ4_n;                     /* '<S153>/Integ4' */
  real_T Delay_f;                      /* '<S153>/Gain' */
  real_T SFunction_j;                  /* '<S157>/S-Function ' */
  real_T Integ4_j;                     /* '<S168>/Integ4' */
  real_T Delay_n;                      /* '<S168>/Gain' */
  real_T SFunction_f;                  /* '<S173>/S-Function ' */
  real_T Integ4_m;                     /* '<S167>/Integ4' */
  real_T Delay_h;                      /* '<S167>/Gain' */
  real_T SFunction_m;                  /* '<S171>/S-Function ' */
  real_T Integ4_o;                     /* '<S251>/Integ4' */
  real_T Delay_o;                      /* '<S251>/Gain' */
  real_T SFunction_k;                  /* '<S253>/S-Function' */
  real_T Switch;                       /* '<S124>/Switch' */
  real_T RateLimiter1;                 /* '<S125>/Rate Limiter   1' */
  real_T dw;                           /* '<S329>/Rotor speed deviation (dw)' */
  real_T Product[6];                   /* '<S354>/Product' */
  real_T Linv[25];                     /* '<S341>/inversion' */
  real_T RLinv[25];                    /* '<S341>/Product1' */
  real_T Lmsatq;                       /* '<S344>/Lmq' */
  real_T Integ4_e;                     /* '<S235>/Integ4' */
  real_T Delay_e;                      /* '<S235>/Gain' */
  real_T SFunction_h;                  /* '<S237>/S-Function' */
  real_T Integ4_h;                     /* '<S238>/Integ4' */
  real_T Delay_l;                      /* '<S238>/Gain' */
  real_T SFunction_i;                  /* '<S240>/S-Function' */
  real_T MathFunction;                 /* '<S226>/Math Function' */
  real_T Fcn;                          /* '<S141>/Fcn' */
  real_T Fcn1;                         /* '<S141>/Fcn1' */
  real_T Fcn_c;                        /* '<S140>/Fcn' */
  real_T Fcn1_b;                       /* '<S140>/Fcn1' */
  real_T LookUpTable;                  /* '<S15>/Look-Up Table' */
  real_T DiscreteTimeIntegrator;       /* '<S12>/Discrete-Time Integrator' */
  real_T Switch2;                      /* '<S12>/Switch2' */
  real_T Switch3;                      /* '<S12>/Switch3' */
  real_T Product1[3];                  /* '<S13>/Product1' */
  real_T Product1_k[3];                /* '<S14>/Product1' */
  boolean_T LogicalOperator1;          /* '<S12>/Logical Operator1' */
  B_Subsystem1_Wind_songweiwei_T Subsystem1;/* '<S274>/Subsystem1' */
  B_Subsystempi2delay_Wind_song_T Subsystempi2delay;/* '<S274>/Subsystem - pi//2 delay' */
  B_Subsystem1_Wind_songweiwei_T Subsystem1_p;/* '<S268>/Subsystem1' */
  B_Subsystempi2delay_Wind_song_T Subsystempi2delay_jr;/* '<S268>/Subsystem - pi//2 delay' */
  B_Subsystem1_Wind_songweiwei_T Subsystem1_l;/* '<S262>/Subsystem1' */
  B_Subsystempi2delay_Wind_song_T Subsystempi2delay_j;/* '<S262>/Subsystem - pi//2 delay' */
  B_Subsystem1_Wind_songweiwei_T Subsystem1_l1;/* '<S254>/Subsystem1' */
  B_Subsystempi2delay_Wind_song_T Subsystempi2delay_p;/* '<S254>/Subsystem - pi//2 delay' */
  B_Subsystem1_Wind_songweiwei_T Subsystem1_lp;/* '<S241>/Subsystem1' */
  B_Subsystempi2delay_Wind_song_T Subsystempi2delay_dg;/* '<S241>/Subsystem - pi//2 delay' */
  B_Subsystem1_Wind_songweiwei_T Subsystem1_oh;/* '<S174>/Subsystem1' */
  B_Subsystempi2delay_Wind_song_T Subsystempi2delay_h;/* '<S174>/Subsystem - pi//2 delay' */
  B_Subsystem1_Wind_songweiwei_T Subsystem1_o;/* '<S160>/Subsystem1' */
  B_Subsystempi2delay_Wind_song_T Subsystempi2delay_k;/* '<S160>/Subsystem - pi//2 delay' */
} B_Wind_songweiwei_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T UnitDelay_DSTATE[6];          /* '<S353>/Unit Delay' */
  real_T UnitDelay6_DSTATE[3];         /* '<S67>/Unit Delay6' */
  real_T Rotorangledthetae_DSTATE;     /* '<S329>/Rotor angle dthetae' */
  real_T fluxes_DSTATE[5];             /* '<S342>/fluxes' */
  real_T StateSpace_DSTATE[45];        /* '<S357>/State-Space' */
  real_T UnitDelay3_DSTATE;            /* '<S67>/Unit Delay3' */
  real_T UnitDelay7_DSTATE;            /* '<S67>/Unit Delay7' */
  real_T dw_delay_DSTATE;              /* '<S351>/dw_delay' */
  real_T dw_predict_DSTATE;            /* '<S351>/dw_predict' */
  real_T DiscreteTimeIntegrator_DSTATE;/* '<S126>/Discrete-Time Integrator' */
  real_T Delay_x1_DSTATE;              /* '<S221>/Delay_x1' */
  real_T Delay_x2_DSTATE;              /* '<S221>/Delay_x2' */
  real_T Delay_x1_DSTATE_c[3];         /* '<S197>/Delay_x1' */
  real_T Delay_x2_DSTATE_e[3];         /* '<S197>/Delay_x2' */
  real_T theta_DSTATE;                 /* '<S329>/theta' */
  real_T Delay_x1_DSTATE_a;            /* '<S213>/Delay_x1' */
  real_T Delay_x2_DSTATE_p;            /* '<S213>/Delay_x2' */
  real_T Delay_x1_DSTATE_i[3];         /* '<S201>/Delay_x1' */
  real_T Delay_x2_DSTATE_m[3];         /* '<S201>/Delay_x2' */
  real_T Delay_x_DSTATE;               /* '<S261>/Delay_x' */
  real_T Delay_x1_DSTATE_ck;           /* '<S217>/Delay_x1' */
  real_T Delay_x2_DSTATE_ee;           /* '<S217>/Delay_x2' */
  real_T Delay_x1_DSTATE_m[3];         /* '<S189>/Delay_x1' */
  real_T Delay_x2_DSTATE_o[3];         /* '<S189>/Delay_x2' */
  real_T DiscreteTimeIntegrator_DSTATE_h;/* '<S225>/Discrete-Time Integrator' */
  real_T UnitDelay_DSTATE_p;           /* '<S225>/Unit Delay' */
  real_T ICi_ic_DSTATE;                /* '<S131>/IC = i_ic' */
  real_T DiscreteTimeIntegrator_DSTATE_g;/* '<S128>/Discrete-Time Integrator' */
  real_T Delay_x_DSTATE_m;             /* '<S134>/Delay_x' */
  real_T UnitDelay1_DSTATE;            /* '<S121>/Unit Delay1' */
  real_T DiscreteTimeIntegrator1_DSTATE;/* '<S121>/Discrete-Time Integrator1' */
  real_T Integ4_DSTATE;                /* '<S154>/Integ4' */
  real_T UnitDelay_DSTATE_f;           /* '<S158>/Unit Delay' */
  real_T UnitDelay_DSTATE_m;           /* '<S154>/Unit Delay' */
  real_T Integ4_DSTATE_n;              /* '<S153>/Integ4' */
  real_T UnitDelay_DSTATE_fd;          /* '<S156>/Unit Delay' */
  real_T UnitDelay_DSTATE_l;           /* '<S153>/Unit Delay' */
  real_T IC0_DSTATE;                   /* '<S131>/IC = 0' */
  real_T Delay_x1_DSTATE_b[3];         /* '<S193>/Delay_x1' */
  real_T Delay_x2_DSTATE_pn[3];        /* '<S193>/Delay_x2' */
  real_T DiscreteTimeIntegrator_DSTATE_n[2];/* '<S129>/Discrete-Time Integrator' */
  real_T Delay_x1_DSTATE_f;            /* '<S209>/Delay_x1' */
  real_T Delay_x2_DSTATE_d;            /* '<S209>/Delay_x2' */
  real_T Integ4_DSTATE_nw;             /* '<S168>/Integ4' */
  real_T UnitDelay_DSTATE_fz;          /* '<S172>/Unit Delay' */
  real_T UnitDelay_DSTATE_b;           /* '<S168>/Unit Delay' */
  real_T Integ4_DSTATE_l;              /* '<S167>/Integ4' */
  real_T UnitDelay_DSTATE_e;           /* '<S170>/Unit Delay' */
  real_T UnitDelay_DSTATE_d;           /* '<S167>/Unit Delay' */
  real_T Delay_x1_DSTATE_in;           /* '<S205>/Delay_x1' */
  real_T Delay_x2_DSTATE_a;            /* '<S205>/Delay_x2' */
  real_T Integ4_DSTATE_a;              /* '<S251>/Integ4' */
  real_T UnitDelay_DSTATE_n;           /* '<S252>/Unit Delay' */
  real_T UnitDelay1_DSTATE_p;          /* '<S251>/Unit Delay1' */
  real_T DiscreteDerivative_states;    /* '<S227>/Discrete Derivative ' */
  real_T DiscreteTimeIntegrator_DSTATE_k;/* '<S227>/Discrete-Time Integrator' */
  real_T Delay_x1_DSTATE_bs;           /* '<S247>/Delay_x1' */
  real_T Delay_x2_DSTATE_oc;           /* '<S247>/Delay_x2' */
  real_T UnitDelay2_DSTATE;            /* '<S67>/Unit Delay2' */
  real_T Delay_x_DSTATE_o;             /* '<S310>/Delay_x' */
  real_T DiscreteTimeIntegrator_DSTATE_o;/* '<S305>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTAT_kk;/* '<S304>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_b;/* '<S306>/Discrete-Time Integrator' */
  real_T Delay_x_DSTATE_i;             /* '<S309>/Delay_x' */
  real_T UnitDelay1_DSTATE_n;          /* '<S67>/Unit Delay1' */
  real_T UnitDelay4_DSTATE;            /* '<S67>/Unit Delay4' */
  real_T DiscreteTimeIntegrator_DSTAT_bh;/* '<S69>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator1_DSTAT_l;/* '<S69>/Discrete-Time Integrator1' */
  real_T voltages_DSTATE[5];           /* '<S342>/voltages' */
  real_T Rotorspeeddeviationdw_DSTATE; /* '<S329>/Rotor speed deviation (dw)' */
  real_T DiscreteTimeIntegrator_DSTAT_b4[6];/* '<S354>/Discrete-Time Integrator' */
  real_T UnitDelay_DSTATE_fv[6];       /* '<S354>/Unit Delay' */
  real_T Lmd_sat_DSTATE;               /* '<S343>/Lmd_sat' */
  real_T Lmq_sat_DSTATE;               /* '<S344>/Lmq_sat' */
  real_T Integ4_DSTATE_nn;             /* '<S235>/Integ4' */
  real_T UnitDelay_DSTATE_i;           /* '<S236>/Unit Delay' */
  real_T UnitDelay1_DSTATE_f;          /* '<S235>/Unit Delay1' */
  real_T Integ4_DSTATE_c;              /* '<S238>/Integ4' */
  real_T UnitDelay_DSTATE_g;           /* '<S239>/Unit Delay' */
  real_T UnitDelay1_DSTATE_c;          /* '<S238>/Unit Delay1' */
  real_T DiscreteTimeIntegrator_DSTATE_l;/* '<S12>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator1_DSTAT_a;/* '<S16>/Discrete-Time Integrator1' */
  real_T UnitDelay_DSTATE_gd;          /* '<S16>/Unit Delay' */
  real_T DiscreteDerivative_tmp;       /* '<S227>/Discrete Derivative ' */
  real_T PrevY;                        /* '<S225>/Rate Limiter' */
  real_T PrevY_m;                      /* '<S125>/Rate Limiter   1' */
  real_T LastMajorTime;                /* '<S125>/Rate Limiter   1' */
  real_T inversion_DWORK4[25];         /* '<S350>/inversion' */
  real_T inversion_DWORK4_h[25];       /* '<S341>/inversion' */
  real_T SFunction_RWORK;              /* '<S253>/S-Function' */
  real_T SFunction_RWORK_e;            /* '<S237>/S-Function' */
  real_T SFunction_RWORK_g;            /* '<S240>/S-Function' */
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
  } StateSpace_PWORK;                  /* '<S357>/State-Space' */

  struct {
    void *uBuffers;
  } SFunction_PWORK;                   /* '<S159>/S-Function ' */

  struct {
    void *uBuffers;
  } SFunction_PWORK_c;                 /* '<S157>/S-Function ' */

  struct {
    void *uBuffers;
  } SFunction_PWORK_l;                 /* '<S173>/S-Function ' */

  struct {
    void *uBuffers;
  } SFunction_PWORK_p;                 /* '<S171>/S-Function ' */

  void *SFunction_PWORK_g;             /* '<S253>/S-Function' */
  void *SFunction_PWORK_ge;            /* '<S237>/S-Function' */
  void *SFunction_PWORK_cs;            /* '<S240>/S-Function' */
  int_T StateSpace_IWORK[5];           /* '<S357>/State-Space' */
  struct {
    int_T indBeg;
    int_T indEnd;
    int_T bufSz;
    int_T maxDiscrDelay;
  } SFunction_IWORK;                   /* '<S159>/S-Function ' */

  struct {
    int_T indBeg;
    int_T indEnd;
    int_T bufSz;
    int_T maxDiscrDelay;
  } SFunction_IWORK_f;                 /* '<S157>/S-Function ' */

  struct {
    int_T indBeg;
    int_T indEnd;
    int_T bufSz;
    int_T maxDiscrDelay;
  } SFunction_IWORK_c;                 /* '<S173>/S-Function ' */

  struct {
    int_T indBeg;
    int_T indEnd;
    int_T bufSz;
    int_T maxDiscrDelay;
  } SFunction_IWORK_k;                 /* '<S171>/S-Function ' */

  int_T SFunction_IWORK_fk;            /* '<S253>/S-Function' */
  int_T SFunction_IWORK_h;             /* '<S237>/S-Function' */
  int_T SFunction_IWORK_n;             /* '<S240>/S-Function' */
  int8_T DiscreteTimeIntegrator_PrevRese[6];/* '<S354>/Discrete-Time Integrator' */
  uint8_T DiscreteTimeIntegrator1_IC_LOAD;/* '<S121>/Discrete-Time Integrator1' */
  uint8_T Integ4_SYSTEM_ENABLE;        /* '<S154>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_f;      /* '<S153>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_l;      /* '<S168>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_h;      /* '<S167>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_a;      /* '<S251>/Integ4' */
  uint8_T Rotorspeeddeviationdw_SYSTEM_EN;/* '<S329>/Rotor speed deviation (dw)' */
  uint8_T Integ4_SYSTEM_ENABLE_e;      /* '<S235>/Integ4' */
  uint8_T Integ4_SYSTEM_ENABLE_m;      /* '<S238>/Integ4' */
  boolean_T Tail_MODE;                 /* '<S353>/Tail' */
  boolean_T Signalgenerator_MODE;      /* '<S10>/Signal generator' */
  boolean_T HarmonicGenerator_MODE;    /* '<S10>/Harmonic Generator' */
  boolean_T AutomaticGainControl_MODE; /* '<S225>/Automatic Gain Control' */
} DW_Wind_songweiwei_T;

/* External outputs (root outports fed by signals with auto storage) */
typedef struct {
  real_T Out_Vabc_B575[3];             /* '<Root>/Out_Vabc_B575' */
  real_T Out_labc_B575[3];             /* '<Root>/Out_labc_B575' */
  real_T Out_P;                        /* '<Root>/Out_P' */
  real_T Out_Q;                        /* '<Root>/Out_Q' */
  real_T Out_Vdc;                      /* '<Root>/Out_Vdc' */
  real_T Out_wr;                       /* '<Root>/Out_wr' */
  real_T Out_Vabc_B25[3];              /* '<Root>/Out_Vabc_B25' */
  real_T Out_labc_B25[3];              /* '<Root>/Out_labc_B25' */
} ExtY_Wind_songweiwei_T;

/* Parameters for system: '<S160>/Subsystem - pi//2 delay' */
struct P_Subsystempi2delay_Wind_song_T_ {
  real_T dq_Y0[2];                     /* Expression: [0,0]
                                        * Referenced by: '<S164>/dq'
                                        */
};

/* Parameters for system: '<S160>/Subsystem1' */
struct P_Subsystem1_Wind_songweiwei_T_ {
  real_T dq_Y0[2];                     /* Expression: [0,0]
                                        * Referenced by: '<S165>/dq'
                                        */
};

/* Parameters (auto storage) */
struct P_Wind_songweiwei_T_ {
  real_T Ts;                           /* Variable: Ts
                                        * Referenced by:
                                        *   '<S124>/Nb_of_samples_by_cycle'
                                        *   '<S124>/Ts'
                                        *   '<S131>/Gain'
                                        *   '<S153>/Gain'
                                        *   '<S154>/Gain'
                                        *   '<S167>/Gain'
                                        *   '<S168>/Gain'
                                        *   '<S251>/Gain'
                                        *   '<S235>/Gain'
                                        *   '<S238>/Gain'
                                        */
  real_T AlphaBetaZerotodq0_Alignment; /* Mask Parameter: AlphaBetaZerotodq0_Alignment
                                        * Referenced by: '<S241>/Constant'
                                        */
  real_T AlphaBetaZerotodq0_Alignment_i;/* Mask Parameter: AlphaBetaZerotodq0_Alignment_i
                                         * Referenced by: '<S274>/Constant'
                                         */
  real_T AlphaBetaZerotodq0_Alignment_n;/* Mask Parameter: AlphaBetaZerotodq0_Alignment_n
                                         * Referenced by: '<S262>/Constant'
                                         */
  real_T AlphaBetaZerotodq0_Alignment_iw;/* Mask Parameter: AlphaBetaZerotodq0_Alignment_iw
                                          * Referenced by: '<S160>/Constant'
                                          */
  real_T AlphaBetaZerotodq0_Alignment_m;/* Mask Parameter: AlphaBetaZerotodq0_Alignment_m
                                         * Referenced by: '<S268>/Constant'
                                         */
  real_T AlphaBetaZerotodq0_Alignment_l;/* Mask Parameter: AlphaBetaZerotodq0_Alignment_l
                                         * Referenced by: '<S174>/Constant'
                                         */
  real_T dq0toAlphaBetaZero_Alignment; /* Mask Parameter: dq0toAlphaBetaZero_Alignment
                                        * Referenced by: '<S137>/Constant'
                                        */
  real_T AlphaBetaZerotodq0_Alignment_f;/* Mask Parameter: AlphaBetaZerotodq0_Alignment_f
                                         * Referenced by: '<S254>/Constant'
                                         */
  real_T DriveTrain_D_mutual;          /* Mask Parameter: DriveTrain_D_mutual
                                        * Referenced by: '<S69>/Mutual damping'
                                        */
  real_T Discrete3phasePLLDrivenPositive;/* Mask Parameter: Discrete3phasePLLDrivenPositive
                                          * Referenced by:
                                          *   '<S153>/Step'
                                          *   '<S154>/Step'
                                          *   '<S167>/Step'
                                          *   '<S168>/Step'
                                          */
  real_T PulseGenerator_Freq_sawtooth; /* Mask Parameter: PulseGenerator_Freq_sawtooth
                                        * Referenced by: '<S124>/Nb_of_samples_by_cycle'
                                        */
  real_T DriveTrain_H_WT;              /* Mask Parameter: DriveTrain_H_WT
                                        * Referenced by: '<S69>/1_2H_WT'
                                        */
  real_T u0kV_HarmonicGeneration;      /* Mask Parameter: u0kV_HarmonicGeneration
                                        * Referenced by: '<S10>/valp_nom7'
                                        */
  real_T DiscretePIController_Init;    /* Mask Parameter: DiscretePIController_Init
                                        * Referenced by: '<S126>/Discrete-Time Integrator'
                                        */
  real_T DiscretePIController_Init_n;  /* Mask Parameter: DiscretePIController_Init_n
                                        * Referenced by: '<S128>/Discrete-Time Integrator'
                                        */
  real_T DiscretePIController1_Init;   /* Mask Parameter: DiscretePIController1_Init
                                        * Referenced by: '<S129>/Discrete-Time Integrator'
                                        */
  real_T Discrete_Init;                /* Mask Parameter: Discrete_Init
                                        * Referenced by: '<S227>/Discrete-Time Integrator'
                                        */
  real_T DiscretePIController1_Init_i; /* Mask Parameter: DiscretePIController1_Init_i
                                        * Referenced by: '<S305>/Discrete-Time Integrator'
                                        */
  real_T DiscretePIController_Init_a;  /* Mask Parameter: DiscretePIController_Init_a
                                        * Referenced by: '<S304>/Discrete-Time Integrator'
                                        */
  real_T DiscretePIController2_Init;   /* Mask Parameter: DiscretePIController2_Init
                                        * Referenced by: '<S306>/Discrete-Time Integrator'
                                        */
  real_T Discrete_Kd;                  /* Mask Parameter: Discrete_Kd
                                        * Referenced by: '<S227>/Discrete Derivative '
                                        */
  real_T DiscretePIController_Ki;      /* Mask Parameter: DiscretePIController_Ki
                                        * Referenced by: '<S126>/Kp5'
                                        */
  real_T DiscretePIController_Ki_f;    /* Mask Parameter: DiscretePIController_Ki_f
                                        * Referenced by: '<S128>/Kp5'
                                        */
  real_T DiscretePIController1_Ki;     /* Mask Parameter: DiscretePIController1_Ki
                                        * Referenced by: '<S129>/Kp5'
                                        */
  real_T DiscretePIController_Ki_p;    /* Mask Parameter: DiscretePIController_Ki_p
                                        * Referenced by: '<S304>/Kp5'
                                        */
  real_T DiscretePIController1_Ki_j;   /* Mask Parameter: DiscretePIController1_Ki_j
                                        * Referenced by: '<S305>/Kp5'
                                        */
  real_T DiscretePIController2_Ki;     /* Mask Parameter: DiscretePIController2_Ki
                                        * Referenced by: '<S306>/Kp5'
                                        */
  real_T DiscretePIController_Kp;      /* Mask Parameter: DiscretePIController_Kp
                                        * Referenced by: '<S126>/Kp4'
                                        */
  real_T DiscretePIController_Kp_h;    /* Mask Parameter: DiscretePIController_Kp_h
                                        * Referenced by: '<S128>/Kp4'
                                        */
  real_T DiscretePIController1_Kp;     /* Mask Parameter: DiscretePIController1_Kp
                                        * Referenced by: '<S129>/Kp4'
                                        */
  real_T Discrete_Kp;                  /* Mask Parameter: Discrete_Kp
                                        * Referenced by: '<S227>/Kp4'
                                        */
  real_T DiscretePIController1_Kp_p;   /* Mask Parameter: DiscretePIController1_Kp_p
                                        * Referenced by: '<S305>/Kp4'
                                        */
  real_T DiscretePIController_Kp_c;    /* Mask Parameter: DiscretePIController_Kp_c
                                        * Referenced by: '<S304>/Kp4'
                                        */
  real_T DiscretePIController2_Kp;     /* Mask Parameter: DiscretePIController2_Kp
                                        * Referenced by: '<S306>/Kp4'
                                        */
  real_T DriveTrain_Ksh;               /* Mask Parameter: DriveTrain_Ksh
                                        * Referenced by:
                                        *   '<S69>/Discrete-Time Integrator1'
                                        *   '<S69>/Stiffness'
                                        */
  real_T HarmonicAgeneration_Mag_Harmo;/* Mask Parameter: HarmonicAgeneration_Mag_Harmo
                                        * Referenced by: '<S13>/Phase_Harmo1'
                                        */
  real_T HarmonicBgeneration_Mag_Harmo;/* Mask Parameter: HarmonicBgeneration_Mag_Harmo
                                        * Referenced by: '<S14>/Phase_Harmo1'
                                        */
  real_T PWMGenerator2Level_MinMax[2]; /* Mask Parameter: PWMGenerator2Level_MinMax
                                        * Referenced by: '<S123>/Constant10'
                                        */
  real_T HarmonicAgeneration_Phase_Harmo;/* Mask Parameter: HarmonicAgeneration_Phase_Harmo
                                          * Referenced by: '<S13>/Phase_Harmo'
                                          */
  real_T HarmonicBgeneration_Phase_Harmo;/* Mask Parameter: HarmonicBgeneration_Phase_Harmo
                                          * Referenced by: '<S14>/Phase_Harmo'
                                          */
  real_T WindTurbine_Prated;           /* Mask Parameter: WindTurbine_Prated
                                        * Referenced by: '<S80>/->pu'
                                        */
  real_T WindTurbineType4_Resistance;  /* Mask Parameter: WindTurbineType4_Resistance
                                        * Referenced by:
                                        *   '<S145>/Rs2'
                                        *   '<S145>/Rs4'
                                        */
  real_T HarmonicAgeneration_Seq_Harmo;/* Mask Parameter: HarmonicAgeneration_Seq_Harmo
                                        * Referenced by: '<S13>/Phase_Harmo2'
                                        */
  real_T HarmonicBgeneration_Seq_Harmo;/* Mask Parameter: HarmonicBgeneration_Seq_Harmo
                                        * Referenced by: '<S14>/Phase_Harmo2'
                                        */
  real_T Tail_Tf;                      /* Mask Parameter: Tail_Tf
                                        * Referenced by:
                                        *   '<S354>/2'
                                        *   '<S356>/Constant2'
                                        *   '<S356>/-0.9//Tf'
                                        */
  real_T VariationSubSystem_Toff_Variati;/* Mask Parameter: VariationSubSystem_Toff_Variati
                                          * Referenced by: '<S16>/Step1'
                                          */
  real_T VariationSubSystem_Ton_Variatio;/* Mask Parameter: VariationSubSystem_Ton_Variatio
                                          * Referenced by: '<S16>/Step'
                                          */
  real_T Tail_Tt;                      /* Mask Parameter: Tail_Tt
                                        * Referenced by:
                                        *   '<S354>/2'
                                        *   '<S356>/Constant2'
                                        *   '<S356>/0.1//Tt'
                                        */
  real_T u0kV_VariationEntity;         /* Mask Parameter: u0kV_VariationEntity
                                        * Referenced by:
                                        *   '<S10>/valp_nom3'
                                        *   '<S12>/valp_nom3'
                                        */
  real_T VariationSubSystem_VariationFre;/* Mask Parameter: VariationSubSystem_VariationFre
                                          * Referenced by: '<S16>/valp_nom9'
                                          */
  real_T VariationSubSystem_VariationMag;/* Mask Parameter: VariationSubSystem_VariationMag
                                          * Referenced by: '<S16>/valp_nom8'
                                          */
  real_T u0kV_VariationPhaseA;         /* Mask Parameter: u0kV_VariationPhaseA
                                        * Referenced by: '<S10>/SinglePhase'
                                        */
  real_T VariationSubSystem_VariationRat;/* Mask Parameter: VariationSubSystem_VariationRat
                                          * Referenced by: '<S16>/valp_nom7'
                                          */
  real_T VariationSubSystem_VariationSte;/* Mask Parameter: VariationSubSystem_VariationSte
                                          * Referenced by: '<S16>/valp_nom6'
                                          */
  real_T u0kV_VariationType;           /* Mask Parameter: u0kV_VariationType
                                        * Referenced by: '<S12>/valp_nom5'
                                        */
  real_T DiscreteVariableFrequencyMeanva;/* Mask Parameter: DiscreteVariableFrequencyMeanva
                                          * Referenced by: '<S154>/Unit Delay'
                                          */
  real_T DiscreteVariableFrequencyMean_m;/* Mask Parameter: DiscreteVariableFrequencyMean_m
                                          * Referenced by: '<S153>/Unit Delay'
                                          */
  real_T DiscreteVariableFrequencyMean_j;/* Mask Parameter: DiscreteVariableFrequencyMean_j
                                          * Referenced by: '<S168>/Unit Delay'
                                          */
  real_T DiscreteVariableFrequencyMean_l;/* Mask Parameter: DiscreteVariableFrequencyMean_l
                                          * Referenced by: '<S167>/Unit Delay'
                                          */
  real_T CompareToConstant_const;      /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S243>/Constant'
                                        */
  real_T CompareToConstant1_const;     /* Mask Parameter: CompareToConstant1_const
                                        * Referenced by: '<S244>/Constant'
                                        */
  real_T CompareToConstant_const_o;    /* Mask Parameter: CompareToConstant_const_o
                                        * Referenced by: '<S276>/Constant'
                                        */
  real_T CompareToConstant1_const_o;   /* Mask Parameter: CompareToConstant1_const_o
                                        * Referenced by: '<S277>/Constant'
                                        */
  real_T CompareToConstant_const_c;    /* Mask Parameter: CompareToConstant_const_c
                                        * Referenced by: '<S264>/Constant'
                                        */
  real_T CompareToConstant1_const_b;   /* Mask Parameter: CompareToConstant1_const_b
                                        * Referenced by: '<S265>/Constant'
                                        */
  real_T CompareToConstant_const_g;    /* Mask Parameter: CompareToConstant_const_g
                                        * Referenced by: '<S162>/Constant'
                                        */
  real_T CompareToConstant1_const_g;   /* Mask Parameter: CompareToConstant1_const_g
                                        * Referenced by: '<S163>/Constant'
                                        */
  real_T CompareToConstant_const_i;    /* Mask Parameter: CompareToConstant_const_i
                                        * Referenced by: '<S270>/Constant'
                                        */
  real_T CompareToConstant1_const_k;   /* Mask Parameter: CompareToConstant1_const_k
                                        * Referenced by: '<S271>/Constant'
                                        */
  real_T CompareToConstant_const_gf;   /* Mask Parameter: CompareToConstant_const_gf
                                        * Referenced by: '<S176>/Constant'
                                        */
  real_T CompareToConstant1_const_gp;  /* Mask Parameter: CompareToConstant1_const_gp
                                        * Referenced by: '<S177>/Constant'
                                        */
  real_T CompareToConstant_const_ia;   /* Mask Parameter: CompareToConstant_const_ia
                                        * Referenced by: '<S138>/Constant'
                                        */
  real_T CompareToConstant1_const_o0;  /* Mask Parameter: CompareToConstant1_const_o0
                                        * Referenced by: '<S139>/Constant'
                                        */
  real_T CompareToConstant_const_oc;   /* Mask Parameter: CompareToConstant_const_oc
                                        * Referenced by: '<S256>/Constant'
                                        */
  real_T CompareToConstant1_const_j;   /* Mask Parameter: CompareToConstant1_const_j
                                        * Referenced by: '<S257>/Constant'
                                        */
  real_T HarmonicAgeneration_n_Harmo;  /* Mask Parameter: HarmonicAgeneration_n_Harmo
                                        * Referenced by: '<S13>/Gain1'
                                        */
  real_T HarmonicBgeneration_n_Harmo;  /* Mask Parameter: HarmonicBgeneration_n_Harmo
                                        * Referenced by: '<S14>/Gain1'
                                        */
  real_T DriveTrain_torque0;           /* Mask Parameter: DriveTrain_torque0
                                        * Referenced by: '<S69>/Discrete-Time Integrator1'
                                        */
  real_T DriveTrain_w_wt0;             /* Mask Parameter: DriveTrain_w_wt0
                                        * Referenced by: '<S69>/Discrete-Time Integrator'
                                        */
  real_T DriveTrain_wbase;             /* Mask Parameter: DriveTrain_wbase
                                        * Referenced by: '<S69>/wbase'
                                        */
  real_T Negativesequence_Value[3];    /* Expression: [0 2*pi/3 -2*pi/3]
                                        * Referenced by: '<S13>/Negative-sequence'
                                        */
  real_T Positivesequence_Value[3];    /* Expression: [0 -2*pi/3 2*pi/3]
                                        * Referenced by: '<S13>/Positive-sequence'
                                        */
  real_T Zerosequence_Value[3];        /* Expression: [0 0 0]
                                        * Referenced by: '<S13>/Zero-sequence'
                                        */
  real_T Negativesequence_Value_m[3];  /* Expression: [0 2*pi/3 -2*pi/3]
                                        * Referenced by: '<S14>/Negative-sequence'
                                        */
  real_T Positivesequence_Value_k[3];  /* Expression: [0 -2*pi/3 2*pi/3]
                                        * Referenced by: '<S14>/Positive-sequence'
                                        */
  real_T Zerosequence_Value_b[3];      /* Expression: [0 0 0]
                                        * Referenced by: '<S14>/Zero-sequence'
                                        */
  real_T Out1_Y0;                      /* Expression: 0
                                        * Referenced by: '<S11>/Out1'
                                        */
  real_T Out2_Y0;                      /* Expression: 0
                                        * Referenced by: '<S11>/Out2'
                                        */
  real_T Gain3_Gain;                   /* Expression: pi/180
                                        * Referenced by: '<S13>/Gain3'
                                        */
  real_T valp_nom2_Value;              /* Expression: 1
                                        * Referenced by: '<S13>/valp_nom2'
                                        */
  real_T Step_Time;                    /* Expression: Ton_Harmo
                                        * Referenced by: '<S11>/Step'
                                        */
  real_T Step_Y0;                      /* Expression: 0
                                        * Referenced by: '<S11>/Step'
                                        */
  real_T Step_YFinal;                  /* Expression: 1
                                        * Referenced by: '<S11>/Step'
                                        */
  real_T Step1_Time;                   /* Expression: Toff_Harmo
                                        * Referenced by: '<S11>/Step1'
                                        */
  real_T Step1_Y0;                     /* Expression: 0
                                        * Referenced by: '<S11>/Step1'
                                        */
  real_T Step1_YFinal;                 /* Expression: -1
                                        * Referenced by: '<S11>/Step1'
                                        */
  real_T Gain3_Gain_p;                 /* Expression: pi/180
                                        * Referenced by: '<S14>/Gain3'
                                        */
  real_T valp_nom2_Value_a;            /* Expression: 1
                                        * Referenced by: '<S14>/valp_nom2'
                                        */
  real_T Constant5_Value;              /* Expression: 0
                                        * Referenced by: '<S12>/Constant5'
                                        */
  real_T Constant4_Value;              /* Expression: 0
                                        * Referenced by: '<S12>/Constant4'
                                        */
  real_T Gain4_Gain;                   /* Expression: pi/180
                                        * Referenced by: '<S12>/Gain4'
                                        */
  real_T Constant1_Value;              /* Expression: 0
                                        * Referenced by: '<S12>/Constant1'
                                        */
  real_T Constant1_Value_h;            /* Expression: 0
                                        * Referenced by: '<S16>/Constant1'
                                        */
  real_T Constant5_Value_h;            /* Expression: 0
                                        * Referenced by: '<S16>/Constant5'
                                        */
  real_T Gain1_Gain;                   /* Expression: 2*pi
                                        * Referenced by: '<S16>/Gain1'
                                        */
  real_T Constant4_Value_d;            /* Expression: 0
                                        * Referenced by: '<S16>/Constant4'
                                        */
  real_T Constant2_Value;              /* Expression: 1
                                        * Referenced by: '<S16>/Constant2'
                                        */
  real_T Constant_Value;               /* Expression: 2
                                        * Referenced by: '<S16>/Constant'
                                        */
  real_T timer_Y0;                     /* Expression: 0
                                        * Referenced by: '<S12>/timer'
                                        */
  real_T magnitude_Y0;                 /* Expression: 0
                                        * Referenced by: '<S12>/magnitude'
                                        */
  real_T frequency_Y0;                 /* Expression: 0
                                        * Referenced by: '<S12>/frequency'
                                        */
  real_T phase_Y0;                     /* Expression: 0
                                        * Referenced by: '<S12>/phase'
                                        */
  real_T LookUpTable_XData[7];         /* Expression: tv
                                        * Referenced by: '<S15>/Look-Up Table'
                                        */
  real_T LookUpTable_YData[7];         /* Expression: opv
                                        * Referenced by: '<S15>/Look-Up Table'
                                        */
  real_T Constant_Value_f;             /* Expression: 2
                                        * Referenced by: '<S12>/Constant'
                                        */
  real_T Constant2_Value_m;            /* Expression: 3
                                        * Referenced by: '<S12>/Constant2'
                                        */
  real_T Constant3_Value;              /* Expression: 4
                                        * Referenced by: '<S12>/Constant3'
                                        */
  real_T Constant6_Value;              /* Expression: 4
                                        * Referenced by: '<S12>/Constant6'
                                        */
  real_T DiscreteTimeIntegrator_gainval;/* Computed Parameter: DiscreteTimeIntegrator_gainval
                                         * Referenced by: '<S12>/Discrete-Time Integrator'
                                         */
  real_T DiscreteTimeIntegrator_IC;    /* Expression: 0
                                        * Referenced by: '<S12>/Discrete-Time Integrator'
                                        */
  real_T Step1_Y0_f;                   /* Expression: 1
                                        * Referenced by: '<S16>/Step1'
                                        */
  real_T Step1_YFinal_m;               /* Expression: 0
                                        * Referenced by: '<S16>/Step1'
                                        */
  real_T Constant3_Value_f;            /* Expression: 1
                                        * Referenced by: '<S16>/Constant3'
                                        */
  real_T Step_Y0_i;                    /* Expression: 0
                                        * Referenced by: '<S16>/Step'
                                        */
  real_T Step_YFinal_b;                /* Expression: 1
                                        * Referenced by: '<S16>/Step'
                                        */
  real_T DiscreteTimeIntegrator1_gainval;/* Computed Parameter: DiscreteTimeIntegrator1_gainval
                                          * Referenced by: '<S16>/Discrete-Time Integrator1'
                                          */
  real_T DiscreteTimeIntegrator1_IC;   /* Expression: 0
                                        * Referenced by: '<S16>/Discrete-Time Integrator1'
                                        */
  real_T UnitDelay_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S16>/Unit Delay'
                                        */
  real_T Switch2_Threshold;            /* Expression: 0.5
                                        * Referenced by: '<S16>/Switch2'
                                        */
  real_T Gain2_Gain;                   /* Expression: 2*pi
                                        * Referenced by: '<S12>/Gain2'
                                        */
  real_T Switch_Threshold;             /* Expression: 0.5
                                        * Referenced by: '<S16>/Switch'
                                        */
  real_T Flux_ref_Value;               /* Expression: 1
                                        * Referenced by: '<S120>/Flux_ref'
                                        */
  real_T avoiddivisionby0_UpperSat;    /* Expression: inf
                                        * Referenced by: '<S120>/avoid division by 0'
                                        */
  real_T avoiddivisionby0_LowerSat;    /* Expression: 1e-6
                                        * Referenced by: '<S120>/avoid division by 0'
                                        */
  real_T Flux_ref_Value_c;             /* Expression: 1
                                        * Referenced by: '<S120>/Flux_ref '
                                        */
  real_T alpha_beta_Y0[2];             /* Expression: [0,0]
                                        * Referenced by: '<S140>/alpha_beta'
                                        */
  real_T alpha_beta_Y0_o[2];           /* Expression: [0,0]
                                        * Referenced by: '<S141>/alpha_beta'
                                        */
  real_T Gain1_Gain_k;                 /* Expression: 0.5
                                        * Referenced by: '<S156>/Gain1'
                                        */
  real_T Gain1_Gain_g;                 /* Expression: 0.5
                                        * Referenced by: '<S158>/Gain1'
                                        */
  real_T Gain1_Gain_m;                 /* Expression: 0.5
                                        * Referenced by: '<S170>/Gain1'
                                        */
  real_T Gain1_Gain_e;                 /* Expression: 0.5
                                        * Referenced by: '<S172>/Gain1'
                                        */
  real_T Gain1_Gain_me;                /* Expression: 0.5
                                        * Referenced by: '<S236>/Gain1'
                                        */
  real_T Gain1_Gain_o;                 /* Expression: 0.5
                                        * Referenced by: '<S239>/Gain1'
                                        */
  real_T Gain_Y0;                      /* Expression: [1]
                                        * Referenced by: '<S226>/Gain'
                                        */
  real_T Gain3_Gain_l[9];              /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
                                        * Referenced by: '<S242>/Gain3'
                                        */
  real_T Gain1_Gain_ms;                /* Expression: 2/3
                                        * Referenced by: '<S242>/Gain1'
                                        */
  real_T Integ4_gainval;               /* Computed Parameter: Integ4_gainval
                                        * Referenced by: '<S235>/Integ4'
                                        */
  real_T Integ4_IC;                    /* Expression: 0
                                        * Referenced by: '<S235>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSat;/* Expression: 1e6
                                         * Referenced by: '<S235>/To avoid division  by zero'
                                         */
  real_T Toavoiddivisionbyzero_LowerSat;/* Expression: eps
                                         * Referenced by: '<S235>/To avoid division  by zero'
                                         */
  real_T SFunction_P1_Size[2];         /* Computed Parameter: SFunction_P1_Size
                                        * Referenced by: '<S237>/S-Function'
                                        */
  real_T SFunction_P1;                 /* Expression: MaxDelay
                                        * Referenced by: '<S237>/S-Function'
                                        */
  real_T SFunction_P2_Size[2];         /* Computed Parameter: SFunction_P2_Size
                                        * Referenced by: '<S237>/S-Function'
                                        */
  real_T SFunction_P2;                 /* Expression: Ts
                                        * Referenced by: '<S237>/S-Function'
                                        */
  real_T SFunction_P3_Size[2];         /* Computed Parameter: SFunction_P3_Size
                                        * Referenced by: '<S237>/S-Function'
                                        */
  real_T SFunction_P3;                 /* Expression: InitialValue
                                        * Referenced by: '<S237>/S-Function'
                                        */
  real_T SFunction_P4_Size[2];         /* Computed Parameter: SFunction_P4_Size
                                        * Referenced by: '<S237>/S-Function'
                                        */
  real_T SFunction_P4;                 /* Expression: DFT
                                        * Referenced by: '<S237>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_b; /* Expression: 0
                                        * Referenced by: '<S236>/Unit Delay'
                                        */
  real_T Constant_Value_g;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S235>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition;  /* Expression: sps.Vinit
                                        * Referenced by: '<S235>/Unit Delay1'
                                        */
  real_T Integ4_gainval_j;             /* Computed Parameter: Integ4_gainval_j
                                        * Referenced by: '<S238>/Integ4'
                                        */
  real_T Integ4_IC_k;                  /* Expression: 0
                                        * Referenced by: '<S238>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_f;/* Expression: 1e6
                                          * Referenced by: '<S238>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_m;/* Expression: eps
                                          * Referenced by: '<S238>/To avoid division  by zero'
                                          */
  real_T SFunction_P1_Size_j[2];       /* Computed Parameter: SFunction_P1_Size_j
                                        * Referenced by: '<S240>/S-Function'
                                        */
  real_T SFunction_P1_h;               /* Expression: MaxDelay
                                        * Referenced by: '<S240>/S-Function'
                                        */
  real_T SFunction_P2_Size_j[2];       /* Computed Parameter: SFunction_P2_Size_j
                                        * Referenced by: '<S240>/S-Function'
                                        */
  real_T SFunction_P2_n;               /* Expression: Ts
                                        * Referenced by: '<S240>/S-Function'
                                        */
  real_T SFunction_P3_Size_e[2];       /* Computed Parameter: SFunction_P3_Size_e
                                        * Referenced by: '<S240>/S-Function'
                                        */
  real_T SFunction_P3_d;               /* Expression: InitialValue
                                        * Referenced by: '<S240>/S-Function'
                                        */
  real_T SFunction_P4_Size_c[2];       /* Computed Parameter: SFunction_P4_Size_c
                                        * Referenced by: '<S240>/S-Function'
                                        */
  real_T SFunction_P4_e;               /* Expression: DFT
                                        * Referenced by: '<S240>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_n; /* Expression: 0
                                        * Referenced by: '<S239>/Unit Delay'
                                        */
  real_T Constant_Value_d;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S238>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_f;/* Expression: sps.Vinit
                                        * Referenced by: '<S238>/Unit Delay1'
                                        */
  real_T Saturation_UpperSat;          /* Expression: inf
                                        * Referenced by: '<S226>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: eps
                                        * Referenced by: '<S226>/Saturation'
                                        */
  real_T Gain1_Gain_et;                /* Expression: 0.5
                                        * Referenced by: '<S252>/Gain1'
                                        */
  real_T Constant4_Value_c;            /* Expression: 0
                                        * Referenced by: '<S124>/Constant4'
                                        */
  real_T Constant3_Value_m;            /* Expression: 1
                                        * Referenced by: '<S124>/Constant3'
                                        */
  real_T wref_Value;                   /* Expression: 1
                                        * Referenced by: '<S125>/wref'
                                        */
  real_T Constant4_Value_k[25];        /* Expression: SM.Linv
                                        * Referenced by: '<S331>/Constant4'
                                        */
  real_T u3_Value;                     /* Expression: SM.Lmq
                                        * Referenced by: '<S341>/u3'
                                        */
  real_T Constant1_Value_n;            /* Expression: SM.Lmsatd(1)
                                        * Referenced by: '<S343>/Constant1'
                                        */
  real_T Constant1_Value_e;            /* Expression: SM.Lmsatq(1)
                                        * Referenced by: '<S344>/Constant1'
                                        */
  real_T Ll_q_Gain[2];                 /* Expression: SM.One_Llq
                                        * Referenced by: '<S349>/1//Ll_q'
                                        */
  real_T u2_Value[3];                  /* Expression: [ 1/SM.Ll 1/SM.Llkq1 1/SM.Llkq2]
                                        * Referenced by: '<S348>/u2'
                                        */
  real_T Lmq_sat_InitialCondition;     /* Expression: SM.Lmsatq(1)
                                        * Referenced by: '<S344>/Lmq_sat'
                                        */
  real_T LookupTable_XData[2];         /* Expression: SM.Phisat
                                        * Referenced by: '<S344>/Lookup Table'
                                        */
  real_T LookupTable_YData[2];         /* Expression: [ 0 SM.Phisat(2:end)./SM.Lmsatq(2:end)*SM.Lmq ]
                                        * Referenced by: '<S344>/Lookup Table'
                                        */
  real_T Lmq_Gain;                     /* Expression: SM.Lmq
                                        * Referenced by: '<S344>/Lmq'
                                        */
  real_T Ll_d_Gain[3];                 /* Expression: [ 1/SM.Ll   1/SM.Llfd   1/SM.Llkd ]
                                        * Referenced by: '<S347>/1//Ll_d'
                                        */
  real_T u1_Value[3];                  /* Expression: [1/SM.Ll 1/SM.Llkd 1/SM.Llfd]
                                        * Referenced by: '<S346>/u1'
                                        */
  real_T Lmd_sat_InitialCondition;     /* Expression: SM.Lmsatd(1)
                                        * Referenced by: '<S343>/Lmd_sat'
                                        */
  real_T LookupTable_XData_k[2];       /* Expression: SM.Phisat
                                        * Referenced by: '<S343>/Lookup Table'
                                        */
  real_T LookupTable_YData_j[2];       /* Expression: [ 0 SM.Phisat(2:end)./SM.Lmsatd(2:end)*SM.Lmd ]
                                        * Referenced by: '<S343>/Lookup Table'
                                        */
  real_T Lmd_Gain;                     /* Expression: SM.Lmd
                                        * Referenced by: '<S343>/Lmd'
                                        */
  real_T u1_Value_p[25];               /* Expression: SM.R
                                        * Referenced by: '<S341>/u1'
                                        */
  real_T u1_Value_m[25];               /* Expression: zeros(SM.nState,SM.nState)
                                        * Referenced by: '<S345>/u1'
                                        */
  real_T u5_Value[25];                 /* Expression: SM.Llqd
                                        * Referenced by: '<S345>/u5'
                                        */
  real_T Constant6_Value_d[25];        /* Expression: SM.RLinv
                                        * Referenced by: '<S331>/Constant6'
                                        */
  real_T Constant2_Value_l;            /* Expression: SM.Sat
                                        * Referenced by: '<S331>/Constant2'
                                        */
  real_T Switch1_Threshold;            /* Expression: 0.5
                                        * Referenced by: '<S331>/Switch1'
                                        */
  real_T u1_Value_f[25];               /* Expression: zeros(SM.nState,SM.nState)
                                        * Referenced by: '<S339>/u1'
                                        */
  real_T Gain1_Gain_l;                 /* Expression: -1
                                        * Referenced by: '<S339>/Gain1'
                                        */
  real_T wbaseTs2_Gain;                /* Expression: SM.web*Ts/2
                                        * Referenced by: '<S350>/wbase*Ts//2'
                                        */
  real_T u5_Value_j[25];               /* Expression: eye(SM.nState,SM.nState)
                                        * Referenced by: '<S350>/u5'
                                        */
  real_T wbaseTs2_Gain_f;              /* Expression: SM.web*Ts/2
                                        * Referenced by: '<S350>/wbase*Ts//2 '
                                        */
  real_T itail_Y0;                     /* Expression: 0
                                        * Referenced by: '<S354>/itail'
                                        */
  real_T _Value;                       /* Expression: 1
                                        * Referenced by: '<S354>/1'
                                        */
  real_T DiscreteTimeIntegrator_gainva_c;/* Computed Parameter: DiscreteTimeIntegrator_gainva_c
                                          * Referenced by: '<S354>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_IC_j;  /* Expression: 0
                                        * Referenced by: '<S354>/Discrete-Time Integrator'
                                        */
  real_T Constant_Value_c;             /* Expression: 0.9
                                        * Referenced by: '<S356>/Constant'
                                        */
  real_T Saturation1_UpperSat;         /* Expression: 0.9
                                        * Referenced by: '<S356>/Saturation1'
                                        */
  real_T Saturation1_LowerSat;         /* Expression: 0
                                        * Referenced by: '<S356>/Saturation1'
                                        */
  real_T Saturation2_UpperSat;         /* Expression: 0.1
                                        * Referenced by: '<S356>/Saturation2'
                                        */
  real_T Saturation2_LowerSat;         /* Expression: 0
                                        * Referenced by: '<S356>/Saturation2'
                                        */
  real_T UnitDelay_InitialCondition_nz;/* Expression: 0
                                        * Referenced by: '<S354>/Unit Delay'
                                        */
  real_T Switch_Threshold_b;           /* Expression: 0.5
                                        * Referenced by: '<S354>/Switch'
                                        */
  real_T SwitchCurrents_Value[13];     /* Expression: zeros(13,1)
                                        * Referenced by: '<S359>/SwitchCurrents'
                                        */
  real_T UnitDelay_InitialCondition_o; /* Expression: 0
                                        * Referenced by: '<S353>/Unit Delay'
                                        */
  real_T UnitDelay6_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S67>/Unit Delay6'
                                        */
  real_T Constant3_Value_e;            /* Expression: sps.Delay
                                        * Referenced by: '<S302>/Constant3'
                                        */
  real_T Constant1_Value_m;            /* Expression: sps.Period
                                        * Referenced by: '<S302>/Constant1'
                                        */
  real_T ib1_Gain;                     /* Expression: sps.Freq
                                        * Referenced by: '<S302>/1\ib1'
                                        */
  real_T LookupTable_XData_f[3];       /* Expression: [0 .5 1]
                                        * Referenced by: '<S302>/Lookup Table'
                                        */
  real_T LookupTable_YData_k[3];       /* Expression: [0 2 0]
                                        * Referenced by: '<S302>/Lookup Table'
                                        */
  real_T Constant2_Value_a;            /* Expression: 1
                                        * Referenced by: '<S302>/Constant2'
                                        */
  real_T Gain1_Gain_d;                 /* Expression: 0.5
                                        * Referenced by: '<S280>/Gain1'
                                        */
  real_T Rotorangledthetae_gainval;    /* Computed Parameter: Rotorangledthetae_gainval
                                        * Referenced by: '<S329>/Rotor angle dthetae'
                                        */
  real_T Rotorangledthetae_IC;         /* Expression: SM.tho
                                        * Referenced by: '<S329>/Rotor angle dthetae'
                                        */
  real_T web2_Gain;                    /* Expression: SM.web
                                        * Referenced by: '<S329>/web2'
                                        */
  real_T fluxes_InitialCondition[5];   /* Expression: SM.phiqd0
                                        * Referenced by: '<S342>/fluxes'
                                        */
  real_T Constant1_Value_f;            /* Expression: SM.Sat
                                        * Referenced by: '<S331>/Constant1'
                                        */
  real_T Constant3_Value_er;           /* Expression: SM.Sat
                                        * Referenced by: '<S331>/Constant3'
                                        */
  real_T Switch_Threshold_j;           /* Expression: 0.5
                                        * Referenced by: '<S331>/Switch'
                                        */
  real_T changeIqIdcurrentsigns_Gain[5];/* Expression: SM.IqdSign
                                         * Referenced by: '<S331>/change Iq Id  current signs'
                                         */
  real_T ib_Gain;                      /* Expression: SM.ib
                                        * Referenced by: '<S333>/ib'
                                        */
  real_T valp_nom2_Value_o;            /* Expression: MagnitudeVps
                                        * Referenced by: '<S10>/valp_nom2'
                                        */
  real_T Switch5_Threshold;            /* Expression: 0.5
                                        * Referenced by: '<S10>/Switch5'
                                        */
  real_T valp_nom1_Value;              /* Expression: FreqVps
                                        * Referenced by: '<S10>/valp_nom1'
                                        */
  real_T Gain_Gain;                    /* Expression: 2*pi
                                        * Referenced by: '<S10>/Gain'
                                        */
  real_T valp_nom_Value;               /* Expression: PhaseVps
                                        * Referenced by: '<S10>/valp_nom'
                                        */
  real_T Gain3_Gain_n;                 /* Expression: pi/180
                                        * Referenced by: '<S10>/Gain3'
                                        */
  real_T P1_Value[3];                  /* Expression: [0  -2*pi/3  2*pi/3]
                                        * Referenced by: '<S10>/P1'
                                        */
  real_T donotdeletethisgain_Gain;     /* Expression: 1
                                        * Referenced by: '<S49>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_k;   /* Expression: 1
                                        * Referenced by: '<S50>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_g;   /* Expression: 1
                                        * Referenced by: '<S51>/do not delete this gain'
                                        */
  real_T Kv1_Gain;                     /* Expression: Kv
                                        * Referenced by: '<S6>/Kv1'
                                        */
  real_T donotdeletethisgain_Gain_p;   /* Expression: 1
                                        * Referenced by: '<S46>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_c;   /* Expression: 1
                                        * Referenced by: '<S47>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_e;   /* Expression: 1
                                        * Referenced by: '<S48>/do not delete this gain'
                                        */
  real_T Kv_Gain;                      /* Expression: Ki
                                        * Referenced by: '<S6>/Kv'
                                        */
  real_T UnitDelay3_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S67>/Unit Delay3'
                                        */
  real_T UnitDelay7_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S67>/Unit Delay7'
                                        */
  real_T MW_Gain;                      /* Expression: 5*2/0.9
                                        * Referenced by: '<Root>/MW'
                                        */
  real_T Gain2_Gain_i[2];              /* Expression: Gain
                                        * Referenced by: '<S76>/Gain2'
                                        */
  real_T nominalspeed_Value;           /* Expression: 1
                                        * Referenced by: '<S329>/nominal speed'
                                        */
  real_T dw_delay_InitialCondition;    /* Expression: SM.dwo
                                        * Referenced by: '<S351>/dw_delay'
                                        */
  real_T F2_Gain;                      /* Expression: 2
                                        * Referenced by: '<S351>/F2'
                                        */
  real_T dw_predict_InitialCondition;  /* Expression: SM.dwo
                                        * Referenced by: '<S351>/dw_predict'
                                        */
  real_T units_Gain;                   /* Expression: SM.Nb
                                        * Referenced by: '<S329>/units'
                                        */
  real_T Constant_Value_c5;            /* Expression: 15
                                        * Referenced by: '<Root>/Constant'
                                        */
  real_T Constant1_Value_fy;           /* Expression: 0
                                        * Referenced by: '<Root>/Constant1'
                                        */
  real_T u_Value;                      /* Expression: 0
                                        * Referenced by: '<S353>/0 4'
                                        */
  real_T Ron_Gain;                     /* Expression: 1./Ron
                                        * Referenced by: '<S353>/1//Ron'
                                        */
  real_T Switch_Threshold_l;           /* Expression: 0.5
                                        * Referenced by: '<S353>/Switch'
                                        */
  real_T Saturation_UpperSat_c;        /* Expression: inf
                                        * Referenced by: '<S353>/Saturation'
                                        */
  real_T Saturation_LowerSat_i;        /* Expression: 0
                                        * Referenced by: '<S353>/Saturation'
                                        */
  real_T pu_Gain;                      /* Expression: 1/(Pnom/(3*sqrt(2)/pi*Vnom_gen))
                                        * Referenced by: '<S8>/-> pu'
                                        */
  real_T donotdeletethisgain_Gain_ga;  /* Expression: 1
                                        * Referenced by: '<S84>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_i;   /* Expression: 1
                                        * Referenced by: '<S85>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_d;   /* Expression: 1
                                        * Referenced by: '<S86>/do not delete this gain'
                                        */
  real_T Kv_Gain_n;                    /* Expression: Ki
                                        * Referenced by: '<S64>/Kv'
                                        */
  real_T donotdeletethisgain_Gain_m;   /* Expression: 1
                                        * Referenced by: '<S87>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_eg;  /* Expression: 1
                                        * Referenced by: '<S88>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_cj;  /* Expression: 1
                                        * Referenced by: '<S89>/do not delete this gain'
                                        */
  real_T Kv1_Gain_c;                   /* Expression: Kv
                                        * Referenced by: '<S64>/Kv1'
                                        */
  real_T donotdeletethisgain_Gain_i1;  /* Expression: 1
                                        * Referenced by: '<S99>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_mu;  /* Expression: 1
                                        * Referenced by: '<S100>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_d5;  /* Expression: 1
                                        * Referenced by: '<S101>/do not delete this gain'
                                        */
  real_T Kv_Gain_h;                    /* Expression: Ki
                                        * Referenced by: '<S65>/Kv'
                                        */
  real_T donotdeletethisgain_Gain_m4;  /* Expression: 1
                                        * Referenced by: '<S102>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_h;   /* Expression: 1
                                        * Referenced by: '<S103>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_kc;  /* Expression: 1
                                        * Referenced by: '<S104>/do not delete this gain'
                                        */
  real_T Kv1_Gain_cr;                  /* Expression: Kv
                                        * Referenced by: '<S65>/Kv1'
                                        */
  real_T DiscreteTimeIntegrator_gainva_a;/* Computed Parameter: DiscreteTimeIntegrator_gainva_a
                                          * Referenced by: '<S126>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_UpperSat;/* Expression: UpperLimit
                                          * Referenced by: '<S126>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_LowerSat;/* Expression: LowerLimit
                                          * Referenced by: '<S126>/Discrete-Time Integrator'
                                          */
  real_T Duk_Gain;                     /* Expression: sps.D
                                        * Referenced by: '<S221>/D*u(k)'
                                        */
  real_T Delay_x1_InitialCondition;    /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S221>/Delay_x1'
                                        */
  real_T C11_Gain;                     /* Expression: sps.C11
                                        * Referenced by: '<S224>/C11'
                                        */
  real_T Delay_x2_InitialCondition;    /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S221>/Delay_x2'
                                        */
  real_T C12_Gain;                     /* Expression: sps.C12
                                        * Referenced by: '<S224>/C12'
                                        */
  real_T Switch_Threshold_p;           /* Expression: 1
                                        * Referenced by: '<S120>/Switch'
                                        */
  real_T Duk_Gain_l;                   /* Expression: sps.D
                                        * Referenced by: '<S197>/D*u(k)'
                                        */
  real_T Delay_x1_InitialCondition_d;  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S197>/Delay_x1'
                                        */
  real_T C11_Gain_m;                   /* Expression: sps.C11
                                        * Referenced by: '<S200>/C11'
                                        */
  real_T Delay_x2_InitialCondition_d;  /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S197>/Delay_x2'
                                        */
  real_T C12_Gain_b;                   /* Expression: sps.C12
                                        * Referenced by: '<S200>/C12'
                                        */
  real_T theta_gainval;                /* Computed Parameter: theta_gainval
                                        * Referenced by: '<S329>/theta'
                                        */
  real_T theta_IC;                     /* Expression: SM.tho
                                        * Referenced by: '<S329>/theta'
                                        */
  real_T t_Gain;                       /* Expression: 180/pi
                                        * Referenced by: '<S329>/t'
                                        */
  real_T degrd_Gain;                   /* Expression: pi/180
                                        * Referenced by: '<S8>/deg->rd'
                                        */
  real_T Duk_Gain_g;                   /* Expression: sps.D
                                        * Referenced by: '<S213>/D*u(k)'
                                        */
  real_T Delay_x1_InitialCondition_h;  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S213>/Delay_x1'
                                        */
  real_T C11_Gain_k;                   /* Expression: sps.C11
                                        * Referenced by: '<S216>/C11'
                                        */
  real_T Delay_x2_InitialCondition_n;  /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S213>/Delay_x2'
                                        */
  real_T C12_Gain_h;                   /* Expression: sps.C12
                                        * Referenced by: '<S216>/C12'
                                        */
  real_T pairsofpoles_Gain;            /* Expression: p
                                        * Referenced by: '<S122>/# pairs of poles'
                                        */
  real_T Constant4_Value_m;            /* Expression: 2*pi
                                        * Referenced by: '<S122>/Constant4'
                                        */
  real_T Duk_Gain_m;                   /* Expression: sps.D
                                        * Referenced by: '<S201>/D*u(k)'
                                        */
  real_T Delay_x1_InitialCondition_e;  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S201>/Delay_x1'
                                        */
  real_T C11_Gain_o;                   /* Expression: sps.C11
                                        * Referenced by: '<S204>/C11'
                                        */
  real_T Delay_x2_InitialCondition_m;  /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S201>/Delay_x2'
                                        */
  real_T C12_Gain_j;                   /* Expression: sps.C12
                                        * Referenced by: '<S204>/C12'
                                        */
  real_T Gain3_Gain_m[9];              /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
                                        * Referenced by: '<S275>/Gain3'
                                        */
  real_T Gain1_Gain_h;                 /* Expression: 2/3
                                        * Referenced by: '<S275>/Gain1'
                                        */
  real_T avoiddivisionby0_UpperSat_d;  /* Expression: inf
                                        * Referenced by: '<S145>/avoid division by 0'
                                        */
  real_T avoiddivisionby0_LowerSat_n;  /* Expression: 1e-6
                                        * Referenced by: '<S145>/avoid division by 0'
                                        */
  real_T Gain_Gain_b;                  /* Expression: -1
                                        * Referenced by: '<S145>/Gain'
                                        */
  real_T D_Gain;                       /* Expression: sps.D
                                        * Referenced by: '<S261>/D'
                                        */
  real_T Delay_x_InitialCondition;     /* Expression: sps.x0
                                        * Referenced by: '<S261>/Delay_x'
                                        */
  real_T C_Gain;                       /* Expression: sps.C
                                        * Referenced by: '<S261>/C'
                                        */
  real_T Saturation2_UpperSat_o;       /* Expression: UpperLimit
                                        * Referenced by: '<S126>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_f;       /* Expression: LowerLimit
                                        * Referenced by: '<S126>/Saturation2'
                                        */
  real_T Vdc_refV_Value;               /* Expression: Vdc_nom
                                        * Referenced by: '<S121>/Vdc_ref (V)'
                                        */
  real_T Duk_Gain_a;                   /* Expression: sps.D
                                        * Referenced by: '<S217>/D*u(k)'
                                        */
  real_T Delay_x1_InitialCondition_f;  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S217>/Delay_x1'
                                        */
  real_T C11_Gain_e;                   /* Expression: sps.C11
                                        * Referenced by: '<S220>/C11'
                                        */
  real_T Delay_x2_InitialCondition_e;  /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S217>/Delay_x2'
                                        */
  real_T C12_Gain_a;                   /* Expression: sps.C12
                                        * Referenced by: '<S220>/C12'
                                        */
  real_T pu_Gain_b;                    /* Expression: 1/Vdc_nom
                                        * Referenced by: '<S121>/->pu'
                                        */
  real_T Duk_Gain_j;                   /* Expression: sps.D
                                        * Referenced by: '<S189>/D*u(k)'
                                        */
  real_T Delay_x1_InitialCondition_l;  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S189>/Delay_x1'
                                        */
  real_T C11_Gain_a;                   /* Expression: sps.C11
                                        * Referenced by: '<S192>/C11'
                                        */
  real_T Delay_x2_InitialCondition_di; /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S189>/Delay_x2'
                                        */
  real_T C12_Gain_p;                   /* Expression: sps.C12
                                        * Referenced by: '<S192>/C12'
                                        */
  real_T Gain3_Gain_a[9];              /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
                                        * Referenced by: '<S263>/Gain3'
                                        */
  real_T Gain1_Gain_k0;                /* Expression: 2/3
                                        * Referenced by: '<S263>/Gain1'
                                        */
  real_T DiscreteTimeIntegrator_gainva_l;/* Computed Parameter: DiscreteTimeIntegrator_gainva_l
                                          * Referenced by: '<S225>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_IC_n;  /* Expression: sps.Phase_Init*pi/180
                                        * Referenced by: '<S225>/Discrete-Time Integrator'
                                        */
  real_T Constant4_Value_mc;           /* Expression: 2*pi
                                        * Referenced by: '<S225>/Constant4'
                                        */
  real_T Constant1_Value_f5;           /* Expression: L_RL
                                        * Referenced by: '<S121>/Constant1'
                                        */
  real_T UnitDelay_InitialCondition_d; /* Expression: sps.Finit
                                        * Referenced by: '<S225>/Unit Delay'
                                        */
  real_T Fnom_Value;                   /* Expression: Fnom
                                        * Referenced by: '<S122>/Fnom'
                                        */
  real_T Constant2_Value_i;            /* Expression: 0
                                        * Referenced by: '<S121>/Constant2'
                                        */
  real_T ICi_ic_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S131>/IC = i_ic'
                                        */
  real_T DiscreteTimeIntegrator_gainva_p;/* Computed Parameter: DiscreteTimeIntegrator_gainva_p
                                          * Referenced by: '<S128>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_UpperS_b;/* Expression: UpperLimit
                                          * Referenced by: '<S128>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_LowerS_p;/* Expression: LowerLimit
                                          * Referenced by: '<S128>/Discrete-Time Integrator'
                                          */
  real_T Saturation2_UpperSat_g;       /* Expression: UpperLimit
                                        * Referenced by: '<S128>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_m;       /* Expression: LowerLimit
                                        * Referenced by: '<S128>/Saturation2'
                                        */
  real_T Imax2_Value;                  /* Expression: Imax_grid_conv^2
                                        * Referenced by: '<S121>/Imax^2'
                                        */
  real_T D_Gain_b;                     /* Expression: sps.D
                                        * Referenced by: '<S134>/D'
                                        */
  real_T Delay_x_InitialCondition_c;   /* Expression: sps.x0
                                        * Referenced by: '<S134>/Delay_x'
                                        */
  real_T C_Gain_d;                     /* Expression: sps.C
                                        * Referenced by: '<S134>/C'
                                        */
  real_T Gain_Gain_j;                  /* Expression: -1
                                        * Referenced by: '<S121>/Gain'
                                        */
  real_T Constant3_Value_i;            /* Expression: Ki_volt_reg
                                        * Referenced by: '<S121>/Constant3'
                                        */
  real_T UnitDelay1_InitialCondition_a;/* Expression: 1
                                        * Referenced by: '<S121>/Unit Delay1'
                                        */
  real_T DiscreteTimeIntegrator1_gainv_d;/* Computed Parameter: DiscreteTimeIntegrator1_gainv_d
                                          * Referenced by: '<S121>/Discrete-Time Integrator1'
                                          */
  real_T DiscreteTimeIntegrator1_UpperSa;/* Expression: Var_reg_output_hi_limit
                                          * Referenced by: '<S121>/Discrete-Time Integrator1'
                                          */
  real_T DiscreteTimeIntegrator1_LowerSa;/* Expression: Var_reg_output_low_limit
                                          * Referenced by: '<S121>/Discrete-Time Integrator1'
                                          */
  real_T puV_Gain;                     /* Expression: Vnom/sqrt(3)*sqrt(2)
                                        * Referenced by: '<S122>/pu->V'
                                        */
  real_T Gain3_Gain_nj[9];             /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
                                        * Referenced by: '<S161>/Gain3'
                                        */
  real_T Gain1_Gain_dl;                /* Expression: 2/3
                                        * Referenced by: '<S161>/Gain1'
                                        */
  real_T Integ4_gainval_g;             /* Computed Parameter: Integ4_gainval_g
                                        * Referenced by: '<S154>/Integ4'
                                        */
  real_T Integ4_IC_c;                  /* Expression: 0
                                        * Referenced by: '<S154>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperS_ft;/* Expression: 1e6
                                          * Referenced by: '<S154>/To avoid division by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_b;/* Expression: 1e-6
                                          * Referenced by: '<S154>/To avoid division by zero'
                                          */
  real_T UnitDelay_InitialCondition_g; /* Expression: 0
                                        * Referenced by: '<S158>/Unit Delay'
                                        */
  real_T Step_Y0_c;                    /* Expression: 0
                                        * Referenced by: '<S154>/Step'
                                        */
  real_T Step_YFinal_g;                /* Expression: 1
                                        * Referenced by: '<S154>/Step'
                                        */
  real_T Switch_Threshold_n;           /* Expression: 0.5
                                        * Referenced by: '<S154>/Switch'
                                        */
  real_T Integ4_gainval_e;             /* Computed Parameter: Integ4_gainval_e
                                        * Referenced by: '<S153>/Integ4'
                                        */
  real_T Integ4_IC_m;                  /* Expression: 0
                                        * Referenced by: '<S153>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_p;/* Expression: 1e6
                                          * Referenced by: '<S153>/To avoid division by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_d;/* Expression: 1e-6
                                          * Referenced by: '<S153>/To avoid division by zero'
                                          */
  real_T UnitDelay_InitialCondition_k; /* Expression: 0
                                        * Referenced by: '<S156>/Unit Delay'
                                        */
  real_T Step_Y0_a;                    /* Expression: 0
                                        * Referenced by: '<S153>/Step'
                                        */
  real_T Step_YFinal_f;                /* Expression: 1
                                        * Referenced by: '<S153>/Step'
                                        */
  real_T Switch_Threshold_m;           /* Expression: 0.5
                                        * Referenced by: '<S153>/Switch'
                                        */
  real_T Vpu_Gain;                     /* Expression: 1/(Vnom/sqrt(3)*sqrt(2))
                                        * Referenced by: '<S122>/V->pu'
                                        */
  real_T IC0_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<S131>/IC = 0'
                                        */
  real_T Constant5_Value_c;            /* Expression: R_RL
                                        * Referenced by: '<S121>/Constant5'
                                        */
  real_T Duk_Gain_f;                   /* Expression: sps.D
                                        * Referenced by: '<S193>/D*u(k)'
                                        */
  real_T Delay_x1_InitialCondition_i;  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S193>/Delay_x1'
                                        */
  real_T C11_Gain_h;                   /* Expression: sps.C11
                                        * Referenced by: '<S196>/C11'
                                        */
  real_T Delay_x2_InitialCondition_nr; /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S193>/Delay_x2'
                                        */
  real_T C12_Gain_h3;                  /* Expression: sps.C12
                                        * Referenced by: '<S196>/C12'
                                        */
  real_T Gain3_Gain_a3[9];             /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
                                        * Referenced by: '<S269>/Gain3'
                                        */
  real_T Gain1_Gain_db;                /* Expression: 2/3
                                        * Referenced by: '<S269>/Gain1'
                                        */
  real_T DiscreteTimeIntegrator_gainv_lq;/* Computed Parameter: DiscreteTimeIntegrator_gainv_lq
                                          * Referenced by: '<S129>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_UpperS_d;/* Expression: UpperLimit
                                          * Referenced by: '<S129>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_LowerS_k;/* Expression: LowerLimit
                                          * Referenced by: '<S129>/Discrete-Time Integrator'
                                          */
  real_T Saturation2_UpperSat_gz;      /* Expression: UpperLimit
                                        * Referenced by: '<S129>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_a;       /* Expression: LowerLimit
                                        * Referenced by: '<S129>/Saturation2'
                                        */
  real_T Constant6_Value_g;            /* Expression: R_RL
                                        * Referenced by: '<S121>/Constant6'
                                        */
  real_T Constant4_Value_de;           /* Expression: L_RL
                                        * Referenced by: '<S121>/Constant4'
                                        */
  real_T K_Value;                      /* Expression: Vnom*2*sqrt(2/3)
                                        * Referenced by: '<S121>/K'
                                        */
  real_T Avoiddivisionbyzero_UpperSat; /* Expression: 1e6
                                        * Referenced by: '<S121>/Avoid division by zero'
                                        */
  real_T Avoiddivisionbyzero_LowerSat; /* Expression: 1e-6
                                        * Referenced by: '<S121>/Avoid division by zero'
                                        */
  real_T Mod_index_max_UpperSat;       /* Expression: Mod_index_max
                                        * Referenced by: '<S121>/0-Mod_index_max'
                                        */
  real_T Mod_index_max_LowerSat;       /* Expression: 0
                                        * Referenced by: '<S121>/0-Mod_index_max'
                                        */
  real_T A_Gain;                       /* Expression: sps.A
                                        * Referenced by: '<S134>/A'
                                        */
  real_T B_Gain;                       /* Expression: sps.B
                                        * Referenced by: '<S134>/B'
                                        */
  real_T Duk_Gain_c;                   /* Expression: sps.D
                                        * Referenced by: '<S209>/D*u(k)'
                                        */
  real_T Delay_x1_InitialCondition_p;  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S209>/Delay_x1'
                                        */
  real_T C11_Gain_j;                   /* Expression: sps.C11
                                        * Referenced by: '<S212>/C11'
                                        */
  real_T Delay_x2_InitialCondition_g;  /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S209>/Delay_x2'
                                        */
  real_T C12_Gain_ax;                  /* Expression: sps.C12
                                        * Referenced by: '<S212>/C12'
                                        */
  real_T Constant1_Value_g;            /* Expression: Vnom/sqrt(3)*sqrt(2)
                                        * Referenced by: '<S122>/Constant1'
                                        */
  real_T C_var_filter_Gain;            /* Expression: C_var_filter
                                        * Referenced by: '<S122>/C_var_filter'
                                        */
  real_T puA_Gain;                     /* Expression: Pnom/sqrt(3)/Vnom*sqrt(2)
                                        * Referenced by: '<S122>/pu->A'
                                        */
  real_T Gain3_Gain_o[9];              /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
                                        * Referenced by: '<S175>/Gain3'
                                        */
  real_T Gain1_Gain_ep;                /* Expression: 2/3
                                        * Referenced by: '<S175>/Gain1'
                                        */
  real_T Integ4_gainval_l;             /* Computed Parameter: Integ4_gainval_l
                                        * Referenced by: '<S168>/Integ4'
                                        */
  real_T Integ4_IC_n;                  /* Expression: 0
                                        * Referenced by: '<S168>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_e;/* Expression: 1e6
                                          * Referenced by: '<S168>/To avoid division by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_a;/* Expression: 1e-6
                                          * Referenced by: '<S168>/To avoid division by zero'
                                          */
  real_T UnitDelay_InitialCondition_f; /* Expression: 0
                                        * Referenced by: '<S172>/Unit Delay'
                                        */
  real_T Step_Y0_l;                    /* Expression: 0
                                        * Referenced by: '<S168>/Step'
                                        */
  real_T Step_YFinal_j;                /* Expression: 1
                                        * Referenced by: '<S168>/Step'
                                        */
  real_T Switch_Threshold_nb;          /* Expression: 0.5
                                        * Referenced by: '<S168>/Switch'
                                        */
  real_T Integ4_gainval_b;             /* Computed Parameter: Integ4_gainval_b
                                        * Referenced by: '<S167>/Integ4'
                                        */
  real_T Integ4_IC_f;                  /* Expression: 0
                                        * Referenced by: '<S167>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperSa_a;/* Expression: 1e6
                                          * Referenced by: '<S167>/To avoid division by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerS_dg;/* Expression: 1e-6
                                          * Referenced by: '<S167>/To avoid division by zero'
                                          */
  real_T UnitDelay_InitialCondition_h; /* Expression: 0
                                        * Referenced by: '<S170>/Unit Delay'
                                        */
  real_T Step_Y0_o;                    /* Expression: 0
                                        * Referenced by: '<S167>/Step'
                                        */
  real_T Step_YFinal_gq;               /* Expression: 1
                                        * Referenced by: '<S167>/Step'
                                        */
  real_T Switch_Threshold_jy;          /* Expression: 0.5
                                        * Referenced by: '<S167>/Switch'
                                        */
  real_T RadDeg_Gain;                  /* Expression: 180/pi
                                        * Referenced by: '<S150>/Rad->Deg.'
                                        */
  real_T RadDeg_Gain_j;                /* Expression: 180/pi
                                        * Referenced by: '<S151>/Rad->Deg.'
                                        */
  real_T DegRad_Gain;                  /* Expression: pi/180
                                        * Referenced by: '<S142>/Deg->Rad'
                                        */
  real_T Gain1_Gain_ki;                /* Expression: 3/2
                                        * Referenced by: '<S142>/Gain1'
                                        */
  real_T varpu_Gain;                   /* Expression: -1/Pnom
                                        * Referenced by: '<S122>/var->pu '
                                        */
  real_T V0_Value;                     /* Expression: 0
                                        * Referenced by: '<S121>/V0'
                                        */
  real_T Gain3_Gain_g[9];              /* Expression: [ 1   0   1; -1/2  sqrt(3)/2   1; -1/2  -sqrt(3)/2  1 ]
                                        * Referenced by: '<S136>/Gain3'
                                        */
  real_T C_var_filterQ_Gain;           /* Expression: C_var_filter/Q_filter
                                        * Referenced by: '<S122>/C_var_filter//Q'
                                        */
  real_T A11_Gain;                     /* Expression: sps.A11
                                        * Referenced by: '<S190>/A11'
                                        */
  real_T A12_Gain;                     /* Expression: sps.A12
                                        * Referenced by: '<S190>/A12'
                                        */
  real_T A21_Gain;                     /* Expression: sps.A21
                                        * Referenced by: '<S190>/A21'
                                        */
  real_T A22_Gain;                     /* Expression: sps.A22
                                        * Referenced by: '<S190>/A22'
                                        */
  real_T B11_Gain;                     /* Expression: sps.B11
                                        * Referenced by: '<S191>/B11'
                                        */
  real_T B21_Gain;                     /* Expression: sps.B21
                                        * Referenced by: '<S191>/B21'
                                        */
  real_T A11_Gain_k;                   /* Expression: sps.A11
                                        * Referenced by: '<S194>/A11'
                                        */
  real_T A12_Gain_g;                   /* Expression: sps.A12
                                        * Referenced by: '<S194>/A12'
                                        */
  real_T A21_Gain_l;                   /* Expression: sps.A21
                                        * Referenced by: '<S194>/A21'
                                        */
  real_T A22_Gain_k;                   /* Expression: sps.A22
                                        * Referenced by: '<S194>/A22'
                                        */
  real_T B11_Gain_k;                   /* Expression: sps.B11
                                        * Referenced by: '<S195>/B11'
                                        */
  real_T B21_Gain_m;                   /* Expression: sps.B21
                                        * Referenced by: '<S195>/B21'
                                        */
  real_T A11_Gain_j;                   /* Expression: sps.A11
                                        * Referenced by: '<S198>/A11'
                                        */
  real_T A12_Gain_o;                   /* Expression: sps.A12
                                        * Referenced by: '<S198>/A12'
                                        */
  real_T A21_Gain_j;                   /* Expression: sps.A21
                                        * Referenced by: '<S198>/A21'
                                        */
  real_T A22_Gain_i;                   /* Expression: sps.A22
                                        * Referenced by: '<S198>/A22'
                                        */
  real_T B11_Gain_n;                   /* Expression: sps.B11
                                        * Referenced by: '<S199>/B11'
                                        */
  real_T B21_Gain_n;                   /* Expression: sps.B21
                                        * Referenced by: '<S199>/B21'
                                        */
  real_T A11_Gain_d;                   /* Expression: sps.A11
                                        * Referenced by: '<S202>/A11'
                                        */
  real_T A12_Gain_b;                   /* Expression: sps.A12
                                        * Referenced by: '<S202>/A12'
                                        */
  real_T A21_Gain_a;                   /* Expression: sps.A21
                                        * Referenced by: '<S202>/A21'
                                        */
  real_T A22_Gain_k5;                  /* Expression: sps.A22
                                        * Referenced by: '<S202>/A22'
                                        */
  real_T B11_Gain_m;                   /* Expression: sps.B11
                                        * Referenced by: '<S203>/B11'
                                        */
  real_T B21_Gain_o;                   /* Expression: sps.B21
                                        * Referenced by: '<S203>/B21'
                                        */
  real_T Delay_x1_InitialCondition_j;  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S205>/Delay_x1'
                                        */
  real_T A11_Gain_b;                   /* Expression: sps.A11
                                        * Referenced by: '<S206>/A11'
                                        */
  real_T Delay_x2_InitialCondition_b;  /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S205>/Delay_x2'
                                        */
  real_T A12_Gain_gm;                  /* Expression: sps.A12
                                        * Referenced by: '<S206>/A12'
                                        */
  real_T A21_Gain_lt;                  /* Expression: sps.A21
                                        * Referenced by: '<S206>/A21'
                                        */
  real_T A22_Gain_o;                   /* Expression: sps.A22
                                        * Referenced by: '<S206>/A22'
                                        */
  real_T B11_Gain_me;                  /* Expression: sps.B11
                                        * Referenced by: '<S207>/B11'
                                        */
  real_T B21_Gain_f;                   /* Expression: sps.B21
                                        * Referenced by: '<S207>/B21'
                                        */
  real_T Duk_Gain_k;                   /* Expression: sps.D
                                        * Referenced by: '<S205>/D*u(k)'
                                        */
  real_T C11_Gain_ai;                  /* Expression: sps.C11
                                        * Referenced by: '<S208>/C11'
                                        */
  real_T C12_Gain_e;                   /* Expression: sps.C12
                                        * Referenced by: '<S208>/C12'
                                        */
  real_T A11_Gain_jo;                  /* Expression: sps.A11
                                        * Referenced by: '<S218>/A11'
                                        */
  real_T A12_Gain_c;                   /* Expression: sps.A12
                                        * Referenced by: '<S218>/A12'
                                        */
  real_T A21_Gain_e;                   /* Expression: sps.A21
                                        * Referenced by: '<S218>/A21'
                                        */
  real_T A22_Gain_kz;                  /* Expression: sps.A22
                                        * Referenced by: '<S218>/A22'
                                        */
  real_T B11_Gain_h;                   /* Expression: sps.B11
                                        * Referenced by: '<S219>/B11'
                                        */
  real_T B21_Gain_np;                  /* Expression: sps.B21
                                        * Referenced by: '<S219>/B21'
                                        */
  real_T A11_Gain_h;                   /* Expression: sps.A11
                                        * Referenced by: '<S210>/A11'
                                        */
  real_T A12_Gain_j;                   /* Expression: sps.A12
                                        * Referenced by: '<S210>/A12'
                                        */
  real_T A21_Gain_h;                   /* Expression: sps.A21
                                        * Referenced by: '<S210>/A21'
                                        */
  real_T A22_Gain_l;                   /* Expression: sps.A22
                                        * Referenced by: '<S210>/A22'
                                        */
  real_T B11_Gain_f;                   /* Expression: sps.B11
                                        * Referenced by: '<S211>/B11'
                                        */
  real_T B21_Gain_nt;                  /* Expression: sps.B21
                                        * Referenced by: '<S211>/B21'
                                        */
  real_T A11_Gain_jk;                  /* Expression: sps.A11
                                        * Referenced by: '<S222>/A11'
                                        */
  real_T A12_Gain_e;                   /* Expression: sps.A12
                                        * Referenced by: '<S222>/A12'
                                        */
  real_T A21_Gain_i;                   /* Expression: sps.A21
                                        * Referenced by: '<S222>/A21'
                                        */
  real_T A22_Gain_ir;                  /* Expression: sps.A22
                                        * Referenced by: '<S222>/A22'
                                        */
  real_T B11_Gain_g;                   /* Expression: sps.B11
                                        * Referenced by: '<S223>/B11'
                                        */
  real_T B21_Gain_fv;                  /* Expression: sps.B21
                                        * Referenced by: '<S223>/B21'
                                        */
  real_T A11_Gain_f;                   /* Expression: sps.A11
                                        * Referenced by: '<S214>/A11'
                                        */
  real_T A12_Gain_i;                   /* Expression: sps.A12
                                        * Referenced by: '<S214>/A12'
                                        */
  real_T A21_Gain_p;                   /* Expression: sps.A21
                                        * Referenced by: '<S214>/A21'
                                        */
  real_T A22_Gain_c;                   /* Expression: sps.A22
                                        * Referenced by: '<S214>/A22'
                                        */
  real_T B11_Gain_kp;                  /* Expression: sps.B11
                                        * Referenced by: '<S215>/B11'
                                        */
  real_T B21_Gain_l;                   /* Expression: sps.B21
                                        * Referenced by: '<S215>/B21'
                                        */
  real_T Constant1_Value_o;            /* Expression: sps.AGC
                                        * Referenced by: '<S225>/Constant1'
                                        */
  real_T Gain3_Gain_c[9];              /* Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
                                        * Referenced by: '<S255>/Gain3'
                                        */
  real_T Gain1_Gain_b;                 /* Expression: 2/3
                                        * Referenced by: '<S255>/Gain1'
                                        */
  real_T Integ4_gainval_d;             /* Computed Parameter: Integ4_gainval_d
                                        * Referenced by: '<S251>/Integ4'
                                        */
  real_T Integ4_IC_ka;                 /* Expression: 0
                                        * Referenced by: '<S251>/Integ4'
                                        */
  real_T Toavoiddivisionbyzero_UpperS_en;/* Expression: 1e6
                                          * Referenced by: '<S251>/To avoid division  by zero'
                                          */
  real_T Toavoiddivisionbyzero_LowerSa_n;/* Expression: eps
                                          * Referenced by: '<S251>/To avoid division  by zero'
                                          */
  real_T SFunction_P1_Size_c[2];       /* Computed Parameter: SFunction_P1_Size_c
                                        * Referenced by: '<S253>/S-Function'
                                        */
  real_T SFunction_P1_i;               /* Expression: MaxDelay
                                        * Referenced by: '<S253>/S-Function'
                                        */
  real_T SFunction_P2_Size_e[2];       /* Computed Parameter: SFunction_P2_Size_e
                                        * Referenced by: '<S253>/S-Function'
                                        */
  real_T SFunction_P2_f;               /* Expression: Ts
                                        * Referenced by: '<S253>/S-Function'
                                        */
  real_T SFunction_P3_Size_m[2];       /* Computed Parameter: SFunction_P3_Size_m
                                        * Referenced by: '<S253>/S-Function'
                                        */
  real_T SFunction_P3_g;               /* Expression: InitialValue
                                        * Referenced by: '<S253>/S-Function'
                                        */
  real_T SFunction_P4_Size_n[2];       /* Computed Parameter: SFunction_P4_Size_n
                                        * Referenced by: '<S253>/S-Function'
                                        */
  real_T SFunction_P4_b;               /* Expression: DFT
                                        * Referenced by: '<S253>/S-Function'
                                        */
  real_T UnitDelay_InitialCondition_i; /* Expression: 0
                                        * Referenced by: '<S252>/Unit Delay'
                                        */
  real_T Constant_Value_e;             /* Expression: 1/sps.Finit
                                        * Referenced by: '<S251>/Constant'
                                        */
  real_T UnitDelay1_InitialCondition_c;/* Expression: sps.Vinit
                                        * Referenced by: '<S251>/Unit Delay1'
                                        */
  real_T DiscreteDerivative_DenCoef[2];/* Expression: [ TcD  Ts-TcD ]
                                        * Referenced by: '<S227>/Discrete Derivative '
                                        */
  real_T DiscreteDerivative_InitialState;/* Expression: 0
                                          * Referenced by: '<S227>/Discrete Derivative '
                                          */
  real_T DiscreteTimeIntegrator_gainva_o;/* Computed Parameter: DiscreteTimeIntegrator_gainva_o
                                          * Referenced by: '<S227>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_UpperS_e;/* Expression: Par_Limits(1)
                                          * Referenced by: '<S227>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_LowerS_b;/* Expression: Par_Limits(2)
                                          * Referenced by: '<S227>/Discrete-Time Integrator'
                                          */
  real_T Saturation1_UpperSat_h;       /* Expression: Par_Limits(1)
                                        * Referenced by: '<S227>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_a;       /* Expression: Par_Limits(2)
                                        * Referenced by: '<S227>/Saturation1'
                                        */
  real_T Gain10_Gain;                  /* Expression: 1/2/pi
                                        * Referenced by: '<S225>/Gain10'
                                        */
  real_T RateLimiter_RisingLim;        /* Computed Parameter: RateLimiter_RisingLim
                                        * Referenced by: '<S225>/Rate Limiter'
                                        */
  real_T RateLimiter_FallingLim;       /* Computed Parameter: RateLimiter_FallingLim
                                        * Referenced by: '<S225>/Rate Limiter'
                                        */
  real_T RateLimiter_IC;               /* Expression: sps.Finit
                                        * Referenced by: '<S225>/Rate Limiter'
                                        */
  real_T Delay_x1_InitialCondition_m;  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S247>/Delay_x1'
                                        */
  real_T A11_Gain_n;                   /* Expression: sps.A11
                                        * Referenced by: '<S248>/A11'
                                        */
  real_T Delay_x2_InitialCondition_o;  /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S247>/Delay_x2'
                                        */
  real_T A12_Gain_a;                   /* Expression: sps.A12
                                        * Referenced by: '<S248>/A12'
                                        */
  real_T A21_Gain_c;                   /* Expression: sps.A21
                                        * Referenced by: '<S248>/A21'
                                        */
  real_T A22_Gain_g;                   /* Expression: sps.A22
                                        * Referenced by: '<S248>/A22'
                                        */
  real_T B11_Gain_e;                   /* Expression: sps.B11
                                        * Referenced by: '<S249>/B11'
                                        */
  real_T B21_Gain_a;                   /* Expression: sps.B21
                                        * Referenced by: '<S249>/B21'
                                        */
  real_T Duk_Gain_p;                   /* Expression: sps.D
                                        * Referenced by: '<S247>/D*u(k)'
                                        */
  real_T C11_Gain_g;                   /* Expression: sps.C11
                                        * Referenced by: '<S250>/C11'
                                        */
  real_T C12_Gain_jf;                  /* Expression: sps.C12
                                        * Referenced by: '<S250>/C12'
                                        */
  real_T A_Gain_h;                     /* Expression: sps.A
                                        * Referenced by: '<S261>/A'
                                        */
  real_T B_Gain_k;                     /* Expression: sps.B
                                        * Referenced by: '<S261>/B'
                                        */
  real_T Wpu_Gain;                     /* Expression: -1/Pnom
                                        * Referenced by: '<S122>/W->pu'
                                        */
  real_T UnitDelay2_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S67>/Unit Delay2'
                                        */
  real_T inf_UpperSat;                 /* Expression: inf
                                        * Referenced by: '<S125>/0-inf'
                                        */
  real_T inf_LowerSat;                 /* Expression: 0
                                        * Referenced by: '<S125>/0-inf'
                                        */
  real_T pu_elecpu_mec_Gain;           /* Expression: Pnom/Pmec
                                        * Referenced by: '<S125>/pu_elec->pu_mec'
                                        */
  real_T Switch_Threshold_ny;          /* Expression: 0.75
                                        * Referenced by: '<S125>/Switch'
                                        */
  real_T D_Gain_i;                     /* Expression: sps.D
                                        * Referenced by: '<S310>/D'
                                        */
  real_T Delay_x_InitialCondition_b;   /* Expression: sps.x0
                                        * Referenced by: '<S310>/Delay_x'
                                        */
  real_T C_Gain_g;                     /* Expression: sps.C
                                        * Referenced by: '<S310>/C'
                                        */
  real_T Gain2_Gain_n;                 /* Expression: Kp_pitch
                                        * Referenced by: '<S125>/Gain2'
                                        */
  real_T pitch_max_UpperSat;           /* Expression: pitch_max
                                        * Referenced by: '<S125>/0-pitch_max'
                                        */
  real_T pitch_max_LowerSat;           /* Expression: 0
                                        * Referenced by: '<S125>/0-pitch_max'
                                        */
  real_T Constant2_Value_e;            /* Expression: 1
                                        * Referenced by: '<S125>/Constant2'
                                        */
  real_T DiscreteTimeIntegrator_gainva_f;/* Computed Parameter: DiscreteTimeIntegrator_gainva_f
                                          * Referenced by: '<S305>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_UpperS_f;/* Expression: UpperLimit
                                          * Referenced by: '<S305>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_Lower_pa;/* Expression: LowerLimit
                                          * Referenced by: '<S305>/Discrete-Time Integrator'
                                          */
  real_T Saturation2_UpperSat_f;       /* Expression: UpperLimit
                                        * Referenced by: '<S305>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_m2;      /* Expression: LowerLimit
                                        * Referenced by: '<S305>/Saturation2'
                                        */
  real_T RateLimiter1_RisingLim;       /* Expression: pitch_rate
                                        * Referenced by: '<S125>/Rate Limiter   1'
                                        */
  real_T RateLimiter1_FallingLim;      /* Expression: -pitch_rate
                                        * Referenced by: '<S125>/Rate Limiter   1'
                                        */
  real_T pitch_max_UpperSat_a;         /* Expression: pitch_max
                                        * Referenced by: '<S125>/0-pitch_max '
                                        */
  real_T pitch_max_LowerSat_n;         /* Expression: 0
                                        * Referenced by: '<S125>/0-pitch_max '
                                        */
  real_T DiscreteTimeIntegrator_gainva_e;/* Computed Parameter: DiscreteTimeIntegrator_gainva_e
                                          * Referenced by: '<S304>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_UpperS_h;/* Expression: UpperLimit
                                          * Referenced by: '<S304>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_LowerS_l;/* Expression: LowerLimit
                                          * Referenced by: '<S304>/Discrete-Time Integrator'
                                          */
  real_T Saturation2_UpperSat_b;       /* Expression: UpperLimit
                                        * Referenced by: '<S304>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_d;       /* Expression: LowerLimit
                                        * Referenced by: '<S304>/Saturation2'
                                        */
  real_T DiscreteTimeIntegrator_gainv_cb;/* Computed Parameter: DiscreteTimeIntegrator_gainv_cb
                                          * Referenced by: '<S306>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_UpperS_i;/* Expression: UpperLimit
                                          * Referenced by: '<S306>/Discrete-Time Integrator'
                                          */
  real_T DiscreteTimeIntegrator_Lower_lu;/* Expression: LowerLimit
                                          * Referenced by: '<S306>/Discrete-Time Integrator'
                                          */
  real_T Saturation2_UpperSat_j;       /* Expression: UpperLimit
                                        * Referenced by: '<S306>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_aj;      /* Expression: LowerLimit
                                        * Referenced by: '<S306>/Saturation2'
                                        */
  real_T Delay_x_InitialCondition_d;   /* Expression: sps.x0
                                        * Referenced by: '<S309>/Delay_x'
                                        */
  real_T A_Gain_m;                     /* Expression: sps.A
                                        * Referenced by: '<S309>/A'
                                        */
  real_T B_Gain_b;                     /* Expression: sps.B
                                        * Referenced by: '<S309>/B'
                                        */
  real_T C_Gain_a;                     /* Expression: sps.C
                                        * Referenced by: '<S309>/C'
                                        */
  real_T D_Gain_p;                     /* Expression: sps.D
                                        * Referenced by: '<S309>/D'
                                        */
  real_T A_Gain_e;                     /* Expression: sps.A
                                        * Referenced by: '<S310>/A'
                                        */
  real_T B_Gain_h;                     /* Expression: sps.B
                                        * Referenced by: '<S310>/B'
                                        */
  real_T UnitDelay1_InitialCondition_i;/* Expression: 0
                                        * Referenced by: '<S67>/Unit Delay1'
                                        */
  real_T UnitDelay4_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S67>/Unit Delay4'
                                        */
  real_T Avoiddivbyzero_UpperSat;      /* Expression: inf
                                        * Referenced by: '<S80>/Avoid div. by zero'
                                        */
  real_T Avoiddivbyzero_LowerSat;      /* Expression: eps
                                        * Referenced by: '<S80>/Avoid div. by zero'
                                        */
  real_T DiscreteTimeIntegrator_gainva_m;/* Computed Parameter: DiscreteTimeIntegrator_gainva_m
                                          * Referenced by: '<S69>/Discrete-Time Integrator'
                                          */
  real_T Avoiddivbyzero_UpperSat_l;    /* Expression: inf
                                        * Referenced by: '<S80>/Avoid div. by zero '
                                        */
  real_T Avoiddivbyzero_LowerSat_j;    /* Expression: eps
                                        * Referenced by: '<S80>/Avoid div. by zero '
                                        */
  real_T DiscreteTimeIntegrator1_gainv_i;/* Computed Parameter: DiscreteTimeIntegrator1_gainv_i
                                          * Referenced by: '<S69>/Discrete-Time Integrator1'
                                          */
  real_T Mode_Value;                   /* Expression: Display
                                        * Referenced by: '<S76>/Mode'
                                        */
  real_T PowerbasefortheGenerator_Gain;/* Expression: Pmec/Pnom
                                        * Referenced by: '<S8>/Power base for the Generator'
                                        */
  real_T N_Gain;                       /* Expression: SM.N
                                        * Referenced by: '<S327>/N'
                                        */
  real_T u_Gain[2];                    /* Expression: [1 -1]
                                        * Referenced by: '<S338>/1-1'
                                        */
  real_T _Vb_Gain;                     /* Expression: 1/SM.Vb
                                        * Referenced by: '<S332>/1_Vb'
                                        */
  real_T Vkd0Vkq10Vkq20_Value[2];      /* Expression: zeros(1, SM.nState-3)
                                        * Referenced by: '<S327>/[ Vkd =0 Vkq1=0  Vkq2=0 ]'
                                        */
  real_T voltages_InitialCondition;    /* Expression: 0
                                        * Referenced by: '<S342>/voltages'
                                        */
  real_T IC_Threshold;                 /* Expression: Ts
                                        * Referenced by: '<S342>/IC'
                                        */
  real_T _Pb_Gain;                     /* Expression: 1/SM.Pb
                                        * Referenced by: '<S329>/1_Pb'
                                        */
  real_T F_Gain;                       /* Expression: SM.F
                                        * Referenced by: '<S329>/F'
                                        */
  real_T uH_Gain;                      /* Expression: 1/(2*SM.H)
                                        * Referenced by: '<S329>/1 ----- 2H'
                                        */
  real_T Rotorspeeddeviationdw_gainval;/* Computed Parameter: Rotorspeeddeviationdw_gainval
                                        * Referenced by: '<S329>/Rotor speed deviation (dw)'
                                        */
  real_T Rotorspeeddeviationdw_IC;     /* Expression: SM.dwo
                                        * Referenced by: '<S329>/Rotor speed deviation (dw)'
                                        */
  real_T webase_Gain;                  /* Expression: SM.web
                                        * Referenced by: '<S329>/we base'
                                        */
  real_T web3_Gain;                    /* Expression: SM.web
                                        * Referenced by: '<S329>/web3'
                                        */
  real_T In_ac_switch_Value;           /* Expression: 1
                                        * Referenced by: '<S8>/In_ac_switch '
                                        */
  real_T In_dc_switch1_Value;          /* Expression: 1
                                        * Referenced by: '<S8>/In_dc_switch 1'
                                        */
  real_T g_Value[6];                   /* Expression: zeros(1,6)
                                        * Referenced by: '<S78>/g'
                                        */
  real_T donotdeletethisgain_Gain_a;   /* Expression: 1
                                        * Referenced by: '<S34>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_l;   /* Expression: 1
                                        * Referenced by: '<S35>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_h0;  /* Expression: 1
                                        * Referenced by: '<S36>/do not delete this gain'
                                        */
  real_T Kv1_Gain_p;                   /* Expression: Kv
                                        * Referenced by: '<S5>/Kv1'
                                        */
  real_T donotdeletethisgain_Gain_mh;  /* Expression: 1
                                        * Referenced by: '<S31>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_gl;  /* Expression: 1
                                        * Referenced by: '<S32>/do not delete this gain'
                                        */
  real_T donotdeletethisgain_Gain_b;   /* Expression: 1
                                        * Referenced by: '<S33>/do not delete this gain'
                                        */
  real_T Kv_Gain_b;                    /* Expression: Ki
                                        * Referenced by: '<S5>/Kv'
                                        */
  boolean_T selector_Y0;               /* Computed Parameter: selector_Y0
                                        * Referenced by: '<S12>/selector'
                                        */
  boolean_T Constant1_Value_k;         /* Expression: SM.nState==6
                                        * Referenced by: '<S341>/Constant1'
                                        */
  boolean_T Constant2_Value_j;         /* Expression: SM.nState==6
                                        * Referenced by: '<S341>/Constant2'
                                        */
  boolean_T _Value_e;                  /* Expression: Tf_sps>0 | Tt_sps>0
                                        * Referenced by: '<S353>/2'
                                        */
  P_Subsystem1_Wind_songweiwei_T Subsystem1;/* '<S274>/Subsystem1' */
  P_Subsystempi2delay_Wind_song_T Subsystempi2delay;/* '<S274>/Subsystem - pi//2 delay' */
  P_Subsystem1_Wind_songweiwei_T Subsystem1_p;/* '<S268>/Subsystem1' */
  P_Subsystempi2delay_Wind_song_T Subsystempi2delay_jr;/* '<S268>/Subsystem - pi//2 delay' */
  P_Subsystem1_Wind_songweiwei_T Subsystem1_l;/* '<S262>/Subsystem1' */
  P_Subsystempi2delay_Wind_song_T Subsystempi2delay_j;/* '<S262>/Subsystem - pi//2 delay' */
  P_Subsystem1_Wind_songweiwei_T Subsystem1_l1;/* '<S254>/Subsystem1' */
  P_Subsystempi2delay_Wind_song_T Subsystempi2delay_p;/* '<S254>/Subsystem - pi//2 delay' */
  P_Subsystem1_Wind_songweiwei_T Subsystem1_lp;/* '<S241>/Subsystem1' */
  P_Subsystempi2delay_Wind_song_T Subsystempi2delay_dg;/* '<S241>/Subsystem - pi//2 delay' */
  P_Subsystem1_Wind_songweiwei_T Subsystem1_oh;/* '<S174>/Subsystem1' */
  P_Subsystempi2delay_Wind_song_T Subsystempi2delay_h;/* '<S174>/Subsystem - pi//2 delay' */
  P_Subsystem1_Wind_songweiwei_T Subsystem1_o;/* '<S160>/Subsystem1' */
  P_Subsystempi2delay_Wind_song_T Subsystempi2delay_k;/* '<S160>/Subsystem - pi//2 delay' */
};

/* Real-time Model Data Structure */
struct tag_RTM_Wind_songweiwei_T {
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
    SimStruct childSFunctions[3];
    SimStruct *childSFunctionPtrs[3];
    struct _ssBlkInfo2 blkInfo2[3];
    struct _ssSFcnModelMethods2 methods2[3];
    struct _ssSFcnModelMethods3 methods3[3];
    struct _ssStatesInfo2 statesInfo2[3];
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
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
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
extern P_Wind_songweiwei_T Wind_songweiwei_P;

/* Block signals (auto storage) */
extern B_Wind_songweiwei_T Wind_songweiwei_B;

/* Block states (auto storage) */
extern DW_Wind_songweiwei_T Wind_songweiwei_DW;

/* External outputs (root outports fed by signals with auto storage) */
extern ExtY_Wind_songweiwei_T Wind_songweiwei_Y;

/* Model entry point functions */
extern void Wind_songweiwei_initialize(void);
extern void Wind_songweiwei_step(void);
extern void Wind_songweiwei_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Wind_songweiwei_T *const Wind_songweiwei_M;

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
 * '<Root>' : 'Wind_songweiwei'
 * '<S1>'   : 'Wind_songweiwei/120 kV'
 * '<S2>'   : 'Wind_songweiwei/120 kV//25 kV 47 MVA'
 * '<S3>'   : 'Wind_songweiwei/25 kV// 575 V 5*2.5 MVA'
 * '<S4>'   : 'Wind_songweiwei/B120 (120 kV)'
 * '<S5>'   : 'Wind_songweiwei/B25 (25 kV)'
 * '<S6>'   : 'Wind_songweiwei/B575 (575 V)'
 * '<S7>'   : 'Wind_songweiwei/Grounding Transformer '
 * '<S8>'   : 'Wind_songweiwei/Wind Turbine Type 4'
 * '<S9>'   : 'Wind_songweiwei/powergui'
 * '<S10>'  : 'Wind_songweiwei/120 kV/Model'
 * '<S11>'  : 'Wind_songweiwei/120 kV/Model/Harmonic Generator'
 * '<S12>'  : 'Wind_songweiwei/120 kV/Model/Signal generator'
 * '<S13>'  : 'Wind_songweiwei/120 kV/Model/Harmonic Generator/Harmonic A generation'
 * '<S14>'  : 'Wind_songweiwei/120 kV/Model/Harmonic Generator/Harmonic B generation'
 * '<S15>'  : 'Wind_songweiwei/120 kV/Model/Signal generator/Timer'
 * '<S16>'  : 'Wind_songweiwei/120 kV/Model/Signal generator/Variation SubSystem'
 * '<S17>'  : 'Wind_songweiwei/120 kV//25 kV 47 MVA/Model'
 * '<S18>'  : 'Wind_songweiwei/120 kV//25 kV 47 MVA/Model/Linear'
 * '<S19>'  : 'Wind_songweiwei/120 kV//25 kV 47 MVA/Model/Linear1'
 * '<S20>'  : 'Wind_songweiwei/120 kV//25 kV 47 MVA/Model/Linear2'
 * '<S21>'  : 'Wind_songweiwei/25 kV// 575 V 5*2.5 MVA/Model'
 * '<S22>'  : 'Wind_songweiwei/25 kV// 575 V 5*2.5 MVA/Model/Linear'
 * '<S23>'  : 'Wind_songweiwei/25 kV// 575 V 5*2.5 MVA/Model/Linear1'
 * '<S24>'  : 'Wind_songweiwei/25 kV// 575 V 5*2.5 MVA/Model/Linear2'
 * '<S25>'  : 'Wind_songweiwei/B120 (120 kV)/Mode I'
 * '<S26>'  : 'Wind_songweiwei/B120 (120 kV)/Mode V'
 * '<S27>'  : 'Wind_songweiwei/B120 (120 kV)/Model'
 * '<S28>'  : 'Wind_songweiwei/B25 (25 kV)/Mode I'
 * '<S29>'  : 'Wind_songweiwei/B25 (25 kV)/Mode V'
 * '<S30>'  : 'Wind_songweiwei/B25 (25 kV)/Model'
 * '<S31>'  : 'Wind_songweiwei/B25 (25 kV)/Model/I A:'
 * '<S32>'  : 'Wind_songweiwei/B25 (25 kV)/Model/I B:'
 * '<S33>'  : 'Wind_songweiwei/B25 (25 kV)/Model/I C:'
 * '<S34>'  : 'Wind_songweiwei/B25 (25 kV)/Model/U A:'
 * '<S35>'  : 'Wind_songweiwei/B25 (25 kV)/Model/U B:'
 * '<S36>'  : 'Wind_songweiwei/B25 (25 kV)/Model/U C:'
 * '<S37>'  : 'Wind_songweiwei/B25 (25 kV)/Model/I A:/Model'
 * '<S38>'  : 'Wind_songweiwei/B25 (25 kV)/Model/I B:/Model'
 * '<S39>'  : 'Wind_songweiwei/B25 (25 kV)/Model/I C:/Model'
 * '<S40>'  : 'Wind_songweiwei/B25 (25 kV)/Model/U A:/Model'
 * '<S41>'  : 'Wind_songweiwei/B25 (25 kV)/Model/U B:/Model'
 * '<S42>'  : 'Wind_songweiwei/B25 (25 kV)/Model/U C:/Model'
 * '<S43>'  : 'Wind_songweiwei/B575 (575 V)/Mode I'
 * '<S44>'  : 'Wind_songweiwei/B575 (575 V)/Mode V'
 * '<S45>'  : 'Wind_songweiwei/B575 (575 V)/Model'
 * '<S46>'  : 'Wind_songweiwei/B575 (575 V)/Model/I A:'
 * '<S47>'  : 'Wind_songweiwei/B575 (575 V)/Model/I B:'
 * '<S48>'  : 'Wind_songweiwei/B575 (575 V)/Model/I C:'
 * '<S49>'  : 'Wind_songweiwei/B575 (575 V)/Model/U A:'
 * '<S50>'  : 'Wind_songweiwei/B575 (575 V)/Model/U B:'
 * '<S51>'  : 'Wind_songweiwei/B575 (575 V)/Model/U C:'
 * '<S52>'  : 'Wind_songweiwei/B575 (575 V)/Model/I A:/Model'
 * '<S53>'  : 'Wind_songweiwei/B575 (575 V)/Model/I B:/Model'
 * '<S54>'  : 'Wind_songweiwei/B575 (575 V)/Model/I C:/Model'
 * '<S55>'  : 'Wind_songweiwei/B575 (575 V)/Model/U A:/Model'
 * '<S56>'  : 'Wind_songweiwei/B575 (575 V)/Model/U B:/Model'
 * '<S57>'  : 'Wind_songweiwei/B575 (575 V)/Model/U C:/Model'
 * '<S58>'  : 'Wind_songweiwei/Grounding Transformer /T1'
 * '<S59>'  : 'Wind_songweiwei/Grounding Transformer /T2'
 * '<S60>'  : 'Wind_songweiwei/Grounding Transformer /T3'
 * '<S61>'  : 'Wind_songweiwei/Grounding Transformer /T1/Model'
 * '<S62>'  : 'Wind_songweiwei/Grounding Transformer /T2/Model'
 * '<S63>'  : 'Wind_songweiwei/Grounding Transformer /T3/Model'
 * '<S64>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen'
 * '<S65>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_grid'
 * '<S66>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_grid_conv'
 * '<S67>'  : 'Wind_songweiwei/Wind Turbine Type 4/Control'
 * '<S68>'  : 'Wind_songweiwei/Wind Turbine Type 4/Diode'
 * '<S69>'  : 'Wind_songweiwei/Wind Turbine Type 4/Drive Train'
 * '<S70>'  : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch'
 * '<S71>'  : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch1'
 * '<S72>'  : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch2'
 * '<S73>'  : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch3'
 * '<S74>'  : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch4'
 * '<S75>'  : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch5'
 * '<S76>'  : 'Wind_songweiwei/Wind Turbine Type 4/Multimeter'
 * '<S77>'  : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard'
 * '<S78>'  : 'Wind_songweiwei/Wind Turbine Type 4/Universal Bridge'
 * '<S79>'  : 'Wind_songweiwei/Wind Turbine Type 4/Universal Bridge1'
 * '<S80>'  : 'Wind_songweiwei/Wind Turbine Type 4/Wind Turbine'
 * '<S81>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Mode I'
 * '<S82>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Mode V'
 * '<S83>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Model'
 * '<S84>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Model/I A:'
 * '<S85>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Model/I B:'
 * '<S86>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Model/I C:'
 * '<S87>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Model/U AB:'
 * '<S88>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Model/U BC:'
 * '<S89>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Model/U CA:'
 * '<S90>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Model/I A:/Model'
 * '<S91>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Model/I B:/Model'
 * '<S92>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Model/I C:/Model'
 * '<S93>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Model/U AB:/Model'
 * '<S94>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Model/U BC:/Model'
 * '<S95>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_gen/Model/U CA:/Model'
 * '<S96>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Mode I'
 * '<S97>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Mode V'
 * '<S98>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Model'
 * '<S99>'  : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Model/I A:'
 * '<S100>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Model/I B:'
 * '<S101>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Model/I C:'
 * '<S102>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Model/U A:'
 * '<S103>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Model/U B:'
 * '<S104>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Model/U C:'
 * '<S105>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Model/I A:/Model'
 * '<S106>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Model/I B:/Model'
 * '<S107>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Model/I C:/Model'
 * '<S108>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Model/U A:/Model'
 * '<S109>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Model/U B:/Model'
 * '<S110>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid/Model/U C:/Model'
 * '<S111>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid_conv/Mode I'
 * '<S112>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid_conv/Mode V'
 * '<S113>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid_conv/Model'
 * '<S114>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid_conv/Model/U A:'
 * '<S115>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid_conv/Model/U B:'
 * '<S116>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid_conv/Model/U C:'
 * '<S117>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid_conv/Model/U A:/Model'
 * '<S118>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid_conv/Model/U B:/Model'
 * '<S119>' : 'Wind_songweiwei/Wind Turbine Type 4/B_grid_conv/Model/U C:/Model'
 * '<S120>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Excitation Control system'
 * '<S121>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system'
 * '<S122>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation'
 * '<S123>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)'
 * '<S124>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Pulse Generator'
 * '<S125>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Speed & Pitch Control'
 * '<S126>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Excitation Control system/Discrete PI Controller'
 * '<S127>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/Cartesian to Polar'
 * '<S128>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/Discrete PI Controller'
 * '<S129>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/Discrete PI Controller1'
 * '<S130>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/First-Order Filter'
 * '<S131>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/INT'
 * '<S132>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/Polar to Cartesian'
 * '<S133>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/dq0 to abc'
 * '<S134>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/First-Order Filter/Model'
 * '<S135>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/INT/S6'
 * '<S136>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/dq0 to abc/Alpha-Beta-Zero to abc'
 * '<S137>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/dq0 to abc/dq0 to Alpha-Beta-Zero'
 * '<S138>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/dq0 to abc/dq0 to Alpha-Beta-Zero/Compare To Constant'
 * '<S139>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/dq0 to abc/dq0 to Alpha-Beta-Zero/Compare To Constant1'
 * '<S140>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/dq0 to abc/dq0 to Alpha-Beta-Zero/Subsystem - pi//2 delay'
 * '<S141>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Grid-side converter Control system/dq0 to abc/dq0 to Alpha-Beta-Zero/Subsystem1'
 * '<S142>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power'
 * '<S143>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter'
 * '<S144>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)'
 * '<S145>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Stator flux estimator'
 * '<S146>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq0'
 * '<S147>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq1'
 * '<S148>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq2'
 * '<S149>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc2qd'
 * '<S150>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value'
 * '<S151>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2'
 * '<S152>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value/Cartesian to Polar'
 * '<S153>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value/Discrete Variable Frequency Mean value1'
 * '<S154>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value/Discrete Variable Frequency Mean value2'
 * '<S155>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value/abc to dq0'
 * '<S156>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value/Discrete Variable Frequency Mean value1/Correction subsystem'
 * '<S157>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value/Discrete Variable Frequency Mean value1/Discrete Variable Transport Delay'
 * '<S158>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value/Discrete Variable Frequency Mean value2/Correction subsystem'
 * '<S159>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value/Discrete Variable Frequency Mean value2/Discrete Variable Transport Delay'
 * '<S160>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value/abc to dq0/Alpha-Beta-Zero to dq0'
 * '<S161>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value/abc to dq0/abc to Alpha-Beta-Zero'
 * '<S162>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S163>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S164>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S165>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S166>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2/Cartesian to Polar'
 * '<S167>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2/Discrete Variable Frequency Mean value1'
 * '<S168>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2/Discrete Variable Frequency Mean value2'
 * '<S169>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2/abc to dq0'
 * '<S170>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2/Discrete Variable Frequency Mean value1/Correction subsystem'
 * '<S171>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2/Discrete Variable Frequency Mean value1/Discrete Variable Transport Delay'
 * '<S172>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2/Discrete Variable Frequency Mean value2/Correction subsystem'
 * '<S173>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2/Discrete Variable Frequency Mean value2/Discrete Variable Transport Delay'
 * '<S174>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2/abc to dq0/Alpha-Beta-Zero to dq0'
 * '<S175>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2/abc to dq0/abc to Alpha-Beta-Zero'
 * '<S176>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S177>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S178>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S179>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Discrete 3-phase PLL-Driven  Positive-Sequence Active & Reactive Power/Discrete 3-phase PLL-Driven Positive-Sequence Fundamental Value2/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S180>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter'
 * '<S181>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter1'
 * '<S182>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter2'
 * '<S183>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter3'
 * '<S184>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter4'
 * '<S185>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter5'
 * '<S186>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter6'
 * '<S187>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter7'
 * '<S188>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter8'
 * '<S189>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter/Model'
 * '<S190>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter/Model/A*k(k-1)'
 * '<S191>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter/Model/B*(u(k)+u(k-1))'
 * '<S192>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter/Model/C*x(k)'
 * '<S193>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter1/Model'
 * '<S194>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter1/Model/A*k(k-1)'
 * '<S195>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter1/Model/B*(u(k)+u(k-1))'
 * '<S196>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter1/Model/C*x(k)'
 * '<S197>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter2/Model'
 * '<S198>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter2/Model/A*k(k-1)'
 * '<S199>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter2/Model/B*(u(k)+u(k-1))'
 * '<S200>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter2/Model/C*x(k)'
 * '<S201>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter3/Model'
 * '<S202>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter3/Model/A*k(k-1)'
 * '<S203>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter3/Model/B*(u(k)+u(k-1))'
 * '<S204>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter3/Model/C*x(k)'
 * '<S205>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter4/Model'
 * '<S206>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter4/Model/A*k(k-1)'
 * '<S207>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter4/Model/B*(u(k)+u(k-1))'
 * '<S208>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter4/Model/C*x(k)'
 * '<S209>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter5/Model'
 * '<S210>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter5/Model/A*k(k-1)'
 * '<S211>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter5/Model/B*(u(k)+u(k-1))'
 * '<S212>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter5/Model/C*x(k)'
 * '<S213>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter6/Model'
 * '<S214>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter6/Model/A*k(k-1)'
 * '<S215>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter6/Model/B*(u(k)+u(k-1))'
 * '<S216>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter6/Model/C*x(k)'
 * '<S217>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter7/Model'
 * '<S218>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter7/Model/A*k(k-1)'
 * '<S219>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter7/Model/B*(u(k)+u(k-1))'
 * '<S220>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter7/Model/C*x(k)'
 * '<S221>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter8/Model'
 * '<S222>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter8/Model/A*k(k-1)'
 * '<S223>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter8/Model/B*(u(k)+u(k-1))'
 * '<S224>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Filter/Second-Order Filter8/Model/C*x(k)'
 * '<S225>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model'
 * '<S226>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control'
 * '<S227>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Discrete'
 * '<S228>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Second-Order Filter'
 * '<S229>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Variable Frequency Mean value'
 * '<S230>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/abc to dq0'
 * '<S231>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)'
 * '<S232>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/Mean (Variable Frequency)1'
 * '<S233>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/Mean (Variable Frequency)2'
 * '<S234>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/abc ->dq0'
 * '<S235>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/Mean (Variable Frequency)1/Model'
 * '<S236>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/Mean (Variable Frequency)1/Model/Correction subsystem'
 * '<S237>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/Mean (Variable Frequency)1/Model/Discrete Variable Time Delay'
 * '<S238>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/Mean (Variable Frequency)2/Model'
 * '<S239>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/Mean (Variable Frequency)2/Model/Correction subsystem'
 * '<S240>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/Mean (Variable Frequency)2/Model/Discrete Variable Time Delay'
 * '<S241>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/abc ->dq0/Alpha-Beta-Zero to dq0'
 * '<S242>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/abc ->dq0/abc to Alpha-Beta-Zero'
 * '<S243>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/abc ->dq0/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S244>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/abc ->dq0/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S245>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/abc ->dq0/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S246>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Automatic Gain Control/Positive-Sequence (PLL-Driven)/abc ->dq0/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S247>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Second-Order Filter/Model'
 * '<S248>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Second-Order Filter/Model/A*k(k-1)'
 * '<S249>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Second-Order Filter/Model/B*(u(k)+u(k-1))'
 * '<S250>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Second-Order Filter/Model/C*x(k)'
 * '<S251>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Variable Frequency Mean value/Model'
 * '<S252>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Variable Frequency Mean value/Model/Correction subsystem'
 * '<S253>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/Variable Frequency Mean value/Model/Discrete Variable Time Delay'
 * '<S254>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/abc to dq0/Alpha-Beta-Zero to dq0'
 * '<S255>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/abc to dq0/abc to Alpha-Beta-Zero'
 * '<S256>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S257>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S258>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S259>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/PLL (3ph)/Model/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S260>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Stator flux estimator/First-Order Filter'
 * '<S261>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/Stator flux estimator/First-Order Filter/Model'
 * '<S262>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq0/Alpha-Beta-Zero to dq0'
 * '<S263>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq0/abc to Alpha-Beta-Zero'
 * '<S264>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S265>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq0/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S266>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S267>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq0/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S268>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq1/Alpha-Beta-Zero to dq0'
 * '<S269>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq1/abc to Alpha-Beta-Zero'
 * '<S270>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq1/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S271>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq1/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S272>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq1/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S273>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq1/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S274>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq2/Alpha-Beta-Zero to dq0'
 * '<S275>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq2/abc to Alpha-Beta-Zero'
 * '<S276>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq2/Alpha-Beta-Zero to dq0/Compare To Constant'
 * '<S277>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq2/Alpha-Beta-Zero to dq0/Compare To Constant1'
 * '<S278>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq2/Alpha-Beta-Zero to dq0/Subsystem - pi//2 delay'
 * '<S279>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Measurement and Transformation/abc to dq2/Alpha-Beta-Zero to dq0/Subsystem1'
 * '<S280>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Cr_MinMax'
 * '<S281>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Modulator type'
 * '<S282>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Reference signal'
 * '<S283>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling'
 * '<S284>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Modulator type/Full Bridge Bipolar'
 * '<S285>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Modulator type/Full Bridge Unipolar'
 * '<S286>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Modulator type/One Three Phase Bridge'
 * '<S287>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Reference signal/External'
 * '<S288>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Reference signal/Internal'
 * '<S289>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Sync Asymmetrical'
 * '<S290>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Sync Natural'
 * '<S291>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Sync Symmetrical'
 * '<S292>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Unsync Asymmetrical'
 * '<S293>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Unsync Natural'
 * '<S294>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Unsync Symmetrical'
 * '<S295>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Sync Asymmetrical/Sample & Hold'
 * '<S296>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Sync Natural/Sync_NaturalSampling'
 * '<S297>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Sync Symmetrical/Sync_SymmetricalSampling'
 * '<S298>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Sync Symmetrical/Sync_SymmetricalSampling/Sample & Hold'
 * '<S299>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Unsync Asymmetrical/Unsync_AsymmetricalSampling'
 * '<S300>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Unsync Natural/Unsync_NaturalSampling'
 * '<S301>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Unsync Natural/Unsync_NaturalSampling/Triangle Generator'
 * '<S302>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Unsync Natural/Unsync_NaturalSampling/Triangle Generator/Model'
 * '<S303>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/PWM Generator (2-Level)/Sampling/Unsync Symmetrical/Unsync_SymmetricalSampling'
 * '<S304>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Speed & Pitch Control/Discrete PI Controller'
 * '<S305>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Speed & Pitch Control/Discrete PI Controller1'
 * '<S306>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Speed & Pitch Control/Discrete PI Controller2'
 * '<S307>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Speed & Pitch Control/First-Order Filter'
 * '<S308>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Speed & Pitch Control/First-Order Filter1'
 * '<S309>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Speed & Pitch Control/First-Order Filter/Model'
 * '<S310>' : 'Wind_songweiwei/Wind Turbine Type 4/Control/Speed & Pitch Control/First-Order Filter1/Model'
 * '<S311>' : 'Wind_songweiwei/Wind Turbine Type 4/Diode/Model'
 * '<S312>' : 'Wind_songweiwei/Wind Turbine Type 4/Diode/Model/Measurement list'
 * '<S313>' : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch/Model'
 * '<S314>' : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch/Model/Measurement list'
 * '<S315>' : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch1/Model'
 * '<S316>' : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch1/Model/Measurement list'
 * '<S317>' : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch2/Model'
 * '<S318>' : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch2/Model/Measurement list'
 * '<S319>' : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch3/Model'
 * '<S320>' : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch3/Model/Measurement list'
 * '<S321>' : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch4/Model'
 * '<S322>' : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch4/Model/Measurement list'
 * '<S323>' : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch5/Model'
 * '<S324>' : 'Wind_songweiwei/Wind Turbine Type 4/Ideal Switch5/Model/Measurement list'
 * '<S325>' : 'Wind_songweiwei/Wind Turbine Type 4/Multimeter/Model'
 * '<S326>' : 'Wind_songweiwei/Wind Turbine Type 4/Multimeter/StoreData'
 * '<S327>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model'
 * '<S328>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Measurements'
 * '<S329>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Mechanical model'
 * '<S330>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Machine measurements'
 * '<S331>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Synchronous Machine Discrete Model'
 * '<S332>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/abc to qd transformation'
 * '<S333>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/qd to abc transformation'
 * '<S334>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Machine measurements/Delta angle'
 * '<S335>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Machine measurements/PQ'
 * '<S336>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Machine measurements/Delta angle/Cartesian to Polar'
 * '<S337>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Machine measurements/Delta angle/Radians to Degrees'
 * '<S338>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Synchronous Machine Discrete Model/Electromagnetic Torque'
 * '<S339>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Synchronous Machine Discrete Model/Matrix W'
 * '<S340>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Synchronous Machine Discrete Model/Phimqd'
 * '<S341>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Synchronous Machine Discrete Model/Saturation'
 * '<S342>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Synchronous Machine Discrete Model/phi'
 * '<S343>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmd_sat'
 * '<S344>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmq_sat'
 * '<S345>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Synchronous Machine Discrete Model/Saturation/Update Matrix L'
 * '<S346>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmd_sat/Lad'
 * '<S347>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmd_sat/phimd'
 * '<S348>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmq_sat/Laq'
 * '<S349>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Synchronous Machine Discrete Model/Saturation/Lmq_sat/phimq'
 * '<S350>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Electrical model/Synchronous Machine Discrete Model/phi/Subsystem'
 * '<S351>' : 'Wind_songweiwei/Wind Turbine Type 4/Synchronous Machine pu Standard/Mechanical model/Delay Prediction'
 * '<S352>' : 'Wind_songweiwei/Wind Turbine Type 4/Universal Bridge/Model'
 * '<S353>' : 'Wind_songweiwei/Wind Turbine Type 4/Universal Bridge1/Model'
 * '<S354>' : 'Wind_songweiwei/Wind Turbine Type 4/Universal Bridge1/Model/Tail'
 * '<S355>' : 'Wind_songweiwei/Wind Turbine Type 4/Universal Bridge1/Model/Vf 1'
 * '<S356>' : 'Wind_songweiwei/Wind Turbine Type 4/Universal Bridge1/Model/Tail/y=f(t)'
 * '<S357>' : 'Wind_songweiwei/powergui/EquivalentModel1'
 * '<S358>' : 'Wind_songweiwei/powergui/EquivalentModel1/Gates'
 * '<S359>' : 'Wind_songweiwei/powergui/EquivalentModel1/Sources'
 * '<S360>' : 'Wind_songweiwei/powergui/EquivalentModel1/Status'
 * '<S361>' : 'Wind_songweiwei/powergui/EquivalentModel1/Yout'
 */
#endif                                 /* RTW_HEADER_Wind_songweiwei_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
