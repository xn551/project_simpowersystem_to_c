/*
 * File: Wind_songweiwei_data.c
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

/* Block parameters (auto storage) */
P_Wind_songweiwei_T Wind_songweiwei_P = {
  2.0E-6,                              /* Variable: Ts
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
  2.0,                                 /* Mask Parameter: AlphaBetaZerotodq0_Alignment
                                        * Referenced by: '<S241>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: AlphaBetaZerotodq0_Alignment_i
                                        * Referenced by: '<S274>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: AlphaBetaZerotodq0_Alignment_n
                                        * Referenced by: '<S262>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: AlphaBetaZerotodq0_Alignment_iw
                                        * Referenced by: '<S160>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: AlphaBetaZerotodq0_Alignment_m
                                        * Referenced by: '<S268>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: AlphaBetaZerotodq0_Alignment_l
                                        * Referenced by: '<S174>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: dq0toAlphaBetaZero_Alignment
                                        * Referenced by: '<S137>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: AlphaBetaZerotodq0_Alignment_f
                                        * Referenced by: '<S254>/Constant'
                                        */
  1.5,                                 /* Mask Parameter: DriveTrain_D_mutual
                                        * Referenced by: '<S69>/Mutual damping'
                                        */
  60.0,                                /* Mask Parameter: Discrete3phasePLLDrivenPositive
                                        * Referenced by:
                                        *   '<S153>/Step'
                                        *   '<S154>/Step'
                                        *   '<S167>/Step'
                                        *   '<S168>/Step'
                                        */
  2000.0,                              /* Mask Parameter: PulseGenerator_Freq_sawtooth
                                        * Referenced by: '<S124>/Nb_of_samples_by_cycle'
                                        */
  4.32,                                /* Mask Parameter: DriveTrain_H_WT
                                        * Referenced by: '<S69>/1_2H_WT'
                                        */
  0.0,                                 /* Mask Parameter: u0kV_HarmonicGeneration
                                        * Referenced by: '<S10>/valp_nom7'
                                        */
  1.0,                                 /* Mask Parameter: DiscretePIController_Init
                                        * Referenced by: '<S126>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Mask Parameter: DiscretePIController_Init_n
                                        * Referenced by: '<S128>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Mask Parameter: DiscretePIController1_Init
                                        * Referenced by: '<S129>/Discrete-Time Integrator'
                                        */
  376.99111843077515,                  /* Mask Parameter: Discrete_Init
                                        * Referenced by: '<S227>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Mask Parameter: DiscretePIController1_Init_i
                                        * Referenced by: '<S305>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Mask Parameter: DiscretePIController_Init_a
                                        * Referenced by: '<S304>/Discrete-Time Integrator'
                                        */
  0.5,                                 /* Mask Parameter: DiscretePIController2_Init
                                        * Referenced by: '<S306>/Discrete-Time Integrator'
                                        */
  1.0,                                 /* Mask Parameter: Discrete_Kd
                                        * Referenced by: '<S227>/Discrete Derivative '
                                        */
  20.0,                                /* Mask Parameter: DiscretePIController_Ki
                                        * Referenced by: '<S126>/Kp5'
                                        */
  27.5,                                /* Mask Parameter: DiscretePIController_Ki_f
                                        * Referenced by: '<S128>/Kp5'
                                        */
  50.0,                                /* Mask Parameter: DiscretePIController1_Ki
                                        * Referenced by: '<S129>/Kp5'
                                        */
  1.0,                                 /* Mask Parameter: DiscretePIController_Ki_p
                                        * Referenced by: '<S304>/Kp5'
                                        */
  6.0,                                 /* Mask Parameter: DiscretePIController1_Ki_j
                                        * Referenced by: '<S305>/Kp5'
                                        */
  100.0,                               /* Mask Parameter: DiscretePIController2_Ki
                                        * Referenced by: '<S306>/Kp5'
                                        */
  10.0,                                /* Mask Parameter: DiscretePIController_Kp
                                        * Referenced by: '<S126>/Kp4'
                                        */
  1.1,                                 /* Mask Parameter: DiscretePIController_Kp_h
                                        * Referenced by: '<S128>/Kp4'
                                        */
  1.0,                                 /* Mask Parameter: DiscretePIController1_Kp
                                        * Referenced by: '<S129>/Kp4'
                                        */
  60.0,                                /* Mask Parameter: Discrete_Kp
                                        * Referenced by: '<S227>/Kp4'
                                        */
  1.5,                                 /* Mask Parameter: DiscretePIController1_Kp_p
                                        * Referenced by: '<S305>/Kp4'
                                        */
  5.0,                                 /* Mask Parameter: DiscretePIController_Kp_c
                                        * Referenced by: '<S304>/Kp4'
                                        */
  0.025,                               /* Mask Parameter: DiscretePIController2_Kp
                                        * Referenced by: '<S306>/Kp4'
                                        */
  0.3,                                 /* Mask Parameter: DriveTrain_Ksh
                                        * Referenced by:
                                        *   '<S69>/Discrete-Time Integrator1'
                                        *   '<S69>/Stiffness'
                                        */
  -9797.9589711327117,                 /* Mask Parameter: HarmonicAgeneration_Mag_Harmo
                                        * Referenced by: '<S13>/Phase_Harmo1'
                                        */
  0.0,                                 /* Mask Parameter: HarmonicBgeneration_Mag_Harmo
                                        * Referenced by: '<S14>/Phase_Harmo1'
                                        */

  /*  Mask Parameter: PWMGenerator2Level_MinMax
   * Referenced by: '<S123>/Constant10'
   */
  { -1.0, 1.0 },
  0.0,                                 /* Mask Parameter: HarmonicAgeneration_Phase_Harmo
                                        * Referenced by: '<S13>/Phase_Harmo'
                                        */
  90.0,                                /* Mask Parameter: HarmonicBgeneration_Phase_Harmo
                                        * Referenced by: '<S14>/Phase_Harmo'
                                        */
  1.0E+7,                              /* Mask Parameter: WindTurbine_Prated
                                        * Referenced by: '<S80>/->pu'
                                        */
  0.006,                               /* Mask Parameter: WindTurbineType4_Resistance
                                        * Referenced by:
                                        *   '<S145>/Rs2'
                                        *   '<S145>/Rs4'
                                        */
  2.0,                                 /* Mask Parameter: HarmonicAgeneration_Seq_Harmo
                                        * Referenced by: '<S13>/Phase_Harmo2'
                                        */
  0.0,                                 /* Mask Parameter: HarmonicBgeneration_Seq_Harmo
                                        * Referenced by: '<S14>/Phase_Harmo2'
                                        */
  0.0,                                 /* Mask Parameter: Tail_Tf
                                        * Referenced by:
                                        *   '<S354>/2'
                                        *   '<S356>/Constant2'
                                        *   '<S356>/-0.9//Tf'
                                        */
  0.13333333333333333,                 /* Mask Parameter: VariationSubSystem_Toff_Variati
                                        * Referenced by: '<S16>/Step1'
                                        */
  0.05,                                /* Mask Parameter: VariationSubSystem_Ton_Variatio
                                        * Referenced by: '<S16>/Step'
                                        */
  0.0,                                 /* Mask Parameter: Tail_Tt
                                        * Referenced by:
                                        *   '<S354>/2'
                                        *   '<S356>/Constant2'
                                        *   '<S356>/0.1//Tt'
                                        */
  2.0,                                 /* Mask Parameter: u0kV_VariationEntity
                                        * Referenced by:
                                        *   '<S10>/valp_nom3'
                                        *   '<S12>/valp_nom3'
                                        */
  2.0,                                 /* Mask Parameter: VariationSubSystem_VariationFre
                                        * Referenced by: '<S16>/valp_nom9'
                                        */
  29393.876913398137,                  /* Mask Parameter: VariationSubSystem_VariationMag
                                        * Referenced by: '<S16>/valp_nom8'
                                        */
  0.0,                                 /* Mask Parameter: u0kV_VariationPhaseA
                                        * Referenced by: '<S10>/SinglePhase'
                                        */
  979795.89711327117,                  /* Mask Parameter: VariationSubSystem_VariationRat
                                        * Referenced by: '<S16>/valp_nom7'
                                        */
  -19595.917942265423,                 /* Mask Parameter: VariationSubSystem_VariationSte
                                        * Referenced by: '<S16>/valp_nom6'
                                        */
  4.0,                                 /* Mask Parameter: u0kV_VariationType
                                        * Referenced by: '<S12>/valp_nom5'
                                        */
  406.58639918226493,                  /* Mask Parameter: DiscreteVariableFrequencyMeanva
                                        * Referenced by: '<S154>/Unit Delay'
                                        */
  234.74276701672125,                  /* Mask Parameter: DiscreteVariableFrequencyMean_m
                                        * Referenced by: '<S153>/Unit Delay'
                                        */
  15538.012815566848,                  /* Mask Parameter: DiscreteVariableFrequencyMean_j
                                        * Referenced by: '<S168>/Unit Delay'
                                        */
  -2739.770885978337,                  /* Mask Parameter: DiscreteVariableFrequencyMean_l
                                        * Referenced by: '<S167>/Unit Delay'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const
                                        * Referenced by: '<S243>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: CompareToConstant1_const
                                        * Referenced by: '<S244>/Constant'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const_o
                                        * Referenced by: '<S276>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: CompareToConstant1_const_o
                                        * Referenced by: '<S277>/Constant'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const_c
                                        * Referenced by: '<S264>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: CompareToConstant1_const_b
                                        * Referenced by: '<S265>/Constant'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const_g
                                        * Referenced by: '<S162>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: CompareToConstant1_const_g
                                        * Referenced by: '<S163>/Constant'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const_i
                                        * Referenced by: '<S270>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: CompareToConstant1_const_k
                                        * Referenced by: '<S271>/Constant'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const_gf
                                        * Referenced by: '<S176>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: CompareToConstant1_const_gp
                                        * Referenced by: '<S177>/Constant'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const_ia
                                        * Referenced by: '<S138>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: CompareToConstant1_const_o0
                                        * Referenced by: '<S139>/Constant'
                                        */
  1.0,                                 /* Mask Parameter: CompareToConstant_const_oc
                                        * Referenced by: '<S256>/Constant'
                                        */
  2.0,                                 /* Mask Parameter: CompareToConstant1_const_j
                                        * Referenced by: '<S257>/Constant'
                                        */
  1.0,                                 /* Mask Parameter: HarmonicAgeneration_n_Harmo
                                        * Referenced by: '<S13>/Gain1'
                                        */
  0.0,                                 /* Mask Parameter: HarmonicBgeneration_n_Harmo
                                        * Referenced by: '<S14>/Gain1'
                                        */
  1.0,                                 /* Mask Parameter: DriveTrain_torque0
                                        * Referenced by: '<S69>/Discrete-Time Integrator1'
                                        */
  1.0,                                 /* Mask Parameter: DriveTrain_w_wt0
                                        * Referenced by: '<S69>/Discrete-Time Integrator'
                                        */
  376.99111843077515,                  /* Mask Parameter: DriveTrain_wbase
                                        * Referenced by: '<S69>/wbase'
                                        */

  /*  Expression: [0 2*pi/3 -2*pi/3]
   * Referenced by: '<S13>/Negative-sequence'
   */
  { 0.0, 2.0943951023931953, -2.0943951023931953 },

  /*  Expression: [0 -2*pi/3 2*pi/3]
   * Referenced by: '<S13>/Positive-sequence'
   */
  { 0.0, -2.0943951023931953, 2.0943951023931953 },

  /*  Expression: [0 0 0]
   * Referenced by: '<S13>/Zero-sequence'
   */
  { 0.0, 0.0, 0.0 },

  /*  Expression: [0 2*pi/3 -2*pi/3]
   * Referenced by: '<S14>/Negative-sequence'
   */
  { 0.0, 2.0943951023931953, -2.0943951023931953 },

  /*  Expression: [0 -2*pi/3 2*pi/3]
   * Referenced by: '<S14>/Positive-sequence'
   */
  { 0.0, -2.0943951023931953, 2.0943951023931953 },

  /*  Expression: [0 0 0]
   * Referenced by: '<S14>/Zero-sequence'
   */
  { 0.0, 0.0, 0.0 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S11>/Out1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S11>/Out2'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S13>/Gain3'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S13>/valp_nom2'
                                        */
  0.1,                                 /* Expression: Ton_Harmo
                                        * Referenced by: '<S11>/Step'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S11>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S11>/Step'
                                        */
  2.0,                                 /* Expression: Toff_Harmo
                                        * Referenced by: '<S11>/Step1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S11>/Step1'
                                        */
  -1.0,                                /* Expression: -1
                                        * Referenced by: '<S11>/Step1'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S14>/Gain3'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S14>/valp_nom2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/Constant5'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/Constant4'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S12>/Gain4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/Constant1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S16>/Constant1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S16>/Constant5'
                                        */
  6.2831853071795862,                  /* Expression: 2*pi
                                        * Referenced by: '<S16>/Gain1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S16>/Constant4'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S16>/Constant2'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S16>/Constant'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/timer'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/magnitude'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/frequency'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/phase'
                                        */

  /*  Expression: tv
   * Referenced by: '<S15>/Look-Up Table'
   */
  { -1.0, 0.0, 0.03, 0.03, 0.13, 0.13, 1.13 },

  /*  Expression: opv
   * Referenced by: '<S15>/Look-Up Table'
   */
  { 97979.58971132712, 97979.58971132712, 97979.58971132712, 73484.692283495344,
    73484.692283495344, 97979.58971132712, 97979.58971132712 },
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S12>/Constant'
                                        */
  3.0,                                 /* Expression: 3
                                        * Referenced by: '<S12>/Constant2'
                                        */
  4.0,                                 /* Expression: 4
                                        * Referenced by: '<S12>/Constant3'
                                        */
  4.0,                                 /* Expression: 4
                                        * Referenced by: '<S12>/Constant6'
                                        */
  2.0E-6,                              /* Computed Parameter: DiscreteTimeIntegrator_gainval
                                        * Referenced by: '<S12>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S12>/Discrete-Time Integrator'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S16>/Step1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S16>/Step1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S16>/Constant3'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S16>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S16>/Step'
                                        */
  2.0E-6,                              /* Computed Parameter: DiscreteTimeIntegrator1_gainval
                                        * Referenced by: '<S16>/Discrete-Time Integrator1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S16>/Discrete-Time Integrator1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S16>/Unit Delay'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S16>/Switch2'
                                        */
  6.2831853071795862,                  /* Expression: 2*pi
                                        * Referenced by: '<S12>/Gain2'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S16>/Switch'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S120>/Flux_ref'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S120>/avoid division by 0'
                                        */
  1.0E-6,                              /* Expression: 1e-6
                                        * Referenced by: '<S120>/avoid division by 0'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S120>/Flux_ref '
                                        */

  /*  Expression: [0,0]
   * Referenced by: '<S140>/alpha_beta'
   */
  { 0.0, 0.0 },

  /*  Expression: [0,0]
   * Referenced by: '<S141>/alpha_beta'
   */
  { 0.0, 0.0 },
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S156>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S158>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S170>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S172>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S236>/Gain1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S239>/Gain1'
                                        */
  1.0,                                 /* Expression: [1]
                                        * Referenced by: '<S226>/Gain'
                                        */

  /*  Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S242>/Gain3'
   */
  { 1.0, 0.0, 0.5, -0.5, 0.8660254037844386, 0.5, -0.5, -0.8660254037844386, 0.5
  },
  0.66666666666666663,                 /* Expression: 2/3
                                        * Referenced by: '<S242>/Gain1'
                                        */
  1.0E-6,                              /* Computed Parameter: Integ4_gainval
                                        * Referenced by: '<S235>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S235>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S235>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S235>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size
   * Referenced by: '<S237>/S-Function'
   */
  { 1.0, 1.0 },
  0.022224222222222222,                /* Expression: MaxDelay
                                        * Referenced by: '<S237>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size
   * Referenced by: '<S237>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-6,                              /* Expression: Ts
                                        * Referenced by: '<S237>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size
   * Referenced by: '<S237>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S237>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size
   * Referenced by: '<S237>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S237>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S236>/Unit Delay'
                                        */
  0.016666666666666666,                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S235>/Constant'
                                        */
  1.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S235>/Unit Delay1'
                                        */
  1.0E-6,                              /* Computed Parameter: Integ4_gainval_j
                                        * Referenced by: '<S238>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S238>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S238>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S238>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_j
   * Referenced by: '<S240>/S-Function'
   */
  { 1.0, 1.0 },
  0.022224222222222222,                /* Expression: MaxDelay
                                        * Referenced by: '<S240>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_j
   * Referenced by: '<S240>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-6,                              /* Expression: Ts
                                        * Referenced by: '<S240>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_e
   * Referenced by: '<S240>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S240>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_c
   * Referenced by: '<S240>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S240>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S239>/Unit Delay'
                                        */
  0.016666666666666666,                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S238>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S238>/Unit Delay1'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S226>/Saturation'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S226>/Saturation'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S252>/Gain1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S124>/Constant4'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S124>/Constant3'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S125>/wref'
                                        */

  /*  Expression: SM.Linv
   * Referenced by: '<S331>/Constant4'
   */
  { 4.11522633744856, 0.0, 0.0, 0.0, -3.0481171290895483, 0.0,
    3.9682539682539653, -2.2303134801555142, -1.4258943432976816, 0.0, 0.0,
    -2.2303134801555142, 3.992482413187421, -1.4457725858578245, 0.0, 0.0,
    -1.4258943432976816, -1.4457725858578245, 3.07394690714502, 0.0,
    -3.0481171290895483, 0.0, 0.0, 0.0, 4.63271878370818 },
  0.31187032639869727,                 /* Expression: SM.Lmq
                                        * Referenced by: '<S341>/u3'
                                        */
  1.1428703263986972,                  /* Expression: SM.Lmsatd(1)
                                        * Referenced by: '<S343>/Constant1'
                                        */
  0.31187032639869727,                 /* Expression: SM.Lmsatq(1)
                                        * Referenced by: '<S344>/Constant1'
                                        */

  /*  Expression: SM.One_Llq
   * Referenced by: '<S349>/1//Ll_q'
   */
  { 6.1679023820101309, 9.1590141322923113 },

  /*  Expression: [ 1/SM.Ll 1/SM.Llkq1 1/SM.Llkq2]
   * Referenced by: '<S348>/u2'
   */
  { 6.1679023820101309, 9.1590141322923113, 0.0 },
  0.31187032639869727,                 /* Expression: SM.Lmsatq(1)
                                        * Referenced by: '<S344>/Lmq_sat'
                                        */

  /*  Expression: SM.Phisat
   * Referenced by: '<S344>/Lookup Table'
   */
  { 0.0, 1.0 },

  /*  Expression: [ 0 SM.Phisat(2:end)./SM.Lmsatq(2:end)*SM.Lmq ]
   * Referenced by: '<S344>/Lookup Table'
   */
  { 0.0, 1.0 },
  0.31187032639869727,                 /* Expression: SM.Lmq
                                        * Referenced by: '<S344>/Lmq'
                                        */

  /*  Expression: [ 1/SM.Ll   1/SM.Llfd   1/SM.Llkd ]
   * Referenced by: '<S347>/1//Ll_d'
   */
  { 6.1679023820101309, 6.2538884581967578, 3.9982649324863662 },

  /*  Expression: [1/SM.Ll 1/SM.Llkd 1/SM.Llfd]
   * Referenced by: '<S346>/u1'
   */
  { 6.1679023820101309, 3.9982649324863662, 6.2538884581967578 },
  1.1428703263986972,                  /* Expression: SM.Lmsatd(1)
                                        * Referenced by: '<S343>/Lmd_sat'
                                        */

  /*  Expression: SM.Phisat
   * Referenced by: '<S343>/Lookup Table'
   */
  { 0.0, 1.0 },

  /*  Expression: [ 0 SM.Phisat(2:end)./SM.Lmsatd(2:end)*SM.Lmd ]
   * Referenced by: '<S343>/Lookup Table'
   */
  { 0.0, 1.0 },
  1.1428703263986972,                  /* Expression: SM.Lmd
                                        * Referenced by: '<S343>/Lmd'
                                        */

  /*  Expression: SM.R
   * Referenced by: '<S341>/u1'
   */
  { 0.006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.006, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.00080278661729558, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01457817382743372, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.011161320463928598 },

  /*  Expression: zeros(SM.nState,SM.nState)
   * Referenced by: '<S345>/u1'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /*  Expression: SM.Llqd
   * Referenced by: '<S345>/u5'
   */
  { 0.1621296736013027, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1621296736013027, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.15990051736361466, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.25010848877844089, 0.0, 0.0, 0.0, 0.0, 0.0, 0.10918205666636752 },

  /*  Expression: SM.RLinv
   * Referenced by: '<S331>/Constant6'
   */
  { 0.024691358024691357, 0.0, 0.0, 0.0, -0.034021012089358461, 0.0,
    0.023809523809523791, -0.001790465814242778, -0.020786935596148054, 0.0, 0.0,
    -0.013381880880933086, 0.0032051114510948239, -0.02107672407157371, 0.0, 0.0,
    -0.00855536605978609, -0.0011606468835794865, 0.044812532348662358, 0.0,
    -0.018288702774537289, 0.0, 0.0, 0.0, 0.051707258964228509 },
  0.0,                                 /* Expression: SM.Sat
                                        * Referenced by: '<S331>/Constant2'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S331>/Switch1'
                                        */

  /*  Expression: zeros(SM.nState,SM.nState)
   * Referenced by: '<S339>/u1'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
  -1.0,                                /* Expression: -1
                                        * Referenced by: '<S339>/Gain1'
                                        */
  0.00037699111843077514,              /* Expression: SM.web*Ts/2
                                        * Referenced by: '<S350>/wbase*Ts//2'
                                        */

  /*  Expression: eye(SM.nState,SM.nState)
   * Referenced by: '<S350>/u5'
   */
  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },
  0.00037699111843077514,              /* Expression: SM.web*Ts/2
                                        * Referenced by: '<S350>/wbase*Ts//2 '
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S354>/itail'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S354>/1'
                                        */
  2.0E-6,                              /* Computed Parameter: DiscreteTimeIntegrator_gainva_c
                                        * Referenced by: '<S354>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S354>/Discrete-Time Integrator'
                                        */
  0.9,                                 /* Expression: 0.9
                                        * Referenced by: '<S356>/Constant'
                                        */
  0.9,                                 /* Expression: 0.9
                                        * Referenced by: '<S356>/Saturation1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S356>/Saturation1'
                                        */
  0.1,                                 /* Expression: 0.1
                                        * Referenced by: '<S356>/Saturation2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S356>/Saturation2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S354>/Unit Delay'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S354>/Switch'
                                        */

  /*  Expression: zeros(13,1)
   * Referenced by: '<S359>/SwitchCurrents'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S353>/Unit Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S67>/Unit Delay6'
                                        */
  8.3333333333333331E-5,               /* Expression: sps.Delay
                                        * Referenced by: '<S302>/Constant3'
                                        */
  0.00033333333333333332,              /* Expression: sps.Period
                                        * Referenced by: '<S302>/Constant1'
                                        */
  3000.0,                              /* Expression: sps.Freq
                                        * Referenced by: '<S302>/1\ib1'
                                        */

  /*  Expression: [0 .5 1]
   * Referenced by: '<S302>/Lookup Table'
   */
  { 0.0, 0.5, 1.0 },

  /*  Expression: [0 2 0]
   * Referenced by: '<S302>/Lookup Table'
   */
  { 0.0, 2.0, 0.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S302>/Constant2'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S280>/Gain1'
                                        */
  2.0E-6,                              /* Computed Parameter: Rotorangledthetae_gainval
                                        * Referenced by: '<S329>/Rotor angle dthetae'
                                        */
  0.0,                                 /* Expression: SM.tho
                                        * Referenced by: '<S329>/Rotor angle dthetae'
                                        */
  376.99111843077515,                  /* Expression: SM.web
                                        * Referenced by: '<S329>/web2'
                                        */

  /*  Expression: SM.phiqd0
   * Referenced by: '<S342>/fluxes'
   */
  { -0.0, 1.0, 1.1399113387320832, 1.0, -0.0 },
  0.0,                                 /* Expression: SM.Sat
                                        * Referenced by: '<S331>/Constant1'
                                        */
  0.0,                                 /* Expression: SM.Sat
                                        * Referenced by: '<S331>/Constant3'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S331>/Switch'
                                        */

  /*  Expression: SM.IqdSign
   * Referenced by: '<S331>/change Iq Id  current signs'
   */
  { -1.0, -1.0, 1.0, 1.0, 1.0 },
  12427.649633603134,                  /* Expression: SM.ib
                                        * Referenced by: '<S333>/ib'
                                        */
  97979.58971132712,                   /* Expression: MagnitudeVps
                                        * Referenced by: '<S10>/valp_nom2'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S10>/Switch5'
                                        */
  60.0,                                /* Expression: FreqVps
                                        * Referenced by: '<S10>/valp_nom1'
                                        */
  6.2831853071795862,                  /* Expression: 2*pi
                                        * Referenced by: '<S10>/Gain'
                                        */
  0.0,                                 /* Expression: PhaseVps
                                        * Referenced by: '<S10>/valp_nom'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S10>/Gain3'
                                        */

  /*  Expression: [0  -2*pi/3  2*pi/3]
   * Referenced by: '<S10>/P1'
   */
  { 0.0, -2.0943951023931953, 2.0943951023931953 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S49>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S50>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S51>/do not delete this gain'
                                        */
  0.0021299910806810243,               /* Expression: Kv
                                        * Referenced by: '<S6>/Kv1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S46>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S47>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S48>/do not delete this gain'
                                        */
  6.3380547094514731E-5,               /* Expression: Ki
                                        * Referenced by: '<S6>/Kv'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S67>/Unit Delay3'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S67>/Unit Delay7'
                                        */
  11.111111111111111,                  /* Expression: 5*2/0.9
                                        * Referenced by: '<Root>/MW'
                                        */

  /*  Expression: Gain
   * Referenced by: '<S76>/Gain2'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S329>/nominal speed'
                                        */
  0.0,                                 /* Expression: SM.dwo
                                        * Referenced by: '<S351>/dw_delay'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S351>/F2'
                                        */
  0.0,                                 /* Expression: SM.dwo
                                        * Referenced by: '<S351>/dw_predict'
                                        */
  1.0,                                 /* Expression: SM.Nb
                                        * Referenced by: '<S329>/units'
                                        */
  15.0,                                /* Expression: 15
                                        * Referenced by: '<Root>/Constant'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<Root>/Constant1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S353>/0 4'
                                        */
  5.0E+6,                              /* Expression: 1./Ron
                                        * Referenced by: '<S353>/1//Ron'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S353>/Switch'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S353>/Saturation'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S353>/Saturation'
                                        */
  8.8726172957282831E-5,               /* Expression: 1/(Pnom/(3*sqrt(2)/pi*Vnom_gen))
                                        * Referenced by: '<S8>/-> pu'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S84>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S85>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S86>/do not delete this gain'
                                        */
  8.04657380504274E-5,                 /* Expression: Ki
                                        * Referenced by: '<S64>/Kv'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S87>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S88>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S89>/do not delete this gain'
                                        */
  0.0016777327005364231,               /* Expression: Kv
                                        * Referenced by: '<S64>/Kv1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S99>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S100>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S101>/do not delete this gain'
                                        */
  6.3380547094514731E-5,               /* Expression: Ki
                                        * Referenced by: '<S65>/Kv'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S102>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S103>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S104>/do not delete this gain'
                                        */
  0.0021299910806810243,               /* Expression: Kv
                                        * Referenced by: '<S65>/Kv1'
                                        */
  2.0E-6,                              /* Computed Parameter: DiscreteTimeIntegrator_gainva_a
                                        * Referenced by: '<S126>/Discrete-Time Integrator'
                                        */
  10.0,                                /* Expression: UpperLimit
                                        * Referenced by: '<S126>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: LowerLimit
                                        * Referenced by: '<S126>/Discrete-Time Integrator'
                                        */
  0.0036217260393912988,               /* Expression: sps.D
                                        * Referenced by: '<S221>/D*u(k)'
                                        */
  0.00012631827789682254,              /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S221>/Delay_x1'
                                        */
  7887.839278290794,                   /* Expression: sps.C11
                                        * Referenced by: '<S224>/C11'
                                        */
  0.0,                                 /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S221>/Delay_x2'
                                        */
  0.0072434520787825977,               /* Expression: sps.C12
                                        * Referenced by: '<S224>/C12'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S120>/Switch'
                                        */
  0.0036217260393912988,               /* Expression: sps.D
                                        * Referenced by: '<S197>/D*u(k)'
                                        */
  0.0,                                 /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S197>/Delay_x1'
                                        */
  7887.839278290794,                   /* Expression: sps.C11
                                        * Referenced by: '<S200>/C11'
                                        */
  0.0,                                 /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S197>/Delay_x2'
                                        */
  0.0072434520787825977,               /* Expression: sps.C12
                                        * Referenced by: '<S200>/C12'
                                        */
  2.0E-6,                              /* Computed Parameter: theta_gainval
                                        * Referenced by: '<S329>/theta'
                                        */
  0.0,                                 /* Expression: SM.tho
                                        * Referenced by: '<S329>/theta'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S329>/t'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S8>/deg->rd'
                                        */
  0.0036217260393912988,               /* Expression: sps.D
                                        * Referenced by: '<S213>/D*u(k)'
                                        */
  0.0,                                 /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S213>/Delay_x1'
                                        */
  7887.839278290794,                   /* Expression: sps.C11
                                        * Referenced by: '<S216>/C11'
                                        */
  0.0,                                 /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S213>/Delay_x2'
                                        */
  0.0072434520787825977,               /* Expression: sps.C12
                                        * Referenced by: '<S216>/C12'
                                        */
  1.0,                                 /* Expression: p
                                        * Referenced by: '<S122>/# pairs of poles'
                                        */
  6.2831853071795862,                  /* Expression: 2*pi
                                        * Referenced by: '<S122>/Constant4'
                                        */
  0.0036217260393912988,               /* Expression: sps.D
                                        * Referenced by: '<S201>/D*u(k)'
                                        */
  0.0,                                 /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S201>/Delay_x1'
                                        */
  7887.839278290794,                   /* Expression: sps.C11
                                        * Referenced by: '<S204>/C11'
                                        */
  0.0,                                 /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S201>/Delay_x2'
                                        */
  0.0072434520787825977,               /* Expression: sps.C12
                                        * Referenced by: '<S204>/C12'
                                        */

  /*  Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S275>/Gain3'
   */
  { 1.0, 0.0, 0.5, -0.5, 0.8660254037844386, 0.5, -0.5, -0.8660254037844386, 0.5
  },
  0.66666666666666663,                 /* Expression: 2/3
                                        * Referenced by: '<S275>/Gain1'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S145>/avoid division by 0'
                                        */
  1.0E-6,                              /* Expression: 1e-6
                                        * Referenced by: '<S145>/avoid division by 0'
                                        */
  -1.0,                                /* Expression: -1
                                        * Referenced by: '<S145>/Gain'
                                        */
  5.99964002879784E-5,                 /* Expression: sps.D
                                        * Referenced by: '<S261>/D'
                                        */
  500000.0,                            /* Expression: sps.x0
                                        * Referenced by: '<S261>/Delay_x'
                                        */
  1.9998800071994243E-6,               /* Expression: sps.C
                                        * Referenced by: '<S261>/C'
                                        */
  10.0,                                /* Expression: UpperLimit
                                        * Referenced by: '<S126>/Saturation2'
                                        */
  0.0,                                 /* Expression: LowerLimit
                                        * Referenced by: '<S126>/Saturation2'
                                        */
  1100.0,                              /* Expression: Vdc_nom
                                        * Referenced by: '<S121>/Vdc_ref (V)'
                                        */
  0.0036217260393912988,               /* Expression: sps.D
                                        * Referenced by: '<S217>/D*u(k)'
                                        */
  0.13895010568650479,                 /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S217>/Delay_x1'
                                        */
  7887.839278290794,                   /* Expression: sps.C11
                                        * Referenced by: '<S220>/C11'
                                        */
  1.1368683772161603E-13,              /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S217>/Delay_x2'
                                        */
  0.0072434520787825977,               /* Expression: sps.C12
                                        * Referenced by: '<S220>/C12'
                                        */
  0.00090909090909090909,              /* Expression: 1/Vdc_nom
                                        * Referenced by: '<S121>/->pu'
                                        */
  0.0036217260393912988,               /* Expression: sps.D
                                        * Referenced by: '<S189>/D*u(k)'
                                        */
  0.0,                                 /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S189>/Delay_x1'
                                        */
  7887.839278290794,                   /* Expression: sps.C11
                                        * Referenced by: '<S192>/C11'
                                        */
  0.0,                                 /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S189>/Delay_x2'
                                        */
  0.0072434520787825977,               /* Expression: sps.C12
                                        * Referenced by: '<S192>/C12'
                                        */

  /*  Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S263>/Gain3'
   */
  { 1.0, 0.0, 0.5, -0.5, 0.8660254037844386, 0.5, -0.5, -0.8660254037844386, 0.5
  },
  0.66666666666666663,                 /* Expression: 2/3
                                        * Referenced by: '<S263>/Gain1'
                                        */
  2.0E-6,                              /* Computed Parameter: DiscreteTimeIntegrator_gainva_l
                                        * Referenced by: '<S225>/Discrete-Time Integrator'
                                        */
  -0.78539816339744828,                /* Expression: sps.Phase_Init*pi/180
                                        * Referenced by: '<S225>/Discrete-Time Integrator'
                                        */
  6.2831853071795862,                  /* Expression: 2*pi
                                        * Referenced by: '<S225>/Constant4'
                                        */
  0.15,                                /* Expression: L_RL
                                        * Referenced by: '<S121>/Constant1'
                                        */
  60.0,                                /* Expression: sps.Finit
                                        * Referenced by: '<S225>/Unit Delay'
                                        */
  60.0,                                /* Expression: Fnom
                                        * Referenced by: '<S122>/Fnom'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S121>/Constant2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S131>/IC = i_ic'
                                        */
  2.0E-6,                              /* Computed Parameter: DiscreteTimeIntegrator_gainva_p
                                        * Referenced by: '<S128>/Discrete-Time Integrator'
                                        */
  1.1,                                 /* Expression: UpperLimit
                                        * Referenced by: '<S128>/Discrete-Time Integrator'
                                        */
  -1.1,                                /* Expression: LowerLimit
                                        * Referenced by: '<S128>/Discrete-Time Integrator'
                                        */
  1.1,                                 /* Expression: UpperLimit
                                        * Referenced by: '<S128>/Saturation2'
                                        */
  -1.1,                                /* Expression: LowerLimit
                                        * Referenced by: '<S128>/Saturation2'
                                        */
  1.2100000000000002,                  /* Expression: Imax_grid_conv^2
                                        * Referenced by: '<S121>/Imax^2'
                                        */
  3.3333222222716048E-6,               /* Expression: sps.D
                                        * Referenced by: '<S134>/D'
                                        */
  550000.00000000012,                  /* Expression: sps.x0
                                        * Referenced by: '<S134>/Delay_x'
                                        */
  1.9999933333555554E-6,               /* Expression: sps.C
                                        * Referenced by: '<S134>/C'
                                        */
  -1.0,                                /* Expression: -1
                                        * Referenced by: '<S121>/Gain'
                                        */
  2.0,                                 /* Expression: Ki_volt_reg
                                        * Referenced by: '<S121>/Constant3'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S121>/Unit Delay1'
                                        */
  1.0E-7,                              /* Computed Parameter: DiscreteTimeIntegrator1_gainv_d
                                        * Referenced by: '<S121>/Discrete-Time Integrator1'
                                        */
  1.1,                                 /* Expression: Var_reg_output_hi_limit
                                        * Referenced by: '<S121>/Discrete-Time Integrator1'
                                        */
  0.9,                                 /* Expression: Var_reg_output_low_limit
                                        * Referenced by: '<S121>/Discrete-Time Integrator1'
                                        */
  469.48553403344255,                  /* Expression: Vnom/sqrt(3)*sqrt(2)
                                        * Referenced by: '<S122>/pu->V'
                                        */

  /*  Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S161>/Gain3'
   */
  { 1.0, 0.0, 0.5, -0.5, 0.8660254037844386, 0.5, -0.5, -0.8660254037844386, 0.5
  },
  0.66666666666666663,                 /* Expression: 2/3
                                        * Referenced by: '<S161>/Gain1'
                                        */
  1.0E-6,                              /* Computed Parameter: Integ4_gainval_g
                                        * Referenced by: '<S154>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S154>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S154>/To avoid division by zero'
                                        */
  1.0E-6,                              /* Expression: 1e-6
                                        * Referenced by: '<S154>/To avoid division by zero'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S158>/Unit Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S154>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S154>/Step'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S154>/Switch'
                                        */
  1.0E-6,                              /* Computed Parameter: Integ4_gainval_e
                                        * Referenced by: '<S153>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S153>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S153>/To avoid division by zero'
                                        */
  1.0E-6,                              /* Expression: 1e-6
                                        * Referenced by: '<S153>/To avoid division by zero'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S156>/Unit Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S153>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S153>/Step'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S153>/Switch'
                                        */
  0.0021299910806810243,               /* Expression: 1/(Vnom/sqrt(3)*sqrt(2))
                                        * Referenced by: '<S122>/V->pu'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S131>/IC = 0'
                                        */
  0.003,                               /* Expression: R_RL
                                        * Referenced by: '<S121>/Constant5'
                                        */
  0.0036217260393912988,               /* Expression: sps.D
                                        * Referenced by: '<S193>/D*u(k)'
                                        */
  0.0,                                 /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S193>/Delay_x1'
                                        */
  7887.839278290794,                   /* Expression: sps.C11
                                        * Referenced by: '<S196>/C11'
                                        */
  0.0,                                 /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S193>/Delay_x2'
                                        */
  0.0072434520787825977,               /* Expression: sps.C12
                                        * Referenced by: '<S196>/C12'
                                        */

  /*  Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S269>/Gain3'
   */
  { 1.0, 0.0, 0.5, -0.5, 0.8660254037844386, 0.5, -0.5, -0.8660254037844386, 0.5
  },
  0.66666666666666663,                 /* Expression: 2/3
                                        * Referenced by: '<S269>/Gain1'
                                        */
  2.0E-6,                              /* Computed Parameter: DiscreteTimeIntegrator_gainv_lq
                                        * Referenced by: '<S129>/Discrete-Time Integrator'
                                        */
  1.1,                                 /* Expression: UpperLimit
                                        * Referenced by: '<S129>/Discrete-Time Integrator'
                                        */
  -1.1,                                /* Expression: LowerLimit
                                        * Referenced by: '<S129>/Discrete-Time Integrator'
                                        */
  1.1,                                 /* Expression: UpperLimit
                                        * Referenced by: '<S129>/Saturation2'
                                        */
  -1.1,                                /* Expression: LowerLimit
                                        * Referenced by: '<S129>/Saturation2'
                                        */
  0.003,                               /* Expression: R_RL
                                        * Referenced by: '<S121>/Constant6'
                                        */
  0.15,                                /* Expression: L_RL
                                        * Referenced by: '<S121>/Constant4'
                                        */
  938.971068066885,                    /* Expression: Vnom*2*sqrt(2/3)
                                        * Referenced by: '<S121>/K'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S121>/Avoid division by zero'
                                        */
  1.0E-6,                              /* Expression: 1e-6
                                        * Referenced by: '<S121>/Avoid division by zero'
                                        */
  1.1,                                 /* Expression: Mod_index_max
                                        * Referenced by: '<S121>/0-Mod_index_max'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S121>/0-Mod_index_max'
                                        */
  0.99999333335555551,                 /* Expression: sps.A
                                        * Referenced by: '<S134>/A'
                                        */
  3.333322222271605,                   /* Expression: sps.B
                                        * Referenced by: '<S134>/B'
                                        */
  0.0036217260393912988,               /* Expression: sps.D
                                        * Referenced by: '<S209>/D*u(k)'
                                        */
  0.0,                                 /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S209>/Delay_x1'
                                        */
  7887.839278290794,                   /* Expression: sps.C11
                                        * Referenced by: '<S212>/C11'
                                        */
  0.0,                                 /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S209>/Delay_x2'
                                        */
  0.0072434520787825977,               /* Expression: sps.C12
                                        * Referenced by: '<S212>/C12'
                                        */
  469.48553403344255,                  /* Expression: Vnom/sqrt(3)*sqrt(2)
                                        * Referenced by: '<S122>/Constant1'
                                        */
  750000.0,                            /* Expression: C_var_filter
                                        * Referenced by: '<S122>/C_var_filter'
                                        */
  15777.711708748329,                  /* Expression: Pnom/sqrt(3)/Vnom*sqrt(2)
                                        * Referenced by: '<S122>/pu->A'
                                        */

  /*  Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S175>/Gain3'
   */
  { 1.0, 0.0, 0.5, -0.5, 0.8660254037844386, 0.5, -0.5, -0.8660254037844386, 0.5
  },
  0.66666666666666663,                 /* Expression: 2/3
                                        * Referenced by: '<S175>/Gain1'
                                        */
  1.0E-6,                              /* Computed Parameter: Integ4_gainval_l
                                        * Referenced by: '<S168>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S168>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S168>/To avoid division by zero'
                                        */
  1.0E-6,                              /* Expression: 1e-6
                                        * Referenced by: '<S168>/To avoid division by zero'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S172>/Unit Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S168>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S168>/Step'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S168>/Switch'
                                        */
  1.0E-6,                              /* Computed Parameter: Integ4_gainval_b
                                        * Referenced by: '<S167>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S167>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S167>/To avoid division by zero'
                                        */
  1.0E-6,                              /* Expression: 1e-6
                                        * Referenced by: '<S167>/To avoid division by zero'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S170>/Unit Delay'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S167>/Step'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S167>/Step'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S167>/Switch'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S150>/Rad->Deg.'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S151>/Rad->Deg.'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S142>/Deg->Rad'
                                        */
  1.5,                                 /* Expression: 3/2
                                        * Referenced by: '<S142>/Gain1'
                                        */
  -9.0000000000000012E-8,              /* Expression: -1/Pnom
                                        * Referenced by: '<S122>/var->pu '
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S121>/V0'
                                        */

  /*  Expression: [ 1   0   1; -1/2  sqrt(3)/2   1; -1/2  -sqrt(3)/2  1 ]
   * Referenced by: '<S136>/Gain3'
   */
  { 1.0, -0.5, -0.5, 0.0, 0.8660254037844386, -0.8660254037844386, 1.0, 1.0, 1.0
  },
  15000.0,                             /* Expression: C_var_filter/Q_filter
                                        * Referenced by: '<S122>/C_var_filter//Q'
                                        */
  0.99275654792121737,                 /* Expression: sps.A11
                                        * Referenced by: '<S190>/A11'
                                        */
  1.8299607852399318E-6,               /* Expression: sps.A12
                                        * Referenced by: '<S190>/A12'
                                        */
  -7243.4520787824358,                 /* Expression: sps.A21
                                        * Referenced by: '<S190>/A21'
                                        */
  0.82996078523991357,                 /* Expression: sps.A22
                                        * Referenced by: '<S190>/A22'
                                        */
  9.14980392619977E-7,                 /* Expression: sps.B11
                                        * Referenced by: '<S191>/B11'
                                        */
  0.91498039261995667,                 /* Expression: sps.B21
                                        * Referenced by: '<S191>/B21'
                                        */
  0.99275654792121737,                 /* Expression: sps.A11
                                        * Referenced by: '<S194>/A11'
                                        */
  1.8299607852399318E-6,               /* Expression: sps.A12
                                        * Referenced by: '<S194>/A12'
                                        */
  -7243.4520787824358,                 /* Expression: sps.A21
                                        * Referenced by: '<S194>/A21'
                                        */
  0.82996078523991357,                 /* Expression: sps.A22
                                        * Referenced by: '<S194>/A22'
                                        */
  9.14980392619977E-7,                 /* Expression: sps.B11
                                        * Referenced by: '<S195>/B11'
                                        */
  0.91498039261995667,                 /* Expression: sps.B21
                                        * Referenced by: '<S195>/B21'
                                        */
  0.99275654792121737,                 /* Expression: sps.A11
                                        * Referenced by: '<S198>/A11'
                                        */
  1.8299607852399318E-6,               /* Expression: sps.A12
                                        * Referenced by: '<S198>/A12'
                                        */
  -7243.4520787824358,                 /* Expression: sps.A21
                                        * Referenced by: '<S198>/A21'
                                        */
  0.82996078523991357,                 /* Expression: sps.A22
                                        * Referenced by: '<S198>/A22'
                                        */
  9.14980392619977E-7,                 /* Expression: sps.B11
                                        * Referenced by: '<S199>/B11'
                                        */
  0.91498039261995667,                 /* Expression: sps.B21
                                        * Referenced by: '<S199>/B21'
                                        */
  0.99275654792121737,                 /* Expression: sps.A11
                                        * Referenced by: '<S202>/A11'
                                        */
  1.8299607852399318E-6,               /* Expression: sps.A12
                                        * Referenced by: '<S202>/A12'
                                        */
  -7243.4520787824358,                 /* Expression: sps.A21
                                        * Referenced by: '<S202>/A21'
                                        */
  0.82996078523991357,                 /* Expression: sps.A22
                                        * Referenced by: '<S202>/A22'
                                        */
  9.14980392619977E-7,                 /* Expression: sps.B11
                                        * Referenced by: '<S203>/B11'
                                        */
  0.91498039261995667,                 /* Expression: sps.B21
                                        * Referenced by: '<S203>/B21'
                                        */
  0.0,                                 /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S205>/Delay_x1'
                                        */
  0.99275654792121737,                 /* Expression: sps.A11
                                        * Referenced by: '<S206>/A11'
                                        */
  0.0,                                 /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S205>/Delay_x2'
                                        */
  1.8299607852399318E-6,               /* Expression: sps.A12
                                        * Referenced by: '<S206>/A12'
                                        */
  -7243.4520787824358,                 /* Expression: sps.A21
                                        * Referenced by: '<S206>/A21'
                                        */
  0.82996078523991357,                 /* Expression: sps.A22
                                        * Referenced by: '<S206>/A22'
                                        */
  9.14980392619977E-7,                 /* Expression: sps.B11
                                        * Referenced by: '<S207>/B11'
                                        */
  0.91498039261995667,                 /* Expression: sps.B21
                                        * Referenced by: '<S207>/B21'
                                        */
  0.0036217260393912988,               /* Expression: sps.D
                                        * Referenced by: '<S205>/D*u(k)'
                                        */
  7887.839278290794,                   /* Expression: sps.C11
                                        * Referenced by: '<S208>/C11'
                                        */
  0.0072434520787825977,               /* Expression: sps.C12
                                        * Referenced by: '<S208>/C12'
                                        */
  0.99275654792121737,                 /* Expression: sps.A11
                                        * Referenced by: '<S218>/A11'
                                        */
  1.8299607852399318E-6,               /* Expression: sps.A12
                                        * Referenced by: '<S218>/A12'
                                        */
  -7243.4520787824358,                 /* Expression: sps.A21
                                        * Referenced by: '<S218>/A21'
                                        */
  0.82996078523991357,                 /* Expression: sps.A22
                                        * Referenced by: '<S218>/A22'
                                        */
  9.14980392619977E-7,                 /* Expression: sps.B11
                                        * Referenced by: '<S219>/B11'
                                        */
  0.91498039261995667,                 /* Expression: sps.B21
                                        * Referenced by: '<S219>/B21'
                                        */
  0.99275654792121737,                 /* Expression: sps.A11
                                        * Referenced by: '<S210>/A11'
                                        */
  1.8299607852399318E-6,               /* Expression: sps.A12
                                        * Referenced by: '<S210>/A12'
                                        */
  -7243.4520787824358,                 /* Expression: sps.A21
                                        * Referenced by: '<S210>/A21'
                                        */
  0.82996078523991357,                 /* Expression: sps.A22
                                        * Referenced by: '<S210>/A22'
                                        */
  9.14980392619977E-7,                 /* Expression: sps.B11
                                        * Referenced by: '<S211>/B11'
                                        */
  0.91498039261995667,                 /* Expression: sps.B21
                                        * Referenced by: '<S211>/B21'
                                        */
  0.99275654792121737,                 /* Expression: sps.A11
                                        * Referenced by: '<S222>/A11'
                                        */
  1.8299607852399318E-6,               /* Expression: sps.A12
                                        * Referenced by: '<S222>/A12'
                                        */
  -7243.4520787824358,                 /* Expression: sps.A21
                                        * Referenced by: '<S222>/A21'
                                        */
  0.82996078523991357,                 /* Expression: sps.A22
                                        * Referenced by: '<S222>/A22'
                                        */
  9.14980392619977E-7,                 /* Expression: sps.B11
                                        * Referenced by: '<S223>/B11'
                                        */
  0.91498039261995667,                 /* Expression: sps.B21
                                        * Referenced by: '<S223>/B21'
                                        */
  0.99275654792121737,                 /* Expression: sps.A11
                                        * Referenced by: '<S214>/A11'
                                        */
  1.8299607852399318E-6,               /* Expression: sps.A12
                                        * Referenced by: '<S214>/A12'
                                        */
  -7243.4520787824358,                 /* Expression: sps.A21
                                        * Referenced by: '<S214>/A21'
                                        */
  0.82996078523991357,                 /* Expression: sps.A22
                                        * Referenced by: '<S214>/A22'
                                        */
  9.14980392619977E-7,                 /* Expression: sps.B11
                                        * Referenced by: '<S215>/B11'
                                        */
  0.91498039261995667,                 /* Expression: sps.B21
                                        * Referenced by: '<S215>/B21'
                                        */
  1.0,                                 /* Expression: sps.AGC
                                        * Referenced by: '<S225>/Constant1'
                                        */

  /*  Expression: [ 1   -1/2   -1/2; 0   sqrt(3)/2   -sqrt(3)/2; 1/2  1/2  1/2 ]
   * Referenced by: '<S255>/Gain3'
   */
  { 1.0, 0.0, 0.5, -0.5, 0.8660254037844386, 0.5, -0.5, -0.8660254037844386, 0.5
  },
  0.66666666666666663,                 /* Expression: 2/3
                                        * Referenced by: '<S255>/Gain1'
                                        */
  1.0E-6,                              /* Computed Parameter: Integ4_gainval_d
                                        * Referenced by: '<S251>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S251>/Integ4'
                                        */
  1.0E+6,                              /* Expression: 1e6
                                        * Referenced by: '<S251>/To avoid division  by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S251>/To avoid division  by zero'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_c
   * Referenced by: '<S253>/S-Function'
   */
  { 1.0, 1.0 },
  0.022224222222222222,                /* Expression: MaxDelay
                                        * Referenced by: '<S253>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_e
   * Referenced by: '<S253>/S-Function'
   */
  { 1.0, 1.0 },
  2.0E-6,                              /* Expression: Ts
                                        * Referenced by: '<S253>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_m
   * Referenced by: '<S253>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S253>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_n
   * Referenced by: '<S253>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S253>/S-Function'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S252>/Unit Delay'
                                        */
  0.016666666666666666,                /* Expression: 1/sps.Finit
                                        * Referenced by: '<S251>/Constant'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S251>/Unit Delay1'
                                        */

  /*  Expression: [ TcD  Ts-TcD ]
   * Referenced by: '<S227>/Discrete Derivative '
   */
  { 0.0001, -9.800000000000001E-5 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S227>/Discrete Derivative '
                                        */
  0.0028,                              /* Computed Parameter: DiscreteTimeIntegrator_gainva_o
                                        * Referenced by: '<S227>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(1)
                                        * Referenced by: '<S227>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(2)
                                        * Referenced by: '<S227>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Par_Limits(1)
                                        * Referenced by: '<S227>/Saturation1'
                                        */
  0.0,                                 /* Expression: Par_Limits(2)
                                        * Referenced by: '<S227>/Saturation1'
                                        */
  0.15915494309189535,                 /* Expression: 1/2/pi
                                        * Referenced by: '<S225>/Gain10'
                                        */
  2.4E-5,                              /* Computed Parameter: RateLimiter_RisingLim
                                        * Referenced by: '<S225>/Rate Limiter'
                                        */
  -2.4E-5,                             /* Computed Parameter: RateLimiter_FallingLim
                                        * Referenced by: '<S225>/Rate Limiter'
                                        */
  60.0,                                /* Expression: sps.Finit
                                        * Referenced by: '<S225>/Rate Limiter'
                                        */
  1215.8541837080536,                  /* Expression: sps.x0(1,:)
                                        * Referenced by: '<S247>/Delay_x1'
                                        */
  0.99999995066293668,                 /* Expression: sps.A11
                                        * Referenced by: '<S248>/A11'
                                        */
  -3.5527136788005009E-15,             /* Expression: sps.x0(2,:)
                                        * Referenced by: '<S247>/Delay_x2'
                                        */
  1.9995558281133505E-6,               /* Expression: sps.A12
                                        * Referenced by: '<S248>/A12'
                                        */
  -0.049337063314990645,               /* Expression: sps.A21
                                        * Referenced by: '<S248>/A21'
                                        */
  0.99955582811335064,                 /* Expression: sps.A22
                                        * Referenced by: '<S248>/A22'
                                        */
  9.9977791405667527E-7,               /* Expression: sps.B11
                                        * Referenced by: '<S249>/B11'
                                        */
  0.99977791405667527,                 /* Expression: sps.B21
                                        * Referenced by: '<S249>/B21'
                                        */
  2.4668531657495326E-8,               /* Expression: sps.D
                                        * Referenced by: '<S247>/D*u(k)'
                                        */
  0.04934802159984595,                 /* Expression: sps.C11
                                        * Referenced by: '<S250>/C11'
                                        */
  4.9337063314990652E-8,               /* Expression: sps.C12
                                        * Referenced by: '<S250>/C12'
                                        */
  0.99988000719942416,                 /* Expression: sps.A
                                        * Referenced by: '<S261>/A'
                                        */
  59.9964002879784,                    /* Expression: sps.B
                                        * Referenced by: '<S261>/B'
                                        */
  -9.0000000000000012E-8,              /* Expression: -1/Pnom
                                        * Referenced by: '<S122>/W->pu'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S67>/Unit Delay2'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S125>/0-inf'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S125>/0-inf'
                                        */
  1.1111111111111109,                  /* Expression: Pnom/Pmec
                                        * Referenced by: '<S125>/pu_elec->pu_mec'
                                        */
  0.75,                                /* Expression: 0.75
                                        * Referenced by: '<S125>/Switch'
                                        */
  1.9999996000001067E-7,               /* Expression: sps.D
                                        * Referenced by: '<S310>/D'
                                        */
  500000.0,                            /* Expression: sps.x0
                                        * Referenced by: '<S310>/Delay_x'
                                        */
  1.99999960000008E-6,                 /* Expression: sps.C
                                        * Referenced by: '<S310>/C'
                                        */
  15.0,                                /* Expression: Kp_pitch
                                        * Referenced by: '<S125>/Gain2'
                                        */
  27.0,                                /* Expression: pitch_max
                                        * Referenced by: '<S125>/0-pitch_max'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S125>/0-pitch_max'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S125>/Constant2'
                                        */
  2.0E-6,                              /* Computed Parameter: DiscreteTimeIntegrator_gainva_f
                                        * Referenced by: '<S305>/Discrete-Time Integrator'
                                        */
  27.0,                                /* Expression: UpperLimit
                                        * Referenced by: '<S305>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: LowerLimit
                                        * Referenced by: '<S305>/Discrete-Time Integrator'
                                        */
  27.0,                                /* Expression: UpperLimit
                                        * Referenced by: '<S305>/Saturation2'
                                        */
  0.0,                                 /* Expression: LowerLimit
                                        * Referenced by: '<S305>/Saturation2'
                                        */
  10.0,                                /* Expression: pitch_rate
                                        * Referenced by: '<S125>/Rate Limiter   1'
                                        */
  -10.0,                               /* Expression: -pitch_rate
                                        * Referenced by: '<S125>/Rate Limiter   1'
                                        */
  27.0,                                /* Expression: pitch_max
                                        * Referenced by: '<S125>/0-pitch_max '
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S125>/0-pitch_max '
                                        */
  2.0E-6,                              /* Computed Parameter: DiscreteTimeIntegrator_gainva_e
                                        * Referenced by: '<S304>/Discrete-Time Integrator'
                                        */
  1.1,                                 /* Expression: UpperLimit
                                        * Referenced by: '<S304>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: LowerLimit
                                        * Referenced by: '<S304>/Discrete-Time Integrator'
                                        */
  1.1,                                 /* Expression: UpperLimit
                                        * Referenced by: '<S304>/Saturation2'
                                        */
  0.0,                                 /* Expression: LowerLimit
                                        * Referenced by: '<S304>/Saturation2'
                                        */
  2.0E-6,                              /* Computed Parameter: DiscreteTimeIntegrator_gainv_cb
                                        * Referenced by: '<S306>/Discrete-Time Integrator'
                                        */
  1.0,                                 /* Expression: UpperLimit
                                        * Referenced by: '<S306>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: LowerLimit
                                        * Referenced by: '<S306>/Discrete-Time Integrator'
                                        */
  1.0,                                 /* Expression: UpperLimit
                                        * Referenced by: '<S306>/Saturation2'
                                        */
  0.0,                                 /* Expression: LowerLimit
                                        * Referenced by: '<S306>/Saturation2'
                                        */
  0.0,                                 /* Expression: sps.x0
                                        * Referenced by: '<S309>/Delay_x'
                                        */
  0.99980001999733381,                 /* Expression: sps.A
                                        * Referenced by: '<S309>/A'
                                        */
  99.990001333166674,                  /* Expression: sps.B
                                        * Referenced by: '<S309>/B'
                                        */
  1.9998000199973339E-6,               /* Expression: sps.C
                                        * Referenced by: '<S309>/C'
                                        */
  9.9990001333166674E-5,               /* Expression: sps.D
                                        * Referenced by: '<S309>/D'
                                        */
  0.99999960000008,                    /* Expression: sps.A
                                        * Referenced by: '<S310>/A'
                                        */
  0.19999996000001069,                 /* Expression: sps.B
                                        * Referenced by: '<S310>/B'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S67>/Unit Delay1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S67>/Unit Delay4'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S80>/Avoid div. by zero'
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S80>/Avoid div. by zero'
                                        */
  2.0E-6,                              /* Computed Parameter: DiscreteTimeIntegrator_gainva_m
                                        * Referenced by: '<S69>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: inf
                                        * Referenced by: '<S80>/Avoid div. by zero '
                                        */
  2.2204460492503131E-16,              /* Expression: eps
                                        * Referenced by: '<S80>/Avoid div. by zero '
                                        */
  2.0E-6,                              /* Computed Parameter: DiscreteTimeIntegrator1_gainv_i
                                        * Referenced by: '<S69>/Discrete-Time Integrator1'
                                        */
  0.0,                                 /* Expression: Display
                                        * Referenced by: '<S76>/Mode'
                                        */
  0.90000000000000013,                 /* Expression: Pmec/Pnom
                                        * Referenced by: '<S8>/Power base for the Generator'
                                        */
  0.000702430187180766,                /* Expression: SM.N
                                        * Referenced by: '<S327>/N'
                                        */

  /*  Expression: [1 -1]
   * Referenced by: '<S338>/1-1'
   */
  { 1.0, -1.0 },
  0.0016777327005364235,               /* Expression: 1/SM.Vb
                                        * Referenced by: '<S332>/1_Vb'
                                        */

  /*  Expression: zeros(1, SM.nState-3)
   * Referenced by: '<S327>/[ Vkd =0 Vkq1=0  Vkq2=0 ]'
   */
  { 0.0, 0.0 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S342>/voltages'
                                        */
  2.0E-6,                              /* Expression: Ts
                                        * Referenced by: '<S342>/IC'
                                        */
  1.0,                                 /* Expression: 1/SM.Pb
                                        * Referenced by: '<S329>/1_Pb'
                                        */
  0.01,                                /* Expression: SM.F
                                        * Referenced by: '<S329>/F'
                                        */
  0.80645161290322587,                 /* Expression: 1/(2*SM.H)
                                        * Referenced by: '<S329>/1 ----- 2H'
                                        */
  1.0E-6,                              /* Computed Parameter: Rotorspeeddeviationdw_gainval
                                        * Referenced by: '<S329>/Rotor speed deviation (dw)'
                                        */
  0.0,                                 /* Expression: SM.dwo
                                        * Referenced by: '<S329>/Rotor speed deviation (dw)'
                                        */
  376.99111843077515,                  /* Expression: SM.web
                                        * Referenced by: '<S329>/we base'
                                        */
  376.99111843077515,                  /* Expression: SM.web
                                        * Referenced by: '<S329>/web3'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S8>/In_ac_switch '
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S8>/In_dc_switch 1'
                                        */

  /*  Expression: zeros(1,6)
   * Referenced by: '<S78>/g'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S34>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S35>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S36>/do not delete this gain'
                                        */
  4.8989794855663556E-5,               /* Expression: Kv
                                        * Referenced by: '<S5>/Kv1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S31>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S32>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S33>/do not delete this gain'
                                        */
  0.0030618621784789723,               /* Expression: Ki
                                        * Referenced by: '<S5>/Kv'
                                        */
  0,                                   /* Computed Parameter: selector_Y0
                                        * Referenced by: '<S12>/selector'
                                        */
  0,                                   /* Expression: SM.nState==6
                                        * Referenced by: '<S341>/Constant1'
                                        */
  0,                                   /* Expression: SM.nState==6
                                        * Referenced by: '<S341>/Constant2'
                                        */
  0,                                   /* Expression: Tf_sps>0 | Tt_sps>0
                                        * Referenced by: '<S353>/2'
                                        */

  /* Start of '<S274>/Subsystem1' */
  {
    /*  Expression: [0,0]
     * Referenced by: '<S279>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S274>/Subsystem1' */
  ,

  /* Start of '<S274>/Subsystem - pi//2 delay' */
  {
    /*  Expression: [0,0]
     * Referenced by: '<S278>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S274>/Subsystem - pi//2 delay' */
  ,

  /* Start of '<S268>/Subsystem1' */
  {
    /*  Expression: [0,0]
     * Referenced by: '<S273>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S268>/Subsystem1' */
  ,

  /* Start of '<S268>/Subsystem - pi//2 delay' */
  {
    /*  Expression: [0,0]
     * Referenced by: '<S272>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S268>/Subsystem - pi//2 delay' */
  ,

  /* Start of '<S262>/Subsystem1' */
  {
    /*  Expression: [0,0]
     * Referenced by: '<S267>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S262>/Subsystem1' */
  ,

  /* Start of '<S262>/Subsystem - pi//2 delay' */
  {
    /*  Expression: [0,0]
     * Referenced by: '<S266>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S262>/Subsystem - pi//2 delay' */
  ,

  /* Start of '<S254>/Subsystem1' */
  {
    /*  Expression: [0,0]
     * Referenced by: '<S259>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S254>/Subsystem1' */
  ,

  /* Start of '<S254>/Subsystem - pi//2 delay' */
  {
    /*  Expression: [0,0]
     * Referenced by: '<S258>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S254>/Subsystem - pi//2 delay' */
  ,

  /* Start of '<S241>/Subsystem1' */
  {
    /*  Expression: [0,0]
     * Referenced by: '<S246>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S241>/Subsystem1' */
  ,

  /* Start of '<S241>/Subsystem - pi//2 delay' */
  {
    /*  Expression: [0,0]
     * Referenced by: '<S245>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S241>/Subsystem - pi//2 delay' */
  ,

  /* Start of '<S174>/Subsystem1' */
  {
    /*  Expression: [0,0]
     * Referenced by: '<S179>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S174>/Subsystem1' */
  ,

  /* Start of '<S174>/Subsystem - pi//2 delay' */
  {
    /*  Expression: [0,0]
     * Referenced by: '<S178>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S174>/Subsystem - pi//2 delay' */
  ,

  /* Start of '<S160>/Subsystem1' */
  {
    /*  Expression: [0,0]
     * Referenced by: '<S165>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S160>/Subsystem1' */
  ,

  /* Start of '<S160>/Subsystem - pi//2 delay' */
  {
    /*  Expression: [0,0]
     * Referenced by: '<S164>/dq'
     */
    { 0.0, 0.0 }
  }
  /* End of '<S160>/Subsystem - pi//2 delay' */
};

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
