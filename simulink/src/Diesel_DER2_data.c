/*
 * File: Diesel_DER2_data.c
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

/* Block parameters (auto storage) */
P_Diesel_DER2_T Diesel_DER2_P = {
  50.0,                                /* Mask Parameter: Fourier_A_Freq
                                        * Referenced by:
                                        *   '<S87>/cos(wt)'
                                        *   '<S87>/sin(wt)'
                                        */
  50.0,                                /* Mask Parameter: Fourier_B_Freq
                                        * Referenced by:
                                        *   '<S88>/cos(wt)'
                                        *   '<S88>/sin(wt)'
                                        */
  50.0,                                /* Mask Parameter: Fourier_C_Freq
                                        * Referenced by:
                                        *   '<S89>/cos(wt)'
                                        *   '<S89>/sin(wt)'
                                        */
  1.0090570934256056,                  /* Mask Parameter: DiscretePIController_Init
                                        * Referenced by: '<S85>/Discrete-Time Integrator'
                                        */
  40.0,                                /* Mask Parameter: DieselEngineSpeedRegulator_K
                                        * Referenced by: '<S2>/Gain K'
                                        */
  8.0,                                 /* Mask Parameter: DiscretePIController_Ki
                                        * Referenced by: '<S85>/Kp5'
                                        */
  10.0,                                /* Mask Parameter: DiscretePIController_Kp
                                        * Referenced by: '<S85>/Kp4'
                                        */
  0.25,                                /* Mask Parameter: DieselEngineSpeedRegulator_Pm0
                                        * Referenced by:
                                        *   '<S2>/Discrete-Time Integrator'
                                        *   '<S2>/Delay Td'
                                        */
  1.0,                                 /* Mask Parameter: ThreePhaseBreaker_SwitchA
                                        * Referenced by: '<S8>/Constant1'
                                        */
  1.0,                                 /* Mask Parameter: ThreePhaseBreaker_SwitchB
                                        * Referenced by: '<S8>/Constant2'
                                        */
  1.0,                                 /* Mask Parameter: ThreePhaseBreaker_SwitchC
                                        * Referenced by: '<S8>/Constant3'
                                        */
  0.024,                               /* Mask Parameter: DieselEngineSpeedRegulator_Td
                                        * Referenced by: '<S2>/Delay Td'
                                        */

  /*  Expression: SM.Linv
   * Referenced by: '<S56>/Constant4'
   */
  { 4.8319110392301683, 0.0, 0.0, 0.0, -4.03960252471989, 0.0,
    7.1430894071971167, -2.6978145895866126, -4.1640661513663018, 0.0, 0.0,
    -2.6978145895866126, 4.4755318245091189, -1.6652586101797939, 0.0, 0.0,
    -4.1640661513663018, -1.6652586101797939, 6.0029042041541549, 0.0,
    -4.03960252471989, 0.0, 0.0, 0.0, 4.3929799281866906 },
  0.82304818785531353,                 /* Expression: SM.Lmq
                                        * Referenced by: '<S66>/u3'
                                        */
  1.7272083710354986,                  /* Expression: SM.Lmsatd(1)
                                        * Referenced by: '<S68>/Constant1'
                                        */

  /*  Expression: [ 1/SM.Ll   1/SM.Llfd   1/SM.Llkd ]
   * Referenced by: '<S72>/1//Ll_d'
   */
  { 13.889115727965949, 5.554419337500164, 8.5732242842244144 },

  /*  Expression: [1/SM.Ll 1/SM.Llkd 1/SM.Llfd]
   * Referenced by: '<S71>/u1'
   */
  { 13.889115727965949, 8.5732242842244144, 5.554419337500164 },
  1.7272083710354986,                  /* Expression: SM.Lmsatd(1)
                                        * Referenced by: '<S68>/Lmd_sat'
                                        */

  /*  Expression: SM.Phisat
   * Referenced by: '<S68>/Lookup Table'
   */
  { 0.0, 1.0 },

  /*  Expression: [ 0 SM.Phisat(2:end)./SM.Lmsatd(2:end)*SM.Lmd ]
   * Referenced by: '<S68>/Lookup Table'
   */
  { 0.0, 1.0 },
  1.7272083710354986,                  /* Expression: SM.Lmd
                                        * Referenced by: '<S68>/Lmd'
                                        */

  /*  Expression: SM.R
   * Referenced by: '<S66>/u1'
   */
  { 0.082012500000000016, 0.0, 0.0, 0.0, 0.0, 0.0, 0.082012500000000016, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.061155, 0.0, 0.0, 0.0, 0.0, 0.0, 0.15906375, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.24158250000000003 },

  /*  Expression: zeros(SM.nState,SM.nState)
   * Referenced by: '<S70>/u1'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /*  Expression: SM.Llqd
   * Referenced by: '<S70>/u5'
   */
  { 0.071998824085430041, 0.0, 0.0, 0.0, 0.0, 0.0, 0.071998824085430041, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.18003682099559709, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.11664223013972695, 0.0, 0.0, 0.0, 0.0, 0.0, 0.16142877500930303 },

  /*  Expression: SM.RLinv
   * Referenced by: '<S56>/Constant6'
   */
  { 0.39627710410486428, 0.0, 0.0, 0.0, -0.975897276928143, 0.0,
    0.58582262000775365, -0.16498485122616929, -0.66235197728439155, 0.0, 0.0,
    -0.2212545190284721, 0.27370114872785517, -0.2648822792549862, 0.0, 0.0,
    -0.34150547523892888, -0.1018388903055453, 0.95484445360352543, 0.0,
    -0.33129790205859005, 0.0, 0.0, 0.0, 1.0612670735011613 },
  0.0,                                 /* Expression: SM.Sat
                                        * Referenced by: '<S56>/Constant2'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S56>/Switch1'
                                        */

  /*  Expression: zeros(SM.nState,SM.nState)
   * Referenced by: '<S64>/u1'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
  -1.0,                                /* Expression: -1
                                        * Referenced by: '<S64>/Gain1'
                                        */
  0.0078539816339744835,               /* Expression: SM.web*Ts/2
                                        * Referenced by: '<S75>/wbase*Ts//2'
                                        */

  /*  Expression: eye(SM.nState,SM.nState)
   * Referenced by: '<S75>/u5'
   */
  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },
  0.0078539816339744835,               /* Expression: SM.web*Ts/2
                                        * Referenced by: '<S75>/wbase*Ts//2 '
                                        */
  0.0,                                 /* Expression: sps.K1
                                        * Referenced by: '<S95>/Gain'
                                        */
  0.0,                                 /* Expression: sps.K2
                                        * Referenced by: '<S95>/Gain1'
                                        */
  0.0,                                 /* Expression: sps.K1
                                        * Referenced by: '<S97>/Gain'
                                        */
  0.0,                                 /* Expression: sps.K2
                                        * Referenced by: '<S97>/Gain1'
                                        */
  0.0,                                 /* Expression: sps.K1
                                        * Referenced by: '<S101>/Gain'
                                        */
  0.0,                                 /* Expression: sps.K2
                                        * Referenced by: '<S101>/Gain1'
                                        */
  0.0,                                 /* Expression: sps.K1
                                        * Referenced by: '<S103>/Gain'
                                        */
  0.0,                                 /* Expression: sps.K2
                                        * Referenced by: '<S103>/Gain1'
                                        */
  0.0,                                 /* Expression: sps.K1
                                        * Referenced by: '<S107>/Gain'
                                        */
  0.0,                                 /* Expression: sps.K2
                                        * Referenced by: '<S107>/Gain1'
                                        */
  0.0,                                 /* Expression: sps.K1
                                        * Referenced by: '<S109>/Gain'
                                        */
  0.0,                                 /* Expression: sps.K2
                                        * Referenced by: '<S109>/Gain1'
                                        */
  0.33333333333333331,                 /* Expression: 1/3
                                        * Referenced by: '<S92>/Gain3'
                                        */
  0.0,                                 /* Expression: sps.K1
                                        * Referenced by: '<S111>/Gain'
                                        */
  0.0,                                 /* Expression: sps.K2
                                        * Referenced by: '<S111>/Gain1'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S120>/Constant2'
                                        */

  /*  Expression: [SM.Laq  SM.Lad]
   * Referenced by: '<S121>/Laqd_nosat'
   */
  { 0.039124613682218694, 0.041871388702726067 },

  /*  Expression: [SM.Lmq SM.Lmd]
   * Referenced by: '<S121>/Lmqd_nosat'
   */
  { 1.5099279691315945, 2.0600993825915066 },

  /*  Expression: SM.Linv
   * Referenced by: '<S121>/Constant4'
   */
  { 4.3480276717263306, 0.0, 0.0, 0.0, -3.8296986296302156, 0.0,
    3.2489041286284333, -2.4507567435320281, -0.39157216704533232, 0.0, 0.0,
    -2.4507567435320281, 2.5675028794883334, -0.057275809335712975, 0.0, 0.0,
    -0.39157216704533232, -0.057275809335712975, 0.45834990868468062, 0.0,
    -3.8296986296302156, 0.0, 0.0, 0.0, 3.9564920328946958 },
  0.039124613682218694,                /* Expression: SM.Laq
                                        * Referenced by: '<S131>/u2'
                                        */
  1.5099279691315945,                  /* Expression: SM.Lmq
                                        * Referenced by: '<S131>/u3'
                                        */
  2.0600993825915066,                  /* Expression: SM.Lmsatd(1)
                                        * Referenced by: '<S133>/Constant1'
                                        */

  /*  Expression: [ 1/SM.Ll   1/SM.Llfd   1/SM.Llkd ]
   * Referenced by: '<S137>/1//Ll_d'
   */
  { 20.003763467952279, 2.925978501057481, 0.46750120974303755 },

  /*  Expression: [1/SM.Ll 1/SM.Llkd 1/SM.Llfd]
   * Referenced by: '<S136>/u1'
   */
  { 20.003763467952279, 0.46750120974303755, 2.925978501057481 },
  2.0600993825915066,                  /* Expression: SM.Lmsatd(1)
                                        * Referenced by: '<S133>/Lmd_sat'
                                        */

  /*  Expression: SM.Phisat
   * Referenced by: '<S133>/Lookup Table'
   */
  { 0.0, 1.0 },

  /*  Expression: [ 0 SM.Phisat(2:end)./SM.Lmsatd(2:end)*SM.Lmd ]
   * Referenced by: '<S133>/Lookup Table'
   */
  { 0.0, 1.0 },
  2.0600993825915066,                  /* Expression: SM.Lmd
                                        * Referenced by: '<S133>/Lmd'
                                        */

  /*  Expression: SM.R
   * Referenced by: '<S131>/u1'
   */
  { 0.0095, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0095, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00197,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.20124999999999998, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.026812500000000003 },

  /*  Expression: zeros(SM.nState,SM.nState)
   * Referenced by: '<S135>/u1'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  /*  Expression: SM.Llqd
   * Referenced by: '<S135>/u5'
   */
  { 0.049990593100247588, 0.0, 0.0, 0.0, 0.0, 0.0, 0.049990593100247588, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.3417660108023996, 0.0, 0.0, 0.0, 0.0, 0.0,
    2.1390318980129504, 0.0, 0.0, 0.0, 0.0, 0.0, 0.20436060211601606 },

  /*  Expression: SM.RLinv
   * Referenced by: '<S121>/Constant6'
   */
  { 0.041306262881400137, 0.0, 0.0, 0.0, -0.10268379450696016, 0.0,
    0.030864589221970114, -0.0048279907847580952, -0.078803898617873125, 0.0,
    0.0, -0.023282189063554268, 0.0050579806725920169, -0.011526756628812236,
    0.0, 0.0, -0.003719935586930657, -0.00011283334439135455,
    0.092242919122791964, 0.0, -0.036382136981487047, 0.0, 0.0, 0.0,
    0.10608344263198904 },
  0.0,                                 /* Expression: SM.Sat
                                        * Referenced by: '<S121>/Constant2'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S121>/Switch1'
                                        */

  /*  Expression: zeros(SM.nState,SM.nState)
   * Referenced by: '<S129>/u1'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
  -1.0,                                /* Expression: -1
                                        * Referenced by: '<S129>/Gain1'
                                        */
  0.0078539816339744835,               /* Expression: SM.web*Ts/2
                                        * Referenced by: '<S140>/wbase*Ts//2'
                                        */

  /*  Expression: eye(SM.nState,SM.nState)
   * Referenced by: '<S140>/u5'
   */
  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },
  0.0078539816339744835,               /* Expression: SM.web*Ts/2
                                        * Referenced by: '<S140>/wbase*Ts//2 '
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S8>/com'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S8>/Switch'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S8>/Switch1'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S8>/Switch2'
                                        */
  5.0E-5,                              /* Computed Parameter: Rotorangledthetae_gainval
                                        * Referenced by: '<S119>/Rotor angle dthetae'
                                        */
  -1.1997602717719262,                 /* Expression: SM.tho
                                        * Referenced by: '<S119>/Rotor angle dthetae'
                                        */
  314.15926535897933,                  /* Expression: SM.web
                                        * Referenced by: '<S119>/web2'
                                        */

  /*  Expression: SM.phiqd0
   * Referenced by: '<S132>/fluxes'
   */
  { -0.36344212219181876, 0.9341684712488878, 1.125407609158289,
    0.9386998800853984, -0.35179491977633692 },
  0.0,                                 /* Expression: SM.Sat
                                        * Referenced by: '<S121>/Constant1'
                                        */
  0.0,                                 /* Expression: SM.Sat
                                        * Referenced by: '<S121>/Constant3'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S121>/Switch'
                                        */

  /*  Expression: SM.IqdSign
   * Referenced by: '<S121>/change Iq Id  current signs'
   */
  { -1.0, -1.0, 1.0, 1.0, 1.0 },

  /*  Expression: SM.One_Llq
   * Referenced by: '<S130>/1//Ll_q'
   */
  { 20.003763467952279, 4.8933110866070813 },
  0.0,                                 /* Expression: SM.Sat
                                        * Referenced by: '<S121>/Constant8'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S121>/Switch3'
                                        */

  /*  Expression: [ 1/SM.Ll   1/SM.Llfd   1/SM.Llkd ]
   * Referenced by: '<S130>/1//Ll_d'
   */
  { 20.003763467952279, 2.925978501057481, 0.46750120974303755 },

  /*  Expression: zeros(9,1)
   * Referenced by: '<S158>/SwitchCurrents'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
  4082.48290463863,                    /* Expression: SM.ib
                                        * Referenced by: '<S123>/ib'
                                        */
  5.0E-5,                              /* Computed Parameter: Rotorangledtheta_gainval
                                        * Referenced by: '<S54>/Rotor angle dtheta'
                                        */
  -1.5155008054332124,                 /* Expression: SM.tho
                                        * Referenced by: '<S54>/Rotor angle dtheta'
                                        */
  314.15926535897933,                  /* Expression: SM.web
                                        * Referenced by: '<S54>/web2'
                                        */

  /*  Expression: SM.phiqd0
   * Referenced by: '<S67>/fluxes'
   */
  { -0.055502684037450858, 1.0035582019973848, 1.1089098748034592,
    1.003764886440033, -0.051037971088332391 },
  0.0,                                 /* Expression: SM.Sat
                                        * Referenced by: '<S56>/Constant1'
                                        */
  0.0,                                 /* Expression: SM.Sat
                                        * Referenced by: '<S56>/Constant3'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S56>/Switch'
                                        */

  /*  Expression: SM.IqdSign
   * Referenced by: '<S56>/change Iq Id  current signs'
   */
  { -1.0, -1.0, 1.0, 1.0, 1.0 },
  16.534055763786451,                  /* Expression: SM.ib
                                        * Referenced by: '<S58>/ib'
                                        */
  112.53999999999999,                  /* Expression: 1.1254*100
                                        * Referenced by: '<S33>/Unit Delay2'
                                        */
  0.0030618621784789728,               /* Expression: 1/SM.Vb
                                        * Referenced by: '<S122>/1_Vb'
                                        */
  0.0,                                 /* Expression: SM.Sat
                                        * Referenced by: '<S121>/Constant5'
                                        */
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S121>/Switch2'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S127>/Gain'
                                        */
  2.0E+6,                              /* Expression: SM.Gain1
                                        * Referenced by: '<S125>/Gain'
                                        */
  2.0E+6,                              /* Expression: SM.Gain1
                                        * Referenced by: '<S125>/Gain1'
                                        */

  /*  Expression: [SM.Ib2*ones(1,5),SM.N2,SM.Ib2,SM.Ib2,SM.Ib2,SM.phib,SM.phib,SM.Vb2,SM.Vb2,(SM.phib/SM.Ib2)*ones(1,2),1,1,1]
   * Referenced by: '<S120>/output formatting'
   */
  { 4082.48290463863, 4082.48290463863, 4082.48290463863, 4082.48290463863,
    4082.48290463863, 206.00993825915066, 4082.48290463863, 4082.48290463863,
    4082.48290463863, 1.0395957349782348, 1.0395957349782348, 326.59863237109039,
    326.59863237109039, 0.00025464790894703254, 0.00025464790894703254, 1.0, 1.0,
    1.0 },
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S119>/nominal speed'
                                        */
  0.0,                                 /* Expression: SM.dwo
                                        * Referenced by: '<S141>/dw_delay'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S141>/F2'
                                        */
  0.0,                                 /* Expression: SM.dwo
                                        * Referenced by: '<S141>/dw_predict'
                                        */
  157.07963267948966,                  /* Expression: SM.Nb
                                        * Referenced by: '<S119>/units'
                                        */
  0.0063661977236758134,               /* Expression: 1/(pi*50)
                                        * Referenced by: '<S2>/rad//s to pu'
                                        */
  2.0E+6,                              /* Expression: 2e6
                                        * Referenced by: '<S2>/pu2Watts'
                                        */
  0.0063661977236758134,               /* Expression: 1/SM.Nb
                                        * Referenced by: '<S54>/1_Pb'
                                        */

  /*  Expression: [1 -1]
   * Referenced by: '<S63>/1-1'
   */
  { 1.0, -1.0 },
  8100.0,                              /* Expression: SM.Pb
                                        * Referenced by: '<S76>/units1'
                                        */
  -9025.313417326428,                  /* Expression: sps.x0
                                        * Referenced by: '<S28>/Delay_x'
                                        */
  0.99869876364190868,                 /* Expression: sps.A
                                        * Referenced by: '<S28>/A'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<Root>/Reference speed (pu)'
                                        */
  0.50062421972534332,                 /* Expression: sps.D
                                        * Referenced by: '<S29>/D'
                                        */
  0.0,                                 /* Expression: sps.x0
                                        * Referenced by: '<S29>/Delay_x'
                                        */
  0.0012484394506866417,               /* Expression: sps.C
                                        * Referenced by: '<S29>/C'
                                        */
  27.70360110803324,                   /* Expression: sps.D
                                        * Referenced by: '<S30>/D'
                                        */
  0.0,                                 /* Expression: sps.x0
                                        * Referenced by: '<S30>/Delay_x'
                                        */
  -0.14835333948907356,                /* Expression: sps.C
                                        * Referenced by: '<S30>/C'
                                        */
  26.024727161826554,                  /* Expression: sps.B
                                        * Referenced by: '<S28>/B'
                                        */
  4.9967469091047721E-5,               /* Expression: sps.C
                                        * Referenced by: '<S28>/C'
                                        */
  0.00065061817904566387,              /* Expression: sps.D
                                        * Referenced by: '<S28>/D'
                                        */
  5.0E-5,                              /* Computed Parameter: DiscreteTimeIntegrator_gainval
                                        * Referenced by: '<S2>/Discrete-Time Integrator'
                                        */
  1.1,                                 /* Expression: Tlim(2)
                                        * Referenced by: '<S2>/Discrete-Time Integrator'
                                        */
  0.0,                                 /* Expression: Tlim(1)
                                        * Referenced by: '<S2>/Discrete-Time Integrator'
                                        */
  0.99750312109862671,                 /* Expression: sps.A
                                        * Referenced by: '<S29>/A'
                                        */
  0.99875156054931336,                 /* Expression: sps.B
                                        * Referenced by: '<S29>/B'
                                        */
  0.9944598337950139,                  /* Expression: sps.A
                                        * Referenced by: '<S30>/A'
                                        */
  0.997229916897507,                   /* Expression: sps.B
                                        * Referenced by: '<S30>/B'
                                        */
  10.4482,                             /* Expression: 10.4482
                                        * Referenced by: '<S33>/Unit Delay4'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S78>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S3>/Reference voltage (pu)'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S36>/Unit Delay1'
                                        */
  5.0E-5,                              /* Computed Parameter: DiscreteTimeIntegrator_gainva_n
                                        * Referenced by: '<S85>/Discrete-Time Integrator'
                                        */
  5.0,                                 /* Expression: UpperLimit
                                        * Referenced by: '<S85>/Discrete-Time Integrator'
                                        */
  -5.0,                                /* Expression: LowerLimit
                                        * Referenced by: '<S85>/Discrete-Time Integrator'
                                        */
  5.0,                                 /* Expression: UpperLimit
                                        * Referenced by: '<S85>/Saturation2'
                                        */
  -5.0,                                /* Expression: LowerLimit
                                        * Referenced by: '<S85>/Saturation2'
                                        */
  11.56,                               /* Expression: 11.56
                                        * Referenced by: '<S36>/pu to volts'
                                        */
  0.0030618621784789728,               /* Expression: SM.N
                                        * Referenced by: '<S52>/N'
                                        */
  0.0030618621784789728,               /* Expression: 1/SM.Vb
                                        * Referenced by: '<S57>/1_Vb'
                                        */

  /*  Expression: zeros(1, SM.nState-3)
   * Referenced by: '<S52>/[ Vkd =0 Vkq1=0  Vkq2=0 ]'
   */
  { 0.0, 0.0 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S67>/voltages'
                                        */
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S67>/IC'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S54>/Nominal speed'
                                        */
  314.15926535897933,                  /* Expression: SM.web
                                        * Referenced by: '<S54>/web1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S16>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S17>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S18>/do not delete this gain'
                                        */
  0.0030618621784789723,               /* Expression: Kv
                                        * Referenced by: '<S1>/Kv1'
                                        */
  2.0,                                 /* Expression: sps.k
                                        * Referenced by: '<S87>/sin(wt)'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S87>/sin(wt)'
                                        */
  0.015707317311820675,                /* Computed Parameter: sinwt_Hsin
                                        * Referenced by: '<S87>/sin(wt)'
                                        */
  0.99987663248166059,                 /* Computed Parameter: sinwt_HCos
                                        * Referenced by: '<S87>/sin(wt)'
                                        */
  -0.015707317311820675,               /* Computed Parameter: sinwt_PSin
                                        * Referenced by: '<S87>/sin(wt)'
                                        */
  0.99987663248166059,                 /* Computed Parameter: sinwt_PCos
                                        * Referenced by: '<S87>/sin(wt)'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval
                                        * Referenced by: '<S97>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S97>/Integ4'
                                        */
  0.02,                                /* Expression: sps.Delay
                                        * Referenced by: '<S97>/K1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size
   * Referenced by: '<S98>/S-Function'
   */
  { 1.0, 1.0 },
  0.0201,                              /* Expression: MaxDelay
                                        * Referenced by: '<S98>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size
   * Referenced by: '<S98>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S98>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size
   * Referenced by: '<S98>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S98>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size
   * Referenced by: '<S98>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S98>/S-Function'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S97>/K2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S97>/Unit Delay'
                                        */
  1.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S97>/Unit Delay1'
                                        */
  2.0,                                 /* Expression: sps.k
                                        * Referenced by: '<S87>/cos(wt)'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S87>/cos(wt)'
                                        */
  0.015707317311820675,                /* Computed Parameter: coswt_Hsin
                                        * Referenced by: '<S87>/cos(wt)'
                                        */
  0.99987663248166059,                 /* Computed Parameter: coswt_HCos
                                        * Referenced by: '<S87>/cos(wt)'
                                        */
  0.99987663248166059,                 /* Computed Parameter: coswt_PSin
                                        * Referenced by: '<S87>/cos(wt)'
                                        */
  0.015707317311820648,                /* Computed Parameter: coswt_PCos
                                        * Referenced by: '<S87>/cos(wt)'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_o
                                        * Referenced by: '<S95>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S95>/Integ4'
                                        */
  0.02,                                /* Expression: sps.Delay
                                        * Referenced by: '<S95>/K1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_j
   * Referenced by: '<S96>/S-Function'
   */
  { 1.0, 1.0 },
  0.0201,                              /* Expression: MaxDelay
                                        * Referenced by: '<S96>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_l
   * Referenced by: '<S96>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S96>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_k
   * Referenced by: '<S96>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S96>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_f
   * Referenced by: '<S96>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S96>/S-Function'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S95>/K2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S95>/Unit Delay'
                                        */
  0.0,                                 /* Expression: sps.Vinit
                                        * Referenced by: '<S95>/Unit Delay1'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S87>/Rad->Deg.'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S86>/deg->rad'
                                        */
  2.0,                                 /* Expression: sps.k
                                        * Referenced by: '<S88>/sin(wt)'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S88>/sin(wt)'
                                        */
  0.015707317311820675,                /* Computed Parameter: sinwt_Hsin_j
                                        * Referenced by: '<S88>/sin(wt)'
                                        */
  0.99987663248166059,                 /* Computed Parameter: sinwt_HCos_p
                                        * Referenced by: '<S88>/sin(wt)'
                                        */
  -0.015707317311820675,               /* Computed Parameter: sinwt_PSin_l
                                        * Referenced by: '<S88>/sin(wt)'
                                        */
  0.99987663248166059,                 /* Computed Parameter: sinwt_PCos_d
                                        * Referenced by: '<S88>/sin(wt)'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_d
                                        * Referenced by: '<S103>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S103>/Integ4'
                                        */
  0.02,                                /* Expression: sps.Delay
                                        * Referenced by: '<S103>/K1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_jj
   * Referenced by: '<S104>/S-Function'
   */
  { 1.0, 1.0 },
  0.0201,                              /* Expression: MaxDelay
                                        * Referenced by: '<S104>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_d
   * Referenced by: '<S104>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S104>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_a
   * Referenced by: '<S104>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S104>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_c
   * Referenced by: '<S104>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S104>/S-Function'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S103>/K2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S103>/Unit Delay'
                                        */
  -0.49999999999999978,                /* Expression: sps.Vinit
                                        * Referenced by: '<S103>/Unit Delay1'
                                        */
  2.0,                                 /* Expression: sps.k
                                        * Referenced by: '<S88>/cos(wt)'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S88>/cos(wt)'
                                        */
  0.015707317311820675,                /* Computed Parameter: coswt_Hsin_p
                                        * Referenced by: '<S88>/cos(wt)'
                                        */
  0.99987663248166059,                 /* Computed Parameter: coswt_HCos_m
                                        * Referenced by: '<S88>/cos(wt)'
                                        */
  0.99987663248166059,                 /* Computed Parameter: coswt_PSin_k
                                        * Referenced by: '<S88>/cos(wt)'
                                        */
  0.015707317311820648,                /* Computed Parameter: coswt_PCos_n
                                        * Referenced by: '<S88>/cos(wt)'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_e
                                        * Referenced by: '<S101>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S101>/Integ4'
                                        */
  0.02,                                /* Expression: sps.Delay
                                        * Referenced by: '<S101>/K1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_a
   * Referenced by: '<S102>/S-Function'
   */
  { 1.0, 1.0 },
  0.0201,                              /* Expression: MaxDelay
                                        * Referenced by: '<S102>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_i
   * Referenced by: '<S102>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S102>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_b
   * Referenced by: '<S102>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S102>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_cp
   * Referenced by: '<S102>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S102>/S-Function'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S101>/K2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S101>/Unit Delay'
                                        */
  -0.86602540378443871,                /* Expression: sps.Vinit
                                        * Referenced by: '<S101>/Unit Delay1'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S88>/Rad->Deg.'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S86>/deg->rad1'
                                        */
  2.0,                                 /* Expression: sps.k
                                        * Referenced by: '<S89>/sin(wt)'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S89>/sin(wt)'
                                        */
  0.015707317311820675,                /* Computed Parameter: sinwt_Hsin_d
                                        * Referenced by: '<S89>/sin(wt)'
                                        */
  0.99987663248166059,                 /* Computed Parameter: sinwt_HCos_k
                                        * Referenced by: '<S89>/sin(wt)'
                                        */
  -0.015707317311820675,               /* Computed Parameter: sinwt_PSin_k
                                        * Referenced by: '<S89>/sin(wt)'
                                        */
  0.99987663248166059,                 /* Computed Parameter: sinwt_PCos_c
                                        * Referenced by: '<S89>/sin(wt)'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_l
                                        * Referenced by: '<S109>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S109>/Integ4'
                                        */
  0.02,                                /* Expression: sps.Delay
                                        * Referenced by: '<S109>/K1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_l
   * Referenced by: '<S110>/S-Function'
   */
  { 1.0, 1.0 },
  0.0201,                              /* Expression: MaxDelay
                                        * Referenced by: '<S110>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_g
   * Referenced by: '<S110>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S110>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_d
   * Referenced by: '<S110>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S110>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_a
   * Referenced by: '<S110>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S110>/S-Function'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S109>/K2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S109>/Unit Delay'
                                        */
  -0.49999999999999978,                /* Expression: sps.Vinit
                                        * Referenced by: '<S109>/Unit Delay1'
                                        */
  2.0,                                 /* Expression: sps.k
                                        * Referenced by: '<S89>/cos(wt)'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S89>/cos(wt)'
                                        */
  0.015707317311820675,                /* Computed Parameter: coswt_Hsin_ph
                                        * Referenced by: '<S89>/cos(wt)'
                                        */
  0.99987663248166059,                 /* Computed Parameter: coswt_HCos_b
                                        * Referenced by: '<S89>/cos(wt)'
                                        */
  0.99987663248166059,                 /* Computed Parameter: coswt_PSin_j
                                        * Referenced by: '<S89>/cos(wt)'
                                        */
  0.015707317311820648,                /* Computed Parameter: coswt_PCos_d
                                        * Referenced by: '<S89>/cos(wt)'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_a
                                        * Referenced by: '<S107>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S107>/Integ4'
                                        */
  0.02,                                /* Expression: sps.Delay
                                        * Referenced by: '<S107>/K1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_b
   * Referenced by: '<S108>/S-Function'
   */
  { 1.0, 1.0 },
  0.0201,                              /* Expression: MaxDelay
                                        * Referenced by: '<S108>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_la
   * Referenced by: '<S108>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S108>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_a3
   * Referenced by: '<S108>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S108>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_n
   * Referenced by: '<S108>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S108>/S-Function'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S107>/K2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S107>/Unit Delay'
                                        */
  0.86602540378443871,                 /* Expression: sps.Vinit
                                        * Referenced by: '<S107>/Unit Delay1'
                                        */
  57.295779513082323,                  /* Expression: 180/pi
                                        * Referenced by: '<S89>/Rad->Deg.'
                                        */
  0.017453292519943295,                /* Expression: pi/180
                                        * Referenced by: '<S86>/deg->rad2'
                                        */
  1.0,                                 /* Expression: sps.PosOn
                                        * Referenced by: '<S86>/Constant'
                                        */
  0.0,                                 /* Expression: sps.NegOn
                                        * Referenced by: '<S86>/Constant1'
                                        */
  0.0,                                 /* Expression: sps.ZeroOn
                                        * Referenced by: '<S86>/Constant2'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S40>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S41>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S42>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: Ki
                                        * Referenced by: '<S31>/Kv'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S43>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S44>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S45>/do not delete this gain'
                                        */
  0.0030618621784789723,               /* Expression: Kv
                                        * Referenced by: '<S31>/Kv1'
                                        */

  /*  Expression: zeros(1,6)
   * Referenced by: '<S34>/g'
   */
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
  0.00010300496912957533,              /* Expression: SM.N
                                        * Referenced by: '<S117>/N'
                                        */

  /*  Expression: [1 -1]
   * Referenced by: '<S128>/1-1'
   */
  { 1.0, -1.0 },

  /*  Expression: zeros(1, SM.nState-3)
   * Referenced by: '<S117>/[ Vkd =0 Vkq1=0  Vkq2=0 ]'
   */
  { 0.0, 0.0 },
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S132>/voltages'
                                        */
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S132>/IC'
                                        */
  5.0E-7,                              /* Expression: 1/SM.Pb
                                        * Referenced by: '<S119>/1_Pb'
                                        */
  0.01973920880217872,                 /* Expression: SM.F
                                        * Referenced by: '<S119>/F'
                                        */
  0.71855682994995107,                 /* Expression: 1/(2*SM.H)
                                        * Referenced by: '<S119>/1 ----- 2H'
                                        */
  2.5E-5,                              /* Computed Parameter: Rotorspeeddeviationdw_gainval
                                        * Referenced by: '<S119>/Rotor speed deviation (dw)'
                                        */
  0.0,                                 /* Expression: SM.dwo
                                        * Referenced by: '<S119>/Rotor speed deviation (dw)'
                                        */
  314.15926535897933,                  /* Expression: SM.web
                                        * Referenced by: '<S119>/we base'
                                        */
  2.5E-5,                              /* Computed Parameter: Integ4_gainval_j
                                        * Referenced by: '<S111>/Integ4'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S111>/Integ4'
                                        */
  0.02,                                /* Expression: sps.Delay
                                        * Referenced by: '<S111>/K1'
                                        */

  /*  Computed Parameter: SFunction_P1_Size_jv
   * Referenced by: '<S112>/S-Function'
   */
  { 1.0, 1.0 },
  0.0201,                              /* Expression: MaxDelay
                                        * Referenced by: '<S112>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P2_Size_k
   * Referenced by: '<S112>/S-Function'
   */
  { 1.0, 1.0 },
  5.0E-5,                              /* Expression: Ts
                                        * Referenced by: '<S112>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P3_Size_i
   * Referenced by: '<S112>/S-Function'
   */
  { 1.0, 1.0 },
  0.0,                                 /* Expression: InitialValue
                                        * Referenced by: '<S112>/S-Function'
                                        */

  /*  Computed Parameter: SFunction_P4_Size_p
   * Referenced by: '<S112>/S-Function'
   */
  { 1.0, 1.0 },
  1.0,                                 /* Expression: DFT
                                        * Referenced by: '<S112>/S-Function'
                                        */
  50.0,                                /* Expression: sps.Freq
                                        * Referenced by: '<S111>/K2'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S111>/Unit Delay'
                                        */
  10.5,                                /* Expression: sps.Vinit
                                        * Referenced by: '<S111>/Unit Delay1'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S13>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S14>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S15>/do not delete this gain'
                                        */
  1.0,                                 /* Expression: Ki
                                        * Referenced by: '<S1>/Kv'
                                        */
  0.0,                                 /* Expression: Gain
                                        * Referenced by: '<S5>/Gain2'
                                        */
  0.0,                                 /* Expression: Display
                                        * Referenced by: '<S5>/Mode'
                                        */
  0.0,                                 /* Expression: Gain
                                        * Referenced by: '<S6>/Gain2'
                                        */
  0.0,                                 /* Expression: Display
                                        * Referenced by: '<S6>/Mode'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S146>/0 1'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S146>/Gain'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S146>/Delay Ts'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S146>/Delay Ts '
                                        */
  0.0,                                 /* Expression: External
                                        * Referenced by: '<S8>/C4'
                                        */

  /*  Expression: sps.tv
   * Referenced by: '<S155>/Look-Up Table'
   */
  { 0.0, 2.9999500000000006, 2.9999500000000006, 3.9999500000000006 },

  /*  Expression: sps.opv
   * Referenced by: '<S155>/Look-Up Table'
   */
  { 0.0, 0.0, 1.0, 1.0 },
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S8>/Switch3'
                                        */
  0.0,                                 /* Expression: InitialState
                                        * Referenced by: '<S8>/Constant5'
                                        */
  1.0,                                 /* Expression: BR.com
                                        * Referenced by: '<S146>/C4'
                                        */

  /*  Expression: sps.tv
   * Referenced by: '<S148>/Look-Up Table'
   */
  { 0.0, 999999.99995000008, 999999.99995000008, 1.0000009999500001E+6 },

  /*  Expression: sps.opv
   * Referenced by: '<S148>/Look-Up Table'
   */
  { 0.0, 0.0, 0.0, 0.0 },
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S146>/Switch3'
                                        */
  1.5,                                 /* Expression: 1.5
                                        * Referenced by: '<S146>/>1.5'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S146>/>1.5'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S146>/>1.5'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S146>/>1.5'
                                        */
  0.0,                                 /* Expression: double(InitialState)
                                        * Referenced by: '<S146>/IC'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S149>/0 1'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S149>/Gain'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S149>/Delay Ts'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S149>/Delay Ts '
                                        */
  1.0,                                 /* Expression: BR.com
                                        * Referenced by: '<S149>/C4'
                                        */

  /*  Expression: sps.tv
   * Referenced by: '<S151>/Look-Up Table'
   */
  { 0.0, 999999.99995000008, 999999.99995000008, 1.0000009999500001E+6 },

  /*  Expression: sps.opv
   * Referenced by: '<S151>/Look-Up Table'
   */
  { 0.0, 0.0, 0.0, 0.0 },
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S149>/Switch3'
                                        */
  1.5,                                 /* Expression: 1.5
                                        * Referenced by: '<S149>/>1.5'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S149>/>1.5'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S149>/>1.5'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S149>/>1.5'
                                        */
  0.0,                                 /* Expression: double(InitialState)
                                        * Referenced by: '<S149>/IC'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S152>/0 1'
                                        */
  2.0,                                 /* Expression: 2
                                        * Referenced by: '<S152>/Gain'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S152>/Delay Ts'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S152>/Delay Ts '
                                        */
  1.0,                                 /* Expression: BR.com
                                        * Referenced by: '<S152>/C4'
                                        */

  /*  Expression: sps.tv
   * Referenced by: '<S154>/Look-Up Table'
   */
  { 0.0, 999999.99995000008, 999999.99995000008, 1.0000009999500001E+6 },

  /*  Expression: sps.opv
   * Referenced by: '<S154>/Look-Up Table'
   */
  { 0.0, 0.0, 0.0, 0.0 },
  0.5,                                 /* Expression: 0.5
                                        * Referenced by: '<S152>/Switch3'
                                        */
  1.5,                                 /* Expression: 1.5
                                        * Referenced by: '<S152>/>1.5'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S152>/>1.5'
                                        */
  1.0,                                 /* Expression: 1
                                        * Referenced by: '<S152>/>1.5'
                                        */
  0.0,                                 /* Expression: 0
                                        * Referenced by: '<S152>/>1.5'
                                        */
  0.0,                                 /* Expression: double(InitialState)
                                        * Referenced by: '<S152>/IC'
                                        */
  0,                                   /* Expression: SM.nState==6
                                        * Referenced by: '<S66>/Constant1'
                                        */
  0,                                   /* Expression: SM.nState==6
                                        * Referenced by: '<S66>/Constant2'
                                        */
  0,                                   /* Expression: SM.nState==6
                                        * Referenced by: '<S131>/Constant1'
                                        */
  0,                                   /* Expression: SM.nState==6
                                        * Referenced by: '<S131>/Constant2'
                                        */
  0,                                   /* Expression: SM.nState==6
                                        * Referenced by: '<S131>/Constant3'
                                        */
  0,                                   /* Expression: SM.nState==6
                                        * Referenced by: '<S120>/Constant1'
                                        */

  /* Start of '<S131>/Lmq_sat' */
  {
    1.5099279691315945,                /* Expression: SM.Lmsatq(1)
                                        * Referenced by: '<S134>/Constant1'
                                        */

    /*  Expression: SM.One_Llq
     * Referenced by: '<S139>/1//Ll_q'
     */
    { 20.003763467952279, 4.8933110866070813 },

    /*  Expression: [ 1/SM.Ll 1/SM.Llkq1 1/SM.Llkq2]
     * Referenced by: '<S138>/u2'
     */
    { 20.003763467952279, 4.8933110866070813, 0.0 },
    1.5099279691315945,                /* Expression: SM.Lmsatq(1)
                                        * Referenced by: '<S134>/Lmq_sat'
                                        */

    /*  Expression: SM.Phisat
     * Referenced by: '<S134>/Lookup Table'
     */
    { 0.0, 1.0 },

    /*  Expression: [ 0 SM.Phisat(2:end)./SM.Lmsatq(2:end)*SM.Lmq ]
     * Referenced by: '<S134>/Lookup Table'
     */
    { 0.0, 1.0 },
    1.5099279691315945                 /* Expression: SM.Lmq
                                        * Referenced by: '<S134>/Lmq'
                                        */
  }
  /* End of '<S131>/Lmq_sat' */
  ,

  /* Start of '<S86>/Pos. Seq. Computation' */
  {
    0.33333333333333331,               /* Expression: 1/3
                                        * Referenced by: '<S91>/Gain3'
                                        */

    /*  Expression: [1 sps.a sps.a2]
     * Referenced by: '<S91>/Gain1'
     */
    { { 1.0, 0.0 }, { -0.49999999999999978, 0.86602540378443871 }, { -
        0.49999999999999978, -0.86602540378443871 } }
  }
  /* End of '<S86>/Pos. Seq. Computation' */
  ,

  /* Start of '<S86>/Neg. Seq. Computation' */
  {
    0.33333333333333331,               /* Expression: 1/3
                                        * Referenced by: '<S90>/Gain3'
                                        */

    /*  Expression: [1 sps.a2 sps.a]
     * Referenced by: '<S90>/Gain1'
     */
    { { 1.0, 0.0 }, { -0.49999999999999978, -0.86602540378443871 }, { -
        0.49999999999999978, 0.86602540378443871 } }
  }
  /* End of '<S86>/Neg. Seq. Computation' */
  ,

  /* Start of '<S66>/Lmq_sat' */
  {
    0.82304818785531353,               /* Expression: SM.Lmsatq(1)
                                        * Referenced by: '<S69>/Constant1'
                                        */

    /*  Expression: SM.One_Llq
     * Referenced by: '<S74>/1//Ll_q'
     */
    { 13.889115727965949, 6.1946824532514153 },

    /*  Expression: [ 1/SM.Ll 1/SM.Llkq1 1/SM.Llkq2]
     * Referenced by: '<S73>/u2'
     */
    { 13.889115727965949, 6.1946824532514153, 0.0 },
    0.82304818785531353,               /* Expression: SM.Lmsatq(1)
                                        * Referenced by: '<S69>/Lmq_sat'
                                        */

    /*  Expression: SM.Phisat
     * Referenced by: '<S69>/Lookup Table'
     */
    { 0.0, 1.0 },

    /*  Expression: [ 0 SM.Phisat(2:end)./SM.Lmsatq(2:end)*SM.Lmq ]
     * Referenced by: '<S69>/Lookup Table'
     */
    { 0.0, 1.0 },
    0.82304818785531353                /* Expression: SM.Lmq
                                        * Referenced by: '<S69>/Lmq'
                                        */
  }
  /* End of '<S66>/Lmq_sat' */
};

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
