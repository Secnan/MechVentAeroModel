;; 1. Based on: Simulation_9.0
;; 2. Description: M_DEP_AA
;; x1. Author: user

$SIZES DIMNEW=-5000  
$PROBLEM Semi-Mechanistic Model for Aerosol Delivery in Mechanical Ventilation
; Humidifier 90*2
; 目标：提高 $DES 中沉积效率计算的可读性与维护性

$DATA Simulation_8.csv IGNORE=@
$INPUT ID TIME AMT DV CMT

$SUBROUTINE ADVAN9 TOL=9
$MODEL
; Define compartments for suspended (M) and deposited (A) masses in each room.
; Compartments numbered sequentially.

COMP=(M_DRY) ; 1: Suspended in Humidifier Dry End 
COMP=(A_DRY) ; 2: Deposited in Humidifier Dry End 
COMP=(M_HUM) ; 3: Suspended in Humidifier
COMP=(A_HUM) ; 4: Deposited in Humidifier
COMP=(M_INSP) ; 5: Suspended in Inspiratory Limb
COMP=(A_INSP) ; 6: Deposited in Inspiratory Limb
COMP=(M_Y) ; 7: Suspended in Y-piece
COMP=(A_Y) ; 8: Deposited in Y-piece
COMP=(M_EXT) ; 9: Suspended in Extension tube
COMP=(A_EXT) ; 10: Deposited in Extension tube
COMP=(M_AA) ; 11: Suspended in Artificial Airway
COMP=(A_AA) ; 12: Deposited in Artificial Airway
COMP=(A_LUNG) ; 13: Deposited in Lung (no suspended, fully captured)
COMP=(M_EXP) ; 14: Suspended in Expiratory Limb
COMP=(A_EXP) ; 15: Deposited in Expiratory Limb
COMP=(A_FILT) ; 16: Deposited in Filter (no suspended, fully captured)
COMP=(M_RES) ; 17: Reservoir mass for nebulizer (for R(t) calculation)

$ABBREVIATED DERIV2=NO COMRES=150 COMSAV=1
$PK
; ─────────────── 物理常数 ───────────────
GRAV = 9.81 ; Gravity (m/s^2)
MU = 1.81E-5 ; Air viscosity (Pa*s)
LAMBDA = 6.8E-8; Mean free path (m)
KB = 1.38E-23 ; Boltzmann constant (J/K)
TEMP = 310 ; Temperature (K, body temp ~37C)
PI = 3.14159

; ─────────────── 气溶胶性质 ───────────────
  D_P    = THETA(1)           ; 空气动力学直径 (m)
  RHOP   = THETA(2)           ; 粒子密度 (kg/m³)
  CC     = 1 + (2*LAMBDA/D_P)*(1.257 + 0.4*EXP(-1.1*D_P/(2*LAMBDA)))
  DB     = (KB * TEMP * CC) / (3 * PI * MU * D_P)   ; 扩散系数 (m²/s)

; ─────────────── 雾化器参数 ───────────────
  CSOL   = THETA(3)
  V_NEB  = THETA(4)
  M_RES0 = CSOL * V_NEB
  K_NEB  = THETA(5)
  M_SW   = THETA(18)          ; 切换阈值 (mg)
  R_CONST= K_NEB * M_RES0     ; 初期恒定速率

  ; ─────────────── 通气机参数 ───────────────
  VT     = THETA(6)
  FREQ   = THETA(7)
  IE     = THETA(8)
  TCYCLE = 60 / FREQ
  T_INSP = TCYCLE / (1 + IE)
  T_EXP  = TCYCLE - T_INSP


; ─────────────── 各段体积 & 几何 (L & m) ───────────────
  V_DRY = THETA(16) ; Dry End Volume (L)
  D_DRY = 0.019 ; m
  L_DRY = 0.1

  V_HUM  = THETA(9)
  D_HUM  = 0.019   ; m
  L_HUM  = 0.06

  V_INSP = THETA(10)
  D_INSP = 0.019
  L_INSP = 1.8

  V_Y    = THETA(11)
  D_Y    = 0.019
  L_Y    = 0.05
  THETA_Y   = PI/4

  V_EXT  = THETA(12)
  D_EXT  = 0.014
  L_EXT  = 0.3
  THETA_EXT = PI/3
  ALPHA_EXT = 0

  V_AA   = THETA(13)
  D_AA   = 0.007
  L_AA   = 0.28
  THETA_AA  = PI*40/180
  ALPHA_AA = PI*50/180   ; 使用弯曲角度计算倾角 (假设)

  V_EXP  = THETA(14)
  D_EXP  = 0.019
  L_EXP  = 1.8

  ; 吸气支分段（向上）
  L_SEC1_INSP = 0.45
  ALPHA1_INSP = 0
  THETA1_INSP = 3*PI/4
  L_SEC2_INSP = 0.9
  ALPHA2_INSP = PI/4
  THETA2_INSP = 3*PI/4
  L_SEC3_INSP = 0.45
  ALPHA3_INSP = 0

  ; 呼气支分段（向下）
  L_SEC1_EXP  = 0.45
  ALPHA1_EXP  = 0
  THETA1_EXP  = 3*PI/4
  L_SEC2_EXP  = 0.9
  ALPHA2_EXP  = PI/4
  THETA2_EXP  = 3*PI/4
  L_SEC3_EXP  = 0.45
  ALPHA3_EXP  = 0

  F_HUM  = THETA(15)   ; 湿化器界面沉积常数
  F_DRY = THETA(17) ; 干端界面沉积常数
  
  W_DEP_HUM = 5.4
  W_DEP_INSP = 4.2
  W_DEP_Y = 12
  W_DEP_EXT  = 22
  W_DEP_AA  = 60
  W_DEP_EXP   = 12
  K_V    = 0.2
  K_M    = 0.1
  K_M_AA = 0.3

$DES
  ; ─────────────── 通气相位判断 ───────────────
TCUR   = MOD(T, TCYCLE)
  IT     = 0
  IF (TCUR .LE. T_INSP) IT = 1

  TPHASE = TCUR
  IF (IT .EQ. 0) TPHASE = TCUR - T_INSP

  ; 流量 (L/s) → QM (m³/s)
  Q_INSP = 0
  Q_EXP  = 0
  IF (IT .EQ. 1) THEN
    Q_INSP = (PI * VT / (2 * T_INSP)) * SIN((PI / T_INSP) * TPHASE)
  ELSE
    Q_EXP  = (PI * VT / (2 * T_EXP )) * SIN((PI / T_EXP ) * TPHASE)
  ENDIF

  QM_INSP = Q_INSP * 0.001
  QM_EXP  = Q_EXP  * 0.001

  ; ─────────────── 雾化输出速率 ───────────────
  ; 分段雾化速率
MRES = A(17)

IF (MRES .GT. M_SW) THEN
  RT = R_CONST
ELSE
  RT = K_NEB * MRES
ENDIF

; 防止出现负值
RT = MAX(0.0D0, RT)

  ; ─────────────── 浓度 (mg/L) ───────────────
C_DRY = A(1)/V_DRY
C_HUM = A(3)/V_HUM
C_INSP = A(5)/V_INSP
C_Y = A(7)/V_Y
C_EXT = A(9)/V_EXT
C_AA = A(11)/V_AA
C_EXP = A(14)/V_EXP

; ───────────────────────────────────────────────
; 1. 公共物理量（所有节段共用）
; ───────────────────────────────────────────────
  TAUP = (RHOP * D_P**2 * CC) / (18 * MU)
  VS   = GRAV * TAUP
  
  RR0    = 1           ; daughter/parent 半径比（暂固定）
  RR02 = RR0**2
  RR04 = RR0**4


; ───────────────────────────────────────────────
; 0. 公共物理量（所有节段共用）
; ───────────────────────────────────────────────
TAUP = (RHOP * D_P**2 * CC) / (18 * MU)
VS   = GRAV * TAUP

; ───────────────────────────────────────────────
; 1. 干端湿化器 (Dry End Humidifier) ─ 喷射撞击沉积
; ───────────────────────────────────────────────
QM_DRY = QM_INSP
EINT_DRY = F_DRY * (RT / R_CONST)
ETOT_DRY = EINT_DRY
KDEP_DRY = ETOT_DRY * Q_INSP / V_DRY

; ───────────────────────────────────────────────
; 2. 湿化器 (Humidifier) ─ 两个90°弯 + 沉降 + 扩散 + 界面
; ───────────────────────────────────────────────
QM_HUM   = QM_INSP
AREA_HUM = PI * D_HUM**2 / 4
U_HUM    = VT / T_INSP * 0.001/AREA_HUM
R_HUM    = D_HUM / 2
ST0_HUM  = (RHOP * CC * D_P**2 * U_HUM) / (36 * MU * R_HUM)

ALPHA_HUM = PI/2
COS2A_HUM = COS(ALPHA_HUM)**2
SIN_A_HUM = SIN(ALPHA_HUM)

F0_HUM = PI - (PI/4 + (5/4*PI - 8/3)*COS2A_HUM)*RR02
F1_HUM = 1 + (-1/3 + (PI - 11/3)*COS2A_HUM - SIN_A_HUM/3)*RR02 + ((2/3 - PI/8)*COS2A_HUM + SIN_A_HUM**2/5 + (6 - 15/8*PI)*COS2A_HUM**2 + (7/15 - PI/8)*SIN_A_HUM**2*COS2A_HUM)*RR04

EI_HUM1 = (8 * SIN(ALPHA_HUM) * F1_HUM) / (RR0 * F0_HUM) * ST0_HUM
EI_HUM2 = EI_HUM1   ; 第二个90°弯相同

EINT_HUM = F_HUM * (RT / R_CONST)

ETOT_HUM = (1 - (1-EI_HUM1)*(1-EI_HUM2)*(1-EINT_HUM))
ETOT_HUM = MIN(0.999, ETOT_HUM)
KDEP_HUM =  W_DEP_HUM *ETOT_HUM * Q_INSP / V_HUM

; ───────────────────────────────────────────────
; 3. 吸气肢 (Inspiratory Limb) ─ 三段直管 + 两弯
; ───────────────────────────────────────────────
RR0 = 0.5
RR02 = RR0**2
RR04 = RR0**4
AREA_INSP = PI * D_INSP**2 / 4
U_INSP    = VT / T_INSP * 0.001/AREA_INSP
R_INSP    = D_INSP / 2
ST0_INSP  = (RHOP * CC * D_P**2 * U_INSP) / (36 * MU * R_INSP)

; Sec1 ────────────────────────────────
; Bend1
COS2A1_INSP = COS(THETA1_INSP)**2
SIN_A1_INSP = SIN(THETA1_INSP)
F0_BEND1_INSP = PI - (PI/4 + (5/4*PI - 8/3)*COS2A1_INSP)*RR02
F1_BEND1_INSP = 1 + (-1/3 + (PI - 11/3)*COS2A1_INSP - SIN_A1_INSP/3)*RR02 + ((2/3 - PI/8)*COS2A1_INSP + SIN_A1_INSP**2/5 + (6 - 15/8*PI)*COS2A1_INSP**2 + (7/15 - PI/8)*SIN_A1_INSP**2*COS2A1_INSP)*RR04
EI_BEND1 = (8 * SIN(THETA1_INSP) * F1_BEND1_INSP) / (RR0 * F0_BEND1_INSP) * ST0_INSP

; Sec2 ────────────────────────────────
; 沉降 Sec2
; uphill
GAMMA_SEC2 = ((3*VS*L_SEC2_INSP)/(8*U_INSP*R_INSP) * COS(ALPHA2_INSP))**(2.0/3.0) / (1 - VS/(2*U_INSP)*SIN(ALPHA2_INSP))
SIGMA_SEC2 = (VS/(6*U_INSP)*SIN(ALPHA2_INSP)) / (1 - VS/(2*U_INSP)*SIN(ALPHA2_INSP))
COND1 = PI/2 - (2*R_INSP/(9*L_SEC2_INSP))*SQRT(2*VS/(3*U_INSP))
COND2 = (PI/2 - R_INSP*VS/(8*U_INSP*L_SEC2_INSP))

; 初始化 ES_SEC2 并使用一系列独立的 IF 语句进行赋值
ES_SEC2 = 0 ; 初始化默认值
IF (ALPHA2_INSP .LT. COND1 .AND. ALPHA2_INSP .GT. 0) THEN
; 值1: 对应条件 (ALPHA2_INSP < COND1 AND ALPHA2_INSP > 0)
ES_SEC2_VAL1 = (1/PI) * (3*SQRT(SIGMA_SEC2*(1-SIGMA_SEC2)) + ASIN(SQRT(1-SIGMA_SEC2)) + (1 - 9*SIGMA_SEC2**2)*ASIN(SQRT((1-SIGMA_SEC2)/(1+3*SIGMA_SEC2)))) - (2/PI)*(SQRT(GAMMA_SEC2*(1-GAMMA_SEC2))*(1-2*GAMMA_SEC2) + ASIN(SQRT((1-GAMMA_SEC2)/(1+3*GAMMA_SEC2))))
  ES_SEC2 = ES_SEC2_VAL1
ENDIF

IF (ALPHA2_INSP .LT. COND2  .AND. ALPHA2_INSP .GT. COND1) THEN
; 值2: 对应条件 (ALPHA2_INSP < COND2 AND ALPHA2_INSP > COND1)
ZETA0_SEC2 = (R_INSP*VS*SIN(ALPHA2_INSP)**2)/(8*U_INSP*L_SEC2_INSP*COS(ALPHA2_INSP)) - (1/16)*SQRT(VS/(6*U_INSP))*SIN(ALPHA2_INSP) + (7*L_SEC2_INSP)/(8*R_INSP)*(1/TAN(ALPHA2_INSP))
ES_SEC2_VAL2 = (1 + 3*SIGMA_SEC2)/PI * SQRT(1 - ZETA0_SEC2**2) * (ZETA0_SEC2 + 4*GAMMA_SEC2**(3.0/2.0)/SQRT(1 + 3*SIGMA_SEC2) - SQRT(ZETA0_SEC2**2 - 3*SIGMA_SEC2/(1 + 3*SIGMA_SEC2))) + (1 - 9*SIGMA_SEC2**2)/PI * ASIN(SQRT(1 - ZETA0_SEC2**2)) - (1/PI) * ASIN(SQRT((1 + 3*SIGMA_SEC2)*(1 - ZETA0_SEC2**2)))
  ES_SEC2 = ES_SEC2_VAL2
ENDIF

IF (ALPHA2_INSP .GT. COND2 ) THEN
; 值3: 对应条件 (ALPHA2_INSP > COND2)
ES_SEC2_VAL3 = 0
  ES_SEC2 = ES_SEC2_VAL3
ENDIF
; *** 代码修正结束 ***

; Bend2
EI_BEND2 = (8 * SIN(THETA2_INSP) * F1_BEND1_INSP) / (RR0 * F0_BEND1_INSP) * ST0_INSP   ; 角度相同，使用相同 F0 F1

; Sec3 ────────────────────────────────

; 总效率
ETOT_INSP = (1 - (1-EI_BEND1)*(1-ES_SEC2)*(1-EI_BEND2))
ETOT_INSP = MIN(0.999, ETOT_INSP)

KDEP_INSP =  W_DEP_INSP *ETOT_INSP * Q_INSP / V_INSP

; ───────────────────────────────────────────────
; 4. Y-piece ─ 主要撞击 + 少量沉降扩散
; ───────────────────────────────────────────────
RR0 = 1
RR02 = RR0**2
RR04 = RR0**4
AREA_Y = PI * D_Y**2 / 4
U_Y    = 2*VT / (T_INSP+ T_EXP) * 0.001/AREA_Y
R_Y    = D_Y / 2
ST0_Y  = (RHOP * CC * D_P**2 * U_Y) / (36 * MU * R_Y)

COS2A_Y = COS(THETA_Y)**2
SIN_A_Y = SIN(THETA_Y)

F0_Y = PI - (PI/4 + (5/4*PI - 8/3)*COS2A_Y)*RR02
F1_Y = 1 + (-1/3 + (PI - 11/3)*COS2A_Y - SIN_A_Y/3)*RR02 + ((2/3 - PI/8)*COS2A_Y + SIN_A_Y**2/5 + (6 - 15/8*PI)*COS2A_Y**2 + (7/15 - PI/8)*SIN_A_Y**2*COS2A_Y)*RR04

EI_Y = (8 * SIN(THETA_Y) * F1_Y) / (RR0 * F0_Y) * ST0_Y

ETOT_Y = EI_Y
KDEP_Y = W_DEP_Y * (IT *ETOT_Y * Q_INSP / V_Y + (1-IT) *ETOT_Y * Q_EXP / V_Y)

; ───────────────────────────────────────────────
; 5. Extension tube ─ 直管 + 一弯
; ───────────────────────────────────────────────
IF (IT .EQ. 1) THEN
RR0 = 14/19
ELSE
RR0 =2
ENDIF

AREA_EXT = PI * D_EXT**2 / 4
U_EXT    = 2*VT / (T_INSP+ T_EXP) * 0.001/AREA_EXT
R_EXT    = D_EXT / 2
ST0_EXT  = (RHOP * CC * D_P**2 * U_EXT) / (36 * MU * R_EXT)

COS2A_EXT = COS(THETA_EXT)**2
SIN_A_EXT = SIN(THETA_EXT)

F0_EXT = PI - (PI/4 + (5/4*PI - 8/3)*COS2A_EXT)*RR02
F1_EXT = 1 + (-1/3 + (PI - 11/3)*COS2A_EXT - SIN_A_EXT/3)*RR02 + ((2/3 - PI/8)*COS2A_EXT + SIN_A_EXT**2/5 + (6 - 15/8*PI)*COS2A_EXT**2 + (7/15 - PI/8)*SIN_A_EXT**2*COS2A_EXT)*RR04

EI_EXT = (8 * SIN(THETA_EXT) * F1_EXT) / (RR0 * F0_EXT) * ST0_EXT

ETOT_EXT =  EI_EXT
ETOT_EXT = MIN(0.999, ETOT_EXT)
KDEP_EXT = W_DEP_EXT * (IT *ETOT_EXT * Q_INSP / V_EXT + (1-IT) *ETOT_EXT * Q_EXP / V_EXT)

; ───────────────────────────────────────────────
; 6. Artificial Airway (气管导管) ─ 根据相位动态计算
; ───────────────────────────────────────────────
IF (IT .EQ. 1) THEN
RR0 = 0.5
ELSE
RR0 = 1
ENDIF
RR02 = RR0**2
RR04 = RR0**4

AREA_AA = PI * D_AA**2 / 4
U_AA    = 2*VT / (T_INSP+ T_EXP) * 0.001/AREA_AA
R_AA    = D_AA / 2
ST0_AA  = (RHOP * CC * D_P**2 * U_AA) / (36 * MU * R_AA)

COS2A_AA = COS(THETA_AA)**2   ; THETA_AA 用于 impaction，不变
SIN_A_AA = SIN(THETA_AA)

F0_AA = PI - (PI/4 + (5/4*PI - 8/3)*COS2A_AA)*RR02
F1_AA = 1 + (-1/3 + (PI - 11/3)*COS2A_AA - SIN_A_AA/3)*RR02 + ((2/3 - PI/8)*COS2A_AA + SIN_A_AA**2/5 + (6 - 15/8*PI)*COS2A_AA**2 + (7/15 - PI/8)*SIN_A_AA**2*COS2A_AA)*RR04

EI_AA = (8 * SIN(THETA_AA) * F1_AA) / (RR0 * F0_AA) * ST0_AA

; 沉降（根据相位：吸气 downhill，呼气 uphill）
; 参数计算
GAMMA_AA = ((3*VS*L_AA)/(8*U_AA*R_AA) * COS(ALPHA_AA))**(2.0/3.0) / (1 - VS/(2*U_AA)*SIN(ALPHA_AA))
SIGMA_AA = (VS/(6*U_AA)*SIN(ALPHA_AA)) / (1 - VS/(2*U_AA)*SIN(ALPHA_AA))
COND3 = PI/2 - (2*R_AA/(9*L_AA))*SQRT(2*VS/(3*U_AA))
COND4 = PI/2 - R_AA*VS/(8*U_AA*L_AA)

; 初始化 ES_AA 并使用一系列独立的 IF 语句进行赋值
ES_AA = 0 ; 初始化默认值
IF (IT .EQ. 1) THEN
; 值1: 吸气 (downhill)
ZETA1_AA = ((3*VS*L_AA)/(8*U_AA*R_AA) * COS(ALPHA_AA) / (1 + 3*VS/(4*U_AA)*SIN(ALPHA_AA))) ** (1.0/3.0)
ES_AA_INSP = - (2/PI)*ASIN(SQRT(1 - ZETA1_AA**2)) - (SQRT(1 - ZETA1_AA**2) / (PI*(1 + (VS/U_AA)*SIN(ALPHA_AA)))) * ((3*VS*L_AA)/(2*U_AA*R_AA)*COS(ALPHA_AA) - (2 + (VS/U_AA)*SIN(ALPHA_AA))*ZETA1_AA)
  ES_AA = ES_AA_INSP
ENDIF

IF (IT .EQ. 0 .AND. ALPHA_AA .LT. COND3 .AND. ALPHA_AA .GT. 0) THEN
; 值2: 呼气 (uphill), 条件1
ES_AA_EXP1 = (1/PI) * (3*SQRT(SIGMA_AA*(1-SIGMA_AA)) + ASIN(SQRT(1-SIGMA_AA)) + (1 - 9*SIGMA_AA**2)*ASIN(SQRT((1-SIGMA_AA)/(1+3*SIGMA_AA)))) - (2/PI)*(SQRT(GAMMA_AA*(1-GAMMA_AA))*(1-2*GAMMA_AA) + ASIN(SQRT((1-GAMMA_AA)/(1+3*GAMMA_AA))))
    ES_AA = ES_AA_EXP1
ENDIF

IF (IT .EQ. 0 .AND. ALPHA_AA .LT. COND4 .AND. ALPHA_AA .GT. COND3) THEN
; 值3: 呼气 (uphill), 条件2
ZETA0_AA = (R_AA*VS*SIN(ALPHA_AA)**2)/(8*U_AA*L_AA*COS(ALPHA_AA)) - (1/16)*SQRT(VS/(6*U_AA))*SIN(ALPHA_AA) + (7*L_AA)/(8*R_AA)*(1/TAN(ALPHA_AA))
ES_AA_EXP2 = (1 + 3*SIGMA_AA)/PI * SQRT(1 - ZETA0_AA**2) * (ZETA0_AA + 4*GAMMA_AA**(3.0/2.0)/SQRT(1 + 3*SIGMA_AA) - SQRT(ZETA0_AA**2 - 3*SIGMA_AA/(1 + 3*SIGMA_AA))) + (1 - 9*SIGMA_AA**2)/PI * ASIN(SQRT(1 - ZETA0_AA**2)) - (1/PI) * ASIN(SQRT((1 + 3*SIGMA_AA)*(1 - ZETA0_AA**2)))
    ES_AA = ES_AA_EXP2
ENDIF

IF (IT .EQ. 0 .AND. ALPHA_AA .GT. COND4) THEN
; 值4: 呼气 (uphill), 条件3
ES_AA_EXP3 = 0
  ES_AA = ES_AA_EXP3
ENDIF

ES_AA = MAX(0.0D0, ES_AA)
ETOT_AA = (1 - (1 - EI_AA)*(1 - ES_AA))

KDEP_AA = W_DEP_AA * (IT *ETOT_AA * Q_INSP / V_AA + (1-IT) *ETOT_AA * Q_EXP / V_AA)


; ───────────────────────────────────────────────
; 7. 呼气肢 (Expiratory Limb) ─ 三段直管 + 两弯（向下）
; ───────────────────────────────────────────────
AREA_EXP = PI * D_EXP**2 / 4
U_EXP    = VT / T_EXP * 0.001/AREA_EXP
R_EXP    = D_EXP / 2
ST0_EXP  = (RHOP * CC * D_P**2 * U_EXP) / (36 * MU * R_EXP)

; Sec1_EXP ────────────────────────────────
; Bend1_EXP
COS2A1_EXP = COS(THETA1_EXP)**2
SIN_A1_EXP = SIN(THETA1_EXP)
F0_BEND1_EXP = PI - (PI/4 + (5/4*PI - 8/3)*COS2A1_EXP)*RR02
F1_BEND1_EXP = 1 + (-1/3 + (PI - 11/3)*COS2A1_EXP - SIN_A1_EXP/3)*RR02 + ((2/3 - PI/8)*COS2A1_EXP + SIN_A1_EXP**2/5 + (6 - 15/8*PI)*COS2A1_EXP**2 + (7/15 - PI/8)*SIN_A1_EXP**2*COS2A1_EXP)*RR04
EI_BEND1_EXP = (8 * SIN(THETA1_EXP) * F1_BEND1_EXP) / (RR0 * F0_BEND1_EXP) * ST0_EXP

; Sec2_EXP ────────────────────────────────
; 沉降 Sec2_EXP (下坡 downhill)
; downhill or horizontal
  ZETA1_SEC2E = ((3*VS*L_SEC2_EXP)/(8*U_EXP*R_EXP) * COS(ALPHA2_EXP) / (1 + 3*VS/(4*U_EXP)*SIN(ALPHA2_EXP))) ** (1.0/3.0)
  ES_SEC2_EXP = - (2/PI)*ASIN(SQRT(1 - ZETA1_SEC2E**2)) - (SQRT(1 - ZETA1_SEC2E**2) / (PI*(1 + (VS/U_EXP)*SIN(ALPHA2_EXP)))) * ((3*VS*L_SEC2_EXP)/(2*U_EXP*R_EXP)*COS(ALPHA2_EXP) - (2 + (VS/U_EXP)*SIN(ALPHA2_EXP))*ZETA1_SEC2E)

  ES_SEC2_EXP = MAX(0.0D0, ES_SEC2_EXP)

; Bend2_EXP
EI_BEND2_EXP = (8 * SIN(THETA2_EXP) * F1_BEND1_EXP) / (RR0 * F0_BEND1_EXP) * ST0_EXP   ; 角度相同

; Sec3_EXP ────────────────────────────────

; 总效率
ETOT_EXP =  1 - (1-EI_BEND1_EXP)*(1-ES_SEC2_EXP)*(1-EI_BEND2_EXP)
KDEP_EXP = W_DEP_EXP *ETOT_EXP * Q_EXP / V_EXP

; ───────────────────────────────────────────────
; 肺和滤器 ─ 完全捕集（无悬浮相）
; ───────────────────────────────────────────────
; Lung: DADT(13) = IT * Q_INSP * C_AA
; Filter: DADT(16) = (1-IT) * Q_EXP * C_EXP


; ─────────────── 微分方程（保持原样，仅微调可读性） ───────────────
; Dry
DADT(1) = RT - IT * Q_INSP * C_DRY - KDEP_DRY * C_DRY * V_DRY
DADT(2) = KDEP_DRY * C_DRY * V_DRY

; Hum
DADT(3) = IT * Q_INSP * C_DRY - IT * Q_INSP * C_HUM - KDEP_HUM * C_HUM * V_HUM
DADT(4) = KDEP_HUM * C_HUM * V_HUM

; Insp
DADT(5) = IT * Q_INSP * C_HUM - IT * Q_INSP * C_INSP - KDEP_INSP * C_INSP * V_INSP
DADT(6) = KDEP_INSP * C_INSP * V_INSP

; Y (bidirectional)
DADT(7) = IT * (Q_INSP * C_INSP - Q_INSP * C_Y) + (1-IT) * (Q_EXP * C_EXT - Q_EXP * C_Y) - KDEP_Y * C_Y * V_Y
DADT(8) = KDEP_Y * C_Y * V_Y

; Ext (bidirectional)
DADT(9) = IT * (Q_INSP * C_Y - Q_INSP * C_EXT) + (1-IT) * (Q_EXP * C_AA - Q_EXP * C_EXT) - KDEP_EXT * C_EXT * V_EXT
DADT(10) = KDEP_EXT * C_EXT * V_EXT

; AA (bidirectional, but exp outflow without inflow from lung)
DADT(11) = IT * (Q_INSP * C_EXT - Q_INSP * C_AA) - (1-IT) * Q_EXP * C_AA - KDEP_AA * C_AA * V_AA
DADT(12)= KDEP_AA * C_AA * V_AA

; Lung (only deposition, no M)
DADT(13)= IT * Q_INSP * C_AA

; Exp
DADT(14)= (1-IT) * Q_EXP * C_Y - (1-IT) * Q_EXP * C_EXP - KDEP_EXP * C_EXP * V_EXP
DADT(15)= KDEP_EXP * C_EXP * V_EXP

; Filt
DADT(16)= (1-IT) * Q_EXP * C_EXP

; Reservoir dM_RES/dt = -RT (for completeness, though not strictly needed if RT direct)
DADT(17) = -RT

$THETA
; Initial guesses for THETAs
(3E-6) ;1 D_P (m)
(1000) ;2 RHOP
(5) ;3 CSOL (mg/mL)
(5) ;4 V_NEB (mL)
(0.0023) ;5 K_NEB (1/s)
(0.425) ;6 VT (L)
(20) ;7 FREQUENCY (bpm)
(2) ;8 IE (E in I:E=1:E)
(0.394) ;9 V_HUM
(0.684) ;10 V_INSP
(0.015) ;11 V_Y
(0.060) ;12 V_EXT
(0.0108);13 V_AA
(0.684) ;14 V_EXP
(0.01) ;15 F_HUM
(0.038) ;16 V_DRY
(0.9) ;17 F_DRY
(3) ;18 M_SW (mg) 切换阈值


$ERROR
Y = F
; ─────────────── 计算未收集分数和质量 ───────────────
; 限制 f_un 在 [0,1]，所有常数写成实数形式
F_UN_DRY = MIN(1.0, MAX(0.0, K_V*V_DRY + K_M))
M_UN_DRY = F_UN_DRY * A(2)

F_UN_HUM = MIN(1.0, MAX(0.0, K_V*V_HUM + K_M))
M_UN_HUM = F_UN_HUM * A(4)

F_UN_INSP = MIN(1.0, MAX(0.0, K_V*V_INSP + K_M))
M_UN_INSP = F_UN_INSP * A(6)

F_UN_Y = MIN(1.0, MAX(0.0, K_V*V_Y + K_M))
M_UN_Y = F_UN_Y * A(8)

F_UN_EXT = MIN(1.0, MAX(0.0, K_V*V_EXT + K_M))
M_UN_EXT = F_UN_EXT * A(10)

F_UN_AA = MIN(1.0, MAX(0.0, K_V*V_AA + K_M_AA)) ; 使用特殊 k_m_AA
M_UN_AA = F_UN_AA * A(12)

F_UN_EXP = MIN(1.0, MAX(0.0, K_V*V_EXP + K_M))
M_UN_EXP = F_UN_EXP * A(15)

; 总未收集质量
M_UN_TOTAL = M_UN_DRY + M_UN_HUM + M_UN_INSP + M_UN_Y + M_UN_EXT + M_UN_AA + M_UN_EXP

; ─────────────── 调整输出沉积质量（扣除未收集） ───────────────
A_DRY = A(2) - M_UN_DRY
A_HUM = A(4) - M_UN_HUM
A_INSP = A(6) - M_UN_INSP
A_Y = A(8) - M_UN_Y
A_EXT = A(10) - M_UN_EXT
A_AA = A(12) - M_UN_AA
A_EXP = A(15) - M_UN_EXP

; 未调整部分保持原样
A_LUNG = A(13)
A_FILT = A(16)
M_RES = A(17)

;$ESTIMATION METHOD=1 INTER MAXEVAL=50 PRINT=5
$SIMULATION (123456) ONLYSIM NSUB=1 ; For simulation only (no estimation)

$TABLE ID TIME
; --- 通用变量 ---
IT Q_INSP Q_EXP RT
; --- 各区段药物质量 (悬浮+沉积) ---
A_DRY A_HUM A_INSP A_Y A_EXT A_AA A_LUNG A_EXP A_FILT M_RES
; --- 新增: 未收集质量 ---
M_UN_DRY M_UN_HUM M_UN_INSP M_UN_Y M_UN_EXT M_UN_AA M_UN_EXP M_UN_TOTAL
; --- 湿化器 (Humidifier) 参数 ---
ST0_HUM F0_HUM F1_HUM EI_HUM1 EINT_HUM ETOT_HUM KDEP_HUM
; --- 吸气支 (Inspiratory) 参数 ---
ST0_INSP EI_BEND1 ES_SEC2 EI_BEND2 ETOT_INSP KDEP_INSP
; --- Y形管 (Y-piece) 参数 ---
ST0_Y EI_Y ETOT_Y KDEP_Y
; --- 延长管 (Extension) 参数 ---
ST0_EXT EI_EXT ETOT_EXT KDEP_EXT
; --- 人工气道 (Artificial Airway) 参数 ---
ST0_AA EI_AA COND3 COND4 ALPHA_AA ZETA1_AA ZETA0_AA ES_AA ETOT_AA KDEP_AA
; --- 呼气支 (Expiratory) 参数 ---
ST0_EXP EI_BEND1_EXP ZETA1_SEC2E ES_SEC2_EXP EI_BEND2_EXP ETOT_EXP KDEP_EXP
NOPRINT ONEHEADER FILE=sim_9-1.tab
