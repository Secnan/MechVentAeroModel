;; 1. Based on: Simulation_9.2
;; 2. Description: M_DEP_Neb_at_Insp_15cm
;; x1. Author: user

$SIZES DIMNEW=-5000
$PROBLEM Semi-Mechanistic Model for Aerosol Delivery in Mechanical Ventilation
; Nebulization at Inspiratory Limb 15cm from Y-piece
$DATA Simulation_9.csv IGNORE=@
$INPUT ID TIME AMT DV CMT
$SUBROUTINE ADVAN9 TOL=9

$MODEL
; Define compartments for suspended (M) and deposited (A) masses in each room.
; Compartments numbered sequentially.
COMP=(M_INSP) ; 1: Suspended in Inspiratory Limb (nebulization site)
COMP=(A_INSP) ; 2: Deposited in Inspiratory Limb
COMP=(M_Y) ; 3: Suspended in Y-piece
COMP=(A_Y) ; 4: Deposited in Y-piece
COMP=(M_EXT) ; 5: Suspended in Extension tube
COMP=(A_EXT) ; 6: Deposited in Extension tube
COMP=(M_AA) ; 7: Suspended in Artificial Airway
COMP=(A_AA) ; 8: Deposited in Artificial Airway
COMP=(A_LUNG) ; 9: Deposited in Lung (no suspended, fully captured)
COMP=(M_EXP) ; 10: Suspended in Expiratory Limb
COMP=(A_EXP) ; 11: Deposited in Expiratory Limb
COMP=(A_FILT) ; 12: Deposited in Filter (no suspended, fully captured)
COMP=(M_RES) ; 13: Reservoir mass for nebulizer (for R(t) calculation)

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
  M_SW   = THETA(15)          ; 切换阈值 (mg)
  R_CONST= K_NEB * M_RES0     ; 初期恒定速率

; ─────────────── 通气参数 ───────────────
  VT     = THETA(6)
  FREQ   = THETA(7)
  IE     = THETA(8)
  TCYCLE = 60 / FREQ
  T_INSP = TCYCLE / (1 + IE)
  T_EXP  = TCYCLE - T_INSP

; ─────────────── 各段体积 & 几何 (L & m) ───────────────
V_INSP = THETA(9)
D_INSP = 0.019
L_INSP = 0.15

V_Y = THETA(10)
D_Y = 0.019
L_Y = 0.05
THETA_Y = PI/4

V_EXT = THETA(11)
D_EXT = 0.014
L_EXT = 0.3
THETA_EXT = PI/3
ALPHA_EXT = 0

V_AA = THETA(12)
D_AA = 0.007
L_AA = 0.28
THETA_AA = PI*40/180
ALPHA_AA = PI*50/180

V_EXP = THETA(13)
D_EXP = 0.019
L_EXP = 1.8

; 呼气支分段（向下）
L_SEC1_EXP = 0.45
ALPHA1_EXP = 0
THETA1_EXP = 3*PI/4
L_SEC2_EXP = 0.9
ALPHA2_EXP = PI/4
THETA2_EXP = 3*PI/4
L_SEC3_EXP = 0.45
ALPHA3_EXP = 0

F_INSP = THETA(14) ; 
W_DEP_Y = 12
W_DEP_EXT = 22
W_DEP_AA = 80
W_DEP_EXP = 12
K_V = 0.2
K_M = 0.1
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
MRES = A(13)
IF (MRES .GT. M_SW) THEN
  RT = R_CONST
ELSE
  RT = K_NEB * MRES
ENDIF

RT = MAX(0.0D0, RT)

; ─────────────── 浓度 (mg/L) ───────────────
C_INSP = A(1)/V_INSP
C_Y = A(3)/V_Y
C_EXT = A(5)/V_EXT
C_AA = A(7)/V_AA
C_EXP = A(10)/V_EXP

; ───────────────────────────────────────────────
; 0. 公共物理量
; ───────────────────────────────────────────────
TAUP = (RHOP * D_P**2 * CC) / (18 * MU)
VS = GRAV * TAUP
RR0 = 1 ; daughter/parent 半径比
RR02 = RR0**2
RR04 = RR0**4

; ───────────────────────────────────────────────
; 1. 吸气支 (Inspiratory Limb)
; ───────────────────────────────────────────────
QM_INSP = QM_INSP
EINT_INSP = F_INSP * (RT / R_CONST)
ETOT_INSP = EINT_INSP
KDEP_INSP = ETOT_INSP * Q_INSP / V_INSP

; ───────────────────────────────────────────────
; 2. Y-piece
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
; 3. Extension tube 
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
; 4. Artificial Airway
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

; 沉降
GAMMA_AA = ((3*VS*L_AA)/(8*U_AA*R_AA) * COS(ALPHA_AA))**(2.0/3.0) / (1 - VS/(2*U_AA)*SIN(ALPHA_AA))
SIGMA_AA = (VS/(6*U_AA)*SIN(ALPHA_AA)) / (1 - VS/(2*U_AA)*SIN(ALPHA_AA))
COND3 = PI/2 - (2*R_AA/(9*L_AA))*SQRT(2*VS/(3*U_AA))
COND4 = PI/2 - R_AA*VS/(8*U_AA*L_AA)

; 初始化
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
; 5. 呼气肢 (Expiratory Limb)
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
EI_BEND2_EXP = (8 * SIN(THETA2_EXP) * F1_BEND1_EXP) / (RR0 * F0_BEND1_EXP) * ST0_EXP

; Sec3_EXP ────────────────────────────────

; 总效率
ETOT_EXP =  1 - (1-EI_BEND1_EXP)*(1-ES_SEC2_EXP)*(1-EI_BEND2_EXP)
KDEP_EXP = W_DEP_EXP *ETOT_EXP * Q_EXP / V_EXP

; ───────────────────────────────────────────────
; 肺和滤器 
; ───────────────────────────────────────────────
; Lung: DADT(9) = IT * Q_INSP * C_AA
; Filter: DADT(12) = (1-IT) * Q_EXP * C_EXP

; ─────────────── 微分方程 ───────────────
; Insp (nebulization site)
DADT(1) = RT - IT * Q_INSP * C_INSP - KDEP_INSP * C_INSP * V_INSP
DADT(2) = KDEP_INSP * C_INSP * V_INSP

; Y (bidirectional)
DADT(3) = IT * (Q_INSP * C_INSP - Q_INSP * C_Y) + (1-IT) * (Q_EXP * C_EXT - Q_EXP * C_Y) - KDEP_Y * C_Y * V_Y
DADT(4) = KDEP_Y * C_Y * V_Y

; Ext (bidirectional)
DADT(5) = IT * (Q_INSP * C_Y - Q_INSP * C_EXT) + (1-IT) * (Q_EXP * C_AA - Q_EXP * C_EXT) - KDEP_EXT * C_EXT * V_EXT
DADT(6) = KDEP_EXT * C_EXT * V_EXT

; AA (bidirectional, but exp outflow without inflow from lung)
DADT(7) = IT * (Q_INSP * C_EXT - Q_INSP * C_AA) - (1-IT) * Q_EXP * C_AA - KDEP_AA * C_AA * V_AA
DADT(8)= KDEP_AA * C_AA * V_AA

; Lung (only deposition, no M)
DADT(9)= IT * Q_INSP * C_AA

; Exp
DADT(10)= (1-IT) * Q_EXP * C_Y - (1-IT) * Q_EXP * C_EXP - KDEP_EXP * C_EXP * V_EXP
DADT(11)= KDEP_EXP * C_EXP * V_EXP

; Filt
DADT(12)= (1-IT) * Q_EXP * C_EXP
; Reservoir dM_RES/dt = -RT (for completeness, though not strictly needed if RT direct)

DADT(13) = -RT

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
(0.057) ;9 V_INSP
(0.015) ;10 V_Y
(0.060) ;11 V_EXT
(0.0108);12 V_AA
(0.684) ;13 V_EXP
(0.9) ;14 F_INSP
(3) ;15 M_SW (mg) 

$ERROR
Y = F
; ─────────────── 质量平衡 ───────────────
F_UN_INSP = MIN(1.0, MAX(0.0, K_V*V_INSP + K_M))
M_UN_INSP = F_UN_INSP * A(2)
F_UN_Y = MIN(1.0, MAX(0.0, K_V*V_Y + K_M))
M_UN_Y = F_UN_Y * A(4)
F_UN_EXT = MIN(1.0, MAX(0.0, K_V*V_EXT + K_M))
M_UN_EXT = F_UN_EXT * A(6)
F_UN_AA = MIN(1.0, MAX(0.0, K_V*V_AA + K_M_AA))
M_UN_AA = F_UN_AA * A(8)
F_UN_EXP = MIN(1.0, MAX(0.0, K_V*V_EXP + K_M))
M_UN_EXP = F_UN_EXP * A(11)
M_UN_TOTAL = M_UN_INSP + M_UN_Y + M_UN_EXT + M_UN_AA + M_UN_EXP

; ─────────────── 沉积质量校正 ───────────────
A_INSP = A(2) - M_UN_INSP
A_Y = A(4) - M_UN_Y
A_EXT = A(6) - M_UN_EXT
A_AA = A(8) - M_UN_AA
A_EXP = A(11) - M_UN_EXP
A_LUNG = A(9)
A_FILT = A(12)
M_RES = A(13)

;$ESTIMATION METHOD=1 INTER MAXEVAL=50 PRINT=5
$SIMULATION (123456) ONLYSIM NSUB=1

$TABLE ID TIME
; --- 通用变量 ---
IT Q_INSP Q_EXP RT
; --- 各区段药物质量 (悬浮+沉积) ---
A_INSP A_Y A_EXT A_AA A_LUNG A_EXP A_FILT M_RES
; --- 未收集质量 ---
M_UN_INSP M_UN_Y M_UN_EXT M_UN_AA M_UN_EXP M_UN_TOTAL
; --- 吸气支 (Inspiratory) 参数 ---
EINT_INSP ETOT_INSP KDEP_INSP
; --- Y形管 (Y-piece) 参数 ---
ST0_Y EI_Y ETOT_Y KDEP_Y
; --- 延长管 (Extension) 参数 ---
ST0_EXT EI_EXT ETOT_EXT KDEP_EXT
; --- 人工气道 (Artificial Airway) 参数 ---
ST0_AA EI_AA COND3 COND4 ALPHA_AA ZETA1_AA ZETA0_AA ES_AA ETOT_AA KDEP_AA
; --- 呼气支 (Expiratory) 参数 ---
ST0_EXP EI_BEND1_EXP ZETA1_SEC2E ES_SEC2_EXP EI_BEND2_EXP ETOT_EXP KDEP_EXP
NOPRINT ONEHEADER FILE=VMN_Pos2_ETT-1.tab
