function R = noon_charging(P_bar, Tcell_C, TBch_C, th)
% NOON_CHARGING  Noon Energy 탄소-산소 저장 시스템의 Charging mode 정상상태 모델
%                (Term Project 2026, Figure 2(a), Stream 1~13)  -- 담당자 A
%
%   R = noon_charging(P_bar, Tcell_C, TBch_C)
%   R = noon_charging(P_bar, Tcell_C, TBch_C, th)   % Cantera 객체 재사용(스윕용)
%
% 입력
%   P_bar    : 시스템 압력 [bar]                    (기본 5)
%   Tcell_C  : 셀(SOEC) 온도 [degC]                 (기본 850)
%   TBch_C   : 충전 시 Boudouard 반응기 온도 [degC] (기본 650)
%   th       : (선택) 이전 호출의 R.th. 넘겨주면 YAML 을 다시 읽지 않는다.
%
% 출력 (구조체 R)
%   R.S        : 스트림 구조체 배열 (1~13). n=[CO CO2 O2 N2 C(s)] [kmol/s], T [K], P [Pa]
%   R.streams  : 스트림 결과표 (table)                         -> 과제 (a), 역할 A1
%   R.eq       : 화학평형/재순환 상세 (Kp, 순환량 등)          -> 역할 A3
%   R.cell     : 셀 전압/전류/전력
%   R.duty     : 장치별 열부하/회수열/방출열 [W]               -> 과제 (b), 역할 A4
%   R.power    : 전력 항목별 내역 [W], 총 충전 입력 P_in
%   R.eta_ch   : 충전 효율 = nC*LHV_C / P_in                   -> 과제 (c)
%   R.balance  : 장치별/전체 원소수지, 에너지수지 잔차 (table)  -> 역할 A5
%   R.th       : Cantera 객체 묶음 (재사용용)
%
% 부호 규약
%   Q > 0 : 장치로 들어가는 열 (가열 필요),  Q < 0 : 장치에서 빼내는 열 (냉각)
%   W > 0 : 시스템으로 들어가는 전력
%   엔탈피는 Cantera 의 절대 엔탈피(생성엔탈피 포함)를 쓰므로 반응열이 자동으로 포함된다.
%
% 단위: 몰유량 kmol/s, 엔탈피 J/kmol, 열/전력 W(=J/s), 온도 K (표는 degC), 압력 Pa (표는 bar)
%
% 기준: MATLAB R2026a + Cantera 3.0.0 MATLAB toolbox (legacy 문법: Solution/set/enthalpy_mole)

    % ------------------------------------------------------------------
    % 0. 입력 기본값
    % ------------------------------------------------------------------
    if nargin < 1 || isempty(P_bar),   P_bar   = 5;   end
    if nargin < 2 || isempty(Tcell_C), Tcell_C = 850; end
    if nargin < 3 || isempty(TBch_C),  TBch_C  = 650; end
    if nargin < 4 || isempty(th),      th = load_thermo(); end

    % ------------------------------------------------------------------
    % 1. 상수 및 과제 지정값
    % ------------------------------------------------------------------
    Ru     = 8314;          % 일반기체상수 [J/(kmol K)]  (과제 지정 R = 8.314 J/mol/K)
    F      = 96485;         % Faraday 상수 [C/mol]
    p0     = 1e5;           % 표준압력 p° = 1 bar [Pa]
    T0     = 298.15;        % 주위/저장 온도 25 degC [K]
    P      = P_bar*1e5;     % 시스템 압력 [Pa]
    Tcell  = Tcell_C + 273.15;   % 셀 온도 [K]
    TB     = TBch_C  + 273.15;   % 충전 반응기 온도 [K]

    nCO2_feed = 1.0;        % Stream 1 의 CO2 공급량 [kmol/s]
    X_cell    = 0.70;       % SOEC 에서 입구 CO2 의 CO 전환율 [-]
    eff_hex   = 0.85;       % 열교환기 유용도(effectiveness) [-]
    f_link    = 0.85;       % thermal link 로 회수 가능한 donor 열의 비율 [-]
    dT_link   = 20;         % thermal link 최소 온도차 [K]
    dV_loss   = 0.10;       % 셀 과전압 가정 [V]  (V_ch = E_rev + 0.10)
    LHV_C     = 393.5e6;    % 탄소 LHV [J/kmol] (= 393.5 kJ/mol, 25 degC)

    % 성분 순서(팀 공통으로 쓰기 쉽게 N2 포함 5성분): 1 CO, 2 CO2, 3 O2, 4 N2, 5 C(s)
    iCO = 1; iCO2 = 2; iO2 = 3; iC = 5;

    % ------------------------------------------------------------------
    % 2. Boudouard 평형상수  2CO <-> C(s) + CO2  (T = TB,ch)
    % ------------------------------------------------------------------
    % 표준 Gibbs 에너지 g°_i(T) = h_i(T) - T*s°_i(T, 1 bar)  [J/kmol]  (Cantera)
    gCO_B  = g0_sp(th, iCO,  TB);
    gCO2_B = g0_sp(th, iCO2, TB);
    gC_B   = g0_sp(th, iC,   TB);      % graphite, activity = 1

    % dG°_B(T) = g°_C + g°_CO2 - 2 g°_CO   (탄소 석출 방향)
    dG_B = gC_B + gCO2_B - 2*gCO_B;                    % [J/kmol]
    % K_B(T) = exp(-dG°_B / (R T)) = (p_CO2/p°) / (p_CO/p°)^2
    K_B  = exp(-dG_B/(Ru*TB));                         % [-]

    % ------------------------------------------------------------------
    % 3. 재순환 루프의 물질수지 (정상상태, purge 없음)
    % ------------------------------------------------------------------
    % 미지수: Stream 2 의 CO 유량 a, CO2 유량 b [kmol/s], 반응기 탄소 석출량 xi [kmol/s]
    %   SOEC   : CO2 -> CO + 1/2 O2, 전환량 = X*b
    %            Stream 5 : CO = a + X*b,  CO2 = (1-X)*b,   Stream 12 : O2 = X*b/2
    %   반응기 : 2CO -> C + CO2, 진행도 xi
    %            Stream 7 : CO = a + X*b - 2*xi,  CO2 = (1-X)*b + xi,   Stream 10 : C = xi
    %   혼합기 : Stream 2 = Stream 1 + Stream 9(= Stream 7 전량 재순환)
    %            CO  : a = a + X*b - 2*xi            ->  xi = X*b/2
    %            CO2 : b = nCO2_feed + (1-X)*b + xi  ->  X*b = nCO2_feed + xi
    %   두 식을 풀면 xi = nCO2_feed, b = 2*nCO2_feed/X
    %   (시스템 전체 탄소 원자수지와 동일: 들어온 CO2 1 kmol 당 C(s) 1 kmol 석출)
    xi = nCO2_feed;                        % 탄소 석출(저장) 속도 [kmol/s]
    b  = 2*xi/X_cell;                      % Stream 2 의 CO2 유량 [kmol/s]  (X*b = 2*xi)
    c  = (1 - X_cell)*b + xi;              % Stream 7(=9) 의 CO2 유량 [kmol/s]

    % 평형 조건으로 재순환 CO 유량 a 결정 (Stream 7: CO = a, CO2 = c, 전체 = a + c)
    %   K_B = (y_CO2 P/p°) / (y_CO P/p°)^2 = c (a + c) / (a^2 (P/p°))
    %   ->  (K_B P/p°) a^2 - c a - c^2 = 0   (a > 0 인 근)
    KP = K_B*(P/p0);                                   % K_B * (P/p°) [-]
    a  = c*(1 + sqrt(1 + 4*KP))/(2*KP);                % 재순환 CO 유량 [kmol/s]

    % 평형 잔차 확인용: ln K_B - ln[(p_CO2/p°)/(p_CO/p°)^2]
    yCO_7  = a/(a + c);                                % Stream 7 CO 몰분율
    yCO2_7 = c/(a + c);                                % Stream 7 CO2 몰분율
    res_eq = log(K_B) - log((yCO2_7*P/p0)/(yCO_7*P/p0)^2);

    % ------------------------------------------------------------------
    % 4. 스트림 조성 [CO CO2 O2 N2 C] (kmol/s)
    % ------------------------------------------------------------------
    n1  = [0, nCO2_feed, 0, 0, 0];                     % CO2 탱크 -> 혼합기
    n9  = [a, c, 0, 0, 0];                             % 재순환 가스 (25 degC)
    n2  = n1 + n9;                                     % 혼합기 출구  (CO = a, CO2 = b)
    n5  = [a + X_cell*b, (1 - X_cell)*b, 0, 0, 0];     % SOEC 연료극 출구
    n12 = [0, 0, X_cell*b/2, 0, 0];                    % SOEC 산소극 출구 (순수 O2)
    n7  = n9;                                          % 반응기 출구 가스 (전량 재순환)
    n10 = [0, 0, 0, 0, xi];                            % 석출 탄소

    % ------------------------------------------------------------------
    % 5. 에너지수지: 혼합기 -> HEX1 -> H1 -> SOEC -> Trim -> 반응기 -> 냉각기
    % ------------------------------------------------------------------
    % (1) 단열 혼합기: H2(T2) = H1(T0) + H9(T0)  -> T2 (이상기체, 동일 온도이므로 25 degC)
    H1s = Hdot(th, n1, T0);                            % Stream 1 엔탈피 유량 [W]
    H9s = Hdot(th, n9, T0);                            % Stream 9 엔탈피 유량 [W]
    T2  = T_from_H(th, n2, H1s + H9s, T0, T0);         % 혼합기 출구 온도 [K]
    H2s = Hdot(th, n2, T2);

    % (2) HEX1 (단열 향류, 유용도 0.85): 고온측 Stream 7 -> 8, 저온측 Stream 2 -> 3
    %     Q_max = min[ H_hot(T_h,in) - H_hot(T_c,in),  H_cold(T_h,in) - H_cold(T_c,in) ]
    %     Q     = eff * Q_max,  출구온도는 각 유체의 에너지수지로 결정
    hx1 = hex_counterflow(th, n7, TB, n2, T2, eff_hex);
    T8  = hx1.Th_out;                                  % Stream 8 온도 [K]
    T3  = hx1.Tc_out;                                  % Stream 3 온도 [K]

    % (3) Heater H1: Stream 3 -> 4 (셀 온도까지 전기 가열)
    H3s  = Hdot(th, n2, T3);
    H4s  = Hdot(th, n2, Tcell);
    Q_H1 = H4s - H3s;                                  % H1 열부하 [W] (>0)

    % (4) SOEC (등온, Tcell): 가역전압, 운전전압, 전류, 전력
    %     E°(T) = -dG°_ox(T)/(2F),  dG°_ox = g°_CO2 - g°_CO - 1/2 g°_O2  (CO + 1/2 O2 -> CO2)
    dG_ox = g0_sp(th, iCO2, Tcell) - g0_sp(th, iCO, Tcell) - 0.5*g0_sp(th, iO2, Tcell); % [J/kmol]
    E0    = -(dG_ox/1000)/(2*F);                       % 표준 가역전압 [V]  (J/kmol -> J/mol)

    %     Nernst: E_rev = E° + (RT/2F) ln[(p_CO/p°)(p_O2/p°)^(1/2)/(p_CO2/p°)]
    %     충전 시 조성은 Stream 5(연료극)와 Stream 12(산소극, 순수 O2 -> p_O2 = P) 로 평가
    yCO_5  = n5(iCO)/sum(n5);                          % Stream 5 CO 몰분율
    yCO2_5 = n5(iCO2)/sum(n5);                         % Stream 5 CO2 몰분율
    pCO  = yCO_5*P/p0;                                 % 무차원 분압 p_CO/p°
    pCO2 = yCO2_5*P/p0;                                % 무차원 분압 p_CO2/p°
    pO2  = 1.0*P/p0;                                   % 무차원 분압 p_O2/p° (Stream 12)
    E_rev = E0 + (Ru/1000)*Tcell/(2*F)*log(pCO*sqrt(pO2)/pCO2);   % 가역전압 [V]

    V_ch   = E_rev + dV_loss;                          % 충전 운전전압 V_ch = E_rev + 0.10 [V]
    nCO_conv = X_cell*b;                               % 셀에서 전환된 CO2(=생성 CO) [kmol/s]
    I_cell = 2*F*(nCO_conv*1000);                      % 전류 I = 2F * n_CO [A]  (kmol/s -> mol/s)
    W_cell = I_cell*V_ch;                              % 셀 전력 = I*V [W] (입력)

    %     등온 셀 에너지수지: H4 + W_cell + Q_cell = H5 + H12  ->  Q_cell
    H5s  = Hdot(th, n5,  Tcell);
    H12s = Hdot(th, n12, Tcell);
    Q_cell = (H5s + H12s) - H4s - W_cell;              % 셀 열부하 [W] (>0 이면 흡열 운전)
    V_tn   = (H5s + H12s - H4s)/I_cell;                % 열중립 전압 V_tn = dH/(2F) [V]

    % (5) Trim unit: Stream 5 (Tcell) -> Stream 6 (TB,ch), 조성 동결
    H6s    = Hdot(th, n5, TB);
    Q_trim = H6s - H5s;                                % Trim 열부하 [W] (Tcell > TB 이면 <0)

    % (6) Boudouard 반응기 + 고체 분리기 (등온, TB,ch): Stream 6 -> 7(가스) + 10(탄소)
    H7s  = Hdot(th, n7,  TB);
    H10s = Hdot(th, n10, TB);
    Q_rx = (H7s + H10s) - H6s;                         % 반응기 열부하 [W] (석출은 발열 -> <0)

    % (7) 냉각기들 (모두 25 degC 까지 냉각, 열은 방출)
    H8s  = Hdot(th, n7,  T8);
    Q_c89 = H9s - H8s;                                 % 재순환 냉각기 8 -> 9 [W] (<0)
    H11s = Hdot(th, n10, T0);
    Q_cC  = H11s - H10s;                               % 탄소 냉각기 10 -> 11 [W] (<0)
    H13s = Hdot(th, n12, T0);
    Q_cO2 = H13s - H12s;                               % O2 냉각기 12 -> 13 [W] (<0)

    % ------------------------------------------------------------------
    % 6. Thermal link (셀 <-> Boudouard 반응기, 같은 모드 안에서만)
    % ------------------------------------------------------------------
    % 음(-)의 열부하 장치(donor)에서 양(+)의 열부하 장치로, donor 온도가 20 K 이상 높을 때만
    %   Q_link = min( 0.85*|Q_donor| , Q_receiver )
    [Q_link, donor] = thermal_link(Q_cell, Tcell, Q_rx, TB, f_link, dT_link);
    link_to_cell = strcmp(donor, 'reactor') * Q_link;  % 반응기 -> 셀 로 전달된 열 [W]
    link_to_rx   = strcmp(donor, 'cell')    * Q_link;  % 셀 -> 반응기 로 전달된 열 [W]

    % ------------------------------------------------------------------
    % 7. 전기 가열(저항 가열 100 %) 및 방출열 집계
    % ------------------------------------------------------------------
    % 회수 후 남은 양(+)의 열부하는 전부 전기 가열
    Qel_H1   = max(Q_H1, 0);                           % Heater H1 [W]
    Qel_cell = max(Q_cell, 0) - link_to_cell;          % 셀 보조 가열 [W]
    Qel_trim = max(Q_trim, 0);                         % Trim unit 가열 [W]
    Qel_rx   = max(Q_rx, 0) - link_to_rx;              % 반응기 가열 [W]
    Qel_tot  = Qel_H1 + Qel_cell + Qel_trim + Qel_rx;  % 전기 가열 합계 [W]

    % 회수되지 않은 음(-)의 열부하는 전부 방출 (양수로 표기)
    Qrej_cell = max(-Q_cell, 0) - link_to_rx;          % 셀 방출열 [W]
    Qrej_trim = max(-Q_trim, 0);                       % Trim unit 방출열 [W]
    Qrej_rx   = max(-Q_rx, 0) - link_to_cell;          % 반응기 방출열 [W]
    Qrej_c89  = -Q_c89;                                % 재순환 냉각기 방출열 [W]
    Qrej_cC   = -Q_cC;                                 % 탄소 냉각기 방출열 [W]
    Qrej_cO2  = -Q_cO2;                                % O2 냉각기 방출열 [W]
    Qrej_tot  = Qrej_cell + Qrej_trim + Qrej_rx + Qrej_c89 + Qrej_cC + Qrej_cO2;

    % 충전 입력 전력: P_in = W_cell + (전기 가열 합계)
    %  (충전 모드에는 공기 압축기와 CO2 분리기가 없다. O2 는 밸브로 방출, 일 회수 없음)
    P_in = W_cell + Qel_tot;                           % [W]

    % 충전 효율: eta_ch = N_C * LHV_C / E_in = (xi * LHV_C) / P_in
    eta_ch = xi*LHV_C/P_in;                            % [-]

    % ------------------------------------------------------------------
    % 8. 스트림 구조체 및 결과표
    % ------------------------------------------------------------------
    S = repmat(struct('n', zeros(1,5), 'T', NaN, 'P', P, 'name', ""), 13, 1);
    S(1)  = mk(n1,  T0,    P, "CO2 tank -> mixer");
    S(2)  = mk(n2,  T2,    P, "mixer -> HEX1 (cold in)");
    S(3)  = mk(n2,  T3,    P, "HEX1 -> heater H1");
    S(4)  = mk(n2,  Tcell, P, "H1 -> SOEC fuel inlet");
    S(5)  = mk(n5,  Tcell, P, "SOEC fuel outlet");
    S(6)  = mk(n5,  TB,    P, "trim unit -> carbon reactor");
    S(7)  = mk(n7,  TB,    P, "reactor gas -> HEX1 (hot in)");
    S(8)  = mk(n7,  T8,    P, "HEX1 -> recycle cooler");
    S(9)  = mk(n9,  T0,    P, "recycle -> mixer");
    S(10) = mk(n10, TB,    P, "deposited carbon");
    S(11) = mk(n10, T0,    P, "carbon -> storage");
    S(12) = mk(n12, Tcell, P, "SOEC O2 outlet");
    S(13) = mk(n12, T0,    P, "O2 -> vent valve");

    R.S       = S;
    R.streams = stream_table(S, th.MW);

    % ------------------------------------------------------------------
    % 9. 원소수지 / 에너지수지 잔차 (장치별 + 전체 시스템)
    % ------------------------------------------------------------------
    % 장치 정의: {이름, 입구 스트림, 출구 스트림, 열 Q [W], 전력 W [W]}
    %   에너지 잔차 = sum(H_in) + Q + W - sum(H_out)
    %   원소 잔차   = sum(원자 in) - sum(원자 out)   [kmol/s]
    U = { ...
        'Mixer',            [1 9],  2,        0,         0; ...
        'HEX 1',            [2 7],  [3 8],    0,         0; ...
        'Heater H1',        3,      4,        Q_H1,      0; ...
        'rSOC (SOEC)',      4,      [5 12],   Q_cell,    W_cell; ...
        'Trim unit',        5,      6,        Q_trim,    0; ...
        'Reactor/separator',6,      [7 10],   Q_rx,      0; ...
        'Recycle cooler',   8,      9,        Q_c89,     0; ...
        'C cooler',         10,     11,       Q_cC,      0; ...
        'O2 cooler',        12,     13,       Q_cO2,     0; ...
        'SYSTEM (overall)', 1,      [11 13],  -Qrej_tot, P_in};
    R.balance = balance_table(th, S, U);

    % ------------------------------------------------------------------
    % 10. 결과 정리
    % ------------------------------------------------------------------
    R.in = struct('P_bar', P_bar, 'Tcell_C', Tcell_C, 'TBch_C', TBch_C);

    R.eq = struct( ...
        'dG_B',        dG_B/1000, ...        % Boudouard(석출) 표준 Gibbs 변화 [J/mol]
        'K_B',         K_B, ...              % 평형상수 [-]
        'yCO_eq',      yCO_7, ...            % 반응기 출구 CO 몰분율
        'yCO2_eq',     yCO2_7, ...           % 반응기 출구 CO2 몰분율
        'res_lnK',     res_eq, ...           % 평형 잔차 ln K - ln Q
        'xi',          xi, ...               % 탄소 석출 속도 [kmol/s]
        'n_recycle',   a + c, ...            % 재순환 유량 Stream 9 [kmol/s]
        'n_cell_in',   a + b, ...            % 셀 입구 유량 Stream 4 [kmol/s]
        'recycle_ratio', (a + c)/nCO2_feed, ...   % 재순환비 n9/n1 [-]
        'CO_conv_reactor', 2*xi/n5(iCO));    % 반응기 1회 통과 CO 전환율 [-]

    R.cell = struct('E0', E0, 'E_rev', E_rev, 'V', V_ch, 'V_tn', V_tn, ...
                    'I', I_cell, 'W', W_cell, 'yCO_5', yCO_5, 'yCO2_5', yCO2_5);

    R.duty = struct( ...
        'Q_H1', Q_H1, 'Q_cell', Q_cell, 'Q_trim', Q_trim, 'Q_rx', Q_rx, ...  % 장치 열부하
        'Q_c89', Q_c89, 'Q_cC', Q_cC, 'Q_cO2', Q_cO2, ...                    % 냉각기 열부하
        'Q_HEX1', hx1.Q, 'Q_HEX1_max', hx1.Qmax, 'dTmin_HEX1', hx1.dTmin, ... % 회수열(recuperator)
        'Q_link', Q_link, 'link_donor', donor, ...                           % thermal link
        'Qrej_cell', Qrej_cell, 'Qrej_trim', Qrej_trim, 'Qrej_rx', Qrej_rx, ...
        'Qrej_c89', Qrej_c89, 'Qrej_cC', Qrej_cC, 'Qrej_cO2', Qrej_cO2, ...
        'Qrej_tot', Qrej_tot);

    R.power = struct('W_cell', W_cell, 'Qel_H1', Qel_H1, 'Qel_cell', Qel_cell, ...
                     'Qel_trim', Qel_trim, 'Qel_rx', Qel_rx, 'Qel_tot', Qel_tot, ...
                     'P_in', P_in);

    R.eta_ch  = eta_ch;
    R.NC_rate = xi;          % 탄소 저장 속도 [kmol/s] (방전측 소비 속도와 매칭용)
    R.LHV_C   = LHV_C;
    R.th      = th;
end

% ======================================================================
%  로컬 함수
% ======================================================================

function th = load_thermo()
% Cantera 물성 객체 생성 (YAML 은 한 번만 읽는다)
%   가스: gri30.yaml 의 CO, CO2, O2, N2 (NASA 7계수 다항식, 이상기체)
%   고체: graphite.yaml 의 C(gr)
    th.gas  = Solution('gri30.yaml');          % 기체 물성 (GRI-Mech 3.0 열역학 데이터)
    try
        th.gr = Solution('graphite.yaml');     % 흑연 물성 (과제 지정 파일)
    catch
        th.gr = ThermoPhase('graphite.yaml');  % 반응 메커니즘이 없는 상이라 Solution 이 실패할 때
    end
    th.p0   = 1e5;                             % 표준압력 1 bar [Pa]
    th.Xstr = {'CO:1', 'CO2:1', 'O2:1', 'N2:1'};   % 순수 성분 조성 문자열
    mwg = molecularWeights(th.gas);            % 기체 분자량 [kg/kmol]
    mwc = molecularWeights(th.gr);             % 흑연 분자량 [kg/kmol]
    th.MW = [mwg(speciesIndex(th.gas, 'CO')), mwg(speciesIndex(th.gas, 'CO2')), ...
             mwg(speciesIndex(th.gas, 'O2')), mwg(speciesIndex(th.gas, 'N2')), mwc(1)];
end

function h = h_sp(th, i, T)
% 순수 성분 i 의 몰 엔탈피 h_i(T) [J/kmol] (생성엔탈피 포함, 이상기체이므로 압력 무관)
    if i == 5
        set(th.gr, 'T', T, 'P', th.p0);
        h = enthalpy_mole(th.gr);
    else
        set(th.gas, 'T', T, 'P', th.p0, 'X', th.Xstr{i});
        h = enthalpy_mole(th.gas);
    end
end

function s = s0_sp(th, i, T)
% 순수 성분 i 의 표준 몰 엔트로피 s°_i(T, 1 bar) [J/(kmol K)]
    if i == 5
        set(th.gr, 'T', T, 'P', th.p0);
        s = entropy_mole(th.gr);
    else
        set(th.gas, 'T', T, 'P', th.p0, 'X', th.Xstr{i});
        s = entropy_mole(th.gas);
    end
end

function g = g0_sp(th, i, T)
% 표준 몰 Gibbs 에너지 g°_i(T) = h_i(T) - T s°_i(T) [J/kmol]
    g = h_sp(th, i, T) - T*s0_sp(th, i, T);
end

function H = Hdot(th, n, T)
% 스트림 엔탈피 유량 H = sum_i n_i h_i(T) [W]   (n [kmol/s] * h [J/kmol])
    H = 0;
    for i = 1:5
        if n(i) ~= 0
            H = H + n(i)*h_sp(th, i, T);
        end
    end
end

function T = T_from_H(th, n, Htarget, Tlo, Thi)
% 에너지수지로 온도 결정: Hdot(n, T) = Htarget 를 만족하는 T [K] (fzero, 구간 [Tlo-1, Thi+1])
    f = @(T) (Hdot(th, n, T) - Htarget)/1e6;   % [MW] 단위로 스케일
    T = fzero(f, [Tlo - 1, Thi + 1]);
end

function hx = hex_counterflow(th, nh, Th_in, nc, Tc_in, eff)
% 단열 향류 열교환기 (유용도 모델)
%   nh, Th_in : 고온측 조성 [kmol/s], 입구 온도 [K]
%   nc, Tc_in : 저온측 조성 [kmol/s], 입구 온도 [K]
%   hx.Q [W], hx.Qmax [W], hx.Th_out [K], hx.Tc_out [K], hx.dTmin [K]
    if Th_in <= Tc_in          % 열이 흐를 수 없는 경우
        hx = struct('Q', 0, 'Qmax', 0, 'Th_out', Th_in, 'Tc_out', Tc_in, 'dTmin', Th_in - Tc_in);
        return
    end
    Hh_in = Hdot(th, nh, Th_in);
    Hc_in = Hdot(th, nc, Tc_in);
    % 최대 전달 가능 열: 고온측이 T_c,in 까지 식는 경우와 저온측이 T_h,in 까지 데워지는 경우 중 작은 값
    Qmax_h = Hh_in - Hdot(th, nh, Tc_in);
    Qmax_c = Hdot(th, nc, Th_in) - Hc_in;
    Qmax   = min(Qmax_h, Qmax_c);
    Q      = eff*Qmax;                          % 실제 전달열 Q = eff * Q_max
    % 출구 온도: 각 유체의 에너지수지
    Th_out = T_from_H(th, nh, Hh_in - Q, Tc_in, Th_in);
    Tc_out = T_from_H(th, nc, Hc_in + Q, Tc_in, Th_in);
    % 내부 온도분포 점검 (고온 -> 저온 방향 열전달이 전 구간에서 성립하는지)
    %   저온 입구 끝에서부터 누적 전달열 q 에 대해 두 유체의 온도를 비교
    Np = 20;
    dT = zeros(Np + 1, 1);
    for j = 0:Np
        q  = Q*j/Np;
        Tc = T_from_H(th, nc, Hc_in + q,       Tc_in, Th_in);
        Th = T_from_H(th, nh, Hh_in - Q + q,   Tc_in, Th_in);
        dT(j + 1) = Th - Tc;
    end
    hx = struct('Q', Q, 'Qmax', Qmax, 'Th_out', Th_out, 'Tc_out', Tc_out, 'dTmin', min(dT));
end

function [Q_link, donor] = thermal_link(Q_cell, T_cell, Q_rx, T_rx, frac, dTmin)
% 셀-반응기 thermal link
%   donor(음의 열부하) 온도가 receiver(양의 열부하) 보다 dTmin 이상 높을 때만
%   Q_link = min(frac*|Q_donor|, Q_receiver),  그 외에는 0
    Q_link = 0;
    donor  = 'none';
    if Q_rx < 0 && Q_cell > 0 && (T_rx - T_cell) >= dTmin
        Q_link = min(frac*(-Q_rx), Q_cell);
        donor  = 'reactor';
    elseif Q_cell < 0 && Q_rx > 0 && (T_cell - T_rx) >= dTmin
        Q_link = min(frac*(-Q_cell), Q_rx);
        donor  = 'cell';
    end
end

function s = mk(n, T, P, name)
% 스트림 구조체 생성
    s = struct('n', n, 'T', T, 'P', P, 'name', name);
end

function tab = stream_table(S, MW)
% 스트림 결과표: 번호, 설명, T [degC], P [bar], 몰유량 [kmol/s], 질량유량 [kg/s], 몰분율
    ns = numel(S);
    No = (1:ns)';  Desc = strings(ns, 1);
    T_C = zeros(ns, 1);  P_bar = zeros(ns, 1);
    n_kmol_s = zeros(ns, 1);  m_kg_s = zeros(ns, 1);  X = zeros(ns, 5);
    for k = 1:ns
        Desc(k)     = S(k).name;
        T_C(k)      = S(k).T - 273.15;
        P_bar(k)    = S(k).P/1e5;
        n_kmol_s(k) = sum(S(k).n);              % 전체 몰유량
        m_kg_s(k)   = sum(S(k).n.*MW);          % 질량유량 = sum n_i M_i
        X(k, :)     = S(k).n/n_kmol_s(k);       % 몰분율
    end
    x_CO = X(:,1); x_CO2 = X(:,2); x_O2 = X(:,3); x_N2 = X(:,4); x_C = X(:,5);
    tab = table(No, Desc, T_C, P_bar, n_kmol_s, m_kg_s, x_CO, x_CO2, x_O2, x_N2, x_C);
end

function tab = balance_table(th, S, U)
% 장치별 원소(C, O, N) 수지와 에너지수지 잔차
    A = [1 1 0 0 1;      % C 원자 수 / 분자  [CO CO2 O2 N2 C]
         1 2 2 0 0;      % O 원자 수
         0 0 0 2 0];     % N 원자 수
    nu = size(U, 1);
    Unit = strings(nu, 1);
    res_C = zeros(nu, 1); res_O = zeros(nu, 1); res_N = zeros(nu, 1);
    res_E_W = zeros(nu, 1); Q_W = zeros(nu, 1); W_W = zeros(nu, 1);
    for u = 1:nu
        Unit(u) = U{u, 1};
        atoms = zeros(3, 1);  Hsum = 0;
        for k = U{u, 2}                          % 입구 스트림
            atoms = atoms + A*S(k).n';
            Hsum  = Hsum + Hdot(th, S(k).n, S(k).T);
        end
        for k = U{u, 3}                          % 출구 스트림
            atoms = atoms - A*S(k).n';
            Hsum  = Hsum - Hdot(th, S(k).n, S(k).T);
        end
        Q_W(u) = U{u, 4};  W_W(u) = U{u, 5};
        res_C(u) = atoms(1); res_O(u) = atoms(2); res_N(u) = atoms(3);   % [kmol/s]
        res_E_W(u) = Hsum + Q_W(u) + W_W(u);                             % [W]
    end
    tab = table(Unit, Q_W, W_W, res_C, res_O, res_N, res_E_W);
end
