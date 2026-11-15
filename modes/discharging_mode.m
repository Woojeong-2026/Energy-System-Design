function R = discharging_mode(cfg, th)
% DISCHARGING_MODE  방전 모드 정상상태 모델 (Figure 2(b), Stream 14~29)
%
%   R = discharging_mode(cfg, th)
%
%   cfg : noon_config() 설정 구조체,  th : thermo_init() 물성 구조체
%
% [통합 인터페이스] 출력 필드는 charging_mode 와 같은 형식이다.
%   R.S, R.streams, R.balance, R.units, R.heat, R.hex, R.link, R.eq, R.cell, R.elec
%   R.P_out   : 방전 순출력 전력 [W],  R.nC : 탄소 소비 속도 [kmol/s]
%   R.eta_dis : 방전 효율 [-]
%   R.Q_exhaust : 배기(Stream 29)가 25 degC 기준으로 가지고 나가는 현열 [W]
%
% 이 파일은 팀원 2(방전 담당) 모듈 자리다. 위 필드 형식만 지키면 내부 구현을 바꿔도
% main.m 은 그대로 동작한다.
%
% 부호 규약: Q > 0 장치로 들어가는 열, W > 0 장치로 들어가는 전력
%            (SOFC 는 전력을 내보내므로 수지표에서 W = -W_cell)

    % ------------------------------------------------------------------
    % 1. 운전 조건
    % ------------------------------------------------------------------
    P     = cfg.P_bar*1e5;               % 시스템 압력 [Pa]
    Pamb  = cfg.P0_bar*1e5;              % 주위 압력 [Pa]
    T0    = cfg.T0_C + 273.15;           % 주위/저장 온도 [K]
    Tcell = cfg.Tcell_C + 273.15;        % 셀 온도 [K]
    TB    = cfg.TBdis_C + 273.15;        % 방전 Boudouard 반응기 온도 [K]
    Uf    = cfg.U_sofc;                  % SOFC CO 산화율 [-]
    xi    = cfg.nC_dis;                  % 탄소 소비 속도 = CO2 분리량 [kmol/s]

    % ------------------------------------------------------------------
    % 2. Boudouard 평형상수 (T = TB,dis). 석출 방향 정의를 그대로 쓴다.
    % ------------------------------------------------------------------
    % K_B = (p_CO2/p°)/(p_CO/p°)^2.  가스화(C + CO2 -> 2CO)는 그 역반응이고 조성 조건은 같다.
    [K_B, dG_B] = K_boudouard(th, TB, cfg.Ru);

    % ------------------------------------------------------------------
    % 3. 재순환 루프 물질수지 (정상상태) - 해석해
    % ------------------------------------------------------------------
    % 미지수: Stream 15 의 CO 유량 a, CO2 유량 b [kmol/s]
    %   반응기 (C + CO2 -> 2CO, 진행도 xi)
    %            Stream 18 : CO = a + 2 xi = m,  CO2 = b - xi = c
    %   SOFC   (CO + 1/2 O2 -> CO2, 산화량 Uf m)
    %            Stream 20 : CO = (1-Uf) m,  CO2 = c + Uf m
    %   분리기 (순수 CO2 xi 를 Stream 23 으로 제거, 나머지는 Stream 15 로 전량 재순환)
    %            CO  수지 : a = (1-Uf) m = (1-Uf)(a + 2 xi)   ->  m = 2 xi / Uf
    %            CO2 수지 : b = c + Uf m - xi = c + xi         (항등식, b 는 평형으로 정해진다)
    m = 2*xi/Uf;                         % Stream 18 CO 유량 [kmol/s]
    a = (1 - Uf)*m;                      % Stream 15 CO 유량 [kmol/s]

    % 반응기 출구(Stream 18: CO = m, CO2 = c)의 평형 (탄소가 충분히 존재, a_C = 1)
    %   K_B = (y_CO2 P/p°)/(y_CO P/p°)^2 = c (m + c) / (m^2 P/p°)
    %   ->  c^2 + m c - (K_B P/p°) m^2 = 0  의 양의 근
    KP = K_B*(P/cfg.p0);                 % K_B * (P/p°) [-]
    c  = m*(-1 + sqrt(1 + 4*KP))/2;      % Stream 18 CO2 유량 [kmol/s]
    b  = c + xi;                         % Stream 15 CO2 유량 [kmol/s]

    yCO_18  = m/(m + c);                 % Stream 18 CO 몰분율 (평형)
    yCO2_18 = c/(m + c);                 % Stream 18 CO2 몰분율 (평형)
    % 평형 잔차: ln K_B - ln[(p_CO2/p°)/(p_CO/p°)^2]  (0 이어야 한다)
    res_eq = log(K_B) - log((yCO2_18*P/cfg.p0)/(yCO_18*P/cfg.p0)^2);

    % 공기: 양론 O2 = (산화된 CO)/2 = Uf m / 2 = xi,  공급 O2 = lambda * 양론 O2
    nO2_st = Uf*m/2;                     % 양론 O2 [kmol/s]
    nO2_in = cfg.lambda_O2*nO2_st;       % 공급 O2 [kmol/s] (3배 = 200 % 과잉공기)
    nN2    = nO2_in*cfg.xN2_air/cfg.xO2_air;   % 동반 N2 [kmol/s]

    % ------------------------------------------------------------------
    % 4. 스트림 조성 [CO CO2 O2 N2 C(s)] [kmol/s]
    % ------------------------------------------------------------------
    n14  = [0, 0, 0, 0, xi];             % 저장 탄소
    n15  = [a, b, 0, 0, 0];              % 재순환 가스
    n18  = [m, c, 0, 0, 0];              % 반응기 출구 (평형)
    n20  = [(1 - Uf)*m, c + Uf*m, 0, 0, 0];    % SOFC 연료극 출구
    n23  = [0, xi, 0, 0, 0];             % 분리된 순수 CO2 -> 탱크
    nair = [0, 0, nO2_in, nN2, 0];       % 공급 공기 (Stream 24~27)
    n28  = [0, 0, nO2_in - nO2_st, nN2, 0];    % O2 가 소모된 공기 (Stream 28, 29)

    % ------------------------------------------------------------------
    % 5. 에너지수지
    % ------------------------------------------------------------------
    % (1) HEX 3 (단열 향류): 고온측 Stream 20 -> 21, 저온측 Stream 15 -> 16
    hx3 = hex_counterflow(th, n20, Tcell, n15, T0, cfg.eff_hex);
    T21 = hx3.Th_out;                    % Stream 21 온도 [K]
    T16 = hx3.Tc_out;                    % Stream 16 온도 [K]

    % (2) Heater H2: Stream 16 -> 17 (반응기 온도까지),  Q_H2 = H17 - H16
    H16  = stream_H(th, n15, T16);
    H17  = stream_H(th, n15, TB);
    Q_H2 = H17 - H16;

    % (3) 가스화 반응기 (등온, TB,dis): Stream 14 (탄소, 25 degC) + 17 -> 18
    %     Q_rx = H18 - H14 - H17   (가스화는 흡열 + 탄소 가열 -> Q_rx > 0)
    H14  = stream_H(th, n14, T0);
    H18  = stream_H(th, n18, TB);
    Q_rx = H18 - H14 - H17;

    % (4) Trim unit: Stream 18 (TB,dis) -> 19 (Tcell), 조성 동결,  Q_trim = H19 - H18
    H19    = stream_H(th, n18, Tcell);
    Q_trim = H19 - H18;

    % (5) 공기 압축기: 25 degC 가역 등온 압축 (Stream 24 -> 25)
    %     W_comp = n_air R T0 ln(P/P_amb),  등온 이상기체이므로 같은 크기의 열을 방출
    n_air  = nO2_in + nN2;               % 공기 몰유량 [kmol/s]
    W_comp = n_air*cfg.Ru*T0*log(P/Pamb);
    Q_comp = -W_comp;

    % (6) HEX 4 (단열 향류): 고온측 Stream 28 -> 29, 저온측 Stream 25 -> 26
    hx4 = hex_counterflow(th, n28, Tcell, nair, T0, cfg.eff_hex);
    T29 = hx4.Th_out;                    % Stream 29 온도 [K]
    T26 = hx4.Tc_out;                    % Stream 26 온도 [K]

    % (7) Heater H3: Stream 26 -> 27 (셀 온도까지),  Q_H3 = H27 - H26
    H26  = stream_H(th, nair, T26);
    H27  = stream_H(th, nair, Tcell);
    Q_H3 = H27 - H26;

    % (8) SOFC (등온, Tcell)
    %     가역전압은 과제 지정 스트림 조성으로 평가: 연료극 Stream 20, 산소극 Stream 28
    yCO_20  = n20(1)/sum(n20);           % Stream 20 CO 몰분율
    yCO2_20 = n20(2)/sum(n20);           % Stream 20 CO2 몰분율
    yO2_28  = n28(3)/sum(n28);           % Stream 28 O2 몰분율
    [E_rev, E0] = E_cell_CO(th, Tcell, yCO_20*P/cfg.p0, yCO2_20*P/cfg.p0, yO2_28*P/cfg.p0, ...
                            cfg.Ru, cfg.F);
    V_dis  = E_rev - cfg.dV_loss;        % 방전 운전전압 V_dis = E_rev - 0.10 [V]
    nCO_el = Uf*m;                       % 셀에서 산화된 CO [kmol/s]
    I_cell = 2*cfg.F*(nCO_el*1000);      % 전류 I = 2 F n_CO [A]
    W_cell = I_cell*V_dis;               % 셀 출력 전력 = I V [W] (양수 = 생산)

    %     등온 셀 에너지수지: H19 + H27 + Q_cell = H20 + H28 + W_cell
    H20    = stream_H(th, n20, Tcell);
    H28    = stream_H(th, n28, Tcell);
    Q_cell = (H20 + H28 + W_cell) - (H19 + H27);   % 셀 열부하 [W] (발열 -> Q_cell < 0)
    V_tn   = (H19 + H27 - H20 - H28)/I_cell;       % 열중립 전압 V_tn = -dH/(2F) [V]

    % (9) 냉각기: Stream 21 -> 22 (25 degC)
    H21   = stream_H(th, n20, T21);
    H22   = stream_H(th, n20, T0);
    Q_c22 = H22 - H21;

    % (10) CO2 분리기: 25 degC 등온 이상적 전기화학 펌프 (Faraday 효율 100 %, 전압 손실 없음)
    %      V_sep,rev = (R T_sep / 2F) ln(p_CO2,23 / p_CO2,15) = (R T_sep / 2F) ln(1 / x_CO2,15)
    %      I_sep = 2 F n_CO2,23
    %      W_sep = I_sep V_sep,rev = n_CO2,23 R T_sep ln(1 / x_CO2,15)
    %      등온 유지를 위해 W_sep 만큼의 열을 방출
    xCO2_15 = b/(a + b);                 % Stream 15 CO2 몰분율 (혼합측은 완전혼합, Stream 15 조성)
    V_sep   = (cfg.Ru/1000)*T0/(2*cfg.F)*log(1/xCO2_15);   % 분리기 가역전압 [V]
    I_sep   = 2*cfg.F*(xi*1000);         % 분리기 전류 [A]
    W_sep   = I_sep*V_sep;               % 분리기 입력 전력 [W]
    Q_sep   = -W_sep;

    % ------------------------------------------------------------------
    % 6. thermal link 적용, 전기 가열/방출열 정산
    % ------------------------------------------------------------------
    units = [ ...
        mkunit("Heater H2",              Q_H2,   NaN); ...
        mkunit("Carbon reactor (dis)",   Q_rx,   TB); ...
        mkunit("Trim unit (18->19)",     Q_trim, NaN); ...
        mkunit("rSOC (SOFC)",            Q_cell, Tcell); ...
        mkunit("Heater H3",              Q_H3,   NaN); ...
        mkunit("Cooler (21->22)",        Q_c22,  NaN); ...
        mkunit("CO2 separator",          Q_sep,  NaN); ...
        mkunit("Air compressor",         Q_comp, NaN)];
    [units, link] = settle_heat(units, 4, 2, cfg);     % 4 = 셀, 2 = 반응기

    Qel_tot  = sum([units.Qel]);         % 전기 가열 합계 [W]
    Qrej_tot = sum([units.Qrej]);        % 방출열 합계 [W]

    % 방전 순출력: P_out = W_cell - W_sep - W_comp - (전기 가열 합계)
    P_out = W_cell - W_sep - W_comp - Qel_tot;

    % 방전 효율: eta_dis = E_out / (N_C LHV_C) = P_out / (xi LHV_C)
    eta_dis = P_out/(xi*cfg.LHV_C);

    % 배기(Stream 29)가 25 degC 기준으로 가지고 나가는 현열 (일/열 회수 없이 버려지는 손실)
    H29 = stream_H(th, n28, T29);
    Q_exhaust = H29 - stream_H(th, n28, T0);

    % ------------------------------------------------------------------
    % 7. 스트림 목록 (Stream 24 만 주위 압력, 나머지는 시스템 압력)
    % ------------------------------------------------------------------
    S = [ ...
        mkstream(14, n14,  T0,    P,    "stored carbon -> reactor"); ...
        mkstream(15, n15,  T0,    P,    "recycle -> HEX 3 (cold in)"); ...
        mkstream(16, n15,  T16,   P,    "HEX 3 -> heater H2"); ...
        mkstream(17, n15,  TB,    P,    "H2 -> carbon reactor"); ...
        mkstream(18, n18,  TB,    P,    "reactor outlet (equilibrium)"); ...
        mkstream(19, n18,  Tcell, P,    "trim unit -> SOFC fuel inlet"); ...
        mkstream(20, n20,  Tcell, P,    "SOFC fuel outlet -> HEX 3 (hot in)"); ...
        mkstream(21, n20,  T21,   P,    "HEX 3 -> cooler"); ...
        mkstream(22, n20,  T0,    P,    "cooler -> CO2 separator"); ...
        mkstream(23, n23,  T0,    P,    "pure CO2 -> tank"); ...
        mkstream(24, nair, T0,    Pamb, "ambient air"); ...
        mkstream(25, nair, T0,    P,    "compressed air -> HEX 4 (cold in)"); ...
        mkstream(26, nair, T26,   P,    "HEX 4 -> heater H3"); ...
        mkstream(27, nair, Tcell, P,    "H3 -> SOFC air inlet"); ...
        mkstream(28, n28,  Tcell, P,    "depleted air -> HEX 4 (hot in)"); ...
        mkstream(29, n28,  T29,   P,    "exhaust")];

    % ------------------------------------------------------------------
    % 8. 원소/에너지수지 잔차: {이름, 입구, 출구, Q [W], W [W]}
    % ------------------------------------------------------------------
    U = { ...
        'HEX 3',              [15 20], [16 21], 0,         0; ...
        'Heater H2',          16,      17,      Q_H2,      0; ...
        'Carbon reactor',     [14 17], 18,      Q_rx,      0; ...
        'Trim unit',          18,      19,      Q_trim,    0; ...
        'rSOC (SOFC)',        [19 27], [20 28], Q_cell,    -W_cell; ...
        'Cooler',             21,      22,      Q_c22,     0; ...
        'CO2 separator',      22,      [23 15], Q_sep,     W_sep; ...
        'Air compressor',     24,      25,      Q_comp,    W_comp; ...
        'HEX 4',              [25 28], [26 29], 0,         0; ...
        'Heater H3',          26,      27,      Q_H3,      0; ...
        'DISCHARGING SYSTEM', [14 24], [23 29], -Qrej_tot, -P_out};
    % 전체 시스템 행: 순수하게 들어오는 전기 = -P_out, 나가는 열 = Qrej_tot

    % ------------------------------------------------------------------
    % 9. 결과 정리
    % ------------------------------------------------------------------
    R.mode    = 'discharging';
    R.cfg     = cfg;
    R.S       = S;
    R.streams = stream_table(S, th.MW);
    R.balance = balance_table(th, S, U);
    R.units   = units;
    R.heat    = heat_table(units);
    R.link    = link;
    R.hex     = [ ...
        struct('name', "HEX 3", 'Q', hx3.Q, 'Qmax', hx3.Qmax, 'dTmin', hx3.dTmin, ...
               'Th_in', Tcell, 'Th_out', T21, 'Tc_in', T0, 'Tc_out', T16); ...
        struct('name', "HEX 4", 'Q', hx4.Q, 'Qmax', hx4.Qmax, 'dTmin', hx4.dTmin, ...
               'Th_in', Tcell, 'Th_out', T29, 'Tc_in', T0, 'Tc_out', T26)];
    R.eq      = struct( ...
        'T_C',        cfg.TBdis_C, ...           % 반응기 온도 [degC]
        'dG',         dG_B/1000, ...             % 석출 방향 dG° [J/mol]
        'K_B',        K_B, ...                   % 평형상수 (석출 방향 정의) [-]
        'yCO',        yCO_18, ...                % 반응기 출구 CO 몰분율
        'yCO2',       yCO2_18, ...               % 반응기 출구 CO2 몰분율
        'res_lnK',    res_eq, ...                % 평형 잔차
        'n_recycle',  a + b, ...                 % 재순환 유량 (Stream 15) [kmol/s]
        'n_cell_in',  m + c, ...                 % 셀 연료극 입구 유량 (Stream 19) [kmol/s]
        'X_pass',     xi/b);                     % 반응기 1회 통과 CO2 전환율 [-]
    R.cell    = struct('E0', E0, 'E_rev', E_rev, 'V', V_dis, 'V_tn', V_tn, ...
                       'I', I_cell, 'W', W_cell);
    R.sep     = struct('xCO2_15', xCO2_15, 'V', V_sep, 'I', I_sep, 'W', W_sep);
    R.elec    = struct('W_cell', W_cell, 'Qel_tot', Qel_tot, 'W_sep', W_sep, 'W_comp', W_comp);
    R.Qel_tot  = Qel_tot;
    R.Qrej_tot = Qrej_tot;
    R.Q_exhaust = Q_exhaust;
    R.P_out   = P_out;
    R.nC      = xi;
    R.eta_dis = eta_dis;
end

% ======================================================================
function u = mkunit(name, Q, T)
% 장치 구조체: 이름, 열부하 Q [W], 온도 T [K]
    u = struct('name', name, 'Q', Q, 'T', T);
end

function s = mkstream(no, n, T, P, name)
% 스트림 구조체: 번호, 몰유량 [kmol/s], 온도 [K], 압력 [Pa], 설명
    s = struct('no', no, 'n', n, 'T', T, 'P', P, 'name', name);
end
