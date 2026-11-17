function R = charging_mode(cfg, th)
% Charging validation extension based on Woojeong-2026/Energy-System-Design
% commit fb28621885af0df18ac6d20412ad4db788992efe. Original equations retained.
% Added input-domain checks; see charging/verify_charging.m for independent checks.
% CHARGING_MODE  충전 모드 정상상태 모델 (Figure 2(a), Stream 1~13)
%
%   R = charging_mode(cfg, th)
%
%   cfg : noon_config() 설정 구조체,  th : thermo_init() 물성 구조체
%
% [통합 인터페이스] main.m / run_roundtrip.m 이 쓰는 출력 필드 (discharging_mode 와 같은 형식)
%   R.S       : 스트림 구조체 배열 (no, n [kmol/s], T [K], P [Pa], name)
%   R.streams : 스트림 결과표                                   -> (a)
%   R.balance : 장치별 + 전체 원소/에너지수지 잔차표            -> (a)
%   R.units   : 장치별 열부하/link/전기가열/방출열 구조체 배열  -> (b)
%   R.heat    : 위 내용을 표로 만든 것 [MW]
%   R.hex     : 열교환기(recuperator) 결과                      -> (b)
%   R.link    : thermal link 전달열과 방향
%   R.eq      : Boudouard 평형과 재순환량
%   R.cell    : 셀 전압/전류/전력
%   R.elec    : 전력 항목별 내역 [W]
%   R.P_in    : 충전 입력 전력 [W],  R.nC : 탄소 저장 속도 [kmol/s]
%   R.eta_ch  : 충전 효율 [-]                                   -> (c)
%
% 이 파일은 충전 담당 모듈이다. 위 필드 형식만 지키면 내부 구현을 바꿔도
% main.m 은 그대로 동작한다.
%
% 부호 규약: Q > 0 장치로 들어가는 열, W > 0 장치로 들어가는 전력
% 단위: 몰유량 kmol/s, 엔탈피 J/kmol, 열/전력 W, 온도 K, 압력 Pa

    % ------------------------------------------------------------------
    % 1. 운전 조건
    % ------------------------------------------------------------------
    validate_charge_inputs(cfg, th);
    P     = cfg.P_bar*1e5;               % 시스템 압력 [Pa]
    T0    = cfg.T0_C + 273.15;           % 주위/저장 온도 [K]
    Tcell = cfg.Tcell_C + 273.15;        % 셀 온도 [K]
    TB    = cfg.TBch_C + 273.15;         % 충전 Boudouard 반응기 온도 [K]
    X     = cfg.X_soec;                  % SOEC CO2 전환율 [-]
    nF    = cfg.nCO2_feed;               % Stream 1 CO2 공급량 [kmol/s]

    % ------------------------------------------------------------------
    % 2. Boudouard 평형상수 (석출: 2CO -> C(s) + CO2, T = TB,ch)
    % ------------------------------------------------------------------
    % K_B = exp(-dG°_B/RT) = (p_CO2/p°)/(p_CO/p°)^2
    [K_B, dG_B] = K_boudouard(th, TB, cfg.Ru);

    % ------------------------------------------------------------------
    % 3. 재순환 루프 물질수지 (정상상태, purge 없음) - 해석해
    % ------------------------------------------------------------------
    % 미지수: Stream 2 의 CO 유량 a, CO2 유량 b, 탄소 석출량 xi [kmol/s]
    %   SOEC   (CO2 -> CO + 1/2 O2, 전환량 X*b)
    %            Stream 5 : CO = a + X b,  CO2 = (1-X) b      Stream 12 : O2 = X b / 2
    %   반응기 (2CO -> C + CO2, 진행도 xi)
    %            Stream 7 : CO = a + X b - 2 xi,  CO2 = (1-X) b + xi      Stream 10 : C = xi
    %   혼합기 (Stream 2 = Stream 1 + Stream 9, Stream 9 = Stream 7 전량)
    %            CO  수지 : a = a + X b - 2 xi        ->  xi  = X b / 2
    %            CO2 수지 : b = nF + (1-X) b + xi     ->  X b = nF + xi
    %   두 식에서  xi = nF,  b = 2 nF / X
    %   (전체 탄소 원자수지와 같다: CO2 1 kmol 이 들어오면 C(s) 1 kmol 이 쌓인다)
    xi = nF;                             % 탄소 석출(저장) 속도 [kmol/s]
    b  = 2*xi/X;                         % Stream 2 CO2 유량 [kmol/s]
    c  = (1 - X)*b + xi;                 % Stream 7(=9) CO2 유량 [kmol/s]

    % 재순환 CO 유량 a 는 반응기 출구(Stream 7: CO = a, CO2 = c)의 평형으로 정해진다.
    %   K_B = (y_CO2 P/p°)/(y_CO P/p°)^2 = c (a + c) / (a^2 P/p°)
    %   ->  (K_B P/p°) a^2 - c a - c^2 = 0  의 양의 근
    KP = K_B*(P/cfg.p0);                 % K_B * (P/p°) [-]
    a  = c*(1 + sqrt(1 + 4*KP))/(2*KP);  % 재순환 CO 유량 [kmol/s]

    yCO_7  = a/(a + c);                  % Stream 7 CO 몰분율 (평형)
    yCO2_7 = c/(a + c);                  % Stream 7 CO2 몰분율 (평형)
    % 평형 잔차: ln K_B - ln[(p_CO2/p°)/(p_CO/p°)^2]  (0 이어야 한다)
    res_eq = log(K_B) - log((yCO2_7*P/cfg.p0)/(yCO_7*P/cfg.p0)^2);

    % ------------------------------------------------------------------
    % 4. 스트림 조성 [CO CO2 O2 N2 C(s)] [kmol/s]
    % ------------------------------------------------------------------
    n1  = [0, nF, 0, 0, 0];              % CO2 탱크 -> 혼합기
    n9  = [a, c, 0, 0, 0];               % 재순환 가스
    n2  = n1 + n9;                       % 혼합기 출구 (CO = a, CO2 = b)
    n5  = [a + X*b, (1 - X)*b, 0, 0, 0]; % SOEC 연료극 출구
    n12 = [0, 0, X*b/2, 0, 0];           % SOEC 산소극 출구 (순수 O2)
    n7  = n9;                            % 반응기 출구 가스 (전량 재순환)
    n10 = [0, 0, 0, 0, xi];              % 석출 탄소

    % ------------------------------------------------------------------
    % 5. 에너지수지
    % ------------------------------------------------------------------
    % (1) 단열 혼합기: H2(T2) = H1(T0) + H9(T0)  ->  T2
    H1 = stream_H(th, n1, T0);
    H9 = stream_H(th, n9, T0);
    T2 = stream_T(th, n2, H1 + H9, T0, T0);
    H2 = stream_H(th, n2, T2);

    % (2) HEX 1 (단열 향류, 유용도 eff): 고온측 Stream 7 -> 8, 저온측 Stream 2 -> 3
    hx1 = hex_counterflow(th, n7, TB, n2, T2, cfg.eff_hex);
    T8  = hx1.Th_out;                    % Stream 8 온도 [K]
    T3  = hx1.Tc_out;                    % Stream 3 온도 [K]

    % (3) Heater H1: Stream 3 -> 4 (셀 온도까지),  Q_H1 = H4 - H3
    H3   = stream_H(th, n2, T3);
    H4   = stream_H(th, n2, Tcell);
    Q_H1 = H4 - H3;

    % (4) SOEC (등온, Tcell)
    %     가역전압은 과제 지정 스트림 조성으로 평가: 연료극 Stream 5, 산소극 Stream 12 (순수 O2 -> p_O2 = P)
    yCO_5  = n5(1)/sum(n5);              % Stream 5 CO 몰분율
    yCO2_5 = n5(2)/sum(n5);              % Stream 5 CO2 몰분율
    [E_rev, E0] = E_cell_CO(th, Tcell, yCO_5*P/cfg.p0, yCO2_5*P/cfg.p0, 1.0*P/cfg.p0, ...
                            cfg.Ru, cfg.F);
    V_ch   = E_rev + cfg.dV_loss;        % 충전 운전전압 V_ch = E_rev + 0.10 [V]
    nCO_el = X*b;                        % 셀에서 전환된 CO2 (= 생성 CO) [kmol/s]
    I_cell = 2*cfg.F*(nCO_el*1000);      % 전류 I = 2 F n_CO [A]  (kmol/s -> mol/s)
    W_cell = I_cell*V_ch;                % 셀 입력 전력 = I V [W]

    %     등온 셀 에너지수지: H4 + W_cell + Q_cell = H5 + H12
    H5     = stream_H(th, n5,  Tcell);
    H12    = stream_H(th, n12, Tcell);
    Q_cell = (H5 + H12) - H4 - W_cell;   % 셀 열부하 [W] (V < V_tn 이면 > 0, 흡열 운전)
    V_tn   = (H5 + H12 - H4)/I_cell;     % 열중립 전압 V_tn = dH/(2F) [V]

    % (5) Trim unit: Stream 5 (Tcell) -> 6 (TB,ch), 조성 동결,  Q_trim = H6 - H5
    H6     = stream_H(th, n5, TB);
    Q_trim = H6 - H5;

    % (6) Boudouard 반응기 + 고체 분리기 (등온, TB,ch): Stream 6 -> 7 (가스) + 10 (탄소)
    %     Q_rx = H7 + H10 - H6   (석출은 발열 -> Q_rx < 0)
    H7   = stream_H(th, n7,  TB);
    H10  = stream_H(th, n10, TB);
    Q_rx = (H7 + H10) - H6;

    % (7) 냉각기: 모두 25 degC 까지 냉각 (Q < 0, 방출)
    H8    = stream_H(th, n7, T8);
    Q_c89 = H9 - H8;                     % 재순환 냉각기 8 -> 9
    H11   = stream_H(th, n10, T0);
    Q_cC  = H11 - H10;                   % 탄소 냉각기 10 -> 11
    H13   = stream_H(th, n12, T0);
    Q_cO2 = H13 - H12;                   % O2 냉각기 12 -> 13

    % ------------------------------------------------------------------
    % 6. thermal link 적용, 전기 가열/방출열 정산
    % ------------------------------------------------------------------
    units = [ ...
        mkunit("Heater H1",              Q_H1,   NaN); ...
        mkunit("rSOC (SOEC)",            Q_cell, Tcell); ...
        mkunit("Trim unit (5->6)",       Q_trim, NaN); ...
        mkunit("Carbon reactor (ch)",    Q_rx,   TB); ...
        mkunit("Recycle cooler (8->9)",  Q_c89,  NaN); ...
        mkunit("C cooler (10->11)",      Q_cC,   NaN); ...
        mkunit("O2 cooler (12->13)",     Q_cO2,  NaN)];
    [units, link] = settle_heat(units, 2, 4, cfg);     % 2 = 셀, 4 = 반응기

    Qel_tot  = sum([units.Qel]);         % 전기 가열 합계 [W]
    Qrej_tot = sum([units.Qrej]);        % 방출열 합계 [W]

    % 충전 입력 전력: P_in = W_cell + (전기 가열 합계)
    %   충전 모드에는 공기 압축기와 CO2 분리기가 없고, O2 는 일 회수 없이 밸브로 방출한다.
    P_in = W_cell + Qel_tot;

    % 충전 효율: eta_ch = N_C LHV_C / E_in = xi LHV_C / P_in
    eta_ch = xi*cfg.LHV_C/P_in;

    % ------------------------------------------------------------------
    % 7. 스트림 목록
    % ------------------------------------------------------------------
    S = [ ...
        mkstream(1,  n1,  T0,    P, "CO2 tank -> mixer"); ...
        mkstream(2,  n2,  T2,    P, "mixer -> HEX 1 (cold in)"); ...
        mkstream(3,  n2,  T3,    P, "HEX 1 -> heater H1"); ...
        mkstream(4,  n2,  Tcell, P, "H1 -> SOEC fuel inlet"); ...
        mkstream(5,  n5,  Tcell, P, "SOEC fuel outlet"); ...
        mkstream(6,  n5,  TB,    P, "trim unit -> carbon reactor"); ...
        mkstream(7,  n7,  TB,    P, "reactor gas -> HEX 1 (hot in)"); ...
        mkstream(8,  n7,  T8,    P, "HEX 1 -> recycle cooler"); ...
        mkstream(9,  n9,  T0,    P, "recycle -> mixer"); ...
        mkstream(10, n10, TB,    P, "deposited carbon"); ...
        mkstream(11, n10, T0,    P, "carbon -> storage"); ...
        mkstream(12, n12, Tcell, P, "SOEC O2 outlet"); ...
        mkstream(13, n12, T0,    P, "O2 -> vent valve")];

    % ------------------------------------------------------------------
    % 8. 원소/에너지수지 잔차: {이름, 입구, 출구, Q [W], W [W]}
    % ------------------------------------------------------------------
    U = { ...
        'Mixer',               [1 9],  2,       0,         0; ...
        'HEX 1',               [2 7],  [3 8],   0,         0; ...
        'Heater H1',           3,      4,       Q_H1,      0; ...
        'rSOC (SOEC)',         4,      [5 12],  Q_cell,    W_cell; ...
        'Trim unit',           5,      6,       Q_trim,    0; ...
        'Reactor / separator', 6,      [7 10],  Q_rx,      0; ...
        'Recycle cooler',      8,      9,       Q_c89,     0; ...
        'C cooler',            10,     11,      Q_cC,      0; ...
        'O2 cooler',           12,     13,      Q_cO2,     0; ...
        'CHARGING SYSTEM',     1,      [11 13], -Qrej_tot, P_in};
    % 전체 시스템 행: 들어오는 전기 P_in (셀 + 전기 가열), 나가는 열 Qrej_tot

    % ------------------------------------------------------------------
    % 9. 결과 정리
    % ------------------------------------------------------------------
    R.mode    = 'charging';
    R.cfg     = cfg;
    R.S       = S;
    R.streams = stream_table(S, th.MW);
    R.balance = balance_table(th, S, U);
    R.units   = units;
    R.heat    = heat_table(units);
    R.link    = link;
    R.hex     = struct('name', "HEX 1", 'Q', hx1.Q, 'Qmax', hx1.Qmax, 'dTmin', hx1.dTmin, ...
                       'Th_in', TB, 'Th_out', T8, 'Tc_in', T2, 'Tc_out', T3);
    R.eq      = struct( ...
        'T_C',        cfg.TBch_C, ...            % 반응기 온도 [degC]
        'dG',         dG_B/1000, ...             % 석출 방향 dG° [J/mol]
        'K_B',        K_B, ...                   % 평형상수 [-]
        'yCO',        yCO_7, ...                 % 반응기 출구 CO 몰분율
        'yCO2',       yCO2_7, ...                % 반응기 출구 CO2 몰분율
        'res_lnK',    res_eq, ...                % 평형 잔차
        'n_recycle',  a + c, ...                 % 재순환 유량 (Stream 9) [kmol/s]
        'n_cell_in',  a + b, ...                 % 셀 연료극 입구 유량 (Stream 4) [kmol/s]
        'X_pass',     2*xi/n5(1));               % 반응기 1회 통과 CO 전환율 [-]
    R.cell    = struct('E0', E0, 'E_rev', E_rev, 'V', V_ch, 'V_tn', V_tn, ...
                       'I', I_cell, 'W', W_cell);
    R.elec    = struct('W_cell', W_cell, 'Qel_tot', Qel_tot, 'W_sep', 0, 'W_comp', 0);
    R.Qel_tot  = Qel_tot;
    R.Qrej_tot = Qrej_tot;
    R.P_in    = P_in;
    R.nC      = xi;
    R.eta_ch  = eta_ch;
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

function validate_charge_inputs(cfg, th)
% Fixed process topology, ideal-gas/graphite model; excludes zero-flow shutdown.
    fields = {'P_bar','T0_C','Tcell_C','TBch_C','nCO2_feed','X_soec', ...
              'eff_hex','f_link','dT_link','dV_loss','Ru','F','p0','LHV_C'};
    for k = 1:numel(fields)
        f = fields{k};
        if ~isfield(cfg,f) || ~isnumeric(cfg.(f)) || ~isscalar(cfg.(f)) || ...
                ~isreal(cfg.(f)) || ~isfinite(cfg.(f))
            error('charging_mode:input','%s must be a finite real scalar.',f);
        end
    end
    if cfg.P_bar <= 0 || cfg.nCO2_feed <= 0 || cfg.X_soec <= 0 || cfg.X_soec >= 1
        error('charging_mode:input','Require P>0, feed>0 and 0<X_soec<1.');
    end
    if cfg.eff_hex < 0 || cfg.eff_hex > 1 || cfg.f_link < 0 || cfg.f_link > 1 || ...
            cfg.dT_link < 0 || cfg.dV_loss < 0 || cfg.Ru <= 0 || cfg.F <= 0 || cfg.LHV_C <= 0
        error('charging_mode:input','Invalid heat-recovery fraction or constant.');
    end
    T = [cfg.T0_C cfg.Tcell_C cfg.TBch_C] + 273.15;
    if any(T < 200 | T > 3500) || cfg.Tcell_C <= cfg.T0_C || cfg.TBch_C <= cfg.T0_C
        error('charging_mode:input','Temperatures outside this model domain (200..3500 K, hot units>T0).');
    end
    if abs(cfg.p0-1e5)>1e-6 || abs(th.p0-cfg.p0)>1e-6
        error('charging_mode:input','Use a consistent 1 bar standard state in cfg and thermo.');
    end
end

