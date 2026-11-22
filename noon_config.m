function cfg = noon_config(varargin)
% NOON_CONFIG  과제 지정 조건을 한곳에 모은 설정 구조체 (팀 공통 입력)
%
%   cfg = noon_config()                      % 과제 기준 조건
%   cfg = noon_config('P_bar', 10)           % 이름-값 쌍으로 일부만 변경 (스윕/민감도용)
%   cfg = noon_config('Tcell_C', 900, 'eff_hex', 0.95)
%
% 모든 모드 함수(run_charging_analysis, run_discharging_analysis)는 이 구조체만 입력으로 받는다.
% 숫자를 코드 안에 직접 쓰지 않고 여기서만 바꾸면 전체 결과가 일관되게 바뀐다.
%
% 출처: Term_Project_2026_Noon_Energy_revV3.pdf "Assumptions and model definitions"

    % ---------------- 운전 조건 ----------------
    cfg.P_bar    = 5;        % 시스템 압력 P [bar] (연료측/산소측 모든 번호 스트림)
    cfg.Tcell_C  = 850;      % 셀(rSOC) 온도 [degC] (충전/방전 공통)
    cfg.TBch_C   = 650;      % 충전 시 Boudouard 반응기 온도 TB,ch [degC]
    cfg.TBdis_C  = 800;      % 방전 시 Boudouard 반응기 온도 TB,dis [degC]
    cfg.T0_C     = 25;       % 주위/저장 온도 [degC]
    cfg.P0_bar   = 1;        % 주위 압력 [bar] (압축 전 공기)

    % ---------------- 유량 기준 ----------------
    cfg.nCO2_feed = 1.0;     % 충전: Stream 1 의 CO2 공급량 [kmol/s]
    cfg.nC_dis    = 1.0;     % 방전: 탄소 소비 속도 (= Stream 23 CO2 분리량) [kmol/s]

    % ---------------- 셀 ----------------
    cfg.X_soec   = 0.70;     % SOEC: 입구 CO2 의 CO 전환율 [-]
    cfg.U_sofc   = 0.80;     % SOFC: 입구 CO 의 산화율(연료 이용률) [-]
    cfg.dV_loss  = 0.10;     % 전압 손실 가정 [V]  (V_ch = E_rev + dV, V_dis = E_rev - dV)

    % ---------------- 공기 ----------------
    cfg.lambda_O2 = 3.0;     % 공급 O2 / 양론 O2 [-]  (3배 = 200 % 과잉공기)
    cfg.xO2_air   = 0.21;    % 건공기 O2 몰분율 [-]
    cfg.xN2_air   = 0.79;    % 건공기 N2 몰분율 [-]

    % ---------------- 열회수 ----------------
    cfg.eff_hex  = 0.85;     % 향류 열교환기 유용도 (HEX 1, 3, 4) [-]
    cfg.f_link   = 0.85;     % thermal link 로 회수 가능한 donor 열의 비율 [-]
    cfg.dT_link  = 20;       % thermal link 가 작동하는 최소 온도차 [K]

    % ---------------- 상수 (과제 지정값) ----------------
    cfg.Ru    = 8314;        % 일반기체상수 [J/(kmol K)]  (= 8.314 J/(mol K))
    cfg.F     = 96485;       % Faraday 상수 [C/mol]
    cfg.p0    = 1e5;         % 표준압력 p° = 1 bar [Pa]
    cfg.LHV_C = 393.5e6;     % 탄소 LHV [J/kmol]  (= 393.5 kJ/mol, 25 degC)

    % ---------------- 이름-값 쌍으로 덮어쓰기 ----------------
    if mod(numel(varargin), 2) ~= 0
        error('noon_config:args', '이름-값 쌍으로 입력해야 합니다.');
    end
    for k = 1:2:numel(varargin)
        name = varargin{k};
        if ~isfield(cfg, name)
            error('noon_config:field', '알 수 없는 설정 이름: %s', name);
        end
        cfg.(name) = varargin{k + 1};
    end
end
