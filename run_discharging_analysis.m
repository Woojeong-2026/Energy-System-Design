function R = run_discharging_analysis(cfg, th)
% RUN_DISCHARGING_ANALYSIS 방전 시스템 계산 (Streams 14~29).
%   R = run_discharging_analysis()        : 기준 조건, 결과 화면 출력
%   R = run_discharging_analysis(cfg)     : 사용자 조건, 출력 없이 반환
%   R = run_discharging_analysis(cfg, th) : Cantera 객체 재사용 (main.m용)
%   cfg: noon_config() 반환값 / th: thermo_init() 반환값
%
% Energy System Design (2026) - Discharging mode, Stream 14~29
%
% [This file belongs to the discharging part; NOT a modification of team 3 helpers.]
%
%
% 사용 공통 함수 
%   noon_config, thermo_init, K_boudouard, E_cell_CO,
%   stream_H, hex_counterflow, settle_heat,
%   stream_table, balance_table, heat_table
%   위 함수 안에서 sp_h, sp_s0, sp_g0, stream_T 등이 호출된다.
%
% 실행법
%   1) 팀 저장소 최상위에 이 파일을 놓는다.
%   2) Cantera MATLAB toolbox를 MATLAB path에 등록한다.
%   3) MATLAB 명령창에서 run_discharging_analysis 실행.
%
% 단위 및 부호
%   n_i [kmol/s], T [K], P [Pa], H 및 Q 및 W [W]
%   Q > 0: 장치로 들어가는 열, Q < 0: 장치에서 나가는 열
%   W_SOFC > 0: 셀이 생산하는 전력, W_comp/W_sep > 0: 소비전력
%   Q_rejected > 0: 외부로 버린 열
%
% 이 코드는 과제 조건 외의 평형상수, 성분 유량, 손실 계수 등을
% 수치로 고정하지 않는다. 물성값은 공통 Cantera 함수가 계산한다.

%% 1. 공통 함수의 경로와 입력 조건
% main.m이나 임시 modes 파일을 읽지 않는다.
projectRoot = fileparts(mfilename('fullpath'));
if isempty(projectRoot), projectRoot = pwd; end
if ~isfile(fullfile(projectRoot, 'noon_config.m'))
    % 이 스크립트를 별도 하위 폴더에 놓은 경우: 상위 저장소 탐색
    projectRoot = fileparts(projectRoot);
end
if ~isfile(fullfile(projectRoot, 'noon_config.m'))
    error(['팀 저장소의 최상위 폴더에 run_discharging_analysis.m을 넣어 주세요. ' ...
           'noon_config.m을 찾지 못했습니다.']);
end

addpath(projectRoot, fullfile(projectRoot,'thermo'), ...
         fullfile(projectRoot,'system'));

required = {'noon_config','thermo_init','K_boudouard','E_cell_CO', ...
    'stream_H','hex_counterflow','settle_heat', ...
    'stream_table','balance_table','heat_table'};
for j = 1:numel(required)
    if exist(required{j},'file') ~= 2
        error('팀원 3 공통 함수를 찾지 못했습니다: %s.m', required{j});
    end
end

% 외부에서 받은 설정/Cantera 객체는 덮어쓰지 않는다.
if nargin < 1 || isempty(cfg), cfg = noon_config(); end
if nargin < 2 || isempty(th),  th  = thermo_init(); end

P = cfg.P_bar*1e5;                 % [Pa] 운전 압력
P_amb = cfg.P0_bar*1e5;           % [Pa] 압축 전 대기압 / 배기 압력
P0 = cfg.p0;                      % [Pa] standard pressure
T_amb = cfg.T0_C + 273.15;       % [K] 저장 및 주위 온도
T_cell = cfg.Tcell_C + 273.15;   % [K] SOFC 온도
T_reactor = cfg.TBdis_C + 273.15; % [K] Boudouard reactor 온도
nC = cfg.nC_dis;                  % [kmol/s] 탄소 소비량 = 분리할 CO2 양
X_SOFC = cfg.U_sofc;              % [-] 셀 CO conversion

assert(P >= P_amb && nC > 0 && X_SOFC > 0 && X_SOFC <= 1, ...
    '압력, 탄소 유량 또는 SOFC 산화율 설정을 확인하세요.');
assert(cfg.xO2_air > 0 && cfg.xN2_air > 0 && ...
       abs(cfg.xO2_air + cfg.xN2_air - 1) < 100*eps, ...
    '공기 조성은 양수이며 합이 1이어야 합니다.');

%% 2. Boudouard 평형 및 Recycle (Stream 15, 18, 20, 23)
% Gasification: C(s) + CO2 -> 2CO
% 팀원 3의 K_boudouard는 역반응이 아닌 '석출 방향' K_B를 돌려준다.
% 2CO -> C(s) + CO2  : K_B = (yCO2*P/P0)/(yCO*P/P0)^2
[Kp_dep, dG_dep] = K_boudouard(th, T_reactor, cfg.Ru);
Kp_gas = 1/Kp_dep;               % gasification 방향의 평형상수

% CO 정상상태 recycle 수지:
%   nCO18 = nCO15 + 2*nC
%   nCO15 = (1-X_SOFC)*nCO18
% 따라서 nCO18 = 2*nC/X_SOFC  (임의 recycle 유량을 가정하지 않음)
n18_CO = 2*nC/X_SOFC;
n15_CO = (1-X_SOFC)*n18_CO;

% 평형식에서 반응기 출구 CO2 유량을 직접 계산한다.
%   Kp_dep = [c*(m+c)/m^2]*(P0/P),  m=nCO18, c=nCO2_18
%   c^2+m*c - Kp_dep*(P/P0)*m^2 = 0
% 양의 근을 안정적인 형태로 평가 (작은 수끼리 뺄셈 방지).
KP = Kp_dep*(P/P0);
n18_CO2 = n18_CO*(2*KP)/(1 + sqrt(1 + 4*KP));
n15_CO2 = n18_CO2 + nC;          % 반응기에 들어와 소비되는 CO2

nCO_reacted = X_SOFC*n18_CO;      % SOFC에서 산화되는 CO
n20_CO = n18_CO - nCO_reacted;
n20_CO2 = n18_CO2 + nCO_reacted;
nCO2_23 = nC;                    % separator에서 순수 CO2를 분리

n15_total = n15_CO + n15_CO2;
n18_total = n18_CO + n18_CO2;
n20_total = n20_CO + n20_CO2;
y15_CO2 = n15_CO2/n15_total;
y18_CO = n18_CO/n18_total;
y18_CO2 = n18_CO2/n18_total;

% 조성 벡터 순서는 thermo_init의 공통 정의: [CO CO2 O2 N2 C(s)]
n14 = [0 0 0 0 nC];
n15 = [n15_CO n15_CO2 0 0 0];
n18 = [n18_CO n18_CO2 0 0 0];
n20 = [n20_CO n20_CO2 0 0 0];
n23 = [0 nCO2_23 0 0 0];

% 평형식 재대입: 0에 가까우면 평형 조건이 만족됨
Kp_dep_check = (y18_CO2*P/P0)/(y18_CO*P/P0)^2;
res_eq = log(Kp_dep_check/Kp_dep);

%% 3. 공기 공급 및 Reversible Isothermal Compressor (24, 25, 27, 28)
% SOFC 반응: CO + (1/2)O2 -> CO2
nO2_stoich = nCO_reacted/2;
nO2_24 = cfg.lambda_O2*nO2_stoich;
n24_total = nO2_24/cfg.xO2_air;
nN2_24 = n24_total*cfg.xN2_air;

nair = [0 0 nO2_24 nN2_24 0];
nO2_28 = nO2_24-nO2_stoich;
n28 = [0 0 nO2_28 nN2_24 0];
n28_total = sum(n28);
yO2_28 = nO2_28/n28_total;

% 공기: 25 C의 가역 등온 압축. 문제의 R=8.314 J/(mol K)를
% 그대로 적용하므로 cfg.Ru=8314 J/(kmol K)를 사용한다.
% 예전 코드에서 cp_mole-cv_mole로 얻은 Cantera 정밀 R과는
% 수 kW의 반올림 차이가 날 수 있다 (아래 README 참고).
W_comp = n24_total*cfg.Ru*T_amb*log(P/P_amb);
Q_comp = -W_comp;                % 이상기체 등온: H25-H24=0

%% 4. SOFC 전압, 전류, 출력 (19, 20, 27, 28)
% 문제는 SOC 전극의 stream 20(연료측), 28(공기측) 조성을 사용한다.
y20_CO = n20_CO/n20_total;
y20_CO2 = n20_CO2/n20_total;
[E_rev_SOFC, E0_SOFC, dG_SOFC] = E_cell_CO(th, T_cell, ...
    y20_CO*P/P0, y20_CO2*P/P0, yO2_28*P/P0, cfg.Ru, cfg.F);
V_SOFC = E_rev_SOFC-cfg.dV_loss;
I_SOFC = 2*cfg.F*nCO_reacted*1000;     % [A], kmol/s -> mol/s
W_SOFC_rev = I_SOFC*E_rev_SOFC;
W_SOFC = I_SOFC*V_SOFC;
W_voltage_loss = W_SOFC_rev-W_SOFC;

%% 5. CO2 Separator (22 -> 23 + 15)
% 과제: 분리기는 25 C 등온, 100% Faradaic, voltage loss 없음.
% recycle 15와 같은 조성의 완전혼합 측에서 pure CO2 23으로 이송한다.
T_sep = T_amb;
E_sep = (cfg.Ru/1000)*T_sep/(2*cfg.F)*log(1/y15_CO2);
I_sep = 2*cfg.F*nCO2_23*1000;
W_sep = I_sep*E_sep;
Q_sep = -W_sep;                    % 등온 유지를 위해 입력 전력만큼 열 방출

%% 6. HEX 3, HEX 4 (팀원 3 공통 hex_counterflow 사용)
% HEX 3: Hot 20->21, Cold 15->16
hx3 = hex_counterflow(th, n20, T_cell, n15, T_amb, cfg.eff_hex);
T16 = hx3.Tc_out;
T21 = hx3.Th_out;
Q_HEX3 = hx3.Q;

% HEX 4: Hot 28->29, Cold 25->26
hx4 = hex_counterflow(th, n28, T_cell, nair, T_amb, cfg.eff_hex);
T26 = hx4.Tc_out;
T29 = hx4.Th_out;
Q_HEX4 = hx4.Q;

% 내부 모든 검사점에서 열 흐름의 방향이 올바른지 확인한다.
assert(hx3.dTmin > 0 && hx4.dTmin > 0, ...
    '열교환기 내부에서 hot/cold stream 온도가 교차했습니다.');

%% 7. H2/H3, Gasifier, Trim, SOFC 열수지
% stream_H는 Cantera에서 온도별 성분 엔탈피를 구해
% 성분별 몰유량을 곱해 합산하므로 반응열이 포함된다.
H14 = stream_H(th, n14, T_amb);
H15 = stream_H(th, n15, T_amb);
H16 = stream_H(th, n15, T16);
H17 = stream_H(th, n15, T_reactor);
H18 = stream_H(th, n18, T_reactor);
H19 = stream_H(th, n18, T_cell);
H20 = stream_H(th, n20, T_cell);
H21 = stream_H(th, n20, T21);
H22 = stream_H(th, n20, T_amb);
H23 = stream_H(th, n23, T_amb);
H24 = stream_H(th, nair, T_amb);
H25 = stream_H(th, nair, T_amb);
H26 = stream_H(th, nair, T26);
H27 = stream_H(th, nair, T_cell);
H28 = stream_H(th, n28, T_cell);
H29 = stream_H(th, n28, T29);

Q_H2 = H17-H16;                    % Heater 2
Q_gasifier = H18-H17-H14;         % 반응기 흡열, carbon inlet 14 포함
Q_conditioner = H19-H18;          % Trim: 18->19 (조성 동결)
Q_H3 = H27-H26;                    % Heater 3
Q_SOFC = H20+H28-H19-H27+W_SOFC; % SOFC: Q-W=Hout-Hin
Q_cooler_22 = H22-H21;           % Cooler 21->22 (음수면 열 방출)

%% 8. Thermal Link 및 Electrical Heating (팀원 3 settle_heat 사용)
% Heat duty를 모두 배열에 등록해 SOFC <-> Reactor 사이에서만 열을 회수한다.
% 나머지 양(+)의 열부하는 전기 가열, 음(-)의 열부하는 열 방출로 처리한다.
unit_names = { ...
    'Heater H2', 'Carbon gasifier', 'Trim 18-19', 'SOFC', ...
    'Heater H3', 'Cooler 21-22', 'CO2 separator', 'Air compressor'};
unit_Q = [Q_H2, Q_gasifier, Q_conditioner, Q_SOFC, ...
          Q_H3, Q_cooler_22, Q_sep, Q_comp];
unit_T = [NaN, T_reactor, NaN, T_cell, NaN, NaN, NaN, NaN];
units = repmat(struct('name','','Q',0,'T',NaN), numel(unit_names),1);
for j = 1:numel(unit_names)
    units(j).name = unit_names{j};
    units(j).Q = unit_Q(j);
    units(j).T = unit_T(j);
end
idxCell = find(strcmp({units.name}, 'SOFC'));
idxRx   = find(strcmp({units.name}, 'Carbon gasifier'));
[units, link] = settle_heat(units, idxCell, idxRx, cfg);

Q_link = link.Q;
W_heating_total = sum([units.Qel]);
Q_rejected_total = sum([units.Qrej]);
W_H2 = units(strcmp({units.name},'Heater H2')).Qel;
W_H3 = units(strcmp({units.name},'Heater H3')).Qel;
W_conditioner = units(strcmp({units.name},'Trim 18-19')).Qel;
W_gasifier_heating = units(strcmp({units.name},'Carbon gasifier')).Qel;
W_SOFC_heating = units(strcmp({units.name},'SOFC')).Qel;
Q_SOFC_rejected = units(strcmp({units.name},'SOFC')).Qrej;
Q_gasifier_rejected = units(strcmp({units.name},'Carbon gasifier')).Qrej;

%% 9. Net Discharging Power, Efficiency, Exhaust Sensible Heat
W_dis_input = W_comp + W_sep + W_heating_total;
W_dis_net = W_SOFC - W_dis_input;
E_C_LHV = nC*cfg.LHV_C;
eta_dis = W_dis_net/E_C_LHV;

% Stream 29가 기준온도 25 C보다 더 높은 온도로 배출되면서 운반하는 현열.
% 대기 배기이며 별도 냉각기가 없으므로 Q_rejected_total에 중복하여 더하지 않는다.
Q_exhaust = H29-stream_H(th, n28, T_amb);

%% 10. Stream 14~29 테이블 (팀원 3 stream_table 이용)
% 공통 함수 표와 기존 사용자 표의 항목을 둘 다 생성한다.
stream_no = (14:29)';
flowN = [n14; n15; n15; n15; n18; n18; n20; n20; n20; n23; ...
         nair; nair; nair; nair; n28; n28];
T_stream = [T_amb; T_amb; T16; T_reactor; T_reactor; T_cell; ...
            T_cell; T21; T_amb; T_sep; T_amb; T_amb; T26; ...
            T_cell; T_cell; T29];
P_stream = P*ones(numel(stream_no),1);
P_stream(stream_no==24) = P_amb; % 압축 전 대기 공기
P_stream(stream_no==29) = P_amb; % 배출된 공기 (vent 후)
stream_description = { ...
    'Stored graphite'; 'Recycle CO/CO2'; 'After HEX 3'; 'After H2'; ...
    'Gasifier outlet'; 'After trim'; 'SOFC fuel outlet'; 'After HEX 3 hot'; ...
    'After fuel cooler'; 'Pure CO2 storage'; 'Ambient air'; 'Compressed air'; ...
    'After HEX 4'; 'After H3'; 'SOFC air outlet'; 'Vented exhaust'};

S = repmat(struct('no',0,'n',zeros(1,5),'T',0,'P',0,'name',''), ...
           numel(stream_no),1);
for j = 1:numel(stream_no)
    S(j).no   = stream_no(j);
    S(j).n    = flowN(j,:);
    S(j).T    = T_stream(j);
    S(j).P    = P_stream(j);
    S(j).name = stream_description{j};
end

% 팀 공통 표: T/P/총몰유량/질량유량/몰분율
TeamStreamTable = stream_table(S, th.MW);

% 이전에 사용한 열 이름을 유지한 표: StreamTable
n_stream = sum(flowN,2);
mdot_stream = flowN*th.MW(:);
y_stream = flowN./n_stream;
H_stream = zeros(numel(S),1);
for j = 1:numel(S)
    H_stream(j) = stream_H(th, S(j).n, S(j).T);
end
StreamTable = table(stream_no,n_stream,mdot_stream, ...
    y_stream(:,5),y_stream(:,1),y_stream(:,2),y_stream(:,3),y_stream(:,4), ...
    T_stream-273.15,P_stream/1e5,H_stream/1e6, ...
    'VariableNames',{'Stream','MolarFlow_kmol_s','MassFlow_kg_s', ...
    'yC','yCO','yCO2','yO2','yN2','Temperature_C', ...
    'Pressure_bar','EnthalpyFlow_MW'});

%% 11. 각 장치와 전체 시스템 Material / Energy Balance
% 공통 함수의 규칙: Q는 장치로 들어가는 열(+), W는 장치로 들어가는 전력(+).
% SOFC 발전: W=-W_SOFC. 전체 방전: W=-W_dis_net.
% Thermal Link는 시스템 내부 열 이동이므로 반응기/셀 자체 Q에는 원래 열부하 사용.
U = { ...
    'HEX 3',       [15 20], [16 21], 0, 0; ...
    'Heater H2',   16, 17, Q_H2, 0; ...
    'Gasifier',    [14 17], 18, Q_gasifier, 0; ...
    'Trim unit',   18, 19, Q_conditioner, 0; ...
    'SOFC',        [19 27], [20 28], Q_SOFC, -W_SOFC; ...
    'Cooler',      21, 22, Q_cooler_22, 0; ...
    'CO2 sep',     22, [23 15], Q_sep, W_sep; ...
    'Compressor',  24, 25, Q_comp, W_comp; ...
    'HEX 4',       [25 28], [26 29], 0, 0; ...
    'Heater H3',   26, 27, Q_H3, 0; ...
    'Overall discharging', [14 24], [23 29], -Q_rejected_total, -W_dis_net};

BalanceTable = balance_table(th,S,U);
HeatTable = heat_table(units);

% 별도 전체 원소 수지 (단위: kmol atoms/s)
atoms_in  = th.A*(n14 + nair)';
atoms_out = th.A*(n23 + n28)';
res_element_total = atoms_in-atoms_out;
res_C_total = res_element_total(1);
res_O_total = res_element_total(2);
res_N_total = res_element_total(3);
res_mass_total = (sum(n14.*th.MW)+sum(nair.*th.MW)) ...
               - (sum(n23.*th.MW)+sum(n28.*th.MW));
res_energy_total = H14+H24-H23-H29-W_dis_net-Q_rejected_total;
res_energy_relative = res_energy_total/(abs(H14+H24) + ...
    abs(H23+H29) + abs(W_dis_net)+abs(Q_rejected_total));
res_net_power = W_SOFC-(W_dis_net+W_comp+W_sep+W_heating_total);

%% 12. Results - 직접 실행할 때만 출력 (main.m 반복 호출 시 출력 없음)
if nargin == 0
fprintf('\n========== DISCHARGING: COMMON HELPER VERSION ==========\n');
fprintf('P = %.3f bar; Tcell = %.2f C; TBdis = %.2f C\n', ...
    cfg.P_bar,cfg.Tcell_C,cfg.TBdis_C);
fprintf('\n--- Boudouard equilibrium and recycle ---\n');
fprintf('Delta G deposition = %.6f kJ/mol\n', dG_dep/1e6);
fprintf('Kp deposition = %.9g; Kp gasification = %.9g\n', Kp_dep,Kp_gas);
fprintf('Recycle Stream 15: CO %.6f, CO2 %.6f kmol/s\n', n15_CO,n15_CO2);
fprintf('Gasifier Stream 18: CO %.6f, CO2 %.6f kmol/s\n',n18_CO,n18_CO2);
fprintf('Equilibrium ln(K) residual = %.3e\n',res_eq);

fprintf('\n--- SOFC and auxiliary electric power ---\n');
fprintf('Air supply = %.6f kmol/s; compressor = %.6f MW\n', ...
    n24_total,W_comp/1e6);
fprintf('E0 / Erev / V = %.8f / %.8f / %.8f V\n',E0_SOFC,E_rev_SOFC,V_SOFC);
fprintf('SOFC output = %.6f MW\n',W_SOFC/1e6);
fprintf('CO2 separator input = %.6f MW\n',W_sep/1e6);

fprintf('\n--- Heat recovery and unit duties ---\n');
fprintf('HEX 3: %.6f MW; dTmin = %.3f K\n',Q_HEX3/1e6,hx3.dTmin);
fprintf('HEX 4: %.6f MW; dTmin = %.3f K\n',Q_HEX4/1e6,hx4.dTmin);
fprintf('Gasifier heat duty = %.6f MW\n',Q_gasifier/1e6);
fprintf('SOFC heat duty = %.6f MW\n',Q_SOFC/1e6);
fprintf('Thermal link recovery = %.6f MW\n',Q_link/1e6);
fprintf('Heater H2 / H3 / Trim = %.6f / %.6f / %.6f MW\n', ...
    W_H2/1e6,W_H3/1e6,W_conditioner/1e6);
fprintf('Total electrical heating = %.6f MW\n',W_heating_total/1e6);
fprintf('Total rejected heat (excluding exhaust) = %.6f MW\n', ...
    Q_rejected_total/1e6);
fprintf('Exhaust stream 29 sensible heat = %.6f MW\n',Q_exhaust/1e6);

fprintf('\n========== DISCHARGING SUMMARY ==========\n');
fprintf('Net discharging power = %.6f MW\n',W_dis_net/1e6);
fprintf('Carbon LHV input = %.6f MW\n',E_C_LHV/1e6);
fprintf('Discharging efficiency = %.6f %%\n',100*eta_dis);

fprintf('\n========== DISCHARGING STREAM TABLE ==========\n');
disp(StreamTable);
fprintf('\n========== UNIT BALANCE TABLE ==========\n');
disp(BalanceTable);
fprintf('\n========== UNIT HEAT TABLE ==========\n');
disp(HeatTable);

fprintf('\n--- Overall residuals ---\n');
fprintf('Carbon = %.3e kmol atoms/s\n',res_C_total);
fprintf('Oxygen = %.3e kmol atoms/s\n',res_O_total);
fprintf('Nitrogen = %.3e kmol atoms/s\n',res_N_total);
fprintf('Mass = %.3e kg/s\n',res_mass_total);
fprintf('Overall energy = %.3e W (relative %.3e)\n', ...
    res_energy_total,res_energy_relative);
fprintf('Electrical power = %.3e W\n',res_net_power);
end % if nargin == 0

% 누락/발산 검사: 음의 net power는 온도 스윕에서 정상 결과일 수 있으므로 금지하지 않음.
assert(all(isfinite([W_dis_net,eta_dis,res_energy_total, ...
    res_element_total(:)'])), 'Non-finite results detected.');
assert(all(flowN(:)>=0), 'Negative species molar flow detected.');


%% 13. 팀 통합 인터페이스: 계산 결과를 구조체 R로 반환
% 팀 main.m, run_roundtrip.m이 요구하는 필드만 전달한다.
% 단위: 유량 kmol/s, 엔탈피/전력 W, 전압 V, 평형 dG J/mol.
V_tn = (H19 + H27 - H20 - H28)/I_SOFC;  % SOFC 열중립 전압 [V]

R.mode    = 'discharging';
R.cfg     = cfg;
R.S       = S;
R.streams = TeamStreamTable;
R.balance = BalanceTable;
R.units   = units;
R.heat    = HeatTable;
R.link    = link;
R.hex     = [ ...
    struct('name',"HEX 3",'Q',hx3.Q,'Qmax',hx3.Qmax,'dTmin',hx3.dTmin, ...
           'Th_in',T_cell,'Th_out',T21,'Tc_in',T_amb,'Tc_out',T16); ...
    struct('name',"HEX 4",'Q',hx4.Q,'Qmax',hx4.Qmax,'dTmin',hx4.dTmin, ...
           'Th_in',T_cell,'Th_out',T29,'Tc_in',T_amb,'Tc_out',T26)];
R.eq   = struct('T_C',cfg.TBdis_C,'dG',dG_dep/1000,'K_B',Kp_dep, ...
                'yCO',y18_CO,'yCO2',y18_CO2,'res_lnK',res_eq, ...
                'n_recycle',n15_total,'n_cell_in',n18_total, ...
                'X_pass',nC/n15_CO2);
R.cell = struct('E0',E0_SOFC,'E_rev',E_rev_SOFC,'V',V_SOFC, ...
                'V_tn',V_tn,'I',I_SOFC,'W',W_SOFC);
R.sep  = struct('xCO2_15',y15_CO2,'V',E_sep,'I',I_sep,'W',W_sep);
R.elec = struct('W_cell',W_SOFC,'Qel_tot',W_heating_total, ...
                'W_sep',W_sep,'W_comp',W_comp);

R.Qel_tot   = W_heating_total;
R.Qrej_tot  = Q_rejected_total;
R.Q_exhaust = Q_exhaust;
R.P_out     = W_dis_net;
R.nC        = nC;
R.eta_dis   = eta_dis;
end
