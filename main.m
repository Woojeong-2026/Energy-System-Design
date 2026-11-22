%% main.m
%  Energy System Design (38609) 2026 Fall - Term Project
%  Carbon Oxygen Energy Storage Inspired by Noon Energy
%
%  팀원 3 (통합 담당): 물성 도우미 함수 + 전체 시스템 통합 + 민감도 분석
%
%  이 스크립트 하나로 보고서의 모든 결과를 재현한다 (과제: "single main MATLAB script").
%    [1] 기준 조건 (P = 5 bar, Tcell = 850 degC)
%        (a) Stream 1~29 표, 원소/에너지수지 잔차
%        (b) 충전 입력 전력, 방전 순출력, 열부하/회수열/방출열
%        (c) 충전/방전/왕복 효율, 운전 시간 맞추기
%    [2] (d) 압력 스윕     P = 1, 2, 5, 10 bar
%    [3] (e) 셀 온도 스윕  Tcell = 750, 800, 850, 900 degC
%    [4] (f) 손실 내역과 민감도 분석 (효율 개선 방향)
%    [5] 그림
%  결과 CSV/PNG 는 results/ 폴더에 저장된다.
%
%  필요: MATLAB R2026a, Cantera 3.0.0 MATLAB toolbox (경로 등록)
%
%  폴더 구성
%    noon_config.m : 과제 지정 조건 (모든 입력값)
%    thermo/       : Cantera 물성 도우미 함수 (팀원 3)
%    system/       : 공통 장치 모델, 수지표, 왕복 효율 통합 (팀원 3)
%    run_charging_analysis.m    : 충전 모드 Stream 1~13 (팀원 1)
%    run_discharging_analysis.m : 방전 모드 Stream 14~29 (팀원 2)
%      두 파일 모두 R = run_xxx_analysis(cfg, th) 로 부르면 화면 출력 없이 결과만 돌려준다.
%      인자 없이 단독 실행하면 각 모드만 따로 분석한다.
%
%  부호 규약 (전체 공통)
%    Q > 0 : 장치로 들어가는 열 (가열 필요),  Q < 0 : 장치에서 빼내는 열
%    W > 0 : 장치로 들어가는 전력
%    보고용 P_in (충전 입력), P_out (방전 순출력) 은 모두 양수로 표기

clear; clc; close all;

here = fileparts(mfilename('fullpath'));         % 이 스크립트가 있는 폴더
if isempty(here), here = pwd; end                % 섹션 단위로 실행할 때는 현재 폴더 사용
addpath(here, fullfile(here, 'thermo'), fullfile(here, 'system'));   % 최상위(모드 파일) + 공통 함수 폴더

outdir = fullfile(here, 'results');              % 결과 저장 폴더
if ~exist(outdir, 'dir'), mkdir(outdir); end

th  = thermo_init();                             % Cantera 물성 객체 (한 번만 생성)
cfg = noon_config();                             % 과제 기준 조건
MW  = 1e-6;                                      % W -> MW 환산 계수

fprintf('Cantera MATLAB toolbox 문법: %s\n', th.api);

%% [1] 기준 조건 ------------------------------------------------------------
fprintf('\n================ [1] 기준 조건: P = %g bar, Tcell = %g degC ================\n', ...
        cfg.P_bar, cfg.Tcell_C);
B   = run_roundtrip(cfg, th);                    % 충전 + 방전 + 왕복 효율
ch  = B.ch;                                      % 충전 모드 결과 (Stream 1~13)
dis = B.dis;                                     % 방전 모드 결과 (Stream 14~29)

% ---- (a) 스트림 표: 몰유량, 질량유량, 조성, 온도, 압력 ----
streams = [ch.streams; dis.streams];             % Stream 1~29 를 한 표로
fprintf('\n--- (a) Stream 1~29 ---\n');
disp(streams);
writetable(streams, fullfile(outdir, 'base_streams.csv'));

% ---- (a) 원소수지 / 에너지수지 잔차 (장치별 + 전체 시스템) ----
balance = [ch.balance; dis.balance];
fprintf('--- (a) 원소수지 [kmol/s] / 에너지수지 [W] 잔차 ---\n');
disp(balance);
writetable(balance, fullfile(outdir, 'base_balance.csv'));
fprintf('  최대 원소 잔차   = %.2e kmol/s\n', ...
        max(abs([balance.res_C; balance.res_O; balance.res_N])));
fprintf('  최대 에너지 잔차 = %.2e W  (충전 입력 전력 대비 %.1e)\n', ...
        max(abs(balance.res_E_W)), max(abs(balance.res_E_W))/ch.P_in);

% ---- 화학평형 / 재순환 ----
fprintf('\n--- Boudouard 평형 (K_B = (p_CO2/p°)/(p_CO/p°)^2) 과 순환량 ---\n');
fprintf('  충전 %g degC: K_B = %.4f, y_CO(7)  = %.4f, 셀 입구 n4  = %.4f kmol/s, 잔차 %.1e\n', ...
        ch.eq.T_C, ch.eq.K_B, ch.eq.yCO, ch.eq.n_cell_in, ch.eq.res_lnK);
fprintf('  방전 %g degC: K_B = %.4f, y_CO(18) = %.4f, 셀 입구 n19 = %.4f kmol/s, 잔차 %.1e\n', ...
        dis.eq.T_C, dis.eq.K_B, dis.eq.yCO, dis.eq.n_cell_in, dis.eq.res_lnK);

% ---- 셀 / 분리기 전기화학 ----
fprintf('\n--- 셀 전압 ---\n');
fprintf('  SOEC: E° = %.4f V, E_rev = %.4f V, V_ch  = %.4f V, V_tn = %.4f V, I = %.4e A\n', ...
        ch.cell.E0, ch.cell.E_rev, ch.cell.V, ch.cell.V_tn, ch.cell.I);
fprintf('  SOFC: E° = %.4f V, E_rev = %.4f V, V_dis = %.4f V, V_tn = %.4f V, I = %.4e A\n', ...
        dis.cell.E0, dis.cell.E_rev, dis.cell.V, dis.cell.V_tn, dis.cell.I);
fprintf('  CO2 분리기: x_CO2,15 = %.4f, V_sep,rev = %.5f V, I_sep = %.4e A\n', ...
        dis.sep.xCO2_15, dis.sep.V, dis.sep.I);

% ---- (b) 열부하 / thermal link / 전기 가열 / 방출열 ----
fprintf('\n--- (b) 충전 모드 열부하 [MW] ---\n');   disp(ch.heat);
fprintf('--- (b) 방전 모드 열부하 [MW] ---\n');     disp(dis.heat);
writetable(ch.heat,  fullfile(outdir, 'base_heat_charging.csv'));
writetable(dis.heat, fullfile(outdir, 'base_heat_discharging.csv'));

fprintf('--- (b) thermal link ---\n');
fprintf('  충전: %.3f MW (%s -> %s)\n', ch.link.Q*MW,  ch.link.from,  ch.link.to);
fprintf('  방전: %.3f MW (%s -> %s)\n', dis.link.Q*MW, dis.link.from, dis.link.to);

fprintf('--- (b) 열교환기 (recuperator) ---\n');
hexAll = [ch.hex; dis.hex];
for k = 1:numel(hexAll)
    h = hexAll(k);
    fprintf('  %s: Q = %8.3f MW (Q_max = %8.3f MW), 고온측 %.1f -> %.1f degC, 저온측 %.1f -> %.1f degC, 최소 온도차 %.1f K\n', ...
            h.name, h.Q*MW, h.Qmax*MW, h.Th_in - 273.15, h.Th_out - 273.15, ...
            h.Tc_in - 273.15, h.Tc_out - 273.15, h.dTmin);
    if h.dTmin <= 0                              % 열이 고온 -> 저온으로 흐르는지 확인
        warning('%s 내부에서 온도가 역전됩니다 (최소 온도차 %.2f K).', h.name, h.dTmin);
    end
end

% ---- (b) 전력 내역 ----
fprintf('\n--- (b) 충전 입력 전력 [MW] ---\n');
fprintf('  셀 전력 W_cell         : %9.3f\n', ch.elec.W_cell*MW);
fprintf('  전기 가열 합계         : %9.3f\n', ch.elec.Qel_tot*MW);
fprintf('  P_in                   : %9.3f\n', ch.P_in*MW);
fprintf('--- (b) 방전 순출력 [MW] ---\n');
fprintf('  셀 전력 W_cell         : %9.3f\n', dis.elec.W_cell*MW);
fprintf('  CO2 분리기 W_sep       : %9.3f\n', -dis.elec.W_sep*MW);
fprintf('  공기 압축 W_comp       : %9.3f\n', -dis.elec.W_comp*MW);
fprintf('  전기 가열 합계         : %9.3f\n', -dis.elec.Qel_tot*MW);
fprintf('  P_out                  : %9.3f\n', dis.P_out*MW);
fprintf('  방출열 합계: 충전 %.3f MW, 방전 %.3f MW, 배기(29) 현열 %.3f MW\n', ...
        ch.Qrej_tot*MW, dis.Qrej_tot*MW, dis.Q_exhaust*MW);

% ---- (c) 효율 ----
% 같은 탄소 재고 N_C 로 맞춘다: N_C = nC_ch t_ch = nC_dis t_dis
%   eta_ch  = N_C LHV_C / E_in  = nC_ch LHV_C / P_in
%   eta_dis = E_out / (N_C LHV_C) = P_out / (nC_dis LHV_C)
%   eta_RT  = E_out / E_in = eta_ch eta_dis
fprintf('\n--- (c) 효율 (LHV_C = %.1f kJ/mol) ---\n', cfg.LHV_C*1e-6);
fprintf('  운전 시간 비 t_dis/t_ch = %.4f  (N_C,stored = N_C,consumed)\n', B.t_ratio);
fprintf('  eta_ch  = %.4f (%.2f %%)\n', B.eta_ch,  100*B.eta_ch);
fprintf('  eta_dis = %.4f (%.2f %%)\n', B.eta_dis, 100*B.eta_dis);
fprintf('  eta_RT  = %.4f (%.2f %%)   [E_out/E_in 직접 계산: %.4f]\n', ...
        B.eta_RT, 100*B.eta_RT, B.eta_RT_energy);

% ---- (c) 맞추지 않은 전력의 비가 왕복 효율이 아닌 이유 (수치 예) ----
% 방전 탄소 소비 속도만 2 kmol/s 로 바꾸면 P_out 이 2배가 되어 P_out/P_in 도 2배가 된다.
% 그러나 같은 탄소 재고는 절반 시간에 소진되므로(t_dis = t_ch/2) 에너지 비 E_out/E_in 은 그대로다.
% 즉 전력의 비는 두 모드의 탄소 처리 속도를 어떻게 잡느냐에 따라 달라지는 값이고,
% 왕복 효율은 같은 N_C 에 대한 에너지의 비로만 정의된다.
B2 = run_roundtrip(noon_config('nC_dis', 2.0), th);
fprintf('  [예] 방전 속도 2 kmol/s: P_out/P_in = %.4f (전력 비),  t_dis/t_ch = %.2f,  E_out/E_in = %.4f (왕복 효율)\n', ...
        B2.dis.P_out/B2.ch.P_in, B2.t_ratio, B2.eta_RT_energy);

%% [2] (d) 압력 스윕 --------------------------------------------------------
P_list = [1 2 5 10];                             % 시스템 압력 [bar]
sweepP = table();
for k = 1:numel(P_list)
    C = run_roundtrip(noon_config('P_bar', P_list(k)), th);
    sweepP = [sweepP; summary_row(C)];                                   %#ok<AGROW>
    save_case(C, outdir, sprintf('P%gbar', P_list(k)));                  % (a)~(c) 표 저장
end
fprintf('\n================ [2] (d) 압력 스윕 (Tcell = %g degC) ================\n', cfg.Tcell_C);
disp(sweepP);
writetable(sweepP, fullfile(outdir, 'sweep_pressure.csv'));

%% [3] (e) 셀 온도 스윕 (반응기 온도는 그대로) ------------------------------
T_list = [750 800 850 900];                      % 셀 온도 [degC]
sweepT = table();
for k = 1:numel(T_list)
    C = run_roundtrip(noon_config('Tcell_C', T_list(k)), th);
    sweepT = [sweepT; summary_row(C)];                                   %#ok<AGROW>
    save_case(C, outdir, sprintf('Tcell%gC', T_list(k)));
end
fprintf('\n================ [3] (e) 셀 온도 스윕 (P = %g bar) ================\n', cfg.P_bar);
disp(sweepT);
writetable(sweepT, fullfile(outdir, 'sweep_Tcell.csv'));
% 셀 온도가 방전 반응기 온도(800 degC) + 20 K 보다 낮으면 thermal link 가 끊겨
% 가스화 열을 전부 전기로 공급해야 하므로 방전 순출력이 음수가 될 수 있다.
if any(sweepT.P_out_MW <= 0)
    fprintf('  주의: P_out <= 0 인 조건이 있습니다 (순발전 불가). Qlink_dis_MW 열을 확인하세요.\n');
end

%% [4] (f) 손실 내역과 민감도 분석 ------------------------------------------
fprintf('\n================ [4] (f) 손실 내역과 민감도 ================\n');

% ---- 손실 내역: 충전 입력 전기 1 단위가 어디로 사라지는가 (탄소 1 kmol 기준) ----
%   e_in - e_out = sum(방출열) + 배기 현열
fprintf('--- 기준 조건 손실 내역 (e_in = %.2f, e_out = %.2f MJ/kmol-C) ---\n', ...
        B.e_in*MW, B.e_out*MW);
disp(B.loss);
fprintf('  사이클 에너지 닫힘 잔차 = %.2e J/kmol-C\n', B.closure);
writetable(B.loss, fullfile(outdir, 'base_loss_breakdown.csv'));

% ---- 민감도: 한 번에 변수 하나만 바꾼다 (나머지는 기준 조건) ----
%   {설정 이름, 시험 값들}
sens_def = { ...
    'eff_hex',   [0.75 0.95];   ...  % 열교환기 유용도
    'lambda_O2', [2 4];         ...  % 공기 과잉률 (양론 대비 O2 배수)
    'dV_loss',   [0.05 0.15];   ...  % 셀 전압 손실 [V]
    'U_sofc',    [0.6 0.9];     ...  % SOFC 연료 이용률
    'X_soec',    [0.5 0.9];     ...  % SOEC CO2 전환율
    'TBdis_C',   [750 830];     ...  % 방전 반응기 온도 [degC]
    'TBch_C',    [600 700 870]; ...  % 충전 반응기 온도 [degC] (870 = Tcell + 20 K, link 작동)
    'f_link',    [0.70 1.00]};       % thermal link 회수 비율
Parameter = strings(0, 1);  Value = zeros(0, 1);
eta_ch_pct = zeros(0, 1);  eta_dis_pct = zeros(0, 1);  eta_RT_pct = zeros(0, 1);
d_eta_RT_pt = zeros(0, 1);
for i = 1:size(sens_def, 1)
    for v = sens_def{i, 2}
        C = run_roundtrip(noon_config(sens_def{i, 1}, v), th);
        Parameter(end + 1, 1)   = sens_def{i, 1};                        %#ok<SAGROW>
        Value(end + 1, 1)       = v;                                     %#ok<SAGROW>
        eta_ch_pct(end + 1, 1)  = 100*C.eta_ch;                          %#ok<SAGROW>
        eta_dis_pct(end + 1, 1) = 100*C.eta_dis;                         %#ok<SAGROW>
        eta_RT_pct(end + 1, 1)  = 100*C.eta_RT;                          %#ok<SAGROW>
        d_eta_RT_pt(end + 1, 1) = 100*(C.eta_RT - B.eta_RT);             %#ok<SAGROW> % 기준 대비 변화 [%p]
    end
end
sens = table(Parameter, Value, eta_ch_pct, eta_dis_pct, eta_RT_pct, d_eta_RT_pt);
fprintf('--- 민감도 (기준 eta_RT = %.2f %%) ---\n', 100*B.eta_RT);
disp(sens);
writetable(sens, fullfile(outdir, 'sensitivity.csv'));

% ---- 개선안 조합 예시: 민감도가 큰 항목을 함께 적용 ----
%   열교환기 유용도 0.95, 공기 과잉률 2, 충전 반응기 870 degC (thermal link 작동)
cfg_imp = noon_config('eff_hex', 0.95, 'lambda_O2', 2, 'TBch_C', 870);
I1 = run_roundtrip(cfg_imp, th);
fprintf('--- 개선안 조합 (eff_hex 0.95, lambda_O2 2, TB,ch 870 degC) ---\n');
fprintf('  eta_ch %.2f -> %.2f %%,  eta_dis %.2f -> %.2f %%,  eta_RT %.2f -> %.2f %%\n', ...
        100*B.eta_ch, 100*I1.eta_ch, 100*B.eta_dis, 100*I1.eta_dis, 100*B.eta_RT, 100*I1.eta_RT);
save_case(I1, outdir, 'improved');

%% [5] 그림 -----------------------------------------------------------------
% (d) 압력 영향: 평형, 순환량, 셀 전압, 압축 전력, 효율
figure('Name', '(d) Pressure sweep', 'Color', 'w', 'Position', [100 100 1100 650]);
subplot(2, 3, 1);
plot(sweepP.P_bar, sweepP.yCO_ch, 'o-', sweepP.P_bar, sweepP.yCO_dis, 's-', 'LineWidth', 1.5); grid on;
xlabel('P [bar]'); ylabel('y_{CO} at reactor outlet [-]'); title('Boudouard equilibrium');
legend('charging (650 °C)', 'discharging (800 °C)', 'Location', 'best');
subplot(2, 3, 2);
plot(sweepP.P_bar, sweepP.n_cell_ch, 'o-', sweepP.P_bar, sweepP.n_cell_dis, 's-', 'LineWidth', 1.5); grid on;
xlabel('P [bar]'); ylabel('Cell fuel inlet flow [kmol/s]'); title('Circulation rate');
legend('stream 4', 'stream 19', 'Location', 'best');
subplot(2, 3, 3);
plot(sweepP.P_bar, sweepP.V_ch, 'o-', sweepP.P_bar, sweepP.V_dis, 's-', ...
     sweepP.P_bar, sweepP.Erev_ch, 'o--', sweepP.P_bar, sweepP.Erev_dis, 's--', 'LineWidth', 1.5); grid on;
xlabel('P [bar]'); ylabel('Voltage [V]'); title('Cell voltage');
legend('V_{ch}', 'V_{dis}', 'E_{rev,ch}', 'E_{rev,dis}', 'Location', 'best');
subplot(2, 3, 4);
plot(sweepP.P_bar, sweepP.Wcomp_MW, 'o-', sweepP.P_bar, sweepP.Wsep_MW, 's-', 'LineWidth', 1.5); grid on;
xlabel('P [bar]'); ylabel('Power [MW]'); title('Compression / separation power');
legend('air compressor', 'CO_2 separator', 'Location', 'best');
subplot(2, 3, 5);
plot(sweepP.P_bar, sweepP.P_in_MW, 'o-', sweepP.P_bar, sweepP.P_out_MW, 's-', 'LineWidth', 1.5); grid on;
xlabel('P [bar]'); ylabel('Power [MW]'); title('Charge input / discharge net output');
legend('P_{in}', 'P_{out}', 'Location', 'best');
subplot(2, 3, 6);
plot(sweepP.P_bar, sweepP.eta_ch_pct, 'o-', sweepP.P_bar, sweepP.eta_dis_pct, 's-', ...
     sweepP.P_bar, sweepP.eta_RT_pct, '^-', 'LineWidth', 1.5); grid on;
xlabel('P [bar]'); ylabel('Efficiency [%]'); title('Efficiency');
legend('\eta_{ch}', '\eta_{dis}', '\eta_{RT}', 'Location', 'best');
saveas(gcf, fullfile(outdir, 'fig_sweep_pressure.png'));

% (e) 셀 온도 영향: 전압, 전기 가열, thermal link, 효율
figure('Name', '(e) Cell temperature sweep', 'Color', 'w', 'Position', [120 120 1100 650]);
subplot(2, 2, 1);
plot(sweepT.Tcell_C, sweepT.V_ch, 'o-', sweepT.Tcell_C, sweepT.V_dis, 's-', ...
     sweepT.Tcell_C, sweepT.Erev_ch, 'o--', sweepT.Tcell_C, sweepT.Erev_dis, 's--', 'LineWidth', 1.5); grid on;
xlabel('T_{cell} [°C]'); ylabel('Voltage [V]'); title('Cell voltage');
legend('V_{ch}', 'V_{dis}', 'E_{rev,ch}', 'E_{rev,dis}', 'Location', 'best');
subplot(2, 2, 2);
plot(sweepT.Tcell_C, sweepT.Qel_ch_MW, 'o-', sweepT.Tcell_C, sweepT.Qel_dis_MW, 's-', ...
     sweepT.Tcell_C, sweepT.Qlink_dis_MW, '^-', 'LineWidth', 1.5); grid on;
xlabel('T_{cell} [°C]'); ylabel('Heat [MW]'); title('Electric heating and thermal link');
legend('Q_{el} charging', 'Q_{el} discharging', 'thermal link (discharging)', 'Location', 'best');
subplot(2, 2, 3);
plot(sweepT.Tcell_C, sweepT.P_in_MW, 'o-', sweepT.Tcell_C, sweepT.P_out_MW, 's-', 'LineWidth', 1.5); grid on;
xlabel('T_{cell} [°C]'); ylabel('Power [MW]'); title('Charge input / discharge net output');
legend('P_{in}', 'P_{out}', 'Location', 'best');
subplot(2, 2, 4);
plot(sweepT.Tcell_C, sweepT.eta_ch_pct, 'o-', sweepT.Tcell_C, sweepT.eta_dis_pct, 's-', ...
     sweepT.Tcell_C, sweepT.eta_RT_pct, '^-', 'LineWidth', 1.5); grid on;
xlabel('T_{cell} [°C]'); ylabel('Efficiency [%]'); title('Efficiency');
legend('\eta_{ch}', '\eta_{dis}', '\eta_{RT}', 'Location', 'best');
saveas(gcf, fullfile(outdir, 'fig_sweep_Tcell.png'));

% (f) 손실 내역 (기준 조건): 충전 입력 전기 대비 비율
figure('Name', '(f) Loss breakdown', 'Color', 'w', 'Position', [140 140 900 500]);
barh(B.loss.Share_of_Ein_pct); grid on;
set(gca, 'YTick', 1:height(B.loss), 'YTickLabel', B.loss.Item + " [" + B.loss.Mode + "]", ...
         'YDir', 'reverse');
xlabel('Loss / charge input electricity [%]');
title(sprintf('Loss breakdown at base case (\\eta_{RT} = %.1f %%)', 100*B.eta_RT));
saveas(gcf, fullfile(outdir, 'fig_loss_breakdown.png'));

% (f) 민감도: 기준 대비 왕복 효율 변화 [%p]
figure('Name', '(f) Sensitivity', 'Color', 'w', 'Position', [160 160 900 600]);
barh(sens.d_eta_RT_pt); grid on;
set(gca, 'YTick', 1:height(sens), 'YTickLabel', sens.Parameter + " = " + string(sens.Value), ...
         'YDir', 'reverse', 'TickLabelInterpreter', 'none');
xlabel('\Delta\eta_{RT} from base case [%-point]');
title(sprintf('One-at-a-time sensitivity (base \\eta_{RT} = %.1f %%)', 100*B.eta_RT));
saveas(gcf, fullfile(outdir, 'fig_sensitivity.png'));

fprintf('\n완료: 결과 파일은 "%s" 폴더에 저장되었습니다.\n', outdir);

%% 로컬 함수 ----------------------------------------------------------------
function row = summary_row(C)
% 한 조건의 핵심 결과를 한 행짜리 표로 만든다 (스윕 요약표용).
%   과제 (d) 가 묻는 항목: Boudouard 평형, 순환량, 셀 전압, 압축 전력, 효율
    ch = C.ch;  dis = C.dis;
    s.P_bar        = C.cfg.P_bar;                % 압력 [bar]
    s.Tcell_C      = C.cfg.Tcell_C;              % 셀 온도 [degC]
    s.K_ch         = ch.eq.K_B;                  % 충전 반응기 평형상수
    s.yCO_ch       = ch.eq.yCO;                  % Stream 7 CO 몰분율
    s.n_cell_ch    = ch.eq.n_cell_in;            % Stream 4 유량 [kmol/s]
    s.Erev_ch      = ch.cell.E_rev;              % SOEC 가역전압 [V]
    s.V_ch         = ch.cell.V;                  % SOEC 운전전압 [V]
    s.Wcell_ch_MW  = ch.elec.W_cell*1e-6;        % SOEC 전력 [MW]
    s.Qel_ch_MW    = ch.Qel_tot*1e-6;            % 충전 전기 가열 [MW]
    s.Qlink_ch_MW  = ch.link.Q*1e-6;             % 충전 thermal link [MW]
    s.P_in_MW      = ch.P_in*1e-6;               % 충전 입력 전력 [MW]
    s.K_dis        = dis.eq.K_B;                 % 방전 반응기 평형상수
    s.yCO_dis      = dis.eq.yCO;                 % Stream 18 CO 몰분율
    s.n_cell_dis   = dis.eq.n_cell_in;           % Stream 19 유량 [kmol/s]
    s.Erev_dis     = dis.cell.E_rev;             % SOFC 가역전압 [V]
    s.V_dis        = dis.cell.V;                 % SOFC 운전전압 [V]
    s.Wcell_dis_MW = dis.elec.W_cell*1e-6;       % SOFC 전력 [MW]
    s.Wcomp_MW     = dis.elec.W_comp*1e-6;       % 공기 압축 전력 [MW]
    s.Wsep_MW      = dis.elec.W_sep*1e-6;        % CO2 분리기 전력 [MW]
    s.Qel_dis_MW   = dis.Qel_tot*1e-6;           % 방전 전기 가열 [MW]
    s.Qlink_dis_MW = dis.link.Q*1e-6;            % 방전 thermal link [MW]
    s.P_out_MW     = dis.P_out*1e-6;             % 방전 순출력 [MW]
    s.eta_ch_pct   = 100*C.eta_ch;               % 충전 효율 [%]
    s.eta_dis_pct  = 100*C.eta_dis;              % 방전 효율 [%]
    s.eta_RT_pct   = 100*C.eta_RT;               % 왕복 효율 [%]
    % 수지 잔차 (모든 조건에서 닫히는지 확인)
    bal = [ch.balance; dis.balance];
    s.max_res_atom = max(abs([bal.res_C; bal.res_O; bal.res_N]));   % [kmol/s]
    s.max_res_E_W  = max(abs(bal.res_E_W));                         % [W]
    row = struct2table(s);
end

function save_case(C, outdir, tag)
% 한 조건의 (a)~(c) 표를 CSV 로 저장한다 (과제 (d), (e): "Repeat parts (a)-(c)").
    writetable([C.ch.streams; C.dis.streams], fullfile(outdir, ['streams_' tag '.csv']));
    writetable([C.ch.balance; C.dis.balance], fullfile(outdir, ['balance_' tag '.csv']));
    writetable([C.ch.heat; C.dis.heat],       fullfile(outdir, ['heat_' tag '.csv']));
    writetable(C.loss,                        fullfile(outdir, ['loss_' tag '.csv']));
end
