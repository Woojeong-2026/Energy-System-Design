%% main_charging_A.m
%  Term Project 2026 - Carbon Oxygen Energy Storage (Noon Energy)
%  담당자 A: Charging System (Stream 1~13)
%
%  이 스크립트 하나로 Charging 관련 결과를 모두 재현한다.
%    [1] 기준 조건 (P = 5 bar, Tcell = 850 degC, TB,ch = 650 degC)
%        - 스트림 결과표 (a), 화학평형/재순환 상세, 열부하/전력 (b), 충전 효율 (c)
%        - 원소수지/에너지수지 잔차
%    [2] 압력 스윕  P = 1, 2, 5, 10 bar           (d)
%    [3] 셀 온도 스윕 Tcell = 750, 800, 850, 900 degC (e)
%  결과 CSV/그림은 results/ 폴더에 저장된다.
%
%  필요: MATLAB R2026a, Cantera 3.0.0 MATLAB toolbox (경로 등록 필요), noon_charging.m

clear; clc; close all;

outdir = 'results';                          % 결과 저장 폴더
if ~exist(outdir, 'dir'), mkdir(outdir); end

MW2W = 1e-6;                                 % W -> MW 환산 계수

%% [1] 기준 조건 -----------------------------------------------------------
P_base     = 5;                              % 시스템 압력 [bar]
Tcell_base = 850;                            % 셀 온도 [degC]
TBch       = 650;                            % 충전 Boudouard 반응기 온도 [degC]

R  = noon_charging(P_base, Tcell_base, TBch);
th = R.th;                                   % Cantera 객체 재사용 (스윕에서 YAML 재로딩 방지)

fprintf('\n===== [1] 기준 조건: P = %g bar, Tcell = %g degC, TB,ch = %g degC =====\n', ...
        P_base, Tcell_base, TBch);

% (a) 스트림 결과표 (Stream 1~13)
fprintf('\n--- (a) 스트림 결과표 ---\n');
disp(R.streams);
writetable(R.streams, fullfile(outdir, 'charging_streams_base.csv'));

% 화학평형 및 재순환 상세
fprintf('--- 화학평형 / 재순환 (2CO <-> C(s) + CO2, T = %g degC) ---\n', TBch);
fprintf('  dG°_B            = %10.2f J/mol\n', R.eq.dG_B);
fprintf('  K_B              = %10.4f\n', R.eq.K_B);
fprintf('  반응기 출구 y_CO  = %10.4f,  y_CO2 = %.4f\n', R.eq.yCO_eq, R.eq.yCO2_eq);
fprintf('  평형 잔차 lnK-lnQ = %10.2e\n', R.eq.res_lnK);
fprintf('  재순환 유량 n9    = %10.4f kmol/s (재순환비 n9/n1 = %.4f)\n', R.eq.n_recycle, R.eq.recycle_ratio);
fprintf('  셀 입구 유량 n4   = %10.4f kmol/s\n', R.eq.n_cell_in);
fprintf('  반응기 1회 통과 CO 전환율 = %.4f\n', R.eq.CO_conv_reactor);

% 셀 전기화학
fprintf('\n--- SOEC ---\n');
fprintf('  E°(Tcell)   = %.4f V\n', R.cell.E0);
fprintf('  E_rev       = %.4f V  (Stream 5: y_CO = %.4f, y_CO2 = %.4f / Stream 12: 순수 O2)\n', ...
        R.cell.E_rev, R.cell.yCO_5, R.cell.yCO2_5);
fprintf('  V_ch        = %.4f V  (= E_rev + 0.10),  열중립 전압 V_tn = %.4f V\n', R.cell.V, R.cell.V_tn);
fprintf('  I           = %.4e A\n', R.cell.I);
fprintf('  W_cell      = %.3f MW\n', R.cell.W*MW2W);

% (b) 열부하 / 회수열 / 방출열 / 전력
d = R.duty;  p = R.power;
fprintf('\n--- (b) 장치 열부하 [MW] (Q > 0: 가열 필요, Q < 0: 냉각) ---\n');
fprintf('  Heater H1          : %9.3f\n', d.Q_H1*MW2W);
fprintf('  rSOC (SOEC)        : %9.3f\n', d.Q_cell*MW2W);
fprintf('  Trim unit          : %9.3f\n', d.Q_trim*MW2W);
fprintf('  Boudouard reactor  : %9.3f\n', d.Q_rx*MW2W);
fprintf('  Recycle cooler     : %9.3f\n', d.Q_c89*MW2W);
fprintf('  C cooler           : %9.3f\n', d.Q_cC*MW2W);
fprintf('  O2 cooler          : %9.3f\n', d.Q_cO2*MW2W);
fprintf('--- 회수열 [MW] ---\n');
fprintf('  HEX 1 (recuperator): %9.3f  (Q_max = %.3f, 최소 온도차 = %.1f K)\n', ...
        d.Q_HEX1*MW2W, d.Q_HEX1_max*MW2W, d.dTmin_HEX1);
fprintf('  Thermal link       : %9.3f  (donor: %s)\n', d.Q_link*MW2W, d.link_donor);
fprintf('--- 방출열 [MW] ---\n');
fprintf('  rSOC %.3f | Trim %.3f | Reactor %.3f | Recycle cooler %.3f | C cooler %.3f | O2 cooler %.3f\n', ...
        d.Qrej_cell*MW2W, d.Qrej_trim*MW2W, d.Qrej_rx*MW2W, d.Qrej_c89*MW2W, d.Qrej_cC*MW2W, d.Qrej_cO2*MW2W);
fprintf('  합계               : %9.3f\n', d.Qrej_tot*MW2W);
fprintf('--- 충전 입력 전력 [MW] ---\n');
fprintf('  셀 전력 W_cell     : %9.3f\n', p.W_cell*MW2W);
fprintf('  전기 가열 H1       : %9.3f\n', p.Qel_H1*MW2W);
fprintf('  전기 가열 셀       : %9.3f\n', p.Qel_cell*MW2W);
fprintf('  전기 가열 Trim     : %9.3f\n', p.Qel_trim*MW2W);
fprintf('  전기 가열 반응기   : %9.3f\n', p.Qel_rx*MW2W);
fprintf('  P_in (합계)        : %9.3f\n', p.P_in*MW2W);

% (c) 충전 효율: eta_ch = N_C * LHV_C / E_in
fprintf('\n--- (c) 충전 효율 ---\n');
fprintf('  eta_ch = (%.1f kmol/s x %.1f MJ/kmol) / %.3f MW = %.4f (%.2f %%)\n', ...
        R.NC_rate, R.LHV_C*1e-6, p.P_in*MW2W, R.eta_ch, 100*R.eta_ch);

% 원소수지 / 에너지수지 잔차
fprintf('\n--- 원소수지 [kmol/s] / 에너지수지 [W] 잔차 ---\n');
disp(R.balance);
writetable(R.balance, fullfile(outdir, 'charging_balance_base.csv'));

%% [2] 압력 스윕 (d) -------------------------------------------------------
P_list = [1 2 5 10];                         % 압력 [bar]
sweepP = run_sweep(P_list, Tcell_base*ones(size(P_list)), TBch, th, outdir, 'P');
fprintf('\n===== [2] 압력 스윕 (Tcell = %g degC) =====\n', Tcell_base);
disp(sweepP);
writetable(sweepP, fullfile(outdir, 'charging_sweep_pressure.csv'));

%% [3] 셀 온도 스윕 (e) ----------------------------------------------------
T_list = [750 800 850 900];                  % 셀 온도 [degC]
sweepT = run_sweep(P_base*ones(size(T_list)), T_list, TBch, th, outdir, 'T');
fprintf('\n===== [3] 셀 온도 스윕 (P = %g bar) =====\n', P_base);
disp(sweepT);
writetable(sweepT, fullfile(outdir, 'charging_sweep_Tcell.csv'));

%% 그림 ---------------------------------------------------------------------
% 압력 영향: 평형 CO 분율, 순환량, 셀 전압, 입력 전력 내역, 충전 효율
figure('Name', 'Charging - pressure sweep', 'Color', 'w');
subplot(2,2,1); plot(sweepP.P_bar, sweepP.yCO_eq, 'o-', 'LineWidth', 1.5); grid on;
xlabel('P [bar]'); ylabel('y_{CO} at reactor outlet [-]'); title('Boudouard equilibrium');
subplot(2,2,2); plot(sweepP.P_bar, sweepP.n_cell_in, 'o-', sweepP.P_bar, sweepP.n_recycle, 's-', 'LineWidth', 1.5); grid on;
xlabel('P [bar]'); ylabel('Molar flow [kmol/s]'); legend('cell inlet (4)', 'recycle (9)', 'Location', 'best'); title('Circulation rate');
subplot(2,2,3); plot(sweepP.P_bar, sweepP.E_rev_V, 'o-', sweepP.P_bar, sweepP.V_ch_V, 's-', 'LineWidth', 1.5); grid on;
xlabel('P [bar]'); ylabel('Voltage [V]'); legend('E_{rev}', 'V_{ch}', 'Location', 'best'); title('Cell voltage');
subplot(2,2,4); yyaxis left; bar(categorical(sweepP.P_bar), [sweepP.W_cell_MW, sweepP.Qel_cell_MW, sweepP.Qel_H1_MW], 'stacked');
ylabel('Input power [MW]'); yyaxis right; plot(categorical(sweepP.P_bar), 100*sweepP.eta_ch, 'ko-', 'LineWidth', 1.5);
ylabel('\eta_{ch} [%]'); xlabel('P [bar]'); legend('W_{cell}', 'Q_{el,cell}', 'Q_{el,H1}', '\eta_{ch}', 'Location', 'southoutside', 'Orientation', 'horizontal');
title('Input power and efficiency');
saveas(gcf, fullfile(outdir, 'charging_sweep_pressure.png'));

% 셀 온도 영향
figure('Name', 'Charging - cell temperature sweep', 'Color', 'w');
subplot(1,3,1); plot(sweepT.Tcell_C, sweepT.E_rev_V, 'o-', sweepT.Tcell_C, sweepT.V_ch_V, 's-', sweepT.Tcell_C, sweepT.V_tn_V, '^-', 'LineWidth', 1.5); grid on;
xlabel('T_{cell} [°C]'); ylabel('Voltage [V]'); legend('E_{rev}', 'V_{ch}', 'V_{tn}', 'Location', 'best'); title('Cell voltage');
subplot(1,3,2); plot(sweepT.Tcell_C, sweepT.Qel_H1_MW, 'o-', sweepT.Tcell_C, sweepT.Qel_cell_MW, 's-', sweepT.Tcell_C, sweepT.Qrej_trim_MW, '^-', 'LineWidth', 1.5); grid on;
xlabel('T_{cell} [°C]'); ylabel('Heat [MW]'); legend('Q_{el,H1}', 'Q_{el,cell}', 'Q_{rej,trim}', 'Location', 'best'); title('Heating / rejection');
subplot(1,3,3); yyaxis left; plot(sweepT.Tcell_C, sweepT.P_in_MW, 'o-', 'LineWidth', 1.5); ylabel('P_{in} [MW]');
yyaxis right; plot(sweepT.Tcell_C, 100*sweepT.eta_ch, 's-', 'LineWidth', 1.5); ylabel('\eta_{ch} [%]');
xlabel('T_{cell} [°C]'); grid on; title('Input power and efficiency');
saveas(gcf, fullfile(outdir, 'charging_sweep_Tcell.png'));

fprintf('\n완료: 결과 파일은 "%s" 폴더에 저장되었습니다.\n', outdir);

%% 로컬 함수 ----------------------------------------------------------------
function tab = run_sweep(P_list, T_list, TBch, th, outdir, tag)
% 조건 목록에 대해 noon_charging 을 반복 실행하고 요약표를 만든다.
% 각 조건의 전체 스트림표와 수지 잔차도 CSV 로 저장한다 (과제 (d), (e): (a)~(c) 반복).
    N = numel(P_list);
    P_bar = zeros(N,1); Tcell_C = zeros(N,1);
    K_B = zeros(N,1); yCO_eq = zeros(N,1); n_recycle = zeros(N,1); n_cell_in = zeros(N,1);
    E_rev_V = zeros(N,1); V_ch_V = zeros(N,1); V_tn_V = zeros(N,1);
    T3_C = zeros(N,1); T8_C = zeros(N,1);
    W_cell_MW = zeros(N,1); Qel_H1_MW = zeros(N,1); Qel_cell_MW = zeros(N,1);
    Q_rx_MW = zeros(N,1); Q_HEX1_MW = zeros(N,1); Q_link_MW = zeros(N,1);
    Qrej_trim_MW = zeros(N,1); Qrej_tot_MW = zeros(N,1);
    P_in_MW = zeros(N,1); eta_ch = zeros(N,1); res_E_sys_W = zeros(N,1);
    for k = 1:N
        r = noon_charging(P_list(k), T_list(k), TBch, th);
        P_bar(k) = P_list(k);           Tcell_C(k) = T_list(k);
        K_B(k) = r.eq.K_B;              yCO_eq(k) = r.eq.yCO_eq;          % 평형
        n_recycle(k) = r.eq.n_recycle;  n_cell_in(k) = r.eq.n_cell_in;    % 순환량 [kmol/s]
        E_rev_V(k) = r.cell.E_rev;      V_ch_V(k) = r.cell.V;  V_tn_V(k) = r.cell.V_tn;
        T3_C(k) = r.S(3).T - 273.15;    T8_C(k) = r.S(8).T - 273.15;      % HEX1 출구 온도
        W_cell_MW(k)   = r.power.W_cell*1e-6;
        Qel_H1_MW(k)   = r.power.Qel_H1*1e-6;
        Qel_cell_MW(k) = r.power.Qel_cell*1e-6;
        Q_rx_MW(k)     = r.duty.Q_rx*1e-6;
        Q_HEX1_MW(k)   = r.duty.Q_HEX1*1e-6;
        Q_link_MW(k)   = r.duty.Q_link*1e-6;
        Qrej_trim_MW(k) = r.duty.Qrej_trim*1e-6;
        Qrej_tot_MW(k) = r.duty.Qrej_tot*1e-6;
        P_in_MW(k)     = r.power.P_in*1e-6;
        eta_ch(k)      = r.eta_ch;
        res_E_sys_W(k) = r.balance.res_E_W(end);                          % 전체 에너지수지 잔차
        % 조건별 전체 스트림표 / 수지 잔차 저장
        name = sprintf('%s_P%gbar_T%gC', tag, P_list(k), T_list(k));
        writetable(r.streams, fullfile(outdir, ['charging_streams_' name '.csv']));
        writetable(r.balance, fullfile(outdir, ['charging_balance_' name '.csv']));
    end
    tab = table(P_bar, Tcell_C, K_B, yCO_eq, n_recycle, n_cell_in, E_rev_V, V_ch_V, V_tn_V, ...
                T3_C, T8_C, W_cell_MW, Qel_H1_MW, Qel_cell_MW, Q_rx_MW, Q_HEX1_MW, Q_link_MW, ...
                Qrej_trim_MW, Qrej_tot_MW, P_in_MW, eta_ch, res_E_sys_W);
end
