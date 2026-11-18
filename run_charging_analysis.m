%% 충전 담당 실행 파일: 과제 (a)~(e) 충전 결과, (f) 개선 비교.
% MATLAB R2026a + Cantera 3.0.0. 저장소 최상위에서 이 파일을 Run.
% charging_student.m + 정우정 팀원의 thermo/, system/, noon_config만 사용.
% 기존 modes/charging_mode.m, Python 참조 계산, 오류 주입 시험은 호출하지 않는다.
clear; clc; close all;
projectRoot=fileparts(mfilename('fullpath'));
addpath(projectRoot,fullfile(projectRoot,'thermo'),fullfile(projectRoot,'system'),'-begin');
clear charging_student
assert(strcmp(which('charging_student'),fullfile(projectRoot,'charging_student.m')), ...
    '다른 폴더의 charging_student가 선택되었습니다. 이 저장소를 Current Folder로 선택하세요.');
fprintf('Charging calculation: %s\n',which('charging_student'));
th=thermo_init();
outdir=fullfile(projectRoot,'results','charging_student');
if ~isfolder(outdir), mkdir(outdir); end

% 기준 조건과 중복되는 P=5, Tcell=850은 base 결과를 재사용한다.
conditions=[5 850;1 850;2 850;10 850;5 750;5 800;5 900];
tags=["base";"P1";"P2";"P10";"T750";"T800";"T900"];
chargingSummary=table();
for k=1:size(conditions,1)
    cfg=noon_config('P_bar',conditions(k,1),'Tcell_C',conditions(k,2));
    R=charging_student(cfg,th);
    chargingSummary=[chargingSummary; result_row(R,tags(k))]; %#ok<AGROW>
    folder=fullfile(outdir,tags(k)); if ~isfolder(folder), mkdir(folder); end
    writetable(R.streams,fullfile(folder,'streams.csv'));       % (a): 1~13번
    writetable(R.balance,fullfile(folder,'balances.csv'));      % (a): 장치/전체 잔차
    writetable(R.heat,fullfile(folder,'heat_duties.csv'));       % (b): 열부하/회수/가열/방출
    hx=R.hex;
    hexTable=table(hx.Q/1e6,hx.Qmax/1e6,hx.Th_in-273.15,hx.Th_out-273.15, ...
        hx.Tc_in-273.15,hx.Tc_out-273.15,hx.dTmin, ...
        'VariableNames',{'Q_MW','Qmax_MW','HotIn_C','HotOut_C','ColdIn_C','ColdOut_C','MinDT_K'});
    writetable(hexTable,fullfile(folder,'heat_exchanger.csv')); % (b): HEX 열량도 별도 표
    writetable(chargingSummary(end,:),fullfile(folder,'performance.csv')); % (b),(c)
    if k==1, baseResults=R; end
end
pressureResults=chargingSummary([2 3 1 4],:);
temperatureResults=chargingSummary([5 6 1 7],:);
writetable(chargingSummary,fullfile(outdir,'summary.csv'));
writetable(pressureResults,fullfile(outdir,'pressure_sweep.csv'));
writetable(temperatureResults,fullfile(outdir,'temperature_sweep.csv'));
fprintf('\n(a) 기준 Stream 1~13 및 수지 잔차\n');
disp(baseResults.streams); disp(baseResults.balance);
fprintf('\n(b) 기준 열부하 [MW]\n'); disp(baseResults.heat);
fprintf('\n(b,c) 전력/효율 및 (d,e) 조건 비교\n'); disp(chargingSummary);

% (f) 아래 두 비교만 기준 가정을 변경한 개선 검토. (a)~(e)와 구분한다.
improvements=result_row(baseResults,"base");
improvements=[improvements;result_row(charging_student(noon_config('eff_hex',0.95),th),"HEX_0.95")];
improvements=[improvements;result_row(charging_student(noon_config('TBch_C',870),th),"Reactor_870C")];
writetable(improvements,fullfile(outdir,'improvement_comparison.csv'));
fprintf('\n(f) 변경된 가정에 따른 충전 개선 비교\n'); disp(improvements);

% 압력과 셀 온도 각각의 경향. 실제 보고 수치는 CSV 표를 사용한다.
make_plot(pressureResults,'P_bar','Pressure [bar]',fullfile(outdir,'pressure_sweep.png'));
make_plot(temperatureResults,'Tcell_C','Cell temperature [degC]',fullfile(outdir,'temperature_sweep.png'));
fprintf('\n결과 폴더: %s\n',outdir);
fprintf('전체 방전/왕복 효율은 팀 통합 계산에서 별도로 구합니다.\n');

function row=result_row(R,label)
% 과제의 전력·효율·평형·순환량과 검산 잔차를 한 행으로 정리.
    s.Case=string(label); s.P_bar=R.cfg.P_bar; s.Tcell_C=R.cfg.Tcell_C;
    s.TBch_C=R.cfg.TBch_C; s.eff_hex=R.cfg.eff_hex;
    s.K_B=R.eq.K_B; s.yCO7=R.eq.yCO; s.Recycle_kmol_s=R.eq.n_recycle;
    s.CellIn_kmol_s=R.eq.n_cell_in; s.Erev_V=R.cell.E_rev; s.Vch_V=R.cell.V;
    s.Current_A=R.cell.I; s.CellPower_MW=R.cell.W/1e6;
    s.ElectricHeat_MW=R.Qel_tot/1e6; s.Pin_MW=R.P_in/1e6;
    s.ChargeEfficiency_pct=100*R.eta_ch;
    s.HEX_MW=R.hex.Q/1e6; s.ThermalLink_MW=R.link.Q/1e6; s.RejectedHeat_MW=R.Qrej_tot/1e6;
    s.MaxAtomResidual_kmol_s=max(abs([R.balance.res_C;R.balance.res_O;R.balance.res_N]));
    s.MaxEnergyResidual_W=max(abs(R.balance.res_E_W));
    row=struct2table(s);
end

function make_plot(t,xname,xlabelText,filename)
    f=figure('Visible','off'); tiledlayout(2,3);
    names={'yCO7','CellIn_kmol_s','Vch_V','ElectricHeat_MW','Pin_MW','ChargeEfficiency_pct'};
    labels={'Reactor outlet CO fraction','Cell inlet [kmol/s]','Cell voltage [V]', ...
        'Electrical heating [MW]','Charging input [MW]','Charging efficiency [%]'};
    for j=1:numel(names)
        nexttile; plot(t.(xname),t.(names{j}),'o-'); grid on;
        xlabel(xlabelText); ylabel(labels{j});
    end
    exportgraphics(f,filename,'Resolution',180); close(f);
end
