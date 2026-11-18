function R = run_charging_analysis(cfg, th)
% RUN_CHARGING_ANALYSIS Energy System Design 2026 - Charging, Stream 1~13.
% MATLAB R2026a + Cantera 3.0.0. 이 파일을 noon_config.m 옆에 놓고 Run.
%
% [과제 범위] (a) 스트림/수지표 (b) 전력/열량 (c) 충전 효율
%             (d) 압력 비교 (e) 셀 온도 비교 (f) 충전 개선안의 수치 근거.
% 방전 및 왕복 효율은 팀 통합 계산에서 다룬다.
% [실행] run_charging_analysis;  -> 기준/조건 비교의 CSV/PNG 저장.
% [통합] ch = run_charging_analysis(cfg, th); -> 한 조건의 결과만 반환.
%
% 사용 공통함수 (정우정 팀원 작성, 수정 없이 재사용)
%   noon_config, thermo_init, K_boudouard, E_cell_CO, stream_H,
%   hex_counterflow, settle_heat, stream_table, balance_table, heat_table.
%   sp_h, sp_s0, sp_g0, stream_T는 위 공통함수 내부에서 호출된다.
% 이 파일의 calculate_charging이 충전 수지와 장치 연결을 계산한다.
% 기존 charging_mode, 별도 student 파일, Python 참조 계산은 필요 없다.
%
% 단위/부호: n [kmol/s], T [K], P [Pa], H/Q/W [W].
% Q>0: 장치 가열, Q<0: 장치 냉각. W>0: 전력 소비.
% 표에서는 전력/열량을 MW, 온도를 degC, 압력을 bar로 변환한다.
% 전제: 정상상태, CO/CO2 완전 재순환, 탄소 완전 분리, 반응기 평형.

    %% 1. 공통함수와 물성 객체 준비
    projectRoot=fileparts(mfilename('fullpath'));
    addpath(projectRoot,fullfile(projectRoot,'thermo'),fullfile(projectRoot,'system'),'-begin');
    if nargin<2, th=thermo_init(); end
    if nargin>0 % 팀 main에서 cfg/th를 넘겼을 때: 표 저장·스윕 없이 계산만 한다.
        R=calculate_charging(cfg,th);
        return
    end
    fprintf('Charging calculation: %s.m (single file)\n',mfilename('fullpath'));
    outdir=fullfile(projectRoot,'results','charging_analysis');
    if ~isfolder(outdir), mkdir(outdir); end

    %% 2. 과제 (a)~(e): 기준/압력/셀 온도별 계산 및 표 저장
    % 기준 조건과 중복되는 P=5, Tcell=850은 base 결과를 재사용한다.
    % 각 행의 값만 변경: 나머지는 noon_config의 PDF 지정 조건을 유지한다.
    conditions=[5 850;1 850;2 850;10 850;5 750;5 800;5 900];
    tags=["base";"P1";"P2";"P10";"T750";"T800";"T900"];
    chargingSummary=table();
    for k=1:size(conditions,1)
        cfg=noon_config('P_bar',conditions(k,1),'Tcell_C',conditions(k,2));
        R=calculate_charging(cfg,th);
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

    %% 3. 과제 (f): 지정 조건 결과와 구분한 개선 검토
    % HEX 0.95: 재순환 기체의 열을 더 회수해 전기 예열을 줄이는 경우.
    % 반응기 870 C: 850 C 셀보다 20 K 높여 열전달 조건을 만족시키는 경우.
    % 반응기 온도는 평형/순환량도 바꾸므로 열회수량만 더해서 효율을 계산하지 않는다.
    % 변경 가정의 수치 비교이며 실제 반응속도·장치 크기·재료의 타당성은 별도다.
    improvements=result_row(baseResults,"base");
    improvements=[improvements;result_row(calculate_charging(noon_config('eff_hex',0.95),th),"HEX_0.95")];
    improvements=[improvements;result_row(calculate_charging(noon_config('TBch_C',870),th),"Reactor_870C")];
    writetable(improvements,fullfile(outdir,'improvement_comparison.csv'));
    fprintf('\n(f) 변경된 가정에 따른 충전 개선 비교\n'); disp(improvements);

    %% 4. 과제 (d),(e) 비교 그림
    % 압력과 셀 온도 각각의 경향. 실제 보고 수치는 CSV 표를 사용한다.
    make_plot(pressureResults,'P_bar','Pressure [bar]',fullfile(outdir,'pressure_sweep.png'));
    make_plot(temperatureResults,'Tcell_C','Cell temperature [degC]',fullfile(outdir,'temperature_sweep.png'));
    fprintf('\n결과 폴더: %s\n',outdir);
    fprintf('전체 방전/왕복 효율은 팀 통합 계산에서 별도로 구합니다.\n');

    R=baseResults; % 입력 없이 호출하고 출력도 요청하면 기준 조건 결과를 반환한다.
end

function R = calculate_charging(cfg, th)
% 한 운전 조건의 실제 충전 계산. 아래 순서대로 읽으면 설계 흐름을 설명할 수 있다.
% 1) 재순환 물질수지/평형 -> 2) HEX/스트림 -> 3) SOEC 전기
% -> 4) 장치 열수지/열회수 -> 5) 검산 -> 6) 결과 반환.
% 반환 구조체 R은 팀 run_roundtrip.m에서 사용하는 인터페이스와 같다.

    fields = {'P_bar','T0_C','Tcell_C','TBch_C','nCO2_feed','X_soec', ...
        'eff_hex','f_link','dT_link','dV_loss','Ru','F','p0','LHV_C'};
    for k = 1:numel(fields)
        validateattributes(cfg.(fields{k}), {'numeric'}, {'real','finite','scalar'});
    end
    assert(cfg.P_bar>0 && cfg.nCO2_feed>0 && cfg.X_soec>0 && cfg.X_soec<1, ...
        'run_charging_analysis:input', 'Require P>0, feed>0, 0<X_soec<1.');
    assert(all([cfg.eff_hex cfg.f_link]>=0) && all([cfg.eff_hex cfg.f_link]<=1) ...
        && cfg.dT_link>=0 && cfg.dV_loss>=0 && cfg.Ru>0 && cfg.F>0 && cfg.LHV_C>0 ...
        && cfg.p0==1e5 && th.p0==cfg.p0, 'run_charging_analysis:input', 'Invalid recovery/constants.');
    P=cfg.P_bar*1e5; T0=cfg.T0_C+273.15; % 운전 압력 [Pa], 저장/주위 온도 [K]
    Tc=cfg.Tcell_C+273.15; Tb=cfg.TBch_C+273.15; % 셀 및 충전 반응기 온도 [K]
    assert(Tb>T0 && Tc>T0 && min([T0 Tc Tb])>=200 && max([T0 Tc Tb])<=3500, ...
        'run_charging_analysis:input', 'Invalid temperature range.');

    % 1. 물질수지: 새 CO2 공급 f = 탄소 저장률, SOEC 전환량 = 2f.
    % SOEC: CO2 -> CO+0.5 O2, 반응기: 2CO -> C+CO2.
    % C 1 mol에 CO 2 mol이 필요하며, 생성 CO2는 재순환한다.
    % X=0.70은 1회 통과 전환율이며 전체 탄소 저장률/에너지 효율이 아니다.
    % 반응기 출구 평형: K=(1-y)/(y^2*(P/p0)), graphite 활동도=1.
    % 따라서 K*(P/p0)*y^2+y-1=0의 양의 근으로 CO 몰분율 y를 구한다.
    f=cfg.nCO2_feed; X=cfg.X_soec; % 신선 CO2 공급 [kmol/s], SOEC 전환율 [-]
    [K,dG]=K_boudouard(th,Tb,cfg.Ru);
    KP=K*P/cfg.p0;
    y=2/(1+sqrt(1+4*KP));             % 0<y<1인 근: 반응기 출구 CO 몰분율
    n4_CO2=2*f/X;                  % SOEC 입구 CO2: X*n4_CO2=2f
    n7_CO2=n4_CO2-f;               % 재순환 CO2: 입구량 - 신선 공급량
    n7_total=n7_CO2/(1-y);         % 평형 조성이 정해지면 재순환 총유량도 정해진다.
    recycle=[y*n7_total n7_CO2 0 0 0];
    feed=[0 f 0 0 0]; fuel=feed+recycle;
    outlet=fuel+2*f*[1 -1 0 0 0];     % SOEC: CO2 -> CO + 0.5 O2
    carbon=[0 0 0 0 f]; oxygen=[0 0 f 0 0];

    % 2. 공통 HEX: 7->8의 열로 2->3을 예열. Q=eff*Qmax.
    % 혼합기 양쪽 입구는 모두 T0: 이상기체·단열 혼합 출구도 T0.
    hx=hex_counterflow(th,recycle,Tb,fuel,T0,cfg.eff_hex);
    % 행 번호=PDF 스트림 번호, 열=[CO CO2 O2 N2 C(s)].
    % 셀/반응기 밖에서는 기체 조성이 변하지 않는다.
    n=[feed;fuel;fuel;fuel;outlet;outlet;recycle;recycle;recycle; ...
       carbon;carbon;oxygen;oxygen];
    T=[T0 T0 hx.Tc_out Tc Tc Tb Tb hx.Th_out T0 Tb T0 Tc T0];
    names={'CO2 feed','Mixer outlet','HEX cold outlet','SOEC inlet', ...
        'SOEC fuel outlet','Reactor inlet','Reactor gas outlet','HEX hot outlet', ...
        'Recycle','Deposited carbon','Stored carbon','SOEC oxygen','O2 before vent valve'};
    % 13번 산소는 vent valve 전이므로 P. 밸브 후에는 대기압이며 일 회수는 없다.
    S=repmat(struct('no',0,'n',zeros(1,5),'T',0,'P',P,'name',""),13,1);
    H=zeros(13,1);
    for k=1:13
        S(k)=struct('no',k,'n',n(k,:),'T',T(k),'P',P,'name',string(names{k}));
        H(k)=stream_H(th,n(k,:),T(k)); % 생성엔탈피 포함: 반응열을 따로 더하지 않음
    end

    % 3. SOEC: PDF가 명시한 5·12번 조성, +0.10 V, CO당 전자 2개.
    % PDF의 'inlet' 표현과 번호가 충돌하므로 명시된 번호를 적용했다.
    y5=outlet/sum(outlet);
    [Er,E0]=E_cell_CO(th,Tc,y5(1)*P/cfg.p0,y5(2)*P/cfg.p0,P/cfg.p0,cfg.Ru,cfg.F);
    V=Er+cfg.dV_loss;               % 충전은 가역전압에 손실을 더해야 한다 [V].
    I=2*cfg.F*(2*f)*1000;           % kmol/s -> mol/s의 1000을 포함한 전류 [A].
    W=I*V;                         % 셀 소비전력 [W]. I는 전체 활성 셀의 합산 등가 전류.
    Vtn=(H(5)+H(12)-H(4))/I;       % 열중립 전압: V<Vtn이면 셀에 추가 열이 필요.

    % 4. 각 장치 Q = Hout-Hin-W. SOEC 이외에는 직접 전력 항이 없다.
    % 순서: H1 / SOEC / 5->6 조절기 / 탄소 반응기 / 재순환·탄소·산소 냉각기.
    q=[H(4)-H(3); H(5)+H(12)-H(4)-W; H(6)-H(5); ...
       H(7)+H(10)-H(6); H(9)-H(8); H(11)-H(10); H(13)-H(12)];
    unitNames={'Heater H1','SOEC','Trim 5-6','Carbon reactor', ...
        'Recycle cooler','Carbon cooler','Oxygen cooler'};
    unitT=[NaN Tc NaN Tb NaN NaN NaN]; % 단일 등온 장치가 아닌 곳의 온도는 NaN.
    units=repmat(struct('name',"",'Q',0,'T',NaN),7,1);
    for k=1:7
        units(k)=struct('name',string(unitNames{k}),'Q',q(k),'T',unitT(k));
    end
    % 음의 Q를 가진 장치가 양의 Q 장치보다 20 K 이상 뜨거울 때만 전달.
    % 회수열=min(0.85*가용열,필요열). 기준 650->850 C에서는 직접 회수가 불가.
    [units,link]=settle_heat(units,2,4,cfg);
    % 남은 가열은 100% 전기 히터. HEX 효과는 H3에 포함되어 있으므로 또 빼지 않는다.
    Qel=sum([units.Qel]); Qrej=sum([units.Qrej]); Pin=W+Qel;

    % 5. 공통 수지표: 원소 in-out=0, 에너지 Hin+Q+W-Hout=0.
    % 전체 경계는 신선 CO2(1) -> 저장 탄소(11)+방출 산소(13).
    % thermal link와 HEX는 내부 열 이동이므로 전체 경계에서 상쇄된다.
    U={'Mixer',[1 9],2,0,0; 'HEX 1',[2 7],[3 8],0,0; ...
       'Heater H1',3,4,q(1),0; 'SOEC',4,[5 12],q(2),W; ...
       'Trim',5,6,q(3),0; 'Reactor',6,[7 10],q(4),0; ...
       'Recycle cooler',8,9,q(5),0; 'Carbon cooler',10,11,q(6),0; ...
       'Oxygen cooler',12,13,q(7),0; 'CHARGING SYSTEM',1,[11 13],-Qrej,Pin};
    balance=balance_table(th,S,U);
    resEq=log((1-y)/(y^2*(P/cfg.p0)))-log(K);
    atomResidual=max(abs([balance.res_C;balance.res_O;balance.res_N]));
    energyResidual=max(abs(balance.res_E_W));
    assert(all(isfinite([n(:);H;q;V;I;Pin])) && all(n(:)>=0) && Pin>0, ...
        'run_charging_analysis:physical', 'Nonfinite or nonphysical charging result.');
    % 공통 HEX의 내부 검사점에서 고온측>저온측인지도 확인한다.
    assert(atomResidual<=1e-8*max(1,max(n(:))) && energyResidual<=0.1+1e-8*Pin ...
        && abs(resEq)<1e-8 && hx.dTmin>0, ...
        'run_charging_analysis:balance', 'Check element/energy/equilibrium/HEX residuals.');

    % 6. 팀 통합 인터페이스 + 제출용 표. 계산의 원본은 이 함수 하나다.
    R=struct('mode','charging','cfg',cfg,'S',S,'units',units,'link',link);
    R.streams=stream_table(S,th.MW); R.balance=balance; R.heat=heat_table(units);
    R.hex=struct('name',"HEX 1",'Q',hx.Q,'Qmax',hx.Qmax,'dTmin',hx.dTmin, ...
        'Th_in',Tb,'Th_out',hx.Th_out,'Tc_in',T0,'Tc_out',hx.Tc_out);
    R.eq=struct('T_C',cfg.TBch_C,'dG',dG/1000,'K_B',K,'yCO',y,'yCO2',1-y, ...
        'res_lnK',-resEq,'n_recycle',sum(recycle),'n_cell_in',sum(fuel),'X_pass',2*f/outlet(1));
    R.cell=struct('E0',E0,'E_rev',Er,'V',V,'V_tn',Vtn,'I',I,'W',W);
    R.elec=struct('W_cell',W,'Qel_tot',Qel,'W_sep',0,'W_comp',0);
    % eta_ch=(저장 탄소의 LHV 에너지)/(소비 전기). 같은 시간으로 나누면 전력 비.
    % 방출 산소와 버린 열에는 효율 credit을 주지 않는다.
    R.Qel_tot=Qel; R.Qrej_tot=Qrej; R.P_in=Pin; R.nC=f; R.eta_ch=f*cfg.LHV_C/Pin;
end

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
