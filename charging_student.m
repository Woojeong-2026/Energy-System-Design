function R = charging_student(cfg, th)
% CHARGING_STUDENT 과제 PDF의 충전 계산: Stream 1~13, 정상상태·완전 재순환.
% 기존 charging_mode를 호출하지 않는다. thermo/와 system/ 공통함수만 사용.
% 단위: kmol/s, K, Pa, W. Q>0은 가열, W>0은 전력 소비.
% 반환 구조체 R은 팀 main.m / run_roundtrip.m의 충전 인터페이스와 같다.

    fields = {'P_bar','T0_C','Tcell_C','TBch_C','nCO2_feed','X_soec', ...
        'eff_hex','f_link','dT_link','dV_loss','Ru','F','p0','LHV_C'};
    for k = 1:numel(fields)
        validateattributes(cfg.(fields{k}), {'numeric'}, {'real','finite','scalar'});
    end
    assert(cfg.P_bar>0 && cfg.nCO2_feed>0 && cfg.X_soec>0 && cfg.X_soec<1, ...
        'charging_student:input', 'Require P>0, feed>0, 0<X_soec<1.');
    assert(all([cfg.eff_hex cfg.f_link]>=0) && all([cfg.eff_hex cfg.f_link]<=1) ...
        && cfg.dT_link>=0 && cfg.dV_loss>=0 && cfg.Ru>0 && cfg.F>0 && cfg.LHV_C>0 ...
        && cfg.p0==1e5 && th.p0==cfg.p0, 'charging_student:input', 'Invalid recovery/constants.');
    P=cfg.P_bar*1e5; T0=cfg.T0_C+273.15;
    Tc=cfg.Tcell_C+273.15; Tb=cfg.TBch_C+273.15;
    assert(Tb>T0 && Tc>T0 && min([T0 Tc Tb])>=200 && max([T0 Tc Tb])<=3500, ...
        'charging_student:input', 'Invalid temperature range.');

    % 1. 물질수지: 새 CO2 공급 f = 탄소 저장률, SOEC 전환량 = 2f.
    % 반응기 출구 평형: K*(P/p0)*yCO^2 + yCO - 1 = 0.
    f=cfg.nCO2_feed; X=cfg.X_soec;
    [K,dG]=K_boudouard(th,Tb,cfg.Ru);
    KP=K*P/cfg.p0;
    y=2/(1+sqrt(1+4*KP));             % 0<y<1인 근: 반응기 출구 CO 몰분율
    n4_CO2=2*f/X;
    n7_CO2=n4_CO2-f;
    n7_total=n7_CO2/(1-y);
    recycle=[y*n7_total n7_CO2 0 0 0];
    feed=[0 f 0 0 0]; fuel=feed+recycle;
    outlet=fuel+2*f*[1 -1 0 0 0];     % SOEC: CO2 -> CO + 0.5 O2
    carbon=[0 0 0 0 f]; oxygen=[0 0 f 0 0];

    % 2. 공통 HEX 함수가 열량과 두 출구 온도를 계산한다.
    % 혼합기 양쪽 입구는 모두 T0: 이상기체·단열 혼합 출구도 T0.
    hx=hex_counterflow(th,recycle,Tb,fuel,T0,cfg.eff_hex);
    n=[feed;fuel;fuel;fuel;outlet;outlet;recycle;recycle;recycle; ...
       carbon;carbon;oxygen;oxygen];
    T=[T0 T0 hx.Tc_out Tc Tc Tb Tb hx.Th_out T0 Tb T0 Tc T0];
    names={'CO2 feed','Mixer outlet','HEX cold outlet','SOEC inlet', ...
        'SOEC fuel outlet','Reactor inlet','Reactor gas outlet','HEX hot outlet', ...
        'Recycle','Deposited carbon','Stored carbon','SOEC oxygen','O2 before vent valve'};
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
    V=Er+cfg.dV_loss; I=2*cfg.F*(2*f)*1000; W=I*V;
    Vtn=(H(5)+H(12)-H(4))/I;

    % 4. 각 장치 Q = Hout-Hin-W. SOEC 이외에는 직접 전력 항이 없다.
    q=[H(4)-H(3); H(5)+H(12)-H(4)-W; H(6)-H(5); ...
       H(7)+H(10)-H(6); H(9)-H(8); H(11)-H(10); H(13)-H(12)];
    unitNames={'Heater H1','SOEC','Trim 5-6','Carbon reactor', ...
        'Recycle cooler','Carbon cooler','Oxygen cooler'};
    unitT=[NaN Tc NaN Tb NaN NaN NaN];
    units=repmat(struct('name',"",'Q',0,'T',NaN),7,1);
    for k=1:7
        units(k)=struct('name',string(unitNames{k}),'Q',q(k),'T',unitT(k));
    end
    [units,link]=settle_heat(units,2,4,cfg); % 셀↔반응기 열회수의 온도·방향 조건
    Qel=sum([units.Qel]); Qrej=sum([units.Qrej]); Pin=W+Qel;

    % 5. 공통 수지표: 열교환기·혼합기 및 전체 충전 경계까지 검산.
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
        'charging_student:physical', 'Nonfinite or nonphysical charging result.');
    assert(atomResidual<=1e-8*max(1,max(n(:))) && energyResidual<=0.1+1e-8*Pin ...
        && abs(resEq)<1e-8 && hx.dTmin>0, ...
        'charging_student:balance', 'Check element/energy/equilibrium/HEX residuals.');

    % 6. 팀 통합 인터페이스 + 제출용 표. 계산의 원본은 이 함수 하나다.
    R=struct('mode','charging','cfg',cfg,'S',S,'units',units,'link',link);
    R.streams=stream_table(S,th.MW); R.balance=balance; R.heat=heat_table(units);
    R.hex=struct('name',"HEX 1",'Q',hx.Q,'Qmax',hx.Qmax,'dTmin',hx.dTmin, ...
        'Th_in',Tb,'Th_out',hx.Th_out,'Tc_in',T0,'Tc_out',hx.Tc_out);
    R.eq=struct('T_C',cfg.TBch_C,'dG',dG/1000,'K_B',K,'yCO',y,'yCO2',1-y, ...
        'res_lnK',-resEq,'n_recycle',sum(recycle),'n_cell_in',sum(fuel),'X_pass',2*f/outlet(1));
    R.cell=struct('E0',E0,'E_rev',Er,'V',V,'V_tn',Vtn,'I',I,'W',W);
    R.elec=struct('W_cell',W,'Qel_tot',Qel,'W_sep',0,'W_comp',0);
    R.Qel_tot=Qel; R.Qrej_tot=Qrej; R.P_in=Pin; R.nC=f; R.eta_ch=f*cfg.LHV_C/Pin;
end
