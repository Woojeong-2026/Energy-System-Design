function V = verify_charging(R, cfg, th)
% VERIFY_CHARGING Recompute checks from raw streams, not R.balance residuals.
% Units: kmol/s, K, Pa, W. Numerical tolerances are not engineering margins.
% A passing result verifies this assignment model, not physical equipment sizing.
    rows = cell(0,4);
    need = {'S','units','cell','hex','link','P_in','nC','eta_ch','Qel_tot','Qrej_tot','eq','elec'};
    assert(all(isfield(R,need)), 'verify_charging:schema','Incomplete R interface.');
    assert(numel(R.S)==13 && isequal([R.S.no],1:13), ...
        'verify_charging:schema','Streams must be ordered 1:13.');
    n = vertcat(R.S.n); T = [R.S.T]'; P = [R.S.P]';
    assert(isequal(size(n),[13 5]) && isreal(n) && all(isfinite(n(:))) && ...
        isreal(T) && all(isfinite(T)) && isreal(P) && all(isfinite(P)), ...
        'verify_charging:schema','Invalid raw stream values.');
    assert(numel(R.units)==7,'verify_charging:schema','Expected seven duty units.');
    q = [R.units.Q]'; ql = [R.units.Qlink]'; qe = [R.units.Qel]'; qr = [R.units.Qrej]';
    assert(isreal([q;ql;qe;qr]) && all(isfinite([q;ql;qe;qr])), ...
        'verify_charging:schema','Invalid heat duty.');
    flowtol = 1e-9 + 1e-8*max(abs(n(:)));
    powtol = 0.1 + 1e-8*max([abs(R.P_in);abs(q);1]);
    T0 = cfg.T0_C+273.15; Tc = cfg.Tcell_C+273.15; Tb = cfg.TBch_C+273.15;
    add('Nonnegative flows',max(0,-min(n(:))),flowtol);
    add('System pressure',max(abs(P-cfg.P_bar*1e5)),1e-4);
    add('Ambient stream temperatures',max(abs(T([1 2 9 11 13])-T0)),1e-5);
    add('Cell stream temperatures',max(abs(T([4 5 12])-Tc)),1e-5);
    add('Reactor stream temperatures',max(abs(T([6 7 10])-Tb)),1e-5);
    expected = zeros(1,5); expected(2) = cfg.nCO2_feed;
    eqrow('Fresh feed',n(1,:),expected,flowtol);
    add('Fuel loop contains only CO CO2',max(abs(n([2:9],3:5)),[],'all'),flowtol);
    add('Carbon outlet purity',max(abs(n([10 11],1:4)),[],'all'),flowtol);
    add('Oxygen outlet purity',max(abs(n([12 13],[1 2 4 5])),[],'all'),flowtol);
    eqrow('Mixer',n(2,:),n(1,:)+n(9,:),flowtol);
    pairs = [2 3;3 4;5 6;7 8;8 9;10 11;12 13];
    for k=1:size(pairs,1)
        a=pairs(k,1); b=pairs(k,2);
        eqrow(sprintf('Frozen composition %d->%d',a,b),n(a,:),n(b,:),flowtol);
    end
    extent = cfg.X_soec*n(4,2);
    eqrow('SOEC fuel stoichiometry',n(5,:),n(4,:)+extent*[1 -1 0 0 0],flowtol);
    eqrow('SOEC oxygen stoichiometry',n(12,:),extent*[0 0 0.5 0 0],flowtol);
    xi = n(10,5);
    eqrow('Reactor stoichiometry',n(7,:)+n(10,:),n(6,:)+xi*[-2 1 0 0 1],flowtol);
    add('Net carbon inventory',abs(xi-cfg.nCO2_feed),flowtol);
    add('Net oxygen generation',abs(n(13,3)-cfg.nCO2_feed),flowtol);
    add('R.nC agrees with stream 11',abs(R.nC-n(11,5)),flowtol);

    % Explicit elemental matrices: independent of balance_table and th.A.
    A = [1 1 0 0 1;1 2 2 0 0;0 0 0 2 0];
    ins = { [1 9], [2 7],3,4,5,6,8,10,12,1 };
    outs = { 2,[3 8],4,[5 12],6,[7 10],9,11,13,[11 13] };
    labels = {'Mixer','HEX1','H1','SOEC','Trim','Reactor','Recycle cooler','C cooler','O2 cooler','SYSTEM'};
    H = zeros(13,1);
    for j=1:13, H(j)=stream_H(th,n(j,:),T(j)); end
    Q = [0;0;q(1);q(2);q(3);q(4);q(5);q(6);q(7);-sum(qr)];
    W = [0;0;0;R.cell.W;0;0;0;0;0;R.P_in];
    for j=1:numel(ins)
        dn=sum(n(ins{j},:),1)-sum(n(outs{j},:),1);
        add([labels{j} ' C/O/N'],max(abs(A*dn')),flowtol);
        res=sum(H(ins{j}))+Q(j)+W(j)-sum(H(outs{j}));
        add([labels{j} ' energy'],abs(res),powtol);
    end
    % Independent equilibrium and voltage calculations from h and s.
    gB=zeros(1,5); gC=zeros(1,5);
    for j=1:5
        gB(j)=sp_h(th,j,Tb)-Tb*sp_s0(th,j,Tb);
        gC(j)=sp_h(th,j,Tc)-Tc*sp_s0(th,j,Tc);
    end
    lnK=-(gB(5)+gB(2)-2*gB(1))/(cfg.Ru*Tb);
    y=n(7,1:2)/sum(n(7,1:2)); pr=P(7)/cfg.p0;
    add('Boudouard ln quotient',abs(log(y(2))-2*log(y(1))-log(pr)-lnK),1e-8);
    add('Reported K_B',abs(log(R.eq.K_B)-lnK),1e-8);
    pCO=n(5,1)/sum(n(5,1:4))*P(5)/cfg.p0;
    pCO2=n(5,2)/sum(n(5,1:4))*P(5)/cfg.p0;
    pO2=n(12,3)/sum(n(12,1:4))*P(12)/cfg.p0;
    E0=-(gC(2)-gC(1)-0.5*gC(3))/(2000*cfg.F);
    Er=E0+cfg.Ru*Tc/(2000*cfg.F)*log(pCO*sqrt(pO2)/pCO2);
    I=2000*cfg.F*extent;
    add('Nernst streams 5 and 12',abs(R.cell.E_rev-Er),1e-8);
    add('SOEC voltage loss sign',abs(R.cell.V-Er-cfg.dV_loss),1e-8);
    add('Faraday current kmol to mol',abs(R.cell.I-I),1e-2+1e-8*abs(I));
    add('Electrical power IV',abs(R.cell.W-I*(Er+cfg.dV_loss)),powtol);
    add('Thermoneutral voltage',abs(R.cell.V_tn-(H(5)+H(12)-H(4))/I),1e-8);
    add('Endo/exothermic cell sign',abs(q(2)-I*(R.cell.V_tn-R.cell.V)),powtol);

    % Reconstruct link without calling settle_heat.
    want=zeros(7,1); L=0;
    if q(4)<0 && q(2)>0 && Tb-Tc>=cfg.dT_link-1e-9
        L=min(-cfg.f_link*q(4),q(2)); want(4)=-L; want(2)=L;
    elseif q(2)<0 && q(4)>0 && Tc-Tb>=cfg.dT_link-1e-9
        L=min(-cfg.f_link*q(2),q(4)); want(2)=-L; want(4)=L;
    end
    eqrow('Thermal link eligibility and amount',ql,want,powtol);
    add('Reported link',abs(R.link.Q-L),powtol);
    eqrow('Residual electric heat',qe,max(q,0)-max(want,0),powtol);
    eqrow('Residual rejected heat',qr,max(-q,0)-max(-want,0),powtol);
    add('Nonnegative electric/rejected heat',max(0,-min([qe;qr])),powtol);
    add('No double counting recovered heat',abs(R.P_in-R.cell.W-sum(qe)),powtol);
    add('Qel total',abs(R.Qel_tot-sum(qe)),powtol);
    add('Qrej total',abs(R.Qrej_tot-sum(qr)),powtol);
    add('Electrical detail totals',max(abs([R.elec.W_cell-R.cell.W,R.elec.Qel_tot-sum(qe),R.elec.W_sep,R.elec.W_comp])),powtol);
    add('Charge efficiency',abs(R.eta_ch-xi*cfg.LHV_C/R.P_in),1e-9);
    add('Positive charge input',double(R.P_in<=0),0);

    % HEX independently recomputed at 201 positions, not just its terminals.
    Qh=H(7)-H(8); Qc=H(3)-H(2);
    Qmax=min(H(7)-stream_H(th,n(7,:),T(2)),stream_H(th,n(2,:),T(7))-H(2));
    add('HEX hot/cold duties',abs(Qh-Qc),powtol);
    add('HEX effectiveness',abs(Qh-cfg.eff_hex*Qmax),powtol);
    add('HEX reported duty',abs(R.hex.Q-Qh),powtol);
    add('HEX reported Qmax',abs(R.hex.Qmax-Qmax),powtol);
    qp=linspace(0,Qh,201)'; Th=zeros(size(qp)); Tl=Th;
    for j=1:numel(qp)
        Th(j)=stream_T(th,n(7,:),H(8)+qp(j),T(2),T(7));
        Tl(j)=stream_T(th,n(2,:),H(2)+qp(j),T(2),T(7));
    end
    dT=Th-Tl;
    add('HEX heat flows hot to cold (201 samples)',double(min(dT)<=0),0);
    % Entropy generation for the whole adiabatic recuperator.
    sgen=entropy_flow(n(3,:),T(3),P(3))+entropy_flow(n(8,:),T(8),P(8)) ...
        -entropy_flow(n(2,:),T(2),P(2))-entropy_flow(n(7,:),T(7),P(7));
    add('HEX nonnegative entropy generation',max(0,-sgen),1e-3);
    V.checks=cell2table(rows,'VariableNames',{'Check','Residual','Limit','Pass'});
    V.all_pass=all(V.checks.Pass);
    V.hex_profile=table(qp,Th,Tl,dT,'VariableNames',{'q_W','T_hot_K','T_cold_K','dT_K'});
    V.dTmin_K=min(dT); V.Sgen_hex_W_K=sgen;
    V.scope='Assignment steady-state model only; hardware size and kinetics unverified.';

    function add(name,res,lim)
        ok=isreal(res) && isscalar(res) && isfinite(res) && res<=lim;
        rows(end+1,:)={name,res,lim,ok};
    end
    function eqrow(name,a,b,tol)
        add(name,max(abs(a(:)-b(:))),tol);
    end
    function sf=entropy_flow(nn,tt,pp)
        sf=0; ng=sum(nn(1:4));
        for ii=1:4
            if nn(ii)>0
                % Cantera gas constant for entropy of mixing, not rounded assignment Ru.
                sf=sf+nn(ii)*(sp_s0(th,ii,tt)-8314.46261815324*log(nn(ii)/ng*pp/cfg.p0));
            end
        end
        sf=sf+nn(5)*sp_s0(th,5,tt);
    end
end
