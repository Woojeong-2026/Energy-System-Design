function test_charging_review(repoRoot)
% MATLAB/Cantera regression + fault injection tests (must be run in MATLAB).
    if nargin<1, repoRoot=fileparts(fileparts(mfilename('fullpath'))); end
    oldpath=path; cleanup=onCleanup(@() path(oldpath)); %#ok<NASGU>
    addpath(repoRoot,fullfile(repoRoot,'thermo'),fullfile(repoRoot,'system'),fullfile(repoRoot,'charging'));
    addpath(fullfile(repoRoot,'modes'),'-begin'); clear charging_mode
    th=thermo_init(); cfg=noon_config(); R=charging_mode(cfg,th);
    V=verify_charging(R,cfg,th); assert(V.all_pass,'Baseline checks failed.');
    assert(abs(R.nC-1)<1e-10 && abs(R.cell.I-385940000)<1e-3);
    assert(R.link.Q==0,'Baseline thermal link must be off.');
    cfg2=noon_config('nCO2_feed',2); R2=charging_mode(cfg2,th);
    assert(abs(R2.P_in/R.P_in-2)<1e-8 && abs(R2.eta_ch-R.eta_ch)<1e-9);
    for name={'X_soec','nCO2_feed','P_bar'}
        c=cfg; c.(name{1})=0; rejected=false;
        try, charging_mode(c,th); catch e, rejected=strcmp(e.identifier,'charging_mode:input'); end
        assert(rejected,['Invalid input accepted: ' name{1}]);
    end
    c=cfg; c.X_soec=1; rejected=false;
    try, charging_mode(c,th); catch e, rejected=strcmp(e.identifier,'charging_mode:input'); end
    assert(rejected,'X=1 accepted despite singular Nernst quotient.');
    bad=R; bad.S(11).n(5)=0.9; V=verify_charging(bad,cfg,th);
    assert(~V.all_pass,'Carbon inventory mutation went undetected.');
    bad=R; bad.cell.I=bad.cell.I/1000; V=verify_charging(bad,cfg,th);
    assert(~V.all_pass,'Faraday unit mutation went undetected.');
    bad=R; bad.link.Q=1e6; V=verify_charging(bad,cfg,th);
    assert(~V.all_pass,'Impossible heat link went undetected.');
    c=noon_config('TBch_C',869.999); r=charging_mode(c,th); assert(r.link.Q==0);
    c=noon_config('TBch_C',870); r=charging_mode(c,th); assert(r.link.Q>0);
    V=verify_charging(r,c,th); assert(V.all_pass,'Thermal link boundary case failed.');
    fprintf('MATLAB charging regression and fault-injection tests passed.\n');
end
