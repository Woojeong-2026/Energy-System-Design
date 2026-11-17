function summary = run_charging_review(repoRoot)
% RUN_CHARGING_REVIEW Run the integrated charging model and independent checks.
%   summary = run_charging_review();  % infer repository root from this file
%   Results: results/charging_review/{base,P1,...,T900}/ and summary.csv.
%   Nine saved cases cover seven distinct operating points.
    if nargin<1, repoRoot=fileparts(fileparts(mfilename('fullpath'))); end
    assert(isfile(fullfile(repoRoot,'noon_config.m')),'Repository root is invalid.');
    oldpath=path; restore=onCleanup(@() path(oldpath)); %#ok<NASGU>
    addpath(repoRoot,fullfile(repoRoot,'thermo'),fullfile(repoRoot,'system'), ...
        fullfile(repoRoot,'modes'),fullfile(repoRoot,'charging'),'-begin');
    clear charging_mode
    modelPath=which('charging_mode'); fprintf('Charging implementation: %s\n',modelPath);
    th=thermo_init();
    outdir=fullfile(repoRoot,'results','charging_review');
    if ~isfolder(outdir), mkdir(outdir); end
    cases={'base','P1','P2','P5','P10','T750','T800','T850','T900'};
    cfgs={noon_config(),noon_config('P_bar',1),noon_config('P_bar',2), ...
        noon_config('P_bar',5),noon_config('P_bar',10),noon_config('Tcell_C',750), ...
        noon_config('Tcell_C',800),noon_config('Tcell_C',850),noon_config('Tcell_C',900)};
    summary=table();
    for k=1:numel(cases)
        cfg=cfgs{k}; R=charging_mode(cfg,th); V=verify_charging(R,cfg,th);
        folder=fullfile(outdir,cases{k}); if ~isfolder(folder), mkdir(folder); end
        writetable(R.streams,fullfile(folder,'streams.csv'));
        writetable(R.balance,fullfile(folder,'reported_balances.csv'));
        writetable(R.heat,fullfile(folder,'heat.csv'));
        writetable(V.checks,fullfile(folder,'independent_checks.csv'));
        writetable(V.hex_profile,fullfile(folder,'hex_profile.csv'));
        export_charging_json(fullfile(folder,'raw.json'),R,cfg);
        save(fullfile(folder,'result.mat'),'R','V','cfg','modelPath');
        row=table(string(cases{k}),cfg.P_bar,cfg.Tcell_C,R.eq.K_B,R.eq.yCO, ...
            R.eq.n_cell_in,R.cell.V,R.cell.W/1e6,R.Qel_tot/1e6,R.Qrej_tot/1e6, ...
            R.P_in/1e6,100*R.eta_ch,R.link.Q/1e6,V.dTmin_K,V.all_pass, ...
            'VariableNames',{'Case','P_bar','Tcell_C','K_B','yCO7','n4_kmol_s','V_ch', ...
            'Wcell_MW','Qel_MW','Qrej_MW','Pin_MW','eta_ch_pct','Qlink_MW','dTmin_K','AllChecksPass'});
        summary=[summary;row]; %#ok<AGROW>
    end
    writetable(summary,fullfile(outdir,'summary.csv')); disp(summary);
    % Required pressure and cell-temperature trends, all other values fixed.
    f=figure('Visible','off'); tiledlayout(2,2);
    nexttile; plot(summary.P_bar(2:5),summary.eta_ch_pct(2:5),'o-'); grid on;
    xlabel('P [bar]'); ylabel('Charge efficiency [%]');
    nexttile; plot(summary.P_bar(2:5),summary.n4_kmol_s(2:5),'o-'); grid on;
    xlabel('P [bar]'); ylabel('Cell fuel inlet [kmol/s]');
    nexttile; plot(summary.Tcell_C(6:9),summary.eta_ch_pct(6:9),'o-'); grid on;
    xlabel('Cell temperature [degC]'); ylabel('Charge efficiency [%]');
    nexttile; plot(summary.Tcell_C(6:9),summary.Pin_MW(6:9),'o-'); grid on;
    xlabel('Cell temperature [degC]'); ylabel('Charging input [MW]');
    exportgraphics(f,fullfile(outdir,'charging_sweeps.png'),'Resolution',180); close(f);
    assert(all(summary.AllChecksPass),'run_charging_review:failed', ...
        'Some model checks failed; inspect independent_checks.csv.');
end
