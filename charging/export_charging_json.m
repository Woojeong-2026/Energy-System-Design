function export_charging_json(filename,R,cfg)
% Raw numerical interface for the independent Python/Cantera cross-check.
    p.schema_version=1; p.cfg=cfg;
    p.n=vertcat(R.S.n); p.T=[R.S.T]; p.P=[R.S.P];
    p.q=[R.units.Q]; p.ql=[R.units.Qlink]; p.qe=[R.units.Qel]; p.qr=[R.units.Qrej];
    p.metrics=struct('K_B',R.eq.K_B,'yCO',R.eq.yCO,'E0',R.cell.E0, ...
        'E_rev',R.cell.E_rev,'V',R.cell.V,'V_tn',R.cell.V_tn,'I',R.cell.I, ...
        'W_cell',R.cell.W,'Q_hex',R.hex.Q,'Qmax_hex',R.hex.Qmax, ...
        'Q_link',R.link.Q,'Qel_tot',R.Qel_tot,'Qrej_tot',R.Qrej_tot, ...
        'P_in',R.P_in,'nC',R.nC,'eta_ch',R.eta_ch);
    fid=fopen(filename,'w');
    assert(fid>=0,'export_charging_json:file','Cannot open output file.');
    cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
    fwrite(fid,jsonencode(p,'PrettyPrint',true),'char');
end
