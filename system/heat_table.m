function tab = heat_table(units)
% HEAT_TABLE  정산된 장치 구조체 배열(settle_heat 출력)을 열부하 표로 만든다  -> 과제 (b)
%
%   열: 장치, 온도 [degC], 열부하 Q [MW], thermal link [MW], 전기 가열 [MW], 방출열 [MW]

    nu = numel(units);
    Unit = strings(nu, 1);
    T_C = zeros(nu, 1); Q_MW = zeros(nu, 1); Qlink_MW = zeros(nu, 1);
    Qel_MW = zeros(nu, 1); Qrej_MW = zeros(nu, 1);
    for u = 1:nu
        Unit(u)     = units(u).name;
        T_C(u)      = units(u).T - 273.15;
        Q_MW(u)     = units(u).Q*1e-6;
        Qlink_MW(u) = units(u).Qlink*1e-6;
        Qel_MW(u)   = units(u).Qel*1e-6;
        Qrej_MW(u)  = units(u).Qrej*1e-6;
    end
    tab = table(Unit, T_C, Q_MW, Qlink_MW, Qel_MW, Qrej_MW);
end
