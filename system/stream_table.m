function tab = stream_table(S, MW)
% STREAM_TABLE  스트림 구조체 배열을 결과표로 만든다  -> 과제 (a)
%
%   tab = stream_table(S, MW)
%
%   S(k).no   : 스트림 번호 (Figure 2 의 번호)
%   S(k).n    : 성분 몰유량 [CO CO2 O2 N2 C(s)] [kmol/s]
%   S(k).T    : 온도 [K],  S(k).P : 압력 [Pa],  S(k).name : 설명
%   MW        : 분자량 [kg/kmol] (1x5)
%
%   표의 열: 번호, 설명, T [degC], P [bar], 몰유량 [kmol/s], 질량유량 [kg/s], 몰분율

    ns = numel(S);
    No = zeros(ns, 1);  Desc = strings(ns, 1);
    T_C = zeros(ns, 1);  P_bar = zeros(ns, 1);
    n_kmol_s = zeros(ns, 1);  m_kg_s = zeros(ns, 1);  X = zeros(ns, 5);
    for k = 1:ns
        No(k)       = S(k).no;
        Desc(k)     = S(k).name;
        T_C(k)      = S(k).T - 273.15;
        P_bar(k)    = S(k).P/1e5;
        n_kmol_s(k) = sum(S(k).n);              % 전체 몰유량 = sum n_i
        m_kg_s(k)   = sum(S(k).n.*MW);          % 질량유량   = sum n_i M_i
        X(k, :)     = S(k).n/n_kmol_s(k);       % 몰분율 x_i = n_i / sum n
    end
    x_CO = X(:, 1); x_CO2 = X(:, 2); x_O2 = X(:, 3); x_N2 = X(:, 4); x_C = X(:, 5);
    tab = table(No, Desc, T_C, P_bar, n_kmol_s, m_kg_s, x_CO, x_CO2, x_O2, x_N2, x_C);
end
