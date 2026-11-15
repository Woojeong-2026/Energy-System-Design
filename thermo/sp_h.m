function h = sp_h(th, i, T)
% SP_H  순수 성분 i 의 몰 엔탈피 h_i(T) [J/kmol]
%
%   h = sp_h(th, i, T)     i: 성분 번호 (1 CO, 2 CO2, 3 O2, 4 N2, 5 C(s)),  T [K]
%
% Cantera 의 절대 엔탈피(25 degC 생성엔탈피 포함)를 그대로 쓴다.
%   h_i(T) = h_f,i(298 K) + int_298^T cp_i dT
% 따라서 반응 전후 엔탈피 차이에 반응열이 자동으로 들어가고, 반응열을 따로 더하지 않는다.
% 이상기체/순수 고체로 다루므로 압력에 무관하다 (항상 p° = 1 bar 에서 평가).

    if i == 5                                         % 흑연
        if strcmp(th.api, 'legacy')
            set(th.gr, 'T', T, 'P', th.p0);
            h = enthalpy_mole(th.gr);
        else
            th.gr.TP = {T, th.p0};
            h = th.gr.H;
        end
    else                                              % 기체 (순수 성분으로 설정)
        if strcmp(th.api, 'legacy')
            set(th.gas, 'T', T, 'P', th.p0, 'X', th.Xstr{i});
            h = enthalpy_mole(th.gas);
        else
            th.gas.TPX = {T, th.p0, th.Xstr{i}};
            h = th.gas.H;
        end
    end
end
