function [K, dG] = K_boudouard(th, T, Ru)
% K_BOUDOUARD  석출 방향 Boudouard 반응 2CO <-> C(s) + CO2 의 평형상수
%
%   [K, dG] = K_boudouard(th, T, Ru)     T [K],  Ru [J/(kmol K)]
%
%   dG°_B(T) = g°_C(s) + g°_CO2 - 2 g°_CO                       [J/kmol]
%   K_B(T)   = exp(-dG°_B/(R T)) = (p_CO2/p°) / (p_CO/p°)^2      (p° = 1 bar, a_C = 1)
%
% 방전(가스화, C + CO2 -> 2CO)에서도 같은 K_B 를 쓴다. 역반응의 평형상수는 1/K_B 이고
% 평형 조성 조건식 자체는 동일하다.

    dG = sp_g0(th, th.iC, T) + sp_g0(th, th.iCO2, T) - 2*sp_g0(th, th.iCO, T);
    K  = exp(-dG/(Ru*T));
end
