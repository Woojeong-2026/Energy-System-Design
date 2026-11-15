function [E_rev, E0, dG] = E_cell_CO(th, T, pCO, pCO2, pO2, Ru, F)
% E_CELL_CO  CO + 1/2 O2 -> CO2 의 가역 셀 전압 (Nernst 식)
%
%   [E_rev, E0, dG] = E_cell_CO(th, T, pCO, pCO2, pO2, Ru, F)
%
%   pCO, pCO2, pO2 : 무차원 분압 p_i/p° (p° = 1 bar)
%   Ru [J/(kmol K)], F [C/mol]
%
%   dG°_ox(T) = g°_CO2 - g°_CO - 1/2 g°_O2                                    [J/kmol]
%   E°(T)     = -dG°_ox / (2F)
%   E_rev     = E° + (R T / 2F) * ln[ (p_CO/p°)(p_O2/p°)^(1/2) / (p_CO2/p°) ]
%
% CO 1 mol 당 전자 2 mol 이 이동한다. J/kmol 을 J/mol 로 바꾸기 위해 1000 으로 나눈다.

    dG    = sp_g0(th, th.iCO2, T) - sp_g0(th, th.iCO, T) - 0.5*sp_g0(th, th.iO2, T);
    E0    = -(dG/1000)/(2*F);                                     % 표준 가역전압 [V]
    E_rev = E0 + (Ru/1000)*T/(2*F)*log(pCO*sqrt(pO2)/pCO2);       % 가역전압 [V]
end
