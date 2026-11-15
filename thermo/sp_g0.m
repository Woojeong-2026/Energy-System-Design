function g = sp_g0(th, i, T)
% SP_G0  순수 성분 i 의 표준 몰 Gibbs 에너지 g°_i(T) [J/kmol]
%
%   g°_i(T) = h_i(T) - T * s°_i(T, 1 bar)

    g = sp_h(th, i, T) - T*sp_s0(th, i, T);
end
