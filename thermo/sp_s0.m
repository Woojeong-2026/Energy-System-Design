function s = sp_s0(th, i, T)
% SP_S0  순수 성분 i 의 표준 몰 엔트로피 s°_i(T, p° = 1 bar) [J/(kmol K)]
%
%   s = sp_s0(th, i, T)
%
% 순수 성분을 p° = 1 bar 로 설정하므로 혼합 엔트로피/압력 항이 없는 표준 상태 값이다.
% 분압의 영향은 평형상수 식과 Nernst 식에서 ln(p_i/p°) 로 따로 넣는다.

    if i == 5                                         % 흑연
        if strcmp(th.api, 'legacy')
            set(th.gr, 'T', T, 'P', th.p0);
            s = entropy_mole(th.gr);
        else
            th.gr.TP = {T, th.p0};
            s = th.gr.S;
        end
    else                                              % 기체
        if strcmp(th.api, 'legacy')
            set(th.gas, 'T', T, 'P', th.p0, 'X', th.Xstr{i});
            s = entropy_mole(th.gas);
        else
            th.gas.TPX = {T, th.p0, th.Xstr{i}};
            s = th.gas.S;
        end
    end
end
