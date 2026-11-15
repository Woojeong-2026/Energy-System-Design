function H = stream_H(th, n, T)
% STREAM_H  스트림의 엔탈피 유량 H [W]
%
%   H = stream_H(th, n, T)     n: 성분 몰유량 [CO CO2 O2 N2 C(s)] [kmol/s],  T [K]
%
%   H = sum_i n_i * h_i(T)      [kmol/s * J/kmol = W]
% 이상기체 혼합물이므로 혼합 엔탈피는 0 이고 압력에 무관하다.

    H = 0;
    for i = 1:5
        if n(i) ~= 0
            H = H + n(i)*sp_h(th, i, T);
        end
    end
end
