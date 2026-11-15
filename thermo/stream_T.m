function T = stream_T(th, n, Htarget, Tlo, Thi)
% STREAM_T  에너지수지로 스트림 온도를 구한다 (엔탈피 -> 온도 역산)
%
%   T = stream_T(th, n, Htarget, Tlo, Thi)
%
%   stream_H(n, T) = Htarget 를 만족하는 T [K] 를 구간 [Tlo, Thi] 안에서 fzero 로 찾는다.
%   cp > 0 이라 H(T) 는 단조증가이므로 근은 하나뿐이다.

    f = @(T) (stream_H(th, n, T) - Htarget)/1e6;      % 잔차를 [MW] 로 스케일
    T = fzero(f, [Tlo - 1, Thi + 1]);                 % 양 끝에 1 K 여유
end
