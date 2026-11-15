function hx = hex_counterflow(th, nh, Th_in, nc, Tc_in, eff)
% HEX_COUNTERFLOW  단열 향류 열교환기 (유용도 모델, HEX 1 / 3 / 4 공통)
%
%   hx = hex_counterflow(th, nh, Th_in, nc, Tc_in, eff)
%
%   nh, Th_in : 고온측 몰유량 [kmol/s] (1x5), 입구 온도 [K]
%   nc, Tc_in : 저온측 몰유량 [kmol/s] (1x5), 입구 온도 [K]
%   eff       : 유용도 [-]
%
%   hx.Q      : 실제 전달열 [W]
%   hx.Qmax   : 최대 전달 가능 열 [W]
%   hx.Th_out, hx.Tc_out : 출구 온도 [K]
%   hx.dTmin  : 열교환기 내부의 최소 (T_hot - T_cold) [K]  (> 0 이어야 타당)
%
% 과제 지시: "Calculate maximum transferable heat from inlet enthalpies.
%            Solve outlet temperatures by energy balance using the effectiveness.
%            Check that heat flows from the hotter stream to the colder stream throughout."

    if Th_in <= Tc_in                    % 온도가 역전되어 열이 흐를 수 없는 경우
        hx = struct('Q', 0, 'Qmax', 0, 'Th_out', Th_in, 'Tc_out', Tc_in, ...
                    'dTmin', Th_in - Tc_in);
        return
    end

    Hh_in = stream_H(th, nh, Th_in);     % 고온측 입구 엔탈피 유량 [W]
    Hc_in = stream_H(th, nc, Tc_in);     % 저온측 입구 엔탈피 유량 [W]

    % 최대 전달 가능 열 (입구 엔탈피 기준)
    %   고온측이 저온측 입구 온도까지 식을 때 내놓는 열  : H_h(T_h,in) - H_h(T_c,in)
    %   저온측이 고온측 입구 온도까지 데워질 때 받는 열  : H_c(T_h,in) - H_c(T_c,in)
    %   Q_max = 둘 중 작은 값 (열용량이 작은 쪽이 한계)
    Qmax_h = Hh_in - stream_H(th, nh, Tc_in);
    Qmax_c = stream_H(th, nc, Th_in) - Hc_in;
    Qmax   = min(Qmax_h, Qmax_c);

    % 실제 전달열 Q = eff * Q_max
    Q = eff*Qmax;

    % 출구 온도: 각 유체의 에너지수지 (단열이므로 고온측이 잃은 열 = 저온측이 얻은 열)
    %   H_h(T_h,out) = H_h(T_h,in) - Q,   H_c(T_c,out) = H_c(T_c,in) + Q
    Th_out = stream_T(th, nh, Hh_in - Q, Tc_in, Th_in);
    Tc_out = stream_T(th, nc, Hc_in + Q, Tc_in, Th_in);

    % 내부 온도분포 점검: 저온측 입구 끝에서 누적 전달열 q (0 -> Q) 에 대해
    %   저온측 온도 T_c(q) : H_c = H_c,in + q
    %   고온측 온도 T_h(q) : H_h = H_h,out + q   (향류이므로 같은 위치에서 맞닿는 유체)
    % 전 구간에서 T_h - T_c > 0 이면 열이 항상 고온 -> 저온으로 흐른다.
    Np = 20;                             % 구간 수
    dT = zeros(Np + 1, 1);
    for j = 0:Np
        q  = Q*j/Np;
        Tc = stream_T(th, nc, Hc_in + q,     Tc_in, Th_in);
        Th = stream_T(th, nh, Hh_in - Q + q, Tc_in, Th_in);
        dT(j + 1) = Th - Tc;
    end

    hx = struct('Q', Q, 'Qmax', Qmax, 'Th_out', Th_out, 'Tc_out', Tc_out, 'dTmin', min(dT));
end
