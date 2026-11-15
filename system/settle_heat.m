function [units, link] = settle_heat(units, iCell, iRx, cfg)
% SETTLE_HEAT  thermal link 를 적용하고 장치별 전기 가열량/방출열을 정산한다
%
%   [units, link] = settle_heat(units, iCell, iRx, cfg)
%
%   units : 장치 구조체 배열. 입력 필드
%             name : 장치 이름
%             Q    : 열부하 [W]  (Q > 0 가열 필요, Q < 0 냉각 필요)
%             T    : 장치 온도 [K] (thermal link 판정용, 해당 없으면 NaN)
%   iCell, iRx : 셀과 Boudouard 반응기의 units 내 번호
%
%   출력에서 units 에 추가되는 필드
%     Qlink : thermal link 로 받은 열(+) / 내준 열(-) [W]
%     Qel   : 전기 저항 가열량 [W] (>= 0)
%     Qrej  : 방출열 [W] (>= 0)
%   link.Q [W], link.from, link.to : thermal link 전달열과 방향
%
% 과제 규칙
%   (1) 음(-)의 열부하 장치(donor)에서 양(+)의 열부하 장치(receiver)로만 전달
%   (2) donor 온도가 receiver 보다 dT_link(20 K) 이상 높을 때만 가능
%   (3) Q_link = min( f_link * |Q_donor| , Q_receiver )          (f_link = 0.85)
%   (4) 남은 양(+)의 열부하는 전부 전기 저항 가열(효율 100 %)
%   (5) 남은 음(-)의 열부하는 전부 방출 (열 수출에 대한 credit 없음)

    Qc = units(iCell).Q;  Tc = units(iCell).T;        % 셀 열부하 / 온도
    Qr = units(iRx).Q;    Tr = units(iRx).T;          % 반응기 열부하 / 온도
    tol = 1e-9;                                       % 온도차 비교용 부동소수 여유 [K]

    link = struct('Q', 0, 'from', 'none', 'to', 'none');
    dQ   = zeros(numel(units), 1);                    % 장치별 link 수수량 (+ 받음, - 내줌)

    if Qc < 0 && Qr > 0 && (Tc - Tr) >= cfg.dT_link - tol
        % 셀(발열) -> 반응기(흡열): 방전 모드에서 일어나는 경우
        link.Q = min(cfg.f_link*(-Qc), Qr);
        link.from = units(iCell).name;  link.to = units(iRx).name;
        dQ(iCell) = -link.Q;  dQ(iRx) = +link.Q;
    elseif Qr < 0 && Qc > 0 && (Tr - Tc) >= cfg.dT_link - tol
        % 반응기(발열) -> 셀(흡열): 충전 모드에서 반응기가 셀보다 20 K 이상 뜨거울 때만
        link.Q = min(cfg.f_link*(-Qr), Qc);
        link.from = units(iRx).name;  link.to = units(iCell).name;
        dQ(iRx) = -link.Q;  dQ(iCell) = +link.Q;
    end

    for u = 1:numel(units)
        units(u).Qlink = dQ(u);
        % 전기 가열 = (필요한 열) - (link 로 받은 열)
        units(u).Qel  = max(units(u).Q, 0)  - max(dQ(u), 0);
        % 방출열   = (내보내야 할 열) - (link 로 내준 열)
        units(u).Qrej = max(-units(u).Q, 0) - max(-dQ(u), 0);
    end
end
