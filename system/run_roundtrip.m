function C = run_roundtrip(cfg, th)
% RUN_ROUNDTRIP  충전 + 방전 모드를 같은 조건으로 풀고 왕복 효율을 계산한다 (통합 함수)
%
%   C = run_roundtrip(cfg, th)
%
%   C.ch, C.dis : charging_mode / discharging_mode 결과 구조체
%   C.t_ratio   : t_dis / t_ch  (같은 탄소 재고 N_C 를 쓰도록 맞춘 운전 시간 비)
%   C.eta_ch, C.eta_dis, C.eta_RT : 충전 / 방전 / 전기 왕복 효율 [-]
%   C.eta_RT_energy : E_out/E_in 으로 직접 구한 왕복 효율 (C.eta_RT 와 같아야 한다)
%   C.loss      : 손실 내역표 (탄소 1 kmol 기준 [MJ/kmol-C])
%   C.closure   : 사이클 전체 에너지 닫힘 잔차 [J/kmol-C]
%
% 과제 (c)
%   eta_ch  = N_C * LHV_C / E_in
%   eta_dis = E_out / (N_C * LHV_C)
%   eta_RT  = E_out / E_in = eta_ch * eta_dis        (N_C,stored = N_C,consumed 일 때)

    ch  = charging_mode(cfg, th);        % 팀원 1 모듈: Stream 1~13
    dis = discharging_mode(cfg, th);     % 팀원 2 모듈: Stream 14~29

    % --- 운전 시간 맞추기: N_C,stored = N_C,consumed ---
    %   N_C = nC_ch * t_ch = nC_dis * t_dis   ->   t_dis/t_ch = nC_ch/nC_dis
    %   과제 기준 유량(둘 다 1 kmol/s)에서는 t_dis = t_ch
    t_ratio = ch.nC/dis.nC;

    % --- 탄소 1 kmol 당 전기 에너지 [J/kmol-C] ---
    e_in  = ch.P_in/ch.nC;               % 충전 입력 전기 E_in/N_C  = P_in / nC_ch
    e_out = dis.P_out/dis.nC;            % 방전 순출력 전기 E_out/N_C = P_out / nC_dis

    % --- 효율 ---
    eta_ch  = cfg.LHV_C/e_in;            % eta_ch  = N_C LHV_C / E_in
    eta_dis = e_out/cfg.LHV_C;           % eta_dis = E_out / (N_C LHV_C)
    eta_RT  = eta_ch*eta_dis;            % eta_RT  = eta_ch * eta_dis
    % 같은 값을 에너지로 직접 계산: E_out/E_in = (P_out t_dis)/(P_in t_ch)
    eta_RT_energy = (dis.P_out*t_ratio)/ch.P_in;

    % --- 손실 내역 (탄소 1 kmol 기준) ---
    % 사이클 전체 에너지수지: 들어온 전기 = 나간 전기 + 방출열 + 배기(Stream 29) 현열
    %   e_in - e_out = sum(q_rej,ch) + sum(q_rej,dis) + q_exhaust
    Item = strings(0, 1);  Mode = strings(0, 1);  q = zeros(0, 1);
    for u = 1:numel(ch.units)
        if ch.units(u).Qrej > 0
            Item(end + 1, 1) = ch.units(u).name;  Mode(end + 1, 1) = "charging"; %#ok<AGROW>
            q(end + 1, 1) = ch.units(u).Qrej/ch.nC;                               %#ok<AGROW>
        end
    end
    for u = 1:numel(dis.units)
        if dis.units(u).Qrej > 0
            Item(end + 1, 1) = dis.units(u).name;  Mode(end + 1, 1) = "discharging"; %#ok<AGROW>
            q(end + 1, 1) = dis.units(u).Qrej/dis.nC;                                 %#ok<AGROW>
        end
    end
    Item(end + 1, 1) = "Exhaust (stream 29) sensible heat";  Mode(end + 1, 1) = "discharging";
    q(end + 1, 1) = dis.Q_exhaust/dis.nC;

    Loss_MJ_per_kmolC = q*1e-6;                        % 손실 [MJ/kmol-C]
    Share_of_Ein_pct  = 100*q/e_in;                    % 충전 입력 전기에 대한 비율 [%]
    loss = table(Item, Mode, Loss_MJ_per_kmolC, Share_of_Ein_pct);
    loss = sortrows(loss, 'Loss_MJ_per_kmolC', 'descend');

    C.cfg = cfg;  C.ch = ch;  C.dis = dis;
    C.t_ratio = t_ratio;
    C.e_in = e_in;  C.e_out = e_out;
    C.eta_ch = eta_ch;  C.eta_dis = eta_dis;  C.eta_RT = eta_RT;
    C.eta_RT_energy = eta_RT_energy;
    C.loss = loss;
    C.closure = e_in - e_out - sum(q);                 % 0 이어야 한다 [J/kmol-C]
end
