function tab = balance_table(th, S, U)
% BALANCE_TABLE  장치별 원소(C, O, N) 수지와 에너지수지 잔차  -> 과제 (a)
%
%   tab = balance_table(th, S, U)
%
%   S : 스트림 구조체 배열 (S(k).no 가 Figure 2 의 스트림 번호)
%   U : 장치 정의 셀 배열, 각 행 = {이름, 입구 스트림 번호들, 출구 스트림 번호들, Q [W], W [W]}
%       Q : 장치로 들어가는 열 (+), W : 장치로 들어가는 전력 (+)
%
%   원소 잔차   = sum(원자 유량 in) - sum(원자 유량 out)          [kmol/s]
%   에너지 잔차 = sum(H_in) + Q + W - sum(H_out)                  [W]
%   모든 잔차가 0 (부동소수 오차 수준) 이면 수지가 닫힌 것이다.

    no = [S.no];                                 % 스트림 번호 목록
    nu = size(U, 1);
    Unit = strings(nu, 1);
    Q_W = zeros(nu, 1); W_W = zeros(nu, 1);
    res_C = zeros(nu, 1); res_O = zeros(nu, 1); res_N = zeros(nu, 1); res_E_W = zeros(nu, 1);
    for u = 1:nu
        Unit(u) = U{u, 1};
        atoms = zeros(3, 1);  Hsum = 0;
        for k = U{u, 2}                          % 입구 스트림: 더한다
            s = S(no == k);
            atoms = atoms + th.A*s.n';
            Hsum  = Hsum + stream_H(th, s.n, s.T);
        end
        for k = U{u, 3}                          % 출구 스트림: 뺀다
            s = S(no == k);
            atoms = atoms - th.A*s.n';
            Hsum  = Hsum - stream_H(th, s.n, s.T);
        end
        Q_W(u) = U{u, 4};  W_W(u) = U{u, 5};
        res_C(u) = atoms(1); res_O(u) = atoms(2); res_N(u) = atoms(3);
        res_E_W(u) = Hsum + Q_W(u) + W_W(u);
    end
    tab = table(Unit, Q_W, W_W, res_C, res_O, res_N, res_E_W);
end
