function th = thermo_init()
% THERMO_INIT  Cantera 물성 객체를 한 번만 만들어 구조체로 돌려준다 (팀 공통 물성 진입점)
%
%   th = thermo_init();
%
% 물성 데이터 출처
%   기체 CO, CO2, O2, N2 : Cantera 'gri30.yaml'  (GRI-Mech 3.0, NASA 7계수 다항식, 이상기체)
%   고체 탄소 C(gr)      : Cantera 'graphite.yaml' (과제 지정 파일), 활동도 = 1
%
% 성분 번호 (팀 공통, 모든 조성 벡터 n 의 순서)
%   1 = CO,  2 = CO2,  3 = O2,  4 = N2,  5 = C(s, graphite)
%
% 출력 th 의 필드
%   th.gas, th.gr : Cantera 객체 (기체, 흑연)
%   th.api        : 'legacy' (Cantera 3.0.0 기본 MATLAB toolbox) 또는 'new' (실험적 toolbox)
%   th.MW         : 분자량 [kg/kmol], 1x5
%   th.A          : 원소 행렬 3x5 (행: C, O, N 원자 수 / 열: 성분)
%   th.names      : 성분 이름
%
% Cantera 를 직접 부르는 곳은 thermo_init, sp_h, sp_s0 세 파일뿐이다.
% toolbox 문법이 다르면 이 세 파일만 고치면 된다.

    th.p0    = 1e5;                                   % 표준압력 1 bar [Pa]
    th.names = {'CO', 'CO2', 'O2', 'N2', 'C(s)'};     % 성분 이름
    th.Xstr  = {'CO:1', 'CO2:1', 'O2:1', 'N2:1'};     % 순수 기체 조성 문자열
    th.iCO = 1; th.iCO2 = 2; th.iO2 = 3; th.iN2 = 4; th.iC = 5;

    th.gas = Solution('gri30.yaml');                  % 기체 물성
    try
        th.gr = Solution('graphite.yaml');            % 흑연 물성
    catch
        th.gr = ThermoPhase('graphite.yaml');         % 반응이 없는 상이라 Solution 이 실패할 때
    end

    % toolbox 문법 판별: legacy 는 set(gas,'T',...) / enthalpy_mole(gas),
    %                    new    는 gas.TPX = {...} / gas.H (basis = 'molar')
    try
        set(th.gas, 'T', 300, 'P', th.p0, 'X', th.Xstr{1});
        enthalpy_mole(th.gas);
        th.api = 'legacy';
    catch
        th.gas.basis = 'molar';
        th.gr.basis  = 'molar';
        th.gas.TPX = {300, th.p0, th.Xstr{1}};
        th.gas.H;
        th.api = 'new';
    end

    % 분자량 [kg/kmol]: 순수 성분 상태의 평균 분자량 = 그 성분의 분자량
    th.MW = zeros(1, 5);
    for i = 1:4
        if strcmp(th.api, 'legacy')
            set(th.gas, 'T', 300, 'P', th.p0, 'X', th.Xstr{i});
            th.MW(i) = meanMolecularWeight(th.gas);
        else
            th.gas.TPX = {300, th.p0, th.Xstr{i}};
            th.MW(i) = th.gas.meanMolecularWeight;
        end
    end
    if strcmp(th.api, 'legacy')
        th.MW(5) = meanMolecularWeight(th.gr);
    else
        th.MW(5) = th.gr.meanMolecularWeight;
    end

    % 원소 행렬: 원자 유량 [kmol/s] = A * n'
    %          CO CO2 O2 N2 C(s)
    th.A = [   1   1   0  0  1;      % C 원자
               1   2   2  0  0;      % O 원자
               0   0   0  2  0];     % N 원자
end
