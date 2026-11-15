# Noon Energy 탄소-산소 에너지 저장 시스템 — Charging System (담당자 A)

Energy System Design (38609) 2026 Fall Term Project 중 **Charging mode (Figure 2(a), Stream 1~13)** 모델입니다.
기준 문서는 `Term_Project_2026_Noon_Energy_revV3.pdf` 입니다.

## 파일

| 파일 | 내용 |
|---|---|
| `main_charging_A.m` | 실행 스크립트. 기준 조건, 압력 스윕, 셀 온도 스윕을 계산하고 표·그림을 `results/`에 저장 |
| `noon_charging.m` | Charging 모델 함수 `R = noon_charging(P_bar, Tcell_C, TBch_C, th)` |

## 실행

1. MATLAB R2026a에서 Cantera 3.0.0 MATLAB toolbox 경로를 등록합니다.
2. 이 폴더에서 `main_charging_A`를 실행합니다.

Cantera 호출은 `noon_charging.m`의 `load_thermo`, `h_sp`, `s0_sp` 세 함수에만 있습니다(legacy 문법: `Solution`, `set`, `enthalpy_mole`, `entropy_mole`, `molecularWeights`, `speciesIndex`).

## 담당 범위와 출력

| 역할 | 출력 |
|---|---|
| A1 Stream 1~13 결과표 | `R.streams`, `results/charging_streams_*.csv` |
| A3 충전 화학평형·재순환 | `R.eq` |
| A4 열부하·회수열·방출열·전력 | `R.duty`, `R.power`, `R.cell` |
| A5 원소·에너지수지 (충전분) | `R.balance`, `results/charging_balance_*.csv` |
| A6 압력별 (충전분) | `results/charging_sweep_pressure.csv/.png` |
| A7 셀 온도별 (충전분) | `results/charging_sweep_Tcell.csv/.png` |

### 통합용 인터페이스 (담당자 B, C)

```matlab
R = noon_charging(5, 850, 650);       % P [bar], Tcell [degC], TB,ch [degC]
R.power.P_in    % 충전 입력 전력 [W]
R.eta_ch        % 충전 효율 [-]
R.NC_rate       % 탄소 저장 속도 [kmol/s]  (방전 소비 속도와 매칭)
R.S(k).n        % Stream k 몰유량 [CO CO2 O2 N2 C(s)] [kmol/s],  R.S(k).T [K],  R.S(k).P [Pa]
```

왕복 효율은 `eta_RT = eta_ch * eta_dis = P_out,dis / P_in,ch` (두 모드 모두 탄소 1 kmol/s 기준)입니다.

## 계산 방법

### 물성

- 기체 CO, CO2, O2, N2: Cantera `gri30.yaml` (NASA 7계수, 이상기체)
- 고체 탄소: Cantera `graphite.yaml`의 C(gr), 활동도 1
- 엔탈피는 생성엔탈피를 포함한 절대값 `h_i(T)`를 쓰므로 반응열이 따로 필요 없습니다.
- 표준 Gibbs 에너지: `g°_i(T) = h_i(T) − T·s°_i(T, 1 bar)`
- 상수: R = 8.314 J/(mol K), F = 96485 C/mol, p° = 1 bar, LHV_C = 393.5 kJ/mol

### 부호 규약

- `Q > 0`: 장치로 들어가는 열(가열), `Q < 0`: 빼내는 열(냉각)
- `W > 0`: 시스템으로 들어가는 전력
- 에너지 잔차 = ΣH_in + Q + W − ΣH_out, 원소 잔차 = Σ(원자 in) − Σ(원자 out)

### 재순환 루프 (해석해)

Stream 2의 CO, CO2 유량을 a, b, 탄소 석출량을 ξ라 하면

- SOEC: Stream 5는 CO = a + 0.7b, CO2 = 0.3b / Stream 12는 O2 = 0.35b
- 반응기: Stream 7은 CO = a + 0.7b − 2ξ, CO2 = 0.3b + ξ
- 혼합기(전량 재순환, purge 없음): CO 수지에서 ξ = 0.35b, CO2 수지에서 0.7b = 1 + ξ

따라서 **ξ = 1 kmol/s, b = 20/7 kmol/s**이고, 반응기 출구 CO2 유량은 c = 13/7 kmol/s로 압력·온도와 무관합니다.
남은 미지수 a는 Boudouard 평형으로 정해집니다.

```
K_B(T) = exp(−ΔG°_B/RT) = (p_CO2/p°)/(p_CO/p°)²,   ΔG°_B = g°_C + g°_CO2 − 2g°_CO
(K_B·P/p°)·a² − c·a − c² = 0   →   a = c·[1 + √(1 + 4·K_B·P/p°)] / (2·K_B·P/p°)
```

반복 계산(tear stream)이 필요 없고, 평형 잔차 `ln K − ln Q`를 `R.eq.res_lnK`로 확인합니다.

### 에너지수지 순서

1. 혼합기(단열): `H2(T2) = H1 + H9` → T2 (`fzero`)
2. HEX 1(단열 향류, ε = 0.85): `Q_max = min[H_h(T_h,in) − H_h(T_c,in), H_c(T_h,in) − H_c(T_c,in)]`, `Q = ε·Q_max`, 출구 온도는 각 유체 에너지수지를 `fzero`로 풀이. 내부 21개 지점에서 T_hot − T_cold > 0을 확인(`R.duty.dTmin_HEX1`)
3. Heater H1: `Q_H1 = H4 − H3`
4. SOEC(등온): `E_rev = −ΔG°_ox/2F + (RT/2F)·ln[(p_CO/p°)(p_O2/p°)^½/(p_CO2/p°)]` (Stream 5, 12 조성), `V_ch = E_rev + 0.10`, `I = 2F·ṅ_CO`, `W_cell = I·V_ch`, `Q_cell = H5 + H12 − H4 − W_cell`
5. Trim unit: `Q_trim = H6 − H5`
6. 반응기(등온): `Q_rx = H7 + H10 − H6`
7. 냉각기 3개(재순환, 탄소, O2): 25 °C까지 냉각, 전량 방출
8. Thermal link: donor(Q<0) 온도가 receiver(Q>0)보다 20 K 이상 높을 때만 `min(0.85·|Q_donor|, Q_receiver)`
9. `P_in = W_cell + Σ(회수 후 남은 양의 열부하)`, `η_ch = ξ·LHV_C / P_in`

## 참고 수치 (MATLAB 실행 전 대조용)

아래는 **MATLAB/Cantera로 얻은 값이 아니라**, 같은 식을 `gri30.yaml`·`graphite.yaml`과 동일한 NASA 계수로 Python에 옮겨 계산한 값입니다. MATLAB 실행 결과와 대조하고, 보고서에는 MATLAB 결과를 쓰세요.

기준 조건 (P = 5 bar, Tcell = 850 °C, TB,ch = 650 °C)

| Stream | T [°C] | CO [kmol/s] | CO2 [kmol/s] | O2 [kmol/s] | C(s) [kmol/s] |
|---|---|---|---|---|---|
| 1 | 25 | 0 | 1.0000 | 0 | 0 |
| 2 | 25 | 0.5312 | 2.8571 | 0 | 0 |
| 3 | 412.4 | 0.5312 | 2.8571 | 0 | 0 |
| 4 | 850 | 0.5312 | 2.8571 | 0 | 0 |
| 5 | 850 | 2.5312 | 0.8571 | 0 | 0 |
| 6 | 650 | 2.5312 | 0.8571 | 0 | 0 |
| 7 | 650 | 0.5312 | 1.8571 | 0 | 0 |
| 8 | 134.1 | 0.5312 | 1.8571 | 0 | 0 |
| 9 | 25 | 0.5312 | 1.8571 | 0 | 0 |
| 10 | 650 | 0 | 0 | 0 | 1.0000 |
| 11 | 25 | 0 | 0 | 0 | 1.0000 |
| 12 | 850 | 0 | 0 | 1.0000 | 0 |
| 13 | 25 | 0 | 0 | 1.0000 | 0 |

| 항목 | 값 |
|---|---|
| K_B (650 °C) | 3.143 |
| E° / E_rev / V_ch | 0.958 / 1.049 / 1.149 V |
| W_cell | 443.6 MW |
| Q_cell (전기 가열) | +120.7 MW |
| Q_H1 (전기 가열) | +73.6 MW |
| Q_trim / Q_rx | −26.2 / −171.3 MW (방출) |
| HEX 1 회수열 | 54.8 MW |
| Thermal link | 0 MW |
| 냉각기 방출 (O2 / C / 재순환) | 27.0 / 10.2 / 9.7 MW |
| **P_in** | **637.9 MW** |
| **η_ch** | **61.7 %** |

스윕 요약

| P [bar] | y_CO (Stream 7) | n4 [kmol/s] | E_rev [V] | W_cell [MW] | Q_cell [MW] | Q_H1 [MW] | P_in [MW] | η_ch [%] |
|---|---|---|---|---|---|---|---|---|
| 1 | 0.427 | 4.241 | 1.024 | 434.0 | 130.3 | 81.8 | 646.0 | 60.9 |
| 2 | 0.327 | 3.760 | 1.034 | 437.6 | 126.7 | 77.2 | 641.4 | 61.3 |
| 5 | 0.222 | 3.388 | 1.049 | 443.6 | 120.7 | 73.6 | 637.9 | 61.7 |
| 10 | 0.163 | 3.219 | 1.063 | 448.8 | 115.5 | 72.0 | 636.3 | 61.8 |

| Tcell [°C] | E_rev [V] | W_cell [MW] | Q_cell [MW] | Q_H1 [MW] | Q_trim [MW] | P_in [MW] | η_ch [%] |
|---|---|---|---|---|---|---|---|
| 750 | 1.086 | 457.8 | 107.2 | 56.1 | −13.0 | 621.2 | 63.3 |
| 800 | 1.068 | 450.7 | 114.0 | 64.8 | −19.6 | 629.5 | 62.5 |
| 850 | 1.049 | 443.6 | 120.7 | 73.6 | −26.2 | 637.9 | 61.7 |
| 900 | 1.031 | 436.5 | 127.3 | 82.5 | −32.9 | 646.3 | 60.9 |

## 발표용 해석 포인트 (Charging)

- **Thermal link는 충전 모드에서 작동하지 않습니다.** 셀은 열중립 전압(약 1.46 V)보다 낮은 1.15 V에서 운전되어 흡열(+121 MW)이고, 발열하는 반응기(−171 MW)는 650 °C로 셀(850 °C)보다 차갑습니다. donor가 20 K 이상 뜨거워야 한다는 조건을 만족하지 못하므로 반응기 열 171 MW는 전부 방출되고, 셀 부족열은 전기로 채웁니다.
- 그 결과 `W_cell + Q_cell`은 셀 반응 엔탈피(약 564 MW)로 고정됩니다. 전압이 올라가면 그만큼 셀 가열이 줄어 서로 상쇄되므로, **P_in ≈ 564 MW + Q_H1**이고 효율은 사실상 H1 열부하가 결정합니다.
- **압력↑**: Boudouard 평형이 석출 쪽으로 이동해 재순환 CO가 줄고(n4: 4.24 → 3.22 kmol/s), 예열 부하가 줄어 효율이 60.9 → 61.8 %로 소폭 상승합니다. E_rev 상승(1.024 → 1.063 V)은 위 상쇄 때문에 P_in에 영향이 없습니다.
- **셀 온도↑**: 반응기 온도가 고정이라 평형·순환량은 그대로이고, H1 부하와 Trim 방출열만 늘어 효율이 63.3 → 60.9 %로 떨어집니다.
- **개선 방향**: 반응기 발열(650 °C, 171 MW)과 Trim 방출열(26 MW)을 Stream 3 예열에 쓰면 H1 부하(74 MW)의 대부분을 대체할 수 있습니다. 셀과 반응기 온도 차를 줄이는 것도 같은 효과입니다.
