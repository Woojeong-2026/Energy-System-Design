# 시스템 구성도

각 그림은 `.mmd`(Mermaid 원본)에서 그렸고, 보고서·발표용 `.png`(2배 해상도)와 `.svg`를 같이 둡니다.
아래 코드 블록은 깃허브에서 바로 그림으로 보입니다.

## 1. 코드 구조도

`main.m` → 통합 함수 → 충전·방전 모델 → 공통 장치 함수(`system/`) → 물성 함수(`thermo/`) → Cantera 순서로 호출합니다.
파란색은 팀원 3, 초록색은 팀원 1, 주황색은 팀원 2 담당입니다. 점선은 아직 연결되지 않은 부분입니다.

```mermaid
flowchart TB
    classDef t3 fill:#E8F0FE,stroke:#3B6FD8,color:#13254d
    classDef t1 fill:#E7F6EC,stroke:#2E8B57,color:#123d25
    classDef t2 fill:#FFF1E0,stroke:#D9822B,color:#4a2a07
    classDef cfg fill:#F3F3F3,stroke:#777,color:#222

    MAIN["main.m<br/>(a)~(f) 전체 재현 · 압력/셀 온도 스윕 · 민감도 · 그림"]:::t3
    CFG["noon_config.m<br/>과제 지정 조건 (P, Tcell, TB, 효율 가정, 상수)"]:::cfg
    RT["system/run_roundtrip.m<br/>충전+방전 실행 · t_dis/t_ch 맞추기 · η_RT · 손실 내역"]:::t3

    subgraph MODES["모드 모델 (저장소 최상위)"]
        direction LR
        CH["run_charging_analysis.m<br/>calculate_charging(cfg, th)<br/>Stream 1~13 · P_in · η_ch"]:::t1
        DIS["run_discharging_analysis.m<br/>Stream 14~29 · P_out · η_dis"]:::t2
    end

    subgraph SYS["system/ — 공통 장치·정산·표 (팀원 3)"]
        direction LR
        HEX["hex_counterflow<br/>향류 열교환기 ε-모델"]:::t3
        SET["settle_heat<br/>thermal link · 전기 가열 · 방출열"]:::t3
        TAB["stream_table · balance_table · heat_table<br/>결과표 · 원소/에너지 잔차"]:::t3
    end

    subgraph TH["thermo/ — Cantera 물성 도우미 (팀원 3)"]
        direction LR
        INIT["thermo_init<br/>gri30.yaml · graphite.yaml"]:::t3
        KB["K_boudouard<br/>K_B = exp(−ΔG°/RT)"]:::t3
        EC["E_cell_CO<br/>Nernst E_rev"]:::t3
        SH["stream_H · stream_T<br/>H = Σ n·h(T), H→T 역산"]:::t3
        SP["sp_h · sp_s0 · sp_g0<br/>순수 성분 h, s°, g° (Cantera 직접 호출)"]:::t3
    end

    CANTERA[("Cantera 3.0.0<br/>MATLAB toolbox")]:::cfg

    MAIN --> CFG
    MAIN --> RT
    RT -. "통합 연결 필요<br/>(현재 charging_mode 호출)" .-> CH
    RT -. "통합 연결 필요<br/>(현재 스크립트 형태)" .-> DIS
    MODES --> SYS
    MODES --> TH
    SYS --> TH
    TH --> CANTERA
```

## 2. 충전 모드 공정 흐름도 (Figure 2(a), Stream 1~13)

노란색은 전력을 쓰거나 내는 장치, 하늘색은 냉각기입니다. 화살표 위 숫자는 스트림 번호입니다.

```mermaid
flowchart LR
    classDef elec fill:#FFF4CC,stroke:#C9A100,color:#3d3000
    classDef unit fill:#E8F0FE,stroke:#3B6FD8,color:#13254d
    classDef store fill:#EDEDED,stroke:#666,color:#222
    classDef cool fill:#E3F4F8,stroke:#2B8FA8,color:#0d3a45

    TANK[("CO2 탱크<br/>25 °C")]:::store
    MIX["혼합기"]:::unit
    HEX1["HEX 1<br/>향류 ε = 0.85"]:::unit
    H1["Heater H1<br/>(전기 가열)"]:::elec
    SOEC["rSOC — SOEC 모드<br/>Tcell = 850 °C<br/>CO2 → CO + ½O2, X = 0.70<br/>V_ch = E_rev + 0.10 V"]:::elec
    TRIM["Trim unit<br/>850 → 650 °C"]:::unit
    RX["Boudouard 반응기 + 고체 분리기<br/>TB,ch = 650 °C<br/>2CO → C(s) + CO2 (평형)"]:::unit
    RC["재순환 냉각기"]:::cool
    CC["탄소 냉각기"]:::cool
    OC["O2 냉각기"]:::cool
    CST[("탄소 저장<br/>25 °C")]:::store
    VENT["O2 방출 밸브"]:::store
    PWR(["전력 입력 P_in"]):::elec

    TANK -- "1  CO2" --> MIX
    MIX -- "2" --> HEX1
    HEX1 -- "3" --> H1
    H1 -- "4  연료극 입구" --> SOEC
    SOEC -- "5  CO/CO2" --> TRIM
    TRIM -- "6" --> RX
    RX -- "7  CO/CO2 (고온측)" --> HEX1
    HEX1 -- "8" --> RC
    RC -- "9  재순환 25 °C" --> MIX
    RX -- "10  C(s)" --> CC
    CC -- "11" --> CST
    SOEC -- "12  순수 O2" --> OC
    OC -- "13" --> VENT
    PWR -.-> SOEC
    PWR -.-> H1
    RX <-. "thermal link<br/>(donor가 20 K 이상 뜨거울 때만)" .-> SOEC
```

## 3. 방전 모드 공정 흐름도 (Figure 2(b), Stream 14~29)

```mermaid
flowchart TB
    classDef elec fill:#FFF4CC,stroke:#C9A100,color:#3d3000
    classDef unit fill:#E8F0FE,stroke:#3B6FD8,color:#13254d
    classDef store fill:#EDEDED,stroke:#666,color:#222
    classDef cool fill:#E3F4F8,stroke:#2B8FA8,color:#0d3a45

    subgraph FUEL["연료측 루프 (CO/CO2, P = 5 bar)"]
        direction LR
        CST[("탄소 저장<br/>25 °C")]:::store
        SEP["CO2 분리기<br/>25 °C 전기화학 펌프"]:::elec
        HEX3["HEX 3<br/>향류 ε = 0.85"]:::unit
        H2["Heater H2"]:::elec
        RX["Boudouard 반응기 (가스화)<br/>TB,dis = 800 °C<br/>C(s) + CO2 → 2CO (평형)"]:::unit
        TRIM["Trim unit<br/>800 → 850 °C"]:::elec
        COOL["냉각기"]:::cool
        TANK[("CO2 탱크")]:::store
    end

    SOFC["rSOC — SOFC 모드<br/>Tcell = 850 °C<br/>CO + ½O2 → CO2, U = 0.80<br/>V_dis = E_rev − 0.10 V"]:::elec
    PWR(["순출력 P_out"]):::elec

    subgraph AIRL["공기측 (O2/N2, 공급 O2 = 양론 × 3)"]
        direction LR
        AIR(["주위 공기<br/>1 bar, 25 °C"]):::store
        COMP["공기 압축기<br/>25 °C 등온 가역"]:::elec
        HEX4["HEX 4<br/>향류 ε = 0.85"]:::unit
        H3["Heater H3"]:::elec
        EXH(["배기"]):::store
    end

    CST -- "14  C(s)" --> RX
    SEP -- "15  재순환" --> HEX3
    HEX3 -- "16" --> H2
    H2 -- "17" --> RX
    RX -- "18" --> TRIM
    TRIM -- "19  연료극 입구" --> SOFC
    SOFC -- "20  고온측" --> HEX3
    HEX3 -- "21" --> COOL
    COOL -- "22" --> SEP
    SEP -- "23  순수 CO2" --> TANK

    AIR -- "24" --> COMP
    COMP -- "25" --> HEX4
    HEX4 -- "26" --> H3
    H3 -- "27  공기극 입구" --> SOFC
    SOFC -- "28  고온측" --> HEX4
    HEX4 -- "29" --> EXH

    SOFC -.-> PWR
    SOFC -. "thermal link<br/>셀 발열 → 가스화 흡열" .-> RX
```
