# 충전 담당 시작 안내

이 폴더는 충전 설계 의도, 코드 연결, 검산 증거와 발표 준비 자료를 모읍니다. 실행 코드는 저장소 루트의 `modes/`, `charging/`, `reference/`, `tests/`에 통합되어 있습니다.

## 1. 가장 먼저 실행할 것

MATLAB R2026a에 Cantera 3.0.0 MATLAB toolbox를 설치·등록한 뒤, 저장소 루트를 Current Folder로 선택합니다.

```matlab
main_charging
```

이 명령은 기준 조건, 압력 1·2·5·10 bar, 셀 온도 750·800·850·900°C를 계산하고 각 결과를 검산합니다. 기준/P5/T850은 같은 조건이므로 9개 저장 케이스에 서로 다른 조건 7개가 있습니다. 이어서 유량 배율·입력 오류·오류 주입·열회수 문턱 회귀 시험을 실행합니다. 실패 시 성공 메시지를 출력하지 않고 중단합니다.

결과 위치:

```text
results/charging_review/
  summary.csv                  조건별 성능과 AllChecksPass
  charging_sweeps.png           압력·온도 비교
  base/                        P1, P2, P5, P10, T750, T800, T850, T900도 동일
    streams.csv                스트림 1–13
    reported_balances.csv       모델이 보고한 장치별 수지
    heat.csv                   열부하 / 회수열 / 전기가열 / 방출열
    independent_checks.csv     원시 스트림으로 다시 계산한 검산 결과
    hex_profile.csv            향류 HEX 내부 온도차 201점
    raw.json                   Python 독립 대조 입력
    result.mat                 R, V, cfg, modelPath
```

`main_charging`은 충전 전용입니다. 과제 전체 (a)~(f)와 방전·왕복 효율을 재현하는 기존 실행 명령은 **`main`**입니다. 두 실행기는 같은 `modes/charging_mode.m`을 호출합니다.

## 2. 기존 팀 코드에 연결되는 방식

```text
main.m → system/run_roundtrip.m → modes/charging_mode.m → 결과 R
                               → modes/discharging_mode.m

main_charging.m → charging/run_charging_review.m → 같은 charging_mode
                                               → verify_charging
                                               → export_charging_json
                → tests/test_charging_review.m
```

`cfg`는 `noon_config.m`의 공통 조건, `th`는 `thermo_init.m`의 물성 객체입니다. 반환값 `R`의 필드·단위·스트림 번호는 기존 인터페이스를 유지합니다. 따라서 기존 main에 별도 복사 코드나 후보 선택 스위치를 넣지 않습니다.

| 파일 | 읽을 때의 질문 |
|---|---|
| [charging_mode.m](../../modes/charging_mode.m) | 수지·평형·전기화학·열부하가 어떤 순서로 상태를 정하는가? |
| [verify_charging.m](../../charging/verify_charging.m) | 원시 스트림이 보존법칙과 지정된 물리식을 만족하는가? |
| [charging_reference.py](../../reference/charging_reference.py) | 다른 해법으로 계산해도 같은 값인가? |
| [test_reference.py](../../tests/test_reference.py) | 단위·열회수·유량을 틀리게 만들면 실제로 잡히는가? |

기존 팀 코드가 제공한 충전 해석식·공통 물성·열회수 모델을 바탕으로 입력 검사, 별도 검산기, 독립 계산, 시험, 지침서를 추가했습니다. 자세한 기여 구분은 [학습 지침서](CHARGING_STUDY_GUIDE_KO.md)에 있습니다.

## 3. 독립 Python 계산과 MATLAB 대조

아래 명령은 **저장소 루트의 터미널**에서 실행합니다. Python 3.11을 권장합니다. Windows에서 `python`이 다른 버전이면 가상환경 생성에 `py -3.11`을 사용하십시오.

```sh
python3.11 -m venv .venv
```

Windows PowerShell은 `.venv\Scripts\Activate.ps1`, macOS/Linux는 `source .venv/bin/activate`로 활성화한 뒤 실행합니다.

```sh
python -m pip install -r reference/requirements.txt
python -m unittest discover -s tests -p 'test_*.py' -v
python reference/charging_reference.py --out results/charging_reference
```

MATLAB에서 `main_charging`을 실행한 다음, 기준 케이스의 실제 출력과 대조합니다.

```sh
python reference/charging_reference.py --compare results/charging_review/base/raw.json --out results/charging_comparison/base
```

`46/46 checks passed`와 종료 코드 0이 통과 조건입니다. 실패하면 `comparison.csv`의 `Pass=False` 행을 확인합니다. 다른 8개 저장 케이스도 `base`를 해당 이름으로 바꿔 대조합니다. Python 자신이 만든 JSON을 대조한 것은 MATLAB 교차 검증으로 세지 않습니다.

선택 사항인 그림·학습 데이터 생성:

```sh
python -m pip install -r reference/requirements-plots.txt
python reference/plot_reference.py --results results/charging_reference
python reference/build_lessons.py --out results/charging_lessons
```

## 4. 이해·발표 준비 순서

1. [설계 가이드](DESIGN_GUIDE_KO.md): 담당 범위, 공정도, 손계산, 열회수 조건.
2. [학습 지침서](CHARGING_STUDY_GUIDE_KO.md): 왜 이 식·코드·인터페이스를 선택했는가.
3. [인터랙티브 학습](interactive/README.md): 조건 변경과 의도적 오류를 비교.
4. [검증 기록](VERIFICATION_KO.md): 실제 확인한 항목과 아직 실행하지 못한 항목.
5. [발표 노트](PRESENTATION_NOTES_KO.md): 8장 구성, 90초 설명, 예상 질문.

기준 독립 계산은 Pin 약 **637.959 MW**, 충전 효율 약 **61.6811%**입니다. [참고 결과](reference_results/)는 Python/Cantera 실행값이며 MATLAB 실행 결과가 아닙니다. 이 작업 환경에는 MATLAB이 없어 MATLAB 런타임과 실제 MATLAB JSON 대조는 미실행입니다. 정상상태 공정 계산의 검증 범위와 실제 장치 크기·반응속도 검증을 구분해 설명하십시오.
