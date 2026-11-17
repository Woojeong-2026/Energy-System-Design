# 조건 변경과 오류 주입 학습

저장소를 본인 PC에 받은 다음 아래 HTML을 브라우저로 여십시오. GitHub 파일 화면에서는 코드가 보이며 인터랙션이 실행되지 않습니다. 파일별 화면의 **Download raw file**로 HTML 하나만 내려받아 열 수도 있습니다. ZIP이나 MATLAB 설치 없이 작동하며, 별도 웹사이트를 배포한 것은 아닙니다.

- [charging-physics.html](charging-physics.html): 압력, 셀 온도, 전압 손실, 반응기 온도, HEX 유용도, SOEC 전환율을 한 번에 하나씩 바꿉니다.
- [charging-verification.html](charging-verification.html): 정상 결과와 열회수 온도 무시, 전류 단위 오류, Nernst 평가 조성 오류, HEX 이중 차감, 탄소 누락을 비교합니다.

첫 자료에서는 전압 손실을 0.05→0.15 V로 바꾸고 셀 전력과 총 입력이 왜 다르게 반응하는지 설명해 보십시오. 다음으로 반응기 온도를 869→870°C로 바꿔 과제의 20 K 열회수 규칙을 확인합니다.

두 번째 자료에서는 ‘온도 조건을 무시한 열회수’를 선택합니다. 전체 에너지 수지는 통과해도 열회수 조건에서 실패하고 효율이 거짓 상승하는 이유가 핵심입니다. 표시된 다섯 항목 외에도 총 46개 검사를 수행한 결과를 판정에 반영합니다.

자료의 수치는 Cantera 3.0.0으로 미리 계산한 22개 조건 및 정상 포함 6개 오류 시나리오입니다. 브라우저에서 새로운 열역학 계산을 수행하거나 여러 조건을 조합하는 도구는 아닙니다. [학습 검증 기록](lesson_verification.json)과 [원시 조건 데이터](physical_lessons.csv)를 제공합니다.

데이터 재생성은 저장소 루트에서 `python reference/build_lessons.py --out results/charging_lessons`로 수행합니다. HTML의 내장 데이터는 게시 시점의 스냅샷이므로 생성 파일이 바뀌어도 자동 갱신되지 않습니다. 편집용 원본은 `fragments/`에 있습니다.
