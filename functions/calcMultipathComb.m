function [obsTable, outlierTable] = calcMultipathComb(obsTable)
    % 수정된 함수: 무조건 Dual Frequency (L1, L2) 사용
    % L2가 없는 경우 outlierTable로 반환
    
    % 1. 데이터 유효성 마스킹
    hasP1 = ~isnan(obsTable.P1);
    hasL1 = ~isnan(obsTable.L1);
    hasL2 = ~isnan(obsTable.L2);
    
    % [Case Main] Dual Frequency 계산 가능 (P1, L1, L2 모두 존재)
    % Triple 가능 여부와 상관없이 무조건 L1, L2만 봅니다.
    mask_Dual = hasP1 & hasL1 & hasL2;
    
    % [Case Outlier] P1, L1은 있는데 L2가 없는 경우 (이상치 확인용)
    mask_NoL2 = hasP1 & hasL1 & ~hasL2;
    
    % -------------------------------------------------------
    % 2. 이상치 데이터 분리 저장
    % -------------------------------------------------------
    outlierTable = obsTable(mask_NoL2, :);
    
    % -------------------------------------------------------
    % 3. Dual Frequency MP 계산 수행
    % -------------------------------------------------------
    % 결과 저장용 벡터 및 플래그 초기화
    raw_MP = nan(height(obsTable), 1);
    method_flag = nan(height(obsTable), 1); % 2: Dual Only
    
    if any(mask_Dual)
        fprintf('Processing Dual-Freq MP (L1-L2) (%d epochs)...\n', sum(mask_Dual));
        
        % Dual 계산 함수 호출
        vals_Dual = calcDualFreqMP(obsTable(mask_Dual, :));
        
        % 결과 할당
        raw_MP(mask_Dual) = vals_Dual;
        method_flag(mask_Dual) = 2; % 2번 방식으로 고정
    end
    
    % -------------------------------------------------------
    % 4. 바이어스 제거 (Detrending / Levelling)
    % -------------------------------------------------------
    % 정상적으로 계산된 데이터에 대해서만 바이어스 제거 수행
    obsTable.MP_raw = raw_MP;
    obsTable.MP_Hybrid = applyDetrending(obsTable, raw_MP, method_flag);
    
    obsTable.MP_Method = method_flag;
    
    fprintf('전체 파이프라인 완료 (Dual Only).\n');
end