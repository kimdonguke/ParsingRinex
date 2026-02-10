function final_MP = applyDetrending(fullTable, raw_MP, method_vec)
    % applyDetrending: 위성별, Arc별 바이어스 제거 (Mean Removal)
    % [추가 기능]: MW(Melbourne-Wübbena) 조합을 이용한 정밀 사이클 슬립 감지
    
    final_MP = nan(size(raw_MP));
    
    % -----------------------------------------------------------
    % 상수 정의 (MW 조합용)
    % -----------------------------------------------------------
    c = 299792458;
    f1 = 1575.42e6;
    f2 = 1227.60e6;
    lambda1 = c / f1;
    lambda2 = c / f2;
    
    % Wide-lane 파장 (~0.86m)
    % MW 조합은 이 파장 단위로 움직이므로, 슬립 감지에 매우 유리함
    lambda_WL = c / (f1 - f2); 
    
    % -----------------------------------------------------------
    % PRN별 루프 수행
    % -----------------------------------------------------------
    if iscategorical(fullTable.PRN) || iscell(fullTable.PRN)
         numPRN = double(string(fullTable.PRN));
    else
         numPRN = fullTable.PRN;
    end
    uPRN = unique(numPRN);

    for i = 1:length(uPRN)
        prn = uPRN(i);
        idx = find(numPRN == prn);
        
        if isempty(idx), continue; end
        
        % 해당 위성의 데이터 추출
        subTable = fullTable(idx, :);
        vals = raw_MP(idx);
        times = subTable.Time;
        methods = method_vec(idx);
        
        % 유효한 데이터만 처리 (MP 값이 있는 경우)
        validMask = ~isnan(vals);
        if sum(validMask) == 0, continue; end
        
        % 마스킹된 서브 데이터셋
        t_sub = times(validMask);
        m_sub = methods(validMask);
        v_sub = vals(validMask);
        
        % -------------------------------------------------------
        % [Step 1] MW 조합 계산 (Cycle Slip 감지용)
        % -------------------------------------------------------
        % L1, L2, P1, P2가 모두 있어야 MW 계산 가능
        % (없으면 MW 검사는 건너뛰고 기존 LLI/TimeGap만 사용)
        hasDual = ~isnan(subTable.P1(validMask)) & ...
                  ~isnan(subTable.P2(validMask)) & ...
                  ~isnan(subTable.L1(validMask)) & ...
                  ~isnan(subTable.L2(validMask));
        
        isMWSlip = false(size(v_sub)); % 초기화
        
        if any(hasDual)
            % 데이터 추출
            L1_sub = subTable.L1(validMask);
            L2_sub = subTable.L2(validMask);
            P1_sub = subTable.P1(validMask);
            P2_sub = subTable.P2(validMask);
            
            % MW 수식: MW = lambda_WL * (L1 - L2) - (f1*P1 + f2*P2)/(f1+f2)
            % (단위: 미터)
            term_Phi = lambda_WL * (L1_sub - L2_sub);
            term_Code = (f1 * P1_sub + f2 * P2_sub) / (f1 + f2);
            MW_val = term_Phi - term_Code;
            
            % 차분(Diff)을 통한 점프 감지
            % 임계값: Wide-lane 파장(0.86m)의 약 1.5배 (1.2 ~ 1.3m) 
            % 코드 노이즈를 고려하여 약간 여유 있게 설정
            mw_diff = [0; diff(MW_val)];
            
            % MW 값 자체도 노이즈가 있으므로, 튀는 값만 잡음
            isMWSlip = abs(mw_diff) > 1.2; 
            
            % 첫 번째 epoch은 차분값이 없으므로 false 처리
            isMWSlip(1) = false;
        end

        % -------------------------------------------------------
        % [Step 2] Arc Segmentation (구간 나누기)
        % -------------------------------------------------------
        
        % 1. 시간 갭 (60초 이상 끊기면 분리)
        isTimeGap = [0; diff(t_sub)] > 60;
        
        % 2. 방식 변경 (Dual <-> Triple 등)
        isMethodChange = [0; diff(m_sub) ~= 0];
        
        % 3. 기존 LLI (Loss of Lock Indicator) 확인
        if ismember('LLI1', subTable.Properties.VariableNames)
            lli1 = subTable.LLI1(validMask);
            isLLI = bitand(uint8(lli1), 1) == 1;
        else
            isLLI = false(size(v_sub));
        end
        
        % 4. (옵션) Geometry-Free 점프 or MP 자체의 거대 점프
        % MW가 있으므로 MP 자체 점프 임계값은 좀 크게 잡습니다 (10m)
        jumpDiff = [0; abs(diff(v_sub))];
        isLargeJump = jumpDiff > 10.0; 
        
        % *** [최종 조건 병합] ***
        % MW Slip이 감지되었거나, 시간이 끊겼거나, LLI가 떴거나
        isNewArc = isTimeGap | isMethodChange | isLLI | isLargeJump | isMWSlip;
        
        % 누적합으로 Arc ID 부여 (1, 1, 1, 2, 2, 3, ...)
        arcIDs = cumsum(isNewArc) + 1;
        
        % -------------------------------------------------------
        % [Step 3] Mean Removal & Short Arc Rejection
        % -------------------------------------------------------
        detrended_sub = v_sub;
        uArcs = unique(arcIDs);
        
        for k = 1:length(uArcs)
            currArcMask = (arcIDs == uArcs(k));
            nEpochs = sum(currArcMask);
            
            % [중요] 데이터 개수가 너무 적으면(10개 미만) 통계적 의미가 없고
            % 평균 제거 시 0으로 수렴하는 문제 발생 -> NaN 처리 후 스킵
            if nEpochs < 10
                detrended_sub(currArcMask) = NaN;
                continue; 
            end
            
            % 해당 구간 평균(Bias) 계산 및 제거
            bias = mean(v_sub(currArcMask), 'omitnan');
            detrended_sub(currArcMask) = v_sub(currArcMask) - bias;
            
            % (선택) 여기서 3-sigma Outlier 제거를 추가할 수도 있음
            % arc_std = std(detrended_sub(currArcMask), 'omitnan');
            % outlier_idx = abs(detrended_sub(currArcMask)) > 3 * arc_std;
            % detrended_sub(currArcMask(outlier_idx)) = NaN; 
        end
        
        % 전체 배열에 결과 매핑
        saveIdx = idx(validMask);
        final_MP(saveIdx) = detrended_sub;
    end
end