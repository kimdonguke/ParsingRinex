function final_MP = applyDetrending(fullTable, raw_MP, method_vec)
    % applyDetrending: 위성별, Arc별 바이어스 제거 (Mean Removal)
    % [업데이트]: MW 조합 Slip 감지 + 상세 로그(fprintf) 출력 추가
    
    fprintf('========================================\n');
    fprintf('[Start] Detrending 및 Arc 분리 작업 시작\n');
    fprintf('========================================\n');
    
    final_MP = nan(size(raw_MP));
    
    % 상수 정의 (MW 조합용)
    c = 299792458;
    f1 = 1575.42e6;
    f2 = 1227.60e6;
    lambda1 = c / f1;
    lambda2 = c / f2;
    lambda_WL = c / (f1 - f2); % Wide-lane 파장 (~0.86m)
    
    % PRN 리스트 추출
    if iscategorical(fullTable.PRN) || iscell(fullTable.PRN)
         numPRN = double(string(fullTable.PRN));
    else
         numPRN = fullTable.PRN;
    end
    uPRN = unique(numPRN);
    
    total_sats = length(uPRN);
    fprintf('총 %d개의 위성 데이터가 감지되었습니다.\n', total_sats);

    for i = 1:total_sats
        prn = uPRN(i);
        idx = find(numPRN == prn);
        
        if isempty(idx), continue; end
        
        % 진행 상황 표시 (줄바꿈 없이 시작)
        fprintf('  > PRN %02d 처리 중... ', prn);
        
        % 데이터 추출
        subTable = fullTable(idx, :);
        vals = raw_MP(idx);
        times = subTable.Time;
        methods = method_vec(idx);
        
        validMask = ~isnan(vals);
        if sum(validMask) == 0
            fprintf('[Skip] 유효한 MP 데이터 없음.\n');
            continue; 
        end
        
        % 서브 데이터셋
        t_sub = times(validMask);
        m_sub = methods(validMask);
        v_sub = vals(validMask);
        
        % -------------------------------------------------------
        % [Step 1] MW 조합 계산 (Cycle Slip 감지)
        % -------------------------------------------------------
        hasDual = ~isnan(subTable.P1(validMask)) & ...
                  ~isnan(subTable.P2(validMask)) & ...
                  ~isnan(subTable.L1(validMask)) & ...
                  ~isnan(subTable.L2(validMask));
        
        isMWSlip = false(size(v_sub));
        mw_check_msg = "Skipped (L2/P2 missing)";
        
        if any(hasDual)
            L1_sub = subTable.L1(validMask);
            L2_sub = subTable.L2(validMask);
            P1_sub = subTable.P1(validMask);
            P2_sub = subTable.P2(validMask);
            
            term_Phi = lambda_WL * (L1_sub - L2_sub);
            term_Code = (f1 * P1_sub + f2 * P2_sub) / (f1 + f2);
            MW_val = term_Phi - term_Code;
            
            mw_diff = [0; diff(MW_val)];
            isMWSlip = abs(mw_diff) > 1.2; % 임계값 1.2m
            isMWSlip(1) = false;
            
            mw_check_msg = "Done";
        end

        % -------------------------------------------------------
        % [Step 2] Arc Segmentation
        % -------------------------------------------------------
        isTimeGap = [0; diff(t_sub)] > 60;
        isMethodChange = [0; diff(m_sub) ~= 0];
        
        if ismember('LLI1', subTable.Properties.VariableNames)
            lli1 = subTable.LLI1(validMask);
            isLLI = bitand(uint8(lli1), 1) == 1;
        else
            isLLI = false(size(v_sub));
        end
        
        jumpDiff = [0; abs(diff(v_sub))];
        isLargeJump = jumpDiff > 10.0; 
        
        % 최종 분리 조건
        isNewArc = isTimeGap | isMethodChange | isLLI | isLargeJump | isMWSlip;
        arcIDs = cumsum(isNewArc) + 1;
        
        % 통계 집계
        count_slips = sum(isMWSlip);
        count_gaps = sum(isTimeGap);
        
        % -------------------------------------------------------
        % [Step 3] Mean Removal & Log
        % -------------------------------------------------------
        detrended_sub = v_sub;
        uArcs = unique(arcIDs);
        valid_arcs = 0;
        dropped_arcs = 0;
        
        for k = 1:length(uArcs)
            currArcMask = (arcIDs == uArcs(k));
            nEpochs = sum(currArcMask);
            
            if nEpochs < 10
                detrended_sub(currArcMask) = NaN;
                dropped_arcs = dropped_arcs + 1;
                continue; 
            end
            
            bias = mean(v_sub(currArcMask), 'omitnan');
            detrended_sub(currArcMask) = v_sub(currArcMask) - bias;
            valid_arcs = valid_arcs + 1;
        end
        
        % 결과 저장
        saveIdx = idx(validMask);
        final_MP(saveIdx) = detrended_sub;
        
        % [완료 로그 출력]
        fprintf('완료. [MW검사: %s] [Slip: %d회] [Arc: %d 생성 / %d 제거]\n', ...
            mw_check_msg, count_slips, valid_arcs, dropped_arcs);
    end
    
    fprintf('========================================\n');
    fprintf('[End] 모든 위성 처리 완료.\n');
    fprintf('========================================\n');
end