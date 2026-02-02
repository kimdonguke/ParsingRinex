%Obs 받고, T_sv는 수렴시켜야되자너 위성에서 나오는 tsv 세 개 다 수렴 시키고 그러고 여기로 넘겨서 각 위성마다 MP combination
%적용 시킨다?
%적용시키고 나오면 계단식(lambda integer ambiguty + reciever noise, multipath error + residual error만 남는데)
%이때 수신기 노이즈만 보고 싶으니까 integer ambiguty 구간 적분으로 없애고, multipath도 미소 시간에 대해선
%바이어슨데 예가 integer ambiguty랑 구분이 안 될거같은데

MP_triple = calcTripleFreqMultipath(return_OBS,'Code');

function obsTable = calcTripleFreqMultipath(obsTable, measType)
    % calcTripleFreqMultipath
    % 3주파수(L1, L2, L5) 선형 결합을 통한 다중경로 및 잡음 추출
    %
    % 입력:
    %   obsTable : P1, P2, P5 (또는 L1, L2, L5) 컬럼이 있는 Flat Table
    %   measType : 'Code' (기본값, 의사거리 P 사용) 또는 'Phase' (반송파 L 사용)
    %
    % 출력:
    %   obsTable : 'MP_Triple' 컬럼이 추가된 테이블
    
    if nargin < 2, measType = 'Code'; end

    % 1. 시스템별 계수 설정 (Table 20.2 기준)
    % 현재는 GPS만 구현 (확장 가능)
    % 계수: betaA(L1), betaB(L2), betaC(L5)
    % 조건: Sum=0 (GF), Ion-Free, Norm=1
    
    coeff.GPS = [0.142, -0.767, 0.626];     % L1, L2, L5
    coeff.GAL = [0.085, -0.746, 0.661];     % E1, E5b, E5a
    coeff.BDS = [0.183, -0.781, 0.597];     % B1, B3, B2
    
    % 여기서는 GPS 계수 사용
    beta = coeff.GPS; 
    
    fprintf('Triple-Freq Combination (%s): GPS L1/L2/L5 계수 [%.3f, %.3f, %.3f] 적용\n', ...
        measType, beta(1), beta(2), beta(3));

    % 2. 입력 데이터 로드 (Step 1)
    % Code(P) 또는 Phase(L) 데이터 추출
    if strcmpi(measType, 'Code')
        % P1, P2, P5가 모두 존재해야 함
        if ~ismember('P5', obsTable.Properties.VariableNames)
            error('테이블에 P5(L5 의사거리) 데이터가 없습니다.');
        end
        DATA_A = obsTable.P1;
        DATA_B = obsTable.P2;
        DATA_C = obsTable.P5;
        
    elseif strcmpi(measType, 'Phase')
        % 위상 조합의 경우 (미터 단위 변환 필요 여부는 논문 수식 정의에 따름)
        % 보통 위상 조합 식은 Cycle 단위 혹은 Meter 단위 둘 다 가능하나, 
        % 여기서는 Code와 동일한 차원(Meter)으로 맞추기 위해 파장 곱셈 권장.
        % 하지만 제시해주신 수식은 phi 자체의 선형결합이므로 입력 그대로 사용합니다.
        DATA_A = obsTable.L1;
        DATA_B = obsTable.L2;
        DATA_C = obsTable.L5;
    else
        error('measType은 "Code" 또는 "Phase"여야 합니다.');
    end
    
    % 3. 선형 결합 계산 (Step 3)
    % oMP = beta1*P1 + beta2*P2 + beta3*P5
    % 행렬 연산(Vectorization)으로 한 번에 계산
    raw_MP = beta(1) * DATA_A + beta(2) * DATA_B + beta(3) * DATA_C;
    
    % 4. 바이어스 제거 (Detrending) (Step 4)
    % 위성별(PRN), 연속 구간(Arc)별로 평균을 구해서 빼줌
    
    obsTable.MP_Triple = nan(height(obsTable), 1);
    
    uniquePRN = unique(obsTable.PRN);
    
    for i = 1:length(uniquePRN)
        prn = uniquePRN(i);
        idx = find(obsTable.PRN == prn);
        
        if isempty(idx), continue; end
        
        % 해당 위성의 데이터
        subData = raw_MP(idx);
        subTime = obsTable.Time(idx);
        
        % ---------------------------------------------------
        % [Arc 분리 로직]
        % 데이터가 끊기는 구간(Time Gap)이나 NaN이 있는 구간 분리
        % ---------------------------------------------------
        % 1) NaN이 아닌 유효 데이터만 처리
        validMask = ~isnan(subData);
        if sum(validMask) == 0, continue; end
        
        % 유효 인덱스 내에서 시간 차이 계산
        validTime = subTime(validMask);
        dt = [0; diff(validTime)];
        
        % 60초 이상 끊기면 새로운 Arc로 간주
        isNewArc = dt > 60; 
        arcIDs = cumsum(isNewArc) + 1;
        
        % 그룹별 바이어스 제거 (Detrending)
        % 여기서는 가장 Robust한 '평균 제거(Mean Removal)' 사용
        % 필요 시 'detrend(..., 1)' 등으로 선형 추세 제거 가능
        
        temp_result = subData(validMask); % 결과 담을 벡터
        
        % 각 Arc 별로 평균 제거
        uArcs = unique(arcIDs);
        for k = 1:length(uArcs)
            aID = uArcs(k);
            inArc = (arcIDs == aID);
            
            arcValues = temp_result(inArc);
            
            % 평균 제거 (Bias 상쇄)
            bias = mean(arcValues, 'omitnan');
            temp_result(inArc) = arcValues - bias;
        end
        
        % 원래 테이블 위치에 저장
        saveIdx = idx(validMask);
        obsTable.MP_Triple(saveIdx) = temp_result;
    end
    
    % 5. 결과 통계 출력 (Step 5)
    rms_val = rms(obsTable.MP_Triple, 'omitnan');
    fprintf('완료: 전체 데이터 RMS = %.4f (m/cycle)\n', rms_val);
end