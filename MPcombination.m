%Obs 받고, T_sv는 수렴시켜야되자너 위성에서 나오는 tsv 세 개 다 수렴 시키고 그러고 여기로 넘겨서 각 위성마다 MP combination
%적용 시킨다?
%적용시키고 나오면 계단식(lambda integer ambiguty + reciever noise, multipath error + residual error만 남는데)
%이때 수신기 노이즈만 보고 싶으니까 integer ambiguty 구간 적분으로 없애고, multipath도 미소 시간에 대해선
%바이어슨데 예가 integer ambiguty랑 구분이 안 될거같은데



%%% main %%%

MPtable = calcMultipathComb(return_OBS);
MPcell = groupingTable(MPtable,'PRN');

%%%
figure
title("whole")
hold on
grid on
% plot(MPcell{1}.MP_Hybrid)
for i = 1:32
    plot(MPcell{i}.Time,MPcell{i}.MP_Hybrid,'o',MarkerSize=2)
end

%%%
figure
title("triple")
hold on
grid on
% plot(MPcell{1}.MP_Hybrid)
for i = 1:32
    plot(MPcell{i}.Time(MPcell{i}.MP_Method == 3),MPcell{i}.MP_Hybrid(MPcell{i}.MP_Method == 3),'o',MarkerSize=2)
end
%%%
figure
title("dual")
hold on
grid on
% plot(MPcell{1}.MP_Hybrid)
for i = 1:32
    plot(MPcell{i}.Time(MPcell{i}.MP_Method == 2),MPcell{i}.MP_Hybrid(MPcell{i}.MP_Method == 2),'o',MarkerSize=2)
end
% fprintf("whole : %f\n dual : %f\ntriple : %f\n",...
%     rms(MP_triple.MP_Hybrid(~isnan(MP_triple.MP_Method)),...
%     rms(MP_triple.MP_Hybrid(MP_triple.MP_Method == 2)),...
%     rms(MP_triple.MP_Hybrid(MP_triple.MP_Method == 3))))



function obsTable = calcMultipathComb(obsTable)
    % 1. 데이터 유효성 마스킹 (어떤 함수를 쓸지 결정)
    % L5(P5), L1(P1,L1), L2(L2) 데이터 존재 여부 확인
    hasP1 = ~isnan(obsTable.P1);
    hasP2 = ~isnan(obsTable.P2);
    hasP5 = ~isnan(obsTable.P5);
    hasL1 = ~isnan(obsTable.L1);
    hasL2 = ~isnan(obsTable.L2);
    
    % [Case A] Triple 가능한 경우: P1, P2, P5 모두 존재
    mask_Tri = hasP1 & hasP2 & hasP5;
    
    % [Case B] Triple 불가능하지만 Dual은 가능한 경우: (Not Triple) & (P1, L1, L2 존재)
    mask_Dual = (~mask_Tri) & (hasP1 & hasL1 & hasL2);
    
    % -------------------------------------------------------
    % 2. 하위 함수 호출 (계산 수행)
    % -------------------------------------------------------
    % 결과 저장용 벡터 및 플래그 초기화
    raw_MP = nan(height(obsTable), 1);
    method_flag = nan(height(obsTable), 1); % 3: Triple, 2: Dual
    
    % [Triple 실행] 해당하는 행만 추출하여 함수로 전달
    if any(mask_Tri)
        fprintf('Processing Triple-Freq MP (%d epochs)...\n', sum(mask_Tri));
        % 부분 테이블 전달 -> 계산 결과 반환
        vals_Tri = calcTripleFreqMP(obsTable(mask_Tri, :));
        
        % 결과 할당
        raw_MP(mask_Tri) = vals_Tri;
        method_flag(mask_Tri) = 3; 
    end
    
    % [Dual 실행] 해당하는 행만 추출하여 함수로 전달
    if any(mask_Dual)
        fprintf('Processing Dual-Freq MP (%d epochs)...\n', sum(mask_Dual));
        vals_Dual = calcDualFreqMP(obsTable(mask_Dual, :));
        
        raw_MP(mask_Dual) = vals_Dual;
        method_flag(mask_Dual) = 2;
    end
    
    % -------------------------------------------------------
    % 3. 바이어스 제거 (Detrending / Levelling)
    % -------------------------------------------------------
    % 계산 방식(Method)이 바뀌거나 시간이 끊기면 Bias가 달라지므로
    % 이를 고려하여 최종 Levelling 수행
    
    obsTable.MP_Hybrid = applyDetrending(obsTable, raw_MP, method_flag);
    obsTable.MP_Method = method_flag;
    
    fprintf('전체 파이프라인 완료.\n');
end

function mp_values = calcDualFreqMP(subTable)
    % 입력: Dual 계산 대상인 행들만 담긴 subTable
    % 출력: 계산된 Raw MP 값 벡터
    
    % 상수 정의
    c = 299792458;
    f1 = 1575.42e6;
    f2 = 1227.60e6;
    
    lambda1 = c / f1;
    lambda2 = c / f2;
    alpha = (f1 / f2)^2;
    
    % 계수 (MP1 수식)
    % MP1 = P1 - Phi1 - 2/(alpha-1)*(Phi1 - Phi2)
    m_factor = 2 / (alpha - 1); 
    
    % 데이터 추출
    P1 = subTable.P1;
    Phi1 = subTable.L1 * lambda1; % Cycle -> Meter
    Phi2 = subTable.L2 * lambda2;
    
    % 벡터 연산 수행
    mp_values = P1 - Phi1 - m_factor * (Phi1 - Phi2);
end

function mp_values = calcTripleFreqMP(subTable)
    % 입력: Triple 계산 대상인 행들만 담긴 subTable
    % 출력: 계산된 Raw MP 값 벡터
    
    % GPS Triple Frequency 계수 (Geometry-free, Ion-free, Noise-normalized)
    % beta = [L1, L2, L5]
    b1 =  0.14150974;
    b2 = -0.76716082;
    b3 =  0.62565108;
    
    % 데이터 추출 (Code Measurement)
    P1 = subTable.P1;
    P2 = subTable.P2;
    P5 = subTable.P5;
    
    % 선형 결합 수행
    mp_values = b1 * P1 + b2 * P2 + b3 * P5;
end

function final_MP = applyDetrending(fullTable, raw_MP, method_vec)
    % 위성별, Arc별, 방식별 바이어스 제거 (Mean Removal) + Cycle Slip 감지 추가
    
    final_MP = nan(size(raw_MP));
    
    % PRN 처리 (categorical 대응)
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
        
        vals = raw_MP(idx);
        times = fullTable.Time(idx);
        methods = method_vec(idx);
        
        % LLI 데이터 가져오기 (없으면 0 처리)
        if ismember('LLI1', fullTable.Properties.VariableNames)
            lli1 = fullTable.LLI1(idx);
        else
            lli1 = zeros(size(idx));
        end
        
        % 유효한 데이터만 처리
        validMask = ~isnan(vals);
        if sum(validMask) == 0, continue; end
        
        % ---------------------------------------------------
        % [Arc Segmentation 강화]
        % ---------------------------------------------------
        
        % 1. 시간 갭 (60초)
        t_sub = times(validMask);
        isTimeGap = [0; diff(t_sub)] > 60;
        
        % 2. 방식 변경 (Dual <-> Triple)
        m_sub = methods(validMask);
        isMethodChange = [0; diff(m_sub) ~= 0];
        
        % 3. LLI (Cycle Slip) 플래그 확인 [추가됨]
        % LLI의 1번 비트가 1이면 슬립 발생
        l_sub = lli1(validMask);
        % bitand는 정수형 입력이 필요하므로 변환
        isLLI = bitand(uint8(l_sub), 1) == 1; 
        
        % 4. Geometry-Free (GF) 점프 감지 [고급 옵션]
        % LLI가 놓친 슬립을 잡기 위해 값 자체가 급격히 튀는지 확인
        % (임계값: 5미터 이상 튀면 슬립으로 간주)
        v_sub = vals(validMask);
        jumpDiff = [0; abs(diff(v_sub))];
        isJump = jumpDiff > 5.0; 
        
        % 모든 조건을 합쳐서 Arc ID 생성
        % (셋 중 하나라도 발생하면 자른다)
        isNewArc = isTimeGap | isMethodChange | isLLI | isJump;
        arcIDs = cumsum(isNewArc) + 1;
        
        % ---------------------------------------------------
        % [Mean Removal]
        % ---------------------------------------------------
        detrended_sub = v_sub;
        uArcs = unique(arcIDs);
        
        for k = 1:length(uArcs)
            currArc = (arcIDs == uArcs(k));
            
            % 구간의 길이가 너무 짧으면(예: 1~2 epoch) 평균이 부정확하므로 제외 가능
            % if sum(currArc) < 5, detrended_sub(currArc) = NaN; continue; end
            
            % 해당 구간 평균 제거
            bias = mean(v_sub(currArc), 'omitnan');
            detrended_sub(currArc) = v_sub(currArc) - bias;
        end
        
        % 결과 저장
        saveIdx = idx(validMask);
        final_MP(saveIdx) = detrended_sub;
    end
end

function [group] = groupingTable(dataTable, varName)
    parameterGroup=findgroups(dataTable.(varName));
    num_group=length(unique(parameterGroup));
    group = cell(num_group,1);
    for i=1:num_group
        group{i}=dataTable(parameterGroup==i,:);
    end
end