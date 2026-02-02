%Obs 받고, T_sv는 수렴시켜야되자너 위성에서 나오는 tsv 세 개 다 수렴 시키고 그러고 여기로 넘겨서 각 위성마다 MP combination
%적용 시킨다?
%적용시키고 나오면 계단식(lambda integer ambiguty + reciever noise, multipath error + residual error만 남는데)
%이때 수신기 노이즈만 보고 싶으니까 integer ambiguty 구간 적분으로 없애고, multipath도 미소 시간에 대해선
%바이어슨데 예가 integer ambiguty랑 구분이 안 될거같은데

MP_triple = calc(return_OBS,'Code');

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
    b1 = 0.142;
    b2 = -0.767;
    b3 = 0.626;
    
    % 데이터 추출 (Code Measurement)
    P1 = subTable.P1;
    P2 = subTable.P2;
    P5 = subTable.P5;
    
    % 선형 결합 수행
    mp_values = b1 * P1 + b2 * P2 + b3 * P5;
end