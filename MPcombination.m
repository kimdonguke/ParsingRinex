%Obs 받고, T_sv는 수렴시켜야되자너 위성에서 나오는 tsv 세 개 다 수렴 시키고 그러고 여기로 넘겨서 각 위성마다 MP combination
%적용 시킨다?
%적용시키고 나오면 계단식(lambda integer ambiguty + reciever noise, multipath error + residual error만 남는데)
%이때 수신기 노이즈만 보고 싶으니까 integer ambiguty 구간 적분으로 없애고, multipath도 미소 시간에 대해선
%바이어슨데 예가 integer ambiguty랑 구분이 안 될거같은데



%%% main %%%

% 함수 호출 결과로 정상 처리된 MPtable과 L2가 없어서 빠진 OutlierTable을 받습니다.
[MPtable, OutlierTable] = calcMultipathComb(return_OBS);

% 그룹핑 (정상 데이터)
MPcell = groupingTable(MPtable, 'PRN');

% (옵션) L2 결측 데이터도 확인하고 싶다면 그룹핑
if ~isempty(OutlierTable)
    OutlierCell = groupingTable(OutlierTable, 'PRN');
    fprintf('L2 측정치 결측 데이터가 %d개 존재하여 별도 저장되었습니다.\n', height(OutlierTable));
else
    fprintf('L2 측정치 결측 데이터가 없습니다.\n');
end




%%% Functions %%%

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

% Triple Freq 함수는 사용하지 않으므로 삭제하거나 유지해도 무방 (호출되지 않음)
% function mp_values = calcTripleFreqMP(subTable) ... end 



function [group] = groupingTable(dataTable, varName)
    parameterGroup=findgroups(dataTable.(varName));
    num_group=length(unique(parameterGroup));
    group = cell(num_group,1);
    for i=1:num_group
        group{i}=dataTable(parameterGroup==i,:);
    end
end