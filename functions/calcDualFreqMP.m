%%% Functions %%%

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