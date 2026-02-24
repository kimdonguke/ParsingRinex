
% ==========================================
% 1. 경로 설정 및 저장 폴더 준비
% ==========================================
addpath(genpath('functions'));

% 파싱된 데이터가 있는 폴더
parsed_dir = fullfile('data', 'PARSED_MAT');

% MP 조합 결과를 저장할 새 폴더 생성
save_dir = fullfile('data', 'MP_RESULTS');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    fprintf('  > [알림] MP 결과 저장용 폴더 생성 완료: %s\n', save_dir);
end

% ==========================================
% 2. MP 조합 일괄 계산 및 그룹핑 시작
% ==========================================
fprintf('\n========================================\n');
fprintf('[Start] 14일 치 MP 조합(Multipath Comb.) 일괄 처리 시작\n');
fprintf('========================================\n');

% '_MO_parsed.mat'로 끝나는 OBS 파싱 파일 목록 불러오기
obs_parsed_files = dir(fullfile(parsed_dir, '*_MO_parsed.mat'));

for i = 1:length(obs_parsed_files)
    filename = obs_parsed_files(i).name;
    filepath = fullfile(parsed_dir, filename);
    
    % 저장될 파일명 미리 계산 (예: YONS00KOR..._MO_MP.mat)
    [~, name_only, ~] = fileparts(filename);
    save_filename = strrep(name_only, '_parsed', '_MP'); 
    save_path = fullfile(save_dir, strcat(save_filename, '.mat'));
    
    % 🔥 파일이 이미 존재하면 건너뛰기
    if exist(save_path, 'file')
        fprintf('  > [%02d/%02d] [PASS] 이미 처리됨: %s\n', i, length(obs_parsed_files), filename);
        continue;
    end
    
    fprintf('  > [%02d/%02d] MP 계산 중: %s ... ', i, length(obs_parsed_files), filename);
    
    % 1) 파싱된 .mat 파일 로드 (return_OBS 변수가 메모리에 올라옴)
    % 주의: rinexread보다 수십 배 빠릅니다!
    load(filepath, 'return_OBS');
    
    % 2) 사용자 정의 MP 조합 함수 실행
    [MPtable, OutlierTable] = calcMultipathComb(return_OBS);
    
    % 3) PRN 기준으로 그룹핑 (정상 데이터)
    MPcell = groupingTable(MPtable, 'PRN');
    
    % 4) 결측치(Outlier) 처리 및 그룹핑
    if ~isempty(OutlierTable)
        OutlierCell = groupingTable(OutlierTable, 'PRN');
        outlier_count = height(OutlierTable);
    else
        OutlierCell = {}; % 빈 셀로 초기화
        outlier_count = 0;
    end
    
    % 5) 결과 저장 (.mat 파일)
    % MPcell, OutlierCell 두 개의 변수를 저장합니다.
    save(save_path, 'MPcell', 'OutlierCell');
    
    % 처리 완료 로그 출력
    fprintf('완료! (결측치: %d개)\n', outlier_count);
end

fprintf('\n========================================\n');
fprintf('[End] 모든 MP 조합 계산 및 저장 완료!\n');
fprintf('저장 위치: %s\n', save_dir);
fprintf('========================================\n');



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