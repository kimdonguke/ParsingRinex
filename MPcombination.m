
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



