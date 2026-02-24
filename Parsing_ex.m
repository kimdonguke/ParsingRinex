clear;
clc;

% ==========================================
% 1. 경로 설정 및 저장 폴더 준비
% ==========================================
addpath(genpath('data'));
addpath(genpath('functions'));

obs_dir = fullfile('data', 'OBS');
nav_dir = fullfile('data', 'NAV');
save_dir = fullfile('data', 'PARSED_MAT');

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
    fprintf('  > [알림] 파싱 데이터 저장용 폴더 생성 완료: %s\n', save_dir);
end

% ==========================================
% 2. OBS 데이터 일괄 파싱 및 저장 (중복 건너뛰기)
% ==========================================
fprintf('\n========================================\n');
fprintf('[Start] OBS 데이터 일괄 파싱 시작\n');
fprintf('========================================\n');

obs_files = dir(fullfile(obs_dir, '*_MO.rnx'));

for i = 1:length(obs_files)
    filename = obs_files(i).name;
    filepath = fullfile(obs_dir, filename);
    
    % 저장될 파일 경로 미리 계산
    [~, name_only, ~] = fileparts(filename);
    save_path = fullfile(save_dir, strcat(name_only, '_parsed.mat'));
    
    % 🔥 핵심 로직: 파일이 이미 존재하면 건너뛰기
    if exist(save_path, 'file')
        fprintf('  > [%02d/%02d] [PASS] 이미 파싱됨: %s\n', i, length(obs_files), filename);
        continue; % 아래 코드를 실행하지 않고 다음 i로 넘어감
    end
    
    fprintf('  > [%02d/%02d] 파싱 중: %s ... ', i, length(obs_files), filename);
    
    % 파싱 및 저장 수행
    [obs_data, obs_header] = rinexread(filepath);
    return_OBS = parsingObsBody(obs_data);
    save(save_path, 'return_OBS', 'obs_header');
    
    fprintf('저장 완료!\n');
end

% ==========================================
% 3. NAV 데이터 일괄 파싱 및 저장 (중복 건너뛰기)
% ==========================================
fprintf('\n========================================\n');
fprintf('[Start] NAV 데이터 일괄 파싱 시작\n');
fprintf('========================================\n');

nav_files = dir(fullfile(nav_dir, '*_MN.rnx'));

for i = 1:length(nav_files)
    filename = nav_files(i).name;
    filepath = fullfile(nav_dir, filename);
    
    % 저장될 파일 경로 미리 계산
    [~, name_only, ~] = fileparts(filename);
    save_path = fullfile(save_dir, strcat(name_only, '_parsed.mat'));
    
    % 🔥 핵심 로직: 파일이 이미 존재하면 건너뛰기
    if exist(save_path, 'file')
        fprintf('  > [%02d/%02d] [PASS] 이미 파싱됨: %s\n', i, length(nav_files), filename);
        continue;
    end
    
    fprintf('  > [%02d/%02d] 파싱 중: %s ... ', i, length(nav_files), filename);
    
    % 파싱 및 저장 수행
    [nav_data, nav_header] = rinexread(filepath);
    return_NAV = parsingNavigationBody(nav_data);
    save(save_path, 'return_NAV', 'nav_header');
    
    fprintf('저장 완료!\n');
end

fprintf('\n========================================\n');
fprintf('[End] 모든 데이터 처리 완료!\n');
fprintf('========================================\n');


%% Navigation Parsing %%
%======================%
%Indexing by sorted PRN(satellite ID) and rematch format
%Input : Navigation Body.GPS(timetable)
%Output : (PRN * 1) cell
%======================%
function [return_NAV] = parsingNavigationBody(nav_data)
    navTable = nav_data.GPS;
    navTable=timetable2table(navTable);
    navTable = sortrows(navTable,1);
    navTable = convertToRtkLibFormat(navTable);
    return_NAV = groupingTable_PRN(navTable);
end


function [return_NAV] = parsingNavigationHeader(nav_header)

end


%% Obs Parsing %%
%======================%
%refactor table
%Input : Obs body.GPS(timetable)
%Output : Obs body(table)
%======================%
function flatObsTable = parsingObsBody(rawObs)
    rawObs = timetable2table(rawObs.GPS);
    % parseObsToFlatTable
    % rinexread 결과를 학습 및 MPC 계산용 Flat Table로 변환
    %
    % Output Columns:
    % Time(Unix), PRN, P1, L1, S1, LLI1, P2, L2, S2, LLI2, P5, L5, S5, LLI5
    
    % 데이터 행 개수
    numRows = height(rawObs);
    
    % -------------------------------------------------------
    % 2. 기본 정보 (Time, PRN) 변환
    % -------------------------------------------------------
    
    % Time: datetime -> Unix Timestamp (double)
    % TimeZone이 없으면 UTC로 가정하여 변환 (로컬 시간 간섭 방지)
    dt = rawObs.Time;
    if isempty(dt.TimeZone)
        dt.TimeZone = 'UTC';
    end
    Time = posixtime(dt);
    
    % PRN: SatelliteID (이미 숫자라고 하셨으므로 그대로 사용)
    % 만약 categorical이나 cell이라면 double로 변환
    if isnumeric(rawObs.SatelliteID)
        PRN = double(rawObs.SatelliteID);
    else
        % 혹시 모를 상황 대비 (문자열 '1' -> 숫자 1)
        PRN = double(string(rawObs.SatelliteID));
    end
    
    % -------------------------------------------------------
    % 3. L1 Frequency Data (C1C, L1C, S1C, LLI)
    % -------------------------------------------------------
    P1   = getCol(rawObs, 'C1C');       % Pseudorange
    L1   = getCol(rawObs, 'L1C');       % Carrier Phase
    SNR1   = getCol(rawObs, 'S1C');       % SNR (Signal Strength)
    LLI1 = getCol(rawObs, 'L1C_LLI');   % Lock Loss Indicator
    
    % -------------------------------------------------------
    % 4. L2 Frequency Data (Priority: W > X)
    % -------------------------------------------------------
    % C2W(Z-tracking)가 있으면 사용, 없으면 C2X(L2C) 사용
    if ismember('C2W', rawObs.Properties.VariableNames)
        P2   = getCol(rawObs, 'C2W');
        L2   = getCol(rawObs, 'L2W');
        SNR2   = getCol(rawObs, 'S2W');
        LLI2 = getCol(rawObs, 'L2W_LLI');
    else
        P2   = getCol(rawObs, 'C2X');
        L2   = getCol(rawObs, 'L2X');
        SNR2   = getCol(rawObs, 'S2X');
        LLI2 = getCol(rawObs, 'L2X_LLI');
    end
    
    % -------------------------------------------------------
    % 5. L5 Frequency Data (C5X) - 필요 시 사용
    % -------------------------------------------------------
    P5   = getCol(rawObs, 'C5X');
    L5   = getCol(rawObs, 'L5X');
    SNR5   = getCol(rawObs, 'S5X');
    LLI5 = getCol(rawObs, 'L5X_LLI');
    
    % -------------------------------------------------------
    % 6. 테이블 생성
    % -------------------------------------------------------
    flatObsTable = table(Time, PRN, ...
        P1, L1, SNR1, LLI1, ...
        P2, L2, SNR2, LLI2, ...
        P5, L5, SNR5, LLI5);
    
    % (옵션) NaN 데이터 제거가 필요하면 여기서 수행
    % flatObsTable = rmmissing(flatObsTable, 'DataVariables', {'P1', 'L1'});
end

%% Obs Paring Helper function %%
%======================%
%check colomn class and is column empty and return NAN vector
%Input : table, variable name
%Output : column vector
%======================%
function val = getCol(tbl, varName)
    if ismember(varName, tbl.Properties.VariableNames)
        val = tbl.(varName);
        % 혹시 데이터가 cell array로 되어있다면 double 변환
        if iscell(val)
            val = str2double(val);
        end
    else
        val = nan(height(tbl), 1); % 컬럼 없으면 NaN으로 채움
    end
    
    % rinexread가 간혹 LLI 등을 빈 값으로 읽을 때 NaN 처리
    if ~isfloat(val) 
        val = double(val);
    end
end

function rtkNav = convertToRtkLibFormat(rawNav)
    % convertToRtkLibFormat
    % rinexread 결과를 RTKLIB 스타일의 포맷으로 변환
    % - toes, ttrs : GPS Week Seconds (Second Domain)
    % - toe, toc, ttr : Unix Timestamp (Absolute Time Domain)
    
    % 1. 기본 상수 설정 (GPS Epoch: 1980년 1월 6일)
    gpsEpoch = datetime(1980, 1, 6, 0, 0, 0); 
    
    % ---------------------------------------------------------
    % 2. 데이터 매핑 및 계산
    % ---------------------------------------------------------
    
    % [1~16번 열]
    PRN  = rawNav.SatelliteID;
    IODE = rawNav.IODE;
    IODC = rawNav.IODC;
    sva  = rawNav.SVAccuracy;
    svh  = rawNav.SVHealth;
    week = rawNav.GPSWeek;
    
    % toes (Second Domain): Ephemeris 기준 시각 (주 내 초)
    toes = rawNav.Toe; 
    
    % ttrs (Second Domain): 전송 시각 (주 내 초)
    ttrs = rawNav.TransmissionTime; 
    
    af0  = rawNav.SVClockBias;
    af1  = rawNav.SVClockDrift;
    af2  = rawNav.SVClockDriftRate;
    
    % A (Semi-major axis): sqrtA를 제곱하여 계산
    if ismember('sqrtA', rawNav.Properties.VariableNames)
        A = rawNav.sqrtA .^ 2;
    else
        error('sqrtA 데이터가 없습니다.');
    end
    
    e    = rawNav.Eccentricity;
    i0   = rawNav.i0;
    OMG0 = rawNav.OMEGA0;
    omg  = rawNav.omega;
    
    % [17~32번 열]
    M0   = rawNav.M0;
    deln = rawNav.Delta_n;
    OMGd = rawNav.OMEGA_DOT;
    idot = rawNav.IDOT;
    crc  = rawNav.Crc;
    crs  = rawNav.Crs;
    cuc  = rawNav.Cuc;
    cus  = rawNav.Cus;
    cic  = rawNav.Cic;
    cis  = rawNav.Cis;
    code = rawNav.L2ChannelCodes;
    flag = rawNav.L2PDataFlag;
    fit  = rawNav.FitInterval;
    tgd  = rawNav.TGD;
    
    % toe (Unix Time Domain) 계산
    % GPS Epoch + Week * 7일 + toes(초) -> Unix Time 변환
    % calweeks는 달력 기준 계산이므로 정확한 초 계산을 위해 hours 사용 권장
    % 여기서는 week * 7 * 24 * 3600 방식이 안전함
    
    % datetime 배열 생성
    toe_dt = gpsEpoch + hours(week * 168) + seconds(toes);
    toe    = posixtime(toe_dt); % datetime -> unix timestamp (double)
    
    % toc (Unix Time Domain) 계산
    % rinexread의 Time은 이미 날짜 정보를 포함한 datetime임
    toc_dt = rawNav.Time;
    toc    = posixtime(toc_dt); 
    
    % [33번 열]
    % ttr (Unix Time Domain) 계산
    ttr_dt = gpsEpoch + hours(week * 168) + seconds(ttrs);
    ttr    = posixtime(ttr_dt);
    
    % ---------------------------------------------------------
    % 3. 테이블 생성 (지정된 순서 엄수)
    % ---------------------------------------------------------
    rtkNav = table(PRN, IODE, IODC, sva, svh, week, toes, ttrs, ...
                   af0, af1, af2, A, e, i0, OMG0, omg, ...
                   M0, deln, OMGd, idot, crc, crs, cuc, cus, cic, cis, ...
                   code, flag, fit, tgd, toe, toc, ttr);
               
    disp('변환 완료: RTKLIB 포맷 (Unix Time / GPS Seconds 분리 적용됨)');
end