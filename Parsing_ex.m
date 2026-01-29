clear;
clc;

addpath("C:\Users\김동욱\OneDrive\Desktop\GNSSLAB\2026-동계방학\Study_1\ParsingRinex\data");
addpath(genpath('data'));
load('Nav.mat');
load('Obs_1hour.mat')

[nav_data, nav_header] = rinexread("BRDC00IGS_R_20253590000_01D_MN.rnx");
[obs_data, obs_header] = rinexread("YONS00KOR_R_20253590000_01D_30S_MO.rnx");

retrun_NAV = parsingNavigationBody(nav_data);

%% Navigation Parsing %%
%======================%
%Indexing by sorted PRN(satellite ID)
%Input : Rinex file
%Output : (PRN * 1) cell
%======================%


function [return_NAV] = parsingNavigationBody(nav_data)
    navTable = nav_data.GPS;
    navTable=timetable2table(navTable);
    navTable = sortrows(navTable,1);
    navTable = convertToRtkLibFormat(navTable);
    return_NAV = groupingTable(navTable,'PRN');
end
function [return_NAV] = parsingNavigationHeader(nav_header)

end

function [group] = groupingTable(dataTable, varName)
    parameterGroup=findgroups(dataTable.(varName));
    num_group=length(unique(parameterGroup));
    group = cell(num_group,1);
    for i=1:num_group
        group{i}=dataTable(parameterGroup==i,:);
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