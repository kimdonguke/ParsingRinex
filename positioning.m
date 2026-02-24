clear; clc;
addpath(genpath('functions')); % positioning.m 파일이 있는 폴더 경로

parsed_dir = fullfile('data', 'PARSED_MAT');
save_dir = fullfile('data', 'MP_RESULTS');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

fprintf('\n========================================\n');
fprintf('[Start] MP 조합 + 위성 고도각 병합 일괄 처리 시작\n');
fprintf('========================================\n');

obs_parsed_files = dir(fullfile(parsed_dir, '*_MO_parsed.mat'));

for i = 1:length(obs_parsed_files)
    filename = obs_parsed_files(i).name;
    filepath = fullfile(parsed_dir, filename);
    
    [~, name_only, ~] = fileparts(filename);
    save_filename = strrep(name_only, '_parsed', '_MP'); 
    save_path = fullfile(save_dir, strcat(save_filename, '.mat'));
    
    % 🔥 수정된 핵심 로직: 파일 존재 여부 + Az/El 컬럼 존재 여부 확인
    need_processing = true; % 기본값: 처리가 필요함
    
    if exist(save_path, 'file')
        % 파일이 있으면 잠시 로드해서 내부 구조 확인 (속도 빠름)
        loaded_data = load(save_path, 'MPcell');
        
        % 데이터가 들어있는 첫 번째 유효한 위성(셀) 찾기
        valid_idx = find(~cellfun(@isempty, loaded_data.MPcell), 1);
        
        if ~isempty(valid_idx)
            % 해당 테이블의 컬럼 이름(VariableNames) 목록 추출
            var_names = loaded_data.MPcell{valid_idx}.Properties.VariableNames;
            
            % 'el'과 'az' 컬럼이 모두 존재하는지 확인
            if ismember('el', var_names) && ismember('az', var_names)
                fprintf('  > [%02d/%02d] [PASS] Az/El 컬럼 확인 완료 (건너뜀): %s\n', i, length(obs_parsed_files), filename);
                need_processing = false; % 처리 불필요, 패스!
            else
                fprintf('  > [%02d/%02d] [UPDATE] 구버전 감지. Az/El 추가를 위해 재처리합니다: %s\n', i, length(obs_parsed_files), filename);
            end
        end
    end
    
    % 처리가 필요 없는 경우 다음 파일로 넘어감
    if ~need_processing
        continue;
    end
    
    fprintf('  > [%02d/%02d] 계산 중: %s ... ', i, length(obs_parsed_files), filename);
    
    fprintf('  > [%02d/%02d] 계산 중: %s ... ', i, length(obs_parsed_files), filename);
    
    % 1. OBS 데이터 로드
    load(filepath, 'return_OBS');
    
   % ==========================================================
    % 2. 짝이 맞는 NAV 데이터 로드 (스마트 매칭 방식 적용)
    % ==========================================================
    % 파일명 예시: YONS00KOR_R_20260151004_01D_30S_MO_parsed.mat
    
    % '_R_' 문자열의 위치를 찾아서 그 뒤의 7글자(연도4 + DOY3)를 추출합니다.
    idx = strfind(filename, '_R_');
    if isempty(idx)
        error('파일명 형식이 예상과 다릅니다: %s', filename);
    end
    
    yyyy_doy = filename(idx+3 : idx+9); % 예: '2026015'
    
    % 추출한 날짜를 BRDC 정규 포맷에 끼워 넣어 완벽한 NAV 파일명을 만듭니다.
    nav_filename = sprintf('BRDC00IGS_R_%s0000_01D_MN_parsed.mat', yyyy_doy);
    nav_filepath = fullfile(parsed_dir, nav_filename); 
    
    % NAV 파일이 존재하는지 안전하게 확인 후 로드
    if ~exist(nav_filepath, 'file')
        fprintf('  > [경고] 짝이 맞는 NAV 파일이 없습니다: %s\n', nav_filename);
        continue; % 에러를 내지 않고 다음 날짜로 넘어감
    end
    
    load(nav_filepath, 'return_NAV');
    % ==========================================================
    
    % 3. MP 조합 계산
    [MPtable, OutlierTable] = calcMultipathComb(return_OBS);
    
    % 4. 위성 좌표 및 고도각(el) 계산 (수정된 함수들 호출)
    [~, ttx_table] = clock_Correction(return_OBS, return_NAV);
    xyz_table = sat_setting2(return_NAV, ttx_table);
    
    % 5. 💎 핵심: MPtable과 xyz_table을 Time과 PRN 기준으로 완벽하게 병합 (Inner Join)
    % 중복되는 컬럼 없이 el, az, x, y, z 가 MPtable 옆에 예쁘게 달라붙습니다.
    MPtable_Joined = innerjoin(MPtable, xyz_table, 'Keys', {'Time', 'PRN'});
    
    % 6. PRN 기준으로 그룹핑
    MPcell = groupingTable(MPtable_Joined, 'PRN');
    
    if ~isempty(OutlierTable)
        OutlierCell = groupingTable(OutlierTable, 'PRN');
        outlier_count = height(OutlierTable);
    else
        OutlierCell = {}; outlier_count = 0;
    end
    
    % 7. 최종 저장
    save(save_path, 'MPcell', 'OutlierCell');
    
    fprintf('완료! (결측치: %d개)\n', outlier_count);
end


%%%% SV Clock Correction %%%%
function [ttx_epoch, ttx_table] = clock_Correction(Obs,Nav)  % ttx_cell{i,1} = time, ttx_cell{i,2}=satellite number
    %% get Obs data %%
    %sObs=table2array(Obs);
    c=299792458;                    % speed of light(meter/sec)
    tolerance = 1e-12;              %convergence threshold
    mu=3.986005e+14;                % earth's gravitational constant meter^3/sec^2
    F=-2*sqrt(mu)/c^2;              % 상대성 효과 보정 힘
    L1 = 1575.42;                   % L1CA 주파수
    L2 = 1227.60;                   % L2 주파수
    mu2 = (L1/L2)^2;                % TGD 계수
    
    Obs = groupingTable(Obs, 'Time');
    
    ttx_list = [];
    ttx_epoch=cell(height(Obs),1);
    
    
    for i=1:height(Obs) % iteration by whole grouped data
        
        obs_i=Obs{i};
        trx = obs_i.Time(1);
        PRN=obs_i.PRN;
        t_k=0;
        index_num=0;
        
        TTX=nan(length(PRN),1);
        P=Obs{i}.P1;     %Pseudorange
        SNR = Obs{i}.SNR1; % 
        % code = Obs{i}.code;
        
        %%% algorithm problem on here
        for j=1:length(PRN) % iteration by each satellite, calculate t_tx
             sNav=Nav{PRN(j)};
             index_num=get_NavIndexbyTime(sNav,trx);
    
             M0=sNav.M0(index_num);
             A=sNav.A(index_num);
             toe=sNav.toe(index_num);
             toc=sNav.toc(index_num);
             deln=sNav.deln(index_num);
             e=sNav.e(index_num); 
             af0=sNav.af0(index_num);
             af1=sNav.af1(index_num);
             af2=sNav.af2(index_num);
             TGD=sNav.tgd(index_num);
             % if obs_i.code == 23
             %     TGD = mu2*TGD;
             % end
             
             ttx_new=obs_i.Time(j)-obs_i.P1(j)/c;   % define ttx_new to first iteration set ttx_new
             t_SV=ttx_new;
             ttx=0;
             % fprintf("ttx_new = ")
             % disp(fff(ttx_new))
             max_iter=1;
            
             while abs(ttx_new-ttx)>tolerance && max_iter<200
                max_iter=max_iter+1;

                ttx=ttx_new;
                t_k=ttx_new-toe; %  time from reference epoch
                M_k=M0+(sqrt(mu/(A^3))+deln)*t_k; % M_k = M_0 + n*t_k /// n=n_0+del n
                E_k=M_k;
                for k=1:10
                    E_kOld=E_k;
                    E_k=M_k+e*sin(E_k);
                    if norm(E_k-E_kOld)<tolerance
                        break
                    end
                end
                delt_r=F*sqrt(A)*e*sin(E_k);
                delt_SV=af0+af1*(ttx-toc)+af2*(ttx-toc)^2+delt_r-TGD;
                ttx_new=t_SV-delt_SV;
             end
             P(j)=P(j)+c*delt_SV;
             TTX(j)=ttx_new;
        end

        ttx_epoch{i}=table(PRN,TTX,P,SNR);
        ttx_list=[ttx_list; [obs_i.Time, PRN, TTX, P, SNR]];
    end
    Time = ttx_list(:,1); % Time 추가!
    PRN  = ttx_list(:,2);
    TTX  = ttx_list(:,3);
    P    = ttx_list(:,4);
    SNR  = ttx_list(:,5);
    ttx_table = table(Time, PRN, TTX, P, SNR); % Time 포함해서 테이블 생성
end

%%%%% sattelite Positioning function %%%%%
function [x y z] = sat_Positioning(Nav, TTX)
    c=299792458;                    % speed of light(meter/sec)
    tolerance = 1e-12;              %convergence threshold
    mu=3.986005e+14;                % earth's gravitational constant meter^3/sec^2
    OMGd_earth=7.2921151467e-5;     % earth rotation rate (rad/sec)
    
    ttx = TTX.TTX;
    P = TTX.P;
    sat_numbering = TTX.PRN;
    
    del_t = P./c;
    
    fixed_time = ttx; %fixed time;
    
    % variable %
    [sat, IODE, IODC, sva, svh, week, toes, ttrs, af0, af1, af2, A, e, i0,...
    OMG0, omg, M0, deln, OMGd, idot, crc, crs, cuc, cus, cic, cis, code,...
    flag, fit, tgd, toe, toc, ttr]=define_Sparameter(Nav,length(sat_numbering));
    
    t_k=zeros(length(sat_numbering),1);                   %time start on toe
    
    for i=1:length(sat_numbering)
        sNav=Nav{sat_numbering(i)};
        % if i==1
        %     E = datetime(1970,1,1);
        %     fprintf("fixed time = %f\n",fixed_time(1))
        %     disp(fixed_time(1))
        %     disp(datetime(fixed_time(1),'ConvertFrom','epochtime','Epoch',E));
        % end
        ttx_time=fixed_time(i);

        j=get_NavIndexbyTime(sNav, ttx_time);
        
        toe=sNav.toe(j);
        t_k(i)=ttx_time-toe;
        toes(i)=sNav.toes(j);
        toes(i) = sNav.toes(j);
        A(i) = sNav.A(j);                       % get semi major axis
        e(i) = sNav.e(j);                       % get eccentricity
        i0(i) = sNav.i0(j);                     % get inclination angle
        OMG0(i) = sNav.OMG0(j);                 % get longitude of ascending node
        omg(i) = sNav.omg(j);                   % get argument of perigee
        M0(i) = sNav.M0(j);                     % get mean anomaly
        deln(i) = sNav.deln(j);                 % get mean motion difference
        OMGd(i) = sNav.OMGd(j);                 % get rate of right ascension
        idot(i) = sNav.idot(j);                 % get rate of inclination angle
        crc(i) = sNav.crc(j);                   % get harmonic cos radius
        crs(i) = sNav.crs(j);                   % get harmonic sine radius
        cuc(i) = sNav.cuc(j);                   % get harmonic cos latitude
        cus(i) = sNav.cus(j);                   % get harmonic sine latitude
        cic(i) = sNav.cic(j);                   % get harmonic cos inclination
        cis(i) = sNav.cis(j);                   % get harmonic sine inclination
    end
    
    %%%%%   Computation of GPS Coordinate   %%%%%%%
    M_k=M0+(deln+sqrt(mu./A.^3)).*t_k;
    E_k = M_k;  % Initialize eccentric anomaly
    for k = 1:10  % Iterate to solve Kepler's equation
        E_kOld=E_k;
        E_k = M_k + e .* sin(E_k);
        if norm(E_k-E_kOld) < tolerance
            break
        end
    end
    nu_k=2.*atan2(sqrt(1.+e).*sin(E_k./2),sqrt(1.-e).*cos(E_k./2));
    ffi=nu_k+omg;
    u_k=ffi+cus.*sin(2.*(ffi))+cuc.*cos(2.*(ffi));
    r_k=A.*(1.-e.*cos(E_k))+crs.*sin(2.*(ffi))+crc.*cos(2.*(ffi));
    i_k=i0+idot.*t_k+cis.*sin(2.*(ffi))+cic.*cos(2.*(ffi));
    OMG_k=OMG0+(OMGd-OMGd_earth).*t_k-OMGd_earth.*toes;
    
    % Convert coordinates
    x_ko = r_k .* cos(u_k);
    y_ko = r_k .* sin(u_k);
    z = y_ko .* sin(i_k);
    
    x=x_ko.*cos(OMG_k)-y_ko.*cos(i_k).*sin(OMG_k);
    y=x_ko.*sin(OMG_k)+y_ko.*cos(i_k).*cos(OMG_k);
end


function [returnSatList] = sat_setting(Nav, ttx)

    Re = 6378137;      % 지구 반지름
    hm = 350000;       % 전리층 고도 (350km)
   
        
    sat_xyzlist=cell(height(ttx),1);
    
    for i=1:height(ttx)
        
        [x, y, z]=sat_Positioning(Nav,ttx{i}); %ttx_list, %grouped sat Obs Data
        PRN=ttx{i}.PRN;
        
        [az, el] = getAzEl(x,y,z);
        sat_xyzlist{i} = table(PRN,x,y,z,az,el);
    end
    returnSatList=sat_xyzlist;
end

% === 수정할 코드 ===
function returnSatlist = sat_setting2(Nav, ttx)
    Time = ttx.Time; % Time 가져오기!
    PRN = ttx.PRN;
    [x, y, z] = sat_Positioning(Nav,ttx);
    [az, el] = getAzEl(x,y,z);
    returnSatlist = table(Time, PRN, x, y, z, az, el); % Time 포함해서 반환
end


function [az, el] = getAzEl(x, y, z)
    ref_position=[-3047507.380, 4043980.305, 3865242.828];
    ref_position=ecef2lla(ref_position);
    ref_positionRAD=deg2rad(ref_position);
    
    [enu]=lla2enu(ecef2lla([x y z]),ref_position,'ellipsoid');
    e=enu(:,1);
    n=enu(:,2);
    u=enu(:,3);
    az=atan2(e,n);
    el = atan2(u,sqrt(e.^2+n.^2));     %calculate azimuth, elevation angle
end