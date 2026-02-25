fprintf("F5 또 잘못 눌렀어용~/n");

clc;

function a = aaaa()
    % =========================================================================
    % [사전 준비] F9 실행 전, 혹은 aaaa() 함수 내에서 가장 먼저 실행할 블록
    % MPcell만 로드된 상태에서 통계 분석용 임시 MPtable을 만듭니다.
    % =========================================================================
    valid_cells = MPcell(~cellfun(@isempty, MPcell)); % 빈 셀 제거
    MPtable = vertcat(valid_cells{:});                % 셀을 하나의 테이블로 세로 병합
    
    numSats = length(MPcell);
    
    % =========================================================================
    % 1. MP hybrid 출력 (전체 위성)
    % =========================================================================
    figure
    title("Dual Frequency MP (L1-L2) - All Satellites")
    hold on; grid on;
    for i = 1:numSats
        if isempty(MPcell{i}), continue; end
        valid_idx = MPcell{i}.MP_Method == 2;
        if any(valid_idx)
            plot(MPcell{i}.Time(valid_idx), MPcell{i}.MP_Hybrid(valid_idx), 'o', 'MarkerSize', 2);
        end
    end
    xlabel('Time'); ylabel('MP (m)');
    
    % 이상치 데이터 출력 예시
    % 1. 데이터 추출
    
   
    
    x_raw = rad2deg(MPtable.el);
    y_raw = MPtable.SNR1;
    z_raw = abs(MPtable.MP_Hybrid);
    
    % 2. 데이터 격자화 (Scattered 데이터일 경우)
    % 만약 x, y가 이미 규칙적인 격자의 벡터라면 바로 reshape을 쓰셔도 되지만,
    % 일반적인 테이블 데이터라면 아래처럼 그리드를 생성하는 것이 정확합니다.
    [X, Y] = meshgrid(unique(x_raw), unique(y_raw));
    Z = griddata(x_raw, y_raw, z_raw, X, Y);
    
    % 3. 시각화
    figure('Color', 'w');
    surf(X, Y, Z);
    
    % 스타일링
    shading interp;      % 면을 부드럽게 보정
    colorbar;            % 색상 바 추가
    colormap('jet');     % 색상 맵 설정 (데이터 변화가 잘 보임)
    
    % 라벨 설정
    xlabel('Elevation Angle (deg)');
    ylabel('SNR1');
    zlabel('MP Hybrid');
    title('Multi-path Hybrid Analysis');
    
    % 보기 편한 각도로 회전
    view(45, 30);
    grid on;
    
    
    
    % =========================================================================
    % 2. MP raw 출력 (전체 위성, bias 포함)
    % =========================================================================
    figure
    title("Dual Frequency MP raw(bias O) - All Satellites")
    hold on; grid on;
    for i = 1:numSats
        if isempty(MPcell{i}), continue; end
        valid_idx = MPcell{i}.MP_Method == 2;
        if any(valid_idx)
            plot(MPcell{i}.Time(valid_idx), MPcell{i}.MP_raw(valid_idx), 'o', 'MarkerSize', 2);
        end
    end
    xlabel('Time'); ylabel('MP (m)');
    
    % =========================================================================
    % 3. MP hybrid 출력 (위성별 개별 Figure 생성)
    % =========================================================================
    for i = 1:numSats
        if isempty(MPcell{i}), continue; end
        figure
        title(sprintf("Dual Frequency MP (L1-L2) - PRN %d", i))
        hold on; grid on;
        valid_idx = MPcell{i}.MP_Method == 2;
        if any(valid_idx)
            plot(MPcell{i}.Time(valid_idx), MPcell{i}.MP_Hybrid(valid_idx), 'o', 'MarkerSize', 2);
        end
        ylim([-2.5 2.5])
        xlabel('Time'); ylabel('MP (m)');
    end
    
    % =========================================================================
    % 4. 통계: Elevation - SNR
    % (참고: 기존 xyz_table 대신 통합된 MPtable의 고도각 컬럼 'el' 사용 가정)
    % 만약 컬럼명이 'Elevation'이라면 MPtable.Elevation 으로 변경하세요.
    % =========================================================================
    % 1. 위성 번호(PRN)별로 나누어 그리기 (MPtable에 PRN 열이 있다고 가정)
    figure; hold on; grid on;
    unique_prns = unique(MPtable.PRN); % PRN 열 이름을 확인하세요
    
    for i = 1:length(unique_prns)
        idx = (MPtable.PRN == unique_prns(i));
        plot(rad2deg(MPtable.el(idx)), MPtable.SNR1(idx), '.', 'MarkerSize', 5);
    end
    
    xlabel('Elevation [degree]'); ylabel('SNR');
    title('Elevation vs SNR (Color by PRN)');
    figure
    plot3(rad2deg(MPtable.el), MPtable.SNR1, abs(MPtable.MP_Hybrid),'o','MarkerSize',2);
    grid on
    xlabel("el")
    ylabel("SNR")
    zlabel("MP")
    
    % =========================================================================
    % 5. 통계: Elevation - MP-combination
    % =========================================================================
    figure
    title("Elevation vs MP-combination")
    hold on; grid on;
    for i = 1:numSats
        if isempty(MPcell{i}), continue; end
        valid_idx = MPcell{i}.MP_Method == 2;
        if any(valid_idx)
            plot(rad2deg(MPcell{i}.el(valid_idx)), MPcell{i}.MP_Hybrid(valid_idx), 'o', 'MarkerSize', 2);
        end
    end
    ylim([-2.5 2.5]);
    xlabel('Elevation [degree]'); ylabel('MP (m)');
    

    
    % =========================================================================
    % 6. 통계: SNR - MP combination (전체 위성 산점도 + 표준편차 오버레이)
    % =========================================================================
    figure
    title("SNR vs MP (L1-L2) with Standard Deviation")
    hold on; grid on;
    
    % 6-1. 산점도 출력
    for i = 1:numSats
        if isempty(MPcell{i}), continue; end
        valid_idx = MPcell{i}.MP_Method == 2;
        if any(valid_idx)
            plot(MPcell{i}.SNR1(valid_idx), MPcell{i}.MP_Hybrid(valid_idx), 'o', 'MarkerSize', 2);
        end
    end
    ylim([-2.5 2.5]);
    xlabel('SNR [dB-Hz]'); ylabel('MP (m)');
    
    % 6-2. 표준편차(Std) 오버레이 출력
    groupingSNR = groupingTable(MPtable, 'SNR1');
    ex_cell = [];
    for i = 1:height(groupingSNR)
        % NaN이 아닌 값들만 남기기
        valid_MP = groupingSNR{i}.MP_Hybrid(~isnan(groupingSNR{i}.MP_Hybrid));
        if ~isempty(valid_MP)
            % std2 대신 기본 내장 함수인 std 사용 (1차원 배열)
            ex_cell = [ex_cell; [groupingSNR{i}.SNR1(1), std(valid_MP)]];
        end
    end
    if ~isempty(ex_cell)
        plot(ex_cell(:,1), ex_cell(:,2), 'ro', 'MarkerSize', 3, 'MarkerFaceColor','red')
    end
    
    % =========================================================================
    % 7. SNR - MP-combination (위성별 개별 Figure 생성)
    % =========================================================================
    for i = 1:numSats
        if isempty(MPcell{i}), continue; end
        figure
        title(sprintf("SNR vs MP-combination - PRN %d", i))
        hold on; grid on;
        valid_idx = MPcell{i}.MP_Method == 2;
        if any(valid_idx)
            plot(MPcell{i}.SNR1(valid_idx), MPcell{i}.MP_Hybrid(valid_idx), 'o', 'MarkerSize', 2);
        end
        ylim([-2.5 2.5]);
        xlabel('SNR [dB-Hz]'); ylabel('MP (m)');
    end
    
    % =========================================================================
    % 8. SNR 구간별 데이터 개수 확인 (SNR - number of data)
    % =========================================================================
    ex_table = table(MPtable.SNR1, MPtable.MP_Hybrid, 'VariableNames', {'SNR1', 'MP_Hybrid'});
    ex_grouped = groupingTable(ex_table, 'SNR1');
    
    figure
    title("SNR vs Number of Data Points")
    hold on; grid on;
    for i = 1:height(ex_grouped)
        % NaN이 아닌 유효한 MP 데이터의 개수 카운트
        valid_num = sum(~isnan(ex_grouped{i}.MP_Hybrid));
        plot(ex_grouped{i}.SNR1(1), valid_num, 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'red');
    end
    xlabel("SNR [dB-Hz]")
    ylabel("Number of Observations")

end