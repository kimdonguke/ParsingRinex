fprintf("F5 또 잘못 눌렀어용~/n");

function a=aaaa()
    %%% Plotting (Dual Frequency 결과만 출력)
    figure
    title("Dual Frequency MP (L1-L2)")
    hold on
    grid on
    
    % 전체 결과 출력
    numSats = length(MPcell);
    for i = 1:numSats
        if isempty(MPcell{i}), continue; end
        
        valid_idx = MPcell{i}.MP_Method == 2;
        if any(valid_idx)
            plot(MPcell{i}.Time(valid_idx), MPcell{i}.MP_Hybrid(valid_idx), 'o', 'MarkerSize', 2);
        end
    end
    xlabel('Time'); ylabel('MP (m)');
    
    
    
    
    % 위성별 결과 출력
    numSats = length(MPcell);
    for i = 1:numSats
        figure
        title("Dual Frequency MP (L1-L2)")
        grid on
        if isempty(MPcell{i}), continue; end
        
        valid_idx = MPcell{i}.MP_Method == 2;
        if any(valid_idx)
            plot(MPcell{i}.Time(valid_idx), MPcell{i}.MP_Hybrid(valid_idx), 'o', 'MarkerSize', 2);
        end
        xlim([1766620800 1766707140])
    end
    xlabel('Time'); ylabel('MP (m)');
    
    
    % 통계 분석(SNR, Elevation angle, MP comb)
    % 1. Elevation - SNR
    % 2. Elevation - MP-combination
    % 3. SNR - MP-combination
    
    % 1. Elevation - SNR
    figure
    hold on
    grid on
    plot(rad2deg(xyz_table.el), ttx_table.SNR, 'o', 'MarkerSize', 2)
    title("Elevation - SNR")
    xlabel('Elevation[degree]'); ylabel('SNR');
    
    % 2. Elevation - MP-combination
    figure
    title("Elevation - MP-combination")
    hold on
    grid on
    groupingSatPos = groupingTable(xyz_table, 'PRN');
    numSats = length(MPcell);
    for i = 1:numSats
        if isempty(MPcell{i}), continue; end
        
        valid_idx = MPcell{i}.MP_Method == 2;
        if any(valid_idx)
            plot(rad2deg(groupingSatPos{i}.el(valid_idx)), MPcell{i}.MP_Hybrid(valid_idx), 'o', 'MarkerSize', 2);
        end
    end
    ylim([-2.5 2.5]);
    xlabel('Elevation(degree)'); ylabel('MP (m)');
    
    
    % 3. SNR - MP-combination 개별 보기
    numSats = length(MPcell);
    for i = 1:numSats
        figure
        title("SNR - MP-combination")
    
        if isempty(MPcell{i}), continue; end
        
        valid_idx = MPcell{i}.MP_Method == 2;
        if any(valid_idx)
            plot(MPcell{i}.SNR1(valid_idx), MPcell{i}.MP_Hybrid(valid_idx), 'o', 'MarkerSize', 2);
        end
        grid on
    end
    xlabel('SNR'); ylabel('MP (m)');
    
    
    % 3-2. SNR - MP combination 결과 전체 보기
    figure
    hold on
    grid on
    
    numSats = length(MPcell);
    for i = 1:numSats
        title("Dual Frequency MP (L1-L2)")
    
        if isempty(MPcell{i}), continue; end
        
        valid_idx = MPcell{i}.MP_Method == 2;
        if any(valid_idx)
            plot(MPcell{i}.SNR1(valid_idx), MPcell{i}.MP_Hybrid(valid_idx), 'o', 'MarkerSize', 2);
        end
    end
    ylim([-2.5 2.5]);
    xlabel('SNR'); ylabel('MP (m)');
    
    % 표준편차 확인 위에 코드랑 같이 실행
    groupingSNR = groupingTable(MPtable,'SNR1');
    ex_cell = []
    for i = 1:height(groupingSNR)
        % NaN이 아닌 값들만 남기기
        delNAN = groupingSNR{i}.MP_Hybrid(~isnan(groupingSNR{i}.MP_Hybrid));
        ex_cell = [ex_cell; [groupingSNR{i}.SNR1(1), std2(delNAN)]];
    end
    plot(ex_cell(:,1),ex_cell(:,2), 'ro', 'MarkerSize', 3, 'MarkerFaceColor','red')

    
    
    
end