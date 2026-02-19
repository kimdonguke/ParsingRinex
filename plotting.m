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

numSats = length(MPcell);
for i = 1:numSats
    figure
    title("Dual Frequency MP (L1-L2)")

    if isempty(MPcell{i}), continue; end
    
    valid_idx = MPcell{i}.MP_Method == 2;
    if any(valid_idx)
        plot(MPcell{i}.SNR1(valid_idx), MPcell{i}.MP_Hybrid(valid_idx), 'o', 'MarkerSize', 2);
    end
    grid on
end
xlabel('SNR'); ylabel('MP (m)');


%whole SNR
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
xlabel('SNR'); ylabel('MP (m)');



% test
groupingSNR = groupingTable(MPtable,'SNR1');
ex_cell = []
for i = 1:height(groupingSNR)
    % NaN이 아닌 값들만 남기기
    delNAN = groupingSNR{i}.MP_Hybrid(~isnan(groupingSNR{i}.MP_Hybrid));
    ex_cell = [ex_cell; [groupingSNR{i}.SNR1(1), std2(delNAN)]];
end
plot(ex_cell(:,1),ex_cell(:,2), 'ro', 'MarkerSize', 3, 'MarkerFaceColor','red')
