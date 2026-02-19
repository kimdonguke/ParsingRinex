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
