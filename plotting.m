%%% Plotting (Dual Frequency 결과만 출력)
figure
title("Dual Frequency MP (L1-L2)")
hold on
grid on

% PRN 개수만큼 반복 (MPcell 크기에 맞춰 자동 조절)
numSats = length(MPcell);
for i = 1:numSats
    if isempty(MPcell{i}), continue; end
    
    % MP_Method == 2 (Dual)인 것만 플롯 (사실상 전 데이터가 2번임)
    valid_idx = MPcell{i}.MP_Method == 2;
    if any(valid_idx)
        plot(MPcell{i}.Time(valid_idx), MPcell{i}.MP_Hybrid(valid_idx), 'o', 'MarkerSize', 2);
    end
end
xlabel('Time'); ylabel('MP (m)');

