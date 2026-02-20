%%% supervised train Loss Function %%%
function [theta_day, theta_night, history] = trainGnssCovariance(sat_xyzlist, ttx, ref_pos_ecef)
    tolerance = 1e-12;              %convergence threshold
    % 1. 하이퍼파라미터 및 초기값 설정
    total_epochs = size(sat_xyzlist, 1);
    num_train = round(total_epochs * 0.8); % 80% 훈련 데이터
    
    theta_day = [5, 5, 5];   
    theta_night = [5, 5, 5]; 
    
    lr = 0.0001;      % 학습률 (Learning Rate)
    max_iter = 100;    % 전체 데이터 반복 횟수 (Epochs)
    
    % 데이터 셔플링을 위한 인덱스 생성
    train_idx = 1:num_train;
    history = struct('loss', []);

    fprintf('최적화 시작 (훈련 데이터: %d epochs)...\n', num_train);
    ref_positionlla=ecef2lla(ref_pos_ecef);
    ref_positionRAD=deg2rad(ref_positionlla);

    for iter = 1:max_iter
        % 에포크마다 순서 섞기 (Stochastic Gradient Descent 효과)
        shuffled_idx = train_idx(randperm(num_train));
        epoch_loss = 0;
        
        grad_day_total = [0, 0, 0];
        grad_night_total = [0, 0, 0];

        for t = shuffled_idx
            % 현재 에포크 데이터 추출
            H_mat = []; % 관측 행렬
            rho = ttx{t, 1}.P;
            el_list = sat_xyzlist{t, 3}(:, 2);
            isDay_list = sat_xyzlist{t, 4}(:, 4);
            num_sat = length(rho);
            s_sat_xyzlist=sat_xyzlist(t,:);
            s_ttx=ttx(t,:);
            
            % --- Forward Pass: WLS 수행 ---
            W_diag = zeros(num_sat, 1);
            for i = 1:num_sat
                el = el_list(i);
                if isDay_list(i) == 1
                    sig2 = theta_day(1)^2 + (theta_day(2)^2 / sin(el)^theta_day(3));
                else
                    sig2 = theta_night(1)^2 + (theta_night(2)^2 / sin(el)^theta_night(3));
                end
                W_diag(i) = 1 / sig2;
            end
            W = diag(W_diag);
            
            % WLS 위치 계산 (단순화를 위해 1회 루프로 가정하거나 함수 호출)
            % H_mat 구성을 위해 위성-사용자 벡터 계산
            % 여기서는 기존 코드의 H 행렬 구성 로직을 간략화하여 적용
            [Ion_delay, Trop_delay]=get_IonTropDelay(s_ttx,s_sat_xyzlist,ref_positionRAD,1);
            [ux, uy, uz, cb, H, resid] = solve_wls_internal(sat_xyzlist{t, 1}, rho-Ion_delay-Trop_delay, W);
            x_est = [ux; uy; uz];
            
            % Loss 계산
            err_vec = x_est - ref_pos_ecef';
            epoch_loss = epoch_loss + norm(err_vec)^2;

            % --- Backward Pass: Gradient 계산 ---
            % inv 대신 직접 역행렬 계산 (더 효율적)
            M = (H' * W * H) \ eye(4);
            err_vec_full=[err_vec;0];
            G_common = err_vec_full' * M * H'; % 1 x N 행렬
            
            % 벡터화된 그라디언트 계산 (for 루프 제거)
            Gi_vec = G_common(:) .* resid; % 벡터화된 공통 기여도
            
            % 낮/밤 케이스별로 분리하여 벡터 연산
            isDay_logical = (isDay_list == 1);
            sin_el = sin(el_list);
            
            % 낮 케이스 그라디언트 계산
            day_indices = find(isDay_logical);
            if ~isempty(day_indices)
                w_i_day = W_diag(day_indices);
                el_day = el_list(day_indices);
                sin_el_day = sin_el(day_indices);
                Gi_day = Gi_vec(day_indices);
                
                dw_da_day = -w_i_day.^2 * 2 * theta_day(1);
                dw_db_day = -w_i_day.^2 .* (2 * theta_day(2) ./ sin_el_day.^theta_day(3));
                dw_dn_day = w_i_day.^2 .* (theta_day(2)^2 ./ sin_el_day.^theta_day(3)) .* log(sin_el_day);
                
                % 각 그라디언트 요소별로 sum 계산
                grad_day_update = [sum(Gi_day .* dw_da_day), ...
                                  sum(Gi_day .* dw_db_day), ...
                                  sum(Gi_day .* dw_dn_day)];
                grad_day_total = grad_day_total + grad_day_update;
            end
            
            % 밤 케이스 그라디언트 계산
            night_indices = find(~isDay_logical);
            if ~isempty(night_indices)
                w_i_night = W_diag(night_indices);
                el_night = el_list(night_indices);
                sin_el_night = sin_el(night_indices);
                Gi_night = Gi_vec(night_indices);
                
                dw_da_night = -w_i_night.^2 * 2 * theta_night(1);
                dw_db_night = -w_i_night.^2 .* (2 * theta_night(2) ./ sin_el_night.^theta_night(3));
                dw_dn_night = w_i_night.^2 .* (theta_night(2)^2 ./ sin_el_night.^theta_night(3)) .* log(sin_el_night);
                
                % 각 그라디언트 요소별로 sum 계산
                grad_night_update = [sum(Gi_night .* dw_da_night), ...
                                    sum(Gi_night .* dw_db_night), ...
                                    sum(Gi_night .* dw_dn_night)];
                grad_night_total = grad_night_total + grad_night_update;
            end
        end
        
        % 파라미터 업데이트 (Batch 평균 적용)
        theta_day = theta_day - lr * (grad_day_total / num_train);
        theta_night = theta_night - lr * (grad_night_total / num_train);
        
        % 물리적 제약 조건 적용 (양수 유지)
        theta_day(1:2) = max(theta_day(1:2), 0.01); theta_day(3) = max(theta_day(3), 1.0);
        theta_night(1:2) = max(theta_night(1:2), 0.01); theta_night(3) = max(theta_night(3), 1.0);
        
        history.loss(iter) = epoch_loss / num_train;
        if iter~=1 && abs(history.loss(iter)-history.loss(iter-1))<tolerance
            fprintf("==================================\n Loss converged in %d iteration\n=============================================\n",iter);
            break
        end
        fprintf('Iter %d: Avg Loss = %.4f | Day_n = %.2f | Night_n = %.2f\n', ...
            iter, history.loss(iter), theta_day(3), theta_night(3));
    end
end