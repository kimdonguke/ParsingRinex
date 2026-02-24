% ---------------- get_NavIndexbyTime ----------------
% input : (Table, variable Name)
% output : index (double)
% ------------------------------------------------
function [index] = get_NavIndexbyTime(Nav, time)
    % 1. 데이터가 비어있는지 안전 검사
    if isempty(Nav)
        error('해당 위성의 항법 데이터가 존재하지 않습니다.');
    end

    toe_all_raw = Nav.toe(:);
    
    % 2. toe 데이터 타입 통일 (Unix Time 기반)
    if isdatetime(toe_all_raw)
        toe_all = seconds(toe_all_raw);
    else
        toe_all = double(toe_all_raw);
    end
    
    % 3. 핵심 로직: ttr 조건 삭제 & toe와의 절대 시간 차이 계산
    % 후처리(Post-processing)이므로 ttr(송신시간)은 무시하고 toe(기준시간)만 봅니다.
    dt = abs(toe_all - time);
    
    % 4. 시간 차이가 '가장 작은(최소)' 인덱스 추출
    [min_dt, best_idx] = min(dt);
    
    % 5. 유효성 검사 및 반환 (경계 완화)
    % 보통 유효기간은 2시간(7200초)이지만, 파일 경계선(자정)이나 
    % 데이터 누락을 고려해 4시간(14400초)까지는 허용해줍니다.
    if min_dt <= 14400
        index = best_idx;
    else
        % 14일치 루프가 에러로 터지는 것을 막기 위해 error 대신 warning 사용
        warning('PRN %d: 최적의 궤도력을 찾지 못해 가장 가까운 궤도력(시간차 %d초)을 임시로 사용합니다.', Nav.PRN(1), round(min_dt));
        index = best_idx; % 일단 가장 가까운 걸 던져줘서 코드는 살림
    end
end