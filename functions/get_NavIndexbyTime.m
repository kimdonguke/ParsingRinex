% ---------------- get_NavIndexbyTime ----------------
% input : (Table, variable Name)
% output : cell
% ------------------------------------------------
function [index] = get_NavIndexbyTime(Nav, time)        %Nav, double class  
    % 벡터화된 검색 (for 루프 제거)
    ttr_all = Nav.ttr(:);
    toe_all_raw = Nav.toe(:);
    
    % toe가 datetime인지 확인하고 처리
    if isdatetime(toe_all_raw)
        toe_all = seconds(toe_all_raw);
    else
        toe_all = double(toe_all_raw);
    end
    
    % 조건 만족하는 인덱스 찾기 (원래 코드 로직 유지: abs(seconds(toe-time))<=7200)
    valid_mask = (ttr_all < time) & (abs(toe_all - time) <= 7200);
    
    % 첫 번째 유효한 인덱스 반환
    idx_found = find(valid_mask, 1, 'first');
    if ~isempty(idx_found)
        index = idx_found;
    else
        % 원래 코드는 break로 루프를 빠져나갔지만, 여기서는 기본값 반환
        error(sprintf("Can't find Index"));
    end
end
