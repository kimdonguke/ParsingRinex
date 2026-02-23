%% Table grouping %%
%======================%
% grouping Table by specific(parameter) column(satellite ID)
%Input : Table class, colum name(string)
%Output : (total Group Number * 1) cell
%======================%
function [group] = groupingTable_PRN(dataTable)
    % groupingTable_PRN
    % PRN 번호와 Cell의 인덱스를 1:1로 일치시켜 그룹화합니다.
    % (결측된 PRN의 인덱스에는 빈 테이블이 들어갑니다.)
    
    % GPS 위성은 1~32번까지 있으므로 Cell 크기를 32로 고정
    MAX_PRN = 32; 
    group = cell(MAX_PRN, 1);
    
    % 원본 테이블과 형태(Column)가 똑같은 '빈 테이블' 템플릿 생성
    % 이렇게 해야 나중에 빈 셀을 참조했을 때 에러가 나지 않고 깔끔하게 처리됩니다.
    emptyTemplate = dataTable([], :);
    
    for prn = 1:MAX_PRN
        % 현재 PRN과 일치하는 데이터 인덱스 찾기
        idx = (dataTable.PRN == prn);
        
        if any(idx)
            % 데이터가 존재하면 해당 인덱스 방에 할당
            group{prn} = dataTable(idx, :);
        else
            % 위성이 unhealthy 이거나 데이터가 아예 없으면 빈 테이블 할당
            group{prn} = emptyTemplate;
        end
    end
end
