% ---------------- grouping Table ----------------
% input : (Table, variable Name)
% output : cell
% ------------------------------------------------
function [group] = groupingTable(dataTable, varName)
    parameterGroup=findgroups(dataTable.(varName));
    num_group=length(unique(parameterGroup));
    group = cell(num_group,1);
    for i=1:num_group
        group{i}=dataTable(parameterGroup==i,:);
    end
end



