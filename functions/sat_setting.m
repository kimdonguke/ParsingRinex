function [returnSatCell] = sat_setting(Nav, ttx)

    Re = 6378137;      % 지구 반지름
    hm = 350000;       % 전리층 고도 (350km)
   
        
    sat_xyzlist=cell(height(ttx),1);
    
    for i=1:height(ttx)
        
        [x, y, z]=sat_Positioning(Nav,ttx{i}); %ttx_list, %grouped sat Obs Data
        PRN=ttx{i}.PRN;
        
        [az, el] = getAzEl(x,y,z);
        sat_xyzlist{i} = table(PRN,x,y,z,az,el);
    end
    returnSatCell=sat_xyzlist;
end