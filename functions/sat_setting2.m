function returnSatlist = sat_setting2(Nav, ttx)
    PRN = ttx.PRN;
    [x, y, z] = sat_Positioning(Nav,ttx);
    [az, el] = getAzEl(x,y,z);
    returnSatlist = table(PRN, x,y,z,az,el);
end