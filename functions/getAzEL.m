function [az, el] = getAzEl(x, y, z)
    ref_position=[-3047507.380, 4043980.305, 3865242.828];
    ref_position=ecef2lla(ref_position);
    ref_positionRAD=deg2rad(ref_position);
    
    [enu]=lla2enu(ecef2lla([x y z]),ref_position,'ellipsoid');
    e=enu(:,1);
    n=enu(:,2);
    u=enu(:,3);
    az=atan2(e,n);
    el = atan2(u,sqrt(e.^2+n.^2));     %calculate azimuth, elevation angle
end