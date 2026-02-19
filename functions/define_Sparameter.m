% ---------------- define_Sparameter ----------------
% input : (Navigation, )
% output : each NAN column vector
% ------------------------------------------------
function [sat, IODE, IODC, sva, svh, week, toes, ttrs, af0, af1, af2, A, e, i0, OMG0, ...
omg, M0, deln, OMGd, idot, crc, crs, cuc, cus, cic, cis, code, flag, fit, tgd, toe, toc, ttr]=define_Sparameter(Nav, Nav_length)
    % variable %
    n_rows = Nav_length;
    
    % --- 1. Basic SV & Data Information ---
    sat  = nan(n_rows, 1);  % Satellite PRN number
    IODE = nan(n_rows, 1);  % Issue of Data, Ephemeris
    IODC = nan(n_rows, 1);  % Issue of Data, Clock
    sva  = nan(n_rows, 1);  % SV accuracy (User Range Accuracy)
    svh  = nan(n_rows, 1);  % SV health
    
    % --- 2. Clock Correction Parameters ---
    af0  = nan(n_rows, 1);  % SV clock bias (sec)
    af1  = nan(n_rows, 1);  % SV clock drift (sec/sec)
    af2  = nan(n_rows, 1);  % SV clock drift rate (sec/sec^2)
    tgd  = nan(n_rows, 1);  % Total Group Delay (sec)
    
    % --- 3. Time Parameters ---
    week = nan(n_rows, 1);  % GPS week number
    toe  = nan(n_rows, 1);  % Reference time ephemeris (sec of week)
    toc  = nan(n_rows, 1);  % Reference time clock (sec of week)
    ttr  = nan(n_rows, 1);  % Transmission time of message
    toes = nan(n_rows, 1);  % Toe from ephemeris (for subframe)
    ttrs = nan(n_rows, 1);  % Transmission time of subframe
    
    % --- 4. Keplerian Orbital Elements ---
    A    = nan(n_rows, 1);  % Semi-major axis (meters)
    e    = nan(n_rows, 1);  % Eccentricity
    i0   = nan(n_rows, 1);  % Inclination angle at reference epoch (rad)
    OMG0 = nan(n_rows, 1);  % Longitude of ascending node (rad)
    omg  = nan(n_rows, 1);  % Argument of perigee (rad)
    M0   = nan(n_rows, 1);  % Mean anomaly at reference epoch (rad)
    
    % --- 5. Perturbations & Rate Parameters ---
    deln = nan(n_rows, 1);  % Mean motion difference (rad/sec)
    OMGd = nan(n_rows, 1);  % Rate of right ascension (rad/sec)
    idot = nan(n_rows, 1);  % Rate of inclination angle (rad/sec)
    crc  = nan(n_rows, 1);  % Harmonic correction (cos, radius, m)
    crs  = nan(n_rows, 1);  % Harmonic correction (sine, radius, m)
    cuc  = nan(n_rows, 1);  % Harmonic correction (cos, lat, rad)
    cus  = nan(n_rows, 1);  % Harmonic correction (sine, lat, rad)
    cic  = nan(n_rows, 1);  % Harmonic correction (cos, incl, rad)
    cis  = nan(n_rows, 1);  % Harmonic correction (sine, incl, rad)
    
    % --- 6. Miscellaneous Flags ---
    code = nan(n_rows, 1);  % Codes on L2 channel
    flag = nan(n_rows, 1);  % L2 P data flag
    fit  = nan(n_rows, 1);  % Curve fit interval
end
