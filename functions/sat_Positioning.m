%%%%% sattelite Positioning function %%%%%
function [x y z] = sat_Positioning(Nav, TTX)
    c=299792458;                    % speed of light(meter/sec)
    tolerance = 1e-12;              %convergence threshold
    mu=3.986005e+14;                % earth's gravitational constant meter^3/sec^2
    OMGd_earth=7.2921151467e-5;     % earth rotation rate (rad/sec)
    
    ttx = TTX.TTX;
    P = TTX.P;
    sat_numbering = TTX.PRN;
    
    del_t = P./c;
    
    fixed_time = ttx; %fixed time;
    
    % variable %
    [sat, IODE, IODC, sva, svh, week, toes, ttrs, af0, af1, af2, A, e, i0,...
    OMG0, omg, M0, deln, OMGd, idot, crc, crs, cuc, cus, cic, cis, code,...
    flag, fit, tgd, toe, toc, ttr]=define_Sparameter(Nav,length(sat_numbering));
    
    t_k=zeros(length(sat_numbering),1);                   %time start on toe
    
    for i=1:length(sat_numbering)
        sNav=Nav{sat_numbering(i)};
        % if i==1
        %     E = datetime(1970,1,1);
        %     fprintf("fixed time = %f\n",fixed_time(1))
        %     disp(fixed_time(1))
        %     disp(datetime(fixed_time(1),'ConvertFrom','epochtime','Epoch',E));
        % end
        ttx_time=fixed_time(i);

        j=get_NavIndexbyTime(sNav, ttx_time);
        
        toe=sNav.toe(j);
        t_k(i)=ttx_time-toe;
        toes(i)=sNav.toes(j);
        toes(i) = sNav.toes(j);
        A(i) = sNav.A(j);                       % get semi major axis
        e(i) = sNav.e(j);                       % get eccentricity
        i0(i) = sNav.i0(j);                     % get inclination angle
        OMG0(i) = sNav.OMG0(j);                 % get longitude of ascending node
        omg(i) = sNav.omg(j);                   % get argument of perigee
        M0(i) = sNav.M0(j);                     % get mean anomaly
        deln(i) = sNav.deln(j);                 % get mean motion difference
        OMGd(i) = sNav.OMGd(j);                 % get rate of right ascension
        idot(i) = sNav.idot(j);                 % get rate of inclination angle
        crc(i) = sNav.crc(j);                   % get harmonic cos radius
        crs(i) = sNav.crs(j);                   % get harmonic sine radius
        cuc(i) = sNav.cuc(j);                   % get harmonic cos latitude
        cus(i) = sNav.cus(j);                   % get harmonic sine latitude
        cic(i) = sNav.cic(j);                   % get harmonic cos inclination
        cis(i) = sNav.cis(j);                   % get harmonic sine inclination
    end
    
    %%%%%   Computation of GPS Coordinate   %%%%%%%
    M_k=M0+(deln+sqrt(mu./A.^3)).*t_k;
    E_k = M_k;  % Initialize eccentric anomaly
    for k = 1:10  % Iterate to solve Kepler's equation
        E_kOld=E_k;
        E_k = M_k + e .* sin(E_k);
        if norm(E_k-E_kOld) < tolerance
            break
        end
    end
    nu_k=2.*atan2(sqrt(1.+e).*sin(E_k./2),sqrt(1.-e).*cos(E_k./2));
    ffi=nu_k+omg;
    u_k=ffi+cus.*sin(2.*(ffi))+cuc.*cos(2.*(ffi));
    r_k=A.*(1.-e.*cos(E_k))+crs.*sin(2.*(ffi))+crc.*cos(2.*(ffi));
    i_k=i0+idot.*t_k+cis.*sin(2.*(ffi))+cic.*cos(2.*(ffi));
    OMG_k=OMG0+(OMGd-OMGd_earth).*t_k-OMGd_earth.*toes;
    
    % Convert coordinates
    x_ko = r_k .* cos(u_k);
    y_ko = r_k .* sin(u_k);
    z = y_ko .* sin(i_k);
    
    x=x_ko.*cos(OMG_k)-y_ko.*cos(i_k).*sin(OMG_k);
    y=x_ko.*sin(OMG_k)+y_ko.*cos(i_k).*cos(OMG_k);
end