tic
fprintf("satellite positioning start....\n")
[ttx_epoch, ttx_table] = clock_Correction(return_OBS,return_NAV);
xyz_epoch = sat_setting(return_NAV, ttx_epoch);   % Nav(table), ttx(cell grouped by epoch)
xyz_table = sat_setting2(return_NAV, ttx_table);  % Nav(table), ttx(table ungrouped)
fprintf("satellite positioning finish. \n")
toc


%%%% SV Clock Correction %%%%
function [ttx_epoch, ttx_table] = clock_Correction(Obs,Nav)  % ttx_cell{i,1} = time, ttx_cell{i,2}=satellite number
    %% get Obs data %%
    %sObs=table2array(Obs);
    c=299792458;                    % speed of light(meter/sec)
    tolerance = 1e-12;              %convergence threshold
    mu=3.986005e+14;                % earth's gravitational constant meter^3/sec^2
    F=-2*sqrt(mu)/c^2;              % 상대성 효과 보정 힘
    L1 = 1575.42;                   % L1CA 주파수
    L2 = 1227.60;                   % L2 주파수
    mu2 = (L1/L2)^2;                % TGD 계수
    
    Obs = groupingTable(Obs, 'Time');
    
    ttx_list = [];
    ttx_epoch=cell(height(Obs),1);
    
    
    for i=1:height(Obs) % iteration by whole grouped data
        
        obs_i=Obs{i};
        trx = obs_i.Time(1);
        PRN=obs_i.PRN;
        t_k=0;
        index_num=0;
        
        TTX=nan(length(PRN),1);
        P=Obs{i}.P1;     %Pseudorange
        SNR = Obs{i}.SNR1; % 
        % code = Obs{i}.code;
        
        %%% algorithm problem on here
        for j=1:length(PRN) % iteration by each satellite, calculate t_tx
             sNav=Nav{PRN(j)};
             index_num=get_NavIndexbyTime(sNav,trx);
    
             M0=sNav.M0(index_num);
             A=sNav.A(index_num);
             toe=sNav.toe(index_num);
             toc=sNav.toc(index_num);
             deln=sNav.deln(index_num);
             e=sNav.e(index_num); 
             af0=sNav.af0(index_num);
             af1=sNav.af1(index_num);
             af2=sNav.af2(index_num);
             TGD=sNav.tgd(index_num);
             % if obs_i.code == 23
             %     TGD = mu2*TGD;
             % end
             
             ttx_new=obs_i.Time(j)-obs_i.P1(j)/c;   % define ttx_new to first iteration set ttx_new
             t_SV=ttx_new;
             ttx=0;
             % fprintf("ttx_new = ")
             % disp(fff(ttx_new))
             max_iter=1;
            
             while abs(ttx_new-ttx)>tolerance && max_iter<200
                max_iter=max_iter+1;

                ttx=ttx_new;
                t_k=ttx_new-toe; %  time from reference epoch
                M_k=M0+(sqrt(mu/(A^3))+deln)*t_k; % M_k = M_0 + n*t_k /// n=n_0+del n
                E_k=M_k;
                for k=1:10
                    E_kOld=E_k;
                    E_k=M_k+e*sin(E_k);
                    if norm(E_k-E_kOld)<tolerance
                        break
                    end
                end
                delt_r=F*sqrt(A)*e*sin(E_k);
                delt_SV=af0+af1*(ttx-toc)+af2*(ttx-toc)^2+delt_r-TGD;
                ttx_new=t_SV-delt_SV;
             end
             P(j)=P(j)+c*delt_SV;
             TTX(j)=ttx_new;
        end

        ttx_epoch{i}=table(PRN,TTX,P,SNR);
        ttx_list=[ttx_list; [PRN,TTX,P,SNR]];
    end
    PRN = ttx_list(:,1);
    TTX = ttx_list(:,2);
    P = ttx_list(:,3);
    SNR = ttx_list(:,4);
    ttx_table = table(PRN,TTX,P,SNR);
end

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


function [returnSatList] = sat_setting(Nav, ttx)

    Re = 6378137;      % 지구 반지름
    hm = 350000;       % 전리층 고도 (350km)
   
        
    sat_xyzlist=cell(height(ttx),1);
    
    for i=1:height(ttx)
        
        [x, y, z]=sat_Positioning(Nav,ttx{i}); %ttx_list, %grouped sat Obs Data
        PRN=ttx{i}.PRN;
        
        [az, el] = getAzEl(x,y,z);
        sat_xyzlist{i} = table(PRN,x,y,z,az,el);
    end
    returnSatList=sat_xyzlist;
end

function returnSatlist = sat_setting2(Nav, ttx)
    PRN = ttx.PRN;
    [x, y, z] = sat_Positioning(Nav,ttx);
    [az, el] = getAzEl(x,y,z);
    returnSatlist = table(PRN, x,y,z,az,el);
end


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