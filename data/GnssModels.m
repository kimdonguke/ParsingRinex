% ==============================================================================
% Class for GNSS models
%
% - Ionospheric model:
%   Broadcast GPS ionospheric model (Klobuchar)
%
% - Tropospheric model: 
%   Saastamoinen model with standard atmoshpere pres. & temp. + Niell mapping 
%   function
% 
% ------------------------------------------------------------------------------
% Yongrae Jo, 0727ggame@sju.ac.kr
% ==============================================================================
classdef GnssModels

    % Constants
    properties (Constant, Hidden)

        % Speed of light
        C_LIGHT = 299792458.0; 

        % Factor of ionospehric model standard deviation
        FACTOR_STD_KLOB = 0.5;

        % Standard deviation of Saastamoinen model
        STD_SAAS = 0.3;

        % Niell coefficients
        NIELL_HGTS = [2.53E-5, 5.49E-3, 1.14E-3];
        NIELL_LATS = [0, 15, 30, 45, 60, 75, 90];
        NIELL_COEF = [...
            1.2769934E-3, 1.2769934E-3, 1.2683230E-3, 1.2465397E-3, 1.2196049E-3, 1.2045996E-3, 1.2045996E-3; ...
            2.9153695E-3, 2.9153695E-3, 2.9152299E-3, 2.9288445E-3, 2.9022565E-3, 2.9024912E-3, 2.9024912E-3; ...
            62.610505E-3, 62.610505E-3, 62.837393E-3, 63.721774E-3, 63.824265E-3, 64.258455E-3, 64.258455E-3; ...
            0.0000000E-0, 0.0000000E-0, 1.2709626E-5, 2.6523662E-5, 3.4000452E-5, 4.1202191E-5, 4.1202191E-5; ...
            0.0000000E-0, 0.0000000E-0, 2.1414979E-5, 3.0160779E-5, 7.2562722E-5, 11.723375E-5, 11.723375E-5; ...
            0.0000000E-0, 0.0000000E-0, 9.0128400E-5, 4.3497037E-5, 84.795348E-5, 170.37206E-5, 170.37206E-5; ...
            5.8021897E-4, 5.8021897E-4, 5.6794847E-4, 5.8118019E-4, 5.9727542E-4, 6.1641693E-4, 6.1641693E-4; ...
            1.4275268E-3, 1.4275268E-3, 1.5138625E-3, 1.4572752E-3, 1.5007428E-3, 1.7599082E-3, 1.7599082E-3; ...
            4.3472961E-2, 4.3472961E-2, 4.6729510E-2, 4.3908931E-2, 4.4626982E-2, 5.4736038E-2, 5.4736038E-2; ...
            ];
    end

    % Methods
    methods (Access = public, Static)

        % ======================================================================
        % Klobuchar ionospheric model
        % 
        % args: 
        %   time : GPS time
        %   llh  : user geodetic position [rad, rad, m]
        %   azel : satellite azimuth and elevation angle [rad, rad]
        %   ceof : GPS broadcast ionospheric model parameters
        % 
        % returns:
        %   iono     : ionospheric delay [m]
        %   var_iono : variance of ionospheric delay [m^2]
        %
        % ----------------------------------------------------------------------
        % Yongrae Jo, 0727ggame@sju.ac.kr
        % ======================================================================
        function [iono, var_iono] = IonoModel(time, llh, azel, coef)
            arguments
                time (1,1) datetime
                llh  (1,3) double
                azel (:,2) double
                coef (1,8) double = [0.1118E-7, -0.7451E-8, -0.5961E-7, 0.1192E-6, ...
                                     0.1167E+6, -0.2294E+6, -0.1311E+6, 0.1049E+7]; % 2004/01/01
            end

            % Number of satellites
            nsat = size(azel, 1);

            % Check satellites
            if ~all(azel(:,2) >= 0.0 & azel(:,2) <= pi * 0.5)
                error("Invalid satellite elevation angle");
            end

            % Check coefficients
            if any(isnan(coef))
                error("Invalid GPS broadcast ionospheric parameters");
            end

            % Initialize
            iono     = zeros(nsat,1);
            var_iono = zeros(nsat,1);

            % Check receiver altitude
            if llh(3) < -1E3
                return
            end

            % Earth centered angle (semi-circle)
            psi = 0.0137 ./ (azel(:,2)./pi + 0.11) - 0.022;

            % Subionospheric latitude and longitude (semi-circle)
            phi = llh(1)./pi + psi .* cos(azel(:,1));
            phi(phi >  0.416) =  0.416;
            phi(phi < -0.416) = -0.416;

            lam = llh(2)./pi + psi .* sin(azel(:,1)) ./ cos(phi * pi);

            % Geometric latitude (semi-circle)
            phi = phi + 0.064 .* cos((lam - 1.617) * pi);

            % Local time
            tt = mod(43200 * lam + hour(time) * 3600 + minute(time) * 60 + second(time), 86400);

            % Slant factor
            f = 1 + 16 .* (0.53 - azel(:,2) ./ pi).^3;

            % Ionospheric delay
            amp = coef(1) + phi .* (coef(2) + phi .* (coef(3) + phi .* coef(4)));
            per = coef(5) + phi .* (coef(6) + phi .* (coef(7) + phi .* coef(8)));

            amp(amp < 0, 1) = 0;
            per(per < 72000, 1) = 72000;

            x = 2 .* pi .* (tt - 50400) ./ per;

            i = abs(x) < 1.57;

            iono( i,1) = GnssModels.C_LIGHT .* f( i) .* (5E-9 + amp(i) .* (1 + x(i).^2.*(-0.5 + x(i).^2./24)));
            iono(~i,1) = GnssModels.C_LIGHT .* f(~i) .* (5E-9);

            var_iono = (GnssModels.FACTOR_STD_KLOB * iono).^2;
        end

        % ======================================================================
        % Compute tropospheric delay using standard atmosphere and Saastamoinen 
        % model
        % 
        % args: 
        %   time : GPS time
        %   llh  : user geodetic position [rad, rad, m]
        %   azel : satellite azimuth and elevation angle [rad, rad]
        %   humi : relative humidity (optional for wet delay)
        % 
        % returns:
        %   trop     : tropospheric delay [m]
        %   var_trop : variance of tropospheric delay [m^2]
        %
        % ----------------------------------------------------------------------
        % Yongrae Jo, 0727ggame@sju.ac.kr
        % ======================================================================
        function [trop, var_trop] = TropModel(time, llh, azel, humi)

            arguments
                time (1,1) datetime
                llh  (1,3) double
                azel (:,2) double
                humi (1,1) double = 0;
            end

            % Number of satellites
            nsat = size(azel, 1);

            % Check satellites
            if ~all(azel(:,2) >= 0.0 & azel(:,2) <= pi * 0.5)
                error("Invalid satellite elevation angle");
            end

            % Initialize
            trop     = zeros(nsat,1);
            var_trop = zeros(nsat,1);

            % Check receiver altitude
            if llh(3) < -1E2 || llh(3) > 1E4
                return
            end

            % Standard atmosphere
            if llh(3) < 0
                hgt = 0;
            else
                hgt = llh(3);
            end

            pres = 1013.25 * (1.0 - 2.2557E-5 * hgt)^5.2568;
            temp = 15 - 6.5E-3 * hgt + 273.16;
            e    = 6.108 * humi * exp((17.15 * temp - 4684) / (temp - 38.45));

            % Mapping function
            [mapf_dry, mapf_wet] = GnssModels.TropMapFunc(time, llh, azel);

            % Dry and wet delays
            dry = 0.0022768 * pres / (1 - 0.00266 * cos(2 * llh(1)) - 0.00028 * hgt / 1E3) .* mapf_dry;
            wet = 0.002277 * (1225 / temp + 0.05) * e .* mapf_wet;

            % Total delay and error variance
            trop     = dry + wet;
            var_trop = ( GnssModels.STD_SAAS ./ (sin(azel(:,2)) + 0.1) ).^2;
        end

        % ======================================================================
        % Compute tropospheric delay mapping function by Niell mapping function
        % 
        % args: 
        %   time : GPS time
        %   llh  : user geodetic position [rad, rad, m]
        %   azel : satellite azimuth and elevation angle [rad, rad]
        % 
        % returns:
        %   mapf_dry : dry mapping function
        %   mapf_wet : wet mapping function
        %
        % ----------------------------------------------------------------------
        % Yongrae Jo, 0727ggame@sju.ac.kr
        % ======================================================================
        function [mapf_dry, mapf_wet] = TropMapFunc(time, llh, azel)

            arguments
                time (1,1) datetime
                llh  (1,3) double
                azel (:,2) double
            end

            GM = GnssModels;
            
            % Check satellites
            if ~all(azel(:,2) >= 0.0 & azel(:,2) <= pi * 0.5)
                error("Invalid satellite elevation angle");
            end

            % Latitude and height
            lat = rad2deg(llh(1));
            hgt = llh(3);

            % Year from DOY 28, added half year for southern latitudes
            y = (day(time, "dayofyear") - 28) / 365.25;
            if lat < 0, y = y + 0.5; end

            % Mapping function
            el   = azel(:,2);
            lat  = abs(lat);
            cosy = cos(2 * pi * y);

            % Coefficients
            ah = interp1(GM.NIELL_LATS, GM.NIELL_COEF(1,:), lat) ...
               - interp1(GM.NIELL_LATS, GM.NIELL_COEF(4,:), lat) * cosy;
            bh = interp1(GM.NIELL_LATS, GM.NIELL_COEF(2,:), lat) ...
               - interp1(GM.NIELL_LATS, GM.NIELL_COEF(5,:), lat) * cosy;
            ch = interp1(GM.NIELL_LATS, GM.NIELL_COEF(3,:), lat) ...
               - interp1(GM.NIELL_LATS, GM.NIELL_COEF(6,:), lat) * cosy;

            aw = interp1(GM.NIELL_LATS, GM.NIELL_COEF(7,:), lat);
            bw = interp1(GM.NIELL_LATS, GM.NIELL_COEF(8,:), lat);
            cw = interp1(GM.NIELL_LATS, GM.NIELL_COEF(9,:), lat);

            % Ellipsoidal height is used instead of height above sea level
            dm = (1 ./ sin(el) - GM.MapFunc(el, GM.NIELL_HGTS(1), GM.NIELL_HGTS(2), GM.NIELL_HGTS(3))) .* hgt ./ 1E3;

            mapf_dry = GM.MapFunc(el, ah, bh, ch) + dm;      % Dry mapping function
            mapf_wet = GM.MapFunc(el, aw, bw, cw);           % Wet mapping function
        end
    end

    % Methods (private)
    methods (Access = private, Static)

        % ======================================================================
        % Mapping function
        % 
        % ----------------------------------------------------------------------
        % Yongrae Jo, 0727ggame
        % ======================================================================
        function ret = MapFunc(el, a, b, c)

            sinel = sin(el);

            ret = (1 + a./(1 + b./(1 + c))) ./ (sinel + a./(sinel + b./(sinel + c)));
        end
    end
end