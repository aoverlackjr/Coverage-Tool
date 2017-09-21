classdef CTS
% Object to maintain coordinate frame and time systems and the conversions
% thereof.
    
    properties
    end
    
    properties (Constant)
        
        % The J2000 to Geocentric Celestial Reference frame transformation
        % matrix.
        J2K2GCRF = [ 1              0.000273E-8    -9.740996E-8;...
                    -0.000273E-8    1              -1.324146E-8;...
                     9.740996E-8    1.324146E-8     1];
        
        % Arcsecond expressed in radians
        as_rad = pi/648000;
        
    end
    
    properties (Constant) % Constants for Vallado code
       
        % iar80       - integers for fk5 1980
        % rar80       - reals for fk5 1980  
        iar80 = [0 0 0 0 1;0 0 2 -2 2;0 0 2 0 2;0 0 0 0 2;0 1 0 0 0;1 0 0 0 0;0 1 2 -2 2;0 0 2 0 1;1 0 2 0 2;0 -1 2 -2 2;1 0 0 -2 0;0 0 2 -2 1;-1 0 2 0 2;1 0 0 0 1;0 0 0 2 0;-1 0 2 2 2;-1 0 0 0 1;1 0 2 0 1;2 0 0 -2 0;-2 0 2 0 1;0 0 2 2 2;2 0 2 0 2;2 0 0 0 0;1 0 2 -2 2;0 0 2 0 0;0 0 2 -2 0;-1 0 2 0 1;0 2 0 0 0;0 2 2 -2 2;-1 0 0 2 1;0 1 0 0 1;1 0 0 -2 1;0 -1 0 0 1;2 0 -2 0 0;-1 0 2 2 1;1 0 2 2 2;0 -1 2 0 2;0 0 2 2 1;1 1 0 -2 0;0 1 2 0 2;-2 0 0 2 1;0 0 0 2 1;2 0 2 -2 2;1 0 0 2 0;1 0 2 -2 1;0 0 0 -2 1;0 -1 2 -2 1;2 0 2 0 1;1 -1 0 0 0;1 0 0 -1 0;0 0 0 1 0;0 1 0 -2 0;1 0 -2 0 0;2 0 0 -2 1;0 1 2 -2 1;1 1 0 0 0;1 -1 0 -1 0;-1 -1 2 2 2;0 -1 2 2 2;1 -1 2 0 2;3 0 2 0 2;-2 0 2 0 2;1 0 2 0 0;-1 0 2 4 2;1 0 0 0 2;-1 0 2 -2 1;0 -2 2 -2 1;-2 0 0 0 1;2 0 0 0 1;3 0 0 0 0;1 1 2 0 2;0 0 2 1 2;1 0 0 2 1;1 0 2 2 1;1 1 0 -2 1;0 1 0 2 0;0 1 2 -2 0;0 1 -2 2 0;1 0 -2 2 0;1 0 -2 -2 0;1 0 2 -2 0;1 0 0 -4 0;2 0 0 -4 0;0 0 2 4 2;0 0 2 -1 2;-2 0 2 4 2;2 0 2 2 2;0 -1 2 0 1;0 0 -2 0 1;0 0 4 -2 2;0 1 0 0 2;1 1 2 -2 2;3 0 2 -2 2;-2 0 2 2 2;-1 0 0 0 2;0 0 -2 2 1;0 1 2 0 1;-1 0 4 0 2;2 1 0 -2 0;2 0 0 2 0;2 0 2 -2 1;2 0 -2 0 1;1 -1 0 -2 0;-1 0 0 1 1;-1 -1 0 2 1;0 1 0 1 0];
        rar80 = [-8.33860138961158e-05 -8.44545432492812e-08 4.46149790041051e-05 4.31484176187487e-09;-6.39323801279145e-06 -7.75701889775258e-10 2.78089127484430e-06 -1.50292241143956e-09;-1.10246631084308e-06 -9.69627362219072e-11 4.73662966444017e-07 -2.42406840554768e-10;9.99685810447863e-07 9.69627362219072e-11 -4.33908244593035e-07 2.42406840554768e-10;6.91344309262198e-07 -1.64836651577242e-09 2.61799387799149e-08 -4.84813681109536e-11;3.45187340949990e-07 4.84813681109536e-11 -3.39369576776675e-09 0;-2.50648673133630e-07 5.81776417331443e-10 1.08598264568536e-07 -2.90888208665722e-10;-1.87138080908281e-07 -1.93925472443814e-10 9.69627362219072e-08 0;-1.45928918013970e-07 0 6.25409648631301e-08 -4.84813681109536e-11;1.05204568800769e-07 -2.42406840554768e-10 -4.60572997054059e-08 1.45444104332861e-10;-7.66005616153067e-08 0 -4.84813681109536e-10 0;6.25409648631301e-08 4.84813681109536e-11 -3.39369576776675e-08 0;5.96320827764729e-08 0 -2.56951250988054e-08 0;3.05432619099008e-08 4.84813681109536e-11 -1.59988514766147e-08 0;3.05432619099008e-08 0 -9.69627362219072e-10 0;-2.86040071854626e-08 0 1.26051557088479e-08 0;-2.81191935043531e-08 -4.84813681109536e-11 1.55140377955052e-08 0;-2.47254977365863e-08 0 1.30899693899575e-08 0;2.32710566932577e-08 0 4.84813681109536e-10 0;2.23014293310387e-08 0 -1.16355283466289e-08 0;-1.84229198821624e-08 0 7.75701889775258e-09 0;-1.50292241143956e-08 0 6.30257785442397e-09 0;1.40595967521765e-08 0 -4.84813681109536e-10 0;1.40595967521765e-08 0 -5.81776417331443e-09 0;1.26051557088479e-08 0 -4.84813681109536e-10 0;-1.06659009844098e-08 0 0 0;1.01810873033003e-08 0 -4.84813681109536e-09 0;8.24183257886211e-09 -4.84813681109536e-11 0 0;-7.75701889775258e-09 4.84813681109536e-11 3.39369576776675e-09 0;7.75701889775258e-09 0 -3.87850944887629e-09 0;-7.27220521664304e-09 0 4.36332312998582e-09 0;-6.30257785442397e-09 0 3.39369576776675e-09 0;-5.81776417331443e-09 0 2.90888208665722e-09 0;5.33295049220490e-09 0 0 0;-4.84813681109536e-09 0 2.42406840554768e-09 0;-3.87850944887629e-09 0 1.45444104332861e-09 0;-3.39369576776675e-09 0 1.45444104332861e-09 0;-3.39369576776675e-09 0 1.45444104332861e-09 0;-3.39369576776675e-09 0 0 0;3.39369576776675e-09 0 -1.45444104332861e-09 0;-2.90888208665722e-09 0 1.45444104332861e-09 0;-2.90888208665722e-09 0 1.45444104332861e-09 0;2.90888208665722e-09 0 -1.45444104332861e-09 0;2.90888208665722e-09 0 0 0;2.90888208665722e-09 0 -1.45444104332861e-09 0;-2.42406840554768e-09 0 1.45444104332861e-09 0;-2.42406840554768e-09 0 1.45444104332861e-09 0;-2.42406840554768e-09 0 1.45444104332861e-09 0;2.42406840554768e-09 0 0 0;-1.93925472443814e-09 0 0 0;-1.93925472443814e-09 0 0 0;-1.93925472443814e-09 0 0 0;1.93925472443814e-09 0 0 0;1.93925472443814e-09 0 -9.69627362219072e-10 0;1.93925472443814e-09 0 -9.69627362219072e-10 0;-1.45444104332861e-09 0 0 0;-1.45444104332861e-09 0 0 0;-1.45444104332861e-09 0 4.84813681109536e-10 0;-1.45444104332861e-09 0 4.84813681109536e-10 0;-1.45444104332861e-09 0 4.84813681109536e-10 0;-1.45444104332861e-09 0 4.84813681109536e-10 0;-1.45444104332861e-09 0 4.84813681109536e-10 0;1.45444104332861e-09 0 0 0;-9.69627362219072e-10 0 4.84813681109536e-10 0;-9.69627362219072e-10 0 4.84813681109536e-10 0;-9.69627362219072e-10 0 4.84813681109536e-10 0;-9.69627362219072e-10 0 4.84813681109536e-10 0;-9.69627362219072e-10 0 4.84813681109536e-10 0;9.69627362219072e-10 0 -4.84813681109536e-10 0;9.69627362219072e-10 0 0 0;9.69627362219072e-10 0 -4.84813681109536e-10 0;9.69627362219072e-10 0 -4.84813681109536e-10 0;-4.84813681109536e-10 0 0 0;-4.84813681109536e-10 0 4.84813681109536e-10 0;-4.84813681109536e-10 0 0 0;-4.84813681109536e-10 0 0 0;-4.84813681109536e-10 0 0 0;-4.84813681109536e-10 0 0 0;-4.84813681109536e-10 0 0 0;-4.84813681109536e-10 0 0 0;-4.84813681109536e-10 0 0 0;-4.84813681109536e-10 0 0 0;-4.84813681109536e-10 0 0 0;-4.84813681109536e-10 0 0 0;-4.84813681109536e-10 0 0 0;-4.84813681109536e-10 0 4.84813681109536e-10 0;-4.84813681109536e-10 0 0 0;-4.84813681109536e-10 0 0 0;-4.84813681109536e-10 0 0 0;4.84813681109536e-10 0 0 0;4.84813681109536e-10 0 0 0;4.84813681109536e-10 0 -4.84813681109536e-10 0;4.84813681109536e-10 0 0 0;4.84813681109536e-10 0 -4.84813681109536e-10 0;4.84813681109536e-10 0 -4.84813681109536e-10 0;4.84813681109536e-10 0 0 0;4.84813681109536e-10 0 0 0;4.84813681109536e-10 0 0 0;4.84813681109536e-10 0 0 0;4.84813681109536e-10 0 0 0;4.84813681109536e-10 0 -4.84813681109536e-10 0;4.84813681109536e-10 0 0 0;4.84813681109536e-10 0 0 0;4.84813681109536e-10 0 0 0;4.84813681109536e-10 0 0 0;4.84813681109536e-10 0 0 0];
    
        % Default options for the true equator to mean equator
        % transformation (truemean)
        TRUEMEAN_ORDER    = 106;
        TRUEMEAN_EQETERMS = 2;
        TRUEMEAN_OPT      = 'a';
        
    end
    
    methods
        
        function cts = CTS()
            
        end

    end
    
    methods (Static)
       
        function R = R1(angle)
           
            R = [1  0          0;...
                 0  cos(angle) sin(angle);...
                 0 -sin(angle) cos(angle)];
            
        end
        
        function R = R2(angle)
           
            R = [cos(angle)  0 -sin(angle);...
                 0           1  0;...
                 sin(angle)  0  cos(angle)];
            
        end
        
        function R = R3(angle)
           
            R = [ cos(angle) sin(angle) 0;...
                 -sin(angle) cos(angle) 0;...
                  0          0          1];
            
        end
        
        function R = precession(UTC_date)
           
            % Earth precession matrix as part of transformation between ECI
            % and ECEF and TEME type reference frames.
            % McCarthy, 1996.
            
            % Get the julian centuries since JD2000.0
            [~,T] = CTS.juliandate(UTC_date);
            
            % Calculate the precession parameters
            z = 2306.2181 * CTS.as_rad * T + ...
                1.0946800 * CTS.as_rad * T^2 + ...
                0.0182030 * CTS.as_rad * T^3;
            
            th = 2004.3109 * CTS.as_rad * T - ...
                 0.42665   * CTS.as_rad * T^2 - ...
                 0.041833  * CTS.as_rad * T^3;
            
            ze = 2306.2181 * CTS.as_rad * T + ...
                 0.30188   * CTS.as_rad * T^2 + ...
                 0.017998  * CTS.as_rad * T^3;
            
            R = CTS.R3(-z) * CTS.R2(th) * CTS.R3(-ze);
            
        end
        
        function vec_out = j2k_to_gcrf(vec_in)
           
            vec_out = CTS.J2K2GCRF*[vec_in(1); vec_in(2); vec_in(3)];
            
        end
        
        function vec_out = gcrf_to_j2k(vec_in)
           
            vec_out = CTS.J2K2GCRF\[vec_in(1); vec_in(2); vec_in(3)];
            
        end
        
        function vec_out = wgs84_to_itrf90(vec_in)
           
            % Function to convert World Geodetic System 84 coordinates to
            % International Terrestrial Reference Frame 90.
            % McCarthy, 1996

            R_WGS84_to_ITRF90 = eye(3,3) + [0      -0.0070 -0.0003;...
                                            0.0070  0      -0.0183;...   
                                            0.0003  0.0183  0]*CTS.as_rad;
            
            vec_out = [0.060; -0.517; -0.223] + ...
                       0.999999989*R_WGS84_to_ITRF90*[vec_in(1);...
                                                      vec_in(2);...
                                                      vec_in(3)];
                                        
        end    
        
        function [JD, T] = juliandate(varargin) 
            % Returns the julian date for a UTC date and time, as well as
            % the julian centuries since J2000.0, meaning Jan 1, 12:00:00
            % for the year 2000.
            
            [year, month, day, hour, min, sec] = datevec(datenum(varargin{:}));
            
            for k = length(month):-1:1
                if ( month(k) <= 2 ) % january & february
                    year(k)  = year(k) - 1.0;
                    month(k) = month(k) + 12.0;
                end
            end

            JD = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
                floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5 + ...
                (hour + min/60 + sec/3600)/24;
            
            % Counting from JD2000.0
            dJD = JD - 2451545; 
            
            % Julian centuries since JD2000.0
            T = dJD/36525;
            
        end
        
        function [yr, mo, dy, hr, mn, sc] = jd_to_date(jd)
           
            b = floor(jd + 0.5) + 1537;
            c = floor((b - 122.1)/365.25);
            d = floor(365.25*c);
            e = floor((b-d)/30.6001);
            dy_fr = jd + 0.5 - floor(jd + 0.5);
            dy = b-d-floor(30.6001*e);
            mo = e-1-12*floor(e/14);
            yr = c-4715-floor((7+mo)/10);
            hr_fr = dy_fr*24;
            hr = floor(hr_fr);
            min_fr = (hr_fr - hr)*60;
            mn = floor(min_fr);
            sc = floor((min_fr - mn)*60);
            
        end
        
        function [DCM] = nadir_fixed_z_down_frame(pos,vel)
            % Z axis along geocentric position, towards nadir
            % X axis in plane and positve towards velocity vector
            % Y axis completes right handed frame
            z = -[pos(1) pos(2) pos(3)]/norm([pos(1) pos(2) pos(3)]);
            v =  [vel(1) vel(2) vel(3)]/norm([vel(1) vel(2) vel(3)]);
            
            y = cross(z,v);
            y = y/norm(y);
            
            x = cross(y,z);
            x = x/norm(x);
            
            DCM = [x;y;z];
        end
        
        function [DCM] = nadir_fixed_z_up_frame(pos,vel)
            % Z axis along geocentric position, towards zenith
            % X axis in plane and positve towards velocity vector
            % Y axis completes right handed frame
            z =  [pos(1) pos(2) pos(3)]/norm([pos(1) pos(2) pos(3)]);
            v =  [vel(1) vel(2) vel(3)]/norm([vel(1) vel(2) vel(3)]);
            
            y = cross(z,v);
            y = y/norm(y);
            
            x = cross(y,z);
            x = x/norm(x);
            
            DCM = [x;y;z];
        end
        
        function [DCM] = velocity_fixed_z_down_frame(pos, vel)
            % Z in plane positive towards nadir vector
            % X axis aligned with velocity vector
            % Y axis completes right handed frame
            r =  [pos(1) pos(2) pos(3)]/norm([pos(1) pos(2) pos(3)]);
            x =  [vel(1) vel(2) vel(3)]/norm([vel(1) vel(2) vel(3)]);
            
            y = cross(x,r);
            y = y/norm(y);
            
            z = cross(x,y);
            z = z/norm(z);
            
            DCM = [x;y;z];
        end
        
        function [DCM] = velocity_fixed_z_up_frame(pos,vel)
            % Z in plane positive towards zenith vector
            % X axis aligned with velocity vector
            % Y axis completes right handed frame
            r =  [pos(1) pos(2) pos(3)]/norm([pos(1) pos(2) pos(3)]);
            x =  [vel(1) vel(2) vel(3)]/norm([vel(1) vel(2) vel(3)]);
            
            y = cross(r,x);
            y = y/norm(y);
            
            z = cross(x,y);
            z = z/norm(z);
            
            DCM = [x;y;z];
        end
        
        function [DCM] = yaw_pitch_roll_frame(yaw_rad, pitch_rad, roll_rad)
           % Provides the frame of a body under yaw, pitch and roll
           % rotations, with the convention that the rotations are applied
           % in the following order: yaw, pitch, roll. Roll is around the x
           % axis, pitch around the y-axis and yaw around the z-axis.
           R_roll  = CTS.R1(roll_rad);
           R_pitch = CTS.R2(pitch_rad);
           R_yaw   = CTS.R3(yaw_rad);
           
           DCM = (R_roll*R_pitch*R_yaw);
           
        end
        
        function [DCM] = Rbi_from_orbit_frame(Roi, yaw_rad, pitch_rad, roll_rad)
            % Provides the inertial orientation frame of an object with its
            % yaw, pitch and roll orientations with respect to the orbital
            % frame definition of choice.
           Rbo = CTS.yaw_pitch_roll_frame(yaw_rad, pitch_rad, roll_rad);
           DCM = Rbo*Roi;
        end
        
        function [DCM] = relative_pointing_frame(first_pointing_vector, first_axis_str, second_pointing_vector, second_axis_str)
           
            % Returns right handed DCM based on a primary vector which is
            % aligned with the first axis given, where the second vector
            % further determines the orientation of the frame
            
            pv1 = first_pointing_vector/norm(first_pointing_vector);
            pv2 = second_pointing_vector/norm(second_pointing_vector);
            
            x = [NaN;NaN;NaN];
            y = [NaN;NaN;NaN];
            z = [NaN;NaN;NaN];
            
            switch first_axis_str
                
                case 'x'
                    
                    x = [pv1(1); pv1(2); pv1(3)];
                    
                    switch second_axis_str
                        
                        case 'y'
                            
                            z = cross(x,pv2);
                            z = z/norm(z);
                            
                            y = cross(z,x);
                            y = y/norm(y);
                            
                        case 'z'
                            
                            y = cross(pv2,x);
                            y = y/norm(y);
                            
                            z = cross(x,y);
                            z = z/norm(z);
                            
                        otherwise
                            
                            disp('CTS: Second axis label cannot be the same as first axis label')
                            
                    end
                    
                case 'y'
                    
                    y = [pv1(1); pv1(2); pv1(3)];
                    
                    switch second_axis_str
                        
                        case 'x'
                            
                            z = cross(pv2,y);
                            z = z/norm(z);
                            
                            x = cross(y,z);
                            x = x/norm(x);
                            
                        case 'z'
                            
                            x = cross(y,pv2);
                            x = x/norm(x);
                            
                            z = cross(x,y);
                            z = z/norm(z);
                            
                        otherwise
                            
                            disp('CTS: Second axis label cannot be the same as first axis label')
                            
                    end
                    
                case 'z'
                    
                    z = [pv1(1); pv1(2); pv1(3)];
                    
                    switch second_axis_str
                        
                        case 'x'
                            
                            y = cross(z,pv2);
                            y = y/norm(y);
                            
                            x = cross(y,z);
                            x = x/norm(x);
                            
                        case 'y'
                            
                            x = cross(pv2,z);
                            x = x/norm(x);
                            
                            y = cross(z,x);
                            y = y/norm(y);

                            
                        otherwise
                            
                            disp('CTS: Second axis label cannot be the same as first axis label')
                            
                    end
                    
                otherwise
                    
                    disp('CTS: Did not recognize first axis label')
                    
            end
            
            DCM = [x,y,z]';
            
            
        end
        
        function [RI2T] = ECI2TEME(time_UTC)

            [~, T] = CTS.juliandate(time_UTC); 
            [prec,~,~,~,~] = CTS.precess(T, '80');
            [~,~,~,~,~,nutteme] = CTS.truemean(T, CTS.TRUEMEAN_ORDER, CTS.TRUEMEAN_EQETERMS, CTS.TRUEMEAN_OPT);
            RI2T = nutteme'*prec';
        
        end
        
        function [q] = DCM2Q(DCM)
            
            dcm11 = DCM(1,1);
            dcm12 = DCM(1,2);
            dcm13 = DCM(1,3);
            
            dcm21 = DCM(2,1);
            dcm22 = DCM(2,2);
            dcm23 = DCM(2,3);
            
            dcm31 = DCM(3,1);
            dcm32 = DCM(3,2);
            dcm33 = DCM(3,3);
            
            divisor(1) = 0.5*sqrt( 1 + dcm11 + dcm22 + dcm33);
            divisor(2) = 0.5*sqrt( 1 + dcm11 - dcm22 - dcm33);
            divisor(3) = 0.5*sqrt( 1 - dcm11 - dcm22 + dcm33);
            divisor(4) = 0.5*sqrt( 1 - dcm11 + dcm22 - dcm33);

            absdivisor = abs(divisor);

            [~,indexd] = max(absdivisor);

            if      indexd == 1

                q4 = divisor(1);
                q1 = 0.25*(dcm23-dcm32)/q4;
                q2 = 0.25*(dcm31-dcm13)/q4;
                q3 = 0.25*(dcm12-dcm21)/q4;

            elseif  indexd == 2

                q1 = divisor(2);
                q2 = 0.25*(dcm12+dcm21)/q1;
                q3 = 0.25*(dcm13+dcm31)/q1;
                q4 = 0.25*(dcm23-dcm32)/q1;

            elseif  indexd == 3

                q3 = divisor(3);
                q1 = 0.25*(dcm13+dcm31)/q3;
                q2 = 0.25*(dcm23+dcm32)/q3;
                q4 = 0.25*(dcm12-dcm21)/q3;

            elseif  indexd == 4

                q2 = divisor(4);
                q1 = 0.25*(dcm12+dcm21)/q2;
                q3 = 0.25*(dcm23+dcm32)/q2;
                q4 = 0.25*(dcm31-dcm13)/q2;

            else
                display('quaternion error')
            end

            q = [q1;q2;q3;q4];
            
        end
        
        function [DCM] = Q2DCM(q)
            
            q = q/norm(q);
            
            q1 = q(1);
            q2 = q(2);
            q3 = q(3);
            q4 = q(4);
            
            dcm11 = q1^2-q2^2-q3^2+q4^2;
            dcm12 = 2*(q1*q2+q3*q4);
            dcm13 = 2*(q1*q3-q2*q4);
            dcm21 = 2*(q1*q2-q3*q4);
            dcm22 = -q1^2+q2^2-q3^2+q4^2;
            dcm23 = 2*(q2*q3+q1*q4);
            dcm31 = 2*(q1*q3+q2*q4);
            dcm32 = 2*(q2*q3-q1*q4);
            dcm33 = -q1^2-q2^2+q3^2+q4^2;

            DCM =  [dcm11 dcm12 dcm13;
                    dcm21 dcm22 dcm23;
                    dcm31 dcm32 dcm33];
            
        end
        
    end
    
    methods (Static) % Functions by Vallado
        
        function [deltapsi,trueeps,meaneps,omega,eqe,nutteme] = truemean(ttt, order, eqeterms, opt)
        % ----------------------------------------------------------------------------
        %
        %                           function truemean
        %
        %  this function forms the transformation matrix to go between the
        %    norad true equator mean equinox of date and the mean equator mean equinox
        %    of date (eci).  the results approximate the effects of nutation and
        %    precession.
        %
        %  author        : david vallado                  719-573-2600   25 jun 2002
        %
        %  revisions
        %    vallado     - fixes to order                                29 sep 2002
        %    vallado     - fixes to all options                           6 may 2003
        %
        %  inputs          description                    range / units
        %    ttt         - julian centuries of tt
        %    order       - number of terms for nutation   4, 50, 106, ...
        %    eqeterms    - number of terms for eqe        0, 2
        %    opt         - option for processing          a - complete nutation
        %                                                 b - truncated nutation
        %                                                 c - truncated transf matrix
        %
        %  outputs       :
        %    nutteme     - matrix for mod - teme - an approximation for nutation
        %
        %  locals        :
        %    prec        - matrix for mod - j2000
        %    tm          - combined matrix for teme
        %    ttt2        - ttt squared
        %    ttt3        - ttt cubed
        %    l           -                                rad
        %    ll          -                                rad
        %    f           -                                rad
        %    d           -                                rad
        %    omega       -                                rad
        %    deltapsi    - nutation angle                 rad
        %    deltaeps    - change in obliquity            rad
        %    eps         - mean obliquity of the ecliptic rad
        %    trueeps     - true obliquity of the ecliptic rad
        %    meaneps     - mean obliquity of the ecliptic rad
        %
        %  coupling      :
        %
        %
        %  references    :
        %    vallado       2004, 230
        %
        % [deltapsi,trueeps,meaneps,omega,eqe,nutteme] = truemean ( ttt,order,eqeterms,opt );
        % ----------------------------------------------------------------------------
            deg2rad = pi/180.0;

            % ---- determine coefficients for iau 1980 nutation theory ----
            ttt2= ttt*ttt;
            ttt3= ttt2*ttt;
            ttt4= ttt2*ttt2;

            meaneps = -46.8150 *ttt - 0.00059 *ttt2 + 0.001813 *ttt3 + 84381.448;
            meaneps = rem( meaneps/3600.0 ,360.0  );
            meaneps = meaneps * deg2rad;

            l    =  134.96340251  + ( 1717915923.2178 *ttt + ...
                    31.8792 *ttt2 + 0.051635 *ttt3 - 0.00024470 *ttt4 ) / 3600.0;
            l1   =  357.52910918  + (  129596581.0481 *ttt - ...
                     0.5532 *ttt2 - 0.000136 *ttt3 - 0.00001149*ttt4 )  / 3600.0;
            f    =   93.27209062  + ( 1739527262.8478 *ttt - ...
                    12.7512 *ttt2 + 0.001037 *ttt3 + 0.00000417*ttt4 )  / 3600.0;
            d    =  297.85019547  + ( 1602961601.2090 *ttt - ...
                     6.3706 *ttt2 + 0.006593 *ttt3 - 0.00003169*ttt4 )  / 3600.0;
            omega=  125.04455501  + (   -6962890.2665 *ttt + ...
                     7.4722 *ttt2 + 0.007702 *ttt3 - 0.00005939*ttt4 )  / 3600.0;

            l    = rem( l,360.0  )     * deg2rad;
            l1   = rem( l1,360.0  )    * deg2rad;
            f    = rem( f,360.0  )     * deg2rad;
            d    = rem( d,360.0  )     * deg2rad;
            omega= rem( omega,360.0  ) * deg2rad;

            deltapsi= 0.0;
            deltaeps= 0.0;

            for i = 1:order   % the eqeterms in nut80.dat are already sorted
                tempval= CTS.iar80(i,1)*l + CTS.iar80(i,2)*l1 + CTS.iar80(i,3)*f + ...
                         CTS.iar80(i,4)*d + CTS.iar80(i,5)*omega;
                deltapsi= deltapsi + (CTS.rar80(i,1)+CTS.rar80(i,2)*ttt) * sin( tempval );
                deltaeps= deltaeps + (CTS.rar80(i,3)+CTS.rar80(i,4)*ttt) * cos( tempval );
            end

            % --------------- find nutation parameters --------------------
            deltapsi = rem( deltapsi,360.0  ) * deg2rad;
            deltaeps = rem( deltaeps,360.0  ) * deg2rad;
            trueeps  = meaneps + deltaeps;
            % fprintf(1,'dpsi %16.9f   deps %16.9f  trueps %16.8f meaneps %16.8f degpre \n',deltapsi/deg2rad,deltaeps/deg2rad,trueeps/deg2rad, (trueeps-deltaeps)/deg2rad);
            cospsi  = cos(deltapsi);
            sinpsi  = sin(deltapsi);
            coseps  = cos(meaneps);
            sineps  = sin(meaneps);
            costrueeps = cos(trueeps);
            sintrueeps = sin(trueeps);

            jdttt = ttt*36525.0 + 2451545.0;
            % small disconnect with ttt instead of ut1
            if (jdttt > 2450449.5 ) && (eqeterms > 0)
                eqe= deltapsi* cos(meaneps) ...
                    + 0.00264*pi /(3600*180)*sin(omega) ...
                    + 0.000063*pi /(3600*180)*sin(2.0 *omega);
            else
                eqe= deltapsi* cos(meaneps);
            end

            nut(1,1) =  cospsi;
            nut(1,2) =  costrueeps * sinpsi;
            
            if (opt=='b')
                nut(1,2) = 0.0;
            end;
            
            nut(1,3) =  sintrueeps * sinpsi;
            nut(2,1) = -coseps * sinpsi;
            
            if (opt=='b')
                nut(2,1) = 0.0;
            end;
            nut(2,2) =  costrueeps * coseps * cospsi + sintrueeps * sineps;
            nut(2,3) =  sintrueeps * coseps * cospsi - sineps * costrueeps;
            nut(3,1) = -sineps * sinpsi;
            nut(3,2) =  costrueeps * sineps * cospsi - sintrueeps * coseps;
            nut(3,3) =  sintrueeps * sineps * cospsi + costrueeps * coseps;

            st(1,1) =  cos(eqe);
            st(1,2) = -sin(eqe);
            st(1,3) =  0.0;
            st(2,1) =  sin(eqe);
            st(2,2) =  cos(eqe);
            st(2,3) =  0.0;
            st(3,1) =  0.0;
            st(3,2) =  0.0;
            st(3,3) =  1.0;

            nutteme = nut*st;

            if (opt=='c')
                nutteme(1,1) =  1.0;
                nutteme(1,2) =  0.0;
                nutteme(1,3) =  deltapsi * sineps;
                nutteme(2,1) =  0.0;
                nutteme(2,2) =  1.0;
                nutteme(2,3) =  deltaeps;
                nutteme(3,1) = -deltapsi * sineps;
                nutteme(3,2) = -deltaeps;
                nutteme(3,3) =  1.0;
            end
        end
        
        function [prec,psia,wa,ea,xa] = precess(ttt, opt)
        %
        % ----------------------------------------------------------------------------
        %
        %                           function precess
        %
        %  this function calulates the transformation matrix that accounts for the effects
        %    of precession. both the 1980 and 2000 theories are handled. note that the 
        %    required parameters differ a little. 
        %
        %  author        : david vallado                  719-573-2600   25 jun 2002
        %
        %  revisions
        %    vallado     - consolidate with iau 2000                     14 feb 2005
        %
        %  inputs          description                    range / units
        %    ttt         - julian centuries of tt
        %    opt         - method option                  '01', '02', '96', '80'
        %
        %  outputs       :
        %    prec        - transformation matrix for mod - j2000 (80 only)
        %    psia        - cannonical precession angle    rad    (00 only)
        %    wa          - cannonical precession angle    rad    (00 only)
        %    ea          - cannonical precession angle    rad    (00 only)
        %    xa          - cannonical precession angle    rad    (00 only)
        %
        %  locals        :
        %    ttt2        - ttt squared
        %    ttt3        - ttt cubed
        %    zeta        - precession angle               rad
        %    z           - precession angle               rad
        %    theta       - precession angle               rad
        %    oblo        - obliquity value at j2000 epoch "%
        %
        %  coupling      :
        %    none        -
        %
        %  references    :
        %    vallado       2004, 214-216, 219-221
        %
        % [prec,psia,wa,ea,xa] = precess ( ttt, opt );
        % ----------------------------------------------------------------------------
            %sethelp;

            % " to rad
            convrt = pi / (180.0*3600.0);
            ttt2= ttt * ttt;
            ttt3= ttt2 * ttt;

            % ------------------- fk4 b1950 precession angles --------------------
            if (strcmp(opt,'50')==1)

                % ---- Seidelmann pg 107             
                % for these calls, ttt will come in with the current jd
    %            t1= 0.0; %(ttt - 2433282.42345905)/365242.198782;  % set start as B1850, 0.0 to simplify
    %            t2= (ttt - 2396758.203)/365242.198782;  % uses B1850
    %            ttt = t2 - t1;
    %            ttt2 = ttt * ttt;
    %            ttt3 = ttt * ttt2;
    %            fprintf(1,'50prec %15.9f  \n',ttt);
                % exp supp 61 pg 38
    %            psia = 50.3708 + 0.0050 * ttt;
    %            wa = 0.0;
    %            ea = 0.0;
    %            xa = 0.1247 - 0.0188 * ttt;
    %            zeta =  (23035.545 + 139.720*t1 + 0.060 *t1*t1)*ttt + (30.240-0.270*t1)*ttt2 + 17.995*ttt3; % "
    %            theta=  (20051.12 - 85.29*t1 - 0.37 *t1*t1)*ttt + (-42.65-0.37*t1)*ttt2 - 41.80*ttt3;
    %            z    =  (23035.545 + 139.720*t1 + 0.060 *t1*t1)*ttt + (109.480+0.390*t1)*ttt2 + 18.325*ttt3;

                % ---- Newcomb Exp Supp 61 approach, but see GTDS pg 3-17
                % Exp Supp 61 says use 1900? but gtds says use 1950. 
                % for these calls, ttt will come in with the current jd
                t1= 0.0; %(ttt - 2415020.31352)/36524.2198782;  % set start as B1900, 0.0 to simplify
    %            t2= (ttt - 2415020.31352)/36525;  % uses B1900   
                t2= (ttt - 2433282.42345905)/36525;  % uses B1950  
    %            ttt = t2 - t1;
    %            ttt2 = ttt * ttt;
    %            ttt3 = ttt * ttt2;
                % exp supp 61 pg 38
                psia = 50.3708 + 0.0050 * ttt;
                wa = 0.0; % not sure which one is which...
                ea = 84428.26 -   46.845*ttt - 0.00059*ttt2 + 0.00181*ttt3;
                xa = 0.1247 - 0.0188 * ttt;
               %fprintf(1,'50prec %15.9f  \n',ttt);
               % seems like Exp supp 61 is ok with 1900 as epoch, and Seidlemann is ok with listed measr,  
    %            zeta =  (2304.25 + 1.396*t1)*ttt + 0.302*ttt2 + 0.018*ttt3; % "
    %            theta=  (2004.682 - 0.853*t1)*ttt -0.426*ttt2 - 0.042*ttt3;
    %            z    =  (2304.25 + 1.396*t1)*ttt + 1.093*ttt2 + 0.018*ttt3;
                % GTDS pg 3-17 using days from 1950 - avoids long rpecession
                % constants...
                zeta =  2304.9969*ttt + 0.302*ttt2 + 0.01808*ttt3; % "
                theta=  2004.2980*ttt -0.425936*ttt2 - 0.0416*ttt3;
                z    =  2304.9969*ttt + 1.092999*ttt2 + 0.0192*ttt3;


                % tp-008 36-45
                % ttt is tropical centruies from 1950 36524.22 days
                prec(1,1) =  1.0 - 2.9696e-4 * ttt2 - 1.3e-7 * ttt3;
                prec(1,2) =  2.234941e-2 * ttt + 6.76e-6 * ttt2 - 2.21e-6 * ttt3;
                prec(1,3) =  9.7169e-3 * ttt - 2.07e-6 * ttt2 - 9.6e-7 * ttt3;
                prec(2,1) =  -prec(1,2);
                prec(2,2) =  1.0 - 2.4975e-4 * ttt2 - 1.5e-7 * ttt3;
                prec(2,3) =  - 1.0858e-4 * ttt2;
                prec(3,1) =  -prec(1,3);
                prec(3,2) =  prec(2,3);
                prec(3,3) =  1.0 - 4.721e-5 * ttt2;


               % pass these back out for testing
               psia = zeta;
               wa = theta;
               ea = z;
            % ------------------- iau 76 precession angles --------------------
            else
              if (strcmp(opt,'80')==1)
    %     fprintf(1,'80prec %15.9f  \n',ttt);
                psia =             5038.7784*ttt - 1.07259*ttt2 - 0.001147*ttt3; % "
                wa   = 84381.448                 + 0.05127*ttt2 - 0.007726*ttt3;
                ea   = 84381.448 -   46.8150*ttt - 0.00059*ttt2 + 0.001813*ttt3;
                xa   =               10.5526*ttt - 2.38064*ttt2 - 0.001125*ttt3;
                zeta =             2306.2181*ttt + 0.30188*ttt2 + 0.017998*ttt3; % "
                theta=             2004.3109*ttt - 0.42665*ttt2 - 0.041833*ttt3;
                z    =             2306.2181*ttt + 1.09468*ttt2 + 0.018203*ttt3;
              % ------------------ iau 03 precession angles -------------------
              else
                oblo =  84381.406; % "
                psia =  (((( -0.0000000951 * ttt + 0.000132851 ) * ttt - 0.00114045 ) * ttt - 1.0790069 ) * ttt + 5038.481507 ) * ttt; % "
                wa   =  ((((  0.0000003337 * ttt - 0.000000467 ) * ttt - 0.00772503 ) * ttt + 0.0512623 ) * ttt -    0.025754 ) * ttt + oblo;
                ea   =  (((( -0.0000000434 * ttt - 0.000000576 ) * ttt + 0.00200340 ) * ttt - 0.0001831 ) * ttt -   46.836769 ) * ttt + oblo;
                xa   =  (((( - 0.0000000560 * ttt + 0.000170663 ) * ttt - 0.00121197 ) * ttt - 2.3814292 ) * ttt +   10.556403 ) * ttt;

                zeta =  (((( - 0.0000003173 * ttt - 0.000005971 ) * ttt + 0.01801828 ) * ttt + 0.2988499 ) * ttt + 2306.083227 ) * ttt + 2.650545; % "
                theta=  (((( - 0.0000001274 * ttt - 0.000007089 ) * ttt - 0.04182264 ) * ttt - 0.4294934 ) * ttt + 2004.191903 ) * ttt;
                z    =  ((((   0.0000002904 * ttt - 0.000028596 ) * ttt + 0.01826837 ) * ttt + 1.0927348 ) * ttt + 2306.077181 ) * ttt - 2.650545;
              end
            end

            % convert units to rad
            psia = psia  * convrt; % rad
            wa   = wa    * convrt;
            ea   = ea    * convrt;
            xa   = xa    * convrt;
            zeta = zeta  * convrt; 
            theta= theta * convrt;
            z    = z     * convrt;
            %iauhelp='y';         
            %fprintf(1,'pr %11.7f  %11.7f  %11.7fdeg \n',zeta*180/pi,theta*180/pi,z*180/pi );
%             if (iauhelp == 'y')
%                 fprintf(1,'pr %11.7f  %11.7f  %11.7f %11.7fdeg \n',psia*180/pi,wa*180/pi,ea*180/pi,xa*180/pi );
%                 fprintf(1,'pr %11.7f  %11.7f  %11.7fdeg \n',zeta*180/pi,theta*180/pi,z*180/pi );
%             end

            if (strcmp(opt,'80')==1) % || (strcmp(opt,'50') == 1)
                coszeta  = cos(zeta);
                sinzeta  = sin(zeta);
                costheta = cos(theta);
                sintheta = sin(theta);
                cosz     = cos(z);
                sinz     = sin(z);

                % ----------------- form matrix  mod to j2000 -----------------
                prec(1,1) =  coszeta * costheta * cosz - sinzeta * sinz;
                prec(1,2) =  coszeta * costheta * sinz + sinzeta * cosz;
                prec(1,3) =  coszeta * sintheta;
                prec(2,1) = -sinzeta * costheta * cosz - coszeta * sinz;
                prec(2,2) = -sinzeta * costheta * sinz + coszeta * cosz;
                prec(2,3) = -sinzeta * sintheta;
                prec(3,1) = -sintheta * cosz;
                prec(3,2) = -sintheta * sinz;
                prec(3,3) =  costheta;

                % ----------------- do rotations instead ----------------------
                % p1 = rot3mat( z );
                % p2 = rot2mat( -theta );
                % p3 = rot3mat( zeta );
                % prec = p3*p2*p1;
            elseif (strcmp(opt,'50') ~= 1)
                oblo = oblo * convrt; % " to rad
                a4  = rot3mat(-xa);
                a5  = rot1mat(wa);
                a6  = rot3mat(psia);
                a7  = rot1mat(-oblo);
                prec = a7*a6*a5*a4;
               % prec = zeros(3);
            end
        
        end
    end
    
end

