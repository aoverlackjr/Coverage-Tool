classdef Kepler < handle

    % Kepler orbit object, written to be part of a more elaborate space
    % object for simulation. The principle of operation is to 'insert' an
    % object into a Kepler orbit at a given simulation epoch. From this
    % point the object can be propagated in for instance a payload cluster
    % deployment framework.
    
    % The kepler state is kept in memory, where the
    % (first) kepler orbit is initiated by the insert function, wich is
    % based on cartesian position and velocity data, at a reference time in
    % seconds. The reference time (as input) should on a whole be
    % consistent throughout the simulation.
    % The Kepler orbit can be propagated by simply polling the get_pos_vel
    % function at increasing times from the epoch/init_time.
    % If a manuever is issued, the internal kepler state is altered,
    % including the init time, to make sure the propagation on the basis of
    % the general timeline is consistent.
    
    % Currently, this object is not written (yet) to allow orbit creation
    % with kepler elements directly, although one can manually set the
    % public Kepler elements and init time and mean anomaly, after which 
    % the other functions will work nominally.
    
    % There are three ways to configure the Kepler orbit:
    % 1) To insert an object into orbit directly with position and velocity values
    % 2) To configure SSO values
    % 3) To manually set the six Kepler elements
    % 
    % After configuration the functions maneuver and get_pos_vel will work.
    
    properties
        
        % The Kepler state vector
        semi_major_axis;
        eccentricity;
        inclination;
        argument_of_perigee;
        right_ascension;
        mean_anomaly;
        
        init_time;
        init_mean_anomaly;
        
        VERBOSE = 0;
        notifications;
        history;
    end
    
    properties (Access = 'private')
        mu = 398600.4418E9;
        mu_alt;
        Re = 6378135;
        J2 = 0.00108263;
        last_update_time = 0;
    end
    
    methods
        
        function obj = Kepler()
           obj.notify('Kepler orbit created. use .insert or .set_SSO functions to create Kepler elements');
        end
        
        function set_keps(obj,a_km,e,i_deg,argper_deg,raan_deg, ma_deg, init_time)
           
           obj.init_time = init_time;
           obj.last_update_time = init_time;
           toRad = pi/180;
           
           obj.semi_major_axis     = a_km*1000;
           obj.eccentricity        = e;
           obj.inclination         = i_deg*toRad;
           obj.argument_of_perigee = argper_deg*toRad;
           obj.right_ascension     = raan_deg*toRad;
           obj.mean_anomaly        = ma_deg*toRad;
           obj.init_mean_anomaly   = ma_deg*toRad;
            
        end
        
        function insert(obj, pos, vel, init_time)
           % Generate the kepler vector in this object through inertial
           % reference of velocity and position.
           
           obj.init_time = init_time;
           obj.last_update_time = init_time;
           
           keps = obj.CART2KEP(pos,vel);
           
           obj.semi_major_axis     = keps(1);
           obj.eccentricity        = keps(2);
           obj.inclination         = keps(3);
           obj.argument_of_perigee = keps(4);
           obj.right_ascension     = keps(5);
           obj.mean_anomaly        = keps(6);
           obj.init_mean_anomaly   = keps(6);
            
        end
        
        function [pos, vel] = get_pos_vel(obj, at_time)
           
            n = sqrt(obj.mu / obj.semi_major_axis^3);            
            obj.mean_anomaly = (n*(at_time - obj.init_time)) +  obj.init_mean_anomaly ;                  
            
            obj.mean_anomaly = mod(obj.mean_anomaly,2*pi);
            
            [pos, vel] = obj.KEP2CART();
            
        end
        
        function maneuver(obj, DV_eci, at_time)
           
            % First get the current state
            [pos, vel] = get_pos_vel(obj, at_time);
            
            % Add DV (as impulsive shot) to current velocity
            vel_new = vel + [DV_eci(1) DV_eci(2) DV_eci(3)];
            
            % 'Insert' into new orbit, and re-set the new init_time
            insert(obj, pos, vel_new, at_time);
            
            obj.notify(['Maneuvered at: ',num2str(at_time),' seconds.']);

        end
        
        function set_mu(obj, mu_alt)
           
            % Function to set the gravitational parameter to another value
            % than the one of Earth, in case other celestial bodies are
            % orbited
            
            obj.mu_alt = mu_alt;
            
        end
        
        function [pos, vel] = KEP2CART(obj)
       
            a      = obj.semi_major_axis;
            e      = obj.eccentricity;
            i      = obj.inclination;
            argper = obj.argument_of_perigee;
            raan   = obj.right_ascension;
            M      = obj.mean_anomaly;

            % First get the eccentric anomaly:
            if e < 0.8
                Ei = M;
            else
               Ei = pi;
            end
            
            % Use Newtons method to determine the eccentic anomaly;
            Ei_1 = Ei - (Ei - e*sin(Ei) - M) / (1 - e*cos(Ei));
            while abs(Ei/Ei_1 - 1) > 1E-10
                Ei = Ei_1;
                Ei_1 = Ei - (Ei - e*sin(Ei) - M) / (1 - e*cos(Ei));
            end
            E = Ei_1;

            % Get the P and Q vectors
            % P is unit vector that points to perigee, and Q corresponds to
            % the point where the true anomaly is 90 degrees.
            
            cos_arg = cos(argper);
            cos_ra  = cos(raan);
            cos_i   = cos(i);
            
            sin_arg = sin(argper);
            sin_ra  = sin(raan);
            sin_i   = sin(i);
            
            P = [cos_arg * cos_ra - sin_arg * cos_i * sin_ra;...
                 cos_arg * sin_ra + sin_arg * cos_i * cos_ra;...
                 sin_arg * sin_i];
            
            Q = [-sin_arg * cos_ra - cos_arg * cos_i * sin_ra;...
                 -sin_arg * sin_ra + cos_arg * cos_i * cos_ra;...
                  cos_arg * sin_i];
            
            pos = a*(cos(E) - e)*P + a*sqrt(1 - e^2)*sin(E)*Q;
            r = norm(pos);
            
            vel = (sqrt(obj.mu*a)/r) * (-sin(E)*P + sqrt(1 - e^2)*cos(E)*Q);
     

%             Eold=M;
%             Enew=Eold+(M-Eold+e*sin(Eold))/(1-e*cos(Eold));
%             while abs(Eold/Enew-1)>1E-10
%                 Eold=Enew;
%                 Enew=Eold+(M-Eold+e*sin(Eold))/(1-e*cos(Eold));
%             end
%             E=Enew;
% 
%               % True anomaly [rad]
%             theta=2*atan(((1+e)/(1-e))^0.5*tan(E/2));
% 
%             % Radius to center of Earth [m]
%             r=a*(1-e*cos(E));
% 
%             % Transform to Xi-Eta-Zeta [m]
%             XiEtaZeta(1)=r*cos(theta);
%             XiEtaZeta(2)=r*sin(theta);
%             XiEtaZeta(3)=0;
% 
%             % Calculation l1,l2,m1,m2,n1,n2 [-]
%             l1=cos(raan)*cos(argper)-sin(raan)*sin(argper)*cos(i);
%             l2=-cos(raan)*sin(argper)-sin(raan)*cos(argper)*cos(i);
%             m1=sin(raan)*cos(argper)+cos(raan)*sin(argper)*cos(i);
%             m2=-sin(raan)*sin(argper)+cos(raan)*cos(argper)*cos(i);
%             n1=sin(argper)*sin(i);
%             n2=cos(argper)*sin(i);
% 
%             % Calculation XYZ ECI [m]
%             pos=([l1,l2;m1,m2;n1,n2]*[XiEtaZeta(1);XiEtaZeta(2)])';
% 
%             % Calculation Vx,Vy,Vz ECI [m]
%             H=(obj.mu*a*(1-e^2))^0.5;
%             vel(1)=obj.mu/H*(-l1*sin(theta)+l2*(e+cos(theta)));
%             vel(2)=obj.mu/H*(-m1*sin(theta)+m2*(e+cos(theta)));
%             vel(3)=obj.mu/H*(-n1*sin(theta)+n2*(e+cos(theta)));

        end
        
        function [keps] = CART2KEP(obj, pos, vel)

            % Calculations 
            r = [pos(1); pos(2); pos(3)];
            v = [vel(1); vel(2); vel(3)];

            % Get the areal velocity vector.
            h = cross(r,v);
           
            h_norm = norm(h);
            
            % Check for velocity colinear with the position:
            if h_norm == 0
                DCM = get_pointing_frame(r);
                h_n = DCM(:,1); 
            else
                h_n = h/h_norm;
            end
            
            % Get the inclination.
            i = atan2(sqrt(h_n(1)^2 + h_n(2)^2) , h_n(3) );
            
            % Get the right ascension.
            raan = atan2(h_n(1) , -h_n(2) );
            if raan < 0
                raan = 2*pi + raan;
            end
            
            % Get the semi-major axis:
            a  = 1/( (2/norm(r)) - (norm(v)^2)/obj.mu);
            
            % Mean motion:
            n = sqrt(obj.mu/a^3);
            
            % Get the eccentricity:
            % Semi latus rectum;
            p = h_norm^2 / obj.mu;
            e = sqrt(1 - p/a);
            
            % Eccentric anomaly:
            E = atan2( (r'*v)/(n*a^2) , (1 - norm(r)/a) );
            
            % Mean anomaly:
            M = E - e*sin(E);
            
            % Get the argument of perigee:
            % Argument of latitude:
            den = (-r(1)*h_n(2) + r(2)*h_n(1));
            %if den == 0
                %u = 0.0;
            %else
                u = atan2( r(3) , den );
                if u < 0
                    u = 2*pi + u;
                end
            %end
            
            % True anomaly:
            nu = atan2( sqrt(1 - e^2) * sin(E), cos(E) - e );
            if nu < 0
                nu = 2*pi + nu;
            end
            argper = u - nu;

            keps(1) = a;               
            keps(2) = e;                
            keps(3) = i;
            keps(4) = argper;
            keps(5) = raan;
            keps(6) = M;
            
        end
        
        function set_SSO(obj, LTAN_str, h_apo_km, h_peri_km, argper_deg, ma_deg, init_time, vec_to_sun_eci)
           
            % Assuming LTAN is given as string eg. '22:30'.
            
            obj.init_time = init_time;
            obj.last_update_time = init_time;
            
            obj.semi_major_axis = (h_apo_km*1000 + h_peri_km*1000 + 2*obj.Re)/2;
            obj.eccentricity = (obj.Re+h_apo_km*1000)/obj.semi_major_axis - 1;
            obj.argument_of_perigee = argper_deg*pi/180;
            
            OHMdot = 360/365;
            T = get_orbital_period(obj);
            n = 360*(24*60*60)/T;
            
            obj.inclination = acos(OHMdot / (-1.5*n*obj.J2*(obj.Re/obj.semi_major_axis)^2 * (1-obj.eccentricity^2)^-2));

            if ischar(LTAN_str)
                hr  = str2num(LTAN_str(1:2));
                min = str2num(LTAN_str(4:5));
            else
                obj.notify('Input LTAN value not valid string, using default value: 22:30');
                hr  = 22;
                min = 30;
            end
            
            % the 24:00 vector:
            v = -vec_to_sun_eci/norm(vec_to_sun_eci);
            
            coords = obj.cart2sph(v(1),v(2),v(3));
            
            theta_raan = (pi/12)*hr + (pi/720)*min;
            
            obj.right_ascension = coords(2,1) + theta_raan;
            
            
            obj.mean_anomaly        = ma_deg*pi/180;
            obj.init_mean_anomaly   = ma_deg*pi/180;
            
            obj.notify('Kepler orbit (re)configured as SSO orbit');
            
        end
        
        function T = get_orbital_period(obj)
            
            T = 2*pi*sqrt(obj.semi_major_axis^3/obj.mu);
            
        end
        
        function notify(obj, notificationStr)
        % Function that adds to notifications, which keep track of changes
        % and events in the filter object. This notifications list can be
        % later examined during analysis of a simulation.
            [~, nNotifications] = size(obj.notifications);
            noteString = [num2str(nNotifications),': ', notificationStr];
            obj.notifications{nNotifications+1} = noteString;
            % If the string of the notification begins with 'ERROR', pop up
            % a message box for the user.
            if strcmp(notificationStr(1:5), 'ERROR')
                msgbox(notificationStr);
            end
            
            if obj.VERBOSE == 1
                disp(notificationStr)
            end
        end % End of notify
        
        function showNotifications(obj)
        % Function that allows the user to use the Matlab console to 
        % quickly show all notifications in the object.
            [~, nNotifications] = size(obj.notifications);
            for i = 1:nNotifications
                disp(obj.notifications{i});
            end
        end % End of showNotifications
        
        function verbose(obj, bool)
            if bool == 1
                obj.VERBOSE = 1;
            elseif bool == 0
                obj.VERBOSE = 0;
            end
        end % End of verbose
        
        function update_orbit(obj, at_time)
           
            % Update nodal precession and argument of perigee drift:
            dt = at_time - obj.last_update_time;
            
            T = get_orbital_period(obj);
            n = 360*(24*60*60)/T;
            
            % Radians per second.
            raan_dot = -1.5*n*obj.J2*((obj.Re/obj.semi_major_axis)^2) * ((1-obj.eccentricity^2)^-2) * cos(obj.inclination)*(pi/180) /(24*60*60);
            
            delta_raan = raan_dot*dt;
            
            argper_dot = 0.75*n*obj.J2*((obj.Re/obj.semi_major_axis)^2) * (4-5*sin(obj.inclination)^2) * ((1-obj.eccentricity^2)^-2)*(pi/180) /(24*60*60);
            
            delta_argper = argper_dot*dt;
            
            obj.right_ascension = obj.right_ascension + delta_raan;
            obj.argument_of_perigee = obj.argument_of_perigee + delta_argper;
            
            obj.last_update_time = at_time;
            
        end
        
        function pre_allocate_history(obj,nr_of_steps)
           
            obj.history.nr_of_steps = nr_of_steps;
            obj.history.positions   = zeros(nr_of_steps,3);
            obj.history.velocities  = zeros(nr_of_steps,3);
            
        end
        
        function [pos, vel] = get_pos_vel_h(obj,stepNr, at_time)
        
            [pos, vel] = get_pos_vel(obj, at_time);
            obj.history.positions(stepNr,:)  = pos;
            obj.history.velocities(stepNr,:) = vel;
        
        end
 
        function keps = get_keps(obj)
            
            keps = [obj.semi_major_axis,...
                    obj.eccentricity,...
                    obj.inclination,...
                    obj.right_ascension,...
                    obj.argument_of_perigee,...
                    obj.mean_anomaly];

        end
        
        function show_orbit(obj,fig,res)
           
            figure(fig);
            
            [orb] = obj.get_orbit(res);
            
            hold on;
            
            plot3(orb(:,1),orb(:,2),orb(:,3));
            
            [xd,yd,zd] = sphere(24);
            
            xd = xd*obj.Re;
            yd = yd*obj.Re;
            zd = zd*obj.Re;
            
            hEarth = surf(xd,yd,zd, 'FaceColor','white', 'FaceAlpha',0.5);
            
            x_l = line([0,obj.Re*1.2],[0,0],[0,0], 'color', 'red');
            y_l = line([0,0],[0,obj.Re*1.2],[0,0], 'color', 'green');
            z_l = line([0,0],[0,0],[0,obj.Re*1.2], 'color', 'blue');
            
            
            axis equal
        end
        
        function [h_apo, h_peri] = get_apoperi(obj)
           
            h_peri = obj.semi_major_axis*(1-obj.eccentricity) - obj.Re;
            h_apo  = obj.semi_major_axis*(1+obj.eccentricity) - obj.Re;
            
        end
        
        function [points] = get_orbit(obj, nr_of_segments)
           
            points = zeros(nr_of_segments+1,3);
            
            ma_store = obj.mean_anomaly;
            
            d_ma = 2*pi/nr_of_segments;
            
            for i = 1:nr_of_segments+1
            
                obj.mean_anomaly = 0.0 + (i-1)*d_ma;
                [pos, ~] = obj.KEP2CART();
                points(i,:) = [pos(1) pos(2) pos(3)];
                
            end
            
            obj.mean_anomaly = ma_store;
        end
            
        
    end
    
    methods (Static)
       
        function [coordVec] = cart2sph(x,y,z)
            hypotxy = hypot(x,y);
            r = hypot(hypotxy,z);
            elev = atan2(z,hypotxy);
            az = atan2(y,x);
            coordVec = [elev;az;r];
        end
        
    end
    
end