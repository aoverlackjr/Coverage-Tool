classdef contact_time_sim < handle
    % Simulation to calculate contact times and gaps between contacts
    % between a set of given satellites and ground points.
    
    properties
        sat_list = {}       % List holding satellite orbit objects.
        gp_list  = []       % Array holding points on the earth.
        dt       = 10       % Simulation timestep of 10 seconds (default)
        start_date;         % A starting date of the simulation.
        minimum_elevation;  % Single number with minimum elevation for making contact between ground and sat.
        stop_date;          % A stopping date of the simulation.
        output_data;        % Data structure to hold the simulation results.
        nr_of_steps;        % Size of simulation data
        nr_of_gps;          % Amount of points used over the Earth
        mesh;               % Structure with ground mesh
        gp_option = 'none'
    end
    
    methods
        
        function sim = contact_time_sim()
            % Empty constructor. Configuration is done with below functions.   
        end
        
        function add_satellite(sim, name_str, LTAN_str, h_apo_km, h_peri_km, argper_deg, ma_deg, timevector)
            % Create another satellite and append to the list:
            [~,n] = size(sim.sat_list);
            ind = n + 1;
            sim.sat_list{ind} = Satellite(name_str);
            % Set the orbit according to the SSO params given.                    
            sim.sat_list{ind}.set_orbit_SSO(LTAN_str, h_apo_km, h_peri_km, argper_deg, ma_deg, timevector);
        end
        
        function add_satellite_kepler(sim, name_str, a_km, e, i_deg, argper_deg, raan_deg, ma_deg)
            % Create another satellite and append to the list:
            [~,n] = size(sim.sat_list);
            ind = n + 1;
            sim.sat_list{ind} = Satellite(name_str);
            % Set the orbit according to the Kepler params given. 
            sim.sat_list{ind}.set_orbit_keps(a_km, e, i_deg, argper_deg, raan_deg, ma_deg)
        end
        
        function add_ground_points(sim, array_with_points)
            % Simple copy function to input a nx3 list of Earth-centered-Earth-Fixed
            % cartesian coordinates.
            sim.gp_list = array_with_points;
            sim.gp_option = 'points';
        end
        
        function add_mesh(sim, mesh_x, mesh_y, mesh_z)
            sim.mesh.x = mesh_x;
            sim.mesh.y = mesh_y;
            sim.mesh.z = mesh_z;
            [n,m] = size(mesh_x);
            sim.mesh.n = n;
            sim.mesh.m = m;
            sim.gp_option = 'mesh';
        end
        
        function run(sim, start_date, stop_date, dt_seconds, min_elev_deg)
           
            % Calculate how many calculation steps are required.
            sim.start_date = start_date;
            sim.stop_date  = stop_date;
            nr_of_days = datenum(stop_date) - datenum(start_date);
            sim.dt = dt_seconds;
            n_steps = floor(nr_of_days*24*60*60/dt_seconds);
            sim.nr_of_steps = n_steps;
            % Set the minimum elevation for the sim
            sim.minimum_elevation = min_elev_deg*pi/180;
            
            if strcmp(sim.gp_option, 'points')
                n_sats = length(sim.sat_list);
                m_gps  = length(sim.gp_list);
                sim.nr_of_gps = m_gps;
                if (n_sats < 1) || (m_gps < 1)
                    disp('Cannot run sim. Define at least one satellite and one ground point.');
                else

                    % The variable where the satellites' instantaneous
                    % positions are placed.
                    sat_positions = zeros(3,n_sats);

                    % Pre-allocate output data.
                    % The output data consists of rows per each timestep, and
                    % each column indicates (with a 1 or 0) whether a ground point (associated
                    % with the column index) is 'seen' by a satellite.
                    sim.output_data.contact_data = zeros(n_steps,m_gps);

                    % Start the loop of all timesteps between the start and
                    % stop date.
                    disp('Starting simulation loop.');
                    for stepNr = 1:n_steps
                        % Calculate the time parameters; the time passed in
                        % seconds, and the actual date for frame
                        % transformations.
                        time_seconds = stepNr*dt_seconds;
                        time_vector  = datevec(datenum(start_date + [0 0 0 0 0 time_seconds]));
                        % Rotation of the Earth
                        ROT = CTS.R3(sim.JD2GAST(sim.juliandate(time_vector))*pi/180)';
                        % Propagate all sats for the time step.
                        for satNr = 1:n_sats
                            [pos,~] = sim.sat_list{satNr}.Kepler.get_pos_vel(time_seconds);
                            sat_positions(:,satNr) = pos';
                        end
                        % Find all GS positions in the inertial frame. and
                        % check their connections to the sats.
                        for gpNr = 1:m_gps
                            % get ground point position in inertial frame.
                            gp_position = ROT * sim.gp_list(gpNr,:)';
                            % Start looping satellites to see if any is
                            % in range of the GS at hand.
                            for satNr = 1:n_sats
                                % Check if the sat is in range:
                                [in_range, ~, ~] = sim.is_in_range(gp_position, sat_positions(:,satNr));
                                % If the satellite is in range, put in contact
                                % data:
                                if in_range
                                    % Mark score
                                    sim.output_data.contact_data(stepNr,gpNr) = 1;
                                end
                            end
                        end
                    end
                end
            elseif strcmp(sim.gp_option, 'mesh')
                
                n_sats = length(sim.sat_list);
                m_gps  = sim.mesh.n*sim.mesh.m;
                sim.nr_of_gps = m_gps;
                if (n_sats < 1) || (m_gps < 1)
                    disp('Cannot run sim. Define at least one satellite and one ground point.');
                else

                    % The variable where the satellites' instantaneous
                    % positions are placed.
                    sat_positions = zeros(3,n_sats);

                    % Pre-allocate output data.
                    % The output data consists of rows per each timestep, and
                    % each column indicates (with a 1 or 0) whether a ground point (associated
                    % with the column index) is 'seen' by a satellite.
                    sim.output_data.contact_data = zeros(n_steps,sim.mesh.n,sim.mesh.m);

                    % Start the loop of all timesteps between the start and
                    % stop date.
                    disp('Starting simulation loop.');
                    for stepNr = 1:n_steps
                        % Calculate the time parameters; the time passed in
                        % seconds, and the actual date for frame
                        % transformations.
                        time_seconds = stepNr*dt_seconds;
                        time_vector  = datevec(datenum(start_date + [0 0 0 0 0 time_seconds]));
                        % Rotation of the Earth
                        ROT = CTS.R3(sim.JD2GAST(sim.juliandate(time_vector))*pi/180)';
                        % Propagate all sats for the time step.
                        for satNr = 1:n_sats
                            [pos,~] = sim.sat_list{satNr}.Kepler.get_pos_vel(time_seconds);
                            sat_positions(:,satNr) = pos';
                        end
                        % Find all GS positions in the inertial frame. and
                        % check their connections to the sats.
                        for gpNr_n = 1:sim.mesh.n
                            for gpNr_m = 1:sim.mesh.m
                                % get ground point position in inertial frame.
                                vec = [sim.mesh.x(gpNr_n, gpNr_m); sim.mesh.y(gpNr_n, gpNr_m); sim.mesh.z(gpNr_n, gpNr_m)];
                                gp_position = ROT * vec;
                                % Start looping satellites to see if any is
                                % in range of the GS at hand.
                                for satNr = 1:n_sats
                                    % Check if the sat is in range:
                                    [in_range, ~, ~] = sim.is_in_range(gp_position, sat_positions(:,satNr));
                                    % If the satellite is in range, put in contact
                                    % data:
                                    if in_range
                                        % Mark score
                                        sim.output_data.contact_data(stepNr,gpNr_n, gpNr_m) = 1;
                                    end
                                end
                            end
                        end
                    end
                end

                disp('Simulation loop finished.');
            end
            
        end
        
        function process(sim, max_contact_gap)
            
            if strcmp(sim.gp_option, 'points')
                % Process some statistics for the ground contacts:
                sim.output_data.contact_percentages = zeros(1, sim.nr_of_gps);
                sim.output_data.contact_gaps        = zeros(1, sim.nr_of_gps);
                sim.output_data.contact_gaps_binary = zeros(1, sim.nr_of_gps);
                
                % Get percentage of time 'seen' by a satellite per point.
                for pt = 1:sim.nr_of_gps
                    sim.output_data.contact_percentages(pt) = sum(sim.output_data.contact_data(:,pt))/sim.nr_of_steps;
                end
                % Get total mean contact time percentage
                sim.output_data.average_contact_percentage = sum(sim.output_data.contact_percentages)/sim.nr_of_gps;

                % Find largest amount of time that is uncontacted for an point:
                current_void = 0;
                largest_void = 0;
                for pt = 1:sim.nr_of_gps
                   for stepNr = 1:sim.nr_of_steps
                        % If a data point is zero add to the total gap
                        if sim.output_data.contact_data(stepNr, pt) == 0
                            current_void = current_void + 1;
                        else % The data point is not zero, check if this was the biggest gap
                            if current_void > largest_void
                                largest_void = current_void;
                            end
                            current_void = 0;
                        end
                   end
                   % A new point is being checked so place largest gap in data and zero the counter:
                   sim.output_data.contact_gaps(pt) = largest_void;
                   current_void = 0;
                   largest_void = 0;
                end

                % The largest gap in time experienced by a point on the ground;
                sim.output_data.contact_gaps = sim.output_data.contact_gaps*sim.dt;
                sim.output_data.largest_void = max(sim.output_data.contact_gaps);
                for pt = 1:sim.nr_of_gps
                    if sim.output_data.contact_gaps(pt) > max_contact_gap
                           sim.output_data.contact_gaps_binary(pt) = 1;
                    end
                end
                
            elseif strcmp(sim.gp_option, 'mesh')
               
                % Process some statistics for the ground contacts:
                sim.output_data.contact_percentages = zeros(sim.mesh.n, sim.mesh.m);
                sim.output_data.contact_gaps        = zeros(sim.mesh.n, sim.mesh.m);
                sim.output_data.contact_gaps_binary = zeros(sim.mesh.n, sim.mesh.m);

                % Get percentage of time 'seen' by a satellite per point.
                for ptn = 1:sim.mesh.n
                    for ptm = 1:sim.mesh.m
                        sim.output_data.contact_percentages(ptn,ptm) = sum(sim.output_data.contact_data(:,ptn,ptm))/sim.nr_of_steps;
                    end
                end
                % Get total mean contact time percentage
                sim.output_data.average_contact_percentage = sum(sim.output_data.contact_percentages)/sim.nr_of_gps;

                % Find largest amount of time that is uncontacted for an point:
                current_void = 0;
                largest_void = 0;
                for ptn = 1:sim.mesh.n
                    for ptm = 1:sim.mesh.m
                       for stepNr = 1:sim.nr_of_steps
                            % If a data point is zero add to the total gap
                            if sim.output_data.contact_data(stepNr, ptn, ptm) == 0
                                current_void = current_void + 1;
                            else % The data point is not zero, check if this was the biggest gap
                                if current_void > largest_void
                                    largest_void = current_void;
                                end
                                current_void = 0;
                            end
                       end
                       % A new point is being checked so place largest gap in data and zero the counter:
                       sim.output_data.contact_gaps(ptn, ptm) = largest_void;
                       current_void = 0;
                       largest_void = 0;
                    end
                end

                % The largest gap in time experienced by a point on the ground;
                sim.output_data.contact_gaps = sim.output_data.contact_gaps*sim.dt;
                sim.output_data.largest_void = max(max(sim.output_data.contact_gaps));
                for ptn = 1:sim.mesh.n
                    for ptm = 1:sim.mesh.m
                        if sim.output_data.contact_gaps(ptn, ptm) > max_contact_gap
                           sim.output_data.contact_gaps_binary(ptn, ptm) = 1;
                        end
                    end
                end

            end
        end
        
        function [in_range, slant_vector, dist] = is_in_range(sim, vec_to_gs_eci, vec_to_sat_eci)
            slant_vector = vec_to_sat_eci - vec_to_gs_eci;
            dist = norm(slant_vector);
            s_n = slant_vector/dist;
            gs_eci_dir = vec_to_gs_eci/norm(vec_to_gs_eci);
            angle = acos(gs_eci_dir'*s_n);
            if angle < (0.5*pi - sim.minimum_elevation);
                in_range = 1;
            else
                in_range = 0;
            end
        end
        
        function show_configuration(sim, hFig)
           
            figure(hFig);
            if strcmp(sim.gp_option, 'points')
                plot_points(hFig, sim.gp_list);
            end
            
            n_sats = length(sim.sat_list);
            
            for satNr = 1:n_sats
                sim.sat_list{satNr}.Kepler.show_orbit(hFig, 60)
            end
            
        end
       
        function jd = juliandate(sim, varargin ) 

            [year month day hour min sec] = datevec(datenum(varargin{:}));

            for k = length(month):-1:1
                if ( month(k) <= 2 ) % january & february
                    year(k)  = year(k) - 1.0;
                    month(k) = month(k) + 12.0;
                end
            end

            jd = floor( 365.25*(year + 4716.0)) + floor( 30.6001*( month + 1.0)) + 2.0 - ...
                floor( year/100.0 ) + floor( floor( year/100.0 )/4.0 ) + day - 1524.5 + ...
                (hour + min/60 + sec/3600)/24;
        end
        
        function GAST = JD2GAST(sim, JD)
            %THETAm is the mean siderial time in degrees
            THETAm = sim.JD2GMST(JD);

            %Compute the number of centuries since J2000
            T = (JD - 2451545.0)./36525;

            %Mean obliquity of the ecliptic (EPSILONm)
            % see http://www.cdeagle.com/ccnum/pdf/demogast.pdf equation 3
            % also see Vallado, Fundamentals of Astrodynamics and Applications, second edition.
            %pg. 214 EQ 3-53
            EPSILONm = 23.439291-0.0130111.*T - 1.64E-07.*(T.^2) + 5.04E-07.*(T.^3);

            %Nutations in obliquity and longitude (degrees)
            % see http://www.cdeagle.com/ccnum/pdf/demogast.pdf equation 4
            L = 280.4665 + 36000.7698.*T;
            dL = 218.3165 + 481267.8813.*T;
            OMEGA = 125.04452 - 1934.136261.*T;

            %Calculate nutations using the following two equations:
            % see http://www.cdeagle.com/ccnum/pdf/demogast.pdf equation 5
            dPSI = -17.20.*sind(OMEGA) - 1.32.*sind(2.*L) - .23.*sind(2.*dL) ...
                + .21.*sind(2.*OMEGA);
            dEPSILON = 9.20.*cosd(OMEGA) + .57.*cosd(2.*L) + .10.*cosd(2.*dL) - ...
                .09.*cosd(2.*OMEGA);

            %Convert the units from arc-seconds to degrees
            dPSI = dPSI.*(1/3600);
            dEPSILON = dEPSILON.*(1/3600);

            %(GAST) Greenwhich apparent sidereal time expression in degrees
            % see http://www.cdeagle.com/ccnum/pdf/demogast.pdf equation 1
            GAST = mod(THETAm + dPSI.*cosd(EPSILONm+dEPSILON),360);
        end
        
        function GMST = JD2GMST(sim, JD)
            %Find the Julian Date of the previous midnight, JD0
            JD0 = NaN(size(JD));
            JDmin = floor(JD)-.5;
            JDmax = floor(JD)+.5;
            JD0(JD > JDmin) = JDmin(JD > JDmin);
            JD0(JD > JDmax) = JDmax(JD > JDmax);
            H = (JD-JD0).*24;       %Time in hours past previous midnight
            D = JD - 2451545.0;     %Compute the number of days since J2000
            D0 = JD0 - 2451545.0;   %Compute the number of days since J2000
            T = D./36525;           %Compute the number of centuries since J2000
            %Calculate GMST in hours (0h to 24h) ... then convert to degrees
            GMST = mod(6.697374558 + 0.06570982441908.*D0  + 1.00273790935.*H + ...
            0.000026.*(T.^2),24).*15;
        end
        
    end
    
end

