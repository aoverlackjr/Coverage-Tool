classdef Satellite
    % This satellite class is used within the ITU coordination tool. The
    % class stores the orbit model and some basic RF properties
    
    properties
        Kepler; % The object that holds the orbit parameters and propagation
        SGP4;   % The object that holds the alternative orbit model
        name; % Satellite name
    end
    
    methods
        
        function sat = Satellite(name_str)
            % Set the satellites name
            sat.name = name_str;
            % Create the Kepler object
            sat.Kepler = Kepler();
        end
        
        function set_orbit_SSO(sat, LTAN_str, h_apo_km, h_peri_km, argper_deg, ma_deg, timevector)
            % One option to set the orbit as SSO, where the inclcination
            % and RAAN are calculted accoridng to theory and the given
            % LTAN.
            vec_to_sun_eci = Sun.solar_direction(timevector);
            sat.Kepler.set_SSO(LTAN_str, h_apo_km, h_peri_km, argper_deg, ma_deg, 0.0 , vec_to_sun_eci)
        end
        
        function set_orbit_keps(sat, a_km,e,i_deg,argper_deg,raan_deg, ma_deg)
            % One oprion to directly set all Kepler elements.
            sat.Kepler.set_keps(a_km,e,i_deg,argper_deg,raan_deg, ma_deg, 0.0)
        end
        
    end
    
end

