classdef Sun < handle
    
    properties (Constant)
        radius              = 1392000000;
        OHM_plus_ohm        = 282.9400;
        obliquityOfEcliptic = 23.43929111;
        earthAngularRate    = 0.25068447*pi/180/60;
    end
    
    methods
        
        function sun = Sun()

        end
    end
    
    methods (Static)

        function [sol_eci] = solar_direction(timevector)
            [~,T] = CTS.juliandate(timevector);
            M = 357.5256 + 35999.049 * T;
            M_rad = M * pi/180;
            lambda = (Sun.OHM_plus_ohm + M + (6892 / 3600)*sin(M_rad) + 0.02 * sin(2*M_rad)) * pi/180;
            sol_eci = [cos(lambda); sin(lambda) * cos(Sun.obliquityOfEcliptic * pi/180); sin(lambda)* sin(Sun.obliquityOfEcliptic * pi/180)];
        end

    end

end

