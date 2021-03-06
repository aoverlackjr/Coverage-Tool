% Simulation configuration script

% The add satellites according to their definition:
% Random joint mean anomaly, but with 90 degree spacing per orbit
ang = rand()*360;
sim.add_satellite('LPGAN', '06:00', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '06:00', 600, 600, 0, 90+ang, start_date);
sim.add_satellite('LPGAN', '06:00', 600, 600, 0, 180+ang, start_date);
sim.add_satellite('LPGAN', '06:00', 600, 600, 0, 270+ang, start_date);

% Random joint mean anomaly, but with 90 degree spacing per orbit
ang = rand()*360;
sim.add_satellite('LPGAN', '07:30', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '07:30', 600, 600, 0, 90+ang, start_date);
sim.add_satellite('LPGAN', '07:30', 600, 600, 0, 180+ang, start_date);
sim.add_satellite('LPGAN', '07:30', 600, 600, 0, 270+ang, start_date);

% Random joint mean anomaly, but with 90 degree spacing per orbit
ang = rand()*360;
sim.add_satellite('LPGAN', '09:00', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '09:00', 600, 600, 0, 90+ang, start_date);
sim.add_satellite('LPGAN', '09:00', 600, 600, 0, 180+ang, start_date);
sim.add_satellite('LPGAN', '09:00', 600, 600, 0, 270+ang, start_date);

% Random joint mean anomaly, but with 90 degree spacing per orbit
ang = rand()*360;
sim.add_satellite('LPGAN', '10:30', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '10:30', 600, 600, 0, 90+ang, start_date);
sim.add_satellite('LPGAN', '10:30', 600, 600, 0, 180+ang, start_date);
sim.add_satellite('LPGAN', '10:30', 600, 600, 0, 270+ang, start_date);

% Random joint mean anomaly, but with 90 degree spacing per orbit
ang = rand()*360;
sim.add_satellite('LPGAN', '12:00', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '12:00', 600, 600, 0, 90+ang, start_date);
sim.add_satellite('LPGAN', '12:00', 600, 600, 0, 180+ang, start_date);
sim.add_satellite('LPGAN', '12:00', 600, 600, 0, 270+ang, start_date);

% Random joint mean anomaly, but with 90 degree spacing per orbit
ang = rand()*360;
sim.add_satellite('LPGAN', '13:30', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '13:30', 600, 600, 0, 90+ang, start_date);
sim.add_satellite('LPGAN', '13:30', 600, 600, 0, 180+ang, start_date);
sim.add_satellite('LPGAN', '13:30', 600, 600, 0, 270+ang, start_date);

% Random joint mean anomaly, but with 90 degree spacing per orbit
ang = rand()*360;
sim.add_satellite('LPGAN', '15:00', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '15:00', 600, 600, 0, 90+ang, start_date);
sim.add_satellite('LPGAN', '15:00', 600, 600, 0, 180+ang, start_date);
sim.add_satellite('LPGAN', '15:00', 600, 600, 0, 270+ang, start_date);

% Random joint mean anomaly, but with 90 degree spacing per orbit
ang = rand()*360;
sim.add_satellite('LPGAN', '16:30', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '16:30', 600, 600, 0, 90+ang, start_date);
sim.add_satellite('LPGAN', '16:30', 600, 600, 0, 180+ang, start_date);
sim.add_satellite('LPGAN', '16:30', 600, 600, 0, 270+ang, start_date);