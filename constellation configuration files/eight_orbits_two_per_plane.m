% Simulation configuration script

% The add satellites according to their definition:
ang = rand()*360;
sim.add_satellite('LPGAN', '06:00', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '06:00', 600, 600, 0, 180+ang, start_date);

ang = rand()*360;
sim.add_satellite('LPGAN', '07:30', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '07:30', 600, 600, 0, 180+ang, start_date);

ang = rand()*360;
sim.add_satellite('LPGAN', '09:00', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '09:00', 600, 600, 0, 180+ang, start_date);

ang = rand()*360;
sim.add_satellite('LPGAN', '10:30', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '10:30', 600, 600, 0, 180+ang, start_date);

ang = rand()*360;
sim.add_satellite('LPGAN', '12:00', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '12:00', 600, 600, 0, 180+ang, start_date);

ang = rand()*360;
sim.add_satellite('LPGAN', '13:30', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '13:30', 600, 600, 0, 180+ang, start_date);

ang = rand()*360;
sim.add_satellite('LPGAN', '15:00', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '15:00', 600, 600, 0, 180+ang, start_date);

ang = rand()*360;
sim.add_satellite('LPGAN', '16:30', 600, 600, 0, 0+ang, start_date);
sim.add_satellite('LPGAN', '16:30', 600, 600, 0, 180+ang, start_date);