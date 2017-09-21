addpath(genpath(pwd));
clc;
clear;

% Example of running the contact time tool:
sim = contact_time_sim();

% Define the start date of the simulation:
start_date = [2016 9 6 12 0 0];
% Define the stop date of the simulation:
stop_date  = [2016 9 8 12 0 0];
% Set the simulation iteration time step (in seconds):
dt = 10;

% Set radius of Earth in meters:
Re = 6378135;

% Run sat-vs-orbit configuration script of choice
eight_orbits_two_per_plane;
             
% Obtain a mesh spread over earth:
% Number of resolution steps of Earth-covering mesh;
res = 36;
[x,y,z] = sphere(res);
sim.add_mesh(-x*Re, -y*Re, -z*Re);

% Show what the constellation looks like in figure 1:
sim.show_configuration(1);

% The minumum elevation needed for a ground point to be 'in contact' with a
% satellite
minimum_elevation = 18;

% Run the simulation
sim.run(start_date, stop_date, dt, minimum_elevation);

% Process the data:
max_contact_gap = 3600;
sim.process(max_contact_gap);

% Obtain an image of the maximum contact time gaps:
figure(2);
%img = imread('earth.jpg');
img = imread('bw_lines.jpg');
imagesc([1 res-1], [1 res-1], flipud(img)); %flipdim(img,1);
hold on;
surf(sim.output_data.contact_gaps_binary, 'EdgeColor', 'none', 'FaceAlpha', 0.3 );
colorbar;
shading interp; 
set(gca,'ydir','normal');

