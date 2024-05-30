close all;clc;clear;
addpath('module')
Directional_coupler = ChenYi_BPM();                        % Initial the Object parameter
Directional_coupler.xmin = -1.5;
Directional_coupler.xmax = 1.5;
Directional_coupler.xmesh = 0.01;                         % height : mesh size (um)
Directional_coupler.zmesh = 0.01;                         % height : mesh size (um)
Directional_coupler.zsample = 20000;                       % zmesh : Number of samples
Directional_coupler.Waveguide.gap = 0.25;                  % Gap of two waveguides (um)
Directional_coupler.Waveguide.Width = 0.45;                % Width of waveguide (um)
Directional_coupler.Waveguide.cladding_index = 1.444;      % Refractive index of Silicon dioxide
Directional_coupler.Waveguide.core_index = 3.4777;         % Refractive index of Silicon 
Directional_coupler.Source.postion = 0.45/2+0.25/2;        % Position of gaussian source (um)
Directional_coupler.Source.bandwidth = 0.05;                % Bandwidth of gaussian source (um)
Directional_coupler = Directional_coupler.StructureView(); % Show the waveguide & source
Directional_coupler = Directional_coupler.BPM_Run();       % Run the BPM
Directional_coupler = Directional_coupler.Field_SideView();% Show the Intensity distribion(side view)
Directional_coupler = Directional_coupler.Field_topView(); % Show the Intensity distribion(top view)
Directional_coupler = Directional_coupler.Coupling_efficeincy();% Show the Coupling efficeincy