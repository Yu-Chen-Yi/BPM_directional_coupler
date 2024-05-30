# BPM_directional_coupler
beam propagation method directional coupler waveguide
# Directional Coupler Simulation

This code simulates the behavior of a directional coupler using the Beam Propagation Method (BPM).

## Prerequisites

Ensure you have the required module in your working directory:
- `module`

## Usage

#### 1. Initialize the simulation environment:

```matlab
close all; clc; clear;
addpath('module');
```
#### 2. Create and configure the directional coupler object:
```matlab
Directional_coupler = ChenYi_BPM();                    % Initialize the object parameter
Directional_coupler.xmin = -1.5;
Directional_coupler.xmax = 1.5;
Directional_coupler.xmesh = 0.01;                      % Mesh size in the x-direction (um)
Directional_coupler.zmesh = 0.01;                      % Mesh size in the z-direction (um)
Directional_coupler.zsample = 20000;                   % Number of samples in the z-direction
Directional_coupler.Waveguide.gap = 0.25;              % Gap between the two waveguides (um)
Directional_coupler.Waveguide.Width = 0.45;            % Width of the waveguide (um)
Directional_coupler.Waveguide.cladding_index = 1.444;  % Refractive index of silicon dioxide
Directional_coupler.Waveguide.core_index = 3.4777;     % Refractive index of silicon
Directional_coupler.Source.position = 0.45/2 + 0.25/2; % Position of the Gaussian source (um)
Directional_coupler.Source.bandwidth = 0.05;           % Bandwidth of the Gaussian source (um)

```
#### 3. Visualize the structure:
```matlab
Directional_coupler = Directional_coupler.StructureView(); % Display the waveguide and source
```
#### 4. Run the BPM simulation:
```matlab
Directional_coupler = Directional_coupler.BPM_Run();  % Execute the BPM
```
#### 5. Visualize the intensity distribution:
```matlab
Directional_coupler = Directional_coupler.Field_SideView(); % Display side view of intensity distribution
Directional_coupler = Directional_coupler.Field_topView();  % Display top view of intensity distribution
```
#### 6. Calculate and display the coupling efficiency:
```matlab
Directional_coupler = Directional_coupler.Coupling_efficiency(); % Display the coupling efficiency
```
## Demonstrate
![image](https://github.com/Yu-Chen-Yi/BPM_directional_coupler/assets/64921305/412bc71d-315d-4469-a412-bd83da9c0fef)
![image](https://github.com/Yu-Chen-Yi/BPM_directional_coupler/assets/64921305/d33f73ec-b49f-48fe-a794-63d382f7d70a)
![image](https://github.com/Yu-Chen-Yi/BPM_directional_coupler/assets/64921305/534b1759-7d5f-4fc6-b1ff-f4c0d2104d84)
![image](https://github.com/Yu-Chen-Yi/BPM_directional_coupler/assets/64921305/e14020dd-7866-4fb4-b1b6-1ed1dff64833)


