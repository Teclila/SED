%% main.m
% @author: Agnello Hugo, Borbouse Maxime, Camus Sarah, Lucarelli Lorenzo, 
%          Miny Héloïse, Trifilò Tecla
% Date: May 2022
% Description: main code for the space experiment development project

%% Initialization
Orbit.rho = 0;
Orbit.R = 35000e3; % [km]
A8U = 2*0.1*0.125 + 2*0.1*0.64 + 2*0.125*0.64;
A8Ubis = 2*0.43*0.1 + 2*0.43*0.186 + 2*0.1*0.186; 
A6U = 2*0.3*0.2 + 2*0.1*0.2 + 2*0.3*0.1;
A3U = 2*0.3*0.1 + 2*0.3*0.1 + 2*0.3*0.1;
% Batteries, Solar panels, Optics, ADCS, Antenna
minRequirements = [-20 -40 0 -20 -40] + 273.15;
maxRequirements = [45 80 30 70 85] + 273.15;
Q = 0;
Mission.tau_eclipse = 3500;
Mission.tau_sun = 8800 - 3500;
Mission.tau = 8800;
Mission.m = 12;
% ADCS devices, Optical system, OBDH, Transmitter-receiver
powRequirements = [10.85 9 0.88 0.65; 10.85 0 0.88 5; 0 0 0.88 0.65];  

FOV = 9;
height = 500;               % must change with hugo
%% Run codes
addpath(genpath('Mars_Landform'))
n_days = 20;

addpath(genpath('Mars_Landform'))
n_days = 20;
Orbit_parameters_polar = Orbit_func (FOV, height, 'polar', n_days);
Orbit_parameters_sunsyncnoecl = Orbit_func (FOV, height, 'sun-sync-noecl', n_days);
Orbit_parameters_sunsyncecl = Orbit_func (FOV, height, 'sun-sync-ecl', n_days);


visibility_deg = FOV/(Orbit_parameters_sunsyncecl.kep_el(1)-astroConstants(24))*astroConstants(24);

% run ('Link/telecommunication_strategy')
%only run if you have the antenna toolbox, can't commit to github :/

%% Thermal and power budget
Orbit.rho = 0;
Orbit.R = parameters.Rpolar + 500; % [km]
Names = {'Polished beryllium'; 'Goldized kapton (gold outside)'; ...
    'Gold'; 'Aluminium tape'; 'Polished aluminium'; 'Aluminized kapton (aluminium outside)'; ...
    'Polished titanium'; 'Black paint (epoxy)'; 'Black paint (polyurethane)'; ...
    'Silver paint (electrically conducting)'; 'White paint (silicone)'; ...
    'White paint (silicate)'; 'Solar cells, GaAs (typical values)'; ...
    'Solar cells, silicon (typical values)'; 'Aluminized kapton (kapton outside)'; ...
    'Aluminized FEP'; 'Silver coated FEP (SSM)'; 'OSR'};
data = [0.44 0.01 44.00; ...
        0.25 0.02 12.5; ...
        0.25 0.04 6.25; ...
        0.21 0.04 5.25; ...
        0.24 0.08 3.00; ...
        0.14 0.05 2.80; ...
        0.60 0.60 1.00; ...
        0.95 0.85 1.12; ...
        0.95 0.90 1.06; ...
        0.37 0.44 0.84; ...
        0.26 0.83 0.31; ...
        0.12 0.90 0.13; ...
        0.88 0.80 1.10; ...
        0.75 0.82 0.91; ...
        0.40 0.63 0.63; ...
        0.16 0.47 0.34; ...
        0.08 0.78 0.10; ...
        0.07 0.74 0.09];
T = zeros(length(data(:,1)), 2);
TRange = [max(minRequirements) min(maxRequirements)];
for i = 1:length(data(:,1))
    alpha = data(i, 1);
    epsilon = data(i, 2);
    
    [TCold, THot] = ThermalDesign(Orbit, A8Ubis, alpha, epsilon, Q, 0);
    T(i, 2) = THot;
    T(i, 1) = TCold;
    if TCold > max(minRequirements) && THot < min(maxRequirements)
        [TCold, THot] = ThermalDesign(Orbit, A8Ubis, alpha, epsilon, Q, 1);
        title(Names{i}, 'Interpreter', 'Latex')
        fprintf('<strong>- Thermal design - </strong>\n')
        fprintf(['<strong>Coating: </strong>                          ' Names{i} '\n'])
        fprintf(['<strong>Absorptivity:</strong>                      alpha = ' num2str(alpha) ' [-]\n'])
        fprintf(['<strong>Emissivity:</strong>                        epsilon = ' num2str(epsilon) ' [-]\n'])
        fprintf(['<strong>Cold case temperature:</strong>             TCold = ' num2str(TCold) ' K\n'])
        fprintf(['<strong>Hot case temperature:</strong>              THot = ' num2str(THot) ' K\n\n'])
    end
end

for i = 1:length(data(:,1))
    fprintf([Names{i} '  ' num2str(data(i,1)) '  ' num2str(data(i,2)) '  ' num2str(T(i,1)-273.15) '  ' num2str(T(i,2)-273.15) '\n'])
end

Array = PowerBudget(powRequirements);
fprintf('<strong>- Power budget - </strong>\n')
fprintf(['<strong>Needed power: </strong>                    ' num2str(Array.P) ' W\n'])
fprintf(['<strong>Battery stored energy: </strong>           ' num2str(Array.E_B) ' W-hrs\n'])
fprintf(['<strong>Total energy required from array: </strong>' num2str(Array.epsilon) ' W-hrs\n'])
fprintf(['<strong>Array surface area: </strong>              ' num2str(Array.A) ' m²\n'])