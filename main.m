%% main.m
% @author: Agnello Hugo, Borbouse Maxime, Camus Sarah, Lucarelli Lorenzo, 
%          Miny Héloïse, Trifilò Tecla
% Date: May 2022
% Description: main code for the space experiment development project

%% Initialization
Orbit.rho = 0;
Orbit.R = 35000e3; % [km]
A6U = 2*0.3*0.2 + 2*0.1*0.2 + 2*0.3*0.1;
A3U = 2*0.3*0.1 + 2*0.3*0.1 + 2*0.3*0.1;
minRequirements = [10 -40 -30 10 -30] + 273.15; % Optic, Batt, 
maxRequirements = [30 80 65 30 65] + 273.15;
Q = 0;
Mission.tau_eclipse = 3500;
Mission.tau_sun = 8800 - 3500;
Mission.tau = 8800;
Mission.m = 12;
powRequirements = [3.5 0.01 0.06 0.65 3.22; 0 0 0.06 6.5 3.22]; 

FOV = 9;
height = 500;               % must change with hugo
%% Run codes
addpath(genpath('Mars_Landform'))
n_days = 20;

addpath(genpath('Mars_Landform'))
n_days = 20;
Orbit_parameters_polar = Orbit_func (FOV, height, 'polar', n_days)
Orbit_parameters_sunsyncnoecl = Orbit_func (FOV, height, 'sun-sync-noecl', n_days)
Orbit_parameters_sunsyncecl = Orbit_func (FOV, height, 'sun-sync-ecl', n_days)


visibility_deg = FOV/(Orbit_parameters_sunsyncecl.kep_el(1)-astroConstants(24))*astroConstants(24);

% run ('Link/telecommunication_strategy')
%only run if you have the antenna toolbox, can't commit to github :/
%% Thermal and power budget
Names = ["Polished beryllium"; "Goldized kapton (gold outside)"; ...
    "Gold"; "Aluminium tape"; "Polished aluminium"; "Aluminized kapton (aluminium outside)"; ...
    "Polished titanium"; "Black paint (epoxy)"; "Black paint (polyurethane)"; ...
    "Silver paint (electrically conducting)"; "White paint (silicone)"; ...
    "White paint (silicate)"; "Solar cells, GaAs (typical values)"; ...
    "Solar cells, silicon (typical values)"; "Aluminized kapton (kapton outside)"; ...
    "Aluminized FEP"; "Silver coated FEP (SSM)"; "OSR"];
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
Orbit.rho = 0;
Orbit.R = 0;

for i = 1:length(data(:,1))
    alpha = data(i, 1);
    epsilon = data(i, 2);
    
    [TCold, THot, alphaEpsilon] = ThermalDesign(Orbit, A, alpha, epsilon, minRequirements, maxRequirements, Q, showplot);
    
    if TCold > max(minRequirements) && THot < min(maxRequirements)
        fprintf(['<strong> Absorptivity:</strong>                      alpha = ' num2str(alpha) ' [-]\n'])
        fprintf(['<strong> Emissivity:</strong>                        epsilon = ' num2str(epsilon) ' [-]\n'])
        fprintf(['<strong> Cold case temperature:</strong>             TCold = ' num2str(TCold) ' K\n'])
        fprintf(['<strong> Hot case temperature:</strong>              THot = ' num2str(THot) ' K\n'])
    end
end