%% main.m
% @author: Agnello Hugo, Borbouse Maxime, Camus Sarah, Lucarelli Lorenzo, 
%          Miny Héloïse, Trifilò Tecla
% Date: May 2022
% Description: main code for the space experiment development project

%% Initialization
Orbit.rho = 0;
Orbit.R = 35000e3; % [km]
A = [0.3*0.2 0.3*0.2 0.3*0.2 0.3*0.2];
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
Orbit_parameters_polar = Orbit_func (FOV, height, 'polar', n_days)
Orbit_parameters_sunsyncnoecl = Orbit_func (FOV, height, 'sun-sync-noecl', n_days)
Orbit_parameters_sunsyncecl = Orbit_func (FOV, height, 'sun-sync-ecl', n_days)

visibility_deg = FOV/(Orbit_parameters_sunsyncecl.kep_el(1)-astroConstants(24))*astroConstants(24);

%%
[TCold, THot, alphaEpsilon] = ThermalDesign(Orbit, Mission, A, minRequirements, maxRequirements, Q, 1);
Array = PowerBudget(Mission, powRequirements);