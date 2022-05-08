%% PowerBudget.m
% @author: Maxime Borbouse
% Date: May 2022
% Description: computes the dimensions of solar arrays.
% Inputs: * Mission: structure containing mission parameters
%         * Requirements [W]: array containing power requirements for each
%                         component
% Outputs: * Array: structure containing array parameters

function Array = PowerBudget(Mission, Requirements)
%%%%%%%%%%%%%%%%%%%%%%%%%%%% To be updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
etaAR = 0.9; % [-] Efficiency of array regulator
etaCharge = 0.9; % [-] Charge efficiency
DOD = 0.1; % Depth of discharge
D = 0.1; % Degradation factor (see 10.3.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Efficiencies
etaBCR = 0.91; % [-] Efficiency of the battery charge regulator
etaBDR = 0.89; % [-] Efficiency of the battery discharge regulator
etaCell = 0.9; % [-] Solar cell efficiency
etaPacking = 0.9; % [-] Cell packing efficiency

% Pointing error
deltaTheta = 1; % [°] Array pointing error with respect to the sun

% Radiation
S = parameters.P/(4*pi*(parameters.D*1000)^2); % [W/m^2] Solar radiation intensity
eta = etaBCR*etaBDR*etaAR;

% Power and energy
PCharge = 1/eta*sum(Requirements(:,1))*Mission.tau_eclipse/Mission.tau_sun;
Array.P = sum(Requirements(:,2)) + PCharge;
Array.E_B = sum(Requirements(:,1))*(Mission.tau - Mission.tau_sun)/(etaCharge*DOD);
Array.epsilon = Array.P*Mission.tau;

% Dimension of the arrays
Array.A = Array.P/(S*cosd(deltaTheta)*etaCell*etaPacking*(1-D));
end