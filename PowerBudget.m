%% PowerBudget.m
% @author: Maxime Borbouse
% Date: May 2022
% Description: computes the dimensions of solar arrays.
% Inputs: * Mission: structure containing mission parameters
%         * Requirements [W]: array containing power requirements for each
%                         component
% Outputs: * Array: structure containing array parameters

function Array = PowerBudget(Requirements)
% Efficiencies
etaBCR = 0.91; % [-] Efficiency of the battery charge regulator
etaBDR = 0.89; % [-] Efficiency of the battery discharge regulator
etaCell = 0.3; % [-] Solar cell efficiency
etaPacking = 0.9; % [-] Cell packing efficiency
etaAR = 0.9; % [-] Efficiency of array regulator, assumed
etaCharge = 0.9; % [-] Charge efficiency, assumed
DOD = 0.8; % Depth of discharge (see Eurostar 3000, Li-ion)
D = 0.05; % Degradation factor

% Pointing error
deltaTheta = 1.34; % [°] Array pointing error with respect to the sun

% Radiation
S = parameters.P/(4*pi*(parameters.D*1000)^2); % [W/m^2] Solar radiation intensity
eta = etaBCR*etaBDR*etaAR;

% Power and energy
PCharge = 1/eta*max([sum(Requirements(1,:)) sum(Requirements(2,:)) sum(Requirements(3,:))])*(1/0.8 - 1);
Array.P = max([sum(Requirements(1,:)) sum(Requirements(2,:)) sum(Requirements(3,:))]) + PCharge;
Array.E_B = max([sum(Requirements(1,:)) sum(Requirements(2,:)) sum(Requirements(3,:))])*(28*24)*(1 - 0.8)/(etaCharge*DOD);
Array.epsilon = Array.P*(28*24);

% Dimension of the arrays
Array.A = Array.P/(S*cosd(deltaTheta)*etaCell*etaPacking*(1-D));
end