%% ThermalDesign.m
% @author: Maxime Borbouse
% Date: May 2022
% Description: computes the temperature of the satellite in both cold and
%   hot cases by performing a power balance (1 node).
% Inputs: * Orbit: structure containing orbit parameters
%         * A [m^2] : array containing surface areas (1. projected area 
%              receiving solar radiation, 2. projected area receiving 
%              albedo radiation, 3. projected area receiving planetary 
%              radiation, 4. total area)
%         * minRequirements: array containing minimum temperature
%                            requirements for each component
%         * maxRequirements: array containing maximum temperature
%                            requirements for each component
%         * showplot: boolean (1 to show the plots, 0 to hide them)
% Outputs: * TCold [K]: Satellite temperature in cold case
%          * THot [K]: Satellite temperature in hot case

function [TCold, THot] = ThermalDesign(Orbit, A, alpha, epsilon, minRequirements, maxRequirements, Q, showplot)
%% THERMAL ENVIRONMENT
% Solar radiation
JsHot = parameters.P/(4*pi*(parameters.perihelion*1000)^2); % [W/m^2] Solar radiation intensity
JsCold = parameters.P/(4*pi*(parameters.aphelion*1000)^2); % [W/m^2] Solar radiation intensity

% Albedo radiation
d = Orbit.R/parameters.R;
x = sqrt(d^2 - 1);
y = -x*cot(Orbit.rho);
if abs(Orbit.rho) < pi/2 - asin(1/d)
    F = cos(Orbit.rho)/d^2; % [-] Visibility factor
else
    F = 1/(pi*d^2)*(cos(Orbit.rho)*acos(y) - c*sin(Orbit.rho*sqrt(1-y^2)) + 1/pi*atan(sin(Orbit.rho)*sqrt(1 - y^2)/x));
end
JaHot = JsHot*parameters.albedo*F; % Intensity of albedo radiation
JaCold = JsCold*parameters.albedo*F; % Intensity of albedo radiation

% Planetary radiation
J_M = parameters.eps*parameters.sigma*parameters.Tb; % [W/m^2] Mars thermal radiation
Qir = J_M;
Rrad = parameters.R; % [m] Radius of Mars effective radiating surface
Jp = Qir*(Rrad/Orbit.R)^2;

%% THERMAL BALANCE
Asurf = A; % [m^2] Total area
Asolar = Asurf/4; % [m^2] Projected area receiving solar radiation
Aalbedo = Asolar; % [m^2] Projected area receiving albedo radiation
Aplanetary = Asolar; % [m^2] Projected area receiving planetary radiation
%CAl = 960; % [J/(kg*K)] Aluminium specific heat

THot = (Aplanetary*Jp/(Asurface*parameters.sigma) + Q/(Asurface*parameters.sigma*epsilon) + ...
    (Asolar*JsHot + Aalbedo*JaHot)/(Asurface*parameters.sigma)*(alpha/epsilon));
TCold = (Aplanetary*Jp/(Asurface*parameters.sigma) + Q/(Asurface*parameters.sigma*epsilon) + ...
    (Asolar*JsCold + Aalbedo*JaCold)/(Asurface*parameters.sigma)*(alpha/epsilon));

% Temperature with respect to absorptance and emissivity
if showplot
    alphas = 0:0.001:1;
    epsilons = 1;
    T_ae = (Aplanetary*Jp/(Asurf*parameters.sigma) + Q./(Asurf*parameters.sigma*epsilons) + ...
        (Asolar*JsHot + Aalbedo*JaHot)/(Asurf*parameters.sigma)*(alphas./epsilons)).^(1/4); % [K] Spacecraft temperature

    figure(1)
    plot(alphas, T_ae, 'LineWidth', 1);
    xlabel('$\alpha = \epsilon$ [-]', 'Interpreter', 'Latex')
    ylabel('$T$ [K]', 'Interpreter', 'Latex')
end
end