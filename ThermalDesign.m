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

function [TCold, THot] = ThermalDesign(Orbit, A, alpha, epsilon, Q, showplot)
%% THERMAL ENVIRONMENT
% Solar radiation
JsHot = parameters.P/(4*pi*(parameters.D*1000)^2); % [W/m^2] Solar radiation intensity
JsCold = parameters.P/(4*pi*(parameters.D*1000)^2); % [W/m^2] Solar radiation intensity

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
J_M = parameters.eps*parameters.sigma*parameters.Tb^4; % [W/m^2] Mars thermal radiation
Qir = J_M;
Rrad = parameters.R; % [m] Radius of Mars effective radiating surface
Jp = Qir*(Rrad/Orbit.R)^2;

%% THERMAL BALANCE
Asurface = A; % [m^2] Total area
Asolar = Asurface/2; % [m^2] Projected area receiving solar radiation
Aalbedo = Asolar; % [m^2] Projected area receiving albedo radiation
Aplanetary = Asolar; % [m^2] Projected area receiving planetary radiation

THot = (Aplanetary*Jp/(Asurface*parameters.sigma) + Q/(Asurface*parameters.sigma*epsilon) + ...
    (Asolar*JsHot + Aalbedo*JaHot)/(Asurface*parameters.sigma)*(alpha/epsilon))^(1/4);
TCold = (Aplanetary*Jp/(Asurface*parameters.sigma) + Q/(Asurface*parameters.sigma*epsilon) + ...
    (Asolar*JsCold + Aalbedo*JaCold)/(Asurface*parameters.sigma)*(alpha/epsilon)*0.8)^(1/4);

% Temperature with respect to absorptance and emissivity
if showplot
    alphaEpsilon = 0:0.001:1.5;
    THot_plot = (Aplanetary*Jp/(Asurface*parameters.sigma) + Q./(Asurface*parameters.sigma) + ...
        (Asolar*JsHot + Aalbedo*JaHot)/(Asurface*parameters.sigma)*(alphaEpsilon)).^(1/4); % [K] Spacecraft temperature
    TCold_plot = (Aplanetary*Jp/(Asurface*parameters.sigma) + Q./(Asurface*parameters.sigma) + ...
        (Asolar*JsHot + Aalbedo*JaHot)/(Asurface*parameters.sigma)*(alphaEpsilon)*0.8).^(1/4); % [K] Spacecraft temperature
    
    set(0,'defaultTextInterpreter','latex')

    color1 = '#FF6600';
    color2 = '#04194E';
    color3 = '#FFD700';
    color4 = '#759AAB';
    color5 = '#9E2A2B';
    
    pt = 14;
    
    figure
    hold on;
    p1 = plot(alphaEpsilon, THot_plot, 'LineWidth', 1, 'Color', color1);
    plot(alpha/epsilon, THot, 'd', 'LineWidth', 1, 'Color', color1);
    p2 = plot(alphaEpsilon, TCold_plot, 'LineWidth', 1, 'Color', color2);
    plot(alpha/epsilon, TCold, 'd', 'LineWidth', 1, 'Color', color2);
    hold off
    xlabel('$\alpha/\epsilon$ [-]', 'Interpreter', 'Latex')
    ylabel('$T$ [K]', 'Interpreter', 'Latex')
    legend([p1 p2], {'Hot case', 'Hot case'}, 'Location', 'SouthEast', 'Interpreter', 'Latex')
    set(gca, 'FontSize', pt, 'FontName', 'Times', 'LineWidth', 0.5)
end
end