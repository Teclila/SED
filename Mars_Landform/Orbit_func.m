function Orbit = Orbit_func (FOV, height, orbittype, n_d)
addpath(genpath('../functions'))


% 1- NOMINAL ORBIT DATA 

R_mars = astroConstants(24); %Radius of Mars
mu=astroConstants(14);


h0 = height;            %altitude [km] must be set as function of GSD
a0= h0 + R_mars;        %semi major axis [km]
e0= 0;                  %eccentricity
om0=deg2rad(0);         %choosen perigee anomaly [rad]
th0=0;                  %choosen initial real anomaly [rad]
lon_0 = 0;              %Initial longitude

% A_M=0.03;               %area to mass ratio [m^2/kg]
% 
% Cd=0.5;                 %drag coefficient
% Mars parameters
j2 = 0.00196045;

run('AssignConstants_Orbiter_mro')
a_mro = myConst.a_mro;
e_mro = myConst.e_mro;
i_mro = myConst.i_mro;
OM_mro = myConst.OM_mro;
om_mro = myConst.om_mro;
th_mro = myConst.trueAnom_mro;
greenwich0_mro = myConst.greenwich0_mro;


switch orbittype
    case 'polar'
OM0=deg2rad(0);         %choosen RAAN [rad]    
i0=deg2rad(90);         %inclination [rad]
Orbit.kep_el = [a0 e0 i0 OM0 om0 th0];
FOV = deg2rad(FOV);
p = a0*(1-e0^2);
Orbit.DELTA_OM = -3*pi*j2*R_mars^2/p^2*cos(i0);

% 2- GROUND TRACKS
th_G0=deg2rad(-170); %"greenwich" LST at t=0 [rad];


% NOMINAL ORBIT------------------------------------------------------------

% n_d = 0.5;
N = 1;
Orbit.t_orbit= N*(2*pi*sqrt(a0^3/mu));          %final time for N nominal orbit g-t
Orbit.t_fin= n_d *86400;                     %final time for n days nominal orbit g-t
N_points = 1000;
tspan = linspace(0, Orbit.t_fin, N_points);
omega_m = 14.62;  % Planet's angular velocity
t0 = 0;      % Initial time

% Repeating parameters
k = 1;
m = 1;




% Computation of the ground track

[lon1,lat1,~,~] = perturbedgroundtracksmars(a0, e0, i0, OM0,...
   om0, th0, Orbit.t_fin, th_G0, mu, omega_m, ...
   t0, R_mars,j2, 'none', 'unperturbed');

% Plot of the ground track
for ii=2:length(lon1)
    if abs(lon1(ii))>0.5
        if lon1(ii-1)*lon1(ii)<0
            lon1(ii-1)=NaN;
        end
    end
end

lon1=lon1*180/pi;
lat1=lat1*180/pi;
figure
A=imread('MarsTexture.jpg');
image('XData',[-180 180],'YData',[90 -90],'CData',A);
hold on
plot(lon1,lat1,'g','linewidth',1.2);
plot(lon1(1),lat1(1),'go','linewidth',2)
plot(lon1(length(lon1)),lat1(length(lat1)),'gs','linewidth',2)
xlim([-180 180])
ylim([-90 90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('Ground track','Start','End','Location','best','Orientation','horizontal')
title('Ground track plot')

%% Modification of the semi-major axis

a_perturbed = j2_groundTrack(omega_m, k, m, mu, j2, R_mars, a0, e0, i0);
T_perturbed = 2*pi*sqrt(a_perturbed^3/mu);
tspan_perturbed = linspace(0,Orbit.t_fin,N_points);

% Computation of the perturbed ground track
% [~, ~, lon_per, lat_per] = perturbed_groundTrack(kep_el, tspan_perturbed, mu, omega_m, t0, th_G0, R_mars, j2);
[lon_per,lat_per,~,~] = perturbedgroundtracksmars(a0, e0, i0, OM0,...
   om0, th0, Orbit.t_fin, th_G0, mu, omega_m, ...
   t0, R_mars,j2, 'orbit', 'j2');

[lon1,lat1,~,~] = perturbedgroundtracksmars(a0, e0, i0, OM0,...
   om0, th0, Orbit.t_fin, th_G0, mu, omega_m, ...
   t0, R_mars,j2, 'none', 'unperturbed');

for ii=2:length(lon1)
    if abs(lon1(ii))>0.5
        if lon1(ii-1)*lon1(ii)<0
            lon1(ii-1)=NaN;
        end
    end
end

lon1=lon1*180/pi;
lat1=lat1*180/pi;

figure
A=imread('MarsTexture.jpg');
image('XData',[-180 180],'YData',[90 -90],'CData',A);
hold on
plot(lon1,lat1,'g','linewidth',1.2);
xlim([-180 180])
ylim([-90 90])

% Plot of the repeating ground track for secular j2

for ii=2:length(lon_per)
    if abs(lon_per(ii))>0.5
        if lon_per(ii-1)*lon_per(ii)<0
            lon_per(ii-1)=NaN;
        end
    end
end

lon_per=lon_per*180/pi;
lat_per=lat_per*180/pi;
plot(lon_per,lat_per,'y','linewidth',1.2);
plot(lon1(1),lat1(1),'go','linewidth',2)
plot(lon1(length(lon1)),lat1(length(lat1)),'gs','linewidth',2)
plot(lon_per(1),lat_per(1),'yo','linewidth',2)
plot(lon_per(length(lon_per)),lat_per(length(lat_per)),'ys','linewidth',2)
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('Original Orbit','Perturbed Orbit','Start','End','Start','End','Location','best')
title('Ground track effect of perturbation')

%colouredplot(lon_per, lat_per, R_mars, h0, FOV, i0)


%% plot with MRO distance

% Computation of the perturbed ground track
% [~, ~, lon_per, lat_per] = perturbed_groundTrack(kep_el, tspan_perturbed, mu, omega_m, t0, th_G0, R_mars, j2);
[lon_mro,lat_mro,~,~] = perturbedgroundtracksmars(a_mro, e_mro, i_mro, OM_mro,...
   om_mro, th_mro, Orbit.t_fin, greenwich0_mro, mu, omega_m, ...
   t0, R_mars,j2, 'none', 'j2');

% [lon1mro,lat1mro,~,~] = perturbedgroundtracksmars(a0, e0, i0, OM0,...
%    om0, th0, Orbit.t_fin, myConst.greenwich0_mro, mu, omega_m, ...
%    t0, R_mars,j2, 'orbit', 'unperturbed');
% 
% for ii=2:length(lon1mro)
%     if abs(lon1mro(ii))>0.5
%         if lon1mro(ii-1)*lon1mro(ii)<0
%             lon1mro(ii-1)=NaN;
%         end
%     end
% end
% 
% lon1mro=lon1mro*180/pi;
% lat1mro=lat1mro*180/pi;
% 

% Plot of the repeating ground track for secular j2

for ii=2:length(lon_mro)
    if abs(lon_mro(ii))>0.5
        if lon_mro(ii-1)*lon_mro(ii)<0
            lon_mro(ii-1)=NaN;
        end
    end
end

lon_mro=lon_mro*180/pi;
lat_mro=lat_mro*180/pi;

figure
A=imread('MarsTexture.jpg');
image('XData',[-180 180],'YData',[90 -90],'CData',A);
hold on
plot(lon_mro,lat_mro,'g','linewidth',1.2);
xlim([-180 180])
ylim([-90 90])

plot(lon_per,lat_per,'y','linewidth',1.2);
plot(lon_mro(1),lat_mro(1),'go','linewidth',2)
plot(lon_mro(length(lon_mro)),lat_mro(length(lat_mro)),'gs','linewidth',2)
plot(lon_per(1),lat_per(1),'yo','linewidth',2)
plot(lon_per(length(lon_per)),lat_per(length(lat_per)),'ys','linewidth',2)
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('MRO Orbit','Satellite Orbit','Start','End','Start','End','Location','best')
title('Ground track confrontation')
% save('latlonpol_mro.mat', 'lon_mro', 'lat_mro')   
% colouredplot(lon_per, lat_per, R_mars, h0, FOV, i0)

delta_lon = abs(lon_per - lon_mro);
visible2_lon = delta_lon > 308.55;
visible1_lon = delta_lon < 51.45;
visible_lon = visible1_lon + visible2_lon;
visibility_lonpol = find(diff([0, visible_lon', 0]));
save('visibilitypol.mat', 'visibility_lonpol')

save('latlon_pol.mat', 'lon_per', 'lat_per')






    case 'sun-sync-ecl',

        
% apo =  R_mars + 320;
% peri = R_mars + 255;
% a0 = (apo + peri) / 2;
% e0 = (apo-peri)/(apo+peri);
p = a0*(1-e0^2);

OM0=-0.256563400043166;         %choosen RAAN [rad]

T_sidereal_mars = 686.980 * 24*3600;
rho = 2*pi/T_sidereal_mars;
i0=acos(-2*rho/(3*j2*R_mars^2*sqrt(mu))*(a0)^(3.5));         %inclination [rad]
Orbit.DELTA_OM = -3*pi*j2*R_mars^2/p^2*cos(i0);
% height max 5.4931e+03

Orbit.kep_el = [a0 e0 i0 OM0 om0 th0];
FOV = deg2rad(FOV);

% 2- GROUND TRACKS
th_G0=deg2rad(-170); %"greenwich" LST at t=0 [rad];


% NOMINAL ORBIT------------------------------------------------------------

% n_d = 0.5;
N = 1;
Orbit.t_orbit= N*(2*pi*sqrt(a0^3/mu));          %final time for N nominal orbit g-t
Orbit.t_fin= n_d *86400;                     %final time for n days nominal orbit g-t
N_points = 1000;
tspan = linspace(0, Orbit.t_fin, N_points);
omega_m = 14.62;  % Planet's angular velocity
t0 = 0;      % Initial time

% Repeating parameters
k = 1;
m = 1;


% Mars parameters
j2 = 0.00196045;

% Computation of the ground track

[lon1,lat1,~,~] = perturbedgroundtracksmars(a0, e0, i0, OM0,...
   om0, th0, Orbit.t_fin, th_G0, mu, omega_m, ...
   t0, R_mars,j2, 'none', 'unperturbed');

% Plot of the ground track
for ii=2:length(lon1)
    if abs(lon1(ii))>0.5
        if lon1(ii-1)*lon1(ii)<0
            lon1(ii-1)=NaN;
        end
    end
end

lon1=lon1*180/pi;
lat1=lat1*180/pi;
figure
A=imread('MarsTexture.jpg');
image('XData',[-180 180],'YData',[90 -90],'CData',A);
hold on
plot(lon1,lat1,'g','linewidth',1.2);
plot(lon1(1),lat1(1),'go','linewidth',2)
plot(lon1(length(lon1)),lat1(length(lat1)),'gs','linewidth',2)
xlim([-180 180])
ylim([-90 90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('Ground track','Start','End','Location','best','Orientation','horizontal')
title('Ground track plot')

%% Modification of the semi-major axis

a_perturbed = j2_groundTrack(omega_m, k, m, mu, j2, R_mars, a0, e0, i0);
T_perturbed = 2*pi*sqrt(a_perturbed^3/mu);
tspan_perturbed = linspace(0,Orbit.t_fin,N_points);

% Computation of the perturbed ground track
% [~, ~, lon_per, lat_per] = perturbed_groundTrack(kep_el, tspan_perturbed, mu, omega_m, t0, th_G0, R_mars, j2);
[lon_per,lat_per,~,~] = perturbedgroundtracksmars(a0, e0, i0, OM0,...
   om0, th0, Orbit.t_fin, th_G0, mu, omega_m, ...
   t0, R_mars,j2, 'orbit', 'j2');

[lon1,lat1,~,~] = perturbedgroundtracksmars(a0, e0, i0, OM0,...
   om0, th0, Orbit.t_fin, th_G0, mu, omega_m, ...
   t0, R_mars,j2, 'none', 'unperturbed');

for ii=2:length(lon1)
    if abs(lon1(ii))>0.5
        if lon1(ii-1)*lon1(ii)<0
            lon1(ii-1)=NaN;
        end
    end
end

lon1=lon1*180/pi;
lat1=lat1*180/pi;

figure
A=imread('MarsTexture.jpg');
image('XData',[-180 180],'YData',[90 -90],'CData',A);
hold on
plot(lon1,lat1,'y','linewidth',1.2);
xlim([-180 180])
ylim([-90 90])

% Plot of the repeating ground track for secular j2

for ii=2:length(lon_per)
    if abs(lon_per(ii))>0.5
        if lon_per(ii-1)*lon_per(ii)<0
            lon_per(ii-1)=NaN;
        end
    end
end

lon_per=lon_per*180/pi;
lat_per=lat_per*180/pi;
plot(lon_per,lat_per,'g','linewidth',1.2);
plot(lon1(1),lat1(1),'yo','linewidth',2)
plot(lon1(length(lon1)),lat1(length(lat1)),'ys','linewidth',2)
plot(lon_per(1),lat_per(1),'go','linewidth',2)
plot(lon_per(length(lon_per)),lat_per(length(lat_per)),'gs','linewidth',2)
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('Original Orbit','Perturbed Orbit','Start','End','Start','End','Location','best')
title('Ground track effect of perturbation')
    
%colouredplot(lon_per, lat_per, R_mars, h0, FOV, i0)

%% plot with MRO distance

% Computation of the perturbed ground track
% [~, ~, lon_per, lat_per] = perturbed_groundTrack(kep_el, tspan_perturbed, mu, omega_m, t0, th_G0, R_mars, j2);
[lon_mro,lat_mro,~,~] = perturbedgroundtracksmars(a_mro, e_mro, i_mro, OM_mro,...
   om_mro, th_mro, Orbit.t_fin, greenwich0_mro, mu, omega_m, ...
   t0, R_mars,j2, 'none', 'j2');

% [lon1mro,lat1mro,~,~] = perturbedgroundtracksmars(a0, e0, i0, OM0,...
%    om0, th0, Orbit.t_fin, myConst.greenwich0_mro, mu, omega_m, ...
%    t0, R_mars,j2, 'orbit', 'unperturbed');
% 
% for ii=2:length(lon1mro)
%     if abs(lon1mro(ii))>0.5
%         if lon1mro(ii-1)*lon1mro(ii)<0
%             lon1mro(ii-1)=NaN;
%         end
%     end
% end
% 
% lon1mro=lon1mro*180/pi;
% lat1mro=lat1mro*180/pi;
% 

% Plot of the repeating ground track for secular j2

for ii=2:length(lon_mro)
    if abs(lon_mro(ii))>0.5
        if lon_mro(ii-1)*lon_mro(ii)<0
            lon_mro(ii-1)=NaN;
        end
    end
end

lon_mro=lon_mro*180/pi;
lat_mro=lat_mro*180/pi;

figure
A=imread('MarsTexture.jpg');
image('XData',[-180 180],'YData',[90 -90],'CData',A);
hold on
plot(lon_mro,lat_mro,'g','linewidth',1.2);
xlim([-180 180])
ylim([-90 90])

plot(lon_per,lat_per,'y','linewidth',1.2);
plot(lon_mro(1),lat_mro(1),'go','linewidth',2)
plot(lon_mro(length(lon_mro)),lat_mro(length(lat_mro)),'gs','linewidth',2)
plot(lon_per(1),lat_per(1),'yo','linewidth',2)
plot(lon_per(length(lon_per)),lat_per(length(lat_per)),'ys','linewidth',2)
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('MRO Orbit','Satellite Orbit','Start','End','Start','End','Location','best')
title('Ground track confrontation')
% save('latlonsse_mro.mat', 'lon_mro', 'lat_mro')   
%colouredplot(lon_per, lat_per, R_mars, h0, FOV, i0)


delta_lon = abs(lon_per - lon_mro);
visible2_lon = delta_lon > 308.55;
visible1_lon = delta_lon < 51.45;
visible_lon = visible1_lon + visible2_lon;
visibility_lonsse = find(diff([0, visible_lon', 0]));
save('visibilitysse.mat', 'visibility_lonsse')
save('latlon_sse.mat', 'lon_per', 'lat_per')









    case 'sun-sync-noecl',

        
% apo =  R_mars + 320;
% peri = R_mars + 255;
% a0 = (apo + peri) / 2;
% e0 = (apo-peri)/(apo+peri);
p = a0*(1-e0^2);

OM0=deg2rad(270);         %choosen RAAN [rad]

T_sidereal_mars = 686.980 * 24*3600;
rho = 2*pi/T_sidereal_mars;
i0=acos(-2*rho/(3*j2*R_mars^2*sqrt(mu))*(a0)^(3.5));         %inclination [rad]
Orbit.DELTA_OM = -3*pi*j2*R_mars^2/p^2*cos(i0);
% height max 5.4931e+03

Orbit.kep_el = [a0 e0 i0 OM0 om0 th0];
FOV = deg2rad(FOV);

% 2- GROUND TRACKS
th_G0=deg2rad(-170); %"greenwich" LST at t=0 [rad];


% NOMINAL ORBIT------------------------------------------------------------

% n_d = 0.5;
N = 1;
Orbit.t_orbit= N*(2*pi*sqrt(a0^3/mu));          %final time for N nominal orbit g-t
Orbit.t_fin= n_d *86400;                     %final time for n days nominal orbit g-t
N_points = 1000;
tspan = linspace(0, Orbit.t_fin, N_points);
omega_m = 14.62;  % Planet's angular velocity
t0 = 0;      % Initial time

% Repeating parameters
k = 1;
m = 1;


% Mars parameters
j2 = 0.00196045;

% Computation of the ground track

[lon1,lat1,~,~] = perturbedgroundtracksmars(a0, e0, i0, OM0,...
   om0, th0, Orbit.t_fin, th_G0, mu, omega_m, ...
   t0, R_mars,j2, 'none', 'unperturbed');

% Plot of the ground track
for ii=2:length(lon1)
    if abs(lon1(ii))>0.5
        if lon1(ii-1)*lon1(ii)<0
            lon1(ii-1)=NaN;
        end
    end
end

lon1=lon1*180/pi;
lat1=lat1*180/pi;
figure
A=imread('MarsTexture.jpg');
image('XData',[-180 180],'YData',[90 -90],'CData',A);
hold on
plot(lon1,lat1,'g','linewidth',1.2);
plot(lon1(1),lat1(1),'go','linewidth',2)
plot(lon1(length(lon1)),lat1(length(lat1)),'gs','linewidth',2)
xlim([-180 180])
ylim([-90 90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('Ground track','Start','End','Location','best','Orientation','horizontal')
title('Ground track plot')

%% Modification of the semi-major axis

a_perturbed = j2_groundTrack(omega_m, k, m, mu, j2, R_mars, a0, e0, i0);
T_perturbed = 2*pi*sqrt(a_perturbed^3/mu);
tspan_perturbed = linspace(0,Orbit.t_fin,N_points);

% Computation of the perturbed ground track
% [~, ~, lon_per, lat_per] = perturbed_groundTrack(kep_el, tspan_perturbed, mu, omega_m, t0, th_G0, R_mars, j2);
[lon_per,lat_per,~,~] = perturbedgroundtracksmars(a0, e0, i0, OM0,...
   om0, th0, Orbit.t_fin, th_G0, mu, omega_m, ...
   t0, R_mars,j2, 'orbit', 'j2');

[lon1,lat1,~,~] = perturbedgroundtracksmars(a0, e0, i0, OM0,...
   om0, th0, Orbit.t_fin, th_G0, mu, omega_m, ...
   t0, R_mars,j2, 'none', 'unperturbed');

for ii=2:length(lon1)
    if abs(lon1(ii))>0.5
        if lon1(ii-1)*lon1(ii)<0
            lon1(ii-1)=NaN;
        end
    end
end

lon1=lon1*180/pi;
lat1=lat1*180/pi;

figure
A=imread('MarsTexture.jpg');
image('XData',[-180 180],'YData',[90 -90],'CData',A);
hold on
plot(lon1,lat1,'y','linewidth',1.2);
xlim([-180 180])
ylim([-90 90])

% Plot of the repeating ground track for secular j2

for ii=2:length(lon_per)
    if abs(lon_per(ii))>0.5
        if lon_per(ii-1)*lon_per(ii)<0
            lon_per(ii-1)=NaN;
        end
    end
end

lon_per=lon_per*180/pi;
lat_per=lat_per*180/pi;
plot(lon_per,lat_per,'g','linewidth',1.2);
plot(lon1(1),lat1(1),'yo','linewidth',2)
plot(lon1(length(lon1)),lat1(length(lat1)),'ys','linewidth',2)
plot(lon_per(1),lat_per(1),'go','linewidth',2)
plot(lon_per(length(lon_per)),lat_per(length(lat_per)),'gs','linewidth',2)
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('Original Orbit','Perturbed Orbit','Start','End','Start','End','Location','best')
title('Ground track effect of perturbation')
    
%colouredplot(lon_per, lat_per, R_mars, h0, FOV, i0)

%% plot with MRO distance

% Computation of the perturbed ground track
% [~, ~, lon_per, lat_per] = perturbed_groundTrack(kep_el, tspan_perturbed, mu, omega_m, t0, th_G0, R_mars, j2);
[lon_mro,lat_mro,~,~] = perturbedgroundtracksmars(a_mro, e_mro, i_mro, OM_mro,...
   om_mro, th_mro, Orbit.t_fin, greenwich0_mro, mu, omega_m, ...
   t0, R_mars,j2, 'none', 'j2');

% [lon1mro,lat1mro,~,~] = perturbedgroundtracksmars(a0, e0, i0, OM0,...
%    om0, th0, Orbit.t_fin, myConst.greenwich0_mro, mu, omega_m, ...
%    t0, R_mars,j2, 'orbit', 'unperturbed');
% 
% for ii=2:length(lon1mro)
%     if abs(lon1mro(ii))>0.5
%         if lon1mro(ii-1)*lon1mro(ii)<0
%             lon1mro(ii-1)=NaN;
%         end
%     end
% end
% 
% lon1mro=lon1mro*180/pi;
% lat1mro=lat1mro*180/pi;
% 

% Plot of the repeating ground track for secular j2

for ii=2:length(lon_mro)
    if abs(lon_mro(ii))>0.5
        if lon_mro(ii-1)*lon_mro(ii)<0
            lon_mro(ii-1)=NaN;
        end
    end
end

lon_mro=lon_mro*180/pi;
lat_mro=lat_mro*180/pi;

figure
A=imread('MarsTexture.jpg');
image('XData',[-180 180],'YData',[90 -90],'CData',A);
hold on
plot(lon_mro,lat_mro,'g','linewidth',1.2);
xlim([-180 180])
ylim([-90 90])

plot(lon_per,lat_per,'y','linewidth',1.2);
plot(lon_mro(1),lat_mro(1),'go','linewidth',2)
plot(lon_mro(length(lon_mro)),lat_mro(length(lat_mro)),'gs','linewidth',2)
plot(lon_per(1),lat_per(1),'yo','linewidth',2)
plot(lon_per(length(lon_per)),lat_per(length(lat_per)),'ys','linewidth',2)
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
legend('MRO Orbit','Satellite Orbit','Start','End','Start','End','Location','best')
title('Ground track confrontation')
    
%colouredplot(lon_per, lat_per, R_mars, h0, FOV, i0)
%save('latlonssne_mro.mat', 'lon_mro', 'lat_mro')

delta_lon = abs(lon_per - lon_mro);
visible2_lon = delta_lon > 308.55;
visible1_lon = delta_lon < 51.45;
visible_lon = visible1_lon + visible2_lon;
visibility_lonssne = find(diff([0, visible_lon', 0]));
save('visibilityssne.mat', 'visibility_lonssne')
save('latlon_ssne.mat', 'lon_per', 'lat_per')



end

Orbit.visibility_lon = find(visible_lon);


rho = asin(R_mars/(R_mars+h0));
eps = acos(sin(FOV)/sin(rho));
Orbit.lambda = pi/2-FOV-eps;

%colouredplot(lon_per, lat_per, R_mars, h0, FOV, i0);
Orbit.lon = lon_per;
Orbit.lat = lat_per;

save('latlon_mro.mat', 'lon_mro', 'lat_mro')
end


%% other

% % nominal orbit ground tracks (unperturbed) 
% [lon_u,lat_u,~,~] =...
%     groundtracks(a0,e0,i0,OM0,om0,th0,t_fin,'unperturbed',th_gw0,'all');
% 
% % nominal orbit ground tracks (j2) ---> must find J2 mars
% [lon_p,lat_p,~,~] =...
%     groundtracks(a0,e0,i0,OM0,om0,th0,t_fin,'j2',th_gw0,'all');

% %% plot results
% figure
% I =imread('mars.jpg');
% image([-180,180],[90 , -90],I)
% hold on
% plot(rad2deg(lon_u),rad2deg(lat_u),'g','linewidth',0.8,'handlevisibility','off')
% plot(rad2deg(lon_u(1)),rad2deg(lat_u(1)),'yo','markersize',5,'linewidth',4)
% plot(rad2deg(lon_u(end)),rad2deg(lat_u(end)),'gx','markersize',13,'linewidth',4)
% %plot(rad2deg(lon_p),rad2deg(lat_p),'r','linewidth',0.8,'handlevisibility','off')
% %plot(rad2deg(lon_p(end)),rad2deg(lat_p(end)),'rx','markersize',13,'linewidth',4)
% grid on, xlabel('Longitude [°]'), ylabel('Latitude [°]')
% leg=legend('Inital point','Final point (unperturbed)');             %, 'Final point (perturbed)'
% set(leg, 'TextColor','w', 'Color','none', 'EdgeColor','w')
% ax = gca; ax.YDir = 'normal'; ax.GridColor = 'white';
% 
% 
% [Rxo,Ryo,Rzo] = plot_orbit_dt(a0,e0,i0,OM0,om0,th0,2*pi, 0.001, mu); 
% [Rxp_i,Ryp_i,Rzp_i]= plot_point(a0,e0,i0,OM0,om0,th0);
% [Rxp_f,Ryp_f,Rzp_f]= plot_point(a0,e0,i0,OM0,om0,2*pi);
% 
% 
% figure
% plot3 (Rxo,Ryo,Rzo,'b')                                     %orbit
% hold on
% grid on
% scatter3 (Rxp_i,Ryp_i,Rzp_i,'c','filled');                  %initial point
% theta_peri = 0;
% theta_ap = pi;
% scatter3(Rxp_f,Ryp_f,Rzp_f,'r','*');                        %final point
% 
% %Stars
% S = 5.*10e3 - 10e4.*rand(3,100);
% scatter3(S(1,:),S(2,:),S(3,:),0.2,'w', 'handleVisibility', 'off')
% 
% % Plot Mars and black background
% Celestial_body('mars',1); 
% whitebg('k')
% title 'One Orbit'
% axis equal
% legend('One Orbit Propagation', 'Initial Point', 'Final Point')



