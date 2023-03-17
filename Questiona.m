%% Fluid mechanics coursework (Part 1)
% Created by Kenneth Yap
% 21st of November 2021
% Part 1

clear
clc
close all

%% Constants

% Quanta emission rate (approximate)
G_breath_rest = 1; % /hr
G_breath_heavy = 8; % /hr
G_speak_light = 15; % /hr
G_speak_loud = 98; % /hr

% Turbulent diffusion coefficient
D_T = 1; % m^2/s

% Average height of man and women
h = 1.68; % m

%% Question a
% Conditions of no wind means no advection, only diffusion
% Assumptions:
% - The quanta has unit of Lˆ3 which makes equivalent to M
% - Continuous release with the floor as a perfect reflector
% - 3D diffusion
% - Audience is speaking during light activity and breathing during heavy
%   activity, so an average quanta is taken (15+8)/2 = 11.5

% Convert quanta from hours to second
G = G_speak_loud*4/(60*60); % account for 4 persons

% Time taken for 3 hours in terms of seconds
t = 0:3*60*60;

% Calculation for Platform 1
x = -8; y = 8; z = 0; D = D_T;
r = sqrt(x^2 + y^2 + (z)^2); r_reflect = sqrt(x^2 + y^2 + (z + 2*h)^2);
C_1 = G/(4*pi*D*r)*erfc(r./sqrt(4*D.*t)) + G/(4*pi*D*r_reflect)*erfc(r_reflect./sqrt(4*D.*t));

% Calculation for Platform 3
x = 0; y = 8; z = 0; D = D_T;
r = sqrt(x^2 + y^2 + (z)^2); r_reflect = sqrt(x^2 + y^2 + (z + 2*h)^2);
C_3 = G/(4*pi*D*r)*erfc(r./sqrt(4*D.*t)) + G/(4*pi*D*r_reflect)*erfc(r_reflect./sqrt(4*D.*t));

% Calculation for Platform 7
x = -4; y = 4; z = 0; D = D_T;
r = sqrt(x^2 + y^2 + (z)^2); r_reflect = sqrt(x^2 + y^2 + (z + 2*h)^2);
C_7 = G/(4*pi*D*r)*erfc(r./sqrt(4*D.*t)) + G/(4*pi*D*r_reflect)*erfc(r_reflect./sqrt(4*D.*t));

% Calculation for Platform 13
x = 0; y = 1; z = 0; D = D_T; % Account for infinitely large value
r = sqrt(x^2 + y^2 + (z)^2); r_reflect = sqrt(x^2 + y^2 + (z + 2*h)^2);
C_13 = G/(4*pi*D*r)*erfc(r./sqrt(4*D.*t)) + G/(4*pi*D*r_reflect)*erfc(r_reflect./sqrt(4*D.*t));

figure
hold on
plot(t,C_1)
plot(t,C_3)
plot(t,C_7)
plot(t,C_13)
xlabel('Time (s)')
ylabel('Concentration of virus (m^{-3})')
title('Temporal variation of virus')
legend('Platform 1', 'Platform 3', 'Platform 7', 'Platform 13')
hold off

%% Question b
% similarities and differences in the spatial distribution of the virus concentration within the public area at head
% height at the end of the concert in the following conditions:
% (i) no wind (calculated in question a)
% (ii) Spatially uniform wind of 3m/s from south west

% Split velocity componenets into x and y directions
U = 3; % m/s
z = 0; % m
% Subsitute into continuous release advection-diffusion equation

% Calculation of the maximum at each location
% Platform 1
x = -8; y = 8; D = D_T;
r = sqrt(x^2 + y^2 + (z)^2); r_reflect = sqrt(x^2 + y^2 + (z + 2*h)^2);
C_nowind_1 = G/(4*pi*D*r) + G/(4*pi*D*r_reflect);
C_wind_1 = G/(4*pi*D*r)*exp(-U/(2*D)*(r-y)) + G/(4*pi*D*r_reflect)*exp(-U/(2*D)*(r_reflect-y));

% Platform 3
x = 0; y = 8; D = D_T;
r = sqrt(x^2 + y^2 + (z)^2); r_reflect = sqrt(x^2 + y^2 + (z + 2*h)^2);
C_nowind_3 = G/(4*pi*D*r) + G/(4*pi*D*r_reflect);
C_wind_3 = G/(4*pi*D*r)*exp(-U/(2*D)*(r-y)) + G/(4*pi*D*r_reflect)*exp(-U/(2*D)*(r_reflect-y));

% Platform 7
x = -4; y = 4; D = D_T;
r = sqrt(x^2 + y^2 + (z)^2); r_reflect = sqrt(x^2 + y^2 + (z + 2*h)^2);
C_nowind_7 = G/(4*pi*D*r) + G/(4*pi*D*r_reflect);
C_wind_7 = G/(4*pi*D*r)*exp(-U/(2*D)*(r-y)) + G/(4*pi*D*r_reflect)*exp(-U/(2*D)*(r_reflect-y));

% Platform 13
x = 0; y = 1; D = D_T;
r = sqrt(x^2 + y^2 + (z)^2); r_reflect = sqrt(x^2 + y^2 + (z + 2*h)^2);
C_nowind_13 = G/(4*pi*D*r) + G/(4*pi*D*r_reflect);
C_wind_13 = G/(4*pi*D*r)*exp(-U/(2*D)*(r-y)) + G/(4*pi*D*r_reflect)*exp(-U/(2*D)*(r_reflect-y));

% Plot figures
figure
hold on
platform = categorical({'Platform 1','Platform 3','Platform 7','Edge of platform 13'});
platform = reordercats(platform,{'Platform 1','Platform 3','Platform 7','Edge of platform 13'});
wind_vs_nowind = [C_nowind_1 C_nowind_3 C_nowind_7 C_nowind_13; C_wind_1 C_wind_3 C_wind_7 C_wind_13];
title('Concentration outdoor')
ylabel('Concentration (m^-3)')
b = bar(platform, wind_vs_nowind);
set(b, {'DisplayName'}, {'No wind','Wind'}')
legend('Location','northwest')
hold off

% Plot out surface plots
% This is for not sanitized floor
x = linspace(-9,9,25);
y = linspace(-9,9,25);
z = 0;

C_nowind_gen = zeros(length(x),length(y));
C_wind_gen = zeros(length(x),length(y));

for i = 1:length(x)
    for j = 1:length(y)
        r = sqrt(x(i).^2 + y(j).^2 + (z).^2); 
        r_reflect = sqrt(x(i).^2 + y(j).^2 + (z + 2*h).^2);
        C_nowind_gen(i,j) = G./(4*pi*D*r) + G./(4*pi*D.*r_reflect);
        C_wind_gen(i,j) = G./(4*pi*D*r).*exp(-U./(2*D).*(r-y(j))) + G./(4*pi*D.*r_reflect).*exp(-U./(2*D).*(r_reflect-y(j)));
    end
end

figure
hold on
subplot(1,2,1)
contourf(x,y,C_nowind_gen')
xlabel('x coordinates [m]')
ylabel('y coordinates [m]')
title({'Concentration without applying disinfectant'}, {'without wind'})
colorbar
axis equal
hold off

subplot(1,2,2)
hold on
contourf(x,y,C_wind_gen')
xlabel('x coordinates [m]')
ylabel('y coordinates [m]')
title({'Concentration without applying disinfectant'}, {'with wind'})
colorbar
hold off
axis equal

% figure
% subplot(2,1,1)
% surf(x,y,C_wind_gen')
% xlabel('x coordinates [m]')
% ylabel('y coordinates [m]')
% zlabel('Concentration (m^-3)')
% title('Concentration with sanitizer with wind')
% caxis([0,0.01])
% colorbar
% 
% subplot(2,1,2)
% surf(x,y,C_nowind_gen')
% xlabel('x coordinates [m]')
% ylabel('y coordinates [m]')
% zlabel('Concentration (m^-3)')
% title('Concentration with sanitizer without wind')
% caxis([0,0.01])
% colorbar


%% Question c
% Assuming the floor is appied with disinfectant, the floor becomes a
% perfect absorber
x = linspace(-9,9,25);
y = linspace(-9,9,25);
z = 0;

C_nowind_san = zeros(length(x),length(y));
C_wind_san = zeros(length(x),length(y));

for i = 1:length(x)
    for j = 1:length(y)
        r = sqrt(x(i).^2 + y(j).^2 + (z).^2); 
        r_reflect = sqrt(x(i).^2 + y(j).^2 + (z + 2*h).^2);
        C_nowind_san(i,j) = G./(4*pi*D*r) - G./(4*pi*D.*r_reflect);
        C_wind_san(i,j) = G./(4*pi*D*r).*exp(-U./(2*D).*(r-y(j))) - G./(4*pi*D.*r_reflect).*exp(-U./(2*D).*(r_reflect-y(j)));
    end
end

figure
hold on
subplot(2,2,1)
contourf(x,y,C_nowind_gen')
xlabel('x coordinates [m]')
ylabel('y coordinates [m]')
title({'Concentration without applying disinfectant'}, {'without wind'})
caxis([0 0.012])
colorbar
axis equal
hold off

subplot(2,2,2)
hold on
contourf(x,y,C_wind_gen')
xlabel('x coordinates [m]')
ylabel('y coordinates [m]')
title({'Concentration without applying disinfectant'}, {'with wind'})
caxis([0 0.012])
colorbar
hold off
axis equal

hold on
subplot(2,2,3)
contourf(x,y,C_nowind_san')
xlabel('x coordinates [m]')
ylabel('y coordinates [m]')
title({'Concentration after applying disinfectant'}, {'without wind'})
caxis([0 0.012])
colorbar
axis equal
hold off

subplot(2,2,4)
hold on
contourf(x,y,C_wind_san')
xlabel('x coordinates [m]')
ylabel('y coordinates [m]')
title({'Concentration after applying disinfectant'}, {'with wind'})
caxis([0 0.012])
colorbar
axis equal
hold off

% figure
% subplot(2,1,1)
% surf(x,y,C_wind_san')
% xlabel('x coordinates [m]')
% ylabel('y coordinates [m]')
% zlabel('Concentration (m^-3)')
% title('Concentration with sanitizer with wind')
% caxis([0,0.01])
% colorbar
% 
% subplot(2,1,2)
% surf(x,y,C_nowind_san')
% xlabel('x coordinates [m]')
% ylabel('y coordinates [m]')
% zlabel('Concentration (m^-3)')
% title('Concentration with sanitizer without wind')
% caxis([0,0.01])
% colorbar


%% Question d
% Plot the variation in height at the stage where h = 3+1.68

x = 0;
y = 15;
h = 1.68;
z = linspace(-(h+3),10,100);

C_nowind_san_ver = zeros(1,length(z));
C_wind_san_ver = zeros(1,length(z));
C_nowind_gen_ver = zeros(1,length(z));
C_wind_gen_ver = zeros(1,length(z));

for k = 1:length(z)

    r = sqrt(x.^2 + y.^2 + (z(k)).^2); 
    r_reflect = sqrt(x.^2 + y.^2 + (z(k) + 2*h).^2);

    C_nowind_san_ver(k) = G/(4*pi*D*r)*erfc(r./sqrt(4*D.*t(end))) - G/(4*pi*D*r_reflect)*erfc(r_reflect./sqrt(4*D.*t(end)));
    C_nowind_gen_ver(k) = G/(4*pi*D*r)*erfc(r./sqrt(4*D.*t(end))) + G/(4*pi*D*r_reflect)*erfc(r_reflect./sqrt(4*D.*t(end)));
    C_wind_san_ver(k) = G/(4*pi*D*r).*exp(-U./(2*D).*(r-y)) - G./(4*pi*D.*r_reflect).*exp(-U./(2*D).*(r_reflect-y));
    C_wind_gen_ver(k) = G/(4*pi*D*r).*exp(-U./(2*D).*(r-y)) + G./(4*pi*D.*r_reflect).*exp(-U./(2*D).*(r_reflect-y));

end

figure
hold on 
plot(C_wind_gen_ver,z+h,'--r')
plot(C_nowind_gen_ver,z+h, 'r')
plot(C_wind_san_ver,z+h, '--b')
plot(C_nowind_san_ver,z+h, 'b')
legend('Wind, not sanitized','No wind, not sanitized', 'Wind, sanitized', 'No wind, sanitized')
legend('Location','NorthEastOutside')
xlim([0,1.2*10^-3])
xlabel('Concentration[mˆ{-3}]')
ylabel('Height [m]')
title('Variation of concentration at stage')
ylim([0,10])
hold off

%% Question e
% Plot out the concentration at the concert place given that the platform is
% at a height of 2m

h = 2+1.68; % m distance from floor
x = linspace(-9,9,25);
y = linspace(-9,9,25);
z = 0; % head height increase by 2m

C_nowind_gen_inc = zeros(length(x),length(y));
C_wind_gen_inc = zeros(length(x),length(y));
C_nowind_san_inc = zeros(length(x),length(y));
C_wind_san_inc = zeros(length(x),length(y));

C_nowind_san_ver_inc = zeros(1,length(z));
C_wind_san_ver_inc = zeros(1,length(z));
C_nowind_gen_ver_inc = zeros(1,length(z));
C_wind_gen_ver_inc = zeros(1,length(z));

for i = 1:length(x)
    for j = 1:length(y)
        r = sqrt(x(i).^2 + y(j).^2 + (z).^2); 
        r_reflect = sqrt(x(i).^2 + y(j).^2 + (z + 2*h).^2);
        C_nowind_gen_inc(i,j) = G./(4*pi*D*r) + G./(4*pi*D.*r_reflect);
        C_wind_gen_inc(i,j) = G./(4*pi*D*r).*exp(-U./(2*D).*(r-y(j))) + G./(4*pi*D.*r_reflect).*exp(-U./(2*D).*(r_reflect-y(j)));
        C_nowind_san_inc(i,j) = G./(4*pi*D*r) - G./(4*pi*D.*r_reflect);
        C_wind_san_inc(i,j) = G./(4*pi*D*r).*exp(-U./(2*D).*(r-y(j))) - G./(4*pi*D.*r_reflect).*exp(-U./(2*D).*(r_reflect-y(j)));
    end
end

figure
subplot(2,2,1)
contourf(x,y,C_nowind_gen_inc')
xlabel('x coordinates [m]')
ylabel('y coordinates [m]')
title({'Concentration without applying disinfectant'}, {'without wind at h = 3.68m'})
ylim([-3,3])
caxis([0 0.0045])
axis equal
colorbar

subplot(2,2,2)
contourf(x,y,C_wind_gen_inc')
xlabel('x coordinates [m]')
ylabel('y coordinates [m]')
title({'Concentration without applying disinfectant'}, {'with wind at h = 3.68m'})
ylim([-3,3])
caxis([0 0.0045])
axis equal
colorbar

subplot(2,2,3)
contourf(x,y,C_nowind_san_inc')
xlabel('x coordinates [m]')
ylabel('y coordinates [m]')
title({'Concentration after applying disinfectant'}, {'without wind at h = 3.68m'})
ylim([-3,3])
caxis([0 0.0045])
axis equal
colorbar

subplot(2,2,4)
contourf(x,y,C_wind_san_inc')
xlabel('x coordinates [m]')
ylabel('y coordinates [m]')
title({'Concentration after applying disinfectant'}, {'without wind at h = 3.68m'})
ylim([-3,3])
caxis([0 0.0045])
axis equal
colorbar

% figure
% subplot(2,1,1)
% surf(x,y,C_wind_san_inc')
% xlabel('x coordinates [m]')
% ylabel('y coordinates [m]')
% zlabel('Concentration (m^-3)')
% title('Concentration sanitized with wind at h = 3.68m')
% caxis([0, 0.01])
% colorbar
% 
% subplot(2,1,2)
% surf(x,y,C_nowind_san_inc')
% xlabel('x coordinates [m]')
% ylabel('y coordinates [m]')
% zlabel('Concentration (m^-3)')
% title('Concentration sanitized without wind at h = 3.68m')
% caxis([0, 0.01])
% colorbar

% figure
% subplot(2,1,1)
% surf(x,y,C_wind_gen_inc')
% xlabel('x coordinates [m]')
% ylabel('y coordinates [m]')
% zlabel('Concentration (m^-3)')
% title('Concentration not sanitized with wind at h = 3.68m')
% caxis([0, 0.01])
% colorbar
% 
% subplot(2,1,2)
% surf(x,y,C_nowind_gen_inc')
% xlabel('x coordinates [m]')
% ylabel('y coordinates [m]')
% zlabel('Concentration (m^-3)')
% title('Concentration not sanitized without wind at h = 3.68m')
% caxis([0, 0.01])
% colorbar

x = 0;
y = 13;
z = linspace(-h,2,100);

for k = 1:length(z)
    r = sqrt(x.^2 + y.^2 + (z(k)).^2); 
    r_reflect = sqrt(x.^2 + y.^2 + (z(k) + 2*h).^2);
    
    C_nowind_san_ver_inc(k) = G/(4*pi*D*r)*erfc(r./sqrt(4*D.*t(end))) - G/(4*pi*D*r_reflect)*erfc(r_reflect./sqrt(4*D.*t(end)));
    C_wind_san_ver_inc(k) = G./(4*pi*D*r).*exp(- U./(2*D).*(r-y)) - G./(4*pi*D.*r_reflect).*exp(-U./(2*D).*(r_reflect-y));
    C_nowind_gen_ver_inc(k) = G/(4*pi*D*r)*erfc(r./sqrt(4*D.*t(end))) + G/(4*pi*D*r_reflect)*erfc(r_reflect./sqrt(4*D.*t(end)));
    C_wind_gen_ver_inc(k) = G./(4*pi*D*r).*exp(-U./(2*D).*(r-y)) + G./(4*pi*D.*r_reflect).*exp(-U./(2*D).*(r_reflect-y));

end

figure
hold on 
plot(C_wind_gen_ver_inc,z+h,'--r')
plot(C_nowind_gen_ver_inc,z+h, 'r')
plot(C_wind_san_ver_inc,z+h, '--b')
plot(C_nowind_san_ver_inc,z+h, 'b')
legend('Wind, not sanitized','No wind, not sanitized', 'Wind, sanitized', 'No wind, sanitized')
legend('Location','NorthEastOutside')
xlabel('Concentration[mˆ-3]')
ylabel('Height [m]')
title('Variation of concentration with height at h = 3.68 [m]')
hold off