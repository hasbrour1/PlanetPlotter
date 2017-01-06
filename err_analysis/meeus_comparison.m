%Name: meeus_comparison
%Date: 12/5/16
%Purpose: Compare Meeus Ephemerides to JPL DE405 Numerical Integration
%Model

%Initialize Console:
clear; close all; clc;

%JPL Astronomical Constants Required Include File:
%Include JPL DE405 directory:
addpath('./JPL_DE405');
AST_Const

%Constants:
mu_Sun = 1.3271244018e11; %(km^3/s^2)

%Select Comparison Epoch:
count = 0;
yr_vec = [2016:3:2056];
for yr =  yr_vec
count = count +1;
[ JD ] = get_JD( 12,12,11,yr ); %Julian Date at 12:00 UT, December 11, 2016 

%% Generate Meeus Ephemeris:

%Planetary Ephemerides at Mercury:
ephem = meeus_ephem(JD);
%Compute Earth position:
[r_Mercury_M,~] = COE2RV(ephem.mercury_a,ephem.mercury_e,rad2deg(ephem.mercury_i),...
                        rad2deg(ephem.mercury_w),rad2deg(ephem.mercury_Omega),rad2deg(ephem.mercury_nu),...
                        mu_Sun); %Mercury Position, Helio-Ecliptic [x, y, z] (km)
r_Mercury_M = r_Mercury_M*1e3; %Mercury Position, Helio-Ecliptic [x, y, z] (m)
  
                                      
%Planetary Ephemerides at Venus:
%Compute Venus position:
[r_Venus_M,~] = COE2RV(ephem.venus_a,ephem.venus_e,rad2deg(ephem.venus_i),...
                        rad2deg(ephem.venus_w),rad2deg(ephem.venus_Omega),rad2deg(ephem.venus_nu),...
                        mu_Sun); %Venus Position, Helio-Ecliptic [x y z] (km)
r_Venus_M = r_Venus_M*1e3; %Venus Position, Helio-Ecliptic [x y z] (m)

%Compute Earth position:
[r_Earth_M,~] = COE2RV(ephem.earth_a,ephem.earth_e,rad2deg(ephem.earth_i),...
                        rad2deg(ephem.earth_w),rad2deg(ephem.earth_Omega),rad2deg(ephem.earth_nu),...
                        mu_Sun); %Earth Position, Helio-Ecliptic [x, y, z] (km)
r_Earth_M = r_Earth_M*1e3; %Earth Position, Helio-Ecliptic [x, y, z] (m)
  

%Planetary Ephemerides at Mars:
%Compute Mars position:
[r_Mars_M,~] = COE2RV(ephem.mars_a,ephem.mars_e,rad2deg(ephem.mars_i),...
                        rad2deg(ephem.mars_w),rad2deg(ephem.mars_Omega),rad2deg(ephem.mars_nu),...
                        mu_Sun); %Mars Position, Helio-Ecliptic [x y z] (km)
r_Mars_M = r_Mars_M*1e3; %Mars Position, Helio-Ecliptic [x y z] (m)

%Planetary Ephemerides at Jupiter:
%Compute Jupiter position:
[r_Jupiter_M,~] = COE2RV(ephem.jupiter_a,ephem.jupiter_e,rad2deg(ephem.jupiter_i),...
                        rad2deg(ephem.jupiter_w),rad2deg(ephem.jupiter_Omega),rad2deg(ephem.jupiter_nu),...
                        mu_Sun); %Jupiter Position, Helio-Ecliptic [x y z] (km)
r_Jupiter_M = r_Jupiter_M*1e3; %Jupiter Position, Helio-Ecliptic [x y z] (m)
                    
%Planetary Ephemerides at Saturn:
%Compute Saturn position:
[r_Saturn_M,~] = COE2RV(ephem.saturn_a,ephem.saturn_e,rad2deg(ephem.saturn_i),...
                        rad2deg(ephem.saturn_w),rad2deg(ephem.saturn_Omega),rad2deg(ephem.saturn_nu),...
                        mu_Sun); %Saturn Position, Helio-Ecliptic [x y z] (km)
r_Saturn_M = r_Saturn_M*1e3; %Saturn Position, Helio-Ecliptic [x y z] (m)   

%Planetary Ephemerides at Uranus:
%Compute Uranus position:
[r_Uranus_M,~] = COE2RV(ephem.uranus_a,ephem.uranus_e,rad2deg(ephem.uranus_i),...
                        rad2deg(ephem.uranus_w),rad2deg(ephem.uranus_Omega),rad2deg(ephem.uranus_nu),...
                        mu_Sun); %Uranus Position, Helio-Ecliptic [x y z] (km)
r_Uranus_M = r_Uranus_M*1e3; %Uranus Position, Helio-Ecliptic [x y z] (m) 

%Planetary Ephemerides at Neptune:
%Compute Neptune position:
[r_Neptune_M,~] = COE2RV(ephem.neptune_a,ephem.neptune_e,rad2deg(ephem.neptune_i),...
                        rad2deg(ephem.neptune_w),rad2deg(ephem.neptune_Omega),rad2deg(ephem.neptune_nu),...
                        mu_Sun); %Neptune Position, Helio-Ecliptic [x y z] (km)
r_Neptune_M = r_Neptune_M*1e3; %Neptune Position, Helio-Ecliptic [x y z] (m) 
 
 %% Generate DE405 Ephemeris:                   

%Add DE405 Required Chebychev Polynomial Coefficients:
global PC
load DE405Coeff.mat
PC = DE405Coeff;

%Generate JPL DE405 Ephemeris:
Mjd_UTC = Mjday(yr, 12, 11, 12, 0, 0); %Modified Julian Date
[r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
          r_Neptune,r_Pluto,r_Moon,r_Sun] = JPL_Eph_DE405(Mjd_UTC);
      
%% Comparison:

phi = 23.4;
R1 = [1 0 0; 0 cosd(phi) sind(phi); 0 -sind(phi) cosd(phi)];
R2 = [cosd(phi) 0 -sind(phi); 0 1 0; sind(phi) 0 cosd(phi)];
R3 = [cosd(phi) sind(phi) 0; -sind(phi) cosd(phi) 0; 0 0 1];

%Mercury:
err_Mercury(count) = norm(abs((r_Mercury_M-r_Mercury)./r_Mercury)*1e2);
%Venus:
err_Venus(count) = norm(abs((r_Venus_M-r_Venus)./r_Venus)*1e2);
%Earth:
err_Earth(count) = norm(abs((r_Earth_M-r_Earth)./r_Earth)*1e2);
%Mars:
err_Mars(count) = norm(abs((r_Mars_M-r_Mars)./r_Mars)*1e2);
%Jupiter:
err_Jupiter(count) = norm(abs((r_Jupiter_M-r_Jupiter)./r_Jupiter)*1e2);
%Saturn:
err_Saturn(count) = norm(abs((r_Saturn_M-r_Saturn)./r_Saturn)*1e2);
%Uranus:
err_Uranus(count) = norm(abs((r_Uranus_M-r_Uranus)./r_Uranus)*1e2);
%Neptune:
err_Neptune(count) = norm(abs((r_Neptune_M-r_Neptune)./r_Neptune)*1e2);

%Mercury:
pos_Mercury(count) = norm(abs(r_Mercury_M-r_Mercury))/1e3;
%Venus:
pos_Venus(count) = norm(abs(r_Venus_M-r_Venus))/1e3;
%Earth:
pos_Earth(count) = norm(abs(r_Earth_M-r_Earth))/1e3;
%Mars:
pos_Mars(count) = norm(abs(r_Mars_M-r_Mars))/1e3;
%Jupiter:
pos_Jupiter(count) = norm(abs(r_Jupiter_M-r_Jupiter))/1e3;
%Saturn:
pos_Saturn(count) = norm(abs(r_Saturn_M-r_Saturn))/1e3;
%Uranus:
pos_Uranus(count) = norm(abs(r_Uranus_M-r_Uranus))/1e3;
%Neptune:
pos_Neptune(count) = norm(abs(r_Neptune_M-r_Neptune))/1e3;

end
%% Plot:

%Upload Planet Info:
global moon mercury venus earth mars jupiter saturn uranus neptune;
planetary_database;

%{
%Plot error increase over time:
figure()
whitebg('k')
scatter(yr_vec,err_Mercury,'m')
hold on
fit = polyfit(yr_vec-2016,err_Mercury,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'m')
scatter(yr_vec,err_Venus,'g')
fit = polyfit(yr_vec-2016,err_Venus,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'g')
scatter(yr_vec,err_Earth,'b')
fit = polyfit(yr_vec-2016,err_Earth,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'b')
scatter(yr_vec,err_Mars,'r')
fit = polyfit(yr_vec-2016,err_Mars,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'r')
xlabel('Year')
ylabel('Percent Error')
title({'Percent Error Magnitude Between Meeus and DE405 Models', 'Inner Planets'})
legend('Mercury Sample','Mercury Fit','Venus Sample','Venus Fit','Earth Sample','Earth Fit','Mars Sample','Mars Fit')
grid on

figure()
whitebg('k')
scatter(yr_vec,err_Jupiter,'y')
hold on
fit = polyfit(yr_vec-2016,err_Jupiter,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'y')
or = [1.0,0.687,0.387];
scatter(yr_vec,err_Saturn,[],or)
fit = polyfit(yr_vec-2016,err_Saturn,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'Color',or)
pu = [0.5,0,0.5];
scatter(yr_vec,err_Uranus,[],pu)
fit = polyfit(yr_vec-2016,err_Uranus,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'Color',pu)
scatter(yr_vec,err_Neptune,'c')
fit = polyfit(yr_vec-2016,err_Neptune,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'c')
xlabel('Year')
ylabel('Percent Error')
title({'Percent Error Magnitude Between Meeus and DE405 Models','Outer Planets'})
legend('Jupiter Sample','Jupiter Fit','Saturn Sample','Saturn Fit','Uranus Sample','Uranus Fit','Neptune Sample','Neptune Fit')
grid on
%}

%Plot position error over time:
figure()
whitebg('k')
scatter(yr_vec,pos_Mercury,'m')
hold on
fit = polyfit(yr_vec-2016,pos_Mercury,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'m')
fit = polyfit(yr_vec-2016,pos_Mercury,1);
plot(yr_vec,polyval(fit,yr_vec-2016),'m--')
scatter(yr_vec,pos_Venus,'g')
fit = polyfit(yr_vec-2016,pos_Venus,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'g')
fit = polyfit(yr_vec-2016,pos_Venus,1);
plot(yr_vec,polyval(fit,yr_vec-2016),'g--')
scatter(yr_vec,pos_Earth,'b')
fit = polyfit(yr_vec-2016,pos_Earth,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'b')
fit = polyfit(yr_vec-2016,pos_Earth,1);
plot(yr_vec,polyval(fit,yr_vec-2016),'b--')
scatter(yr_vec,pos_Mars,'r')
fit = polyfit(yr_vec-2016,pos_Mars,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'r')
fit = polyfit(yr_vec-2016,pos_Mars,1);
plot(yr_vec,polyval(fit,yr_vec-2016),'r--')
xlabel('Year')
ylabel('Position Magnitude Error (km)')
title({'Position Error Magnitude Between Meeus and DE405 Models', 'Inner Planets'})
legend('Mercury Sample','Mercury Fit','Mars Linear Trend','Venus Sample','Venus Fit','Venus Linear Trend','Earth Sample','Earth Fit','Earth Linear Trend','Mars Sample','Mars Fit','Linear Trend')
grid on

figure()
whitebg('k')
scatter(yr_vec,pos_Jupiter,'y')
hold on
fit = polyfit(yr_vec-2016,pos_Jupiter,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'y')
or = [1.0,0.687,0.387];
scatter(yr_vec,pos_Saturn,[],or)
fit = polyfit(yr_vec-2016,pos_Saturn,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'Color',or)
pu = [0.5,0,0.5];
scatter(yr_vec,pos_Uranus,[],pu)
fit = polyfit(yr_vec-2016,pos_Uranus,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'Color',pu)
scatter(yr_vec,pos_Neptune,'c')
fit = polyfit(yr_vec-2016,pos_Neptune,5);
plot(yr_vec,polyval(fit,yr_vec-2016),'c')
xlabel('Year')
ylabel('Position Error Magnitude (km)')
title({'Position Error Magnitude Between Meeus and DE405 Models','Outer Planets'})
legend('Jupiter Sample','Jupiter Fit','Saturn Sample','Saturn Fit','Uranus Sample','Uranus Fit','Neptune Sample','Neptune Fit')
grid on

%{
%Plot Mercury:
whitebg('k')
grid on
h = plot_sphere(7e6,[0 0 0],1,1,'y' )
hold on
planet_ellipse(mercury,1,'m','-' )
plot_sphere( 3e6,r_Mercury_M/1e3,1,1,'m' )
plot_sphere( 3e6,r_Mercury/1e3,1,1,'r' )
xlabel('X-axis (km)','Color','k')
ylabel('Y-axis (km)','Color','k')
zlabel('Z-axis (km)','Color','k')
title('Mercury Position Comparison, December 11, 2016','Color','k')
legend('Sun','Mercury Orbit','Meeus','JPL DE405')
axis equal

figure
%Plot Venus:
whitebg('k')
grid on
h = plot_sphere(7e6,[0 0 0],1,1,'y' )
hold on
planet_ellipse(venus,1,'g','-' )
plot_sphere( 3e6,r_Venus_M/1e3,1,1,'g' )
plot_sphere( 3e6,r_Venus/1e3,1,1,'b' )
xlabel('X-axis (km)','Color','k')
ylabel('Y-axis (km)','Color','k')
zlabel('Z-axis (km)','Color','k')
title('Venus Position Comparison, December 11, 2016','Color','k')
legend('Sun','Venus Orbit','Meeus','JPL DE405')
axis equal

figure
%Plot Earth:
whitebg('k')
grid on
h = plot_sphere(7e6,[0 0 0],1,1,'y' )
hold on
planet_ellipse(earth,1,'b','-' )
plot_sphere( 3e6,r_Earth_M/1e3,1,1,'b' )
plot_sphere( 3e6,r_Earth/1e3,1,1,'w' )
xlabel('X-axis (km)','Color','k')
ylabel('Y-axis (km)','Color','k')
zlabel('Z-axis (km)','Color','k')
title('Earth Position Comparison, December 11, 2016','Color','k')
legend('Sun','Earth Orbit','Meeus','JPL DE405')
axis equal

figure
%Plot Mars:
whitebg('k')
grid on
h = plot_sphere(14e6,[0 0 0],1,1,'y' )
hold on
planet_ellipse(mars,1,'r','-' )
plot_sphere( 6e6,r_Mars_M/1e3,1,1,'r' )
plot_sphere( 6e6,r_Mars/1e3,1,1,'c' )
xlabel('X-axis (km)','Color','k')
ylabel('Y-axis (km)','Color','k')
zlabel('Z-axis (km)','Color','k')
title('Mars Position Comparison, December 11, 2016','Color','k')
legend('Sun','Mars Orbit','Meeus','JPL DE405')
axis equal


figure
%Plot Jupiter:
whitebg('k')
grid on
h = plot_sphere(7e7,[0 0 0],1,1,'y' )
hold on
planet_ellipse(jupiter,1,'r','-' )
plot_sphere( 3e7,r_Jupiter_M/1e3,1,1,'r' )
plot_sphere( 3e7,r_Jupiter/1e3,1,1,'y' )
xlabel('X-axis (km)','Color','k')
ylabel('Y-axis (km)','Color','k')
zlabel('Z-axis (km)','Color','k')
title('Jupiter Position Comparison, December 11, 2016','Color','k')
legend('Sun','Jupiter Orbit','Meeus','JPL DE405')
axis equal

figure
%Plot Saturn:
whitebg('k')
grid on
h = plot_sphere(9e7,[0 0 0],1,1,'y' )
hold on
planet_ellipse(saturn,1,'m','-' )
plot_sphere( 3e7,r_Saturn_M/1e3,1,1,'m' )
plot_sphere( 3e7,r_Saturn/1e3,1,1,'g' )
xlabel('X-axis (km)','Color','k')
ylabel('Y-axis (km)','Color','k')
zlabel('Z-axis (km)','Color','k')
title('Saturn Position Comparison, December 11, 2016','Color','k')
legend('Sun','Saturn Orbit','Meeus','JPL DE405')
axis equal

figure
%Plot Uranus:
whitebg('k')
grid on
h = plot_sphere(18e7,[0 0 0],1,1,'y' )
hold on
planet_ellipse(uranus,1,'b','-' )
plot_sphere( 6e7,r_Uranus_M/1e3,1,1,'b' )
plot_sphere( 6e7,r_Uranus/1e3,1,1,'w' )
xlabel('X-axis (km)','Color','k')
ylabel('Y-axis (km)','Color','k')
zlabel('Z-axis (km)','Color','k')
title('Uranus Position Comparison, December 11, 2016','Color','k')
legend('Sun','Uranus Orbit','Meeus','JPL DE405')
axis equal

figure
%Plot Neptune:
whitebg('k')
grid on
h = plot_sphere(54e7,[0 0 0],1,1,'y' )
hold on
planet_ellipse(neptune,1,'c','-' )
plot_sphere( 18e7,r_Neptune_M/1e3,1,1,'c' )
plot_sphere( 18e7,r_Neptune/1e3,1,1,'b' )
xlabel('X-axis (km)','Color','k')
ylabel('Y-axis (km)','Color','k')
zlabel('Z-axis (km)','Color','k')
title('Neptune Position Comparison, December 11, 2016','Color','k')
legend('Sun','Neptune Orbit','Meeus','JPL DE405')
axis equal
%}
