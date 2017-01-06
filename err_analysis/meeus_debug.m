%Name: meeus_debug
%Date: 12/5/16
%Purpose: Compare Kayla's Meeus to known Solution

%Initialize Console:
clear; close all; clc;

%JPL Astronomical Constants Required Include File:
mu_Sun = 1.3271244018e11; %(km^3/s^2)

%Select Comparison Epoch:
[ JD ] = get_JD( 12,12,11,2016 ); %Julian Date at 12:00 UT, December 11, 2016 

%% Generate Meeus Ephemeris Solution, Michael:
ephem = meeus_ephem(JD);
                    
 %% Generate Meeus Ephemeris, Kayla:
 
 %File location:
 addpath('./Kayla');
 
 %Planetary Ephem:
 [ephem2] = meeus( JD ) 
 
%% Comparison:
%Venus:
dMe_a = ephem.mercury_a - cell2mat(ephem2(2,2))
dMe_e = ephem.mercury_e - cell2mat(ephem2(3,2))
dMe_i = rad2deg(ephem.mercury_i) - cell2mat(ephem2(4,2))
dMe_Omega = rad2deg(ephem.mercury_Omega) - cell2mat(ephem2(5,2))
dMe_LOP = ephem.mercury_LOP - cell2mat(ephem2(6,2))
dMe_L = ephem.mercury_L - mod(cell2mat(ephem2(7,2)),360)
dMe_M = ephem.mercury_M - mod(cell2mat(ephem2(8,2)),360)
dMe_w = rad2deg(ephem.mercury_w) - mod(cell2mat(ephem2(9,2)),360)


%Venus:
dV_a = ephem.venus_a - cell2mat(ephem2(2,3))
dV_e = ephem.venus_e - cell2mat(ephem2(3,3))
dV_i = rad2deg(ephem.venus_i) - cell2mat(ephem2(4,3))
dV_Omega = rad2deg(ephem.venus_Omega) - cell2mat(ephem2(5,3))
dV_LOP = ephem.venus_LOP - cell2mat(ephem2(6,3))
dV_L = ephem.venus_L - mod(cell2mat(ephem2(7,3)),360)
dV_M = ephem.venus_M - mod(cell2mat(ephem2(8,3)),360)
dV_w = rad2deg(ephem.venus_w) - mod(cell2mat(ephem2(9,3)),360)

%Earth:
dE_a = ephem.earth_a - cell2mat(ephem2(2,4))
dE_e = ephem.earth_e - cell2mat(ephem2(3,4))
dE_i = rad2deg(ephem.earth_i) - cell2mat(ephem2(4,4))
dE_Omega = rad2deg(ephem.earth_Omega) - cell2mat(ephem2(5,4))
dE_LOP = ephem.earth_LOP - cell2mat(ephem2(6,4))
dE_L = ephem.earth_L - mod(cell2mat(ephem2(7,4)),360)
dE_M = ephem.earth_M - mod(cell2mat(ephem2(8,4)),360)
dE_w = rad2deg(ephem.earth_w) - mod(cell2mat(ephem2(9,4)),360)

%Mars:
dE_a = ephem.mars_a - cell2mat(ephem2(2,5))
dE_e = ephem.mars_e - cell2mat(ephem2(3,5))
dE_i = rad2deg(ephem.mars_i) - cell2mat(ephem2(4,5))
dE_Omega = rad2deg(ephem.mars_Omega) - cell2mat(ephem2(5,5))
dE_LOP = ephem.mars_LOP - cell2mat(ephem2(6,5))
dE_L = ephem.mars_L - mod(cell2mat(ephem2(7,5)),360)
dE_M = ephem.mars_M - mod(cell2mat(ephem2(8,5)),360)
dE_w = rad2deg(ephem.mars_w) - mod(cell2mat(ephem2(9,5)),360)

%Jupiter:
dE_a = ephem.jupiter_a - cell2mat(ephem2(2,6))
dE_e = ephem.jupiter_e - cell2mat(ephem2(3,6))
dE_i = rad2deg(ephem.jupiter_i) - cell2mat(ephem2(4,6))
dE_Omega = rad2deg(ephem.jupiter_Omega) - cell2mat(ephem2(5,6))
dE_LOP = ephem.jupiter_LOP - cell2mat(ephem2(6,6))
dE_L = ephem.jupiter_L - mod(cell2mat(ephem2(7,6)),360)
dE_M = ephem.jupiter_M - mod(cell2mat(ephem2(8,6)),360)
dE_w = rad2deg(ephem.jupiter_w) - mod(cell2mat(ephem2(9,6)),360)

%Saturn:
dE_a = ephem.saturn_a - cell2mat(ephem2(2,7))
dE_e = ephem.saturn_e - cell2mat(ephem2(3,7))
dE_i = rad2deg(ephem.saturn_i) - cell2mat(ephem2(4,7))
dE_Omega = rad2deg(ephem.saturn_Omega) - cell2mat(ephem2(5,7))
dE_LOP = ephem.saturn_LOP - cell2mat(ephem2(6,7))
dE_L = ephem.saturn_L - mod(cell2mat(ephem2(7,7)),360)
dE_M = ephem.saturn_M - mod(cell2mat(ephem2(8,7)),360)
dE_w = rad2deg(ephem.saturn_w) - mod(cell2mat(ephem2(9,7)),360)


